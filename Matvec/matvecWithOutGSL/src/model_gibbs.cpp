/*
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Library General Public
 License as published by the Free Software Foundation; either
 version 2 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Library General Public License for more details.
 
 You should have received a copy of the GNU Library General Public
 License along with this library; if not, write to the Free
 Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 MA 02111-1307, USA 
 */

#include <sstream>
#include <iomanip>

#include "model.h"
#include "stat.h"
#include "population.h"

namespace matvec {

void Model::sampleW(Vector<double>& bu,Vector<double>& newrhs)
{
   double sd_error = std::sqrt(residual_var[0][0]);
   Vector<double> xbzu_trait(numtrait);
   double vy, catW;

   unsigned i,ii,k,rec,t1;
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);
   if (!modelfile) throw exception(" Model::sampleW(): cannot open binary file");
   newrhs.assign(0.0);
   double *rhs = newrhs.begin() - 1;  // now rhs starts from 1
   double *sol = bu.begin() - 1;      // now sol starts from 1

   for (rec=0; rec<numrec; rec++) {
      input_pos_val(modelfile);
      k = pos_term[numterm];
      for (t1=0; t1<numtrait; t1++) {
         xbzu_trait[t1] = 0.0;
         for (i=0; i<numterm; i++) {
            ii = term[i].addr[t1];          // only one trait allowed now
            if (ii) xbzu_trait[t1] += sol[pos_term[i]]*xval_term[i];
         }
      }
      catW = snorm();
      if (trait_vec[0] == 0.0) while (catW > 0.0) catW = snorm();
      else while (catW < 0.0) catW = snorm();

      catW = sd_error*catW + xbzu_trait[0];
      vy = rve[k][0][0]*catW;
      for (i=0; i<numterm; i++) {
         ii = term[i].addr[0];          // only one trait allowed now
         if (ii) rhs[ii] += vy*xval_term[i];
      }
   }
   //BRS modelfile.close();
}

double Model::log_likelihood_gibbs(const unsigned wup,const unsigned glen)
{
   warning("Model.log_likelihood_gibbs(): not inplemented yet");
   return 0;
}

Vector<double>* Model::genotypic_value_cat(const unsigned warmup,
                                  const unsigned gibbslen)
{
   if (!data_prepared) prepare_data();

   blupsol.resize(hmmesize,0.0);
   Vector<double> tmpsol(hmmesize);
   Vector<double> wspace; wspace.reserve(hmmesize);
   setup_mme(&rellrhs);            // relrhs is used here

   /////////////////////////////////////////////////////
   //    Gibbs Sampler warm up
   /////////////////////////////////////////////
   unsigned it;
   for (it=0; it<warmup; it++) {
      sampleW(tmpsol,rellrhs);
      hmmec.gibbs_iterate(tmpsol.begin(),rellrhs.begin(),wspace.begin());
   }

   /////////////////////////////////////////////////////
   //    samples from Gibbs Sampler now will be used
   ////////////////////////////////////////////////////
   for (it=0; it<gibbslen; it++) {
      sampleW(tmpsol,rellrhs);
      hmmec.gibbs_iterate(tmpsol.begin(),rellrhs.begin(),wspace.begin());
      blupsol += tmpsol;
   }
   blupsol /= gibbslen;
   hmmec.resize(0,0);
   return &blupsol;
}

double Model::genotype_dist_gibbs(const unsigned warmup,const unsigned gibbslen)
{
   if (ntermGdist==0) {
      warning(" Model::genotype_dist_gibbs(): no GeneticDist term");
      return 0.0;
   }
   if (!data_prepared) prepare_data();
   unsigned i,it;
   blupsol.resize(hmmesize,0.0);
   Vector<double> wspace; wspace.reserve(hmmesize);
   pop->resize_genotype_counter();
   pop->residual_var = residual_var;
   pop->input_data(data);
   pop->cond_genotype_config();   //  an initial configuration

   if (numterm) setup_mme(&rellrhs);

   /////////////////////////////////////////////////////
   //    Gibbs Sampler warm up
   /////////////////////////////////////////////
   for (it=0; it<warmup; it++) {
      if (numterm) hmmec.gibbs_iterate(blupsol.begin(),rellrhs.begin(),wspace.begin());
      for (i=0; i<ntermGdist; i++) {
         compute_xbzu(blupsol);                   // to adjust trait
         pop->gibbs_iterate();
      }
   }

   /////////////////////////////////////////////////////
   //    samples from Gibbs Sampler now will be used
   ////////////////////////////////////////////////////

   for (it=0; it<gibbslen; it++) {
      if (numterm) hmmec.gibbs_iterate(blupsol.begin(),rellrhs.begin(),wspace.begin());
      for (i=0; i<ntermGdist; i++) {
         compute_xbzu(blupsol);                       // to adjust trait
         pop->gibbs_iterate(1);
      }
   }
   hmmec.resize(0,0);
   ntermGdist = 0;
   return 0.0;     // no useful value could be returned so far
}

Vector<double>* Model::genotypic_value_gibbs(const unsigned warmup,
                           const unsigned gibbslen)
{
   if (categorical) return genotypic_value_cat(warmup,gibbslen);
   if (!data_prepared) prepare_data();
   unsigned i,it;
   if (numterm) setup_mme(&rellrhs);
   blupsol.resize(hmmesize + ntermGdist*popsize, 0.0);
   Vector<double> tmpsol(hmmesize);
          tmpsol.assign(0.0);
   Vector<double> wspace; wspace.reserve(hmmesize);
   double * ve = &(blupsol[hmmesize]);
   if (ntermGdist) {
      pop->residual_var = residual_var;
      pop->input_data(data);
      pop->cond_genotype_config();
   }

   /////////////////////////////////////////////////////
   //    Gibbs Sampler warm up
   /////////////////////////////////////////////
   for (it=0; it<warmup; it++) {
      if (numterm) hmmec.gibbs_iterate(tmpsol.begin(),rellrhs.begin(),wspace.begin());
      for (i=0; i<ntermGdist; i++) {
         compute_xbzu(tmpsol);                   // to adjust trait
         pop->gibbs_iterate();
      }
   }

   /////////////////////////////////////////////////////
   //    samples from Gibbs Sampler now will be used
   ////////////////////////////////////////////////////

   for (it=0; it<gibbslen; it++) {
      if (numterm) {
         hmmec.gibbs_iterate(tmpsol.begin(),rellrhs.begin(),wspace.begin());
         for (i=0; i<hmmesize; i++) blupsol[i] += tmpsol[i];
      }
      for (i=0; i<ntermGdist; i++) {
         compute_xbzu(tmpsol);                       // to adjust trait
         pop->gibbs_iterate();
         for (i=0; i<popsize; i++) {
            ve[i] += pop->popmember[i]->genotypic_val();
         }
      }
   }
   blupsol /= gibbslen;
   hmmec.resize(0,0);
   return &(blupsol);
}

void Model::compute_xbzu(Vector<double>& bu)      // xbzu will be used to adjust trait
{
   unsigned i,ii,rec;
   double val;
   double *sol = bu.begin()-1;
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);
   if (!modelfile) throw exception("Model::compute_xbzu(): cannot open file");
   for (rec=0; rec<numrec; rec++) {
      input_pos_val(modelfile);
      val = 0.0;
      for (i=0; i<numterm; i++) {
         ii = term[i].addr[0];          // only one trait allowed now
         if (ii) val += sol[pos_term[i]]*xval_term[i];
      }
      ii = rec_indid[rec];
      if (ii > 0) pop->popmember[ii-1]->xbzu(val);
   }
   //BRS modelfile.close();
}
} //////// end of namespace matvec

