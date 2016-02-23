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

#include "population.h"
#include "stat.h"

namespace matvec {

unsigned sampling(double* dist_prob,const unsigned size)
{
  //////////////////////////////////////////////////////////////////////
  // sampling a genotype id from the distribution given in dist_prob
  // note that dist_prob is destroyed in this routine
  ///////////////////////////////////////////////////////////////////////

   unsigned retval,i;
   for (i=1; i<size; i++) dist_prob[i] += dist_prob[i-1];
   double u = ranf();
   for (retval=size, i=0; i<size; i++) if (u<=dist_prob[i]) retval=i;
   return retval;
}

void G_to_gamete(const unsigned G_id,unsigned k1,unsigned k2,unsigned size)
 // convert genotype id into two gamete ids. size is the size of gamete pool
 // For example, size = 3
 //  k1    0 0 0  1 1 1  2 2 2
 //  k2    0 1 2  0 1 2  0 1 2
 // ---------------------------
 // G_id   0 1 2  3 4 5  6 7 8
{
   k1 = G_id/size;
   k2 = G_id - k1*size;
}

void Population::count_genotype(Individual* I)
{
   unsigned i,j,k,t;
   for (t=0; t<numchrom; t++) {
      i = I->paternal_chrom_id(t);
      j = I->maternal_chrom_id(t);
      if (i >= j) k = i*(i+1)/2+j;
      else k = j*(j+1)/2 + i;
      I->genotype_counter[t][k] += 1.0;
   }
}

unsigned Population::ind_gamete(const Individual* I, const unsigned jc,
                                      unsigned g_id[], double fq[])
{
   unsigned i,k,retval,ii,jj;
   if (I) {
      ii = I->paternal_chrom_id(jc);
      jj = I->maternal_chrom_id(jc);
      if (ii == jj) {
         retval = 1;
         g_id[0] = ii;
         fq[0] = 1.0;
      }
      else {
         if (ii > jj)  k = ii*(ii-1)/2 + jj;
         else  k = jj*(jj-1)/2 + ii;
         Chromosome* C = &(gametebase[jc].chromosome[k]);
         retval = C->nloci();
         for (i=0; i<retval; i++) {
            g_id[i] = C->locus[i].allele;
            fq[i] = C->locus[i].effect;
         }
      }
   }
   else {
      retval = pop_gamete[jc].nchrom();
      for (i=0; i<retval; i++) {
         g_id[i] = i;  fq[i] = pop_gamete[jc].chromosome[i].freq();
      }
   }
   return retval;
}

void Population::build_pop_gamete()
{
   /////////////////////////////////////////////////////////////
   //  It is assumed that population is in equilibrium status
   /////////////////////////////////////////////////////////////

   if (pop_gamete) return;
   if (!stdid) renum();
   if (pop_gamete) {
     delete [] pop_gamete;
     pop_gamete=0;
   }
   pop_gamete = new_Genome_vec(numchrom);
   unsigned chr,nl,nc,num,allele_id;
   unsigned i,j,k,h;
   double gamete_freq;
   ChromStruct *Chrom = prior->chrom();
   Genome *G, *GB;

   for (tn_gamete=1,chr = 0; chr< numchrom; chr++) {
      nl = Chrom[chr].nloci();
      Vector<unsigned> ngtvec(nl);
      Vector<unsigned> gtvec(nl);
      for (num=1,k=0; k<nl; k++) {
         num *= ngtvec[k] = Chrom[chr].locus[k].nallele();
      }
      tn_gamete *= num;
      G = &(pop_gamete[chr]);
      G->resize(num,nl);
      for (i=0; i<num; i++) {
         gtindex(i,nl,ngtvec,gtvec);
         for (gamete_freq=1.0, k=0; k<nl; k++) {
            allele_id = gtvec[k];
            G->chromosome[i].locus[k].allele = allele_id;
            gamete_freq *= Chrom[chr].locus[k].allele_freq[allele_id];
         }
         G->chromosome[i].id(i);
         G->chromosome[i].freq(gamete_freq);
      }
   }

   /////////////////////  build gamete base information ////////////////////
   // gamete id 0, 1, 2,.... num-1. There are num*(num-1)/2 combinations:
   // (1,0)
   // (2,0)  (2,1)
   // (3,0)  (3,1)  (3,2)
   //    . . .
   ///////////////////////////////////////////////////////////////////////
   if(gametebase){
     if (gametebase) delete [] gametebase; 
     gametebase=0;
   }
   gametebase = new_Genome_vec(numchrom);
   Chromosome gamete0, *gamete1, *gamete2;
   double r_k2_kk, freq;
   unsigned ngt,k2,kk,gamete_id;
   for (chr = 0; chr<numchrom; chr++) {
      nl = Chrom[chr].nloci();
      gamete0.resize(nl);
      Vector<unsigned> ngtvec(nl);
      Vector<unsigned> gtvec(nl);
      for (num=1,i=0; i<nl; i++) num *= Chrom[chr].locus[i].nallele();
      GB = &(gametebase[chr]);
      GB->resize(num*(num-1)/2, 0);
      G = &(pop_gamete[chr]);

      for (nc =0, i=1; i<num; i++) {
         gamete1 = &(G->chromosome[i]);
         for (j=0; j<i; j++) {
            gamete2 = &(G->chromosome[j]);
            for (ngt=1, k=0; k<nl; k++) {
               if (gamete1->locus[k].allele == gamete2->locus[k].allele)
                  ngtvec[k] = 1;
               else {
                  ngtvec[k] = 2;
                  ngt *= 2;
               }
            }
            GB->chromosome[nc].resize(ngt);
            for (k=0; k<ngt; k++) {
               gtindex(k,nl,ngtvec,gtvec);
               freq = 0.5;
               k2 = 0;     // the first 2-different allele case hasn't found
               for (kk=0; kk<nl; kk++) {
                  if (ngtvec[kk] == 2) {
                     if (k2 > 0 ) {  // not the first 2-different allele case
                        r_k2_kk = prior->recomb_rate(chr+1,k2,kk+1);
                        if (gtvec[kk] == gtvec[k2-1]) freq *= 1.0 -  r_k2_kk;
                        else freq *= r_k2_kk;
                     }
                     k2 = kk+1;
                  }
               }
               ///////////////////////////////////////////////////////////////
               // GB->chromosome[nc].locus[k].allele = gamete_id in pop_gamete
               // GB->chromosome[nc].locus[k].effect = gamete freq
               ///////////////////////////////////////////////////////////////
               GB->chromosome[nc].locus[k].effect = freq;
               for (kk=0; kk<nl; kk++) {
                  if (gtvec[kk] == 0) {
                    gamete0.locus[kk].allele = gamete1->locus[kk].allele;
                  }
                  else gamete0.locus[kk].allele = gamete2->locus[kk].allele;
               }
               gamete_id = G->chrom_id(gamete0);
               GB->chromosome[nc].locus[k].allele = gamete_id;
            }
            nc++;                            // nc = i*(i-1)/2 + j
         }
      }
   }
   tn_genotype = tn_gamete*(tn_gamete+1)/2;

   //////////////////////////////////////////////////////////
   //  build a vector of genotypic value: mean_for_genotype
   //////////////////////////////////////////////////////////
   mean_for_genotype.resize(tn_genotype);
   unsigned nchrom = prior->nchrom();
  std::cout << "pop nchrom " << nchrom << std::endl;
   Vector<unsigned> ngtvec(nchrom);
   Vector<unsigned> sire_gtvec(nchrom);
   Vector<unsigned> dam_gtvec(nchrom);
   for (i=0; i<nchrom; i++) ngtvec[i] = prior->chrom()[i].locus[0].nallele();

   for (kk=0,i=0; i<tn_gamete; i++) {
      gtindex(i,nchrom,ngtvec,sire_gtvec);
      for (j=0; j<=i; j++) {
         gtindex(j,nchrom,ngtvec,dam_gtvec);
         for (freq=0.0,k=0; k<nchrom; k++) {
             for (h=0; h<nl; h++) {     // nl = 1
                freq +=  prior->genotypic_val(k,h)[sire_gtvec[k]][dam_gtvec[k]];
             }
         }
         mean_for_genotype[kk++] = freq;
      }
   }
}

double  Population::cprob_children(const Individual* I, const unsigned jc,
                       unsigned sg[],unsigned dg[],double sgf[],double dgf[])
{
 ////////////////////////////////////////////////////
 //   note that if no children, return 1.0
 //////////////////////////////////////////////////

   unsigned i,j,t,k,nc,ns,gid,ssize,dsize,noffs;
   Individual *child, **offspring;
   ns = I->nspouse();
   double prob,tmp1,tmp2;
   noffs = I->noffs();
   offspring = I->offspring();

   prob=1.0;
   if (I->sex()== 'F') {  // I is a mother
      dsize = ind_gamete(I, jc,dg,dgf);
      for (k=0; k<noffs; k++) {
         child = offspring[k];
         gid = child->maternal_chrom_id(jc);
         for (tmp1=0.0, t=0; t < dsize; t++) {
            if (gid == dg[t]) {
               tmp1 = dgf[t];
               break;
            }
         }
         if  (tmp1 == 0.0) return tmp1;
         prob *= tmp1;
      }
      for (k=0, i=0; i<ns; i++) {
         ssize = ind_gamete(I->myoffspring[k]->father(), jc, sg,sgf);
         nc = I->numoffs_spouse[i];
         for (j=0; j<nc; j++) {
            child = offspring[k++];
            gid = child->paternal_chrom_id(jc);
            for (tmp2=0.0, t=0; t < ssize; t++) {
               if (gid == sg[t]) {
                  tmp2 = sgf[t];
                  break;
               }
            }
            if  (tmp2 == 0.0) return tmp2;
            prob *= tmp2;
         }
      }
   }
   else {   // I is a father
      ssize = ind_gamete(I, jc,sg,sgf);
      for (k=0; k<noffs; k++) {
         child = offspring[k];
         gid = child->paternal_chrom_id(jc);
         for (tmp1=0.0, t=0; t < ssize; t++) {
            if (gid == sg[t]) {
               tmp1 = sgf[t];
               break;
            }
         }
         if  (tmp1 == 0.0) return tmp1;
         prob *= tmp1;
      }

      for (k=0, i=0; i<ns; i++) {
         dsize = ind_gamete(I->myoffspring[k]->mother(), jc, dg,dgf);
         nc = I->numoffs_spouse[i];
         for (j=0; j<nc; j++) {
            child = offspring[k++];
            gid = child->maternal_chrom_id(jc);
            for (tmp2=0.0, t=0; t < dsize; t++) {
               if (gid == dg[t]) {
                  tmp2 = dgf[t];
                  break;
               }
            }
            if  (tmp2 == 0.0) return tmp2;
            prob *= tmp2;
         }
      }
   }
   return prob;
}

unsigned Population::partial_cdist(Individual *I,const unsigned jc,
                                   unsigned** cdist_value, double* cdist_prob,
                                   unsigned** gid_mat, double** freq_mat)
{
   unsigned *s_gamete_id = gid_mat[0];
   unsigned *d_gamete_id = gid_mat[1];
   double   *s_gamete_freq = freq_mat[0];
   double   *d_gamete_freq = freq_mat[1];
   unsigned ss = ind_gamete(I->father(), jc, s_gamete_id, s_gamete_freq);
   unsigned ds = ind_gamete(I->mother(), jc, d_gamete_id, d_gamete_freq);

   unsigned i,j,k,sg,dg;
   double prob;
   for (k=0,i=0; i<ss; i++) {
      sg = s_gamete_id[i];
      I->genome0.chromosome[jc] = pop_gamete[jc].chromosome[sg];
      for (j=0; j<ds; j++) {
         dg = d_gamete_id[j];
         I->genome1.chromosome[jc] = pop_gamete[jc].chromosome[dg];
         prob = (*penetrance_f)(I,(const double**)residual_var.begin());
         prob *= s_gamete_freq[i]*d_gamete_freq[j];
         if (prob > 1.0e-10) {
            cdist_value[0][k] = sg;
            cdist_value[1][k] = dg;
            cdist_prob[k++] = prob;
         }
      }
   }
   if (k == 0) throw exception(" Population::partial_cdist(): you have probably found a bug!");
   for (prob=0.0,i=0; i<k; i++) prob += cdist_prob[i];
   for (i=0; i<k; i++) cdist_prob[i] /= prob;
   return k;
}

unsigned Population::full_cdist(Individual *I,const unsigned jc,
                                unsigned** cdist_value, double* cdist_prob,
                                unsigned** gid_mat,     double** freq_mat)
{
   ////////////////////////////////////////////////////////////////////
   //   working space gid_mat,freq_mat must have at least four rows
   ////////////////////////////////////////////////////////////////////

   unsigned *s_gamete_id = gid_mat[0];
   unsigned *d_gamete_id = gid_mat[1];
   double   *s_gamete_freq = freq_mat[0];
   double   *d_gamete_freq = freq_mat[1];
   unsigned ss = ind_gamete(I->father(), jc, s_gamete_id, s_gamete_freq);
   unsigned ds = ind_gamete(I->mother(), jc, d_gamete_id, d_gamete_freq);

   unsigned i,j,k,sg,dg;
   double prob;
   for (k=0,i=0; i<ss; i++) for (j=0; j<ds; j++) {
      sg = s_gamete_id[i];
      dg = d_gamete_id[j];
      I->genome0.chromosome[jc] = pop_gamete[jc].chromosome[sg];
      I->genome1.chromosome[jc] = pop_gamete[jc].chromosome[dg];
      prob = 1.0;
      if (I->noffs() > 0) {
         prob = cprob_children(I,jc,gid_mat[2],gid_mat[3],freq_mat[2],
                               freq_mat[3]);
      }
      if (prob != 0.0) {
         prob *= (*penetrance_f)(I,(const double**)residual_var.begin());
         prob *= s_gamete_freq[i]*d_gamete_freq[j];
         if (prob > 1.0e-10) {
            cdist_value[0][k] = sg;
            cdist_value[1][k] = dg;
            cdist_prob[k++]   = prob;
         }
      }
   }
   if (k == 0) throw exception(" Population.full_cdist(...): you have probably found a bug!");
   for (prob=0.0,i=0; i<k; i++) prob += cdist_prob[i];
   for (i=0; i<k; i++) cdist_prob[i] /= prob;
   return k;
}

double Population::llhood_phenotype(const unsigned nrep)
{
   if (!pop_gamete) build_pop_gamete();

   Individual *I;
   unsigned i,j,it,s,d,nsg,ndg;
   Chromosome *C0, *C1;

   unsigned maxg = pop_gamete[0].nchrom(); // maxg = maximum number of gametes
   for (i=1; i<numchrom; i++) {           // among chromosomes
      j = pop_gamete[i].nchrom();
      if (j > maxg) maxg = j;
   }

   /////////////////////////////////////////////////////////////
   //   get working space  to hold parents gametes information
   /////////////////////////////////////////////////////////////

   Vector<unsigned> s_gamete_id(maxg);
   Vector<unsigned> d_gamete_id(maxg);
   Vector<double> s_gamete_freq(maxg);
   Vector<double> d_gamete_freq(maxg);

   Vector<double> work(nrep);
   double prob,scaled_coef;
   for (it=0; it<nrep; it++) {
      for (prob=0.0, i=0; i<popsize; i++) {
         I = popmember[i];
         C0 = I->paternal_chrom();
         C1 = I->maternal_chrom();
         for (j=0; j<numchrom; j++) {
            nsg = ind_gamete(I->father(), j, s_gamete_id.begin(),s_gamete_freq.begin());
            ndg = ind_gamete(I->mother(), j, d_gamete_id.begin(),d_gamete_freq.begin());
            s = sampling(s_gamete_freq.begin(), nsg);
            d = sampling(d_gamete_freq.begin(), ndg);
            C0[j] = pop_gamete[j].chromosome[s_gamete_id[s]];
            C1[j] = pop_gamete[j].chromosome[d_gamete_id[d]];
         }
         prob += std::log((*penetrance_f)(I,(const double**)residual_var.begin()));
      }
      work[it] = prob;
   }
   scaled_coef = work.sum()/work.size();
   for (it=0; it<nrep; it++) work[it] = std::exp(work[it]-scaled_coef);
   return  std::log(work.sum()/work.size()) + scaled_coef;
}

void Population::genotype_config(const char type[])

{
  ////////////////////////////////////////////////////////////////
  //   a possible genotype configuration G over a given pedigree,
  /////////////////////////////////////////////////////////////////

   if (!pop_gamete) build_pop_gamete();

   Individual *I;
   unsigned i,j, s,d,nsg,ndg;
   Chromosome *C0, *C1;

   unsigned maxg = pop_gamete[0].nchrom(); // ng = maxximum number of gametes
   for (i=1; i<numchrom; i++) {           // among chromosomes
      j = pop_gamete[i].nchrom();
      if (j > maxg) maxg = j;
   }
   /////////////////////////////////////////////////////////////
   //   get working space  to hold parents gametes information
   /////////////////////////////////////////////////////////////

   Vector<unsigned> s_gamete_id(maxg);
   Vector<unsigned> d_gamete_id(maxg);
   Vector<double> s_gamete_freq(maxg);
   Vector<double> d_gamete_freq(maxg);

   if (strcmp(type,"qtl") == 0) {

      for (i=0; i<popsize; i++) {
         I = popmember[i];
         C0 = I->paternal_chrom();
         C1 = I->maternal_chrom();
         for (j=0; j<numchrom; j++) {
            nsg = ind_gamete(I->father(), j, s_gamete_id.begin(),s_gamete_freq.begin());
            ndg = ind_gamete(I->mother(), j, d_gamete_id.begin(),d_gamete_freq.begin());
            s = sampling(s_gamete_freq.begin(), nsg);
            d = sampling(d_gamete_freq.begin(), ndg);
            C0[j] = pop_gamete[j].chromosome[s_gamete_id[s]];
            C1[j] = pop_gamete[j].chromosome[d_gamete_id[d]];
         }
      }
   }
   else if (strcmp(type,"ml") == 0) {
      ;
   }
   else if (strcmp(type,"ml_qtl") == 0) {
     ;
   }
   else {
     ;
   }
}

void Population::cond_genotype_config()
{
  ////////////////////////////////////////////////////////////////////
  //   a conditional genotype configuration G over a given pedigree,
  //  given y.
  ///////////////////////////////////////////////////////////////////

   if (!pop_gamete) build_pop_gamete();

   Individual *I;
   unsigned i,j,k,k0,k1,size;
   Chromosome *C0, *C1;

   unsigned maxg = pop_gamete[0].nchrom(); // ng = maxximum number of gametes
   for (i=1; i<numchrom; i++) {            // among chromosomes
      k = pop_gamete[i].nchrom();
      if (k > maxg) maxg = k;
   }
   //////////////////////////////////////////////////////
   //   get working space ready for partial_cdist()
   //   freq_mat, gid_mat, cdist_value, cdist_prob
   ////////////////////////////////////////////////////

   Matrix<double> freq_mat(2,maxg);
   Matrix<unsigned> gid_mat(2,maxg);
   Matrix<unsigned> cdist_value(2,maxg*maxg);
   Vector<double> cdist_prob(maxg*maxg);

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      C0 = I->paternal_chrom();
      C1 = I->maternal_chrom();

      ///////////////////////////////////////////////////////////////
      // first randomly assigned genotypes for each chromosome except
      // the first chromosome.
      // there is other alternatives such as using calculations of
      // f(y|Q1), f(y|Q1,Q2) ...,
      ///////////////////////////////////////////////////////////////
      for (j=1; j<numchrom; j++) {
         maxg = pop_gamete[j].nchrom();
         for (k=0; k<maxg; k++) {
            freq_mat[0][k] = pop_gamete[j].chromosome[k].freq();
         }
         k0 = sampling(freq_mat[0], maxg);
         k1 = sampling(freq_mat[0], maxg);
         C0[j] = pop_gamete[j].chromosome[k0];
         C1[j] = pop_gamete[j].chromosome[k1];
      }
      for (j=0; j<numchrom; j++) {
         size = partial_cdist(I,j,cdist_value.begin(),cdist_prob.begin(),gid_mat.begin(),freq_mat.begin());
         k = sampling(cdist_prob.begin(),size);
         k0 = cdist_value[0][k];               // paternal gamete id
         k1 = cdist_value[1][k];                // maternal gamete id
         C0[j] = pop_gamete[j].chromosome[k0];
         C1[j] = pop_gamete[j].chromosome[k1];
      }
   }
}

void Population::gibbs_iterate(const int count_gt)
{
   Individual *I;
   Chromosome *C0, *C1;
   unsigned i,j,k,k0,k1,size;

   unsigned maxg = pop_gamete[0].nchrom();  // maxg = maximum number of gametes
   for (i=1; i<numchrom; i++) {             // among chromosomes
      k = pop_gamete[i].nchrom();
      if (k > maxg) maxg = k;
   }

   //////////////////////////////////////////////////////
   //   get working space ready for full_cdist()
   //   freq_mat, gid_mat must have at least four rows
   //  cdist_value, cdist_prob
   ////////////////////////////////////////////////////

   Matrix<double> freq_mat(4,maxg);
   Matrix<unsigned> gid_mat(4,maxg);
   Matrix<unsigned> cdist_value(2,maxg*maxg);
   Vector<double> cdist_prob(maxg*maxg);

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      C0 = I->paternal_chrom();
      C1 = I->maternal_chrom();
      for (j=0; j<numchrom; j++) {
         size = full_cdist(I,j,cdist_value.begin(),cdist_prob.begin(),gid_mat.begin(),freq_mat.begin());
         k = sampling(cdist_prob.begin(),size);
         k0 = cdist_value[0][k];                    // paternal gamete id
         k1 = cdist_value[1][k];                    // maternal gamete id
         C0[j] = pop_gamete[j].chromosome[k0];
         C1[j] = pop_gamete[j].chromosome[k1];
      }
      if (count_gt) count_genotype(I);
   }
}

void Population::release_genotype_counter(void)
{
   Vector<double> *vv;
   for (unsigned i=0; i<popsize; i++) {
      vv = popmember[i]->genotype_counter;
      if (vv) {
         delete [] vv;
         popmember[i]->genotype_counter = 0;
      }
   }
}

void Population::resize_genotype_counter(void)
{
   if (!pop_gamete) build_pop_gamete();

   unsigned i,t,n;
   Individual *I;
   for (i=0; i<popsize; i++) {
      I = popmember[i];
      if(numchrom>0){
	I->genotype_counter = new Vector<double> [numchrom];
      }
      else {
	I->genotype_counter = 0;
      }
      for (t=0; t<numchrom; t++) {
         n = pop_gamete[t].size();
         I->genotype_counter[t].resize(n*(n+1)/2, 0.0);
      }
   }
}
} ///////// end of namespace matvec

