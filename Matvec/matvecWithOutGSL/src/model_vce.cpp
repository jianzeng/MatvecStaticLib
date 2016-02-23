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

#include "session.h"
#include "model.h"
#include "stat.h"
#include "population.h"
#include "individual.h"


namespace matvec {

extern void ABCt(double **A, const unsigned ma, const unsigned na,
          double **B, const unsigned mb, const unsigned nb,
          double **C, const unsigned mc, const unsigned nc,
          double **OUT,const unsigned mo, const unsigned no);  // OUT = A*B*C'

extern void getlambda(double **lambda, const int n);
extern double nr_powell(Data* D,Model* M,Vector<double> &p,Matrix<double> &xi,int nvc,
                        int& iter,double ftol);

unsigned Model::totalnvc(void) const
{
   ////////////////////////////////////////////////////
   //  total number of distinct variance components
   /////////////////////////////////////////////////////

   unsigned nvc,i;
   nvc = numtrait*(numtrait + 1)/2;
   for (i=0; i<numterm; i++) {
      if (term(i).classi() == 'R' || term(i).classi() =='P') {
	int nt=nt_vec[i];
         nvc += nt*(nt + 1)/2;
      }
   }
   return nvc;
}

void Model::var2vec(Vector<double>&x)
 /*******************************************************************
 * put variance components into a double vector
 * user is totally responsible to have an appropriate size for vector x
 * and matrix residual_var
 *********************************************************************/
{
   int i,k,t1,t2;
   doubleMatrix *var;
   for (k=0,i=0; i<numterm; i++) {
      if (term[i].classi() == 'R' || term[i].classi() == 'P') {
	int nt=nt_vec[i];
         var = term[i].prior->var_matrix();
         for (t1=0; t1<nt; t1++) for (t2=t1; t2<nt; t2++) {
             x[k++] = (*var)[t1][t2];
         }
      }
   }
   for (t1=0; t1<numtrait; t1++) {
      for (t2=t1; t2<numtrait; t2++) x[k++] = residual_var[t1][t2];
   }
}

void  Model::vec2var(const Vector<double>&x)
 ///////////////////////////////////////////////////////////////////////
 // put variance components into a double vector
 // user is totally responsible to have an appropriate size for vector x
 ///////////////////////////////////////////////////////////////////////
 {
   int i,k,t1,t2;
   doubleMatrix *var;
   for (k=0,i=0; i<numterm; i++) {
      if (term[i].classi() == 'R' || term[i].classi() == 'P') {
	int nt=nt_vec[i];
         var = term[i].prior->var_matrix();
         for (t1=0; t1<nt; t1++) {
            for (t2=t1; t2<nt; t2++,k++) {
               (*var)[t1][t2] = x[k];
               if ( t2>t1 ) (*var)[t2][t1] = x[k];
            }
         }
      }
   }
   for (t1=0; t1<numtrait; t1++) {
      for (t2=t1; t2<numtrait; t2++,k++) {
         residual_var[t1][t2] = x[k];
         if ( t2>t1 ) residual_var[t2][t1] = x[k];
      }
   }
}

double Model::restricted_log_likelihood(void)
{
  ////////////////////////////////////////////////////////////////////
  // it returns restricted log likelihood under this model with data D
  //  Important Notice:
  //  because  restricted_log_likelihood will be called sequentially
  //  under the same MME structure.
  //  hmmec is kept after restricted_log_likelihood(D) call.
  //  Call release_mme() to release hmmec from memory
  ////////////////////////////////////////////////////////////////////

  setup_mme(&rellrhs);
   double rellhood = 0.0;
   unsigned i;
   SESSION.warning = 0;
   if (!hmmec.solve(blupsol,rellrhs,"ysmp1")) return rellhood;
   rellhood=hmmec.info_vec[0];
   //rellhood = hmmec.logdet();
   SESSION.warning = 1;
   rellhood += yry - blupsol.inner_product(rellrhs);
   for (i=0; i<npattern; i++) rellhood += lnr0vec[i]*kvec[i];
   for (i=0; i<numterm; i++) {
      if (term[i].classi() == 'R'){
         rellhood += lng0vec[i]*term[i].nlevel();
      }
      else if (term[i].classi() =='P') {
         rellhood += lng0vec[i]*popsize;
      }
   }
   return (-0.5*rellhood);
}

double Model::minfun_vce(const Vector<double> &x,const int n)
{
  ////////////////////////////////////////////////////////////////////////
  // it returns negative restricted log likelihood.
  //   (it is for minimization program)
  //
  // x contains upper triangular part of variance matrices for residual
  // and other random terms in the model.
  // the hiden order in x is :
  //    residual--randtermlist[0]--randtermlist[1] .....
  // the order of randtermlist[i] is based on "random = A B A*B" given by user
  // user is totally responsible for size of vector x.
  // size of x == (1+nrandterm)*numtrait*(numtrait+1)/2
  //
  //  Important Notice:
  //  because restricted_log_likelihood(D) will be called sequentially
  //  under the same MME structure.
  //  hmmec is kept after restricted_log_likelihood(D) call.
  //  Call release_mme() to release hmmec from memory
  ///////////////////////////////////////////////////////////////////////////

   double mlhood;
   int i;
   vec2var(x);
   if (!residual_var.psd()) return (1.0e30);
   for (mlhood=0.0,i=0; i<numterm; i++) {
      if (term[i].classi() == 'R' || term[i].classi() == 'P') {
         if (!(term[i].prior->var_matrix()->psd())) {
            mlhood =1.0e30 + 1.0e10*(i+1);
            break;
         }
      }
   }
   if (mlhood == 0.0) mlhood = -1.0*restricted_log_likelihood();
   return ( mlhood );
}

void Model::uAu_trCA(Vector<double>& tsol,Vector<double> &uAu, Vector<double> &trca,Vector<double> &ivect,
                     Vector<double> &invect,const Vector<double> &ratio,const Vector<double> &zz) {

  non_zero = hmmec.close();
  unsigned i,j,k,t,ii,jj,je,startaddr,endaddr,nl;
  double rr,ee,xval,val;
  int *m_ia = hmmec.ia();
  int *m_ja = hmmec.ja();
  double *m_a = hmmec.a();
  double *sol;

  for (k=0,t=0; t<numterm; t++) {
    startaddr =  term[t].start + 1;
    sol = &(tsol[term[t].start]);
    nl = term[t].nlevel();
    switch (term[t].classi()) {
    case 'P':
      endaddr = startaddr+nl-1;
      for (ee=0.0,rr=0.0,ii=startaddr,i=0; i<nl; i++,ii++) {
	rr += (*sol * zz[i] * *sol);  sol++;
	memset(ivect.begin(),'\0',sizeof(double)*hmmesize);
	ivect[ii-1] = 1.0;
	if (!hmmec.solve(invect,ivect,"ysmp1")) break;
	je = m_ia[ii+1];
	for (j=m_ia[ii]; j<je; j++) {
	  jj = m_ja[j];
	  if (jj>=startaddr && jj<=endaddr) {
	    val = invect[jj-1];
	    xval = m_a[j]*val;
	    if (jj != ii) {ee += (xval + xval);}
	    else          {ee += (xval - val*zz[i]);}
	  }
	}
      }
      uAu[k] = hmmec.qTLW(tsol.begin(),tsol.begin(),startaddr,endaddr) - rr;
      trca[k++] = ee;
      break;
    case 'R':
      for (val=0.0,ii=startaddr-1,i=0; i<nl; i++,ii++) {
	 memset(ivect.begin(),'\0',sizeof(double)*hmmesize);
	 ivect[ii] = 1.0;
	 if (!hmmec.solve(invect,ivect,"ysmp1")) break;
	 val += invect[ii];
      }
      for (ee=0.0,i=0; i<nl; ++i) ee += *sol * *sol++;
      uAu[k] = ee*ratio[k];
      trca[k] = val*ratio[k];
      k++;
      break;
    }
  }
}

void Model::update_mme(const Vector<double> &ratio,const Vector<double> &oldratio,
                       const Vector<double> &zz)
{
   unsigned i,j,t,k,ii,jj,je,nl,startaddr,endaddr;
   double rr,val;
   int *m_ia = hmmec.ia();
   int *m_ja = hmmec.ja();
   double *m_a = hmmec.a();

   for (k=0,t=0; t<numterm; t++) {
      nl = term[t].nlevel();
      startaddr = term[t].start + 1;
      endaddr = startaddr + nl - 1;
      if (term[t].classi() == 'P') {
         rr = ratio[k]/oldratio[k];
         k++;
         for (ii=startaddr,i=0; i<nl; i++,ii++) {
            je = m_ia[ii+1];
            for (j=m_ia[ii]; j<je; j++) {
               jj = m_ja[j];
               if (jj >=startaddr && jj <=endaddr) {
                  val = m_a[j];
                  if (jj != ii) {val *= rr;}
                  else          {val = (val - zz[i]) * rr + zz[i];}
                  m_a[j] = val;
               }
            }
         }
      }
      else if (term[t].classi() == 'R') {
         val = ratio[k]-oldratio[k];
         k++;
         hmmec.add_diag(startaddr,endaddr,val);
      }
   }
}

void Model::vce_gibbs_sampler(const unsigned nvc,Vector<double> &vc,
                              const Vector<double> &ratio,const Vector<unsigned> &nlevel,
                              const double yy,double& ee,Vector<double> &uAu,
                              Vector<double> &bu,Vector<double> &wspace,const Vector<double> &zz)
{
   ////////////////////////////////////////////////////////
   // Gibbs sampling for b,u, and variance components
   // REQUIREMENTS:
   //     (1)  hmme must be already built
   //     (2)  if you want blupsol later, set it to zero
   ///////////////////////////////////////////////////////
   double *sol,rr,mean,se = std::sqrt(vc[nvc]);
   double *bsol = blupsol.begin() - 1;
   unsigned i,j,t,k,je,nl,startaddr,endaddr;
   int *m_ia = hmmec.ia();
   int *m_ja = hmmec.ja();
   double *m_a = hmmec.a();

   /////////////////////////////////////
   //  Gibbs sampling over mme
   /////////////////////////////////////
   memcpy(wspace.begin(),rellrhs.begin(),sizeof(double)*hmmesize);
   for (i=1; i<=hmmesize; i++) {
      mean = wspace[i - 1];
      k = m_ia[i];  je = m_ia[i+1];
      rr = m_a[k];                   // first element in each is the diagonal
      for (j=k+1; j<je; j++) mean -= m_a[j]*bu[m_ja[j] - 1];
      mean /= rr;
      bsol[i] += mean;
      mean += se*snorm()/std::sqrt(rr);     // variance = 1/diag
      bu[i - 1] = mean;
      for (j=k; j<je; j++) wspace[m_ja[j] - 1] -= m_a[j]*mean;
   }

   for (k=0,t=0; t<numterm; t++) {
      sol = &(bu[term[t].start]);
      switch (term[t].classi()) {
         case 'P':
            nl = nlevel[k];
            startaddr =  term[t].start + 1;
            endaddr = startaddr+nl-1;
            uAu[k] = hmmec.q(bu.begin(),bu.begin(),startaddr,endaddr);
            for (rr=0.0,i=0; i<nl; i++)  rr += (*sol * zz[i] * *sol++);
            uAu[k] = (uAu[k] - rr)/ratio[k];
            k++;
            break;
         case 'R':
            for (rr=0.0,j=0; j<nlevel[k]; ++j) rr += *sol * *sol++;
            uAu[k] = rr;
            k++;
            break;
      }
   }
   ee = yy - 2.0*bu.inner_product(rellrhs) + hmmec.q(bu.begin(),bu.begin());
   for (k=0; k<nvc; k++) {
      ee -=  uAu[k]*ratio[k];
      vc[k] = uAu[k]/genchi(static_cast<int>(nlevel[k]-2));     // sampling for vc
   }
   if (ee < 0.0) throw exception("Model::vce_gibbs_sampler(): you have probably found a bug!");
   vc[nvc] = ee/genchi(static_cast<int>(numobs-2));
}

void Model::add_G_1_single_trait(const Vector<double> &ratio)
{
   unsigned i,j,ii,jj,t,k,rec,nl,startaddr;
   double val,rr,xval;
   for (k=0,t=0; t<numterm; t++) {
      startaddr = term[t].startaddr();
      nl = term[t].nlevel();
      if (term[t].classi() == 'P') {
         Individual *I;
         unsigned asd[3];
         val = ratio[k++];
         Matrix<double> lambda(3,3);
         getlambda(lambda.begin(),3);
         for (rec=0; rec<nl; rec++) {
            I = pop->popmember[rec];
            asd[0]=I->id(); asd[1]=I->father_id(); asd[2] = I->mother_id();
            rr = val*4.0/(2.0-I->father_inbcoef()- I->mother_inbcoef());
            for (i=0; i<3; i++) {
               if (asd[i] != 0) {
                  ii = startaddr + asd[i];
                  for (j=0; j<3; j++) {
                     if (asd[j] != 0) {
                        jj = startaddr + asd[j];
                        if (jj >= ii) {
                           xval = lambda[i][j]*rr;
                           hmmec.insert(ii,jj,xval);
                        }
                     }
                  }
               }
            }
         }
      }
      else if (term[t].classi() == 'R') {
         val = ratio[k++];
         for (ii=startaddr+1,i=0; i<nl; i++,ii++) hmmec.insert(ii,ii,val);
      }
   }
}

double Model::setup_ww_single_trait(const Vector<double> &ratio,const unsigned pt,
                                    const unsigned ibeg,Vector<double> *zz)
{
   unsigned i,j,ii,jj,rec;
   double wt,xval,yy;
   hmmec.resize(hmmesize,max_nz);
   if (zz) zz->assign(0.0);
   rellrhs.assign(0.0);
   double *rhstmp = rellrhs.begin() -1;
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with    
   std::istringstream modelfile(modelstringstr);
   if (!modelfile) throw exception(" Model::setup_ww_single_trait(): cannot file");
   for (yy=0.0,rec=0; rec<numrec; rec++) {
      input_pos_val(modelfile);
      if (ntermGdist) {
         ii = rec_indid[rec];
         if (ii > 0) trait_vec[0] -= pop->popmember[ii-1]->xbzu();
      }
      wt = xval_term[numterm];             // wt = weight variable
      yy += trait_vec[0]*trait_vec[0]*wt;
      if (zz) (*zz)[term[pt].addr[0] - ibeg] += xval_term[pt]*wt;
      for (i=0; i<numterm; i++) {
         ii = term[i].addr[0];
         xval = xval_term[i]*wt;
         for (j=i; j<numterm; j++) {
            jj = term[j].addr[0];
            hmmec.insert(ii,jj,xval*xval_term[j]);
         }
         rhstmp[ii] += (xval*trait_vec[0]);
      }
   }
   //BRS modelfile.close();

   if (type == mixed_model) add_G_1_single_trait(ratio);
   // the new sparse routines don't work if closed.
   // non_zero = hmmec.close();
   //std::cout << hmmec << std::endl;
   return yy;
}

void Model::vce_emreml_single_trait(const int miter,const int intive,
                                    const double stoptol)
{
   rellrhs.resize(hmmesize);
   blupsol.resize(hmmesize);

   int maxiter = miter;
   unsigned i,k,t;
   unsigned pedigree_term;
   unsigned rank_mme,startaddr,iter = 0;
   unsigned nvc = totalnvc() - 1;              // number of ratio's
   double coef = 0.2;
   double sumtrace,ee,reldiff,yy,val;
   Vector<double> zz;
   Vector<double> ivect;
   Vector<double> invect;

   Vector<double> alfa1(nvc+1);
   Vector<double> alfa2(nvc+1);
   Vector<double> vc(nvc+1);
   Vector<double> oldvc(nvc+1);
   Vector<double> ratio(nvc);
   Vector<double> oldratio(nvc);
   Vector<double> uAu(nvc);     // including ratio
   Vector<double> trca(nvc);     // including ratio
   Vector<double> d1(nvc);
   Vector<double> d2(nvc);

   var2vec(vc);
   for (i=0; i<nvc; i++) alfa1[i] = ratio[i] = vc[nvc]/vc[i];
   for (t=0; t<numterm; t++) {
      if (term[t].classi() == 'P') {
         startaddr = term[t].start + 1;
         pedigree_term = t;
         zz.resize(popsize);
         ivect.reserve(hmmesize);
         invect.reserve(hmmesize);
      }
      else if (term[t].classi() == 'R') {
         ivect.reserve(hmmesize);
         invect.reserve(hmmesize);
      }
   }
   SESSION.warning = 0;
   int W = SESSION.output_precision + 6;            // 6 = +.e+00
   std::cout.precision(SESSION.output_precision);
   std::cout << " iteration  sigma_1 ..... sigma_e  restricted_log_likelihood\n";
   std::cout << std::setw(5) << iter << ": ";
   for (i=0; i<=nvc; i++) std::cout << " " << std::setw(W) << vc[i];
      do {
         iter++;
      yy = setup_ww_single_trait(ratio,pedigree_term,startaddr,&zz);
      if (!hmmec.solve(blupsol,rellrhs)) break;
      reml_value = hmmec.logdet(); // log determinant
      rank_mme = hmmec.irank();
      uAu_trCA(blupsol,uAu, trca,ivect,invect,ratio,zz);
      ee = yy - blupsol.inner_product(rellrhs);
      reml_value += ee/vc[nvc];
      for (sumtrace=0.0,i=0; i<nvc; i++) {
         sumtrace += trca[i];
         ee -= uAu[i];
      }
      if (ee < 0.0) throw exception("Model::vce_emreml_single_trait(): a bug!");

      /////////////////////////////////////////////////
      // estimate of residual error variance,ve
      //////////////////////////////////////////////////
      reml_value += std::log(vc[nvc])*static_cast<double>(double(numobs)-double(hmmesize));
      vc[nvc] = ee/(sumtrace + (double(numobs) - double(rank_mme))); 
      oldvc = vc;
      for (k=0,t=0; t<numterm; t++) {
         if ( term[t].classi() == 'P' || term[t].classi() == 'R') {
            val = static_cast<double>(term[t].nlevel());
            reml_value += val*std::log(oldvc[k]);
            vc[k] = uAu[k]/(ratio[k]*(val - trca[k]));
            k++;
         }
      }
      reml_value *= -0.5;
      //////////////////////////////////////////////////////////////////
      //  try to compare the difference between estimates at (t-1)
      //  and t round
      ////////////////////////////////////////////////////////////////
      reldiff = fabs((vc[nvc]-oldvc[nvc])/oldvc[nvc]);
      for (i=0; i<nvc; i++) {
         reldiff = std::max(reldiff,fabs((vc[i]-oldvc[i])/oldvc[i]));
      }
      for (i=0; i<nvc; i++) {
         oldratio[i] = ratio[i];
         ratio[i] =  vc[nvc]/vc[i];
      }
      if (iter == 1) { // extrapolation stuff
         for (i=0; i<nvc; i++) {
            d1[i] = ratio[i] - oldratio[i];
            alfa2[i] = ratio[i] = oldratio[i] + coef;
         }
      }
      else if (iter == 2) {
         for (i=0; i<nvc; i++)  d2[i] =  ratio[i] - oldratio[i];
         for (i=0; i<nvc; i++) {
            if (fabs(d2[i]-d1[i]) > SESSION.epsilon) {
               ratio[i] = (alfa1[i]*d2[i]-alfa2[i]*d1[i])/(d2[i]-d1[i]);
            }
            else {
               ratio[i] = 0.0;
            }
            if (ratio[i] <= stoptol) {
               coef += 1.0;
               iter--;
               for (k=0; k<nvc; k++) ratio[k] = oldratio[k] + coef;
               break;
            }
         }
      }
      if (iter > 3 || iter==1){
	std::cout.precision(SESSION.output_precision + 4);
	std::cout << " " << std::setw(W+4) << reml_value << std::endl;
	std::cout.precision(SESSION.output_precision);
      }
      else {
	std::cout << std::endl;
      }

      if (iter>=maxiter ) {
         std::cout << std::endl;
         if (intive) {
            std::cout << " how many more iterations do you want?" << std::endl;
            std::cin >> k;
            if (k==0) break;
            maxiter += k;
         }
         else {
            std::cout<<"maximun iteration has been exceeded without convergence\n";
            break;
         }
      }
      std::cout << std::setw(5) << iter << ": ";
      for (i=0; i<=nvc; i++) std::cout << " " << std::setw(W) << vc[i];
   } while (reldiff > stoptol);

   //////////////////////////////////////////////////////////////////
   // get blup solution and likelihood for vce from last iteration
   ////////////////////////////////////////////////////////////////////
   yy = setup_ww_single_trait(ratio,pedigree_term,startaddr,&zz);
   if (!hmmec.solve(blupsol,rellrhs,"ysmp1")) throw exception(" Model::vce_emreml_single_trait(): returned with error");
   reml_value = hmmec.logdet();
   for (k=0,t=0; t<numterm; t++) {
      if (term[t].classi() == 'P' || term[t].classi() == 'R') {
         reml_value += static_cast<double>(term[t].nlevel()) * std::log(vc[k++]);
      }
   }
   reml_value += std::log(vc[nvc])*(numobs-hmmesize);
   reml_value += (yy - blupsol.inner_product(rellrhs))/vc[nvc];
   reml_value *= -0.5;
   std::cout.precision(SESSION.output_precision + 4);
   std::cout << " " << std::setw(W+4) << reml_value << std::endl;
   std::cout.precision(SESSION.output_precision);

   vec2var(vc);

   SESSION.warning = 1;
   hmmec.close();
   info(std::cout);
   hmmec.resize(0,0);
}

void Model::vce_emreml_multi_trait(const int miter,const int intive,
                                   const double stoptol)
{
   if (npattern != 1) throw exception(" Model::vce_emreml_multi_trait(): missing data isn't allowed");
   int maxiter = miter;
   int iter = 0;
   unsigned i,j,ii,jj,t,t1,t2,rec,pt,nlevel;
   unsigned rank_mme,startaddr,endaddr = 0;
   double *rh,*rhs,*Lii;
   double val,xval,rr,reldiff,ee,u1Au2;
   Vector<double> *solmat; 
   if(numtrait>0){
     solmat = new Vector<double> [numtrait];
   }
   else {
     solmat = 0;
   }

   hmmesize /= numtrait;
   max_nz = est_nze(hmmesize,1,popsize);
   hmmec.resize(hmmesize,max_nz);
   rellrhs.resize(hmmesize);
   blupsol.resize(0);
   for (term[0].start = 0, t=1; t<numterm; t++) {
      term[t].start = term[t-1].startaddr() + term[t-1].nlevel();
   }

   int nvc = 1;
   Vector<double> zz;
   Vector<double> ivect;
   Vector<double> invect;

   Vector<double> oldratio(nvc);
   Matrix<double> rhsmat(hmmesize,numtrait);

   Vector<Vector<double> > ratio(numtrait);
   Vector<Vector<double> > uAu(numtrait);
   Vector<Vector<double> > trca(numtrait);
   for (i=0; i<nvc; ++i) {
      ratio[i].reserve(nvc);
      uAu[i].reserve(nvc);
      trca[i].reserve(nvc);
   }

   doubleMatrix yy(numtrait,numtrait);
   doubleMatrix tyy(numtrait,numtrait);
   doubleMatrix R(numtrait,numtrait);
   doubleMatrix G(numtrait,numtrait);
   doubleMatrix oldR(numtrait,numtrait);
   doubleMatrix oldG(numtrait,numtrait);
   doubleMatrix tR(numtrait,numtrait);
   doubleMatrix tG(numtrait,numtrait);
   doubleMatrix Ti(numtrait,numtrait);
   doubleMatrix Tti(numtrait,numtrait);
   doubleMatrix P(numtrait,numtrait);
   doubleMatrix L(numtrait,numtrait);
   doubleMatrix Li(numtrait,numtrait);
   Vector<double> D(numtrait);

   for (i=0,t=0; t<numterm; t++) {
      if (term[t].classi() == 'P' || term[t].classi() == 'R') {
         i++;
         pt = t;
         nlevel = term[t].nlevel();
         startaddr = term[t].startaddr();
         endaddr = startaddr + nlevel;
         ivect.reserve(hmmesize);
         invect.reserve(hmmesize);
         if (term[t].classi() == 'P')  zz.resize(nlevel);
      }
   }
   if (i != 1) throw exception(" Model::vce_emreml_multi_trait(): one and only one random effect is allowed");
   for (t1=0; t1<numtrait; t1++) {
      for (t2=0; t2<numtrait; t2++) yy[t1][t2] = 0.0;
      for (i=0; i<hmmesize; i++) rhsmat[i][t1] = 0.0;
   }
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);
   if (!modelfile) throw exception(" Model::vce_emreml_multi_trait():cannot open");
   for (rec=0; rec<numrec; rec++) {
      modelfile.read((char *)pos_term, sizeof(unsigned)*(numterm+1));
      modelfile.read((char *)xval_term, sizeof(double)*numterm);
      modelfile.read((char *)trait_vec, sizeof(double)*numtrait);
      for (i=0; i<numterm; i++) {
         term[i].addr[0] = term[i].start + pos_term[i];
      }
      if (!zz.empty()) zz[term[pt].addr[0] - startaddr -1] += xval_term[pt];
      for (i=0; i<numterm; i++) {
         ii = term[i].addr[0];
         xval = xval_term[i];
         for (j=i; j<numterm; j++) {
            jj = term[j].addr[0];
            hmmec.insert(ii,jj,xval*xval_term[j]);
         }
         rh = rhsmat[ii-1];
         for (t1=0; t1<numtrait; t1++) rh[t1] += xval*trait_vec[t1];
      }
      for (t1=0; t1<numtrait; t1++) for (t2=0; t2<numtrait; t2++) {
         yy[t1][t2] +=  trait_vec[t1]*trait_vec[t2];
      }
   }
   //BRS modelfile.close();
   oldratio[0] = 1.0;  // ratio has only one element
   add_G_1_single_trait(oldratio);
   non_zero = hmmec.close();
   hmmec.reorder();
   rank_mme = hmmec.factorization(5);

   int W = SESSION.output_precision + 6;            // 6 = +.e+00
   std::cout.precision(SESSION.output_precision);

   R = residual_var;
   G = *(term[pt].prior->var_matrix());
   std::cout << " iteration:    G0      and      R0\n";
   std::cout << std::setw(5) << iter << ": ";
   for (t1=0; t1<numtrait; t1++) for (t2=t1; t2<numtrait; t2++) {
      std::cout << " " << std::setw(W) << G[t1][t2];
   }
   std::cout << "\n" << std::setw(7) << " ";
   for (t1=0; t1<numtrait; t1++) for (t2=t1; t2<numtrait; t2++) {
      std::cout << " " << std::setw(W) << R[t1][t2];
   }
   std::cout << "\n\n";
   SESSION.warning = 0;
   double *sol1,*sol2;
   do {
      iter++;
      oldR = R;
      oldG = G;
      R.sqrtm();          // R turns into a lower triangular T
      Tti = R.transpose();    Tti.inv();
      Ti = R;         Ti.inv();
      P = Ti*G*Tti;
      D = P.eigen();
      L = R*P;
      Li = L;        Li.inv();
      tyy = Li*yy*Li.transpose();
      for (t1=0; t1<numtrait; t1++) {
         ratio[t1][0] = 1.0/D[t1];
         Lii = Li[t1];
         rhs = rellrhs.begin();
         for (i=0; i<hmmesize; i++) {
            rh = rhsmat[i];
            for (rr=0.0,t=0; t<numtrait; t++) rr += rh[t]*Lii[t];
            *rhs++ = rr;
         }
         update_mme(ratio[t1],oldratio,zz);
         oldratio[0] = ratio[t1][0];
         if (hmmec.solve(solmat[t1],rellrhs,"ysmp1")) break;
         uAu_trCA(solmat[t1],uAu[t1],trca[t1],ivect,invect,ratio[t1],zz);
         xval = static_cast<double>(nlevel) - trca[t1][0];
         tG[t1][t1] = uAu[t1][0]/(ratio[t1][0]*xval);
         ee = tyy[t1][t1] - solmat[t1].inner_product(rellrhs) - uAu[t1][0];
         val = trca[t1][0] + (numobs - rank_mme);
         tR[t1][t1] = ee/val;
         for (t2=0; t2<t1; t2++) {
            sol1 = &(solmat[t1][startaddr]);
            sol2 = &(solmat[t2][startaddr]);
            if (!zz.empty()) {
               for (val=0.0,i=0; i<nlevel; i++) val += (*sol1++ * zz[i] * *sol2++);
               u1Au2 = hmmec.q(solmat[t1].begin(),solmat[t2].begin(),startaddr+1,endaddr);
               u1Au2 = (u1Au2-val)/ratio[t1][0];
            }
            else {
               for (u1Au2=0.0,i=0; i<nlevel; i++) u1Au2 += *sol1++ * *sol2++;
            }
            tG[t2][t1] = tG[t1][t2] = u1Au2/nlevel;
            ee = tyy[t1][t2] - solmat[t2].inner_product(rellrhs) - ratio[t2][0]*u1Au2;
            tR[t2][t1] = tR[t1][t2] = ee/numobs;
         }
      }
      R = L*tR*L.transpose();
      G = L*tG*L.transpose();
      reldiff = 0.0;
      for (t1=0; t1<numtrait; t1++) for (t2=t1; t2<numtrait; t2++) {
         reldiff = std::max(reldiff,fabs((R[t1][t2]-oldR[t1][t2])/oldR[t1][t2]));
         reldiff = std::max(reldiff,fabs((G[t1][t2]-oldG[t1][t2])/oldR[t1][t2]));
      }
      std::cout << std::setw(5) << iter << ": ";
      for (t1=0; t1<numtrait; t1++) for (t2=t1; t2<numtrait; t2++) {
         std::cout << " " << std::setw(W) << G[t1][t2];
      }
      std::cout << "\n" << std::setw(7) << " ";
      for (t1=0; t1<numtrait; t1++) for (t2=t1; t2<numtrait; t2++) {
         std::cout << " " << std::setw(W) << R[t1][t2];
      }
      std::cout << "\n\n";
      if (iter >= maxiter) {
         std::cout << std::flush;
         if (intive) {
            std::cout << "how many more iterations do you want?  0 for end" << std::endl;
            std::cin >> t;
            if (t == 0) break;
            maxiter += t;
         }
         else {
            std::cout<<"maximun iteration has been exceeded without convergence\n";
            break;
         }
      }
   } while(reldiff > stoptol);
   std::cout << "variance components for random effect\n";
   G.print();
   std::cout << "variance components for residual error\n";
   R.print();
   residual_var = R;
   term[pt].prior->var_matrix()->copy(G);
   if (solmat){
     delete [] solmat;
     solmat = 0;
   }
   SESSION.warning = 1;
   info(std::cout);
   hmmec.resize(0,0);
}

/*!
 * Estimate variance component using EM-REML.
 *
 *
 * \sa vce_dfreml()
 */
Vector<double>* Model::vce_emreml(const int miter,const int intive,
                          const double stoptol)
{
   if (type == bad_model) throw exception("Model::vce_emreml(): bad model");
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return 0;
      }
   }

   if (type == mixed_model) {
      if (residual_var.num_rows() != numtrait) {
         type = bad_model;
         throw exception("Model::vce_emreml(): residual variance has not yet assigned");
      }
   }
   else if (type == fixed_model) {
      if (residual_var.num_rows() != numtrait) {
         residual_var.identity(numtrait,numtrait);
      }
   }
   else {
      throw exception("Model::vce_emreml(): inappropriate model");
   }

   if (numtrait == 1) {
      vce_emreml_single_trait(miter,intive,stoptol);
   }
   else if (numtrait > 1) {
      vce_emreml_multi_trait(miter,intive,stoptol);
   }
   else {
      throw exception("Model::vce_emreml(): no trait available");
   }
   rellrhs.resize(0);
   hmmec.release();
   return &blupsol;
}

/*!
 * Estimate variance component using DF-REML.
 *
 * \sa vce_emreml()
 */
Vector<double>* Model::vce_dfreml(const std::string &method,int maxiter,double funtol)
{
   if (type == bad_model) throw exception("Model::vce_dfreml(): bad model");
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return 0;
      }
   }

   if (type == mixed_model) {
      if (residual_var.num_rows() != numtrait) {
         type = bad_model;
         throw exception("Model::vce_dfreml(): residual variance has not yet assigned");
      }
   }
   else if (type == fixed_model) {
      if (residual_var.num_rows() != numtrait) {
         residual_var.identity(numtrait,numtrait);
      }
   }
   else {
      throw exception("Model::vce_dfreml(): inappropriate model");
   }

   minfun_indx = 1;
   int nvc,i,iter = maxiter;
   double mlhood;
   double ftol = funtol;

   nvc = totalnvc();
   Vector<double> x(nvc);
   if (method == "powell") {
      var2vec(x);
      for (i=0; i<nvc; i++) {
         if (fabs(x[i]) <= 10.0*SESSION.epsilon) x[i] = 0.0001;
      }
      if (ftol<1.0e-8) ftol = 1.0e-8;
      mlhood = praxis(x,nvc,iter,ftol*ftol,ftol);
      vec2var(x);
   }
   else if (method == "nr_powell") {
      Matrix<double> p(nvc+1,nvc);
      for (i=0; i<nvc; i++) {
         p[i][i] = 1.0;
      }
      var2vec(x);
      for (i=0; i<nvc; i++) {
         if (fabs(x[i]) <= 10.0*SESSION.epsilon) x[i] = 0.0001;
      }
      mlhood = nr_powell(data,this,x,p,nvc,iter,ftol);
      vec2var(x);
   }
   else {               // I will be adding more methods later
      throw exception(" Model::dfreml(): bad arg");
   }
   minfun_indx = 0;
   //hmmec.release();
   rellrhs.resize(0);

   nfunk_in_dfreml = iter;
   reml_value = -mlhood;
   dfreml_method = method;
   dfreml_called = 1;
   info(std::cout);
   return &blupsol;
}

Vector<double>* Model::vce_gibbs(const unsigned warmup,const unsigned gibbslen)
{
   if (type == bad_model) throw exception("Model::vce_gibbs(): model is too bad");
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return 0;
      }
   }

   if (numtrait != 1) throw exception("Model::vce_gibbs(): only for single trait model");
   if (numobs < 7) throw exception("Model::vce_gibbs(): number of records must > 6");
   hmmec.resize(hmmesize,max_nz);
   rellrhs.resize(hmmesize);
   blupsol.resize(hmmesize);

   unsigned nvc = totalnvc() - 1;                     // number of ratio's
   unsigned i,k,t,pt,nl,startaddr;
   double yy,ee;
   Vector<double> zz;
   Vector<double> vc(nvc+1);
   Vector<double> bu; bu.reserve(hmmesize);
   Vector<double> uAu(nvc);
   Vector<double> vchat(nvc+1);
   Vector<double> vcsquare(nvc+1);
   Vector<double> ratio(nvc);
   Vector<unsigned> nlevel(nvc);
   Vector<double> oldratio(nvc);
   Vector<double> wspace; wspace.reserve(hmmesize);

   int W = SESSION.output_precision + 6;            // 6 = +.e+00
   std::cout.precision(SESSION.output_precision);
   var2vec(vc);
   for (k=0; k<nvc; k++) ratio[k] = vc[nvc]/vc[k];
   for (k=0,t=0; t<numterm; t++) {
      if (term[t].classi() == 'P') {
         nlevel[k] = term[t].nlevel();
         startaddr = term[t].start + 1;
         pt = t;
         zz.reserve(nlevel[k]);
         k++;
      }
      else if (term[t].classi() == 'R') {
         nlevel[k++] = term[t].nlevel();
      }
   }
   for (k=0; k<nvc; k++) {
      if (nlevel[k] < 7) throw exception(" Model::vce_gibbs(): random_effects levels > 6");
   }
   yy = setup_ww_single_trait(ratio,pt,startaddr,&zz);
   ///////////////////////////////////////////////////////////////
   //  Gibbs Sampling:  WARM UP
   ////////////////////////////////////////////////////////////////
   for (i=0; i<hmmesize; i++) bu[i] = snorm();   // initial values
   vce_gibbs_sampler(nvc,vc,ratio,nlevel,yy,ee,uAu,bu,wspace,zz);
   for (k=0; k<nvc; k++) {
      oldratio[k] = ratio[k]; ratio[k] = vc[nvc]/vc[k];
   }
   for (i=1; i<warmup; i++) {
      update_mme(ratio,oldratio,zz);
      vce_gibbs_sampler(nvc,vc,ratio,nlevel,yy,ee,uAu,bu,wspace,zz);
      for (k=0; k<nvc; k++) {oldratio[k] = ratio[k]; ratio[k] = vc[nvc]/vc[k];}
   }
   ///////////////////////////////////////////////////////////////
   //   Gibbs Sampling
   ////////////////////////////////////////////////////////////////
   blupsol.assign(0.0);
   vchat.assign(0.0);
   vcsquare.assign(0.0);
   for (i=0; i<gibbslen; i++) {
      update_mme(ratio,oldratio,zz);
      vce_gibbs_sampler(nvc,vc,ratio,nlevel,yy,ee,uAu,bu,wspace,zz);
      for (k=0; k<nvc; k++) {
         nl = nlevel[k];
         vchat[k] += uAu[k]/(nl-4);
         vcsquare[k] += uAu[k]*uAu[k]/((nl-6)*(nl-4));
         oldratio[k] = ratio[k];
         ratio[k] = vc[nvc]/vc[k];
      }
      vchat[nvc] += ee/(numobs-4);
      vcsquare[nvc] += ee*ee/((numobs-6)*(numobs-4));
   }
   for (k=0; k<=nvc; k++) {
      vchat[k] /= gibbslen;
      vcsquare[k] = vcsquare[k]/gibbslen - vchat[k]*vchat[k];
   }
   blupsol /= gibbslen;
   std::cout.precision(SESSION.output_precision);
   std::cout << " posterior mean and variance of variance components \n";
   for (k=0; k<=nvc; k++) {
      std::cout << " " << std::setw(8) << k;
      std::cout << " " << std::setw(W) << vchat[k] << " " << std::setw(W) << vcsquare[k]<<"\n";
   }
   vec2var(vchat);
   hmmec.release();
   rellrhs.resize(hmmesize);
   info(std::cout);
   return &blupsol;
}

void Model::info(std::ostream& stream) const
{
   if (type == bad_model) throw exception("Model::info(stream): model is too bad");
   stream.precision(SESSION.output_precision);
   stream <<"\n             some extra information in the model\n";
   stream << " --------------------------------------------------------\n";
   if (dfreml_called) {
      stream << "optimizor used in dfreml = " << dfreml_method <<"\n";
      stream <<"# of function evaluations used in dfreml = "
             << nfunk_in_dfreml<<"\n";
      stream << "maximum log restricted likelihood in dfreml = "
             <<  reml_value <<"\n";
   }
   else {
      ;
   }
   const ModelTerm *T;
   int i;
   for (i=0; i<numterm; i++) {
      T = &(term(i));
      if (T->classi() == 'R' || T->classi() == 'P') {
         stream << " variance for "
              << factor_struct[T->factorindx[0]].name() << " = \n";
         T->prior->var_matrix()->print(stream);
      }
   }
   stream << " residual variance =\n" << residual_var;

   if (data_prepared) {
      stream << "\n";
      stream << " MME dimension   : " << hmmesize << "\n";
      stream << " non-zeros in MME: " << hmmec.nz() << "\n\n";

      stream << "       basic statistics for dependent variables\n";
      stream << "  -----------------------------------------------------\n";
      stream << " trait-name              n          mean         std\n";
      unsigned n;
      for (i=0; i<numtrait; i++) {
         n = trait_struct[i]->nlevel();
         stream << " " << std::setw(16) << trait_struct[i]->name()
                << " " << std::setw(10) << n
                << " " << std::setw(12) << trait_struct[i]->mean();
         if (n==1) {
            stream << " " << std::setw(12) << "." << "\n";
         }
         else {
            stream<< " " << std::setw(12) << trait_struct[i]->std() << "\n";
         }
      }
   }
   stream << " --------------------------------------------------------\n";
}
} ///////// end of namespace matvec

