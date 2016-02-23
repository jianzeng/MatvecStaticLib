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

#include <iostream>
#include <iomanip>
#include <sstream>

#include "session.h"
#include "util.h"
#include "doublematrix.h"
#include "model.h"
#include "population.h"
#include "statdist.h"
#include "glmm.h"

namespace matvec {
#undef DO_THREADS
#ifdef SKIP_THREADS
#undef  DO_THREADS
#endif

#undef DO_PMAT
#define DO_LOG
#undef DO_CHOL

#ifdef DO_CHOL
#undef DO_PMAT
#undef DO_LOG
#endif

#ifdef DO_LOG
#undef DO_PMAT
#undef DO_CHOL
#endif

#define DO_HINV

#undef DEBUG

#define EPSILON  SESSION.epsilon
   //#define act_numtrait numtrait // Temporary fix

 extern void getlambda(double **lambda, const int n);

 void GLMM::Vec2Var(const double *x)
   ///////////////////////////////////////////////////////////////////////
   // put variance components into a double vector
   // user is totally responsible to have an appropriate size for vector x
   ///////////////////////////////////////////////////////////////////////
 {
   int i,k,t1,t2,ped_done;
   ped_done=0;
   doubleMatrix *var;
 #ifdef DO_CHOL
   doubleMatrix L;
 #endif
   for (k=0,i=0; i<numterm; i++) {
     if (term[i].classi() == 'R' || (!ped_done && term[i].classi() == 'P')) {
       if(var_link[i]==i){
	 int nt=nt_vec[i];
#ifdef DO_CHOL
	 L.resize(numtrait,numtrait);
	 L.assign(0.0);
	 var=&L;
#else
	 var = term[i].prior->var_matrix();
#endif
	 for (t1=0; t1<nt; t1++) {
	   for (t2=t1; t2<nt; t2++,k++) {
	     (*var)[t1][t2] = x[k];
#ifndef DO_CHOL
	     if ( t2>t1 ) (*var)[t2][t1] = x[k];
#endif
	   }
	 }
#ifdef DO_CHOL
	 var= term[i].prior->var_matrix();
	 *var=L.transpose()*L;
#endif
       }
       else{
	 int ii=var_link[i];
	 var= term[i].prior->var_matrix();
	 doubleMatrix *var2;
	 var2=term[ii].prior->var_matrix();
	 *var=*var2;
       }
     }       
       //if(corrmap[i]) corrvar[i][i]=*var;
   }
   /*
     if(ncorr) {
     for(int ict=0;ict<numterm;ict++) {
     if(corrmap[ict]) {
     for(int jct=(ict+1);jct<numterm;jct++) {
     if(corrmap[jct]) {
     for (t1=0; t1<numtrait; t1++) for (t2=0; t2<numtrait; t2++,k++) {
     corrvar[ict][jct].me[t1][t2]=x[k];
     }
     corrvar[jct][ict]=corrvar[ict][jct];
     }
     }
     }
     }
     }
   */
   if(link_function){
     for(t1=0;t1<res_nvc;t1++) varcomp[t1]=x[k++];
   }
   else{
     for (t1=0; t1<numtrait; t1++) {
       for (t2=t1; t2<numtrait; t2++,k++) {
 	residual_var[t1][t2] = x[k];
 	if ( t2>t1 ) residual_var[t2][t1] = x[k];
       }
     }
   }
 }


 Vector<double> * GLMM::glim(int iterations,double **ke, int nr, int nc, double *mean)
 {
   int iter,i,j,k;
   doubleMatrix kvk;
   Vector<double> kpm,wrkvec;
   double old_like,new_like,*krow,*m,*tmpsol,kpd;
   double *kp,*sol;
   Vector<double> old_sol;
   Vector<double> *retval;
   if(!link_function) {
     blup("ysmp1");
     if(ke)
       {
	 tmpsol=new double [hmmesize];
	 krow=new double [hmmesize];
       }
     if(ke)
       {
 	kpm.resize(nr);
 	kvk.resize(nr,nr);
 	if(mean) memcpy(kpm.begin(),mean,sizeof(double)*nr);
 	for(i=0;i<nr;i++){
 	  for(j=0;j<nc;j++)  kpm[i]-=ke[i][j]*blupsol[j];
 	  memset(krow,'\0',sizeof(double)*hmmesize);
 	  memcpy(krow,ke[i],sizeof(double)*nc);
 	  hmmec.solve(tmpsol,krow,"ysmp");
 	  for (k=0; k<nr; k++) {
 	    kp = ke[k]; sol = tmpsol;
 	    for (kpd=0.0,j=0; j<nc; j++) kpd += *kp++ * *sol++;
 	    kvk[k][i] = kpd;
 	  }
 	}
 	kvk.ginv1();
 	kpm=kvk*kpm;
 	memset(krow,'\0',sizeof(double)*hmmesize);
 	for(i=0;i<nr;i++) {
 	  for(j=0;j<nc;j++) { 
 	    krow[j]+=ke[i][j]*kpm[i];
 	  }
 	}
 	hmmec.solve(tmpsol,krow,"ysmp");
 	for(i=0;i<hmmesize;i++) blupsol[i]+=tmpsol[i];
	
 	delete [] krow;
 	delete [] tmpsol;	
       }



     return(&blupsol);
   }
   for(iter=0;iter< iterations;iter++) {
     setup_mme(&rellrhs); 
#ifdef DEBUG
     doubleMatrix LHS;
     LHS=mmec()->dense();
     std::cout << "\nRHS\n"<<rellrhs.subvec()<<"\n";
     std::cout << "\nLHS\n"<<LHS <<"\n";
#endif
     blup("ysmp1");

     if(ke && iter==0)
       {
	 tmpsol=new double [hmmesize];
	 krow=new double [hmmesize];
       }
     if(ke)
       {
	 kpm.resize(nr);
	 kvk.resize(nr,nr);
	 if(mean) memcpy(kpm.begin(),mean,sizeof(double)*nr);
	 for(i=0;i<nr;i++){
	   for(j=0;j<nc;j++)  kpm[i]-=ke[i][j]*blupsol[j];
	   memset(krow,'\0',sizeof(double)*hmmesize);
	   memcpy(krow,ke[i],sizeof(double)*nc);
	   hmmec.solve(tmpsol,krow,"ysmp");
	   for (k=0; k<nr; k++) {
	     kp = ke[k]; sol = tmpsol;
	     kpd=0.0;
	     for (j=0; j<nc; j++) kpd += ke[k][j] * tmpsol[j];
	     kvk[k][i] = kpd;
	   }
	 }
	 kvk.ginv1();
	 kpm=kvk*kpm;
	 memset(krow,'\0',sizeof(double)*hmmesize);
	 for(i=0;i<nr;i++) {
	   for(j=0;j<nc;j++) {
	     krow[j]+=ke[i][j]*kpm[i];
	   }
	 }
	 hmmec.solve(tmpsol,krow,"ysmp");
	 for(i=0;i<hmmesize;i++) blupsol[i]+=tmpsol[i]; 
       }
     
   }  
   
   
   retval=&blupsol;
   if(ke)
     {
       delete [] krow;
       delete [] tmpsol;
     }
   return retval;
 }
  

 double GLMM::contrast(const double **kpme, const unsigned nr,
 		      const unsigned nc,const double *M,const int prt_flag,double **result)
 {
   ////////////////////////////////////////////////////
   // H0:  K'*b = M
   // trailing zeros in K' can be omitted
   // REMINDER: blupsol & rellrhs remain intact
   //
   // meanp must be declared as double meanp[nr][4]
   /////////////////////////////////////////////////////


   if(!link_function){
     return(Model::contrast(kpme,nr,nc,M, prt_flag,result));
   }
   int do_varcomp=Info.num_rows();
   int Print_Flag=1;
   if(Print_Level < 0 || !prt_flag) Print_Flag=0;

   if (!kpme) {
     warning("GLMM::contrast(kp): kp is null, thus it is ignored");
     return 0.0;
   }
   if (!get_blupsol()) throw exception("GLMM::contrast(): bad model");

   if (nc > hmmesize) throw exception("GLMM::contrast(): size not conformable");
   int est=1;
   doubleMatrix kvk(nr,nr);
   doubleMatrix kvkori(nr,nr);
   doubleMatrix kv(nr,hmmesize);
   doubleMatrix kGk(nr,nr);
   doubleMatrix kRk(nr,nr);
   Vector<double> kpdiff(nr);
   Vector<double> kpb_m(nr);
   Vector<double> xy(hmmesize);
   Vector<double> tmpsol(hmmesize);
   double *sol, kpd;
   double log_LR;
   log_LR=2.*like_val;
   const double *kp;
   unsigned i,j,k;
   int t1,t2,ii,*ainv_ia,*ainv_ja,ainv_nrow;
   int ipos,jpos,je,jj,jaddr;
   double *ainv_a;
   if(LR==2) glim(10,(double **)kpme,nr,nc,(double *)M);//Lagrange Multiplier uses KVK under Ho

   for (i=0; i<nr; i++) {
     if (nc < hmmesize) {
       memset(xy.begin(),'\0',sizeof(double)*hmmesize);
       memcpy(xy.begin(),kpme[i],sizeof(double)*nc);
       hmmec.solve(tmpsol,xy,"ysmp");
     }
     else {
       hmmec.solve(tmpsol.begin(),kpme[i],"ysmp");
     }
     hmmec.mv(tmpsol,xy);
     memcpy(kv[i],tmpsol.begin(),sizeof(double)*hmmesize);
     kp = kpme[i]; 
     sol = xy.begin();
     for (kpd=0.0,j=0; j<nc; j++) {
       kpd = std::max(kpd,fabs(kpme[i][j] - xy[j]));
       kp++;
       sol++;
     }
     for (j=nc; j<hmmesize; j++) {
       kpd = std::max(kpd,fabs(xy[j]));
       sol++;
     }
     kpdiff[i] = kpd;
     for (k=0; k<nr; k++) {
       kp = kpme[k]; sol = tmpsol.begin();
       for (kpd=0.0,j=0; j<nc; j++) kpd += kpme[k][j] * tmpsol[j];
       kvk[k][i] = kpd;
     }
   }
   for(k=1;k<nr;k++) for(i=0;i<k;i++) kvk[k][i]=kvk[i][k];
   kvkori=kvk;
   kvk = kvkori;
   kvk.ginv1(&k);
   if (k < nr) {
     warning("GLMM::contrast(K'): K'V{-1}K is singular. Possible reason:\n"
 	    " (1) K is non-estimable, and (2) K isn't of full column rank");
     //return 0.0;
   }
   // Get Lambda's for vc calculations
   int nvc=TotalNvc();
   unsigned startaddr,nlevels,iaddr,t;
   Vector<double> Lambda(nvc);
   double addj=0.;
   doubleMatrix *smat;
   int k_vc,i_vc;
   doubleMatrix wmat(numtrait,numtrait);
   doubleMatrix rmat(numtrait,numtrait);
   for(k_vc=0,i_vc=0;do_varcomp && i_vc<numterm;i_vc++) {
     if(term[i_vc].classi() == 'P' || term[i_vc].classi() == 'R'){
       startaddr=term[i_vc].startaddr()+1;
       int nt=nt_vec[i_vc];
       for(t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,k_vc++) {
	
 	rmat=*(term[i_vc].prior->var_matrix());
 	rmat.ginv1();
 	wmat=SinMat[k_vc];
 	kGk.assign(0.0);
 	kRk.assign(0.0);
 	if(term[i_vc].classi() == 'R') {
 	  nlevels=term[i_vc].nlevel();
 	  for(ii=1;ii<=nlevels;ii++) {
 	    iaddr=startaddr+(ii-1)*numtrait;
 	    kGk+=kv.submat(0,iaddr-1,nr,nt)*wmat
 	      *(kv.submat(0,iaddr-1,nr,nt)).transpose();
 	    if(t1==0 && t2==0) kRk+=kv.submat(0,iaddr-1,nr,nt)*rmat
 				 *(kv.submat(0,iaddr-1,nr,nt)).transpose();
 	  }
 	}
 	if(term[i_vc].classi() == 'P') {
 	  ainv();
 	  ainv_ia=hainv.ia();
 	  ainv_ja=hainv.ja();
 	  ainv_a=hainv.a();
 	  ainv_nrow=hainv.num_rows();
 	  for(ii=1;ii<=ainv_nrow;ii++) {
 	    ipos=ii;
 	    iaddr=startaddr+(ipos-1)*nt;
 	    je=ainv_ia[ii+1];
 	    for(jj=ainv_ia[ii];jj<je;jj++) {
 	      jpos=ainv_ja[jj];
	      double offdiag;
	      offdiag=2.;
	      if(ipos==jpos) offdiag=1.;
 	      jaddr=startaddr+(jpos-1)*nt;
 	      kGk+=ainv_a[jj]*offdiag*
		kv.submat(0,iaddr-1,nr,nt)*wmat
 		*(kv.submat(0,iaddr-1,nr,nt)).transpose();
 	      if(t1==0 && t2==0)  kRk+=kv.submat(0,iaddr-1,nr,nt)*rmat
 				    *(kv.submat(0,iaddr-1,nr,nt)).transpose();
 	    }
 	  }
 	}
 	Lambda[k_vc]=(kvk*kGk).trace();
 	if(t1 == 0 && t2 == 0)  addj+=(kvk*kRk).trace();
       }
     }
   }

   //   ****
   // Build Residual SSs and CP
   //   ****
   unsigned kvc;

   if(link_function){

     for(kvc=0;do_varcomp && kvc<res_nvc;kvc++,k_vc++){
       wmat.assign(0.0);
       rmat.assign(0.0);
       for(ii=0;ii<numrec;ii++) {
 	wmat+=observation[ii].H.transpose()*observation[ii].Resid_Sin_Mat[kvc]*
 	  observation[ii].H*observation[ii].weight;
 	rmat+=observation[ii].rve*observation[ii].weight;
 	//
 	// At this point I am not sure if it is better to adjust
 	// for prediction errors or leave it alone.
 	// For "most" of the cases I typically run accross it is
 	// probably not much of an issue.
 	//
	
       }
       Lambda[k_vc]=(static_cast<double>(nr) - addj)*(wmat*rmat.inv()).trace();

     }


   }
   else{
     rmat.assign(0.0);
     for(ii=0;ii<numrec;ii++) {
       rmat+=observation[ii].rve*observation[ii].weight;
     }
     for(t1=0;t1<numtrait;t1++) for(t2=t1;t2<numtrait;t2++,k_vc++) {
       wmat.assign(0.0);
       for(ii=0;ii<numrec;ii++) {
 	wmat+=SinMat[k_vc]*observation[ii].weight;
	
       }
       Lambda[k_vc]=(static_cast<double>(nr) - addj)*(wmat*rmat.inv()).trace();
     }
   }

   if(LR==2) glim();
   double dfe;


   // K'b-M
   for (i=0; i<nr; i++) kpb_m[i] = blupsol.inner_product(kpme[i]);
   if (M) for (i=0; i<nr; i++) kpb_m[i] -= M[i];
   double quq;
   if(LR!=1) {
     quq = kvk.quadratic(kpb_m,kpb_m);       // (K'b-m)'(K'(X'V-1X)-K)-1(K'b-m)
   }
   else{
     //    tmpsol=blupsol-((kvk*kpb_m)*kv); //tmpsol now contains the restricted solutions
     {
       //double *blupve;
       //blupve=blupsol.ve;
       //blupsol.ve=tmpsol.ve;
       glim(10,(double **)kpme,nr,nc,(double *)M);
       residual();
       quq= -like_val;
       //blupsol.ve=blupve;
       glim();
       residual();
       quq+=like_val;
       quq*=2.;
     }
   }
   double ssm,sse;
   unsigned rank_x;
   int degen = 0;
   int errcode = 0;
   double f_stat=0.0, prob=0.0,vart=0.0;
   double sigma_e = residual_var[0][0];
   int not_chisq = Info.num_rows();
   if(Asym) not_chisq=0;

   if(not_chisq) vart=Info.ginv0().quadratic(Lambda,Lambda);
   if (!not_chisq || vart <= EPSILON ) {
     f_stat = quq;

     prob = ChiSquare_cdf(f_stat,static_cast<double>(nr),0.0,errcode);
     f_stat /= static_cast<double>(nr);
     dfe=2000.;
   }
   else {

     dfe=(2.0*nr*nr)/vart;
     if(dfe < 1.) dfe=1.;
     if(dfe > 2000.) {
       dfe=2000.;
     }
     f_stat = (quq)/nr;
     // if(f_stat > 1000.) f_stat=1000.;
     prob = F_cdf(f_stat,static_cast<double>(nr),static_cast<double>(dfe),0.0,errcode);

     //   in univariate cases only: the estimated variance of estimator K'B
   }
   std::string raw_data_code;
   double p_value,var,tcal=0.0;
   int W = SESSION.output_precision;
   std::cout.precision(W);
   if(Print_Level < 0) Print_Flag=0;
   if (Print_Flag) {
     std::cout << "\n            RESULTS FROM CONTRAST(S)\n";
     std::cout << " ----------------------------------------------------------\n";
     std::cout << "  Contrast MME_addr    K_coef   Raw_data_code\n";
     std::cout << "  ---------------------------------------------\n";
   }
   for (i=0; i<nr; i++) {
     if (Print_Flag) {
       kp = kpme[i];
       for (j=0; j<nc; j++) {
 	if (kp[j] == 0.0) continue;
 	trait_effect_level(j,raw_data_code);
 	std::cout << std::setw(4) << i+1 <<std::setw(11) << j+1 << std::setw(13) << kp[j]
 	     << "   " << raw_data_code << "\n";
       }
     }
     kpd = kpdiff[i];
     if (result) {
       result[i][0] = kpb_m[i];
       result[i][3] = kpd;
     }
     var = kvkori[i][i];
     if (var <= 0.0) throw exception("GLMM::contrast(): you have probably found a bug!");
     var = std::sqrt(var);
     //    std::cout << var << "                              \n";
     //    std::cout << kpb_m[i] << "      xxxxxx              \n";
     tcal = fabs(kpb_m[i]/var);
     //  std::cout << tcal << "  "<< dfe << "                               \n";
     p_value = 2.0*(1.0-t_cdf(tcal,static_cast<double>(dfe),0.0));

     if (result) {
       result[i][1] = var;
       result[i][2] = p_value;
     }
     if (!Print_Flag) continue;
     if (kpd > 1.0e-5) {
       est = 0;
       std::cout << "            **** NON-ESTIMABLE ****\n\n";
     }
     else {
       if (kpd > 1.0e-10) {
 	std::cout << " ***WARNING***: it may or mayn't be estimable\n";
 	std::cout << "       max(k'*ginv(mme)*mme - k') = " << kpd << "\n";
       }
       std::cout << "          estimated value (K'b-M) = " << kpb_m[i]
 	   << " +- " << var << "\n";
       std::cout << "             Prob(|t| > " << tcal << ") = "
 	   <<  p_value << " (p_value)\n";
       if(nr == 1 && not_chisq)
 	std::cout << "             " << dfe << " (Error degrees of freedom)\n";
       if (i != nr-1) std::cout << "\n";
     }
   }
   if (Print_Flag || Print_Summary) {
     std::cout << " ----------------------------------------------------------\n";
     if ((nr > 1 || Print_Summary) && est) {
       std::cout << "   joint hypothesis test H: K'b = M\n";
       std::cout << "   Prob(F > " <<  f_stat << ") = "
 	   <<  1.0-prob << " (p_value)\n";
       if(not_chisq) std::cout << "  " << dfe << " (Error degrees of freedom)\n";
     }

   }
   if(Return_Stat) return f_stat;
   return 1.0-prob;   // p_value
 }




 //
 //
 //

 double GLMM::contrast(const Vector<double>& Kp,const double m)
 {
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////
	
   double retval, **kpme = new double *[1];
   kpme[0] = Kp.begin();
   unsigned nc = Kp.size();
   double* M = (double *)NULL;
   if (m != 0.0) {
     M = new double [1];
     M[0] = m;
   }
   retval = contrast((const double **)kpme,1,nc,M);
   delete [] kpme;
   if (M) delete [] M;
   return retval;
 }

 double GLMM::contrast(const doubleMatrix& Kp)
 {
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////
	
   return contrast((const double **)Kp.begin(), Kp.num_rows(), Kp.num_cols());
 }

 double GLMM::contrast(const doubleMatrix& Kp,const Vector<double>& M)
 {
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted 
   ///////////////////////////////////////////
	

   double retval = 0.0;
   int nr = Kp.num_rows();
   unsigned nc = Kp.num_cols();
   if (M.size() != nr ) throw exception(" GLMM::contrast: size is incomformable");
   Vector<double> fullsize_M(hmmesize);
   Matrix<double>fullsize_Kp(nr,hmmesize);
   for (unsigned i=0; i<nr; i++) {
     memcpy(fullsize_Kp[i],Kp[i], sizeof(double)*nc);
     }
     memcpy(fullsize_M.begin(),M.begin(), sizeof(double)*nc);
     retval=GLMM::contrast((const double **)fullsize_Kp.begin(),nr,hmmesize,fullsize_M.begin());
     return retval;
 }

 double GLMM::contrast(const std::string &termname,const doubleMatrix& Kp, const Vector<double>& M)
 {
   if (!(Kp.begin())) {
     warning("GLMM::contrast(kp): kp is null, thus it is ignored");
     return 0.0;
   }

   double retval = 0.0;
   if (!get_blupsol()) throw exception("GLMM::contrastranspose(): bad model");

   int k = term.index(termname,factor_struct);
   if (k < 0) {
     warning("GLMM::contrast(): no such term in the model");
     return retval;
   }

   unsigned nr = Kp.num_rows();
   unsigned nc = Kp.num_cols();
   if (M.size() != nr) {
     warning("GLMM::contrast(): unconformable size");
     return retval;
   }
   unsigned i,startaddr;
   int nt=nt_vec[k];
   if (nc <= term[k].nlevel()*nt) {
     startaddr = term[k].startaddr();
     Vector<double> fullsize_M(hmmesize);
     Matrix<double>fullsize_Kp(nr,hmmesize);
     for (i=0; i<nr; i++) {
       memcpy(&(fullsize_Kp[i][startaddr]),Kp[i], sizeof(double)*nc);
     }
     memcpy(&(fullsize_M[startaddr]),M.begin(), sizeof(double)*nc);
     retval=contrast((const double **)fullsize_Kp.begin(),nr,hmmesize,fullsize_M.begin());
   }
   else {
     warning("GLMM::contrast(%s,Kp): too many columns in Kp",termname.c_str());
   }
   return retval;
 }

 double GLMM::restricted_log_likelihood(int solve)
 {
   ////////////////////////////////////////////////////////////////////
   // it returns restricted log likelihood under this model with data D
   //  Important Notice:
   //  because  restricted_log_likelihood will be called sequentially
   //  under the same MME structure.
   //  hmmec is kept after restricted_log_likelihood(D) call.
   //  Call release_mme() to release hmmec from memory
   ////////////////////////////////////////////////////////////////////
		
   if(!link_function) {
		
     return(Model::restricted_log_likelihood());
   }
   int Need_PEV_orig=Need_PEV;
   Need_PEV=0;
   Need_Residual=1;
   setup_mme(&rellrhs);
   double rellhood = 0.0;

   SESSION.warning = 0;
   if(solve){
     hmmec.solve(blupsol,rellrhs,"ysmp1");
     //     Need_Residual=1;
     setup_mme(&rellrhs);
   }
   hmmec.factor();
   //hmmec.display_order();
   //rellhood=0.0;
   rellhood = hmmec.info_vec[0];
   //   SESSION.warning = 1;
   // std::cout << "\n relhood " << rellhood << " rank " <<hmmec.info_vec[1];
   SESSION.warning = 0;


   if(Use_Like)rellhood +=yry+like_adj;
   //std::cout << "\n yry " << yry << " like_adj " << like_adj;
   if(!Use_Like)rellhood += yry - blupsol.inner_product(rellrhs);
   unsigned i;
   int done_ped;
   done_ped=0;
   // for (i=0; i<npattern; i++) rellhood += double(kvec[i])*lnr0vec[i];
   for (i=0; i<numterm; i++) {
     if (term[i].classi() == 'R'){
       rellhood += lng0vec[i]*term[i].nlevel();
       //std::cout << "\n[" << i << "] n "<<term[i].nlevel() << " log|D|" << lng0vec[i];
     }
     else if (term[i].classi() =='P' && !done_ped) {
       // ainv();
       done_ped=1;
       //       build_CorrVar();
       rellhood += lng0vec[i]*double(popsize-numgroup);
       //       std::cout << "\n[" << i << "] n "<< (popsize-numgroup)<< " log|G|" << lng0vec[i] ;
       rellhood -= ldet*double(nt_vec[i]);
       //        std::cout << "|A| " << -double(nt_vec[i])*ldet <<"\n";;
     }
   }
   Need_PEV=Need_PEV_orig;
   return (-0.5*rellhood);
 }




 double GLMM::log_likelihood(int solve)
 {
   ////////////////////////////////////////////////////////////////////
   // it returns restricted log likelihood under this model with data D
   //  Important Notice:
   //  because  restricted_log_likelihood will be called sequentially
   //  under the same MME structure.
   //  hmmec is kept after restricted_log_likelihood(D) call.
   //  Call release_mme() to release hmmec from memory
   // Still need -.5ln|D^{-1}+Z'H'R^{-1}HZ|
   ////////////////////////////////////////////////////////////////////

   if(!link_function) {

     return(Model::restricted_log_likelihood());
   }

   Need_Residual=1;
   setup_mme(&rellrhs);
   double rellhood = 0.0;

   SESSION.warning = 0;
   if(solve){
     hmmec.solve(blupsol,rellrhs,"ysmp1");
     Need_Residual=1;
     setup_mme(&rellrhs);
   }
   //   hmmec.factor();
   //hmmec.display_order();
   rellhood=0.0;
   //rellhood = hmmec.info_vec[0];
   //   SESSION.warning = 1;
   //       std::cout << "\n relhood " << rellhood << " rank " <<hmmec.info_vec[1];
   SESSION.warning = 0;


   if(Use_Like)rellhood +=yry+like_adj;
   //std::cout << "\n yry " << yry << " like_adj " << like_adj;
   if(!Use_Like)rellhood += yry - blupsol.inner_product(rellrhs);
   unsigned i;
   int done_ped;
   done_ped=0;
   // for (i=0; i<npattern; i++) rellhood += double(kvec[i])*lnr0vec[i];
   for (i=0; i<numterm; i++) {
     if (term[i].classi() == 'R'){
       rellhood += lng0vec[i]*term[i].nlevel();
            //std::cout << "\n[" << i << "] n "<<term[i].nlevel() << " log|D|" << lng0vec[i];
     }
     else if (term[i].classi() =='P' && !done_ped) {
       done_ped=1;
       //build_CorrVar();
       rellhood += lng0vec[i]*double(popsize-numgroup);
            //std::cout << "\n[" << i << "] n "<< (popsize-numgroup)<< " log|G|" << CorrVar.logdet() << " CorrVar " << CorrVar;
       rellhood -= hainv.logdet()*double(act_numtrait);
           //std::cout << "|A| " << -double(act_numtrait)*hainv.logdet() <<"\n";;
     }
   }
   return (-0.5*rellhood);
 }

 #ifdef GLMM_QUAD
 extern int ysmp_solve(int n,int *iiv,int *jjv,double *aav,int *p,int *ip,int nsp,
 		      double *rsp,const double *b, double *z,double tol,int& irank,
 		      double& lgdet,int path1,int path2);

 double GLMM::quad(double *x, double *y, double *xsol, double *ysol,  int path){
   ////
   //
   // calculates x`Cy
   // path = 3 return x`Cy,  Du`y in ysol, and ux in xsol
   // path = 1 return Du`y in ysol
   // path = 2 given Du`y in ysol, return x`Cy and ux in xsol
   //
   ///
		
   double xCy,ldet;
   int path1,path2,ir;
   xCy=0;
   if(!hmmec.perm) hmmec.factorization();
   path1=0;
   if(path & 1){
     path2=13;
     ysmp_solve(hmmec.dim,hmmec.iiv,hmmec.jjv,hmmec.aav,
 	       hmmec.perm,hmmec.iperm,hmmec.nsp,hmmec.rsp,y,ysol,
 	       EPSILON,ir,ldet,path1,path2);
   }
   if(path & 2){
     path2=23;
     ysmp_solve(hmmec.dim,hmmec.iiv,hmmec.jjv,hmmec.aav,
 	       hmmec.perm,hmmec.iperm,hmmec.nsp,hmmec.rsp,x,xsol,
 	       EPSILON,ir,ldet,path1,path2);
     for(int i=0;i< hmmec.dim;i++) xCy += xsol[i]*ysol[i];
   }
   return(xCy);

 }
 #endif

  ////////////////////////////////////////////////////////////
  //KP function
  //March 1998, University of Nebraska
  //Copyright (C) 1998 Roger Collins
  /////////////////////////////////////////////////////////////

 doubleMatrix *KP( Data *D,doubleMatrix &SKM_ans,const std::string &lpl,const std::string &censor){

   // 0=censored 1=uncensored;
   doubleMatrix ret,*retpt;
   doubleMatrix SKM;
   int ny;
   ny=D->num_rows();
   doubleMatrix lpl_mat;
   lpl_mat.resize(ny,2);
   int lpl_col;
   if(-1 ==(lpl_col=D->field_index(lpl))) {
     warning("KP: Can't find lpl %s \n",lpl.c_str());
     return(NULL);
   }
   for(int i=0;i<ny;i++){
     lpl_mat[i][0]=D->datasheet[lpl_col][i].double_val();
   }
   int status_col;
   if(-1 ==(status_col=D->field_index(censor))) {
     warning("KP: Can't find status %s \n",censor.c_str());
     return(NULL);
   };
   for(int i=0;i<ny;i++){
     lpl_mat[i][1]=D->datasheet[status_col][i].double_val();
     // std::cout << i << " " << lpl_mat[i][0] << " " << lpl_mat[i][1] << endl;
   }
   lpl_mat.sortby(Matrix<double>::COLUMN,0);
   // std::cout << "\n  Input doubleMatrix col 1 = time  col 2: 0=censored 1=uncensored " << lpl_mat;
   int j;
   float numdead;
   float numcensored;
   float lagdead;
   float lagcensored;
   doubleMatrix skm_dead;
   skm_dead.resize(ny,3);
   j=0;
   lagdead=0.;
   lagcensored=0.;
   for(int i=0;i<ny;i++){
     if(i < ny-1 && lpl_mat[i][0] == lpl_mat[i+1][0]){
       if(lpl_mat[i][1]==1){
 	lagdead++;
       }
       if(lpl_mat[i][1]==0){
 	lagcensored++;
       }
     }
     else{
       if(lpl_mat[i][1]==1){
 	numdead=1.+lagdead;
 	numcensored=lagcensored;
 	lagdead=0.;
 	lagcensored=0.;
 	skm_dead[j][0]=lpl_mat[i][0];
 	skm_dead[j][1]=numdead;
 	skm_dead[j][2]=numcensored;
 	j++;
       }
       if(lpl_mat[i][1]==0.){
 	numdead=lagdead;
 	numcensored=1.+lagcensored;
 	lagdead=0.;
 	lagcensored=0.;
 	skm_dead[j][0]=lpl_mat[i][0];
 	skm_dead[j][1]=numdead;
 	skm_dead[j][2]=numcensored;
 	j++;
       }
     }
   }
   //std::cout << "\n SKM dead vector = " << skm_dead;
   skm_dead=skm_dead.submat(0,0,j,3);
   //std::cout << "\n SKM dead vector = " << skm_dead;


   //doubleMatrix SKM;
   SKM.resize(j+1,3);
   float tot_alive;
   tot_alive=ny;
   int k;
   k=0;
   for(int i=0;i<j+1;i++){
     if(i == 0){
       SKM[k][0]=0.;
       SKM[k][1]=1.;
       k++;
     }
     if(i > 0 && skm_dead[i-1][1] != 0.){
       SKM[k][0]=skm_dead[i-1][0];
       SKM[k][1]=SKM[k-1][1]*(1.-(skm_dead[i-1][1]/tot_alive));
       SKM[k-1][2]=SKM[k][1]/SKM[k-1][1];
       SKM[k-1][2]/=SKM[k][0]-SKM[k-1][0];
       tot_alive=tot_alive-(skm_dead[i-1][1]+skm_dead[i-1][2]);
       k++;
     }
     if(i > 0 && skm_dead[i-1][1] == 0.){
       tot_alive=tot_alive-(skm_dead[i-1][1]+skm_dead[i-1][2]);
     }
   }
   SKM[k-1][2]=SKM[k-2][2];
   //for(int ii=0;ii<k;ii++) cout << SKM[ii][0] << " " << SKM[ii][1] << " " << SKM[ii][2] <<endl;
   SKM_ans=SKM.submat(0,0,k,3);
   return(&SKM_ans);


 }




 void GLMM::add_G_Sand(void)
 {
   for (int t=0; t<numterm; t++) {
     if (term[t].classi() == 'P') {
       add_AgSand(t);
     }
     else if (term[t].classi() == 'R') {
       add_IgSand(t);
     }
   }
 }

 void GLMM::add_IgSand(int t)
 {
   unsigned k,i,ii,jj;
   int t1,t2;
   double **v = (double **)NULL;
   doubleMatrix *tmp = term[t].prior->var_matrix();
   if (!tmp) throw exception("Model::add_Ig():  you probably forgot to use M.variance(...)\n"
 	  "  to set variance for each random effect");
   v = tmp->begin();
   Matrix<double> rv(numtrait,numtrait);
   for (t1=0; t1<numtrait; t1++) {
     for(t2=0; t2<numtrait; t2++) rv[t1][t2]=v[t1][t2];
   }
   ginverse1(rv.begin(),numtrait,lng0vec[t],1,EPSILON);

   unsigned na = term[t].nlevel();
   unsigned startaddr = term[t].start + 1;
   for (k=0; k<na; k++) {
     i = startaddr + k*numtrait;
     for (t1=0; t1<numtrait; t1++) {
       ii = i + t1;
       for (t2=0; t2<numtrait; t2++) {
 	jj = i + t2;
 	if (jj>=ii) hSand.insert(ii,jj,rv[t1][t2]);
       }
     }
   }
 }



 void GLMM::add_Ag(int t,int ct)
 {
   unsigned t1, t2, k, i, j, ii, jj, iii, jjj;
   double dii,val,value;
   double **var = (double **)NULL;
   var = corrvar[t][ct].begin();
   Matrix<double> rvarg(numtrait,numtrait);
   for (i=0; i<numtrait;i++) for (j=0; j<numtrait; j++) rvarg[i][j]= var[i][j];
   ginverse1(rvarg.begin(),numtrait,lng0vec[t],1,EPSILON);

   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   unsigned startaddr = term[t].start + 1;
   unsigned startaddrct = term[ct].start + 1;
   unsigned asd[3];
   Individual *I;
   if(startaddr > startaddrct){
     jjj=startaddr;
     startaddr =startaddrct;
     startaddrct=jjj;
   }
   for (k=0; k<popsize; k++) {
     I = pop->member(k);
     asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
     dii = 4.0/(2.0 - I->father_inbcoef() - I->mother_inbcoef());

     for (i=0; i<3; i++) {
       if (asd[i] != 0) {
 	ii = startaddr + (asd[i]-1)*numtrait;
 	for (j=0; j<3; j++) {
 	  if (asd[j] != 0) {
 	    jj = startaddrct + (asd[j]-1)*numtrait;
 	    val = lambda[i][j]* dii;
 	    for (t1=0; t1<numtrait; t1++) {
 	      iii = ii + t1;
 	      for (t2=0; t2<numtrait; t2++) {
 		jjj = jj + t2;
 		if (jjj>=iii) {
 		  value = val*rvarg[t1][t2];
 		  hmmec.insert(iii,jjj,value);
 		}
 	      }
 	    }
 	  }
 	}
       }
     }
   }
 }


 void GLMM::add_AgSand(int t)
 {
   unsigned t1, t2, k, i, j, ii, jj, iii, jjj;
   double dii,val,value;
   double **var = (double **)NULL;
   doubleMatrix *tmp = term[t].prior->var_matrix();
   if (!tmp) throw exception("Model::add_Ag():  you probably forgot to use M.variance(...)\n"
 	  "  to set variance for each random effect");
   var = tmp->begin();
   Matrix<double> rvarg(numtrait,numtrait);
   for (i=0; i<numtrait;i++) for (j=0; j<numtrait; j++) rvarg[i][j]= var[i][j];
   ginverse1(rvarg.begin(),numtrait,lng0vec[t],1,EPSILON);

   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   unsigned startaddr = term[t].start + 1;
   unsigned asd[3];
   Individual *I;

   for (k=0; k<popsize; k++) {
     I = pop->member(k);
     asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
     dii = 4.0/(2.0 - I->father_inbcoef() - I->mother_inbcoef());

     for (i=0; i<3; i++) {
       if (asd[i] != 0) {
 	ii = startaddr + (asd[i]-1)*numtrait;
 	for (j=0; j<3; j++) {
 	  if (asd[j] != 0) {
 	    jj = startaddr + (asd[j]-1)*numtrait;
 	    val = lambda[i][j]* dii;
 	    for (t1=0; t1<numtrait; t1++) {
 	      iii = ii + t1;
 	      for (t2=0; t2<numtrait; t2++) {
 		jjj = jj + t2;
 		if (jjj>=iii) {
 		  value = val*rvarg[t1][t2];
 		  hSand.insert(iii,jjj,value);
 		}
 	      }
 	    }
 	  }
 	}
       }
     }
   }
 }

 doubleMatrix operator*(Vector<double> &a,Vector<double> &b)
 {
   int row,col;
   row=a.size();
   col=b.size();
   doubleMatrix temp(row,col);
   for(int i=0;i<row;i++) {
     for(int j=0;j<col;j++) {
       temp[i][j]=a[i]*b[j];
     }
   }
   return temp;
 }

 void GLMM::covariance_old(const std::string &termname,const std::string &termname2,const doubleMatrix& v)
   // **********************************************************************
   // * covariance is the method to get variance components for each random
   // * effect.
   // *********************************************************************
 {
   if (!factor_struct) throw exception("Model::variance(args): no equation(s) in the model");

   int k = term.index(termname,factor_struct);
   int k2 = term.index(termname2,factor_struct);
   if(k2<k) {
     int kt=k;
     k=k2;
     k2=kt;
   }
   if ( k < 0) {
     warning("GLMM::covariance(): %s: not in the model",termname.c_str());
     type = bad_model;
     return;
   }
   if ( k2< 0) {
     warning("GLMM::covariance(): %s: not in the model",termname2.c_str());
     type = bad_model;
     return;
   }
   if (term[k].classi() != 'P')   {
     warning("GLMM::covariance(): %s: no associated Pedigree",termname.c_str());
     type = bad_model;
     return;
   }          // factor classi will override
   if (term[k2].classi() != 'P')  {
     warning("GLMM::covariance(): %s: no associated Pedigree",termname2.c_str());
     type = bad_model;
     return;
   }                  // factor classi will override
   if(!corrmap[k]) {
     corrmap[k]=1;
     ncorr++;
   }
   if(!corrmap[k2]) {
     corrmap[k2]=1;
     ncorr++;
   }
   if (v.num_rows() == v.num_cols() && v.num_rows() == numtrait) {
       corrvar[k][k2]=v;
   }
   else {
     warning("Model::variance(): invalid variance matrix size");
     type = bad_model;
   }
 }


 void GLMM::build_CorrVar(void) {
   CorrVar.resize(act_numtrait,act_numtrait);
   int ipos,jpos,i,j;
   doubleMatrix *var;
   Matrix<int> posmat(numterm,act_numtrait);
   for(i=0,ipos=0;i<numterm;i++) {
    if(term[i].classi() == 'P' && corrmap[i]==0) {
      corrmap[i]=1;
      ncorr++;
    }
     if(corrmap[i]) {
       for(int t=0;t<numtrait;t++) {
 	posmat[i][t]=ipos++;
       }
     }
   }
   for(i=0;i<numterm;i++) {
     if(corrmap[i]) {
       var = term[i].prior->var_matrix();
       for(int t1=0;t1<numtrait;t1++) {
 	for(int t2=0;t2<numtrait;t2++) {
 	  CorrVar[posmat[i][t1]][posmat[i][t2]]=(*var)[t1][t2];
 	}
       }
       for(j=(i+1);j<numterm;j++) {
 	if(corrmap[j]) {
 	  var = &corrvar[i][j];
 	  for(int t1=0;t1<numtrait;t1++) {
 	    for(int t2=0;t2<numtrait;t2++) {
 	      CorrVar[posmat[i][t1]][posmat[j][t2]]=(*var)[t1][t2];
 	      CorrVar[posmat[j][t2]][posmat[i][t1]]=(*var)[t1][t2];
 	    }
 	  }
 	}
       }
     }
   }
 }

 void GLMM::add_G_1_old(void)
 {
   int done_ped=0;
    for (int t=0; t<numterm; t++) {
       if (term[t].classi() == 'P' && !done_ped) {
 	done_ped=1;
 	Model::add_Ag(t);
       }
       else if (term[t].classi() == 'R') {
          Model::add_Ig(t);
       }
    }
 }


 void GLMM::add_Ag_old(int t)
 {
    unsigned t1, t2, k, i, j, ii, jj, iii, jjj,icnt;
    double dii,val,value;
    double **var = (double **)NULL;
    Vector<int> tmap;
    build_CorrVar();
    tmap.resize(act_numtrait);
    icnt=0;
    for(i=t;i<numterm;i++){
      if(term[i].classi()=='P') {
	for(j=0;j<numtrait;j++,icnt++){
	  tmap[icnt]=i;
	}
      }
    }
    //    std::cout <<tmap;
    doubleMatrix *tmp = &CorrVar;
    if (!tmp) throw exception("Model::add_Ag():  you probably forgot to use M.variance(...)\n"
             "  to set variance for each random effect");
    var = tmp->begin();
    Matrix<double> rvarg(act_numtrait,act_numtrait);
    for (i=0; i<act_numtrait;i++) for (j=0; j<act_numtrait; j++) rvarg[i][j]= var[i][j];

    ginverse1(rvarg.begin(),act_numtrait,lng0vec[t],1,EPSILON);


    double f1,f2;
    Matrix<double> lambda(3,3);
    getlambda(lambda.begin(),3);
    unsigned startaddr = term[t].start + 1;
    unsigned asd[3];
    Individual *I;
    unsigned nanim=popsize-numgroup;
    for (k=0; k<nanim; k++) {
       I = pop->member(k);
       asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
       f1= I->father_inbcoef() ;
       f2= I->mother_inbcoef();
       if(asd[1]> nanim) f1=-1.;
       if(asd[2]> nanim) f2=-1.;
       dii = 4.0/(2.0 -f1-f2);

       for (i=0; i<3; i++) {
          if (asd[i] != 0) {
             ii = startaddr + (asd[i]-1)*act_numtrait;
             for (j=0; j<3; j++) {
                if (asd[j] != 0) {
                   jj = startaddr + (asd[j]-1)*act_numtrait;
                   val = lambda[i][j]* dii;
                   for (t1=0; t1<act_numtrait; t1++) {
                      iii = ii + t1;
                      for (t2=0; t2<act_numtrait; t2++) {
                         jjj = jj + t2;
                         if (jjj>=iii) {
                            value = val*rvarg[t1][t2];
			    //if(term[tmap[t1]].trait[t1%numtrait] && term[tmap[t2]].trait[t2%numtrait])
			      hmmec.insert(iii,jjj,value);
                         }
                      }
                   }
                }
             }
          }
       }
    }
 }



 void
 bound_pd (double *varnew, int numtrait, doubleMatrix & P)
 {
   int i, j, k, t1, t2;
   doubleMatrix varbound (numtrait, numtrait);
   Vector<double> eig, flg;
   flg.resize (numtrait);
   for (k = 0, i = 0; i < numtrait; i++)
     for (j = i; j < numtrait; j++, k++)
       {
 	varbound[i][j] = varnew[k];
 	varbound[j][i] = varnew[k];
       }
   P = varbound;

   eig = P.eigen ();

   for (t1 = 0; t1 < numtrait; t1++)
     {
       flg[t1] = 0;
       if (eig[t1] < 1.e-3)
 	{
 	  eig[t1] = .99e-3;
 	  flg[t1] = 1;
 	}
       if (eig[t1] > 1.e10)
 	{
 	  eig[t1] = 1.01e10;
 	  flg[t1] = 0;
 	}
     }
   varbound = P * diag(eig) * P.transpose();
   for (t1 = 0; t1 < numtrait; t1++)
     {
       if (flg[t1])
 	{
 	  for (t2 = 0; t2 < numtrait; t2++)
 	    {
 	      P[t1][t2] = 0;
 	    }
 	}
     }

   for (k = 0, t1 = 0; t1 < numtrait; t1++)
     for (t2 = t1; t2 < numtrait; t2++, k++)
       {
 	varnew[k] = varbound[t1][t2];

       }
 }

 void GLMM::info(std::ostream& stream)
 {
    if (type == bad_model) {
       warning("Model::info(stream): model is too bad");
       return;
    }
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
    if(aireml_called){
      stream << "AI REML ";
      if(aireml_called < 0) stream << "failed to converge\n";
      else stream << "converged\n";
      stream << "maximum log restricted likelihood = "
              <<  restricted_log_likelihood(0) <<"\n";
    }
    const ModelTerm *T;
    int i,done_ped,nt;
    doubleMatrix *var;
    done_ped=0;
    for (i=0; i<numterm; i++) {
       T = &(term(i));
       if (T->classi() == 'R' || (T->classi() == 'P' && !done_ped)) {
 	nt=nt_vec[i];
 	var=T->prior->var_matrix();
 	//if(T->classi() == 'P') {
 	//  nt=act_numtrait;
 	//  done_ped=1;
 	//  build_CorrVar();
 	//  var=&CorrVar;
 	//}
          stream << " variance for "
                 << factor_struct[T->factorindx[0]].name() << " = \n"
                 << (*var);
       }
    }
    stream << " residual variance =\n" << residual_var;

    if (data_prepared) {
       stream << "\n";
       stream << " MME dimension   : " << hmmesize << "\n";
       stream << " non-zeros in MME: " << hmmec.nz() << "\n\n";

       stream << "       basic statistics for dependent variables\n";
       stream << "  -----------------------------------------------------\n";
       stream << "      trait-name         n          mean         std\n";
       unsigned n;
       for (i=0; i<numtrait; i++) {
 	// std::cout <<"\n i= " << i << "\n";
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



 void GLMM::save(const std::string &fname, const int io_mode)
 {
    if (type == bad_model) {
       warning("Model::save(): model is too bad, nothing is saved");
       return;
    }
    if((numterm+ntermGdist) == 0) return;
    if (blupsol.size() != hmmesize + ntermGdist*popsize) return;
    if (!data) {
       warning("Model::save(): no data has been fitted");
       return;
    }
    std::ofstream ofs;
    ofs.open(fname.c_str(),(OpenModeType)io_mode);
    if (!ofs) throw exception("Model::save(): cann't open or already exit");
    ofs.setf(std::ios::unitbuf | std::ios::showpoint);
    ofs.precision(SESSION.output_precision);
    int W = SESSION.output_precision + 6;        // 6 = +.e+00
    char S = ' ';                               // one blank space

    if (!data->in_memory()) data->input_datasheet();
    this->info(ofs);

    double *ve = blupsol.begin();

    char* tempstr;
    int nord, *tempval;
    unsigned *tempuns;
    unsigned i,j,t,k,ii,kk,jj,nlevel,total_numterm;
    char CL;

    ofs << "\n            BLUP (BLUE, Posterior_Mean) \n";
    ofs << "    ----------------------------------------------- \n\n";
    total_numterm = numterm+ntermGdist;
    for (i=0; i<total_numterm; i++) {
       ofs << S << std::setw(W)<< "Term";
       for (t=0; t<numtrait; t++) {
          ofs << S << std::setw(W-2) << "Trait" << S << std::setw(2) << t+1;
       }
       ofs <<"\n";

       ofs << S << std::setw(W) << factor_struct[term[i].factorindx[0]].name();
       for (t=1; t<term[i].order(); t++) {
          ofs << "*" << factor_struct[term[i].factorindx[t]].name();
       }
       for (t=0; t<numtrait; t++) ofs <<S <<std::setw(W)<<trait_struct[t]->name();
       ofs << "\n";
       nlevel = term[i].nlevel();
       k = term[i].startaddr();
       nord = term[i].order();
       if (nord == 1) {      // single_factor_term
          kk = factor_struct[term[i].factorindx[0]].index();
          CL = data->datasheet[kk].classi();
          switch(CL) {
             case 'F':
                for (j=0; j<nlevel; j++) {
                   tempval = (int *)(idlist[kk]->find(j+1));
                   ofs << S << std::setw(W) << *tempval;
                   for (t=0; t<numtrait; t++,k++) {
                      ofs << S << std::setw(W);
                      if (term[i].trait[t]) {
                         ofs << ve[k];
                      }
                      else {
                         ofs << S;
                      }
                   }
		   k+=(nt_vec[base_effect[i]]-numtrait);
                   ofs << "\n";
                }
                break;
             case 'U':
                for (j=0; j<nlevel; j++) {
                   tempstr = (char *)(idlist[kk]->find(j+1));
                   ofs << S << std::setw(W) << tempstr;
                   for (t=0; t<numtrait; t++,k++) {
                      ofs << S << std::setw(W);
                      if (term[i].trait[t]) {
                         ofs << ve[k];
                      }
                      else {
                         ofs << S;
                      }
                   }
                   ofs << "\n";
                }
                break;
             case 'P':
                for (j=0; j<nlevel; j++) {
                   ofs << S;
                   if (pop->maxnamelen > 0) {
                      ofs << std::setw(W) << (char *)pop->ind_name(j+1);
                   }
                   else {
                      ofs << std::setw(W) << j+1;
                   }
                   for (t=0; t<numtrait; t++,k++) {
                      ofs << S;
                      if (term[i].trait[t]) {
                         ofs << std::setw(W) << ve[k];
                      }
                      else {
                         ofs << std::setw(W) << S;
                      }
                   }
 		  k+=(nt_vec[base_effect[i]]-numtrait);
                   ofs << "\n";
                }
                break;
             case 'C':
                ofs << S << std::setw(W) << data->datasheet[kk].name();
                for (t=0; t<numtrait; t++,k++) {
                   ofs << S << std::setw(W);
                   if (term[i].trait[t]) {
                      ofs << ve[k];
                   }
                   else {
                      ofs << S;
                   }
                }
                ofs << "\n";
                break;
             case 'I':
                ofs << S << std::setw(W) << "1";
                for (t=0; t<numtrait; t++,k++) {
                   ofs << S << std::setw(W);
                   if (term[i].trait[t]) {
                      ofs << ve[k];
                   }
                   else {
                      ofs << S;
                   }
                }
                ofs << "\n";
                break;
             default:
                warning("Model.save(): unknown column type: %c",CL);
                break;
          }
          if (i+1 < total_numterm) ofs << "\n";
       }
       else {     // there are interactions
          for (j=0; j<nlevel; j++) {
             tempuns = (unsigned *)(xact_htable[i].find(j+1));
             kk = factor_struct[term[i].factorindx[0]].index();
             CL = data->datasheet[kk].classi();
             jj = tempuns[0];

             if (CL == 'F') {
                tempval = (int *)(idlist[kk]->find(jj));
                ofs << S << std::setw(W) << *tempval;
             }
             else if (CL =='U') {
                tempstr = (char *)(idlist[kk]->find(jj));
                ofs << S << std::setw(W) << tempstr;
             }
             else if ( CL =='C') {
                ofs << S << std::setw(W) << data->datasheet[kk].name();
             }
             else if (CL =='P') {
                ofs << S << std::setw(W) << (char *)pop->ind_name(jj);
             }
             else if (CL =='I') {
                ofs << S << std::setw(W) << "1";
             }
             else {
                warning("Model::save(): unknown column type: %c",CL);
                break;
             }

             for (ii=1; ii<nord; ii++) {
                kk = factor_struct[term[i].factorindx[ii]].index();
                CL = data->datasheet[kk].classi();
                jj = tempuns[ii];

                if (CL == 'F') {
                   tempval = (int *)(idlist[kk]->find(jj));
                   ofs << "*" << *tempval;
                }
                else if (CL =='U') {
                   tempstr = (char *)(idlist[kk]->find(jj));
                   ofs << "*" << tempstr;
                }
                else if (CL =='C' || CL =='I') {
                   ofs << "*" << data->datasheet[kk].name();
                }
                else if (CL =='P') {
                   ofs << "*" << (char *)pop->ind_name(jj);
                }
                else {
                  warning("Model::save(): unknown column type: %c",CL);
                }
             }
             for (t=0; t<numtrait; t++,k++) {
                ofs << S << std::setw(W);
                if (term[i].trait[t]) ofs << ve[k];
                else ofs << S;
             }
	     //if(CL == 'P' || CL=='F') 
	      k+=(nt_vec[base_effect[i]]-numtrait);
             ofs << "\n";
          }
          if (i+1 < total_numterm) ofs << "\n";
       }
    }
    ofs.close();
    if(!data_static) data->release_datasheet();
 }

 static int squeeze(unsigned len,int *ia,int *ja,double *a)
 /*******************************************************************
 	squeeze out zeros from (ia,ja,a)	Tianlin Wang at UIUC
 	ia(1...len),jjv(1..len),a(1..len)
 	return the actural number of non-zeros in a
         "actually this moves the nonzeros to the begining"
 *******************************************************************/
 {
    unsigned lmar = 1;
    unsigned rmar = len;
    for (;;) {
       if (lmar >= rmar) {
          if (ia[lmar] == 0 && rmar == lmar) return(rmar-1);
          else return(rmar);
       }
       if ( ia[lmar] != 0) {
          lmar++;
       }
       else if (ia[rmar] == 0) {
          rmar--;
       }
       else {
          ia[lmar] = ia[rmar];
          ja[lmar] = ja[rmar];
          a[lmar]  = a[rmar];
          lmar++;
          rmar--;
       }
    }
 }

 static void hsort_ija(unsigned n,int *ia,int *ja,double *a)
 /******************************************************************
 	this is a Heap sort, sorting (ia,ja) ascendingly
 	ia(1...n), ja(1...n) ,a(1..n)
 	where n is the actual number of non-zeros in A
 **********************************************************************/
 {
    unsigned iia, jja,i,j,l,ir;
    double aa;
    if (n == 1) return;
    l = (n >> 1)+1;		/* l = n/2+1; */
    ir=n;
    for (;;) {
       if (l > 1) {
          l = l-1;  iia = ia[l];  jja = ja[l];  aa = a[l];

       }
       else {
          iia = ia[ir];  jja = ja[ir];  aa = a[ir];
          ia[ir] = ia[1]; ja[ir] = ja[1];  a[ir] = a[1];
          ir--;
          if(ir == 1) {
             ia[1] = iia;  ja[1] = jja;  a[1] = aa;
             return;
          }
       }
       i=l;
       j=l+l;
       while ( j <= ir) {
          if (j < ir) {
             if (ia[j]<ia[j+1]  || (ia[j]==ia[j+1] && ja[j] < ja[j+1])) j++;
          }
          if (iia < ia[j] || (iia == ia[j] && jja < ja[j])) {
             ia[i] = ia[j];  ja[i] = ja[j];  a[i] = a[j];
             i=j;
             j=j+j;
          }
          else {
             j=ir+1;
          }
       }
       ia[i] = iia;  ja[i] = jja;  a[i] = aa;
    }
 }

 #ifdef DO_THREADS
 #include "pt.h"

 struct thread_arg {
   SparseMatrix *A,*Ainv;
   unsigned *rowdone,*order,hmmesize;
   int *ia,*ja,*lhs_ia,*lhs_ja;
   pthread_mutex_t *mdone_mutex;
 };


 int nextrow(int *lhs_ia,int *lhs_ja,unsigned *rowdone,pthread_mutex_t *mdone_mutex,int hmmesize, int i,unsigned *order,SparseMatrix *hInv,Vector<double> &InvCol) {
   int res;

   pt_mutex_lock(mdone_mutex,"mutex lck");     /* lock the mutex      */
   if(i) {
     for(int jj=lhs_ia[i];jj<lhs_ia[i+1];jj++) {
       // std::cout << i<< " " << jj << " "<< lhs_ja[jj] << "\n";
       hInv->insert(order[i-1],order[lhs_ja[jj]-1],InvCol[order[lhs_ja[jj]-1]-1]);
     }
   }
   *rowdone=(*rowdone)+1;
   res=*rowdone;                                 /* get & incr row cntr */
   pt_mutex_unlock(mdone_mutex,"mutex unlck"); /* unlock the mutex    */
   return(res<=hmmesize ? res : -1);                /* return -1 if done   */

 }


 pt_addr_t thread_code(pt_arg_t *t_arg)
 {
   thread_arg *arg=(thread_arg *)t_arg->data;

   //  std::cout << "thread_code " << arg <<  "\n";
   int n,p;                              /* n and p counters         */
   int m;                                /* row we're working on     */
   double sum;                           /* temporary var            */
   int *lhs_ia=arg->lhs_ia;
   int *lhs_ja=arg->lhs_ja;
   int *ia=arg->ia;
   int *ja=arg->ja;
   unsigned *rowdone=arg->rowdone;
   m=0;
   pthread_mutex_t *mdone_mutex=arg->mdone_mutex;
   unsigned *order=arg->order;
   SparseMatrix *hInv=arg->Ainv;
   SparseMatrix *hmmec=arg->A;
   unsigned hmmesize=arg->hmmesize;
   Vector<double> InvCol,wrkVec;
   wrkVec.resize(hmmesize);
   InvCol.resize(hmmesize);


   while ((m=nextrow(lhs_ia,lhs_ja,rowdone,mdone_mutex,hmmesize,m,order,hInv,InvCol))!=-1) {           /* m=row to do              */
     wrkVec(order[m-1])=1;
     hmmec->solvrow(InvCol.begin(),wrkVec.begin(),m);
     wrkVec(order[m-1])=0;
   }
   return(NULL);
 }
 #endif


 void GLMM::build_hInv()
 {

   unsigned row=0;
   Vector<double> InvCol,wrkVec;
   int *lhs_ia,*lhs_ja,*lhs_sh,nz;
   double *lhs_a;
   unsigned hsize=hmmesize;
   InvCol.resize(hmmesize);
   wrkVec.resize(hmmesize);
   hInv.resize(hmmesize,max_nz);
#ifdef DEBUG
     doubleMatrix LHS;
     LHS=hmmec.dense();
     std::cout << "\nbuild_hInv LHS\n"<<LHS <<"\n";
 #endif
   if (!hmmec.factor_done) hmmec.factor();
   unsigned int *order=hmmec.order.begin();

 #ifdef DO_THREADS
   pthread_mutex_t mdone_mutex;
   pt_mutex_init(&mdone_mutex,"init mutex"); /* initialize the mutex */

   thread_arg mult_arg;
   mult_arg.mdone_mutex=&mdone_mutex;
   mult_arg.A=&hmmec;
   mult_arg.Ainv=&hInv;
   mult_arg.hmmesize=hmmesize;
   mult_arg.rowdone=&row;
   mult_arg.order=hmmec.order.begin();
 #endif
   int *ia   = hmmec.iiv-1;
   int *ja   = hmmec.jjv-1;
   double *a = hmmec.aav-1;
   unsigned nrow   = 0;		
   unsigned oldrow = 0;
   nz=squeeze(hmmec.hsize,ia,ja,a);
   for(int i=1;i<=nz;i++) {
     ia[i]=hmmec.inv_order[ia[i]-1];
     ja[i]=hmmec.inv_order[ja[i]-1];
     if(ja[i] < ia[i]) {
       int itmp=ia[i];
       ia[i]=ja[i];
       ja[i]=itmp;
     }
   }
   hsort_ija(nz,ia,ja,a);
   for (int k=1; k<=nz ; k++) {
     if (ia[k] != oldrow) {
       nrow++;
       ia[nrow] = k;
       oldrow = ia[k];    
     }
   }
   ia[hmmec.dim+1] = nz+1;

   lhs_ia=hmmec.ia();
   lhs_ja=hmmec.ja();
   lhs_a=hmmec.a();
 #ifdef DO_THREADS
   mult_arg.ia=ia;
   mult_arg.ja=ja;
   mult_arg.lhs_ja=lhs_ja;
   mult_arg.lhs_ia=lhs_ia;
   int NTHREADS;
   char *MP_NUM;
   NTHREADS=1;
   MP_NUM=getenv("MP_NUM_THREADS");
   if(MP_NUM) NTHREADS=atol(MP_NUM);
   if(NTHREADS >1) {
     NTHREADS--;
     std::cout << "Nthreads "<< NTHREADS << "\n";
     //thread_code(&mult_arg[0]);
     //  std::cout << &mult_arg << "\n";
 //TW    pt_fork(NTHREADS,thread_code,(pt_addr_t) &mult_arg,NULL); /* run code in parallel */
    
    
 //TW    pt_mutex_destroy(&mdone_mutex,"del mutex");
     //  hInv.close();
   }
   else{
 #endif
     Vector<double> solwrk,temp_rhswrk;
     Vector<unsigned int> temp_col,temp_row; 
     int m,jj;
#pragma omp parallel private(jj,m,wrkVec,solwrk,temp_rhswrk,InvCol,temp_col,temp_row)
     {
#pragma omp critical
       {
       wrkVec.resize(hmmesize);
       InvCol.resize(hmmesize); 
       solwrk.resize(hmmesize);
       temp_rhswrk.resize(hmmesize);
       temp_row.resize(hmmesize);
       temp_col.resize(hmmesize);
       }
#pragma omp for schedule(dynamic,50) nowait
       for(m=1;m<=hmmesize;m++) {
	 wrkVec(order[m-1])=1;
	 hmmec.solvrow(InvCol.begin(),wrkVec.begin(),m,solwrk,temp_rhswrk,temp_row.begin(),temp_col.begin());
	 wrkVec(order[m-1])=0;
	 {
	   for(jj=lhs_ia[m];jj<lhs_ia[m+1];jj++) {
	     //	     hInv.insert(order[m-1],order[lhs_ja[jj]-1],InvCol[order[lhs_ja[jj]-1]-1]);
	     //	     hInv.insert(order[m-1],order[lhs_ja[jj]-1],solwrk[order[lhs_ja[jj]-1]-1]);
	     hInv.insert(order[m-1],order[lhs_ja[jj]-1],InvCol[lhs_ja[jj]-1]);
	   }
	 }
       }
     }
 #ifdef DO_THREADS
     pt_mutex_destroy(&mdone_mutex,"del mutex");
   }
 #endif
 }

doubleMatrix GLMM::getVarEstimates(const std::string &termname)
{
  // Authors: Liviu R. Totir and Rohan L. Fernando (March, 2004) 
  // Contributors: 
  doubleMatrix* retval = 0;
  if (type == bad_model) {
    warning("Model::info(stream): model is too bad");
    return *retval;
  }
  int k = term.index(termname,factor_struct);
  if (k < 0) throw exception("GLMM::getVarEstimates: no such term in the model");
  const ModelTerm *T;
  T = &(term(k));
  if (T->classi() == 'R' || (T->classi() == 'P')) {
    retval=T->prior->var_matrix();
  }
  return *retval;
}

} ///////// end of namespace matvec


