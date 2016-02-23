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



#define EPSILON  SESSION.epsilon
   //#define act_numtrait numtrait // Temporary fix

 doubleMatrix *P_Mat=NULL;
 ;
 ;
 ;

 extern void getlambda(double **lambda, const int n);
 
 void normal_loginClasses(matvec::Observe &observe);


 Observe& GLMM::new_obs()
 {
   unsigned i,j;
   if(observation) release_obs();
   observation=new Observe [numrec+1] ;
   observation[numrec].rec=-1;
   observation[numrec].numtrait=numtrait;
   observation[numrec].param=param;
   observation[numrec].data=data;
   observation[numrec].resid_var=&residual_var;
   if(link_function) {
     observation[numrec].nvc=res_nvc;
     observation[numrec].varcomp=&varcomp;
     (*link_function)(observation[numrec]);
     if(res_nvc)observation[numrec].Adj_Cov.identity(res_nvc,res_nvc);
   }
   for(i=0;i<numrec;i++) {
     observation[i].rec=i;
     observation[i].need_data=1;
     observation[i].pattern.reserve(numtrait);
     observation[i].numtrait=numtrait;
     observation[i].residual.resize(numtrait);
     observation[i].estimate.resize(numtrait);
     if(Initial.begin()) observation[i].estimate=Initial;
     observation[i].y.resize(numtrait);
     observation[i].Resid_Var.identity(numtrait);
     observation[i].H.identity(numtrait);
     observation[i].Hinv.identity(numtrait);
     observation[i].pev.resize(numtrait,numtrait);
     observation[i].weight=1;
     observation[i].pev_current=0;
     observation[i].rve.resize(numtrait,numtrait);
     //    observation[i].rve_sand.zeros(numtrait,numtrait);
     observation[i].rpy.resize(numtrait);
     observation[i].nvc=res_nvc;
     observation[i].data=data;
     observation[i].param=param;

     observation[i].resid_var=&residual_var;
     if(res_nvc) {
       observation[i].Resid_Sin_Mat=new doubleMatrix [res_nvc];
       observation[i].varcomp=&varcomp;
     }
     for(j=0;j<res_nvc;j++) {
       observation[i].Resid_Sin_Mat[j].resize(numtrait,numtrait);
     }

   }
   return *observation;
 }




 void GLMM::release_obs(void){
   unsigned i,j;

   if(observation) {
     delete [] observation;
     observation=(Observe *)NULL;
   }
 }

 void GLMM::residual(int get_pev)
 {
   unsigned i,j,k,l,rec,t1,t2,ii,jj;
   double xval;
   //GLMMModelTerm *localterm;
   UniformDist U;

   double *k_new_vec;
   k_new_vec=new double [hmmesize];
   double *k_new_vec2;
   k_new_vec2=new double [hmmesize];
   double *svec;
   svec=new double [hmmesize];
   double *xtmp;
   xtmp=new double [hmmesize];
   double *ytmp;
   ytmp=new double [hmmesize];
   double *koff,*koff2;
   koff=k_new_vec-1;
   koff2=k_new_vec2-1;
   double *soff;
   soff=svec-1;
   like_val=0;
   doubleMatrix var;
   int startaddr,nlevels,iaddr,*ainv_ia,*ainv_ja,ipos,ainv_nrow;
   int ped_cnt=0;
   int je,jpos,jaddr;
   double *ainv_a;

 #ifdef DO_HINV
   if(Need_PEV)
     build_hInv();
 #endif
     like_adj=0;
   for(i=0;i<numterm;i++){
     if(term[i].classi() == 'P' || term[i].classi() == 'R'){
       unsigned nanim=popsize-numgroup;
       startaddr=term[i].startaddr()+1;
       var=*term[i].prior->var_matrix();
       var.ginv1();
       int nt=nt_vec[i];
       if(term[i].classi() == 'R') {
 	nlevels=term[i].nlevel();
 	for(ii=1;ii<=nlevels;ii++) {
 	  iaddr=startaddr+(ii-1)*nt;
 	  like_val+=
 	    -.5*var.quadratic(blupsol.subvec(iaddr-1,nt),
 				blupsol.subvec(iaddr-1,nt));
 	  if(Use_Like)like_adj+=
			var.quadratic(blupsol.subvec(iaddr-1,nt),
 				blupsol.subvec(iaddr-1,nt));
 	}
       }
       if(term[i].classi() == 'P') {;
				    // 	build_CorrVar();
	doubleMatrix Varinv=term[i].prior->var_matrix()->ginv0();
	// 	corrmap[i]=1;
	// 	if(!ncorr) ncorr=1;
	//	ped_cnt++;
 	double mult_corr;
 	unsigned startaddr2;
	int jcorr=i;
 	//for(int jcorr=i;jcorr<=i;jcorr++) 
	  {
	    //if(corrmap[jcorr]) 
	    {
	    mult_corr=1.;
 	    //if(jcorr!=i && corrmap[jcorr]) {
	    // mult_corr=2.;
	    // var=corrvar[i][jcorr];
 	    //}
 	    startaddr2=term[jcorr].startaddr()+1;
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
 		jaddr=startaddr2+(jpos-1)*nt;
 		like_val+=-.5*ainv_a[jj]*offdiag*
 		  mult_corr*Varinv.quadratic(blupsol.subvec(jaddr-1,nt),
 						blupsol.subvec(iaddr-1,nt));
		if(i==jcorr &&Use_Like){
 		like_adj+=ainv_a[jj]*offdiag*
 		  mult_corr*Varinv.quadratic(blupsol.subvec(jaddr-1,nt),
 						blupsol.subvec(iaddr-1,nt));
		}
 	      }
 	    }
 	  }
 	}
       }
     }
   }

   if(!link_function) inverse_residual_var();
   if(!observation) new_obs();
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);
   // std::cout << "numrec = " << numrec << "\n";
   // std::cout << "numobs = " << numobs << "\n";


   double pselect,scalvar;
   int numrem,numneed,numsel;
   int selected;
   numneed=numrec;
   numrem=numrec;
   numsel=0;
   if(AI_sample && AI_sample < numrec) numneed=AI_sample;
   scalvar=static_cast<double>(numrem)/static_cast<double>(numneed);


   for(rec=0;rec<numrec;rec++) {

     selected=0;
     pselect=static_cast<double>(numneed)/static_cast<double>(numrem);
     if(U.sample() <= pselect) selected=1;
     numrem--;
     numsel+=selected;
     numneed-=selected;

     input_pos_val(modelfile);
     observation[rec].pattern = pattern[pos_term[numterm]];
     //std::cout << "\n rec " << rec << " " << pos_term[numterm] << " " << pattern[pos_term[numterm]];

     if(get_pev) observation[rec].pev.resize(numtrait,numtrait);
     for(t1=0;t1< numtrait;t1++) {
       memset(k_new_vec,'\0',sizeof(double)*hmmesize);
       //      observation[rec].residual.ve[t1]=trait_vec[t1];
       if(Update_estimate)
 	observation[rec].estimate[t1]=0;
       if(weightname.size())
 	observation[rec].weight = xval_term[numterm];
       for(i=0;i<numterm;i++) {
 	//	localterm=(GLMMModelTerm *)&term[i];
 	ii=term[i].addr[t1];
 	if(ii){
   	  xval=xval_term[i];
	
 	  //  observation[rec].residual[t1]-=xval*(blupsol(ii));
	
 	  if(Update_estimate)
 	    observation[rec].estimate[t1]+=xval*(blupsol(ii));
 	  koff[ii]=xval;
 	}
 #ifdef DO_HINV
 	if(get_pev) {
 	  for(t2=0;t2<numtrait;t2++) {
 	    for(j=0;j<numterm;j++) {
 	      int jj=term[j].addr[t2];
 	      if(jj){
 		(observation[rec].pev)[t1][t2] += xval_term[i]*xval_term[j]*hInv.getaij(ii,jj);
 	      }
 	    }
 	  }
 	  observation[rec].pev_current=1;
 	}
 #endif
	
       }
 #ifndef DO_HINV
       if(get_pev ) {
 	if(selected) {
 	  hmmec.solve(svec,k_new_vec,"ysmp");
 	  //quad(k_new_vec2,k_new_vec,xtmp,ytmp,1);
 	  for(t2=0;t2<numtrait;t2++) {
 	    memset(k_new_vec2,'\0',sizeof(double)*hmmesize);
	
 	    for(i=0;i<numterm;i++) {
 	      ii=term[i].addr[t2];
 	      if(ii){
 		koff2[ii]=xval_term[i];
 		(observation[rec].pev)[t1][t2] += xval_term[i]*soff[ii];
 	      }
 	    }
 	    //observation[rec].pev.me[t1][t2]=quad(k_new_vec2,k_new_vec,xtmp,ytmp,2);
 	  }
 	  observation[rec].pev_current=1;
 	}
 	else{
 	  observation[rec].pev_current=0;
 	}
       }
 #endif

     }
     if(link_function){
       link_function(observation[rec]);
       like_val+=observation[rec].like_val*observation[rec].weight;
       if(get_pev &&
 	  observation[rec].pev_current)
 	observation[rec].pev=
 	  observation[rec].H*
 	  observation[rec].pev*observation[rec].H.transpose();
     }
     else{
       observation[rec].Resid_Var.assign(0.0);
       for(t1=0;t1<numtrait;t1++) {
 	if(observation[rec].pattern[t1] != '0') {
 	  for(t2=0;t2<numtrait;t2++) {
 	    if(observation[rec].pattern[t2] != '0') {
 	      observation[rec].Resid_Var[t1][t2]=residual_var[t1][t2];
 	      observation[rec].rve[t1][t2]=rve[pos_term[numterm]][t1][t2];
 	    }
 	  }
 	}
       }

       for(t1=0;t1<numtrait;t1++) {
 	observation[rec].rpy[t1]=0.;
 	for(t2=0;t2<numtrait;t2++) {
 	  observation[rec].rpy[t1]+=
 	    observation[rec].rve[t1][t2]*trait_vec[t2];
 	}
       }
     }
   }
   delete [] ytmp;
   delete [] xtmp;
   delete [] svec;
   delete [] k_new_vec2;
   delete [] k_new_vec;
   //  modelfile.close();
   Need_Residual=0;
 }



 void GLMM::new_SinMat(void)
 {

   int i;
   int nvc;
   nvc=TotalNvc();
   if(SinMat) release_SinMat();
   SinMat = new doubleMatrix [nvc];
   for(i=0;i<nvc;i++) {
     
     SinMat[i].resize(numtrait,numtrait);
   }
 }



 void GLMM::release_SinMat(void)
 {
   if(SinMat) {
     int i,nvc;
     nvc=TotalNvc();

     for(i=0;i<nvc;i++) SinMat[i].resize(0,0);
     delete [] SinMat;
     SinMat=NULL;
   }

 }

 void GLMM::Build_SinMat()
 {
   int i,k,t1,t2,raneff;
#ifdef DO_CHOL
   doubleMatrix L;
#endif
#ifdef DO_LOG
   doubleMatrix P,L;
   Vector<double> D;
#endif
   if(!SinMat)new_SinMat();
   if(!varinv)varinv = new doubleMatrix [numterm];
   if(!Var)Var = new doubleMatrix [numterm];
   int done_ped=0;
   doubleMatrix var,part;
   part.resize(numtrait,numtrait);
   kvec.resize(numterm);
   int kk=0;
   for(k=0,i=0,raneff=0;i<numterm;i++) {
     if((term[i].classi() == 'P'&& !done_ped) || term[i].classi() == 'R'){
       int nt=nt_vec[i];
       kvec[i]=k;
       if(var_link[i]!=i) {
	 kvec[i]=kvec[var_link[i]];
       }
       else{
	 raneff++;
	 kk+=(nt*(nt+1))/2;
       }
       k=kvec[i];
       Var[i]=*term[var_link[i]].prior->var_matrix();
       
       
#ifdef DO_CHOL
       L=Var[i];
       L.chol();
       L=L.transpose();
#endif
#ifdef DO_LOG
       P=Var[i];
       D=P.eigen();
       double dmax=D.max();
       if(dmax < 1) dmax=1;
       double dmin=10.*dmax*SESSION.epsilon;
       for(int ii=0;ii<nt;ii++) {
	 if(D[ii] < dmin) D[ii]=.5*dmin;
       }
       D=log(D);
#endif
       var=Var[i].ginv0();
       varinv[i]=var;
#ifdef DO_PMAT
       var*=P_Mat[raneff]*P_Mat[raneff].transpose();
#endif
       for(t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,k++) {
	   part.resize(nt,nt);
	   //	std::cout << "PART\n";
	   part[t1][t2]=1.;
	   part[t2][t1]=1.;
#ifdef DO_CHOL
	   part[t1][t2]=0.;
	   part[t2][t1]=0.;
	   for(int j=0;j<nt;j++) {
	     part[t2][j]+=L[t1][j];
	     if(j != t2) part[j][t2]+=L[t1][j];
	   }
#endif
#ifdef DO_LOG
	   part[t1][t2]=0.;
	   part[t2][t1]=0.;
	   //std::cout << "DO LOG\n";
	   double D_avg,D_delt;
	   D_avg=.5*(D[t1]+D[t2]);
	   D_delt=.5*(D[t1]-D[t2]);
	   if(fabs(D_delt) > 1.e-8) {
	     
	     // Sinh(x)=(e^x-e^{-x})/2
	     part[t1][t2]=std::exp(D_avg)*sinh(D_delt)/(D_delt);
	     part[t2][t1]=std::exp(D_avg)*sinh(D_delt)/(D_delt);
	   }
	   else{
	     part[t1][t2]=std::exp(D_avg);
	     part[t2][t1]=std::exp(D_avg);
	   }
	 
	   /**/
	   if(D[t1]<-11 && D[t2]<-11) {
	     part[t2][t1]=0;
	     part[t1][t2]=0;
	   }
           /**/	 
	 
       
   
	 //std::cout << "part" << part;
	 part=P*part*P.transpose();
	 //	 if(std::exp(D[t1]) <= dmin || std::exp(D[t2]) <= dmin)  part.assign(0.0);
#endif
	 SinMat[k]=var*part*var.transpose();
	 //std::cout << "SinMat" << SinMat[k];
       }
       
     }
     
   }
   k=kk;
   if(!link_function) {
     var=residual_var.ginv0();
     var*=P_Mat[raneff]*P_Mat[raneff].transpose();
     for(t1=0;t1<numtrait;t1++) for(t2=t1;t2<numtrait;t2++,k++) {
       part.assign(0.0);
       part[t1][t2]=1.;
       part[t2][t1]=1.;
       SinMat[k]=var*part*var.transpose();
     }
   }
 }



 #define STRTEND(start,length) (start,start+length-1)	


 void GLMM::SSQCP(void)
 {
   int i,k,j,t1,t2,ii,jj,a_i,a_j,je,ipos,jpos,startaddr,iaddr,jaddr;
   int t_1,t_2,t;
   int *ainv_ia,*ainv_ja,ainv_nrow,nvc;
   //  GLMMModelTerm *localterm;
   UniformDist U;
   unsigned nlevels;
   double *ainv_a,xval;
   doubleMatrix *smat,wmat,*var_pt;
   wmat.resize(numtrait,numtrait);
   Vector<double> wvec;
   wvec.resize(numtrait);
   doubleMatrix vc_mat;
   Vector<double> k_new_vec,*svec,wrhs,wsol;
   k_new_vec.resize(numtrait);
   double *kvecx,*kvecy,*svecx,*svecy;
   kvecx=new double [hmmesize];
   kvecy=new double [hmmesize];
   svecx=new double [hmmesize];
   svecy=new double [hmmesize];
   wrhs.resize(hmmesize);
   wsol.resize(hmmesize);
   svec=new Vector<double> [act_numtrait];
   for(t1=0;t1<numtrait;t1++) {
     svec[t1].resize(hmmesize);
   }
   k_new_vec.resize(hmmesize);
   Need_PEV=1;
   glim(1);
   Need_PEV=0;
   //residual(res_nvc);
   //setup_mme(&rellrhs);
   //blup("ysmp1");
   nvc=TotalNvc();
   vc_mat.resize(nvc,numtrait);
   Build_SinMat();
   Info.resize(nvc,nvc);
   Score.resize(nvc);
   // Bulild SSs and CP for Random Effects;
   like_adj=0;
   int done_ped=0;
   int kk;
   for(kk=0,i=0;i<numterm;i++) {
     doubleMatrix tmpvar;
     if((term[i].classi() == 'P' && !done_ped) || term[i].classi() == 'R'){
       int nt;
       nt=nt_vec[i];
       kvec[i]=kk;
       if(var_link[i]!=i) {
	 kvec[i]=kvec[var_link[i]];
       }
       else{
	 kk+=(nt*(nt+1))/2;
       }
       k=kvec[i];
       var_pt=term[var_link[i]].prior->var_matrix();
       doubleMatrix Varinv;
       Varinv=*var_pt;
       Varinv.ginv1();
       startaddr=term[i].startaddr()+1;
       for(t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,k++) {
	 int t1_i,t2_i;
	 t1_i=trait_map[i][t1];
	 t2_i=trait_map[i][t2];
	 
	 if(term[t1_i].trait[t1%numtrait] && term[t2_i].trait[t2%numtrait]) {
	   smat= &SinMat[k];
	   if(term[i].classi() == 'R') {
	     nlevels=term[i].nlevel();
	     double pselect,scalvar;
	     int numrem,numneed,numsel;
	     int selected;
	     numneed=nlevels;
	     numrem=nlevels;
	     numsel=0;
	     if(AI_sample && AI_sample < nlevels) numneed=AI_sample;
	     scalvar=static_cast<double>(numrem)/static_cast<double>(numneed);
	     for(ii=1;ii<=nlevels;ii++) {
	       selected=0;
	       pselect=static_cast<double>(numneed)/static_cast<double>(numrem);
	       if(U.sample() <= pselect) selected=1;
	       numrem--;
	       numsel+=selected;
	       numneed-=selected;
	       iaddr=startaddr+(ii-1)*nt;
	       for(t=0;t<nt;t++) {
#ifndef DO_HINV
		 if(selected) {
		   k_new_vec.assign(0.0);
		   k_new_vec(iaddr+t)=1.;
		   hmmec.solve(svec[t],k_new_vec,"ysmp");
		 }
#endif
		 //Mem_zero(kvecx,hmmesize*sizeof(double));
		 //kvecx[iaddr+t-1]=1.;
		 //quad(kvecy,kvecx,svecy,svec[t].ve,1);
	       }
	       if(t1==0 && t2==0 && Use_Like) like_adj+=Varinv.quadratic(
									 blupsol.subvec(iaddr-1,nt),
									 blupsol.subvec(iaddr-1,nt));
	       Score[k]+=
		 smat->quadratic(
				 blupsol.subvec(iaddr-1,nt),
				 blupsol.subvec(iaddr-1,nt));
	       wmat.resize(nt,nt);
	       for(t_1=0;t_1<nt;t_1++) {
		 for(t_2=0;t_2<nt;t_2++) {
#ifdef DO_HINV
		   wmat[t_2][t_1]= -hInv.getaij(iaddr+t_1,iaddr+t_2);
#else
		   if(selected){
		     wmat[t_2][t_1]= -scalvar*svec[t_1][iaddr+t_2-1];
		   }
		   else{
		     wmat[t_2][t_1]=0;
		   }
#endif
		 }
	       }
	       wmat+=*(term[i].prior->var_matrix());
	       Score[k]-=
		 ((*smat)*wmat).trace();
	       
	     }
	   }
	   if(term[i].classi() == 'P') {
	     //   Score[k]-=((*smat)*CorrVar).trace();
	     
	     unsigned nanim=popsize-numgroup;
	     
	     double pselect,scalvar;
	     unsigned numrem,numneed,numsel;
	     int selected;
	     numneed=popsize;
	     numrem=popsize;
	     numsel=0;
	     unsigned *sel_vec;
	     sel_vec=new unsigned [popsize+1];
	     if(AI_sample && AI_sample < popsize) numneed=AI_sample;
	     scalvar=static_cast<double>(numrem)/static_cast<double>(numneed);
	     for(ii=1;ii<=popsize;ii++) {
	       selected=0;
	       if(U.sample() <= static_cast<double>(numneed)/static_cast<double>(numrem)) selected=1;
	       sel_vec[ii]=selected;
	       numrem--;
	       numsel+=selected;
	       numneed-=selected;
	     }
	     ainv();
	     ainv_ia=hainv.ia();
	     ainv_ja=hainv.ja();
	     ainv_a=hainv.a();
	     ainv_nrow=hainv.num_rows();
	     for(ii=1;ii<=popsize;ii++) {
	       
	       
	       iaddr=startaddr+(ii-1)*numtrait;
	       ipos=ii;
	       iaddr=startaddr+(ipos-1)*nt;
	       for(t=0;t<nt;t++) {
#ifndef DO_HINV
		 if(sel_vec[ii]) {
		   k_new_vec.assign(0.0);
		   k_new_vec(iaddr+t)=1.;
		   hmmec.solve(svec[t],k_new_vec,"ysmp");
		 }
#endif
		 //Mem_zero(kvecx,hmmesize*sizeof(double));
		 //kvecx[iaddr+t-1]=1.;
		 //quad(kvecy,kvecx,svecy,svec[t].ve,1);
	       }
	       je=ainv_ia[ii+1];
	       for(jj=ainv_ia[ii];jj<je;jj++) {
		 jpos=ainv_ja[jj];
		 double offdiag;
		 offdiag=2.;
		 if(ipos==jpos) offdiag=1.;
		 jaddr=startaddr+(jpos-1)*nt;
		 if(t1==0 && t2==0 && Use_Like) like_adj+=ainv_a[jj]*offdiag*
						  Varinv.quadratic(
								   blupsol.subvec(jaddr-1,nt),
								   blupsol.subvec(iaddr-1,nt));
		 
		 Score[k]+=ainv_a[jj]*offdiag*
		   smat->quadratic(
				   blupsol.subvec(jaddr-1,nt),
				   blupsol.subvec(iaddr-1,nt));
		 //	std::cout << "\nScore " << Score[k] << " k " << k << " jj " << jj << " iaddr " << iaddr << " jaddr " << jaddr;
		 wmat.resize(nt,nt);
		 for(t_1=0;t_1<nt;t_1++) for(t_2=0;t_2<nt;t_2++) {
		   //  Mem_zero(kvecy,hmmesize*sizeof(double));
		   //  kvecy[jaddr+t_2-1]=-1.;
		   //  wmat.me[t_2][t_1]=quad(kvecy,kvecx,svecx,svec[t_1].ve,2);
#ifdef DO_HINV
		   wmat[t_2][t_1]= -hInv.getaij(iaddr+t_1,jaddr+t_2);
		   //std::cout << "\n wmat " << wmat[t_2][t_1] << "t_2" << t_2 << "t_1" << t_1 << " iaddr " << iaddr << " jaddr " << jaddr;
#else
		   if(sel_vec[ii] ){// && sel_vec[jpos]) {
		     wmat[t_2][t_1]= -scalvar*svec[t_1][jaddr+t_2-1];
		   }
		   else{
		     //		    std::cout << ii << " " << jpos << " " << sel_vec[ii] << " " << sel_vec[jpos] << "\n";
		   }
#endif
		 }
		 //	wmat+=CorrVar;
		 //std::cout << "\nScore " << Score[k] << " k " << k;
		 Score[k]-=ainv_a[jj]*offdiag*
		   ((*smat)*wmat).trace();
		 //std::cout << "\nScore " << Score[k] << " k " << k  << " jj " <<jj << " ainv  " <<ainv_a[jj] << " offdiag "<< offdiag<< " wmat " <<  wmat<< " smat " <<  *smat << "prod" <<((*smat)*wmat) ;
	       }
	       //	std::cout << "\nScore " << Score[k] << " k " << k;
	       if(ii<=nanim )Score[k]-=((*smat)*(Var[i])).trace();
	       //std::cout << "\nScore " << Score[k] << " k " << k;
	     }
	     delete [] sel_vec;
	   }
	 }
       }
     }
   }

   //   ****
   // Build Residual SSs and CP
   //   ****
   k=kk;
   unsigned kvc;
   if(link_function){
     double scalvar=1.;
     if(AI_sample && AI_sample < numrec) scalvar=static_cast<double>(numrec)/static_cast<double>(AI_sample);	
     
     for(kvc=0;kvc<res_nvc;kvc++,k++){
       
       for(ii=0;ii<numrec;ii++) {
	 smat=&(observation[ii].Resid_Sin_Mat[kvc]);
	 Score[k]+=
	   smat->quadratic(observation[ii].residual,observation[ii].residual)*
	   observation[ii].weight;
	 wmat=observation[ii].Resid_Var/observation[ii].weight;

	 if(PEV_Scale){
	   //	   std::cout << "\nBefore " << Score[k]<< " " << wmat << observation[ii].pev;	 
	   
 	   wmat-=observation[ii].pev;
	 }
	 
	 //
	 // At this point I am not sure if it is better to adjust
	 // for prediction errors or leave it alone.
	 // For most of the cases I typically run accross it is
	 // probably not much of an issue.
	 //		
	 Score[k]-=
	   observation[ii].weight*
	   ((*smat)*wmat).trace();
	 

	 //	 std::cout << "\nAfter " << Score[k]<< " " << wmat;        
	 //	    if(k == 0) std::cout << "  ii " << Score[k];
       }
     }
     
     //std::cout << "\n";
   }
   else{
     for(t1=0;t1<numtrait;t1++) for(t2=t1;t2<numtrait;t2++,k++) {
       smat = &SinMat[k];
       for(ii=0;ii<numrec;ii++) {
	 Score[k]+=
	   observation[ii].weight*
	   smat->quadratic(observation[ii].residual,observation[ii].residual);
	 
	 wmat=observation[ii].Resid_Var/observation[ii].weight-observation[ii].pev;
	 Score[k]-=
	   observation[ii].weight*
	   ((*smat)*wmat).trace();
	 
       }
     }
     
     
   }
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);  
   FRHS.resize(nvc,hmmesize);
   FSol.resize(nvc,hmmesize);
   Fmat.resize(nvc,numtrait);

   int rec,Fmatrow;
   Fmatrow=numtrait;
   for(rec=0;rec<numrec;rec++) {
     if(link_function) {
       Fmatrow=observation[rec].H.num_rows();
       Fmat.resize(nvc,Fmatrow);
     }
     Fmat.assign(0.0);
     input_pos_val(modelfile);
     done_ped=0;
     for(kk=0,i=0;i<numterm;i++) {
       if((term[i].classi() == 'P' && !done_ped) || term[i].classi() == 'R'){
	 int nt;
	 nt=nt_vec[i];
	 if(var_link[i]==i) {
	   kk+=(nt*(nt+1))/2;
	 }
	 k=kvec[i];
	 int i_link=var_link[i];
	 // 	if(term[i].classi() == 'P') {
	 // 	  done_ped=1;
	 // 	  nt=act_numtrait;
	 // 	}
	 
 	doubleMatrix H,one;
 	H=observation[rec].H;
 	//H=H.kron(one.ones(1,nt/numtrait));
	
 	for(t1=0;t1<nt;t1++)
 	  for(t2=t1;t2<nt;t2++,k++){
 	    wvec.resize(numtrait);
 	    int t1_i,t2_i,t1_a,t2_a;
 	    t1_i=trait_map[i][t1];
 	    t2_i=trait_map[i][t2];
 	    t1_a=t1%numtrait;
 	    t2_a=t2%numtrait;
 	      int jt;ii=0;
 	      Vector<double> uscaled;
 	      uscaled.resize(nt);
 	      for(jt=0;(jt<nt) && (ii<=0);jt++)
 		{
 		  int t3_i;
 		  t3_i=i;
		  t3_i=trait_map[i][jt];
 		  ii=term[t3_i].addr[jt%numtrait]-jt%numtrait;
 		}
 	      uscaled = Var[i] * SinMat[k] * blupsol.subvec(ii-1, nt);
	      //std::cout << "\nuscaled "<<uscaled<< "Var" << Var[i] <<"SinMat " << SinMat[k] << "u" << blupsol.subvec(ii-1, nt);
 	      int t3;
 	      for(t3=0;t3<nt;t3++) {
 		int t3_i;
		t3_i=trait_map[i][t3];
 		ii=term[t3_i].addr[t3%numtrait];
 		if(ii) {
 		  int ioff=ii-t3;
 		  uscaled=Var[i]*SinMat[k]*blupsol.subvec(ioff-1,nt);
 		  xval=xval_term[t3_i];
 		  wvec[t3%numtrait] += xval*uscaled[t3];
 		}
 	      }
 	      if(link_function) wvec=H*wvec;
	      for(int iii=0;iii<Fmatrow;iii++) Fmat[k][iii]+=wvec[iii];
 	      //memcpy(Fmat[k],wvec.begin(),sizeof(double)*Fmatrow);
	      //}
 	  }
       }
     }
     k=kk;

     if(link_function) {

       for(kvc=0;kvc<res_nvc;kvc++,k++) {
	
 	wvec=//observation[rec].Hinv*
 	  observation[rec].Resid_Var*
 	    observation[rec].Resid_Sin_Mat[kvc]*
 	      observation[rec].residual;
	
 	memcpy(Fmat[k],wvec.begin(),sizeof(double)*Fmatrow);
       }

     }
     else{
       for(t1=0;t1<numtrait;t1++) for(t2=t1;t2<numtrait;t2++,k++){
	
 	wvec=observation[rec].Resid_Var*SinMat[k]*observation[rec].residual;
	
 	memcpy(Fmat[k],wvec.begin(),sizeof(double)*Fmatrow);
       }

     }
     doubleMatrix vinv;
     if(observation[rec].Resid_Var.empty()) {
       vinv=observation[rec].rve;
     }
     else {
       vinv=observation[rec].Resid_Var;
       vinv.ginv1();
     }
     //std::cout << "Info B"<< Info << "vinv" << vinv << "Fmat"<<Fmat;
     Info+=Fmat*vinv*Fmat.transpose()*
       observation[rec].weight;
     //std::cout << "Info A"<< Info;
     //    if(Info[20][20] < -1.e-6) {
     //      std::cout << "Wow!\n";
     //    }

     vc_mat=observation[rec].weight*Fmat*vinv*observation[rec].H;


     for(t1=0;t1<numtrait;t1++) {
       for(i=0;i<numterm;i++) {
 	// localterm=(GLMMModelTerm *)&term[i];
 	ii=term[i].addr[t1];
 	if(ii) {
 	  ii--;
 	  xval=xval_term[i];
 	  for(j=0;j<nvc;j++){
 	    FRHS[j][ii]+=vc_mat[j][t1]*xval;
 	    //    if(fabs(vc_mat[j][t1]*xval) > 1.e4){
 	    //  std::cout << "Wow2 " ;
 	    //}
 	  }
 	}
       }
     }
   }
   for(k=0;k<nvc;k++) {
     hmmec.solve(FSol[k],FRHS[k],"ysmp");
   }
   Info-=
     FRHS*FSol.transpose();
   Info*=.5;
   Score*=.5;
   // modelfile.close();

   delete [] svec;
   delete [] svecy;
   delete [] svecx;
   delete [] kvecy;
   delete[] kvecx;

 }



 SparseMatrix& GLMM::ainv(void)
 {
   if(!Need_Ainv) return(hainv);
   Need_Ainv=0;
   //TW  int do_west;
   //TW Vector<double> West;
   //TW  do_west=pop->do_west;
   //TW  West=pop->West;
   unsigned nanim=popsize-numgroup;
   hainv.resize(popsize,4*popsize);
   double dii,val;

   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   unsigned asd[3];
   int  ii,jj,j;
   Individual *I;
   double f1,f2;
   ldet=0.;
   for (int k=0; k<nanim; k++) {
     I = pop->member(k);
     asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
     dii = 4.0/(2.0 - I->father_inbcoef() - I->mother_inbcoef());;
     f1= I->father_inbcoef() ;
     f2= I->mother_inbcoef();
     if(asd[1]> nanim) f1=-1.;
     if(asd[2]> nanim) f2=-1.;
     dii = 4.0/(2.0 -f1-f2);
     ldet+=std::log(dii);
     for (int i=0; i<3; i++) {
       if (asd[i] != 0) {
 	ii = asd[i];
 	for (j=0; j<3; j++) {
 	  if (asd[j] != 0) {
 	    jj = asd[j];
 	    val = lambda[i][j]* dii;
	    // std::cout << "\n ii " << ii << " jj " << jj << " val " << val << " dii " << dii;
 	    if(ii <= jj ) hainv.insert(ii,jj,val);
 	  }
 	}
       }
     }
   }
   ldet=-2.*ldet;
   //std::cout << "ldet " << ldet << "\n";
   //std::cout << hainv.dense();
   hainv.logdet();
   hainv.close();

   return(hainv);
 }
doubleMatrix GLMM::AI_REML(int numiter,double tol,double info_scale)
{
  int Converged,Need_Like;
  aireml_called=-1;
  Converged=0;
  Need_Like=1;
  Vector<double> varold,varnew;
  int Update_estimate_AI;
  int nvc,iteration,numran,k,t1,t2,raneff,kidx,kend;
  Update_estimate=1;
  // info_scale=1.1;
  glim(num_glmm);
  //  glim(5)->display();
  Update_estimate=1;
  nvc=TotalNvc();
  doubleMatrix varbound;
  varbound.resize(numtrait,numtrait);
  Vector<double> eig,flg,fix;
  eig.resize(numtrait);
  flg.resize(numtrait);
  fix.resize(nvc);
  varold.resize(nvc);
  varnew.resize(nvc);
  numran=Var2Vec(varold.begin()); //Have Numran return Number af random
  // effects plus one for P if it exists
  
  doubleMatrix Adj_Cov;
  if(numran) P_Mat = new doubleMatrix [numran];
  doubleMatrix Ip,Info_orig;
  double log_like_old,log_like_new;
  if(Need_Like){ 
    log_like_old=restricted_log_likelihood(0);
    std::cout << "\nOriginal  Residual Log Likelihood:" << std::setprecision(SESSION.output_precision) << log_like_old <<"\n\n";
    Need_Like=0;
  }
  for(k=0;k<numran;k++) P_Mat[k].identity(nt_vec2[k]);
  //if(Print_Level == 0) std::cout << "Initial Estimates " << varold;
  long subit;
  for(iteration=0;(iteration<numiter) && !Converged;iteration++){
    
    Update_estimate=1;
    //log_like_old=restricted_log_likelihood(1);
    //std::cout << "\nOld Residual Log Likelihood:" << log_like_old <<"\n";
    SSQCP();
    if(Need_Like) {
      //glim(2);
      log_like_old=restricted_log_likelihood(0);
      Need_Like=0;
    }
    if(Print_Level > 0)  std::cout << "\n Old Residual Log Likelihood:" << std::setprecision(SESSION.output_precision)<< log_like_old <<"\n";
    Var2Vec(varold.begin());
    subit=0;
    
    info_scale=0.;
    Info_orig=Info;
    do{
      Info=Info_orig;
      //      std::cout << "\nOrig Info" << Info_orig;
    if(info_scale ){
      double info_mult=std::exp(info_scale);
      //std::cout << "\n Info Mult " << info_mult;
      for(k=0;k<nvc;k++) Info[k][k]*=info_mult;
    }
    //std::cout <<"\n Fix " << fix;
    for(k=0;k<nvc;k++) {
      if(Info[k][k] < EPSILON || fix[k]) {
 	for(int k2=0;k2<nvc;k2++) {
 	  Info[k][k2]=0;
 	  Info[k2][k]=0;
 	}
 	Info[k][k]=10.*EPSILON;
 	//fix[k]=1;
      }
    }
#ifdef DO_CHOL_NOT
    //varbound=Info;
    //    std::cout << "Eigen" << varbound.eigen();
    /*
    for(k=0;k<nvc;k++) {
    if(fix[k] ) {
    Score[k]=0;
    for(int k2=0;k2<nvc;k2++) {
	  Info[k2][k]=0;
	  Info[k][k2]=0;
	}
	Info[k][k]=1.e-10
      }
    }
    */
    //    Info_orig=Info;
    //std::cout << "Score " << Score;
    //std::cout << "Info " << Info;
    //    std::cout << "ascov 1" << Info_orig.ginv1();
    //    std::cout << "ascov 0" << Info_orig.ginv0();
    /*Info_orig=Info;*/
#endif
    //    std::cout << "Score " << Score;
    //    std::cout << "Info " << Info;
    Vector<double> Delt;
    Delt=Info.ginv0()*Score;
    varnew=varold+Delt;
    //std::cout << "\nVar " << varnew << varold << Delt;
    //    if(Print_Level >= 0) {
    //  std::cout << " New Estimates before Adjustment" << varnew;
    //}
#ifdef DO_LOG
    doubleMatrix P,L,D;
    
    Adj_Cov.identity(nvc);
    int offset=nvc-res_nvc;
    for(k=0;k<res_nvc;k++)
      for(int k2=0;k2<res_nvc;k2++)
	Adj_Cov[offset+k][offset+k2]=observation[numrec].Adj_Cov[k][k2];
    

    Vector<double> d; 	
    int nt;
    for(k=0,raneff=0;raneff<numran;raneff++) {
      nt=nt_vec2[raneff];
      varbound.resize(nt,nt);
      for(kidx=k,t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,kidx++) {
 	varbound[t1][t2]=varold[kidx];
 	varbound[t2][t1]=varold[kidx];
      }
      //      std::cout << "varbound "<< varbound;
      P=varbound;
      //std::cout << "\nP " << P;
      d=P.eigen();
      D.resize(nt,nt);
      for(t1=0;t1<nt;t1++){
       d[t1]=std::log(d[t1]);
       D[t1][t1]=d[t1];
      }
	  

      doubleMatrix part;
      
      for(kidx=k,t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,kidx++) {
 	part.resize(nt,nt);
 	double D_avg,D_delt;
 	D_avg=.5*(d[t1]+d[t2]);
 	D_delt=.5*(d[t1]-d[t2]);
 	if(fabs(D_delt) > 1.e-8) {
	  
	  //std::cout << "how did I get here!!!\n";
	  part[t1][t2]=std::exp(D_avg)*(sinh(D_delt))/(D_delt);
	  part[t2][t1]=std::exp(D_avg)*(sinh(D_delt))/(D_delt);
	}
	else{
	  part[t1][t2]=std::exp(D_avg);
	  part[t2][t1]=std::exp(D_avg);
	}
	//std::cout << "part\n" << part <<"d\n"<< d;
	part=P*part*P.transpose();
	if(d[t1] <= -9. && d[t2] <= -9.)  part.assign(0.0);
	int kkidx,tt1,tt2;
	for( kkidx=k,tt1=0;tt1<nt;tt1++) for(tt2=tt1;tt2<nt;tt2++,kkidx++) {
	  //	  std::cout <<"\n" << kkidx << " " << kidx << " " << tt1 << " " << tt2 << " " << nt << " " << raneff << " " << k <<"\n",
	  Adj_Cov[kkidx][kidx]=part[tt1][tt2];
	}
      }
      
      /*
      double dmax=0.;
      for(kidx=k,t1=0;t1<nt;t1++) {
	for(t2=t1;t2<nt;t2++,kidx++) {
	  if(fabs(Delt[kidx])> dmax) dmax=fabs(Delt[kidx]);
	}
      }

      
#define	     MAXCHG 1.
	    if(dmax > MAXCHG){
	std::cout << "\ndmax " << dmax <<"\n";
	for(kidx=k,t1=0;t1<nt;t1++) {
	  for(t2=t1;t2<nt;t2++,kidx++) {
	  Delt[kidx]/=dmax/MAXCHG;
	  }
	}
      }
      */
      //std::cout << "D" << D << "P" << P;
      for(kidx=k,t1=0;t1<nt;t1++) {
	for(t2=t1;t2<nt;t2++,kidx++) {
	D[t1][t2]+=Delt[kidx];
	D[t2][t1]=D[t1][t2];
	}
      }
      
      //std::cout << "D new" << D;

      D=P*D*P.transpose();
      doubleMatrix Q;
      doubleMatrix DE;
      Q=D;
      d=Q.eigen();
      //std::cout << " d " << d << " Q " << Q << " D" << D;
      double dmax=d.max();
      //if(dmax < -8.) std::cerr << "\nD\n" << D << "\nP\n" << P<< "\n";
      if(dmax > 15.) dmax=15.;//
      if(dmax < 0.) dmax=0;
      double dmin=std::log(10.*SESSION.epsilon)+dmax;
      for(t1=0;t1<nt;t1++) {
	if(d[t1] < dmin) d[t1]=dmin+std::log(.5);
	if(d[t1] > dmax) d[t1]=15;
      }
      DE.resize(nt,nt);
      for(t1=0;t1<nt;t1++)DE[t1][t1]=std::exp(d[t1]);
      D=Q*DE*Q.transpose();
      Q=D;
      d=Q.eigen(); 
      //std::cout << " d exp " << d << " Q " << Q;
      varbound=D;//P*D*P.t();
      //std::cout << " d " << d;
      //   std::cout << " New Estimates before Adjustment" << varbound;
      for(t1=0;t1<nt;t1++) {
	//	varbound[t1][t1]+=1.0e-4;
	if(varbound[t1][t1] <  std::exp(-11.)) {
	  for(t2=0;t2<nt;t2++){
	    varbound[t1][t2]=0.;
	    varbound[t2][t1]=0.;
	  }
	  varbound[t1][t1]=std::exp(-11.5);
	}
      }
      //std::cout << "varbound	 new" << varbound;
      for(t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,k++) {
	varnew[k]=varbound[t1][t2];
      }
    }
    
#endif

    
    //	 ****
    //	  If the variance covariance matrix is not pd make it slighty pd
    //	 ****
    //Print_Level=2;
#ifdef DO_CHOL_2
    int nt;
    for(k=0,raneff=0;raneff<numran;raneff++) {
      nt=nt_vec2[raneff];
      flg.zeros(nt);
      
      for(t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,k++) {
	if(t1==t2 && varnew[k]<1.e-8 && !fix[k]) {
	  Score[k]=0;;
	  varnew[k]=1.e-4;
	  flg[t1]=1;
	  fix[k]=1;
	}
	else if(flg[t1] && t1!=t2) {
	  Score[k]=0.;
	  varnew[k]=0.;
	  fix[k]=1;
	}
       }
     }
	
    //	  for(k=0;k<nvc;k++) {
    //	    if(fix[k]){
    //	     for(int k2=0;k2<nvc;k2++) {
    //	      Info[k][k2]=0;
    //	      Info[k2][k]=0;
    //	     }
    //	     Info[k][k]=1.e-4;
    //	    }
    //	  }
    // varnew=varold+Info.ginv0()*Score;
#endif

#ifdef DO_PMAT	 
    int nt;
    for(k=0,raneff=0;raneff<numran;raneff++) {
      nt=nt_vec2[raneff];
      varbound.resize(nt,nt);
      eig.resize(nt);
      flg.resize(nt);
      for(kidx=k,t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,kidx++) {
	varbound[t1][t2]=varnew[kidx];
	varbound[t2][t1]=varnew[kidx];
      }

      for(t1=0;t1<nt;t1++){
	if(varbound[t1][t1] < 1.e-8){
	  for(t2=0;t2<nt;t2++) {
	    varbound[t1][t2]=0;
	    varbound[t2][t1]=0;
	  }
	}
       }
	
       //	    P_Mat[raneff]=P_Mat[raneff]*P_Mat[raneff].transpose()*varbound*P_Mat[raneff]*P_Mat[raneff].transpose();
	
	
       P_Mat[raneff]=varbound;
	
       eig=P_Mat[raneff].eigen();
       if(Print_Level > 1) {
	 //	std::cout << "varbound" << varbound << "\n";
	 //	std::cout << "Eigen Values " << eig << "\n";
	 //	std::cout << "P_Mat" << P_Mat[raneff] << "\n";
      }
	    
      for(t1=0;t1<nt;t1++){
	flg[t1]=0;
	if(eig[t1]<1.e-8){
	  eig[t1]=.99e-3;
	  flg[t1]=1;
	}
	if(eig[t1]>1.e10){
	  eig[t1]=1.01e10;
	  flg[t1]=0;
	}
      }
      if(Print_Level > 1) {
	//	std::cout << "Adjusted Eigen Values " << eig << "\n";
       }
       varbound=P_Mat[raneff]*eig.diag()*P_Mat[raneff].transpose();
       for(t1=0;t1<nt;t1++){
	if(flg[t1]){
	  for(t2=0;t2<nt;t2++) {
	    P_Mat[raneff][t1][t2]=0;
	  }
	}
      }
	    
      for(t1=0;t1<nt;t1++) for(t2=t1;t2<nt;t2++,k++) {
	varnew[k]=varbound[t1][t2];
	      
      }
      eig=varbound.eigen();
      if(Print_Level > 1) {
	//	std::cout << " Eigen Values bounded Matrix " << eig << "\n";
      }
	 
    }
#endif /* DO_PMAT */
    Vec2Var(varnew.begin());
    //std::cout << "Varnew " <<varnew;
    glim(num_glmm); 
    log_like_new=restricted_log_likelihood(0);
    std::cout << " Iteration " << (iteration+1) <<"."<< subit << " Res Log Like "<< std::setprecision(SESSION.output_precision)<< log_like_new << " Change " << std::setprecision(SESSION.output_precision)<<(log_like_new-log_like_old) << "\n";
    info_scale+=1.5*ranf();
  }while(tol && log_like_new-log_like_old < -tol && subit++<= 10);
    if(subit == 0 && std::abs(log_like_new-log_like_old) < std::max(tol,EPSILON) ) Converged=1;
    if(Converged) aireml_called=1;
    log_like_old=log_like_new;
    info_scale-=.75;
    varold=varnew;

    //	  if(iteration%3 == 2) {
    //	    Update_estimate=1; 
    //	    glim(2);
    //	    Update_estimate=1;
    //	  }

    //        Print_Level=2;
    
    if(Print_Level > 0) {
      std::cout << " Iteration " << iteration << "  Subit " << subit << "\n";
      std::cout << " Residual log likelihood " 
		<< std::setprecision(SESSION.output_precision)
	   << log_like_new << "\n";
      std::cout << " New Estimates \n" << varnew;
      if(Print_Level > 1) {
	std::cout << " Score \n" << Score;
	std::cout << " Info	\n" << Info;
	std::cout << " Asy Cov \n" << Info.ginv0();
      }
    }  


     Vec2Var(varnew.begin());
     varold=varnew;
   }
   if(Print_Level >= 0) {
     if(Converged){
       std::cout << "\n Iteration " << iteration  << " Converged\n";
     }
     else{
       std::cout << "\n Iteration " << iteration  << "  Subit " << subit << " Failed to Converge\n";
     }
     std::cout << " Final Estimates \n" << varnew;
     std::cout << " Last Change \n" << Info.ginv0()*Score;
     Info=Info_orig;
     if(info_scale) std::cout << " Unscaled Last Change \n" << Info.ginv0()*Score;
     std::cout << " Residual log likelihood "
	       << std::setprecision(SESSION.output_precision)
	 << restricted_log_likelihood(0) << "\n" ;
     Info=Info_orig;
     //  std::cout << "    Adj Matrix " << Adj_Cov;
     std::cout << " Asy Covariance Matrix\n" <<
 #ifdef DO_LOG
       Adj_Cov*
 #endif
       Info.ginv0()
 #ifdef DO_LOG
     *Adj_Cov.transpose()
 #endif
       ;
   }
   Update_estimate=1;
   //glim();
   //  for(k=0;k<numran;k++) P_Mat[k].resize(0,0);
   if(P_Mat)delete [] P_Mat;
   P_Mat=NULL;
   varest=varnew;
   doubleMatrix ans;
   ans = hadjoin(varest,
 #ifdef DO_LOG
     Adj_Cov*
 #endif
		 Info.ginv0()
 #ifdef DO_LOG
     *Adj_Cov.transpose()
 #endif
   );
   return ans;
 }


 SparseMatrix* GLMM::setup_mme(Vector<double>* rhs)
 {
   /////////////////////////////////////////////////
   // build up MME: hmmec, the coefficient matrix,
   //               rellrhs, right-hand-side
   ///////////////////////////////////////////////////

   if(link_function == NULL) {
     return Model::setup_mme( rhs);
   }


   if (type == bad_model) {
     warning("GLMM::setup_mme(): bad model");
     return 0;
   }
   if (!data_prepared) {
     if (!prepare_data()) {
       type = bad_model;
       return 0;
     }
   }
   if (hmmesize==0) return &hmmec;

   if (type == mixed_model) {
      ;
   }
   else if (type == fixed_model) {
      ;
   }
   else {
     warning("Model::setup_mme(): inappropriate model");
     return 0;
   }
   if(res_nvc) {
     varcomp.reserve(res_nvc);
   }

   if(!observation)new_obs();
   if(hmmec.nz() != 0 && Need_Residual) {
     residual(res_nvc*Need_PEV);
   }
   //hmmec.release();
   //hmmec.resize(0,0);
   //std::cout << "\n hmmesize " << hmmesize << " max_nz " << max_nz << "\n";
   hmmec.resize(hmmesize,max_nz);
   if(Sand_Calc) hSand.resize(hmmesize,max_nz);
   Vector<double> localrhs;
   Vector<double> *rhsvec = &localrhs;
   if (rhs) rhsvec = rhs;

   rhsvec->resize(hmmesize);
   setup_ww(rhsvec);
   Model::add_G_1();
   if(Sand_Calc) {
     add_G_Sand();
     hSand.close();
   }
 #ifdef OLD_SPARSE
   non_zero = hmmec.close();
 #else
   non_zero=hmmec.nz();
 #endif
   //std::cout << "\nNon Zero  " << non_zero << "\n";
   Need_Residual=1;
   return &hmmec;
 }

 int GLMM::equation(const std::string &modelspecs){

  
   int retval;
   retval=Model::equation(modelspecs);


   return(retval);
 }


 void GLMM::setup_ww(Vector<double>* rhsvec)
 {
   if(!link_function) {
     Model::setup_ww(rhsvec);
     return;
   }
   //  std::cout << " setup_ww ";
   unsigned i,j,ii,jj,t1,t2,rec;
   double vy,xval,val,otherval;;
   memset(rhsvec->begin(),'\0', sizeof(double)*hmmesize);
   double *rhstmp = rhsvec->begin() -1;
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);  
   if (!modelfile) throw exception(" GLMM::setup_ww(): cannot open binary file");
   double **vep;
   Matrix<double> ve(numtrait,numtrait);
   Vector<double> scale_res(numtrait);
   doubleMatrix ves(numtrait,numtrait);
   for (yry=0, rec=0; rec<numrec; rec++) {
     input_pos_val(modelfile);
     observation[rec].pattern=pattern[pos_term[numterm]];
     if (ntermGdist) {
       trait_vec[0] -= pop->popmember[rec_indid[rec]-1]->xbzu();
     }
     if(observation[rec].need_data) {
       memcpy(observation[rec].y.begin(),trait_vec,sizeof(double)*numtrait);
       observation[rec].need_data=0;
     }
     (*link_function)(observation[rec]);
     if(Sand_Calc) {
       scale_res=observation[rec].rpy-observation[rec].rve*observation[rec].estimate;
       ves=scale_res*scale_res;
     }
     vep = rve[pos_term[numterm]].begin();
     val = xval_term[numterm];        // val = weight  variable
     if(Use_Like) yry -= 2.*observation[rec].like_val*val;
     if(!Use_Like) yry+=observation[rec].Resid_Var.logdet()*val;
     for (t1=0; t1<numtrait; t1++) {
       if(!Use_Like) yry+=observation[rec].y[t1]*observation[rec].rpy[t1]*val;
       for (t2=0; t2<numtrait; t2++)
 	ve[t1][t2] = val*observation[rec].rve[t1][t2];

     }
     for (t1=0; t1<numtrait; t1++) {
       vy = val*observation[rec].rpy[t1];

       for (i=0; i<numterm; i++) {
 	ii = term[i].addr[t1];
 	if (ii) {
 	  xval = xval_term[i];
 	  for (t2=t1; t2<numtrait; t2++) {
 	    otherval = xval* ve[t1][t2];
 	    jj = term[i].addr[t2];
 	    if (jj) {
 	      hmmec.insert(ii,jj,otherval*xval);
 	      if(Sand_Calc) hSand.insert(ii,jj, xval* (ves[t1][t2])*xval);
 	    }
 	  }
 	  for (t2=0; t2<numtrait; t2++) {
 	    otherval = xval * ve[t1][t2];
 	    for (j=i+1; j<numterm; j++) {
 	      jj = term[j].addr[t2];
 	      if (jj) {
 		hmmec.insert(ii,jj,otherval*xval_term[j]);
 		if(Sand_Calc) hSand.insert(ii,jj, xval* (ves[t1][t2])*xval_term[j]);
 	      }
 	    }
 	  }
	  //std::cout << "\nBuilding "<< rec <<" trait " << t1  << " pos " << ii << '\n' << vy <<"\nrpy" << observation[rec].rpy << "\nrve" <<observation[rec].rve<<'\n' <<ve;
 	  rhstmp[ii] += vy*xval_term[i];
 	}
       }
     }
   }

   //  modelfile.close();
 }


 unsigned GLMM::TotalNvc(void) const
 {
   ////////////////////////////////////////////////////
   //  total number of distinct variance components
   /////////////////////////////////////////////////////
	
   unsigned nvc,i,np;
   if(link_function) {
     nvc=res_nvc;
   }
   else{
     nvc = (numtrait*(numtrait + 1))/2; //
   }
   np=0;
   for (i=0; i<numterm; i++) {

     if (term(i).classi() == 'R' || term(i).classi() =='P') {
       if(var_link[i]==i) nvc += nt_vec[i]*(nt_vec[i] + 1)/2;
     }
   }
   if(np) {
     nvc +=(np)*(np-1)*numtrait*numtrait/2;
     //ncorr=np;
   }

   return nvc;
 }

 int GLMM::Var2Vec(double *x)
   /*******************************************************************
       * put variance components into a double vector
       * user is totally responsible to have an appropriate size for vector x
       * and matrix residual_var
       *********************************************************************/
 {
   int i,k,t1,t2,numran,ped_done,nt;
   doubleMatrix *var;
 #ifdef DO_CHOL
   doubleMatrix L;
 #endif
   numran=0;
   ped_done=0;
   if(nrandom) {
     nrandom=0;
   }
   if(nt_vec2) {
     delete [] nt_vec2;
   }
   nt_vec2=new int [numterm+1];
   for (k=0,i=0; i<numterm; i++) {
     if ((term[i].classi() == 'R' ||  (!ped_done && term[i].classi() == 'P'))&&(var_link[i]==i)) {
       numran++;
       nt=nt_vec[i];
       nt_vec2[nrandom++]=nt;
       var = term[i].prior->var_matrix();
       if(term[i].classi()=='R') {
 #ifdef DO_CHOL
 	L=*var;
 	L.chol();
 	L=L.transpose();
 	var=&L;
 #endif
 	for (t1=0; t1<nt; t1++) for (t2=t1; t2<nt; t2++) {
 	  x[k++] = (*var)[t1][t2];
 	}
       }
       if(term[i].classi()=='P'){
 #ifdef DO_CHOL
	 // 	build_CorrVar();
 	L=var;
 	L.chol();
 	L=L.transpose();
 	var=&L;
 	for (t1=0; t1<nt; t1++) for (t2=t1; t2<nt; t2++) {
 	  x[k++] = (*var)[t1][t2];
 	}
 #else
	// 	for(int ii=i;ii<numterm;ii++) {
	// 	  for(t1=0;t1<numtrait;t1++){
	// 	    if(corrmap[ii]) {
	// 	      var=term[ii].prior->var_matrix();
	// 	      for(int jcorr=ii;jcorr<numterm;jcorr++) {
	// 		if(corrmap[jcorr]) {
	// 		  int imult=1;
	// 		  if(ii != jcorr) {
	// 		    var=&corrvar[ii][jcorr];
	// 		    imult=0;
	// 		  }
	// 		  for (t2=imult*t1; t2<numtrait; t2++) {
	// 		    x[k++] = (*var)[t1][t2];
	// 		  }
	// 		}
	// 	      }
	// 	    }
	// 	  }
	// 	}

 	for (t1=0; t1<nt; t1++) for (t2=t1; t2<nt; t2++) {
 	  x[k++] = (*var)[t1][t2];
 	}
 #endif
       }
     }
   }
   /*  if(ncorr) {
       for(int ict=0;ict<numterm;ict++) {
       if(corrmap[ict]) {
       for(int jct=(ict+1);jct<numterm;jct++) {
       if(corrmap[jct]) {
       for (t1=0; t1<numtrait; t1++) for (t2=t1; t2<numtrait; t2++) {
       x[k++] = corrvar[ict][jct].me[t1][t2];
       }
       }
       }
       }
       }
       }*/
   if(link_function) {
     for(t1=0;t1<res_nvc;t1++) x[k++]=varcomp[t1];
   }
   else{
     numran++;
     nt_vec2[nrandom++]=numtrait;
     for (t1=0; t1<numtrait; t1++) {
       for (t2=t1; t2<numtrait; t2++) x[k++] = residual_var[t1][t2];
     }
   }
   return(numran);
   }

// This is a function that I am adding to make the normal_log link function
// work from C++. (March 2, 2004; RLF)
void GLMM::normalLog(void){
	int n_nvc= (ntrait()*(ntrait()+1))/2;
	link(&normal_loginClasses,n_nvc);
	num_glmm=1;
	Like_LR(1);
}
	

// I am moving the following stuff from the link.cc file that is in the interface directory to make the normal_log link function
// work from C++. (March 2, 2004; RLF)

#define twopi 6.283195307179587
#define invsqrt2pi .3989422804014326779399461 
#define BOUND(est,L,U) (est=(est > U ? U:(est < L ? L:est) ));
#define BOUND_FLG(est,L,U,flg) (flg=0,(if (est < L || est > U )flg=1),est=(est > U ? U:(est < L ? L:est) ));
#define UPPER_BOUND(est,U) (est=(est > U ? U:est ));
#define LOWER_BOUND(est,L) (est=(est < L ? L:est ));
	
void normal_loginClasses(matvec::Observe &observe) {
	int numtrait=observe.numtrait;
	
	doubleMatrix H,variance;
	doubleMatrix Alog;
	doubleMatrix P,R;
	Vector<double> D;
	// double final=0;
	//if(observe.param) final=*((double *) observe.param[0]);
	int k,t1,t2;
	if(observe.rec < 0) {
		if(observe.nvc) { 
			if(observe.nvc != (numtrait*(numtrait+1))/2) {
				throw matvec::exception("Error: Should be normal(,%d)");
				return;
			} 
			
			
			//cout <<  variance << observe.resid_var[0];
			*observe.resid_var=observe.resid_var[0];
			Alog=(observe.resid_var[0]).mat_log();
			for(k=0,t1=0;t1<numtrait;t1++) 
				for(t2=t1;t2<numtrait;t2++,k++) 
					observe.varcomp[0][k]=Alog[t1][t2];
		}
		return;
	}
	Vector<double> mean,y;
	Alog.resize(numtrait,numtrait);
	
	Vector<doubleMatrix> part(observe.nvc);
	for(k=0;k<observe.nvc;k++) part[k].resize(numtrait,numtrait);
	if(observe.nvc) {
		
        for(k=0,t1=0;t1<numtrait;t1++) 
			for(t2=t1;t2<numtrait;t2++,k++) {
				Alog[t1][t2]=observe.varcomp[0][k];
				Alog[t2][t1]=observe.varcomp[0][k];
			}
				variance=Alog.mat_exp();
		//cout << variance;
		Alog=variance.mat_log();
		
		for(k=0,t1=0;t1<numtrait;t1++) 
			for(t2=t1;t2<numtrait;t2++,k++) 
				observe.varcomp[0][k]=Alog[t1][t2];
		*observe.resid_var=variance;
		variance.mat_exp_der(part);
		doubleMatrix Atmp;
		/*
		 for(int ii=0,t1=0;t1 <numtrait;t1++)  
		 for(t2=t1;t2<numtrait;t2++,ii++) {
			 Atmp=Alog;
			 Atmp[t1][t2]+=.0001;
			 Atmp[t2][t1]=Atmp[t1][t2];
			 Atmp=(Atmp.mat_exp()-Alog.mat_exp())/.0001;
			 cout << ii << part[ii]<< endl;
			 cout << ii << Atmp << endl;
		 }
		 */
		
	}
	else{
		variance=observe.resid_var[0];
	}
	//cout << "\n Variance " << variance << *observe.resid_var;
	doubleMatrix Miss;
	Miss.resize(numtrait,numtrait);
	//cout << "\n Pattern" ;
	for(t1=0;t1<numtrait;t1++) {
		//cout <<' ' <<observe.pattern[t1];
		if(observe.pattern[t1]=='1') Miss[t1][t1]=1;
	}
	//Miss.identity();
	
	//static double mean,H,variance,scale,y;
	y=observe.y;
	y=Miss*y;
	//cout << "\n y\n" << y;
	// if(observe.estimate[0] < -5.) observe.estimate[0]=-5.;
	//  if(observe.estimate[0] > 5.) observe.estimate[0]=5.;
	mean=Miss*observe.estimate;
	
	variance=*observe.resid_var;
	variance=Miss*variance*Miss;
	//cout << "\n Variance \n" << variance; 
	observe.Resid_Var=variance;
	variance=variance.ginv1();
	//cout << "\n Variance inv \n" << variance; 
	observe.like_val=.5*(variance.logdet()-variance.quadratic((y-mean),(y-mean))) ;
	
	
	observe.H=Miss;
	H=observe.H;
	//observe.Hinv[0][0]=1./H;
	observe.rve=H.transpose()*variance*H;
	
	observe.rpy=H.transpose()*variance*y;
	//cout <<"\n rve rpy "<<observe.rve << " " << observe.rpy<< " H "<< H << "pattern" << observe.pattern <<" y"<<  y <<" var" << variance;
	observe.residual=y-mean;
	if(observe.nvc) {
		for(k=0,t1=0;t1<numtrait;t1++)
			for(t2=t1;t2<numtrait;t2++,k++) {
				
				observe.Resid_Sin_Mat[k]=variance*Miss*part[k]*Miss*variance;
				//	cout << k << " " << part[k] << observe.Resid_Sin_Mat[k] << endl;
			}
				
	}
}

	



 } ///////// end of namespace matvec


