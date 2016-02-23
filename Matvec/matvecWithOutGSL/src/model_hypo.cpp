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

#include <string>
#include <vector>
#include <iomanip>

#include "session.h"
#include "util.h"
#include "model.h"
#include "data.h"
#include "stat.h"

namespace matvec {




double Model::estimate(const Vector<double>& Kp)
{
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////

   double **kpme = new double *[1];
   kpme[0] = Kp.begin();
   unsigned nc = Kp.size();
   double kpb[1];
   estimate((const double **)kpme,1,nc,kpb);
   if(kpme) {
     delete [] kpme;
     kpme = 0;
   }
   return kpb[0];
}

Vector<double> Model::estimate(const doubleMatrix& Kp)
{
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////

   unsigned nr = Kp.num_rows();
   unsigned nc = Kp.num_cols();
   Vector<double> kpb(nr);
   estimate((const double **)Kp.begin(),nr,nc,kpb.begin());
   return kpb;
}


doubleMatrix Model::CovMat(doubleMatrix &Kp){
  unsigned nr=Kp.num_rows();
  unsigned nc = Kp.num_cols();
  double *kp,*sol;
  doubleMatrix KvK(nr,nr);
  if (!get_blupsol("ysmp")) throw exception("Model::CovMat(): bad model");
  if (nc > hmmesize) throw exception("Model::CovMat(): size not conformable");
  
  Vector<double> xy(hmmesize);
  Vector<double> tmpsol(hmmesize);
  double kpd,**kpme=Kp.begin();
  unsigned j;
  for (unsigned i=0; i<nr; i++) {
    if (nc < hmmesize) {
      memset(xy.begin(),'\0',sizeof(double)*hmmesize);
      memcpy(xy.begin(),kpme[i],sizeof(double)*nc);
      if (!hmmec.solve(tmpsol,xy,"ysmp"))  throw exception("Model::CovMat(): hmmec.solve failed");
    }
    else {
      if (!hmmec.solve(tmpsol.begin(),kpme[i],"ysmp"))  throw exception("Model::CovMat(): hmmec.solve failed");
    }
    for (unsigned k=0; k<nr; k++) {
      kp = Kp[k]; sol = tmpsol.begin();
      for (kpd=0.0,j=0; j<nc; j++) kpd += *kp++ * *sol++;
      KvK[k][i] = kpd;
    }
  }
  return(KvK);
}

Vector<double> Model::estimate(const std::string &termname,const doubleMatrix& Kp)
{
   Vector<double> retval;
   if (!get_blupsol()) throw exception("Model::estimate(): bad model");

   int k = term.index(termname,factor_struct);
   if (k < 0) throw exception("Model::estimate(): no such term in the model");

   unsigned nr = Kp.num_rows();
   unsigned nc = Kp.num_cols();
   unsigned i,startaddr;
   if (nc <= term[k].nlevel()*numtrait) {
      startaddr = term[k].startaddr();
      Matrix<double> fullsize_Kp(nr,hmmesize);
      for (i=0; i<nr; i++) {
         memcpy(&(fullsize_Kp[i][startaddr]),Kp[i], sizeof(double)*nc);
      }
      retval.resize(nr);
      estimate((const double **)fullsize_Kp.begin(),nr,hmmesize,retval.begin());
   }
   else {
      throw exception("Model::estimate(): too many columns in Kp");
   }
   return retval;
}




doubleMatrix Model::CovMat(const std::string &termname,doubleMatrix &Kp){
   doubleMatrix retval;
   if (!get_blupsol()) throw exception("Model::CovMat(): bad model");

   int k = term.index(termname,factor_struct);
   if (k < 0) throw exception("Model::CovMat(): no such term in the model");

   unsigned nr = Kp.num_rows();
   unsigned nc = Kp.num_cols();
   unsigned i,startaddr;
   if (nc <= term[k].nlevel()*numtrait) {
      startaddr = term[k].startaddr();
      doubleMatrix fullsize_Kp(nr,hmmesize);
      for (i=0; i<nr; i++) {
         memcpy(&(fullsize_Kp[i][startaddr]),Kp[i], sizeof(double)*nc);
      }
      retval=CovMat(fullsize_Kp);
   }
   else {
      throw exception("Model::CovMat(): too many columns in Kp");
   }
   return retval;
}



void Model::estimate(const double **kpme, const unsigned nr,const unsigned nc, double *kpb)
{
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////
   if (!kpme) {
      warning("Model::estimate(kp): kp is null, thus it is ignored");
      return;
   }
   if (!get_blupsol("ysmp")) throw exception("Model::estimate(): bad model");

   Vector<double> xy;
   if (nc < hmmesize) xy.resize(hmmesize);
   for (unsigned i=0; i<nr; i++) {
      if (nc < hmmesize) {
         memcpy(xy.begin(),kpme[i],sizeof(double)*nc);
         kpb[i] = blupsol.inner_product(xy);
      }
      else {
         kpb[i] = blupsol.inner_product(kpme[i]);
      }
   }
}




double Model::contrast(const Vector<double>& Kp,const double m)
{
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////

   double retval, **kpme = new double *[1];
   kpme[0] = Kp.begin();
   unsigned nc = Kp.size();
   double* M = 0;
   if (m != 0.0) {
      M = new double [1];
      M[0] = m;
   }
   retval = contrast((const double **)kpme,1,nc,M);
   if(kpme){
     delete [] kpme;
     kpme = 0;
   }
   if (M) {
     delete [] M;
     M = 0;
   }
   return retval;
}

double Model::contrast(const doubleMatrix& Kp)
{
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////

   return contrast((const double **)Kp.begin(), Kp.num_rows(), Kp.num_cols());
}

double Model::contrast(const doubleMatrix& Kp,const Vector<double>& M)
{
   ///////////////////////////////////////////
   // trailing zeros in K' can be omitted
   ///////////////////////////////////////////

   int nr = Kp.num_rows();
   if (M.size() != nr ) throw exception(" Model::contrast(): size is incomformable");
   return contrast((const double **)Kp.begin(),nr,Kp.num_cols(),M.begin());
}

double Model::contrast(const std::string &termname,const doubleMatrix& Kp, const Vector<double>& M)
{
   if (!(Kp.begin())) {
      warning("Model::contrast(kp): kp is null, thus it is ignored");
      return 0.0;
   }

   double retval = 0.0;
   if (!get_blupsol()) throw exception("Model::contrast(): bad model");

   int k = term.index(termname,factor_struct);
   if (k < 0) throw exception("Model::contrast(): no such term in the model");

   unsigned nr = Kp.num_rows();
   unsigned nc = Kp.num_cols();
   if (M.size() != nr) throw exception("Model::contrast(): unconformable size");
   unsigned i,startaddr;
   if (nc <= term[k].nlevel()*numtrait) {
      startaddr = term[k].startaddr();
      Matrix<double> fullsize_Kp(nr,hmmesize);
      for (i=0; i<nr; i++) {
         memcpy(&(fullsize_Kp[i][startaddr]),Kp[i],sizeof(double)*nc);
      }
      contrast((const double **)fullsize_Kp.begin(),nr,hmmesize,(const double*)M.begin());
   }
   else {
      throw exception("Model::contrast(): too many columns in Kp");
   }
   return retval;
}



std::string Model::label(const std::string &termname,const unsigned i){
  std::string rawcode;
  if(termname==""){
    trait_effect_level(i,rawcode);
  }
  else{
    int k = term.index(termname,factor_struct);
    if (k < 0) throw exception("Model::label(): no such term in the model");
    unsigned startaddr = term[k].startaddr();
    // cout << startaddr+i << std::endl;
    trait_effect_level(startaddr+i,rawcode);
  }
  return(rawcode);
}


double Model::contrast(const double **kpme, const int nr,
          const int nc,const double *M,const int prt_flag,double **result)
{
   ////////////////////////////////////////////////////
   // H0:  K'*b = M
   // trailing zeros in K' can be omitted
   // REMINDER: blupsol & rellrhs remain intact
   //
   // if M=NULL, then it acts as a vector of zeros.
   // meanp must be declared as double meanp[nr][4]
   /////////////////////////////////////////////////////
   if (!kpme) {
      warning("Model::contrast(kp): kp is null, thus it is ignored");
      return 0.0;
   }
   if (!get_blupsol()) throw exception("Model::contrast(): bad model");
   if (nc > hmmesize) throw exception("Model::contrast(): size not conformable");
   int est=1;
   doubleMatrix kvk(nr,nr);
   doubleMatrix kvkori(nr,nr);
   Vector<double> kpdiff(nr);
   Vector<double> kpb_m(nr);
   Vector<double> xy(hmmesize);
   Vector<double> tmpsol(hmmesize);
   double *sol, kpd;
   const double *kp;
   unsigned i,j,k;
   for (i=0; i<nr; i++) {
      if (nc < hmmesize) {
         memset(xy.begin(),'\0',sizeof(double)*hmmesize);
         memcpy(xy.begin(),kpme[i],sizeof(double)*nc);
         if (!hmmec.solve(tmpsol,xy,"ysmp")) return 0.0;
      }
      else {
         if (!hmmec.solve(tmpsol.begin(),kpme[i],"ysmp")) return 0.0;
      }
      hmmec.mv(tmpsol,xy);
      kp = kpme[i]; sol = xy.begin();
      for (kpd=0.0,j=0; j<nc; j++) {
         kpd = std::max(kpd,fabs(*kp - *sol));
         kp++;
         sol++;
      }
      for (j=nc; j<hmmesize; j++) {
         kpd = std::max(kpd,fabs(*sol));
         sol++;
      }
      kpdiff[i] = kpd;
      for (k=0; k<nr; k++) {
         kp = kpme[k]; sol = tmpsol.begin();
         for (kpd=0.0,j=0; j<nc; j++) kpd += *kp++ * *sol++;
         kvk[k][i] = kpd;
      }
   }
   kvkori = kvk;
   kvk.ginv1(&k);
   if (k < nr) throw exception("Model::contrast(K'): K'V{-1}K is singular. Possible reason:\n"
            " (1) K is non-estimable, and (2) K isn't of full column rank");
   // K'b-M
   for (i=0; i<nr; i++) kpb_m[i] = blupsol.inner_product(kpme[i]);
   if (M) for (i=0; i<nr; i++) kpb_m[i] -= M[i];
   double quq = kvk.quadratic(kpb_m,kpb_m);       // (K'b-m)'(K'(X'V-1X)-K)-1(K'b-m)
   double ssm,sse;
   unsigned rank_x;
   int dfe = 0;
   anova(ssm,sse,rank_x,dfe);
   double res_var_est = sse/dfe;
   int degen = 0;
   int errcode = 0;
   double f_stat=0.0, prob=0.0;
   double sigma_e = residual_var[0][0];
   int chisq = 0;
   if (chisq) {
      f_stat = quq;
      prob = ChiSquare_cdf(f_stat,static_cast<double>(nr),0.0,errcode);
   }
   else {
      if (dfe <= 0 || sse <= SESSION.epsilon) {
         warning("Model::contrast(): either dfe or sse is zero");
         degen = 1;
      }
      else {
         f_stat = (quq*sigma_e/res_var_est)/nr;
         prob = F_cdf(f_stat,static_cast<double>(nr),static_cast<double>(dfe),0.0,errcode);
      }
      //   in univariate cases only: the estimated variance of estimator K'B
      if (!degen)  kvkori *= res_var_est/sigma_e;
   }
   std::string raw_data_code;
   double p_value,var,tcal=0.0;
   int W = SESSION.output_precision;
   std::cout.precision(W);
   if (prt_flag) {
      std::cout << "\n            RESULTS FROM CONTRAST(S)\n";
      std::cout << " ----------------------------------------------------------\n";
      std::cout << "  Contrast MME_addr    K_coef   Raw_data_code\n";
      std::cout << "  ---------------------------------------------\n";
   }
   for (i=0; i<nr; i++) {
      if (prt_flag) {
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
      if (var <= 0.0) throw exception(" Model::contrast(): you have probably found a bug!");
      var = std::sqrt(var);
      tcal = fabs(kpb_m[i]/var);
      p_value = 2.0*(1.0-t_cdf(tcal,static_cast<double>(dfe),0.0));
      if (result) {
         result[i][1] = var;
         result[i][2] = p_value;
      }
      if (!prt_flag) continue;
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
         if (i != nr-1) std::cout << "\n";
      }
   }
   if (prt_flag) {
      std::cout << " ----------------------------------------------------------\n";
      if (nr > 1 && est) {
         std::cout << "   joint hypothesis test H: K'b = M\n";
         std::cout << "   Prob(F > " <<  f_stat << ") = "
              <<  1.0-prob << " (p_value)\n";
      }
   }
   return 1.0-prob;   // p_value
}

void Model::anova(double& ssm,double& sse,unsigned& rank_x,int& dfe)
{
   ////////////////////////////////////////////////////////
   // calculates sums of squares for single trait model
   //  the rank of the random part:
   ////////////////////////////////////////////////////////
   unsigned rank_mme=0,rank_z=0;
   char cl;
   for (unsigned i=0; i<numterm; i++) {
      cl = term[i].classi();
      if (cl == 'P' || cl == 'R') rank_z += term[i].nlevel();
   }
   rank_z = (rank_z - numgroup)*numtrait;
   hmmec.logdet(&rank_mme);
   double ve = residual_var[0][0];
   double sst = yry*ve;
   ssm = blupsol.inner_product(rellrhs)*ve;
   sse = sst - ssm;
   rank_x = rank_mme - rank_z;
   dfe = static_cast<int>(numobs) - static_cast<int>(rank_x);
}

void Model::get_lms_kp(const unsigned termid, const unsigned lvlth,
                       const unsigned yth,Vector<double>& kp, const int *lsmtablevec)
{
   int i,j,k,h;
   unsigned l,t,jb,n,nh,nl,term_start,term_nlevel,endaddr;
   int mcell;

   k = term[termid].startaddr() + lvlth*numtrait + yth;
   int *ia = hmmec.ia();
   int *ja = hmmec.ja();
   double *a = hmmec.a();
   kp.assign(0.0);
   double *kpve = kp.begin();
   kpve--;
   mcell = 0;
   k++;
   endaddr = ia[k+1];
   for (j=ia[k]; j<endaddr; j++) {
      if (ja[j] != k) continue;
      if (a[j] <= SESSION.epsilon) {
         mcell = 1;
         return;
      }
      else {
          break;
      }
   }

   ////////////////////////////////////////////////////////////////
   // set  element kp[i] = 1.0 if MME[k][i] is non-zero element
   //
   ////////////////////////////////////////////////////////////////
   for (j=ia[k]; j<endaddr; j++) kpve[ja[j]] = 1.0;
   for (i=1; i<=hmmesize; i++) {
      endaddr = ia[i+1];
      for (j=ia[i]; j<endaddr; j++) if (ja[j] == k) kpve[i] = 1.0;
   }

   for (t=0; t<numterm; t++) {
      if (!(term[t].classi() == 'F' || term[t].classi() == 'C')) continue;
      if (lsmtablevec[t] >= 1) continue;
      term_start = term[t].startaddr() + 1;
      term_nlevel = term[t].nlevel();
      l = term_start + yth;
      for (i=0; i<term_nlevel; i++) {
         for (n=ia[l+1], j=ia[l]; j<n; j++) {
            if (ja[j] == l) {
               if (fabs(a[j]) > SESSION.epsilon) kpve[l] = 1.0;
            }
         }
         l += numtrait;
      }
   }

   ////////////////////////////////////////////////////////////////////
   // if term[i] is random, then set all elements kp[j] for term[i]
   //    to be zeros,
   // else set elements kp[j] for non-relevent trait to be zero
   ////////////////////////////////////////////////////////////////
   for (i=0; i<numterm; i++) {
      term_start = term[i].startaddr() + 1;
      term_nlevel = term[i].nlevel();
      if (term[i].classi() == 'F' || term[i].classi() == 'C') {
         for (t=0; t<numtrait; t++) {
            if (t != yth) {
               l = term_start;
               for ( j=0; j<term_nlevel; j++) {
                  kpve[l+t] = 0.0;
                  l +=  numtrait;
               }
            }
         }
      }
      else {
         memset(&(kpve[term_start]),'\0',numtrait*sizeof(double)*term_nlevel);
      }
   }

   double dx,xn;
   DataNode dnode;
   for (j=0; j<numterm; j++) {
      if (lsmtablevec[j] == 2) continue;   // j = lsm_term
      term_nlevel = term[j].nlevel();
      term_start = term[j].startaddr() + 1;
      dx = 1.0;
      if (term[j].classi() == 'C') {
         k = factor_struct[term[j].factorindx[0]].index();
         dnode = data->datasheet[k].mean(0);
         if (dnode.missing) throw exception("get_lms_kp(): no valid values in the covariate column");
         dx = dnode.double_val();
      }
      if (!(term[j].classi() == 'F' || term[j].classi() == 'C')) continue;
      if (lsmtablevec[j] == 1) {        // j is related to lsm_term
         n = term_start + term_nlevel*numtrait;
         for (xn=0.0,k=term_start; k<n; k++) xn += kpve[k];
         if (xn >= SESSION.epsilon*10) {
            for (k=term_start; k<n; k++) kpve[k] *= dx/xn;
         }
      }
      else {
         if (term[j].order() == 1) {
            n = term_start + term_nlevel*numtrait;
            for (xn=0.0, k=term_start; k<n; k++) xn += kpve[k];
            if (xn >= SESSION.epsilon*10) {
               for (k=term_start; k<n; k++) kpve[k] *= dx/xn;
            }
         }
         else {
            nl = factor_struct[term[j].factorindx[0]].nlevel();
            nh = term_nlevel/nl;
            for (k=0; k<nh; k++) {
               jb = term_start + k;
               for (xn=0.0,h=0; h<nl; h++) {
                  xn += kpve[jb];
                  jb += nh*numtrait;
               }
               if (xn < SESSION.epsilon*10) continue;
               jb = term_start + k;
               for (h=0; h<nl; h++) {
                  kpve[jb] *= dx/(xn*nh);
                  jb += nh*numtrait;
               }
            }
         }
      }
   }
}

void Model::lsmeans(const std::string &termnames,const std::string &filename,
                    const int io_mode,const int savekp_flag)
{
   ///////////////////////////////////////////////////////////////////////
   // currently not working together with group-model
   // A*B works, but not A * B, we will fix it later
   //
   // REMINDER: blupsol & rellrhs remain intact
   ///////////////////////////////////////////////////////////////////////

   //   lsm_terms.split(nlsm," ");
   std::string sep(" ");
   std::vector<std::string> tmpstr;
   std::string::size_type begidx,endidx;
   begidx = termnames.find_first_not_of(sep);
   while(begidx != std::string::npos) {
      endidx = termnames.find_first_of(sep,begidx);
      if (endidx == std::string::npos) endidx = termnames.length();
      tmpstr.push_back(termnames.substr(begidx,endidx - begidx));
      begidx = termnames.find_first_not_of(sep,endidx);
   }

   unsigned nlsm = tmpstr.size();
   if (nlsm == 0)  {
      warning("Model::lsmeans(termnames): termnames is empty");
      return;
   }
   if (!blup("ysmp")) throw exception("Model::lsmeans(): bad model");

   Matrix<int> lsmtable(nlsm,numterm+1);
   unsigned term_nlevel,i,j,t;
   unsigned nord = 0;
   int k,kk;
   sep = "*";
   std::vector<std::string> effectnames;
   for (term_nlevel = 0, i=0; i<nlsm; i++) {
      k = term.index(tmpstr[i],factor_struct);
      if (k >= 0) {
         if (term[k].classi() != 'F') {
            warning("Model::lsmeans(%s): only fixed(discrete) term allowed",
                    tmpstr[i].c_str());
            lsmtable[i][numterm] = -1; // lsm_term is discarded
            continue;
         }
      }
      else {
         warning("Model::lsmeans(%s): no such term in the model, it's ignored",
               tmpstr[i].c_str());
         lsmtable[i][numterm] = -1;   // lsm_term can't be found in the model
         continue;
      }
      lsmtable[i][k] = 2;          // term k is the lsmeans term
      lsmtable[i][numterm] = k;
      term_nlevel += term[k].nlevel();
      effectnames.clear();
      ////////  effectnames = tmpstr[i].split(nord,"*");
      begidx = tmpstr[i].find_first_not_of(sep);
      while(begidx != std::string::npos) {
         endidx = tmpstr[i].find_first_of(sep,begidx);
         if (endidx == std::string::npos) endidx = tmpstr[i].length();
         effectnames.push_back(tmpstr[i].substr(begidx,endidx - begidx));
         begidx = tmpstr[i].find_first_not_of(sep,endidx);
      }
      nord = effectnames.size();
      for (t=0; t<numterm; t++) {
         if (term[t].order() > 2) {
            lsmtable[i][numterm] = -1;
            throw exception("Model::lsmeans(args): can't handle high order interaction");
         }
         if (lsmtable[i][t] == 2) continue;
         for (j=0; j<nord; j++) {
            kk = term[t].partial_match(effectnames[j],factor_struct);
            if (kk>=0) break;
         }
         if (j<nord) lsmtable[i][t] = 1;
      }
   }
   if (filename == "") {
     std::ofstream ofs;
     ofs.open(filename.c_str(),(OpenModeType)io_mode);
     if (!ofs) throw exception("Model::lsmeans(): cann't open or already exist");
     out_lsmeans_to_stream(ofs,nlsm,lsmtable,savekp_flag);
     info(ofs);
     ofs.close();
   }
   else {
        out_lsmeans_to_stream(std::cout,nlsm,lsmtable,savekp_flag);
   }
}

void Model::out_lsmeans_to_stream(std::ostream& ofs, const unsigned nlsm,
                                  const Matrix<int> &lsmtable,const int savekp)
{
   std::string tmpfname;
   std::fstream tmpfile;
   if (savekp) {

      tmpfname = SESSION.mktemp();
      tmpfile.open(tmpfname.c_str(),std::ios::out);
      if (!tmpfile) throw exception("Model::out_lsmeans_to_stream(): can't open");
   }
   int kk;
   unsigned term_nlevel,startaddr,i,j,k,t,tc,nord;
   int max_levelwidth=0;
   int *levelwidth = new int [1];
   std::string lsm_terms;
   std::string level_name;
   doubleMatrix result(1,4);
   Vector<double> kp(hmmesize);
   double **kpme = new double *[1];
   kpme[0] = kp.begin();
   std::cout.precision(SESSION.output_precision);

   ofs <<  "\n                       LEAST SQUARES MEANS\n";
   ofs << " ------------------------------------------------------------"
           "--------------\n";
   for (i=0; i<nlsm; i++) {
      kk = lsmtable[i][numterm];
      if (kk == -1) continue;
      lsm_terms = "";
      nord = term[kk].order();
      for (t=0; t<nord; t++) {
         lsm_terms.append(factor_struct[term[kk].factorindx[t]].name());
         if (t+1 < nord) lsm_terms.append("*");
      }
      term_nlevel = term[kk].nlevel();
      startaddr = term[kk].startaddr();
      for (t=0; t<numtrait; t++) {
         if (numtrait > 1) {
            ofs << "  trait: "<< trait_struct[t]->name()<< "\n";
            ofs << "  ------\n";
         }
         ofs << " " << std::setw(20) << lsm_terms;
         ofs << "      lsmean  sqrt(Var(kb))     p_value  estimability\n";
         tc = startaddr + t;
         for (j=0; j<term_nlevel; j++) {
            get_lms_kp(kk,j,t,kp,lsmtable[i]);
            contrast((const double **)kpme,(unsigned)1,hmmesize,
                     0,0,result.begin());
            trait_effect_level(tc,level_name,3);
            if (savekp) {
               levelwidth[0] = level_name.size()+1;
               if (*levelwidth > max_levelwidth) max_levelwidth = *levelwidth;
               tmpfile.write((char *)levelwidth,sizeof(int));
               tmpfile.write(level_name.c_str(),(size_t)levelwidth[0]);
               tmpfile.write((char *)kp.begin(),sizeof(double)*hmmesize);
            }

            tc += numtrait;
            ofs << " " << std::setw(20) << level_name;
            if (result[0][3] >= 1.0e-5) {
               ofs << "          *** NON-ESTIMABLE ***";
            }
            else {
               ofs << " " << std::setw(12) << result[0][0]
                    << " " << std::setw(12) << result[0][1]
                    << " " << std::setw(12) << result[0][2];
               if (result[0][3] > 1.0e-12) {
                  ofs  << " " << std::setw(12) << result[0][3];
               }
            }
            ofs << "\n";
         }
      }
      if (i+1 < nlsm) ofs << "\n";
   }
   ofs << " ------------------------------------------------------------"
           "--------------\n";
   ofs << " note 1: estimability = Max(mme*ginv(mme)*k - k);\n";
   ofs << "         if it is non-zero, lsmean may or may not be estimable.\n";
   ofs << " note 2: sqrt(Var(kb)) is also called STD ERR (see SAS)\n";
   if (kpme) {
     delete [] kpme;
     kpme = 0;
   }

   if (savekp == 0) {
     if (levelwidth){
       delete [] levelwidth;
       levelwidth = 0;
     }
      tmpfile.close();
      return;
   }
   char *levelname_holder; 
   if(max_levelwidth>0){
     levelname_holder = new char [max_levelwidth];
   }
   else {
     levelname_holder = 0;
   }
   tmpfile.close();
   tmpfile.open(tmpfname.c_str(),std::ios::in);
   if (!tmpfile) throw exception("Model::out_lsmeans_to_stream(): can't open");
   ofs << "\n";
   ofs << "             COEFFICIENTS FOR LEAST SQUARES MEANS\n";
   ofs << " -----------------------------------------------------------\n";
   for (i=0; i<nlsm; i++) {
      kk = lsmtable[i][numterm];
      if (kk == -1) continue;
      lsm_terms = "";
      nord = term[kk].order();
      for (t=0; t<nord; t++) {
         lsm_terms.append(factor_struct[term[kk].factorindx[t]].name());
         if (t+1 < nord) lsm_terms.append("*");
      }
      term_nlevel = term[kk].nlevel();
      startaddr = term[kk].startaddr();
      for (t=0; t<numtrait; t++) {
         tc = startaddr + t;
         for (j=0; j<term_nlevel; j++) {
            tmpfile.read((char *)levelwidth,sizeof(int));
            tmpfile.read(levelname_holder,(size_t)levelwidth[0]);
            tmpfile.read((char*)kp.begin(),sizeof(double)*hmmesize);
            if (numtrait > 1) {
               ofs << "   for " << trait_struct[t]->name() << ":"
                   << lsm_terms << ": " << levelname_holder << "\n";
            }
            else {
               ofs << "   for "<< lsm_terms << ": " << levelname_holder << "\n";

            }
            ofs << "   -----------------------------------------------------\n";
            ofs << "      mme_addr   coefficient               raw_data_code\n";
            for (k=0; k<hmmesize; k++) {
               if (fabs(kp[k]) > 0.0) {
                  trait_effect_level(k,level_name,2);
                  ofs << " " << std::setw(12) << k
                      << " " << std::setw(12) << kp[k]
                      << " " << std::setw(29) << level_name.c_str() << "\n";
               }
            }
            ofs << "\n";
         }
      }
      if (i+1 < nlsm) ofs << "\n";
   }
   if(levelwidth) {
     delete [] levelwidth;
     levelwidth = 0;
   }
   if (levelname_holder) {
     delete [] levelname_holder;
     levelname_holder = 0;
   }
   tmpfile.close();
}
} ///////////// end of namespace matvec
