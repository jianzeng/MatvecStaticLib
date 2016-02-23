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

#ifndef MATVEC_GLMM__H
#define MATVEC_GLMM__H

namespace matvec {
class Observe {
public: 
  void **param;
  int numtrait,rec,pev_current;
  int need_data;
  std::string pattern;
  Vector<double> residual;
  Vector<double> estimate;
  doubleMatrix pev;
  doubleMatrix rve;
  Vector<double> rpy;
  doubleMatrix W,Adj_Cov; // {sum_j (y_{ij}-u_{ij})_j {\partial R^{-1}_iH_i\over \partial \eta_k}
  Vector<double> q; //q is the q in \ell = y'q+....
  Vector<double> est_null; // Estimate under HO
  double weight;
  double like_val;
  Vector<double> y;
  doubleMatrix Resid_Var;
  doubleMatrix H;
  doubleMatrix Hinv;
  int nvc;
  doubleMatrix *Resid_Sin_Mat;
  doubleMatrix *resid_var;
  Vector<double> *varcomp;
  Data *data;
  Observe() {pattern="";nvc=0;Resid_Sin_Mat=0;varcomp=0;
  pev_current=0;like_val=0;data=0;param=0;};
  ~Observe() {
        if(Resid_Sin_Mat) delete [] Resid_Sin_Mat;Resid_Sin_Mat=0;
    };
  int display(void) const {std::cout<<"\tan object of Observe\n";return 1;};
} ;






class GLMM : public Model{
public:
  virtual        double   contrast(const Vector<double>& Kp, const double m=0.0)
  ;
  virtual        double   contrast(const doubleMatrix& Kp, const Vector<double>& M)
  ;
  virtual        double   contrast(const doubleMatrix& Kp)
  ;
  virtual        double   contrast(const std::string &termname,const doubleMatrix& Kp,const Vector<double>& M)
  ;
protected:
   virtual void     setup_ww(Vector<double>* rhs);

public:
     int      equation(const std::string &modelspecs);
      virtual double   contrast(const double **kpme, const unsigned nr,
                        const unsigned nc,const double *M=0,
                        const int prt_flag=1,double **result=0);

protected:
  double Marq_lambda;
  void **param;
  double numy,ldet;
  int aireml_called;
  int Need_Ainv, Need_Residual,res_nvc,Return_Stat,do_west,num_west,AI_sample;
  int Print_Level,PEV_Scale,Use_Like,Print_Summary,Update_estimate,Sand_Calc,LM,Need_PEV;
  Observe * observation;
  doubleMatrix  * SinMat;
  doubleMatrix  * varinv, *Var;
  doubleMatrix  **corrvar;
  doubleMatrix  CorrVar;
  int     ncorr,*corrmap;
  doubleMatrix  Info;
  Vector<double> Score;
  Vector<double> varcomp,varest;
  Vector<int> kvec;
  double like_val,like_adj;
  doubleMatrix Fmat,FRHS,FSol;
  SparseMatrix hainv;
  SparseMatrix hSand,hInv;
  SparseMatrix itilde,ihat,jtilde,jhat,S;
  Vector<double> sol_null,qhat;
  void (*link_function)(Observe &);
  Vector<double> Initial;
  int Asym,LR;
public:
  int nrandom;  ////// act_numtrait added by TW
  int num_glmm;  //// Number of iterations of glim during each round of ai_reml
  int *nt_vec2; 
  GLMM() {observation = 0;SinMat=0;Need_Ainv=1;varinv=0;like_val=0;Need_PEV=0;
  Var=0;link_function=0;Need_Residual=0;res_nvc=0;Print_Level=0;Print_Summary=0;Asym=0;LR=0;PEV_Scale=1;Use_Like=0;Sand_Calc=0;LM=0;Return_Stat=0;do_west=0;num_west=0;numy=0;
  param=0;Marq_lambda=1.;Update_estimate=1;ncorr=0;corrmap=0,corrvar=0;nrandom=0;nt_vec2=0;AI_sample=0;aireml_called=0;num_glmm=3;};
  ~GLMM()
    {release_obs();
    release_SinMat();
    int i;if(varinv){
      for(i=0;i<numterm;i++)varinv[i].resize(0,0);delete [] varinv; varinv=0;
    }
    if (Var){
      for(i=0;i<numterm;i++)Var[i].resize(0,0); delete [] Var;Var=0;
    }
    if(corrvar) {
      for(int i=0;i<numterm;i++) {
	delete [] corrvar[i];
      }
      delete corrvar;corrvar=0;;
    }
    if(corrmap) delete [] corrmap;corrmap=0;
    if(nt_vec2) delete [] nt_vec2;nt_vec2=0;
    };
  void SSQCP(void);
  void **Param(void **par=0) {if(par) param=par;return(param);}
  int Large_Sample(int Sample=-1) {if(Sample != -1) Asym=Sample;return(Asym);}
  int Pev_scale(int Sample=-1) {if(Sample != -1) PEV_Scale=Sample;return(PEV_Scale);}
  int LR_Test(int LR_set=-1) {if(LR_set != -1) LR=LR_set;return(LR);} //LR=0 WALD LR=1 LR LR=2 LM
  int LM_Test(int LM_set=-1) {if(LM_set != -1) LM=LM_set;return(LM);}
  int Data_static(int DS_set=-1) {if(DS_set != -1) data_static=DS_set;return(data_static);}
  int Like_LR(int L_set=-1) {if(L_set != -1) Use_Like=L_set;return(Use_Like);}
  int Westell(int W_set=-1) {if(W_set != -1) do_west=W_set;return(do_west);}
  int AI_Sample(int S_set=-1) {if(S_set != -1) AI_sample=S_set;return(AI_sample);}
  void Initial_eta(Vector<double> initial){Initial=initial;}
  void residual(int get_pev=0);
  void new_SinMat(void);
  double log_likelihood(int solve=1);
  double restricted_log_likelihood(int solve=1);
  void release_SinMat(void);
  void Build_SinMat(void);
  int print_level(int level=-99) {if(level !=-99) Print_Level=level; return(Print_Level);}
  int return_stat(int level=-99) {if(level !=-99) Return_Stat=level; return(Return_Stat);}
  int print_summary(int level=-99) {if(level !=-99) Print_Summary=level; return(Print_Summary);}
  Observe* get_obs(void) {return observation;}
  doubleMatrix* get_SinMat(void) {return SinMat;}
  Vector<double>& resid(int i){return observation[i].residual;}
  Vector<double>& Linear_Predictor(int i){return observation[i].estimate;}
  doubleMatrix& pev(int i) {return observation[i].pev;}
  double log_like(void) {return like_val;}
  doubleMatrix INFO(void) {return Info;}
  Vector<double> SCORE(void) {return Score;}
  Observe& new_obs();
  void release_obs();
  SparseMatrix& ainv(void);
  unsigned Numrec(void) {return numrec;};
  doubleMatrix AI_REML(int numiter=10,double tol=1.e-4,double info_scale=0);
  virtual      SparseMatrix*  setup_mme(Vector<double> *rhs=0);
  void link(void (*link_fn)(Observe &),int rvc=0) {link_function=link_fn;if(res_nvc!=rvc){res_nvc=rvc,varcomp.resize(res_nvc);}}
  unsigned TotalNvc(void) const;
  int Var2Vec(double *x);
  void Vec2Var(const double *x);
  void covariance_old(const std::string &termname,const std::string &termname2,const doubleMatrix& v);
  Vector<double>* glim(int iterations=10,double **ke=0, int nr=0,int nc=0,double *mean=0);
  void add_G_Sand(void);
  void add_IgSand(int t);
  void add_AgSand(int t);
  void add_Ag(int t,int ct);
  void add_Ag(int t, SparseMatrix &h);
  void add_Ag_old(int t);
  void add_Ig(int t, SparseMatrix &h);
  void add_G_1(SparseMatrix &h);
  void add_G_1_old();
  void build_CorrVar(void);
  void build_hInv(void);
  void info(std::ostream &stream);
  void save(const std::string &fname, const int io_mode=std::ios::out);//SDK |std::ios::noreplace);
#ifdef GLMM_QUAD
  double quad(double *x, double *y, double *xsol, double *ysol, int path=3);
#endif
  void normalLog(void);
  doubleMatrix getVarEstimates(const std::string &termname);
};

doubleMatrix *KP(Data *D,doubleMatrix &SKM,const std::string &lpl,const std::string &censor);

doubleMatrix operator*(Vector<double> &a,Vector<double> &b);
//Vector<double> operator*(Vector<double> &a,Vector<double> &b);

void bound_pd (double *varnew, int ntrait, doubleMatrix & P);
} ///////// end of namespace matvec
#endif
