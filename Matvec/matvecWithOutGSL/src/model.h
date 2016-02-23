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

#ifndef MATVEC_MODEL_H
#define MATVEC_MODEL_H

#include <string>
#include <set>
#include <vector>
#include "session.h"
#include "data.h"
#include "geneticdist.h"
#include "termlist.h"
#include "minimizer.h"
#include "fpmm.h"
#include "sparsematrix.h"
#include "rsamplerparms.h"

namespace matvec {
class Population;
class Pedigree;

/*!
 * \struct pos_val_node
 * \brief used to store mme positions, x-values, and trait values for one obs.
 */

 struct pos_val_node{
   Vector<unsigned> pos_term;
   Vector<double> xval_term, trait_vec;
 };


/*!
 * \struct Parm_Elem mvmodel.h
 * \brief It is used in Model
 */
struct   Parm_Elem {
         double*  Value;
         double   Max;
         double   Min;
         std::string Name;
 };
//BRS

/*!
 *  \class Model  mvmodel.h
 *  \brief a generalized linear mixed model
 *
 * term[t].classi =
 *  <TABLE>
 *  <TR><TD>'F'</TD><TD>a fixed term, including intercept</TD></TR>
 *  <TR><TD>'C'</TD><TD>a covariate term</TD></TR>
 *  <TR><TD>'R'</TD><TD>a regular random term</TD></TR>
 *  <TR><TD>'P'</TD><TD>a random term associated with a pedigree</TD></TR>
 *  </TABLE>
 * \verbatim
  factor_struct[i].classi = 'F', a fixed factor
                          = 'I', intercept
                          = 'C', covariate factor
                          = 'P', a factor associated with a pedigree

  weightinverse  = 0, no inverse is needed (default)
                 = 1, inverse is needed
 * \endverbatim
 *
 *  \sa Data
 */

class Model : public Minimizer {
   protected:
  std::vector<pos_val_node> pos_val_vector;
      std::string   modelfname,weightname,modelstring, *traitname;
      std::string   modelstringstr;
      char         *modelstr; // Change to keep model information in memory
      int          modelpcount; // Change to keep model information in memory

      unsigned     numtrait,numterm,numfact,maxorder,numrec,numcol,popsize;
      unsigned     numobs,hmmesize,non_zero,max_nz,npattern,ntermGdist;
      unsigned     *rawpos,*kvec,numgroup;

      int          weightinverse,data_prepared,pop_created, categorical;
      std::vector<std::string> pattern;
      std::set<std::string>    pattern_tree;          // pattern points to patterrn_tree
      UnknownDist  *unknown_prior;

      FieldStruct   *factor_struct,**trait_struct;
      HashTable    *xact_htable, **idlist;

      Vector<double>       lnr0vec, blupsol, rellrhs, diag_mme, mme_times_res;
      double       yry;

      unsigned     nfunk_in_dfreml,dfreml_called;
      double       reml_value;
      std::string   dfreml_method;

      SparseMatrix hmmec;
  //////////////////// Stuff for  Correlated random effects //////////
  Vector<int>  base_effect;
     // defaults to [0...numterm-1]
     // changed first correlated effect  for later correlated effects
  Vector<int> corr_link;
     // defaults to -1
     // points next correlated term
  Vector<int> pos_vec;
     // Position in link list;
     // defaults to 0
  Vector<Vector<int> > trait_map;
     // Maps traits within correlated effects
     // defaults to Vectors of length numtrait with element i
  Vector<int>  nt_vec;
     // defaults to numtrait for uncorrelated traits
     // numtrait*ncorr for first correlated trait
     // 0 for later correlated traits
  Vector<int> var_link;
     // defaults to [0..numterm-1]

     //////////////////  working spaces //////////////////
      unsigned *pos_term;
      double *xval_term, *trait_vec;

     ///////////////// for gibbs /////////////////////////
      unsigned     *rec_indid;

      virtual void  setup_ww(Vector<double> *rhs);
      void     add_Ig(int t, Vector<double>* x_p=0);
      void     add_Ig_diag(int t, double** M);
      void     add_Ag(int t, Vector<double>* x_p=0);
      void     add_Ag_diag(int t, double** M);
      virtual void add_G_1(void);
      virtual void add_G_1_diag(double** M);
      void     add_G_1_single_trait(const Vector<double> &ratio);

      double   setup_ww_single_trait(const Vector<double> &ratio,const unsigned pt,
                                     const unsigned ibeg,Vector<double> *zz);
      void     compute_rhs_diag_mme(double** M);
      void     mme_times(Vector<double> &x);
      int      build_model_term(const std::string &modelspecs);
      void     inverse_residual_var(void);
      void     assign_id_xact(const Vector<bool>&);
      void     hashxact(const Vector<bool>&);
      void     re_hash_data(Vector<bool>&);
      virtual void save_pos_val(Vector<bool>&, const std::string &solver="ysmp");
      void     input_pos_val(std::istream& modelfile);
      void     input_pos_val_iod(std::vector<pos_val_node>::iterator it);
      void     get_lms_kp(const unsigned termth, const unsigned lvlth,
                       const unsigned yth,Vector<double>& kp, const int *lsmtablevec);
      void     out_lsmeans_to_stream(std::ostream& ofs,const unsigned nlsm,
                                     const Matrix<int> &lsmtable,const int savekp=0);

   public:
      unsigned     num_ped,act_numtrait;
      unsigned     data_static;
      RSamplerParms myRSamplerParms;
      enum ModelType {bad_model, fixed_model, mixed_model, reg_model};
      ModelType    type; // This changed from model_type

      Population   *pop;
      Data         *data;

      Matrix<double>  *rve;
      doubleMatrix     residual_var;
      Vector<double>   lng0vec;
      TermList  term;

      Model(void) {initialize();}
      Model(const std::string &modelspecs) {initialize(); equation(modelspecs);}
      Model(const Model& M):Minimizer() {initialize(); copyfrom(M);}
      virtual ~Model(void){release();}

      int      display(void) const;
      virtual int prepare_data(const std::string &solver="ysmp");
      unsigned nonzero(void) const {return non_zero;}
      unsigned mmesize(void) const {return hmmesize;}
      unsigned ntrait(void) const {return numtrait;}
      unsigned nterm(void) const {return numterm;}
      unsigned totalnvc(void) const;

      void     trait_effect_level(const unsigned mmeaddr, std::string& retval,
                                  int info_flag = 1);
      void     initialize(void);
      void     copyfrom(const Model& M);
      void     release(void);
      void     nonzero(const unsigned nz) {max_nz = nz;}
      void     prior_dist(const std::string &termname,Pedigree& P,GeneticDist *D);
      void     prior_dist(const std::string &termname,GeneticDist *D);
      void     DGSamplerSetup(unsigned numLoci,Data *D);
      void     DGSampler(unsigned PQTL, bool &initSamperToOrderFounders, unsigned numOfSamples,unsigned numOfSL,unsigned numOfHaplo,unsigned numOfSLCas,char* fname); //adapted by CS, Dec. 04/Jan 05
      void     RSamplerPrior_dist(const std::string &termname,Pedigree& P, Data *D,GeneticDist *G);
	  void     RSampler(RSamplerParms rsParms);
	  void     RSampler(std::string fileName);
      void     RSamplerInitialDG(const std::string &samplerType="genotypic",const std::string &whatToCompute="probDescHaplo");
      void     RSampler(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType,const std::string &howToSample, const std::string &samplerUsed, const std::string &whatToCompute); 
      void     RSamplerGibbs(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType="genotypic",const std::string &whatToCompute="genotypeProb");
      void     RSamplerMH(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType="genotypic",const std::string &whatToCompute="genotypeProb");
      void     RSamplerMH(unsigned numLoci,unsigned numOfSamples,const std::string &samplerType="allelic",const std::string &whatToCompute="genotypeProb");
      void     RSamplerGibbsMH(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,unsigned stepGibbsMH,const std::string &samplerType="genotypic",const std::string &whatToCompute="genotypeProb");
	  double   logOldProposal, logNewProposal, logOldTarget, logNewTarget;
	  double   acceptedPosition;
      void     variance(const std::string &termname,const doubleMatrix& v);
      void     variance(const std::string &termname,Pedigree& P,const doubleMatrix& v);

      void     variance(const std::string &termname,const double v00,...);
      void     variance(const std::string &termname,Pedigree& P,const double v00,...);

      void     VarLink(const std::string &termname,const std::string &termname2);
      void     covariance(const std::string &termname,const std::string &termname2,const doubleMatrix& v);
      void     covariance(const std::string &termname,const std::string &termname2,const double v00,...);



      void     weight(const std::string &wtname,const int r=0)
                                      {weightname= wtname; weightinverse=r;}
      void     covariate(const std::string &factorname);
      int      equation(const std::string &modelspecs);
      void     fitdata(Data& D);
      void     info(std::ostream& stream) const;
      void     var2vec(Vector<double> &x);
      void     vec2var(const Vector<double> &x);
      void     release_mme(void){hmmec.resize(0,0);}
      void     vce_gibbs_sampler(const unsigned nvc,Vector<double> &vc,
                                 const Vector<double> &ratio,const Vector<unsigned> &nlevel,
                                 const double yy,double& ee,Vector<double> &uAu,
                                 Vector<double> &bu,Vector<double> &wspace,const Vector<double> &zz);
      void     update_mme(const Vector<double> &ratio,const Vector<double> &oldratio,const Vector<double> &zz);
      void     compute_xbzu(Vector<double>& bu);
      Vector<double> get_solutions(const std::string &termname);
      void     save(const std::string &fname,
                    const int io_mode=std::ios::out/*|std::ios::noreplace*/);
      void     sampleW(Vector<double> &sol,Vector<double> &newrhs);
      void     vce_emreml_single_trait(const int miter=1000,const int intive=1,
                                       const double stoptol=1.0e-5);
      void     vce_emreml_multi_trait(const int miter=1000,const int intive=1,
                                      const double stoptol=1.0e-5);
      void     uAu_trCA(Vector<double> &tsol,Vector<double> &uAu, Vector<double> &trca,Vector<double> &ivect,
                        Vector<double> &invect,const Vector<double> &ratio,const Vector<double> &zz);
      Vector<double>*  blup(const std::string &method = "ysmp",double stopval=0.001,
                            double relax=1.0, unsigned mxiter=100000);
      void             blup_pccg(double stopval=1.0e-8,unsigned int mxiter=100000);
      Vector<double>*  blue(const std::string &method = "ysmp",double stopval=0.001,
                            double relax=1.0, unsigned mxiter=100000);
      Vector<double>*  get_blupsol(const std::string &method = "ysmp",double stopval=0.001,
                                   double relax=1.0, unsigned mxiter=100000);
      Vector<double>*  vce_emreml(const int miter=1000,const int intive=1,
                                  const double stoptol=1.0e-4);
      Vector<double>*  vce_gibbs(const unsigned warmup=2000,
                                 const unsigned gibbslen=5000);
      Vector<double>*  vce_dfreml(const std::string &method = "powell",int maxiter=500000,
                                  double funtol=1.0e-5);
      Vector<double>*  genotypic_value_peeling(void);
      Vector<double>*  genotypic_value_gibbs(const unsigned warmup=30,
                                     const unsigned gibbslen=1000);
      Vector<double>*  genotypic_value_cat(const unsigned warmup,
                           const unsigned gibbslen);
      double   genotype_dist_gibbs(const unsigned warmup=30,
                                   const unsigned gibbslen=1000);
      void     genotype_dist_peeling(const int prtlevel=1,const int estifreq=1);

      double   restricted_log_likelihood(void);
      double   log_likelihood_peeling(const unsigned maxit=3);
      double   log_likelihood_gibbs(const unsigned wup=30,
                                    const unsigned glen=1000);
      double   minfun(const Vector<double> &x, const int n);
      double   minfun_peeling(const Vector<double> &x, const int n);
      double   minfun_vce(const Vector<double> &x,const int n);
      double   estimate_peeling(void); // new version of ml_estimate_peeling?
      double   ml_estimate_peeling(void); // old version
      double   estimate(const Vector<double>& Kp);
      Vector<double>   estimate(const doubleMatrix& Kp);
      Vector<double>   estimate(const std::string &termname,const doubleMatrix& Kp);
      void     estimate(const double **kpme, const unsigned nr,const unsigned nc, double *kpb);
  doubleMatrix CovMat(doubleMatrix& Kp);
  doubleMatrix CovMat(const std::string &termname,doubleMatrix& Kp);
      virtual double   contrast(const Vector<double>& Kp, const double m=0.0);
      virtual double   contrast(const doubleMatrix& Kp, const Vector<double>& M);
      virtual double   contrast(const doubleMatrix& Kp);
      virtual double   contrast(const std::string &termname,const doubleMatrix& Kp,const Vector<double>& M);
      virtual double   contrast(const double **kpme, const int nr,
                        const int nc,const double *M = 0,
                        const int prt_flag=1,double **result = 0);
      void     anova(double& ssm,double& sse,unsigned& rank_x,int& dfe);
      void     lsmeans(const std::string &termnames,const std::string &filename = "",
                       const int io_mode = std::ios::out /*| std::ios::noreplace*/,
                       const int savekp_flag=0);

      virtual SparseMatrix*  setup_mme(Vector<double> *rhs = 0);
      SparseMatrix*  mmec(void){return &hmmec;}
      Vector<double>*        rhs(void) {return &rellrhs;}

//      ModelType type(void) const {return type;} // check this is not found elsewhere
      //BRS
  int Q_pos, maxit, map_parm; //BRS
  char map_function; //BRS
  double new_distance;
  doubleMatrix RecoTable;
  double multipoint(int maxit);
  void multipoint_setup(void);
  void multipoint_setup(Fpmm& F, char map_fun='H', int map_par=2); //BRS
  void multipoint_tables(void);
  void multipoint_QTL_table(void);
  void multipoint_Rec_table(void);
  void multipoint_Rec_recalc(double distance);
//  double MapD(double rate); //BRS
//  double MapR(double dist); //BRS
   double MapF(double dist, int qr=0); //BRS
  double minfun_multipoint(const Vector<double> &x, const int n);
  double ProbGamete (Vector<int> value, int nloci);
//  void multipoint_profile(int method, int maxit,double Min, double Max, double step_size, int Print);
  void multipoint_profile_distance(int method, int maxit,double Min, double Max, double step_size, int Print, const std::string &fname); 
  double multipoint_ml_estimate(const int method=0, int niter=0, double tol=1.0e-16, double epsilon=1.0e-8, int printlevel=0);
  int N_Parm_List;
  void multipoint_estimate_Ve(double Initial, double Max, double Min);
  void  multipoint_estimate_Distance(double Initial, double Max, double Min);
  void multipoint_estimate_Qmeans(Vector<double> Initial, Vector<double> Max, Vector<double> Min , char  mtype='G');
  void multipoint_estimate_Qfreq(double Initial, double Max, double Min); 
  Vector<Parm_Elem> Parm_List;
  void multipoint_genot_probs(void); 
  void multipoint_compute_Rec_table(void);
  double multipoint_mixture_approx(int maxit=0, double P_var=1.0);
  void multipoint_display_parms(void);
  void multipoint_estimate_PolyV(double Initial, double Max, double Min);
  unsigned   get_popsize(void) const {return popsize;};
  //BRS
  void descent_graph_peeling_initialisation(void);
  Vector <int> get_gamete_vector(unsigned index);
  void calc_gamete_prob(void);
  Vector <double> map_distance;
  doubleMatrix gamete_prob_table;
  void map(void);
  doubleMatrix Cr_table;
  double graph_log_likelihood_peeling(const unsigned maxiter);
  std::string label(const std::string &termname,const unsigned i);
  void DGSamplerIBDMatrix(unsigned numOfSamples, unsigned numSL, unsigned numHaplo, unsigned numSLCas);
};
}
#endif
