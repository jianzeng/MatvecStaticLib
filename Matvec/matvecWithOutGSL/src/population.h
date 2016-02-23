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

#ifndef Population_H
#define Population_H

///////////////////////////////////////////////////////////////////
// (1) pm_storage remains original table-index all the time
// (2) hashtable.get_id() returns the table-index of pm_storage
// (3) popmember is changed along with individual id such that
//     popmember[i]->id == i+1.
// April, 1993, University of Illinois
///////////////////////////////////////////////////////////////////
#include <string>
#include <sstream>
#include <iostream>
#include "hashtable.h"
#include "pedigree.h"
#include "rpedigree.h"
#include "data.h"
#include "individual.h"
#include "mim.h"
#include "nufamily.h"
#include "chromosome.h"
#include "gnodederived.h"
#include "gnodesetderived.h"
#include "sparsematrix.h"
#include "fpmm.h"
#include "matrix.h"
#include "safe_vectors.h"

namespace matvec {

  //  typedef double (*penetranceType)(const Individual*, const double **);
  /*!
    \class Population  Population.h
    \brief a genetic population consisting of individuals

    \sa Individual Pedigree
  */

  class IndexVector{
  public:
    IndexVector(){
		maxElementVector.name = "IndexVector::maxElementVector";
		multiplicationCode.name = "IndexVector::multiplicationCode";
	}	
    unsigned size;
    SafeSTLVector<unsigned> maxElementVector;
    SafeSTLVector<unsigned> multiplicationCode;
    void setupVectors(SafeSTLVector<unsigned> aVector);
    SafeSTLVector<unsigned> getVector(unsigned index);
    unsigned getIndex(std::vector<unsigned> stateVector);
  };

  class Model;
  class Population {
    friend class Model;
    friend class GLMM;
    friend class Individual;
    friend class NuFamily;
   public:
   
   Pedigree  *myPedPtr;
   RPedigree *myRPedPtr;

    Matrix<float> ibdMatrix;

    int          stdid,mark_need_remove,sex_confirmed;
    int          offs_info_built,spouse_info_built,fdone;
	bool         nuFamiliesDone;

    double       scaling_val;
    NuFamily     **nufamily_vec;
    unsigned     num_nufamily,num_t_nufamily,nonloop_t_nufamily;
    unsigned     tn_genotype, tn_gamete, maxnamelen, tn_qtl,  n_markerLoci;
    doubleMatrix       *trans_mat;

    unsigned     popsize, maxpopsize, numchrom, numtrait;
    Individual   *pm_storage, **popmember;
    Vector<double>       blupsol;
    std::string     evaluate_method;
    Genome       *pop_gamete, *gametebase;
    static GeneticDist *prior;
    int          nmarker;
    std::string*    markersymbol;
    penetranceType penetrance_f;

    SparseMatrix hmmec;

    void     mark_descendant_recursive(const Individual* ind,const unsigned
				       counter=1);
    void     mark_ancestor_recursive(const Individual* ind,const unsigned
				     counter=1);
    void     mark_relative_recursive(Individual* ind,const unsigned
				     counter=1);
    void     partial_iterative_peeling(Individual *I,doubleMatrix& wspace);
    void     setPenetrance(penetranceType p){penetrance_f=p;}

    HashTable hashtable;
    doubleMatrix    residual_var;
    Vector<double>    mean_for_genotype;

    Population(void);
    Population(GeneticDist *D);
    Population(const int maxsize,GeneticDist *D);
    Population(const Population& A);
    ~Population(void){release();}

    const Population&  operator=(const Population& A);

    void        copyfrom(const Population& A);
    void        initialize(void);

    void        inbcoef_quaas(void);
    void        inbcoef_meuwissen(void);
    void        confirm_sex(const int col_sex);
    doubleMatrix      rela(void);
    doubleMatrix      relv(const double er);
    doubleMatrix      relv_nomissing(const double er);
    double      cprob_op(const unsigned goff[],const unsigned gparent[],
			 const int ori);
    double      llhood_genotype(const unsigned c_id,const unsigned l_id);
    double      llhood_phenotype(const unsigned nrep);
    doubleMatrix      relv_inv(const double er);
    void        allele_freq(const char type[], const unsigned nal,...);
    void        penetrance_ml(const char fname[]);
    void        release(void);
    void        resize(const unsigned maxn,GeneticDist *D);

    void        build_offs_info(void);
    void        build_spouse_info(void);
    void        build_trans_mat(void);
    void        renum(void);
    void        display(const char key[]);
    //    unsigned    input(const char fname[],const char recfmt[]);
    unsigned    input_ped(const char fname[],const char recfmt[],
			  GeneticDist *D);
    unsigned    input_ped0(matvec::Pedigree& P, GeneticDist *D);
	void        input_ped0(GeneticDist& D);
    void        input_ped(matvec::RPedigree& P, GeneticDist& D);
    unsigned    input_ped(matvec::Pedigree& P, GeneticDist *D);
    unsigned    input_data(Data* D);
    unsigned    input_data(const char fname[], GeneticDist* G);
    unsigned    size(void) const {return popsize;};
    unsigned    nbase(void) const;
    unsigned    get_id(const char name[]) {return hashtable.get_id(name);}
    unsigned   input_markerData(Data* D);
	void       input_markerData(MIM& M);
    unsigned   input_descentGraph(char *dgfile);
    void       output_descentGraph(char *dgfile);
	void       output_descentGraphJP(string dgfile);

    const char* ind_name(const unsigned id) const;
    void        compute_pdm(const Individual* ind,Matrix<double> &pdm);
    Vector<double>      inbcoef(void);
    void        setup_m_ww(SparseMatrix& mme,Vector<double>& rhs,Vector<int> &start_addr,const double ev_e);
    void        add_Ga_inv(SparseMatrix& mme, const int start_addr,
			   const double ev_a);
    void        add_Gv_inv(SparseMatrix& mme,const int start_addr,
			   const double ev_v, const double er);
    SparseMatrix*  setup_m_mme1(Vector<double>& rhs,const double ev_v,
				const double ev_u, const double ev_e,
				const double er);
    SparseMatrix*  setup_m_mme(Vector<double>& rhs,const double ev_v,
			       const double ev_u,const double ev_e,
			       const double er);
    Individual*  member(unsigned k);
    Vector<double>*      mblup1(const double ev_v,const double ev_u,const double ev_e,
				const double er);
    Vector<double>*      mblup(const double ev_v,const double ev_u,const double ev_e,
			       const double er);

    void         remove_mark(void);
    unsigned     mark_descendant_of(const unsigned id,const unsigned
				    counter=1);

    unsigned     mark_ancestor_of(const unsigned id,const unsigned counter=1);
    unsigned     mark_relative_of(const unsigned id,const unsigned counter=1);
    unsigned     mark_families(const unsigned start_family_id=1);
    Population*  sub(const unsigned subsize);
    Population*  ancestor_of(const unsigned id);
    Population*  descendant_of(const unsigned id);
    Population*  relative_of(const unsigned id);
    Population*  family_of(const unsigned id){return relative_of(id);}
    Population*  ancestor_of(const char name[]);
    Population*  descendant_of(const char name[]);
    Population*  relative_of(const char name[]);
    Population*  family_of(const char name[]);
    unsigned     fetch_families(const char filename[]);
    unsigned     fetch_families(std::ostream& stream);
    unsigned     n_nufamily(void) const { return num_nufamily; }
    unsigned     ind_gamete(const Individual* I, const unsigned jc,
			    unsigned g_id[], double fq[]);
    double       cprob_children(const Individual* I,const unsigned jc,
				unsigned sg[],unsigned dg[],double sgf[],double dgf[]);
    unsigned     full_cdist(Individual *I,const unsigned jc,
			    unsigned** cdist_value, double* cdist_prob,
			    unsigned** gid_mat, double** freq_mat);
    unsigned     partial_cdist(Individual *I,const unsigned jc,
			       unsigned** cdist_value, double* cdist_prob,
			       unsigned** gid_mat, double** freq_mat);
    int          peeling_sequence(void);
    int          detect_loop(void);
    void         maxant_maxpost(void);
    void         break_loop(void);
    void         cond_genotype_config();
    void         genotype_config(const char type[]);
    void         build_pop_gamete();
    void         gibbs_iterate(const int count_genotype=0);
    void         count_genotype(Individual* I);
    void         release_genotype_counter(void);
    void         resize_genotype_counter(void);
    void         build_nufamily(void);
    void         anterior(Individual* I,doubleMatrix& wspace);
    void         posterior(Individual* I,Individual* J,doubleMatrix& wspace,
			   const unsigned pj);
    void         anterior_posterior(Individual* I, doubleMatrix& wspace);
    void         genotype_dist_peeling(const int prtlevel=1,
				       const int estifreq=1);
    double       log_likelihood_peeling(const unsigned maxit=3);
    double       get_posterior(Individual* I,Individual* exclJ,Vector<double> &vec);
    double       fullsibs_prob(Individual* dam,Individual* sire,
			       Individual* excludeI,doubleMatrix& wspace);
    Vector<double>       get_genotype_freq(void);

    void        maxant_maxpost_old(void); // Old version
    static void gtindex(const unsigned num,const unsigned ni,const Vector<unsigned> &ngt,Vector<unsigned> &gtvec);

    NuFamily**      nufamily(void)   { return nufamily_vec; }
    const DataNode* record_ind(const unsigned i) const {return popmember[i]->record();}
    // BRS
    void        set_switches(void);
    double      multipoint_likelihood(int maxit);
    double      multipoint_init_parm(Fpmm& Farg);
    Matrix<int>     switch_table, switch_table_gmt;
    doubleMatrix      switch_table_prob;
    Vector<double>      Q_freq;
    Vector<Vector<double> > marker_freq;
    Vector<Vector<double> > locusFreq; // LRT (pop_graph modification) 
    void multi_geno_dist_peeling(const unsigned prtlevel);
    std::string    analysis_type;
    Vector<double>    mean_for_pgenotype;
    Vector<double>    P_freq;
    int P_ndim;
    Fpmm *F;

    // from Multi mixture stuffd
    double multi_m_log_likelihood_peeling(int maxit);
    void multi_m_geno_dist_peeling(const unsigned prtlevel);
    // BRS

  //MJS descent graph stuff 
  // moved to matvec-1.0.2d by RLF
  //protected:
  int founder_allele_counter;
  SafeMVVector <SafeMVVector<unsigned> > connected_groups; //CS changed it to vector of vectors to store connected groups for all markerLoci in a sample
  SafeMVVector<SafeMVVector <int> > allele_vector1, allele_vector2; //CS changed it to vector of vectors to store allele vectors for all markerLoci in a sample
  
  SafeMVVector <unsigned> connect_counter; //CS made it a vector to store the counter for all loci in a sample
  SafeMVVector <double> RecoVector;
  void   build_connected_groups(unsigned);
  int    update_allele_vectors(unsigned, Individual *ind);
  void build_founder_allele_neighbors(unsigned lcs);
  int build_allele_vector(unsigned lcs, unsigned connected_group);
  void process_alleles_neighbors(unsigned lcs, unsigned founder_allele);
  void   show_descent_graph_stuff(unsigned lcs);
  double calc_prior_descent_graph(unsigned lcs, bool &initSamplerToOrderFounders);
  void   descent_graph_init_parm(void);
  double descent_graph_log_lhood(bool &initSamplerToOrderFounders); //CS July 04
  void   descent_graph_setup(Data*);
  double MapR(double dist);
  double M_H_sample_ext(unsigned QTL, double l_hood, unsigned &option, bool &initSamplerToOrderFounders); //CS July 04
  double MH_ibd_sample(unsigned l, double log_lhood, unsigned &option, bool &initSamplerToOrderFounders);
  unsigned SLlocus;
  double calc_log_q_ratio(void);
  void   show_change(const int, const int);
  Model  *model;
  void   initPDQs(void);
  void   sum_descentState(unsigned);
  void   sum_genotype_freq1(unsigned);
  void   sum_genotype_freq2(unsigned);
  void   sum_genotype_freq3(unsigned);
  int    valid_graph; // used to indicate result of "process_alleles_neighbors"
  SafeSTLVector< SafeSTLVector< SafeSTLVector<int> > > founder_allele_neighbors;  // CS changed this vector of vectors to a vector of vectors of vectors 
                                                                                  // in order to store the founder allele neighbors for all loci in a sample,
                                                                                  // 1st dimension = locus, 2nd dim = founder allele, 3rd dim = animals
  SafeSTLVector<SafeSTLVector<SafeSTLVector<int> > >founder_allele_neighbors_all; // see above, this stores the founder allele neighbors for typed and untyped 
												                          // ind, after an allele vector has been sampled

  void output_pdq(int n_samples, char *fname);
  void   Gibbs_sample(void);
  //MJS
  void update_ibdMatrix(unsigned locus);
  void output_ibdMatrix(unsigned n_samples,char *fname);

  // CS begin
  void   output_pdq(unsigned locus, int n_samples,char *fname);
  bool orderHeterozygotes(unsigned lcs);
  void displayFounderHaplotypeOriginProbs(unsigned nsamples);
  void accumulateFounderHaplotypeOrigin(unsigned lcs);
  void sampleAlleleVectors(unsigned lcs);
  void initFounderHaplotypeOriginProbs(unsigned QTL);
  void outputFounderHaplotypeOriginProbs(unsigned nsamples, char* fname);
  void build_founder_allele_neighbors_all(unsigned lcs);
  void countHeterozygotes(unsigned lcs);
  void pickFounderLocusToOrder(unsigned orderSamples, unsigned PQTL);
  double findOptimumFounderLocusToOrder(int orderSamples, int QTL, bool & initSamplerToOrderFounders);
  bool orderFoundersGenotypes;
  SafeMVVector<SafeMVVector <int> > previousAlleleVector1, previousAlleleVector2;
  SafeMVVector <unsigned> previousConnectCounter; //CS made it a vector to store the counter for all loci in a sample
  SafeMVVector<SafeMVVector <int> > previousConnectedGroups;
  SafeSTLVector<SafeSTLVector<SafeSTLVector<int> > > previousFounderAlleleNeighbors;
  void accumulateHaplotypes(void);
  void displayHaplotypes(unsigned nsamples);
  template <class T> std::string toString(const T t, std::ios_base & (*f)(std::ios_base&)) const{
	  std::ostringstream oss;
	  oss<< f << t;
	  return oss.str();
  }
  double pop_graph_scaling;	  
  // CS end
  
  //FVP
  Matrix<unsigned> re_pdq;
  Vector<unsigned> temp;
  double log_lhood;
  unsigned initial; 
  ///////6/18/2003 RLF/LRT /////////
    unsigned numFounders;
    void countFounders(void);
    Matrix <double> recombinationMatrix;
    static double disPenetranceTable[2][3];
//    SafeSTLVector<GNodeSet*> vectorOfGNsts;
    GNodeList gNodeList;
    void initJointAlleleNodeList(unsigned howManyLoci);
    void initAlleleNodeList(unsigned whichLocus);
    void initGenotypeNodeList(unsigned whichLocus);
    void getInitialGNodeListSample(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock,string samplerType);
	void getInitialSampleForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock);
    void getInitialGNodeListSample(unsigned maxsize, unsigned numLoci, string samplerType);
	void getInitialGNodeListSampleMIM(void);
	void getInitialSampleForLocusMIM(unsigned j);
	void getGNodeListSampleGibbsMIM(void);
	void getGNodeListSampleSeqInMIM(void);
	void getSampleForLocusMIM(unsigned j);
    void getOldGNodeListProbability(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock,string samplerType);
    void getOldProbabilityForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock);
    void getGNodeListSample(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock,string samplerType);
    void getSampleForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock);
    void getGNodeListSample(unsigned maxsize, unsigned numLoci, string samplerType);
    void setSegregationIndex(unsigned atLocus,string samplerType);
    void setupGNodeSampler(RPedigree& P,MIM& M,GeneticDist& G);
    void setupRSampler(Pedigree& P,Data *D,GeneticDist *G);
    void GNodeStructureSetup(void);
    void SimpleGenotypeElimination(unsigned locus);
    void LGGenotypeElimination(unsigned locus);
    void setAlleleStateVectors(unsigned locus);
    void buildRecombinationMatrix(void);
    void finishUpAlleleStateInit(unsigned locus);
    void switchStateBack(string samplerType);
    void copyGNodeStatesToCandidateStates(string samplerType);
    void copyCandidateToAccepted(string samplerType);
	void copyCandidateGamete(void);
    void storeSampledGametes(void);
    void retreiveSampledGametes(void);
    bool lookToYourLeft;
    bool lookToYourRight;
    IndexVector haplotypeCoder;
    unsigned numHaplotypes;
    void initGenotypeFreq(void);
	void initCounters(void);
	void updateCounters(void);
	void updateCountersSimple(void);
	void outputCounters(unsigned numOfSamples);
    void countHaplotypes(string samplerType);
	void countHaplotypesSimple(string samplerType);
    void countGenotypes (string samplerType); 
	void countGenotypesSimple (string samplerType); 
    void displayHaplotypeFrequencies(unsigned numSamples, std::ostream &os = cout);
	void displayGenotypeFrequencies (unsigned numSamples, std::ostream &outfile = cout);
	void outputQ2Probs              (unsigned numSamples, std::ostream &outfile = cout);
	void displayPDQ(char *fname);
	void displaySegregationIndicators(std::ostream &outfile = cout);
	unsigned calcTotalNumberOfGenotypes(unsigned locus);
	void identifyLociToSamplePosFor(void);
	void samplePositionOnChromosome(void);
	SafeSTLVector<double> sampledPosition;
	double logOldProposalForPosition;
	double logNewProposalForPosition;
	unsigned long samplePosForMe;
    void sampleSegregationIndicators(void);
	void sampleSegregationIndicatorsSimple(void);
    void PrintSegregationIndicators(void);
    void resizeSegregationIndicators(void);
	double MH_ibd_sample_map(unsigned l, double log_lhood, unsigned &option, bool &initSamplerToOrderFounders); //adapted by CS Dec, 2004
    void sum_descentState_map(void);
    bool markerDataIn;
    void setHaploSegregationIndex(unsigned atLocus, unsigned atMarker);
    void ListAlleleFounders(void);
    void SetPossibleHaplotypes(void);
    void SetFreqHaploFounders(void);
    void UpdateFreqHaploFounders(void);
	void UpdateFreqHaploFoundersSimple(void);
    void CalcFreqHaploFounders(unsigned numOfSamples);
    void DisplayFreqHaploFounders(void);
	void getOldGNodeListProbabilityMIM(void);
	void getOldProbabilityForLocusMIM(unsigned j);
	matvec::Matrix<double> getIBDMatrix(void);
	matvec::Matrix<double> getIBDMatrixSimple(void);
	matvec::Matrix<double> getMeans(void);
	matvec::Matrix<double> getMeansSimple(void);
	void setMissingToHet(unsigned locus);
	GNodeSet getSirePivots(unsigned locus);
	void countParentOffspring(string samplerType);
	void countParentOffspringSimple(string samplerType);
    void initParentOffspring(void);
	void displayParentOffspringProbs(std::ostream &outfile);
	void copyGNodeStatesToAcceptedStates(string typeOfGNode);
	void initJointNodeList(unsigned howManyLoci);
	void displayGenotypeProbs(std::ostream &outfile);
	void descent_graph_setup(MIM &mim);
	void displayGNodeProbs(std::ostream &outfile);
	void calcAndDisplayOriginProbsForLocus(unsigned j, std::ostream &outfile);
	void calcAndDisplayAlleleProbsForLocus(unsigned j, std::ostream &outfile);
	void calcAndDisplayGenotypeProbsForLocus(unsigned j, std::ostream &outfile);	
    void calcDistanceToIndividual(std::string strId);
  }; // end of Population class
  
  inline const char* Population::ind_name(const unsigned id) const
  {
	  if (myRPedPtr==0){
		  if (hashtable.size()==0) {
			  string stringName = toString<int>(id, std::dec);
			  return (const char*) stringName.c_str();
		  }
		  else return (const char *)hashtable.find(id);
	  }
	  else {
		return myRPedPtr->pedVector[id-1]->ind_str.c_str();
	  }
  }
} ////// end of matvec namespace

#endif
