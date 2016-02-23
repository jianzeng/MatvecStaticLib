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

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <set>
#include <list>

#include "doublematrix.h"
#include "datanode.h"
#include "genome.h"
#include "geneticdist.h"
#include "gnodederived.h"
#include "dblock.h"
#include "safe_vectors.h"
#include "bg.h"
#include "matrix.h"

namespace matvec {
class Individual;
typedef double (*penetranceType)(const Individual*, const double **);
extern int compare_ind_pt(const void* x, const void* y);
extern int compare_ind_id(const void* a, const void* b);
extern int compare_ind_gid(const void* a, const void* b);
extern int compare_mother_id(const void* x, const void* y);
extern int compare_father_id(const void* x, const void* y);

class Population;
class Chromosome;
class NuFamily;

/*!
 * \struct RNormal Individual.h
 */
struct RNormal
{
  double k, a, b, c;
  unsigned done;
};

/*!
 * struct UNormal Individual.h
 */
struct UNormal
{
  double nu,tsq,k;
};

/*!
 * \class Anterior_c
 * \brief something here
 */
class Anterior_c{
public:
  Vector<UNormal> for_gen;

  Anterior_c(void){;};
  Anterior_c(unsigned n){for_gen.resize(n);};
};

/*!
 *  \class Posterior_c
 *  \brief something here
 */
class Posterior_c{
public:
  Vector<UNormal> for_gen;
  unsigned done;
  Posterior_c(void){;};
  Posterior_c(unsigned n){for_gen.resize(n); done=0;};
};

/*!
 * \class Mixed_post_vec Individual.h
 */
class Mixed_post_vec{
 public:
  Vector<Vector<UNormal> > postvec;
  int done;
};

extern double global_mu;           //TEMP
extern double global_vg;           //TEMP
extern double global_ve;           //TEMP

//RLF
/*!
 *  \class Individual Individual.h
 *  \brief an individual consists of genomes
 *
 *   \sa Genome
 */
class Individual {
   friend class Population;
   friend class NuFamily;
   protected:
      unsigned      myid, numchrom, numoffs, numspouse, numtrait;
      unsigned      *numoffs_spouse;
      Individual    **myoffspring;
      Individual    **spouselist;
      penetranceType  penetrance_f;

      char          mysex; // changed from an int
      int           connect_keeper;

//BRS
  int index_sw, bet_sw, eps_sw; // individual switch index, beta (maternal gamete)  and epsilon (paternal)
    int family;
      Vector<double>        gprobs;
//BRS
      double        inbc;
      Genome        genome0, genome1;
      Vector<double>*       genotype_counter;

      std::set<Individual *>  offs_tree;
      std::list<NuFamily *>   related_family;

      void copyfrom(const Individual& A);

   public:
	  Individual    *myfather, *mymother;
      int           loop_connector,p_origin;
      unsigned      group_id;  // an auxilary group id as working variable
      int           genotype_id;
      double        RBV,MBV, est_GV, true_GV, xbzu_val;
      static Population   *population; 
      DataNode     *myrecord;
      Vector<double>       anterior, *posterior,posterior_iw;
      double             anterior_iw;
      doubleMatrix       *residual_var;

      Individual(void);
      Individual(Population *P);
      Individual(const Individual& A);
      ~Individual(void){release();}

      const Individual&  operator=(const Individual& A);

      int     terminal() const {return !numoffs;}
      int     base() const {return !(myfather || mymother);}
      int     isolated() const {return !(myfather || mymother || numoffs);}

      void    reset_id(const unsigned newid) {myid = newid;}
      void    release(void);
      void    remodel(Population *P);
      void    gamete(Chromosome *C, const int n,const double r) const;
      void    display() {save(std::cout);}
      void    save(std::ostream& out);
      void    xbzu(const double x) {xbzu_val = x;}
      void    initial_anterior(const double *gfreq,const unsigned tng,
                               const int cond=-1);
      void    initial_posterior(const unsigned tng,const int cond=-1);
      void    pretend_missing(int on,const double *gfreq,const unsigned tng);
      void    setPenetrance(penetranceType p){penetrance_f = p;}

      double      get_penetrance(Vector<double> &pen);
      double      xbzu(void) const {return xbzu_val;}
      double      genotypic_val(void) const;
      double      inbcoef(void) const {return inbc;}
      double      father_inbcoef(void) const;
      double      mother_inbcoef(void) const;

      unsigned    nchrom(void) const {return numchrom;}
      unsigned    family_id(void) {return group_id;}
      unsigned    marked(void) const {return group_id;}

      unsigned    id(void) const {return myid;}
      unsigned    nspouse(void) const {return numspouse;}
      unsigned    noffs(void) const {return numoffs;}
      const unsigned* noffs_spouse(void) const {return numoffs_spouse;}
      char  sex(void) const {return mysex;}

      unsigned    gid(void) const {return group_id;}
      unsigned    paternal_chrom_id(const unsigned c) const
                       {return genome0.chromosome[c].id();}
      unsigned    maternal_chrom_id(const unsigned c) const
                       {return genome1.chromosome[c].id();}

      unsigned    father_id(void) const;
      unsigned    mother_id(void) const;
      unsigned    father_gid(void) const;
      unsigned    mother_gid(void) const;

      DataNode*    record(void) const {return myrecord;}
      Chromosome*  paternal_chrom(void) {return genome0.chromosome;}
      Chromosome*  maternal_chrom(void) {return genome1.chromosome;}
      Individual*  father(void) const {return myfather;}
      Individual*  mother(void) const {return mymother;}
      Individual** offspring(void) const {return myoffspring;}
      Individual** spouses(void) const {return spouselist;}

      void   genotype(const unsigned c,const unsigned l,unsigned gtype[]) const;
      void   set_genotype(const unsigned c,const unsigned l,const unsigned a0,
                          const unsigned a1);
//BRS
      void   set_switch(int Nloci);
      void   initial_multi_anterior(doubleMatrix& penetrance);
      void   initial_multi_posterior(const int cond);
      void   collapse_antpost(void);
      int    m_anterior_iw;
      Dblock m_anterior;
      Vector<Dblock> m_posterior;
      double m_anterior_scale;
      Vector<double> m_posterior_scale;
      unsigned n_switches(void);
      double get_m_penetrance(doubleMatrix& pen);
      void pretend_multi_missing(int on, doubleMatrix& pen);
      void put_family(int i){family = i;}
      int  get_family(void ){return family;}
      void cal_gprobs(Dblock& post_mat_f);
      Vector<double> get_gprobs(void){return gprobs;}
      double get_m_posterior(int excJ ,Dblock& post_mat);
      Vector<Vector<UNormal> > mix_anterior;
      Vector<Mixed_post_vec > mix_posterior;
      void initial_multi_m_posterior(const unsigned tn_qtl,const int cond );
      void initial_multi_m_anterior(const unsigned tn_qtl);
      void pretend_multi_m_missing(int on, int tng);
      void cal_m_gprobs(int tng);
      void collapse_mix_antpost(void);
      double get_mix_penetrance(Vector<double> &pen);
      static doubleMatrix g_weight;
//BRS

//MJS
      SafeMVVector <unsigned> m_gamete, m_gamete1, p_gamete, p_gamete1, 
	                          m_gameteAccepted,    p_gameteAccepted;
      SafeMVVector <unsigned> m_founder, p_founder;
      static std::vector <std::vector<double> > vec_prob; //SDK was vector
      static std::vector <std::vector<double> > vec_cutsetval; //SDK was vector
      unsigned p_counter;
      unsigned m_counter;

      void put_gametes(Vector <unsigned> m_g,
		       Vector <unsigned> p_g);
      void  set_founder_alleles(bool);
      Vector <unsigned> &get_founder(std::string, bool);
      void update_offsprings_founder_alleles(void);
      void t0(unsigned SLlocus, unsigned gender_parent);
      void t1(unsigned SLlocus);
      void t2a(unsigned SLlocus);
      void t2b(unsigned SLlocus);
      void apply_SL_transition(const unsigned rule, const unsigned lcs);
      void apply_SL_cascade(const unsigned rule, const unsigned SLlocus);
      double graph_trans_prob(unsigned, unsigned, unsigned);
      int marker_index;
      double calc_q(void);
      void sample_haplotypes(const unsigned);
      double q_prob;
      doubleMatrix genotype_freq1;
      doubleMatrix genotype_freq2;
      int ***mode_matrix;      //this is bad stile
      int assigned_founder_allele;
      void sample_self(int parent_indicator);
//MJS

//CS
	  unsigned ord_heter;
	  std::vector<std::vector<int> > orderHetGenotype1; // vector of vectors, first dimension is to store the locus, second dimension to store the allele
	  std::vector<std::vector<int> > orderHetGenotype2;
	  Genome sampledMaternalGenome;
	  Genome sampledPaternalGenome;
	  Genome previousSampledMaternalGenome;
	  Genome previousSampledPaternalGenome;
	  doubleMatrix paternalFounderHaplotypeProbs, maternalFounderHaplotypeProbs;
	  Vector<unsigned> hetCounter;
	  std::map<std::string, int> matHaplotypeMap;
	  std::map<std::string, int> patHaplotypeMap;
	  void pdq_grid(int samples, std::ofstream &outfile);
	  double DMPDQ, SMPDQ;

//CS
	  
//FVP
      
      Vector <unsigned> id_pdq, id_pdq1; //it stores the pdq index (1,2,3 or 4) for maternal and paternal allele
      Vector<unsigned> get_id_pdq(void);
//      unsigned ord_heter; replaced by declaration above under CS
      unsigned initial_order1, initial_order2;
      void search_heteroz(unsigned qtl);
      void get_allele_v1(unsigned lcs, Vector<int> allele_vector1);
      // HG
      void count_haplotype(void);
      void display_freq_haplotype(unsigned);
      Matrix<unsigned> m_haplotype, p_haplotype;
      void map_pdq(int samples, int interval, double r_aq, double r_bq);
      void map_pdq(int samples, int interval, BG r_aq, BG r_bq);
      double mm_pdq, mp_pdq, pm_pdq, pp_pdq;
      BG bgDMPDQ, bgSMPDQ;
      Vector<unsigned> m_counter_map, p_counter_map;
      void parents(void);
      //HG

      //FVP
      
      ///////6/16/2003 RLF/LRT /////////
	std::map<unsigned, unsigned> matHaplotypeCount;
	std::map<unsigned, unsigned> patHaplotypeCount;
	Vector <unsigned> m_gameteOld, p_gameteOld;
	SafeSTLVector<GenotypeNode>     genotNodeVector;
	SafeSTLVector<AlleleStateNode>  malleleStateNodeVector;
	SafeSTLVector<AlleleStateNode>  palleleStateNodeVector;
	SafeSTLVector<AlleleOriginNode> malleleOriginNodeVector;
	SafeSTLVector<AlleleOriginNode> palleleOriginNodeVector;	
	static SafeSTLVector<unsigned> QTLPosVector;
	static unsigned numLoci;
	static unsigned currentLocus;
	static unsigned currentPosition;
	static unsigned nQTL;
	static double isqrt2pi;
	unsigned allAllelesForLocus(unsigned j);
	unsigned allGenosForLocus(unsigned j);
	double getAllelePenetrance(void);
	double getDisAllelePenetrance(void);
	double getAllelePenetrance(unsigned forLocus);
	double getDisAllelePenetrance(unsigned forLocus);
	double getGenoPenetrance(void);
	double getDisGenoPenetrance(void);
	double getAlleleMTransmissionProb(void);
	double getAllelePTransmissionProb(void);
	double getAlleleMRTransmissionProb(unsigned forLocus);
	double getAllelePRTransmissionProb(unsigned forLocus);	
	double getTransitionProb(void);
	double getMTransmissionProb(void);
	double getPTransmissionProb(void);
	double getTransmissionProb(unsigned iA, unsigned mA, unsigned pA,
                                   Vector <unsigned> &gamete);
	double getMatthiasTransitionProb(void);
	double getMatthiasMTransmissionProb(void);
	double getMatthiasPTransmissionProb(void);
	double getMatthiasTransmissionProb(unsigned iA, unsigned mA, unsigned pA,
					   Vector <unsigned> &gamete,
					   Vector <unsigned> &gameteOld);
	int findLeftLocus(Vector <unsigned> &gamete);
	int findRightLocus(Vector <unsigned> &gamete);
	void setFounderHaplotype(void);
	friend ostream& operator<<(ostream& stream, Individual& i);
	void displayGenotypes(unsigned locus);
	void displayAlleleVectors(unsigned locus);
	void setAlleleStateVectors(unsigned locus);
	void setAlleleOriginVectors(unsigned locus);
	//HG
	Vector<Vector<double> > mHaplotypeFreq , pHaplotypeFreq ;
	//HG
	void setOwnerGNodes(void);
	void setSegregationIndex(unsigned atLocus, string samplerType);
	void sampleSegregationIndicators(void);
	Matrix<unsigned> motherOffspring, fatherOffspring;
	double getGenotypeMRTransmissionProb(unsigned forLocus);
	double getGenotypePRTransmissionProb(unsigned forLocus);
	double getDisGenoPenetrance(unsigned forLocus);
	int distanceToPivot0;
	SafeSTLVector<Individual*> pedBlock, neighborBlock;
	void calcDistanceToNeighbors(SafeSTLVector<unsigned>& individualsToBeProcessed);
	void getPedBlockUpto(unsigned endDist);
	bool pedBlockAdjacent(SafeSTLVector<Individual*> jPedBlock);
	bool isNeighborOf(Individual* indj);
}; // end of Individual class

//RLF
extern RNormal* new_n_posterior_vec(const int m);
//RLF

} // end of namespace matvec
#endif
