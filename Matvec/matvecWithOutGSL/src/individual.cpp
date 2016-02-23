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

#include <iomanip>
#include "population.h"
#include "chromosome.h"
#include "individual.h"
#include "model.h"
namespace matvec {
	//BRS: Declaration of static variables
	Population *Individual::population;
	doubleMatrix Individual::g_weight;
	std::vector <std::vector<double> > Individual::vec_prob;
	std::vector <std::vector<double> > Individual::vec_cutsetval;
	SafeSTLVector<unsigned> Individual::QTLPosVector;
	unsigned Individual::numLoci;
	unsigned Individual::currentLocus;
	unsigned Individual::currentPosition;
	unsigned Individual::nQTL;
	double Individual::isqrt2pi = 1.0/std::sqrt(6.283185308);
	extern double ranf(void);
	
	int compare_ind_pt(const void* x, const void* y)
	{
		Individual **U1,**U2;
		U1 = (Individual **)x;
		U2 = (Individual **)y;
		int retval = 1;
		if (*U1 < *U2)
			retval = -1;
		else if (*U1 == *U2) {
			retval = 0;
		}
		return retval;
	}
	
	int compare_ind_id(const void* a, const void* b)
		//  sort the popmember table which contains just address of each individal
		// in popmember table, it isn't time consuming
	{
		Individual **x,**y;
		x = (Individual **)a; y = (Individual **)b;
		return ((*x)->id() - (*y)->id());
	}
	
	int compare_ind_gid(const void* a, const void* b)
		//  sort the popmember table which contains just address of each individal
		// in popmember table, it isn't time consuming
	{
		Individual **x,**y;
		x = (Individual **)a; y = (Individual **)b;
		return ((*x)->gid() - (*y)->gid());
	}
	
	int compare_mother_id(const void* x, const void* y)
	{
		Individual **U1,**U2;
		U1 = (Individual **)x;
		U2 = (Individual **)y;
		return ((*U1)->mother_id() - (*U2)->mother_id());
	}
	
	int compare_father_id(const void* x, const void* y)
	{
		Individual **U1,**U2;
		U1 = (Individual **)x;
		U2 = (Individual **)y;
		return ((*U1)->father_id() - (*U2)->father_id());
	}
	
	Individual::Individual(void)
	{
		population       = 0;
		genotype_id      = -1;
		myid             = 0;
		mysex            = '.';  // Changed from zero
		numtrait         = 0;
		numchrom         = 0;
		myrecord         = 0;
		posterior        = 0;
		p_origin         = 1;
		RBV              = 0.0;
		MBV              = 0.0;
		inbc             = 0.0;
		est_GV           = 0.0;
		true_GV          = 0.0;
		xbzu_val         = 0.0;
		group_id         = 0;
		spouselist       = 0;
		myoffspring      = 0;
		residual_var     = 0;
		loop_connector   = 0;
		connect_keeper   = 0;
		numoffs_spouse   = 0;
		genotype_counter = 0;
		marker_index     = 0;
		assigned_founder_allele = 0;
	}
	
	Individual::Individual(Population *P)
	{
		if (population == P) return;
		population   = P;
		genotype_id  = -1;
		myid         = 0;
		mysex        = '.';
		myrecord     = 0;
		posterior    = 0;
		p_origin     = 1;
		RBV          = 0.0;
		MBV          = 0.0;
		inbc         = 0.0;
		est_GV       = 0.0;
		true_GV      = 0.0;
		xbzu_val     = 0.0;
		group_id     = 0;
		
		spouselist   = 0;
		myoffspring  = 0;
		residual_var = 0;
		loop_connector = 0;
		connect_keeper = 0;
		numoffs_spouse = 0;
		genotype_counter = 0;
		marker_index     = 0;
		assigned_founder_allele = 0;
		remodel(P);
	}
	
	Individual::Individual(const Individual& A)
	{
		numoffs        = 0;
		myrecord       = 0;
		posterior      = 0;
		myoffspring    = 0;
		spouselist     = 0;
		numoffs_spouse = 0;
		copyfrom(A);
		assigned_founder_allele = 0;
	}
	
	void Individual::copyfrom(const Individual& A)
	{
		if (this == &A) return;
		population  = A.population;
		genotype_id = A.genotype_id;
		myid      = A.myid;
		mysex       = A.mysex;                inbc      = A.inbc;
		myfather    = A.myfather;             mymother  = A.mymother;
		MBV         = A.MBV;                  RBV       = A.RBV;
		true_GV     = A.true_GV;              est_GV    = A.est_GV;
		xbzu_val    = A.xbzu_val;             numtrait  = A.numtrait;
		numspouse   = A.numspouse;            numoffs   = A.numoffs;
		group_id    = A.group_id;             offs_tree = A.offs_tree;
		p_origin    = A.p_origin;             numchrom  = A.numchrom;
		genome0     = A.genome0;              genome1   = A.genome1;
		loop_connector = A.loop_connector;  anterior  = A.anterior;
		residual_var = A.residual_var;      related_family = A.related_family;
		offs_tree = A.offs_tree;            connect_keeper = A.connect_keeper;
		marker_index = A.marker_index;
		assigned_founder_allele = A.assigned_founder_allele;
		//BRS
		index_sw=A.index_sw;
		bet_sw=A.bet_sw;
		eps_sw=A.eps_sw;
		family=A.family;
		gprobs=A.gprobs;
		m_anterior=A.m_anterior;
		m_posterior=A.m_posterior;
		m_anterior_scale=A.m_anterior_scale;
		m_posterior_scale=A.m_posterior_scale;
		//BRS
		
		unsigned i;
		if (posterior) { delete [] posterior; posterior=0;}
		if (A.posterior) {
			if (numspouse>0){
				posterior = new Vector<double> [numspouse];
			}
			else {
				posterior = 0;
			}
			for (i=0; i<numspouse; i++) posterior[i] = A.posterior[i];
		}
		if (myrecord) {delete [] myrecord; myrecord = 0;}
		if (A.myrecord) {
			if(numtrait>0){
				myrecord = new DataNode [numtrait];
			}
			else {
				myrecord = 0;
			}
			for (i=0; i<numtrait; i++) myrecord[i] = A.myrecord[i];
		}
		if (spouselist) {delete [] spouselist; spouselist= 0;}
		if (A.spouselist) {
			if (numspouse>0){
				spouselist = new Individual* [numspouse];
			}
			else {
				spouselist = 0;
			}
			for (i=0; i<numspouse; i++) spouselist[i] = A.spouselist[i];
		}
		if (myoffspring) {delete [] myoffspring; myoffspring=0;}
		if (A.myoffspring) {
			if(numoffs>0){
				myoffspring = new Individual* [numoffs];
			}
			else {
				myoffspring = 0;
			}
			for (i=0; i<numoffs; i++) myoffspring[i] = A.myoffspring[i];
		}
		if (numoffs_spouse) {
			delete [] numoffs_spouse; numoffs_spouse = 0;
		}
		if (A.numoffs_spouse) {
			if(numspouse>0){
				numoffs_spouse = new unsigned [numspouse];
			}
			else {
				numoffs_spouse = 0;
			}
			for (i=0; i<numspouse; i++) numoffs_spouse[i] = A.numoffs_spouse[i];
		}
		if (A.genotype_counter) {
			if(numchrom>0){
				genotype_counter = new Vector<double> [numchrom];
			}
			else {
				genotype_counter = 0;
			}
			for (i=0; i<numchrom; i++) genotype_counter[i] = A.genotype_counter[i];
		}
		else {
			genotype_counter = 0;
		}
	}
	
	void Individual::remodel(Population *P)
	{
		population   = P;
		genotype_id  = -1;
		numtrait = population->prior->ntrait();
		numchrom  = population->prior->nchrom();
		
		genome0.remodel(population->prior);
		genome1.remodel(population->prior);
   		sampledMaternalGenome.remodel(population->prior);
  		sampledPaternalGenome.remodel(population->prior);
		previousSampledMaternalGenome.remodel(population->prior);
   		previousSampledPaternalGenome.remodel(population->prior);
		
		if (myrecord) {
			delete [] myrecord;
			myrecord=0;
		}
		if(numtrait>0){
			myrecord = new DataNode [numtrait];
		}
		else {
			myrecord = 0;
		}
		
		myfather = 0;
		mymother = 0;
		if (numoffs_spouse) {
			delete [] numoffs_spouse; numoffs_spouse = 0;
		}
		numspouse = 0;
		if (myoffspring) { delete [] myoffspring; myoffspring = 0;}
		numoffs  = 0;
		RBV          = 0.0;
		MBV          = 0.0;
		inbc         = 0.0;
		mysex        = '.';
		est_GV       = 0.0;
		true_GV      = 0.0;
		group_id     = 0;
		xbzu_val     = 0.0;
		genotype_counter = 0;
		if (posterior) {delete [] posterior; posterior=0;}
		anterior_iw = 0;
		anterior.resize(0);
	}
	
	const Individual& Individual::operator=(const Individual& A)
	{
		copyfrom(A);
		return *this;
	}
	
	unsigned Individual::father_id(void) const
	{
		if (myfather) return myfather->myid;
		else  return 0;
	}
	
	unsigned Individual::mother_id(void) const
	{
		if (mymother) return mymother->myid;
		else return 0;
	}
	
	unsigned Individual::father_gid(void) const
	{
		if (myfather) return myfather->group_id;
		else  return 0;
	}
	
	unsigned Individual::mother_gid(void) const
	{
		if (mymother) return mymother->group_id;
		else  return 0;
	}
	
	double Individual::father_inbcoef(void) const
	{
		if (myfather) return myfather->inbcoef();
		else  return -1.0;
	}
	
	double Individual::mother_inbcoef(void) const
	{
		if (mymother) return mymother->inbcoef();
		else  return -1.0;
	}
	
	void Individual::initial_anterior(const double *gfreq,const unsigned tng,
									  const int cond)
	{
		/////////////////////////////////////////////////////////////////////////
		// Conditional initialization for anterior
		// cond = -1  is the default, means initialize anterior un-conditionally
		// cond = 1, associated with extra data resulting from iterating
		/////////////////////////////////////////////////////////////////////
		
		if (cond == anterior_iw) return;      // anterior_iw is the condition
		unsigned j;
		if (anterior.size() != tng+1) anterior.resize(tng+1);
		if (myrecord[0].missing) {
			for (j=0; j<tng; j++) anterior[j] = gfreq[j];
			anterior[tng] = 0.0;
			return;
		}
		double s;
		Vector<double> pen_vec(tng);
		anterior[tng] = get_penetrance(pen_vec);
		for (s=0.0,j=0; j<tng; j++) s += (anterior[j] = pen_vec[j]*gfreq[j]);
		for (j=0; j<tng; j++) anterior[j] /= s;
		anterior[tng] += std::log(s);
	}
	
	void Individual::initial_posterior(const unsigned tng, const int cond)
	{
		/////////////////////////////////////////////////////////////////////////
		// Conditional initialization for posterior
		// cond = -1  is the default, means initialize anterior un-conditionally
		// cond = 1, associated with extra data resulting from iterating
		/////////////////////////////////////////////////////////////////////
		int i;
		if (posterior) {
			for (i=0; i<numspouse; i++) {
				if (posterior_iw[i] != cond) {
					posterior[i].resize(tng+1,1.0);
					posterior[i][tng] = 0.0;
				}
			}
		}
		else {
			if(numspouse>0){
				posterior = new Vector<double> [numspouse];
				posterior_iw.resize(numspouse,0.0); //BRS
			}
			else {
				posterior = 0;
			}
			for (i=0; i<numspouse; i++) {
				posterior[i].resize(tng+1,1.0);
				posterior[i][tng] = 0.0;
			}
		}
	}
	
	void Individual::pretend_missing(int on,const double *gfreq,const unsigned tng)
	{
		if (on) {
			if (p_origin) {
				record()[0].missing += 1;
				initial_anterior(gfreq,tng);
				if (connect_keeper == 1) {
					initial_posterior(tng,1);
				}
				else {
					initial_posterior(tng);
				}
			}
			else {
				initial_posterior(tng,1);
			}
		}
		else {
			if (p_origin) record()[0].missing -= 1;
		}
	}
	
	//RLF modified to work with pop->mean_for_genotype[i]
	
	double Individual::get_penetrance(Vector<double> &pen)
	{
		////////////////////////////////////////////////////////////////
		// pen should have enough space to hold all penetrance values
		// this works for multi-chromosomes, but not for multi-loci.
		// For two chromosomes with two alleles each, there are 10 genotypes
		//     AB Ab aB ab
		// AB  0
		// Ab  1  2
		// aB  3  4  5
		// ab  6  7  8  9
		////////////////////////////////////////////////////////////////
		unsigned i,j,t;
		unsigned tn_genotype = pen.size();
		if (myrecord[0].missing) {
			for (i=0; i<tn_genotype; i++) pen[i] = 1.0;
			return 0.0;
		}
		
		//*
		
		 unsigned tn_gamete = static_cast<unsigned>(0.5*std::sqrt(static_cast<double>(8*tn_genotype + 1)) - 0.5);
		 for (genotype_id=0,i=0; i<tn_gamete; i++) {
			 for (j=0; j<=i; j++) {
				 pen[genotype_id] = (*penetrance_f)(this,(const double**)residual_var->begin());
				 genotype_id++;
			 }
		 }
		 
		 //RLF
		 double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));   // pi =2.0*std::asin(1.0)
		 double mu, v = *residual_var[0][0];
		 double y = record()[0].double_val();
		 
		 for (i=0; i < tn_genotype; i++) {
			 mu = population->mean_for_genotype[i];
			 pen[i] = inverse_sqrt2pi/std::sqrt(v) * std::exp(-(y-mu)*(y-mu)/(2*v));
			 
		 }
		 //RLF
		 
		 double scale = pen.sum();
		 if (fabs(scale) < 1.0e-8) {
			 warning("this record is not conformable:");
			 this->display();
			 scale = 3.0;
			 for (i=0; i<tn_genotype; i++) pen[i] = 1.0;
		 }
		 else {
			 for (i=0; i<tn_genotype; i++) pen[i] /= scale;
		 }
		 return std::log(scale);
	}
		 //RLF
		 
		 void Individual::gamete(Chromosome *C,const int n,const double r) const
		 {
			 /************************************************
			 * C is non marked genome,  D is marked chromosome
			 ****************************************************/
			 if (numchrom != n) throw exception("Individual::gamete(): bad arg");
			 const Chromosome *A1, *A2;
			 A1 = genome0.chromosome;
			 A2 = genome1.chromosome;
			 for (unsigned i=0; i<numchrom; i++) {
				 /*
				  if (ranf() <= 0.5) {
					  C[i].qtl = A1[i].qtl;
					  if (ranf() <= r) {
						  C[i].marker = A2[i].marker;
					  }
					  else {
						  C[i].marker = A1[i].marker;
					  }
				  }
				  else {
					  C[i].qtl = A2[i].qtl;
					  if (ranf() <= r) {
						  C[i].marker = A1[i].marker;
					  }
					  else {
						  C[i].marker = A2[i].marker;
					  }
				  }
				  */
			 }
		 }
		 
		 void Individual::release(void)
		 {
			 if (myrecord)    {delete [] myrecord; myrecord = 0;}
			 if (spouselist)  {delete [] spouselist; spouselist = 0;}
			 if (posterior)   {delete [] posterior; posterior = 0; }
			 if (myoffspring) {delete [] myoffspring; myoffspring = 0;}
			 if (numoffs_spouse) {
				 delete [] numoffs_spouse; numoffs_spouse = 0;
			 }
			 if (genotype_counter) {
				 delete [] genotype_counter; genotype_counter = 0;
			 }
		 }
		 
		 double Individual::genotypic_val(void) const
		 {
			 /////////////////////////////////////////////////////////////////////
			 //  assuming there are two alleles [0, 1],
			 //  then there are a total of 4 genotypes
			 //    ----------------------
			 //      |    0        1
			 //   -----------------------
			 //    0 |  (0 0)   (0 1)
			 //    1 |  (1 0)   (1 1)
			 //  ------------------------
			 //    allele id must be an integer starting from 0
			 ///////////////////////////////////////////////////////////////////////
			 double retval = 0;
			 if (genotype_id >= 0) {
				 retval = population->mean_for_genotype[genotype_id];
			 }
			 else {
				 warning("Individual::genotypic_val(): no genotype available");
			 }
			 return retval;
			 /*
			  unsigned i,j,nc, nl,sd[2];
			  nc = prior->nchrom();
			  ChromStruct *Chrom = prior->chrom();
			  const double** gv;
			  double retval = 0;
			  for (i=0; i<nc; i++) {
				  nl = Chrom[i].nloci();
				  for (j=0; j<nl; j++) {
					  genotype(i,j,sd);
					  gv = prior->genotypic_val(i,j);
					  retval +=  gv[sd[0]][sd[1]];
				  }
			  }
			  return retval;
			  */
		 }
		 
		 void Individual::set_genotype(const unsigned c,const unsigned l,
									   const unsigned a0,const unsigned a1)
		 {
			 genome0.chromosome[c].locus[l].allele = a0;
			 genome1.chromosome[c].locus[l].allele = a1;
		 }
		 
		 void Individual::genotype(const unsigned c, const unsigned l,unsigned gtype[])
		 const
		 {
			 gtype[0] = genome0.chromosome[c].locus[l].allele;
			 gtype[1] = genome1.chromosome[c].locus[l].allele;
		 }
		 
		 void Individual::save(std::ostream& out)
		 {
			 out << "  I " << myid << "   F " << father_id() << "   M " << mother_id();
			 out << " Marker Info: ";
			 out << genome0.chromosome[0].locus[1].allele << " " ;
			 out << genome1.chromosome[0].locus[1].allele << " " ;
			 out << genome0.chromosome[0].locus[2].allele << " " ;
			 out << genome1.chromosome[0].locus[2].allele;
			 //out << "   Switch " << index_sw << "   Beta " << bet_sw << "   Epl " << eps_sw;
			 out << std::endl;
		 }
		 
		 //BRS
		 void Individual::set_switch(int Nloci)
		 
		 {
			 /* This assigns the switches and gametes to each individual in the population.
			 By definition, the maternal gamete gets the lower allele number unless it is known
			 which source the allele came from.
			 
			 An locus is switchable only in two circumstances:
			 1)  founder animal is heterozgyous at that locus.
2)  both parents are heterozygous for the same alleles at that locus
Under all other cases, the source of the allele can be tracked.
Thus, each locus needs to be checked to determine the source of the allele.

*/
			 int i,j,k,l,allele1,allele2,sire,dam,sireF,sireM,damF,damM,levelS=1,levelEB=1,Ixswitch=0,epl=0,beta=0,temp,tallele;
			 if (myfather && (mymother == 0)) {
				 std::cout << "error father known, mother not! for " << population->ind_name(id()) << std::endl;
				 std::cout << "This program is not designed to handle this yet!\n Please supply a mother with records." << std::endl;
				 exit(1);
			 }
			 if (mymother && (myfather == 0)) {
				 std::cout << "error mother known, father not! for " << population->ind_name(id()) << std::endl;
				 std::cout << "This program is not designed to handle this yet! Please supply a father with records." << std::endl;
				 exit(1);
			 }
			 
			 if (myfather == 0)  { // animal is a founder
				 for (j=(Nloci); j > 0; j--)  {
					 if (genome0.chromosome[0].locus[j].allele > genome1.chromosome[0].locus[j].allele) {  // swap alleles
						 temp=genome0.chromosome[0].locus[j].allele;
						 genome0.chromosome[0].locus[j].allele=genome1.chromosome[0].locus[j].allele;
						 genome1.chromosome[0].locus[j].allele =temp;
					 }
					 if (genome0.chromosome[0].locus[j].allele != genome1.chromosome[0].locus[j].allele) {   // switchable
						 Ixswitch += levelS;
					 }
					 levelS *=2;
				 }
				 index_sw=Ixswitch;
			 }
			 else { // not a founder so find all the alleles from individual and parents
				 for (j=(Nloci); j>0; j--)  {
					 allele1=genome0.chromosome[0].locus[j].allele;
					 allele2=genome1.chromosome[0].locus[j].allele;
					 sireM=myfather->genome0.chromosome[0].locus[j].allele;
					 sireF=myfather->genome1.chromosome[0].locus[j].allele;
					 damM=mymother->genome0.chromosome[0].locus[j].allele;
					 damF=mymother->genome1.chromosome[0].locus[j].allele;
					 if ((allele1 != allele2) && (((sireF == damF) && (sireM == damM)) || ((sireF == damM) && (sireM == damF)))) {
						 Ixswitch += levelS; // Locus is switchable !!!
											 // Check if alleles in correct order if not swap them
						 if (allele1 > allele2) {
							 genome1.chromosome[0].locus[j].allele=allele1;
							 genome0.chromosome[0].locus[j].allele=allele2;
						 }
					 }
					 // Not switchable so determine source of each allele
					 else if ((allele1 == sireM) && (allele2 == damM)) {
						 genome0.chromosome[0].locus[j].allele=allele2;
						 genome1.chromosome[0].locus[j].allele=allele1;
					 }
					 else if ((allele1 == sireM) && (allele2 == damF))  {
						 genome0.chromosome[0].locus[j].allele=allele2;
						 genome1.chromosome[0].locus[j].allele=allele1;
					 }
					 else if ((allele1 == sireF) && (allele2 == damM)) {
						 genome0.chromosome[0].locus[j].allele=allele2;
						 genome1.chromosome[0].locus[j].allele=allele1;
					 }
					 else if ((allele1 == sireF) && (allele2 == damF))  {
						 genome0.chromosome[0].locus[j].allele=allele2;
						 genome1.chromosome[0].locus[j].allele=allele1;
					 }
					 else if ((allele2 == sireM) && (allele1 == damM))  {
						 genome0.chromosome[0].locus[j].allele=allele1;
						 genome1.chromosome[0].locus[j].allele=allele2;
					 }
					 else if ((allele2 == sireM) && (allele1 == damF))  {
						 genome0.chromosome[0].locus[j].allele=allele1;
						 genome1.chromosome[0].locus[j].allele=allele2;
					 }
					 else if ((allele2 == sireF) && (allele1 == damM))  {
						 genome0.chromosome[0].locus[j].allele=allele1;
						 genome1.chromosome[0].locus[j].allele=allele2;
					 }
					 else if ((allele2 == sireF) && (allele1 == damF))  {
						 genome0.chromosome[0].locus[j].allele=allele1;
						 genome1.chromosome[0].locus[j].allele=allele2;
					 }
					 else    {
						 std::cerr << "ERROR ident alleles different from parents - perhaps parent not in pedigree?\n"
						 << "Id " <<  allele1 << " " << allele2 << " "  << sireM  << " " << sireF  << " " << damM << " "  << damF << std::endl;
						 throw exception("ERROR ident alleles different from parents - perhaps parent not in pedigree?");
					 }
					 // Now determine epsilon
					 if (sireM == sireF)
						 epl += (levelEB*2);
					 else if ( genome1.chromosome[0].locus[j].allele == sireF)
						 epl += levelEB; // Paternal allele from father's father
										 // Now determine  beta
					 if (damM == damF)
						 beta += (levelEB*2);
					 else if (genome0.chromosome[0].locus[j].allele == damF)
						 beta += levelEB; // maternal allele from dam's father
					 
					 levelS *=2;
					 levelEB *=3;
				 }
			 }
			 index_sw=Ixswitch;
			 eps_sw=epl;
			 bet_sw=beta;
			 //std::cout << "Switch " << Ixswitch << " Epli " << epl << " beta " << beta << std::endl;
}

unsigned Individual::n_switches(void){
	return (population->switch_table[index_sw][0]);
}


void Individual::initial_multi_anterior(doubleMatrix& penetrance)
{
	unsigned i,j,k,a1,a2,ndim;
	double freq, scale, sum=0.0, freq_mark=1.0;
	unsigned n_switch = population->switch_table[index_sw][0];
	ndim=population->P_ndim;
	if (m_anterior.get_nrow() != n_switch) m_anterior.resize(ndim,n_switch,4,0.0);
	m_anterior_iw = 0;
	
	for (k=1;k<=population->n_markerLoci; k++){
		a1 = genome0.chromosome[0].locus[k].allele;
		a2 = genome1.chromosome[0].locus[k].allele;
		freq_mark *= population->marker_freq(k)(a1)*population->marker_freq(k)(a2);
	}
	if (myrecord[0].missing) {
		m_anterior_scale =0.0;
		for (k=0;k<ndim; k++){
			for (j=0; j<4; j++) {
				freq = (population->Q_freq[j])*freq_mark*(population->P_freq[k]);
				for (i=0; i<n_switch; i++) {
					m_anterior[k][i][j] = freq;
					sum += freq; // accumulating scaling factor
				}
			}
		}
	}
	else {
		//    std::cout << " before call Pen " <<  std::endl;
		m_anterior_scale = get_m_penetrance(penetrance);
		//   std::cout << "scale in intial_anterior "<< m_anterior_scale << std::endl;
		//    std::cout << " Pen " << penetrance << std::endl;
		for (k=0;k<ndim; k++){
			for (j=0; j<4; j++) {
				freq = (population->Q_freq[j])*freq_mark*(population->P_freq[k])*penetrance[k][j];
				//  std::cout << "k " << k << " j " << j << " Q_f " <<  population->Q_freq[j];
				//   std::cout << " M_f " << freq_mark << " Pen " << penetrance[k][j] << " freq " << freq << std::endl;
				for (i=0; i<n_switch; i++) {
					m_anterior[k][i][j] = freq;
					sum += freq; // accumulating scaling factor
				}
			}
		}
	}
	//  std::cout << "sum " << sum << std::endl;
	
	// Need to rescale due to probability of observing the marker genotypes
	for (k=0;k<ndim; k++){
		for (i=0; i<n_switch; i++) {
			for (j=0; j<4; j++) {
				m_anterior[k][i][j] /= sum;
			}
		}
	}
	m_anterior_scale += std::log(sum);
	//  std::cout << "scale " << m_anterior_scale << std::endl;
}

void Individual::initial_multi_posterior(const int cond) {
	/////////////////////////////////////////////////////////////////////////
	// Conditional initialization for posterior
	// cond = -1  is the default, means initialize anterior un-conditionally
	// cond = 1, associated with extra data resulting from iterating
	/////////////////////////////////////////////////////////////////////
	int i,ndim;
	unsigned n_switch = population->switch_table[index_sw][0];
	ndim=population->P_ndim;
	if (numspouse) {
		if (m_posterior.size() != numspouse ) {
			m_posterior.resize(numspouse);
			m_posterior_scale.resize(numspouse,0.0);
			posterior_iw.resize(numspouse,0.0); //BRS
		}
		if (posterior){
			m_posterior_scale.resize(numspouse,0.0);
			posterior_iw.resize(numspouse,0.0); //BRS
			for (i=0; i<numspouse; i++) {
				if (posterior_iw[i] != cond) {
					m_posterior[i].resize(ndim,n_switch,4,1.0);
				}
			}
		}
		else {
			for (i=0; i<numspouse; i++) {
				m_posterior[i].resize(ndim,n_switch,4,1.0);
			}
			m_posterior_scale.resize(numspouse,0.0);
			posterior_iw.resize(numspouse,0.0); //BRS
		}
	}
}

double Individual::get_m_penetrance(doubleMatrix& pen) {
	
	unsigned QTL,PGN;
	
	if (myrecord[0].missing) {
		for (PGN=0; PGN<population->P_ndim; PGN++) {
			for (QTL=0; QTL<4; QTL++) {
				pen[PGN][QTL] = 1.0;
			}
		}
		return 0.0;
	}
	
	double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));  // pi = 2.0*std::asin(1.0)
	double mu, v = *residual_var[0][0],scale;
	double y = record()[0].double_val();
	
	for (PGN=0; PGN<population->P_ndim; PGN++) {
		for (QTL=0; QTL<4; QTL++) {
			mu = population->mean_for_genotype[QTL]+population->mean_for_pgenotype[PGN];
			pen[PGN][QTL] = inverse_sqrt2pi/std::sqrt(v) * std::exp(-(y-mu)*(y-mu)/(2*v));
			//  std::cout << QTL << " Q_mu " << population->mean_for_genotype[QTL] << " P_mu " << population->mean_for_pgenotype[PGN] <<" mu " << mu << " y " << y << " v " << v << " pen is " << setprecision(12) << pen[PGN][QTL] << std::endl;
		}
	}
	//   std::cout << "pen in get_m_pen " << pen;
	scale = (pen.sum()).sum();
	//  std::cout << "Pen scale " << scale << std::endl;
	if (fabs(scale) < 1.0e-8) {
		std::cout << id() << " has record " << y << std::endl;
		warning("this record is not conformable");
		scale = 3.0;
		for (PGN=0; PGN<population->P_ndim; PGN++) {
			for (QTL=0; QTL<4; QTL++)
				pen[PGN][QTL] = 1.0;
		}
	}
	else {
		for (PGN=0; PGN<population->P_ndim; PGN++) {
			for (QTL=0; QTL<4; QTL++) {
				pen[PGN][QTL] /= scale;
			}
		}
		//    std::cout << "pen in get pen " << pen << std::endl;
	}
	return std::log(scale);
}

void Individual::pretend_multi_missing(int on, doubleMatrix& pen)
{
	if (on) {
		if (p_origin) {
			record()[0].missing += 1;
			initial_multi_anterior(pen);
			//     initial_anterior(gfreq,tng);
			if (connect_keeper == 1) {
				initial_multi_posterior(1);
				// initial_posterior(tng,1);
			}
			else {
				initial_multi_posterior(-1);
				// initial_posterior(tng);
			}
		}
		else {
			initial_multi_posterior(1);
			//   initial_posterior(tng,1);
		}
	}
	else {
		if (p_origin) record()[0].missing -= 1;
	}
}
//RLF had  modified to work with pop->mean_for_genotype[i]

double Individual::get_m_posterior(int excJ , Dblock& post_mat) {
	
	/////////////////////////////////////////////////////////
	//  if excJ < 0, say -1, means taking all posteriors
	/////////////////////////////////////////////////////////
	
	// make sure post_mat is initialized properly !!!!!!!!!!!!!!!!
	
	unsigned nswitches = n_switches();
	int spouse,i,j,k,ndim;
	double retval=0.0;
	ndim=post_mat.get_ndim();
	
	for (spouse=0; spouse<numspouse; spouse++) {
		if (spouse != excJ) {
			for (k=0; k < ndim; k++){
				for (i=0;i<nswitches; i++){
					for (j=0;j<4;j++){
						post_mat[k][i][j] *= (m_posterior[spouse])[k][i][j];
					}
				}
			}
			retval += m_posterior_scale[spouse];
		}
	}
	return retval;
}


void Individual::cal_gprobs(Dblock& post_mat_f) {
	double val;
	int ii,jj,kk;
	
	get_m_posterior(-1,post_mat_f);
	gprobs.resize(4);
	val = 0;
	for (jj=0;jj<4; jj++) {
		for (kk=0; kk < post_mat_f.get_ndim() ; kk++) {
			for (ii=0;ii< n_switches(); ii++) {
				gprobs[jj] += post_mat_f[kk][ii][jj]*m_anterior[kk][ii][jj];
			}
		}
		val += gprobs[jj];
	}
	for (jj=0; jj<4; jj++) {
		gprobs[jj] /= val;
	}
}

void Individual::collapse_antpost() {
	int i,j, k, nswitch,ndim;
	double sum=0.0;
	anterior.resize(4,0.0);
	initial_posterior(3,-1);
	nswitch=n_switches();
	ndim=m_anterior.get_ndim();
	// std::cout << nswitch << m_anterior << anterior << std::endl;
	for (k=0; k<ndim; k++) {
		for (i=0;i<nswitch;i++){
			//  std::cout << i << " " << nswitch << " " << (m_anterior[i][0]) << " " << anterior[i] << std::endl;
			anterior[0] += m_anterior[k][i][0];
			anterior[1] += m_anterior[k][i][1] + m_anterior[k][i][2];
			anterior[2] += m_anterior[k][i][3];
			for (j=0;j<numspouse;j++){
				(posterior[j])[0] += (m_posterior[j])[k][i][0];
				(posterior[j])[1] += (m_posterior[j])[k][i][1] + (m_posterior[j])[k][i][2];
				(posterior[j])[2] += (m_posterior[j])[k][i][3];
			}
		}
	}
}

void Individual::pretend_multi_m_missing(int on, int tng) {
	if (on) {
		if (p_origin) {
			record()[0].missing += 1;
			initial_multi_m_anterior(tng);
			//     initial_anterior(gfreq,tng);
			if (connect_keeper == 1) {
				initial_multi_m_posterior(tng,1);
				// initial_posterior(tng,1);
			}
			else {
				initial_multi_m_posterior(tng,-1);
				// initial_posterior(tng);
			}
		}
		else {
			initial_multi_m_posterior(tng,1);
			//   initial_posterior(tng,1);
		}
	}
	else {
		if (p_origin) record()[0].missing -= 1;
	}
}

void Individual::initial_multi_m_anterior(const unsigned tn_qtl) {
	int i_q, i_sw, i_switches,k, a1, a2;
	double freq_mark=1.0;
	double P_var = population->F->var; // polygenic varaince
	i_switches=n_switches();
	
	// compute frequency of marker information
	for (k=1;k<=population->n_markerLoci; k++){
		a1 = genome0.chromosome[0].locus[k].allele;
		a2 = genome1.chromosome[0].locus[k].allele;
		freq_mark *= population->marker_freq(k)(a1)*population->marker_freq(k)(a2);
	}
	
	if (tn_qtl != mix_anterior.size())
		mix_anterior.resize(tn_qtl);
	for (i_q=0; i_q < tn_qtl; i_q++){
		mix_anterior[i_q].resize(i_switches);
		for (i_sw=0; i_sw < i_switches; i_sw++) {
			mix_anterior[i_q][i_sw].nu  = 0.0;
			mix_anterior[i_q][i_sw].tsq = 1.0/P_var;
			mix_anterior[i_q][i_sw].k   = std::log(1.0/std::sqrt(4.0*std::asin(1.0)*P_var)*population->Q_freq[i_q]*freq_mark);
			//       std::cout << id() << " " << i_q << " " << i_sw << " nu " << mix_anterior[i_q][i_sw].nu << " tsq " << mix_anterior[i_q][i_sw].tsq << " k " << mix_anterior[i_q][i_sw].k << std::endl;
		}
	}
	//  std::cout << "anterior initialized for " << id() << std::endl;
}

void Individual::initial_multi_m_posterior(const unsigned tn_qtl, const int cond) {
	
	int i_q, i_sw,spouse,i_switches, i_s;
	i_switches=n_switches();
	if (numspouse) { 
		if (mix_posterior.size() != numspouse) { 
			mix_posterior.resize(numspouse);
			posterior_iw.resize(numspouse,0.0);
			for (spouse=0; spouse<numspouse; spouse++) {
			    mix_posterior[spouse].done = 0;
			    mix_posterior[spouse].postvec.resize(tn_qtl);
			    for (i_q=0; i_q < tn_qtl; i_q++) {
					mix_posterior[spouse].postvec[i_q].resize(i_switches); 
			    } 
			} 
		} 
		else { 
			if (posterior){ 
				mix_posterior.resize(numspouse); 
				posterior_iw.resize(numspouse,0.0); 
				for (spouse=0; spouse<numspouse; spouse++) { 
					if (posterior_iw[spouse] != cond) { 
						mix_posterior[spouse].done = 0; 
						mix_posterior[spouse].postvec.resize(tn_qtl); 
						for (i_q=0; i_q < tn_qtl; i_q++) { 
							mix_posterior[spouse].postvec[i_q].resize(i_switches); 
						} 
					} 
				} 
			} 
		} 
	} 
	else { 
		mix_posterior.resize(numspouse); 
		for (spouse=0; spouse < numspouse; spouse++){ 
			mix_posterior[spouse].done = 0; 
			mix_posterior[spouse].postvec.resize(tn_qtl); 
			for (i_q=0; i_q < tn_qtl; i_q++) { 
				mix_posterior[spouse].postvec[i_q].resize(i_switches); 
			} 
		} 
	} 
}




void Individual::cal_m_gprobs(int tng) {
	double val,y,v_y,k_y,nu_y, a11, a12, a22, k, v, nu, inverse_sqrt2pi, sum1, sum2;
	int i_q, i_sw, i_sp, i_switches,j;
	i_switches=n_switches();
	inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));  // pi = 2.0*std::asin(1.0)
	gprobs.resize(tng);
	val = 0.0;
	sum1=0.0;
	for (i_q=0;i_q<tng;i_q++){
		// contributions from penetrance function
		// skip if phenotype if missing
		
		if (!(record()[0].missing)) {
			y    = record()[0].double_val();
			v_y    = 1.0/(*residual_var[0][0]);
			k_y    = std::log(inverse_sqrt2pi*std::sqrt(v_y));
			nu_y   = population->mean_for_genotype[i_q] - y;
		}
		else{
			nu_y = 0.0;
			v_y = 0.0;
			k_y = 0.0;
		}
		
		for (i_sw=0; i_sw < i_switches; i_sw++){   // loop for genotypes of J
			
			a11 = nu_y*nu_y*v_y;
			a12 = nu_y*v_y;
			a22 = v_y;
			k   = k_y;
			// contributions from anterior
			
			v   = mix_anterior[i_q][i_sw].tsq;
			nu  = mix_anterior[i_q][i_sw].nu;
			k   += mix_anterior[i_q][i_sw].k;
			a11 += nu*nu*v;
			a12 += -nu*v;
			a22 += v;
			
			// Contributions from posteriors.
			// Skip over any posteriors not yet available
			unsigned nspouses =nspouse();
			for (j=0; j<nspouses; j++){
				if (mix_posterior[j].done)  {
					v     = mix_posterior[j].postvec[i_q][i_sw].tsq;
					nu    = mix_posterior[j].postvec[i_q][i_sw].nu;
					k    += mix_posterior[j].postvec[i_q][i_sw].k;
					a11  += nu*nu*v;
					a12  += -nu*v;
					a22  += v;
				}
			}
			
			// compute "log likelihood for gi"
			
			v          = 1.0/a22;
			g_weight[i_q][i_sw]   = -0.5*(a11 - a12*a12*v)
				+ k + std::log(std::sqrt(4.0*std::asin(1.0)*v));  // pi = 2.0*std::asin(1.0)
			sum1 += g_weight[i_q][i_sw];
		}
	}
	
	// compute sum of the "likelihoods"
	
	sum1 /= (tng*i_switches);
	for (i_q=0;i_q<tng;i_q++){
		sum2=0.0;
		for (i_sw=0; i_sw < i_switches; i_sw++){   // loop for genotypes of I
			gprobs[i_q] += (std::exp(g_weight[i_q][i_sw] - sum1));
		}
		sum2 += gprobs[i_q];
	}
	for (i_q=0;i_q<tng;i_q++){
		gprobs[i_q] /= sum2;
	}
}

void Individual::collapse_mix_antpost(){
	int tng=4;
	anterior.resize(tng,0.0);
	initial_posterior(tng,-1);
	double val,y,v_y,k_y,nu_y, a11, a12, a22, k, v, nu, inverse_sqrt2pi, sum1, sum2;
	int i_q, i_sw, i_sp, i_switches,j;
	i_switches=n_switches();
	inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));   // pi = 2.0*std::asin(1.0)
	gprobs.resize(3);
	val = 0.0;
	sum1=0.0;
	for (i_q=0;i_q<tng;i_q++){
		// contributions from penetrance function, skip if phenotype is missing
		if (!(record()[0].missing)) {
			y    = record()[0].double_val();
			v_y    = 1.0/(*residual_var[0][0]);
			k_y    = std::log(inverse_sqrt2pi*std::sqrt(v_y));
			nu_y   = population->mean_for_genotype[i_q] - y;
		}
		else{
			nu_y = 0.0;
			v_y = 0.0;
			k_y = 0.0;
		}
		
		for (i_sw=0; i_sw < i_switches; i_sw++){   // loop for genotypes of J
			
			a11 = nu_y*nu_y*v_y;
			a12 = nu_y*v_y;
			a22 = v_y;
			k   = k_y;
			// contributions from Father's anterior
			
			v   = mix_anterior[i_q][i_sw].tsq;
			nu  = mix_anterior[i_q][i_sw].nu;
			k   += mix_anterior[i_q][i_sw].k;
			a11 += nu*nu*v;
			a12 += -nu*v;
			a22 += v;
			
			// compute "log likelihood for gi"
			v    = 1.0/a22;
			g_weight[i_q][i_sw]   = -0.5*(a11 - a12*a12*v)
				+ k + std::log(std::sqrt(4.0*std::asin(1.0)*v));    // pi = 2.0*std::asin(1.0)
			sum1 += g_weight[i_q][i_sw];
		}
	}
	
	// compute sum of the "likelihoods"
	
	sum1 /= (tng*i_switches);
	for (i_q=0;i_q<tng;i_q++){
		sum2=0.0;
		for (i_sw=0; i_sw < i_switches; i_sw++){   // loop for genotypes of I
			anterior[i_q] += (std::exp(g_weight[i_q][i_sw] - sum1));
		}
		sum2 += anterior[i_q];
	}
    anterior[0] /= sum2;
    anterior[1] /= sum2;
    anterior[2] /= sum2;
    anterior[3] /= sum2;
    anterior[1] += anterior[2];
    anterior[2]= anterior[3];
	
	// Now do the posteriors:
	unsigned nspouses =nspouse();
	for (j=0; j<nspouses; j++){
		sum1=0.0;
		for (i_q=0;i_q<tng;i_q++){
			for (i_sw=0; i_sw < i_switches; i_sw++){   // loop for genotypes of J
													   // Contributions from posteriors
				
				if (mix_posterior[j].done) {
					v     = mix_posterior[j].postvec[i_q][i_sw].tsq;
					nu    = mix_posterior[j].postvec[i_q][i_sw].nu;
					k    += mix_posterior[j].postvec[i_q][i_sw].k;
					a11  += nu*nu*v;
					a12  += -nu*v;
					a22  += v;
				}
				// compute "log likelihood for gi"
				v    = 1.0/a22;
				g_weight[i_q][i_sw] = -0.5*(a11 - a12*a12*v)
					+ k + std::log(std::sqrt(4.0*std::asin(1.0)*v));    // pi = 2.0*std::asin(1.0)
				sum1 += g_weight[i_q][i_sw];
			}
		}
		
		// compute sum of the "likelihoods"
		
		sum1 /= (tng*i_switches);
		for (i_q=0;i_q<tng;i_q++){
			sum2=0.0;
			for (i_sw=0; i_sw < i_switches; i_sw++){   // loop for genotypes of I
				posterior[j][i_q] += (std::exp(g_weight[i_q][i_sw] - sum1));
			}
			sum2 += posterior[j][i_q];
		}
		posterior[j][0] /= sum2;
		posterior[j][1] /= sum2;
		posterior[j][2] /= sum2;
		posterior[j][3] /= sum2;
		posterior[j][1]  +=posterior[j][2];
		posterior[j][2]=posterior[j][3];
	}
}


double Individual::get_mix_penetrance(Vector<double> &pen) {
	
	unsigned QTL,PGN;
	
	if (myrecord[0].missing) {
		for (QTL=0; QTL<3; QTL++) {
			pen[QTL] = 1.0;
		}
		return 0.0;
	}
	
	double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));       // pi = 2.0*std::asin(1.0)
	double mu, v = *residual_var[0][0],scale;
	double y = record()[0].double_val();
	
    mu = population->mean_for_genotype[0];
    pen[0] = inverse_sqrt2pi/std::sqrt(v) * std::exp(-(y-mu)*(y-mu)/(2*v));
    mu = population->mean_for_genotype[1];
    pen[1] = inverse_sqrt2pi/std::sqrt(v) * std::exp(-(y-mu)*(y-mu)/(2*v));
    mu = population->mean_for_genotype[2];
    pen[1] += inverse_sqrt2pi/std::sqrt(v) * std::exp(-(y-mu)*(y-mu)/(2*v));
    mu = population->mean_for_genotype[3];
    pen[2] = inverse_sqrt2pi/std::sqrt(v) * std::exp(-(y-mu)*(y-mu)/(2*v));
	
	scale = pen[0]+ pen[1] + pen[2];
	
	if (fabs(scale) < 1.0e-8) {
		warning("this record is not conformable");
		scale = 3.0;
		for (QTL=0; QTL<3; QTL++)
			pen[QTL] = 1.0;
	}
	else {
		for (QTL=0; QTL<3; QTL++) {
			pen[QTL] /= scale;
		}
		//    std::cout << "pen in get pen " << pen << std::endl;
	}
	return std::log(scale);
}

//BRS

//MJS
// moved to matvec 1.0.2d RLF
//descent graph stuff
void Individual::put_gametes(Vector <unsigned> m_g, Vector <unsigned> p_g) {
	// Authors: Matthias Schelling and Rohan L. Fernando 
	// (1999) 
	// Contributors:
	m_gamete = m_g;
	p_gamete = p_g;
	if(m_gamete.size()!= p_gamete.size()){
		throw exception("Individual::put_gametes: m_gamete.size()!= p_gamete.size()");
	}
	//init q
	q_prob = 1.0;
	calc_q();
}

double Individual::calc_q(void) {
	// Authors: Matthias Schelling and Rohan L. Fernando 
	// (1999) 
	// Contributors: L. Radu Totir 
	if ((m_gamete.size() == 0) || (p_gamete.size() == 0)) {
		throw exception ("Individual::calc_q: no graph present.");
	}
	double old_q_prob = q_prob; //store old q
	if (mymother) { // if not founder
		q_prob = 0.25;
		for (int j=0; j<numLoci-1; j++) {
			if (m_gamete[j] == m_gamete[j+1]) {
				q_prob *= 1 - population->RecoVector[j];
			}
			else {
				q_prob *= population->RecoVector[j];
			}
			if (p_gamete[j] == p_gamete[j+1]) {
				q_prob *= 1 - population->RecoVector[j];
			}
			else {
				q_prob *= population->RecoVector[j];
			}
		}
	}
	else {
		q_prob = 1.0;
	}
	return std::log(old_q_prob) - std::log(q_prob);
}
/*! \fn double Individual::calc_q(void)
*  \brief method to calc q prob for individual; q is stored in 1 variable
*/

void Individual::sample_haplotypes(const unsigned origin) {
	//sample haplotypes for an SLlocus,
	//origin 0 = maternal, 1 = paternal
	Vector <unsigned>* gamete = &m_gamete;
	if (origin) {
		gamete = &p_gamete;
	}
	//first, go from SLlocus down to locus 1
	for (unsigned i = population->SLlocus - 1; i > 0; i--) {
		if (ranf() < population->RecoVector(i)) { //do crossover
			if ((*gamete)(i+1) == 0) { (*gamete)(i) = 1; }
			else { (*gamete)(i) = 0; }
		}
		else { //no crossover
			if ((*gamete)(i+1) == 0) { (*gamete)(i) = 0; }
			else { (*gamete)(i) = 1; }
		}
	}
	//then, go from SLlocus up to last locus
	for (unsigned i = population->SLlocus; i < (*gamete).size(); i++) {
		if (ranf() < population->RecoVector(i)) { //do crossover
			if ((*gamete)(i) == 0) { (*gamete)(i+1) = 1; }
			else { (*gamete)(i+1) = 0; }
		}
		else { //no crossover
			if ((*gamete)(i) == 0) { (*gamete)(i+1) = 0; }
			else { (*gamete)(i+1) = 1; }
		}
	}
}

void Individual::set_founder_alleles(bool update) {
	// Authors: Matthias Schelling and Rohan L. Fernando 
	// (1999) 
	// Contributors: L. Radu Totir 
	
	//  update == 0  ->  first pass
	//  update == 1  ->  already done
	//		     Done
	//  	          T        F
	//  	       ________________
	//  	     T | yes      yes |
	// Update    |              |
	//  	     F | no       yes |
	//           |______________|
	//  			      
	// This means that we set the founder_alleles from
	// the founder_allele_counter only if not already done
	// When we are updating, we change only offsprings
	
	// get space if needed
	if (m_founder.size() == 0) { 
		int zero=0;
		m_founder.resize(numLoci, zero);
		p_founder.resize(numLoci, zero);
	}
	
	// if not already done
	//if (!m_founder(1)) {
	if (!(!update && m_founder(1))) {
		if (!m_founder(1) && update) {
			std::cout << "should not happen: update = 1, m_founder(1) set" << std::endl;
		}
		// if founder set m/p_founder to founder_allele_counter
		if (mymother == NULL) {
			if (!update) {
				std::cout << "founder_allele_counter++" << std::endl;
				population->founder_allele_counter++;
				for (int i=1; i<= numLoci; i++){
					m_founder(i) = population->founder_allele_counter;
					p_founder(i) = population->founder_allele_counter+1;
				}
				population->founder_allele_counter++;
			}
		}
		// if not founder 
		else {
			Vector <unsigned> &mm_founder = mymother->get_founder("mother",update);
			Vector <unsigned> &mp_founder = mymother->get_founder("father",update);
			Vector <unsigned> &pm_founder = myfather->get_founder("mother",update);
			Vector <unsigned> &pp_founder = myfather->get_founder("father",update);  
			for (int i=1; i <= numLoci; i++){
				if (m_gamete(i) == 0) {
					m_founder(i) = mm_founder(i);
				}
				else {
					m_founder(i) = mp_founder(i);
				}
				if (p_gamete(i) == 0) {
					p_founder(i) = pm_founder(i);
				}
				else {
					p_founder(i) = pp_founder(i);
				}
			}
		} 
	}
	}

Vector <unsigned> &Individual::get_founder(std::string parent, bool update) {
	set_founder_alleles(update);
	
	if (parent == "mother") {
		return m_founder;
	}
	else {
		return p_founder;
	}
}

void Individual::update_offsprings_founder_alleles(void) {
	set_founder_alleles(1); //1 invokes updating
							//do all ofsprings
	if (numoffs) {
		for (int i=0; i<numoffs; i++) {
			myoffspring[i]->update_offsprings_founder_alleles();
		}
	}
}

double Individual::graph_trans_prob(unsigned uj,
									unsigned um,
									unsigned uf){
    double m_prob=0, f_prob=0;
    unsigned FQTable[4][4]={2,3,2,3,0,1,0,1,1,0,1,0,3,2,3,2};
    unsigned MQTable[4][4]={2,2,3,3,0,0,1,1,1,1,0,0,3,3,2,2};
    unsigned FQindex = FQTable[uf][uj];
    unsigned MQindex = MQTable[um][uj];
    unsigned MMindex = 0, FMindex = 0, size = m_gamete.size();
    
    //calc marker indexes (bin -> dec)
    for(int k = size; k > 0; k--){
		/*std::cout << "m_marker" << k << '\t' << m_gamete(k)
		<< "\tp_marker" << k << '\t' << p_gamete(k)
		<< std::endl;
		*/
		MMindex+=(unsigned)((1 << (size-k)))*m_gamete(k);
		FMindex+=(unsigned)((1 << (size-k)))*p_gamete(k);
    }
	
    if(FQindex==3 || MQindex==3){return 0.0;}
    if(FQindex==2){f_prob=1;}
    else{
		f_prob=population->model->gamete_prob_table[FMindex][FQindex];
		//std::cout << "FMindex: " << FMindex << "\tFQindex: " << FQindex << std::endl;
    }
    if(MQindex==2){m_prob=1;}
    else{
		m_prob=population->model->gamete_prob_table[MMindex][MQindex];
		//std::cout << "MMindex: " << MMindex << "\tMQindex: " << MQindex << std::endl;
    }
    //return transition (transmission*transmission) probability
    //std::cerr << "return transition " << f_prob*m_prob << std::endl;
    return f_prob*m_prob;
}


void Individual::t0(unsigned SLlocus, unsigned gender_parent) {
	//Sobel & Lange T0 rule
	if (gender_parent) {
		if (p_gamete(SLlocus)){p_gamete(SLlocus)=0;}
		else {p_gamete(SLlocus)=1;}
	}
	else {
		if (m_gamete(SLlocus)){m_gamete(SLlocus)=0;}
		else {m_gamete(SLlocus)=1;}
	}
	update_offsprings_founder_alleles();
	return;
}

void Individual::t1(unsigned SLlocus) {
	//Sobel & Lange T1 rule
	//get offsprings
	for (int i=0; i<numoffs; i++) {
		if (myoffspring[i]->mymother == this) {
			myoffspring[i]->t0(SLlocus, 0);
		}
		else {
			myoffspring[i]->t0(SLlocus, 1);
		}
	}
	update_offsprings_founder_alleles();
	return;
}

void Individual::t2a(unsigned SLlocus) {
	//Sobel & Lange T2a rule
	//choose spouse
	Individual *my_spouse = spouselist[unsigned(ranf() * numspouse)];
	//get offsprings
	for (int i=0; i<numoffs; i++) {
		//check if related with myspouse
		if (myoffspring[i]->myfather == my_spouse
			|| myoffspring[i]->mymother == my_spouse) {
			
			//do t0 for myoffspring for both alleles if they originate in the opposite gender
			if (myoffspring[i]->m_gamete(SLlocus) != myoffspring[i]->p_gamete(SLlocus)) {
				myoffspring[i]->t0(SLlocus,0);
				myoffspring[i]->t0(SLlocus,1);
			}
			
			//do t1 for myoffspring
			myoffspring[i]->t1(SLlocus);
		}
	}
	return;
}

void Individual::t2b(unsigned SLlocus) {
	//Sobel & Lange T2b rule
	//choose spouse
	Individual *my_spouse = spouselist[unsigned(ranf() * numspouse)];
	//get offsprings
	for (int i=0; i<numoffs; i++) {
		//check if related with myspouse
		if (myoffspring[i]->myfather == my_spouse
			|| myoffspring[i]->mymother == my_spouse) {
			
			//do t0 for myoffspring for both alleles if they originate in the same gender
			if (myoffspring[i]->m_gamete(SLlocus) == myoffspring[i]->p_gamete(SLlocus)) {
				myoffspring[i]->t0(SLlocus,0);
				myoffspring[i]->t0(SLlocus,1);
			}
			
			//do t1 for myoffspring
			myoffspring[i]->t1(SLlocus);
		}
	}
	return;
}   

void Individual::apply_SL_transition(const unsigned rule, const unsigned locus) {
	if (rule > 3) {
		throw exception("Population::apply_SL_transition: nonexistent rule!, abort");
	}
	switch (rule) {
		case 0: {	//T0 rule
			t0(locus,unsigned(ranf()*2)); //choose gender at random
			break;
		}
		case 1: {	//T1 rule
			t1(locus);
			break;
		}
		case 2: {	//T2a rule
			t2a(locus);
			break;
		}
		case 3: {	//T2b rule
			t2b(locus);
			break;
		}
	}
	return;
}

void Individual::apply_SL_cascade(const unsigned rule, const unsigned SLlocus) {
	//going from SL locus, we apply SL to the other loci depending on distance
	//first, go from SLlocus down to locus 1
	unsigned rule_applied = 1;
	for (unsigned locus = SLlocus - 1; locus > 0; locus--) {
		double u = ranf();
		if ((u > population->RecoVector(locus) && rule_applied) ||
			(u < population->RecoVector(locus) && !rule_applied)) { 
			//apply T rule
			apply_SL_transition(rule, locus);
			rule_applied = 1;
		}
		else {
			rule_applied = 0;
		}
	}
	
	//then go from SLlocus up to last locus
	rule_applied = 1;
	for (unsigned locus = SLlocus; locus < m_gamete.size(); locus++) {
		double u = ranf();
		if ((u > population->RecoVector(locus) && rule_applied) ||
			(u < population->RecoVector(locus) && !rule_applied)) { 
			//apply T rule
			apply_SL_transition(rule, locus + 1);
			rule_applied = 1;
		}
		else {
			rule_applied = 0;
		}
	}
	return;
}

//E. Thompson meiosis sampler

void Individual::sample_self(int parent_indicator){
	// Authors: Matthias Schelling and Rohan L. Fernando 
	// (1999) 
	// Contributors: L. Radu Totir 
	
	// parent_indicator = 0 -> maternal
	// parent_indicator = 1 -> paternal
	Vector <unsigned> *gamete_p;
	
	if(parent_indicator==0){
		gamete_p=&m_gamete;
	}
	else {
		gamete_p=&p_gamete;
	}
	
	// Computing the priors
	bool initSamplerToOrderFounders = false;
	for ( int j=1; j<=numLoci; j++){
		(*gamete_p)(j)=0;
		update_offsprings_founder_alleles();
		vec_prob[j-1][0]=population->calc_prior_descent_graph(j, initSamplerToOrderFounders);
		(*gamete_p)(j)=1;
		update_offsprings_founder_alleles();
		vec_prob[j-1][1]=population->calc_prior_descent_graph(j, initSamplerToOrderFounders); 
	}
	// Computing the "cutset" values
	int j=0;
	vec_cutsetval[j][0]=vec_prob[j][0];
	vec_cutsetval[j][1]=vec_prob[j][1];
	for ( int j=1; j<numLoci; j++){ 
		vec_cutsetval[j][0]= vec_prob[j][0]*
		(vec_cutsetval[j-1][0]*(1-population-> RecoVector(j))
		 + vec_cutsetval[j-1][1]*(population-> RecoVector(j)));
		vec_cutsetval[j][1]= vec_prob[j][1]*
			(vec_cutsetval[j-1][1]*(1-population-> RecoVector(j))
			 + vec_cutsetval[j-1][0]*(population-> RecoVector(j)));
	}
	
	// Sample the last locus from the marginal
	
	j=numLoci-1;
	double u = ranf();
	double p1 = vec_cutsetval[j][0];
	double p2 = vec_cutsetval[j][1];
	if (u<(p1/(p1+p2))){
		(*gamete_p)(j+1)=0;
	}
	else{
		(*gamete_p)(j+1)=1;
	}
	update_offsprings_founder_alleles();
	
	// sample the other loci
	
	for ( int j=numLoci-2; j>=0; j--){ 
		double u = ranf();
		double recomb=0.0;
		if((*gamete_p)(j+2)==0){
			recomb=(1-population-> RecoVector(j+1));
		}
		else{
			recomb=population-> RecoVector(j+1);
		}
		double p1= recomb*vec_cutsetval[j][0];
		double p2=(1-recomb)*vec_cutsetval[j][1];
		if (u<(p1/(p1+p2))){
			(*gamete_p)(j+1)=0;
		}
		else{
			(*gamete_p)(j+1)=1;
		}
		update_offsprings_founder_alleles();
	}
}

Vector<unsigned>  Individual::get_id_pdq(void){
	// Authors: Fabiano V. Pita and Rohan L. Fernando 
	// (2003) 
	// Contributors:
	return id_pdq;
}

void Individual::search_heteroz(unsigned qtl){
	// Authors: Fabiano V. Pita and Rohan L. Fernando 
	// (2003) 
	// Contributors: L. Radu Totir, Chris Stricker
	// This method is not optimal, use sampleAlleleVectors(), then findOptimumLocusToOrdder() instead.
	
	int mult = -1;           //  __1___2___3___qtl___5____6__
	unsigned step = 1;
	unsigned endl = 0;
	unsigned endr = 0;
	unsigned lcs = qtl;
	
	ord_heter = 0;  // 0 means individual is not an ordered heterozygous
	
	if (mymother != NULL){
		return;
	}
	
	// first goes to the closest left marker
	lcs += step * mult;
	mult *= -1;
	
	for (;;) {
		if (genome0.chromosome[0].locus[lcs-1].allele != genome1.chromosome[0].locus[lcs-1].allele){
			ord_heter = lcs;
			break;
		}
		
		if (lcs == 1){  //last locus of left side
			endl = 1;
		}
		
		if (lcs == numLoci){  //last locus of right side
			endr = 1;
		}
		
		step += 1;  
		lcs += step*mult;
		mult *= -1;
		
		if (lcs > numLoci){
			break;
		}
		
		//working just with the right side
		if (endl == 1){
			for (int i=lcs; i<= numLoci; i++){
				if (genome0.chromosome[0].locus[lcs-1].allele != genome1.chromosome[0].locus[lcs-1].allele){
					ord_heter = lcs;
					break;
				}
			}
			break;
		}
		
		//working just with the left side
		if (endr == 1){
			for (int i=lcs; i>=1; i--){
				if (genome0.chromosome[0].locus[lcs-1].allele != genome1.chromosome[0].locus[lcs-1].allele){
					ord_heter = lcs;
					break;
				}
			}
			break;
		}
	}
	
}

void Individual::get_allele_v1(unsigned lcs, Vector<int> allele_vector1){
  // Authors: Fabiano V. Pita and Rohan L. Fernando 
  // (2003) 
  // Contributors: 
  // used in conjuction with search_heteroz(), superseeded by sampleAlleleVector() and findOptimumLocusToOrder()
	if (ord_heter == lcs){
		initial_order1 = allele_vector1(m_founder(lcs));
		initial_order2 = allele_vector1(p_founder(lcs));
		std::cout<< id() <<" "<<lcs<<" "<<initial_order1<<" "<<initial_order2<< std::endl;
	}
}

void Individual::count_haplotype(void){
	// Authors: Fabiano V. Pita and Rohan L. Fernando 
	// (2003) 
	// Contributors: L. Radu Totir  
	if (m_haplotype.empty()){
		int zero = 0;
		m_haplotype.resize(4,numLoci-1,zero);
		p_haplotype.resize(4,numLoci-1,zero);
	}
	
	if (mymother == NULL) { // if founder
		for (int j=1; j<=numLoci-1;j++){
			m_haplotype(1,j) = 0;
			m_haplotype(2,j) = 0;
			m_haplotype(3,j) = 0;
			m_haplotype(4,j) = 0;
			p_haplotype(1,j) = 0;
			p_haplotype(2,j) = 0;
			p_haplotype(3,j) = 0;
			p_haplotype(4,j) = 0;
		}
	}
	else {
		for (int j=2; j<=numLoci; j++) {
			if ((m_gamete(j-1) == 0) && (m_gamete(j) == 0)) {
				m_haplotype(1,j-1)+=1;
			}      
			if ((m_gamete(j-1) == 0) && (m_gamete(j) == 1)) {
				m_haplotype(2,j-1)+=1;
			}
			if ((m_gamete(j-1) == 1) && (m_gamete(j) == 0)) {
				m_haplotype(3,j-1)+=1;
			}
			if ((m_gamete(j-1) == 1) && (m_gamete(j) == 1)) {
				m_haplotype(4,j-1)+=1;
			}
			if ((p_gamete(j-1) == 0) && (p_gamete(j) == 0)) {  
				p_haplotype(1,j-1)+=1;
			}      
			if ((p_gamete(j-1) == 0) && (p_gamete(j) == 1)) {
				p_haplotype(2,j-1)+=1;
			}
			if ((p_gamete(j-1) == 1) && (p_gamete(j) == 0)) {
				p_haplotype(3,j-1)+=1;
			}
			if ((p_gamete(j-1) == 1) && (p_gamete(j) == 1)) {
				p_haplotype(4,j-1)+=1;
			}
		}
	}
}  

// HG : to print the haplotype frequencies used in additive.cpp
void Individual::display_freq_haplotype(unsigned n)
{ 
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (2004) 
	// Contributors: L. Radu Totir, Chris Stricker
	double doublen = n;
	cout << std::setw(8)<<population->ind_name(myid);
	for (int j=1; j<numLoci; j++) {
		double m_haplo1 = m_haplotype(1,j);
		double m_haplo2 = m_haplotype(2,j);
		double m_haplo3 = m_haplotype(3,j);
		double m_haplo4 = m_haplotype(4,j);
		double p_haplo1 = p_haplotype(1,j);
		double p_haplo2 = p_haplotype(2,j);
		double p_haplo3 = p_haplotype(3,j);
		double p_haplo4 = p_haplotype(4,j);
		cout << std::setw(15)<<m_haplo1/doublen<< std::setw(8)<<m_haplo2/doublen<< std::setw(8)<<m_haplo3/doublen<< std::setw(8)<<m_haplo4/doublen<< std::setw(8)
			<< std::setw(8)<<p_haplo1/doublen<<std::setw(8)<<p_haplo2/doublen<< std::setw(8)<<p_haplo3/doublen<< std::setw(8)<<p_haplo4/doublen;
	}
	cout<<endl;
}

// Fabiano
// method for calculating the pdq from individual's haplotype probabilities (2 markers interval)
// mm_pdq = considering the individual maternal haplotype, pdq from the mother's maternal haplotype
// mp_pdq = considering the individual maternal haplotype, pdq from the mother's paternal haplotype 

void Individual::map_pdq(int samples, int interval, BG r_aq, BG r_bq){
	double doubleSample = samples;
	BG r_ab = (1-r_aq)*r_bq + r_aq*(1-r_bq);
	if (mymother == NULL) { // if founder
		bgDMPDQ = 0.0;
		bgSMPDQ = 0.0;    
	}
	else {  
		bgDMPDQ = (1-r_aq)*(1-r_bq) / (1-r_ab)     * m_haplotype(1,interval)/doubleSample 
		+ (1-r_aq)*(r_bq)   / r_ab         * m_haplotype(2,interval)/doubleSample
		+ (r_aq)  *(1-r_bq) / r_ab         * m_haplotype(3,interval)/doubleSample 
		+ (r_aq)*(r_bq)     / (1-r_ab)     * m_haplotype(4,interval)/doubleSample; 
		
		bgSMPDQ = (1-r_aq)*(1-r_bq) / (1-r_ab)     * p_haplotype(1,interval)/doubleSample 
            + (1-r_aq)*(r_bq)   / r_ab         * p_haplotype(2,interval)/doubleSample
            + (r_aq)  *(1-r_bq) / r_ab         * p_haplotype(3,interval)/doubleSample 
            + (r_aq)*(r_bq)     / (1-r_ab)     * p_haplotype(4,interval)/doubleSample; 
	}
}
	
void Individual::map_pdq(int samples, int interval, double r_aq, double r_bq){
	double doubleSample = samples;
	double r_ab = (1-r_aq)*r_bq + r_aq*(1-r_bq);
	if (mymother == NULL) { // if founder
		DMPDQ = 0.0;
		SMPDQ = 0.0;    
	}
	else if(r_ab < 0.00000000001){
		DMPDQ = m_haplotype(1,interval)/doubleSample;
		SMPDQ = p_haplotype(1,interval)/doubleSample;
	}
	else if(r_ab >.9999999999){  
		cout<<"ERROR: recombination rate between flanking markers of interval "<<interval<<" is virtually 1. Error in marker map. Please correct."<<endl;
		exit(123);
	}
	else {  
		DMPDQ = (1-r_aq)*(1-r_bq) / (1-r_ab)     * m_haplotype(1,interval)/doubleSample 
		+ (1-r_aq)*(r_bq)   / r_ab         * m_haplotype(2,interval)/doubleSample
		+ (r_aq)  *(1-r_bq) / r_ab         * m_haplotype(3,interval)/doubleSample 
		+ (r_aq)*(r_bq)     / (1-r_ab)     * m_haplotype(4,interval)/doubleSample; 
		
		SMPDQ = (1-r_aq)*(1-r_bq) / (1-r_ab)     * p_haplotype(1,interval)/doubleSample 
            + (1-r_aq)*(r_bq)   / r_ab         * p_haplotype(2,interval)/doubleSample
            + (r_aq)  *(1-r_bq) / r_ab         * p_haplotype(3,interval)/doubleSample 
            + (r_aq)*(r_bq)     / (1-r_ab)     * p_haplotype(4,interval)/doubleSample; 
//			cout<<"m_haplotype(1,interval)="<<m_haplotype(1,interval)<<" m_haplotype(2,interval)="<<m_haplotype(2,interval)<<" m_haplotype(3,interval)="<<m_haplotype(3,interval)<<" m_haplotype(4,interval)="<<m_haplotype(4,interval)<<endl;
//			cout<<"p_haplotype(1,interval)="<<p_haplotype(1,interval)<<" p_haplotype(2,interval)="<<p_haplotype(2,interval)<<" p_haplotype(3,interval)="<<p_haplotype(3,interval)<<" p_haplotype(4,interval)="<<p_haplotype(4,interval)<<endl;
	}
}

void Individual::pdq_grid(int samples, std::ofstream &outfile){
	
	double r_ab,r_aq,r_bq, ab, aq, bq;
	outfile<<population->ind_name(this->myid)<<"\t";
	for (int j=1; j<=numLoci-1; j++) {
		r_ab = population->RecoVector(j);
		ab = population->prior->get_distance(1,j+1) - population->prior->get_distance(1,j);
		aq = bq = ab/2;
		r_aq = population->model->MapF(aq);
		r_bq = r_aq;
		map_pdq(samples,j,r_aq, r_bq);
		if(mymother == NULL) {
			outfile<<-1<<"\t"<<-1<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<-1<<"\t"<<-1;
		}
		else {
			outfile<<DMPDQ<<"\t"<<1-DMPDQ<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<SMPDQ<<"\t"<<1-SMPDQ;
		}
		if(j!=numLoci-1) outfile<<"\t";
	}
	outfile<<endl;
}


void Individual::setFounderHaplotype(void)
{
	// Authors: L. Radu Totir 
	// (August, 2004) 
	// Contributors:
	if(mymother){
		throw exception("Individual::setFounderHaplotype(void) called for a non-founder");
	}
	else {
		for (unsigned j=0;j<numLoci;j++){
			unsigned patAllele = genome0.chromosome[0].locus[j].allele;
			unsigned matAllele = genome1.chromosome[0].locus[j].allele;
			if (matAllele > patAllele){
				genome0.chromosome[0].locus[j].allele=matAllele;
				genome1.chromosome[0].locus[j].allele=patAllele;
			}
		}
	}
}
/*! \fn double Individual::setFounderHaplotype(void)
*  \brief makes sure that in founders the maternal allele 
*  is asigned the lower allele state 
*/

double Individual::getAllelePenetrance(void)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2003) 
	// Contributors: 
	unsigned mState, pState;
	if(population->prior->chrom()[0].locus[currentLocus].qtl_ml=='q'){
		double mu=0.0,d,pen;
		double v = *residual_var[0][0];
		unsigned genoState, QLocus;
		if (myrecord[0].missing) {
			return 1.0;
		}
		for (unsigned i=0;i<nQTL;i++) {
			QLocus = QTLPosVector[i];
			mState    = malleleStateNodeVector[QLocus].getState();
			pState    = palleleStateNodeVector[QLocus].getState();
			mu+=population->prior->chrom()[0].locus[QLocus].genotypic_val_mat(pState+1,mState+1);
		}
		double y = record()[0].double_val();
		pen = isqrt2pi/std::sqrt(v)*std::exp(-(y-mu)*(y-mu));
		return pen;
	}
	else{
		unsigned allelePat = genome0.chromosome[0].locus[currentLocus].allele;
		unsigned alleleMat = genome1.chromosome[0].locus[currentLocus].allele;
		if(allelePat!=0 && allelePat!=alleleMat){
			mState    = malleleStateNodeVector[currentLocus].getMyAlleleState();
			pState    = palleleStateNodeVector[currentLocus].getMyAlleleState();
			if(mState != pState){
				return 1.0;
			}
			else {
				return 0.0;
			}
		}
		else {
			return 1.0;
		}
	}
}
/*! \fn double Individual::getAllelePenetrance(void)
*  \brief returns the allele penetrance for the RSampler
*/

double Individual::getAllelePenetrance(unsigned forLocus)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: 
	unsigned mState, pState;
	if(population->prior->chrom()[0].locus[forLocus].qtl_ml=='q'){
		double mu=0.0,d,pen;
		double v = *residual_var[0][0];
		unsigned genoState, QLocus;
		if (myrecord[0].missing) {
			return 1.0;
		}
		for (unsigned i=0;i<nQTL;i++) {
			QLocus = QTLPosVector[i];
			mState    = malleleStateNodeVector[QLocus].getState();
			pState    = palleleStateNodeVector[QLocus].getState();
			mu+=population->prior->chrom()[0].locus[QLocus].genotypic_val_mat(pState+1,mState+1);
		}
		double y = record()[0].double_val();
		pen = isqrt2pi/std::sqrt(v)*std::exp(-(y-mu)*(y-mu));
		return pen;
	}
	else{
		unsigned allelePat = genome0.chromosome[0].locus[forLocus].allele;
		unsigned alleleMat = genome1.chromosome[0].locus[forLocus].allele;
		if(allelePat!=0 && allelePat!=alleleMat){
			mState    = malleleStateNodeVector[forLocus].getMyAlleleState();
			pState    = palleleStateNodeVector[forLocus].getMyAlleleState();
			if(mState != pState){
				return 1.0;
			}
			else {
				return 0.0;
			}
		}
		else {
			return 1.0;
		}
	}
}
/*! \fn double Individual::getAllelePenetrance(void)
*  \brief returns the allele penetrance for the RSampler
*/

double Individual::getDisAllelePenetrance(void)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: 
	double mu=0.0,d,pen;
	unsigned mState, pState;
	if (myrecord[0].missing) {
		return 1.0;
	}
	mState    = malleleStateNodeVector[currentLocus].getState();
	pState    = palleleStateNodeVector[currentLocus].getState();
	int phenotype = int(record()[0].double_val());
	// cout << phenotype << "   " << mState+pState << "  " 
	// << population->disPenetranceTable[phenotype][mState+pState] 
	// << endl;
	return population->disPenetranceTable[phenotype][mState+pState];
}

double Individual::getDisAllelePenetrance(unsigned forLocus)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: 
	double mu=0.0,d,pen;
	unsigned mState, pState;
	if (myrecord[0].missing) {
		return 1.0;
	}
	mState    = malleleStateNodeVector[forLocus].getState();
	pState    = palleleStateNodeVector[forLocus].getState();
	int phenotype = int(record()[0].double_val());
	// cout << phenotype << "   " << mState+pState << "  " 
	// << population->disPenetranceTable[phenotype][mState+pState] 
	// << endl;
	return population->disPenetranceTable[phenotype][mState+pState];
}

double Individual::getGenoPenetrance(void)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors: 
	double mu=0.0,d,pen;
	double v = *residual_var[0][0];
	unsigned mState, pState, QLocus;
	if (myrecord[0].missing) {
		return 1.0;
	}
	for (unsigned i=0;i<nQTL;i++) {
		QLocus = QTLPosVector[i];
		mState    = genotNodeVector[QLocus].getmState();
		pState    = genotNodeVector[QLocus].getpState();
		mu+=population->prior->chrom()[0].locus[QLocus].genotypic_val_mat(pState+1,mState+1);
	}
	double y = record()[0].double_val();
	pen = isqrt2pi/std::sqrt(v)*std::exp(-(y-mu)*(y-mu));
	return pen;
}
/*! \fn double Individual::getGenoPenetrance(void)
*  \brief returns the genotype penetrance for the RSampler
*/

double Individual::getDisGenoPenetrance(void)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors: 
	double mu=0.0,d,pen;
	unsigned mState, pState;
	if (myrecord[0].missing) {
		return 1.0;
	}
	mState    = genotNodeVector[currentLocus].getmState();
	pState    = genotNodeVector[currentLocus].getpState();
	int phenotype = int(record()[0].double_val());
	// cout << phenotype << "   " << mState+pState << "  " 
	// << population->disPenetranceTable[phenotype][mState+pState] 
	// << endl;
	return population->disPenetranceTable[phenotype][mState+pState];
}

double Individual::getAllelePTransmissionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	// cout << currentLocus << endl;
	unsigned ip = palleleStateNodeVector[currentLocus].getMyAlleleState();
	unsigned fm = myfather->malleleStateNodeVector[currentLocus].getMyAlleleState();
	unsigned fp = myfather->palleleStateNodeVector[currentLocus].getMyAlleleState();
	// cout << ip << " " << fm << " " << fp <<  endl;
	return getTransmissionProb(ip,fm,fp,p_gamete);
}
/*! \fn double Individual::getAllelePTransmissionProb()
*  \brief returns the transmission probability of a paternal allele
in the RSampler
*/

double Individual::getAlleleMTransmissionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	unsigned im = malleleStateNodeVector[currentLocus].getMyAlleleState();
	// cout << "Maternal allele for individual: " << myid << " is " << im << endl;
	unsigned mm = mymother->malleleStateNodeVector[currentLocus].getMyAlleleState();
	unsigned mp = mymother->palleleStateNodeVector[currentLocus].getMyAlleleState();
	return getTransmissionProb(im,mm,mp,m_gamete);
}
/*! \fn double Individual::getAlleleMTransmissionProb()
*  \brief returns the transmission probability of a maternal allele
in the RSampler
*/

double Individual::getAllelePRTransmissionProb(unsigned forLocus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors:
	unsigned ipState = palleleStateNodeVector[forLocus].getMyAlleleState();
	unsigned fmState = myfather->malleleStateNodeVector[forLocus].getMyAlleleState();
	unsigned fpState = myfather->palleleStateNodeVector[forLocus].getMyAlleleState();
	unsigned ipOrigin = palleleOriginNodeVector[forLocus].getMyAlleleOrigin();
	// cout << ipState << " " << fmState << " " << fpState 
	// << " " << ipOrigin << endl;
	if(ipOrigin==0){
		if(ipState==fmState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	if(ipOrigin==1){
		if(ipState==fpState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
    else{
        cerr << "Error: Individual::getAllelePRTransmissionProb(): ipOrigin=" << ipOrigin << endl;
        exit(1);
    }
}
/*! \fn double Individual::getAllelePTransmissionProb()
*  \brief returns the transmission probability of a paternal allele 
when the allele origin is assumed known. Used when we peel across 
pedigree and loci jointly.
*/

double Individual::getAlleleMRTransmissionProb(unsigned forLocus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors:
	unsigned imState = malleleStateNodeVector[forLocus].getMyAlleleState();
	unsigned mmState = mymother->malleleStateNodeVector[forLocus].getMyAlleleState();
	unsigned mpState = mymother->palleleStateNodeVector[forLocus].getMyAlleleState();
	unsigned imOrigin = malleleOriginNodeVector[forLocus].getMyAlleleOrigin();
	//cout << imState << " " << mmState << " " << mpState 
	//     << " " << imOrigin << endl;
	if(imOrigin==0){
		if(imState==mmState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	if(imOrigin==1){
		if(imState==mpState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
    else{
        cerr << "Error: Individual::getAlleleMRTransmissionProb(): imOrigin=" << imOrigin << endl;
        exit(1);
    }
}
/*! \fn double Individual::getAlleleMTransmissionProb()
*  \brief returns the transmission probability of a maternal allele
when the allele origin is assumed known. Used when we peel across 
pedigree and loci jointly.
*/

double Individual::getTransitionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	double result;
	result = getMTransmissionProb()*getPTransmissionProb();
	if (result==0){
		result = -9999.0;
	}
	else {
		result = std::log(result);
	}
	return result;
}
/*! \fn double Individual::getTransitionProb()
*  \brief returns the genotype transition probability in the
RSampler
*/

double Individual::getPTransmissionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	unsigned ip = genotNodeVector[currentLocus].getpState();
	unsigned fm = myfather->genotNodeVector[currentLocus].getmState();
	unsigned fp = myfather->genotNodeVector[currentLocus].getpState();
	return getTransmissionProb(ip,fm,fp,p_gamete);
}
/*! \fn double Individual::getPTransmissionProb(void)
*  \brief returns the paternal allele transmission probability for
the RSampler for genotype peeling
*/

double Individual::getMTransmissionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	unsigned im = genotNodeVector[currentLocus].getmState();
	unsigned mm = mymother->genotNodeVector[currentLocus].getmState();
	unsigned mp = mymother->genotNodeVector[currentLocus].getpState();
	return getTransmissionProb(im,mm,mp,m_gamete);
}
/*! \fn double Individual::getMTransmissionProb(void)
*  \brief returns the maternal allele transmission probability for
the RSampler for genotype peeling
*/

double Individual::getTransmissionProb(unsigned iA, unsigned mA, unsigned pA,
									   Vector <unsigned> &gamete){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// new version,  November, 2004 RLF
	// Contributors:
	double temp = 1.0;
	if(mA==pA){
		if(iA==mA){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	int segregationIndex = gamete[currentLocus]; //LRT 6/28/04
	int	leftLocus  = findLeftLocus(gamete);
	int	rightLocus = findRightLocus(gamete);
	if(leftLocus==-1 && rightLocus==-1){
		if(iA==mA || iA==pA){
			return 0.5;
		}
		else {
			return 0.0;
		}
	}
	else {
		gamete[currentLocus] = 0;
		if (leftLocus!=-1) {
			if(gamete[leftLocus]==gamete[currentLocus]){
				temp*=1 - population->recombinationMatrix[leftLocus][currentLocus];
			}
			else {
				temp*=    population->recombinationMatrix[leftLocus][currentLocus];
			}
		}
		if(rightLocus!=-1){
			if(gamete[currentLocus]==gamete[rightLocus]){
				temp*=1 - population->recombinationMatrix[currentLocus][rightLocus];
			}
			else {
				temp*=    population->recombinationMatrix[currentLocus][rightLocus];
			}
		}
		double temp0 = temp;
		temp = 1.0;
		
		gamete[currentLocus] = 1;
		if (leftLocus!=-1) {
			if(gamete[leftLocus]==gamete[currentLocus]){
				temp*=1 - population->recombinationMatrix[leftLocus][currentLocus];
			}
			else {
				temp*=    population->recombinationMatrix[leftLocus][currentLocus];
			}
		}
		if(rightLocus!=-1){
			if(gamete[currentLocus]==gamete[rightLocus]){
				temp*=1 - population->recombinationMatrix[currentLocus][rightLocus];
			}
			else {
				temp*=    population->recombinationMatrix[currentLocus][rightLocus];
			}
		}
		double temp1 = temp;
		gamete[currentLocus] = segregationIndex; //LRT 6/28/2004
		double iSum  = 1.0/(temp0+temp1);
		if(iA==mA){
			return temp0*iSum;
		}
		else if(iA==pA){
			return temp1*iSum;
		}
		else {
			return 0.0;
		}
	}
}

	/*! \fn double Individual::getTransmissionProb(unsigned im, unsigned mm, unsigned mp, Vector<unsigned> &gamete)
*   \brief computes and returns the transmission probability for a given trio of alleles taking in account flanking loci in the
RSampler
*/


double Individual::getMatthiasTransitionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors:
	double result;
	result = getMatthiasMTransmissionProb()*getMatthiasPTransmissionProb();
	return result;
}
/*! \fn double Individual::getMatthiasTransitionProb()
*  \brief returns the genotype transition probability for the proposal
*  in the RSampler. The cascading origin trick of Matthias Schelling.
*/

double Individual::getMatthiasPTransmissionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors:
	unsigned ip = genotNodeVector[currentLocus].getpState();
	unsigned fm = myfather->genotNodeVector[currentLocus].getmState();
	unsigned fp = myfather->genotNodeVector[currentLocus].getpState();
	return getMatthiasTransmissionProb(ip,fm,fp,p_gamete,p_gameteOld);
}
/*! \fn double Individual::getMatthiasPTransmissionProb(void)
*  \brief returns the paternal allele transmission probability for
*   the proposal in the RSampler for genotype peeling. 
*/

double Individual::getMatthiasMTransmissionProb(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors:
	unsigned im = genotNodeVector[currentLocus].getmState();
	unsigned mm = mymother->genotNodeVector[currentLocus].getmState();
	unsigned mp = mymother->genotNodeVector[currentLocus].getpState();
	return getMatthiasTransmissionProb(im,mm,mp,m_gamete,m_gameteOld);
}
/*! \fn double Individual::getMatthiasMTransmissionProb(void)
*  \brief returns the maternal allele transmission probability for
*   the proposal in the RSampler for genotype peeling.
*/

double Individual::getMatthiasTransmissionProb(unsigned iA, unsigned mA, unsigned pA,
											   Vector <unsigned> &gamete, 
											   Vector <unsigned> &gameteOld){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors:
	double temp = 1.0;
	double rCOLeft,rCORight;
	if(mA==pA){
		if(iA==mA){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	int segregationIndex = gamete[currentLocus]; //LRT 6/28/2004
	int leftLocus,rightLocus=-1;
	if(gamete.size()>1){
		leftLocus  = findLeftLocus(gamete);
		rightLocus = findRightLocus(gamete);
	}
	else{
		leftLocus = -1;
	}
	if(leftLocus==-1){
		if(iA==mA || iA==pA){
			return 0.5;
		}
		else {
			return 0.0;
		}
	}
	else {
		if(gameteOld[leftLocus]==gameteOld[currentLocus]){
			rCOLeft = population->recombinationMatrix[leftLocus][currentLocus];
		}
		else {
			rCOLeft = 1 - population->recombinationMatrix[leftLocus][currentLocus];
		}
		gamete[currentLocus] = 0;
		if(currentLocus>0){
			if(gamete[leftLocus]==gamete[currentLocus]){
				temp *= 1 - rCOLeft;
			}
			else {
				temp *= rCOLeft;
			}
		}
		if(currentLocus<gamete.size()){
			if(!(population->lookToYourLeft) && rightLocus!=-1){
				if(gamete[currentLocus]==gamete[rightLocus]){
					temp*=1 - population->recombinationMatrix[currentLocus][rightLocus];
				}
				else {
					temp*=population->recombinationMatrix[currentLocus][rightLocus];
				}
			}
		}
		double temp0 = temp;
		temp = 1.0;
		gamete[currentLocus] = 1;
		if(currentLocus>0){
			if(gamete[leftLocus]==gamete[currentLocus]){
				temp *= 1 - rCOLeft;
			}
			else {
				temp *= rCOLeft;
			}
		}
		if(currentLocus<gamete.size()){
			if(!(population->lookToYourLeft) && rightLocus!=-1){
				if(gamete[currentLocus]==gamete[rightLocus]){
					temp*=1 - population->recombinationMatrix[currentLocus][rightLocus];
				}
				else {
					temp*=population->recombinationMatrix[currentLocus][rightLocus];
				}
			}
		}
		gamete[currentLocus] = segregationIndex; //LRT 6/28/2004
		double temp1 = temp;
		double iSum  = 1.0/(temp0+temp1);
		if(iA==mA){
			return temp0*iSum;
		}
		else if(iA==pA){
			return temp1*iSum;
		}
		else {
			return 0.0;
		}
	}
	return temp;
}
/*! \fn double Individual::getMatthiasTransmissionProb(unsigned im, 
unsigned mm, unsigned mp, Vector <unsigned> &gamete, 
Vector <unsigned> &gameteOld)
*  \brief computes and returns the transmission probability for a
given trio of alleles taking in account flanking loci in the
RSampler
*/

int Individual::findLeftLocus(Vector <unsigned> &gamete){ 
	if (!population->lookToYourLeft) return -1;
	for (int jj = currentPosition-1; jj>=0; jj--){
		unsigned i = population->prior->chrom()[0].locusMask[jj].index; 
		if(population->prior->chrom()[0].locus[i].qtl_ml!='q'){
			if(gamete[i]!=9999){
				return i;
			}
		}
	}
	return -1;
}
/*! \fn inline int Individual::findLeftLocus(Vector <unsigned> &gamete)
*  \brief finds the first heterozygous marker locus to the left
of the current locus 
*/

int Individual::findRightLocus(Vector <unsigned> &gamete){ 
	if (!population->lookToYourRight) return -1;
	for (unsigned jj = currentPosition+1; jj < numLoci ; jj++){
		unsigned i = population->prior->chrom()[0].locusMask[jj].index; 
		if(population->prior->chrom()[0].locus[i].qtl_ml!='q'){
			if(gamete[i]!=9999){
				return i;
			}
		}
	}
	return -1;
}
/*! \fn inline int Individual::findRightLocus(Vector <unsigned> &gamete)
*  \brief finds the first heterozygous marker locus to the right 
of the current locus
*/

ostream& operator<<(ostream& stream, Individual& i) {
	std::vector<unsigned> matHaplotype, patHaplotype;
	matHaplotype.resize(i.numLoci);
	patHaplotype.resize(i.numLoci);
	for (unsigned j=0; j<i.numLoci; j++){  
		//     if (i.population->samplerType=="genotypic"){
		matHaplotype[j] = i.genotNodeVector[j].getmState()+1;
		patHaplotype[j] = i.genotNodeVector[j].getpState()+1;
		//     }
		//     else if (i.population->samplerType=="allelic"){
		//       matHaplotype[j] = i.malleleStateNodeVector[j].getMyAlleleState();
		//       patHaplotype[j] = i.palleleStateNodeVector[j].getMyAlleleState();
		//     }
	}
	for (unsigned j=0; j<i.numLoci; j++){  
		stream << matHaplotype[j] << "  ";
	}
	stream << endl;
	for (unsigned j=0; j<i.numLoci; j++){  
		stream << "-" << "  ";
	}
	stream << endl;
	for (unsigned j=0; j<i.numLoci; j++){  
		stream << patHaplotype[j] << "  ";
	}
	stream << endl;
	return stream;
}

void Individual::displayGenotypes(unsigned locus){
	unsigned numGenos = genotNodeVector[locus].genotypeVector.size();
	cout << "Individual: " << myid << endl;
	for (unsigned i=0;i<numGenos;i++){
		cout << genotNodeVector[locus].genotypeVector[i].maternal+1;
		cout << "/";
		cout << genotNodeVector[locus].genotypeVector[i].paternal+1 << endl;
	}
	cout << "genotypeState: " << genotNodeVector[locus].genotypeState << endl;
}

void Individual::displayAlleleVectors(unsigned locus){
	cout << "Individual: " << myid << " Locus: " << locus << endl;	
	cout << "maternal allele states" << endl;
	unsigned sizeOfVector = malleleStateNodeVector[locus].alleleStateVector.size();
	for (unsigned i=0; i<sizeOfVector; i++){
		cout << malleleStateNodeVector[locus].alleleStateVector[i] << endl;
	}
	if (mymother){
		cout << "maternal allele origin" << endl;
		sizeOfVector = malleleOriginNodeVector[locus].alleleOriginVector.size();
		for (unsigned i=0; i<sizeOfVector; i++){
			cout << malleleOriginNodeVector[locus].alleleOriginVector[i] << endl;
		}
	}
	cout << "paternal allele states" << endl;
	sizeOfVector = palleleStateNodeVector[locus].alleleStateVector.size();
	for (unsigned i=0; i<sizeOfVector; i++){
		cout << palleleStateNodeVector[locus].alleleStateVector[i] << endl;
	}
	if(mymother){
		cout << "paternal allele origin" << endl;
		sizeOfVector = palleleOriginNodeVector[locus].alleleOriginVector.size();
		for (unsigned i=0; i<sizeOfVector; i++){
			cout << palleleOriginNodeVector[locus].alleleOriginVector[i] << endl;
		}
	}
}

void Individual::setAlleleStateVectors(unsigned locus){
	set<unsigned> matAlleles, patAlleles;
	malleleStateNodeVector[locus].alleleStateVector.clear();
	palleleStateNodeVector[locus].alleleStateVector.clear();	
	unsigned numGenos = genotNodeVector[locus].genotypeVector.size();
	for (unsigned i=0;i<numGenos;i++){
		matAlleles.insert(genotNodeVector[locus].genotypeVector[i].maternal+1);
		patAlleles.insert(genotNodeVector[locus].genotypeVector[i].paternal+1);
	}
	set<unsigned>::iterator it;
	for (it=matAlleles.begin();it!=matAlleles.end();it++){
		malleleStateNodeVector[locus].alleleStateVector.push_back(*it);
	}
	for (it=patAlleles.begin();it!=patAlleles.end();it++){
		palleleStateNodeVector[locus].alleleStateVector.push_back(*it);
	}		
}

void Individual::setAlleleOriginVectors(unsigned locus){
	SafeSTLVector<unsigned> intersectionVector;
	unsigned msize, psize;
	
	set_intersection(myfather->malleleStateNodeVector[locus].alleleStateVector.begin(),
					 myfather->malleleStateNodeVector[locus].alleleStateVector.end(),
					 palleleStateNodeVector[locus].alleleStateVector.begin(),
					 palleleStateNodeVector[locus].alleleStateVector.end(),
					 inserter(intersectionVector,intersectionVector.begin())  );
	msize = intersectionVector.size();
	intersectionVector.clear();
	set_intersection(myfather->palleleStateNodeVector[locus].alleleStateVector.begin(),
					 myfather->palleleStateNodeVector[locus].alleleStateVector.end(),
					 palleleStateNodeVector[locus].alleleStateVector.begin(),
					 palleleStateNodeVector[locus].alleleStateVector.end(),
					 inserter(intersectionVector,intersectionVector.begin())  );
	psize = intersectionVector.size();
	
	if (msize==0 && psize>0) {
		palleleOriginNodeVector[locus].alleleOriginVector.clear();
		palleleOriginNodeVector[locus].alleleOriginVector.push_back(1);
	}
	if(msize>0 && psize==0) {
		palleleOriginNodeVector[locus].alleleOriginVector.clear();
		palleleOriginNodeVector[locus].alleleOriginVector.push_back(0);
	}
	if (msize == 0 && psize ==0){cout << "locus : " << locus; throw exception(" Possible error in Genotypes");}
	
	intersectionVector.clear();
	set_intersection(mymother->malleleStateNodeVector[locus].alleleStateVector.begin(),
					 mymother->malleleStateNodeVector[locus].alleleStateVector.end(),
					 malleleStateNodeVector[locus].alleleStateVector.begin(),
					 malleleStateNodeVector[locus].alleleStateVector.end(),
					 inserter(intersectionVector,intersectionVector.begin())  );
	msize = intersectionVector.size();
	intersectionVector.clear();
	set_intersection(mymother->palleleStateNodeVector[locus].alleleStateVector.begin(),
					 mymother->palleleStateNodeVector[locus].alleleStateVector.end(),
					 malleleStateNodeVector[locus].alleleStateVector.begin(),
					 malleleStateNodeVector[locus].alleleStateVector.end(),
					 inserter(intersectionVector,intersectionVector.begin())  );
	psize = intersectionVector.size();
	if (msize==0 && psize>0) {
		malleleOriginNodeVector[locus].alleleOriginVector.clear();
		malleleOriginNodeVector[locus].alleleOriginVector.push_back(1);
	}
	if(msize>0 && psize==0) {
		malleleOriginNodeVector[locus].alleleOriginVector.clear();
		malleleOriginNodeVector[locus].alleleOriginVector.push_back(0);
	}
	if (msize == 0 && psize ==0){cout << "locus: " << locus; throw exception(" Possible error in Genotypes");}
}

void Individual::setOwnerGNodes(void){
	// Authors: Rohan L. Fernando 
	// (November 2005) 
	// Contributors: 
	for (unsigned i=0;i<Individual::numLoci;i++){
		genotNodeVector[i].owner         = this;
		malleleStateNodeVector[i].owner  = this;
		palleleStateNodeVector[i].owner  = this;
		if(mymother){
			malleleOriginNodeVector[i].owner = this;
			palleleOriginNodeVector[i].owner = this;
		}

	}
}

void Individual::setSegregationIndex(unsigned atLocus,string samplerType){
	// Authors: Rohan L. Fernando 
	// (November 2005) 
	// Contributors: 
	Individual *ind, *mom, *dad;
	unsigned im,ip,mm,mp,fm,fp;
	ind=this;
	mom=ind->mymother;
	dad=ind->myfather;
	if(mom){
		if (samplerType=="genotypic"){
			im = ind->genotNodeVector[atLocus].getmState();
			ip = ind->genotNodeVector[atLocus].getpState();
			mm = mom->genotNodeVector[atLocus].getmState();
			mp = mom->genotNodeVector[atLocus].getpState();
		}
		else if (samplerType=="allelic"){
			im = ind->malleleStateNodeVector[atLocus].getMyAlleleState();
			ip = ind->palleleStateNodeVector[atLocus].getMyAlleleState();
			mm = mom->malleleStateNodeVector[atLocus].getMyAlleleState();
			mp = mom->palleleStateNodeVector[atLocus].getMyAlleleState();	
		}
		if (mm==mp){
			ind->m_gamete[atLocus] =9999;
		}
		else if (im==mm){
			ind->m_gamete[atLocus] = 0;
		}
		else if (im==mp){
			ind->m_gamete[atLocus] = 1;
		}
		else {
			cerr << "Error in: Individual::setSegregationIndex(unsigned atLocus,string samplerType)" << endl;
			exit(1);
		}
		if (samplerType=="genotypic"){
			fm = dad->genotNodeVector[atLocus].getmState();
			fp = dad->genotNodeVector[atLocus].getpState(); 
		}
		else if (samplerType=="allelic"){
			fm = dad->malleleStateNodeVector[atLocus].getMyAlleleState();
			fp = dad->palleleStateNodeVector[atLocus].getMyAlleleState(); 
		}
		if (fm==fp){
			ind->p_gamete[atLocus] =9999;
		}
		else if (ip==fm){
			ind->p_gamete[atLocus] = 0;
		}
		else if (ip==fp){
			ind->p_gamete[atLocus] = 1;
		}
		else {
			cerr << "Error in: Individual::setSegregationIndex(unsigned atLocus,string samplerType)" << endl;
			exit(1);
		}
	}
}
void Individual::sampleSegregationIndicators(void){
	// Authors: Rohan L. Fernando 
	// (November 2005) 
	// Contributors: 
	Individual *ind, *mom;
	ind=this;
	mom=ind->mymother;
	if(mom){
		for (unsigned atMarker = 0; atMarker < Individual::numLoci; atMarker++){
		    Individual::currentPosition = atMarker;
			if (ind->m_gamete[atMarker] == 9999){
				double u = ranf();
				int RightLocus = ind->findRightLocus(ind->m_gamete);
				int LeftLocus  = ind->findLeftLocus(ind->m_gamete);
				if ( (RightLocus != -1) && (LeftLocus != -1) ) {
					double recLeft  = population->recombinationMatrix[LeftLocus][atMarker];
					double recRight = population->recombinationMatrix[atMarker][RightLocus];
					double recFlank = population->recombinationMatrix[LeftLocus][RightLocus];
					if (ind->m_gamete[LeftLocus] == ind->m_gamete[RightLocus]){
						ind->m_gamete[atMarker] = 
						u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->m_gamete[LeftLocus] : 1 - ind->m_gamete[LeftLocus];
					}
					else{
						ind->m_gamete[atMarker] = 
						u <((  recLeft)*(1-recRight))/(  recFlank) ? ind->m_gamete[RightLocus] : ind->m_gamete[LeftLocus];
					}
					
				}
				else if (RightLocus != -1) {	 
					double recRight = population->recombinationMatrix[atMarker][RightLocus];
					ind->m_gamete[atMarker] = u < (1-recRight) ? ind->m_gamete[RightLocus] : 1-ind->m_gamete[RightLocus] ;
				}
				else if (LeftLocus != -1) {	
					double recLeft = population->recombinationMatrix[LeftLocus][atMarker];	
					ind->m_gamete[atMarker] = u < (1-recLeft) ? ind->m_gamete[LeftLocus] : 1-ind->m_gamete[LeftLocus];
				}
				else {
					ind->m_gamete[atMarker] = u < 0.5 ? 1 : 0 ;
				}
			}
			if (ind->p_gamete[atMarker] == 9999){
				double u = ranf();
				int RightLocus = ind->findRightLocus(ind->p_gamete);
				int LeftLocus  = ind->findLeftLocus(ind->p_gamete);
				if ( (RightLocus != -1) && (LeftLocus != -1) ) {
					double recLeft  = population->recombinationMatrix[LeftLocus][atMarker];
					double recRight = population->recombinationMatrix[atMarker][RightLocus];
					double recFlank = population->recombinationMatrix[LeftLocus][RightLocus];
					if (ind->p_gamete[LeftLocus] == ind->p_gamete[RightLocus]){
						ind->p_gamete[atMarker] = 
						u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->p_gamete[LeftLocus] : 1 - ind->p_gamete[LeftLocus];
					}
					else{
						ind->p_gamete[atMarker] = 
						u <((  recLeft)*(1-recRight))/(  recFlank) ? ind->p_gamete[RightLocus] : ind->p_gamete[LeftLocus];
					}
					
				}
				else if (RightLocus != -1) {	 
					double recRight = population->recombinationMatrix[atMarker][RightLocus];
					ind->p_gamete[atMarker] = u < (1-recRight) ? ind->p_gamete[RightLocus] : 1-ind->p_gamete[RightLocus] ;
				}
				else if (LeftLocus != -1) {	
					double recLeft = population->recombinationMatrix[LeftLocus][atMarker];	
					ind->p_gamete[atMarker] = u < (1-recLeft) ? ind->p_gamete[LeftLocus] : 1-ind->p_gamete[LeftLocus];
				}
				else {
					ind->p_gamete[atMarker] = u < 0.5 ? 1 : 0 ;
				}
			}
		}
	}
}

double Individual::getGenotypeMRTransmissionProb(unsigned forLocus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (July, 2006) 
	// Contributors:
	
	unsigned imState = genotNodeVector[forLocus].getmState();
	unsigned mmState = mymother->genotNodeVector[forLocus].getmState();
	unsigned mpState = mymother->genotNodeVector[forLocus].getpState();
	unsigned imOrigin = malleleOriginNodeVector[forLocus].getMyAlleleOrigin();
	
	//cout << imState << " " << mmState << " " << mpState 
	//     << " " << imOrigin << endl;
	if(imOrigin==0){
		if(imState==mmState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	if(imOrigin==1){
		if(imState==mpState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
    else{
        cerr << "Error: Individual::getGenotypeMRTransmissionProb(): imOrigin=" << imOrigin << endl;
        exit(1);
    }
}
/*! \fn double Individual::getGenotypeMTransmissionProb()
*  \brief returns the transmission probability of a maternal allele
when the allele origin is assumed known. Used when we peel across 
pedigree and loci jointly.
*/

double Individual::getGenotypePRTransmissionProb(unsigned forLocus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (July, 2006) 
	// Contributors:
	unsigned ipState = genotNodeVector[forLocus].getpState();
	unsigned fmState = myfather->genotNodeVector[forLocus].getmState();
	unsigned fpState = myfather->genotNodeVector[forLocus].getpState();
	unsigned ipOrigin = palleleOriginNodeVector[forLocus].getMyAlleleOrigin();
	// cout << ipState << " " << fmState << " " << fpState 
	// << " " << ipOrigin << endl;
	if(ipOrigin==0){
		if(ipState==fmState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
	if(ipOrigin==1){
		if(ipState==fpState){
			return 1.0;
		}
		else {
			return 0.0;
		}
	}
    else{
        cerr << "Error: Individual::getGenotypePRTransmissionProb(): ipOrigin=" << ipOrigin << endl;
        exit(1);
    }
}
/*! \fn double Individual::getGenotypePTransmissionProb()
*  \brief returns the transmission probability of a paternal allele 
when the allele origin is assumed known. Used when we peel across 
pedigree and loci jointly.
*/

double Individual::getDisGenoPenetrance(unsigned forLocus)
{
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors: 
	double mu=0.0,d,pen;
	unsigned mState, pState;
	if (myrecord[0].missing) {
		return 1.0;
	}
	mState    = genotNodeVector[forLocus].getmState();
	pState    = genotNodeVector[forLocus].getpState();
	int phenotype = int(record()[0].double_val());

	return population->disPenetranceTable[phenotype][mState+pState];
}

void Individual::calcDistanceToNeighbors(SafeSTLVector<unsigned>& individualsToBeProcessed){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
    Individual* ind;
	if (mymother){
		if (mymother->distanceToPivot0==-1){
			neighborBlock.push_back(mymother);
			mymother->distanceToPivot0 = distanceToPivot0 + 1;
			individualsToBeProcessed.push_back(mymother->myid - 1);
		}
	}
	if (myfather){
		if (myfather->distanceToPivot0==-1){
			neighborBlock.push_back(myfather);
			myfather->distanceToPivot0 = distanceToPivot0 + 1;
			individualsToBeProcessed.push_back(myfather->myid - 1);
		}
	}
	for (unsigned i=0;i<numspouse;i++){
		ind = spouselist[i];
		if (ind->distanceToPivot0==-1){
			neighborBlock.push_back(ind);
			ind->distanceToPivot0 = distanceToPivot0 + 1;
			individualsToBeProcessed.push_back(ind->myid - 1);
		}
	}
	for (unsigned i=0;i<numoffs;i++){
		ind = myoffspring[i];
		if (ind->distanceToPivot0==-1){
			neighborBlock.push_back(ind);
			ind->distanceToPivot0 = distanceToPivot0 + 1;
			individualsToBeProcessed.push_back(ind->myid - 1);
		}
	}
}

void Individual::getPedBlockUpto(unsigned endDist){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
	SafeSTLVector<Individual*>::iterator first, last;
	if (neighborBlock.size()==0) return;
	pedBlock = neighborBlock;
	if (!(distanceToPivot0 < endDist-1)) return;
	for (unsigned i=0;i<neighborBlock.size();i++){
		neighborBlock[i]->getPedBlockUpto(endDist);
	}
	for (unsigned i=0;i<neighborBlock.size();i++){
		if (neighborBlock[i]->pedBlock.size()){
			first = neighborBlock[i]->pedBlock.begin();
			last  = neighborBlock[i]->pedBlock.end();
			pedBlock.insert(pedBlock.end(),first,last);
			neighborBlock[i]->pedBlock.clear();
		}
	}
}

bool Individual::pedBlockAdjacent(SafeSTLVector<Individual*> jPedBlock){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
	for (unsigned i=0;i<pedBlock.size();i++){
		Individual* indi = pedBlock[i];
		for (unsigned j=0;j<jPedBlock.size();j++){
			Individual* indj = jPedBlock[j];
			if (indi->isNeighborOf(indj)) return true;
		}
	}
	return false;
}

bool Individual::isNeighborOf(Individual* indj){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
	if (indj == mymother || indj == myfather) return true;
	for (unsigned i=0;i<numspouse;i++){
		if (indj==spouselist[i]) return true;
	}
	for (unsigned i=0;i<numoffs;i++){
		if (indj==myoffspring[i]) return true;
	}
	return false;
}

} /// end of matvec namespace
