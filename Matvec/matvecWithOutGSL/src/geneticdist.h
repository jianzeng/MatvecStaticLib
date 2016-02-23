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

#ifndef GeneticDist_H
#define GeneticDist_H

#include <iostream>
#include <string>
#include "util.h"
#include "doublematrix.h"
#include "safe_vectors.h"

namespace matvec {
extern int DISCRETE_TRAIT;  // was called  DISCRETE_TRAIT;
extern doubleMatrix PENETRANCE_MATRIX;

class LocusMask {
	public:
	unsigned index;
	double position;
	bool operator < (const LocusMask X) const{
		if(position < X.position) {
		   return 1;
		}
		else return 0;
	  };
};

/////////////////////////// LocusStruct class ///////////////////////

class ChromStruct;
class GeneticDist;

/*!
   \class LocusStruct GeneticDist.h
   \brief a structure for a locus

   \sa Locus
*/
class LocusStruct {

   protected:
      unsigned   numallele;
      void       copyfrom(const LocusStruct& A);

   public:
      string myName;
      char       qtl_ml;
	  string gnodeType;
      Vector<double>     allele_freq;
      doubleMatrix     genotypic_val_mat;
      double     distance; // distance from origin to this locus in Morgans
      std::string   nameOfcol1, nameOfcol2;              // for marker data
      unsigned nextMarker;
      unsigned nHaplotypes ;
      Vector<unsigned> listOfAllele;            // setup in population::ListAlleleFounders
      Vector<unsigned> al1Haplo, al2Haplo;      // setup in population::SetPossibleHaplotypes : related to the next marker on the right
	  bool originProbs, stateProbs;
	  bool samplePositionForThisLocus;
      LocusStruct(void) { 
		numallele = 0;
		samplePositionForThisLocus=false;
		originProbs = stateProbs = false;
		myName = "unspecified ";
	  }
      LocusStruct(const unsigned na) {numallele = 0; resize(na);}
//      ~LocusStruct(void) {numallele =0;}

      const LocusStruct& operator=(const LocusStruct& A)
           {numallele=0; copyfrom(A); return *this;}

      void     resize(const unsigned na);
      void     nallele(const unsigned na) {numallele = na;}
      unsigned nallele(void) const {return numallele;}
      // LRT for the RSampler
      SafeSTLVector<int> peelOrder;
	  Vector<unsigned> gameteMask;
	  bool cutLoops;
	  bool operator < (const LocusStruct L) const{
		if(distance < L.distance) {
		   return 1;
		}
		else return 0;
	  };
};

/////////////////////////// ChromStruct class ///////////////////////
/*!
   \class ChromStruct  GeneticDist.h
   \brief chromosome structure

   \sa Chromosome
*/
class ChromStruct {
   protected:
      unsigned numloci;
      void     copyfrom(const ChromStruct& A);

   public:
      SafeSTLVector<LocusStruct> locus;
      doubleMatrix       recomb_rate_mat;

      ChromStruct(void) {numloci=0;}
      ChromStruct(const unsigned nl) {numloci=0; resize(nl);}
      ~ChromStruct(void) {release();}

      const ChromStruct& operator=(const ChromStruct& A)
           {numloci=0; copyfrom(A); return *this;}

      void      resize(const unsigned nl);
      unsigned  nloci(void) const {return numloci;}

      void      release(void);
      // LRT for the RSampler
      SafeSTLVector<int> peelOrder;
	  double myLength;
	  SafeSTLVector<LocusMask> locusMask;
	  void calcLocusMask();
};

class MaternalPaternalRQTLAlleles {
 public:
  unsigned maternal;
  unsigned paternal;
};
/*! \class MaternalPaternalRQTLAlleles geneticdist.h inc/geneticdist.h
 * \brief This is the building block class for dealing with RQTL.
 *
 * MaternalPaternalRQTLAlleles class defines the maternal and paternal 
 * alleles of a RQTL genotype.
 */ 

/////////////////////////// GeneticDist class ///////////////////////

/*!
   \class GeneticDist  GeneticDist.h
   \brief  a base for Genetic distributions

   \sa StatDist
*/

class GeneticDist  {
   public:
      char         distname[25];
      unsigned     numchrom, numtrait;     // default numtrait = 1
	  static double HaldaneMToR(double M){
			 return 0.5*(1.0 - std::exp(-2.0*M));
	  }		 
	  ChromStruct  *chromosome;
      unsigned numMarkerLoci; //LRT 
      GeneticDist(void);
      GeneticDist(const unsigned nc,...);
      GeneticDist(GeneticDist& G ){copyfrom(G);}
      virtual ~GeneticDist(void) {release();}

      const GeneticDist& operator=(const GeneticDist& A) {numchrom=0;
                  chromosome=0; copyfrom(A); return *this;}


      virtual void     resize(const unsigned nc) {nchrom(nc); }
      virtual unsigned ntrait(void) const {return numtrait;}
      virtual doubleMatrix*  var_matrix(void) { return 0;}

      const char* name(void) const {return distname;}
      ChromStruct  *chrom(void) {return chromosome;}
      double  recomb_rate(const unsigned c,const unsigned li,const unsigned lj);
      unsigned     nchrom(void) const {return numchrom;}
      unsigned     nloci_chrom(const unsigned c) const;
      unsigned     nallele(const unsigned c,const unsigned l) const
                           {return chromosome[c].locus[l].nallele();}

      //int     display(void)const{cout<<"  "<<distname<<"\n"; return 1;}
      int     display(void);

      void    copyfrom(const GeneticDist& A);
      void    nchrom(const unsigned nc);
      void    nloci(const unsigned nl0,...);
      void    release(void);
      void    locus(const unsigned c,const unsigned l, string type,
                     const unsigned na, ...);
     void putColmNames(const unsigned c,const unsigned l,
                                char nm1[],char nm2[]);
	     void putColmNames(const unsigned c,const unsigned l,
                                string nm1,string nm2);							
      void    recomb_rate(const unsigned c,const unsigned li,const unsigned lj,
                           const double r);
      void    genotypic_val(const unsigned c,const unsigned l,const double* v);
      void    genotypic_val(const unsigned c,const unsigned l,
                             const double v0,...);
      double  sample(void) const {return 0.0;}
      const double** genotypic_val(const unsigned c,const unsigned l) const;
      void multi_loci(int num_loci);
      void put_distance(const unsigned c,const unsigned l, double distance);
	  void putGnodeType(const unsigned c,const unsigned l, string type="allelic");
	  void calcProbs(const unsigned c,const unsigned l, string type="origin");
      double get_distance(const unsigned c,const unsigned l);
      void    locus(const unsigned c,const unsigned l, string type, const unsigned na, Vector<double> allele_freq);
	  std::vector<unsigned> markerVec;
	  void createMarkerVec(const unsigned c, const unsigned l);
      SafeSTLVector<MaternalPaternalRQTLAlleles> rqtlVector; 
	  void putChromLength(const unsigned chromID, double length){chromosome[chromID-1].myLength=length;} //LRT
	  void setFlagForSamplingPosition(const unsigned chromID, const unsigned locusID, bool flag)
		{chromosome[chromID-1].locus[locusID-1].samplePositionForThisLocus=flag;}
	  void putLocusName(const unsigned chromID, const unsigned locusID, string name){
		chromosome[chromID-1].locus[locusID-1].myName = name;
	  }	
};

inline const double** GeneticDist::genotypic_val(const unsigned c,
                                                 const unsigned l) const
{
   return (const double**)chromosome[c].locus[l].genotypic_val_mat.begin();
}

inline GeneticDist::GeneticDist(void)
{
   numchrom = 0; chromosome = 0; numtrait = 1;
    strcpy(distname,"GeneticDist");
//   cout << "GeneticDist::GeneticDist(void): numtrait = " << numtrait << endl;
}

inline unsigned  GeneticDist::nloci_chrom(const unsigned c) const
{
   if (c > 0) {
      return chromosome[c-1].nloci();
   }
   else {
      std::cerr << "GeneticDist::nloci(c): arg value out of range";
      exit(1);
   }
}


///////////////////  UnknownDist  class /////////////////////////////
/*!
   \class UnknownDist
   \brief a dummy genetic distribution
*/
class UnknownDist: public GeneticDist {
   protected:
      unsigned dim;
      Vector<double> mu_vec;
      doubleMatrix var_mat;

   public:
      UnknownDist(void) {strcpy(distname,"UnknownDist"); dim=0;}
      UnknownDist(const UnknownDist& u):GeneticDist() {copyfrom(u);}
      UnknownDist(const Vector<double> mu, const doubleMatrix sigma);
      UnknownDist(const double mu, const double sigma);
      ~UnknownDist(void) {;}

      void     copyfrom(const UnknownDist& A);
      void     resize(const unsigned n);
      unsigned ntrait(void) const {return dim;}

      Vector<double>*  mean_vector(void) {return &mu_vec;}
      doubleMatrix*  var_matrix(void) {return &var_mat;}
};
}
#endif


