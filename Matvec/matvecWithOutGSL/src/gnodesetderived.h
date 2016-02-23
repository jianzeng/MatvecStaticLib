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

#ifndef gnodesetderived_h
#define gnodesetderived_h
#include <ext/hash_map>
#include <map>
#include <set>
#include "exception.h"
#include "gnodestuff.h"
#include "individual.h"
#include "safe_vectors.h"

namespace matvec {
  
  class CutSet:public GNodeSet{
  public:
    SafeSTLVector<GNode*> elementVector;
    SafeSTLVector<unsigned> keyMultiplicationCode;
	SafeSTLVector<double> valueVector;
    void setKeyMultiplicationCode(void);
	void setupElementVector(void);
    inline int getKey(void);
    inline int getOldKey(void);
    double getValue(void);
    double getOldValue(void);
    void   putValue(double x);
    void   displayValues(void);
    GNode* getLocation();
    void   normalize();
    void   replaceMeWith(CutSet& A);
    CutSet& operator = (GNodeSet* A);
    CutSet& operator+= (GNodeSet* A);
    CutSet& operator*= (GNodeSet* A);
    ~CutSet(void) {;}
  };
  /*! \class CutSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "cut set" class.
   *
   * CutSet class has methods relevant to a "cut set"
   */
  
  class AlleleFounderSet:public GNodeSet{
  public:
    double getValue(void){
      AlleleFounderSet::iterator it;
      it = begin();
      unsigned allele = (*it)->getState();
      double prob = prior->chrom()[0].locus[currentLocus].allele_freq[allele];
	  double logProb;
	  if (prob==0){
		logProb = -9999.0;
	  }
	  else{
		logProb = std::log(prob);
	  }
	  return logProb;
    }
    ~AlleleFounderSet(void) {;} 
  };
  /*! \class AlleleFounderSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Allele Founder set" class.
   *
   * AlleleFounderSet class has the method that returns founder 
   probabilities for alleles
   */
  
  class GenoFounderSet:public GNodeSet{
  public:
    double getValue(void){
      GenoFounderSet::iterator it;
      it = begin();
      unsigned mallele = (*it)->getmState();
      unsigned pallele = (*it)->getpState();
      double pm = prior->chrom()[0].locus[currentLocus].allele_freq[mallele];
      double pp = prior->chrom()[0].locus[currentLocus].allele_freq[pallele];
      double prob = pm*pp;
	  double logProb;
	  if (prob==0){
		logProb = -9999.0;
	  }
	  else{
		logProb = std::log(prob);
	  }
	  return logProb;
    }
    ~GenoFounderSet(void) {;} 
  };
  /*! \class GenoFounderSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Genotype Founder set" class.
   *
   * GenoFounderSet class has the method that returns founder 
   probabilities for genotypes
   */ 
  
  class TransmissionSet:public GNodeSet{
public:
	  Individual *offspring;
	  bool paternal;
	  double prob;
	  double getValue(void){
		  if (paternal){
			  prob = offspring->getAllelePTransmissionProb();
		  }
		  else {
			  prob = offspring->getAlleleMTransmissionProb();
		  }
		  double logProb;
		  if (prob==0){
			  logProb = -9999.0;
		  }
		  else{
			  logProb = std::log(prob);
		  }
		  return logProb;
	  }
	  ~TransmissionSet(void) {;} 
  };
  /*! \class TransmissionSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Transmission set" class.
   *
   * TransmissionSet class has the method that returns transmission 
   probabilities (alleles)
   */
  
  class TransitionSet:public GNodeSet{
  // log is done in individual.getTransitionProb 
    public:
	  static int transmissionType;
	  Individual *offspring;
	  double getValue(void){
		  switch(transmissionType){
			  case 1:
				  return offspring->getTransitionProb();
				  break;
			  case 2:
				  return offspring->getMatthiasTransitionProb();
				  break;
			  default:
				  throw exception("Proposal type in TransitionSet getValue(void) is non-existent");
		  }
		  return 1.0;
	  }
	  double getTargetValue(void){
		  return offspring->getTransitionProb();
	  }
	  ~TransitionSet(void) {;}
  };
  /*! \class TransitionSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Transition set" class.
   *
   * TransitionSet class has the methods that return transition 
   probabilities (genotypes) for the proposal (the cascading origin 
   trick of Matthias Schelling) and for the target
   */
  
  class GenoPenetranceSet:public GNodeSet{
  public:
    Individual *owner;
	  double getValue(void){
		  double prob = owner->getGenoPenetrance();
		  double logProb;
		  if (prob==0){
			  logProb = -9999.0;
		  }
		  else{
			  logProb = std::log(prob);
		  }
		  return logProb;
	  }
    ~GenoPenetranceSet(void) {;}
  };
  /*! \class GenoPenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Genotype penetrance set" class.
   *
   *  GenoPenetranceSet class has the method that calculates the 
   penetrance function for a quantitative trait in genotype peeling 
   */ 
  
  class DisGenoPenetranceSet:public GNodeSet{
  public:
    Individual *owner;
	  double getValue(void){
		  double prob = owner->getDisGenoPenetrance();
		  double logProb;
		  if (prob==0){
			  logProb = -9999.0;
		  }
		  else{
			  logProb = std::log(prob);
		  }
		  return logProb;
	  }
    ~DisGenoPenetranceSet(void) {;}
  };
  /*! \class DisGenoPenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Disease Genotype penetrance set" class.
   *
   *  DisGenoPenetranceSet class has the method that calculates the 
   penetrance function for a recessive trait in genotype peeling 
   */  

  class AllelePenetranceSet:public GNodeSet{
  public:
    Individual *owner;
	  double getValue(void){
		  double prob = owner->getAllelePenetrance();
		  double logProb;
		  if (prob==0){
			  logProb = -9999.0;
		  }
		  else{
			  logProb = std::log(prob);
		  }
		  return logProb;
	  }
    ~AllelePenetranceSet(void) {;}
  };
  /*! \class AllelePenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Allele penetrance set" class.
   *
   *  AllelePenetranceSet class has the method that calculates the 
   penetrance function for allelic peeling
   */ 

  class DisAllelePenetranceSet:public GNodeSet{
  public:
    Individual *owner;
	  double getValue(void){
		  double prob = owner->getDisAllelePenetrance();
		  double logProb;
		  if (prob==0){
			  logProb = -9999.0;
		  }
		  else{
			  logProb = std::log(prob);
		  }
		  return logProb;
	  }
    ~DisAllelePenetranceSet(void) {;}
  };
  /*! \class DisAllelePenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Disease Allele penetrance set" class.
   *
   *  DisAllelePenetranceSet class has the method that calculates the 
   penetrance function for a recessive disease trait in allelic peeling
  */ 

  class RAlleleFounderSet:public GNodeSet{
  public:
    unsigned forLocus;
    double getValue(void){
      RAlleleFounderSet::iterator it;
      it = begin();
      unsigned allele = (*it)->getState();
	  double prob = prior->chrom()[0].locus[forLocus].allele_freq[allele];
	  double logProb;
	  if (prob==0){
		logProb = -9999.0;
	  }
	  else{
		logProb = std::log(prob);
	  }
	  return logProb;
    }
    ~RAlleleFounderSet(void) {;} 
  };
  /*! \class RAlleleFounderSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "RAllele Founder set" class.
   *
   * AlleleFounderSet class has the method that returns founder 
   probabilities for alleles when the locus passed as an argument.
  */
  
  class RAllelePenetranceSet:public GNodeSet{
  public:
    unsigned forLocus;
    Individual *owner;
    double getValue(void){
		double prob = owner->getAllelePenetrance(forLocus);
		double logProb;
		if (prob==0){
			logProb = -9999.0;
		}
		else{
			logProb = std::log(prob);
		}
		return logProb;
		
    }
    ~RAllelePenetranceSet(void) {;}
  };
  /*! \class RAllelePenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "RAllele penetrance set" class.
   *
   *  RAllelePenetranceSet class has the method that calculates the 
   penetrance function for allelic peeling for the locus passed as an 
   argument.
   */ 

  class RDisAllelePenetranceSet:public GNodeSet{
  public:
    unsigned forLocus;
    Individual *owner;
    double getValue(void){
		double prob = owner->getDisAllelePenetrance(forLocus);
		double logProb;
		if (prob==0){
			logProb = -9999.0;
		}
		else{
			logProb = std::log(prob);
		}
		return logProb;
    }
    ~RDisAllelePenetranceSet(void) {;}
  };
  /*! \class RDisAllelePenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "RDisAllele penetrance set" class.
   *
   *  RDisAllelePenetranceSet class has the method that calculates the
   *  penetrance function a recessive disease trait for allelic
   *  peeling with the disease locus passed as an argument.
   */ 

  class RecombinationSet:public GNodeSet{
public:
	  double r;
	  double getValue(void){
	      double prob;
		  RecombinationSet::iterator it;
		  it = begin();
		  //cout << "Size of recombination set is = " << size() << endl;
		  //cout << "Elements of recombination set are: ";
		  //for(it=begin();it!=end();it++){
		  //cout << (*it)->id << " " << (*it)->getState() << " ";
		  //}
		  //cout << endl;
		  it = begin();
		  if(size()==1){
			  prob = 0.5;
		  }
		  else if(size()==2){
			  unsigned segOne = (*it)->getMyAlleleOrigin();
			  it++;
			  unsigned segTwo = (*it)->getMyAlleleOrigin();
			  if(segOne==segTwo){
				  prob = 1-r;
			  }
			  else{
				  prob = r;
			  }
		  }
		  else{
			  cerr << "Trouble: I have found more than two GNodes in a RecombinationSet" << endl;
			  exit(1);
		  }
		  double logProb;
		  if (prob==0){
			  logProb = -9999.0;
		  }
		  else{
			  logProb = std::log(prob);
		  }
		  return logProb;
	  }
	  ~RecombinationSet(void) {;} 
  }; 
  /*! \class RecombinationSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Recombination set" class.
   *
   * RecombinationSet class is used when we peel across pedigree and loci 
   jointly.
  */   
  
  class RTransmissionSet:public GNodeSet{
  public:
    unsigned forLocus;
    Individual *offspring;
    bool paternal;
    double getValue(void){
		double prob;
		if (paternal){
			prob = offspring->getAllelePRTransmissionProb(forLocus);
		}
		else {
			prob = offspring->getAlleleMRTransmissionProb(forLocus);
		}
		double logProb;
		if (prob==0){
			logProb = -9999.0;
		}
		else{
			logProb = std::log(prob);
		}
		return logProb;
    }
    ~RTransmissionSet(void) {;} 
  };
  /*! \class RTransmissionSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "RTransmission set" class.
   *
   * RTransmissionSet class has the method that returns transmission 
   probabilities (alleles) for known parental origin. Used when we peel 
   across pedigree and loci jointly.
   */
   
  class RTransitionSet:public GNodeSet{
  public:
    unsigned forLocus;
    Individual *offspring;
    double getValue(void){
		double patTrPr = offspring->getGenotypePRTransmissionProb(forLocus);
		double matTrPr = offspring->getGenotypeMRTransmissionProb(forLocus);
		double prob = patTrPr*matTrPr;
		double logProb;
		if (prob==0){
			logProb = -9999.0;
		}
		else{
			logProb = std::log(prob);
		}
		return logProb;
	}
    ~RTransitionSet(void) {;} 
  };
  /*! \class RTransitionSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "RTransition set" class.
   *
   * RTransition class has the method that returns genotype probability of an 
   * individual given the genotypes of the parents and the allele origins of the
   * individual.  Used when we peel across pedigree and loci jointly.
   */
   
  
  class RGenoFounderSet:public GNodeSet{
  public:
    double getValue(void){
      RGenoFounderSet::iterator it;
      it = begin();
      unsigned mallele = (*it)->getmState();
      unsigned pallele = (*it)->getpState();
      double pm = (*it)->locus->allele_freq[mallele];
      double pp = (*it)->locus->allele_freq[pallele];
	  double prob = pm*pp;
	  double logProb;
	  if (prob==0){
		logProb = -9999.0;
	  }
	  else{
		logProb = std::log(prob);
	  }
	  return logProb;
    }
    ~RGenoFounderSet(void) {;} 
  };
  /*! \class RGenoFounderSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Genotype Founder set" class.
   *
   * RGenoFounderSet class has the method that returns founder 
   * probabilities for genotypes.  Used when we peel across pedigree and loci jointly.
   */ 

  class RDisGenoPenetranceSet:public GNodeSet{
  public:
    Individual *owner;
	unsigned forLocus;
    double getValue(void){
		double prob = owner->getDisGenoPenetrance(forLocus);
		double logProb;
		if (prob==0){
			logProb = -9999.0;
		}
		else{
			logProb = std::log(prob);
		}
		return logProb;
    }
    ~RDisGenoPenetranceSet(void) {;}
  };
  /*! \class RDisGenoPenetranceSet gnodesetderived.h inc/gnodesetderived.h
   * \brief This is the "Disease Genotype penetrance set" class.
   *
   *  RDisGenoPenetranceSet class has the method that calculates the 
   *  penetrance function for a recessive trait in genotype peeling 
   *  Used when we peel across pedigree and loci jointly.
   */  



}////// end of namespace matvec 
#endif 
