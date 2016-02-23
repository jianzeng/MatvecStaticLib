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

#ifndef gnodederived_h
#define gnodederived_h
#include <iostream>
#include <ext/hash_map>
#include <map>
#include <set>
#include "exception.h"
#include "gnodestuff.h"
#include "safe_vectors.h"

namespace matvec {
	class AlleleStateNode:public GNode{
public:
		unsigned alleleState,acceptedAlleleState,candidateAlleleState;
		SafeSTLVector<unsigned> alleleStateVector;
		SafeSTLVector<unsigned> acceptedAlleleStateVector;  
		AlleleStateNode(void){
			alleleStateVector.name = "allelStateVector";
			acceptedAlleleStateVector.name = "acceptedAlleleStateVector";
			alleleState = 0;
		};
		unsigned getState(void) {
			return alleleState;
		};
		unsigned getOldState(void) {
			return acceptedAlleleState;
		};
		unsigned getAcceptedAlleleState(void){
			return alleleStateVector[acceptedAlleleState]-1;
		};
		unsigned getMyAlleleState(void){
			return alleleStateVector[alleleState]-1;
		};
		unsigned getWeight(void) {
			return alleleStateVector.size();
		};
		void reset(int i){
			alleleState = i;
		};
		void resetAndSwitch(int i){
			acceptedAlleleState = alleleState; 
			alleleState = i;
		};
		void saveState(void){
			acceptedAlleleState = alleleState; 
		};
		void recoverSavedState(void){
			alleleState = acceptedAlleleState;
		};
		bool incr();
		void updateCounter(void){
			throw exception("AlleleStateNode:upddateCounter() not implimented yet \n");
		};
		void EEUpdateCounter(void){
			throw exception("AlleleStateNode:EEUpddateCounter() not implimented yet \n");
		};
		void updateProbs(CutSet& myGenotProb){
			throw exception("AlleleStateNode:myGenotProbs(CutSet&) not implimented yet \n");
		};
		void resetStateandVector(int i){
			unsigned element = alleleStateVector[i];
			alleleStateVector.clear();
			alleleStateVector.push_back(element);
			//		  alleleStateVector.push_back(i+1);
			alleleState = 0;
		};
		void display(void);
		~AlleleStateNode(void) {;}
	};
	/*! \class AlleleStateNode gnodederived.h inc/gnodederived.h
   * \brief This is the allele state node class.
   * 
   * AlleleStateNode class has methods relevant to an allele state node
   */ 
  class AlleleOriginNode:public GNode{
  public:
    unsigned alleleOrigin,acceptedAlleleOrigin;
    SafeSTLVector<unsigned> alleleOriginVector;
	AlleleOriginNode(void){
		alleleOriginVector.name = "allelOriginVector";
		alleleOrigin = 0;
	};

    unsigned getState(void) {
      return alleleOrigin;
    };
    unsigned getOldState(void) {
      return acceptedAlleleOrigin;
    };
    unsigned getMyAlleleOrigin(void){
      return alleleOriginVector[alleleOrigin];
    };
	unsigned getAcceptedAlleleOrigin(void){
      return alleleOriginVector[acceptedAlleleOrigin];
    };

    unsigned getWeight(void) {
      return alleleOriginVector.size();
    };
    void reset(int i){
      alleleOrigin = i;
    };
    void resetAndSwitch(int i){
      acceptedAlleleOrigin = alleleOrigin; 
      alleleOrigin = i;
    };
	void updateCounter(void){
		throw exception("AlleleOriginNode:upddateCounter() not implimented yet \n");
	};
	void EEUpdateCounter(void){
		throw exception("AlleleOriginNode:EEUpddateCounter() not implimented yet \n");
	};
	void resetStateandVector(int i){
		unsigned element = alleleOriginVector[i];
	   alleleOriginVector.clear();
//	   alleleOriginVector.push_back(i);
       alleleOriginVector.push_back(element);
	   alleleOrigin = 0;
	};
	void updateProbs(CutSet& myGenotProb){
		throw exception("AlleleOriginNode:myGenotProbs(CutSet&) not implimented yet \n");
	};
	void saveState(void){
		acceptedAlleleOrigin = alleleOrigin; 
	};
	void recoverSavedState(void){
		alleleOrigin = acceptedAlleleOrigin;
	};
    bool incr();
    ~AlleleOriginNode(void) {;}
  };
  /*! \class AlleleOriginNode gnodederived.h inc/gnodederived.h
   * \brief This is the allele origin node class.
   * 
   * AlleleOriginNode class has methods relevant to an allele origin node
   */ 
  class MaternalPaternalAlleles {
  public:
    unsigned maternal;
    unsigned paternal;
	bool operator< (const MaternalPaternalAlleles y) const;
  };
  /*! \class MaternalPaternalAlleles gnodederived.h inc/gnodederived.h
   * \brief This is the building block class for the GenotypeNode class.
   *
   * MaternalPaternalAlleles class defines the maternal and paternal 
   * alleles of a genotype.
   */ 
  class GenotypeNode:public GNode{
public:
	  unsigned genotypeState,candidateGenotypeState,acceptedGenotypeState;
	  GenotypeNode(void){
		  genotypeVector.name = "genotypeVector";
		  acceptedGenotypeVector.name = "acceptedGenotypeVector";
	  }
	  SafeSTLVector<MaternalPaternalAlleles> genotypeVector;
	  SafeSTLVector<float> genotypeCount;
	  unsigned sampleCount;
	  SafeSTLVector<MaternalPaternalAlleles> acceptedGenotypeVector;
	  
	  inline unsigned getState(void) {
		  return genotypeState;
	  };
	  unsigned getOldState(void) {
		  return acceptedGenotypeState;
	  };
	  unsigned getAcceptedMatState(void){
		  return acceptedGenotypeVector[acceptedGenotypeState].maternal;
	  };
	  unsigned getAcceptedPatState(void){
		  return acceptedGenotypeVector[acceptedGenotypeState].paternal;
	  };
	  unsigned getmState(void){
		  return genotypeVector[genotypeState].maternal;
	  };
	  unsigned getpState(void){
		  return genotypeVector[genotypeState].paternal;
	  };
	  unsigned getWeight(void) {
		  return genotypeVector.size();
	  };
	  void reset(int i){
		  genotypeState = i;
	  };
	  void saveState(void){
		acceptedGenotypeState = genotypeState; 
	  };
	  void recoverSavedState(void){
		genotypeState = acceptedGenotypeState;
	  };
	  void resetAndSwitch(int i){
		  acceptedGenotypeState = genotypeState; 
		  acceptedGenotypeVector = genotypeVector;
		  genotypeState = i;
	  };
	  void updateCounter(void);
	  void EEUpdateCounter(void);
	  void updateProbs(CutSet& myGenotProb);
	  bool incr();
	  ~GenotypeNode(void) {;}
	  void resetStateandVector(int i){
		  MaternalPaternalAlleles element = genotypeVector[i];
		  genotypeVector.clear();
		  genotypeVector.push_back(element);
		  genotypeState = 0;
	  };
	  void calcGenotypeProbs(void);
	  void display(void);  
  };
  /*! \class GenotypeNode gnodederived.h inc/gnodederived.h
   * \brief This is the genotype node class.
   *
   * GenotypeNode class has methods relevant to a genotype node
   */ 
}////// end of namespace matvec 
#endif
