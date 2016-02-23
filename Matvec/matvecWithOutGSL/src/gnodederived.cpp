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

#include <algorithm>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ext/hash_map>
#include <map>
#include <set>
#include "safe_vectors.h"
#include "gnodederived.h"
#include "gnodesetderived.h"
#include "gnodestuff.h"

namespace matvec {
  
inline int CutSet::getOldKey(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors: Robert Benson
  CutSet::iterator it;
  int key = 0;
  CutSet::iterator beginIt = begin();
  CutSet::iterator endIt = end();
  unsigned i=0;
  for (it=beginIt;it!=endIt;it++){
    key += ((*it)->getOldState())*keyMultiplicationCode[i];
    i++;
  }
  return key;
}
/*! \fn int CutSet::getKey(void) 
 *  \brief returns the string key to be used to put or retrieve
           values from the valueVector
*/

double CutSet::getOldValue(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors: 
  int key = getOldKey();
  return valueVector[key];
}
/*! \fn double CutSet::getOldValue(void)
 *  \brief returns the relevant probability from the valueVector
 */  

inline int CutSet::getKey(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors: Robert Benson
  int key = 0;
  
  vector<GNode*>::iterator it;
  vector<GNode*>::iterator beginIt = elementVector.begin();
  vector<GNode*>::iterator endIt = elementVector.end();
  
//  CutSet::iterator it;
//  CutSet::iterator beginIt = begin();
//  CutSet::iterator endIt   = end();
  
  unsigned i=0;
  for (it=beginIt;it!=endIt;it++){
    key += ((*it)->getState())*keyMultiplicationCode[i];
    i++;
  }
  return key;
}
/*! \fn int CutSet::getKey(void) 
 *  \brief returns the string key to be used to put or retrieve
 values from the valueVector
*/

double CutSet::getValue(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
		int key = getKey();
		if(key > valueVector.size()){
			cout << "key > size in CutSet::getValue\n";
			exit(1);
		}
		return valueVector[key];
}
/*! \fn double CutSet::getValue(void)
 *  \brief returns the relevant probability from the valueVector
 */  

void CutSet::putValue(double x) {
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	int key = getKey();
	if(key > valueVector.size()){
		cout << "key > size in CutSet::putValue\n";
		exit(1);
	}
	valueVector[key] = x;
}
/*! \fn void CutSet::putValue(double x)
 *  \brief puts the value x in the valueVector
 */ 
 
  void CutSet::setupElementVector(void){
	CutSet::iterator it;
	elementVector.clear();
	for (it=begin();it!=end();it++){
		elementVector.push_back(*it);
	}
} 
/*! \fn void CutSet::setupElementVector
 *  \brief copy set elements to a vector for speed
 */ 
  
void CutSet::setKeyMultiplicationCode(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (July, 2004) 
  // Contributors:
  setupElementVector();
  unsigned mySize = (size()) ? size():1;
  keyMultiplicationCode.resize(mySize);
  CutSet::iterator it;
  CutSet::iterator beginIt = begin();
  CutSet::iterator endIt = end();
  unsigned state = 1;
  unsigned i = 0;
  for (it=beginIt;it!=endIt;it++){
	keyMultiplicationCode[i] = (*it)->sampled ? 0 : state;
    i++;
    state *= (*it)->sampled ? 1 : (*it)->getWeight();
  }
  valueVector.resize(state,0.0);
}
/*! \fn void CutSet::setKeyMultiplicationCode(double x) 
 * \brief sets the multiplication code used to create the key to store 
 *        and retreive CutSet probabilities from valueVector. The code 
 *        for the first GNode in the CutSet is set to 1, and as a general 
 *        rule, the code for the n'th GNode in the cutset is set to the 
 *        product between the code of (n -1)'th GNode multiplied by the 
 *        number of possible GNode (genotype or allele) states of GNode n-1. 
 *        Also we resize valueVector and setup elementVector
 */  

GNode* CutSet::getLocation(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
  CutSet::iterator it;
  set<GNodeSet*>::iterator itFound;
  CutSet::iterator beginIt = begin();
  CutSet::iterator endIt = end();
  for (it=beginIt;it!=endIt;it++){
    itFound = (*it)->SetofGNsts.find(this);
    if(itFound!=(*it)->SetofGNsts.end()){
      return *it;
    }
  }
  display();
  throw exception("CutSet::getLocation failed for CutSet displayed \n");
}
  
void CutSet::normalize(){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (February, 2004) 
	// Contributors: Chris Stricker   
	// We are switching from the log scale to unlog'ed scale.
	double sum = 0.0;
	double temp= 1.0, bigValue=-9999.0;
	
	for(unsigned i=0;i!=valueVector.size();i++){
//		cout<<"valueVector["<<i<<"]="<<valueVector[i]<<endl;

		if(valueVector[i] > bigValue) bigValue = valueVector[i];
	}
	for(unsigned i=0;i!=valueVector.size();i++){
		valueVector[i] = std::exp(valueVector[i] - bigValue);		
		sum += valueVector[i];
//		cout << "sum after summing " << i+1 << " value = " << sum << endl;
	}
	if (sum < 0.00000000001) cout << "Normalize: sum = " << sum << endl;
	if(sum!=0.0){
		temp = 1.0/sum;
	}
//	cout << "Sum in normalize() = " << sum << endl;
	if(sum>1e-40){ // original value -14, set to -20 by me (LRT 9/24/04)
				   //       cout << "Unscaled probabilities: " << endl;
				   //       displayValues();
		for(unsigned i=0;i!=valueVector.size();i++){
			valueVector[i]*=temp;
//			cout<<"valueVector["<<i<<"]="<<valueVector[i]<<endl;
		} 
	}
	else if(sum==0.0){
		throw matvec::InvalidSample();
	}
	else{
		cout << "Sum used for scaling = " << sum << endl;
		CutSet::iterator itC;
		cout << "CutSet Members: [";
		for (itC=begin();itC!=end();itC++){
			cout << (*itC)->id << " ";
		} 
		cout << "]" << endl;
		displayValues();
		throw matvec::NormalizeFailed();
		//cerr << " Need scaling in CutSet::normalize(), sum = " << sum << endl;
		//exit(1);
	}
}

void CutSet::displayValues(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (February, 2004) 
  // Contributors:
  for(unsigned i=0;i!=valueVector.size();i++){
    cout << valueVector[i] << endl;
  } 
}
/*! \fn void CutSet::displayValues(void)
 *  \brief displays the values stored in the valueVector
*/  

void CutSet::replaceMeWith(CutSet& A){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  GNodeSet::iterator it;
  connectFlag = A.connectFlag;
  numberOfCuts = A.numberOfCuts;
  prior = A.prior;
  currentLocus = A.currentLocus;
  clear();
  for (it=A.begin();it!=A.end();it++){
    insert(*it);
  }
  keyMultiplicationCode = A.keyMultiplicationCode;
  valueVector = A.valueVector;
  elementVector = A.elementVector;
}
  
CutSet& CutSet::operator= (GNodeSet* A){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:

//  reset();
//  do{
//    double Atemp = A->getValue();
//    (*this).putValue(Atemp);
//  } while(incr());
//  return *this;
//

//	// Authors: L. Radu Totir and Rohan L. Fernando 
//	// (June, 2003) 
//	// Contributors:
	reset();
	unsigned key = 0;
	if(dynamic_cast<CutSet*>( A)){
		CutSet* cutsetPtr = reinterpret_cast<CutSet*>(A);
		do{
			//unsigned key = getKey();
			unsigned keyA = cutsetPtr->getKey();
			double Atemp = cutsetPtr->valueVector[keyA];
			valueVector[key] = Atemp;
			key++;
		} while(incr());
	} 
	else {
		do{
			//unsigned key = getKey();
			double Atemp = A->getValue();
			valueVector[key] = Atemp;
			key++;
		} while(incr());
	}
//	cout << "after = \n";
//	displayValues();	
	return *this;
}
/*! \fn CutSet& CutSet::operator= (GNodeSet* A)
 *  \brief overloaded = operator, puts the values contained by the
     larger set (the neighborhood set) into the new set.
*/  

CutSet& CutSet::operator+= (GNodeSet* A){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:

	for(unsigned i=0; i<valueVector.size(); i++){
		valueVector[i] = 0.0;
	}
	A->reset();
	unsigned keyA = 0;
  	if(dynamic_cast<CutSet*> (A)){
		CutSet* cutsetPtr = reinterpret_cast<CutSet*>(A);
		double bigValue = -9999.0, temp;
		for(unsigned i=0; i<cutsetPtr->valueVector.size();i++){
			temp = cutsetPtr->valueVector[i];
			if(temp > bigValue){
				bigValue = temp;
			}
		}
		do{
			unsigned key = getKey();
			//unsigned keyA = cutsetPtr->getKey();
			double Atemp = std::exp(cutsetPtr->valueVector[keyA] - bigValue);  //Scaling
			valueVector[key] += Atemp;
			keyA++;
		} while(A->incr());
		
		for(unsigned i=0; i<valueVector.size();i++){
			double temp = valueVector[i];
			if(temp < 0.000000000000001) valueVector[i] = -9999.0;
			else valueVector[i] = std::log(temp);
		}
	} 
	else {					// no need to scale GNodesets
		do{
			unsigned key = getKey();
			double Atemp = A->getValue();
			valueVector[key] += std::exp(Atemp);
		}  while(A->incr());
		for(unsigned i=0; i<valueVector.size();i++){
			double temp = valueVector[i];
			if(temp < 0.000000000000001) valueVector[i] = -9999.0;
			else valueVector[i] = std::log(temp);
		}
	}
//	cout << "after += \n";
//	displayValues();	
	
	return *this;
}
/*! \fn CutSet& CutSet::operator+= (LogGNodeSet* A)
 *  \brief overloaded += operator, scales the logs, delogs, adds them and takes logs agian the values contained by the
 larger set (the neighborhood set) into the new set.
*/ 
  
CutSet& CutSet::operator*= (GNodeSet* A){
	//  // Authors: L. Radu Totir and Rohan L. Fernando 
	//  // (June, 2003) 
	//  // Contributors:
	//  reset();
	//  do{
	//    double temp = getValue();
	//    double Atemp = A->getValue();
	//    putValue(temp*Atemp);
	//  } while(incr());
	//  normalize();
	//  return *this;
	
	
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	
//	CutSet::iterator itC;
//	cout << "CutSet Members in *= : [";
//	for (itC=begin();itC!=end();itC++){
//		cout << (*itC)->id << " ";
//	} 
//	cout << "]" << endl;
//	cout << "cutset values before * \n";
//	displayValues();	
	
	reset();
	unsigned key = 0;
   	if(dynamic_cast<CutSet*> (A)){
		CutSet* cutsetPtr = reinterpret_cast<CutSet*>(A);
		do{
			//unsigned key = getKey();
			unsigned keyA = cutsetPtr->getKey();
			double Atemp = cutsetPtr->valueVector[keyA];
			valueVector[key] += Atemp;
			key++;
		} while(incr());
	} 
	else {
		do{
			//unsigned key = getKey();
			double Atemp = A->getValue();
			valueVector[key] += Atemp;
			key++;
		}  while(incr());
	}
//	cout << "after *= \n";
//	displayValues();	
	return *this;
}
/*! \fn CutSet& CutSet::operator*= (GNodeSet* A)
 *  \brief overloaded *= operator, multiplies the values contained
 by the various cutsets.
*/   

bool AlleleStateNode::incr(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (based on the incr() of Bernt Guldbrantdsen)
  // (June, 2003) 
  // Contributors: 
  alleleState++;
  if (alleleState==alleleStateVector.size()){
    alleleState=0;
    return 0;
  }
  else {
    return 1;
  }
}
/*! \fn bool AlleleStateNode::incr(void)
 *  \brief the method used to increment the alleleState of an
 AlleleStateNode
*/
void AlleleStateNode::display(void){
	SafeSTLVector<unsigned>::iterator it;
	for (it=alleleStateVector.begin();it!=alleleStateVector.end();it++){
		cout << *it << endl;
	}
	cout << endl;
}
bool AlleleOriginNode::incr(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (based on the incr() of Bernt Guldbrantdsen)
  // (September, 2004) 
  // Contributors: 
  alleleOrigin++;
  if (alleleOrigin==alleleOriginVector.size()){
    alleleOrigin=0;
    return 0;
  }
  else {
    return 1;
  }
}
/*! \fn bool AlleleOriginNode::incr(void)
 *  \brief the method used to increment the alleleOrigin of an
 AlleleOriginNode
*/
  
bool GenotypeNode::incr(){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (based on the incr() of Bernt Guldbrantdsen)
  // (June, 2003) 
  // Contributors:
  genotypeState++;
  if (genotypeState==genotypeVector.size()){
    genotypeState=0;
    return 0;
  }
  else {
    return 1;
  }
}
/*! \fn bool GenotypeNode::incr(void)
 *  \brief the method used to increment the genotypeState of a
 *  GenotypeNode
*/
void GenotypeNode::updateCounter(void) {
		  unsigned mAllele = genotypeVector[genotypeState].maternal;
		  unsigned pAllele = genotypeVector[genotypeState].paternal;
		  unsigned geno = mAllele*GNodeSet::prior->chromosome[0].locus[GNodeSet::currentLocus].nallele() + pAllele;
		  genotypeCount[geno]++;
		  sampleCount++;
		  // this is temp stuff for aviagen {
		  Individual *ind = owner;
		  if (ind->mother()) {
			  ind->m_counter += ind->m_gamete[Individual::currentLocus];
			  ind->p_counter += ind->p_gamete[Individual::currentLocus];
		  }
		  // ============================== }
}
/*! \fn void GenotypeNode::UpdateCounters(void)
 *  \brief update GenotypeNode counters 
*/	  
void GenotypeNode::EEUpdateCounter(void){
	CutSet myGenotProb = calcGNodeProbs();
	myGenotProb.reset();
	do {
		double prob = myGenotProb.getValue();
		unsigned mAllele = genotypeVector[genotypeState].maternal;
		unsigned pAllele = genotypeVector[genotypeState].paternal;
		unsigned geno = mAllele*GNodeSet::prior->chromosome[0].locus[GNodeSet::currentLocus].nallele() + pAllele;
		
		genotypeCount[geno] += prob;
	} while (myGenotProb.incr());
	sampleCount++;
}
/*! \fn bool GenotypeNode::EEUpdateCounters(void)
 *  \brief update GenotypeNode counters by forward and backward Elston-Stewart algorithm 
*/

void GenotypeNode::updateProbs(CutSet& myGenotProb){
	myGenotProb.reset();
	do {
		double prob = myGenotProb.getValue();
		unsigned mAllele = genotypeVector[genotypeState].maternal;
		unsigned pAllele = genotypeVector[genotypeState].paternal;
		unsigned geno = mAllele*GNodeSet::prior->chromosome[0].locus[GNodeSet::currentLocus].nallele() + pAllele;
		
		genotypeCount[geno] += prob;
	} while (myGenotProb.incr());
	sampleCount++;
}

void GenotypeNode::calcGenotypeProbs(void){
	CutSet myGenotProb = calcGNodeProbs();
	myGenotProb.reset();
	do {
		double prob = myGenotProb.getValue();
		unsigned mAllele = genotypeVector[genotypeState].maternal;
		unsigned pAllele = genotypeVector[genotypeState].paternal;
		unsigned geno = mAllele*GNodeSet::prior->chromosome[0].locus[GNodeSet::currentLocus].nallele() + pAllele;
		genotypeCount[geno] = prob;
	} while (myGenotProb.incr());
}
/*! \fn bool GenotypeNode::calcGenotypeProbs(void)
 *  \brief calculate GenotypeNode probabilities by forward and backward Elston-Stewart algorithm 
*/

void GenotypeNode::display(void){
	SafeSTLVector<MaternalPaternalAlleles>::iterator it;
	for (it=genotypeVector.begin();it!=genotypeVector.end();it++){
		cout << it->maternal << ":"<< it->paternal << endl;
	}
	cout << endl;
}



} // end of namespace matvec


