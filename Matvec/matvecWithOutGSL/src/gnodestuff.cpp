/*
 This program is free software; you can redistribute it and/r
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
#include <map>
#include <set>
#include "safe_vectors.h"
#include "stat.h"
//#include "gnodederived.h"
//#include "gnodesetderived.h"
//#include "gnodestuff.h"
#include "population.h"
#include "model.h"
//#include "util.h"

using namespace std; 

namespace matvec {
GNodeList* GNode::gNodeListPtr;
GeneticDist *GNodeSet::prior;
unsigned GNodeSet::currentLocus;
double GNodeList::logTarget = 0.0;
double GNodeList::logProposal = 0.0;
double GNodeList::logOldProposal = 0.0;
Population *GNodeList::popPtr;
int TransitionSet::transmissionType = 1;
CutSet* GNode::magSet = 0;

////////////////// GNodeList methods /////////////////////

void GNodeList::makealleleGNodeSets(char* infilename,char* outfilename){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  ifstream infile(infilename);
  ofstream outfile(outfilename);
  unsigned ind,mom,dad;
  while (infile >> ind >> dad >> mom) {
    outfile << "2   " << (ind-1)*2 +1 << setw(5) << (ind-1)*2 +2 << endl;
    if(!dad) continue;
    outfile << "3   " << (dad-1)*2 +1 << setw(5) << (dad-1)*2 +2 << setw(5) << (ind-1)*2 +1 << endl;
    outfile << "3   " << (mom-1)*2 +1 << setw(5) << (mom-1)*2 +2 << setw(5) << (ind-1)*2 +2 << endl;
  }
}
/*! \fn void GNodeList::makealleleGNodeSets(char* infilename,char*
    outfilename)
 *  \brief creates the AlleleNodeSets from the pedigree
*/
void GNodeList::makegenotGNodeSets(char* infilename,char* outfilename){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  ifstream infile(infilename);
  ofstream outfile(outfilename);
  unsigned ind,mom,dad;
  while (infile >> ind >> dad >> mom) {
    if(!dad) continue;
    outfile << "3   " << dad << setw(5) << mom << setw(5) << ind << endl;
  }
}
/*! \fn void GNodeList::makegenotGNodeSets(char* infilename,char*
    outfilename)
 *  \brief creates the GenotypeNodeSets from the pedigree
*/ 
void GNodeList::fill(unsigned n){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  resize(n);
  for (unsigned i=0;i<n;i++){
    GNode* elm_ptr = new GNode;
    elm_ptr->peelorder = i+1;
    elm_ptr->weight = 1;
    (*this)[i] = elm_ptr;
  }
}
/*! \fn void GNodeList::fill(unsigned n)
 *  \brief method used in determining peeling order for a pedigree
    (used only if make....GNodeSets is not called)
*/
void GNodeList::inputGNodeSets(char* fname){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:  
  ifstream setfile(fname);
  unsigned n,member;
  lastmember = 0;
  while (setfile >> n ) {
    for (unsigned i=0;i<n;i++){
      setfile >> member;
      if (member>lastmember){
	lastmember=member;
      }
    }
  }
  setfile.close();
  ifstream setfile1(fname);
  fill(lastmember);
  while (setfile1 >> n ) {
    GNodeSet* myset = new GNodeSet;
    for (unsigned i=0;i<n;i++){
      setfile1 >> member;
      GNode* elm_ptr = (*this)[member-1];
      myset->insert(elm_ptr);
      elm_ptr->SetofGNsts.insert(myset);
    }
  }
}
/*! \fn void GNodeList::inputGNodeSets(char* fname)
 *  \brief method used in determining peeling order for a pedigree
    (used only if none of the make....GNodeSets methods is called)
*/
void GNodeList::displayGNodeSets(){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	GNodeList::iterator it;
	set< GNodeSet*>::iterator GNodeSetit;
	GNodeSet::iterator setit;
	for (it=begin();it!=end();it++){
		std::cout << (*it)->id << "   ";
		for (GNodeSetit=(*it)->SetofGNsts.begin();GNodeSetit!=(*it)->SetofGNsts.end();GNodeSetit++){
			(*GNodeSetit)->display();
		}
		std::cout << endl;
	}
	std::cout << "=========== " << endl;
}
/*! \fn void GNodeList::displayGNodeSets()
 *  \brief displays the GNodeSets of each Node (pentrance, founder,
    transition sets)
*/
void GNodeList::getPeelOrder(char* fname){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  ofstream ordfile(fname);
  GNodeList::iterator it1,it2,peel_me;
  GNode* temp;
  float mag, min;
  sum =0.0;
  for (it1=begin();it1!=end();it1++){ 
    min = 1e+100;
    for (it2=it1;it2!=end();it2++){
      mag = (*it2)->getCutsetMagnitude();
      if (mag<min){
	peel_me = it2;
	min=mag;
      }
    }    
    sum+=pow(4.0,(double)min);
    ordfile << (*peel_me)->peelorder <<"     "<< min << endl;
    //(*peel_me)->peelord = ++ord;
    (*peel_me)->updateMysets();
    temp = *it1;
    *it1 = *peel_me;
    *peel_me = temp;  
  }
  ordfile << "Network size = " << sum << endl; 
  std::cout << sum << endl; 
}
/*! \fn void GNodeList::getPeelOrder(char* fname)
 *  \brief returns the peeling order for a pedigree provided in the
    file given as an argument
*/
void GNodeList::peelCutAndCompute(unsigned maxCutSetSize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (October, 2003) 
  // Contributors: 
  unsigned sampledCount = 0;
  GNodeList::reverse_iterator it;
  setGNodeSampleFlags(startBlock, stopBlock, sizeBlock);  	
  if (GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].peelOrder.size()==size()){
    sampleGNodePtr = peelAndCutFast();
  }
  else{
    // cout << "We determine the peeling order again (inefficient!!)" << endl;
    sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  // cout << "I am calculating the probability of the previous sample" << endl;
  // cout << "Log Old Proposal before processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logOldProposal << endl;
  //cout << "Log Target before = " << GNodeList::logTarget << endl;
  while (sampleGNodePtr) {
    sampledCount++;
    //displayGNodeSets();
    cout << "GNode sampled after " << sampledCount << " 'th peeling: " << sampleGNodePtr->id << endl;
    sampleGNodePtr->reComputeGNode();
    reinitGNodeList();
    //displayGNodeSets();
    sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  //cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
  //cout << "GNodeList before sampling all GNodes" << endl;
  for (it=rbegin();it!=rend();it++){
    if(!(*it)->sampled){
      (*it)->reComputeGNode();
    }
  }
  checkSample();
  // cout << "Log Old Proposal after processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logOldProposal << endl;
  //calculateTargetProb();
  //cout << "Log Target after = " << GNodeList::logTarget << endl;
  //cout << "Final gNodeList (after peeling, cutting and sampling):" << endl;
  //displayGNodeSets();
  clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelCutAndSample()
 *  \brief sample genotypes after peeling and cutting(if necesary)
*/
void GNodeList::peelCutAndSample(unsigned maxCutSetSize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (October, 2003) 
  // Contributors: 
  unsigned printFlag = popPtr->model->myRSamplerParms.printFlag; 		
  unsigned sampledCount = 0;
  Model *modelPtr = popPtr->model;	
  GNodeList::reverse_iterator it;
  setGNodeSampleFlags(startBlock, stopBlock, sizeBlock);  
  bool cutLoops = GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops;
  if (GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].peelOrder.size()==size() && !cutLoops){
    sampleGNodePtr = peelAndCutFast();
  }
  else{
	  // cout << "We determine the peeling order again (inefficient!!)" << endl;
	  popPtr->SimpleGenotypeElimination(Individual::currentLocus);			
	  if (sizeBlock == size() && GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].qtl_ml=='m') {
		  popPtr->LGGenotypeElimination(Individual::currentLocus); 
	  }
	  popPtr->setAlleleStateVectors(Individual::currentLocus);
	  sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  //cout << "I am actually sampling a new sample" << endl;
  //cout << "Log Proposal before processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logProposal << endl;
  //cout << "Log Target before = " << GNodeList::logTarget << endl;
  while (sampleGNodePtr) {
    sampledCount++;
    //displayGNodeSets();
    if(printFlag>0) cout << "GNode sampled after " << sampledCount << " 'th peeling: " << sampleGNodePtr->id << endl;
    sampleGNodePtr->sampleGNode();
    reinitGNodeList();
    //displayGNodeSets();
	if (modelPtr->myRSamplerParms.samplerType=="genotypic"){
		popPtr->LGGenotypeElimination(Individual::currentLocus);
		popPtr->setAlleleStateVectors(Individual::currentLocus);
	}
    sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  if(printFlag>0) cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
  //cout << "GNodeList before sampling all GNodes" << endl;
  for (it=rbegin();it!=rend();it++){
    if(!(*it)->sampled){
      (*it)->reverseSampleGNode();
    }
  }
  checkSample();
  //cout << "Log Proposal after processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logProposal << endl;
  calculateTargetProb();
  //cout << "Log Target after = " << GNodeList::logTarget << endl;
  //cout << "Final gNodeList (after peeling, cutting and sampling):" << endl;
  //displayGNodeSets();
  clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelCutAndSample()
 *  \brief sample genotypes after peeling and cutting(if necesary)
*/
void GNodeList::peelOrderCutAndSample(unsigned maxCutSetSize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (February, 2004) 
  // Contributors: 
  unsigned printFlag = popPtr->model->myRSamplerParms.printFlag; 
  Model *modelPtr = popPtr->model;
  unsigned sampledCount = 0;
  GNodeList::reverse_iterator it;
  setGNodeSampleFlags(startBlock, stopBlock, sizeBlock); 
  sampleGNodePtr = peelAndCut(maxCutSetSize);
  GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops = false;
  //cout << "Log Proposal before = " << GNodeList::logProposal << endl;
  //cout << "Log Target before = " << GNodeList::logTarget << endl;
  while (sampleGNodePtr) {
	GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops = true;  
    sampledCount++;
    //displayGNodeSets();
    if(printFlag>0) cout << "GNode sampled after " << sampledCount << " 'th peeling: " << sampleGNodePtr->id << endl;
    sampleGNodePtr->sampleGNode();
    reinitGNodeList();
	if (modelPtr->myRSamplerParms.samplerType=="genotypic"){
		popPtr->LGGenotypeElimination(Individual::currentLocus);
		popPtr->setAlleleStateVectors(Individual::currentLocus);
	}
    //displayGNodeSets();
    sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  if(printFlag>0) cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
  //cout << "GNodeList before sampling all GNodes" << endl;
  for (it=rbegin();it!=rend();it++){
    if(!(*it)->sampled){
      (*it)->reverseSampleGNode();
    }
  }
  checkSample();
  //cout << "Log Proposal after = " << GNodeList::logProposal << endl;
  calculateTargetProb();
  //cout << "Log Target after = " << GNodeList::logTarget << endl;
  //cout << "Final gNodeList (after peeling, cutting and sampling):" << endl;
  //displayGNodeSets();
  clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelOrderCutAndSample()
 *  \brief sample genotypes after determining peeling order, peeling
    and cutting(if necesary)
*/
GNode* GNodeList::peelAndCut(unsigned maxCutSetSize){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2003) 
	// Contributors: Joseph Abraham
	unsigned printFlag;
	string samplerUsed; 	
	if (popPtr->model){
		printFlag = popPtr->model->myRSamplerParms.printFlag;
		samplerUsed = popPtr->model->myRSamplerParms.samplerUsed;
	}
	else {
		printFlag = MIM::parmMap.getUnsignedValue("printFlag");
		samplerUsed = MIM::parmMap["GNodeSampler"];
	}	 		
	if (!distanceMatDone) makeDistanceMatrix(); 
	cutGNodesVector.clear();
	loopSetGNodes.clear();  
	searchStart = begin();
	resetConnectFlags();
	sampleGNodePtr=NULL;
	maxDist = 0;
	GNodeList::iterator it1,it2,peel_me;
	GNode* temp;
	double mag, min;
	sum =0.0;
	unsigned count=0;
	unsigned loopCount=0;
	it1=begin();
	while (it1!=end()){
		min = 1e+100;
		if (samplerUsed != "LG"){
		for (it2=it1;it2!=end();it2++){
			mag = (*it2)->getCutsetMagnitude();
			if (mag<min){
				peel_me = it2;
				min=mag;
				//if (min<17) break;
			}
		} 
		}
		if (samplerUsed == "LG"){
			peel_me = it1;
			min = (*peel_me)->getCutsetMagnitude();
		}
		if(min<=maxCutSetSize){
			sum+=min;
			(*peel_me)->peelorder = ++count;
			(*peel_me)->updateMysets();
			if(printFlag>0) cout << "Peeling node " << count << " with  id  " <<  (*peel_me)->id  << " & cost  " << min << endl;
			(*peel_me)->peel();
			if(howToSample=="joint"){
				GNodeSet::prior->chrom()[0].peelOrder[count-1] = (*peel_me)->id;
			}
			else{
				GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].peelOrder[count-1] = (*peel_me)->id;
			}
			if (samplerUsed != "LG"){
			temp = *it1;
			*it1 = *peel_me;
			*peel_me = temp;
			}
			it1++;
		}
		else {
			if (samplerUsed == "P" || samplerUsed == "LG") throw exception(" Not possible to peel without larger maximum cutset size");
			if(printFlag>0) cout << "Peeled " << count << " gNodes, now I need to cut " << endl;
			searchStart = it1;
			sort(searchStart,end(),compareGNodesWeight());
			findLoopAndCut();
			loopCount++;
		}
		if(printFlag>2) {
			displayGNodeSets();
		}
	}
	if (cutGNodesVector.size()>0){
		sampleGNodePtr = findSampleGNode();
		cout << "The GNode selected for sampling is  " << sampleGNodePtr->id << endl;
	}		
	if(printFlag>0) cout << loopCount  << " loops were cut" << endl;
	cout << "The total cost was  " << sum << endl;
//	cout << " The likelihood is  " << CutSet::NormalizingConstant << endl;
	return sampleGNodePtr;
}
/*! \fn void GNodeList::peelAndCut()
 *  \brief determines peeling order, peels and cuts (if it needs to) 
*/
bool GNodeList::peelNoCut(unsigned maxCutSetSize){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2003) 
	// Contributors: 
	unsigned printFlag, count=0;
	if (popPtr->model){
		printFlag = popPtr->model->myRSamplerParms.printFlag;
	}
	else {
		printFlag = MIM::parmMap.getUnsignedValue("printFlag");
	}	 		
	GNodeList::iterator it1,it2,peel_me;
	GNode* temp;
	float mag, min;
	it1=begin();
	while (it1!=end()){
		min = 1e+100;
		for (it2=it1;it2!=end();it2++){
			mag = (*it2)->getCutsetMagnitude();
			if (mag<min){
				peel_me = it2;
				min=mag;
				//if (min<17) break;
			}
		} 
		if(min<=maxCutSetSize){
			count++;
			(*peel_me)->updateMysets();
			if(printFlag>1) cout << "Peeling node " << count<<"("<<(*peel_me)->id<<")"<< " cutset magnitude " << min << endl;
			(*peel_me)->peel();
			temp = *it1;
			*it1 = *peel_me;
			*peel_me = temp;
			it1++;
		}
		else {
		    cerr << "cutset magnitude = " << mag << endl;
			cerr << "maxCutSetSize    = " << maxCutSetSize << endl;
			return false;
		}
		if(printFlag>2) {
			displayGNodeSets();
		}
	}
	return true;
}
/*! \fn void GNodeList::peelAndCut()
 *  \brief determines peeling order, peels and cuts (if it needs to) 
*/ 

GNode* GNodeList::findSampleGNode(void){
	// Authors: Joseph Abraham and Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	unsigned maxDistance = 0;
	
	GNode* sampleGuy = 0;
	set<GNode*>::iterator it;
	SafeSTLVector<GNode*>::iterator cutListIter;
	for(it=loopSetGNodes.begin();it!=loopSetGNodes.end();it++){
		if((*it)->sampled) continue;
		unsigned rowId = (*it)->id;
		unsigned sumDist = 0;
		for(cutListIter=cutGNodesVector.begin();cutListIter!=cutGNodesVector.end();cutListIter++){
			unsigned colId = (*cutListIter)->id;
			sumDist += distanceMat(rowId,colId);
		}
		if (sumDist>maxDistance){
			maxDistance = sumDist;
			sampleGuy = (*it);
		}
	}
	if (sampleGuy){
		return sampleGuy;
	}
	else {
		throw exception("GNodeList::findSampleGNode: Error- possible bug in code\n");
	}
}
/*! \fn GNode* GNodeList::findSampleGNode()
*  \brief returns pointer to GNode with largest sum of distances to the "cuts"
*/
void GNodeList::makeDistanceMatrix(void){
	// Authors: Joseph Abraham and Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	unsigned numGnodes = size();
	distanceMat.resize(numGnodes,numGnodes);
	GNodeList::iterator it;
	for (it=begin();it!=end();it++){
		calcDistancefrom(*it);
	}
	distanceMatDone = true;
	//cout << distanceMat;
}
/*! \fn void GNodeList::makeDistanceMatrix
*  \brief fills in distanceMatrix for all pairs of nodes
*/

void GNodeList::calcDistancefrom(GNode* refGNode){
    // Authors: Joseph Abraham and Rohan L. Fernando
    // (October, 2004)
    // Contributors:

    SafeSTLVector<unsigned> dist;
    unsigned long UHUGE = 2*ULONG_MAX;
    GNode *gnodePtr;
    GNodeList::iterator GNodeListiter;
    set<GNodeSet*>::iterator GNodeSetiter;
    set<GNode*>::iterator GNodeiter;
    list<GNode*>::iterator GreyListiter;
    list<GNode*> GreyList;

    GreyList.clear();
    unsigned rowId = refGNode->id;
    for(GNodeListiter = begin(); GNodeListiter != end(); GNodeListiter++){
        (*GNodeListiter)->connectFlag =  ULONG_MAX; // every one gets white color
    }
    refGNode->connectFlag = 0; // this is colored gray
    GreyList.push_back(refGNode);
    distanceMat(rowId,rowId) = 0;
    while (GreyList.size() > 0){
        gnodePtr = GreyList.back();
        GreyList.pop_back();
        unsigned uId = gnodePtr->id;
        for(GNodeSetiter= gnodePtr->SetofGNsts.begin();GNodeSetiter !=gnodePtr->SetofGNsts.end();GNodeSetiter++){
            for(GNodeiter = (*GNodeSetiter)->begin(); GNodeiter!= (*GNodeSetiter)->end(); GNodeiter++){
                if ( (*GNodeiter)->connectFlag == ULONG_MAX) { // if color white
                    (*GNodeiter)->connectFlag = 0;             // make it gray
                    unsigned colId = (*GNodeiter)->id;
                    distanceMat(rowId,colId) = distanceMat(rowId,uId) + 1;
                    GreyList.push_front(*GNodeiter);
                }
            }
        }
        gnodePtr->connectFlag = UHUGE;
    }
}
/*! \fn void GNodeList::calcDistancefrom((GNode* refGNode)
*  \brief calculate distance from argument to all the nodes using the breadth first search algorithm (got it from  Boost)
*/
GNodeList GNodeList::findNeighbors(GNode* refGNode, unsigned largestdist){
    // Authors: Joseph Abraham and Rohan L. Fernando
    // (October, 2005)
    // Contributors:

    static SafeSTLVector<unsigned long> distvec;
	SafeSTLVector<unsigned> indexVec;
	GNodeList::iterator GNodeListiter;
    if (distvec.size() == 0) {
		distvec.resize(size(),0);
		for(GNodeListiter = begin(); GNodeListiter != end(); GNodeListiter++){
			(*GNodeListiter)->connectFlag =  ULONG_MAX; // every one gets white color
		}
	}
    GNode *gnodePtr;
    GNodeList Neighbors;
    set<GNodeSet*>::iterator GNodeSetiter;
    set<GNode*>::iterator GNodeiter;
    list<GNode*>::iterator GreyListiter;
    list<GNode*> GreyList;

    GreyList.clear();
    unsigned rowId = refGNode->id -1;
    distvec[rowId] = 0;
	indexVec.push_back(rowId);
//    for(GNodeListiter = begin(); GNodeListiter != end(); GNodeListiter++){
//        (*GNodeListiter)->connectFlag =  ULONG_MAX; // every one gets white color
//    }
    refGNode->connectFlag = 0; // this is colored gray
    GreyList.push_back(refGNode);
	Neighbors.push_back(refGNode);
    while (GreyList.size() > 0 ){
        gnodePtr = GreyList.back();
        GreyList.pop_back();
        for(GNodeSetiter= gnodePtr->SetofGNsts.begin();GNodeSetiter !=gnodePtr->SetofGNsts.end();GNodeSetiter++){
            for(GNodeiter = (*GNodeSetiter)->begin(); GNodeiter!= (*GNodeSetiter)->end(); GNodeiter++){
                if ( (*GNodeiter)->connectFlag == ULONG_MAX) { // if color white
                    (*GNodeiter)->connectFlag = 0;             // make it gray
                    distvec[(*GNodeiter)->id -1] = distvec[gnodePtr->id -1] + 1;
					indexVec.push_back((*GNodeiter)->id -1);
                    if (distvec[(*GNodeiter)->id -1] <= largestdist) {
                        Neighbors.push_back(*GNodeiter);
                        GreyList.push_front(*GNodeiter);
                    }
                }
            }
        }
    }
	// now I will set all back to white
	unsigned index;
	for (unsigned i=0;i<indexVec.size();i++){
		index = indexVec[i];
		distvec[index] = ULONG_MAX;
		(*this)[index]->connectFlag = ULONG_MAX;
	}
	indexVec.clear();
    return Neighbors;
}


GNode* GNodeList::peelAndCutFast(){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (February, 2004) 
  // Contributors: 
  // cout << "Before sorting and peeling the GNodeList looks like this: " << endl;
  // displayGNodeSets();
  //sort(begin(),end(),compareGNodesId());
  // cout << "After sorting but before peeling the GNodeList looks like this: " << endl;
  // displayGNodeSets();
  //cout << "Peeling locus: " << Individual::currentLocus << endl;
  searchStart = begin();
  resetConnectFlags();
  sampleGNodePtr=NULL;
  maxDist = 0;
  GNode* peel_me;
  unsigned idOfPeelGuy=0;
  for (int i=0; i<size();i++){
    if(howToSample=="joint"){
      idOfPeelGuy = GNodeSet::prior->chrom()[0].peelOrder[i]-1;
    }
    else{
      idOfPeelGuy = GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].peelOrder[i]-1;
    }
    peel_me = (*this)[idOfPeelGuy];
    // if(peel_me->id==76){
    // cout << "Ready to peel GNode: " << peel_me->id << endl;
    // }  
    peel_me->updateMysets();
    peel_me->peel();
  }
  sort(begin(),end(),compareGNodesPeelId());
  // cout << "After peeling the GNodeList looks like this: " << endl;
  // displayGNodeSets();
//   }
//   else {
//     cout << "Peeled " << count << " gNodes, now I need to cut" << endl;
//     searchStart = it1;
//     sampleGNodePtr = findLoopAndCut();
//     loopCount++;
//   }
  //displayGNodeSets();
  //cout << "Number of loops that were cut =" << loopCount << endl;; 
  return sampleGNodePtr;
  //std::cout << "Network size = " << sum << endl; 
}
/*! \fn void GNodeList::peelAndCut()
 *  \brief peels and cut (if it needs to) with known peeling order
*/
void GNodeList::clearGNodeListForNextLocus(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (September, 2003) 
  // Contributors:
  set< CutSet*>::iterator GNodeSetit;
  for (GNodeSetit=completeSetofCutSets.begin();GNodeSetit!=completeSetofCutSets.end();GNodeSetit++){ 
    delete *GNodeSetit;
  }
  completeSetofCutSets.clear();
  GNodeList::iterator it;
  for (it=begin();it!=end();it++){ 
    (*it)->SetofGNsts.clear();
  }
}

void GNodeList::releaseGNsts(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (February, 2004) 
  // Contributors:     
  set< GNodeSet*>::iterator GNodeSetit;
  for (GNodeSetit=completeSetofGNsts.begin();GNodeSetit!=completeSetofGNsts.end();GNodeSetit++){ 
    delete *GNodeSetit;
  }
  completeSetofGNsts.clear();
}
/*! \fn void GNodeList::releaseGNsts(void)
 *  \brief delete GNodeSet (penetrance, founder, transition) pointers
 *  and clear the completeSetofGNsts
*/
void GNodeList::reinitGNodeList(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (September, 2003) 
  // Contributors:
  releaseCutSets();
  resetBackSets();
  resetGNodeList();
}
/*! \fn void GNodeList::reinitGNodeList(void)
 *  \brief reinitialize the GNodeList for a new round of peeling 
*/
void GNodeList::releaseCutSets(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (September, 2003) 
  // Contributors:     
  set< CutSet*>::iterator GNodeSetit;
  for (GNodeSetit=completeSetofCutSets.begin();GNodeSetit!=completeSetofCutSets.end();GNodeSetit++){ 
    delete *GNodeSetit;
  }
  completeSetofCutSets.clear();
}
/*! \fn void GNodeList::releaseCutSets(void)
 *  \brief release CutSets before reseting the initial GNodeList
 *  for a new round of peeling 
*/
void GNodeList::resetGNodeList(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
  GNodeList::iterator it;
  for (it=begin();it!=end();it++){ 
    (*it)->SetofGNsts.clear();
  }
  set< GNodeSet*>::iterator GNodeSetit;
  GNodeSet::iterator setit;
  for (GNodeSetit=completeSetofGNsts.begin();GNodeSetit!=completeSetofGNsts.end();GNodeSetit++){  
    for (setit=(*GNodeSetit)->begin();setit!=(*GNodeSetit)->end();setit++){
      (*setit)->SetofGNsts.insert(*GNodeSetit);
    }
  }
}
/*! \fn void GNodeList::resetGNodeList(void)
 *  \brief reset the initial GNodeList for a new round of peeling 
*/
void GNodeList::setGNodeSampleFlags(unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (October, 2003) 
  // Contributors: 
  GNodeList::iterator it;
  if ((stopBlock-startBlock)==size()){
    for (it=begin();it!=end();it++){ 
		(*it)->sampled = false;
    }
    return;
  }
  else {
    for (it=begin();it!=end();it++){
      if(it>=begin()+startBlock && it<begin()+stopBlock){
		  (*it)->sampled = false;
      }
    }
  }
}
/*! \fn void GNodeList::setGNodeSampleFlags()
 *  \brief resets the sampled flags to false
*/      
void GNodeList::findLoopAndCut(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors: 
  unsigned lowFlag;
  unsigned highFlag;
  GNodeList::iterator it;
  set< GNodeSet*>::iterator GNodeSetit, cutGNodeSet;
  GNodeSet::iterator setit;
  resetConnectFlags();
  for (it=searchStart;it!=end();it++){ 
	  //std::cout << (*it)->id << "   ";
	  for (GNodeSetit=(*it)->SetofGNsts.begin();GNodeSetit!=(*it)->SetofGNsts.end();GNodeSetit++){
		  if((*it)->connectFlag < (*GNodeSetit)->connectFlag){
			  lowFlag = (*it)->connectFlag;
			  highFlag= (*GNodeSetit)->connectFlag;
			  (*GNodeSetit)->connectFlag = (*it)->connectFlag;
			  propagateFlags(lowFlag,highFlag);
			  //cout << (*GNodeSetit)->connectFlag;
		  }
		  else if((*it)->connectFlag > (*GNodeSetit)->connectFlag){
			  lowFlag = (*GNodeSetit)->connectFlag;
			  highFlag= (*it)->connectFlag;
			  (*it)->connectFlag = (*GNodeSetit)->connectFlag;
			  propagateFlags(lowFlag,highFlag);
		  }
		  else{
			  // we have found a subset of nodes containing at least one loop 
			  // and possibly many "dangling branches"
			  GNodeList::iterator loopCloser = it;
			  GNodeList loopList = isolatePureLoop(loopCloser);
			  cutGNode =  loopList.getMinWeightGNode();
			  cutGNodesVector.push_back(*cutGNode);
			  cutGNodeSet = loopList.getCutGNodeSet(cutGNode);
			  // cut the loop
			  cutLoop(cutGNode,cutGNodeSet);
			  return;
		  }
	  }
  }
}
/*! \fn void GNodeList::findLoopAndCut()
 *  \brief finds a loop in the unpeeled part of GNodeList and cuts it 
*/
set< GNodeSet*>::iterator GNodeList::getCutGNodeSet(GNodeList::iterator cutGNode){
	// Authors: Joseph Abraham and Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	GNodeSet loopGNodesSet;
	GNodeList intersectionVector;
	GNodeList::iterator it;
	for (it=begin();it!=end();it++){
		if(it!=cutGNode) loopGNodesSet.insert(*it);	
	}
	set< GNodeSet*>::iterator GNodeSetit;
	for (GNodeSetit=(*cutGNode)->SetofGNsts.begin();GNodeSetit!=(*cutGNode)->SetofGNsts.end();GNodeSetit++){
		intersectionVector.clear();
		set_intersection(loopGNodesSet.begin(),  loopGNodesSet.end(),
						 (*GNodeSetit)->begin(), (*GNodeSetit)->end(),
						 inserter(intersectionVector,intersectionVector.begin())  );
		if (intersectionVector.size()>0) return GNodeSetit;
	}
	throw exception("GNodeList::getCutGNodeSet: possible error in program");
}
/*! \fn set< GNodeSet*>::iterator GNodeList::getCutGNodeSet(GNodeList::iterator cutGNode)
*  \brief This method is called for a list of nodes in a loop. 
          It returns a GNodeSet that is in the loop 
*/
GNodeList::iterator GNodeList::getMinWeightGNode(void){
	// Authors: Joseph Abraham and Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	GNodeList::iterator it, minIt;
	unsigned long minWeight = ULONG_LONG_MAX;
	for(it=begin();it!=end();it++){
		unsigned wt = (*it)->sampled ? 1 : (*it)->getWeight();
		if(wt < minWeight){
			minWeight = wt;
			minIt = it;
		}
	}
	return minIt;
}
/*! \fn void GNodeList::getMinWeightGNode
*  \brief returns GNode with smallest weight in this List

*/void GNodeList::propagateFlags(unsigned lowFlag, unsigned highFlag){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors: 
  if(highFlag==ULONG_MAX){
    return;
  }
	else{
		GNodeList::iterator it;
		set< GNodeSet*>::iterator GNodeSetit;  
		for (it=searchStart;it!=end();it++){ 
			if((*it)->connectFlag==highFlag){
				(*it)->connectFlag = lowFlag;
			}
			for (GNodeSetit=(*it)->SetofGNsts.begin();GNodeSetit!=(*it)->SetofGNsts.end();GNodeSetit++){ 
				if((*GNodeSetit)->connectFlag==highFlag){
					(*GNodeSetit)->connectFlag=lowFlag;
				}
			}
		}
	}
}
/*! \fn void GNodeList::propagateFlags(void)
 *  \brief reassign the flag of a GNode or GNodeSet after the
    comparison with an adjacent GNode or GNodeSet
*/
void GNodeList::resetConnectFlags(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors: 
  // resseting the flags for each loop
  GNodeList::iterator it;
  set< GNodeSet*>::iterator GNodeSetit;
  for (it=searchStart;it!=end();it++){
    (*it)->connectFlag = (*it)->id;
    for (GNodeSetit=(*it)->SetofGNsts.begin();GNodeSetit!=(*it)->SetofGNsts.end();GNodeSetit++){
      (*GNodeSetit)->connectFlag = ULONG_MAX;
    }
  }
}
/*! \fn void GNodeList::resetConnectFlags(void)
 *  \brief reset the flag of a GNode and GNodeSet respectively before
    starting to look for a place to cut
*/
GNodeList GNodeList::isolatePureLoop(GNodeList::iterator loopCloser){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2003) 
	// Contributors: 
	GNodeList LoopList;
	GNodeList::iterator it;
	set< GNodeSet*>::iterator GNodeSetit; 
	GNodeSet::iterator setit;
	//cout << "Loop members: " << endl; 
	unsigned numberEliminated;
    do {
		numberEliminated = 0;
		for (it=searchStart;it!=loopCloser;it++){
			if((*loopCloser)->connectFlag==(*it)->connectFlag){
				unsigned GNodeSetcount = 0;
				for (GNodeSetit=(*it)->SetofGNsts.begin();GNodeSetit!=(*it)->SetofGNsts.end();GNodeSetit++){
					unsigned GNodecount = 0;
					for (setit=(*GNodeSetit)->begin();setit!=(*GNodeSetit)->end();setit++){
						if((*setit)->id!=(*it)->id && (*setit)->connectFlag==(*loopCloser)->connectFlag){
							GNodecount++;
						}
					}
					if(GNodecount!=0) GNodeSetcount++;
				}
				if(GNodeSetcount<=1){
					(*it)->connectFlag = ULONG_MAX;
					numberEliminated++;
				}
			}
		}
	} while(numberEliminated>0);
	for (it=searchStart;it!=end();it++){
		if((*loopCloser)->connectFlag==(*it)->connectFlag){
			LoopList.push_back(*it);
		}
	}
	for (it=LoopList.begin();it!=LoopList.end();it++){
		loopSetGNodes.insert(*it);
	}
	//LoopList.displayGNodeSets();
	return LoopList;
}
/*! \fn GNodeList GNodeList::isolatePureLoop (void)
*   \brief find the members of a loop by removing all terminal GNodes;
     returns a GNodeList containing loop members
   
*/

void  GNodeList::cutLoop(GNodeList::iterator cutGNode,set< GNodeSet*>::iterator cutGNodeSet){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2003) 
	// Contributors:  
	unsigned printFlag;
	if (popPtr->model){
		printFlag = popPtr->model->myRSamplerParms.printFlag;
	}
	else {
		printFlag = MIM::parmMap.getUnsignedValue("printFlag");
	}		
	if(printFlag>1){
		cout << "cutGNode: "<< (*cutGNode)->id << " cutGNodeSet: ";
		(*cutGNodeSet)->display();
		cout << endl;
	}
	GNodeSet::iterator setit;
	CutSet* newSet = new CutSet;
	completeSetofCutSets.insert(newSet);
	for (setit=(*cutGNodeSet)->begin();setit!=(*cutGNodeSet)->end();setit++){
		if((*setit)!=(*cutGNode)){
			newSet->insert(*setit);
		}
	}
	newSet->setKeyMultiplicationCode();
	newSet->connectFlag = (*cutGNode)->connectFlag;
	newSet->numberOfCuts++;
	(*cutGNode)->numberOfCuts++;
	*newSet+=(*cutGNodeSet);
	if(dynamic_cast<CutSet*>(*cutGNodeSet)){
		((CutSet*) *cutGNodeSet)->replaceMeWith(*newSet);
		(*cutGNode)->SetofGNsts.erase(*cutGNodeSet);
		delete newSet;
		completeSetofCutSets.erase(newSet);
	}
	else {
		for (setit=(*cutGNodeSet)->begin();setit!=(*cutGNodeSet)->end();setit++){
			(*setit)->SetofGNsts.erase(*cutGNodeSet);
			if((*setit)!=(*cutGNode)){
				(*setit)->SetofGNsts.insert(newSet);
			}
		}
	}
}
/*! \fn void GNodeList::cutLoop()
*  \brief cut the loop 
*/
void GNodeList::checkSample(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2003) 
	// Contributors:
	set< GNodeSet*>::iterator GNodeSetit;
	for (GNodeSetit=completeSetofGNsts.begin();GNodeSetit!=completeSetofGNsts.end();GNodeSetit++){ 
		if(dynamic_cast<TransitionSet*>(*GNodeSetit)){
			if((*GNodeSetit)->getValue()==-9999.0){
				throw matvec::InvalidSample();
			}
		}
		else if(dynamic_cast<TransmissionSet*>(*GNodeSetit)){
			if((*GNodeSetit)->getValue()==-9999.0){
				throw matvec::InvalidSample();
			}
		}
	}
}

void GNodeList::calculateTargetProb(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (September, 2003) 
  // Contributors:
  set< GNodeSet*>::iterator GNodeSetit;
  for (GNodeSetit=completeSetofGNsts.begin();GNodeSetit!=completeSetofGNsts.end();GNodeSetit++){ 
    if(dynamic_cast<TransitionSet*>(*GNodeSetit)){
      //std::cout << "Target = " << (*GNodeSetit)->getTargetValue() << endl;
      logTarget += (*GNodeSetit)->getTargetValue();
    }
    else {
      //std::cout << "Target = " << (*GNodeSetit)->getValue() << endl;
      logTarget += (*GNodeSetit)->getValue();
      //std::cout << "Intermediate logTarget = " << logTarget << endl;
    }
    //std::cout << "Cumulative logTarget = " << logTarget << endl;
  }
}
  
////////////////// GNodeSet methods /////////////////////

GNodeSet::GNodeSet(){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:  
  connectFlag = ULONG_MAX;
  numberOfCuts= 0;
}
/*! \fn GNodeSet::GNodeSet()
 *  \brief GNodeSet constructor
*/

bool GNodeSet::incr(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (based on the incr() of Bernt Guldbrantdsen)
  // (June, 2003) 
  // Contributors:
  GNodeSet::iterator iter;
  for(iter=begin();iter!=end();iter++) {
    if((*iter)->sampled) continue;
    if((*iter)->incr()) return 1;
  }
  return 0;
}
  
void GNodeSet::reset(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  GNodeSet::iterator iter;
  if(empty())
    return;
  for(iter =  begin();iter != end();iter++) {
    if((*iter)->sampled) continue;
    (*iter)->reset(0);
  }
}
/*! \fn void GNodeSet::reset(void)
 *  \brief sets all GNodes to 0
*/
void GNodeSet::attachMeToMyGnodes(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (September, 2004) 
  // Contributors:
  GNodeSet::iterator iter;
  for(iter =  begin();iter != end();iter++) {
    (*iter)->SetofGNsts.insert(this);
  }
}
/*! \fn void GNodeSet::attachMeToMyGnodes(void)
 *  \brief attaches the current GNodeSet to all its GNode members
*/
void GNodeSet::display(void){
	GNodeSet::iterator it;
	cout << "(";
	for (it=begin();it!=end();it++){
		cout << (*it)->id <<" ";
	}
	cout << ") ";
}

////////////////// GNode methods /////////////////////

bool GNode::isMyNeighbor(GNode* refGNode){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
  set< GNodeSet*>::iterator myGNodeSetit;
  GNodeSet::iterator setit; 
  for (myGNodeSetit=SetofGNsts.begin();myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
    for (setit=(*myGNodeSetit)->begin();setit!=(*myGNodeSetit)->end();setit++){ 
      if((*setit)->id==refGNode->id){
	return true;
      }
    }
  }
  return false;
}
/*! \fn bool GNode::isMyNeighbor(GNode* refGNode)
 *  \brief determine the neighbours of the GNode where we will cut;
 *  this is done, to determine the GNode that is farthest away from 
 *  the cut
*/
CutSet* GNode::makeNeighSet(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  set< GNodeSet*>::iterator myGNodeSetit;
  GNodeSet::iterator setit;
  CutSet* neighSet = new CutSet;
  for (myGNodeSetit=SetofGNsts.begin();myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
    for (setit=(*myGNodeSetit)->begin();setit!=(*myGNodeSetit)->end();setit++){
      neighSet->insert(*setit);
    }
  }
  // cout << "Neigh size= " << neighSet->size() << endl;
  return neighSet;
}
/*!/*! \fn CutSet* GNode::makeNeighSet(void)
 *  \brief creates the neighborhood cutset (the resulting cutset
    includes the individual that is going to be peeled as well)
*/  
double GNode::getCutsetMagnitude(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	
	if (magSet == 0){
		magSet = new CutSet;
	}
	else {
		magSet->clear();
	}
	
	set< GNodeSet*>::iterator myGNodeSetit;
	GNodeSet::iterator setit;
	for (myGNodeSetit=SetofGNsts.begin();myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
		for (setit=(*myGNodeSetit)->begin();setit!=(*myGNodeSetit)->end();setit++){
			magSet->insert(*setit);
		}
	}
	magSet->erase(this);
	double mag=1.0;
	for (setit=magSet->begin();setit!=magSet->end();setit++){
		mag *= (*setit)->sampled ? 1.0 : (*setit)->getWeight();
	}
	return mag;
}
/*! \fn double GNode::getCutsetMagnitude(void)
 *  \brief returns the magnitude of the resulting cutset if the
    individual to be peeled is removed from its neighborhood cutset
*/

void GNode::peel() {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors:
  CutSet* tempSet = makeNeighSet();
  tempSet->setKeyMultiplicationCode();
  generatedSet->setKeyMultiplicationCode();
  gNodeListPtr->completeSetofCutSets.insert(tempSet);
  // cout << "Resulting CutSet size = " << tempSet->size() - 1 << endl;
  set< GNodeSet*>::iterator myGNodeSetit;
  myGNodeSetit = SetofGNsts.begin();
  *tempSet = *myGNodeSetit;
  for (myGNodeSetit++;myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
    *tempSet *= *myGNodeSetit;
  }
  *generatedSet+= tempSet;
//   if(generatedSet->size()==0){
//     cout << "Likelihood = ";
//     generatedSet->displayValues();
//   }
  delete tempSet;
  gNodeListPtr->completeSetofCutSets.erase(tempSet);
}
/*! \fn void GNode::peel(void)
 *  \brief the GNode peeling method (peeling the GNode in question)
*/
void GNode::reComputeGNode(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
  //if(GNodeSet::currentLocus==1){
  //cout << "Recompute old sample probability for " << id << " = ";
  //}
  CutSet* tempSet = makeNeighSet();
  tempSet->setKeyMultiplicationCode();
  // next line is used to prevent a memory leak
  gNodeListPtr->completeSetofCutSets.insert(tempSet);
  GNodeSet::iterator setit;
  CutSet::iterator cutSetit;
  CutSet myGNodeProbs;
  myGNodeProbs.insert(this);
  myGNodeProbs.setKeyMultiplicationCode();
  set< GNodeSet*>::iterator myGNodeSetit;
  myGNodeSetit = SetofGNsts.begin();
  *tempSet = *myGNodeSetit;
  for (myGNodeSetit++;myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
    *tempSet *= *myGNodeSetit;
  }
  myGNodeProbs += tempSet;
  delete tempSet;
  gNodeListPtr->completeSetofCutSets.erase(tempSet);
  myGNodeProbs.normalize();
  myGNodeProbs.reset();
  double oldProb = myGNodeProbs.getOldValue();
  if(oldProb==0) cout << "Zero prob when recalculating value for: "<< id << endl; 
  GNodeList::logOldProposal += std::log(oldProb);
  int oldState = getOldState();
  reset(oldState);
  sampled = true; 
  //if(GNodeSet::currentLocus==1){
  //std::cout << getOldState() << "  " << getState() << " (" << getmState()+1 << "|" << getpState()+1 << ")" << endl;
  //}
//   std::cout << "oldProb = " << oldProb << endl;
//   cout << "GNodeList::logOldProposal = " << GNodeList::logOldProposal << endl;
}

void GNode::reverseSampleGNode(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
	
	try{
		CutSet* tempSet = makeNeighSet();
		tempSet->setKeyMultiplicationCode();
		// next line is used to prevent a memory leak
		gNodeListPtr->completeSetofCutSets.insert(tempSet);
		GNodeSet::iterator setit;
		CutSet::iterator cutSetit;
		if(myGNodeProbs->size()==0){
			myGNodeProbs->insert(this);
			myGNodeProbs->setKeyMultiplicationCode();
		}
		set< GNodeSet*>::iterator myGNodeSetit;
		myGNodeSetit = SetofGNsts.begin();
		*tempSet = *myGNodeSetit;
		for (myGNodeSetit++;myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
			*tempSet *= *myGNodeSetit;
		}
		*myGNodeProbs += tempSet;
		delete tempSet;
		gNodeListPtr->completeSetofCutSets.erase(tempSet);
		myGNodeProbs->normalize();
		double u = ranf();
		double sum = 0.0;
		myGNodeProbs->reset();
		unsigned sampledState=0;
		do{
			double tempProb = myGNodeProbs->getValue();
			sum+= tempProb;
			if (u<sum){
				reset(sampledState);
				sampled = true;
				GNodeList::logProposal += std::log(tempProb);
				return;
			}
			sampledState++;
		} while(incr());
	}
	catch (matvec::exception &ex) {    
		cout << ex.what() << "\n"; 
		cout << "in reverse sample\n";
	}	
	catch (...){
		cout << "in reverse sample\n";
	}
}

/*! \fn void GNode::sampleGNode(void)
 *  \brief the GNode sampling method (sample the GNode in question)
*/
void GNode::sampleGNode(void) {
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2003) 
	// Contributors:
	// cout << "Sampled genotype for " << id << " = " ;
	calcMyBCutSet();
	CutSet* tempSet = makeNeighSet();
	tempSet->setKeyMultiplicationCode();
	gNodeListPtr->completeSetofCutSets.insert(tempSet);
	GNodeSet::iterator setit;
	CutSet::iterator cutSetit;
	if(myGNodeProbs->size() == 0){
		myGNodeProbs->insert(this);
		myGNodeProbs->setKeyMultiplicationCode();
	}
	set< GNodeSet*>::iterator myGNodeSetit;
	myGNodeSetit = SetofGNsts.begin();
	*tempSet = *myGNodeSetit;
	for (myGNodeSetit++;myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
		*tempSet *= *myGNodeSetit;
	}
	*tempSet *= backSet;
	*myGNodeProbs += tempSet;
	delete tempSet;
	gNodeListPtr->completeSetofCutSets.erase(tempSet);
	myGNodeProbs->normalize();
	double u = ranf();
	double sum = 0.0;
	myGNodeProbs->reset();
	unsigned sampledState=0;
	do{
		double tempProb = myGNodeProbs->getValue();
		sum+= tempProb;
		if (u<sum){
			reset(sampledState);
			resetStateandVector(sampledState);
			sampled = true;
			GNodeList::logProposal += std::log(tempProb);
			//cout << sampledState << endl;
			return;
		}
		sampledState++;
	} while(incr());
}
/*! \fn void GNode::sampleGNode(void)
*  \brief the GNode sampling method (sample the GNode in question)
*/
void GNode::sampleUninformativeTerminalGNode(void){
	// Authors: Rohan L. Fernando 
	// (November, 2005) 
	// Contributors:
	
	if (sampled) return;
	set<GNodeSet*>::iterator setOfGNSetsIt;
	GNodeSet::iterator gnSetIt;
	if (SetofGNsts.size() > 1){
		throw exception("GNode::sampleUninformativeTerminalGNode(void): SetofGNsts.size > 1 \n");
	 } 
	setOfGNSetsIt = SetofGNsts.begin();
	for (gnSetIt=(*setOfGNSetsIt)->begin(); gnSetIt!=(*setOfGNSetsIt)->end(); gnSetIt++){
		if((*gnSetIt)->sampled==false){
			if((*gnSetIt)->uninformativeGNode>0){
				if (this != *gnSetIt) (*gnSetIt)->sampleUninformativeTerminalGNode();
			}
		}
	}
	if(myGNodeProbs->size() == 0){
		myGNodeProbs->insert(this);
		myGNodeProbs->setKeyMultiplicationCode();
	}
	*myGNodeProbs = *setOfGNSetsIt;
	myGNodeProbs->normalize();
	//updateProbs(myGNodeProbs);
	myGNodeProbs->reset();
	double u = ranf();
	double sum = 0.0;
	unsigned sampledState=0;
	do{
		double tempProb = myGNodeProbs->getValue();
		sum+= tempProb;
		if (u<sum){
			reset(sampledState);
			return;
		}
		sampledState++;
	} while(incr());
}

CutSet GNode::calcGNodeProbs(void) {
	// Authors: Rohan L. Fernando 
	// (October, 2005) 
	// Contributors:
	calcMyBCutSet();
	CutSet* tempSet = makeNeighSet();
	tempSet->setKeyMultiplicationCode();
	gNodeListPtr->completeSetofCutSets.insert(tempSet);
	GNodeSet::iterator setit;
	CutSet::iterator cutSetit;
	if(myGNodeProbs->size() == 0) {
		myGNodeProbs->insert(this);
		myGNodeProbs->setKeyMultiplicationCode();
	}
	set< GNodeSet*>::iterator myGNodeSetit;
	myGNodeSetit = SetofGNsts.begin();
	*tempSet = *myGNodeSetit;
	for (myGNodeSetit++;myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
		*tempSet *= *myGNodeSetit;
	}
	*tempSet *= backSet;
	*myGNodeProbs += tempSet;
	delete tempSet;
	gNodeListPtr->completeSetofCutSets.erase(tempSet);
	myGNodeProbs->normalize();
    return *myGNodeProbs;
}
/*! \fn void GNode::calcGNode(void)
*   \brief returns a CutSet with GNode Probabilities
*/

void GNode::calcMyBCutSet(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
	if (backSet) return;
	backSet = new CutSet;
	gNodeListPtr->completeSetofCutSets.insert(backSet);
	if(generatedSet->size()==0){
		backSet->setKeyMultiplicationCode();  
		backSet->putValue(1.0);
	}
	else {
		GNode* nextGNodePtr = generatedSet->getLocation();
		*backSet = nextGNodePtr->calcBCutSetForPrevGNode(generatedSet); 
	}
}
/*! \fn CutSet GNode::calcMyBackCutSet(void)
 *  \brief method used to calculate the backward cutset contribution 
 *  used in calculating the genotype probabilities of a given GNode
*/
CutSet GNode::calcBCutSetForPrevGNode(CutSet* FCutSetPrevGNode) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors:
  calcMyBCutSet(); // result in *backSet
  CutSet* tempSet = makeNeighSet();
  tempSet->setKeyMultiplicationCode();
  gNodeListPtr->completeSetofCutSets.insert(tempSet);
  *tempSet = backSet;
  GNodeSet::iterator setit;
  CutSet BCutSet;
  for (setit=FCutSetPrevGNode->begin();setit!=FCutSetPrevGNode->end();setit++){
    BCutSet.insert(*setit);
  }
  BCutSet.setKeyMultiplicationCode();
  set< GNodeSet*>::iterator myGNodeSetit;
  for (myGNodeSetit=SetofGNsts.begin();myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
    if((*myGNodeSetit)!=FCutSetPrevGNode){
      *tempSet *= *myGNodeSetit;
    }
  }
  BCutSet += tempSet;
  delete tempSet;
  gNodeListPtr->completeSetofCutSets.erase(tempSet);
  return BCutSet;
}
/*! \fn CutSet GNode::calcPreviousGNodeBackCutSet(CutSet* forwardCutSet)
 *  \brief 
*/
GNode::GNode(){
  generatedSet = 0;
  backSet = 0;
  owner = 0;
  myGNodeProbs = new CutSet;
}

GNode::GNode(const GNode& G){
	generatedSet = 0;
	backSet = 0;
	owner = 0;
	if (G.generatedSet) {
		generatedSet = new CutSet;
		generatedSet->replaceMeWith(*G.generatedSet);
	}
	if (G.backSet) {
		backSet = new CutSet;
		backSet->replaceMeWith(*G.backSet);
	}	
	if (G.myGNodeProbs) {
		myGNodeProbs = new CutSet;
		myGNodeProbs->replaceMeWith(*G.myGNodeProbs);
	}
}

void GNode::release(void){

}

void GNodeList::peelOrderCutAndSample(unsigned maxCutSetSize){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: Joseph Abraham
	if(howToSample != "joint") throw exception("Should be sampling jointly !");
	unsigned sampledCount = 0;
	Model *modelPtr = popPtr->model;
	string samplerUsed  = popPtr->model->myRSamplerParms.samplerUsed; 
	string samplerType  = popPtr->model->myRSamplerParms.samplerType;
	string whatToCompute = popPtr->model->myRSamplerParms.whatToCompute;		
	unsigned nsamples = modelPtr->myRSamplerParms.numOfSamples;
	GNodeList::reverse_iterator it;
	setGNodeSampleFlags();
	loopsCut = false;
	GNode* sampleGNodePtr;
	if(samplerUsed != "SG"){
		sampleGNodePtr = peelAndCut(maxCutSetSize);
//	} else{
//		sampleGNodePtr = peelAndCutSG(maxCutSetSize);/ Method to be added later 
	}
	while (sampleGNodePtr) {
	    cout << "GNode to be sampled is  " << (sampleGNodePtr)->id << endl;
		loopsCut = true;
		if (samplerUsed == "P") throw exception("Cannot peel without cutting");
		sampledCount++;
		sampleGNodePtr->sampleGNode();
		cout << "GNode sampled with id  " << sampleGNodePtr->id  << endl;
		if (dynamic_cast<AlleleOriginNode*>(sampleGNodePtr)) AfterSamplingAlleleOriginNode(sampleGNodePtr);
		if (dynamic_cast<AlleleStateNode*>(sampleGNodePtr))  AfterSamplingAlleleStateNode(sampleGNodePtr);
		reinitGNodeList();
		for (unsigned j = 0;j != Individual::numLoci;j++){
			popPtr->LGGenotypeElimination(j);
			popPtr->setAlleleStateVectors(j);
		}
		sampleGNodePtr = peelAndCut(maxCutSetSize);
	}
	cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
	if (samplerUsed != "P" && samplerUsed != "LG"){
		for (it=rbegin();it!=rend();it++){
			if(!(*it)->sampled){
				(*it)->reverseSampleGNode();
			}
		}
		checkSample();
		calculateTargetProb();
		clearGNodeListForNextLocus();
	}else{
	// this part should be moved to another method 
		for(unsigned i = 0;i != nsamples;i++){
			for (it=rbegin();it!=rend();it++){
				if(!(*it)->sampled){
					(*it)->reverseSampleGNode();
				}
			}
			popPtr->copyGNodeStatesToCandidateStates(samplerType);
			popPtr->storeSampledGametes();
			popPtr->copyCandidateToAccepted(samplerType);
			if (whatToCompute=="haplotypeFreq"){
				popPtr->countHaplotypes(samplerType);
			}
			else if (whatToCompute=="genotypeFreq"){
				popPtr->countGenotypes(samplerType);
			}
			setGNodeSampleFlags();
			cout << "Sample No.  " <<  i << "  generated " << endl;
		}
	}
}
/*! \fn void GNodeList::peelOrderCutAndSample()
 *  \brief sample genotypes after determining peeling order, peeling
    and cutting(if necesary)
*/

void GNodeList::setGNodeSampleFlags(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors: 
  GNodeList::iterator it;
  for (it=begin();it!=end();it++){ 
	  (*it)->sampled = false;
  }
}
/*! \fn void GNodeList::setGNodeSampleFlags()
 *  \brief reset sampled flag to false
*/
void GNodeList::setGNodeSampleFlags(bool state){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (August, 2003) 
  // Contributors: 
  GNodeList::iterator it;
  for (it=begin();it!=end();it++){ 
	  (*it)->sampled = state;
  }
}
/*! \fn void GNodeList::setGNodeSampleFlags()
 *  \brief reset sampled flag to false
*/
void GNodeList::peelCutAndSample(unsigned maxCutSetSize){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: Joseph Abraham
	unsigned sampledCount = 0;
    GNodeList::reverse_iterator it;
	setGNodeSampleFlags();
	GNode* SampleGNodePtr;
    sampleGNodePtr = peelAndCut(maxCutSetSize);
	while (sampleGNodePtr){
	    loopsCut = true;
		sampledCount++;
		sampleGNodePtr->sampleGNode();
		cout << " GNode with id  " << sampleGNodePtr->id << " sampled " << endl;
		if (dynamic_cast<AlleleOriginNode*>(sampleGNodePtr)) AfterSamplingAlleleOriginNode(sampleGNodePtr);
		if (dynamic_cast<AlleleStateNode*>(sampleGNodePtr)) AfterSamplingAlleleStateNode(sampleGNodePtr);
		reinitGNodeList();
		for (unsigned j = 0;j != Individual::numLoci;j++){
			popPtr->LGGenotypeElimination(j);
			popPtr->setAlleleStateVectors(j);
		}
		sampleGNodePtr = peelAndCut(maxCutSetSize);
	}
	cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
	for (it=rbegin();it!=rend();it++){
		if(!(*it)->sampled){
			(*it)->reverseSampleGNode();
		}
	}
	checkSample();
	calculateTargetProb();
	clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelCutAndSample()
 *  \brief sample genotypes after peeling and cutting(if necesary)
*/


void GNodeList::peelOrderCutAndSampleMIM(void){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (September, 2005) 
  // Contributors: 
  unsigned printFlag       = MIM::parmMap.getUnsignedValue("printFlag");
  unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
  std::string typeOfGNode  = MIM::parmMap["typeOfGNodes"];
  unsigned sampledCount = 0;
  GNodeList::reverse_iterator it;
  setGNodeSampleFlags(); 
  sampleGNodePtr = peelAndCut(maxCutSetSize);
  GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops = false;
  while (sampleGNodePtr) {
	GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops = true;  
    sampledCount++;
    if(printFlag>0) cout << "GNode sampled after " << sampledCount << " 'th peeling: " << sampleGNodePtr->id << endl;
    sampleGNodePtr->sampleGNode();
    reinitGNodeList();
	if (typeOfGNode=="genotypic"){
		popPtr->LGGenotypeElimination(Individual::currentLocus);
		popPtr->setAlleleStateVectors(Individual::currentLocus);
	}
    sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  if(printFlag>0) cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
  for (it=rbegin();it!=rend();it++){
    if(!(*it)->sampled){
      (*it)->reverseSampleGNode();
    }
  }
  checkSample();
  //cout << "Log Proposal after = " << GNodeList::logProposal << endl;
  calculateTargetProb();
  //cout << "Log Target after = " << GNodeList::logTarget << endl;
  //cout << "Final gNodeList (after peeling, cutting and sampling):" << endl;
  //displayGNodeSets();
  clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelOrderCutAndSampleMIM()
 *  \brief sample genotypes after determining peeling order, peeling
    and cutting(if necesary)
*/
void GNodeList::peelCutAndSampleMIM(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2005) 
	// Contributors: 
	unsigned printFlag       = MIM::parmMap.getUnsignedValue("printFlag");
    unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	std::string typeOfGNode  = MIM::parmMap["typeOfGNodes"];
	unsigned sampledCount = 0;
	GNodeList::reverse_iterator it;
	setGNodeSampleFlags();
    bool cutLoops = GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops;
    if (GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].peelOrder.size()==size() && !cutLoops){
		sampleGNodePtr = peelAndCutFast();
	}
	else{
		// cout << "We determine the peeling order again (inefficient!!)" << endl;
		sampleGNodePtr = peelAndCut(maxCutSetSize);
	}
	//cout << "I am actually sampling a new sample" << endl;
	//cout << "Log Proposal before processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logProposal << endl;
	//cout << "Log Target before = " << GNodeList::logTarget << endl;
	while (sampleGNodePtr) {
		sampledCount++;
		//displayGNodeSets();
		//cout << "GNode sampled after " << sampledCount << " 'th peeling: " << sampleGNodePtr->id << endl;
		sampleGNodePtr->sampleGNode();
		reinitGNodeList();
		if (typeOfGNode=="genotypic"){
			popPtr->LGGenotypeElimination(Individual::currentLocus);
			popPtr->setAlleleStateVectors(Individual::currentLocus);
		}
		sampleGNodePtr = peelAndCut(maxCutSetSize);
	}
	if(printFlag>0) cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
	//cout << "GNodeList before sampling all GNodes" << endl;
	for (it=rbegin();it!=rend();it++){
		if(!(*it)->sampled){
			(*it)->reverseSampleGNode();
		}
	}
	checkSample();
	//cout << "Log Proposal after processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logProposal << endl;
	calculateTargetProb();
	//cout << "Log Target after = " << GNodeList::logTarget << endl;
	//cout << "Final gNodeList (after peeling, cutting and sampling):" << endl;
	//displayGNodeSets();
	clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelCutAndSampleMIM()
*  \brief sample genotypes after peeling and cutting(if necesary)
*/
void GNodeList::peelCutAndComputeMIM(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2005) 
	// Contributors: 

    unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	std::string typeOfGNode  = MIM::parmMap["typeOfGNodes"];
	unsigned sampledCount = 0;
	GNodeList::reverse_iterator it;
	setGNodeSampleFlags();  	
    bool cutLoops = GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].cutLoops;
    if (GNodeSet::prior->chrom()[0].locus[Individual::currentLocus].peelOrder.size()==size() && !cutLoops){
		sampleGNodePtr = peelAndCutFast();
	}
	else{
		// cout << "We determine the peeling order again (inefficient!!)" << endl;
		sampleGNodePtr = peelAndCut(maxCutSetSize);
	}
	// cout << "I am calculating the probability of the previous sample" << endl;
	// cout << "Log Old Proposal before processing locus " << Individual::currentLocus+1 << " = " << GNodeList::logOldProposal << endl;
	//cout << "Log Target before = " << GNodeList::logTarget << endl;
	while (sampleGNodePtr) {
    sampledCount++;
    //displayGNodeSets();
    //cout << "GNode sampled after " << sampledCount << " 'th peeling: " << sampleGNodePtr->id << endl;
    sampleGNodePtr->reComputeGNode();
    reinitGNodeList();
	if (typeOfGNode=="genotypic"){
		popPtr->LGGenotypeElimination(Individual::currentLocus);
		popPtr->setAlleleStateVectors(Individual::currentLocus);
	}
    //displayGNodeSets();
    sampleGNodePtr = peelAndCut(maxCutSetSize);
  }
  //cout << "Number of GNodes sampled to break loops = " << sampledCount << endl;
  //cout << "GNodeList before sampling all GNodes" << endl;
  for (it=rbegin();it!=rend();it++){
    if(!(*it)->sampled){
      (*it)->reComputeGNode();
    }
  }
  checkSample();
  clearGNodeListForNextLocus();
}
/*! \fn void GNodeList::peelCutAndComputeMIM()
 *  \brief
*/
void GNodeList::saveGNodeSets(void){
	GNodeList::iterator it;
	set<GNodeSet*>::iterator GNodeSetiter;
	for (it=begin();it!=end();it++){
		for(GNodeSetiter=(*it)->SetofGNsts.begin();GNodeSetiter!=(*it)->SetofGNsts.end();GNodeSetiter++){
			completeSetofGNsts.insert(*GNodeSetiter);
		}
	}	
}
/*! \fn void GNodeList::saveGNodeSets()
 *  \brief saves the all GNodeSets of the GNodes in completeSetofGNsts
*/
GNodeList GNodeList::getBigList(void){
    GNodeSet bigSet;
	GNodeList bigList;
	GNodeList::iterator it;
	set<GNodeSet*>::iterator GNodeSetiter;
	GNodeSet::iterator gNodeIt;
	for (it=begin();it!=end();it++){
		for(GNodeSetiter=(*it)->SetofGNsts.begin();GNodeSetiter!=(*it)->SetofGNsts.end();GNodeSetiter++){
			for(gNodeIt=(*GNodeSetiter)->begin();gNodeIt!=(*GNodeSetiter)->end();gNodeIt++){
				bigSet.insert(*gNodeIt);
			}
		}
	}
	for (gNodeIt=bigSet.begin();gNodeIt!=bigSet.end();gNodeIt++){
			bigList.push_back(*gNodeIt);
	}
	return bigList;
}
/*! \fn void GNodeList::getBigList()
 *  \brief BigList contains all the GNodes in the GNodeList + all the neighboring GNodes
*/
void GNodeList::resetBackSets(void){
	GNodeList::iterator it;
	for (it=begin();it!=end();it++){
		(*it)->backSet = 0;
	}
}

void GNodeList::updateCounters(void){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	GNodeList::iterator it;
	for(it=begin();it!=end();it++){
		(*it)->updateCounter();
	}
}
/*! \fn void GNodeList::updateCounters(void)
*  \brief update the GNode counters 
*/
void GNodeList::EEUpdateCounters(void){
	// Authors: Rohan L. Fernando 
	// (October, 2005) 
	// Contributors: 
	// I am not calling calcMyCutSet as the "boundry condition" will not work 
	// with blocking. Note that calcBCutSetForPrevGNode(generatedSet) will call
	// calcMyCutSet, but it should return immediately as backset for this GNode 
	// should have already been calculated. 
	
	GNodeList::reverse_iterator it;
	for (it=rbegin();it!=rend();it++){
		(*it)->backSet = new CutSet;
		(*it)->gNodeListPtr->completeSetofCutSets.insert((*it)->backSet);
		if(it==rbegin()){
			(*it)->backSet->setKeyMultiplicationCode();  
			(*it)->backSet->putValue(1.0);
		}
		else {
		    GNode *forwardGNode = findForwardGNode(it,(*it)->generatedSet);
			*(*it)->backSet = forwardGNode->calcBCutSetForPrevGNode((*it)->generatedSet);
		}
		(*it)->EEUpdateCounter(); 	
	}
}
/*! \fn void GNodeList::EEUpdateCounters(void)
*   \brief update the GNode counters with probabilities computed by forward and 
*   backward application of Elston-Stewart Algorithm
*/

GNode* GNodeList::findForwardGNode(GNodeList::reverse_iterator here, CutSet* generatedSet){
	GNodeList::reverse_iterator it, startIt;
	set<GNodeSet*>::iterator itFound;
	startIt = here-1;
	for (it=startIt;it!=(rbegin()-1);it--){
		itFound = (*it)->SetofGNsts.find(generatedSet); 
		if(itFound!=(*it)->SetofGNsts.end()){
			return *it;
		}
	}
	cout << "failed to find the Forward GNode for GNode" << (*here)->id << endl;
	throw exception ("in GNodeList::findForwardGNode \n");
}void GNodeList::reverseSample(void){
	// Authors: Rohan L. Fernando 
	// (October, 2005) 
	// Contributors: 
	GNodeList::reverse_iterator it;
	for (it=rbegin();it!=rend();it++){
		if(!(*it)->sampled){
			if((*it)->getWeight() == 1){
				(*it)->sampled = true;
			}
			else {
				(*it)->reverseSampleGNode();
			}
		}
	}	
}
void GNodeList::sampleSegregationIndicators(void){
	// Authors: Rohan L. Fernando 
	// (October, 2005) 
	// Contributors: 
	GNodeList::iterator it;
	for (it=begin();it!=end();it++){
		Individual *ind = (*it)->owner;
		for (unsigned locus=0;locus<Individual::numLoci;locus++) ind->setSegregationIndex(locus,"genotypic");
		ind->sampleSegregationIndicators();
	}	
}

void GNodeList::pruneUninformativeTerminalGenotypeNodes(void){
	// Authors: Rohan L. Fernando 
	// (November, 2005) 
	// Contributors: 
	GNodeList::iterator it;
	set<GNodeSet*>::iterator setOfGNSetsIt;
	GNodeSet::iterator gnSetIt;
	bool pruned;
	do {
		pruned = false;
		//displayGNodeSets();
		for (it=begin();it!=end();it++){
			if( (*it)->SetofGNsts.size()==1 && (*it)->uninformativeGNode==1 ){
				pruned = true;
				(*it)->uninformativeGNode=2; // pruned
				(*it)->sampled=false;
				setOfGNSetsIt = (*it)->SetofGNsts.begin();
				for (gnSetIt=(*setOfGNSetsIt)->begin(); gnSetIt!=(*setOfGNSetsIt)->end(); gnSetIt++){
					if( *it != *gnSetIt ){
						(*gnSetIt)->SetofGNsts.erase(*setOfGNSetsIt);
					}
				}
			}
		}
	} while (pruned);
}
/*! \fn void GNodeList::pruneUninformativeTerminalNodes(void)
*   \brief A genotype GNode g for a "terminal" offspring has only one element in its SetOfGNsts. 
*    This member is removed from all GNodes except g in that GNodeSet. 
*/

void GNodeList::saveGNodeStates(void){
	// Authors: Rohan L. Fernando 
	// (November, 2005) 
	// Contributors: 
	GNodeList::iterator it;
		for (it=begin();it!=end();it++){
			(*it)->saveState();
		}
}	

void GNodeList::recoverSavedGNodeStates(void){
	// Authors: Rohan L. Fernando 
	// (November, 2005) 
	// Contributors: 
	GNodeList::iterator it;
		for (it=begin();it!=end();it++){
			(*it)->recoverSavedState();
		}
}	

void GNodeList::setUninformativeGNodesNotSampled(){
	// Authors: Rohan L. Fernando 
	// (November, 2005) 
	// Contributors: 
	GNodeList::iterator it;
	for (it=begin();it!=end();it++){
		if ((*it)->uninformativeGNode==2){
			(*it)->sampled = false;
		} 
	}
}

void GNodeList::sampleUpdateUninformativeTerminalGNodes(void){
	// Authors: Rohan L. Fernando 
	// (November, 2005) 
	// Contributors: 
	GNodeList::iterator it;
	for (it=begin();it!=end();it++){
		if ((*it)->sampled==false){
			if((*it)->uninformativeGNode==2){
				(*it)->sampleUninformativeTerminalGNode(); // sample recursively
			} 
		}
	}
	for (it=begin();it!=end();it++){
		if ((*it)->sampled==false){
			if((*it)->uninformativeGNode==2){
				Individual *ind = (*it)->owner;
				for (unsigned locus=0;locus<Individual::numLoci;locus++) ind->setSegregationIndex(locus,"genotypic");
				ind->sampleSegregationIndicators();
				(*it)->updateCounter();
				(*it)->sampled = true;
			}
		}
	}
}	

bool GNodeList::peelNoOrder(unsigned maxCutSetSize){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (May 2006) 
	// Contributors: 
	unsigned printFlag, count=0;
	if (popPtr->model){
		printFlag = popPtr->model->myRSamplerParms.printFlag;
	}
	else {
		printFlag = MIM::parmMap.getUnsignedValue("printFlag");
	}	 		
	GNodeList::iterator it;
	double mag;
	for (it=begin(); it!=end(); it++){
        if((*it)->getWeight() > 1) {
		   mag = (*it)->getCutsetMagnitude();
		   if(mag<=maxCutSetSize){
			   (*it)->updateMysets();
			   if(printFlag>1) cout << "Peeling node " << ++count<<"("<<(*it)->id<<")"<< " cutset magnitude " << mag << endl;
			   (*it)->peel();
		   }
		   else {
		    cerr << "cutset magnitude = " << mag << endl;
			cerr << "maxCutSetSize    = " << maxCutSetSize << endl;
			return false;
		   }
        }
		else {
			if(printFlag>1) cout << "Skipping node " << count<<"("<<(*it)->id<<")"<< endl;
		}
		if(printFlag>2) {
			displayGNodeSets();
		}
	}
	return true;
}

/*! \fn void GNodeList::peelNoOrder()
 *  \brief peels in GNodeList order 
*/

void GNodeList::calcAndDisplayGNodeProbs(void){
	// Authors: Rohan L. Fernando 
	// (August, 2006) 
	// Contributors: 
	
	GNodeList::reverse_iterator it;
	for (it=rbegin();it!=rend();it++){
		LocusStruct* locusPtr = (*it)->locus;
		if (locusPtr->originProbs || locusPtr->stateProbs){
			CutSet myGenotProb = (*it)->calcGNodeProbs();
			cout << (*it)->type << " probs for " << locusPtr->myName <<" locus of ind: " << (*it)->owner->id() << endl;
			myGenotProb.displayValues();
			cout << endl;
		}
	}
}
/*! \fn void GNodeList::calcAndDisplayGNodeProbs(void)
*   \brief calculates and displays GNode probabilities computed by forward and 
*   backward application of Elston-Stewart Algorithm
*/

void GNodeList::calcGNodeProbs(void){
	// Authors: Rohan L. Fernando 
	// (August, 2006) 
	// Contributors: 
	
	GNodeList::reverse_iterator it;
	for (it=rbegin();it!=rend();it++){
		LocusStruct* locusPtr = (*it)->locus;
		if (locusPtr->originProbs || locusPtr->stateProbs){
			(*it)->calcGNodeProbs();
		}
	}
}
/*! \fn void GNodeList::calcGNodeProbs(void)
*   \brief calculates GNode probabilities computed by forward and 
*   backward application of Elston-Stewart Algorithm
*/

bool GNodeList::peelNoCut2(unsigned maxCutSetSize){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2003) 
	// Contributors: 
	unsigned printFlag, count=0;
	if (popPtr->model){
		printFlag = popPtr->model->myRSamplerParms.printFlag;
	}
	else {
		printFlag = MIM::parmMap.getUnsignedValue("printFlag");
	}	 		
	GNodeList::iterator it1,it2,peel_me;
	GNode* temp;
	float mag, min;
	it1=begin();
	while (it1!=end()){
		min = 1e+100;
		for (it2=it1;it2!=end();it2++){
			if((*it2)->getWeight() ==1){
				mag = 1;
			}
			else {
				mag = (*it2)->getCutsetMagnitude();
			}
			if (mag<min){
				peel_me = it2;
				min=mag;
				if (min==1) break;
			}
		} 
		if(min<=maxCutSetSize){
			count++;
			if((*peel_me)->getWeight() > 1) {
				(*peel_me)->updateMysets2();
				(*peel_me)->generatedSet->valueVector.name = "generatedSetFor"+toString((*peel_me)->id,std::dec);
				if(printFlag>1) cout << "Peeling node " << count<<"("<<(*peel_me)->id<<")"<< " cutset magnitude " << min << endl;
				(*peel_me)->peel2();
			}
			else {
				if(printFlag>1) cout << "Skipping node " << count<<"("<<(*peel_me)->id<<")"<< " cutset magnitude " << min << endl;
			}
			temp = *it1;
			*it1 = *peel_me;
			*peel_me = temp;
			it1++;
		}
		else {
		    cerr << "cutset magnitude = " << mag << endl;
			cerr << "maxCutSetSize    = " << maxCutSetSize << endl;
			return false;
		}
		if(printFlag>2) {
			displayGNodeSets();
		}
	}
	return true;
}
/*! \fn void GNodeList::peelNoCut2()
 *  \brief determines peeling order, peels using peel2() 
*/

void GNode::peel2() {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors: Christian Stricker
  neighborSet = makeNeighSet();
  neighborSet->setKeyMultiplicationCode();
  generatedSet->setKeyMultiplicationCode();
  gNodeListPtr->completeSetofCutSets.insert(neighborSet);
  set< GNodeSet*>::iterator myGNodeSetit;
  myGNodeSetit = SetofGNsts.begin();
  *neighborSet = *myGNodeSetit;
  for (myGNodeSetit++;myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
    *neighborSet *= *myGNodeSetit;
  }
  *generatedSet+= neighborSet;
}
/*! \fn void GNode::peel2(void)
 *  \brief the GNode peeling method without releasing cutsets (peeling the GNode in question)
*/

void GNode::updateMysets(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	set< GNodeSet*>::iterator myGNodeSetit;
	GNodeSet::iterator setit;
	CutSet* newSet = makeNeighSet();
	newSet->erase(this);
	gNodeListPtr->completeSetofCutSets.insert(newSet);
	generatedSet = newSet;
	for (myGNodeSetit=SetofGNsts.begin();myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
		for (setit=(*myGNodeSetit)->begin();setit!=(*myGNodeSetit)->end();setit++){
			if(*setit!=this){
				(*setit)->SetofGNsts.erase(*myGNodeSetit);
			}
		}
	}
	for (setit=newSet->begin();setit!=newSet->end();setit++){
		(*setit)->SetofGNsts.insert(newSet);
	}
}
/*! \fn void GNode::updateMysets(void)
 *  \brief creates the new cutset obtained by removing the node peeled
    from the neighborhood cutset. Also sets the pointers of the
    affected nodes to the generated cutset
*/


void GNode::updateMysets2(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	set< GNodeSet*>::iterator myGNodeSetit;
	GNodeSet::iterator setit;
	CutSet* newSet = makeNeighSet();
	newSet->erase(this);
	gNodeListPtr->completeSetofCutSets.insert(newSet);
	generatedSet = newSet;
	for (myGNodeSetit=SetofGNsts.begin();myGNodeSetit!=SetofGNsts.end();myGNodeSetit++){
		for (setit=(*myGNodeSetit)->begin();setit!=(*myGNodeSetit)->end();setit++){
			if(*setit!=this && ((*setit)->getWeight() > 1)){
				(*setit)->SetofGNsts.erase(*myGNodeSetit);
			}
		}
	}
	for (setit=newSet->begin();setit!=newSet->end();setit++){
		if((*setit)->getWeight() > 1) (*setit)->SetofGNsts.insert(newSet);
	}
}
/*! \fn void GNode::updateMysets(void)
 *  \brief creates the new cutset obtained by removing the node peeled
    from the neighborhood cutset. Also sets the pointers of the
    affected nodes to the generated cutset.
*/

}////// end of namespace matvec 
