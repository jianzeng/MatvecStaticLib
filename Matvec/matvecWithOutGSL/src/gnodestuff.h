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

#ifndef gnodestuff_h
#define gnodestuff_h
#include <ext/hash_map>
#include <map>
#include <set>
#include "geneticdist.h"
#include "exception.h"
#include "safe_vectors.h"


namespace matvec {
  class Population;
  class Individual;
  class GNodeSet; // set of graph nodes
  class CutSet;
  class GNodeList;
  class GNode {  
  public:
    GNode(void); 
	GNode(const GNode& G);
	static CutSet* magSet; // used in get cutset magnitude 
	Individual*  owner;
	LocusStruct* locus;
	string type;
    static GNodeList* gNodeListPtr;
	bool sampled;
	unsigned uninformativeGNode;
	unsigned weight;
    set< GNodeSet*> SetofGNsts;
    CutSet *generatedSet, *backSet, *neighborSet, *myGNodeProbs;
    CutSet* myOldGNodeProbs;
    unsigned peelorder;
    unsigned id;
    unsigned long connectFlag; // used in detecting loops
    unsigned numberOfCuts;
    double getCutsetMagnitude(void);
    void updateMysets(void); 
    void updateMysets2(void); 
    CutSet* makeNeighSet(void);
    CutSet calcBCutSetForPrevGNode(CutSet* FCutSetPrevGNode);
    void   calcMyBCutSet(void);
	CutSet calcGNodeProbs(void);
    void peel(void);
    void peel2(void);
    void reComputeGNode();
    void reverseSampleGNode();
    void sampleGNode();
	void sampleUninformativeTerminalGNode();
    void release(void);
    bool isMyNeighbor(GNode* refGNode);
    SafeSTLVector<int> sampledStateCount;
    virtual unsigned getState(void){
      throw exception("GNode::getState(): called for virtual function");
    };
    virtual unsigned getOldState(void){
      throw exception("GNode::getOldState(): called for virtual function");
    };
    virtual unsigned getAcceptedAlleleState(void){
      throw exception("GNode::getAcceptedAlleleState(): called for virtual function");
    };
    virtual unsigned getMyAlleleState(void){
      throw exception("GNode::getmState(): called for virtual function");
    };
    virtual unsigned getMyAlleleOrigin(void){
      throw exception("GNode::getmState(): called for virtual function");
    }; 
    virtual unsigned getAcceptedMatState(void){
      throw exception("GNode::getAcceptedMatState(): called for virtual function");
    };
    virtual unsigned getAcceptedPatState(void){
      throw exception("GNode::getAcceptedPatState(): called for virtual function");
    };
    virtual unsigned getmState(void){
      throw exception("GNode::getmState(): called for virtual function");
    };
    virtual unsigned getpState(void){
      throw exception("GNode::getpState(): called for virtual function");
    };
    virtual unsigned getWeight(void){
      throw exception("GNode::getWeight(): called for virtual function");
    };	
    virtual bool incr(){
      throw exception("GNode::incr(): called for virtual function");
    };
    virtual void reset(int i){
      throw exception("GNode::reset(int i): called for virtual function");
    };
    virtual void resetAndSwitch(int i){
      throw exception("GNode::resetAndSwitch(int i): called for virtual function");
    };
	virtual void saveState(void){
      throw exception("GNode::saveState(void): called for virtual function");
    };
	virtual void recoverSavedState(void){
      throw exception("GNode::recoverSavedState(void): called for virtual function");
    };
	virtual void updateCounter(void){
		throw exception("GNode::updateCounter(void): called for virtual function");
	};
	virtual void EEUpdateCounter(void){
		throw exception("GNode::EEUpdateCounter(void): called for  virtual function");
	};
		virtual void updateProbs(CutSet& myGenotProb){
		throw exception("GNode::updateProbs(CutSet& ): called for  virtual function");
	};
	virtual void  resetStateandVector(int i){
	     throw exception("GNode::resetStateandVector(int i): called for virtual function");
		 };
    virtual ~GNode(void) {release();}
	virtual void display(void){
		throw exception("GNode::display(): called for virtual function");
	}
  };
  /*! \class GNode gnodestuff.h inc/gnodestuff.h
   * \brief This is the base graph node class.
   *
   * GNode class includes methods relevant to a graph node
   */
  class GNodeSet:public set<GNode*>{
  public:
    GNodeSet(void);
    unsigned long connectFlag; // used in detecting loops
    unsigned numberOfCuts;
    static GeneticDist *prior;
    static unsigned currentLocus;
    bool incr(void);
    void reset(void);
    void attachMeToMyGnodes(void);
	void display(void);
    virtual double getValue(void){
      throw exception("GNodeSet::getValue(void): call for virtual function");
    };
    virtual double getTargetValue(void){
      throw exception("GNodeSet::getTargetValue(void): call for virtual function");
    };
    virtual double getOldValue(void){
      throw exception("GNodeSet::getOldValue(void): call for virtual function");
    }; 
    virtual ~GNodeSet(void) {;}
  };
  /*! \class GNodeSet gnodestuff.h inc/gnodestuff.h
   * \brief This is the  base "set of graph nodes" class.
   *
   * GNodeSet class includes methods relevant to a set of graph nodes
   */
  struct compareGNodesWeight: public binary_function<GNode*,GNode*,bool>{
    bool operator()(GNode* a,GNode* b){
		unsigned aw = a->sampled ? 1 : a->getWeight();
		unsigned bw = b->sampled ? 1 : b->getWeight();		
      return aw>bw;
    }
  };
  /*! \struct compareGNodes gnodestuff.h inc/gnodestuff.h
   * \brief This is the comparison criterion to rank Gnodes in GNodeList
   * based on the size of the alleleStateVector or genotypeVector
   */
    struct compareGNodesPeelId: public binary_function<GNode*,GNode*,bool>{
    bool operator()(GNode* a,GNode* b){
      return a->peelorder<b->peelorder;
    }
  };
  /*! \struct compareGNodes gnodestuff.h inc/gnodestuff.h
   * \brief This is the comparison criterion to rank Gnodes in GNodeList
   * based on the peeling order of the GNode
   */
  struct GNode_info {
    GNode* pointerToGNode; 
    unsigned dist;
  };
  /*! \struct GNode_info gnodestuff.h inc/gnodestuff.h
   * \brief Stores the id and distance from the cut of the GNode that will be sampled first
   */  

  class InvalidSample{}; 
  class NormalizeFailed{};

  class GNodeList:public SafeSTLVector<GNode*> {
  public:
    // from here
    static double logTarget;
    static double logProposal;
    static double logOldProposal;
	static Population *popPtr;
    float sum;
	bool loopsCut;
    unsigned lastmember;
	Matrix<int> distanceMat;
	bool distanceMatDone;
	SafeSTLVector<GNode*> cutGNodesVector;
	set<GNode*> loopSetGNodes;
	GNodeList(void){distanceMatDone = false;};
    void  inputGNodeSets(char* fname);
    void  displayGNodeSets();
    void  getPeelOrder(char* fname);
    void  fill(unsigned n);
    void  makealleleGNodeSets(char* infile,char* outfile);
    void  makegenotGNodeSets(char* infile,char* outfile);
    // up to here are methods from the original peeling order program
    unsigned maxDist;
    GNode* sampleGNodePtr;
    GNode* containerGNode;
    set< GNodeSet*> completeSetofGNsts;
    set< CutSet*> completeSetofCutSets;
    GNodeList::iterator searchStart;
    GNodeList::iterator cutGNode;
	void AfterSamplingAlleleOriginNode(GNode* sampleGNodePtr);
	void AfterSamplingAlleleStateNode(GNode* sampleGNodePtr);
	bool peelNoCut(unsigned maxCutSetSize);
	bool peelNoOrder(unsigned maxCutSetSize);
    GNode* peelAndCut(unsigned maxCutSetSize);
    GNode* peelAndCutFast(void);
    void   peelCutAndCompute(unsigned maxCutSetSize, unsigned startBlock,unsigned stopBlock,unsigned sizeBlock);
    void   peelCutAndSample(unsigned maxCutSetSize, unsigned startBlock,unsigned stopBlock,unsigned sizeBlock);
    void   peelCutAndSample(unsigned maxCutSetSize);
    void   peelOrderCutAndSample(unsigned maxCutSetSize,unsigned startBlock,unsigned stopBlock,unsigned sizeBlock);
    void   peelOrderCutAndSample(unsigned maxCutSetSize);
	void   peelOrderCutAndSampleMIM(void);
	void   peelCutAndSampleMIM(void);
	void   peelCutAndComputeMIM(void);
    void   clearGNodeListForNextLocus(void);
    void   reinitGNodeList(void);
    void   releaseGNsts(void);
    void   releaseCutSets(void);
    void   resetGNodeList(void);
    void   setGNodeSampleFlags(unsigned startBlock, unsigned stopBlock, unsigned sizeBlock);
    void   setGNodeSampleFlags(void);
    void   findLoopAndCut();
	GNodeList::iterator getMinWeightGNode(void);
	set< GNodeSet*>::iterator  getCutGNodeSet(GNodeList::iterator cutGNode);
    void  propagateFlags(unsigned lowFlag, unsigned highFlag);
    void  resetConnectFlags(void);
    GNode_info choosenGNodeInfo;
	GNodeList isolatePureLoop(GNodeList::iterator loopCloser);
    void  cutLoop(GNodeList::iterator cutGNode, set< GNodeSet*>::iterator cutGNodeSet);
    void checkSample(void);
    void calculateTargetProb(void);
    string howToSample;
	void makeDistanceMatrix(void);
	void calcDistancefrom(GNode* refGNode);
	GNodeList findNeighbors(GNode* refGNode, unsigned largestdist);
	void saveGNodeSets(void);
	GNode* findSampleGNode(void);
	void reverseSample(void);
	void updateCounters(void);
	void EEUpdateCounters(void);
	GNode* findForwardGNode(GNodeList::reverse_iterator here, CutSet* generatedSet);
	void countGenotypes(string nodeType);
	GNodeList getBigList(void);
	void resetBackSets(void);
	void pruneUninformativeTerminalGenotypeNodes(void);
	void setGNodeSampleFlags(bool state);
	void saveGNodeStates(void);
	void recoverSavedGNodeStates(void);
	void setUninformativeGNodesNotSampled(void);
	void sampleUpdateUninformativeTerminalGNodes(void);
	void sampleSegregationIndicators(void);
	void calcAndDisplayGNodeProbs(void);
	void calcGNodeProbs(void);
	bool peelNoCut2(unsigned maxCutSetSize);
  };
  /*! \class GNodeList gnodestuff.h inc/gnodestuff.h
   * \brief This is the "list of graph nodes" class.
   *
   * GNodeList class includes methods relevant to the list of graph
   * nodes
   */
}////// end of namespace matvec 
#endif
