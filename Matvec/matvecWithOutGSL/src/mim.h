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

#ifndef MIM_H
#define MIM_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <map>
#include "doublematrix.h"
#include "geneticdist.h"
#include "vector.h"
#include "session.h"
#include "parmMap.h"
#include "rutil.h"
#include "rpedigree.h"
#include "mqtl.h"
#include "IBDMatrix.h"
#include "individual.h"
#include "BNodeList.h"

// classes for multiple trait, mixed models with iteration on data

using namespace std; 

namespace matvec{
 
 class TermData{
  public:
    double value;
    unsigned level;
  };
  
  class MimDataNode{
  public:
    SafeSTLVector<TermData > trmVec;
    SafeSTLVector<double>    depVec;
	std::vector<bool>        misVec;
	doubleMatrix      *RiPtr;
  };
  
  class MIM;
  
  class MIModelTerm{
  public:
    unsigned start;
    unsigned trait;
    string name;
    string depVarName;
    static MIM *myMIMPtr;
	bool pedTerm;
    Recoder<string>* myRecoderPtr;
    MQTL* myMQTLPtr;
    bool secondMQTLEffect;
    SafeSTLVector<unsigned> factors;
    
    unsigned code(string str){return myRecoderPtr->code(str);}
    unsigned nLevels(){return (unsigned)myRecoderPtr->size();}
    void putFactors(string str);
    string   getTermString();
    unsigned getTermLevel (){
      return code(getTermString());
    }
    double   getTermValue ();
  };
  
  class CovBlock {
  public:
    SafeSTLVector<MIModelTerm*> modelTrmPtrVec;
    matvec::doubleMatrix Var, Vari;
    RPedigree* pedPtr;
    MQTL* myMQTLPtr;
    CovBlock(void){
		pedPtr = 0; 
		myMQTLPtr = 0;
	}
    CovBlock(string str, matvec::doubleMatrix V){
      Var = V;
	  pedPtr = 0;
	  myMQTLPtr = 0;
      buildModelTrmVec(str);
    }
    CovBlock(string str, matvec::doubleMatrix V, RPedigree &P){
      Var = V;
      pedPtr = &P;
      myMQTLPtr = 0;
      buildModelTrmVec(str);
    }
    CovBlock(MQTL &mQTL){
      pedPtr = 0;
      myMQTLPtr = &mQTL;
    }
    void buildModelTrmVec(string str);
    void addGinv(void);
  };
  
  class Population;
  class MIM {
  private:
    void putModel(string str);
  public:
    MIM(void) {initialize();}
    ~MIM(void){release();}
    void initialize(void);
    void release(void);
    Population *pop;
	bool printResVar;
    int pop_created;
    static string solMethod;
    static matvec::Vector<double> *vec, diag, res;
	static ParmMap parmMap;
	static IBDMatrix ibdMatrix;
    
	matvec::Vector<double> resid;
    string phenDataFileName;
    string markDataFileName;
    Tokenizer colType;
    Tokenizer colName;
    Tokenizer markerColName;
    Tokenizer depVar;
    Tokenizer colData;
    unsigned numCols;
    unsigned markerNumCols;
    SafeSTLVector <MIModelTerm> modelTrmVec;
    SafeSTLVector <CovBlock>  covBlockVec;
    SafeSTLVector <MimDataNode> dataVec;
    SafeSTLVector<MQTL*> MQTLVec;
    unsigned numTerms, numTraits;
    unsigned mmeSize;
    matvec::doubleMatrix lhs, R;
	map<unsigned,doubleMatrix> RMap;
    matvec::Vector<double> rhs, sol, tempSol;
    
    void putColNames(string str);
    void putMarkerColNames(string str);
    void putColTypes(string str);
    void putModels(string str);
	void isMIMReady(void);
    void putVarCovMatrix(string str, matvec::doubleMatrix V, RPedigree &P){
      CovBlock covBlock(str,V,P);
      covBlockVec.push_back(covBlock);
    }
    void putVarCovMatrix(string str, matvec::doubleMatrix V){
      CovBlock covBlock(str,V);
      covBlockVec.push_back(covBlock);
    }
    void putVarCovMatrix(MQTL &Q){
      CovBlock covBlock(Q);
      covBlockVec.push_back(covBlock);
    }
    void initSetup();
    void inputData();
	doubleMatrix getMissR(SafeSTLVector<unsigned> missVector);
    void displayData();
    static double getDouble(string& Str);
	static int    getInteger(string& Str);
    void calcStarts();
    void getDirectSolution();
    void getJacobiSolution(double p = 0.5);
    void getCGSolution  (double eps=0.00001);
    void getPCCGSolution(double eps=0.00001);
    
    void mmeTimes(matvec::Vector<double>& x);
    void calcWPW();
    void addGinv();
    void display(std::string str="");
	void displayGenotypicValues(RPedigree &ped, std::ostream &outfile=cout);
    
    void putMQTL(MQTL& Q){MQTLVec.push_back(&Q);};
    void mergeMQTLLevelsBy(string mergeStr);
    void addColNames(string str);
    void addColTypes(string str);
    void putMQTLStuffInModelTerms(MQTL& mQTL);
    void processMQTL(void);

    // GNodeSampler stuff comes next:
	unsigned individualModelTerm;
	void runAnalysis(std::string inputFileName,RPedigree& P,GeneticDist& G);
	bool monoMendel(GeneticDist& G);
    void setupPopulation(RPedigree& P,GeneticDist& G);
	void transferTraitData();
	void putGNodeSamplerParms(std::string inputFileName);
	static matvec::doubleMatrix IBDCovMatrix;
	static matvec::doubleMatrix meansVector;
	void GNodeSamplerGS1(void);
	void GNodeSamplerOBGSL(void);
	void GNodeSamplerOBGSL1(void);
	void GNodeSamplerEEGSL(void);
	void GNodeSamplerJNC(void);
	void GNodeSamplerJNCIDG(void);
	void GNodeSamplerLG(void);
	void GNodeSamplerSLNC(void);
	void GNodeSamplerJointNC(void);
	void GNodeProbsAPESJointNC(std::string inputFileName,RPedigree& P,GeneticDist& G);
	void DGSampleIBDMatrix(unsigned numOfSamples, unsigned numSL, unsigned numHaplo, unsigned numSLCas, unsigned leftLocus,
                            unsigned rightLocus, string fileName);
    void GNodeSamplerOBGML(void);
	void makePedigreeBlocks(void);	
	void mergeAdjacentPedBlocks(SafeSTLVector<Individual*>& pedBlockPivots);
	void removeRedundantBlocks(void);
	void makeBlockGNodeListVector(void);	
	void addBlocksFor(unsigned startLocus, unsigned endLocus);
	SafeSTLVector<BlockNodeList> blockNodeListVector;
	RPedigree *pedPtr;
					
  };
}
#endif

