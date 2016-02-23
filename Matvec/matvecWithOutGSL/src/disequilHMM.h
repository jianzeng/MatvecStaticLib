#include <iostream> 
#include <iomanip>

#include <vector>
#include "population.h"
#include "safe_vectors.h"
using namespace std; 



class Haplotype{
 public:
  unsigned hSize;
  unsigned nHaplotypes;
  matvec::Vector <double> alleleFreqV;
  matvec::Vector <double>  hapProbV;
  vector < vector<unsigned> > hVector;
  matvec::Vector <double> markerPos;
  matvec::IndexVector iVector;
	
  void setup(unsigned size){
    hSize = size;
    SafeSTLVector<unsigned> maxV;
    maxV.resize(size,2);
    iVector.setupVectors(maxV);
    nHaplotypes = pow(2.0,double(hSize));
    hVector.resize(nHaplotypes);
    hapProbV.resize(nHaplotypes,0.0);
    for (unsigned i=0;i<nHaplotypes;i++){
      hVector[i].resize(hSize);
      hVector[i] = iVector.getVector(i);
    }
  }
	
  double markerProb(unsigned h){
    double res = 1.0;
    for (unsigned i=0; i<hSize; i++){
      res *= (hVector[h][i]==0) ? alleleFreqV[i] : (1 - alleleFreqV[i]); 
    }
    return res;
  }	
};	

class RLHap{
 public:
  static Haplotype *hapPtr;
  static unsigned QTLInterval;
  static unsigned mutantHap;
  static unsigned currentHap;
  static unsigned nGen;
  static double mutantFreq;
	
  unsigned length;
  bool right;
  matvec::Vector<double> recombV;
	
	
  void setupRight(){
    right = true;
    length = hapPtr->hSize - QTLInterval;
    recombV.resize(length);
    double x;
    for (unsigned i=0;i<length; i++){
      if(i==0){
	x = (hapPtr->markerPos[QTLInterval] - hapPtr->markerPos[QTLInterval-1])*0.5;
      }
      else{
	x = hapPtr->markerPos[QTLInterval+i] - hapPtr->markerPos[QTLInterval+i-1];		
      }
      recombV[i] = 0.5*(1 - exp(-2.0*x));
    }
  }
	
  void setupLeft(){
    right = false;
    length = QTLInterval;
    recombV.resize(length);
    double x;
    for (unsigned i=0;i<length; i++){
      if(i==0){
	x = (hapPtr->markerPos[QTLInterval] - hapPtr->markerPos[QTLInterval-1])*0.5;
      }
      else{
	x = hapPtr->markerPos[QTLInterval-i] - hapPtr->markerPos[QTLInterval-i-1];		
      }
      recombV[i] = 0.5*(1 - exp(-2.0*x));
    }
  }
	
  double penProb(unsigned Locus, unsigned S){
    unsigned RLLocus;
    if(right){
      RLLocus = Locus + QTLInterval - 1;
    }
    else {
      RLLocus = QTLInterval - Locus;
    }
    double res;
    if (S==0){
      res = (hapPtr->hVector[currentHap][RLLocus]==hapPtr->hVector[mutantHap][RLLocus]) ? 1.0 : 0.0;
    }
    else {
      res = (hapPtr->hVector[currentHap][RLLocus]==0) ? hapPtr->alleleFreqV[RLLocus] : 1-hapPtr->alleleFreqV[RLLocus];
    }
    return res;
  }
	
  double transProb(unsigned Locus, unsigned S1, unsigned S2){
    double res;
    double NR = pow(1-recombV[Locus],double(nGen));
    if (S1==0 && S2==0){
      res = NR + (1 - NR)*mutantFreq;
    }
    else if(S1==0 && S2==1) {
      res = (1 - NR)*(1-mutantFreq);
    }
    else if(S1==1 && S2==0) {
      res = (1 - NR)*mutantFreq;
    }
    else {
      res = NR + (1 - NR)*(1-mutantFreq);
    }
    return res;
  }
	
  double calcMGS(unsigned Locus, unsigned S){
    if(Locus==length) return 1.0;
    double P0 = calcMGS(Locus+1,0);
    double P1 = calcMGS(Locus+1,1);
    double res = transProb(Locus,S,0)*penProb(Locus+1,0)*P0 + transProb(Locus,S,1)*penProb(Locus+1,1)*P1;
    return res;	
  }
};

Haplotype* RLHap::hapPtr;
unsigned RLHap::QTLInterval;
unsigned RLHap::mutantHap;
unsigned RLHap::currentHap;
unsigned RLHap::nGen;
double RLHap::mutantFreq;

class HapProbV{
 public:
  unsigned nGen;
  double mutantFreq;
  Haplotype hap;

  HapProbV(matvec::Vector <double> alleleFreq, matvec::Vector <double> markerPos){
    hap.alleleFreqV = alleleFreq;
    hap.setup(alleleFreq.size());
    hap.markerPos = markerPos;
  }

  matvec::Vector<double> getHapProbV (unsigned mutantHap, unsigned QTLInterval){
    RLHap::hapPtr = &hap;
    RLHap::QTLInterval = QTLInterval;
    RLHap::mutantHap = mutantHap;
    RLHap::nGen = nGen;
    RLHap::mutantFreq = mutantFreq;
    RLHap right, left;
    right.setupRight();
    left.setupLeft();
    for (unsigned i=0; i<hap.nHaplotypes;i++){
      RLHap::currentHap = i;
      double p0 = left.calcMGS(0,0)*right.calcMGS(0,0)*RLHap::mutantFreq;
      double p1 = left.calcMGS(0,1)*right.calcMGS(0,1)*(1-RLHap::mutantFreq);
      hap.hapProbV[i] = p0/(p0+p1);
    }
    return hap.hapProbV;
  }
};	
	



                            




     

 


