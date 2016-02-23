/*
 *  IBDMatrix.cpp
 *  OBGML
 *
 *  Created by Rohan Fernando on 1/23/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "IBDMatrix.h"
#include "population.h"
#include "stat.h"

using namespace std;

namespace matvec{

	void IBDMatrix::initialize(unsigned ll, unsigned rl){
		Population *popPtr;
		popPtr = GNodeList::popPtr;
		leftLocus  = ll;
		rightLocus = rl;
		double xab = std::abs(popPtr->prior->get_distance(1,leftLocus) 
							  - popPtr->prior->get_distance(1,rightLocus));
		rab = 0.5*(1-std::exp(-2.0*xab)); // Haldane
		double x1 = xab/2;                // QTL at center;
		r1  = r2 = 0.5*(1-std::exp(-2.0*x1)); 					
		count = 0;
		popSize = popPtr->popsize;
		genotypeVec.resize(popSize);
		resize(popSize,popSize,0.0);
	}
	

void IBDMatrix::sample(){
	// assume QTL at midpoint between markers
	unsigned founderCount = 1;
	double noR = (1-r1)*(1-r2)/(1-rab);
	Population *popPtr;
	popPtr = GNodeList::popPtr;
	
	for (unsigned t=0;t<popSize;t++){
		Individual *ind = popPtr->popmember[t];
		Individual *mom = ind->mymother;
		Individual *dad = ind->myfather;
		if (mom) {
			unsigned momID = mom->id();
			unsigned momMatAllele = genotypeVec[momID-1].maternalAllele;
			unsigned momPatAllele = genotypeVec[momID-1].paternalAllele;
			unsigned leftOrigin  = ind->malleleOriginNodeVector[ leftLocus-1].getAcceptedAlleleOrigin();
			unsigned rightOrigin = ind->malleleOriginNodeVector[rightLocus-1].getAcceptedAlleleOrigin();
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].maternalAllele = myOrigin ? momPatAllele : momMatAllele;
		}
		else {
			genotypeVec[t].maternalAllele = founderCount++;
		}
		if (dad) {
			unsigned dadID = dad->id();
			unsigned dadMatAllele = genotypeVec[dadID-1].maternalAllele;
			unsigned dadPatAllele = genotypeVec[dadID-1].paternalAllele;
			unsigned leftOrigin  = ind->palleleOriginNodeVector[ leftLocus-1].getAcceptedAlleleOrigin();
			unsigned rightOrigin = ind->palleleOriginNodeVector[rightLocus-1].getAcceptedAlleleOrigin();
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].paternalAllele = myOrigin ? dadPatAllele : dadMatAllele;
		}
		else {
			genotypeVec[t].paternalAllele = founderCount++;
		}
	}
	for (unsigned it=0;it<popSize;it++){
		unsigned itMatAllele = genotypeVec[it].maternalAllele;
		unsigned itPatAllele = genotypeVec[it].paternalAllele;
		for (unsigned jt=it;jt<popSize;jt++){
			unsigned jtMatAllele = genotypeVec[jt].maternalAllele;
			unsigned jtPatAllele = genotypeVec[jt].paternalAllele;	
			(*this)[it][jt] += (itMatAllele==jtMatAllele) ? 1 : 0;
			(*this)[it][jt] += (itMatAllele==jtPatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtMatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtPatAllele) ? 1 : 0;	
		}
	}
	count++;
}

void IBDMatrix::sampleSimple(){
	// assume QTL at midpoint between markers
	unsigned founderCount = 1;
	double noR = (1-r1)*(1-r2)/(1-rab);
	Population *popPtr;
	popPtr = GNodeList::popPtr;
	
	for (unsigned t=0;t<popSize;t++){
		Individual *ind = popPtr->popmember[t];
		Individual *mom = ind->mymother;
		Individual *dad = ind->myfather;
		if (mom) {
			unsigned momID = mom->id();
			unsigned momMatAllele = genotypeVec[momID-1].maternalAllele;
			unsigned momPatAllele = genotypeVec[momID-1].paternalAllele;
			unsigned leftOrigin  = ind->malleleOriginNodeVector[ leftLocus-1].getMyAlleleOrigin();
			unsigned rightOrigin = ind->malleleOriginNodeVector[rightLocus-1].getMyAlleleOrigin();
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].maternalAllele = myOrigin ? momPatAllele : momMatAllele;
		}
		else {
			genotypeVec[t].maternalAllele = founderCount++;
		}
		if (dad) {
			unsigned dadID = dad->id();
			unsigned dadMatAllele = genotypeVec[dadID-1].maternalAllele;
			unsigned dadPatAllele = genotypeVec[dadID-1].paternalAllele;
			unsigned leftOrigin  = ind->palleleOriginNodeVector[ leftLocus-1].getMyAlleleOrigin();
			unsigned rightOrigin = ind->palleleOriginNodeVector[rightLocus-1].getMyAlleleOrigin();
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].paternalAllele = myOrigin ? dadPatAllele : dadMatAllele;
		}
		else {
			genotypeVec[t].paternalAllele = founderCount++;
		}
	}
	for (unsigned it=0;it<popSize;it++){
		unsigned itMatAllele = genotypeVec[it].maternalAllele;
		unsigned itPatAllele = genotypeVec[it].paternalAllele;
		for (unsigned jt=it;jt<popSize;jt++){
			unsigned jtMatAllele = genotypeVec[jt].maternalAllele;
			unsigned jtPatAllele = genotypeVec[jt].paternalAllele;	
			(*this)[it][jt] += (itMatAllele==jtMatAllele) ? 1 : 0;
			(*this)[it][jt] += (itMatAllele==jtPatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtMatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtPatAllele) ? 1 : 0;	
		}
	}
	count++;
}


void IBDMatrix::sample1(){
	// in this version allele origin is assumed to be stored in m_gameteAccepted and p_gameteAccepted vectors
	// assume QTL at midpoint between markers
	unsigned founderCount = 1;
	double noR = (1-r1)*(1-r2)/(1-rab);
	Population *popPtr;
	popPtr = GNodeList::popPtr;
	for (unsigned t=0;t<popSize;t++){
		Individual *ind = popPtr->popmember[t];
		Individual *mom = ind->mymother;
		Individual *dad = ind->myfather;
		if (mom) {
			unsigned momID = mom->id();
			unsigned momMatAllele = genotypeVec[momID-1].maternalAllele;
			unsigned momPatAllele = genotypeVec[momID-1].paternalAllele;
			unsigned leftOrigin  = ind->m_gameteAccepted[ leftLocus-1];
			unsigned rightOrigin = ind->m_gameteAccepted[rightLocus-1];
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].maternalAllele = myOrigin ? momPatAllele : momMatAllele;
		}
		else {
			genotypeVec[t].maternalAllele = founderCount++;
		}
		if (dad) {
			unsigned dadID = dad->id();
			unsigned dadMatAllele = genotypeVec[dadID-1].maternalAllele;
			unsigned dadPatAllele = genotypeVec[dadID-1].paternalAllele;
			unsigned leftOrigin  = ind->p_gameteAccepted[ leftLocus-1];
			unsigned rightOrigin = ind->p_gameteAccepted[rightLocus-1];
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].paternalAllele = myOrigin ? dadPatAllele : dadMatAllele;
		}
		else {
			genotypeVec[t].paternalAllele = founderCount++;
		}
	}
	for (unsigned it=0;it<popSize;it++){
		unsigned itMatAllele = genotypeVec[it].maternalAllele;
		unsigned itPatAllele = genotypeVec[it].paternalAllele;
		for (unsigned jt=it;jt<popSize;jt++){
			unsigned jtMatAllele = genotypeVec[jt].maternalAllele;
			unsigned jtPatAllele = genotypeVec[jt].paternalAllele;	
			(*this)[it][jt] += (itMatAllele==jtMatAllele) ? 1 : 0;
			(*this)[it][jt] += (itMatAllele==jtPatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtMatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtPatAllele) ? 1 : 0;	
		}
	}
	count++;
}


void IBDMatrix::display(void){
	for (unsigned it=0;it<popSize;it++){
		for (unsigned jt=0;jt<popSize;jt++){
			cout << setw(5) << (*this)[it][jt]/(2.0*count) <<" ";
	    }
		cout << endl;
	}	
}
void IBDMatrix::output(string outFileName){
	std::ofstream outFile(outFileName.c_str());
	for (unsigned it=0;it<popSize;it++){
		for (unsigned jt=it+1;jt<popSize;jt++){
			outFile << setw(5) << it+1 << " "
					<< setw(5) << jt+1 << " "  
					<< setw(5) << (*this)[it][jt]/(2.0*count) <<"\n";
	    }
	}	
}


void IBDMatrix::sample1Simple(){
	// in this version allele origin is assumed to be stored in m_gamete and p_gamete vectors
	// assume QTL at midpoint between markers
	unsigned founderCount = 1;
	double noR = (1-r1)*(1-r2)/(1-rab);
	Population *popPtr;
	popPtr = GNodeList::popPtr;
	for (unsigned t=0;t<popSize;t++){
		Individual *ind = popPtr->popmember[t];
		Individual *mom = ind->mymother;
		Individual *dad = ind->myfather;
		if (mom) {
			unsigned momID = mom->id();
			unsigned momMatAllele = genotypeVec[momID-1].maternalAllele;
			unsigned momPatAllele = genotypeVec[momID-1].paternalAllele;
			unsigned leftOrigin  = ind->m_gamete[ leftLocus-1];
			unsigned rightOrigin = ind->m_gamete[rightLocus-1];
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].maternalAllele = myOrigin ? momPatAllele : momMatAllele;
		}
		else {
			genotypeVec[t].maternalAllele = founderCount++;
		}
		if (dad) {
			unsigned dadID = dad->id();
			unsigned dadMatAllele = genotypeVec[dadID-1].maternalAllele;
			unsigned dadPatAllele = genotypeVec[dadID-1].paternalAllele;
			unsigned leftOrigin  = ind->p_gamete[ leftLocus-1];
			unsigned rightOrigin = ind->p_gamete[rightLocus-1];
			unsigned myOrigin;
			if (leftOrigin==rightOrigin){
				double u = ranf();
				if (u < noR){
					myOrigin = leftOrigin;			
				}
				else {
					myOrigin = leftOrigin ? 0 : 1;
				}
				
			}
			else {
				double u = ranf();
				myOrigin = (u < 0.5) ? 0 : 1;
			}
			genotypeVec[t].paternalAllele = myOrigin ? dadPatAllele : dadMatAllele;
		}
		else {
			genotypeVec[t].paternalAllele = founderCount++;
		}
	}
	for (unsigned it=0;it<popSize;it++){
		unsigned itMatAllele = genotypeVec[it].maternalAllele;
		unsigned itPatAllele = genotypeVec[it].paternalAllele;
		for (unsigned jt=it;jt<popSize;jt++){
			unsigned jtMatAllele = genotypeVec[jt].maternalAllele;
			unsigned jtPatAllele = genotypeVec[jt].paternalAllele;	
			(*this)[it][jt] += (itMatAllele==jtMatAllele) ? 1 : 0;
			(*this)[it][jt] += (itMatAllele==jtPatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtMatAllele) ? 1 : 0;	
			(*this)[it][jt] += (itPatAllele==jtPatAllele) ? 1 : 0;	
		}
	}
	count++;
}



} ////// end of namespace matvec


