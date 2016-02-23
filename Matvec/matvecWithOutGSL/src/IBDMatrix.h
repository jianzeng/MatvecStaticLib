/*
 *  IBDMatrix.h
 *  OBGML
 *
 *  Created by Rohan Fernando on 1/23/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef IBDMATRIX_H
#define IBDMATRIX_H

#include "safe_vectors.h"
#include "matrix.h"
using namespace std; 

namespace matvec{
	class Genotype {
	public:
		unsigned maternalAllele;
		unsigned paternalAllele;
		Genotype(){
			maternalAllele = 0;
			paternalAllele = 0;
		}
	};
	class IBDMatrix:public Matrix<double> {
	public:
	    unsigned count;
	    Vector<Genotype> genotypeVec;
		unsigned leftLocus, rightLocus;
		unsigned popSize;
		double r1,r2,rab;
		void initialize(unsigned ll, unsigned rl);
		void sample(void);
		void sampleSimple(void);
		void sample1(void);
		void sample1Simple();
		void display(void);	
		void output(string outFileName);
	};
	
}////// end of namespace matvec 
#endif
