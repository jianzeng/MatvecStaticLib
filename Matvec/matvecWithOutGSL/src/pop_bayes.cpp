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

#include "model.h"
#include "population.h"
#include "statdist.h"
using namespace std;

namespace matvec {
	
	void Population::identifyLociToSamplePosFor(void){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (December, 2004) 
		// Contributors:
		samplePosForMe = ULONG_MAX;
		for (unsigned k=0;k<Individual::numLoci;k++){
			if (prior->chrom()[0].locus[k].samplePositionForThisLocus){
				cout << "I am sampling a new position for locus: " << k << endl; 
				samplePosForMe = k;
			}
		}
	}
	
	void Population::samplePositionOnChromosome(void){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (November, 2004) 
		// Contributors:
		for (unsigned k=0;k<Individual::numLoci;k++){
			if (prior->chrom()[0].locus[k].samplePositionForThisLocus){
				// cout << "I am sampling a new position for locus: " << k << endl; 
				samplePosForMe = k;
				double r1,r2;
				// get ready to obtain candidate from a Uniform proposal 
				double delta=0.004;
				double leftEdge  = (prior->chrom()[0].locus[samplePosForMe].distance)-delta;
				double rightEdge = (prior->chrom()[0].locus[samplePosForMe].distance)+delta;
				leftEdge  = (leftEdge  <  0.0) ?  0.0 : leftEdge;
				rightEdge = (rightEdge > prior->chrom()[0].myLength) ? prior->chrom()[0].myLength : rightEdge;
				// calculate pdf for going from old to new
				logNewProposalForPosition = 1/(rightEdge-leftEdge);
				// sample new position
				UniformDist Uniform(leftEdge,rightEdge);
				prior->chrom()[0].locus[samplePosForMe].distance=Uniform.sample();
				// reconstruct the recombination matrix for the new position
				prior->chrom()[0].calcLocusMask();
				buildRecombinationMatrix();
				// get ready to calculate old proposal
				leftEdge  = (prior->chrom()[0].locus[samplePosForMe].distance)-delta;
				rightEdge = (prior->chrom()[0].locus[samplePosForMe].distance)+delta;
				leftEdge  = (leftEdge  <  0.0) ?  0.0 : leftEdge;
				rightEdge = (rightEdge > prior->chrom()[0].myLength) ? prior->chrom()[0].myLength : rightEdge;
				// calculate probability of going from new to old
				logOldProposalForPosition = 1/(rightEdge-leftEdge); 
			}
		}
	}
	
} // end of matvec namespace
