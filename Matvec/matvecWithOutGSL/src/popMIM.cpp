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
/*
 *  popMIM.cpp
 *
 *  Created by Rohan Fernando on 9/3/05.
 *
 */

#include "model.h"
#include "population.h"
#include "nufamily.h"
#include "stat.h"
#include "mim.h"
#include <iomanip>
#include <sstream>
using namespace std;

namespace matvec {
	
	void Population::getInitialGNodeListSampleMIM(void){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		unsigned printFlag  = MIM::parmMap.getUnsignedValue("printFlag"); 
		unsigned startLocus = MIM::parmMap.getUnsignedValue("startLocus") - 1;
		// get a sample for the startLocus
		lookToYourLeft  = false;
		lookToYourRight = false;
		Individual::currentPosition = startLocus;
		unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j << endl;
		getInitialSampleForLocusMIM(j);
		// get samples for the other loci 
		// first to the right
		for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
			lookToYourLeft  = true;
			lookToYourRight = false;
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j << endl;
			getInitialSampleForLocusMIM(j);
		}
		// next to the left
		for (int jj=startLocus-1;jj>=0;jj--){
			lookToYourLeft  = false;
			lookToYourRight = true;
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j << endl;
			getInitialSampleForLocusMIM(j);
		}
	}
	/*! \fn void Population::getInitialGNodeListSampleMIM(void)
	*  \brief get the initial sample and also determine and store the peeling and sampling order for each locus starting with an arbitrary locus
	*/
	
	void Population::getInitialSampleForLocusMIM(unsigned j){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		std::string typeOfGNode = MIM::parmMap["typeOfGNodes"];	
		if (typeOfGNode=="genotypic"){
			prior->chrom()[0].locus[j].peelOrder.resize(popsize);
			initGenotypeNodeList(j);
//			set< GNodeSet*>::iterator GNodeSetit;
//			for (GNodeSetit=gNodeList.completeSetofGNsts.begin();GNodeSetit!=gNodeList.completeSetofGNsts.end();GNodeSetit++){ 
//				(*GNodeSetit)->display();
//			}
		} 
		else if (typeOfGNode=="allelic"){
			prior->chrom()[0].locus[j].peelOrder.resize(2*popsize);
			initAlleleNodeList(j);
			//gNodeList.displayGNodeSets();
		}
		else {
			cerr << "Population::getInitialSampleForLocusMIM: typeOfGNodes " << typeOfGNode << " not allowed "<< endl;
			throw exception("typeOfGNodes must be alleleic or genotypic");
		}
		gNodeList.peelOrderCutAndSampleMIM(); 
		setSegregationIndex(j,typeOfGNode);
	}
	/*! \fn void Population::getInitialSampleForLocusMIM(unsigned j)
	*  \brief get the initial sample and also determine and store the peeling and sampling order for a given locus
	*/

	void Population::getGNodeListSampleGibbsMIM(void){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		unsigned printFlag  = MIM::parmMap.getUnsignedValue("printFlag"); 
		unsigned startLocus = MIM::parmMap.getUnsignedValue("startLocus") - 1;
		// get a sample for the startLocus
		lookToYourLeft  = true;
		lookToYourRight = true;
		Individual::currentPosition = startLocus;
		unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j << endl;
		getSampleForLocusMIM(j);
		// get samples for the other loci 
		// first to the right
		for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j << endl;
			getSampleForLocusMIM(j);
		}
		// next to the left
		for (int jj=startLocus-1;jj>=0;jj--){
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j << endl;
			getSampleForLocusMIM(j);
		}
	}
	/*! \fn void Population::getGNodeListSampleMIM(void)
	*  \brief obtain a new candidate sample using peeling and cutting (if necessary)  
	*/
		
	void Population::getGNodeListSampleSeqInMIM(void){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		unsigned printFlag  = MIM::parmMap.getUnsignedValue("printFlag"); 
		unsigned startLocus = MIM::parmMap.getUnsignedValue("startLocus") - 1;
		// get a sample for the startLocus
		lookToYourLeft  = false;
		lookToYourRight = false;
		Individual::currentPosition = startLocus;
		unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
		getSampleForLocusMIM(j);
		// get samples for the other loci 
		// first to the right
		for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
			lookToYourLeft  = true;
			lookToYourRight = false;
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
			getSampleForLocusMIM(j);
		}
		// next to the left
		for (int jj=startLocus-1;jj>=0;jj--){
			lookToYourLeft  = false;
			lookToYourRight = true;
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
			getSampleForLocusMIM(j);
		}
	}
	/*! \fn void Population::getGNodeListSampleMIM(void)
	*  \brief obtain a new candidate sample using peeling and cutting (if necessary)  
	*/
	
	
	void Population::getSampleForLocusMIM(unsigned j){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		std::string typeOfGNode = MIM::parmMap["typeOfGNodes"];
		if(typeOfGNode=="genotypic"){
			initGenotypeNodeList(j);
		} 
		else if(typeOfGNode=="allelic"){
			initAlleleNodeList(j);
		}
		else {
			cerr << "Sampler type should be either genotypic or allelic;" << endl;
			exit(1);
		}
		gNodeList.peelCutAndSampleMIM(); 
		setSegregationIndex(j,typeOfGNode);
	}
	/*! \fn void Population::getSampleForLocusMIM(unsigned j)
	*  \brief get a sample for a given locus
	*/
	
	void Population::getOldGNodeListProbabilityMIM(){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		unsigned printFlag  = MIM::parmMap.getUnsignedValue("printFlag"); 
		unsigned startLocus = MIM::parmMap.getUnsignedValue("startLocus") - 1;
		// get a sample for the startLocus
		lookToYourLeft  = false;
		lookToYourRight = false;
		Individual::currentPosition = startLocus;
		unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
		getOldProbabilityForLocusMIM(j);
		// get samples for the other loci 
		// first to the right
		for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
			lookToYourLeft  = true;
			lookToYourRight = false;
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
			getOldProbabilityForLocusMIM(j);
		}
		// next to the left
		for (int jj=startLocus-1;jj>=0;jj--){
			lookToYourLeft  = false;
			lookToYourRight = true;
			Individual::currentPosition = jj;
			unsigned j = prior->chrom()[0].locusMask[jj].index; 
			Individual::currentLocus = j;
			if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
			getOldProbabilityForLocusMIM(j);
		}
	}	
	/*! \fn void Population::getOldGNodeListProbability(void)
	*  \brief recompute the needed quantities for the existing (old)
	sample (for the MH step)
	*/ 
	
	void Population::getOldProbabilityForLocusMIM(unsigned j){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		std::string typeOfGNode = MIM::parmMap["typeOfGNodes"];	
		if (typeOfGNode=="genotypic"){
			initGenotypeNodeList(j);
		} 
		else if (typeOfGNode=="allelic"){
			initAlleleNodeList(j);
		}
		else {
			cerr << "Sampler type should be either genotypic or allelic;" << endl;
			exit(1);
		}
		gNodeList.peelCutAndComputeMIM(); 
		setSegregationIndex(j,typeOfGNode);
	}
	/*! \fn void Population::getOldProbabilityForLocus(unsigned j)
	*  \brief calculate the probability of an old genotype configuration given the new probability structures used to obtain 
	*         the new sample.
	*/
	
	void Population::initCounters(void){
    // Authors: Rohan L. Fernando 
    // (September, 2005) 
    // Contributors: 
    string whatToCompute = MIM::parmMap["whatToCompute"];
    if (whatToCompute=="genotypeFreq") {
        initGenotypeFreq();
    }
	else if (whatToCompute=="Q2Probs") {
        initGenotypeFreq();
    }
    else if (whatToCompute=="haplotypeFreq"){
        // all that is needed is done in updateCounters()  
    }
    else if (whatToCompute=="ibdCovMatrix") {
        MIM::IBDCovMatrix.resize(popsize,popsize,0.0);
        MIM::meansVector.resize(popsize,1,0.0);
    }
    else if (whatToCompute=="probDescHaplo"){
        SetFreqHaploFounders();
    }
    else if (whatToCompute=="recessiveLocusLocation"){
        identifyLociToSamplePosFor();
    }
	else if (whatToCompute=="ibdMatrix"){
		unsigned leftLocus  = MIM::parmMap.getUnsignedValue("leftLocus");
		unsigned rightLocus = MIM::parmMap.getUnsignedValue("rightLocus");
		MIM::ibdMatrix.initialize(leftLocus,rightLocus);
	}
	else if(whatToCompute=="parentOffspringProbs"){
	    initGenotypeFreq();
		initParentOffspring();
	}
}
	
void Population::updateCounters(){
    // Authors: Rohan L. Fernando 
    // (September, 2005) 
    // Contributors: L. Radu Totir
    std::string whatToCompute = MIM::parmMap["whatToCompute"];
    std::string typeOfGNode   = MIM::parmMap["typeOfGNodes"];
	std::string typeOfSampler = MIM::parmMap["GNodeSampler"];
    if (whatToCompute=="genotypeFreq"){
        countGenotypes(typeOfGNode);
    }
    else if (whatToCompute=="haplotypeFreq"){
        countHaplotypes(typeOfGNode);
    } 
    else if (whatToCompute=="ibdCovMatrix"){
        MIM::IBDCovMatrix += getIBDMatrix();
        MIM::meansVector  += getMeans();
    }
    else if(whatToCompute=="probDescHaplo"){
        UpdateFreqHaploFounders();
        sampleSegregationIndicators();
    }
    else if (whatToCompute=="recessiveLocusLocation"){
        samplePositionOnChromosome();
    }
	else if (whatToCompute=="ibdMatrix"){
		if (typeOfSampler=="JNC" || typeOfSampler=="LG" || typeOfSampler=="JointNC"){
			MIM::ibdMatrix.sample();
		}
		else {
			sampleSegregationIndicators();
			MIM::ibdMatrix.sample1();
		}
	}
	else if (whatToCompute=="parentOffspringProbs"){
		countParentOffspring(typeOfGNode);
	}

}

void Population::updateCountersSimple(){
    // Authors: Rohan L. Fernando 
    // (September, 2005) 
    // Contributors: L. Radu Totir
    std::string whatToCompute = MIM::parmMap["whatToCompute"];
    std::string typeOfGNode   = MIM::parmMap["typeOfGNodes"];
	std::string typeOfSampler = MIM::parmMap["GNodeSampler"];
    if (whatToCompute=="genotypeFreq"){
        countGenotypesSimple(typeOfGNode);
    }
    else if (whatToCompute=="haplotypeFreq"){
        countHaplotypesSimple(typeOfGNode);
    } 
    else if (whatToCompute=="ibdCovMatrix"){
        MIM::IBDCovMatrix += getIBDMatrixSimple();
        MIM::meansVector  += getMeansSimple();
    }
    else if(whatToCompute=="probDescHaplo"){
        UpdateFreqHaploFoundersSimple();
        sampleSegregationIndicatorsSimple();
    }
    else if (whatToCompute=="recessiveLocusLocation"){
        samplePositionOnChromosome();
    }
	else if (whatToCompute=="ibdMatrix"){
		if (typeOfSampler=="JNC" || typeOfSampler=="LG" || typeOfSampler=="JointNC"){
			MIM::ibdMatrix.sampleSimple();
		}
		else {
			sampleSegregationIndicators();
			MIM::ibdMatrix.sample1Simple();
		}
	}
	else if (whatToCompute=="parentOffspringProbs"){
		countParentOffspringSimple(typeOfGNode);
	}

}
	
	
void Population::outputCounters(unsigned numOfSamples){
    // Authors: Rohan L. Fernando 
    // (September, 2005) 
    // Contributors: L. Radu Totir
    string whatToCompute  = MIM::parmMap["whatToCompute"];
    string resultsFile    = MIM::parmMap["resultsFile"];
    // unsigned numOfSamples = MIM::parmMap.getUnsignedValue("chainLength");
    if (whatToCompute=="genotypeFreq"){
        if(resultsFile == ""){
            displayGenotypeFrequencies(numOfSamples);
			displayPDQ("pdq.res");
        }
        else {
            std::ofstream outfile(resultsFile.c_str());
            std::ofstream outfile1((resultsFile+".1").c_str());
            displayGenotypeFrequencies(numOfSamples, outfile);
            //outputGenotypeFrequenciesSL(numOfSamples, outfile1);
			displayPDQ("pdq.res");
        }
    }
	else if (whatToCompute=="Q2Probs"){
        if(resultsFile == ""){
            outputQ2Probs(numOfSamples);
			displayPDQ("pdq.res");
        }
        else {
            std::ofstream outfile(resultsFile.c_str());
			outputQ2Probs(numOfSamples, outfile);
			displayPDQ("pdq.res");
        }
    }
    else if (whatToCompute=="haplotypeFreq"){
        if(resultsFile == ""){
            displayHaplotypeFrequencies(numOfSamples);
        }
        else {
            std::ofstream outfile(resultsFile.c_str());
            displayHaplotypeFrequencies(numOfSamples, outfile);
        }
    }
    else if(whatToCompute=="ibdCovMatrix"){
        matvec::doubleMatrix temp1 = MIM::IBDCovMatrix/numOfSamples;
        matvec::doubleMatrix temp2 = MIM::meansVector /numOfSamples;
        if(resultsFile == ""){
            std::cout << temp1;
            std::cout << endl;
            std::cout << temp2;
            std::cout << endl;
            std::cout << temp1 - temp2*temp2.transpose();
        }
        else {
            std::ofstream outfile(resultsFile.c_str());
            outfile << temp1;
            outfile << endl;
            outfile << temp2;
            outfile << endl;
            outfile << temp1 - temp2*temp2.transpose();
        } 
    }
    else if(whatToCompute=="probDescHaplo"){
        CalcFreqHaploFounders(numOfSamples);
        DisplayFreqHaploFounders();
    }
    else if (whatToCompute=="recessiveLocusLocation"){
        
    }
	else if (whatToCompute=="ibdMatrix"){
	    if(resultsFile == ""){
			MIM::ibdMatrix.display();
		}
		else {
			MIM::ibdMatrix.output(resultsFile);
		}		
    }
	else if(whatToCompute=="parentOffspringProbs"){
		std::ofstream outfile(resultsFile.c_str());
		displayParentOffspringProbs(outfile);
	}
}

	
	matvec::Matrix<double> Population::getIBDMatrix(void){
  // Author: L. Radu Totir
  // (July, 2004) 
  // Contributors: Rohan L. Fernando
  matvec::Matrix<double> geneticValVec;
  matvec::Matrix<double> CovMatrix;
  CovMatrix.resize(popsize,popsize,0.0);
  matvec::Vector<double> qtlMu;
  qtlMu.resize(4,0.0);
  qtlMu[0]=-1.41421,qtlMu[3]=1.41421;
  geneticValVec.resize(popsize,1,0.0);
  for (unsigned t=0;t<popsize;t++){
	unsigned j = 1;
	int mat = (popmember[t])->genotNodeVector[j].getAcceptedMatState();
	int pat = (popmember[t])->genotNodeVector[j].getAcceptedPatState();
	mat = (mat>1) ? 1 : mat;
	pat = (pat>1) ? 1 : pat;
	int geno = mat*2 + pat;
    geneticValVec[t][0] = qtlMu[geno];
  }
  CovMatrix =geneticValVec*geneticValVec.transpose();
  return CovMatrix;
}

	matvec::Matrix<double> Population::getIBDMatrixSimple(void){
  // Author: L. Radu Totir
  // (July, 2004) 
  // Contributors: Rohan L. Fernando
  matvec::Matrix<double> geneticValVec;
  matvec::Matrix<double> CovMatrix;
  CovMatrix.resize(popsize,popsize,0.0);
  matvec::Vector<double> qtlMu;
  qtlMu.resize(4,0.0);
  qtlMu[0]=-1.41421,qtlMu[3]=1.41421;
  geneticValVec.resize(popsize,1,0.0);
  for (unsigned t=0;t<popsize;t++){
	unsigned j = 1;
	int mat = (popmember[t])->genotNodeVector[j].getmState();
	int pat = (popmember[t])->genotNodeVector[j].getpState();
	mat = (mat>1) ? 1 : mat;
	pat = (pat>1) ? 1 : pat;
	int geno = mat*2 + pat;
    geneticValVec[t][0] = qtlMu[geno];
  }
  CovMatrix =geneticValVec*geneticValVec.transpose();
  return CovMatrix;
}


matvec::Matrix<double> Population::getMeans(void){
  // Author: Rohan L. Fernando
  // (September, 2005) 
  // Contributors: 
  matvec::Matrix<double> geneticValVec;
  matvec::Vector<double> qtlMu;
  qtlMu.resize(4,0.0);
  qtlMu[0]=-1.41421,qtlMu[3]=1.41421;
  geneticValVec.resize(popsize,1,0.0);
  for (unsigned t=0;t<popsize;t++){
	unsigned j = 1;
	int mat = (popmember[t])->genotNodeVector[j].getAcceptedMatState();
	int pat = (popmember[t])->genotNodeVector[j].getAcceptedPatState();
	mat = (mat>1) ? 1 : mat;
	pat = (pat>1) ? 1 : pat;
	int geno = mat*2 + pat;
    geneticValVec[t][0] = qtlMu[geno];
  }
  return geneticValVec;
}

matvec::Matrix<double> Population::getMeansSimple(void){
  // Author: Rohan L. Fernando
  // (September, 2005) 
  // Contributors: 
  matvec::Matrix<double> geneticValVec;
  matvec::Vector<double> qtlMu;
  qtlMu.resize(4,0.0);
  qtlMu[0]=-1.41421,qtlMu[3]=1.41421;
  geneticValVec.resize(popsize,1,0.0);
  for (unsigned t=0;t<popsize;t++){
	unsigned j = 1;
	int mat = (popmember[t])->genotNodeVector[j].getmState();
	int pat = (popmember[t])->genotNodeVector[j].getpState();
	mat = (mat>1) ? 1 : mat;
	pat = (pat>1) ? 1 : pat;
	int geno = mat*2 + pat;
    geneticValVec[t][0] = qtlMu[geno];
  }
  return geneticValVec;
}


void Population::setMissingToHet(unsigned locus){
	// Authors: Rohan L. Fernando and  L. Radu Totir
    // (October, 2005) 
    // Contributors: 
	unsigned nAlleles = prior->chrom()[0].locus[locus].nallele();
	if (nAlleles > 2) {
		throw exception ("Population::setMissingToHet: called with nAlleles > 2 \n");
	}
	Individual *ind;
	unsigned hetInd[2],index, setGenotype;
	for(unsigned i=0; i<popsize; i++){
		ind=popmember[i];
		GenotypeNode* indGenotype = &ind->genotNodeVector[locus];
		unsigned nGenotypes = indGenotype->genotypeVector.size();
		if (nGenotypes > 1){
			index = 0;
			for (unsigned j=0;j<nGenotypes;j++){
				indGenotype->reset(j);
				int matState = indGenotype->getmState();
				int patState = indGenotype->getpState();
				if (matState != patState){
					hetInd[index] = j;
					index++;
				}
			}
			if(index==1){
				setGenotype = hetInd[0];
			}
			else{
				double u = ranf();
				index = (u<0.5) ? 0 : 1;
				setGenotype = hetInd[index];
			}
			indGenotype->reset(setGenotype);
		}
		else {
			indGenotype->reset(0);
		}
		indGenotype->sampled = true;
	}
}

GNodeSet Population::getSirePivots(unsigned locus){
	// Authors: Rohan L. Fernando and  L. Radu Totir
    // (October, 2005) 
    // Contributors: 
	Individual *ind, *dad;
	GNodeSet pivotSet;
	for(unsigned i=0; i<popsize; i++){
		ind=popmember[i];
		dad=ind->myfather;
		if (dad){
			GenotypeNode* dadGenotype = &dad->genotNodeVector[locus];
			dadGenotype->uninformativeGNode = 0; // don't want to prune pivots
			pivotSet.insert(dadGenotype);
		}
	}
	return pivotSet;
}


} ////////// end of namespace matvec



