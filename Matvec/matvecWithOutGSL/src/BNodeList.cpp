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

//#include <algorithm>
//#include <cstdarg>
//#include <cstdlib>
//#include <cmath>
//#include <fstream>
//#include <iomanip>
#include <iostream>
#include <set>
#include "stat.h"
#include "gnodederived.h"
#include "gnodesetderived.h"
#include "gnodestuff.h"
#include "population.h"
//#include "model.h"
#include "mim.h"
#include "BNodeList.h"


namespace matvec {

	void BlockNodeList::initNodes(GNodeList &blockList){
		// Authors: Christian Stricker and Rohan L. Fernando 
		// (June, 2006) 
		// Contributors: 
		BlockNode blockNode;
		set<GNodeSet*>::iterator gnsetit;
		clear();
		for (unsigned i=0;i<blockList.size();i++){
			blockNode.myNode = blockList[i];
			blockNode.generatedSet = blockList[i]->generatedSet;
			blockNode.neighborSet = blockList[i]->neighborSet;
			blockNode.vectorGNSets.clear();
			for(gnsetit = blockList[i]->SetofGNsts.begin(); gnsetit != blockList[i]->SetofGNsts.end(); gnsetit++){
				blockNode.vectorGNSets.push_back(*gnsetit);
			}
			push_back(blockNode);
		}
	}
	/*! \fn void BlockNodeList::initNodes(GNodeList &blockList)
	*   after peeling in GNodeList the GNodeSets that are actually used in the peeling
	*   are in the SetofGNSets. These are transfered here to the vectorGNSets of a BlockNode.
	*   
	*/ 
	
	void BlockNodeList::peelNoCutNoOrder(void){
		// Authors: Christian Stricker and Rohan L. Fernando 
		// (June, 2006) 
		// Contributors: 
		unsigned printFlag, count=0;
		printFlag = MIM::parmMap.getUnsignedValue("printFlag");
		BlockNodeList::iterator it;
		
		for (it=begin();it!=end();it++){
			count++;
			if(it->myNode->getWeight() > 1) {
				if(printFlag>1) cout << "Peeling node " << count<<"("<<it->myNode->id<<")"<< endl;
				it->peel();
			}
			else {
				if(printFlag>1) cout << "Skipping node " << count<<"("<<it->myNode->id<<")"<< endl;
			}
			if(printFlag>2) {
				displayGNodeSets();
			}
		}
	}
	/*! \fn void GNodeList::peelAndCutNoOrder()
	*  \brief peels and cuts (if it needs to) according to peeling order established in vector of block nodes
	*/ 
	
	void BlockNodeList::displayGNodeSets(void){
		// Authors: Christian Stricker and Rohan L. Fernando 
		// (June, 2006) 
		// Contributors:
		BlockNodeList::iterator it;
		GNodeSet::iterator setit;
		for (it=begin();it!=end();it++){
			std::cout << it->myNode->id << "   ";
			for (unsigned i=0;i<it->vectorGNSets.size();i++){
				it->vectorGNSets[i]->display();
			}
			std::cout << endl;
		}
		std::cout << "=========== " << endl;
	}
	/*! \fn void BlockNodeList::displayGNodeSets()
	*  \brief displays the GNodeSets of each Node (pentrance, founder,
												   transition sets)
	*/
	
	
	void BlockNode::peel(void) {
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (June, 2006) 
		// Contributors: chris stricker
		*neighborSet = vectorGNSets[0];
		for (unsigned i=1; i<vectorGNSets.size(); i++){
			*neighborSet *= vectorGNSets[i];
		}
		*generatedSet+= neighborSet;
	}
	/*! \fn void BlockNode::peel(void)
	*  \brief the BlockNode peeling method without releasing cutsets and without updatemysets() (peeling the GNode in question)
	*/
	
void BlockNode::reverseSampleGNode(void) {
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2006) 
  // Contributors: Chris Stricker
	
	try{
		if(myNode->myGNodeProbs->size()==0){
			myNode->myGNodeProbs->insert(myNode);
			myNode->myGNodeProbs->setKeyMultiplicationCode();
		}
		
		*neighborSet = vectorGNSets[0];
		for (unsigned i=1; i<vectorGNSets.size(); i++){
			*neighborSet *= vectorGNSets[i];
		}
		*(myNode->myGNodeProbs) += neighborSet;

		myNode->myGNodeProbs->normalize();
		double u = ranf();
		double sum = 0.0;
		myNode->myGNodeProbs->reset();
		unsigned sampledState=0;
		do{
			double tempProb = myNode->myGNodeProbs->getValue();
			sum+= tempProb;
			if (u<sum){
				myNode->reset(sampledState);
				myNode->sampled = true;
				GNodeList::logProposal += std::log(tempProb);
				return;
			}
			sampledState++;
		} while(myNode->incr());
	}
	catch (...){
		cout << "in reverse sample\n";
	}
}
/*! \fn void BlockNode::reverseSampleGNode(void)
 *  \brief the GNode sampling method in reverse order using Blocknodes (sample the GNode in question)
*/

	
	
	void BlockNodeList::setGNodeSampleFlags(bool state){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (August, 2003) 
		// Contributors: 
		BlockNodeList::iterator it;
		for (it=begin();it!=end();it++){ 
			it->myNode->sampled = state;
		}
	}
	/*! \fn void BlockNodeList::setGNodeSampleFlags()
	*  \brief reset sampled flag to false
	*/
	
	void BlockNodeList::reverseSample(void){
		// Authors: Rohan L. Fernando 
		// (October, 2005) 
		// Contributors: 
		BlockNodeList::reverse_iterator it;
		for (it=rbegin();it!=rend();it++){
			if(!(it->myNode->sampled)){
				if(it->myNode->getWeight() == 1){
					it->myNode->sampled = true;
				}
				else {
					it->reverseSampleGNode();
				}
			}
		}	
	}
	
	void BlockNodeList::sampleSegregationIndicators(void){
		// Authors: Rohan L. Fernando 
		// (October, 2005) 
		// Contributors: 
		BlockNodeList::iterator it;
		for (it=begin();it!=end();it++){
			Individual *ind = it->myNode->owner;
			for (unsigned locus=0;locus<Individual::numLoci;locus++) ind->setSegregationIndex(locus,"genotypic");
			ind->sampleSegregationIndicators();
		}	
	}
	
	void BlockNodeList::saveGNodeStates(void){
		// Authors: Rohan L. Fernando 
		// (November, 2005) 
		// Contributors: 
		BlockNodeList::iterator it;
		for (it=begin();it!=end();it++){
			it->myNode->saveState();
		}
	}	

void BlockNodeList::updateCounters(void){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	BlockNodeList::iterator it;
	for(it=begin();it!=end();it++){
		it->myNode->updateCounter();
	}
}
/*! \fn void BlockNodeList::updateCounters(void)
*  \brief update the GNode counters 
*/
	
}
