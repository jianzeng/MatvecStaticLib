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

#include "mqtl.h"
#include "mim.h"
#include "rpedigree.h"
#include "session.h"
using namespace std; 

namespace matvec{
	RPedigree* MQTL::myPedPtr=0;
	
	void MQTL::inputPDQ(char *filename){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		unsigned rec=0, rec1=0;
		QNode *ptr;
		double pPDQ, mPDQ;
		string indStr;
		ifstream pdqfile;
		pdqfile.open(filename);
		if(!pdqfile){
			cerr << "Couldn't open " << filename << endl;
			exit (-1);
		}
		std::cout << "\n reading PDQ file: " << filename << std::endl;
		resize(myPedPtr->size(),0);
		while (pdqfile >> indStr >> pPDQ >> mPDQ ){

			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";				
				cout.flush();
				rec1 += rec;
				rec = 0;
			}
			RPedigree::iterator pedIt = myPedPtr->find(indStr);
			if (pedIt == myPedPtr->end()) {
				cout << indStr << " in record " << (rec1+rec) 
				<< " of " << filename  
				<<" not found in pedigree file \n";
				exit(1);
			}
			ptr = new QNode();
			ptr->pPDQ = pPDQ;
			ptr->mPDQ = mPDQ;
			(*this)[pedIt->second->ind-1] = ptr;
//			cout<<indStr<<"\t"<<pedIt->second->ind-1<<"\t"<<ptr->mPDQ<<"\t"<<ptr->pPDQ<<endl;
		}
		pdqfile.close();
		for (unsigned i=0;i<size();i++){
			if((*this)[i]==0){
				if (myPedPtr->pedVector[i]->sire==0 && myPedPtr->pedVector[i]->dam==0) {
					QNode *ptr = new QNode();
					ptr->pPDQ = -1.0;
					ptr->mPDQ = -1.0;
					(*this)[i] = ptr;
				}
				else {
					cout << myPedPtr->pedVector[i]->ind_str << " not in PDQ file \n";
					exit(-1);
				}
			}
		}
	}
	void MQTL::inputQTLProbs(char *filename){
		// Authors: Rohan L. Fernando
		// (December, 2005) 
		// Contributors:
		unsigned rec=0, rec1=0, count=0;
		QNode *ptr;
		double prQ11, prQ12, prQ21, prQ22;
		string indStr, scratch, headerLine;
		
		ifstream prfile;
		prfile.open(filename);
		if(!prfile){
			cerr << "Couldn't open " << filename << endl;
			exit (-1);
		}
		std::cout << "\n reading QTLProbs file: " << filename << std::endl;
		while(prfile.peek() != 10 && prfile.peek() != 13){ //first line needs to be a header line with a string labelling each column, each string shall NOT contain blanks
			prfile >>scratch; // here we count how many columns there are, 
			count++;
		}
		while((prfile.peek() != EOF) && (prfile.peek() == 10 || prfile.peek() == 13)) prfile.ignore(1); 
		while (prfile >> indStr >> prQ11 >> prQ12 >> prQ21 >> prQ22 ){
			for(unsigned i=0; i<count-5;i++){ //skip all the columns after the fifth
				prfile>>scratch;
			}
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}
			RPedigree::iterator pedIt = myPedPtr->find(indStr);
			if (pedIt == myPedPtr->end()) {
				cout << indStr << " in record " << (rec1+rec) 
				<< " of " << filename  
				<<" not found in pedigree file \n";
				exit(1);
			}
			(*this)[pedIt->second->ind-1]->prQ11 = prQ11;
			(*this)[pedIt->second->ind-1]->prQ12 = prQ12;
			(*this)[pedIt->second->ind-1]->prQ21 = prQ21;
			(*this)[pedIt->second->ind-1]->prQ22 = prQ22;
			(*this)[pedIt->second->ind-1]->patQ2Prob = prQ12 + prQ22;
			(*this)[pedIt->second->ind-1]->matQ2Prob = prQ21 + prQ22;			
		}
		prfile.close();
	}
	
	void MQTL::calcQTLVars(double mu){
	    std::cout << "\n calculating variances of QTL allele effects"<< std::endl;
	    unsigned rec=0, rec1=0;
		double patQ2Prob, matQ2Prob;
		for (unsigned qi=0; qi<size(); qi++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}		
			QNode *ptr = ((*this)[qi]);
			if(ptr->mLevel != pseudoLevel) {
				matQ2Prob = ptr->matQ2Prob;
				(*this)[qi]->varvj = mu*mu*matQ2Prob*(1-matQ2Prob);
			}
			else {
				(*this)[qi]->varvj = 0.000000000001;
			}
			if(ptr->pLevel != pseudoLevel) {
				patQ2Prob = ptr->patQ2Prob;
				(*this)[qi]->varvi = mu*mu*patQ2Prob*(1-patQ2Prob);
			}
			else {
				(*this)[qi]->varvi = 0.000000000001;
			}
			//cout << (*this)[qi]->varvi << "    " << (*this)[qi]->varvj << endl;
		}
	}
	
	void MQTL::generateMQTLLevels(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		std::cout << "\n generating MQTL levels"<< std::endl;
		unsigned rec=0, rec1=0;
		MQTL::iterator it;
		count = 0;
		RPedigree::iterator pedIt;

		for (unsigned i=0;i<size();i++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}		
			if ((*this)[i]->mLevel==0) {
				codeMQTL(i);
			}
		}
		codeLevels();
	}
	
		void MQTL::generateQTLLevels(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		std::cout << "\n generating QTL levels"<< std::endl;
		unsigned rec=0, rec1=0;
		MQTL::iterator it;
		count = 0;
		RPedigree::iterator pedIt;
		for (unsigned i=0;i<size();i++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}		
			if ((*this)[i]->mLevel==0) {
				codeQTL(i);
			}
		}
		codeLevels();
	}

	
	void MQTL::codeMQTL(unsigned i){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		if ((*this)[i]->mLevel) return;       // already coded
		int dam   = myPedPtr->pedVector[i]->dam;
		int sire  = myPedPtr->pedVector[i]->sire;
		
		if (sire==0) {
			(*this)[i]->pLevel = ++count;
		}
		else {
			codeMQTL(sire-1);
			if((*this)[i]->pPDQ==1){
				(*this)[i]->pLevel = (*this)[sire-1]->mLevel;
			}
			else if ((*this)[i]->pPDQ==0){
				(*this)[i]->pLevel = (*this)[sire-1]->pLevel;
			}
			else {
				(*this)[i]->pLevel = ++count;
			}
		}
		
		if (dam==0) {
			(*this)[i]->mLevel = ++count;
		}
		else {
			codeMQTL(dam-1);
			if((*this)[i]->mPDQ==1){
				(*this)[i]->mLevel = (*this)[dam-1]->mLevel;
			}
			else if ((*this)[i]->mPDQ==0){
				(*this)[i]->mLevel = (*this)[dam-1]->pLevel;
			}
			else {
				(*this)[i]->mLevel = ++count;
			}
		}
	}
	
		void MQTL::codeQTL(unsigned i){
		// Authors: Rohan L. Fernando
		// (December, 2005) 
		// Contributors:
		// Note that pedVector and the MQTL vector 
		// are ordered such that parents precede offspring
		
		QNode *QNodePtr ((*this)[i]);
		
		bool observedQG = false;
		observedQG = (QNodePtr->prQ11 > 0.99)                   ? true : observedQG;
		observedQG = ((QNodePtr->prQ12+QNodePtr->prQ21) > 0.99) ? true : observedQG;
		observedQG = (QNodePtr->prQ22 > 0.99)                   ? true : observedQG;
		if (observedQG) {
			QNodePtr->mPDQ = -1.0;
			QNodePtr->pPDQ = -1.0;
			if (pseudoLevel == 0){
				pseudoLevel = ++count;
			}
			QNodePtr->mLevel = pseudoLevel;
			QNodePtr->pLevel = pseudoLevel;				
			return;
		}	
		int dam   = myPedPtr->pedVector[i]->dam;
		int sire  = myPedPtr->pedVector[i]->sire;
		double Q2Prob = QNodePtr->prQ12 + QNodePtr->prQ22;
		if (Q2Prob > .99 || Q2Prob < .01) {
		    if (pseudoLevel == 0){
				pseudoLevel = ++count;
			}
			QNodePtr->pLevel = pseudoLevel;
			QNodePtr->pPDQ = -1.0;
		}
		else if (sire==0) {
			QNodePtr->pLevel = ++count;
		}
		else if ((*this)[sire-1]->mLevel==0) {
			QNodePtr->pLevel = ++count;
		}
		else {
			if(QNodePtr->pPDQ > .99){
				QNodePtr->pLevel = (*this)[sire-1]->mLevel;
			}
			else if (QNodePtr->pPDQ<.01){
				QNodePtr->pLevel = (*this)[sire-1]->pLevel;
			}
			else {
				QNodePtr->pLevel = ++count;
			}
		}
		Q2Prob = QNodePtr->prQ21 + QNodePtr->prQ22;
		if (Q2Prob > .99 || Q2Prob < .01) {
		    if (pseudoLevel == 0){
				pseudoLevel = ++count;
			}		
			QNodePtr->mLevel = pseudoLevel;
			QNodePtr->mPDQ = -1.0;
		}
		else if (dam==0) {
			QNodePtr->mLevel = ++count;
		}
		else if ((*this)[dam-1]->mLevel==0) {
			QNodePtr->mLevel = ++count;
		}
		else {
			if(QNodePtr->mPDQ>.99){
				QNodePtr->mLevel = (*this)[dam-1]->mLevel;
			}
			else if (QNodePtr->mPDQ<.01){
				QNodePtr->mLevel = (*this)[dam-1]->pLevel;
			}
			else {
				QNodePtr->mLevel = ++count;
			}
		}
	}
	
	void MQTL::codeLevels(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		for (unsigned i=1;i<=count;i++){
			string str = getString(i);
			myRecoder.code(str);
		}
	}
	
	string MQTL::getString(unsigned i){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		ostringstream outputStrStream(ostringstream::out);
		outputStrStream << i; 
		return outputStrStream.str();
	}
	
	void MQTL::display(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		for (unsigned i=0;i<size();i++){
			cout << setw(10) <<  myPedPtr->pedVector[i]->ind_str
			<< setw(5)  << (*this)[i]->pLevel 
			<< setw(5)  << (*this)[i]->mLevel 
			<< setw(10)  << (*this)[i]->pPDQ 
			<< setw(10)  << (*this)[i]->mPDQ 
			<< setw(10) << (*this)[i]->covvivj   
			<< endl;
		}
	}
	void MQTL::calcQTLcov(void){
	// It is assumed that in the pedgree file ancestors are 
	// listed before descendants
	std::cout << "\n calculating covariances between maternal and paternal QTL allele effects"<< std::endl;
	SpCij.clear();
	unsigned rec=0, rec1=0;
	for (unsigned qi=0; qi<size(); qi++){
		rec++;
		if(rec==1000){
			cout<<rec+rec1<<"\r";
			cout.flush();
			rec1 += rec;
			rec = 0;
		}	
		if(myPedPtr->pedVector[qi]->sire==0||myPedPtr->pedVector[qi]->dam==0){
			(*this)[qi]->covvivj = 0.0;
		}
		else{
			unsigned ind = qi + 1;
			(*this)[qi]->covvivj = getQTLcovvivj(ind, ind, 0, 1);
		}
	}
}

double MQTL::getQTLcovvivj(unsigned i, unsigned j, unsigned pi, unsigned pj){
	if (i==0 || j==0){
		return 0.0;
	}
	unsigned ai = pi ? (*this)[i-1]->pLevel : (*this)[i-1]->mLevel;
	unsigned aj = pj ? (*this)[j-1]->pLevel : (*this)[j-1]->mLevel;
	if (ai==pseudoLevel || aj==pseudoLevel) {
		return 0.0;
	}
	float x = SpCij.retrieve_cij(ai,aj);
	if(x != -1.0) {
		return x;
	}
	unsigned oldp,youngp;
	unsigned old, young;
	if(i < j){
		old = i;
		oldp = pi;
		young = j;
		youngp = pj;
	}
	else if(j < i){
		old = j;
		oldp = pj;
		young = i;
		youngp = pi;
	}
	else if(pi==pj) {
	        if(pi) return (*this)[i-1]->varvi;
	        else   return (*this)[i-1]->varvj;
	}
	else {
		old = i;
		oldp = pi;
		young = j;
		youngp = pj;
	}
	double PDQym, PDQyp;
	unsigned y_parent;
	if (youngp==1){
		y_parent = myPedPtr->pedVector[young - 1]->sire;
		if(y_parent==0) return 0.0;
		PDQym = (*this)[young - 1]->pPDQ;
	}
	else{
		y_parent = myPedPtr->pedVector[young -1]->dam;
		if(y_parent==0) return 0.0;
		PDQym    = (*this)[young - 1]->mPDQ;
	}
	PDQyp = 1.0 - PDQym;
	x = PDQym*getQTLcovvivj(old,y_parent,oldp,0) + PDQyp*getQTLcovvivj(old,y_parent,oldp,1);
	SpCij.put_cij(ai,aj,x);
	return x;
}

	void MQTL::calcMQTLInbreeding(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		std::cout << "\n calculating MQTL inbreeding " << std::endl;
		SpCij.clear();
		unsigned rec=0, rec1=0;
		for (unsigned i=0;i<size();i++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}		
			if(myPedPtr->pedVector[i]->sire==0 || myPedPtr->pedVector[i]->dam==0){
				(*this)[i]->f= 0.0;
			}
			else{
				unsigned ind = i+1;
				(*this)[i]->f = getMQTLrij(ind,ind,0,1);
			}
		}
	}
	
	double MQTL::getMQTLrij(unsigned i, unsigned j, unsigned pi, unsigned pj){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		if (i==0 || j==0){
			return 0.0;
		}
		unsigned vi = pi ? (*this)[i-1]->pLevel : (*this)[i-1]->mLevel;
		unsigned vj = pj ? (*this)[j-1]->pLevel : (*this)[j-1]->mLevel;
		float x = SpCij.retrieve_cij(vi,vj);
		if(x != -1.0) {
			return x;
		}
		unsigned old, young,oldp,youngp;
		if(i < j){
			old = i;
			oldp = pi;
			young = j;
			youngp = pj;
		}
		else if(j < i){
			old = j;
			oldp = pj;
			young = i;
			youngp = pi;
		}
		else if(pi==pj) {
			return 1.0;
		}
		else {
			old = i;
			oldp = pi;
			young = j;
			youngp = pj;
		}
		unsigned y_parent;
		double PDQym, PDQyp;
		if (youngp==1){
			y_parent = myPedPtr->pedVector[young-1]->sire;
			PDQym = (*this)[young-1]->pPDQ;
		}
		else{
			y_parent = myPedPtr->pedVector[young-1]->dam;
			PDQym    = (*this)[young-1]->mPDQ;
		}
		PDQyp = 1.0 - PDQym;
		x = PDQym*getMQTLrij(old,y_parent,oldp,0) + PDQyp*getMQTLrij(old,y_parent,oldp,1);
		SpCij.put_cij(vi,vj,x);
		return x;
	}
	
	void MQTL::addGinv(matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol, double ratio){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		MQTL::iterator it;
		unsigned maxLevel=0;
		double q[3];
		unsigned pos[3];
		double v;
		for (unsigned i=0;i<size();i++){
			if ((*this)[i]->pLevel > maxLevel){
				if (myPedPtr->pedVector[i]->sire==0){ // founder
					pos[0] = 0;
					pos[1] = 0;
					pos[2] = (*this)[i]->pLevel;
					q[0]   = 0;
					q[1]   = 0;
					q[2]   = 1.0;
					v = 1.0;
					addtoGinv(pos,q,v,lhs,startRow,startCol,ratio);
				}
				else { // not founder
					unsigned dad = myPedPtr->pedVector[i]->sire;
					pos[0] = (*this)[dad-1]->mLevel;
					pos[1] = (*this)[dad-1]->pLevel;
					pos[2] = (*this)[i]->pLevel;
					q[0]   = -(*this)[i]->pPDQ;  
					q[1]   = -(1+q[0]);
					q[2]   = 1.0;    
					v      = 1 - q[0]*q[0] - q[1]*q[1] - 2*q[0]*q[1]*(*this)[dad-1]->f;  
					addtoGinv(pos,q,1.0/v,lhs,startRow,startCol,ratio);
				}
				maxLevel = pos[2]; 
			}
			if ((*this)[i]->mLevel > maxLevel){
				if (myPedPtr->pedVector[i]->dam==0) { // founder
					pos[0] = 0;
					pos[1] = 0;
					pos[2] = (*this)[i]->mLevel;
					q[0]   = 0;
					q[1]   = 0;
					q[2]   = 1.0;
					v = 1.0;
					addtoGinv(pos,q,v,lhs,startRow,startCol,ratio);
				}
				else { // not founder
					unsigned mom = myPedPtr->pedVector[i]->dam;
					pos[0] = (*this)[mom-1]->mLevel;
					pos[1] = (*this)[mom-1]->pLevel;
					pos[2] = (*this)[i]->mLevel;
					q[0]   = -(*this)[i]->mPDQ;  
					q[1]   = -(1+q[0]);
					q[2]   = 1.0;    
					v      = 1 - q[0]*q[0] - q[1]*q[1] - 2*q[0]*q[1]*(*this)[mom-1]->f;  
					addtoGinv(pos,q,1.0/v,lhs,startRow,startCol,ratio);
				}
				maxLevel = pos[2];
			}
		}
	}
	
	void MQTL::addGinvDis(matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol){
	MQTL::iterator it;
	unsigned maxLevel=0;
	double q[3];
	unsigned pos[3];
	double v;
	for (unsigned i=0;i<size();i++){
		if ((*this)[i]->pLevel > maxLevel){
			unsigned mySire = myPedPtr->pedVector[i]->sire;
			bool founder = false;
			if (mySire==0 || (*this)[i]->pLevel==pseudoLevel) founder = true;
			if (mySire){
				founder = ((*this)[mySire-1]->pLevel==pseudoLevel) ? true : founder;
			}
			if (founder){ // founder, pseudo or "founder"
				pos[0] = 0;
				pos[1] = 0;
				pos[2] = (*this)[i]->pLevel;
				q[0]   = 0;
				q[1]   = 0;
				q[2]   = 1.0;
				v = (*this)[i]->varvi;
				addtoGinv(pos,q,1.0/v,lhs,startRow,startCol);
			}
			else { // not founder
				unsigned dad = myPedPtr->pedVector[i]->sire;
				pos[0] = (*this)[dad-1]->pLevel;
				pos[1] = (*this)[dad-1]->mLevel;
				pos[2] = (*this)[i]->pLevel;
				q[0]   = -(1-(*this)[i]->pPDQ);  
				q[1]   = -(*this)[i]->pPDQ;
				q[2]   = 1.0; 
  				double covvivjSire = (*this)[dad-1]->covvivj;
				double varviSire   = (*this)[dad-1]->varvi;
				double varvjSire   = (*this)[dad-1]->varvj;
				v      = (*this)[i]->varvi - (q[0]*varviSire*q[0]) 
                                       - (q[1]*varvjSire*q[1]) - 2*q[0]*q[1]*covvivjSire; 
				addtoGinv(pos,q,1.0/v,lhs,startRow,startCol);
			}
			maxLevel = pos[2]; 
		}
		if ((*this)[i]->mLevel > maxLevel){
			unsigned myDam = myPedPtr->pedVector[i]->dam;
			bool founder = false;
			if (myPedPtr->pedVector[i]->dam==0 || (*this)[i]->mLevel==pseudoLevel) founder = true;
			if (myDam) {
			 founder = ((*this)[myDam]->pLevel==pseudoLevel) ? true : founder;
			}
			if (founder) { // founder, pseudo or "founder"
				pos[0] = 0;
				pos[1] = 0;
				pos[2] = (*this)[i]->mLevel;
				q[0]   = 0;
				q[1]   = 0;
				q[2]   = 1.0;
				v = (*this)[i]->varvj;
				addtoGinv(pos,q,1.0/v,lhs,startRow,startCol);
			}
			else { // not founder
				unsigned mom = myPedPtr->pedVector[i]->dam;
				pos[0] = (*this)[mom-1]->pLevel;
				pos[1] = (*this)[mom-1]->mLevel;
				pos[2] = (*this)[i]->mLevel;
				q[0]   = -(1-(*this)[i]->mPDQ);  
				q[1]   = -   (*this)[i]->mPDQ;
				q[2]   = 1.0; 
   				double covvivjDam  = (*this)[mom-1]->covvivj;
				double varviDam   = (*this)[mom-1]->varvi;
				double varvjDam   = (*this)[mom-1]->varvj;
				v      = (*this)[i]->varvj - (q[0]*varviDam*q[0]) 
                                       - (q[1]*varvjDam*q[1]) - 2*q[0]*q[1]*covvivjDam; 
				addtoGinv(pos,q,1.0/v,lhs,startRow,startCol);  
			}
			maxLevel = pos[2];
		}
	}
}

	
	void MQTL::addtoGinv(unsigned pos[], double q[], double v,
						 matvec::doubleMatrix& lhs, unsigned startRow, 
						 unsigned startCol, double ratio){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		double vi,vj;
		unsigned ii,jj;
		for (unsigned i=0;i<3;i++){
			if(pos[i]){
				ii = startRow + pos[i] - 1;
				vi = q[i]*v;
				for (unsigned j=0;j<3;j++){
					if(pos[j]) { 
						vj = q[j];
						jj = startCol + pos[j] - 1;
						if(MIM::solMethod=="direct"){
							lhs[ii][jj] += ratio*vi*vj;
						}
						else{
							MIM::res[ii] += ratio*vi*vj * (*MIM::vec)[jj];
							if(ii==jj) MIM::diag[ii] += ratio*vi*vj;
						}
					}
				}
			}
		}
	}  
	
	void MQTL::addtoGinv(unsigned pos[], double q[], double v,
					 matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol){
	double vi,vj;
	unsigned ii,jj;
	for (unsigned i=0;i<3;i++){
		if(pos[i]){
			ii = startRow + pos[i] - 1;
			vi = q[i]*v;
			for (unsigned j=0;j<3;j++){
				if(pos[j]) { 
					vj = q[j];
					jj = startCol + pos[j] - 1;
					if(MIM::solMethod=="direct"){
						lhs[ii][jj] += vi*vj;
					}
					else{
						MIM::res[ii] += vi*vj * (*MIM::vec)[jj];
						if(ii==jj) MIM::diag[ii] += vi*vj;
					}
				}
			}
		}
	}
}  
	
} ////// end of namespace matvec




