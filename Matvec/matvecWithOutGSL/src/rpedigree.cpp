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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cctype>
#include <sstream>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>
#include "rpedigree.h"
#include "mim.h"
using namespace std; 

namespace matvec{
	void RPedigree::putColNames(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" \t,");
		colName.getTokens(str,sep);
	}
	
	void RPedigree::inputPed(string fileString){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		const char* fname = fileString.c_str();
		cout << "reading pedigree file \n";
		if (colName.size() < 3) {
			cerr << "RPedigree::input(): colName.size() < 3 \n";
			cerr << "Did you forget the putColNames method? \n";
			throw exception("Error in RPedigree::input()");
		}
		int indIndex = colName.getIndex("individual");
		if (indIndex == -1){
			cerr << "RPedigree::input(): individual column is missing in colName \n";
			throw exception("Error in RPedigree::input()");
		}
		int sireIndex = colName.getIndex("sire");
		if (sireIndex == -1){
			cerr << "RPedigree::input(): sire column is missing in colName \n";
			throw exception("Error in RPedigree::input()");
		}
		int damIndex = colName.getIndex("dam");
		if (damIndex == -1){
			cerr << "RPedigree::input(): dam column is missing in colName \n";
			throw exception("Error in RPedigree::input()");			
		}
		unsigned numCol = colName.size();			
		double rec = 0, rec1 = 0;
		string indstr, sirestr, damstr;
		ifstream datafile(fname);
		if(!datafile){
			cout<< "Cannot open pedigree file: " << fname << endl;
			exit(1);
		}
		datafile.setf(ios::skipws);
		PNode *ptr;
		std::string sep(" \t,\n\r");
		std::string inputStr;
		Tokenizer colData;
		unsigned COUNT = 0;
		while ( getline(datafile,inputStr) ){
			colData.getTokens(inputStr,sep);
			indstr  = colData[indIndex];
			sirestr = colData[sireIndex];
			damstr  = colData[damIndex];
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}
			if (colData.size() != numCol){
				cerr << " Record " << rec1 + rec << " has " << colData.size() << " columns \n";
				cerr << " Expected " << numCol << endl;
				throw exception("Error in RPedigree::input()"); 
			}
			ptr = new PNode(indstr, sirestr, damstr);
			if (orderedPed) ptr->ind = ++COUNT;
			(*this)[indstr] = ptr;
		}
		datafile.close();
		if(orderedPed){
			seqnPed();
		}
		else {
			generateEntriesforParents();
			codePed();
		}	
		makePedVector();
		calc_inbreeding();
		fillCoder(); 
	}
	
	void RPedigree::displayPed(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		RPedigree::iterator it;
		SafeSTLVector<PNode*>::iterator vecit;
		for (vecit=pedVector.begin();vecit!=pedVector.end();vecit++){
			cout << setw(10) << (*vecit)->ind 
			<< setw(10) << (*vecit)->sire
			<< setw(10) << (*vecit)->dam
			<< setw(20) << ((*vecit)->ind_str).c_str()
			<< setw(20) << (*vecit)->f <<     endl;
		}
		cout << endl;
	}
	
	void RPedigree::generateEntriesforParents(void) {
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		// if a parent does not have an entry, we make it a founder
		cout << "\n generating missing entries for parents \n";
		unsigned sire_count = 0;
		unsigned dam_count = 0;
		double rec = 0, rec1 = 0;
		RPedigree::iterator it, parent_it;
		for(it=begin();it!=end();it++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}
			PNode *ptr = (*it).second;
			if(ptr->sire_str!="0"){
				parent_it = (*this).find(ptr->sire_str);
				if(parent_it == end()){ // sire has no entry
					PNode *ptrs = new PNode(ptr->sire_str, "0", "0");
					(*this)[ptr->sire_str] = ptrs;
				}
			}
			if(ptr->dam_str!="0"){
				parent_it = (*this).find(ptr->dam_str);
				if(parent_it == end()){ // dam has no entry
					PNode *ptrd =  new PNode(ptr->dam_str, "0", "0");
					(*this)[ptr->dam_str] = ptrd;
				}
			}
		}
	}
	
	void RPedigree::codePed(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		cout << "\n coding pedigree \n";
		RPedigree::iterator it;
		COUNT = 0;
		unsigned rec = 0, rec1 = 0;
		
		for(it=begin();it!=end();it++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			};
			cout.flush();
			PNode *ptr =(*it).second;
			code(ptr);                                     
		}
	}
	
	void RPedigree::code(PNode *ptr){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		if(ptr->ind != -1) { // already coded 
			return;
		}
		if(ptr->sire_str == "0" && ptr->dam_str == "0"){ // founder
			ptr->ind  = ++COUNT;
			ptr->sire = 0;
			ptr->dam  = 0;
		}
		else if(ptr->sire_str != "0" && ptr->dam_str == "0"){ // dam missing, sire is not missing
			
			PNode* sire_ptr = (*this)[ptr->sire_str];
			if (sire_ptr->ind == -1) { 
				code(sire_ptr);
			}
			ptr->ind  = ++COUNT;
			ptr->sire = sire_ptr->ind;
			ptr->dam  = 0;
		}
		else if(ptr->dam_str != "0" && ptr->sire_str == "0"){ // sire missing, dam is not missing
			PNode* dam_ptr = (*this)[ptr->dam_str];
			if (dam_ptr->ind == -1) {
				code(dam_ptr);
			}
			ptr->ind  = ++COUNT;
			ptr->sire = 0;
			ptr->dam  = dam_ptr->ind;
		}
		else{
			PNode* sire_ptr = (*this)[ptr->sire_str];
			if (sire_ptr->ind == -1) {
				code(sire_ptr);
			}
			PNode* dam_ptr = (*this)[ptr->dam_str];
			if (dam_ptr->ind == -1) { 
				code(dam_ptr); 
			}
			ptr->ind = ++COUNT;
			ptr->sire = sire_ptr->ind;
			ptr->dam  = dam_ptr->ind;
		}
	}
	
	void RPedigree::seqnPed(void){
		// Authors: Rohan L. Fernando
		// (November, 2005) 
		// Contributors:	
		cout << "\n assigning sequential ids \n";
		RPedigree::iterator it;
		unsigned COUNT = 0;
		unsigned rec = 0, rec1 = 0;
		
		for(it=begin();it!=end();it++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			};
			cout.flush();
			PNode *ptr =(*it).second;
			if(ptr->sire_str == "0" && ptr->dam_str == "0"){ // founder
				ptr->sire = 0;
				ptr->dam  = 0;
			}
			else if(ptr->sire_str != "0" && ptr->dam_str == "0"){ // dam missing, sire is not missing
				
				PNode* sire_ptr = (*this)[ptr->sire_str];
				ptr->sire = sire_ptr->ind;
				ptr->dam  = 0;
			}
			else if(ptr->dam_str != "0" && ptr->sire_str == "0"){ // sire missing, dam is not missing
				PNode* dam_ptr = (*this)[ptr->dam_str];
				ptr->sire = 0;
				ptr->dam  = dam_ptr->ind;
			}
			else{
				PNode* sire_ptr = (*this)[ptr->sire_str];
				PNode* dam_ptr = (*this)[ptr->dam_str];
				ptr->sire = sire_ptr->ind;
				ptr->dam  = dam_ptr->ind;
			}
		}
	}
	
	void RPedigree::calc_inbreeding(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		SafeSTLVector <PNode*>::iterator it;
		unsigned rec = 0, rec1 = 0, non_rec = 0;
		cout << "\n calculating inbreeding \n";
		for (it=pedVector.begin();it!=pedVector.end();it++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			};
			(*it)->f = get_rij((*it)->sire,(*it)->dam);
		}	
	}
	
	double RPedigree::get_rij(int i, int j){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		if (i==0||j==0){
			return 0.0;
		}
		double x = SpCij.retrieve_cij(i,j);
		if(x != -1.0) {
			return x;
		}
		int old, young;
		if(i < j){
			old = i;
			young = j;
		}
		else if(j < i){
			old = j;
			young = i;
		}
		else{
			double f = pedVector[i-1]->f;
			x = 0.5*(1 + f);
			SpCij.put_cij(i,j,x);
			return x;
		}
		int y_sire = pedVector[young-1]->sire;
		int y_dam  = pedVector[young-1]->dam;
		x = (get_rij(old,y_sire)+get_rij(old,y_dam))/2.0;
		SpCij.put_cij(i,j,x);
		return x;
	}
	
	void RPedigree:: output(char* ped){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		ofstream pedfile(ped);
		SafeSTLVector<PNode*>::iterator vecit;
		pedfile.setf(ios::fixed | ios::right);
		for (vecit=pedVector.begin();vecit!=pedVector.end();vecit++){
			pedfile << setw(10) << (*vecit)->ind 
			<< setw(10) << (*vecit)->sire
			<< setw(10) << (*vecit)->dam
			<< setw(20) << ((*vecit)->ind_str).c_str()
			<< setw(20) << (*vecit)->f <<     endl;
		}
	}
	
	void RPedigree::makePedVector(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		RPedigree::iterator it;
		pedVector.resize(size());
		for(it=begin();it!=end();it++){
			PNode *ptr = (*it).second;
			unsigned i = ptr->ind - 1;
			pedVector[i] = ptr;
		}
	}
	
	void RPedigree::fillCoder(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		SafeSTLVector<PNode *>::iterator it;
		for(it=pedVector.begin();it!=pedVector.end();it++){
			coder.code((*it)->ind_str);
		}
	}
	
	void RPedigree::addAinv(matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol, double ratio){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		bool direct = (MIM::solMethod=="direct") ? true : false;
		double q[3];
		double d,fs,fd;
		unsigned pos[3];
		vector<PNode*>::iterator it;
		if (coder.size()>pedVector.size()){
			cout << "RPedigree is not complete \n";
			exit(-1);
		}
		for (it=pedVector.begin();it!=pedVector.end();it++){
			pos[0] = (*it)->sire;
			pos[1] = (*it)->dam;
			pos[2] = (*it)->ind;
			if((*it)->sire && (*it)->dam){
				q[0] = -0.5;
				q[1] = -0.5;
				q[2] =  1.0;
				fs = pedVector[pos[0]-1]->f;
				fd = pedVector[pos[1]-1]->f;
				d = 4.0/(2 - fs - fd);
			}
			else if((*it)->sire){
				q[0] = -0.5;
				q[1] =  0.0;
				q[2] =  1.0;
				fs = pedVector[pos[0]-1]->f;
				d = 4.0/(3-fs);
			}
			else if((*it)->dam){
				q[0] =  0.0;
				q[1] = -0.5;
				q[2] =  1.0;
				fd = pedVector[pos[1]-1]->f;
				d = 4.0/(3-fd);
			}
			else{
				q[0] =  0.0;
				q[1] =  0.0;
				q[2] =  1.0;
				d = 1.0;
			}
			for (unsigned i=0;i<3;i++){
				if(pos[i]){
					unsigned ii = startRow + pos[i] - 1;
					for (unsigned j=0;j<3;j++){
						if(pos[j]) { 
							unsigned jj = startCol + pos[j] - 1;
							if (direct){
								lhs[ii][jj] += ratio*q[i]*d*q[j];
							}
							else {
								MIM::res[ii] += ratio*q[i]*d*q[j] * (*MIM::vec)[jj];
								if (ii==jj)  MIM::diag[ii] += ratio*q[i]*d*q[j];
							}
						}
					}
				}
			}
		}
	}
	
} ////// end of namespace matvec


