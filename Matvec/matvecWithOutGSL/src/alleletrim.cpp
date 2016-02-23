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

#include <algorithm>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include "safe_vectors.h"
#include "stat.h"
#include "gnodederived.h"
#include "gnodesetderived.h"
#include "gnodestuff.h"
#include "population.h"
#include "model.h"
using namespace std; 

namespace matvec {
	
	void GNodeList::AfterSamplingAlleleOriginNode(GNode* sampleGNodePtr){
		unsigned locus;
		Individual* ind;
		set <GNodeSet*>::iterator GNstsiter;
		set<GNode*>::iterator setit;
		SafeSTLVector<unsigned> intersectionVector,dadvec,myvec,momvec;
		unsigned intersize =0;
		unsigned element;
		SafeSTLVector<MaternalPaternalAlleles> tempvec,indgenovec;
		SafeSTLVector<MaternalPaternalAlleles>::iterator stlit;
		
		for(GNstsiter = completeSetofGNsts.begin();GNstsiter != completeSetofGNsts.end();GNstsiter++){
			if (intersize > 0)  {intersectionVector.clear();  intersize = 0; break;} 
			if (dynamic_cast<RTransmissionSet*>(*GNstsiter)){ // check Transmission Sets 
				for(setit= (*GNstsiter)->begin();setit != (*GNstsiter)->end();setit++){ // check each element of each transmission set 
					if (intersize > 0) break; 
					if((*setit)->id == sampleGNodePtr->id){ // Sampled GNode found in some set
						ind = (reinterpret_cast<RTransmissionSet*>(*GNstsiter))->offspring;
						locus =  (reinterpret_cast<RTransmissionSet*>(*GNstsiter))->forLocus;
						intersectionVector.clear();
						tempvec.clear();
						indgenovec = ind->genotNodeVector[locus].genotypeVector; // individual genotype vector 
						if (reinterpret_cast<RTransmissionSet*>(*GNstsiter)->paternal && 
							popPtr->prior->chrom()[0].locus[locus].qtl_ml =='m'){ // we have sampled a paternal alleleorigin Node and marker locus
							element = ind->palleleOriginNodeVector[locus].alleleOriginVector[0];
							myvec = ind->palleleStateNodeVector[locus].alleleStateVector;
							if (element == 0) { //paternal allele from father's mother
								dadvec = ind->myfather->malleleStateNodeVector[locus].alleleStateVector;
							}
							else {// paternal allele from father's father
								dadvec = ind->myfather->palleleStateNodeVector[locus].alleleStateVector;
							}
							set_intersection(dadvec.begin(),dadvec.end(),myvec.begin(),myvec.end(),
											 inserter(intersectionVector,intersectionVector.begin())  );
							intersize = intersectionVector.size();
							if (intersize != myvec.size()){// Retain consistent elements of the individual genotype vector. 
								ind->palleleStateNodeVector[locus].alleleStateVector = intersectionVector;
								for(stlit = indgenovec.begin();stlit != indgenovec.end();stlit++){
									for(unsigned int1 = 0;int1 != intersize; int1++){
										if(stlit->paternal != intersectionVector[int1] -1) continue;
										if(stlit->paternal == intersectionVector[int1] -1) tempvec.push_back(*stlit);
									}
								}
							}
						}
						if ( !(reinterpret_cast<RTransmissionSet*>(*GNstsiter)->paternal) && 
							 popPtr->prior->chrom()[0].locus[locus].qtl_ml =='m'){ // we have sampled a maternal alleleorigin Node and marker locus   
							element = ind->malleleOriginNodeVector[locus].alleleOriginVector[0];
							myvec = ind->malleleStateNodeVector[locus].alleleStateVector;
							if (element == 0) { //maternal allele from mother's mother
								momvec = ind->mymother->malleleStateNodeVector[locus].alleleStateVector;
							}
							else {// maternal allele from mother's father
								momvec = ind->mymother->palleleStateNodeVector[locus].alleleStateVector;
							}
							set_intersection(momvec.begin(),momvec.end(),myvec.begin(),myvec.end(),
											 inserter(intersectionVector,intersectionVector.begin())  );
							intersize = intersectionVector.size();
							if (intersize != myvec.size()){
								ind->malleleStateNodeVector[locus].alleleStateVector = intersectionVector;
								for(stlit = indgenovec.begin();stlit != indgenovec.end();stlit++){
									for(unsigned int1 = 0;int1 != intersize; int1++){
										if(stlit->maternal != intersectionVector[int1] -1) continue;
										if(stlit->maternal == intersectionVector[int1] -1) tempvec.push_back(*stlit);
									}
								}
							}
						}
						if (tempvec.size() > 0){
							indgenovec.clear();
							ind->genotNodeVector[locus].genotypeVector.clear();
							ind->genotNodeVector[locus].genotypeVector = tempvec;
							tempvec.clear();
						}
					}					
				}
			}
		}
	}
	
	void GNodeList::AfterSamplingAlleleStateNode(GNode* sampleGNodePtr){
		unsigned locus;
		Individual* owner;
		set <GNodeSet*>::iterator GNstsiter;
		set<GNode*>::iterator setit;
		owner = NULL;
		bool paternal = false;
		bool maternal = false;
		
		for(GNstsiter = completeSetofGNsts.begin();GNstsiter != completeSetofGNsts.end();GNstsiter++){
			if (dynamic_cast<RAllelePenetranceSet*>(*GNstsiter)){ // check Penetrance Sets 
				for(setit= (*GNstsiter)->begin();setit != (*GNstsiter)->end();setit++){
					if((*setit)->id == sampleGNodePtr->id){ // Sampled GNode found in some set
						owner = (reinterpret_cast<RAllelePenetranceSet*>(*GNstsiter))->owner;
						locus = (reinterpret_cast<RAllelePenetranceSet*>(*GNstsiter))->forLocus;
						break;
					}
				}
			}
			if (dynamic_cast<RDisAllelePenetranceSet*>(*GNstsiter)){ // check Disease Penetrance Sets 
				for(setit= (*GNstsiter)->begin();setit != (*GNstsiter)->end();setit++){
					if((*setit)->id == sampleGNodePtr->id){ // Sampled GNode found in some set
						owner = (reinterpret_cast<RDisAllelePenetranceSet*>(*GNstsiter))->owner;
						locus = (reinterpret_cast<RDisAllelePenetranceSet*>(*GNstsiter))->forLocus;
						break;
					}
				}
			}
			if(owner != NULL) break;
		}
		if (owner == NULL) cout << "No owner of GNode found ! " << endl;
		if (owner == NULL) throw exception("Error in trimming alleles !");
		if (popPtr->prior->chrom()[0].locus[locus].qtl_ml =='m') {
			SafeSTLVector<unsigned> momvec = owner->malleleStateNodeVector[locus].alleleStateVector;
			SafeSTLVector<unsigned> dadvec = owner->palleleStateNodeVector[locus].alleleStateVector;
			SafeSTLVector<MaternalPaternalAlleles> genovec = owner->genotNodeVector[locus].genotypeVector;
			SafeSTLVector<MaternalPaternalAlleles> tempvec;
			SafeSTLVector<MaternalPaternalAlleles>::iterator stlit;
			tempvec.clear();
			
			if(momvec.size() == 1 && dadvec.size() != 1) { // maternal allele state vector sampled 
				maternal = true;
				for(stlit = genovec.begin(); stlit != genovec.end(); stlit++){
					if(stlit->maternal != momvec[0]-1) continue;
					if(stlit->maternal == momvec[0]-1) tempvec.push_back(*stlit);
				}
			}
			
			if(momvec.size() != 1 && dadvec.size() == 1) {// paternal allele state vector sampled 
				paternal = true;
				for(stlit = genovec.begin(); stlit != genovec.end(); stlit++){
					if(stlit->paternal != dadvec[0]-1) continue;
					if(stlit->paternal == dadvec[0]-1) tempvec.push_back(*stlit);
				}
			}
			
			if(momvec.size() == 1 && dadvec.size() == 1) {
				for(stlit = genovec.begin(); stlit != genovec.end(); stlit++){
					if(stlit->paternal != (dadvec[0]-1) || stlit->maternal !=  (momvec[0]-1) ) continue;
					if(stlit->paternal == (dadvec[0]-1) && stlit->maternal ==  (momvec[0]-1) ) tempvec.push_back(*stlit);
				}
			}
			genovec.clear();
			if (tempvec.size() == 0) cout << "Problem trimming allele state vector ! " << endl;
			owner->genotNodeVector[locus].genotypeVector.clear();
			owner->genotNodeVector[locus].genotypeVector = tempvec;
			if ( (paternal == true) || (maternal == true)  ){
				if (paternal) owner->malleleStateNodeVector[locus].alleleStateVector.clear();
				if (maternal) owner->palleleStateNodeVector[locus].alleleStateVector.clear();
				set<unsigned>alleleset;
				set<unsigned>::iterator it;
				alleleset.clear();
				for(stlit = tempvec.begin();stlit != tempvec.end(); stlit++){
					if (paternal) alleleset.insert(stlit->maternal +1);
					if (maternal) alleleset.insert(stlit->paternal +1);
				}
				for(it = alleleset.begin(); it != alleleset.end(); it++){
					if (paternal) owner->malleleStateNodeVector[locus].alleleStateVector.push_back(*it);
					if (maternal) owner->palleleStateNodeVector[locus].alleleStateVector.push_back(*it);
				}
			}
			tempvec.clear();
		}
	}
}