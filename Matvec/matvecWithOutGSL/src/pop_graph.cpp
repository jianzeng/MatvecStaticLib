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

#include <iomanip>
#include <cmath>
#include <map>
#include <cstdlib>
#include "session.h"
#include "population.h"
#include "individual.h"
#include "stat.h"
#include "model.h"

namespace matvec {
	
	
	void Population::descent_graph_init_parm(void) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: L. Radu Totir, Chris Stricker

		unsigned nLoci = Individual::numLoci;
		locusFreq.resize(nLoci);
		// Reading in frequencies for locus i into vector locusFreq
		for (unsigned i=0;i<nLoci; i++){
			locusFreq[i] = prior->chrom()[0].locus[i].allele_freq;
		}
		for (int i=0; i<popsize; i++) {
			popmember[i]->orderHetGenotype1.resize(nLoci);
			popmember[i]->orderHetGenotype2.resize(nLoci);
			popmember[i]->hetCounter.resize(nLoci);		
		}
	}
	
	void Population::descent_graph_setup(Data *data) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita, L. Radu Totir
		
		// set up all the necessary tables for descent graph calculations
		input_markerData(data);
		descent_graph_init_parm();
		int nLociInt = Individual::numLoci-1;
		RecoVector.resize(nLociInt);
		for (int i=1; i<=nLociInt; i++) {
			RecoVector(i) = model->MapF(prior->get_distance(1,i+1) - prior->get_distance(1,i));
		}
		std::cout << "RecoVector:\n" << RecoVector;
	}
	
	
	void Population::build_connected_groups(unsigned lcs) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: L. Radu Totir, Chris Stricker
		
		if (!founder_allele_counter) {
			throw exception("Individual::set_founder_alleles() must be called first" ); 
		}
		int zero = 0;
		connected_groups(lcs).resize(0);							//making sure that connected groups are destroyed
		connected_groups(lcs).resize(founder_allele_counter, zero); //making space to store connected_groups for all markerLoci.
		connect_counter.resize(0);									//making sure connect_counter is destroyed.
		connect_counter.resize(Individual::numLoci, zero);			//making space to store connect_counters for all markerLoci.
		
		connect_counter(lcs) = 0;
		for (int i=0; i<popsize; i++) {
			if(popmember[i]->genome0.chromosome[0].locus[lcs-1].allele && popmember[i]->genome1.chromosome[0].locus[lcs-1].allele) {           //genotype known
				popmember[i]->assigned_founder_allele=0;
				int mfounder = popmember[i]->m_founder(lcs);
				int pfounder = popmember[i]->p_founder(lcs);
				if (!connected_groups(lcs)(mfounder) && !connected_groups(lcs)(pfounder)){    //both zero
					connect_counter(lcs)++;
					connected_groups(lcs)(mfounder) = connect_counter(lcs);
					connected_groups(lcs)(pfounder) = connect_counter(lcs);
					//					std::cout <<"both zero i= "<<i<<" "<< connected_groups(lcs)(mfounder) << " " << connected_groups(lcs)(pfounder) << "mfounder = "<<mfounder<<" pfounder = "<<pfounder<<std::endl;
				}
				else if (connected_groups(lcs)(mfounder) && connected_groups(lcs)(pfounder)) {//both non-zero
																							  //set them equal
					int ccgpf = connected_groups(lcs)(pfounder);
					int ccgmf = connected_groups(lcs)(mfounder);
					for (int j=1; j<=founder_allele_counter; j++) {
						if (connected_groups(lcs)(j) == ccgpf ) {
							connected_groups(lcs)(j) = ccgmf;
						}
					}
					//					std::cout <<"both non-zero i= "<<i<<" "<< ccgmf <<std::endl;
				}
				else { //only one non-zero -> let the other inherit it's connected group
					if (connected_groups(lcs)(mfounder)) {
						connected_groups(lcs)(pfounder) = connected_groups(lcs)(mfounder);
					}
					else {
						connected_groups(lcs)(mfounder) = connected_groups(lcs)(pfounder);
					}
					//					std::cout <<"only one non-zero i= "<<i<<" "<< connected_groups(lcs)(mfounder) << " " << connected_groups(lcs)(pfounder) << std::endl;
				}
			}
		}
//		for(int i=0;i<popsize;i++){
//			std::cout<<"i="<<i<<" conn.groups(i, mfounder)="<<connected_groups(lcs)(popmember[i]->m_founder(lcs))<<" conn.groups(i, pfounder)="<<connected_groups(lcs)(popmember[i]->p_founder(lcs))<<std::endl;
//		}
		
	}
	
	
	
	
	
	int Population::update_allele_vectors(unsigned lcs, Individual *ind) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: L. Radu Totir, Chris Stricker
		
		int mfounder = ind->m_founder(lcs);
		int pfounder = ind->p_founder(lcs);
		int mallele =  ind->genome0.chromosome[0].locus[lcs-1].allele;
		int pallele =  ind->genome1.chromosome[0].locus[lcs-1].allele;
//		std::cout<<"in updateallele vec: mallele = "<<mallele<<" pallele= "<<pallele<<std::endl;		
		unsigned id = ind -> id();
//		std::cout << "update allele_vec: ind "<<ind->myid<< " mfounder = "<<mfounder<<" pfounder = "<<pfounder<<std::endl;
//		for (int j =1; j<=founder_allele_counter; j++){
//			if(allele_vector1(lcs)(j)) {
//				std::cout<<"j = "<<j<<" lcs = "<<lcs<<" allelevector1(lcs)(j) = "<<allele_vector1(lcs)(j)<<" allele_vector2(lcs)(j) = "<<allele_vector2(lcs)(j)<<std::endl;
//			}
//		}
		
		if (!allele_vector1(lcs)(mfounder) && !allele_vector1(lcs)(pfounder)) {  //both zero
			allele_vector1(lcs)(mfounder) = mallele; 
			allele_vector1(lcs)(pfounder) = pallele;
			if (mallele == pallele) { // homozygous
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector2(lcs)(j) = -1;
					}
				}
			}  
			else {   // heterozygous   
				allele_vector2(lcs)(mfounder) = pallele;
				allele_vector2(lcs)(pfounder) = mallele;
			}	    
		}
		else if (allele_vector1(lcs)(mfounder) && allele_vector1(lcs)(pfounder)) {// both non-zero
			bool M1 = (allele_vector1(lcs)(mfounder) == mallele && allele_vector1(lcs)(pfounder) == pallele ) ||
			(allele_vector1(lcs)(mfounder) == pallele && allele_vector1(lcs)(pfounder) == mallele );
			bool M2 = (allele_vector2(lcs)(mfounder) == mallele && allele_vector2(lcs)(pfounder) == pallele ) ||
				(allele_vector2(lcs)(mfounder) == pallele && allele_vector2(lcs)(pfounder) == mallele );
			if (M1 && M2) {
				;
			}
			else if (M1) {
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if (M2) {
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector1(lcs)(j) = allele_vector2(lcs)(j);
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else {
				return 1; // descent graph not consistent
			}
		}
		else if (allele_vector1(lcs)(mfounder)) { // only one non-zero
			bool M1 = (allele_vector1(lcs)(mfounder) == mallele);
			bool M2 = (allele_vector1(lcs)(mfounder) == pallele);
			bool M3 = (allele_vector2(lcs)(mfounder) == mallele);
			bool M4 = (allele_vector2(lcs)(mfounder) == pallele);
			if (M1 && M4) {
				allele_vector1(lcs)(pfounder) = pallele;
				allele_vector2(lcs)(pfounder) = mallele;
			}
			else if (M2 && M3) {
				allele_vector1(lcs)(pfounder) = mallele;
				allele_vector2(lcs)(pfounder) = pallele;
			}
			else if(M1) {
				allele_vector1(lcs)(pfounder) = pallele;
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if(M2) {
				allele_vector1(lcs)(pfounder) = mallele;
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if (M3) {
				allele_vector2(lcs)(pfounder) = pallele;
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++) {
					if (connected_groups(lcs)(j) == group) {
						allele_vector1(lcs)(j) = allele_vector2(lcs)(j);
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if (M4) {
				allele_vector2(lcs)(pfounder) = mallele;
				int group = connected_groups(lcs)(mfounder);
				for (int j =1; j<=founder_allele_counter; j++) {
					if (connected_groups(lcs)(j) == group) {
						allele_vector1(lcs)(j) = allele_vector2(lcs)(j);
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else {
				return 1;
			}
		}
		else if (allele_vector1(lcs)(pfounder)) {
			bool M1 = (allele_vector1(lcs)(pfounder) == mallele);
			bool M2 = (allele_vector1(lcs)(pfounder) == pallele);
			bool M3 = (allele_vector2(lcs)(pfounder) == mallele);
			bool M4 = (allele_vector2(lcs)(pfounder) == pallele);
			if (M1 && M4) {
				allele_vector1(lcs)(mfounder) = pallele;
				allele_vector2(lcs)(mfounder) = mallele;
			}
			else if (M2 && M3) {
				allele_vector1(lcs)(mfounder) = mallele;
				allele_vector2(lcs)(mfounder) = pallele;
			}
			else if(M1) {
				allele_vector1(lcs)(mfounder) = pallele;
				int group = connected_groups(lcs)(pfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if(M2) {
				allele_vector1(lcs)(mfounder) = mallele;
				int group = connected_groups(lcs)(pfounder);
				for (int j =1; j<=founder_allele_counter; j++){
					if (connected_groups(lcs)(j) == group) {
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if (M3) {
				allele_vector2(lcs)(mfounder) = pallele;
				int group = connected_groups(lcs)(pfounder);
				for (int j =1; j<=founder_allele_counter; j++) {
					if (connected_groups(lcs)(j) == group) {
						allele_vector1(lcs)(j) = allele_vector2(lcs)(j);
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else if (M4) {
				allele_vector2(lcs)(mfounder) = mallele;
				int group = connected_groups(lcs)(pfounder);
				for (int j =1; j<=founder_allele_counter; j++) {
					if (connected_groups(lcs)(j) == group) {
						allele_vector1(lcs)(j) = allele_vector2(lcs)(j);
						allele_vector2(lcs)(j) = -1;
					}
				}
			}
			else {
				return 1;
			}
		}
		else {
			throw exception("Oops bug in Population::build_allele_vectors");
		}
//		std::cout << "update allele_vec: ind "<<ind->myid<< " mfounder = "<<mfounder<<" pfounder = "<<pfounder<<std::endl;
//		for (int j =1; j<=founder_allele_counter; j++){
//			std::cout<<"j = "<<j<<" lcs = "<<lcs<<" allelevector1(lcs)(j) = "<<allele_vector1(lcs)(j)<<" allele_vector2(lcs)(j) = "<<allele_vector2(lcs)(j)<<std::endl;
//		}
		return 0;
	}
	
	
	
	void Population::show_descent_graph_stuff(unsigned lcs) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: L. Radu Totir, Chris Stricker
		
		std::cout << "\ndescent graph parameters:" << std::endl;
		for (int i=1; i<=founder_allele_counter; i++) {
			std::cout << "founder allele  " << i << "\tis in connected group "
			<< connected_groups(lcs)(i) << "\twith allele(s)\t"
			<< allele_vector1(lcs)(i)<< '\t' << allele_vector2(lcs)(i)<< std::endl;
		}
		std::cout << std::endl;
	}
	
	void Population::build_founder_allele_neighbors(unsigned lcs){
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: L. Radu Totir, Chris Stricker
		
		allele_vector1(lcs).resize(0); allele_vector2(lcs).resize(0); 
		allele_vector1(lcs).resize(founder_allele_counter,0); //need to put this somewhere before building
		allele_vector2(lcs).resize(founder_allele_counter,0); //allele vectors. So just stuck it here
		founder_allele_neighbors[lcs-1].resize(0);            //this sets the vector of vectors founder_allele_neighbors to a null vector for the 'second' dimension
		founder_allele_neighbors[lcs-1].resize(founder_allele_counter);
		// founder_allele_neighbors is, for each lcs, a vector with as many elements as founder_allele_counter. 
		//	Each position in the vector refers to a founder allele. We store at each position  
		//	a vector of genotyped individuals which carry that respective founder allele
		int mfounder, pfounder;
		for (int i=0; i<popsize; i++){
			if(popmember[i]->genome0.chromosome[0].locus[lcs-1].allele && popmember[i]->genome1.chromosome[0].locus[lcs-1].allele){ // genotype known
				mfounder = popmember[i]->m_founder(lcs);
				pfounder = popmember[i]->p_founder(lcs);  
//				std::cout<<"mfounder="<<mfounder<<" pfounder="<<pfounder<<std::endl;

				founder_allele_neighbors[lcs-1][mfounder-1].push_back(i);
				if(mfounder != pfounder) founder_allele_neighbors[lcs-1][pfounder-1].push_back(i);
				//			std::cout<< "i= " << i << " founder_allele_neighbours =  " << i << std::endl;
			}
		}
	}
	
	
	void Population::build_founder_allele_neighbors_all(unsigned lcs){
		// Authors: Chris Stricker 
		// (2004) 
		// Contributors: 
		
		founder_allele_neighbors_all[lcs-1].resize(0);
		founder_allele_neighbors_all[lcs-1].resize(founder_allele_counter);
		// founder_allele_neighbors is a vector with as many elements as founder_allele_counter. 
		//	Each position in the vector refers to a founder allele. We store at each position  
		//	a vector of individuals which carry that respective founder allele
		// irrespective whether that individual is genotyped or not
		// we need also to make sure, that only non-founder individuals are stored, since founders
		// may have singletons (indicated by "0" in the allele vector). The allele vector entries 
		// will be subsequently used in accumulateFounderHaplotypeOriginProbs(lcs) to access rows and
		// colums of a the matrix where these probs are stored. This matrix has the dimension of the 
		// number of alleles at flanking markers and start  with row and column 1.
		int mfounder, pfounder;
		for (int i=0; i<popsize; i++){
				mfounder = popmember[i]->m_founder(lcs);
				pfounder = popmember[i]->p_founder(lcs);  
				founder_allele_neighbors_all[lcs-1][mfounder-1].push_back(i);
				if(mfounder != pfounder) founder_allele_neighbors_all[lcs-1][pfounder-1].push_back(i);
				//			std::cout<< "i= " << i << " founder_allele_neighbours =  " << i << std::endl;
			}
	}
	
	
	int Population::build_allele_vector(unsigned lcs, unsigned connected_group){
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Chris Stricker

		valid_graph = 1; // graph is valid
						 //		std::cout<<"connect_counter(lcs) = "<<connect_counter(lcs)<<std::endl;
		for (int i=1; i<=founder_allele_counter; i++) {
			if(connected_groups(lcs)(i) == connected_group) {
				//		        std::cout<<"connected_groups(lcs)(i) = "<< connected_groups(lcs)(i) <<" i = "<< i<<std::endl;
				process_alleles_neighbors(lcs,i);
				return valid_graph; 
			}
		}
		return valid_graph;
	}
	
	
	void Population::process_alleles_neighbors(unsigned lcs, unsigned founder_allele){
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Chris Stricker
				
		std::vector<int>::iterator vec_iter;
		std::vector<int> temp;
		int neighbor_allele; 
		for (vec_iter =founder_allele_neighbors[lcs-1][founder_allele-1].begin(); 
			 vec_iter!=founder_allele_neighbors[lcs-1][founder_allele-1].end();
			 vec_iter++){ // loop through allele recepients which are stored at the vector position founder_allele-1
//			std::cout<<"vec_iter in process_allele_vector = "<<*vec_iter<<std::endl;
			Individual *ind = popmember[*vec_iter];
//			std::cout<<"Ind in process_allele_vector = "<<ind->myid<<std::endl;
			if (ind->assigned_founder_allele == 0) {
				ind->assigned_founder_allele = 1;
				int mfounder = ind->m_founder(lcs);
				int pfounder = ind->p_founder(lcs);
//				std::cout<<"in process_alleles_neighbors: mfounder="<<mfounder<<" pfounder="<<pfounder<<std::endl;

				if (founder_allele == mfounder) {
					neighbor_allele = pfounder;
				}
				else {
					neighbor_allele = mfounder;
				}
				temp.push_back(neighbor_allele);
//				std::cout<<"ind = "<<ind->myid<<" ind->assigned_founder_allele = "<<ind->assigned_founder_allele<<std::endl;
				int update_allele_vector_done=update_allele_vectors(lcs,ind);
				if (update_allele_vector_done== 1){
					valid_graph = 0; // graph is not valid
					//cout<<"Descent Graph is invalid at Individual "<<ind->myid<<" with its genotype at locus "<<lcs<<" that is incompatible with the allele vector at the relevant position. The graph will be rejected."<<endl;
				}
			}
		}
		for (vec_iter=temp.begin();vec_iter!=temp.end(); vec_iter++){
			process_alleles_neighbors(lcs, *vec_iter);
		}
//		for(unsigned i=1; i<=founder_allele_counter;i++) std::cout<<"in process allele neighbors: allele_vector1("<<lcs<<")("<<i<<") = "<<allele_vector1(lcs)(i)<<std::endl;
	}
	
	   
	void Population::sampleAlleleVectors(unsigned lcs){
		// Author: Chris Stricker
		// (2004) 
		// Contributors: 
		
		std::vector<int>::iterator vec_iter;
		double prior_prob1, prior_prob2;
		
		for (unsigned i=1; i<=connect_counter(lcs); i++) {
			prior_prob1 = 0.0;
			prior_prob2 = 0.0;
			for (unsigned j=1; j<=founder_allele_counter; j++) {
				if(i == connected_groups(lcs)(j)) {
					prior_prob1 += std::log(locusFreq(lcs)(allele_vector1(lcs)(j)));
					if (allele_vector2(lcs)(j) == -1) {
						prior_prob2 = 0.0;
					}
					else {
						prior_prob2 += std::log(locusFreq(lcs)(allele_vector2(lcs)(j)));
					}
				}
			}
			if(prior_prob2){ //i.e. if there is a second allele vector
				if(std::log(ranf()) > prior_prob1/(prior_prob1+prior_prob2)){ // i.e. the second allele vector gets sampled
					for (unsigned j=1; j<=founder_allele_counter; j++) {
						if(i == connected_groups(lcs)(j)) {
							allele_vector1(lcs)(j) = allele_vector2(lcs)(j); // fill the entries from the second allele vector into the first and erase the second allele vector.
							allele_vector2(lcs)(j) = -1;
						}
					}
				}
				else { //i.e. the first allele vector gets sampled
					for (unsigned j=1; j<=founder_allele_counter; j++) {
						if(i == connected_groups(lcs)(j)) {
							allele_vector2(lcs)(j) = -1; // erase the second allele vector.
						}
					}
				}
			}
		} 
		
		unsigned mfounderId, pfounderId, nalleles = prior->chrom()[0].locus[lcs-1].nallele();
		double rannumb, cutoff;
		for (unsigned i=0;i<popsize;i++){
			if(popmember[i]->mymother==0) {
				mfounderId = popmember[i]->m_founder(lcs);
				pfounderId = popmember[i]->p_founder(lcs);
				popmember[i]->previousSampledMaternalGenome.chromosome[0].locus[lcs-1].allele = popmember[i]->sampledMaternalGenome.chromosome[0].locus[lcs-1].allele;
				popmember[i]->previousSampledPaternalGenome.chromosome[0].locus[lcs-1].allele = popmember[i]->sampledPaternalGenome.chromosome[0].locus[lcs-1].allele;
				if(allele_vector1(lcs)(mfounderId)){
					popmember[i]->sampledMaternalGenome.chromosome[0].locus[lcs-1].allele =  allele_vector1(lcs)(mfounderId);
				}
				else{ // we sample the singletons according to marker freq.
					rannumb = ranf();
					cutoff=0;
					for(unsigned ii=1; ii<=nalleles; ii++){
						cutoff+= locusFreq(lcs)(ii);
						if(rannumb <= cutoff){
							popmember[i]->sampledMaternalGenome.chromosome[0].locus[lcs-1].allele =  ii;
							break;
						}
					}
				}
				if(allele_vector1(lcs)(pfounderId)){
					popmember[i]->sampledPaternalGenome.chromosome[0].locus[lcs-1].allele =  allele_vector1(lcs)(pfounderId);
				}
				else{ // sampling the singletons according to marker freq.
					rannumb = ranf();
					cutoff=0;
					for(unsigned ii=1; ii<=nalleles; ii++){
						cutoff+= locusFreq(lcs)(ii);
						if(rannumb <= cutoff){
							popmember[i]->sampledPaternalGenome.chromosome[0].locus[lcs-1].allele =  ii;
							break;
						}
					}
				}
			}
		}
	}
	
	
	
	bool Population::orderHeterozygotes(unsigned lcs){
		// Author: Chris Stricker
		// (2004) 
		// Contributors: 
		
		std::vector<int>::iterator ordallele1_iter, ordallele2_iter;
		bool found;
		unsigned pos =0;
		
		for (int i=0; i<size(); i++) {
			found = false;
			if (member(i)->mymother == 0 && member(i)->ord_heter == lcs) { // if founder and has a locus to be ordered
				pos=0;
				Individual *ind = member(i);
				if(allele_vector1(lcs)(ind->m_founder(lcs)) != allele_vector1(lcs)(ind->p_founder(lcs))){
					pos = 0;
					for(ordallele1_iter = ind->orderHetGenotype1[lcs-1].begin(); ordallele1_iter != ind->orderHetGenotype1[lcs-1].end(); ordallele1_iter++){
						pos++;
						if(*ordallele1_iter == allele_vector1(lcs)(ind->m_founder(lcs))){// ind allele at some position in orderHetGenotype1 vector corresponds to ind maternal allele in allele vector
							if(ind->orderHetGenotype2[lcs-1][pos-1] == allele_vector1(lcs)(ind->p_founder(lcs))){ // ind allele at the same position in orderHetGenotype2 vector corresponds to paternal allele
								found = true;
								//								std::cout<<"ACCEPTING ORDER for ind "<<ind->myid<<"\n";
								//								std::cout<<"ind->orderHetGenotype1["<<lcs-1<<"]["<<pos-1<<"] = "<< ind->orderHetGenotype1[lcs-1][pos-1]<<" ind->orderHetGenotype2["<<lcs-1<<"].["<<pos-1<<"] ="<<ind->orderHetGenotype2[lcs-1][pos-1]<<"\n";
								break; //ordered genotype found, order of het genotype is not violated, we can thus quit the loop
							}
						}
					}
					if(!found){
						pos = 0;
						for(ordallele2_iter = ind->orderHetGenotype2[lcs-1].begin(); ordallele2_iter != ind->orderHetGenotype2[lcs-1].end(); ordallele2_iter++){
							pos++;
							if(*ordallele2_iter == allele_vector1(lcs)(ind->m_founder(lcs))){// ind allele at some position in orderHetGenotype2 vector corresponds to ind maternal allele in allele vector
								if(ind->orderHetGenotype1[lcs-1][pos-1] == allele_vector1(lcs)(ind->p_founder(lcs))){ // ind allele at the same position in orderHetGenotype1 vector corresponds to paternal allele
																													  //									std::cout<<"REJECTING ORDER for ind "<<ind->myid<<"\n";
																													  //									std::cout<<"ind->orderHetGenotype1["<<lcs-1<<"]["<<pos-1<<"] = "<< ind->orderHetGenotype1[lcs-1][pos-1]<<" ind->orderHetGenotype2["<<lcs-1<<"]["<<pos-1<<"] ="<<ind->orderHetGenotype2[lcs-1][pos-1]<<"\n";
									return false; // genotype found, but wrong order quit the loop and reject the DG
								}
							}
						}
						ind->orderHetGenotype1[lcs-1].push_back(allele_vector1(lcs)(ind->m_founder(lcs)));
						ind->orderHetGenotype2[lcs-1].push_back(allele_vector1(lcs)(ind->p_founder(lcs)));
						std::cout<<"ASSIGNING ORDER FOR for ind "<<ind->myid<<"\n";
						std::cout<<"ind->orderHetGenotype1["<<lcs-1<<"]["<<pos<<"] = "<< ind->orderHetGenotype1[lcs-1][pos]<<" ind->orderHetGenotype2["<<lcs-1<<"]["<<pos<<"] ="<<ind->orderHetGenotype2[lcs-1][pos]<<"\n";
					}
				}
			}
		}
		return true;
	}
	
	
	
	void Population::countHeterozygotes(unsigned lcs){
		// Author: Chris Stricker
		// (2004) 
		// Contributors: 
		
		Individual *ind;
		for (int i=0; i<size(); i++) {
			ind = member(i);
			if (!ind->mymother) { // if founder
								  //				std::cout<<allele_vector1(lcs)(ind->m_founder(lcs))<<" "<<allele_vector1(lcs)(ind->p_founder(lcs))<<std::endl;
				if(allele_vector1(lcs)(ind->m_founder(lcs)) != allele_vector1(lcs)(ind->p_founder(lcs))){
					ind->hetCounter(lcs)++;
					//					std::cout<<"ind = "<<ind->myid<<" lcs = "<<lcs<<" hetCounter = "<<ind->hetCounter(lcs)<<std::endl;
				}
			}
		}
	}
	
	double Population::findOptimumFounderLocusToOrder(int orderSamples, int PQTL, bool &initSamplerToOrderFounders){
		// Author: Chris Stricker
		// (2004) 
		// Contributors: 
		
		double l_hood = 0.0; 
		
		unsigned s1 = 10;               //# of rounds of SL sampler
		unsigned s2 = 10;               //# of rounds of haplotype sampler
		unsigned s3 = 10;               //# of rounds of SL cascading sampler
		
		// s1+s2+s3 are not supposed to exceed orderSamples
		
		unsigned s = 0, samplerOption = 0;	
		
		for(int sample=1; sample <= orderSamples; sample++) {
			if (s >= s1 + s2 + s3)                {s = 0;}
			if (s < s1)                           {samplerOption = 1;}
			else if ((s >= s1) && (s < s1 + s2))  {samplerOption = 2;}
			else if (s <= s1 + s2 + s3)           {samplerOption = 3;}
			s++;
			
			l_hood = M_H_sample_ext(PQTL, l_hood, samplerOption, initSamplerToOrderFounders); 
			std::cout << "orderSample = " << sample << "  log lh = "
				<< l_hood << "  samplerOption: " << samplerOption << std::endl;
			
		}
		pickFounderLocusToOrder(orderSamples, PQTL);
		return(l_hood);
	}
	
	
	
	void Population::pickFounderLocusToOrder(unsigned orderSamples, unsigned PQTL){
		// Author: Chris Stricker
		// (2004) 
		// Contributors: 
		
		double recRate, informativityCurrentLocus=0, informativityPreviousLocus=0, dist1;
		unsigned nLoci = Individual::numLoci;
		std::cout<<"These are the founders with their respective loci that are used to order their haplotypes: "<<std::endl;
		double dist2 = prior->get_distance(1,PQTL);
		for (int i=0; i<size(); i++) {
			if(member(i)->mymother == 0){
				informativityPreviousLocus = 0;
				for( unsigned j =1; j<=nLoci; j++){
					if(j!=PQTL){
						dist1 = prior->get_distance(1,j);
						if(dist1>dist2) recRate = model->MapF(dist1-dist2);
						else recRate = model->MapF(dist2-dist1);
						informativityCurrentLocus = member(i)->hetCounter(j)/orderSamples * (0.5-recRate)/0.5; // criterion to pick the best locus to order w.r.t distance and heterozygosity
//						cout<<"ind="<<member(i)->myid<<" recrate="<<recRate<<" hetcounter ="<<member(i)->hetCounter(j)<<" infcurrent="<<informativityCurrentLocus<<" infprev="<<informativityPreviousLocus<<"Locus="<<j<<endl;
						if(informativityCurrentLocus > informativityPreviousLocus){
							informativityPreviousLocus = informativityCurrentLocus;
							member(i)->ord_heter = j;
						}
					}
				}
				std::cout<<"Ind = "<<member(i)->myid<<" locus = "<<member(i)->ord_heter<<"\n";
			}
			else member(i)->ord_heter = 0;
		}
		std::cout<<std::endl;
	}
	
		
	double Population::calc_prior_descent_graph(unsigned lcs, bool &initSamplerToOrderFounders) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita, Chris Stricker
		
		build_connected_groups(lcs);
		build_founder_allele_neighbors(lcs);
		
		//	std::cout<<"connect_counter(lcs) = "<<connect_counter(lcs)<<std::endl;
		for (int i=1; i<=connect_counter(lcs);i++){
			int valid = build_allele_vector(lcs,i);
			if (valid == 0){
				return 0.0; // we should never get here
			}
		}
		//for(unsigned i=1; i<=founder_allele_counter;i++) std::cout<<"in calc_prior_descent_graph: allele_vector1("<<lcs<<")("<<i<<") = "<<allele_vector1(lcs)(i)<<std::endl;

		if(!orderFoundersGenotypes){
			if(initSamplerToOrderFounders){
				throw exception("Something wrong with option 'initSamplerToOrderFounders' and 'orderFoundersGenotyps': Both cannot be set to true. Aborting");
			}
		}
		else if(!initSamplerToOrderFounders){
			if(!orderHeterozygotes(lcs)) return 0.0; // this rejects the DG if order is not respected for het founders
		}
		
		// CS 04/17/2003 added scaling, likelihood of DG will be scaled!
		// Scaling here is tricky, since we cannot simply take log's as prior_prob1 
		// and prior_prob2 are added! So we will scale prior_prob1 and prioir_prob2 
		// with the same scaling factor. The log of this scaling factor is accumulated 
		// into pop_graph_scaling, which in turn is returned for de-scaling in 
		// descent_graph_log_lhood().
		
		double prior_prob = 1.0;
		double prior_prob1;
		double prior_prob2;
		double scale;
		pop_graph_scaling = 0;
		
		for (int i=1; i<=connect_counter(lcs); i++) {
			prior_prob1 = 1.0;
			prior_prob2 = 1.0;
			scale = 0;
			bool found = false;
			for (int j=1; j<=founder_allele_counter; j++){
				if (connected_groups(lcs)(j) == i) {
					found = true;
//					std::cout<<" connected_groups("<<lcs<<")("<<j<<") = "<<connected_groups(lcs)(j)<<std::endl;
//					std::cout<<" allele_vector1("<<lcs<<")("<<j<<") = "<<allele_vector1(lcs)(j)<<std::endl;
//					std::cout<<" locusFreq("<<lcs<<")("<<j<<") = "<<locusFreq(lcs)(allele_vector1(lcs)(j))<<std::endl;
//					std::cout<<"in calc_prior_descent_graph: lcs = "<<lcs<<" j = "<<j<<" allelevector1(lcs)(j) = "<<allele_vector1(lcs)(j)<<" markerfreq = "<<locusFreq(lcs)(allele_vector1(lcs)(j))<<std::endl;
					prior_prob1 *= locusFreq(lcs)(allele_vector1(lcs)(j));
//					std::cout<<"prior_prob1 = "<<prior_prob1<<"lcs = "<<lcs<<"allele = "<<allele_vector1(lcs)(j)<<"locusFreq = "<<locusFreq(lcs)(allele_vector1(lcs)(j))<<"\n";
					if(prior_prob1< 1/std::pow(10.0,20.0)){
						scale += std::log(prior_prob1);
						//std::cout<<"Scaling prior_prob1 (the likelihood) ...scale = "<<scale<< "\n";
						prior_prob1 = 1.0;
					}
					//					std::cout<<"prior_prob1 "<<prior_prob1<<"\n";
				}
			}
			for (int j=1; j<=founder_allele_counter; j++){
				if (connected_groups(lcs)(j) == i) {
					if (allele_vector2(lcs)(j) == -1) {
						prior_prob2 = 0.0;
						break;
					}
					// 	  prior_prob2 *= locusFreq(lcs)(allele_vector2(lcs)(j));
					prior_prob2 += std::log(locusFreq(lcs)(allele_vector2(lcs)(j)));
//					std::cout<<"prior_prob2 "<<prior_prob2<<"\n";
				}
			}
			// the following statement scales prior_prob2 by the same scaling factor
			// as prior_prob1 by subtracting the log of the scaling factor from the 
			//log(prior_prob2). But if there is no second allele vector, i.e.
			// prior_prob2=0 or no connected group can be found ('singleton', i.e.
			// prior_prob2=1), then these values are to be used as they were never scaled.
			
			if(found) {
			if(prior_prob2 < 1 && prior_prob2 != 0) prior_prob2 = std::exp(prior_prob2-scale);
				if(prior_prob < 1/std::pow(10.0,20.0)){
//					std::cout<<"Scaling prior_prob (the likelihood) ...\n";
				scale += std::log(prior_prob);
				prior_prob = 1.0;
			}
			pop_graph_scaling += scale;
				prior_prob *= prior_prob1 + prior_prob2;
			}
			
//			std::cout<<"delog prior_prob = "<<prior_prob<<" connect_counter(lcs) = "<<i<<" scaling = "<<pop_graph_scaling<<" prob1 = "<<prior_prob1<<" prob2 = "<< prior_prob2<<" scale = "<<scale<<"\n\n";
		}
		
		//		std::cout<<"PriorProb= "<<prior_prob<<"\n";
		return prior_prob;
	}
	
	
	double Population::descent_graph_log_lhood(bool &initSamplerToOrderFounders) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita, L. Radu Totir, Chris Stricker
		
		// calc prior probability:
		double prior_prob = 0.0, temp=0.0;  
		unsigned nLoci = Individual::numLoci;
		for (int i=0; i<nLoci; i++) {
			if(prior->chrom()[0].locus[i].qtl_ml=='q') continue;
			double temp = 1.0;
			temp = calc_prior_descent_graph(i+1, initSamplerToOrderFounders);
			// std::cout << "prior marker locus " << i << ": " << temp << std::endl;
			if (temp == 0.0) {  //i.e. graph is not valid, return inpossible log lh and return  
			    //cout << "graph is not valid at locus " << i+1 << endl; 
				return 9999.9;
			}
			prior_prob += std::log(temp) + pop_graph_scaling;
			//pop_graph_scaling is the log of the scaling factor from calc_prior_descent_graph()
			
		}
		//     std::cout << "log prior prob: " << prior_prob<< std::endl;
		//calc transmission probability
		double trans_prob = 0.0;
		for (int i=1; i<size(); i++) {
			if (member(i)->mymother) { // if not founder
				double temp = 0.25;
				for (int j=2; j<=nLoci; j++) {
					if (member(i)->m_gamete(j-1) == member(i)->m_gamete(j)) {
						temp *= 1 - RecoVector(j-1);
					}
					else {
						temp *= RecoVector(j-1);
					}
					if (member(i)->p_gamete(j-1) == member(i)->p_gamete(j)) {
						temp *= 1 - RecoVector(j-1);
					}
					else {
						temp *= RecoVector(j-1);
					}
				}
				trans_prob += std::log(temp);
			}
		}
		//     std::cout << "log transmission prob: " << trans_prob<< std::endl;
		//return log lh
		return prior_prob + trans_prob;
	}

	
	double Population::MH_ibd_sample(unsigned l, double log_lhood, unsigned &option, bool &initSamplerToOrderFounders){
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Chris Stricker
		
		// updates ibd states at locus  l
		// then gets new sample
		update_ibdMatrix(l);
		
		double llh;
		llh = M_H_sample_ext(l,log_lhood, option, initSamplerToOrderFounders);
		//	std::cout<<"llh = "<<llh<<std::endl;
		return llh;
	}
	
	double Population::M_H_sample_ext(unsigned PQTL, double log_lhood, unsigned &option, bool &initSamplerToOrderFounders){
	  	// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita, Chris Stricker
		
		//options are:
		//1 = pure SL sampler
		//2 = SL with haplotype sampler
		//3 = SL sampler with cascading
		unsigned rejected=0, accepted=0;
		
		//		std::cout<<" option = "<<option<<"\n";
		
		if ((option < 1)  || (option > 3)) {
			throw exception("nonexistent option, abort.");
		}
		
		//first, calc init likelihood if not already done
		if(log_lhood >= 0.0) {
			initial = 0;
			log_lhood = descent_graph_log_lhood(initSamplerToOrderFounders);
			std::cout << "initial log likelihood: " << log_lhood << std::endl;
			if(log_lhood > 0.0) {
				throw exception("could not calculate initial log lh for MH-Sampler, aborting.");
			}
		}
		initial = 1;
		
		//start
		//store actual gamete sources
		for (int i=0; i<size(); i++) {
			member(i)->m_gamete1 = member(i)->m_gamete;
			member(i)->p_gamete1 = member(i)->p_gamete;
		}
		
		
		
		//choose popmember
		unsigned node = unsigned(ranf()*size());  // allowing entire pedigree
												  //  unsigned node = unsigned(ranf()*2);   // using individual 1 and 2 only
		
		//choose marker
		unsigned nLoci = Individual::numLoci;
		SLlocus = unsigned(ranf()*nLoci); // samples locus at random
		while(prior->chrom()[0].locus[SLlocus].qtl_ml=='q') SLlocus = unsigned(ranf()*nLoci);
		//continues to pick a locus until its a marker, not a QTL
		SLlocus+=1;
		unsigned more=1, rule=0;
		while (more) { //do transitions
			
			//choose transition rule, avoid nodes without offsprings for T1/T2 rule
			//avoid nodes without parents for T0 rule
			if (member(node)->numoffs && member(node)->mymother) {
				//parents and offsprings -> all rules
				rule=(unsigned(ranf()*4));
			}
			else if (member(node)->numoffs) {
				//no parents, but offsprings -> T1/T2
				rule=(unsigned(ranf()*3))+1;
			}
			else if (member(node)->mymother) {
				//parents, but no offsprings -> T0
				rule=0;
			}
			else {
				throw exception("single unrelated individual encountered while choosing transition rules");
			}
			
			//apply SL transition to SL locus    
			member(node)->apply_SL_transition(rule, SLlocus);
			
			if (option == 3) {
				member(node)->apply_SL_cascade(rule, SLlocus);
			}
			
			//do another transition?
			if (unsigned(ranf()*2)) {
				more = 0;
			}
			else {
				more=1;
				//choose a new node? Factor changes proportion!
				if (!(unsigned(ranf()*2))) {
					node = unsigned(ranf()*size());
				}
			}
		}
		
		double q_ratio = 0;
		
		if (option == 2) {
			//generate haplotypes and calc q
			q_ratio = calc_log_q_ratio();
		}
		
		//calculate new likelihood and calc Metropolis Hastings (if graph invalid, log_lhood > 1)
		double old_log_lhood = log_lhood;
		for (int lcs=1; lcs<=nLoci; lcs++) {
			previousConnectCounter(lcs) = connect_counter(lcs);
			for (unsigned j=1; j<=founder_allele_counter; j++) {
				previousAlleleVector1(lcs)(j) = allele_vector1(lcs)(j);
				previousAlleleVector2(lcs)(j) = allele_vector2(lcs)(j);
				previousConnectedGroups(lcs)(j) = connected_groups(lcs)(j);
				unsigned numinds = founder_allele_neighbors[lcs-1][j-1].size();
				previousFounderAlleleNeighbors[lcs-1][j-1].resize(0);
				previousFounderAlleleNeighbors[lcs-1][j-1].resize(numinds);
				for(unsigned k=1; k<=numinds; k++) {
					previousFounderAlleleNeighbors[lcs-1][j-1][k-1] = founder_allele_neighbors[lcs-1][j-1][k-1];
				}
			}
		}
		log_lhood =  descent_graph_log_lhood(initSamplerToOrderFounders);
		//		std::cout << "new log likelihood: " << log_lhood << std::endl;
		
		
		if ((log_lhood > 1) || !(ranf() <= std::exp(log_lhood - old_log_lhood + q_ratio))) {
			//changes not accepted -> set back changes
			//			std::cout << "rejected\n" << std::endl;
			log_lhood = old_log_lhood;
			//set sampled alleles for founders back
			for (int i=0; i<size(); i++) {
				if(member(i)->mymother == 0){
					for (int lcs=1; lcs<=nLoci; lcs++) {
						member(i)->sampledMaternalGenome.chromosome[0].locus[lcs-1].allele = member(i)->previousSampledMaternalGenome.chromosome[0].locus[lcs-1].allele;
						member(i)->sampledPaternalGenome.chromosome[0].locus[lcs-1].allele = member(i)->previousSampledPaternalGenome.chromosome[0].locus[lcs-1].allele;
					}
				}			
				//set all founder indicators back
				member(i)->m_gamete = member(i)->m_gamete1;
				member(i)->p_gamete = member(i)->p_gamete1;
				
				//reset founder_alleles (updating)
				member(i)->set_founder_alleles(1);
			}
			for (int lcs=1; lcs<=nLoci; lcs++) {
				connect_counter(lcs) = previousConnectCounter(lcs);
				for (unsigned j=1; j<=founder_allele_counter; j++) {
					allele_vector1(lcs)(j) = previousAlleleVector1(lcs)(j);
					allele_vector2(lcs)(j) = previousAlleleVector2(lcs)(j);
					connected_groups(lcs)(j) = previousConnectedGroups(lcs)(j);
					unsigned numinds = previousFounderAlleleNeighbors[lcs-1][j-1].size();
					founder_allele_neighbors[lcs-1][j-1].resize(0);
					founder_allele_neighbors[lcs-1][j-1].resize(numinds);
					for(unsigned k=1; k<=numinds; k++) {
						founder_allele_neighbors[lcs-1][j-1][k-1] = previousFounderAlleleNeighbors[lcs-1][j-1][k-1];
					}
				}
			}
		}
		//    else std::cout << "accepted\n" << std::endl;
		
		if(initSamplerToOrderFounders){
			for(unsigned j=1;j<=nLoci;j++){
				if(j!=PQTL) countHeterozygotes(j);
			}
		}
		//		std::cout<<"log_lhood = "<<log_lhood<<std::endl;
		return log_lhood;
	}
	
	double Population::calc_log_q_ratio(void) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita, Chris Stricker
		
		double accumulator = 0.0;
		for (int i=0; i<size(); i++) {
			unsigned flag = 0;
			if (member(i)->m_gamete(SLlocus) != member(i)->m_gamete1(SLlocus)) {
				flag = 1;
				member(i)->sample_haplotypes(0);
			}
			if (member(i)->p_gamete(SLlocus) != member(i)->p_gamete1(SLlocus)) {
				flag = 1;
				member(i)->sample_haplotypes(1);
			}
			if (flag != 0) {
				accumulator += member(i)->calc_q();
				member(i)->update_offsprings_founder_alleles();
			}
		}
		return accumulator;
	}
	
	
	void Population::show_change(const int ind, const int counter) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita
		unsigned p_flag = 0, m_flag = 0;
		for (int i = 1; i <= member(ind)->m_gamete.size(); i++) {
			if (member(ind)->m_gamete(i) != member(ind)->m_gamete1(i)) {m_flag++;}
			if (member(ind)->p_gamete(i) != member(ind)->p_gamete1(i)) {p_flag++;}
		}
		if (m_flag > 0 || p_flag > 0) {
			std::cout << "Individual "<< member(ind)->myid << ", round: " << counter << std::endl;
			for (int i = 1; i <= member(ind)->m_gamete.size(); i++) {
				std::cout << member(ind)->m_gamete(i) << " || "  << member(ind)->p_gamete(i) << std::endl;
			}
		}
		return;
	}
	
	void Population::initPDQs(void){
		// Authors: Rohan L. Fernando 
		// (September, 2005) 
		// Contributors: 
		for (int i=0; i<size(); i++) {
			if (member(i)->mymother){ 
				member(i)->m_counter = 0;
				member(i)->p_counter = 0;
			}
		}		
	}

	void Population::sum_descentState(unsigned locus) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita
		
		for (int i=0; i<size(); i++) {
			if (member(i)->mymother){ 
				member(i)->m_counter += member(i)->m_gamete(locus);
				member(i)->p_counter += member(i)->p_gamete(locus);
			}
		}
	}
	
	
	void Population::sum_descentState_map(void) {  
		// Authors: Fabiano V. Pita and Rohan L. Fernando 
		// (2003) 
		// Contributors: Chris Stricker
		unsigned nLoci = Individual::numLoci;
		for (int i=0; i<size(); i++) {
			if (member(i)->mymother){ 
				for (int locus=1;locus<=nLoci;locus++){
					member(i)->m_counter_map(locus) += member(i)->m_gamete(locus);
					member(i)->p_counter_map(locus) += member(i)->p_gamete(locus);
//					cout<<"Ind="<<member(i)->myid<<" locus="<<locus<<" m_counter_map="<<member(i)->m_counter_map(locus)<<" p_counter_map="<<member(i)->p_counter_map(locus)<<endl;
				}
			}
		}
	}

	void Population::Gibbs_sample(){
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita
		
		for (int i=0; i<size(); i++){
			member(i)->sample_self(0);//maternal
			member(i)->sample_self(1);//paternal
		}
	}
	
	void Population::update_ibdMatrix(unsigned locus){ 
		// Authors: Rohan L. Fernando and Fabiano V. Pita
		// (2002) 
		// Contributors: 
		
		if (ibdMatrix.empty()){
			ibdMatrix.resize(popsize,popsize,0.0);
		}
		for (unsigned i=1;i<=popsize;i++){
			unsigned im = member(i-1)->m_founder(locus);
			unsigned ip = member(i-1)->p_founder(locus);
			for (unsigned j=i;j<=popsize;j++){
				unsigned jm = member(j-1)->m_founder(locus);
				unsigned jp = member(j-1)->p_founder(locus);
				ibdMatrix(i,j) += (im==jm) + (im==jp) + (ip==jm) + (ip==jp);
			}
		}
	}
	
	void Population::accumulateFounderHaplotypeOrigin(unsigned lcs){
		// Authors: Chris Stricker 
		// (2004) 
		// Contributors: 

		std::vector<int>::iterator vec_iter;
		unsigned nLoci = Individual::numLoci;
		for(unsigned markerlocus=1; markerlocus<=nLoci; markerlocus++){
			sampleAlleleVectors(markerlocus);
			build_founder_allele_neighbors_all(markerlocus);
		}
		for (int i=0; i<popsize; i++){			
			if(popmember[i]->mymother==0){ // i is a founder
										   // this assumes that the loci are entred in the GeneticDist object according to map position
				int mfounder = popmember[i]->m_founder(lcs);
				int pfounder = popmember[i]->p_founder(lcs);  
				int maternalAi = popmember[i]->sampledMaternalGenome.chromosome[0].locus[lcs-2].allele; // these are the sampled allele states of the flanking markers
				int maternalBj = popmember[i]->sampledMaternalGenome.chromosome[0].locus[lcs].allele;
				int paternalAi = popmember[i]->sampledPaternalGenome.chromosome[0].locus[lcs-2].allele; 
				int paternalBj = popmember[i]->sampledPaternalGenome.chromosome[0].locus[lcs].allele;
				//				std::cout<<"ind="<<popmember[i]->myid<<" maternalAi="<<maternalAi<<" maternalBj="<<maternalBj<<" paternalAi="<<paternalAi<<" paternalBj="<<paternalBj<<std::endl;
				for (vec_iter =founder_allele_neighbors_all[lcs-1][mfounder-1].begin(); 
					 vec_iter!=founder_allele_neighbors_all[lcs-1][mfounder-1].end();
					 vec_iter++){ // loop through allele recepients which are stored at the vector position mfounder-1
					Individual *ind = popmember[*vec_iter];
					int indmfounder = ind->m_founder(lcs);
					int indpfounder = ind->p_founder(lcs);
					if (indmfounder == mfounder) {
						ind->maternalFounderHaplotypeProbs(maternalAi,maternalBj)++;
					}
					if(indpfounder == mfounder) {
						ind->paternalFounderHaplotypeProbs(maternalAi,maternalBj)++;
					}
				}
				for (vec_iter =founder_allele_neighbors_all[lcs-1][pfounder-1].begin(); 
					 vec_iter!=founder_allele_neighbors_all[lcs-1][pfounder-1].end();
					 vec_iter++){ // loop through allele recepients which are stored at the vector position pfounder-1
					Individual *ind = popmember[*vec_iter];
					int indmfounder = ind->m_founder(lcs);
					int indpfounder = ind->p_founder(lcs);
					if (indmfounder == pfounder) {
						ind->maternalFounderHaplotypeProbs(paternalAi,paternalBj)++;
					}
					if(indpfounder == pfounder) {
						ind->paternalFounderHaplotypeProbs(paternalAi,paternalBj)++;
					}
				}
			}
		}
	}
	
	
	void Population::initFounderHaplotypeOriginProbs(unsigned QTL){
		// Authors: Rohan L. Fernando and Chris Stricker
		// (2004) 
		// Contributors: 
		
		// number of alleles in locus to the left of the QTL
		unsigned nallelesAlocus = prior->chrom()[0].locus[QTL-2].nallele();
		// number of alleles in locus to the right of the QTL
		unsigned nallelesBlocus = prior->chrom()[0].locus[QTL].nallele();
		for (unsigned i=0;i<popsize;i++){ 
			member(i)->maternalFounderHaplotypeProbs.resize(nallelesAlocus,nallelesBlocus,0.0);
			member(i)->paternalFounderHaplotypeProbs.resize(nallelesAlocus,nallelesBlocus,0.0);
		}
	}
	
	void Population::displayFounderHaplotypeOriginProbs(unsigned nsamples){
		// Authors: Rohan L. Fernando and Chris Stricker
		// (2004) 
		// Contributors: 
		std::cout <<"Matrix FounderHaplotypeProbs:" << std::endl;
		const char *strid;
		for (unsigned i=0;i<popsize;i++){ 
			int id = member(i)->id();
			std::cout << ind_name(id) << std::endl;
			std::cout << member(i)->maternalFounderHaplotypeProbs/nsamples << std::endl;
			std::cout << member(i)->paternalFounderHaplotypeProbs/nsamples << std::endl;
		}
	}

	void Population::accumulateHaplotypes(){
		// Authors: Chris Stricker
		// (2004) 
		// Contributors: 
		
		unsigned nLoci = Individual::numLoci;
		for(unsigned markerlocus=1; markerlocus<=nLoci; markerlocus++){
			sampleAlleleVectors(markerlocus);
			build_founder_allele_neighbors_all(markerlocus);
		}
		unsigned mfounder, pfounder;
		string matHapIndex, patHapIndex;
		int l;
		double ten = 10.0;
		for (int i=0; i<popsize; i++){	
			l = 1;
			matHapIndex = "";
			patHapIndex = "";
			while(l<=nLoci){
				mfounder = popmember[i]->m_founder(l);
				pfounder = popmember[i]->p_founder(l);
				toString<int>(allele_vector1(l)(mfounder), std::dec);
				matHapIndex = matHapIndex + "-" + toString<int>(allele_vector1(l)(mfounder), std::dec) ;
				patHapIndex = patHapIndex + "-" + toString<int>(allele_vector1(l)(pfounder), std::dec) ;
				l++;
			}
			matHapIndex = matHapIndex + "-";
			patHapIndex = patHapIndex + "-";			
			if(popmember[i]->matHaplotypeMap.find(matHapIndex) != popmember[i]->matHaplotypeMap.end() ) {
				popmember[i]->matHaplotypeMap[matHapIndex] +=1;
			}
			else{
				popmember[i]->matHaplotypeMap[matHapIndex] = 1;
			}
			if(popmember[i]->patHaplotypeMap.find(patHapIndex) != popmember[i]->patHaplotypeMap.end() ) {
				popmember[i]->patHaplotypeMap[patHapIndex] +=1;
			}
			else{
				popmember[i]->patHaplotypeMap[patHapIndex] = 1;
			}
		}
	}
	
	
	void Population::displayHaplotypes(unsigned nsamples){
		// Authors: Chris Stricker
		// (2004) 
		// Contributors: 
		std::cout <<"'Reconstructed' Haplotypes with corresponding Probs:" << std::endl;
		const char *strid;
		map<string, int>::iterator m;
		for (unsigned i=0;i<popsize;i++){ 
			int id = member(i)->id();
			std::cout << "Ind. " << ind_name(id) <<std::endl;
			std::cout << "maternal haplotypes" << std::endl;
			for(m = member(i)->matHaplotypeMap.begin(); m != member(i)->matHaplotypeMap.end(); ++m) {
				std::cout << m->first <<"\t\t" << ((double) m->second)/nsamples<< std::endl;
			}
			std::cout << "paternal haplotypes" << std::endl;
			for(m = member(i)->patHaplotypeMap.begin(); m != member(i)->patHaplotypeMap.end(); ++m) {
				std::cout << m->first <<"\t\t" << ((double) m->second)/nsamples << std::endl;
			}
			std::cout << std::endl << std::endl;
		}
		std::cout << std::endl << std::endl;
	}
	
	
	void Population::outputFounderHaplotypeOriginProbs(unsigned nsamples, char* fname){
		// Authors: Rohan L. Fernando and Chris Stricker
		// (2004) 
		// Contributors: 

		const char *strid;
		std::ofstream outfile(fname);
		for (unsigned i=0;i<popsize;i++){ 
			int id = member(i)->id();
			outfile << ind_name(id)<< std::endl;
			outfile << member(i)->maternalFounderHaplotypeProbs/nsamples << std::endl;
			outfile << member(i)->paternalFounderHaplotypeProbs/nsamples << std::endl;
		}
	}
	
	double Population::MH_ibd_sample_map(unsigned l, double log_lhood, unsigned &option, bool &initSamplerToOrderFounders){
		// Authors: Fabiano V. Pita and Rohan L. Fernando 
		// (2003) 
		// Contributors: Chris Stricker
		// gets the current haplotype for each individual 
		// and all pairs of loci -- QTL Mapping stuff
		for (int i=0; i<size(); i++) {
			member(i)->Individual::count_haplotype();
		}
		double llh = M_H_sample_ext(l,log_lhood, option, initSamplerToOrderFounders);
		return llh;
	}
	
	void Population::descent_graph_setup(MIM& mim) {
		// Authors: Mathias Schelling and Rohan L. Fernando 
		// (1999) 
		// Contributors: Fabiano V. Pita, L. Radu Totir
		
		// set up all the necessary tables for descent graph calculations
		input_markerData(mim);
		descent_graph_init_parm();
		int nLociInt = Individual::numLoci-1;
		RecoVector.resize(nLociInt);
		for (int i=1; i<=nLociInt; i++) {
			double xab = prior->get_distance(1,i+1) -  prior->get_distance(1,i);
			RecoVector(i) = 0.5*(1-std::exp(-2.0*xab)); // Haldane
		}
		std::cout << "RecoVector:\n" << RecoVector;
	}
	

	
	
}//end of namespace matvec

