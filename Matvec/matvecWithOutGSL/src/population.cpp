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
#include "nufamily.h"
#include "stat.h"
#include "mim.h"
#include <iomanip>
#include <sstream>
using namespace std;

namespace matvec {

GeneticDist* Population::prior;
	
	
	
	Population::Population(void)              // one of constructors
{
		initialize();
}

Population::Population(GeneticDist *D)    // one of constructors
{
	initialize();
	prior = D;
}

Population::Population(const Population& A)
{
	initialize();
	copyfrom(A);
}

Population::Population(const int maxsize, GeneticDist* D)
{
	resize(maxsize,D); // hashtable needs to resize() afterwards if necessary
}

void Population::copyfrom(const Population& A)
{
	if (this == &A)
		return;
	unsigned i;
	this->resize(A.maxpopsize, A.prior);
	prior              = A.prior;
	fdone              = A.fdone;
	stdid              = A.stdid;
	popsize            = A.popsize;
	numchrom           = A.numchrom;
	numtrait           = A.numtrait;
	nmarker            = A.nmarker;
	maxpopsize         = A.maxpopsize;
	maxnamelen         = A.maxnamelen;
	//RLF
	n_markerLoci       = A.n_markerLoci;
	nuFamiliesDone      = A.nuFamiliesDone;
	//RLF
	trans_mat          = A.trans_mat;
	blupsol            = A.blupsol;
	hashtable          = A.hashtable;
	evaluate_method    = A.evaluate_method;
	mark_need_remove   = A.mark_need_remove;
	
	spouse_info_built  = 0;
	offs_info_built    = 0;
	
	popsize = A.popsize;
	
	for (i=0; i<popsize; i++)
		pm_storage[i] = A.pm_storage[i];
	for (i=0; i<numchrom; i++) {
		pop_gamete[i] = A.pop_gamete[i];
		gametebase[i] = A.gametebase[i];
	}
	
	if (nmarker != A.nmarker) {
		if(markersymbol){
			delete [] markersymbol;
			markersymbol=0;
		}
		if(A.nmarker>0){
			markersymbol = new std::string [A.nmarker];
		}
		else {
			markersymbol = 0;
		}
		nmarker = A.nmarker;
	}
	for (i=0; i<nmarker; i++) markersymbol[i] = A.markersymbol[i];
}

void Population::initialize(void)
{
    model        = 0;
	myPedPtr     = 0;
	myRPedPtr    = 0;
	popsize      = 0;
	maxpopsize   = 0;
	numchrom     = 0;
	numtrait     = 0;
	nmarker      = 0;
	fdone        = 0;
	stdid        = 0;
	num_nufamily = 0;
	maxnamelen   = 0;
	//RLF
	n_markerLoci = 0;
	nuFamiliesDone = false;
	//RLF
	sex_confirmed     = 0;
	mark_need_remove  = 0;
	spouse_info_built = 0;
	offs_info_built   = 0;
	
	prior        = 0;
	trans_mat    = 0;
	pop_gamete   = 0;
	gametebase   = 0;
	pm_storage   = 0;
	popmember    = 0;
	nufamily_vec = 0;
	markersymbol = 0;
	markerDataIn = false;
	// LRT
	lookToYourLeft = true;
	lookToYourRight= true;
}

void Population::resize(const unsigned maxn,GeneticDist *D)
{
	if (maxn == maxpopsize && prior == D) return;
	prior = D;
	if (maxn != maxpopsize) {
		if (pm_storage) {
			delete [] pm_storage;
			pm_storage=0;
		}
		if (popmember) {
			delete [] popmember;
			popmember=0;
		}
		if(maxn>0){
			pm_storage = new Individual [maxn];
		}
		else {
			pm_storage = 0;
		}
		for (int i=0; i<maxn; ++i) pm_storage[i].remodel(this);
		if(maxn>0){
			popmember  = new Individual* [maxn];
		}
		else {
			popmember = 0;
		}
		for (unsigned i=0; i<maxn; i++) popmember[i] = &(pm_storage[i]);
	}
	maxpopsize = maxn;
	if (pop_gamete) {delete [] pop_gamete; pop_gamete = 0;}
	if (gametebase) {delete [] gametebase; gametebase = 0;}
	
	numchrom          = prior->nchrom();
	numtrait          = prior->ntrait();
	popsize           = 0;
	fdone             = 0;
	stdid             = 0;
	
	mark_need_remove  = 0;
	sex_confirmed     = 0;
	offs_info_built   = 0;
	spouse_info_built = 0;
	nuFamiliesDone    = false;
}

unsigned Population::nbase(void) const
{
	unsigned n = 0;
	for (unsigned i=0; i<popsize; i++) if (popmember[i]->base()) n++;
	return n;
}

const Population& Population::operator=(const Population& A)
{
	copyfrom(A);
	return *this;
}

unsigned Population::input_ped(Pedigree& P,GeneticDist *D){
	myPedPtr = &P;
	unsigned ret = input_ped0(P,D);
	if (strcmp(D->name(),"UnknownDist") != 0) {
		if (P.sex_confirmed) {
			confirm_sex(1);
		}
		else {
			confirm_sex(-1);
		}
		countFounders();
		build_offs_info();
		build_spouse_info();
	}
    return ret;
}

void Population::input_ped(RPedigree& P,GeneticDist& D){
	myRPedPtr = &P;
	prior = &D;
	input_ped0(D);
	countFounders();
	build_offs_info();
	build_spouse_info();
}

unsigned Population::input_ped0(Pedigree& P,GeneticDist *D)
{
	//////////////////////////////////////////////////
	// return 0 if everything is successful
	//        1  something wrong
	/////////////////////////////////////////////////
	if (P.size() == 0) {
		warning("Population::input_ped(pedigree,data): pedigree is empty");
		return 1;  // This has been changed from mv12
	}
	resize(P.size()+P.ngroup(),D);
	hashtable.resize(maxpopsize);
	popsize = maxpopsize;
	maxnamelen = P.maxnamelen;
	if (P.in_memory()) P.release_pedsheet();
	unsigned i,s, idssex[4];
	//   unsigned i,s, isd[4];  old declaration
	
	size_t recsize = sizeof(unsigned)*4;
	unsigned nanim;
	nanim=P.size();
	Individual* I;
	std::ifstream pedfile(P.diskfname().c_str(),std::ios::in);
	//std::cout << P.diskfname().c_str() << " " << recsize << "\n";
	if (!pedfile) throw exception("Population::input_ped(): cannot open file");
	popmember--;
	for (i=1; i<=nanim; i++) {
		I = popmember[i];
		pedfile.read((char *)idssex,recsize);
		I->myid = idssex[0];
		if (idssex[1] > 0) {
			I->mymother = popmember[idssex[1]];
		}
		else {
			I->mymother = 0;
		}
		if (idssex[2] > 0) {
			I->myfather = popmember[idssex[2]];
		}
		else {
			I->myfather = 0;
		}
		I->mysex = (char)idssex[3];
		// std::cout << I->myid << " "<< I->mymother << " " << I->myfather << "\n";
	}
	for(i=(nanim+1);i<=popsize;i++) {
		I = popmember[i];
		I->myid=i;
	}
	popmember++;
	if (P.maxnamelen > 0) {
		char *strid = new char [P.maxnamelen+1];
		for (i=0; i<popsize; i++) {
			pedfile.read((char *)&s,sizeof(unsigned));
			pedfile.read(strid,s);
			hashtable.insert(strid);
		}
		if(strid){
			delete [] strid;
			strid=0;
		}
	}
	stdid = 1;
	if(!P.inbcoef_done()) P.inbcoef(); //SDK Guarantee that the
									   //inbreeding coefficients are
									   //available when needed.
	if (P.inbcoef_done()) {
		fdone = 1;
		double *fvec = P.inbcoef()->begin();
		for (i=0; i<popsize; i++) {
			I = popmember[i];
			I->inbc = *fvec++;
		}
	}
	if (strcmp(D->name(),"UnknownDist") != 0) {
		if (P.sex_confirmed) {
			confirm_sex(1);
		}
		else {
			confirm_sex(-1);
		}
	}
	pedfile.close();
	return 0;// This has been changed from mv12
}

void Population::input_ped0(GeneticDist& D)
{
	// Authors: Rohan L. Fernando 
	// (August, 2005) 
	// Contributors: 
	
	if (myRPedPtr->size() == 0) {
		throw exception("Population::input_ped0(Rpedigree,GeneticDist): pedigree is empty");
	}
	unsigned n = myRPedPtr->size();
	cout <<"\n Building Indiviual List   \n";
	resize(n,&D); // for some reason, this sets popsize to 0 
	popsize = n;
	Individual* I;
	PNode* pedPtr;
	unsigned rec = 0,rec1=0;
	for (unsigned i=0; i<popsize; i++) {
		rec++;
		if(rec==1000){
			cout<<rec1+rec <<"\r";			
			cout.flush();
			rec1+=rec;
			rec = 0;
		}
		I = popmember[i];
		pedPtr = myRPedPtr->pedVector[i];
		I->myid = pedPtr->ind; // actually this is equivalent to i+1
		if (pedPtr->dam > 0) {
			I->mymother = popmember[pedPtr->dam - 1];
		}
		else {
			I->mymother = 0;
		}
		if (pedPtr->sire > 0) {
			I->myfather = popmember[pedPtr->sire - 1];
		}
		else {
			I->myfather = 0;
		}
		I->mysex = 'U';
		I->inbc = pedPtr->f;
	}
	confirm_sex(-1);
	stdid = 1;
}


unsigned Population::input_ped(const char fname[],const char recfmt[],
							   GeneticDist *D)
{
	//////////////////////////////////////////////////
	// return 0 if everything is successful
	//        1  something wrong
	/////////////////////////////////////////////////
	Pedigree ped;
	ped.input(fname,recfmt);
	return input_ped(ped,D);
}

unsigned Population::input_data(Data* D)
{
	if (!D || D->num_rows() <= 0) throw exception("Population::input_data(): data is empty");
	unsigned numcol = D->num_cols();
	unsigned numrec = D->num_rows();
	unsigned i,j,t,id,dummy_ind;
	const char *strid;
	
	if (!D->in_memory()) D->input_datasheet();
	for (t=1; t<numcol; t++) {   // 1st column is reserved for intercept
								 //   for (t=0; t<numcol; t++) {  Old version
		if (D->datasheet[t].col_struct.classi()=='G') break;
								 }
	if (t==numcol) throw exception("Population::input_data(): no pedigree exits");
	DataNode  *tcol, *column = D->rawcol(t);
	HashTable* dhtable = D->hashtable[t];
	if (strcmp(prior->name(),"GeneticDist") == 0) {
		///////////////////////////////////////////////////////////
		// need more work for finding marker-locus in prior
		///////////////////////////////////////////////////////////
		Individual *I;
		int gtype[2];
		if (prior->chrom()[0].locus[1].qtl_ml == 'm') {
			for (t=0; t<numcol; t++) {
				if (D->datasheet[t].col_struct.name()== "marker1") break;
			}
			if (t==numcol) throw exception(" Population::input_data(): no marker1 exits");
			tcol =  D->rawcol(t);         // marker1 column
			for (i=0; i<numrec; i++) {
				if (tcol[i].double_val() == 0.0) {
					gtype[0]=0;
					gtype[1]=0;
				}
				else if (tcol[i].double_val() == 1.0) {
					gtype[0]=0;
					gtype[1]=1;
				}
				else if (tcol[i].double_val() == 2.0) {
					gtype[0]=1;
					gtype[1]=1;
				}
				else  {
					throw exception("Population::input_data(): marker1: invalid genotype");
				}
				
				id = column[i].unsigned_val();
				strid = (const char *)dhtable->find(id);
				id  = hashtable.get_id(strid);
				if (id == 0) continue;
				I = popmember[id-1];
				I->genome0.chromosome[0].locus[1].allele = gtype[0];
				I->genome1.chromosome[0].locus[1].allele = gtype[1];
			}
		}
	}
	
	for (j=0, t=0; t<numcol; t++) {
		if (D->datasheet[t].col_struct.classi()=='T') j++;
	}
	if (j != numtrait) throw exception("Population::input_data(): Pop's ntrait != Data's ntrait");
	// this is where the phenotype from the data file is
	// merged with the pedigree informatiojn
	for (j=0, t=1; t<numcol; t++) {   // 1st column is reserved for intercept
		if (D->datasheet[t].col_struct.classi()=='T') {
			tcol =  D->rawcol(t);
			dummy_ind = 0;
			for (i=0; i<numrec; i++) {
				id = column[i].unsigned_val();
				strid = (const char *)dhtable->find(id);
				id  = hashtable.get_id(strid);
				if (id > 0) {
					popmember[id-1]->myrecord[j] = tcol[i];
				}
				else {
					column[i].missing = 1;
					dummy_ind++;
				}
			}
			j++;
		}
	}
	D->release_datasheet();
	if (j==0) throw exception("Population::input_data(): trait(s) not available");
	return (numrec - dummy_ind);
	}

void Population::confirm_sex(const int col_sex)
{
	//////////////////////////////////////////////////////////////////////////
	// if Pedigree has an explicity column for,
	//       then it overwrite 'mother' and 'father' columns, i.e. mother
	//       and father just mean parents, not sexuality
	//  else 'mother' and 'father' columns make their sexuality
	//////////////////////////////////////////////////////////////////////////
	
	Individual *I,*temp,*mother,*father;
	unsigned i;
	if (col_sex >= 0) {             // sex has been already set, now check
		for (i=0; i<popsize; i++) {   //  to see it is confirmable
			I = popmember[i];
			mother = I->mother();
			father = I->father();
			if (mother && father) {
				if (mother->sex() != '.' && mother->sex() == father->sex()) {
					throw exception(std::string("Population::confirm_sex(): ") + ind_name(i+1));
				}
				else if (mother->sex() == 'M') {
					temp = mother;
					I->mymother = father;
					I->myfather = temp;
				}
				else if (father->sex() == 'F') {
					temp = father;
					I->myfather = mother;
					I->mymother = temp;
				}
			}
			else if (father && !mother) {
				if (father->sex() == 'F') {
					temp = father;
					I->myfather = mother;
					I->mymother = temp;
				}
			}
			else if (!father && mother) {
				if (mother->sex() == 'M') {
					temp = mother;
					I->mymother = father;
					I->myfather = temp;
				}
			}
		}
	}
	else {
		for (i=0; i<popsize; i++) {
			I = popmember[i];
			mother = I->mother();
			father = I->father();
			//         if (mother && mother->mysex == 'M') throw exception(std::string("its mother has been father earlier: ")+ind_name(i+1));
			//        if (father && father->mysex == 'F') throw exception(std::string("its father has been mother earlier") + ind_name(i+1));
			if (mother) {
				if (mother->mysex == 'M') {
					throw exception(std::string("its mother has been father earlier: ") + ind_name(i+1));
				}
				else {
					mother->mysex = 'F';
				}
			}
			if (father) {
				if  (father->mysex == 'F') {
					throw exception(std::string("its father has been mother earlier: ") + ind_name(i+1));
				}
				else {
					father->mysex = 'M';
				}
			}
		}
		///////////////////////////////////////////////////////////////////
		// if I has never been parent, then set it as a female.
		// actually whether it is set to be male or female does not matter
		///////////////////////////////////////////////////////////////////
		for (i=0; i<popsize; i++) {
			I = popmember[i];
			if (I->sex() == '.') I->mysex = 'F';
		}
	}
}

void Population::countFounders(void){
	// Authors: L. Radu Totir
	// (September, 2004) 
	// Contributors:
	Individual *ind;
	numFounders=0;
	for(int i=0;i<popsize;i++){
		ind=popmember[i];
		if(!(ind->mymother)){
			numFounders++;
		}
	}
}

void Population::build_offs_info(void)
{
	if (offs_info_built) return;
	if (!stdid) renum();
	Individual *I;
	unsigned i,j,noffs;
	for (i=0; i<popsize; i++) {
		I = popmember[i];
		if (I->mother()) I->mother()->offs_tree.insert(I);
		if (I->father()) I->father()->offs_tree.insert(I);
	}
	for (i=0; i<popsize; i++) {
		I = popmember[i];
		if (I->myoffspring) {
			delete [] I->myoffspring;I->myoffspring=0;
		}
		noffs = I->offs_tree.size();
		I->numoffs = noffs;
		if (noffs > 0) {
			I->myoffspring = new Individual* [noffs];
			std::set<Individual *>::iterator pos;
			for (j=0,pos=I->offs_tree.begin(); pos != I->offs_tree.end(); ++pos,++j) {
				I->myoffspring[j] = *pos;
			}
			I->offs_tree.clear();
		}
	}
	mark_need_remove = 1;
	offs_info_built = 1;
}

void Population::build_spouse_info(void)
{
	if (spouse_info_built) return;
	if (!stdid) renum();
	if (!offs_info_built) build_offs_info();
	
	Individual *I,*spouse,**offspring;
	unsigned i,j,k,s,noffs,nsp;
	for (i=0; i<popsize; i++) {
		I = popmember[i];
		if (I->numoffs_spouse) {
			delete [] I->numoffs_spouse; I->numoffs_spouse = 0;
		}
		if (I->spouselist) {
			delete [] I->spouselist; I->spouselist = 0;
		}
		noffs = I->numoffs;
		offspring = I->offspring();
		if (noffs == 1) {
			I->numspouse = 1;
			I->spouselist = new Individual*[1];
			if (I->mysex == 'F') {
				I->spouselist[0] = offspring[0]->father();
			}
			else {
				I->spouselist[0] = offspring[0]->mother();
			}
			I->numoffs_spouse = new unsigned [1];
			I->numoffs_spouse[0] = 1;
		}
		else if (noffs > 1) {
			j = 1;  nsp = 1;
			if (I->mysex == 'F') {
				qsort(offspring,noffs,sizeof(Individual*),compare_father_id);
				spouse = offspring[0]->father();
				while(j<noffs) {
					if (spouse != offspring[j]->father()) {
						spouse = offspring[j]->father();
						nsp++;
					}
					j++;
				}
				I->numspouse = nsp;
				if(nsp>0){
					I->numoffs_spouse = new unsigned [nsp];
				}
				else {
					I->numoffs_spouse = 0;
				}
				if(nsp>0){
					I->spouselist = new Individual* [nsp];
				}
				else {
					I->spouselist = 0;
				}
				spouse = offspring[0]->father();
				I->spouselist[0] = spouse;
				j = 1; s = 0; k = 1;
				while(j<noffs) {
					if (spouse != offspring[j]->father()) {
						spouse = offspring[j]->father();
						I->numoffs_spouse[s++] = k;
						I->spouselist[s] = spouse;
						k = 1;
					}
					else {
						k++;            // count # of myoffspring of each spouses
					}
					j++;
				}
				I->numoffs_spouse[s] = k;
			}
			else {
				qsort(offspring,noffs,sizeof(Individual*),compare_mother_id);
				spouse = offspring[0]->mother();
				while(j<noffs) {
					if (spouse != offspring[j]->mother()) {
						spouse = offspring[j]->mother();
						nsp++;
					}
					j++;
				}
				I->numspouse = nsp;
				if(nsp>0){
					I->numoffs_spouse = new unsigned [nsp];
				}
				else {
					I->numoffs_spouse = 0;
				}
				if(nsp>0){
					I->spouselist = new Individual* [nsp];
				}
				else {
					I->spouselist = 0;
				}
				spouse = offspring[0]->mother();
				I->spouselist[0] = spouse;
				j = 1; s = 0; k = 1;
				while(j<noffs) {
					if (spouse != offspring[j]->mother()) {
						spouse = offspring[j]->mother();
						I->numoffs_spouse[s++] = k;
						I->spouselist[s] = spouse;
						k = 1;
					}
					else {
						k++;         // count # of smyoffspring of each spouses
					}
					j++;
				}
				I->numoffs_spouse[s] = k;
			}
		}
	}
	mark_need_remove = 1;
	spouse_info_built = 1;
}

void Population::renum(void)
{
	if (stdid) return;
	unsigned i,k,newid,oldid;
	for (i=0; i<popsize; i++) popmember[i]->group_id = 1;
	
	Individual *I;
	unsigned j = 1;
	int g = 0;
	while (j) {
		j = 0;
		for (i=0; i<popsize; i++) {
			I = popmember[i];
			k = std::max(I->mother_gid(),I->father_gid());
			if (I->gid() <= k) {
				I->group_id = k+1;
				j = i+1;
			}
		}
		if (g++ >150) {
			std::cout << " A B C" << std::endl;
			std::cout << " B A D" << std::endl;
			std::cout << "is it strange ?  but it's in your population " <<  std::endl;
			std::cout << "check the families with " << ind_name(j) << std::endl;
			exit(1);
		}
	}
	qsort(popmember,popsize,sizeof(Individual*),compare_ind_gid);
	for (i=0; i<popsize; i++) {
		I = popmember[i];
		newid = i+1;
		oldid = I->id();
		hashtable.change_id(oldid,newid);
		I->reset_id(newid);
	}
	hashtable.reorder();  // this is necessary because change_id(i,j) is used
	stdid = 1;
	remove_mark();
}

void Population::release(void)
{
	if (nufamily_vec) {
		for (unsigned i=0; i<num_nufamily; i++) {
			if(nufamily_vec[i]){
				delete nufamily_vec[i];
				nufamily_vec[i]=0;
			}
		}
		if(nufamily_vec){
			delete [] nufamily_vec;
			nufamily_vec=0;
		}
		nufamily_vec = 0;
	}
	if (popmember)    {delete [] popmember;    popmember = 0;}
	if (trans_mat)    {delete [] trans_mat;    trans_mat = 0; }
	//BRS This is line is trying to delete something that no longer exists.
	//It appears to be destroyed when deleting nufamilies in first IF (nufamily_vec)
	if (pm_storage)   {delete [] pm_storage;   pm_storage = 0;} 
	// I have uncommented the line above it is actually needed LRT (2/12/04)
	if (pop_gamete)   {delete [] pop_gamete;   pop_gamete = 0;}
	if (gametebase)   {delete [] gametebase;   gametebase = 0;}
	if (markersymbol) {delete [] markersymbol; markersymbol = 0;}
	// LRT (2/17/2004)
	gNodeList.releaseGNsts();
}

Individual* Population::member(unsigned k)
{
	if ( k >= popsize ) throw exception(" Population::member(), bad arg");
	if (!stdid) renum();
	return popmember[k];
}

Population* Population::sub(const unsigned subsize)
{
	Population* tmp = new Population;
	check_ptr(tmp);
	tmp->resize(subsize,prior);
	tmp->hashtable.resize(subsize);
	unsigned i,id,tpsize = 0;
	Individual *I;
	for (i=0; i<popsize;i++) {
		I = popmember[i];
		if (I->group_id) {
			tmp->hashtable.insert(ind_name(i+1));
			tmp->pm_storage[tpsize++] = *I;
		}
	}
	if (tpsize != subsize) throw exception(" Population::sub(): subsize is too small");
	tmp->popsize = tpsize;
	
	for (tpsize=0,i=0; i<popsize; i++) {
		I = popmember[i];
		if (I->group_id) {
			id = 0;
			if (I->mymother) id = tmp->hashtable.get_id(ind_name(I->mother_id()));
			if (id>0) {
				tmp->pm_storage[tpsize].mymother = &(tmp->pm_storage[id-1]);
			}
			else {
				tmp->pm_storage[tpsize].mymother = 0;
			}
			id = 0;
			if (I->myfather) id = tmp->hashtable.get_id(ind_name(I->father_id()));
			if (id>0) {
				tmp->pm_storage[tpsize].myfather = &(tmp->pm_storage[id-1]);
			}
			else {
				tmp->pm_storage[tpsize].myfather = 0;
			}
			tpsize++;
		}
	}
	tmp->fdone = fdone;
	tmp->renum();
	tmp->build_offs_info();
	return tmp;
}

void Population::gtindex(const unsigned num,const unsigned ni,const Vector<unsigned> &ngt,Vector<unsigned> &gtvec)
{
	////////////////////////////////////////////////////////////////////////
	// suppose, 3 individuals, the 1st and 2nd can have 2 genotypes,
	// the 3nd can have 3 genotypes. thus ngt = {2,2,3}.
	// I want num be transformed into gtvec system according to ngt:
	//
	// num     gtvec[0] gtvec[1] gtvec[2]
	// 0  ->    0        0        0
	// 1  ->    0        0        1
	// 2  ->    0        0        2
	// 3  ->    0        1        0
	// 4  ->    0        1        1
	// 5  ->    0        1        2
	// 6  ->    1        0        0
	// 7  ->    1        0        1
	// 8  ->    1        0        2
	// 9  ->    1        1        0
	// 10 ->    1        1        1
	// 11 ->    1        1        2
	/////////////////////////////////////////////////////////////////////
{
	long i;
	unsigned ir=num;
	for (i=ni-1; i>=0; i--) {
		gtvec[i] = static_cast<unsigned>(fmod(static_cast<double>(ir),static_cast<double>(ngt[i])));
		ir /= ngt[i];
	}
}
}

void Population::display(const char key[])
{
	Individual *I;
	unsigned i,t,t1,t2,t3,k,nl,num,ii;
	if (strcmp(key,"ped") == 0) {
		for (i=0; i<popsize; i++) {
			I = popmember[i];
			std::cout << ind_name(i+1) << " " << I->id() << " " << I->mother_id();
			std::cout  << " " << I->father_id() << " " << I->sex() << "\n";
		}
	}
	else if(strcmp(key,"genotype") == 0) {
		if (!(popmember[0]->genotype_counter)) return;
		Chromosome *C1, *C2;
		Vector<double> *vv;
		for (i=0; i<popsize; i++) {
			I = popmember[i];
			std::cout << "individual " << ind_name(i+1) << ":\n";
			for (t=0; t<numchrom; t++) {
				vv = &(I->genotype_counter[t]);
				std::cout << "  chromosome " << t+1 << ":\n";
				num = pop_gamete[t].size();
				nl = pop_gamete[t].chromosome[0].nloci();
				for (k=0,t1=0; t1<num; t1++) {
					C1 = &(pop_gamete[t].chromosome[t1]);
					for (t2=0; t2<=t1; t2++) {
						C2 = &(pop_gamete[t].chromosome[t2]);
						std::cout << "    ___ ";
						for (t3=0; t3<nl; t3++) std::cout << C1->locus[t3].allele << " ";
						std::cout << "___\n";
						std::cout << "    ~~~ ";
						for (t3=0; t3<nl; t3++) std::cout << C2->locus[t3].allele << " ";
						std::cout << "~~~: "<< static_cast<unsigned>(vv->begin()[k]) << "\n";
						k++;                          // k = t1*(t1+1)/2 + t2;
					}
				}
			}
		}
	}
	else if(strcmp(key,"marker") == 0) {
		for (ii=0; ii<popsize; ii++) {
			I = popmember[ii];
			std::cout << "individual " << ind_name(ii+1) << ":\n";
			for (unsigned i=0; i<prior->nloci_chrom(1); i++){
				if(prior->chrom()[0].locus[i].qtl_ml=='m'){
					std::cout << I->genome0.chromosome[0].locus[i].allele;
				}
			}
			cout << endl;
			for (unsigned i=0; i<prior->nloci_chrom(1); i++){
				if(prior->chrom()[0].locus[i].qtl_ml=='m'){
					std::cout << I->genome1.chromosome[0].locus[i].allele;
				}
			}
			for (unsigned i=0;i<numtrait;i++){
				std::cout << setw(8) 
				<< setprecision (3) 
				<< setiosflags (ios::right | ios::fixed)
				<< I->myrecord[i].double_val();
			}
			cout << endl;
		}
	}	
	std::cout << "\n";
	release_genotype_counter();
}
//RLF


unsigned Population::input_descentGraph(char *dgfile){
	// Authors: Rohan L. Fernando and Fabiano V. Pita
	// (2003) 
	// Contributors: L. Radu Totir, Chris Stricker
	
	std::string strid;
	int id, intJunk;
	Individual *I;
	Vector <unsigned> m_g , p_g;
	unsigned nLoci = prior->nloci_chrom(1);
	allele_vector1.resize(nLoci); //need to put this somewhere before building
	allele_vector2.resize(nLoci); //allele vectors. Need to do this only once.
	previousAlleleVector1.resize(nLoci);
	previousAlleleVector2.resize(nLoci);
	connected_groups.resize(nLoci);
	previousConnectedGroups.resize(nLoci);
	connect_counter.resize(nLoci);
	previousConnectCounter.resize(nLoci);
	founder_allele_neighbors.resize(nLoci);
	founder_allele_neighbors_all.resize(nLoci);
	previousFounderAlleleNeighbors.resize(nLoci);
	m_g.resize(nLoci);
	p_g.resize(nLoci);
	std::ifstream mfile(dgfile);
	if(!mfile) {
		throw exception("Couldn't open file ");
	}
	for(int i=0; i<popsize;i++){
		//    if (!(mfile >> strid)){
		if (!(mfile >> id)){ //this is for standard pedigrees only
			throw exception(" Not enough rows in dgfile ");
		}
		//    id  = hashtable.get_id(strid.c_str()); //this needs to be uncommented if pedigree is non-standard!
		I = popmember[id-1];
		if (I->mymother==0){
			for(int j=1; j<=nLoci;j++){
				if(!(mfile >> intJunk >> intJunk )){
					std::cout << " Not enough columns in dgfile " << std::endl;
				}
			}
			continue;
		}
		for(int j=1; j<=nLoci;j++){
			if(!(mfile >> m_g(j) )){
				std::cout << " Not enough columns in dgfile " << std::endl;
			}
			if(!(mfile >> p_g(j))){
				throw exception(" Not enough columns in dgfile ");
			}
		}
		I->put_gametes(m_g, p_g);
	}
	//Initialize Pop->founder_allele_counter
	//Determine and set founder alleles, 0 for not updating
	founder_allele_counter = 0;
	for (int i=0;i<popsize;i++) {
		member(i)->set_founder_alleles(0);
		cout << "." << flush;
	}
	cout << endl;
	for(unsigned lcs=1; lcs<=nLoci; lcs++){
		previousAlleleVector1(lcs).resize(founder_allele_counter,0);  // Here we can allocate space for the second dimension of the vector 
		previousAlleleVector2(lcs).resize(founder_allele_counter,0);  // to store the previously sampled allele vector, we want to do this 
																	  // once and not redo it for every DG, as the vector shall not be 
																	  // resized and thus initialized for every new DG as it would loose 
																	  // the values stored in it that we need to restore the previous values.
		previousConnectedGroups(lcs).resize(founder_allele_counter,0);
		previousFounderAlleleNeighbors[lcs-1].resize(founder_allele_counter);
		
  	}
  	return 1;
	}


void Population::output_descentGraph(char *dgfile){
	// Authors: Rohan L. Fernando and Fabiano V. Pita
	// (2003) 
	// Contributors: L. Radu Totir
	std::ofstream mfile(dgfile, std::ios::out);//SDK |std::ios::noreplace);
	if(!mfile) {
		throw exception("Couldn't open file ");
	}
	unsigned nLoci = prior->nloci_chrom(1);
	for(int i=0; i<popsize;i++) {
		mfile << ind_name(i+1);
		for(int j=1; j<=nLoci;j++){
			mfile << " " << member(i)->m_gamete(j)
			      << " " << member(i)->p_gamete(j);
		}
		mfile << std::endl;
	}
	return;
}

void Population::output_descentGraphJP(string dgfile){
	// Authors: Rohan L. Fernando 
	// (2003) 
	// Contributors:
	// this version is for when joint peeling is used and allele origin
	// indicators are directly sampled
	
	std::ofstream mfile(dgfile.c_str(), std::ios::out);//SDK |std::ios::noreplace);
	if(!mfile) {
		throw exception("Couldn't open file ");
	}
	unsigned nLoci = prior->nloci_chrom(1);
	for(int i=0; i<popsize;i++) {
		mfile << ind_name(i+1);
		Individual *mom = member(i)->mymother;
		if (mom){
			for(int j=0; j<nLoci; j++){
				unsigned matOrigin = member(i)->malleleOriginNodeVector[j].getAcceptedAlleleOrigin();
				unsigned patOrigin = member(i)->palleleOriginNodeVector[j].getAcceptedAlleleOrigin();
				mfile << " " << matOrigin
					  << " " << patOrigin;
			}
		}
		else {
			for(int j=0; j<nLoci; j++){
				mfile << " -9999 -9999";
			}
		}
		mfile << std::endl;
	}
	return;
}

unsigned Population::input_markerData(Data* D)
{
	// Authors: Rohan L. Fernando 
	// (???, 1996) 
	// Contributors: LRT, Chris Stricker
	

	if (markerDataIn) return D->num_rows();
	if (!D || D->num_rows() <= 0) throw exception(" Population::input_data(): data is empty");
	unsigned numcol = D->num_cols();
	unsigned numrec = D->num_rows();
	unsigned i,j,t,id;
	const char *strid;
	
	if (!D->in_memory()) D->input_datasheet();
	for (t=0; t<numcol; t++) {
		if (D->datasheet[t].col_struct.classi()=='G') break;
	}
	if (t==numcol) throw exception("Population::input_data(): no pedigree exits");
	char TY = D->datasheet[t].type();
	DataNode  *tcol, *column = D->rawcol(t);
	HashTable* dhtable = D->hashtable[t];
	Vector<DataNode*> markerCols1, markerCols2; 
	std::string nm1,nm2;
	
	// number of marker loci is assumed to be = number of loci - 1
	// locus(0) is the qtl and the remaining loci are all markers
	//    n_markerLoci = prior->nloci_chrom(1) - 1;
	
	// I (LRT) have replaced the line above with the following:
	n_markerLoci = prior->numMarkerLoci; 
	std::cout<<"in input_MarkerData(), numMarkerLoci = "<<prior->numMarkerLoci<<std::endl;
	markerCols1.resize(prior->nloci_chrom(1));
	markerCols2.resize(prior->nloci_chrom(1));
	
	// find marker columns in data sheet for each marker locus
	// pointers to columns in data sheet are stored in
	// markerCols1 and markerCols2
	
	for (i=0; i<prior->nloci_chrom(1); i++){
		if(prior->chrom()[0].locus[i].qtl_ml=='m'){
			nm1 = prior->chrom()[0].locus[i].nameOfcol1;
			nm2 = prior->chrom()[0].locus[i].nameOfcol2;
			
			for (t=0; t<numcol; t++) { // column 1
				if (nm1 == D->datasheet[t].col_struct.name()) break;
			}
			if (t==numcol) throw exception("Population::input_data(): no data for nameOfcol1 locus i");
			markerCols1(i+1) = D->rawcol(t);
			
			for (t=0; t<numcol; t++) { // column 2
				if (nm2 == D->datasheet[t].col_struct.name()) break;
			}
			if (t==numcol) throw exception("Population::input_data(): no data for nameOfcol2 locus i");
			markerCols2(i+1) = D->rawcol(t);
		}
	}
	Individual *I;
	int gtype[2];
	
	// Now we transfer the marker data from the data sheet to Individuals
    
	if ((TY == 'S')  || (Pedigree::type==Pedigree::raw)){ 		
		for (i=0; i<numrec; i++) {
			id = column[i].unsigned_val();
			strid = (const char *)dhtable->find(id);
			id  = hashtable.get_id(strid);
			I = popmember[id-1];
			for (j=0; j<prior->nloci_chrom(1); j++){
				if(prior->chrom()[0].locus[j].qtl_ml=='m'){
					int alleleState1 = (int) markerCols1(j+1)[i].double_val();
					int alleleState2 = (int) markerCols2(j+1)[i].double_val();
					if (alleleState1==0 || alleleState2==0){
						alleleState1 = 0;
						alleleState2 = 0;
					}
					if (alleleState1 <= prior->chrom()[0].locus[j].nallele() ){
						I->genome0.chromosome[0].locus[j].allele = alleleState1;
					}
					else {
						std::cerr << "In record " << i+1 << ", allele state " << alleleState1 << " for locus " << j+1 << " is larger than " << prior->chrom()[0].locus[j].nallele() << ", the total number of alleles for that locus."<< std::endl;	
						throw exception("Population::input_markerData(): Allele state at marker locus exceeds total number of alleles for that locus");
					}
					if (alleleState2 <= prior->chrom()[0].locus[j].nallele() ){
						I->genome1.chromosome[0].locus[j].allele = alleleState2;
					}
					else {
						std::cerr << "In record " << i+1 << ", allele state " << alleleState2 << " for locus " << j+1 << " is larger than " << prior->chrom()[0].locus[j].nallele() << ", the total number of alleles for that locus."<< std::endl;	
						throw exception("Population::input_markerData(): Allele state at marker locus exceeds total number of alleles for that locus");
					}
				}
			}
		}
	}
	else {
		for (i=0; i<numrec; i++) {
			id = (unsigned) column[i].double_val();
			I = popmember[id-1];
			for (j=0; j<prior->nloci_chrom(1); j++){
				if(prior->chrom()[0].locus[j].qtl_ml=='m'){
					int alleleState = (int) markerCols1(j+1)[i].double_val();
					if (alleleState <= prior->chrom()[0].locus[j].nallele()){
						I->genome0.chromosome[0].locus[j].allele = (int) markerCols1(j+1)[i].double_val();
						I->genome1.chromosome[0].locus[j].allele = (int) markerCols2(j+1)[i].double_val();
					}
					else {
						std::cerr << "In record " << i+1 << ", allele state " << alleleState << " for locus " << j+1 << " is larger than " << prior->chrom()[0].locus[j].nallele() << ", the total number of alleles for that locus."<< std::endl;	
						throw exception("Population::input_markerData(): Allele state at marker locus exceeds total number of alleles for that locus");
					}
				}
			}
		}
	}
	if ((TY == 'S')  && (Pedigree::type!=Pedigree::raw)){ 	
		std::cerr<<"Warning: Pedigree type is standard, but pedigree variable is of type string."<<std::endl;
	}
	else if ((TY != 'S')  && (Pedigree::type==Pedigree::raw)){ 		
		std::cerr<<"Warning: Pedigree type is raw, but pedigree variable is of not of type string."<<std::endl;
	}
	
	for (j=0, t=0; t<numcol; t++) {
		if (D->datasheet[t].col_struct.classi()=='T') j++;
	}
	if (j != numtrait) throw exception("Population::input_Markerdata(): Pop's ntrait != Data's ntrait");
	for (j=0, t=1; t<numcol; t++) {   // 1st column is reserved for intercept
		if (D->datasheet[t].col_struct.classi()=='T') {
			tcol =  D->rawcol(t);
			for (i=0; i<numrec; i++) {
				if (TY == 'S') {
					id = column[i].unsigned_val();
					strid = (const char *)dhtable->find(id);
					id  = hashtable.get_id(strid);
				}
				else{
					id = (unsigned) column[i].double_val();
				}
				popmember[id-1]->myrecord[j] = tcol[i];
			}
			j++;
		}
	}
	D->release_datasheet();
	if (j == 0) throw exception("Population::input_data(): trait(s) not available");
	//std::cout << "exit from input marker data \n";
	markerDataIn = true;
	return numrec;
}


void Population::input_markerData(MIM& M)
{
	// Authors: Rohan L. Fernando 
	// (August, 2005) 
	// Contributors: 
	
	if (markerDataIn) return;
	ifstream datafile;
    datafile.open(M.markDataFileName.c_str());
    if(!datafile) {
      cerr << "Couldn't open marker data file: " << M.markDataFileName << endl;
      throw exception(" in Population::inputMarkerData(MIM)");
    }
	Vector<int> colIndex0, colIndex1;
	colIndex0.resize(prior->nloci_chrom(1));
	colIndex1.resize(prior->nloci_chrom(1));
	n_markerLoci = prior->numMarkerLoci; 
	// find the colum indices for the markers
	// and store them in colIndex0 and colIndex1
	string nm0, nm1;
	for (unsigned i=0; i<prior->nloci_chrom(1); i++){
		if(prior->chrom()[0].locus[i].qtl_ml=='m'){
			std::string nm0 = prior->chrom()[0].locus[i].nameOfcol1;
			std::string nm1 = prior->chrom()[0].locus[i].nameOfcol2;
			int index0 = M.markerColName.getIndex(nm0);
			int index1 = M.markerColName.getIndex(nm1);
			if (index0 < 0) {
				cout << nm0 << " not found in marker column names" << endl;
				throw exception("in Population::input_markerData(MIM)");
			}
			colIndex0[i] = index0;
			if (index1 < 0) {
				cout << nm1 << " not found in marker column names" << endl;
				throw exception("in Population::input_markerData(MIM)");
			}
			colIndex1[i] = index1;
		}
	}
	Individual *I;
	string sep(" \t");
	Tokenizer colData;
	unsigned rec=0, rec1=0;
	unsigned numberTokens = M.markerColName.size();
	std::string inputStr;
	// Now we transfer the marker data to Individuals
	cout << "\nReading marker data ... \n";
	while (getline(datafile,inputStr)){		
		rec++;
		if(rec==1000){
			cout << "." ;
			rec1 += rec;
			rec = 0;
		}
		colData.getTokens(inputStr,sep);
		if (numberTokens != colData.size()){
			cerr << "Expected " << numberTokens 
			<< " columns in marker data file, but found "
			<< colData.size() 
			<< " for record "  
			<< (rec+rec1) << endl;
			throw exception("in Population::input_markerData(MIM)");
		}
		// find the individual first
		string indStrId = colData[0];
		RPedigree::iterator pedIt = myRPedPtr->find(indStrId);
		if (pedIt == myRPedPtr->end()) {
			cerr << indStrId << " in record " << (rec1+rec) 
			<< " of " << M.markDataFileName 
			<<" not found in pedigree file \n";
			throw exception("in Population::input_markerData(MIM)");
		}
		unsigned indIntId = pedIt->second->ind;
		I = popmember[indIntId-1]; 
		// now transfer marker data
		int alleleState0, alleleState1; 
		for (unsigned j=0; j<prior->nloci_chrom(1); j++){
			if(prior->chrom()[0].locus[j].qtl_ml=='m'){
			    unsigned markIndex0 = colIndex0[j];
				unsigned markIndex1 = colIndex1[j];
				std::string strAllele0 = colData[markIndex0];
				std::string strAllele1 = colData[markIndex1];
				if (strAllele0=="."){
					alleleState0 = 0;
				}
				else {
					alleleState0 = MIM::getInteger(strAllele0);
				}
				if (strAllele1=="."){
					alleleState1 = 0;
				}
				else {
					alleleState1 = MIM::getInteger(strAllele1);
				}
				if (alleleState0==0 || alleleState1==0){
					alleleState0 = 0;
					alleleState1 = 0;
				}
				if (alleleState0 <= prior->chrom()[0].locus[j].nallele()){
					I->genome0.chromosome[0].locus[j].allele = alleleState0;
				}
				else {
					std::cerr << "In record " 
							  << (rec1+rec) 
							  << ", allele state " 
							  << alleleState0 
					          << " for locus " 
							  << j+1 
							  << " is larger than " 
							  << prior->chrom()[0].locus[j].nallele() 
							  << ", the total number of alleles for that locus."<< std::endl;	
					throw exception("in Population::input_markerData(MIM)");
				}
				if (alleleState1 <= prior->chrom()[0].locus[j].nallele()){
					I->genome1.chromosome[0].locus[j].allele = alleleState1;
				}
				else {
					std::cerr << "In record " 
							  << (rec1+rec) 
							  << ", allele state " 
							  << alleleState1 
					          << " for locus " 
							  << j+1 
							  << " is larger than " 
							  << prior->chrom()[0].locus[j].nallele() 
							  << ", the total number of alleles for that locus."<< std::endl;	
					throw exception("in Population::input_markerData(MIM)");
				}
			}
		}
	}
	markerDataIn = true;
}

void Population::displayPDQ(char *fname){
	std::ofstream outfile(fname);
	Individual *ind;
	for (int i=0; i<popsize; i++) {
		ind=popmember[i];
		if (!ind->base()) {
			outfile <<setw(10)  <<  ind_name(i+1)  <<" ";
			outfile <<setw(10)  << 1.0 - ind->p_counter/double(ind->genotNodeVector[0].sampleCount) << " ";
			outfile <<setw(10)  << 1.0 - ind->m_counter/double(ind->genotNodeVector[0].sampleCount) << endl;
		}
	}
	outfile.close();
}
/*! \fn void Population::displayPDQ(char *fname)
*   \brief method for outputting PDQ's for aviagen project
*/


void Population::output_pdq(int n_samples, char *fname){
	std::ofstream outfile(fname);
	if(Pedigree::type==Pedigree::raw){ 		
		for (int i=0; i<popsize; i++) {
			ostringstream ostr;
			unsigned k = i+1;
			ostr <<  k;
			std::string t = ostr.str();
			int id  = hashtable.get_id(t.c_str());
			if (member(id-1)->base()) {
				outfile <<  ind_name(id)  <<" ";
				outfile << "0.0" << " " << "0.0" << " " << "0.0" << " " << "0.0" << " "
				<< "0.0" << " " << "0.0" << " " << "0.0" << " " << "0.0" << " " << std::endl;
			}
			else {
				outfile <<  ind_name(id)  <<" ";
				outfile << 1.0-(member(id-1)->m_counter/float(n_samples))
				<< " "
				<< member(id-1)->m_counter/float(n_samples)
				<< " "
				<< "0.0"
				<< " "
				<< "0.0"
				<< " "
				<< "0.0"
				<< " "
				<< "0.0"
				<< " "
				<< 1.0-(member(id-1)->p_counter/float(n_samples))
				<< " "
				<< member(id-1)->p_counter/float(n_samples)
				<<std::endl;
			}
		}
	}
	else{
		for (int i=0; i<popsize; i++) {
			if (member(i)->base()) {
				outfile << "0.0" << " " << "0.0" << " " << "0.0" << " " << "0.0" << " "
				<< "0.0" << " " << "0.0" << " " << "0.0" << " " << "0.0" << " " << std::endl;
			}
			else {
				outfile << 1.0-(member(i)->m_counter/float(n_samples))
				<< " "
				<< member(i)->m_counter/float(n_samples)
				<< " "
				<< "0.0"
				<< " "
				<< "0.0"
				<< " "
				<< "0.0"
				<< " "
				<< "0.0"
				<< " "
				<< 1.0-(member(i)->p_counter/float(n_samples))
				<< " "
				<< member(i)->p_counter/float(n_samples)
				<<std::endl;
			}
		}
	}
	outfile.close();
}

void Population::output_pdq(unsigned locus, int n_samples, char *fname){
	// Authors: Rohan L. Fernando (2005)
	// Contributors: Chris Stricker (2005)
	std::ofstream outfile(fname);
	for (int i=0; i<popsize; i++) {
		outfile << ind_name(i+1) << " " 
//		<< member(i)->myid << " "
		;
		if (member(i)->base()) {
			outfile << "-1"  << " " << "-1"  << " " << "0.0" << " " << "0.0" << " "
			        << "0.0" << " " << "0.0" << " " << "-1"  << " " << "-1"  << " " << std::endl;
		}
		else {
			//outfile <<  i+1  <<" "
			outfile << 1.0-(member(i)->m_counter_map(locus)/float(n_samples))
			<< " "
			<< member(i)->m_counter_map(locus)/float(n_samples)
			<< " "
			<< "0.0"
			<< " "
			<< "0.0"
			<< " "
			<< "0.0"
			<< " "
			<< "0.0"
			<< " "
			<< 1.0-(member(i)->p_counter_map(locus)/float(n_samples))
			<< " "
			<< member(i)->p_counter_map(locus)/float(n_samples)
			<<std::endl;
		}
	}
	outfile.close();
}



// 10/11/02 RLF
void Population::output_ibdMatrix(unsigned n_samples, char *fname){
	std::ofstream outfile(fname);
	for (unsigned i=1; i<=popsize; i++) {  
		int id = member(i-1)->id();
		for (unsigned j=i; j<=popsize; j++) {
			int jd = member(j-1)->id();
			outfile << std::setw(8) << ind_name(id) << std::setw(8) << ind_name(jd) << std::setw(15) << ibdMatrix(i,j)/(2.0*n_samples) << endl;
		}
	}
}
//RLF

//BRS

void Population::set_switches(void)
{
	int i;
	cout <<"Switches on "<< n_markerLoci<<endl;
	for (i=0; i<popsize; i++){
		popmember[i]->set_switch(n_markerLoci);
	}
}
//BRS

double Population::disPenetranceTable[2][3] =
{ 1.00, 1.00, 0.00,
	0.00, 0.00, 1.00};

void Population::setupGNodeSampler(RPedigree& P, MIM& M, GeneticDist& G){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2005) 
	// Contributors:
	input_ped(P,G);
	input_markerData(M);
	GNodeStructureSetup();
}

void Population::setupRSampler(Pedigree& P, Data* D, GeneticDist *G){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	input_ped(P,G);
	model->fitdata(*D);
	model->prepare_data();
	input_markerData(D);
	GNodeStructureSetup();
}

void Population::GNodeStructureSetup(void){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors:
	cout << "\nGNode structure setup \n";
	
	Individual *ind;
	GNodeList::popPtr = this;
	GNodeSet::prior=prior;
	// next line is necessary for when we sample locus positions
	prior->chrom()[0].calcLocusMask();
	Individual::numLoci = prior->chrom()[0].nloci();
	Individual::nQTL=0;
	numHaplotypes=1;
	SafeSTLVector<unsigned> alleleCountVector;
	alleleCountVector.name = "Population::GNodeStructureSetup::AlleleCountVector";
	alleleCountVector.resize(Individual::numLoci);
	for(unsigned i=0;i<Individual::numLoci;i++){
		if(prior->chrom()[0].locus[i].qtl_ml=='q'){
			Individual::QTLPosVector.push_back(i);
			Individual::nQTL++;
		}
		alleleCountVector[i]=prior->chrom()[0].locus[i].nallele();
		numHaplotypes*=prior->chrom()[0].locus[i].nallele();
	}
	haplotypeCoder.setupVectors(alleleCountVector);
	buildRecombinationMatrix();
	unsigned rec = 0,rec1=0;
	cout << "\n Initializing GNodes for genotype, allele state and allele origin \n";
	for(unsigned i=0;i<popsize;i++){
		rec++;
		if(rec==1000){
			cout<<rec1+rec <<"\r";			
			cout.flush();
			rec1+=rec;
			rec = 0;
		}
		ind=popmember[i];
		ind->m_gamete.resize(Individual::numLoci,9999);
		ind->p_gamete.resize(Individual::numLoci,9999);
		ind->m_gameteOld.resize(Individual::numLoci,9999);
		ind->p_gameteOld.resize(Individual::numLoci,9999);
		ind->residual_var = &model->residual_var;
		ind->genotNodeVector.resize(Individual::numLoci);
		ind->malleleStateNodeVector.resize(Individual::numLoci);
		ind->palleleStateNodeVector.resize(Individual::numLoci);   
		if(ind->mymother){
			ind->malleleOriginNodeVector.resize(Individual::numLoci);
			ind->palleleOriginNodeVector.resize(Individual::numLoci);   
		}
		ind->setOwnerGNodes(); 
	}
	if (!nuFamiliesDone) {
		build_nufamily();
	}
	for(unsigned j=0;j<Individual::numLoci;j++){
		SimpleGenotypeElimination(j);
		if (prior->chrom()[0].locus[j].qtl_ml=='m' && num_nufamily != 0){
			LGGenotypeElimination(j);
		}
		setAlleleStateVectors(j);
	} 
}
/*! \fn void Population::completeGNodeStructureSetup(void)
*  \brief method to setup the structures needed by the GNode-sampler
*/

void Population::SimpleGenotypeElimination(unsigned locus){
	// Authors: Rohan L. Fernando
	// (October, 2004) 
	// Contributors: 	
	Individual *ind;
	unsigned j = locus;
	unsigned rec = 0,rec1=0;
	cout << "\n Simple genotype elimination \n";
	for(unsigned i=0;i<popsize;i++){
		rec++;
		if(rec==1000){
			cout<<rec1+rec <<"\r";			
			cout.flush();
			rec1+=rec;
			rec = 0;
		}

		ind=popmember[i];
		ind->genotNodeVector[j].sampled = false;
		ind->malleleStateNodeVector[j].sampled = false;
		ind->palleleStateNodeVector[j].sampled = false;
		ind->malleleStateNodeVector[j].alleleStateVector.clear();
		ind->palleleStateNodeVector[j].alleleStateVector.clear();
		
		unsigned allelePat = ind->genome0.chromosome[0].locus[j].allele;
		unsigned alleleMat = ind->genome1.chromosome[0].locus[j].allele;
		if(alleleMat==0){ // maternal allele is missing 
			ind->genotNodeVector[j].uninformativeGNode = 1; // missing genotype; may be pruned
			unsigned numAllele = prior->chrom()[0].locus[j].nallele();
			for (unsigned k=1;k<=numAllele;k++){
				ind->malleleStateNodeVector[j].alleleStateVector.push_back(k);
			}
		}
		else if(alleleMat==allelePat){
			ind->genotNodeVector[j].uninformativeGNode = 0;
			ind->malleleStateNodeVector[j].alleleStateVector.push_back(alleleMat);
			ind->malleleStateNodeVector[j].alleleState = 0; // gives position in alleleStateVector of current state 
			ind->palleleStateNodeVector[j].alleleState = 0;
			ind->genotNodeVector[j].genotypeState = 0;
		}
		else {
			ind->genotNodeVector[j].uninformativeGNode = 0;
			ind->malleleStateNodeVector[j].alleleStateVector.push_back(allelePat);
			ind->malleleStateNodeVector[j].alleleStateVector.push_back(alleleMat);
		}
		if(ind->mymother){
			ind->malleleOriginNodeVector[j].alleleOriginVector.push_back(0);
			ind->malleleOriginNodeVector[j].alleleOriginVector.push_back(1);
			ind->palleleOriginNodeVector[j].alleleOriginVector.push_back(0);
			ind->palleleOriginNodeVector[j].alleleOriginVector.push_back(1);
		}
	}
	finishUpAlleleStateInit(locus);
}

/*! \fn void Population::SimpleGenotypeElimination(unsigned locus)
*  \brief eliminates genotypes using the individuals phenotypes
*/

void Population::LGGenotypeElimination(unsigned locus){
	// Authors: Joseph Abraham and Rohan L. Fernando
	// (September, 2004) 
	// Contributors: 	
	if (!nuFamiliesDone) {
		throw exception ("Population::LGGenotypeElimination: nuFamilies not built");
	}
	unsigned totalNumberOfGenotypes = calcTotalNumberOfGenotypes(locus);
	unsigned oldTotal, rec=0, rec1=0;
	do {
	    rec=rec1=0;
		cout <<"\nIn LGGenotypeElimination: Locus: " <<locus<<" Total number of genotypes = "<<totalNumberOfGenotypes << endl;
		if(!prior->chrom()[0].locus[locus].qtl_ml=='m') break; 
		oldTotal = totalNumberOfGenotypes;
		for (unsigned i=0;i<num_nufamily;i++){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1+=rec;
				rec = 0;
			}
			
			nufamily_vec[i]->eliminateGenotypes(locus);
		}
		totalNumberOfGenotypes = calcTotalNumberOfGenotypes(locus);		
	}
	while (oldTotal>totalNumberOfGenotypes);
}
/*! \fn void Population::LGGenotypeElimination(unsigned locus)
*  \brief method to do Lange and Goradia genotype elimination
*/

void Population::setAlleleStateVectors(unsigned locus){
	// Authors: Joseph Abraham 
	// (2004) 
	// Contributors: 
	for(unsigned i=0;i<popsize;i++){
		popmember[i]->setAlleleStateVectors(locus);
		//popmember[i]->displayAlleleVectors(locus);
	}
	for(unsigned i=0;i<popsize;i++){
		if (popmember[i]->mymother) popmember[i]->setAlleleOriginVectors(locus);
		//popmember[i]->displayAlleleVectors(locus);
	}
}



unsigned Population::calcTotalNumberOfGenotypes(unsigned locus){
	// Authors: Joseph Abraham and Rohan L. Fernando
	// (September, 2004) 
	// Contributors: 
	unsigned count = 0;
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		//ind->displayGenotypes(locus);
		unsigned numGenotypes = ind->genotNodeVector[locus].genotypeVector.size();
		if (numGenotypes==0){
			throw exception ("Population::calcTotalNumberOfGenotypes: Incompatible genotypes in pedigree");
		}
		count += numGenotypes;
	}
	return count;
}
/*! \fn unsigned Population::calcTotalNumberOfGenotypes(unsigned locus)
*  \brief method to do count total number of possible genotypes across 
*         the whole pedigree
*/


void Population::buildRecombinationMatrix(void){
	// Authors: L. Radu Totir and Rohan L. Fernando
	// (September, 2003) 
	// Contributors: 
	unsigned nLoci = prior->chrom()[0].nloci();
	recombinationMatrix.resize(nLoci,nLoci,0.0);
	for(unsigned i=1; i<nLoci; i++) {
		for (unsigned j=1; j<nLoci; j++) {
			if(j>=i){
				//recombinationMatrix[i-1][j] = model->MapF(std::abs(prior->get_distance(1,j+1) - prior->get_distance(1,i)));
				recombinationMatrix[i-1][j] = GeneticDist::HaldaneMToR(std::abs(prior->get_distance(1,j+1) - prior->get_distance(1,i)));
				//recombinationMatrix[i-1][j] = 0.5;
			}
		}
	}
	// std::cout << "Recombination Matrix:\n" << recombinationMatrix;
}
/*! \fn unsigned Population::buildRecombinationMatrix(void)               
*  \brief method to compute the recombintion rate between loci  
*                                             
*/


void Population::finishUpAlleleStateInit(unsigned locus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors: 
	Individual *ind;
	MaternalPaternalAlleles matPat;
	unsigned j = locus;
	for(unsigned i=0;i<popsize;i++) { 
		ind=popmember[i];
		ind->genotNodeVector[j].genotypeVector.clear();
		unsigned allelePat = ind->genome0.chromosome[0].locus[j].allele;
		unsigned alleleMat = ind->genome1.chromosome[0].locus[j].allele;
		unsigned n = ind->malleleStateNodeVector[j].alleleStateVector.size();
		for(unsigned k=0;k<n;k++){
			unsigned allelek = ind->malleleStateNodeVector[j].alleleStateVector[k];
			ind->palleleStateNodeVector[j].alleleStateVector.push_back(allelek);
			matPat.maternal = allelek-1;
			for(unsigned l=0;l<n;l++){
				unsigned allelel = ind->malleleStateNodeVector[j].alleleStateVector[l];
				matPat.paternal = allelel-1;
				// both alleles at a marker assumed present or absent
				if(allelePat!=0 && allelePat!=alleleMat && matPat.paternal==matPat.maternal){
					continue;
				}
				ind->genotNodeVector[j].genotypeVector.push_back(matPat);
			}
		}
	}
}
/*! \fn void Population::finishUpAlleleStateInit(void)
*  \brief inputs all possible maternal and paternal allele states
*/

void Population::initAlleleNodeList(unsigned whichLocus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors: 
	GNodeSet::currentLocus=whichLocus;
	renum();
	Individual *ind, *mom, *dad;
	gNodeList.resize(2*popsize);
	gNodeList.howToSample = "single";
	GNode::gNodeListPtr = &gNodeList;
//	vectorOfGNsts.clear();
	gNodeList.releaseGNsts();
	unsigned i,j;
	for(i=0;i<popsize;i++){
		j=i*2;
		ind=popmember[i];
		mom=ind->mymother;
		dad=ind->myfather;
		AlleleStateNode* imAllele = &ind->malleleStateNodeVector[ind->currentLocus];
		AlleleStateNode* ipAllele = &ind->palleleStateNodeVector[ind->currentLocus];
		ipAllele->id = j+1;
		ipAllele->connectFlag = j+1;
		ipAllele->numberOfCuts= 0;
		imAllele->id = j+2;
		imAllele->connectFlag = j+2;
		imAllele->numberOfCuts= 0;
		gNodeList[j]   = ipAllele;
		gNodeList[j+1] = imAllele;
		
		if(prior->chrom()[0].locus[Individual::currentLocus].qtl_ml=='r'){
			DisAllelePenetranceSet *aPenSet = new DisAllelePenetranceSet;
			aPenSet->owner = ind;
			aPenSet->insert(imAllele);
			aPenSet->insert(ipAllele);
//			vectorOfGNsts.push_back(aPenSet);
			gNodeList.completeSetofGNsts.insert(aPenSet);
			aPenSet->attachMeToMyGnodes();
		}
		else{
			AllelePenetranceSet *aPenSet = new AllelePenetranceSet;
			aPenSet->owner = ind;
			aPenSet->insert(imAllele);
			aPenSet->insert(ipAllele);
//			vectorOfGNsts.push_back(aPenSet);
			gNodeList.completeSetofGNsts.insert(aPenSet);
			aPenSet->attachMeToMyGnodes();
		}
		
		if(mom){
			AlleleStateNode *mMAllele = &mom->malleleStateNodeVector[ind->currentLocus];
			AlleleStateNode *mPAllele = &mom->palleleStateNodeVector[ind->currentLocus];
			AlleleStateNode *pMAllele = &dad->malleleStateNodeVector[ind->currentLocus];
			AlleleStateNode *pPAllele = &dad->palleleStateNodeVector[ind->currentLocus];
			
			TransmissionSet *mTransmSet = new TransmissionSet;
			mTransmSet->offspring = ind;
			mTransmSet->paternal = false;
			mTransmSet->insert(imAllele);
			mTransmSet->insert(mMAllele);
			mTransmSet->insert(mPAllele);
//			vectorOfGNsts.push_back(mTransmSet);
			gNodeList.completeSetofGNsts.insert(mTransmSet);
			mTransmSet->attachMeToMyGnodes();
			
			TransmissionSet *pTransmSet = new TransmissionSet;
			pTransmSet->offspring = ind;
			pTransmSet->paternal = true;
			pTransmSet->insert(ipAllele);
			pTransmSet->insert(pMAllele);
			pTransmSet->insert(pPAllele);
//			vectorOfGNsts.push_back(pTransmSet);
			gNodeList.completeSetofGNsts.insert(pTransmSet);
			pTransmSet->attachMeToMyGnodes();
		}
		else {
			AlleleFounderSet *mFoundSet = new AlleleFounderSet;
			mFoundSet->insert(imAllele);
//			vectorOfGNsts.push_back(mFoundSet);
			gNodeList.completeSetofGNsts.insert(mFoundSet);
			mFoundSet->attachMeToMyGnodes();
			
			AlleleFounderSet *pFoundSet = new AlleleFounderSet;
			pFoundSet->insert(ipAllele);
//			vectorOfGNsts.push_back(pFoundSet);
			gNodeList.completeSetofGNsts.insert(pFoundSet);
			pFoundSet->attachMeToMyGnodes();
		}
	}
}
/*! \fn void Population::initAlleleNodeList(unsigned whichLocus)
*  \brief creates the allele node list, and the founder, penetrance,
and transmission sets for each allele node
*/

void Population::initGenotypeNodeList(unsigned whichLocus){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (June, 2003) 
	// Contributors: 
	GNodeSet::currentLocus = whichLocus;
	Individual::currentLocus = whichLocus;
	renum();
	Individual *ind, *mom, *dad;
	gNodeList.resize(popsize);
	gNodeList.howToSample = "single";
	GNode::gNodeListPtr = &gNodeList;
//	vectorOfGNsts.clear();
	gNodeList.releaseGNsts();
	unsigned i;
	unsigned rec = 0,rec1=0;
	cout << "\n Building GNodeList \n";
	for(i=0;i<popsize;i++){
		rec++;
		if(rec==1000){
			cout<<rec1+rec <<"\r";			
			cout.flush();
			rec1+=rec;
			rec = 0;
		}
		ind=popmember[i];
		mom=ind->mymother;
		dad=ind->myfather;
		GenotypeNode* indGenotype = &ind->genotNodeVector[ind->currentLocus];
		
		indGenotype->id = i+1;
		indGenotype->connectFlag = i+1;
		indGenotype->numberOfCuts= 0;
		gNodeList[i]   = indGenotype;
		
		if(prior->chrom()[0].locus[Individual::currentLocus].qtl_ml=='q'){
			GenoPenetranceSet *gPenSet = new GenoPenetranceSet;
			gPenSet->owner = ind;
			gPenSet->insert(indGenotype);
//			vectorOfGNsts.push_back(gPenSet);
			gNodeList.completeSetofGNsts.insert(gPenSet);
			gPenSet->attachMeToMyGnodes();
		}
		else if(prior->chrom()[0].locus[Individual::currentLocus].qtl_ml=='r'){
			DisGenoPenetranceSet *gPenSet = new DisGenoPenetranceSet;
			gPenSet->owner = ind;
			gPenSet->insert(indGenotype);
//			vectorOfGNsts.push_back(gPenSet);
			gNodeList.completeSetofGNsts.insert(gPenSet);
			gPenSet->attachMeToMyGnodes();
		} 
		if(mom){
			GenotypeNode *momGenotype = &mom->genotNodeVector[ind->currentLocus];
			GenotypeNode *dadGenotype = &dad->genotNodeVector[ind->currentLocus];
			TransitionSet *TransitSet = new TransitionSet;
			TransitSet->offspring = ind;
			TransitSet->insert(indGenotype);
			TransitSet->insert(momGenotype);
			TransitSet->insert(dadGenotype);
//			vectorOfGNsts.push_back(TransitSet);
			gNodeList.completeSetofGNsts.insert(TransitSet);
			TransitSet->attachMeToMyGnodes();
		}
		else {
			GenoFounderSet *FoundSet = new GenoFounderSet;
			FoundSet->insert(indGenotype);
//			vectorOfGNsts.push_back(FoundSet);
			gNodeList.completeSetofGNsts.insert(FoundSet);
			FoundSet->attachMeToMyGnodes();
		}
	}
}
/*! \fn void Population::initGenotypeNodeList(unsigned whichLocus)
*  \brief creates the genotype node list, and the founder, penetrance,
and transmission sets for each genotype node
*/

void Population::getInitialGNodeListSample(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock, string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (February, 2004) 
	// Contributors: 
	unsigned printFlag = model->myRSamplerParms.printFlag; 
	unsigned startLocus= model->myRSamplerParms.startLocus;
	// get a sample for the startLocus
	lookToYourLeft  = false;
	lookToYourRight = false;
	Individual::currentPosition = startLocus;
	unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
	Individual::currentLocus = j;
	if (printFlag>0) cout << "Sample locus: " << j << endl;
	getInitialSampleForLocus(j,startBlock,stopBlock,sizeBlock);
	// get samples for the other loci 
	// first to the right
	for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
		std::cout<<"+";
		std::cout.flush();
		lookToYourLeft  = true;
		lookToYourRight = false;
		Individual::currentPosition = jj;
		unsigned j = prior->chrom()[0].locusMask[jj].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j << endl;
		getInitialSampleForLocus(j,startBlock,stopBlock,sizeBlock);
	}
	// next to the left
	for (int jj=startLocus-1;jj>=0;jj--){
		std::cout<<"-";std::cout.flush();
		lookToYourLeft  = false;
		lookToYourRight = true;
		Individual::currentPosition = jj;
		unsigned j = prior->chrom()[0].locusMask[jj].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j << endl;
		getInitialSampleForLocus(j,startBlock,stopBlock,sizeBlock);
	}
	std::cout<<endl;

}
/*! \fn void Population::getInitialGNodeListSample(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock, string samplerType[])
*  \brief get the initial sample and also determine and store the peeling and sampling order for each locus starting with an arbitrary locus
*/

void Population::getInitialSampleForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (December, 2004) 
	// Contributors: 
	std::string samplerType = model->myRSamplerParms.samplerType;	
	unsigned maxsize = model->myRSamplerParms.maxCutsetSize;
	if (samplerType=="genotypic"){
		prior->chrom()[0].locus[j].peelOrder.resize(popsize);
		initGenotypeNodeList(j);
	} 
	else if (samplerType=="allelic"){
		prior->chrom()[0].locus[j].peelOrder.resize(2*popsize);
		initAlleleNodeList(j);
		//gNodeList.displayGNodeSets();
	}
	else {
		cerr << "Sampler type should be either genotypic or allelic;" << endl;
		exit(1);
	}
	gNodeList.peelOrderCutAndSample(maxsize,startBlock,stopBlock,sizeBlock); 
	setSegregationIndex(j,samplerType);
}
/*! \fn void Population::getInitialSampleForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock)
*  \brief get the initial sample and also determine and store the peeling and sampling order for a given locus
*/

void Population::getGNodeListSample(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock, string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// Contributors: 
	unsigned printFlag = model->myRSamplerParms.printFlag; 
	unsigned startLocus= model->myRSamplerParms.startLocus;
	// get a sample for the startLocus
	lookToYourLeft  = false;
	lookToYourRight = false;
	Individual::currentPosition = startLocus;
	unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
	Individual::currentLocus = j;
	if (printFlag>0) cout << "Sample locus: " << j << endl;
	getSampleForLocus(j,startBlock,stopBlock,sizeBlock);
	// get samples for the other loci 
	// first to the right
	for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
		lookToYourLeft  = true;
		lookToYourRight = false;
		Individual::currentPosition = jj;
		unsigned j = prior->chrom()[0].locusMask[jj].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j << endl;
		getSampleForLocus(j,startBlock,stopBlock,sizeBlock);
	}
	// next to the left
	for (int jj=startLocus-1;jj>=0;jj--){
		lookToYourLeft  = false;
		lookToYourRight = true;
		Individual::currentPosition = jj;
		unsigned j = prior->chrom()[0].locusMask[jj].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j << endl;
		getSampleForLocus(j,startBlock,stopBlock,sizeBlock);
	}
}
/*! \fn void Population::getGNodeListSample(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock, string samplerType)
*  \brief obtain a new candidate sample using peeling and cutting (if necessary)  
*/

void Population::getSampleForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (December, 2004) 
	// Contributors: 
	std::string samplerType = model->myRSamplerParms.samplerType;	
	unsigned maxsize = model->myRSamplerParms.maxCutsetSize;
	if(samplerType=="genotypic"){
		initGenotypeNodeList(j);
	} 
	else if(samplerType=="allelic"){
		initAlleleNodeList(j);
	}
	else {
		cerr << "Sampler type should be either genotypic or allelic;" << endl;
		exit(1);
	}
	gNodeList.peelCutAndSample(maxsize,startBlock,stopBlock,sizeBlock); 
	setSegregationIndex(j,samplerType);
}
/*! \fn void Population::getSampleForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock)
*  \brief get a sample for a given locus
*/

void Population::getOldGNodeListProbability(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock, string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// Contributors: 
	unsigned printFlag = model->myRSamplerParms.printFlag; 
	unsigned startLocus= model->myRSamplerParms.startLocus;
	// get a sample for the startLocus
	lookToYourLeft  = false;
	lookToYourRight = false;
	Individual::currentPosition = startLocus;
	unsigned j = prior->chrom()[0].locusMask[startLocus].index; 
	Individual::currentLocus = j;
	if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
	getOldProbabilityForLocus(j,startBlock,stopBlock,sizeBlock);
	// get samples for the other loci 
	// first to the right
	for (unsigned jj=startLocus+1;jj<prior->chrom()[0].nloci();jj++){
		lookToYourLeft  = true;
		lookToYourRight = false;
		Individual::currentPosition = jj;
		unsigned j = prior->chrom()[0].locusMask[jj].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
		getOldProbabilityForLocus(j,startBlock,stopBlock,sizeBlock);
	}
	// next to the left
	for (int jj=startLocus-1;jj>=0;jj--){
		lookToYourLeft  = false;
		lookToYourRight = true;
		Individual::currentPosition = jj;
		unsigned j = prior->chrom()[0].locusMask[jj].index; 
		Individual::currentLocus = j;
		if (printFlag>0) cout << "Sample locus: " << j+1 << endl;
		getOldProbabilityForLocus(j,startBlock,stopBlock,sizeBlock);
	}
}	
/*! \fn void Population::getOldGNodeListProbability(unsigned maxsize, unsigned startBlock, unsigned stopBlock, unsigned sizeBlock, string samplerType)
*  \brief recompute the needed quantities for the existing (old)
sample (for the MH step)
*/ 

void Population::getOldProbabilityForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (December, 2004) 
	// Contributors: 
	std::string samplerType = model->myRSamplerParms.samplerType;	
	unsigned maxsize = model->myRSamplerParms.maxCutsetSize;
	if (samplerType=="genotypic"){
		initGenotypeNodeList(j);
	} 
	else if (samplerType=="allelic"){
		initAlleleNodeList(j);
	}
	else {
		cerr << "Sampler type should be either genotypic or allelic;" << endl;
		exit(1);
	}
	gNodeList.peelCutAndCompute(maxsize,startBlock,stopBlock,sizeBlock); 
	setSegregationIndex(j,samplerType);
}
/*! \fn void Population::getOldProbabilityForLocus(unsigned j,unsigned startBlock, unsigned stopBlock, unsigned sizeBlock)
*  \brief calculate the probability of an old genotype configuration given the new probability structures used to obtain 
*         the new sample.
*/

void Population::copyCandidateToAccepted(string samplerType){
	// Authors: L. Radu Totir
	// (October, 2003) 
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for(unsigned j=0;j<Individual::numLoci;j++){
			if(samplerType=="genotypic"){
				ind->genotNodeVector[j].acceptedGenotypeState = ind->genotNodeVector[j].candidateGenotypeState;  
				ind->genotNodeVector[j].acceptedGenotypeVector = ind->genotNodeVector[j].genotypeVector; 
			}
			else if(samplerType=="allelic"){
				ind->malleleStateNodeVector[j].acceptedAlleleState       = ind->malleleStateNodeVector[j].candidateAlleleState;
				ind->malleleStateNodeVector[j].acceptedAlleleStateVector = ind->malleleStateNodeVector[j].alleleStateVector; 
				ind->palleleStateNodeVector[j].acceptedAlleleState       = ind->palleleStateNodeVector[j].candidateAlleleState;
				ind->palleleStateNodeVector[j].acceptedAlleleStateVector = ind->palleleStateNodeVector[j].alleleStateVector; 
			}
			else {
				cerr << "Sampler type should be either genotypic or allelic;" << endl;
				exit(1);
			}
		}
		ind->m_gameteAccepted = ind->m_gamete1;
		ind->p_gameteAccepted = ind->p_gamete1;
	}
}
/*! \fn Population::copyCandidateToAccepted(string samplerType)
*  \brief save candidate state as accepted states; allele origins for accepted states are already in (p)m_gameteOld 
*/

void Population::copyGNodeStatesToCandidateStates(string samplerType){
	// Authors: L. Radu Totir
	// (October, 2003) 
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for(unsigned j=0;j<Individual::numLoci;j++){
			if(samplerType=="genotypic"){
				ind->genotNodeVector[j].candidateGenotypeState = ind->genotNodeVector[j].genotypeState;  
			}
			else if(samplerType=="allelic"){
				ind->malleleStateNodeVector[j].candidateAlleleState = ind->malleleStateNodeVector[j].alleleState;
				ind->palleleStateNodeVector[j].candidateAlleleState = ind->palleleStateNodeVector[j].alleleState;
			}
			else {
				cerr << "Sampler type should be either genotypic or allelic;" << endl;
				exit(1);
			}
		}
		ind->m_gamete1 = ind->m_gamete;
		ind->p_gamete1 = ind->p_gamete;
	}
}
/*! \fn Population::copyGNodeStatesToCandidateStates(string samplerType)
*  \brief saves the GNode states (alleles or genotypes) for candidate sample
*/

void Population::copyCandidateGamete(void){
	// Authors: L. Radu Totir
	// (October, 2003) 
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		ind->m_gamete = ind->m_gameteAccepted;
		ind->p_gamete = ind->p_gameteAccepted;
	}
}
/*! \fn Population::copyCandidateStateGamete(void)
*   \brief copies candidate gamete into gamete
*/


void Population::storeSampledGametes(void){
	// Authors: L. Radu Totir
	// (Summer, 2004) 
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for(unsigned j=0;j<Individual::numLoci;j++){ 
			ind->m_gameteOld[j] = ind->m_gamete[j];
			ind->p_gameteOld[j] = ind->p_gamete[j];
		}
	}
}
/*! \fn Population::storeSampledGametes(void)
*  \brief store sampled gametes to be able to retreive them if we  
/*        switch from Gibbs to MH in the gibbsMH sampler 
*/

void Population::retreiveSampledGametes(void){
	// Authors: L. Radu Totir
	// (Summer, 2004) 
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for(unsigned j=0;j<Individual::numLoci;j++){ 
			// I need next two lines for the gibbsMH sampler
			ind->m_gamete[j] = ind->m_gameteOld[j];
			ind->p_gamete[j] = ind->p_gameteOld[j];
		}
	}
}
/*! \fn Population::retreiveSampledGametes(void)
*  \brief retreive the gametes sampled by Gibbs if MH rejects the sample
*/

void Population::setSegregationIndex(unsigned atLocus,string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2003) 
	// Contributors: 
	Individual *ind, *mom, *dad;
	unsigned im,ip,mm,mp,fm,fp;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		dad=ind->myfather;
		if(mom){
			if (samplerType=="genotypic"){
				im = ind->genotNodeVector[atLocus].getmState();
				ip = ind->genotNodeVector[atLocus].getpState();
				mm = mom->genotNodeVector[atLocus].getmState();
				mp = mom->genotNodeVector[atLocus].getpState();
			}
			else if (samplerType=="allelic"){
				im = ind->malleleStateNodeVector[atLocus].getMyAlleleState();
				ip = ind->palleleStateNodeVector[atLocus].getMyAlleleState();
				mm = mom->malleleStateNodeVector[atLocus].getMyAlleleState();
				mp = mom->palleleStateNodeVector[atLocus].getMyAlleleState();	
			}
			if (mm==mp){
				ind->m_gamete[atLocus] =9999;
			}
			else if (im==mm){
				ind->m_gamete[atLocus] = 0;
			}
			else if (im==mp){
				ind->m_gamete[atLocus] = 1;
			}
			else {
				cerr << "Error in: Population::setSegregationIndex(unsigned atLocus,string samplerType)" << endl;
				exit(1);
			}
			if (samplerType=="genotypic"){
				fm = dad->genotNodeVector[atLocus].getmState();
				fp = dad->genotNodeVector[atLocus].getpState(); 
			}
			else if (samplerType=="allelic"){
				fm = dad->malleleStateNodeVector[atLocus].getMyAlleleState();
				fp = dad->palleleStateNodeVector[atLocus].getMyAlleleState(); 
			}
			if (fm==fp){
				ind->p_gamete[atLocus] =9999;
			}
			else if (ip==fm){
				ind->p_gamete[atLocus] = 0;
			}
			else if (ip==fp){
				ind->p_gamete[atLocus] = 1;
			}
			else {
				cerr << "Error in: Population::setSegregationIndex(unsigned atLocus,string samplerType)" << endl;
				exit(1);
			}
		}
	}
}
/*! \fn Population::setSegregationIndex(unsigned atLocus,string samplerType)
*  \brief set the segregation indicators when these can be infered from 
*         the ordered genotypes that have been sampled   
*/

void Population::displaySegregationIndicators(std::ostream &outfile) {
	// Authors: Rohan L. Fernando
	// (November, 2004)
	// Contributors: 
	for(unsigned i=0;i<popsize;i++){
		Individual *ind=popmember[i];
		outfile << setw(5) << ind_name(i+1) << "  ";
		Individual *mom=ind->mymother;
		if (mom){
			for (unsigned locus = 0; locus < Individual::numLoci; locus++){
				outfile << ind->m_gamete[locus]	<< " " << ind->p_gamete[locus] << "  ";
			}
		}
		else {
			for (unsigned locus = 0; locus < Individual::numLoci; locus++){
				outfile << " 9  9  ";
			}
		}
		outfile << endl;
	}
}

void Population::PrintSegregationIndicators(void){
	// Authors: Helene Gilbert, L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// Contributors: 
	Individual *ind,*mom;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(mom){
			cout << "individual....."<< popmember[i]->myid <<endl;
			cout << "m_gamete....."<<endl;
			cout << ind->m_gamete;  
			cout << "p_gamete....."<<endl;
			cout << ind->p_gamete;
		}
	}
}
/*! \fn Population::PrintSegregationIndicators(void)
*  \brief print the segregation indicators of the non-founders   
*/

void Population::sampleSegregationIndicators(void){
	// Authors: Helene Gilbert, Honghua Zhao, L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// modifications November 2004, RLF  
	// Contributors: 
	Individual *ind, *mom;
	lookToYourLeft  = true;
	lookToYourRight = true;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(mom){
			for (unsigned atMarker = 0; atMarker < Individual::numLoci; atMarker++){
				Individual::currentPosition = atMarker;
				if (ind->m_gameteAccepted[atMarker] == 9999){
					double u = ranf();
					int RightLocus = ind->findRightLocus(ind->m_gameteAccepted);
					int LeftLocus  = ind->findLeftLocus(ind->m_gameteAccepted);
					if ( (RightLocus != -1) && (LeftLocus != -1) ) {
						double recLeft  = recombinationMatrix[LeftLocus][atMarker];
						double recRight = recombinationMatrix[atMarker][RightLocus];
						double recFlank = recombinationMatrix[LeftLocus][RightLocus];
						if (ind->m_gameteAccepted[LeftLocus] == ind->m_gameteAccepted[RightLocus]){
							ind->m_gameteAccepted[atMarker] = 
							u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->m_gameteAccepted[LeftLocus] : 1 - ind->m_gameteAccepted[LeftLocus];
						}
						else{
							ind->m_gameteAccepted[atMarker] = 
							u <((  recLeft)*(1-recRight))/(  recFlank) ? ind->m_gameteAccepted[RightLocus] : ind->m_gameteAccepted[LeftLocus];
						}
						
					}
					else if (RightLocus != -1) {	 
						double recRight = recombinationMatrix[atMarker][RightLocus];
						ind->m_gameteAccepted[atMarker] = u < (1-recRight) ? ind->m_gameteAccepted[RightLocus] : 1-ind->m_gameteAccepted[RightLocus] ;
					}
					else if (LeftLocus != -1) {	
						double recLeft = recombinationMatrix[LeftLocus][atMarker];	
						ind->m_gameteAccepted[atMarker] = u < (1-recLeft) ? ind->m_gameteAccepted[LeftLocus] : 1-ind->m_gameteAccepted[LeftLocus];
					}
					else {
						ind->m_gameteAccepted[atMarker] = (u < 0.5) ? 1 : 0 ;
					}
				}
				
				if (ind->p_gameteAccepted[atMarker] == 9999){
					double u = ranf();
					int RightLocus = ind->findRightLocus(ind->p_gameteAccepted);
					int LeftLocus  = ind->findLeftLocus(ind->p_gameteAccepted);
					if ( (RightLocus != -1) && (LeftLocus != -1) ) {
						double recLeft  = recombinationMatrix[LeftLocus][atMarker];
						double recRight = recombinationMatrix[atMarker][RightLocus];
						double recFlank = recombinationMatrix[LeftLocus][RightLocus];
						if (ind->p_gameteAccepted[LeftLocus] == ind->p_gameteAccepted[RightLocus]){
							ind->p_gameteAccepted[atMarker] = 
							u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->p_gameteAccepted[LeftLocus] : 1 - ind->p_gameteAccepted[LeftLocus];
						}
						else{
							ind->p_gameteAccepted[atMarker] = 
							u <((  recLeft)*(1-recRight))/(  recFlank) ? ind->p_gameteAccepted[RightLocus] : ind->p_gameteAccepted[LeftLocus];
						}
						
					}
					else if (RightLocus != -1) {	 
						double recRight = recombinationMatrix[atMarker][RightLocus];
						ind->p_gameteAccepted[atMarker] = u < (1-recRight) ? ind->p_gameteAccepted[RightLocus] : 1-ind->p_gameteAccepted[RightLocus] ;
					}
					else if (LeftLocus != -1) {	
						double recLeft = recombinationMatrix[LeftLocus][atMarker];	
						ind->p_gameteAccepted[atMarker] = u < (1-recLeft) ? ind->p_gameteAccepted[LeftLocus] : 1-ind->p_gameteAccepted[LeftLocus];
					}
					else {
						ind->p_gameteAccepted[atMarker] = (u < 0.5) ? 1 : 0 ;
					}
				}
			}
		}
	}
}

void Population::sampleSegregationIndicatorsSimple(void){
	// Authors: Helene Gilbert, Honghua Zhao, L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// modifications November 2004, RLF  
	// Contributors: 
	Individual *ind, *mom;
	lookToYourLeft  = true;
	lookToYourRight = true;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(mom){
			for (unsigned atMarker = 0; atMarker < Individual::numLoci; atMarker++){
				Individual::currentPosition = atMarker;
				if (ind->m_gamete[atMarker] == 9999){
					double u = ranf();
					int RightLocus = ind->findRightLocus(ind->m_gamete);
					int LeftLocus  = ind->findLeftLocus(ind->m_gamete);
					if ( (RightLocus != -1) && (LeftLocus != -1) ) {
						double recLeft  = recombinationMatrix[LeftLocus][atMarker];
						double recRight = recombinationMatrix[atMarker][RightLocus];
						double recFlank = recombinationMatrix[LeftLocus][RightLocus];
						if (ind->m_gamete[LeftLocus] == ind->m_gamete[RightLocus]){
							ind->m_gamete[atMarker] = 
							u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->m_gamete[LeftLocus] : 1 - ind->m_gamete[LeftLocus];
						}
						else{
							ind->m_gamete[atMarker] = 
							u <((  recLeft)*(1-recRight))/(  recFlank) ? ind->m_gamete[RightLocus] : ind->m_gamete[LeftLocus];
						}
						
					}
					else if (RightLocus != -1) {	 
						double recRight = recombinationMatrix[atMarker][RightLocus];
						ind->m_gamete[atMarker] = u < (1-recRight) ? ind->m_gamete[RightLocus] : 1-ind->m_gamete[RightLocus] ;
					}
					else if (LeftLocus != -1) {	
						double recLeft = recombinationMatrix[LeftLocus][atMarker];	
						ind->m_gamete[atMarker] = u < (1-recLeft) ? ind->m_gamete[LeftLocus] : 1-ind->m_gamete[LeftLocus];
					}
					else {
						ind->m_gamete[atMarker] = (u < 0.5) ? 1 : 0 ;
					}
				}
				
				if (ind->p_gamete[atMarker] == 9999){
					double u = ranf();
					int RightLocus = ind->findRightLocus(ind->p_gamete);
					int LeftLocus  = ind->findLeftLocus(ind->p_gamete);
					if ( (RightLocus != -1) && (LeftLocus != -1) ) {
						double recLeft  = recombinationMatrix[LeftLocus][atMarker];
						double recRight = recombinationMatrix[atMarker][RightLocus];
						double recFlank = recombinationMatrix[LeftLocus][RightLocus];
						if (ind->p_gamete[LeftLocus] == ind->p_gamete[RightLocus]){
							ind->p_gamete[atMarker] = 
							u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->p_gamete[LeftLocus] : 1 - ind->p_gamete[LeftLocus];
						}
						else{
							ind->p_gamete[atMarker] = 
							u <((  recLeft)*(1-recRight))/(  recFlank) ? ind->p_gamete[RightLocus] : ind->p_gamete[LeftLocus];
						}
						
					}
					else if (RightLocus != -1) {	 
						double recRight = recombinationMatrix[atMarker][RightLocus];
						ind->p_gamete[atMarker] = u < (1-recRight) ? ind->p_gamete[RightLocus] : 1-ind->p_gamete[RightLocus] ;
					}
					else if (LeftLocus != -1) {	
						double recLeft = recombinationMatrix[LeftLocus][atMarker];	
						ind->p_gamete[atMarker] = u < (1-recLeft) ? ind->p_gamete[LeftLocus] : 1-ind->p_gamete[LeftLocus];
					}
					else {
						ind->p_gamete[atMarker] = (u < 0.5) ? 1 : 0 ;
					}
				}
			}
		}
	}
}



//			void Population::sampleSegregationIndicators(void){
//	// Authors: Helene Gilbert, Honghua Zhao, L. Radu Totir and Rohan L. Fernando 
//	// (October, 2003) 
//	// Contributors: 
//	Individual *ind, *mom;
//	for(unsigned i=0;i<popsize;i++){
//		ind=popmember[i];
//		mom=ind->mymother;
//		if(mom){
//			// set segregation indicators for the first locus
//			unsigned atMarker = 0;
//			Individual::currentLocus = atMarker;
//			if (ind->m_gamete[atMarker] == 9999){
//				double u = ranf();  
//				int RightLocus = ind->findRightLocus((ind->m_gamete));
//				double recRight = recombinationMatrix[atMarker][RightLocus];
//				if(RightLocus != -1){
//					ind->m_gamete[atMarker] = u < (1-recRight) ? ind->m_gamete[RightLocus] : 1-ind->m_gamete[RightLocus];
//				}
//				else{
//					ind->m_gamete[atMarker] = u < 0.5 ? 1 : 0 ;
//				}
//			}
//			if (ind->p_gamete[atMarker] == 9999){
//				double u = ranf();
//				int RightLocus = ind->findRightLocus(ind->p_gamete);
//				double recRight = recombinationMatrix[atMarker][RightLocus];
//				if(RightLocus != -1){
//					ind->p_gamete[atMarker] = u < (1-recRight) ? ind->p_gamete[RightLocus] : 1-ind->p_gamete[RightLocus] ;
//				}
//				else{
//					ind->p_gamete[atMarker] = u < 0.5 ? 1 : 0 ;
//				}
//			}
//		}
//	}
//			// set segregation indicators for the loci between the first and the last
//			for (unsigned atMarker = 1; atMarker < Individual::numLoci-Individual::nQTL - 1; atMarker++){
//				Individual::currentLocus = atMarker ;
//				if (ind->m_gamete[atMarker] == 9999){
//					double u = ranf();
//					int RightLocus = ind->findRightLocus(ind->m_gamete);
//					int LeftLocus  = ind->findLeftLocus(ind->m_gamete);
//					double recLeft  = recombinationMatrix[LeftLocus][atMarker];
//					double recRight = recombinationMatrix[atMarker][RightLocus];
//					double recFlank = recombinationMatrix[LeftLocus][RightLocus];
//					if (ind->m_gamete[LeftLocus] == ind->m_gamete[RightLocus]){
//						ind->m_gamete[atMarker] = 
//						u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->m_gamete[LeftLocus] : 1 - ind->m_gamete[LeftLocus];
//					}
//					else{
//						ind->m_gamete[atMarker] = 
//						u <((  recLeft)*(1-recRight))/(1-recFlank) ? ind->m_gamete[RightLocus] : ind->m_gamete[LeftLocus];
//					}
//				}
//				if (ind->p_gamete[atMarker] == 9999){ 
//					double u = ranf();
//					int RightLocus = ind->findRightLocus(ind->p_gamete);
//					int LeftLocus  = ind->findLeftLocus(ind->p_gamete);  
//					double recLeft  = recombinationMatrix[LeftLocus][atMarker];
//					double recRight = recombinationMatrix[atMarker][RightLocus];
//					double recFlank = recombinationMatrix[LeftLocus][RightLocus];
//					if (ind->p_gamete[LeftLocus] == ind->p_gamete[RightLocus]){
//						ind->p_gamete[atMarker] = 
//						u <((1-recLeft)*(1-recRight))/(1-recFlank) ? ind->p_gamete[LeftLocus] : 1 - ind->p_gamete[LeftLocus];
//					}
//					else{
//						ind->p_gamete[atMarker] = 
//						u <((recLeft)*(1-recRight))/(1-recFlank) ? ind->p_gamete[RightLocus] : ind->p_gamete[LeftLocus];
//					}
//				}      
//			}
//			// set segregation indicators for the last locus
//			atMarker = Individual::numLoci - Individual::nQTL-1;
//			Individual::currentLocus = atMarker;
//			if (ind->m_gamete[atMarker] == 9999){
//				int  LeftLocus = ind->findLeftLocus(ind->m_gamete);
//				double recLeft = recombinationMatrix[LeftLocus][atMarker];
//				double u = ranf();
//				ind->m_gamete[atMarker] = u < (1-recLeft) ? ind->m_gamete[LeftLocus] : 1-ind->m_gamete[LeftLocus];
//			}
//			if (ind->p_gamete[atMarker] == 9999){
//				int  LeftLocus = ind->findLeftLocus(ind->p_gamete);
//				double recLeft = recombinationMatrix[LeftLocus][atMarker];
//				double u = ranf();
//				ind->p_gamete[atMarker] = u < (1-recLeft) ? ind->p_gamete[LeftLocus] : 1-ind->p_gamete[LeftLocus] ;
//			}
//		}
//	}
//}
/*! \fn Population::sampleSegregationIndicators(void)
*  \brief sample the unknown segregation indicators of non-founders based 
*  on the ordered parental genotypes
*/

// Actually not used !!!
void Population::resizeSegregationIndicators(){
	// Authors: Helene Gilbert, Honghua Zhao, L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// Resize the segregation indicators : 1st position for the QTL to use it in the descent graph
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		Vector <unsigned> m_g, p_g ;
		m_g.resize(Individual::numLoci,9999);    
		p_g.resize(Individual::numLoci,9999);
		m_g=ind->m_gamete;
		p_g=ind->p_gamete;
		ind->m_gamete.resize(Individual::numLoci+1,9999);   // 9999 for the QTL ????
		ind->p_gamete.resize(Individual::numLoci+1,9999);
		for (unsigned i = 0 ; i < Individual::numLoci ; i++){
			ind->m_gamete[i+1]=m_g[i];
			ind->p_gamete[i+1]=p_g[i];
		}
		cout<<"Final indicators .........."<<  endl;
		cout << ind->m_gamete<<  endl;
	}
}

void Population::ListAlleleFounders(void){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir
	Individual *ind, *mom;
	unsigned al ;
	unsigned nLoci = Individual::numLoci; 
	for (unsigned i=0; i<nLoci ; i++){            
		if(prior->chrom()[0].locus[i].qtl_ml != 'q'){
			prior->chrom()[0].locus[i].listOfAllele.resize(prior->chrom()[0].locus[i].nallele());
			prior->chrom()[0].locus[i].nextMarker = 0 ;
			unsigned il=i+1 ;
			while (!prior->chrom()[0].locus[i].nextMarker && il < nLoci){
				if(prior->chrom()[0].locus[i].qtl_ml != 'q') prior->chrom()[0].locus[i].nextMarker = il;
				il++;
			}
			unsigned nrecord = 0, j=0 ;
			while(nrecord+1 <= prior->chrom()[0].locus[i].nallele() && j < popsize){
				ind=popmember[j];
				mom=ind->mymother;
				if(!mom){	  
					for (unsigned ia = 0 ; ia < 2 ; ia++){
						unsigned recorded=0;  
						if(!ia){
							al = ind->genome0.chromosome[0].locus[i].allele ;	 
						}else{
							al = ind->genome1.chromosome[0].locus[i].allele ;	 
						}
						for (unsigned k = 0 ; k <= nrecord ; k++){
							unsigned alComp = prior->chrom()[0].locus[i].listOfAllele[k] ;
							if (al == alComp ) recorded=1;	
						}
						if(!recorded) {
							prior->chrom()[0].locus[i].listOfAllele[nrecord] = al ;
							nrecord++;
						}
					}
				}  
				j++ ;
			}       
		}
	}
}

void Population::SetPossibleHaplotypes(void){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir
	unsigned nLoci = Individual::numLoci; 
	for (unsigned i1=0; i1<nLoci-1 ; i1++){  // defined until (n-1)th locus
		if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
			unsigned i2 = prior->chrom()[0].locus[i1].nextMarker ;
			unsigned size1 = prior->chrom()[0].locus[i1].nallele() ;
			unsigned size2 = prior->chrom()[0].locus[i2].nallele() ;
			prior->chrom()[0].locus[i1].nHaplotypes = size1*size2 ;
			prior->chrom()[0].locus[i1].al1Haplo.resize(size1*size2); // should be # of markers
			prior->chrom()[0].locus[i1].al2Haplo.resize(size1*size2); 
			Vector<unsigned> haplo1 = prior->chrom()[0].locus[i1].al1Haplo ; 
			Vector<unsigned> haplo2 = prior->chrom()[0].locus[i1].al2Haplo ; 
			for (unsigned ia1 = 0 ; ia1 < size1 ; ia1++){
				for (unsigned ia2 = 0 ; ia2 < size2 ; ia2++){
					haplo1[ia2 + ia1*size2 ]=prior->chrom()[0].locus[i1].listOfAllele[ia1] ;
					haplo2[ia2 + ia1*size2 ]=prior->chrom()[0].locus[i2].listOfAllele[ia2] ;
				}
			}
			prior->chrom()[0].locus[i1].al1Haplo = haplo1 ; 
			prior->chrom()[0].locus[i1].al2Haplo = haplo2; 
		}
	}
}

void Population::initGenotypeFreq(void){
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for (unsigned j=0;j<Individual::numLoci;j++){
			unsigned numAlleles = prior->chromosome[0].locus[j].nallele();
			ind->genotNodeVector[j].genotypeCount.resize(numAlleles*numAlleles,0);
			ind->genotNodeVector[j].sampleCount = 0;
		}
	}
}

void Population::SetFreqHaploFounders(void){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir  
	Individual *ind, *mom;
	unsigned nLoci = Individual::numLoci;
	unsigned nMarkerLoci = nLoci - Individual::nQTL;
	unsigned hm, hp ;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(!mom){  
			// cout << "Individual " << ind->myid << endl;
			unsigned ihaplo=0;
			// # of intervals = nMarkerLoci-1
			ind->mHaplotypeFreq.resize(nMarkerLoci-1); 
			ind->pHaplotypeFreq.resize(nMarkerLoci-1);
			for (unsigned i1=0; i1<nLoci-1; i1++){  
				if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
					// cout << "Interval " << i1;
					unsigned i2 = prior->chrom()[0].locus[i1].nextMarker ;
					unsigned im1 = ind->genotNodeVector[i1].getmState()+1;    
					unsigned ip1 = ind->genotNodeVector[i1].getpState()+1;
					unsigned im2 = ind->genotNodeVector[i2].getmState()+1;    
					unsigned ip2 = ind->genotNodeVector[i2].getpState()+1;
					// cout << " Haplo mat " << im1 << "/" << im2 
					//      << " Haplo pat " << ip1 << "/" << ip2 << endl;
					ind->mHaplotypeFreq[ihaplo].resize(prior->chrom()[0].locus[i1].nHaplotypes);
					ind->pHaplotypeFreq[ihaplo].resize(prior->chrom()[0].locus[i1].nHaplotypes);
					for (unsigned ih = 0 ; ih < prior->chrom()[0].locus[i1].nHaplotypes ; ih++){
						ind->mHaplotypeFreq[ihaplo][ih] = 0 ;
						ind->pHaplotypeFreq[ihaplo][ih] = 0 ;
						if(im1 == prior->chrom()[0].locus[i1].al1Haplo[ih] && 
						   im2 == prior->chrom()[0].locus[i1].al2Haplo[ih] ) hm=ih ;
						if(ip1 == prior->chrom()[0].locus[i1].al1Haplo[ih] && 
						   ip2 == prior->chrom()[0].locus[i1].al2Haplo[ih] ) hp=ih ;
					}
					ind->mHaplotypeFreq[ihaplo][hm] ++;
					ind->pHaplotypeFreq[ihaplo][hp] ++; 
					ihaplo++;
				}
			}
			// cout << " / " << ind->mHaplotypeFreq[0] << endl; 
			// cout << " / " << ind->mHaplotypeFreq[1] << endl;
		}
	}
}

void Population::UpdateFreqHaploFounders(void){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir  
	Individual *ind, *mom;
	unsigned nLoci = Individual::numLoci; 
	unsigned hm , hp ;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(!mom){  
			//      cout<<"Individual "<<ind->myid<<endl;
			unsigned ihaplo=0;
			for (unsigned i1=0; i1<nLoci-1 ; i1++){ // defined until (n-1)th locus
				if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
					//	  cout<<"Interval "<<i1;
					//          cout<<ind->mHaplotypeFreq[i1];
					unsigned i2 = prior->chrom()[0].locus[i1].nextMarker ;
					unsigned im1 = ind->genotNodeVector[i1].getAcceptedMatState()+1;    
					unsigned ip1 = ind->genotNodeVector[i1].getAcceptedPatState()+1;
					unsigned im2 = ind->genotNodeVector[i2].getAcceptedMatState()+1;    
					unsigned ip2 = ind->genotNodeVector[i2].getAcceptedPatState()+1;
					//	  cout<<" Haplo mat "<<im1<<"/"<<im2<< " Haplo pat "<<ip1<<"/"<<ip2<<endl;
					for (unsigned ih = 0 ; ih < prior->chrom()[0].locus[i1].nHaplotypes ; ih++){
						if(im1 == prior->chrom()[0].locus[i1].al1Haplo[ih] && 
						   im2 == prior->chrom()[0].locus[i1].al2Haplo[ih] ) hm=ih ;
						if(ip1 == prior->chrom()[0].locus[i1].al1Haplo[ih] && 
						   ip2 == prior->chrom()[0].locus[i1].al2Haplo[ih] ) hp=ih ;
					}
					ind->mHaplotypeFreq[ihaplo][hm] ++;
					ind->pHaplotypeFreq[ihaplo][hp] ++; 
					ihaplo++;
				}
			}
			//      cout<<" / "<<ind->mHaplotypeFreq[0]<<endl<<" /  "<<endl<<ind->mHaplotypeFreq[1]<<endl;
		}
	}
}

void Population::UpdateFreqHaploFoundersSimple(void){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir  
	Individual *ind, *mom;
	unsigned nLoci = Individual::numLoci; 
	unsigned hm , hp ;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(!mom){  
			//      cout<<"Individual "<<ind->myid<<endl;
			unsigned ihaplo=0;
			for (unsigned i1=0; i1<nLoci-1 ; i1++){ // defined until (n-1)th locus
				if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
					//	  cout<<"Interval "<<i1;
					//          cout<<ind->mHaplotypeFreq[i1];
					unsigned i2 = prior->chrom()[0].locus[i1].nextMarker ;
					unsigned im1 = ind->genotNodeVector[i1].getmState()+1;    
					unsigned ip1 = ind->genotNodeVector[i1].getpState()+1;
					unsigned im2 = ind->genotNodeVector[i2].getmState()+1;    
					unsigned ip2 = ind->genotNodeVector[i2].getpState()+1;
					//	  cout<<" Haplo mat "<<im1<<"/"<<im2<< " Haplo pat "<<ip1<<"/"<<ip2<<endl;
					for (unsigned ih = 0 ; ih < prior->chrom()[0].locus[i1].nHaplotypes ; ih++){
						if(im1 == prior->chrom()[0].locus[i1].al1Haplo[ih] && 
						   im2 == prior->chrom()[0].locus[i1].al2Haplo[ih] ) hm=ih ;
						if(ip1 == prior->chrom()[0].locus[i1].al1Haplo[ih] && 
						   ip2 == prior->chrom()[0].locus[i1].al2Haplo[ih] ) hp=ih ;
					}
					ind->mHaplotypeFreq[ihaplo][hm] ++;
					ind->pHaplotypeFreq[ihaplo][hp] ++; 
					ihaplo++;
				}
			}
			//      cout<<" / "<<ind->mHaplotypeFreq[0]<<endl<<" /  "<<endl<<ind->mHaplotypeFreq[1]<<endl;
		}
	}
}
void Population::CalcFreqHaploFounders(unsigned numOfSamples){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir  
	Individual *ind, *mom;
	double doubleSamples = numOfSamples ;
	//  cout<<" nb samples "<<doubleSamples<<endl;
	unsigned nLoci = Individual::numLoci;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;
		if(!mom){  
			//      cout<<"Individual "<<ind->myid<<endl;
			unsigned ihaplo=0;
			for (unsigned i1=0; i1<nLoci-1 ; i1++){ // defined until (n-1)th locus
				if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
					ind->mHaplotypeFreq[ihaplo] /= doubleSamples ;
					ind->pHaplotypeFreq[ihaplo] /= doubleSamples ;
					ihaplo++;
				}
			}
			//       cout<<ind->mHaplotypeFreq[1]<<"  "<<ind->pHaplotypeFreq[1]<<"  "<<ind->mHaplotypeFreq[2]<<"  "<<ind->pHaplotypeFreq[2]<<"  "<<endl;
		}
	}
}  

void Population::DisplayFreqHaploFounders(){
	// Authors: Helene Gilbert and Rohan L. Fernando 
	// (Spring 2004) 
	// Contributors: L. Radu Totir
	const char *s = model->myRSamplerParms.resultsFile.c_str(); 
	ofstream outfile;
	outfile.open(s);
	//  cout<<endl<<" Founders'haplotype Frequencies "<<endl;
	Individual *ind, *mom;
	unsigned nLoci = Individual::numLoci;
	for (unsigned i1=0; i1<nLoci-1 ; i1++){ // defined until (n-1)th locus
		if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
			outfile << "Number of possible haplotypes at locus " << i1+1 << " = " 
			<< prior->chrom()[0].locus[i1].nHaplotypes << " => ";
			for (unsigned ih = 0; ih < prior->chrom()[0].locus[i1].nHaplotypes; ih++){
				outfile << "(" << prior->chrom()[0].locus[i1].al1Haplo[ih] 
				<< "," << prior->chrom()[0].locus[i1].al2Haplo[ih] << ") ";
			}
			outfile << endl;
		}
	}
	outfile << endl;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mom=ind->mymother;  
		if(!mom){ 
			outfile << "Individual " << ind_name(i+1) << endl;
			unsigned ihaplo=0;
			for (unsigned i1=0; i1<nLoci-1 ; i1++){ // defined until (n-1)th locus
				outfile << "For interval " << i1+1 << " => " << endl;
				if(prior->chrom()[0].locus[i1].qtl_ml != 'q' ){
					for (unsigned ih = 0; ih < prior->chrom()[0].locus[i1].nHaplotypes; ih++){
						if(ind->mHaplotypeFreq[ihaplo][ih]||ind->pHaplotypeFreq[ihaplo][ih]){
							outfile << "(" << prior->chrom()[0].locus[i1].al1Haplo[ih] 
							<< "," << prior->chrom()[0].locus[i1].al2Haplo[ih] << ") => ";
							outfile << "maternal = " << ind->mHaplotypeFreq[ihaplo][ih] << "    "
								<< "paternal = " << ind->pHaplotypeFreq[ihaplo][ih] << endl;
						}
					}
					ihaplo++;	
				}
			}
		}
	}
}

void IndexVector::setupVectors(SafeSTLVector<unsigned> aVector){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004)
	// Contributors:
	maxElementVector=aVector;
	size = aVector.size();
	multiplicationCode.resize(size);
	multiplicationCode[size-1]=1;
	unsigned temp = 1;
	if (size > 1){// needed for ppc64 
		for (long i=size-2;i>=0;i--){
			temp *= maxElementVector[i+1];
			multiplicationCode[i]=temp;
		} 
	}
}
/*! \fn void IndexVector::setupVectors(std::vector<unsigned> aVector)
*  \brief Initialize the multiplicationCode vector used to convert an 
*         integer (e.g. individual genotype) into the vector decomposition 
*         (e.g vector of genotypes at each locus). The input vector (aVector) 
*         specifies the number of possible values (genotypes) at each position 
*         in the vector (e.g. the number of alleles at 3 loci [2 2 3]). 
*/ 

SafeSTLVector<unsigned> IndexVector::getVector(unsigned index){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (based on the gtindex(...) of Tianlin Wang)
	// (August, 2004) 
	// Contributors:
	SafeSTLVector<unsigned> resultVector;
	resultVector.name = "IndexVector::getVector()::resultVector";
	resultVector.resize(size);
	unsigned ir=index;
	for (long i=size-1; i>=0; i--) {
		unsigned maxElement = maxElementVector[i];
		resultVector[i] = unsigned(fmod(double(ir),double(maxElement)));
		ir /= maxElement;
	}
	return resultVector;
}
/*! \fn std::vector<unsigned> IndexVector::getVector(unsigned index)
*  \brief Returns the corresponding vector decomposition for index.  
*/
unsigned IndexVector::getIndex(std::vector<unsigned> stateVector){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors:
	unsigned index=0;
	for(unsigned i=0;i<size;i++){
		index+=stateVector[i]*multiplicationCode[i];
	}
	return index;
}
/*! \fn unsigned IndexVector::getIndex(std::vector<unsigned> stateVector)
*  \brief Returns the index for the vector decomposition given in stateVector
*/

void Population::countHaplotypes(string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors: Chris Stricker
	SafeSTLVector<unsigned> matHaplotype, patHaplotype;
	unsigned nLoci = Individual::numLoci;
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		matHaplotype.resize(nLoci);
		patHaplotype.resize(nLoci);
		for (unsigned j=0; j<nLoci; j++){  
			if (samplerType=="genotypic"){
				matHaplotype[j] = ind->genotNodeVector[j].getAcceptedMatState();
				patHaplotype[j] = ind->genotNodeVector[j].getAcceptedPatState();
			}
			else if (samplerType=="allelic"){
				matHaplotype[j] = ind->malleleStateNodeVector[j].getAcceptedAlleleState();
				patHaplotype[j] = ind->palleleStateNodeVector[j].getAcceptedAlleleState();
			}
			else if(samplerType=="joint"){
				if (prior->chrom()[0].locus[j].gnodeType=="allelic"){
					matHaplotype[j]  = ind->malleleStateNodeVector[j].getAcceptedAlleleState();
					patHaplotype[j]  = ind->palleleStateNodeVector[j].getAcceptedAlleleState();
				}
				else {
					matHaplotype[j] = ind->genotNodeVector[j].getAcceptedMatState();
					patHaplotype[j] = ind->genotNodeVector[j].getAcceptedPatState();
				}	
			}
		}
		unsigned matIndex = haplotypeCoder.getIndex(matHaplotype);
		unsigned patIndex = haplotypeCoder.getIndex(patHaplotype);
		if(ind->matHaplotypeCount.find(matIndex) != ind->matHaplotypeCount.end() ) {
			ind->matHaplotypeCount[matIndex] +=1;
		}
		else{
			ind->matHaplotypeCount[matIndex] = 1;
		}
		if(ind->patHaplotypeCount.find(patIndex) != ind->patHaplotypeCount.end() ) {
			ind->patHaplotypeCount[patIndex] +=1;
		}
		else{
			ind->patHaplotypeCount[patIndex] = 1;
		}
	}
}
/*! \fn void Population::countHaplotypes(samplerType)
*  \brief count the number haplotypes sampled for each individual 
*/


void Population::countHaplotypesSimple(string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors: Chris Stricker
	SafeSTLVector<unsigned> matHaplotype, patHaplotype;
	unsigned nLoci = Individual::numLoci;
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		matHaplotype.resize(nLoci);
		patHaplotype.resize(nLoci);
		for (unsigned j=0; j<nLoci; j++){  
			if (samplerType=="genotypic"){
				matHaplotype[j] = ind->genotNodeVector[j].getmState();
				patHaplotype[j] = ind->genotNodeVector[j].getpState();
			}
			else if (samplerType=="allelic"){
				matHaplotype[j] = ind->malleleStateNodeVector[j].getMyAlleleState();
				patHaplotype[j] = ind->palleleStateNodeVector[j].getMyAlleleState();
			}
			else if(samplerType=="joint"){
				if (prior->chrom()[0].locus[j].gnodeType=="allelic"){
					matHaplotype[j]  = ind->malleleStateNodeVector[j].getMyAlleleState();
					patHaplotype[j]  = ind->palleleStateNodeVector[j].getMyAlleleState();
				}
				else {
					matHaplotype[j] = ind->genotNodeVector[j].getmState();
					patHaplotype[j] = ind->genotNodeVector[j].getpState();
				}	
			}
		}
		unsigned matIndex = haplotypeCoder.getIndex(matHaplotype);
		unsigned patIndex = haplotypeCoder.getIndex(patHaplotype);
		if(ind->matHaplotypeCount.find(matIndex) != ind->matHaplotypeCount.end() ) {
			ind->matHaplotypeCount[matIndex] +=1;
		}
		else{
			ind->matHaplotypeCount[matIndex] = 1;
		}
		if(ind->patHaplotypeCount.find(patIndex) != ind->patHaplotypeCount.end() ) {
			ind->patHaplotypeCount[patIndex] +=1;
		}
		else{
			ind->patHaplotypeCount[patIndex] = 1;
		}
	}
}
/*! \fn void Population::countHaplotypes(samplerType)
*  \brief count the number haplotypes sampled for each individual 
*/


void Population::displayHaplotypeFrequencies(unsigned numSamples, std::ostream &os){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors: Chris Stricker
    
	unsigned nLoci = Individual::numLoci;
	std::map<unsigned,unsigned>::iterator m;
	SafeSTLVector<unsigned> outVector;
	for(unsigned i=0;i<popsize;i++){
		os << "\nIndividual "<< ind_name(i+1) <<std::endl<<std::endl;
		os << "maternal haplotypes" << std::endl;
		for(m = member(i)->matHaplotypeCount.begin(); m != member(i)->matHaplotypeCount.end(); ++m) {
			outVector = haplotypeCoder.getVector(m->first);
			os<<"-";
			for (unsigned k=0;k<nLoci;k++){
				os  << outVector[k]+1<<"-";
			} 
			
			os <<"\t\t" << ((double) m->second)/numSamples<< std::endl;
		}
		os << "paternal haplotypes" << std::endl;
		for(m = member(i)->patHaplotypeCount.begin(); m != member(i)->patHaplotypeCount.end(); ++m) {
			outVector = haplotypeCoder.getVector(m->first);
			os<<"-";
			for (unsigned k=0;k<nLoci;k++){
				os << outVector[k]+1<<"-";
			} 
			
			os <<"\t\t" << ((double) m->second)/numSamples<< std::endl;
		}
		
		std::cout << std::endl << std::endl;
		
	}
}
/*! \fn void Population::displayHaplotypeFrequncies(unsigned numSamples)
*  \brief calculate and display haplotype frequencies 
*/

void Population::countGenotypes(string samplerType){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	
	unsigned nLoci = Individual::numLoci;
	unsigned mat, pat, geno;
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for (unsigned j=0; j<nLoci; j++){  
			if (samplerType=="genotypic"){
				mat = ind->genotNodeVector[j].getAcceptedMatState();
				pat = ind->genotNodeVector[j].getAcceptedPatState();
			}
			else if (samplerType=="allelic"){
				mat = ind->malleleStateNodeVector[j].getAcceptedAlleleState();
				pat = ind->palleleStateNodeVector[j].getAcceptedAlleleState();
			}
			else if(samplerType=="joint"){
				if (prior->chrom()[0].locus[j].gnodeType=="allelic"){
					mat  = ind->malleleStateNodeVector[j].getAcceptedAlleleState();
					pat  = ind->palleleStateNodeVector[j].getAcceptedAlleleState();
				}
				else {
					mat = ind->genotNodeVector[j].getAcceptedMatState();
					pat = ind->genotNodeVector[j].getAcceptedPatState();
				}	
			}
			geno = mat*prior->chromosome[0].locus[j].nallele() + pat;
			ind->genotNodeVector[j].genotypeCount[geno]++;
			ind->genotNodeVector[j].sampleCount++;
		}
	}
}
/*! \fn void Population::countHaplotypes(samplerType)
*  \brief count the number genotypes sampled for each individual 
*/

void Population::countGenotypesSimple(string samplerType){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	
	unsigned nLoci = Individual::numLoci;
	unsigned mat, pat, geno;
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for (unsigned j=0; j<nLoci; j++){  
			if (samplerType=="genotypic"){
				mat = ind->genotNodeVector[j].getmState();
				pat = ind->genotNodeVector[j].getpState();
			}
			else if (samplerType=="allelic"){
				mat = ind->malleleStateNodeVector[j].getMyAlleleState();
				pat = ind->palleleStateNodeVector[j].getMyAlleleState();
			}
			else if(samplerType=="joint"){
				if (prior->chrom()[0].locus[j].gnodeType=="allelic"){
					mat  = ind->malleleStateNodeVector[j].getMyAlleleState();
					pat  = ind->palleleStateNodeVector[j].getMyAlleleState();
				}
				else {
					mat = ind->genotNodeVector[j].getmState();
					pat = ind->genotNodeVector[j].getpState();
				}	
			}
			geno = mat*prior->chromosome[0].locus[j].nallele() + pat;
			ind->genotNodeVector[j].genotypeCount[geno]++;
			ind->genotNodeVector[j].sampleCount++;
		}
	}
}
/*! \fn void Population::countHaplotypes(samplerType)
*  \brief count the number genotypes sampled for each individual 
*/


void Population::displayGenotypeFrequencies(unsigned numSamples,  std::ostream &outfile){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors:

	Individual *ind;
	unsigned nLoci = Individual::numLoci;
	for (unsigned j=0;j<nLoci;j++){
		unsigned numAlleles = prior->chromosome[0].locus[j].nallele();
		outfile << "Locus: " << j+1 << endl;
		for(unsigned i=0;i<popsize;i++){
			ind=popmember[i];
			outfile << "Genotype (maternal/paternal) for: " << setw(5) << ind_name(i+1) << endl;
			outfile << " ---------------------------- " << endl;
			for (unsigned mat=0;mat<numAlleles;mat++){
				for (unsigned pat=0;pat<numAlleles;pat++){
					unsigned geno = mat*prior->chromosome[0].locus[j].nallele() + pat;
					if (ind->genotNodeVector[j].genotypeCount[geno]){
						outfile <<setw(2)   << mat+1   << "/"   <<setw(2) << pat+1 << "  "
						<<setw(10)  << ind->genotNodeVector[j].genotypeCount[geno]/
						               (double) ind->genotNodeVector[j].sampleCount 
						<<setw(10)  << ind->genotNodeVector[j].sampleCount << endl;
					}
				}
			}
			outfile << " ---------------------------- " << endl;
		}
	}
}




void Population::outputQ2Probs(unsigned numSamples,  std::ostream &outfile){
	// Authors: Rohan L. Fernando 
	// (December, 2005) 
	// Contributors:
	// For Avaigen Project; Not intended for general use

	Individual *ind;
	double Q00,Q01,Q10,Q11;
	for(unsigned i=0;i<popsize;i++){
			ind=popmember[i];
			Q00  = (ind->genotNodeVector[0].genotypeCount[0]) ?   ind->genotNodeVector[0].genotypeCount[0]/(double)ind->genotNodeVector[0].sampleCount : 0;
			Q01  = (ind->genotNodeVector[0].genotypeCount[1]) ?   ind->genotNodeVector[0].genotypeCount[1]/(double)ind->genotNodeVector[0].sampleCount : 0;
			Q10  = (ind->genotNodeVector[0].genotypeCount[2]) ?   ind->genotNodeVector[0].genotypeCount[2]/(double)ind->genotNodeVector[0].sampleCount : 0;
			Q11  = (ind->genotNodeVector[0].genotypeCount[3]) ?   ind->genotNodeVector[0].genotypeCount[3]/(double)ind->genotNodeVector[0].sampleCount : 0;

			outfile << setw(15) << ind_name(i+1) << " " 
			        << setw(10) << Q00           << " " 
			        << setw(10) << Q01           << " "
			        << setw(10) << Q10           << " " 
			        << setw(10) << Q11           << endl;
	}
}

void Population::initJointAlleleNodeList(unsigned howManyLoci){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (August, 2004) 
	// Contributors: 
	renum();
	Individual *ind, *mom, *dad;
	gNodeList.resize(2*howManyLoci*(popsize + popsize-numFounders));
	gNodeList.howToSample = "joint";
	GNode::gNodeListPtr = &gNodeList;
//	vectorOfGNsts.clear();
	gNodeList.releaseGNsts();
	for (unsigned jj=0;jj<howManyLoci;jj++){
		unsigned locus = prior->chrom()[0].locusMask[jj].index; 
		unsigned i,j,k;
		unsigned startPosition = (popsize + popsize-numFounders)*2*locus;
		unsigned countOffspring = 0;
		for(i=0;i<popsize;i++){
			j=startPosition + i*2;
			ind=popmember[i];
			mom=ind->mymother;
			dad=ind->myfather;
			AlleleStateNode* imAlleleState = &ind->malleleStateNodeVector[locus];
			AlleleStateNode* ipAlleleState = &ind->palleleStateNodeVector[locus];
			ipAlleleState->id = j+1; 
			ipAlleleState->connectFlag = j+1;
			ipAlleleState->numberOfCuts= 0; 
			imAlleleState->id = j+2;
			imAlleleState->connectFlag = j+2;
			imAlleleState->numberOfCuts= 0;
			gNodeList[j]   = ipAlleleState;
			gNodeList[j+1] = imAlleleState;
			
			if(prior->chrom()[0].locus[locus].qtl_ml=='r'){
				RDisAllelePenetranceSet *aPenSet = new RDisAllelePenetranceSet;
				aPenSet->forLocus = locus;
				aPenSet->owner = ind;
//				vectorOfGNsts.push_back(aPenSet);
				aPenSet->insert(imAlleleState);
				aPenSet->insert(ipAlleleState);
				gNodeList.completeSetofGNsts.insert(aPenSet);
				aPenSet->attachMeToMyGnodes();
			}
			else{
				RAllelePenetranceSet *aPenSet = new RAllelePenetranceSet;
				aPenSet->forLocus = locus;
				aPenSet->owner = ind;
//				vectorOfGNsts.push_back(aPenSet);
				aPenSet->insert(imAlleleState);
				aPenSet->insert(ipAlleleState);
				gNodeList.completeSetofGNsts.insert(aPenSet);
				aPenSet->attachMeToMyGnodes();
			}
			
			if(mom){
				k=startPosition + 2*popsize + countOffspring*2;
				countOffspring++;
				AlleleOriginNode *imAlleleOrigin = &ind->malleleOriginNodeVector[locus];
				AlleleOriginNode *ipAlleleOrigin = &ind->palleleOriginNodeVector[locus];
				ipAlleleOrigin->id = k+1;
				ipAlleleOrigin->connectFlag = k+1;
				ipAlleleOrigin->numberOfCuts= 0;
				imAlleleOrigin->id = k+2;
				imAlleleOrigin->connectFlag = k+2;
				imAlleleOrigin->numberOfCuts= 0;
				
				gNodeList[k]   = ipAlleleOrigin;
				gNodeList[k+1] = imAlleleOrigin;
				
				RecombinationSet *mRecomSet = new RecombinationSet;
				mRecomSet->r=0.5;
				mRecomSet->insert(imAlleleOrigin);
				if(locus){
					AlleleOriginNode *imLeftAlleleOrigin = &ind->malleleOriginNodeVector[locus-1];
					mRecomSet->r=recombinationMatrix[locus-1][locus];
					mRecomSet->insert(imLeftAlleleOrigin);
				}
//				vectorOfGNsts.push_back(mRecomSet);
				gNodeList.completeSetofGNsts.insert(mRecomSet);
				mRecomSet->attachMeToMyGnodes();
				
				RecombinationSet *pRecomSet = new RecombinationSet;
				pRecomSet->r=0.5;
				pRecomSet->insert(ipAlleleOrigin);
				if(locus){
					AlleleOriginNode *ipLeftAlleleOrigin = &ind->palleleOriginNodeVector[locus-1];
					pRecomSet->r=recombinationMatrix[locus-1][locus];
					pRecomSet->insert(ipLeftAlleleOrigin);
				}
//				vectorOfGNsts.push_back(pRecomSet);
				gNodeList.completeSetofGNsts.insert(pRecomSet);
				pRecomSet->attachMeToMyGnodes();
				
				AlleleStateNode *mMAlleleState = &mom->malleleStateNodeVector[locus];
				AlleleStateNode *mPAlleleState = &mom->palleleStateNodeVector[locus];
				AlleleStateNode *pMAlleleState = &dad->malleleStateNodeVector[locus];
				AlleleStateNode *pPAlleleState = &dad->palleleStateNodeVector[locus];
				
				RTransmissionSet *mTransmSet = new RTransmissionSet;
				mTransmSet->forLocus = locus;
				mTransmSet->offspring = ind;
				mTransmSet->paternal = false;
				mTransmSet->insert(imAlleleState);
				mTransmSet->insert(imAlleleOrigin);
				mTransmSet->insert(mMAlleleState);
				mTransmSet->insert(mPAlleleState);
//				vectorOfGNsts.push_back(mTransmSet);
				gNodeList.completeSetofGNsts.insert(mTransmSet);
				mTransmSet->attachMeToMyGnodes();
				
				RTransmissionSet *pTransmSet = new RTransmissionSet;
				pTransmSet->forLocus = locus;
				pTransmSet->offspring = ind;
				pTransmSet->paternal = true;
				pTransmSet->insert(ipAlleleState);
				pTransmSet->insert(ipAlleleOrigin);
				pTransmSet->insert(pMAlleleState);
				pTransmSet->insert(pPAlleleState);
//				vectorOfGNsts.push_back(pTransmSet);
				gNodeList.completeSetofGNsts.insert(pTransmSet);
				pTransmSet->attachMeToMyGnodes();
			}
			else {
				RAlleleFounderSet *mFoundSet = new RAlleleFounderSet;
				mFoundSet->forLocus = locus;
				mFoundSet->insert(imAlleleState);
//				vectorOfGNsts.push_back(mFoundSet);
				gNodeList.completeSetofGNsts.insert(mFoundSet);
				mFoundSet->attachMeToMyGnodes();
				
				RAlleleFounderSet *pFoundSet = new RAlleleFounderSet;
				pFoundSet->forLocus = locus;
				pFoundSet->insert(ipAlleleState);
//				vectorOfGNsts.push_back(pFoundSet);
				gNodeList.completeSetofGNsts.insert(pFoundSet);
				pFoundSet->attachMeToMyGnodes();
			}
		}
	}
}
/*! \fn void Population::initjointAlleleNodeList(unsigned howManyLoci)
*  \brief creates the allele node list, and the founder, penetrance,
and transmission sets for each allele node for joint peeling
*/

void Population::getInitialGNodeListSample(unsigned maxsize, unsigned numLoci, string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: 
	if (samplerType=="genotypic"){
		cerr << "For joint peeling, the sampler type must be allelic;" << endl;
		exit(1);
	} 
	else if (samplerType=="allelic"){
		prior->chrom()[0].peelOrder.resize(2*numLoci*(popsize + popsize-numFounders));
		initJointAlleleNodeList(numLoci);
		gNodeList.displayGNodeSets();
	}
	else {
		cerr << "Sampler type should be allelic;" << endl;
		exit(1);
	}
	bool notSampled = true;
	while(notSampled){
		try { 
			gNodeList.peelOrderCutAndSample(maxsize); 
			notSampled = false;
		}
		catch (matvec::InvalidSample) {
			cout << "Inconsistent sample, I will try to sample again" << endl;
			gNodeList.reinitGNodeList();
		} 
	}
}
/*! \fn void Population::getInitialGNodeListSample(unsigned maxsize, unsigned numLoci, string samplerType[])
*  \brief get the initial sample and also determine and store the peeling and sampling order for peeling across
the pedigree as well as across the loci
*/

void Population::getGNodeListSample(unsigned maxsize, unsigned numLoci, string samplerType){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (September, 2004) 
	// Contributors: 
	if(samplerType=="genotypic"){
		cerr << "For joint peeling, the sampler type must be allelic;" << endl;
		exit(1);
	} 
	else if(samplerType=="allelic"){
		initJointAlleleNodeList(numLoci);
	}
	else {
		cerr << "Sampler type should be allelic;" << endl;
		exit(1);
	}
	bool notSampled = true;
	while(notSampled){	
		try { 
			gNodeList.peelCutAndSample(maxsize); 
			notSampled = false;
		}
		catch (matvec::InvalidSample) {
			cout << "Inconsistent sample, I will try to sample again" << endl;
			gNodeList.reinitGNodeList();
		} 
		//setSegregationIndex(j,samplerType);
	}
}
/*! \fn void Population::getGNodeListSample(unsigned maxsize, unsigned numLoci, string samplerType[])
*  \brief get a new joint sample across the pedigree as well as loci
*/

void Population::countParentOffspring(string samplerType){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	// This version is only for single locus, which is assumed to be at position 0
	unsigned nLoci = Individual::numLoci;
	unsigned mat, pat, geno, mmat, mpat, mgeno, pmat, ppat, pgeno;
	Individual *ind, *mother, *father;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mother = ind->mymother;
		if (mother == 0) continue;
		father = ind->myfather;
		unsigned j=0; 
		if (samplerType=="genotypic"){
			mat = ind->genotNodeVector[j].getAcceptedMatState();
			pat = ind->genotNodeVector[j].getAcceptedPatState();
			mmat= mother->genotNodeVector[j].getAcceptedMatState();
			mpat= mother->genotNodeVector[j].getAcceptedPatState();
			pmat= father->genotNodeVector[j].getAcceptedMatState();
			ppat= father->genotNodeVector[j].getAcceptedPatState();
		}
		else if (samplerType=="allelic"){
			mat = ind->malleleStateNodeVector[j].getAcceptedAlleleState();
			pat = ind->palleleStateNodeVector[j].getAcceptedAlleleState();
			mmat= mother->malleleStateNodeVector[j].getAcceptedAlleleState();
			mpat= mother->palleleStateNodeVector[j].getAcceptedAlleleState();
			pmat= father->malleleStateNodeVector[j].getAcceptedAlleleState();
			ppat= father->palleleStateNodeVector[j].getAcceptedAlleleState();
		}
		geno = mat*prior->chromosome[0].locus[j].nallele() + pat;
		ind->genotNodeVector[j].genotypeCount[geno]++;
		ind->genotNodeVector[j].sampleCount++;
		if (mmat<mpat){
			unsigned temp = mmat;
			mmat = mpat;
			mpat = temp;
		}
		if (pmat<ppat){
			unsigned temp = pmat;
			pmat = ppat;
			ppat = temp;
		}
		mgeno = mmat*(mmat+1)/2 + mpat;
		pgeno = pmat*(pmat+1)/2 + ppat;
		ind->motherOffspring[mgeno][geno]++;
		ind->fatherOffspring[pgeno][geno]++;
	}
}
/*! \fn void Population::countHaplotypes(samplerType)
*  \brief count the number genotype combinations for offspring-mother and 
*  offspring-father
*/

void Population::countParentOffspringSimple(string samplerType){
	// Authors: Rohan L. Fernando 
	// (October, 2004) 
	// Contributors: 
	// This version is only for single locus, which is assumed to be at position 0
	unsigned nLoci = Individual::numLoci;
	unsigned mat, pat, geno, mmat, mpat, mgeno, pmat, ppat, pgeno;
	Individual *ind, *mother, *father;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		mother = ind->mymother;
		if (mother == 0) continue;
		father = ind->myfather;
		unsigned j=0; 
		if (samplerType=="genotypic"){
			mat = ind->genotNodeVector[j].getmState();
			pat = ind->genotNodeVector[j].getpState();
			mmat= mother->genotNodeVector[j].getmState();
			mpat= mother->genotNodeVector[j].getpState();
			pmat= father->genotNodeVector[j].getmState();
			ppat= father->genotNodeVector[j].getpState();
		}
		else if (samplerType=="allelic"){
			mat = ind->malleleStateNodeVector[j].getMyAlleleState();
			pat = ind->palleleStateNodeVector[j].getMyAlleleState();
			mmat= mother->malleleStateNodeVector[j].getMyAlleleState();
			mpat= mother->palleleStateNodeVector[j].getMyAlleleState();
			pmat= father->malleleStateNodeVector[j].getMyAlleleState();
			ppat= father->palleleStateNodeVector[j].getMyAlleleState();
		}
		geno = mat*prior->chromosome[0].locus[j].nallele() + pat;
		ind->genotNodeVector[j].genotypeCount[geno]++;
		ind->genotNodeVector[j].sampleCount++;
		if (mmat<mpat){
			unsigned temp = mmat;
			mmat = mpat;
			mpat = temp;
		}
		if (pmat<ppat){
			unsigned temp = pmat;
			pmat = ppat;
			ppat = temp;
		}
		mgeno = mmat*(mmat+1)/2 + mpat;
		pgeno = pmat*(pmat+1)/2 + ppat;
		ind->motherOffspring[mgeno][geno]++;
		ind->fatherOffspring[pgeno][geno]++;
	}
}
/*! \fn void Population::countHaplotypes(samplerType)
*  \brief count the number genotype combinations for offspring-mother and 
*  offspring-father
*/
void Population::initParentOffspring(void){
	Individual *ind;
	unsigned zero = 0;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		unsigned j=0;
		unsigned numAlleles = prior->chromosome[0].locus[j].nallele();
		ind->motherOffspring.resize(numAlleles*(numAlleles+1)/2,numAlleles*numAlleles,zero);
		ind->fatherOffspring.resize(numAlleles*(numAlleles+1)/2,numAlleles*numAlleles,zero);
		ind->genotNodeVector[j].sampleCount = 0;
	}
}

void Population::displayParentOffspringProbs(std::ostream &outfile){
	// Authors: Rohan L. Fernando and David Habier 
	// (July, 2006) 
	// Contributors:
	
	Individual *ind;
	unsigned nLoci = Individual::numLoci;
	unsigned j = 0;
	unsigned numAlleles = prior->chromosome[0].locus[j].nallele();
	outfile << "Locus: " << j+1 << endl;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		if (ind->mymother){
			outfile << " mother-offspring probs for: " << setw(5) << ind_name(i+1) << endl;
			outfile <<   ind->motherOffspring;
			outfile << " father-offspring probs for: " << setw(5) << ind_name(i+1) << endl;
			outfile <<   ind->fatherOffspring;
		}
	}
}

void Population::copyGNodeStatesToAcceptedStates(string typeOfGNode){
	// Authors: L. Radu Totir
	// (October, 2003) 
	// Contributors: 
	Individual *ind;
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		for(unsigned j=0;j<Individual::numLoci;j++){
			if(typeOfGNode=="genotypic"){
				ind->genotNodeVector[j].acceptedGenotypeState  = ind->genotNodeVector[j].genotypeState;
				ind->genotNodeVector[j].acceptedGenotypeVector = ind->genotNodeVector[j].genotypeVector;   
			}
			else if(typeOfGNode=="allelic"){
				ind->malleleStateNodeVector[j].acceptedAlleleState = ind->malleleStateNodeVector[j].alleleState;
				ind->palleleStateNodeVector[j].acceptedAlleleState = ind->palleleStateNodeVector[j].alleleState;
				ind->malleleStateNodeVector[j].acceptedAlleleStateVector = ind->malleleStateNodeVector[j].alleleStateVector; 
				ind->palleleStateNodeVector[j].acceptedAlleleStateVector = ind->palleleStateNodeVector[j].alleleStateVector; 
			}
			else {
				cerr << "Sampler type should be either genotypic or allelic;" << endl;
				exit(1);
			}
		}
		ind->m_gameteAccepted = ind->m_gamete;
		ind->p_gameteAccepted = ind->p_gamete;
	}
}
/*! \fn Population::copyGNodeStatesToAcceptedStates(string samplerType)
*  \brief moves the GNode states (alleles or genotypes) to accepted states
*/

void Population::initJointNodeList(unsigned howManyLoci){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (July, 2006) 
	// Contributors: 
	renum();
	Individual *ind, *mom, *dad;
	unsigned numNodes=0;
	for (unsigned i=0;i<howManyLoci;i++){
		if (prior->chrom()[0].locus[i].gnodeType=="allelic"){
			numNodes += 2*(popsize + popsize - numFounders);
		}
		else if(prior->chrom()[0].locus[i].gnodeType=="genotypic"){
			numNodes += popsize + 2*(popsize - numFounders);
		}
		else {
			std::cout << prior->chrom()[0].locus[i].gnodeType << " not recognized for locus " << i+1 << endl;
			throw exception ("Population::initJointNodeList: Error");
		}
	}
	gNodeList.resize(numNodes);
	gNodeList.name = "Population::gNodeList";
	gNodeList.howToSample = "joint";
	GNode::gNodeListPtr = &gNodeList;
	//vectorOfGNsts.clear();
	gNodeList.releaseGNsts();
	unsigned listPos = 0;
	for (unsigned jj=0;jj<howManyLoci;jj++){
		unsigned locus = prior->chrom()[0].locusMask[jj].index; 
		unsigned i,j,k;
		//unsigned countOffspring = 0;
		for(i=0;i<popsize;i++){
			ind=popmember[i];
			mom=ind->mymother;
			dad=ind->myfather;
			if (prior->chrom()[0].locus[locus].gnodeType=="allelic"){
				AlleleStateNode* imAlleleState = &ind->malleleStateNodeVector[locus];
				AlleleStateNode* ipAlleleState = &ind->palleleStateNodeVector[locus];
				ipAlleleState->id = listPos+1; 
				ipAlleleState->connectFlag = listPos+1;
				ipAlleleState->numberOfCuts = 0; 
				ipAlleleState->type = "paternal alleleState ";
				ipAlleleState->locus = &(prior->chrom()[0].locus[locus]);
				imAlleleState->id = listPos+2;
			    imAlleleState->connectFlag = listPos+2;
				imAlleleState->numberOfCuts= 0;
				imAlleleState->type = "maternal alleleState ";
				imAlleleState->locus = &(prior->chrom()[0].locus[locus]);
				
			    gNodeList[listPos++]  = ipAlleleState;
			    gNodeList[listPos++]  = imAlleleState;
			    if(prior->chrom()[0].locus[locus].qtl_ml=='r'){
					RDisAllelePenetranceSet *aPenSet = new RDisAllelePenetranceSet;
					aPenSet->forLocus = locus;
					aPenSet->owner = ind;
					//vectorOfGNsts.push_back(aPenSet);
					aPenSet->insert(imAlleleState);
					aPenSet->insert(ipAlleleState);
					gNodeList.completeSetofGNsts.insert(aPenSet);
					aPenSet->attachMeToMyGnodes();
				}
				else{
					RAllelePenetranceSet *aPenSet = new RAllelePenetranceSet;
					aPenSet->forLocus = locus;
					aPenSet->owner = ind;
					//vectorOfGNsts.push_back(aPenSet);
					aPenSet->insert(imAlleleState);
					aPenSet->insert(ipAlleleState);
					gNodeList.completeSetofGNsts.insert(aPenSet);
					aPenSet->attachMeToMyGnodes();
				}
				if(mom){
					AlleleOriginNode *imAlleleOrigin = &ind->malleleOriginNodeVector[locus];
					AlleleOriginNode *ipAlleleOrigin = &ind->palleleOriginNodeVector[locus];
					ipAlleleOrigin->id = listPos+1;
					ipAlleleOrigin->connectFlag = listPos+1;
					ipAlleleOrigin->numberOfCuts= 0;
					ipAlleleOrigin->type = "paternal alleleOrigin ";
				    ipAlleleOrigin->locus = &(prior->chrom()[0].locus[locus]);
					imAlleleOrigin->id = listPos+2;
					imAlleleOrigin->connectFlag = listPos+2;
					imAlleleOrigin->numberOfCuts= 0;
					imAlleleOrigin->type = "paternal alleleOrigin ";
				    imAlleleOrigin->locus = &(prior->chrom()[0].locus[locus]);
					
					gNodeList[listPos++] = ipAlleleOrigin;
					gNodeList[listPos++] = imAlleleOrigin;
					
					RecombinationSet *mRecomSet = new RecombinationSet;
					mRecomSet->r=0.5;
					mRecomSet->insert(imAlleleOrigin);
					if(locus){
						AlleleOriginNode *imLeftAlleleOrigin = &ind->malleleOriginNodeVector[locus-1];
						mRecomSet->r=recombinationMatrix[locus-1][locus];
						mRecomSet->insert(imLeftAlleleOrigin);
					}
					//vectorOfGNsts.push_back(mRecomSet);
					gNodeList.completeSetofGNsts.insert(mRecomSet);
					mRecomSet->attachMeToMyGnodes();
					
					RecombinationSet *pRecomSet = new RecombinationSet;
					pRecomSet->r=0.5;
					pRecomSet->insert(ipAlleleOrigin);
					if(locus){
						AlleleOriginNode *ipLeftAlleleOrigin = &ind->palleleOriginNodeVector[locus-1];
						pRecomSet->r=recombinationMatrix[locus-1][locus];
						pRecomSet->insert(ipLeftAlleleOrigin);
					}
					//vectorOfGNsts.push_back(pRecomSet);
					gNodeList.completeSetofGNsts.insert(pRecomSet);
					pRecomSet->attachMeToMyGnodes();
					
					AlleleStateNode *mMAlleleState = &mom->malleleStateNodeVector[locus];
					AlleleStateNode *mPAlleleState = &mom->palleleStateNodeVector[locus];
					AlleleStateNode *pMAlleleState = &dad->malleleStateNodeVector[locus];
					AlleleStateNode *pPAlleleState = &dad->palleleStateNodeVector[locus];
					
					RTransmissionSet *mTransmSet = new RTransmissionSet;
					mTransmSet->forLocus = locus;
					mTransmSet->offspring = ind;
					mTransmSet->paternal = false;
					mTransmSet->insert(imAlleleState);
					mTransmSet->insert(imAlleleOrigin);
					mTransmSet->insert(mMAlleleState);
					mTransmSet->insert(mPAlleleState);
					//vectorOfGNsts.push_back(mTransmSet);
					gNodeList.completeSetofGNsts.insert(mTransmSet);
					mTransmSet->attachMeToMyGnodes();
					
					RTransmissionSet *pTransmSet = new RTransmissionSet;
					pTransmSet->forLocus = locus;
					pTransmSet->offspring = ind;
					pTransmSet->paternal = true;
					pTransmSet->insert(ipAlleleState);
					pTransmSet->insert(ipAlleleOrigin);
					pTransmSet->insert(pMAlleleState);
					pTransmSet->insert(pPAlleleState);
					//vectorOfGNsts.push_back(pTransmSet);
					gNodeList.completeSetofGNsts.insert(pTransmSet);
					pTransmSet->attachMeToMyGnodes();
				}
				else {
					RAlleleFounderSet *mFoundSet = new RAlleleFounderSet;
					mFoundSet->forLocus = locus;
					mFoundSet->insert(imAlleleState);
					//vectorOfGNsts.push_back(mFoundSet);
					gNodeList.completeSetofGNsts.insert(mFoundSet);
					mFoundSet->attachMeToMyGnodes();
					
					RAlleleFounderSet *pFoundSet = new RAlleleFounderSet;
					pFoundSet->forLocus = locus;
					pFoundSet->insert(ipAlleleState);
					//vectorOfGNsts.push_back(pFoundSet);
					gNodeList.completeSetofGNsts.insert(pFoundSet);
					pFoundSet->attachMeToMyGnodes();
				}
			}
			else if (prior->chrom()[0].locus[locus].gnodeType=="genotypic"){
				GenotypeNode* indGenotype = &ind->genotNodeVector[locus];
				indGenotype->id = listPos+1;
				indGenotype->connectFlag = listPos+1;
				indGenotype->numberOfCuts= 0;
				indGenotype->type = "genotypeState";
				indGenotype->locus = &(prior->chrom()[0].locus[locus]);				
				gNodeList[listPos++]   = indGenotype;
				if(prior->chrom()[0].locus[locus].qtl_ml=='q'){
					GenoPenetranceSet *gPenSet = new GenoPenetranceSet;
					gPenSet->owner = ind;
					gPenSet->insert(indGenotype);
					//vectorOfGNsts.push_back(gPenSet);
					gNodeList.completeSetofGNsts.insert(gPenSet);
					gPenSet->attachMeToMyGnodes();
				}
				else if(prior->chrom()[0].locus[locus].qtl_ml=='r'){
					RDisGenoPenetranceSet *gPenSet = new RDisGenoPenetranceSet;
					gPenSet->owner = ind;
					gPenSet->forLocus = locus;
					gPenSet->insert(indGenotype);
					//vectorOfGNsts.push_back(gPenSet);
					gNodeList.completeSetofGNsts.insert(gPenSet);
					gPenSet->attachMeToMyGnodes();
				} 
				if(mom){
					AlleleOriginNode *imAlleleOrigin = &ind->malleleOriginNodeVector[locus];
					AlleleOriginNode *ipAlleleOrigin = &ind->palleleOriginNodeVector[locus];
					ipAlleleOrigin->id = listPos+1;
					ipAlleleOrigin->connectFlag = listPos+1;
					ipAlleleOrigin->numberOfCuts= 0;
					ipAlleleOrigin->type = "paternal alleleOrigin ";
				    ipAlleleOrigin->locus = &(prior->chrom()[0].locus[locus]);					
					imAlleleOrigin->id = listPos+2;
					imAlleleOrigin->connectFlag = listPos+2;
					imAlleleOrigin->numberOfCuts= 0;
					imAlleleOrigin->type = "maternal alleleOrigin ";
				    imAlleleOrigin->locus = &(prior->chrom()[0].locus[locus]);					
					gNodeList[listPos++] = ipAlleleOrigin;
					gNodeList[listPos++] = imAlleleOrigin;
					
					RecombinationSet *mRecomSet = new RecombinationSet;
					mRecomSet->r=0.5;
					mRecomSet->insert(imAlleleOrigin);
					if(locus){
						AlleleOriginNode *imLeftAlleleOrigin = &ind->malleleOriginNodeVector[locus-1];
						mRecomSet->r=recombinationMatrix[locus-1][locus];
						mRecomSet->insert(imLeftAlleleOrigin);
					}
					//vectorOfGNsts.push_back(mRecomSet);
					gNodeList.completeSetofGNsts.insert(mRecomSet);
					mRecomSet->attachMeToMyGnodes();
					
					RecombinationSet *pRecomSet = new RecombinationSet;
					pRecomSet->r=0.5;
					pRecomSet->insert(ipAlleleOrigin);
					if(locus){
						AlleleOriginNode *ipLeftAlleleOrigin = &ind->palleleOriginNodeVector[locus-1];
						pRecomSet->r=recombinationMatrix[locus-1][locus];
						pRecomSet->insert(ipLeftAlleleOrigin);
					}
					//vectorOfGNsts.push_back(pRecomSet);
					gNodeList.completeSetofGNsts.insert(pRecomSet);
					pRecomSet->attachMeToMyGnodes();
					
					GenotypeNode *momGenotype = &mom->genotNodeVector[locus];
					GenotypeNode *dadGenotype = &dad->genotNodeVector[locus];
					RTransitionSet *TransitSet = new RTransitionSet; 
					TransitSet->offspring = ind;
					TransitSet->insert(indGenotype);
					TransitSet->insert(momGenotype);
					TransitSet->insert(dadGenotype);
					TransitSet->insert(imAlleleOrigin);
					TransitSet->insert(ipAlleleOrigin);
					TransitSet->forLocus = locus;
					//vectorOfGNsts.push_back(TransitSet);
					gNodeList.completeSetofGNsts.insert(TransitSet);
					TransitSet->attachMeToMyGnodes();
				}
				else {
					RGenoFounderSet *FoundSet = new RGenoFounderSet;
					FoundSet->insert(indGenotype);
					//vectorOfGNsts.push_back(FoundSet);
					gNodeList.completeSetofGNsts.insert(FoundSet);
					FoundSet->attachMeToMyGnodes();
				}
			}
			//ind->displayAlleleVectors(locus);
		}
	}
}
/*! \fn void Population::initjointAlleleNodeList(unsigned howManyLoci)
*  \brief creates the allele node list, and the founder, penetrance,
and transmission sets for each allele node for joint peeling
*/

void Population::displayGenotypeProbs(std::ostream &outfile){
	// Authors: Rohan L. Fernando and Radu Totir
	// (July 2006) 
	// Contributors:

	Individual *ind;
	unsigned nLoci = Individual::numLoci;
	for (unsigned j=0;j<nLoci;j++){
	    if (prior->chrom()[0].locus[j].gnodeType=="allelic") continue;
		unsigned numAlleles = prior->chromosome[0].locus[j].nallele();
		outfile << "Locus: " << j+1 << endl;
		for(unsigned i=0;i<popsize;i++){
			ind=popmember[i];
			outfile << "Genotype (maternal/paternal) for: " << setw(5) << ind_name(i+1) << endl;
			outfile << " ---------------------------- " << endl;
			for (unsigned mat=0;mat<numAlleles;mat++){
				for (unsigned pat=0;pat<numAlleles;pat++){
					unsigned geno = mat*prior->chromosome[0].locus[j].nallele() + pat;
					if (ind->genotNodeVector[j].genotypeCount[geno]){
						outfile <<setw(2)   << mat+1   << "/"   <<setw(2) << pat+1 << "  "
						<<setw(10)  << ind->genotNodeVector[j].genotypeCount[geno] << endl;
					}
				}
			}
			outfile << " ---------------------------- " << endl;
		}
	}
}

void Population::displayGNodeProbs(std::ostream &outfile){
	// Authors: Rohan L. Fernando
	// (September 2006) 
	// Contributors:

	Individual *ind;
	unsigned nLoci = Individual::numLoci;
	for (unsigned j=0;j<nLoci;j++){
	    if (prior->chrom()[0].locus[j].originProbs){
			calcAndDisplayOriginProbsForLocus(j, outfile);
		}
		if (prior->chrom()[0].locus[j].stateProbs){
			if (prior->chrom()[0].locus[j].gnodeType=="allelic"){
				calcAndDisplayAlleleProbsForLocus(j, outfile);
			}
			else if (prior->chrom()[0].locus[j].gnodeType=="genotypic"){
				calcAndDisplayGenotypeProbsForLocus(j, outfile);
			}
		}
	}
}		

void Population::calcAndDisplayOriginProbsForLocus(unsigned j, std::ostream &outfile){
	// Authors: Rohan L. Fernando 
	// (September 2006) 
	// Contributors:
	
	Individual *ind;
	CutSet probSet;
	double matProb, patProb;
	outfile << "Grand Maternal Origin Probabilities for locus " << j+1 << endl; 
	outfile << "individual                maternal allele paternal allele \n";
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		if(!ind->mother()) continue;
		probSet = ind->malleleOriginNodeVector[j].calcGNodeProbs();
		unsigned origin0 = ind->malleleOriginNodeVector[j].alleleOriginVector[0];
		matProb = origin0 ? 1-probSet.valueVector[0] : probSet.valueVector[0];
		probSet = ind->palleleOriginNodeVector[j].calcGNodeProbs();
		origin0 = ind->palleleOriginNodeVector[j].alleleOriginVector[0];
		patProb = origin0 ? 1-probSet.valueVector[0] : probSet.valueVector[0];
		outfile << setw(25) << setiosflags (ios::left) << ind_name(i+1) << " ";
		outfile << setw(15) << setprecision (5)  << setiosflags (ios::right | ios::fixed) << matProb <<" ";
		outfile << setw(15) << setprecision (5)  << setiosflags (ios::right | ios::fixed) << patProb << endl;		
		outfile << resetiosflags(ios::right | ios::fixed); 
	}
}

void Population::calcAndDisplayAlleleProbsForLocus(unsigned j, std::ostream &outfile){
	// Authors: Rohan L. Fernando 
	// (September 2006) 
	// Contributors:
	
	Individual *ind;
	CutSet probSet;
	outfile << "Allele State Probabilities for locus " << j+1 << endl; 
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		//ind->displayAlleleVectors(2);
		probSet = ind->palleleStateNodeVector[j].calcGNodeProbs();
		outfile << "paternal allele probabilities for individual: " << ind_name(i+1) <<endl;
		outfile << "allele           probability \n";
		for (unsigned allele=0;allele<ind->palleleStateNodeVector[j].alleleStateVector.size();allele++){
		    outfile << setw(8) << setiosflags (ios::left);
			outfile << ind->palleleStateNodeVector[j].alleleStateVector[allele]<< " "; 
			outfile << setw(15) << setprecision (5) << setiosflags (ios::right | ios::fixed);  
			outfile << probSet.valueVector[allele] << endl;
			outfile << resetiosflags(ios::right | ios::fixed); 
		}
		probSet = ind->malleleStateNodeVector[j].calcGNodeProbs();
		outfile << "maternal allele probabilities for individual: " << ind_name(i+1) <<endl;
		outfile << "allele           probability \n";
		for (unsigned allele=0;allele<ind->malleleStateNodeVector[j].alleleStateVector.size();allele++){
		    outfile << setw(8) << setiosflags (ios::left);
			outfile << ind->malleleStateNodeVector[j].alleleStateVector[allele] << " "; 
			outfile << setw(15) << setprecision (5) << setiosflags (ios::right | ios::fixed);
			outfile << probSet.valueVector[allele] << endl;
			outfile << resetiosflags(ios::right | ios::fixed); 
		}
 
	}
}

void Population::calcAndDisplayGenotypeProbsForLocus(unsigned j, std::ostream &outfile){
	// Authors: Rohan L. Fernando 
	// (September 2006) 
	// Contributors:
	
	Individual *ind;
	CutSet probSet;
	outfile << "Genotype Probabilities for locus " << j+1 << endl; 
	outfile << "individual               maternal allele paternal allele \n";
	for(unsigned i=0;i<popsize;i++){
		ind=popmember[i];
		probSet = ind->genotNodeVector[j].calcGNodeProbs();
		outfile << "genotype probabilities for individual: " << ind_name(i+1) <<endl;
		outfile << "maternal allele : paternal allele       probability \n";
 		for (unsigned geno=0;geno<ind->genotNodeVector[j].genotypeVector.size();geno++){
			unsigned matAllele = ind->genotNodeVector[j].genotypeVector[geno].maternal;
			unsigned patAllele = ind->genotNodeVector[j].genotypeVector[geno].paternal;
			outfile << setw(8) << setiosflags (ios::right) <<  matAllele;
			outfile << setw(18) << setiosflags (ios::right) << patAllele;
			outfile << setw(22) << setprecision (5) << setiosflags (ios::right | ios::fixed); 
			outfile << probSet.valueVector[geno] << endl;
		}
	}
}

void Population::calcDistanceToIndividual(std::string strId){
	// Authors: Rohan L. Fernando 
	// (October 2006) 
	// Contributors:
	
    Individual* ind;
	if (!spouse_info_built){
		build_spouse_info();
	}
	if (myRPedPtr==0){
		throw exception("Population::calcDistanceToIndividual: RPedigree is empty");
	}
	for(unsigned i=0;i<popsize;i++){ // initialize distance to "not calculated"
		ind=popmember[i];
		ind->distanceToPivot0 = -1;
	}
	SafeSTLVector<unsigned> individualsToBeProcessed;
	RPedigree::iterator pedIt;
	pedIt = myRPedPtr->find(strId);
	if(pedIt==myRPedPtr->end()){
	    cout << "Could not find " << strId << "in Pedigree" << endl;
		throw exception("Population::calcDistanceToIndividual: Error");
	}
	unsigned pivot0ID = pedIt->second->ind - 1;
	unsigned currentInd = 0;
	ind=popmember[pivot0ID];
	ind->distanceToPivot0 = 0;
	individualsToBeProcessed.push_back(pivot0ID);
	while (currentInd<individualsToBeProcessed.size()) {
		unsigned id = individualsToBeProcessed[currentInd];
		ind=popmember[id];
		ind->calcDistanceToNeighbors(individualsToBeProcessed);
		currentInd++;
	}
//	for(unsigned i=0;i<popsize;i++){ 
//		ind=popmember[i];
//		cout << ind->myid << " " << ind->distanceToPivot0 << endl;
//	}
}


} ////////// end of namespace matvec
