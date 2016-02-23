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

#include <list>
#include <iomanip>
#include "session.h"
#include "population.h"
#include "nufamily.h"

namespace matvec {

extern double penetrance(const Individual* I,const double **res_var);

int compare_seq_id(const void *a, const void *b)
{
   NuFamily **x,**y;
   x = (NuFamily **)a; y = (NuFamily **)b;
   return (x[0]->seq_id - y[0]->seq_id);
}

static int compare_nuf_ptr(const void *a, const void *b)
{
   NuFamily **x,**y;
   x = (NuFamily **)a; y = (NuFamily **)b;
   if (x[0] == y[0]) return 0;
   else return 1;
}
// This is Tianlin's maxant_maxpost
void Population::maxant_maxpost_old(void)
{
   NuFamily *family;
   Individual *I, *father, *mother;
   int i,j;
   doubleMatrix wspace(tn_genotype,tn_genotype);
   Vector<double> genotype_freq = get_genotype_freq();

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      I->initial_anterior(genotype_freq.begin(),tn_genotype);
      I->initial_posterior(tn_genotype);
   }
   genotype_dist_peeling(0);
//   for (j=0; j<6; j++) {
//      for (i=num_t_nufamily; i<num_nufamily; i++) {
//          nufamily_vec[i]->iterative_peeling(wspace);
//      }
//   }
   double Pryi,Pry,Prya,Pryp,val1,val2,val3;
   Vector<double> pen_vec(tn_genotype);
   Vector<double> pos_vec(tn_genotype);
Vector<double> gfreq = genotype_freq;
   for (i=0; i<popsize; i++) {
      I = popmember[i];
      if (I->group_id > 1) {    // I is a connector

get_posterior(I,0,pos_vec);
for (Pry=0.0,j=0; j<tn_genotype; j++) Pry += gfreq[j]=I->anterior[j]*pos_vec[j];
gfreq /= Pry;

         I->get_penetrance(pen_vec);
         for (Prya=0.0,Pryi=0.0,j=0; j<tn_genotype; j++) {
            if (pen_vec[j] > 0.0) {
               Prya += I->anterior[j] /= pen_vec[j];
               Pryi  += pen_vec[j] *= genotype_freq[j];
            }
            else {
               I->anterior[j] = 0.0;
            }
         }
         for (j=0; j<tn_genotype; j++) {
             I->anterior[j] /= Prya;  pen_vec[j] /= Pryi;
         }
// Pry = Pr(y); Prya = Pr(ya); Pryi = Pr(yi); Pryp = Pr(yp)
// I->anterior[j] = Pr(ui|ya)
// pen_vec[j] = Pr(ui|yi)
// gfreq[j]  = Pr(ui|y)

         val1 = I->anterior.max();
         val3 = pen_vec.max();
         I->est_GV = val1/val3;
	 //std::cout << ind_name(I->id()) << ": " << I->est_GV << "\n";
      }
   }
   for (i=num_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      family->in_connector = 0;           // re-initialize
      family->out_connector = 0;

      father = family->father();
      mother = family->mother();
family->workvec[0]=-1.0;
family->workvec[1]=-1.0;
      if (father->group_id > 1) {
         father->get_penetrance(pen_vec);
         get_posterior(father,mother, pos_vec);

for (Pry=0.0,j=0; j<tn_genotype; j++) {
   Pry += gfreq[j] = father->anterior[j]*pos_vec[j];
}
gfreq /= Pry;

         for (Pryi=0.0,Pryp=0.0,j=0; j<tn_genotype; j++) {
            Pryi  += pen_vec[j] *= genotype_freq[j];
            Pryp += pos_vec[j] *= genotype_freq[j];
         }
         for (j=0; j<tn_genotype; j++) {
             pos_vec[j] /= Pryp; pen_vec[j] /= Pryi;
         }
         val2 = pos_vec.max();
         val3 = pen_vec.max();
         family->workvec[0] = val2/val3;

      }
      if (mother->group_id > 1) {
         mother->get_penetrance(pen_vec);
         get_posterior(mother,father, pos_vec);

for (Pry=0.0,j=0; j<tn_genotype; j++) {
   Pry += gfreq[j] = mother->anterior[j]*pos_vec[j];
}
gfreq /= Pry;

         for (Pryi=0.0,Pryp=0.0,j=0; j<tn_genotype; j++) {
            Pryi  += pen_vec[j] *= genotype_freq[j];
            Pryp += pos_vec[j] *= genotype_freq[j];
         }
         for (j=0; j<tn_genotype; j++) {
             pos_vec[j] /= Pryp; pen_vec[j] /= Pryi;
         }
         val2 = pos_vec.max();
         val3 = pen_vec.max();
         family->workvec[1] = val2/val3;
      }
std::cout << father->id() << " " << mother->id() << ": "
     << family->workvec[0] << " " << family->workvec[1] << "\n";
   }

   std::list<Individual *>::iterator it;
   ///////////////////////////////////////////////////////////
   // we will build related_family list for each connectors
   //////////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      for (it = family->connectors.begin(); it != family->connectors.end(); ++it) {
         (*it)->related_family.push_back(family);
      }
   }
}

void Population::build_trans_mat(void)
{
  //This if statement appears in the old version
   if (prior->chrom()[0].nloci() != 1) {
     trans_mat = new doubleMatrix[1];
     return;
   }
   if (trans_mat) return;
   if (!pop_gamete) build_pop_gamete();
   if(tn_genotype>0){
     trans_mat = new doubleMatrix[tn_genotype];
   }
   else {
     trans_mat = 0;
   }
   check_ptr(trans_mat);
   unsigned s,s1,s2,d,d1,d2,p,p1,p2,ss1,ss2;
   d = prior->nchrom();
   for (s=0; s<d; s++) {
      if (prior->chrom()[s].nloci() != 1) {
	if(trans_mat){
	  delete [] trans_mat;
	  trans_mat=0;
	}
         throw exception("Population::build_trans_mat(): one locus/chromosome is alowed");
      }
   }
   for (s=0; s<tn_genotype; s++) trans_mat[s].resize(tn_genotype,tn_genotype,0.0);
   for (s=0,s1=0; s1<tn_gamete; s1++) {
      ss1 = s1*(s1+1)/2;
      for (s2=0; s2<=s1; s2++) {
         ss2 = s2*(s2+1)/2;
         for (d=0,d1=0; d1<tn_gamete; d1++) {
            if (s1 >= d1) {
               p1 = ss1 + d1;
            }
            else {
               p1 = d1*(d1+1)/2 + s1;
            }
            if (s2 >= d1) {
               p2 = ss2 + d1;
            }
            else {
               p2 = d1*(d1+1)/2 + s2;
            }
            for (d2=0; d2<=d1; d2++) {
               trans_mat[p1][s][d] += 0.25; trans_mat[p2][s][d] += 0.25;
               if (s1 >= d2) {
                  p = ss1 + d2;
               }
               else {
                  p = d2*(d2+1)/2 + s1;
               }
               trans_mat[p][s][d] += 0.25;
               if (s2 >= d2) {
                  p = ss2 + d2;
               }
               else {
                  p = d2*(d2+1)/2 + s2;
               }
               trans_mat[p][s][d] += 0.25;
               d++;
            }
         }
         s++;
      }
   }
}

void Population::build_nufamily(void)
{
	if (!spouse_info_built){
		build_spouse_info();
	}
	Individual *I;
	unsigned i,j,nsp;
	num_t_nufamily = 0;
	num_nufamily = 0;
	for (i=0; i<popsize; i++) {
		I = popmember[i];
		if (I->father() && !I->mother() || !I->father() && I->mother()) {
			j = I->id();
			throw exception(std::string("one of parent is missing: ") + ind_name(j));
		}
		nsp = I->nspouse();
		if (I->sex() == 'F' && nsp > 0) num_nufamily += nsp;
	}
	if (num_nufamily == 0) return;
	nufamily_vec = new NuFamily*[num_nufamily];
	check_ptr(nufamily_vec);
	for (i=0; i<num_nufamily; i++) nufamily_vec[i] = new NuFamily;
	unsigned k,ii,jj,msp;
	Individual *J,**spouses,**offspring;
	NuFamily* nufamily;
	unsigned rec = 0, rec1=0;
	cout << "\n Building nuclear family list \n";
	const unsigned *noffs_sp;
	for (ii=0,i=0; i<popsize; i++) {
		rec++;
		if(rec==1000){
			cout<<rec1+rec <<"\r";			
			cout.flush();
			rec1+=rec;
			rec = 0;
		}
		I = popmember[i];
		nsp = I->nspouse();
		if (I->sex() == 'F' && nsp > 0) {
			offspring = I->offspring();
			noffs_sp = I->noffs_spouse();
			for (k=0, j=0; j<nsp; j++) {
				J = offspring[k]->father();
				//            J = offspring[k]->mother(); //old version
				nufamily = nufamily_vec[ii++];
				nufamily->pop = this;
				nufamily->mother(I);
				nufamily->father(J);
				nufamily->offspring(&(offspring[k]),noffs_sp[j]);
				nufamily->father_indx = j;
				spouses = J->spouses();
				msp = J->nspouse();
				for (jj=0; jj<msp; jj++) if (spouses[jj] == I) {
					nufamily->mother_indx = jj;
					break;
				}
					k += noffs_sp[j];
			}
		}
	}
	//This appears only in the old version
	for (i=0; i<num_nufamily; i++) {
		int noffs = nufamily_vec[i]->noffs();
		offspring = nufamily_vec[i]->offspring();
		for (k=0;k<noffs;k++) offspring[k]->put_family(i);
	}
	nuFamiliesDone = true;
}

int Population::detect_loop(void) {
   //////////////////////////////////////////////////////////////////////////
   //  if an individual's group_id > 1, then it is a connector.
   //  we count # of connectors within each family, each family can have
   //  more than 1 connectors, however, a terminal family only have
   //  1 connector. If a family has no connector at all, it is the
   //  kernal family. an isolated individial is not a nuclei family
   //  if we can't clip out terminal to the end of the pedigree, that means
   //  there is at least one loop in the rest part of the pedigree.
   //
   //  REQUIREMENTS:    build_nufamily() has to be called;
   //////////////////////////////////////////////////////////////////////////

   unsigned i,nc,keep_going = 2;
   std::list<Individual *>::iterator it;
   NuFamily *family;
   unsigned i0 = num_t_nufamily;
   while (keep_going) {
      keep_going = 2;                   // it remains 2 until it's changed
      for (i=i0; i<num_nufamily; i++) {
         family = nufamily_vec[i];
         if ( !family->terminal) {
            nc = family->nconnector();
            if (nc == 0) {                   // found a kernal nuclei family
               family->seq_id = num_t_nufamily++;
               family->kernal = 1;
               family->terminal = 1;
               keep_going = 1;               // reset to 1
            }
            else if (nc == 1) {              // found a terminal nuclei family
               family->seq_id = num_t_nufamily++;
               it = family->connectors.begin();
               (*it)->related_family.remove(family);
               (*it)->group_id -= 1;
               family->terminal = 1;
               keep_going = 1;               // reset to 1
            }
            else {
               family->seq_id = num_nufamily+nc;    // nothing found
            }
            if (num_t_nufamily == num_nufamily) {keep_going = 0; break;}
         }
      }
      if (keep_going == 2) break;   // can't clip farther, loops may be found
   }
   if (num_t_nufamily > i0 ) {
     qsort((char *)nufamily_vec,num_nufamily,sizeof(NuFamily *),compare_seq_id);
   }
   return (num_nufamily-num_t_nufamily);
}

void Population::break_loop(void)
{
   //////////////////////////////////////////////////////////////////////
   // search for any circles by visiting neighbors by neighbors
   // REQUIREMENT:  call detect_loop() first, ie. there must be loops
   /////////////////////////////////////////////////////////////////////
   if (num_t_nufamily == num_nufamily) return;
   NuFamily *pre_family, *cur_family, *neighbor;
   Individual *enter_door, *exit_door;
   std::list<Individual *>::iterator connector;
   std::list<NuFamily *> family_loop_connector;
   std::list<NuFamily *>::iterator Fam;
   std::list<NuFamily *>::reverse_iterator rFam;
   pre_family = 0;
   cur_family = nufamily_vec[num_t_nufamily];
//   cur_family = nufamily_vec[num_nufamily-1];
   cur_family->in_connector = 0;
   while (cur_family->seq_id != 0) {   // seq_id should not be 0 at this moment
      cur_family->seq_id = 0;          // 0 indicates the family is visited
      family_loop_connector.push_back(cur_family);
      neighbor = 0;
      connector = cur_family->connectors.begin();
      while ((connector != cur_family->connectors.end()) && ((*connector) == cur_family->in_connector)) {
         ++connector;
      }
      if (connector == cur_family->connectors.end()) throw exception("Population::break_loop(): a bug!");
      while (connector != cur_family->connectors.end()) {
         Fam = (*connector)->related_family.begin();
         while (Fam != (*connector)->related_family.end()) {
            if (*Fam != cur_family && *Fam != pre_family ) {
               neighbor = *Fam;    // this may not a good neighbor, so chech
               if (neighbor->out_connector != (*connector)) break;
            }
            ++Fam;
         }
         if (neighbor) break;
         ++connector;
      }
      if (!neighbor) {   // connector=0 too, a short circle has been found
         connector= cur_family->connectors.begin();
         while ((connector!= cur_family->connectors.end())&& (*connector)==cur_family->in_connector) {
            ++connector;
         }
         if (connector == cur_family->connectors.end()) throw exception("Population::break_loop(): a bug!");
         cur_family->out_connector = *connector;
         pre_family->in_connector = *connector;
         cur_family = pre_family;
         break;
      }
      cur_family->out_connector = *connector;
      pre_family = cur_family;
      cur_family = neighbor;
      cur_family->in_connector = *connector;
   }
   cur_family->seq_id = 1;  // cur_family is the first family in the circle
   /////////////////////////////////////////////////////////////////////
   // a circle has been found, cur_family is in this circle,
   // let's break it at the position of connector
   /////////////////////////////////////////////////////////////////////
/*
   Individual *C = connector[0];
   int perfectcut = 0;
   double  smallest_gp = 1.0e38;
   Fam = (NuFamily **)family_loop_connector.bottom();
   while (Fam) {
      enter_door = Fam[0]->in_connector;
      exit_door = Fam[0]->out_connector;
      if (enter_door != Fam[0]->mother() && enter_door != Fam[0]->father()){
         if (enter_door->est_GV < smallest_gp) {
            smallest_gp = enter_door->est_GV;
            C = enter_door;
            cur_family = *Fam;
            perfectcut = 1;
         }
      }
      if (exit_door != Fam[0]->mother() && exit_door != Fam[0]->father()){
         if (exit_door->est_GV < smallest_gp) {
            smallest_gp = exit_door->est_GV;
            C = exit_door;
            cur_family = *Fam;
            perfectcut = 1;
         }
      }
      if (Fam[0]->seq_id == 1) break;
      Fam = (NuFamily **)family_loop_connector.previous();
   }
   if (perfectcut == 0) {
      Individual *mother,*father;
      Fam = (NuFamily **)family_loop_connector.bottom();
      while (Fam) {
         enter_door = Fam[0]->in_connector;
         exit_door = Fam[0]->out_connector;
         mother = Fam[0]->mother();
         father = Fam[0]->father();
         if (!enter_door->connect_keeper) {
            if (enter_door == mother) {
               if (Fam[0]->workvec[1] < smallest_gp) {
                  smallest_gp = Fam[0]->workvec[1];
                  C = enter_door;
                  cur_family = *Fam;
               }
            }
            else if (enter_door == father) {
               if (Fam[0]->workvec[0] < smallest_gp) {
                  smallest_gp = Fam[0]->workvec[0];
                  C = enter_door;
                  cur_family = *Fam;
               }
            }
         }
         if (!exit_door->connect_keeper) {
            if (exit_door == mother) {
               if (Fam[0]->workvec[1] < smallest_gp) {
                  smallest_gp = Fam[0]->workvec[1];
                  C = exit_door;
                  cur_family = *Fam;
               }
            }
            else if (exit_door == father) {
               if (Fam[0]->workvec[0] < smallest_gp) {
                  smallest_gp = Fam[0]->workvec[0];
                  C = exit_door;
                  cur_family = *Fam;
               }
            }
         }
         if (Fam[0]->seq_id == 1) break;
         Fam = (NuFamily **)family_loop_connector.previous();
      }
   }
   // wait a minute, let's first update related_family in C
   C->related_family.clear(&cur_family,compare_nuf_ptr);
   cur_family->cutting(C);
*/

   Individual *C = *connector;
   int perfectcut = 0;
   double  smallest_gp = 1.0e38;
   rFam = family_loop_connector.rbegin();
if (stdid >= 0) {
   while (rFam != family_loop_connector.rend()) {
      enter_door = (*rFam)->in_connector;
      exit_door = (*rFam)->out_connector;
      if ((*rFam)->numcut == 0) {
         if (enter_door != (*rFam)->mother() && enter_door != (*rFam)->father()) {
            if (enter_door->est_GV < smallest_gp) {
               smallest_gp = enter_door->est_GV;
               C = enter_door;
               cur_family = *rFam;
               perfectcut = 1;
            }
         }
         if (exit_door != (*rFam)->father() && exit_door != (*rFam)->mother()){
            if (exit_door->est_GV < smallest_gp) {
               smallest_gp = exit_door->est_GV;
               C = exit_door;
               cur_family = *rFam;
               perfectcut = 1;
            }
         }
      }
      if ((*rFam)->seq_id == 1) break;
      ++rFam;
   }
}
else {
   while (rFam != family_loop_connector.rend()) {
      if (((*rFam)->mother_id() == 4 && (*rFam)->father_id() == 10)
          && (*rFam)->myoffspring[0]->id() == 11) {
         C = (*rFam)->myoffspring[0];
         cur_family = *rFam;
         perfectcut = 1;
      }
      if ((*rFam)->seq_id == 1) break;
      ++rFam;
   }
}
   if (perfectcut == 0) {
      Individual *father, *mother;
      rFam = family_loop_connector.rbegin();
      while (rFam != family_loop_connector.rend()) {
         enter_door = (*rFam)->in_connector;
         exit_door = (*rFam)->out_connector;
         if ((*rFam)->numcut == 0) {
            father = (*rFam)->father();
            mother = (*rFam)->mother();
            if (!enter_door->connect_keeper) {
               if (enter_door == mother) {
                  if ((*rFam)->workvec[1] < smallest_gp) {
                     smallest_gp = (*rFam)->workvec[1];
                     C = enter_door;
                     cur_family = *rFam;
                  }
               }
               else if (enter_door == father) {
                  if ((*rFam)->workvec[0] < smallest_gp) {
                     smallest_gp = (*rFam)->workvec[0];
                     C = enter_door;
                     cur_family = *rFam;
                  }
               }
            }
            if (!exit_door->connect_keeper) {
               if (exit_door == mother) {
                  if ((*rFam)->workvec[1] < smallest_gp) {
                     smallest_gp = (*rFam)->workvec[1];
                     C = exit_door;
                     cur_family = *rFam;
                  }
               }
               else if (exit_door == father) {
                  if ((*rFam)->workvec[0] < smallest_gp) {
                     smallest_gp = (*rFam)->workvec[0];
                     C = exit_door;
                     cur_family = *rFam;
                  }
               }
            }
         }
         if ((*rFam)->seq_id == 1) break;
         ++rFam;
      }
   }
   // wait a minute, let's first update related_family in C
   C->related_family.remove(cur_family);
   cur_family->cutting(C);
stdid++;
}

int Population::peeling_sequence(void)
{
   ////////////////////////////////////////////////////////////////////////////
   // Memories for anterior and posterior are allocated here for ever, although
   //     contents in anterior and posterior will be changed from time to time
   // We also determine who will have extended ped from iterating peeling:
   //    anterior_iw = 1     if extended ped will be clipped into anterior
   //    posterior_iw[i] = 1 if extended ped will be clipped into posterior[i]
   //    anterior_iw = 2     if dangling parts were clipped into anterior
   //    posterior_iw[i] = 2 if dangling parts were clipped into posterior[i]
   //    anterior_iw = 3     extended parts, but pretend there doesn't exit.
   //    posterior_iw[i] = 3 extended parts, but pretend there doesn't exit.
   //
   // IMPORTANT: never change anterior_iw and posterior_iw[i]
   ////////////////////////////////////////////////////////////////////////////

   build_nufamily();
   build_trans_mat();
   if (!trans_mat) throw exception("Population::peeling_sequence(): fail to build trans matrix");
   unsigned i,j,noffs;
   for (i=0; i<num_nufamily; i++) {
     nufamily_vec[i]->tn_genotype = tn_genotype;
     nufamily_vec[i]->tn_qtl = tn_qtl; //BRS added
   }
   ///////////////////////////////////////////////////////////////////
   //  count # of nuclei families each individual belongs to
   //  store this counter into group_id
   ///////////////////////////////////////////////////////////////////
   Individual *I;
   for (i=0; i<popsize; i++) {
      I= popmember[i];
      I->residual_var = &residual_var;         // this is for time-being
      I->anterior.resize(tn_genotype+1);
      noffs = I->nspouse();                   // noffs = num of spouses
      if (noffs) {
         if (I->posterior) {
	   delete [] I->posterior;
           I->posterior=0;
	 }
	 if(noffs>0){
	   I->posterior = new Vector<double> [noffs];
	   I->posterior_iw.resize(noffs,0.0); //BRS
	 }
	 else {
	   I->posterior = 0;
	 }
         check_ptr(I->posterior);
         for (j=0; j<noffs; j++) I->posterior[j].resize(tn_genotype+1);
      }
      I->p_origin = 1;                       // the phenotype is original
      I->group_id = noffs;                   // noffs = num of spouses
      if (I->father() && I->mother()) I->group_id += 1;
   }
   if (!nufamily_vec) return 0;           // everyone is isolated individual

   NuFamily *family;
   for (i=0; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      family->seq_id = family->build_connectors();
   }
   //////////////////////////////////////////////////////////////////////
   // now break loops in the pedigree one by one if there is any.
   // peeling sequence will be build while we are detecting loops
   // two kinds of terminals: nonloop_terminal and inloop_terminal
   /////////////////////////////////////////////////////////////////////
// stdid = 0;
   noffs = 0;
   i = detect_loop();
   nonloop_t_nufamily = num_t_nufamily;        // 1st time detection
   while (i) {
      if (noffs == 0) maxant_maxpost();
      break_loop();  noffs++;
// std::cout << "****** loop " << noffs << " has been broken ******\n" << std::endl;
      i = detect_loop();
   }
std::cout << "the # of loops detected and then broken = " << noffs << std::endl;

   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      if (family->in_connector) family->in_connector->loop_connector = 1;
      if (family->out_connector) family->out_connector->loop_connector = 1;
   }

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      noffs = I->nspouse();
      if (I->anterior_iw == 3) I->anterior_iw = 1;
      for (j=0; j<noffs; j++) {
         if (I->posterior_iw[j] == 3) I->posterior_iw[j] = 1;
      }
   }

   ///////////////////////////////////////////////////////////////////////
   // release memory for connector->related_family,it is useless anymore,
   // however, we lose the track of connectors, so we try everyone
   ////////////////////////////////////////////////////////////////////////
   for (i=0; i<popsize; i++) popmember[i]->related_family.clear();
std::cout << " peeling sequence has been built" << std::endl;
   return 0;
}

Vector<double> Population::get_genotype_freq(void)
{
   if (!pop_gamete) build_pop_gamete();
   unsigned nchrom = prior->nchrom();
   Vector<double> genotype_freq(tn_genotype);
   Vector<double> gamete_freq(tn_gamete);
   Vector<unsigned> ngtvec(nchrom);
   Vector<unsigned> gtvec(nchrom);
   double val;
   unsigned i,j,k;
   for (i=0; i<nchrom; i++) ngtvec[i] = pop_gamete[i].nchrom();
   for (i=0; i<tn_gamete; i++) {
      gtindex(i,nchrom,ngtvec,gtvec);
      for (val=1.0,j=0; j<nchrom; j++) {
         val *= pop_gamete[j].chromosome[gtvec[j]].freq();
      }
      gamete_freq[i] = val;
   }
   for (k=0,i=0; i<tn_gamete; i++) {
      for (j=0; j<i; j++) {
         genotype_freq[k++] = 2.0*gamete_freq[i]*gamete_freq[j];
      }
      genotype_freq[k++] = gamete_freq[i]*gamete_freq[i];
   }
   return genotype_freq;
}

double Population::fullsibs_prob(Individual* sire,Individual* dam,
                                  Individual* excludeI,doubleMatrix& wspace)

{
   wspace.assign(1.0);
   unsigned j,k,um,uf,uj,nsp,nfs;
   const unsigned *noffs_sp;
   double val, scale;
   Vector<double> post_vec(tn_genotype);
   Vector<double> pen_vec(tn_genotype);
   Individual *J,**offspring, **fullsibs,**spouses;

   spouses = sire->spouses();
   offspring = sire->offspring();
   nsp = sire->nspouse();
   noffs_sp = sire->noffs_spouse();
   for (k=0, j=0; j<nsp; j++) {
      if ( spouses[j] == dam) break;
      k += noffs_sp[j];
   }
   nfs = noffs_sp[j];
   fullsibs = &offspring[k];

   for (scale = 0.0, j=0; j<nfs; j++) {
      J = fullsibs[j];
      if (J != excludeI) {
         nsp = J->nspouse();
         if (nsp ) scale += get_posterior(J,0,post_vec);
         scale += J->get_penetrance(pen_vec);
         for (um=0; um<tn_genotype; um++) {
            for (uf=0; uf<tn_genotype; uf++) {
               val = 0.0;
               if (nsp) {
                  for (uj=0; uj<tn_genotype; uj++) {
                     val += trans_mat[uj][um][uf]*pen_vec[uj]*post_vec[uj];
                  }
               }
               else {
                  for (uj=0; uj<tn_genotype; uj++) {
                     val += trans_mat[uj][um][uf]*pen_vec[uj];
                  }
               }
               wspace[um][uf] *=  val;
            }
         }
      }
   }
   //   std::cout << "scale from fullsibs " << scale << std::endl;
   return scale;
}

double Population::get_posterior(Individual* I,Individual* exclJ,Vector<double> &vec)
{
   unsigned i,j,nsp = I->nspouse();
   Individual **spouses = I->spouses();
   double retval,*ve;
   for (i=0; i<tn_genotype; i++) vec[i]=1.0;
   for (retval=0.0,i=0; i<nsp; i++) {
      if (spouses[i] != exclJ) {
         ve = I->posterior[i].begin();
         for (j=0; j<tn_genotype; j++) {
	   vec[j] *= *ve++;
	 }
         retval += *ve;
	 //	 std::cout << i << " retval " << retval << std::endl;
      }
   }
   //   std::cout << "retval from get_posterior " << retval << std::endl;
   return retval;
}

void Population::anterior(Individual* I,doubleMatrix& wspace)
{
   /////////////////////////////////////////////////////////////////
   // REQUIREMENTS:  memory must be allocated for anterior
   /////////////////////////////////////////////////////////////////
   if (I->base()) throw exception("Population::anterior(): it shoould initialize first");
   double ai,val,scale,**trme;
   unsigned um,uf,ui;
   Vector<double> post_vec(tn_genotype);
   Vector<double> pen_vec(tn_genotype);
   Vector<double> apm(tn_genotype);
   Vector<double> apf(tn_genotype);
   Individual *father, *mother;

   mother = I->mother();
   father = I->father();
   scale = fullsibs_prob(father,mother,I,wspace); //Beware of call changes here
   scale += mother->anterior[tn_genotype];
   for (um=0; um<tn_genotype; um++) apm[um] = mother->anterior[um];
   if (mother->nspouse() > 1) {
      scale += get_posterior(mother,father,post_vec);
      for (um=0; um<tn_genotype; um++) apm[um] *= post_vec[um];
   }

   scale += father->anterior[tn_genotype];
   for (uf=0; uf<tn_genotype; uf++) apf[uf] = father->anterior[uf];
   if (father->nspouse() > 1) {
      scale += get_posterior(father,mother,post_vec);
      for (uf=0; uf<tn_genotype; uf++) apf[uf] *= post_vec[uf];
   }

   scale += I->get_penetrance(pen_vec);
   for (val=0.0,ui=0; ui<tn_genotype; ui++) {
      trme = trans_mat[ui].begin();
      for (ai=0.0,um=0; um<tn_genotype; um++) for(uf=0; uf<tn_genotype; uf++) {
         ai += apm[um]*trme[um][uf]*wspace[um][uf]*apf[uf];
      }
      val += I->anterior[ui] = ai*pen_vec[ui];
   }
   //This is new
   if (val == 0.0) throw exception(std::string("Population::anterior(): genotypes conflict. Related animal:") + ind_name(I->id()) );
   for (ui=0; ui<tn_genotype; ui++) I->anterior[ui] /= val;
   scale += std::log(val);
   I->anterior[tn_genotype] = scale;
}

void Population::posterior(Individual* I,Individual* J,doubleMatrix& wspace,
                           const unsigned pj)
{
   /////////////////////////////////////////////////////////////////
   // REQUIREMENTS:  memory must be allocated for posterior
   /////////////////////////////////////////////////////////////////
   unsigned ui,uj;
   double val, scale, *ws, *post_j;
   post_j = I->posterior[pj].begin();
   ////////////////////////////////////////////////////////////////
   // if pj is not provided, then do the following
   //     nsp = I->nspouse();
   //     for (j=0; j<nsp; j++) if (I->spouselist[j]==J) break;
   //     post_j = I->posterior[j].ve;
   ////////////////////////////////////////////////////////////////

   Vector<double> post_vec(tn_genotype);
   Vector<double> apj(tn_genotype);

   scale = fullsibs_prob(I,J,0,wspace);
   for (uj=0; uj<tn_genotype; uj++) apj[uj] = J->anterior[uj];
   scale += J->anterior[tn_genotype];
   if (J->nspouse() > 1) {
      scale += get_posterior(J,I,post_vec);
      for (uj=0; uj<tn_genotype; uj++) apj[uj] *= post_vec[uj];
   }
   for (ui=0; ui<tn_genotype; ui++) {
      ws = wspace[ui];
      for (val=0.0, uj=0; uj<tn_genotype; uj++) val += *ws++ * apj[uj];
      post_j[ui] = val;
   }
   for (val=0.0, uj=0; uj<tn_genotype; uj++) val += post_j[uj];
   //This is new
   if (val == 0.0) throw exception(std::string("Population::posterior(): genotypes conflict. Related animal:") + ind_name(I->id()));
   for (uj=0; uj<tn_genotype; uj++) post_j[uj] /= val;
   scale += std::log(val);
   post_j[tn_genotype] = scale;
}

void  Population::anterior_posterior(Individual* I, doubleMatrix& wspace)
{
   //////////////////////////////////////////////////////////////////////
   //  iterative peeling over individual I
   // REQUIREMENTS:  memories must be allocated for anterior and posterior
   ///////////////////////////////////////////////////////////////////////
   if (!(I->base() ||  I->anterior_iw == 2 || I->anterior_iw == 3)) {
      anterior(I,wspace);
   }

   Individual **spouses = I->spouses();
   unsigned nsp = I->nspouse();
   for (unsigned j=0; j<nsp; j++) {
      if (!(I->posterior_iw[j] == 2 || I->posterior_iw[j] == 3)) {
         posterior(I,spouses[j],wspace,j);
      }
   }
}

void Population::partial_iterative_peeling(Individual *I, doubleMatrix& wspace)
{
   unsigned j,k,nsp,nfs;
   const unsigned *noffs_sp;
   Individual *mother, *father,  **offspring, **spouses, *J, **fullsibs;
   mother = I->mother();
   father = I->father();

   if (father) {
      if (!father->base() && father->anterior_iw != 2) anterior(father,wspace);
      spouses = father->spouses();
      offspring = father->offspring();
      nsp = father->nspouse();
      for (j=0; j<nsp; j++) {
         if (spouses[j] != mother && father->posterior_iw[j] != 2) {
            posterior(father,spouses[j],wspace,j);
         }
      }
      noffs_sp = father->noffs_spouse();
      for (k=0, j=0; j<nsp; j++) {
         if ( spouses[j] == mother) break;
         k += noffs_sp[j];
      }
      nfs = noffs_sp[j];
      fullsibs = &offspring[k];
      for (k=0; k<nfs; k++) {
         J = fullsibs[k];
         if (J != I) {
            spouses = J->spouses();
            nsp = J->nspouse();
            for (j=0; j<nsp; j++) {
               if (J->posterior_iw[j] != 2) posterior(J,spouses[j],wspace,j);
            }
         }
      }
      if (mother) {
         if (!mother->base() && mother->anterior_iw != 2) {
            anterior(mother,wspace);
         }
         spouses = mother->spouses();
         nsp = mother->nspouse();
         for (j=0; j<nsp; j++) {
            if (spouses[j] != father && mother->posterior_iw[j] !=2) {
               posterior(mother,spouses[j],wspace,j);
            }
         }
      }

      if (!I->base() && I->anterior_iw != 2) anterior(I,wspace);
      spouses = I->spouses();
      nsp = I->nspouse();
      for (j=0; j<nsp; j++) {
         if (I->posterior_iw[j] !=2) posterior(I,spouses[j],wspace,j);
      }
   }
}

double Population::log_likelihood_peeling(const unsigned maxit)
{
   ////////////////////////////////////////////////////////////////////////////
   // Memories for I->anterior and I->posterior are allocated once in
   //   peeling_sequence(), of course, peeling sequence is determined once,too
   // Since log_likelihood_peeling(0 will be called sequentially, don't need
   //   to release memories for I->anterior and I->posterior.
   ////////////////////////////////////////////////////////////////////////////
   if (!nufamily_vec ) {
      if (peeling_sequence()) throw exception("Population::log_likelihood_peeling(): "
               "fail to build the peeling sequence");
   }
   Vector<double> genotype_freq = get_genotype_freq();
   Individual *I;
   double llhood_val = 0.0;
   unsigned i,num_isolated = 0;

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      I->initial_anterior(genotype_freq.begin(),tn_genotype);
      I->initial_posterior(tn_genotype);
      if (I->isolated()) {
         num_isolated++;
         llhood_val += I->anterior[tn_genotype];
      }
   }
std::cout << "number of isolated individuals = " << num_isolated << std::endl;
   if (!nufamily_vec) return llhood_val;     // because everyone is isolated

   doubleMatrix wspace(tn_genotype,tn_genotype);
   ///////////////////////////////////////////////////////////////////
   //  first peeling for nonloop_terminal nuclear families only
   // ie. dangling data will be first clipped out, indicated by iw=1
   //////////////////////////////////////////////////////////////////
   NuFamily* family;
   for (i=0; i<nonloop_t_nufamily; i++) {
      family = nufamily_vec[i];
      if (family->kernal) {
         llhood_val += family->log_likelihood(trans_mat,wspace);
      }
      else {
         family->terminal_peeling(trans_mat,2,wspace);
      }
   }
   if (nonloop_t_nufamily == num_nufamily) return llhood_val;

   ////////////////////////////////////////////////////////////////////
   //  iterating on nuclear families in unbroken loops
   //  IMPORTANT: don't touch anterior and posterior into which
   //             dangling data were clipped
   /////////////////////////////////////////////////////////////////////
   //   std::cout << "here to do iterations!!!!!!!!" << std::endl;
   if (maxit > 0) {
      unsigned j,nsp;
      Individual **spouses;
      unsigned niterate = 0;
      while (niterate++ < maxit) {
	std::cout << "Niterate " << niterate << std::endl;
         for (i=0; i<popsize; i++) {
            I = popmember[i];
            if (I->loop_connector) {
               if (!I->base() && I->anterior_iw != 2) anterior(I,wspace);
	       //	       std::cout << "id " << I->id() << " anterior " << I->anterior << std::endl;
               spouses = I->spouses();
               nsp = I->nspouse();
               for (j=0; j<nsp; j++) {
                  if (I->posterior_iw[j] !=2) posterior(I,spouses[j],wspace,j);
               }
            }
         }
      }   // end of while loop
   }

/*
   if (maxit >0) {
      for (i=0; i<popsize; i++) {
         I = popmember[i];
         if (I->connect_keeper) {
            partial_iterative_peeling(I,wspace);
         }
      }
   }
*/
//std::cout << " nonloop_t_nufamily = " << nonloop_t_nufamily << std::endl;
//std::cout << " num_nufamily = " << num_nufamily << std::endl;
//std::cout << " computing log-likelihood for original data and the extended" << std::endl;

  ///////////////////////////////////////////////////////////////////////
  // if there are nuclear families with two offspring being cut,
  // then only one of them is allowed to have extra anterior information.
  // WATCH OUT: I am not sure this rule is necessary
  ////////////////////////////////////////////////////////////////////////

   unsigned j,noffs,nrel,id;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     //   std::cout << " Family "<< " numcut " << std::endl;
      family = nufamily_vec[i];
      if (family->numcut > 1) {
	//	std::cout << i << " ***** numcut = " << family->numcut << "***********" << std::endl;
         nrel = 1;
         noffs = family->noffs();
         for (j=0; j<noffs; j++) {
            id = family->myoffspring[j]->id();
            if (id > popsize) {
               id -= popsize;
               if (nrel < family->numcut) {
                  popmember[id-1]->initial_anterior(genotype_freq.begin(),
                                                     tn_genotype);
                  nrel++;
               }
            }
         } // end offspring for
      } //end if family->numcut
   }  //end loop i

   ////////////////////////////////////////////////////////////////////////
   //  now calculate the log-likelihood for original data and extended ped
   //  which is hidden in the duplicated individuals' anterior or posterior
   ////////////////////////////////////////////////////////////////////////
//   std::cout << " calculate the log-likelihood for original data and extended ped " << std::endl;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      if (family->kernal) {
         llhood_val += family->log_likelihood(trans_mat,wspace);
      }
      else {
         family->terminal_peeling(trans_mat,wspace);
      }
   }
//std::cout << " log likelihood for cut-extended = " << llhood_val << std::endl;
//std::cout << " removing original data" << std::endl;

   //////////////////////////////////////////////////////////////////////////
   // remove original data from individuals in the looped nuclear families
   // and then calculate the likelihood for the extended ped
   ///////////////////////////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      nufamily_vec[i]->pretend_missing(1,genotype_freq);
   }
//std::cout << "computing log-likelihood for the extended ped" << std::endl;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      if (family->kernal) {
         llhood_val -= family->log_likelihood(trans_mat,wspace);
      }
      else {
         family->terminal_peeling(trans_mat,wspace);
      }
   }

   ///////////////////////////////////////////////////////
   //  now put the original data back
   //////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      nufamily_vec[i]->pretend_missing(0,genotype_freq);
   }
   return llhood_val;
}

void Population::genotype_dist_peeling(const int prtlevel,const int estifreq)

  // void Population::genotype_dist_peeling(const unsigned prtlevel) // old call
{
   //////////////////////////////////////////////////////////////////////////
   // this works for arbitrary pedigrees
   // based on iterative peeling algorithm: Arendonk(1989), Fernando(1993)
   // prtlevel = 0, print nothing
   //            1, default, print genotype distribution and genotypic values
   // estifreq = 0, not estimate allele frequencies
   //            1, default, estimate allele frequencies
   //////////////////////////////////////////////////////////////////////////
   build_trans_mat();
   if (!trans_mat) throw exception("Population::genotype_dist_peeling(): fail to build trans matrix");
   Vector<double> post_vec(tn_genotype);
   Vector<double> genotype_freq = get_genotype_freq();
   double *gfreq = genotype_freq.begin();
   doubleMatrix wspace(tn_genotype,tn_genotype);
   unsigned i,j,k,nallele;
   doubleMatrix I_genotype_freq(popsize,tn_genotype);   // this is working space
   Individual *I;
   double *ve,val,diff,maxchange, maxchange_freq;
   unsigned inner_niterate,outer_niterate;
   unsigned nb = this->nbase();
   double allelefreq;
   outer_niterate = 0;
   do {
      genotype_freq = get_genotype_freq();
      gfreq = genotype_freq.begin();

      for (i=0; i<popsize; i++) {
         I = popmember[i];
         I->residual_var = &residual_var;       // this is for time-being
         ////////////////////////////////////////////////////////////////////
         // the last column in anterior and posterior are for scaling factor
         ////////////////////////////////////////////////////////////////////
         I->initial_anterior(gfreq,tn_genotype);
         k = I->nspouse();
         if (!I->posterior) {
	   if(k>0){
	     I->posterior = new Vector<double> [k];
	   }
	   else {
	     I->posterior = 0;
	   }
         }
         for (j=0; j<k; j++) I->posterior[j].resize(tn_genotype+1,1.0);
      }

      inner_niterate = 0;
      do {
         inner_niterate++;
         maxchange = 0.0;
         for (i=0; i<popsize; i++) {
            I = popmember[i];
            ve = I_genotype_freq[i];
            anterior_posterior(I, wspace);
            get_posterior(I,0,post_vec);
            for (val=0.0,j=0; j<tn_genotype; j++) {
               val += gfreq[j] = I->anterior[j]*post_vec[j];
            }
            if (inner_niterate > 1) {
               for (j=0; j<tn_genotype; j++) {
                  gfreq[j] /= val;
                  diff = fabs(ve[j] - gfreq[j]);
                  if (diff > maxchange) maxchange = diff;
               }
            }
            for (j=0; j<tn_genotype; j++) ve[j] = gfreq[j];
         }
      } while (inner_niterate == 1 || maxchange >= 0.000001);

      if (prtlevel==1) {
         std::cout << "   the # of inner iterations:  " << inner_niterate <<  "\n";
      }

      if (estifreq == 1) {
         //////////////////////////////////////////////////////////////////
         // it only works for the case of two-allele and one chromosome
         // PLEASE FIX ME
         ///////////////////////////////////////////////////////////////////

         if (prior->nchrom() != 1 || prior->chrom()[0].nloci() != 1
             || prior->chrom()[0].locus[0].nallele() != 2) {
            throw exception("genotype_dist_peeling(): works only for the 1-chrom and 2-allele");
         }
         outer_niterate++;
         if (prtlevel==1) {
            std::cout << " outer iteration: " << outer_niterate << "\n";
         }
         for (allelefreq=0.0,i=0; i<popsize; i++) {
            I = popmember[i];
            if (I->base()) {
               ve = I_genotype_freq[i];
               allelefreq += ve[0];
               allelefreq += ve[1] * 0.5;
            }
         }
         allelefreq /= nb;
         val = prior->chrom()[0].locus[0].allele_freq[0];
         maxchange_freq = fabs(allelefreq - val);
         prior->chrom()[0].locus[0].allele_freq[0] = allelefreq;
         prior->chrom()[0].locus[0].allele_freq[1] = 1.0-allelefreq;
	 if(pop_gamete){
	   delete [] pop_gamete;
	   pop_gamete=0;
	 }
         pop_gamete = 0;
      }
      else {
         maxchange_freq = 0.0;
      }
   } while (maxchange_freq >= 0.00001);
   if (prtlevel == 1) {
      std::cout << "allele frequencies\n";
      for (i=0; i<prior->nchrom(); i++) {
         for (j=0; j<prior->chrom()[i].nloci(); j++) {
            for (k=0; k<prior->chrom()[i].locus[j].nallele(); k++) {
               std::cout << prior->chrom()[i].locus[j].allele_freq[k] << " ";
            }
            std::cout << "\n";
         }
      }
   }

   I = popmember[popsize-1];
   double llhood_val = get_posterior(I,0,post_vec);
   val = I->anterior.inner_product(post_vec);
   llhood_val +=  std::log(val) + I->anterior[tn_genotype];
   std::cout << "log_likelihood = " << llhood_val << "\n";

   int W = SESSION.output_precision + 6;           // 6 = +.e+00
   const double **gv;

   unsigned kk,nchrom=prior->nchrom();
   for (kk=0,k=0; k<nchrom; k++) {
      gv = prior->genotypic_val(k,0);
      nallele = prior->chrom()[k].locus[0].nallele();
      for (i=0; i<nallele; i++) for (j=0; j<=i; j++) {
         gfreq[kk++] = gv[i][j];
      }
   }
   for (i=0; i<popsize; i++) {
      I = popmember[i];
      if (prtlevel==1) {
         std::cout << ind_name(I->id())
              << " " << I->sex()
              << " " << I->id()
              << " " << I->mother_id()
              << " " << I->father_id()
              << "\n";
      }
      ve = I_genotype_freq[i];
      for (val=0.0, j=0; j<tn_genotype; j++) {
         val += ve[j] * gfreq[j];
         if (prtlevel==1) std::cout << " " << std::setw(W) << ve[j];
      }
      I->est_GV = val;
      if (prtlevel==1) std::cout << "   " << std::setw(W) << val << "\n";
   }
}


//BRS
double Population::multipoint_init_parm(Fpmm& Farg)
{
  mean_for_pgenotype = Farg.getgv();
  P_freq = Farg.getgenfreq();
  P_ndim = Farg.get_dim();
  F = &Farg;

// Reading in genotypic values for QTL
   std::cout << " Reading in genotypic values for QTL \n";

  mean_for_genotype.resize(4);
  mean_for_genotype(1) = prior->chrom()[0].locus[0].genotypic_val_mat(1,1);
  mean_for_genotype(2) = prior->chrom()[0].locus[0].genotypic_val_mat(1,2);
  mean_for_genotype(3) = prior->chrom()[0].locus[0].genotypic_val_mat(2,1);
  mean_for_genotype(4) = prior->chrom()[0].locus[0].genotypic_val_mat(2,2);

  std::cout << mean_for_genotype << std::endl;

// Reading in QTL allele frequencies into Q_freq

  std::cout << "Reading in QTL allele frequencies into Q_freq \n";

   double QTL_allele_freq = prior->chrom()[0].locus[0].allele_freq[0];
   Q_freq.resize(4);
   Q_freq(1) = QTL_allele_freq*QTL_allele_freq;
   Q_freq(2) = QTL_allele_freq*(1-QTL_allele_freq);
   Q_freq(3) = QTL_allele_freq*(1-QTL_allele_freq);
   Q_freq(4) = (1-QTL_allele_freq)*(1-QTL_allele_freq);

   std::cout << Q_freq;
   marker_freq.resize(n_markerLoci);


// Reading in frequencies for marker i into vector marker_freq(i)

   for (unsigned i=1;i<=n_markerLoci;i++) {
     marker_freq(i) = prior->chrom()[0].locus[i].allele_freq;
   }
return 0;
}


double Population::multipoint_likelihood(int maxit) {
   ////////////////////////////////////////////////////////////////////////////
   // Memories for I->anterior and I->posterior are allocated once in
   //   peeling_sequence(), of course, peeling sequence is determined once,too
   // Since log_likelihood_peeling(0 will be called sequentially, don't need
   //   to release memories for I->anterior and I->posterior.
   ////////////////////////////////////////////////////////////////////////////
 // Vector<double> genotype_freq;

  Individual *I;
  double llhood_val = 0.0, temp;
  unsigned num_isolated = 0,i,j;
  unsigned max_switches = (unsigned) (1 << n_markerLoci); //pow(2,n_markerLoci);
  NuFamily::penetrance.resize(P_ndim,tn_qtl);
  NuFamily::wspace.resize(P_ndim*tn_qtl*max_switches,P_ndim*tn_qtl*max_switches);
  doubleMatrix pen(P_ndim,tn_qtl);  // One should not need this!
  NuFamily* family;
  // std::cout << "entered Population::multipoint_likelihood \n";
  // std::cout << "Mean for geno 0: " << mean_for_genotype[0] << endl;
  // std::cout << "Mean for geno 1: " << mean_for_genotype[1] << endl;
  // std::cout << "Mean for geno 2: " << mean_for_genotype[2] << endl;
  // std::cout << "Mean for geno 3: " << mean_for_genotype[3] << endl;

   if (!nufamily_vec) {
      if (peeling_sequence()) throw exception("Population::log_likelihood_peeling(): "
               "fail to build the peeling sequence");
   }

   // peeling_sequence(); Need to reset flags as well!
   // Need to initialize everyone!


   for (i=0; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     family->multi_initialize(pen);
   }

 for (i=0; i<popsize; i++) {
     I = popmember[i];
     //   std::cout << ind_name(I->id()) << " A scale is " << I->m_anterior_scale << std::endl;
     if (I->isolated()) {
       num_isolated++;
       I->initial_multi_anterior(pen);// these do not appear in nufamily_vec and are not initialized
       I->initial_multi_posterior(-1);
       llhood_val += I->m_anterior_scale;
     }
   }
   // std::cout << "number of isolated individuals = " << num_isolated;
   // std::cout  << " with new llh = " << llhood_val << std::endl;
   if (!nufamily_vec) return llhood_val;     // because everyone is isolated


   Dblock post_mat_m, post_mat_f;
   post_mat_m.resize(P_ndim,max_switches,tn_qtl,1.0);
   post_mat_f.resize(P_ndim,max_switches,tn_qtl,1.0);

   ///////////////////////////////////////////////////////////////////
   //  first peeling for nonloop_terminal nuclear families only
   // ie. dangling data will be first clipped out and iw set to 2
   //////////////////////////////////////////////////////////////////


   for (i=0; i<nonloop_t_nufamily; i++) {
   family = nufamily_vec[i];
   //     family->display();
       // family = nufamily_vec[i];
//      family->display(); // prints out peeling sequence
      if (family->kernal) {
	llhood_val += family->multi_llh(post_mat_m);
      }
      else {
	 family->multi_terminal_peeling(post_mat_m, post_mat_f, 2);
      }
      //  std::cout << " family " << i << " llh " << llhood_val << std::endl;
   }

   if (nonloop_t_nufamily == num_nufamily) return llhood_val;

   //   std::cout << "******** doing iterations@@@@@!!!!!!!! " << std::endl;
   ////////////////////////////////////////////////////////////////////
   //  iterating on nuclear families in unbroken loops
   //  IMPORTANT: don't touch anterior and posterior into which
   //             dangling data were clipped
   /////////////////////////////////////////////////////////////////////
   if (maxit > 0) {
      unsigned niterate = 0;
      while (niterate++ < maxit) {
	std::cout << "Iteration: " << niterate << std::endl;
	for (i=0; i<num_nufamily; i++) {
          // compute anteriors and posteriors that are not set to iw=2 or iw=3
	  nufamily_vec[i]->multi_ant_post(post_mat_m, post_mat_f);
	}
      }   // end of while loop
   }


//std::cout << " nonloop_t_nufamily = " << nonloop_t_nufamily << std::endl;
//std::cout << " num_nufamily = " << num_nufamily << std::endl;
//std::cout << " computing log-likelihood for original data and the extended" << std::endl;
  ///////////////////////////////////////////////////////////////////////
  // if there are nuclear families with two offspring being cut,
  // then only one of them is allowed to have extra anterior information.
  // WATCH OUT: I am not sure this rule is necessary
  ////////////////////////////////////////////////////////////////////////

   unsigned noffs,nrel,id;

   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     //   family->display();
     if (family->numcut > 1) {
       //     std::cout << "****** numcut = " << family->numcut << "***********\n";
       nrel = 1;
       noffs = family->noffs();
       for (j=0; j<noffs; j++) {
	 id = family->myoffspring[j]->id();
	 if (id > popsize) {
	   id -= popsize;
	   if (nrel < family->numcut) {
	     popmember[id-1]->initial_multi_anterior(pen);
	     //initial_anterior(genotype_freq.ve,tn_genotype); old initial
	     nrel++;
	   }
	 }
       }
     }
   }
   ////////////////////////////////////////////////////////////////////////
   //  now calculate the log-likelihood for original data and extended ped
   //  which is hidden in the duplicated individuals' anterior or posterior
   ////////////////////////////////////////////////////////////////////////
// std::cout << "Likelihood before loops " << llhood_val << std::endl;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];

      if (family->kernal) {
	 llhood_val += family->multi_llh(post_mat_m);
      }
      else {
	 family->multi_terminal_peeling(post_mat_m, post_mat_f, 2);
      }
   }
//std::cout << " log likelihood for cut-extended = " << llhood_val << std::endl;
//std::cout << " removing original data" << std::endl;
temp=llhood_val;
//t2=0.0;
   //////////////////////////////////////////////////////////////////////////
   // remove original data from individuals in the looped nuclear families
   // and then calculate the likelihood for the extended ped
   ///////////////////////////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      nufamily_vec[i]->pretend_multi_missing(1);
   }
//std::cout << "computing log-likelihood for the extended ped" << std::endl;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      if (family->kernal) {
	llhood_val -= family->multi_llh(post_mat_m);
      }
      else {
      family->multi_terminal_peeling(post_mat_m, post_mat_f, 2);
      }
   }

   ///////////////////////////////////////////////////////
   //  now put the original data back
   //////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      nufamily_vec[i]->pretend_multi_missing(0);
   }
   std::cout << "cut llh " << temp << " removed " << (temp-llhood_val) << " final " << llhood_val <<std::endl;
   return llhood_val;
}

void Population::multi_geno_dist_peeling(const unsigned prtlevel)
{
   //////////////////////////////////////////////////////////////////////////
   // this works for arbitrary pedigrees
   // based on iterative peeling algorithm: Arendonk(1989), Fernando(1993)
   // prtlevel = 0, print nothing
   //            1, default, print genotype distribution and genotypic values
   //////////////////////////////////////////////////////////////////////////

  // Multipoint_setup must be called before coming here!
   unsigned max_switches = (unsigned) (1 << n_markerLoci); //pow(2,n_markerLoci);
   doubleMatrix wspace(P_ndim*tn_qtl*max_switches,P_ndim*tn_qtl*max_switches);
   Dblock post_mat_m(P_ndim,max_switches,tn_qtl), post_mat_f(P_ndim,max_switches,tn_qtl);
   post_mat_m.init(1.0);
   post_mat_f.init(1.0);
   Vector<double> gfreq(tn_qtl), gprobs_old(tn_qtl),gprobs_new(tn_qtl);
   doubleMatrix pen(P_ndim,tn_qtl);
   unsigned i,j,jj;
   doubleMatrix I_genotype_freq(popsize,tn_qtl);   // to store conditional probs for 4 QTL genot.
   Individual *I;
   for (i=0; i<popsize; i++) {
      I = popmember[i];
      I->residual_var = &residual_var;           // this is for time-being
      ////////////////////////////////////////////////////////////////////////
      // the last column in both anterior and posterior are for scaling factor
      ////////////////////////////////////////////////////////////////////////
      I->initial_multi_anterior(pen);
      I->initial_multi_posterior(-1);
   }
   unsigned niterate = 0;
   double val,diff,maxchange = 1.0;
   while (maxchange >= 0.000001) {
     std::cout << "niterate = " << niterate << " " << maxchange << std::endl;
      if (niterate > 1) maxchange = 0.0;
      niterate++;
      if (niterate > 1) maxchange = 0.0;
      for (i=0; i<num_nufamily; i++) {
         nufamily_vec[i]->multi_ant_post(post_mat_m, post_mat_f);
      }
      for (i=0; i<popsize; i++) {
        gprobs_old = popmember[i]->get_gprobs();
	popmember[i]->cal_gprobs(post_mat_f);
        gprobs_new = popmember[i]->get_gprobs();
        if (niterate > 1) {
            for (jj=0; jj<tn_qtl; jj++) {
               diff = fabs(gprobs_new[jj] - gprobs_old[jj]);
               if (diff > maxchange) maxchange = diff;
            }
         }
      }
   }

   if (prtlevel==1) {
      std::cout << " total number of iterations = " << niterate-1 <<  "\n";
   }
   int W = SESSION.output_precision + 6;           // 6 = +.e+00

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      if (prtlevel==1) std::cout << ind_name(i+1) << ":\n";
      gprobs_new = popmember[i]->get_gprobs();
      for (val=0.0, j=0; j<tn_qtl; j++) {
         val += gprobs_new[j] * mean_for_genotype[j];
         if (prtlevel==1) std::cout << " " << std::setw(W) << gprobs_new[j];
      }
      I->est_GV = val;
      if (prtlevel==1) std::cout << "   " << std::setw(W) << val << "\n";
   }
   if (prtlevel==1) std::cout << std::endl;
}


void Population::maxant_maxpost(void)
{

   NuFamily *family,Ftmp;
   Individual *I, *father, *mother,Itmp;
   int i,j,k;
   int max_switches = (int) (1 << n_markerLoci); //pow(2,n_markerLoci);

   doubleMatrix pen(P_ndim,tn_qtl);


   //  SDK
   // I->g_weight.resize(tn_qtl,max_switches);
   Itmp.g_weight.resize(tn_qtl,max_switches);
   // SDK

   Vector<double> gfreq(tn_genotype);

   if (analysis_type == "multipoint") {
     gfreq[0]=Q_freq[0];
     gfreq[1]=Q_freq[1]+Q_freq[2];
     gfreq[2]=Q_freq[3];
     for (i=0; i<popsize; i++) {
       I = popmember[i];
       I->initial_multi_anterior(pen);
       I->initial_multi_posterior(-1);
     }
     multi_geno_dist_peeling(0);
     for (i=0; i<popsize; i++) {
       I = popmember[i];
       I->collapse_antpost();
     }
   }
   else if (analysis_type == "multipoint_m") {

     // SDK  
     //family->multi_wksp_resize(tn_qtl, max_switches);
     Ftmp.multi_wksp_resize(tn_qtl, max_switches);
     //SDK
     gfreq[0]=Q_freq[0];
     gfreq[1]=Q_freq[1]+Q_freq[2];
     gfreq[2]=Q_freq[3];
     for (i=0; i<popsize; i++) {
       I = popmember[i];
       I->initial_multi_m_anterior(tn_qtl);
       I->initial_multi_m_posterior(tn_qtl,-1);
     }
     multi_m_geno_dist_peeling(0);
     for (i=0; i<popsize; i++) {
       I = popmember[i];
       I->collapse_mix_antpost();
     }
   }
   else {
     gfreq=get_genotype_freq();
     for (i=0; i<popsize; i++) {
       I = popmember[i];
       I->initial_anterior(gfreq.begin(),tn_genotype);
       I->initial_posterior(tn_genotype);
     }
     genotype_dist_peeling(0);
   }

   double Pryi,Prya,Pryp,val1,val2,val3;
   Vector<double> pen_vec(tn_genotype);
   Vector<double> pos_vec(tn_genotype);
   doubleMatrix big_pen(P_ndim,tn_qtl);
   for (k=0; k < tn_genotype; k++){
     pen_vec[k]=0.0;
     pos_vec[k]=0.0;
   }

   for (i=0; i<popsize; i++) {
     I = popmember[i];
     if (I->group_id > 1) {    // I is a connector
       get_posterior(I,0,pos_vec);
       if (analysis_type == "multipoint") {
	 I->get_m_penetrance(big_pen);
	 for (k=0; k<P_ndim; k++){
	   pen_vec[0] += big_pen[k][0];
	   pen_vec[1] += big_pen[k][1] + big_pen[k][2];
	   pen_vec[2] += big_pen[k][3];
	 }
       }
       else if (analysis_type == "multipoint_m") {
 	 I->get_mix_penetrance(pen_vec);
       }
       else {
	 I->get_penetrance(pen_vec);
       }
       for (Prya=0.0,Pryi=0.0,j=0; j<tn_genotype; j++) {
	 if (pen_vec[j] > 0.0) {
               Prya += I->anterior[j] /= pen_vec[j];
               Pryi  += pen_vec[j] *= gfreq[j];
            }
            else {
               I->anterior[j] = 0.0;
            }
         }
         for (j=0; j<tn_genotype; j++) {
             I->anterior[j] /= Prya;  pen_vec[j] /= Pryi;
         }
// Pry = Pr(y); Prya = Pr(ya); Pryi = Pr(yi); Pryp = Pr(yp)
// I->anterior[j] = Pr(ui|ya)
// pen_vec[j] = Pr(ui|yi)
// gfreq[j]  = Pr(ui|y)

         val1 = I->anterior.max();
         val3 = pen_vec.max();
         I->est_GV = val1/val3;
	 //	 std::cout << ind_name(I->id()) << ": " << I->est_GV << "\n";
	 /* BRS trying change peeling sequence
	 int tag=I->id();
	 // std::cout << "here " << tag << std::endl;
	 if (tag == 3)
	     I->est_GV=3;
	 else if (tag == 8)
	   I->est_GV=3;
	 else if (tag == 5)
	   I->est_GV=3;
	 std::cout << ind_name(tag) << ": " << I->est_GV << "\n";
	 BRS */
     }
   }
   for (i=num_t_nufamily; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     family->in_connector = 0;           // re-initialize
     family->out_connector = 0;
     father = family->father();
     mother = family->mother();
     family->workvec[0]=-1.0;
     family->workvec[1]=-1.0;
     if (father->group_id > 1) {
       if (analysis_type == "multipoint") {
	 father->get_m_penetrance(big_pen);
	 for (k=0; k<P_ndim; k++){
	   pen_vec[0] = big_pen[k][0];
	   pen_vec[1] = big_pen[k][1] + big_pen[k][2];
	   pen_vec[2] = big_pen[k][3];
	 }
       }
       else if (analysis_type == "multipoint_m") {
	 father->get_mix_penetrance(pen_vec);
       }
       else {
	 father->get_penetrance(pen_vec);
       }
       get_posterior(father,mother, pos_vec);
       for (Pryi=0.0,Pryp=0.0,j=0; j<tn_genotype; j++) {
	 Pryi  += pen_vec[j] *= gfreq[j];
	 Pryp += pos_vec[j] *= gfreq[j];
       }
       for (j=0; j<tn_genotype; j++) {
	 pos_vec[j] /= Pryp; pen_vec[j] /= Pryi;
       }
       val2 = pos_vec.max();
       val3 = pen_vec.max();
       family->workvec[0] = val2/val3;
     }
     if (mother->group_id > 1) {
       if (analysis_type == "multipoint") {
	 mother->get_m_penetrance(big_pen);
	 for (k=0; k<P_ndim; k++){
	   pen_vec[0] = big_pen[k][0];
	   pen_vec[1] = big_pen[k][1] + big_pen[k][2];
	   pen_vec[2] = big_pen[k][3];
	 }
       }
       else if (analysis_type == "multipoint_m") {
	 mother->get_mix_penetrance(pen_vec);
       }
       else {
	 mother->get_penetrance(pen_vec);
       }
       get_posterior(mother,father, pos_vec);
       for (Pryi=0.0,Pryp=0.0,j=0; j<tn_genotype; j++) {
	 Pryi  += pen_vec[j] *= gfreq[j];
	 Pryp += pos_vec[j] *= gfreq[j];
       }
       for (j=0; j<tn_genotype; j++) {
	 pos_vec[j] /= Pryp; pen_vec[j] /= Pryi;
       }
       val2 = pos_vec.max();
       val3 = pen_vec.max();
       family->workvec[1] = val2/val3;
     }
     //     std::cout << father->id() << " " << mother->id() << ": "
     //	  << family->workvec[0] << " " << family->workvec[1] << "\n";
   }

   std::list<Individual *>::iterator connector;
   ///////////////////////////////////////////////////////////
   // we will build related_family list for each connectors
   //////////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
      family = nufamily_vec[i];
      for (connector = family->connectors.begin();
           connector != family->connectors.end();
           ++connector) {
          (*connector)->related_family.push_back(family);
      }
   }
}

double Population::multi_m_log_likelihood_peeling(int maxit ) {
  ////////////////////////////////////////////////////////////////////////////
  // Memories for I->anterior and I->posterior are allocated once in
  //   peeling_sequence(), of course, peeling sequence is determined once,too
  // Since log_likelihood_peeling(0 will be called sequentially, don't need
  //   to release memories for I->anterior and I->posterior.
  ////////////////////////////////////////////////////////////////////////////

//// Approximation of mixed inheritance.

//// genotype is defined by QTL and Markers
// std::cout << "here  in mixture stuff " << std::endl;
  if (!nufamily_vec) {
    if (peeling_sequence()) throw exception("Population::log_likelihood_peeling(): "
	    "fail to build the peeling sequence");
  }

   Individual *I;
   double llhood_val = 0.0,y,v,nu,temp, freq_mark;
   double P_var = F->var; // polygenic varaince
   unsigned num_isolated = 0,i, i_q, a1, a2,k;
   int max_switches = (int) (1 << n_markerLoci); //pow(2,n_markerLoci);
   // Need to initialize everyone!
   NuFamily* family;

   for (i=0; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     family->multi_m_initialize(tn_qtl);
   }
   for (i=0; i<popsize; i++) {
     I = popmember[i];
     if (I->isolated()) {
       num_isolated++;
      // calculate likelhood for this individual
       if (!(I->record()[0].missing)) {
 	freq_mark=1.0;
 	for (k=1; k <= n_markerLoci; k++){
 	  a1 = I->genome0.chromosome[0].locus[k].allele;
 	  a2 = I->genome1.chromosome[0].locus[k].allele;
 	  freq_mark *= marker_freq(k)(a1)*marker_freq(k)(a2);
 	}
 	temp = 0.0;
 	for (i_q=0; i_q<tn_qtl;i_q++){
 	  y     = I->record()[0].double_val();
 	  v     = (P_var+(residual_var[0][0]));
 	  nu    = y - mean_for_genotype[i_q];
	  temp += std::exp(-0.5*nu*nu/v)/std::sqrt(4*std::asin(1.0)*v)*Q_freq[i_q]*freq_mark;
 	}
 	llhood_val += std::log(temp);
       }
       else {
 	freq_mark=1.0;
 	for (k=1; k <= n_markerLoci; k++){
 	  a1 = I->genome0.chromosome[0].locus[k].allele;
 	  a2 = I->genome1.chromosome[0].locus[k].allele;
 	  freq_mark *= marker_freq(k)(a1)*marker_freq(k)(a2);
 	}
 	llhood_val += std::log(freq_mark);
       }
     }
   }
//  std::cout << "number of isolated individuals = " << num_isolated << std::endl;
   if (!nufamily_vec) return llhood_val;     // because everyone is isolated

   ///////////////////////////////////////////////////////////////////
   //  first peeling for nonloop_terminal nuclear families only
   // ie. dangling data will be first clipped out, indicated by iw=1
   //////////////////////////////////////////////////////////////////

   family->multi_wksp_resize(tn_qtl, max_switches);

   for (i=0; i<nonloop_t_nufamily; i++) {
     family = nufamily_vec[i];
     if (family->kernal) {
       llhood_val += family->multi_m_log_likelihood();
     }
     else {
       family->multi_m_terminal_peeling(2);
     }
   }
   if (nonloop_t_nufamily == num_nufamily) {
     return llhood_val;
   }

   // Incorporate loop cutting
   ////////////////////////////////////////////////////////////////////
   //  iterating on nuclear families in unbroken loops
   //  IMPORTANT: don't touch anterior and posterior into which
   //             dangling data were clipped
   /////////////////////////////////////////////////////////////////////
   if (maxit > 0) {
     Individual **spouses;
     unsigned niterate = 0;
     while (niterate++ < maxit) {
       std::cout << "Iteration: " << niterate << std::endl;
       for (i=0; i<num_nufamily; i++) {
 	// compute anteriors and posteriors that are not set to iw=2 or iw=3
 	nufamily_vec[i]->multi_m_ant_post();
       }
     }   // end of while loop
   }


//   std::cout << " nonloop_t_nufamily = " << nonloop_t_nufamily << std::endl;
//   std::cout << " num_nufamily = " << num_nufamily << std::endl;
//   std::cout << " computing log-likelihood for original data and the extended" << std::endl;

   ///////////////////////////////////////////////////////////////////////
   // if there are nuclear families with two offspring being cut,
   // then only one of them is allowed to have extra anterior information.
   // WATCH OUT: I am not sure this rule is necessary
   ////////////////////////////////////////////////////////////////////////

   unsigned noffs,nrel,id,j;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     if (family->numcut > 1) {
       //   std::cout << "****** numcut = " << family->numcut << "***********\n";
       nrel = 1;
       noffs = family->noffs();
       for (j=0; j<noffs; j++) {
 	id = family->myoffspring[j]->id();
 	if (id > popsize) {
 	  id -= popsize;
 	  if (nrel < family->numcut) {
 	    popmember[id-1]->initial_multi_m_anterior(tn_qtl);
 	    nrel++;
 	  }
 	}
       }
     }
   }
   ////////////////////////////////////////////////////////////////////////
   //  now calculate the log-likelihood for original data and extended ped
   //  which is hidden in the duplicated individuals' anterior or posterior
   ////////////////////////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     if (family->kernal) {
       llhood_val += family->multi_m_log_likelihood();
     }
     else {
       family->multi_m_terminal_peeling(2);
     }
   }
//   std::cout << " log likelihood for cut-extended = " << llhood_val << std::endl;
//   std::cout << " removing original data" << std::endl;
   temp=llhood_val;

    //////////////////////////////////////////////////////////////////////////
    // remove original data from individuals in the looped nuclear families
    // and then calculate the likelihood for the extended ped
    ///////////////////////////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     nufamily_vec[i]->pretend_multi_m_missing(1);
   }
//   std::cout << "computing log-likelihood for the extended ped" << std::endl;
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     family = nufamily_vec[i];
     if (family->kernal) {
       llhood_val -= family->multi_m_log_likelihood();
     }
     else {
       family->multi_m_terminal_peeling(2);
     }
   }

   ///////////////////////////////////////////////////////
   //  now put the original data back
   //////////////////////////////////////////////////////
   for (i=nonloop_t_nufamily; i<num_nufamily; i++) {
     nufamily_vec[i]->pretend_multi_m_missing(0);
   }
   std::cout << "cut llh " << temp << " removed " << (temp-llhood_val) << " final " << llhood_val << std::endl;
  return llhood_val;
}

void Population::multi_m_geno_dist_peeling(const unsigned prtlevel) {
   //////////////////////////////////////////////////////////////////////////
   // this works for arbitrary pedigrees
   // based on iterative peeling algorithm: Arendonk(1989), Fernando(1993)
   // prtlevel = 0, print nothing
   //            1, default, print genotype distribution and genotypic values
   //////////////////////////////////////////////////////////////////////////

  // Multipoint_setup must be called before coming here!
   Vector<double> gfreq(tn_qtl), gprobs_old(tn_qtl),gprobs_new(tn_qtl);
//   int family; BRS not used
//   NuFamily *NF; BRS not used

   unsigned i,j,jj;
   doubleMatrix I_genotype_freq(popsize,tn_qtl);   // to store conditional probs for 4 QTL genot.
   Individual *I;
   for (i=0; i<popsize; i++) {
      I = popmember[i];
      I->residual_var = &residual_var;           // this is for time-being
      ////////////////////////////////////////////////////////////////////////
      // the last column in both anterior and posterior are for scaling factor
      ////////////////////////////////////////////////////////////////////////
      I->initial_multi_m_posterior(tn_qtl,-1);
      I->initial_multi_m_anterior(tn_qtl); 
   }
   unsigned niterate = 0;
   double val,diff,maxchange = 1.0;
   while (maxchange >= 0.000001) {
     std::cout << "niterate = " << niterate << " " << maxchange << std::endl;
      if (niterate > 1) maxchange = 0.0;
      niterate++;
      if (niterate > 1) maxchange = 0.0;
      for (i=0; i<num_nufamily; i++) {
         nufamily_vec[i]->multi_m_ant_post();
      }
      for (i=0; i<popsize; i++) {
        gprobs_old = popmember[i]->get_gprobs();
	popmember[i]->cal_m_gprobs(tn_qtl);
        gprobs_new = popmember[i]->get_gprobs();
        if (niterate > 1) {
            for (jj=0; jj<tn_qtl; jj++) {
               diff = fabs(gprobs_new[jj] - gprobs_old[jj]);
               if (diff > maxchange) maxchange = diff;
            }
         }
      }
   }


   if (prtlevel==1) {
      std::cout << " total number of iterations = " << niterate-1 <<  "\n";
   }
   int W = SESSION.output_precision + 6;           // 6 = +.e+00

   for (i=0; i<popsize; i++) {
      I = popmember[i];
      if (prtlevel==1) std::cout << ind_name(i+1) << ":\n";
      gprobs_new = popmember[i]->get_gprobs();
      for (val=0.0, j=0; j< tn_qtl; j++) {
         val += gprobs_new[j] * mean_for_genotype[j];
         if (prtlevel==1) std::cout << " " << std::setw(W) << gprobs_new[j];
      }
      I->est_GV = val;
      if (prtlevel==1) std::cout << "   " << std::setw(W) << val << "\n";  
   }
   if (prtlevel==1) std::cout << std::endl;
}
/* BRS Superceeded by Model::MapF
double Population::MapR(double dist) {
  return (0.5*(1.0 - std::exp(-2.0*dist)));
}
*/
//BRS
} /////////// end of namespace matvec

