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
#include "population.h"

namespace matvec {

void Population::remove_mark(void)
{
   for (unsigned i=0; i<popsize; i++) popmember[i]->group_id = 0;
   mark_need_remove = 0;
}

void Population::mark_ancestor_recursive(const Individual* ind,
                                         const unsigned counter)
{
   ///////////////////////////////////////////////////////////////
   // set group_id variable of all ancestors of individual id
   // usr is responsible for resetting group_id for each individual
   // before calling
   /////////////////////////////////////////////////////////////////////////
   if (ind->mymother)  {
      if (ind->mymother->group_id == 0) {
         ind->mymother->group_id = counter;
         mark_ancestor_recursive(ind->mymother,counter);
      }
   }

   if (ind->myfather)  {
      if (ind->myfather->group_id == 0) {
         ind->myfather->group_id = counter;
         mark_ancestor_recursive(ind->myfather,counter);
      }
   }
}

unsigned Population::mark_ancestor_of(const unsigned id,const unsigned counter)
        // count all ancestors of individual id
{
   if (!stdid) renum();
   if (id<1 || id>popsize || counter == 0) throw exception(" Population::mark_ancestor_of(): bad args");
   if (mark_need_remove) remove_mark();
   mark_ancestor_recursive(popmember[id-1],counter);
   unsigned mark_index = 0;
   for (unsigned i=0; i<popsize;i++) {
      if (popmember[i]->group_id == counter) mark_index++ ;
   }
   mark_need_remove = 1;
   return mark_index;
}

Population* Population::ancestor_of(const unsigned id)
{
   //////////////////////////////////////////////////////////////
   // put all ancestors of individual id into a new population
   // user is resposible to release this new population
   //       using either delete or resize(0,0,0);
   //////////////////////////////////////////////////////////////

   if (!stdid) renum();
   if (id<1 || id >popsize) throw exception(" Population::mark_ancestor_of(): bad arg");
   unsigned counter = 1;
   unsigned nd = mark_ancestor_of(id,counter);
   return sub(nd);
}

void Population::mark_descendant_recursive(const Individual* ind,
                                           const unsigned counter)
{
   ////////////////////////////////////////////////////
   // count all descendants of individual id
   // usr is responsible for resetting group_id variable for each individual
   ///////////////////////////////////////////////////////////////////////

   unsigned n = ind->noffs();
   for (unsigned i=0; i<n;i++) {
      if (ind->myoffspring[i]->group_id == 0) {
         ind->myoffspring[i]->group_id = counter;
         mark_descendant_recursive(ind->myoffspring[i],counter);
      }
   }
}

unsigned Population::mark_descendant_of(const unsigned id,const unsigned
                                         counter)
        // count all descendants of individual id
{
   if (!stdid) renum();
   if (id<1 || id >popsize || counter == 0) throw exception("Population::mark_descendant_of(): bad args");
   if (mark_need_remove) remove_mark();
   mark_descendant_recursive(popmember[id-1],counter);
   unsigned mark_index = 0;
   for (unsigned i=0; i<popsize;i++) {
       if (popmember[i]->group_id == counter) mark_index++ ;
   }
   mark_need_remove = 1;
   return mark_index;
}
/**
 * put all descendants of individual id into a new population.
 * caller is resposible to release this new population
 * using either delete or resize(0,0,0)
 */
Population* Population::descendant_of(const unsigned id)
{
   if (!stdid) renum();
   if (id < 1 || id > popsize ) throw exception(" Population::mark_descendant_of(): bad arg");
   unsigned counter = 1;
   unsigned nd = mark_descendant_of(id,counter);
   return sub(nd);
}

/**
 * mark all relatives of a specific individual ind (including ind).
 * you have to take care the following in caller program
 *       (1) mark_need_remove variabe
 *       (2) stdid
 * a more user-end function is relative_of(id,counter)
 */
void Population::mark_relative_recursive(Individual* ind,const unsigned counter)
{
   if (ind->isolated()) {
      ind->group_id = counter;
      return;
   }

   if (ind->mymother) {
      if (ind->mymother->group_id == 0) {
         ind->mymother->group_id = counter;
         mark_relative_recursive(ind->mymother,counter);
      }
   }
   if (ind->myfather) {
      if (ind->myfather->group_id == 0) {
         ind->myfather->group_id = counter;
         mark_relative_recursive(ind->myfather,counter);
      }
   }

   unsigned n = ind->noffs();
   for (unsigned i=0; i<n; i++) {
      if (ind->myoffspring[i]->group_id == 0) {
         ind->myoffspring[i]->group_id = counter;
         mark_relative_recursive(ind->myoffspring[i],counter);
      }
   }
}

unsigned Population::mark_relative_of(const unsigned id,const unsigned counter)
{
   ////////////////////////////////////////////////
   // count all relatives of individual id
   ////////////////////////////////////////////////
   if (!stdid) renum();
   if (id<1 || id >popsize || counter == 0) throw exception(" Population::mark_relative_of(): bad args");
   if (mark_need_remove) remove_mark();
   popmember[id-1]->group_id = counter;
   mark_relative_recursive(popmember[id-1],counter);
   unsigned mark_index = 0;
   for (unsigned i=0; i<popsize;i++) {
      if (popmember[i]->group_id == counter) mark_index++ ;
   }
   mark_need_remove = 1;
   return mark_index;
}

/**
 * put all relatives of individual id into a new population.
 * user is responsible for releasing this new population
 *       using either delete or resize(0,0,0);
 */
Population* Population::relative_of(const unsigned id)
{
   if (!stdid) renum();
   unsigned counter = 1;
   unsigned nr = mark_relative_of(id,counter);
   return sub(nr);
}

/**
 * mark all families in the population and return the number of accounted families.
 */
unsigned Population::mark_families(const unsigned start_family_id)
{
   if (!stdid) renum();
   if (popsize == 0) return 0;
   if (mark_need_remove) remove_mark();
   unsigned counter = start_family_id;
   for (unsigned i=0; i<popsize; i++) {
      if (popmember[i]->group_id == 0) {
         mark_relative_recursive(popmember[i],counter++);
      }
   }
   mark_need_remove = 1;
   return (counter - start_family_id);
}

Population* Population::ancestor_of(const char name[])
{
   unsigned id = get_id(name);
   if (id == 0) throw exception("Population::ancestor_of(): cannot found");
   return ancestor_of(id);
}

Population* Population::descendant_of(const char name[])
{
   unsigned id = get_id(name);
   if (id == 0) throw exception(" Population::descendant_of(): cannot found");
   return descendant_of(id);
}

Population* Population::relative_of(const char name[])
{
   unsigned id = get_id(name);
   if (id == 0) throw exception("Population::relative_of(): cannot found");
   return relative_of(id);
}

Population* Population::family_of(const char name[])
{
   unsigned id = get_id(name);
   if (id == 0) throw exception("Population::family_of(): cannot found");
   return relative_of(id);
}

unsigned Population::fetch_families(const char filename[])
{
   std::ofstream outfile(filename);
   if (!outfile) throw exception("Population::fetch_families(): cannot open or already exists");
   unsigned retval = this->fetch_families(outfile);
   outfile.close();
   return retval;
}

unsigned Population::fetch_families(std::ostream& stream)
{
  if (!stdid) renum();
  if (popsize == 0) return 0;
  if (mark_need_remove) remove_mark();

   Individual *I;
   int k;
   unsigned i,familyid = 1;
   for (i=0; i<popsize; i++) {
      I = popmember[i];
      if (I->group_id == 0) {
         mark_relative_recursive(I,familyid++);
      }
   }
   familyid--;

   unsigned nf;
   unsigned fd = 1;
   while (fd <= familyid) {
      for (nf=0,i=0; i<popsize; i++) if (popmember[i]->group_id == fd) nf++;
      stream << "family " << fd << ": " << nf << "\n";
      for (i=0; i<popsize; i++) {
         I = popmember[i];
         if (I->group_id != fd) continue;
         if (maxnamelen == 0) {
            stream << std::setw(10) << I->id();
            if (I->mother()) {
               stream << " " << std::setw(10) << I->mother()->id();
            }
            else {
               stream << " " << std::setw(10) << ".";   // TW
            }
            if (I->father()) {
               stream << " " << std::setw(10) << I->father()->id();
            }
            else {
               stream << " " << std::setw(10) << ".";   // TW
            }
         }
         else {
            stream << std::setw(maxnamelen) << ind_name(I->id());
            if (I->mother()) {
               stream << " " << std::setw(maxnamelen)<<ind_name(I->mother()->id());
            }
            else {
               stream << " " << std::setw(maxnamelen) << ".";   // TW
            }
            if (I->father()) {
                stream << " " << std::setw(maxnamelen)<<ind_name(I->father()->id());
            }
            else {
               stream << " " << std::setw(maxnamelen) << ".";   // TW
            }
         }
         for (k=0; k<numtrait; k++) {
            if (I->record()[k].missing) {
               stream << " .";   // TW
            }
            else {
               stream << " " << I->record()[k].double_val();
            }
         }
         stream << "\n";
      }
      fd++;
   }
   mark_need_remove = 1;
   return familyid;
}

void Population::inbcoef_quaas(void)
{
   if (!stdid) renum();
   unsigned i,j,s,d;
   Individual *I;
   double t;
   Vector<double> fvalue; fvalue.reserve(popsize+1); // f vector starting from 0
   fvalue[0] = 0.0;
   for (j=0; j<popsize; j++) {
      I = popmember[j];
      t = 1.0 - 0.25*(fvalue[I->father_id()] + fvalue[I->mother_id()]);
      fvalue[j+1] += t;
      I->inbc = std::sqrt(t);
      for (i=j+1; i<popsize; i++) {
         I = popmember[i]; s = I->father_id(); d = I->mother_id();
         t = 0.0;
         if (s > j) t += 0.5*(I->father_inbcoef());
         if (d > j) t += 0.5*(I->mother_inbcoef());
         I->inbc = t;
         fvalue[i+1] += t*t;
      }
   }
   for (i=0; i<popsize; i++) popmember[i]->inbc = fvalue[i+1] - 1.0;
   fdone = 1;
}

void Population::inbcoef_meuwissen(void)
{
   if (!stdid) renum();
   Vector<double> Fam; Fam.reserve(popsize+1); // f vector starting from 0
   double *L;
   if(popsize>0){
     L = new double [popsize];  
   }
   else {
     L = 0;
   }
   L--;
   double *D;
   if(popsize>0){
     D = new double [popsize];  
   }
   else {
     D = 0;
   }
   D--;
   unsigned *point;
   if(popsize>0){
     point = new unsigned [popsize];
   }
   else {
     point = 0;
   }
   memset(point,'\0',sizeof(unsigned)*popsize);
   point--;
   unsigned i,j,k,s,d,ks,kd,kk,previous_s,previous_d;
   double r,fi;
   Individual *I;

   previous_s = 0;  previous_d = 0;
   Fam[0] = -1.0;
   for (i=1; i<=popsize; i++) {
      I = popmember[i-1];
      s = I->father_id();
      d = I->mother_id();

      D[i] = 0.5 - 0.25*(Fam[s] + Fam[d]);
      if (s == 0 || d == 0) {
         Fam[i] = 0.0;
      }
      else if ((s == previous_s) && (d == previous_d)) {
         Fam[i] = Fam[i-1];
      }
      else {
         D[i] = 0.5 - 0.25*(Fam[s] + Fam[d]);
         fi = -1.0;
         L[i] = 1.0;
         j = i;
         while (j != 0) {
            k = j;
            r = 0.5*L[k];
            I = popmember[k-1];
            ks = I->father_id();
            kd = I->mother_id();
            if (ks < kd) {kk = ks; ks = kd; kd = kk;}
            if (ks > 0) {
               while (point[k] > ks) k = point[k];
               L[ks] += r;
               if (ks != point[k]) {
                  point[ks] = point[k];
                  point[k] = ks;
               }
               if (kd > 0) {
                  while (point[k] > kd) k = point[k];
                  L[kd] += r;
                  if (kd != point[k]) {
                     point[kd] = point[k];
                     point[k] = kd;
                  }
               }
            }
            fi += L[j]*L[j]*D[j];
            L[j] = 0.0;
            k = j;
            j = point[j];
            point[k] = 0;
         }
         Fam[i] = fi;
      }
      previous_s = s;  previous_d = d;
   }
   L++;     
   if(L){
     delete [] L;
     L=0;
   }
   D++;     
   if(D){
     delete [] D;
     D=0;
   }
   point++; 
   if(point){
     delete [] point;
     point=0;
   }
   for (i=0; i<popsize; i++) popmember[i]->inbc = Fam[i+1];
   fdone = 1;
}

Vector<double> Population::inbcoef()
{
   if (!fdone) inbcoef_meuwissen();
   Vector<double> fvec(popsize);
   for (unsigned i=0; i<popsize; i++) fvec[i] = popmember[i]->inbcoef();
   return  fvec;
}

doubleMatrix Population::rela(void)
{
   if (!stdid) renum();
   doubleMatrix amat(popsize,popsize);
   Individual* I;
   
   unsigned i,j,a,s,d;
   double *srow, *drow;
   for (i=0; i<popsize; i++) {
      memset(amat[i],'\0',sizeof(double)*popsize);
      I = popmember[i];
      a = I->id() - 1; s = I->father_id(); d = I->mother_id();
      amat[a][a] = 1.0;
      if (s && d) {
         --s; --d;
         srow = amat[s]; drow = amat[d];
         for (j=0;j<a; j++) amat[a][j] = amat[j][a] = (*srow++ + *drow++)/2.0;
         amat[a][a] += amat[s][d]/2.0;
      }
      else if (!s && d) {
         drow = amat[d-1];
         for (j=0; j<a; j++) amat[a][j] = amat[j][a] = (*drow++)/2.0;
      }
      else if (s && !d) {
         srow = amat[s-1];
         for (j=0; j<a; j++) amat[a][j] = amat[j][a] = (*srow++)/2.0;
      }
   }
   return amat;
}
} ///////// end of namespace matvec

