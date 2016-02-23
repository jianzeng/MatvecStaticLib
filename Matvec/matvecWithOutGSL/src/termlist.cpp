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

#include <vector>
#include <string>

#include "util.h"
#include "termlist.h"
#include "fieldstruct.h"

namespace matvec {

///////////////////////// ModelTerm class ///////////////////////
ModelTerm::ModelTerm(void)
{
   start      = 0;
   numfact    = 0;
   numlevel   = 0;
   numtrait   = 0;
   myclassi   = 'F',
   addr       = 0;
   vartype    = 0;
   prior      = 0;
}

ModelTerm::ModelTerm(const ModelTerm& A)
{
   start      = 0;
   numfact    = 0;
   numlevel   = 0;
   numtrait   = 0;
   myclassi   ='F',
   addr       = 0;
   vartype    = 0;
   prior      = 0;
   copyfrom(A);
}

const ModelTerm& ModelTerm::operator=(const ModelTerm& A)
{
   if (this != &A) copyfrom(A);
   return *this;
}

void ModelTerm::copyfrom(const ModelTerm& A)
{
   if (this == &A) return;
   resize(A.numfact,A.numtrait);
   base_effect = A.base_effect;
   pos_vec     = A.pos_vec;
   nt_vec      = A.nt_vec;
   corr_link   = A.corr_link;
   trait_map   = A.trait_map;
   start       = A.start;
   numlevel    = A.numlevel;
   myclassi    = A.myclassi;
   prior       = A.prior;
   addr        = A.addr;
   vartype     = A.vartype;
   trait       = A.trait;
   factorindx  = A.factorindx;
}

void  ModelTerm::resize(const unsigned ne, const unsigned nt)
{
   if (numfact != ne) {
      numfact = ne;
      //addr.reserve(numfact);
      factorindx.reserve(numfact);
   }
   if (numtrait != nt) {
      numtrait = nt;
      trait.resize(numtrait);
      addr.resize(numtrait);
   }
}

void ModelTerm::variance(const doubleMatrix& A)
{
  doubleMatrix *var = prior->var_matrix();
  *var = A;
}

void ModelTerm::release(void)
{
   numfact = 0; numtrait = 0;
}

int ModelTerm::match(const std::string &termname,const std::string &sep,FieldStruct* FS)
{
   std::string s(termname);
   std::string::size_type begidx,endidx;
   std::vector<std::string> tmpstr;
   /////  s.split(n,sep);
   begidx = s.find_first_not_of(sep);
   while(begidx != std::string::npos) {
      endidx = s.find_first_of(sep,begidx);
      if (endidx == std::string::npos) endidx = s.length();
      tmpstr.push_back(s.substr(begidx,endidx - begidx));
      begidx = s.find_first_not_of(sep,endidx);
   }
   unsigned n = tmpstr.size();
   if (numfact != n) return 0;
   int i,j, k=1;
   for (i=0; i<numfact; i++) {
      for (j=0; j<numfact; j++) if (tmpstr[i] == FS[factorindx[j]].name()) break;
      if (j == numfact) { k=0; break;}
   }
   return k;
}

int ModelTerm::partial_match(const std::string & effectname,FieldStruct* FS)
{
   int i,k=-1;
   for (i=0; i<numfact; i++) {
      if (effectname == FS[factorindx[i]].name()) break;
   }
   if (i < numfact)  k = i;
   return k;
}

void ModelTerm::input(const std::string &str, FieldStruct* FS,const unsigned nf)
{
   int k;
   std::string T(str),sep("* ");
   std::string::size_type begidx,endidx;
   std::vector<std::string> eff_name;
   /////////// T.split(ne,"* ");
   begidx = T.find_first_not_of(sep);
   while(begidx != std::string::npos) {
      endidx = T.find_first_of(sep,begidx);
      if (endidx == std::string::npos) endidx = T.length();
      eff_name.push_back(T.substr(begidx,endidx - begidx));
      begidx = T.find_first_not_of(sep,endidx);
   }
   unsigned ne = eff_name.size();
   resize(ne,numtrait);
   for (unsigned i=0; i<ne; i++) {
      k = fieldindex(eff_name[i],FS,nf);    // it cannot found
      factorindx[i] = k;
   }
}

///////////////////////// TermList class /////////////////////////////////

TermList::TermList(const TermList& T) {
  termlist=0; 
  numterm=0;
  numtrait=0;
  copyfrom(T);
}

TermList::TermList(const unsigned ntrm, const unsigned ntrt)
{
   numterm = 0;
   numtrait = 0;
   termlist = 0;
   resize(ntrm,ntrt);
}

void TermList::resize(const unsigned ntrm, const unsigned ntrt)
{
   unsigned i;
   if (numterm != ntrm) {
      numterm = ntrm;
      numtrait = ntrt;
      if (termlist) {
	delete [] termlist;
	termlist=0;
      }
      if(numterm>0){
	termlist = new ModelTerm [numterm];
      }
      else {
	termlist = 0;
      }
      check_ptr(termlist);
      for (i=0; i<numterm; i++) termlist[i].resize(0,numtrait);
   }
   else {
      if (numtrait != ntrt) {
         numtrait = ntrt;
         for (i=0; i<numterm; i++) termlist[i].resize(0,numtrait);
      }
   }
}

void TermList::release(void)
{
   numterm = 0; numtrait = 0;
   if (termlist) {
      delete [] termlist;
      termlist = 0;
   }
}

unsigned TermList::maxorder(void) const
{
   unsigned retval = 0;
   for (unsigned i=0; i<numterm; i++) {
      if (termlist[i].order() >retval) retval = termlist[i].order();
   }
   return retval;
}

int TermList::index(const std::string &termname,FieldStruct* FS) const
{
   std::string tname(termname);
   std::string::size_type begidx;
   //////// tname.replace_all(')', ' ');
   begidx = tname.find(")",0);
   while(begidx != std::string::npos) {
      tname.replace(begidx,1," ");
      begidx = tname.find(")",begidx);
   }

   /////// tname.replace_all('(', '*');
   begidx = tname.find("(",0);
   while(begidx != std::string::npos) {
      tname.replace(begidx,1,"*");
      begidx = tname.find("(",begidx);
   }

   for (int i=0; i<numterm; i++) if (termlist[i].match(tname,"* ",FS)) return i;
   return -1;
}

void TermList::copyfrom(const TermList& T)
{
   if (this == &T) return;
   resize(T.numterm,T.numtrait);
   for (int i=0; i<numterm; i++) termlist[i].copyfrom(T(i));
}

void TermList::swap(const unsigned t1, const unsigned t2)
{
   if (t1 == t2) return;
   if (t1 >= numterm || t2 >= numterm) throw exception("TermList::swap(): out of range");
   ModelTerm T = termlist[t1];
   termlist[t1] = termlist[t2];
   termlist[t2] = T;
}

} ///////// end of namespace matvec

