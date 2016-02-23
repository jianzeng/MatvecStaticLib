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

#ifndef MATVEC_TERMLIST_H
#define MATVEC_TERMLIST_H

#include "geneticdist.h"

namespace matvec {

class FieldStruct;

///////////////////  ModelTerm class /////////////////////////////
/*!
   \class ModelTerm  TermList.h
   \brief a class to hold term definition for a model

   \sa Model
*/

class ModelTerm {
   friend class TermList;
   friend class Model;
   friend class GLMM;
 public:
   unsigned numfact, start, numtrait,numlevel;
   // steve's stuff
   unsigned base_effect;
   // defaults to term index
   // changed first correlated effect  for later correlated effects
  unsigned vartype;
  // defaults 0 for Unstructured
  // a value of 1 for Estimating a correlation matrix with fixed variances 
  // Used to designate covariance structure
   unsigned corr_link;
   // defaults to -1
   // points next correlated term
   unsigned pos_vec;
   // Position in link list;
   // defaults to 0
   Vector<int>  trait_map;
   // Maps ?traits? within correlated effects
   // defaults to Vectors of length numtrait with element i
   unsigned nt_vec;
   // defaults to numtrait for uncorrelated traits
   // numtrait*ncorr for first correlated trait
   // 0 for later correlated traits
   char             myclassi;                // fixed 'F' is the default
   Vector<int>      trait,factorindx;
   Vector<unsigned> addr;
   
   GeneticDist *prior;
   ModelTerm(void);
   ModelTerm(const ModelTerm& A);
   ~ModelTerm(void) {release();}
   
   const ModelTerm& operator=(const ModelTerm& A);
   int      partial_match(const std::string &effectname,FieldStruct* FS);
   int      match(const std::string &termname, const std::string &sep,FieldStruct* FS);
   void     copyfrom(const ModelTerm& A);
   void     resize(const unsigned ne,const unsigned nt);
   void     input(const std::string &str, FieldStruct* FS,const unsigned nf);
   void     classi(const char c) {myclassi = c;}
   char     classi(void) const {return myclassi;}
   unsigned order(void) const {return numfact;}
   unsigned ntrait(void) const {return numtrait;}
   unsigned size(void) const {return numfact;}
   unsigned nlevel(void) const {return numlevel;}
   unsigned startaddr(void) const {return start;}
   void     variance(const doubleMatrix& A);
   void     release(void);
};

/////////////////// TermList class /////////////////////////////
/*!
   \class TermList  TermList.h
   \brief a class to hold a list of terms in a model

   \sa Model
*/
class TermList {
   friend class Model;
   protected:
      unsigned numterm, numtrait;
      ModelTerm *termlist;

   public:
      TermList(void){numterm=0; numtrait=0; termlist = 0;}
      TermList(const unsigned ntrm, const unsigned ntrt);
      TermList(const TermList& T); 
      ~TermList(void) {release();}

      ModelTerm& operator[](const unsigned i) {return termlist[i];}
      const ModelTerm& operator()(const unsigned i) const{return termlist[i];}
      const TermList& operator=(const TermList& T) {copyfrom(T);return *this;}

      unsigned nterm(void) const {return numterm;}
      unsigned ntrait(void) const {return numtrait;}
      unsigned maxorder(void) const;

      void     copyfrom(const TermList& T);
      void     swap(const unsigned t1, const unsigned t2);
      void     resize(const unsigned ntrm, const unsigned ntrt);
      void     release(void);

      int      index(const std::string &termname,FieldStruct* FS) const;
};
}
#endif
