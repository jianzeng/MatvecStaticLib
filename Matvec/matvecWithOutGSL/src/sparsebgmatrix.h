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

#ifndef MATVEC_BGSPARSEMATRIX_H
#define MATVEC_BGSPARSEMATRIX_H

#include "doublematrix.h"
#include "dblock.h"
#include "session.h"
#include "bgmatrix.h"
#include "bg.h"

namespace matvec {
/*!
   \struct BGNodeStruct
   \brief This is used in SparseMatrix
*/
struct BGNodeStruct {
  Vector<unsigned>  adj_list;
  Vector<BG>    a;
};

//BRS
// This makes qsort work in C++ from Boland site (www.inprise.com)
typedef int (*cmp_func) (const void *, const void *);
int comp1 (const int *i, const int *j);
// {return *i-*j;};
//BRS

/*!
   \class SparseBGMatrix sparsebgmatrix.h
   \brief a SparseBGMatrix is a sparse matrix of BG objects

   \sa Matrix
*/

  typedef int idxtype; //from struct.h

class SparseBGMatrix {  
   friend class GLMM;
 protected:
   // public:
      Vector<BGNodeStruct> node_list;
      Vector<unsigned> mask,masklist;
      Vector<unsigned> order;
      Vector<unsigned> inv_order;
      Vector<BG> diag;
      unsigned   dim, nonzero, max_nz, hsize, insertend, nsp;
      int        *iiv, *jjv, *perm, *iperm;
      idxtype *xadj,*adjncy;
      BG     *aav, *rsp;
      int        *srt_hash,*hash_srt;
      char       solver_name;
      void       copyfrom(const SparseBGMatrix& A);
      unsigned   factor_done,elim_done,initial_lij_done,initialize_node_list_done;
      unsigned   logdet_done;
      unsigned   rank;

 public:

      BG     logdeterm;
      Vector<BG> info_vec;

      SparseBGMatrix(void);                                   // Constructor 1
      SparseBGMatrix(int n,int max_nz);                       // Constructor 2
      SparseBGMatrix(const SparseBGMatrix& A);                  // Constructor 3
      ~SparseBGMatrix(void){release();}                      // Destructor

      const SparseBGMatrix&   operator=(const SparseBGMatrix& A);
      BG                operator()(const int i,const int j);

      friend int   operator==(const SparseBGMatrix& A, const SparseBGMatrix& B);
/* 	{ */
/* 	if (A.dim != B.dim || A.nonzero != B.nonzero) return 0; */
/* 	unsigned n = A.dim; */
/* 	unsigned m = A.nonzero; */
/* 	if (m*n == 0) return 1; */

/* 	unsigned i; */
/* 	BG *A_a = A.aav-1; */
/* 	BG *B_a = B.aav-1; */

/* 	// one loop instead of the two loops of Tianlin */
/* 	for (i=1; i<=m; i++) { */
/* 	  if (A_a[i] != B_a[i])  return 0; */
/* 	} */
/* 	return 1; */
/*       } */
      friend int   operator!=(const SparseBGMatrix& A, const SparseBGMatrix& B);

      friend SparseBGMatrix   operator*(const SparseBGMatrix& A, const BG s);
      friend SparseBGMatrix   operator*(const BG s, const SparseBGMatrix& A);
      friend Vector<BG>         operator*(const SparseBGMatrix& A, const Vector<BG>& V);

      friend std::ostream&  operator<<(std::ostream& stream, const SparseBGMatrix& A);

      SparseBGMatrix&       resize(const int n,const int maxnz);
      BG     q(const Vector<BG>& v1,const Vector<BG>& v2);
      BG     q(const BG* x,const BG *y);
      BG     q(const BG* x,const BG *y,const int i1,const int i2);
      BG     qTLW(const Vector<BG>& v1,const Vector<BG>& v2);
      BG     qTLW(const BG* x,const BG *y);
      BG     qTLW(const BG* x,const BG *y,const int i1,const int i2);
      unsigned   num_rows(void) const {return dim;};
      unsigned   nz(void) const {return nonzero;};
      unsigned   close(void);
      unsigned   factorization(const int path=5,const double stopval=1.0e-10);
      BG*    a(void) {return aav-1;};
      doubleMatrix     dense(unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;
      doubleMatrix     dense(void){return dense(1,dim,1,dim);};
      Vector<BG>     row(const int ithrow);
      BG     logdet(unsigned &irank);
      BG     logdet(unsigned *irank){return logdet(*irank);};
      BG     logdet(void);
      unsigned   irank(void);
      int        solve(BG* x, const BG* b, const std::string &method = "ysmp",
                       double relax=1.0,double stopval=0.001,int mxiter=1000);
      int        solve(Vector<BG>& x,const Vector<BG>& b,const std::string &method = "ysmp",
                       double relax=1.0,double stopval=0.001,int mxiter=1000);
      int*       ia(void) {return iiv-1;};
      int*       ja(void) {return jjv-1;};
      int        reset(const int i,const int j,const BG value);
      int        insert(const int i,const int j,const BG value);
      void       reorder(const int path=2);
      void       add_diag(const int ibeg,const int iend,const BG r);
      void       display(unsigned r1, unsigned r2, unsigned c1, unsigned c2,
                         const std::string &msg = "");
      void       display(const std::string &msg = ""){display(1,dim,1,dim,msg);};
      void       release_nodestuff(void);
      void       save(const std::string &fname,const int io_mode=std::ios::out) const;//SDK|ios::noreplace) const;
      void       release(void);
      void       mv(const Vector<BG>& v, Vector<BG>& result);
      void       mv(const BG *v, const unsigned n, BG *result);
      BG         getaij(const int i,const int j);
      int        initialize_node_list(void);
      int        minfilsymelm_node_list(void);
      int        elim_node(unsigned i, Vector<unsigned> Aorder);
      void       display_node_list(void);
      void       display_order(void);
      void       display_inv_order(void);
      int        initial_lij(void);
      void       display_lij(void);
      int        factor(void);
      int        solv(BG* x, const BG* rhs);
      int        solvrow(BG* x, const BG* rhs,unsigned row,Vector<BG> &sol,Vector<BG> &temp_rhs,unsigned *row_pt,unsigned *col_pt);
      Vector<BG>     solv(Vector<BG> rhs);
 
};

extern void sparse_dense(unsigned n,int* ia,int* ja,BG* a,BG** matx);
}
#endif
