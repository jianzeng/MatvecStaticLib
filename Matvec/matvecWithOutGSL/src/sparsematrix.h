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

#ifndef MATVEC_SPARSEMATRIX_H
#define MATVEC_SPARSEMATRIX_H

#include "doublematrix.h"
#include "dblock.h"

namespace matvec {
/*!
   \struct NodeStruct
   \brief This is used in SparseMatrix
*/
struct NodeStruct {
  Vector<unsigned>  adj_list;
  Vector<double>    a;
};

//BRS
// This makes qsort work in C++ from Boland site (www.inprise.com)
typedef int (*cmp_func) (const void *, const void *);
int comp(const int *i, const int *j);
// {return *i-*j;};
//BRS

/*!
   \class SparseMatrix mvsparsematrix.h
   \brief a sparse matrix is a matrix with very sparse elements

   \sa Matrix
*/

  typedef int idxtype; //from struct.h

class SparseMatrix {  
   friend class GLMM;
   protected:
      Vector<NodeStruct> node_list;
      Vector<unsigned> mask,masklist;
      Vector<unsigned> order;
      Vector<unsigned> inv_order;
      Vector<double> diag;
      unsigned   dim, nonzero, max_nz, hsize, insertend, nsp;
      int        *iiv, *jjv, *perm, *iperm;
  idxtype *xadj,*adjncy;
      double     *aav, *rsp;
      int        *srt_hash,*hash_srt;
      char       solver_name;
      void       copyfrom(const SparseMatrix& A);
      unsigned   factor_done,elim_done,initial_lij_done,initialize_node_list_done;
      unsigned   logdet_done;
      unsigned   rank;

   public:
      double     logdeterm;
      Vector<double> info_vec;

      SparseMatrix(void);                                   // Constructor 1
      SparseMatrix(int n,int max_nz);                       // Constructor 2
      SparseMatrix(const SparseMatrix& A);                  // Constructor 3
      ~SparseMatrix(void){release();}                      // Destructor

      const SparseMatrix&   operator=(const SparseMatrix& A);
      double                operator()(const int i,const int j);

      friend int   operator==(const SparseMatrix& A, const SparseMatrix& B);
      friend int   operator!=(const SparseMatrix& A, const SparseMatrix& B);

      friend SparseMatrix   operator*(const SparseMatrix& A, const double s);
      friend SparseMatrix   operator*(const double s, const SparseMatrix& A);
      friend Vector<double>         operator*(const SparseMatrix& A, const Vector<double>& V);

      friend std::ostream&  operator<<(std::ostream& stream, const SparseMatrix& A);

      SparseMatrix&       resize(const int n,const int maxnz);
      double     q(const Vector<double>& v1,const Vector<double>& v2);
      double     q(const double* x,const double *y);
      double     q(const double* x,const double *y,const int i1,const int i2);
      double     qTLW(const Vector<double>& v1,const Vector<double>& v2);
      double     qTLW(const double* x,const double *y);
      double     qTLW(const double* x,const double *y,const int i1,const int i2);
      unsigned   num_rows(void) const {return dim;};
      unsigned   nz(void) const {return nonzero;};
      unsigned   close(void);
      unsigned   factorization(const int path=5,const double stopval=1.0e-10);
      double*    a(void) {return aav-1;};
      doubleMatrix     dense(unsigned r1, unsigned r2, unsigned c1, unsigned c2) const;
      doubleMatrix     dense(void){return dense(1,dim,1,dim);};
      Vector<double>     row(const int ithrow);
      double     logdet(unsigned &irank);
      double     logdet(unsigned *irank){return logdet(*irank);};
      double     logdet(void);
      unsigned   irank(void);
      int        solve(double* x, const double* b, const std::string &method = "ysmp",
                       double relax=1.0,double stopval=0.001,int mxiter=1000);
      int        solve(Vector<double>& x,const Vector<double>& b,const std::string &method = "ysmp",
                       double relax=1.0,double stopval=0.001,int mxiter=1000);
      int*       ia(void) {return iiv-1;};
      int*       ja(void) {return jjv-1;};
      int        reset(const int i,const int j,const double value);
      int        insert(const int i,const int j,const double value);
      void       reorder(const int path=2);
      void       add_diag(const int ibeg,const int iend,const double r);
      void       gibbs_iterate(double *x,const double *r,double *w,
                               const double se=1.0);
      void       display(unsigned r1, unsigned r2, unsigned c1, unsigned c2,
                         const std::string &msg = "");
      void       display(const std::string &msg = ""){display(1,dim,1,dim,msg);};
      void       release_nodestuff(void);
  void       save(const std::string &fname,const int io_mode=std::ios::out) const;//SDK|ios::noreplace) const;
      void       release(void);
      void       mv(const Vector<double>& v, Vector<double>& result);
      void       mv(const double *v, const unsigned n, double *result);
      double     getaij(const int i,const int j);
      int        initialize_node_list(void);
      int        minfilsymelm_node_list(void);
      int        elim_node(unsigned i, Vector<unsigned> Aorder);
      void       display_node_list(void);
      void       display_order(void);
      void       display_inv_order(void);
      int        initial_lij(void);
      void       display_lij(void);
      int        factor(void);
      int        solv(double* x, const double* rhs);
      int        solvrow(double* x, const double* rhs,unsigned row,Vector<double> &sol,Vector<double> &temp_rhs,unsigned *row_pt,unsigned *col_pt);
      Vector<double>     solv(Vector<double> rhs);
 
};

extern void sparse_dense(unsigned n,int* ia,int* ja,double* a,double** matx);
}
#endif
