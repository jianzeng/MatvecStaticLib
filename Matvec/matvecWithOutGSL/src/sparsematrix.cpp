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

//DO_MINORD also do the MINORD
#define DO_MINORD
#define NEW_ELIMNODE

#include "session.h"
#include "util.h"
#include "sparsematrix.h"

#define DEBUG
namespace matvec {

int comp(const int *i, const int *j) {
  return *i-*j;
}
/* calls to ysmpack - may need equivalent?
extern void odrvd(int n,int *ia,int *ja,double *a,int *p,int *ip,
           int nsp,int *isp,int path,int& flag);
extern int sdrvd(int n,int *p,int *ip,int *ia,int *ja,double *a,
           const double *b,double *z,int nsp,int *isp,double *rsp,double& esp,
           int path,int& flag,double tol,int& rank,double& lgdet);
extern void detdrv(int n,int *p,int *ip,int *ia,int *ja,double *a,
                int nsp,int *isp,double *rsp,double& esp,int path,
                int& flag,double tol,int& rank,double& lgdet);
*/
extern void gs_ioc(int *ia,int *ja,double *a,const double *rhs,double *sol,
                   const int n,const double relax,const double stopval,
                   const int mxiter,const double tol,unsigned& cov);
extern long next_prime(long n);
// extern double snorm(void);

static int squeeze(unsigned len,int *iiv,int *jjv,double *aav);
static void hsort_ija(unsigned n,int *iiv,int *jjv,double *aav);

void gs_ioc(int *ia,int *ja,double *a,const double *rhs,double *sol,const int n,
            const double relax,const double stopval,const int mxiter,
            const double tol,unsigned& cov)
  //***********************************************************************
  //   12/26/91, University of Illinois
  //   Copyright (C) 1991, 1992, 1993, 1994 Tianlin Wang
  //
  //   This routine solves a linear system by Gauss-Seidal(GS) or Successive
  //   Over-Relaxation(SOR) Iteration on Coefficients (IOC).
  //   the (ia,ja,a) holds only half of coefficient matrix
  //   (either upper or lower, not mixed)
  //    the diagonal element is not necessary to be the first (or last)
  //    one in each row.
  //
  // input
  //      ia,ja,a  the complete coefficient matrix in half
  //               stored format (upper triangular).
  //      rhs      the right-hand-side, rhs[1..n]
  //      n        the dimension of the system of equations
  //      mxiter   maximum number of iterations allowed
  //      stopval  the accuracy at which iteration stops.
  //      relax    relaxation coefficient, relax=1.2 seems good for some case
  //      tol      tolerence
  //
  // input/output
  //      sol      the initial solutions as input, returned
  //               final solutions as output. sol[1..n]
  //      niter    actual number of iterations
  // working space
  //      adj_rhs  adjusted rhs
  // ***********************************************************************
{
   int j,row,col;
   double oldsol,newsol,local,diag, cmax,cval;
   Vector<double> adj_rhs(n);
   int niter = 0;
   cov = 0;

   std::cout << std::endl;
   do {                             // now iteration begins
      cmax = 0.0;
      for (row=1; row<=n; row++) adj_rhs[row-1] = rhs[row];
      for (row=1; row<=n; row++) {
         oldsol = sol[row];
         local = 0.0;
         diag = 0.0;             // subtract all the other solutions*count:
         for (j=ia[row]; j<=ia[row+1]-1; j++) {
            col = ja[j];
            if (col != row) {
               local += a[j]*sol[col];
            }
            else {
               diag = a[j];
            }
         }
         if (diag > tol) {
            newsol = (adj_rhs[row - 1] - local)/diag;   // new solution for GS
            cval = (newsol - oldsol)*relax;         // change vs old solution
            if (fabs(cval) > cmax) cmax=fabs(cval);
            newsol = sol[row] = oldsol + cval;      // new solution for SOR
            newsol = sol[row];
            for (j = ia[row]; j<=ia[row+1]-1; j++) {
               col = ja[j];
               adj_rhs[col - 1] -= a[j]*newsol;         // adjust the adj_rhs
            }
         }
         else if (diag > -tol) {
            throw exception(" zero-diagonal found in gs_ioc.c");
         }
         else {
            throw exception(" in gs_ioc(), matrix is not psd");
         }
      }
      niter += 1;
      if ( (niter % 10) ==0) {
         std::cout << " IOC: # of iter = " << niter << ", max_change = "
              << cmax << std::endl;
      }
   } while (niter <= mxiter && cmax > stopval);
   if (cmax <= stopval) cov=1;
}

/**
 * A constructor to create an SparseMatrix object.
 */
SparseMatrix::SparseMatrix(void)
{
   dim         = 0;
   nonzero     = 0;
   max_nz      = 0;
   hsize       = 0;
   insertend   = 0;
   factor_done = 0;
   elim_done   = 0;
   logdet_done = 0;
   initial_lij_done = 0;
   initialize_node_list_done = 0;
   iiv         = 0;
   jjv         = 0;
   xadj        = 0;
   adjncy      = 0;
   aav         = 0;
   srt_hash    = 0;
   hash_srt    = 0;

   solver_name = 'N';
   info_vec.reserve(2);
}

/**
 * A constructor to create an identical SparseMatrix object to A.
 */
SparseMatrix::SparseMatrix(const SparseMatrix& A)
{
   dim         = 0;
   nonzero     = 0;
   max_nz      = 0;
   hsize       = 0;
   insertend   = 0;
   factor_done = 0;
   elim_done   = 0;
   logdet_done = 0;
   initial_lij_done = 0;  
   initialize_node_list_done = 0;
   iiv         = 0;
   jjv         = 0;
   xadj        = 0;
   adjncy      = 0;
   aav         = 0;
   srt_hash    = 0;
   hash_srt    = 0;

   solver_name  = 'N';
   copyfrom(A);
}

/**
 * A constructor to create an SparseMatrix object of size
 *  n by n and maximum non-zero elemenat maxnz.
 */
SparseMatrix::SparseMatrix(int n,int maxnz)
{
   if (n<0 || maxnz<0) throw exception("SparseMatrix(m,n): requring nonnegative m,n");
   dim = n;
   max_nz = maxnz;
   nonzero = 0;
   factor_done = 0;
   elim_done = 0;
   logdet_done = 0;
   initial_lij_done = 0;
   initialize_node_list_done = 0;
   if (maxnz > 0) {
      hsize   = static_cast<unsigned>(next_prime(static_cast<long>(1.25*maxnz + 50)));
   }
   iiv     = new int [hsize];
   jjv     = new int [hsize];
   aav     = new double [hsize];
   hash_srt = new int [hsize+1];
   srt_hash = new int [hsize+1];
   nonzero=0;
   solver_name  = 'N';
   insertend = 0;
   info_vec.reserve(2);

   memset(iiv,'\0',sizeof(int)*hsize);
   memset(jjv,'\0',sizeof(int)*hsize);
   memset(aav,'\0',sizeof(double)*hsize);
   if(hsize){
     memset(hash_srt,'\0',sizeof(int)*(hsize+1));
     memset(srt_hash,'\0',sizeof(int)*(hsize+1));
   }
   nonzero=0;
   for (unsigned i=1; i<=dim; i++) insert(i,i,0.0);
}

/**
 * Get a copy from A
 */
void SparseMatrix::copyfrom(const SparseMatrix& A)
{
   if (this == &A) return;
   if (dim != A.dim || max_nz != A.max_nz) resize(A.dim,A.max_nz);
   nonzero      = A.nonzero;
   hsize        = A.hsize;
   info_vec     = A.info_vec;
   insertend    = A.insertend;
   solver_name  = A.solver_name;
   factor_done  = A.factor_done;
   elim_done    = A.elim_done;
   initial_lij_done = A.initial_lij_done;
   initialize_node_list_done = A.initialize_node_list_done;

   memcpy(iiv,A.iiv,sizeof(int)*hsize);
   memcpy(jjv,A.jjv,sizeof(int)*hsize);
   memcpy(aav,A.aav,sizeof(double)*hsize);
   if(hsize){
     memcpy(hash_srt,A.hash_srt,sizeof(int)*(hsize+1));
     memcpy(srt_hash,A.srt_hash,sizeof(int)*(hsize+1));
   }
}

/**
 * Retrieve element (i,j).
 */
double SparseMatrix::operator()(const int i,const int j) {

  if (i>=1 && i<=dim && j>=1 && j<=dim) {
    return getaij(i,j);
  }
  else {
    throw exception("SparseMatrix(): subscript out of range");
  }
  return 0.0;
}

/**
 * Assignment operator.
 */
const SparseMatrix& SparseMatrix::operator=(const SparseMatrix& A)
{
  copyfrom(A);
  return *this;
}

/**
 * Resize the object with n by n and maximum  non-zero elements maxnz.
 *
 * if  n and maxnz are the same as the previous ones,
 * then we assume that the matrix has the same pattern
 * of non-zeros. So, we will just change the numerical
 * values using initial_lij and then factor etc.
 * This is useful for DFREML
 */
SparseMatrix& SparseMatrix::resize(const int n,const int maxnz)
{
  // std::cout << "\n n " << n << " dim " << dim << " max_nz " << max_nz << " maxnz " << maxnz <<"\n";
   if ( dim != n || max_nz != maxnz) {
      dim     = n;
      max_nz  = maxnz;
      if (max_nz == 0) {
         hsize = 0;
      }
      else {
         hsize = static_cast<unsigned>(next_prime(static_cast<long>(1.25*max_nz + 50)));
      }
      if (iiv) {
	delete [] iiv;
	iiv=0;
      }
      if (jjv) {
	delete [] jjv;
	jjv=0;
      }
      if (aav) {
	delete [] aav;
	aav=0;
      }
      if (hash_srt) {
	delete [] hash_srt;
	hash_srt=0;
      }
      if (srt_hash) {
	delete [] srt_hash;
	srt_hash=0;
      }

      iiv = new int [hsize];
      jjv = new int [hsize];
      aav = new double [hsize];
      hash_srt = new int [hsize+1];
      srt_hash = new int [hsize+1];
      nonzero=0;
      // std::cout <<"\n release nodes \n";
      release_nodestuff(); 
      solver_name = 'N';
   }  
   memset(iiv,'\0',sizeof(int)*hsize);
   memset(jjv,'\0',sizeof(int)*hsize); 
   if(hsize){
     memset(hash_srt,'\0',sizeof(int)*(hsize+1));
     memset(srt_hash,'\0',sizeof(int)*(hsize+1));
   }
   memset(aav,'\0',sizeof(double)*hsize);
   initial_lij_done=0;
   factor_done=0;
   logdet_done = 0;
   insertend   = 0;
   nonzero     = 0;

   for (unsigned i=1; i<=dim; i++) insert(i,i,0.0);
   return *this;
}

/**
 * Relational operator. Return 1 if A and B are identical,
 * otherwise return 0.
 */
int operator==(const SparseMatrix& A, const SparseMatrix& B)
{
   if (A.dim != B.dim || A.nonzero != B.nonzero) return 0;
   unsigned n = A.dim;
   unsigned m = A.nonzero;
   if (m*n == 0) return 1;

   unsigned i;
   double *A_a = A.aav-1;
   double *B_a = B.aav-1;

   // one loop instead of the two loops of Tianlin
   for (i=1; i<=m; i++) {
     if (A_a[i] != B_a[i])  return 0;
   }
   return 1;
}

/**
 * Not equal operator.
 */
int operator!=(const SparseMatrix& A, const SparseMatrix& B)
{return !operator==(A,B);}

/**
 * Retrieve the i'th row. Note that the first row is 0.
 */
Vector<double> SparseMatrix::row(const int ithrow)
{
   if ( dim == 0 ) throw exception("SparseMatrix is empty");
   if (ithrow<1 || ithrow>dim) throw exception("SparseMatrix.row(k): out of range");
   Vector<double> out(dim);
   int *ia = iiv-1;
   int *ja = jjv-1;
   double *a = aav-1;
   double *o_pt = out.begin()-1;
   for (unsigned j=ia[ithrow]; j<ia[ithrow+1]; j++) o_pt[ja[j]] = a[j];
   return out;
}

/**
 * SparseMatrix times vector operation.
 */
Vector<double> operator*(const SparseMatrix& A, const Vector<double>& v)
{
  if ( A.dim == 0) throw exception("SparseMatrix is empty");
  if (A.dim != v.size()) throw exception("SparseMatrix*Vector<double>, size not conformable");
  Vector<double> vec_out;
  int n = A.dim;
  int *A_ia = A.iiv-1;   /* offset by 1 */
  int *A_ja = A.jjv-1;
  double *A_a = A.aav-1;
  double *v_ve = v.begin()-1;
  double tmp,x;
  unsigned k,i,j;

  vec_out.reserve(n);
  double *out = vec_out.begin();

  out--; // offset by 1
  int ksrt=0;       //SDK
  for(ksrt=0;ksrt<A.nonzero;ksrt++) { //SDK
    k=A.srt_hash[ksrt];  //SDK
  //  for (k=1; k<=A.hsize; k++) {  //SDK
    if (i = A_ia[k]) {
        j = A_ja[k];
        x = A_a[k];
	out[i] = out[i] + x*v_ve[j];
        if (j>i) {
	  out[j] = out[j] + x*v_ve[i];
	}
    }
  }
  return vec_out;
}

void SparseMatrix::mv(const Vector<double> &v, Vector<double> &result)
{
   if (dim != v.size()) throw exception("SparseMatrix::mv(): sizenot conformable");
   if (v.begin() == result.begin()) throw exception("SparseMatrix::mv(): can't be done in situ");
   result.reserve(dim);
   this->mv(v.begin(),dim,result.begin());
}

void SparseMatrix::mv(const double *v, const unsigned n, double *result) {
   /////////////////////////////////////////////////////////////////////
   // REQUIREMENTS:  both v and result must have n(dim) elements (space)
   // modified 8/4/1998
   //////////////////////////////////////////////////////////////////////

  if (dim != n) throw exception("SparseMatrix::mv(): not conformable");
  memset(result,'\0',sizeof(double)*n);
  int *A_ia   = iiv-1;   // offset by 1
  int *A_ja   = jjv-1;
  double *A_a = aav-1;
  const double *v_ve = v-1;
  unsigned i,j,k;
  double x;

  //  result.resize(n,0.0);

  result--;      // offset by 1
  //SDK
  int ksrt;
  for(ksrt=0;ksrt<nonzero;ksrt++) {
    k=srt_hash[ksrt];
    // for (k=1; k<= hsize; k++) {
    //SDK
    if (i = A_ia[k]) {
        j = A_ja[k];
        x = A_a[k];
	result[i] = result[i] + x*v_ve[j];
        if (j>i) {
	  result[j] = result[j] + x*v_ve[i];
	}
    }
  }
  result++; // offset back
}

/**
 * Multiply SparseMatrix object by s.
 */
SparseMatrix operator*(const SparseMatrix& A, const double s)
{
   if ( A.dim == 0 ) throw exception("SparseMatrix is empty");

   unsigned new_dim = A.dim;
   unsigned new_nz = A.nonzero;
   SparseMatrix B(new_dim,A.max_nz);
   B.nonzero = new_nz;
   memcpy(B.iiv,A.iiv,sizeof(int)*A.hsize);
   memcpy(B.jjv,A.jjv,sizeof(int)*A.hsize);
   if(A.hsize){
     memcpy(B.hash_srt,A.hash_srt,sizeof(int)*(A.hsize+1));
     memcpy(B.srt_hash,A.srt_hash,sizeof(int)*(A.hsize+1));
   }
   double *A_pt = A.aav;
   double *New_pt = B.aav;
   for (unsigned i=0; i<A.hsize; i++) {
     if (A.iiv[i] > 0) {
       New_pt[i] = A_pt[i]*s;
     }
   }
   return B;
}

SparseMatrix operator*(const double s, const SparseMatrix& A)
{return operator*(A,s);}


void SparseMatrix::release(void)
{
   if (adjncy) {delete [] adjncy; adjncy = 0;}
   if (xadj) {delete [] xadj; xadj = 0;}
   if (iiv) {delete [] iiv; iiv = 0;}
   if (jjv) {delete [] jjv; jjv = 0;}
   if (aav) {delete [] aav; aav = 0;}
   if (hash_srt) {delete [] hash_srt; hash_srt = 0;}
   if (srt_hash) {delete [] srt_hash; srt_hash = 0;}
   release_nodestuff();
   dim=0;
   hsize=0;
   nonzero=0;
   max_nz=0;
   insertend=0;
}

void SparseMatrix::release_nodestuff(void)
{
   node_list.clear();
   mask.clear();
   order.clear();
   inv_order.clear();
   diag.clear();
   initialize_node_list_done=0;
   elim_done=0;
   factor_done=0;
   initial_lij_done=0;
}

/**
 * Calculates the product v1*A*v2.
 *
 * SparseMatrix A is half storage only(upper,lower or mixed)
 */
double SparseMatrix::q(const Vector<double>& v1,const Vector<double>& v2)
{
   if (dim == 0) throw exception("SparseMatrix is empty");
   if (dim != v1.size() || dim != v2.size()) throw exception("SparseMatrix::q(v1,v2): size not conformable");
   return q(v1.begin(),v2.begin());
}

/**
 * Calculates the quadratic form x*A*y.
 *
 * SparseMatrix A is half storage (upper or lower or mixed)
 * REQUIREMENTS: x and y must have at least this->dim elements
 *               x[0] and y[0] are the first elements
*/
double SparseMatrix::q(const double* x,const double* y)
{
  if (dim == 0) exception("SparseMatrix is empty");

  int   *ia = iiv-1;
  int   *ja = jjv-1;
  double *a = aav-1;
  double tmp,qaq ;
  const double *xx,*yy;
  xx = x-1;  yy = y-1;
  unsigned i,j,k;
  //SDK
  int ksrt;
  for(qaq=0.0,ksrt=0;ksrt<nonzero;ksrt++) {
    k=srt_hash[ksrt];
    //  for (qaq = 0.0, k=1; k<=hsize; k++) {
    //SDK
    if (i = ia[k]) {
        j = ja[k];
	if (i == j) {
	  tmp = xx[i]*a[k]*yy[i];
	  qaq += tmp;
	}
	else {
	  tmp = xx[i]*a[k]*yy[j] + xx[j]*a[k]*yy[i];
	  qaq += tmp;
	}
    }
  }
  return qaq;
}

/**
 * Calculates the quadratic form x*A*y  from i1'th row to i2'th row.
 *
 * SparseMatrix A is half storage (upper or lower or mixed)
 * REQUIREMENTS: x and y must have at least this->dim elements
 *               x[0] and y[0] are the first elements.
 */
double SparseMatrix::q(const double* x,const double* y,const int i1,const int i2)
{
  if (dim == 0) throw exception("SparseMatrix::q(): SparseMatrix is empty");
  if (i1<=0 || i1>dim || i2<=0 || i2>dim || i1>i2) throw exception("SparseMatrix::q(): out of range");

  int   *ia = iiv-1;
  int   *ja = jjv-1;
  double *a = aav-1;
  double tmp,qaq ;
  const double *xx,*yy;
  xx = x-1;  yy = y-1;
  unsigned i,j,k;
  //SDK
  int ksrt;
  for (qaq = 0.0, ksrt=0; ksrt< nonzero; ksrt++) {
    k=srt_hash[ksrt];
    //  for (qaq = 0.0, k=1; k<= hsize; k++) {
    //SDK
    if (i = ia[k]) {
        j = ja[k];
	if (i == j && i >= i1 && i <= i2 ) {
	  tmp = xx[i]*a[k]*yy[i];
	  qaq += tmp;
	}
	else if (i >= i1 && i <= i2 && j >= i1 && j <= i2) {
	  tmp = xx[i]*a[k]*yy[j] + xx[j]*a[k]*yy[i];
	  qaq += tmp;
	}
    }
  }
  return qaq;
}

/**
 *  Calculates the quadratic form v1*A*v2.
 *
 *  SparseMatrix A is half storage only(upper,lower or mixed)
 */
double SparseMatrix::qTLW(const Vector<double>& v1,const Vector<double>& v2)
{
   if (dim == 0) throw exception("SparseMatrix::qTLW(): it is empty");
   if (dim != v1.size() || dim != v2.size()) throw exception("SparseMatrix::qTLW(): size not conformable");
   return q(v1.begin(),v2.begin());
}

/**
 * Calculates the quadratic form x*A*y
 *
 *   SparseMatrix A is half storage (upper or lower or mixed)
 *   REQUIREMENTS: x and y must have at least this->dim elements
 *                 x[0] and y[0] are the first elements
 */
double SparseMatrix::qTLW(const double* x,const double* y)
{
   if (dim ==0) throw exception("SparseMatrix::qTLW(): SparseMatrix is empty");
   if (!insertend) throw exception("SparseMatrix::qTLW(): first call SparseMatrix::close() first");

   int   *ia = iiv-1;
   int   *ja = jjv-1;
   double qaq,tmp, *a = aav-1;
   const double *xx,*yy;
   xx = x-1;  yy = y-1;
   int i,j,je,jj;
   if (x == y) {
      for (qaq=0.0,i=1; i<=dim; i++) {
         je = ia[i+1];
         for (j=ia[i]; j<je; j++) {
            jj = ja[j];
            tmp = xx[i]*a[j]*yy[jj];
            qaq += tmp;
            if (jj != i) qaq += tmp;
         }
      }
   }
   else {
      for (qaq=0.0,i=1; i<=dim; i++) {
         je = ia[i+1];
         for (j=ia[i]; j<je; j++) {
            jj = ja[j];
            tmp = xx[i]*a[j]*yy[jj];
            qaq += tmp;
            if (jj != i) qaq +=  xx[jj]*a[j]*yy[i];
         }
      }
   }
   return qaq;
}

/**
 * Calculates the quadratic form x*A*y  from i1'th row to i2'th row.
 *
 * SparseMatrix A is half storage (upper or lower or mixed)
 *  REQUIREMENTS: x and y must have at least this->dim elements
 *                x[0] and y[0] are the first elements
 */
double SparseMatrix::qTLW(const double* x,const double* y,const int i1,
                  const int i2)
{
   if (dim == 0) throw exception("SparseMatrix::qTLW(): SparseMatrix is empty");
   if (!insertend) throw exception("SparseMatrix::qTLW(): call SparseMatrix::close() first");
   if (i1<=0 || i1>dim || i2<=0 || i2>dim || i1>i2) throw exception("SparseMatrix::qTLW(), out of range");

   int   *ia = iiv-1;
   int   *ja = jjv-1;
   double qaq,tmp, *a = aav-1;
   const double *xx,*yy;
   xx = x-1;  yy = y-1;
   int i,j,jj,je;
   if (x == y) {
      for (qaq=0.0,i=i1; i<=i2; i++) {
         je = ia[i+1];
         for (j=ia[i]; j<je; j++) {
            jj = ja[j];
            if (jj >= i1 && jj <= i2) {     // this is necessary since M may not
               tmp = xx[i]*a[j]*yy[jj];     // be strictly upper or lower
               qaq += tmp;
               if (jj != i) qaq += tmp;
            }
         }
      }
   }
   else {
      for (qaq=0.0,i=i1; i<=i2; i++) {
         je = ia[i+1];
         for (j=ia[i]; j<je; j++) {
            jj = ja[j];
            if (jj >= i1 && jj <= i2) {     // this is necessary since M may not
               tmp = xx[i]*a[j]*yy[jj];     // be strictly upper or lower
               qaq += tmp;
               if (jj != i) qaq += xx[jj]*a[j]*yy[i];
            }
         }
      }
   }
   return qaq;
}

/**
 * Add  r to the diagonals from row ibeg to row iend.
 */
void SparseMatrix::add_diag(const int ibeg,const int iend,const double r)
{
   int   *ia = iiv-1;
   int   *ja = jjv-1;
   double *a = aav-1;
   int i,j;
   for (i=ibeg; i<=iend; i++) {
      for (j=ia[i]; j<ia[i+1]; j++) if (ja[j] == i) a[j] += r;
   }
}

/**
 * Obtain a solution with the right-hand side  b.
 *
 * @param x a solution
 * @param b right hand side
 * @param method solver
 * @param relax relaxation factor
 * @param stopval criteria for stop iteration
 * @param mxiter maximum iteration allowed
 */
int SparseMatrix::solve(Vector<double>& x,const Vector<double>& b,const std::string &method,
                           double relax, double stopval,int mxiter)
{
   if (dim != b.size()) throw exception("SparseMatrix::solve(): size not conformable");
   if (x.size() < dim ) x.reserve(dim);
   if (method == "ysmp" || method == "ysmp1") {
      if (!factor_done && !factor()) throw exception("SparseMatrix::solve- returned with error from factor");
      x = solv(b);
      return 1;
   }
   else if (method == "ioc") {
      close();
      unsigned ir=0;
      memset(x.begin(),'\0',sizeof(double)*dim);
      gs_ioc(iiv-1,jjv-1,aav-1,b.begin()-1,x.begin()-1,dim,relax,stopval,mxiter,SESSION.epsilon,ir);
      if (!ir) throw exception("SparseMatrix::solve(): not converged");
      return 1;
   }
   else if (method == "iod") {
      throw exception("SparseMatrix::solve(): not available yet");
   }
   else {
      throw exception("SparseMatrix::solve(): no such solver");
  }
  return 0;
}

/**
 * Obtain a solution with the right-hand side b.
 *
 * @param x a solution
 * @param b right hand side
 * @param method solver
 * @param relax relaxation factor
 * @param stopval criteria for stop iteration
 * @param mxiter maximum iteration allowed
 *
 * @warning x and b must have at least this->dim elements
 *  x[0] and b[0] are the first elements.
 */
int SparseMatrix::solve(double* x, const double* b,const std::string &method,
                           double relax, double stopval,int mxiter)
{
   if (method == "ysmp" || method == "ysmp1") {
      if (!factor_done && !factor()) throw exception("SparseMatrix::solve() returned with error from factor");
      solv(x,b);
      return 1;
    }
    else if (method == "ioc") {
      close();
      unsigned ir=0;
      memset(x,'\0',sizeof(double)*dim);
      gs_ioc(iiv-1,jjv-1,aav-1,b-1,x-1,dim,relax,stopval,
               mxiter,SESSION.epsilon,ir);
      if (!ir) warning("SparseMatrix::solve(): not converged");
      return 1;
    }
    else if (method == "iod") {
      throw exception("SparseMatrix::solve(): solver not available yet");
   }
   else {
      throw exception("SparseMatrix::solve():no such solver");
   }
}

/**
 * Re-order elements.
 *  @param path  an integer path specification;  values and their meanings are -
 *  <UL>
 *  <LI>       1 find minimum degree ordering only
 *  <LI>       2  find minimum degree ordering and reorder symmetrically
 *              stored matrix (used when only the nonzero entries in
 *              the upper triangle of m are being stored)
 *  <LI>       3  reorder symmetrically stored matrix as specified by
 *              input permutation (used when an ordering has already
 *              been determined and only the nonzero entries in the
 *              upper triangle of m are being stored)
 *  <LI>       4  same as 2 but put diagonal entries at start of each row
 *  <LI>       5  same as 3 but put diagonal entries at start of each row
 * </UL>
 */
void SparseMatrix::reorder(const int path)
{
   if (path ==1) {
     initialize_node_list();
     minfilsymelm_node_list();
   }
   else if (path==2) {
     initialize_node_list();
     minfilsymelm_node_list();
     initial_lij();
   }
   else if (path==3) {
     initial_lij();
   }
   else if (path==4) {
     initialize_node_list();
     minfilsymelm_node_list();
     initial_lij();
  }
   else if (path==5) {
     initial_lij();
   }
}

/**
 * Factorization.
 *
 * @param path an integer path specification;  values and their meanings are -
 * <UL>
 *  <LI> 4  perform SSF
 *  <LI> 5  perform SSF and SNF
 *  <LI> 6  perform SNF only (isp/rsp is assumed to have been set
 *             up in an earlier call to sdrv (for SSF))
 * </UL>
 */
unsigned SparseMatrix::factorization(const int path,const double tol)
{
   if (path < 4) throw exception("SparseMatrix::factorization(path): path >= 4 is required");

   if (path==4) {
     initialize_node_list();
     minfilsymelm_node_list();
     initial_lij();
     factor();
   }
   if (path==6){
     factor();
   }
   return irank();
}

/**
 * Returns log of its determinant.
 */
double SparseMatrix::logdet() {
  if (!logdet_done) {
    logdet(rank);
  }
  return logdeterm;
}

/**
 * Return the rank of the obhect.
 */
unsigned SparseMatrix::irank() {
  if (!logdet_done) {
    logdet(rank);
  }
  return rank;
}

double SparseMatrix::logdet(unsigned &irank)  {
   if (dim == 0) throw exception("SparseMatrix::logdet(): SparseMatrix is empty");
  unsigned i;
  double res=0.0;
   if (!factor_done) {
     if(!factor()) return 0.0;
   }
   irank = 0;
   for (i=1;i<=dim;i++){
     if (diag(i) > 0.0000000001) {
       irank++;
       res += std::log(diag(i));
     }
   }
   logdeterm = 2*res;
   return (logdeterm);
}

/**
 * Add value to element (i,j).
 *
 * Note that Indexing is from 1 instead of 0
 */
int SparseMatrix::insert(const int i,const int j,const double value)
{
  int ii,jj; 
  ii = i;    
  jj = j;    
  
  if (dim == 0) throw exception("SparseMatrix::insert(): SparseMatrix is empty");
  if (i<=0 || j<=0 || j>dim) {
    //std::cout << "insert "<< i<<","<<j;
    throw exception("SparseMatrix::insert(i,j,a): range error"); 
  }
  if (i>j) {
    // symmetry is assumed, and only upper triangular elements saved
    ii = j;
    jj = i;
  }
  //    if (i>j) throw exception("SparseMatrix::insert(i,j,a): i>j, only upper triangular allowed");

  if (fabs(value) <= SESSION.epsilon && ii != jj) return(1);
  int *ia   = iiv-1;
  int *ja   = jjv-1;
  double *a = aav-1;
  unsigned new_ill, ill;
  ill =  1 + (433*ii+53*jj) % hsize;
  for (unsigned k=1; k<=200; k++) {
    if (ia[ill] == 0) {	// first time store
      srt_hash[nonzero]=ill;
      hash_srt[ill]=nonzero;
      nonzero++;  //SDK
      ia[ill]=ii;
      ja[ill]=jj;
      a[ill]=value;
      return(1);
    }
    else if (ia[ill] == ii && ja[ill] == jj) {
      a[ill] += value;
      return(1);
    }
    else {
      new_ill = (ill+640) % hsize + 1;
      ill = (new_ill == ill) ? (ill+1)% hsize + 1 : new_ill;
    }
  }
  throw exception("SparseMatrix::insert(), max_nz is too small");
  return(0);
}

/**
 * Return the element at (i,j) position.
 */
double SparseMatrix::getaij(const int ii,const int jj)
{
  int i,j;
  i=ii;
  j=jj;
  
   if (dim == 0) throw exception("SparseMatrix is empty");
   if (i>j) {
     i=jj;
     j=ii;
   }
//throw exception("SparseMatrix::insert(i,j,a): i>j, only upper triangular allowed");
   if (j>dim) {
     //std::cout << "getaij "<< i<<","<<j;
    throw exception("SparseMatrix::getaij(i,j): range error"); 
  }

   int *ia   = iiv-1;
   int *ja   = jjv-1;
   double *a = aav-1;
   unsigned new_ill, ill;
   ill =  1 + (433*i+53*j) % hsize;
   for (unsigned k=1; k<=200; k++) {
      if (ia[ill] == 0) {	// aij = 0
	return(0.0);
      }
      else if (ia[ill] == i && ja[ill] == j) {
         return(a[ill]);
      }
      else {
         new_ill = (ill+640) % hsize + 1;
         ill = (new_ill == ill) ? (ill+1)% hsize + 1 : new_ill;
      }
   }
   throw exception("SparseMatrix::getaij This should not happen! Probably a bug :-( ");
   return(0.0);
}

/**
 * Reset the element (i,j) to the value.
 */
int SparseMatrix::reset(const int ii,const int jj,const double value)
{
  int i=ii;
  int j=jj;

   if ( dim == 0 ) throw exception("SparseMatrix::reset(): SparseMatrix is empty");
   if (i<=0 || j<=0 ) throw exception("SparseMatrix::reset(i,j,a), range error");
   if ( i>j){
     i=jj;
     j=ii;
     //throw exception("SparseMatrix::reset(i,j,a): i<j only upper triangular allowed");
   }
   
   if (fabs(value) <= SESSION.epsilon && i != j) return(1);
   int *ia   = iiv-1;
   int *ja   = jjv-1;
   double *a = aav-1;
   unsigned new_ill, ill;
   ill =  1 + (433*i+53*j) % hsize;
   for (unsigned k=1; k<=200; k++) {
      if (ia[ill] == 0) {	// first time store
         ia[ill] = i;
         ja[ill] = j;
         a[ill] = value;
         return(1);
      }
      else if (ia[ill] == i && ja[ill] == j) {
         a[ill] = value;      // the previous value has been reset
         return(1);
      }
      else {
         new_ill = (ill+640) % hsize + 1;
         ill = (new_ill == ill) ? (ill+1)% hsize + 1 : new_ill;
      }
   }
   throw exception("SparseMatrix::insert(): max_nz too small");
   return(0);
}

/**
 * After close(), no more element can be added.
 */
unsigned SparseMatrix::close(void)
{                      /* replace iiv by starting elements of rows:*/
   if (insertend) return nonzero;
   unsigned nrow   = 0;		
   unsigned oldrow = 0;
   int *ia   = iiv-1;
   int *ja   = jjv-1;
   double *a = aav-1;
   nonzero = squeeze(hsize,ia,ja,a);      /* squeeze out zeros */
   hsort_ija(nonzero,ia,ja,a);	          /* sort iiv,jjv      */
   unsigned k;
   for (k=1; k<=nonzero ; k++) {
      if (ia[k] != oldrow) {
         nrow++;
         ia[nrow] = k;
         oldrow = ia[k];
      }
   }
   if (nrow != dim) {
      throw exception("SparseMatrix::close(): actual size <. expectec size");
      // filling diagonals with zeros is one way to avoid this ERROR
   }
   unsigned i;
   ia[dim+1] = nonzero+1;
   int zero_pivot = 0;          // zero diagonals ?
   for (k=1; k<=dim; k++) {
      i = ia[k];
      while ( i<ia[k+1] && ja[i] != k ) i++;
      if (i < ia[k+1]  && a[i] == 0.0) {
         a[i] = 1.0;
         zero_pivot++;
      }
   }
   if (zero_pivot) {
      warning("zero-diagonals is found in a SparseMatrix, 1.0 has been added");
   }
   insertend = 1;
   //   std::cout << "\n--> A Sparse Matrix has been Set up" << std::endl;
   return nonzero;
}
// BRS added const as suggested by Rohan

/**
 * Convert a sub sparse matrix(r1:r2,c1:c2) into Matrix object.
 */
doubleMatrix SparseMatrix::dense(unsigned r1, unsigned r2, unsigned c1, unsigned c2)  const
{
  doubleMatrix out(r2-r1+1,c2-c1+1);
  int   *ia = iiv-1;
  int   *ja = jjv-1;
  double *a = aav-1;

  unsigned i,j,k;

  for (k=1; k<= hsize; k++) {
    if (i = ia[k]) {
        j = ja[k];
	if (i!=j) {
	  if (i>=r1 && i<=r2 && j>=c1 && j<=c2 ) {
	    out(i-r1+1,j-c1+1) = a[k];
	  }
          if (j>=r1 && j<=r2 && i>=c1 && i<=c2 ) {
	    out(j-r1+1,i-c1+1) = a[k];
	  }
	}
	else {
	  if (i>=r1 && i<=r2) {
	   out(i-r1+1,i-r1+1) = a[k];
	  }
	}
    }
  }
  return out;
}

/**
 * Display the content.
 */
void SparseMatrix::display(unsigned r1, unsigned r2, unsigned c1, unsigned c2, const std::string &msg)
{
  if (msg != "") std::cout <<"\n     " << msg << "\n";
  if (dim == 0) {
    std::cout << "      null sparse matrix     \n\n";
    return;
  }
  std::cout << dense(r1,r2,c1,c2);
}

/**
 * Save the contents into a disk file.
 */
void SparseMatrix::save(const std::string &fname, const int io_mode) const
{

   std::ofstream ofs;
   ofs.open(fname.c_str(),(OpenModeType)io_mode);
   if (!ofs) throw exception(std::string("SparseMatrix::save(): cannot open file:") + fname);

   ofs.precision(16);
   ofs.setf(std::ios::unitbuf);
   int i,j,ii,jj,jjj,last_j;
   for (i=0; i<dim; i++) {
      ii = i + 1; last_j = iiv[ii] - 1;
      for (j = iiv[i]; j<=last_j; j++) {
         jjj = j-1;
         jj = jjv[jjj];
         ofs << ii << " " << jj << " " << aav[jjj] << "\n";
         if (jj != ii) ofs << jj << " " << ii << " " << aav[jjj] << "\n";

      }
   }
   ofs.close();
   ofs.unsetf(std::ios::unitbuf);
}

void SparseMatrix::gibbs_iterate(double *x,const double *r,double *w,
                                 const double se)
{
   /////////////////////////////////////////////////////////////////
   // REQUIREMENTS: x,r, and w must have enough spaces
   // se = sqrt(residual variance), if MMe is built with variance
   // ratio, then se is necessary, otherwise set 1 which is default
   //////////////////////////////////////////////////////////////////
   double *sample = x-1;
   double *adj_rhs = w-1;
   double *a = aav-1;
   int *ia = iiv-1;
   int *ja = jjv-1;
   int i,j,jj,k,ke;
   double mean,diagon=0.0;

   memcpy(w,r,sizeof(double)*dim);
   for (i=1; i<=dim; i++) {
      mean = adj_rhs[i];
      k = ia[i];  ke = ia[i+1];
      for (j=k; j<ke; j++) {
         jj = ja[j];
         if (jj == i) {diagon = a[j];}  // diagonal may not be the first elements
         else         {mean -= a[j]*sample[jj];}
      }
      mean /= diagon;
//      mean += se*snorm()/std::sqrt(diagon);   // variance = 1/diagon
      sample[i] = mean;                // adjust the rhs, use all coef of row i
      for (j=k; j<ke; j++)  adj_rhs[ja[j]] -= a[j]*mean;
   }
}
/*
//BRS this has error with passing a const and trying to copy from it
ostream& operator<<(ostream& stream, const SparseMatrix& A)
{
  unsigned dim = A.nrow();
   Matrix out;
//   out  = A.dense(1,dim,1,dim);
//  out.display();
  return stream;
}
*/

/**
 * squeeze out zeros from (ia,ja,a).
 *
 * ia(1...len),jjv(1..len),a(1..len)
 *  return the actural number of non-zeros in a
 *  actually this moves the nonzeros to the begining
 */

static int squeeze(unsigned len,int *ia,int *ja,double *a)
{
   unsigned lmar = 1;
   unsigned rmar = len;
   for (;;) {
      if (lmar >= rmar) {
         if (ia[lmar] == 0 && rmar == lmar) return(rmar-1);
         else return(rmar);
      }
      if ( ia[lmar] != 0) {
         lmar++;
      }
      else if (ia[rmar] == 0) {
         rmar--;
      }
      else {
         ia[lmar] = ia[rmar];
         ja[lmar] = ja[rmar];
         a[lmar]  = a[rmar];
         lmar++;
         rmar--;
      }
   }
}

/**
 * this is a Heap sort, sorting (ia,ja) ascendingly.
 *
 * ia(1...n), ja(1...n) ,a(1..n)
 * where n is the actual number of non-zeros in A
 */
static void hsort_ija(unsigned n,int *ia,int *ja,double *a)
{
   unsigned iia, jja,i,j,l,ir;
   double aa;
   if (n == 1) return;
   l = (n >> 1)+1;		/* l = n/2+1; */
   ir=n;
   for (;;) {
      if (l > 1) {
         l = l-1;  iia = ia[l];  jja = ja[l];  aa = a[l];

      }
      else {
         iia = ia[ir];  jja = ja[ir];  aa = a[ir];
         ia[ir] = ia[1]; ja[ir] = ja[1];  a[ir] = a[1];
         ir--;
         if(ir == 1) {
            ia[1] = iia;  ja[1] = jja;  a[1] = aa;
            return;
         }
      }
      i=l;
      j=l+l;
      while ( j <= ir) {
         if (j < ir) {
            if (ia[j]<ia[j+1]  || (ia[j]==ia[j+1] && ja[j] < ja[j+1])) j++;
         }
         if (iia < ia[j] || (iia == ia[j] && jja < ja[j])) {
            ia[i] = ia[j];  ja[i] = ja[j];  a[i] = a[j];
            i=j;
            j=j+j;
         }
         else {
            j=ir+1;
         }
      }
      ia[i] = iia;  ja[i] = jja;  a[i] = aa;
   }
}


int SparseMatrix::initialize_node_list(void) {
  int *ia   = iiv-1;
  int *ja   = jjv-1;
  double *a = aav-1;
  unsigned i,j;
  unsigned n_rows,row,col;
  unsigned count;

  // I am going to count number of nabors for each node
  // I will use "mask" to do the counting!

  mask.resize(dim);

  node_list.reserve(dim);
  // do the counting
  //SDK
  int isrt;
  for(isrt=0;isrt<nonzero;isrt++) {
    i=srt_hash[isrt];
    //  for (i=1;i<=hsize;i++){
    //SDK
    if (!ia[i]) continue; //SDK
    row = ia[i];
    col = ja[i];
    if(row != col){
      mask(row) = mask(row) + 1;
      mask(col) = mask(col) + 1;
    }
  }
    // now resize the adjacent list for each node
  for (i=1;i<=dim;i++){
    node_list(i).adj_list.reserve(mask(i));
  }
  // now I am going to use mask to mark the current position in the adj_list for each node
  // start at position 1
  for (i=1;i<=dim;i++){
    mask(i) = 1;
  }
    // now do the filling
  //SDK
  for(isrt=0;isrt<nonzero;isrt++) {
    i=srt_hash[isrt];
    //  for (i=1;i<=hsize;i++){
    //SDK
    if (!ia[i]) continue; //SDK
    row = ia[i];
    col = ja[i];
    if(row!= col){
      node_list(row).adj_list(mask(row)) = col;
      node_list(col).adj_list(mask(col)) = row;
      mask(row) = mask(row) + 1;
      mask(col) = mask(col) + 1;
    }
  }
  initialize_node_list_done = 1;
  return 1;
}

int SparseMatrix::minfilsymelm_node_list(void)
{
  unsigned i,j;
  unsigned orderj,inor;
  unsigned mindeg=dim,deg,oldmindeg=0;
  Vector<unsigned> same;
  unsigned nsame;
  unsigned minpos,elimpos,temp;
  unsigned nabor;
  unsigned count;
  unsigned *num;
  //std::cout << "entering minfilsymelm_node_list" << std::endl;


  if (!initialize_node_list_done) {
    if (!initialize_node_list()) return 0;
  }

  order.reserve(dim);
  inv_order.reserve(dim);
  same.resize(dim);
  // initialize order
  for (i=1;i<=dim;i++) {
    order(i) = i;
  }


  if(xadj) delete [] xadj;
  if(adjncy) delete [] adjncy;
  xadj=new idxtype [dim+1];
  adjncy=new idxtype [2*(nonzero-dim)];
  xadj[0]=1;
  idxtype *xoff,*aoff;
  xoff=xadj;xoff--;
  aoff=adjncy;aoff--;
  idxtype idcnt=1;
  for(i=1;i<=dim;i++) {
      deg = node_list(i).adj_list.size();
      xoff[i+1]=xoff[i]+deg; 
      for(j=1;j<=deg;j++,idcnt++) {
	aoff[idcnt]= node_list(i).adj_list(j);
      } 
  }
  masklist.reserve(dim);
  mask.resize(dim);
  
#ifdef DO_MINORD
  // loop for elimination
  Vector<unsigned> degpt,degnxt,deglst;
  degpt.resize(dim+1);
  degnxt.resize(dim+1);
  deglst.resize(dim+1);
  unsigned tmp1,tmp2,tmp3,k;
  for(i=1;i<=dim;i++){
    deg=node_list(i).adj_list.size();
    tmp1=degpt[deg];
    degpt[deg]=i;
    degnxt[i]=tmp1;
    if(tmp1) deglst[tmp1]=i;
  }
  //  cout << "degpt\n" << degpt; 
  //  cout << "degnxt\n" << degnxt; 
  //  cout << "deglst\n" << deglst; 
  mindeg=1;
  for (i=1;i<=dim;i++) {
    // loop to look through all remaining nodes to update mindeg
    //mindeg = dim+1;
    nsame=0;
    if(mindeg) mindeg--;
    for(k=mindeg;degpt[k]==0&&k<=dim;k++);
    minpos=degpt[k];
    //if(k<mindeg) cout << "Mindeg " << mindeg <<" " << k << std::endl;
    mindeg=k;
    /*
    for (j=i;j<=dim;j++) {
      orderj = order(j);
      deg = node_list(orderj).adj_list.size();
      if (deg < mindeg) {
	mindeg = deg;
	minpos = j;
      }
    }
    */
    
    //    cout << "Mindeg " << mindeg << std::endl;
    // update order vector
    for(j=0;j<mindeg;j++) {
      //orderj = order(i);
      k=node_list(minpos).adj_list[j];
      same[j]=node_list(k).adj_list.size();
    }
    temp = order(i);
    order(i)=minpos;
    //order(i) = order(minpos);
    //order(minpos) = temp;
    // now eliminate node at minpos
    //    std::cout << i <<" eliminating row at " << minpos<< " " << std::endl;
    if (!elim_node(i,order)) return 0;

    {
      k=minpos;
      unsigned a,b,c;
      a=deglst[k];
      c=degnxt[k];
      deglst[c]=0;
      degpt[mindeg]=c;
      if(a !=0) std::cout << "Err " << a << " " << k << std::endl;

    }
    for(j=0;j<mindeg;j++) {
      k=node_list(minpos).adj_list[j];
      deg=node_list(k).adj_list.size();
      if(same[j]!=deg){
	unsigned a,b,c;
	b=degpt[deg];
	c=degnxt[k];
        a=deglst[k];
	//if(a == minpos)a=0;
	degpt[deg]=k;
	degnxt[k]=b;
	deglst[k]=0;
	if(b) deglst[b]=k;
	if(a) degnxt[a]=c;
	if(c) deglst[c]=a;
	if(a==0) degpt[same[j]]=c;

      }
    }
    //  cout << "degpt\n" << degpt; 
    //  cout << "degnxt\n" << degnxt; 
    //  cout << "deglst\n" << deglst; 
  }
#endif
  for (i = 1;i<=dim;i++) {
    // loop to get inv_order from order
#ifndef DO_MINORD
    if (!elim_node(i,order)) return 0;
#endif
    inor = order(i);
    inv_order(inor) = i;
  }
  
  //display_order();
  //display_inv_order();
  // renumber and rearrange adj_lists according to order of elimination
  for (i=1;i<=dim;i++) {
    deg = node_list(i).adj_list.size();
    // put nabors into mask according to elimination order
    for (j=1;j<=deg;j++){
      nabor = node_list(i).adj_list(j);
      elimpos = inv_order(nabor);
      node_list(i).adj_list(j) = elimpos;
    }
    num = node_list(i).adj_list.begin();
    qsort(num,deg,sizeof(unsigned), (cmp_func) comp);
  }
  // get memory for lij and diagonals
  unsigned totsize=0;
  for (i=1;i<=dim;i++) {
    deg = node_list(i).adj_list.size();	
    node_list(i).a.reserve(deg);
    totsize+=deg;
    //    std::cout << " i " << i << ":";
    //    num = node_list(i).adj_list.begin();
    //    for(j=0;j<deg;j++) std::cout << " " << num[j];
    //    std::cout << "\n";
    if (node_list(i).a.size() < 0) return 0;
  }
  //  std::cout << "Total : " << totsize << "\n";
  diag.reserve(dim);
  elim_done = 1;
  //std::cout << "exiting minfilsymelm_node_list" << std::endl;
  return 1;
}

#ifndef NEW_ELIMNODE
int SparseMatrix::elim_node(unsigned i, Vector<unsigned> Aorder)
{

  // updates elimination graph for eliminating elimnode

  unsigned j,k;
  unsigned naborj,nabork,nabor_nabor;
  unsigned degnode,degnabor;
  unsigned count;
  unsigned *maskptr,nmask;
  Vector<unsigned> masklist(dim);
  //Array<unsigned> mask_base(dim);
  
  NodeStruct *elimnode, *nabornode;
  elimnode = &node_list(Aorder(i));
  maskptr=mask.begin();
  // loop to update adjacent lists for neigbors of elimnode
  degnode = elimnode->adj_list.size();
  memset(mask.begin(),'\0',dim*sizeof(unsigned)); // SDK
  for(j=1;j<=degnode;j++) {
    naborj = elimnode->adj_list(j);
    //  std::cout << " i " << i << " naborj " << naborj << " Aorder  " << Aorder(i) << " Aorder j " << Aorder(j) << "\n";
    mask[naborj-1]=1;
    masklist[j-1]=naborj-1;
  }
  
  for (j=1;j<=degnode;j++) {
    naborj = elimnode->adj_list(j);
    nabornode = &node_list(naborj);
    nmask=degnode;
    // initialize and make mask for adjcent list of naborj
    //for (k=0;k<dim;k++){ // SDK
    //  mask[k] = 0;       // SDK
    //}                    // SDK
    //nmask=0;                                // SDK
    // mark naborj's nabors in mask (except elimnode)
    degnabor = nabornode->adj_list.size();
    for (k=0;k<degnabor;k++){
      nabor_nabor = nabornode->adj_list[k];
            if (nabor_nabor != Aorder(i)) { //elim elimnode
	//if (inv_order(nabor_nabor) >= i) { //elim elimnode BAD IDEA
	//	mask[nabor_nabor-1] = 1; //SDK
	if(!mask[nabor_nabor-1]) {  //SDK
	  // mask[nabor_nabor-1]=1;    //SDK
          masklist[nmask++]=nabor_nabor-1; //SDK
	}                      //SDK
      }
    }
    // mark elimnodes nabors in mask (except naborj)
    //    for (k=0;k<degnode;k++) {
    //      nabork = elimnode->adj_list[k];
    //      if (nabork != naborj) {
    //	        if(!mask[nabork-1]) {  //SDK as above
    //           mask[nabork-1]=1;
    //           masklist[nmask++]=nabork-1;
    //          }
    //       }
    //    }
    // transfer nabors from mask into adjcent list for naborj
    // first get space
    //count = 0;               //SDK 
    //for (k=0;k<dim;k++) {    //SDK
    //  if (mask[k]) count++;  //SDK
    //}                        //SDK
    count=nmask-1; //SDK
    if (nabornode->adj_list.size() != count) {
      nabornode->adj_list.reserve(count);
      if (nabornode->adj_list.size() < count) return 0;
    }
    // now do the transfer
    //count = 0;                           //SDK
    //for (k=0;k<dim;k++) {                //SDK 
    //  if (mask[k]) {                     //SDK
    //	nabornode->adj_list[count] = k+1;  //SDK 
    //	count++;                           //SDK
    //  }                                  //SDK
    //}    
                                //SDK
    int adj_count=0;
    for(count=0;count<nmask;count++){                  //SDK
      if(masklist[count] != (naborj-1)) {
	nabornode->adj_list[adj_count++] = masklist[count]+1;  //SDK
      }
    }
    if(count != (adj_count+1)) {
      std::cout << count << " " << adj_count << "\n" << std::flush;
      throw exception("Elim node: not possible");
      }//SDK
  }
  return 1;
}
#endif


int SparseMatrix::initial_lij(void)
{
  if (insertend) throw exception("SparseMatrix::initial_lij: cannot initial_lij after closing");
  unsigned i,j;
  unsigned row,col;
  unsigned deg,nabor;
  double aij;
  double *paij;
  unsigned *nabors;
  NodeStruct *nodes, *node;

  if (!elim_done) {
    if (!minfilsymelm_node_list()) return 0;
  }

  nodes = node_list.begin();
  for (i=0;i<dim;i++) {
    row = i+1;
    diag[i] = getaij(row,row);
    node = &nodes[i];
    nabors = node->adj_list.begin();
    deg    = node->adj_list.size();

    for (j=0;j<deg;j++){
      nabor = nabors[j];
      col = order(nabor);
      if (col < row) {
	aij = getaij(col,row);
      }
      else{
	aij = getaij(row,col);
      }
      paij = node->a.begin();
      paij[j] = aij;
    }
  }
  initial_lij_done = 1;
  return 1;
}


int SparseMatrix::factor(void){
  unsigned i,j,nextj,k,count,m;
  unsigned coli,colj,rowm,rowk;
  unsigned degi,degj;
  unsigned startj;
  Vector<double> temp;
  Vector<unsigned> start, contr_cols;
  double lij,diagi;
  unsigned *prows;
  double   *plij;
  NodeStruct *nodes, *node;

  if (!initial_lij_done) {
    if(!initial_lij()) return 0;
  }
  temp.reserve(dim);
  start.reserve(dim);
  contr_cols.reserve(dim);


  double   *ptemp        = temp.begin();       ptemp--;
  double   *pdiag        = diag.begin();       pdiag--;
  unsigned *pstart       = start.begin();      pstart--;
  unsigned *porder       = order.begin();      porder--;
  unsigned *pcontri_cols = contr_cols.begin(); pcontri_cols--;
  double    mult;
  unsigned lastpos;
  lastpos=0;

  info_vec[0]=0;
  info_vec[1]=0;
  //for (i=1;i<=dim;i++) {
  //  ptemp[i]        = 0.0;
  //  pstart[i]       = 0;
  //  pcontri_cols[i] = 0;
  // }
  memset(contr_cols.begin(),'\0',dim*sizeof(unsigned));
  memset(start.begin(),'\0',dim*sizeof(unsigned));
  memset(temp.begin(),'\0',dim*sizeof(double));
  // loop for calculating column i of lij
  nodes = node_list.begin(); nodes--;
  for (i=1;i<=dim;i++){
    coli = porder[i];
    j = pcontri_cols[i];
    count = 0;
    // loop for getting contributions from columns contributing to column i
    // contri_cols is used to point to these columns
    while (j) {
      startj = pstart[j];
      colj  = porder[j];
      node  = &nodes[colj];
      degj  = node->adj_list.size();
      nextj = pcontri_cols[j];
      prows = node->adj_list.begin();
      plij  = node->a.begin();
      lij = plij[startj];
      // update pstart[j] for next time
      pstart[j] = startj+1;
      ptemp[prows[startj++]]+=lij*lij;
      // accumulate contributions from col j to col i in temp
      for (k=startj;k <degj ;k++) {
	rowk=prows[k];
	ptemp[rowk]  =ptemp[rowk]+ lij*plij[k] ;
      }
    
      
      k = startj;
      if (k<degj) {
	rowk = prows[k];
	pcontri_cols[j] = pcontri_cols[rowk];
	pcontri_cols[rowk] = j;
       }
      j = nextj;
    }
    
    node  = &nodes[coli];
    prows = node->adj_list.begin();
    plij  = node->a.begin();
    diagi = pdiag[coli] - ptemp[i];
    degi = node->adj_list.size();
    //         std::cout << "No of contributions to row: " << count << std::endl;
    //         std::cout << "No of adj elements       : " << degi   << std::endl;
    //         std::cout << "diagonal element: " << coli << " with value: " << pdiag[coli] << std::endl;
    //         std::cout << "row             : " << i    << " with value: " << ptemp[i] << std::endl;
    //         std::cout << "check psd:" << diagi << std::endl;
    //         std::cout << std::endl;
    if (diagi > std::max(1.e-8,1.e-10*pdiag[coli])) { 
      info_vec[0]+=std::log(diagi);
      info_vec[1]++;
      //std::cout << " rank " << info_vec[1] ;
      diagi = std::sqrt(diagi);
      pdiag[coli] = diagi;
      mult = 1/diagi;
      // use contributions in temp to update col i
      for (k=0;k<degi;k++){
	rowk = prows[k];
	plij[k]-=ptemp[rowk];
	plij[k]*=mult;
	ptemp[rowk] = 0.0;
      }
      
      if (degi) {
	rowk = prows[0];
	pcontri_cols[i] = pcontri_cols[rowk];
	pcontri_cols[rowk] = i;
      }
      
    }
    else if (diagi > -std::max(1.e-10,1.e-10*pdiag[coli])) {
      pdiag[coli] = 0;
      degi = node->adj_list.size();
      for (k=0;k<degi;k++){
	rowk = prows[k];
        plij[k] = 0.0;
	ptemp[rowk] = 0.0;
      }
    }
    else if (diagi < 0) {
      throw exception("SparseMatrix::factor matrix is not positive semi-definite");
    }
  }

  factor_done = 1;
  return 1;
}


Vector<double> SparseMatrix::solv(Vector<double> rhs)
{
  if (!factor_done) throw exception("SparseMatrix::factor must be called first");
  unsigned i,j;
  unsigned coli,rowj,corowj;
  unsigned degi;
  double diagi,soli;
  Vector<double> sol,fsol;
  unsigned *prows;
  double *pa;
  NodeStruct *nodes, *node;

  sol.reserve(dim);
  fsol.reserve(dim);

  double   *pdiag  = diag.begin();  pdiag--;
  double   *psol   = sol.begin();    psol--;
  double   *pfsol  = fsol.begin();   pfsol--;
  double   *prhs   = rhs.begin();    prhs--;
  unsigned *porder = order.begin(); porder--;


  // forward elimination
  // loop for getting solution 1..dim
  nodes = node_list.begin(); nodes--;
  for (i=1;i<=dim;i++) {
    coli  = porder[i];
    diagi = pdiag[coli];
    node  = &nodes[coli];
    prows = node->adj_list.begin();
    pa    = node->a.begin();

    if (diagi == 0) { //lineraly dependent equation
      psol[coli] = 0.0;
    }
    else {
      soli = psol[coli] = prhs[coli]/diagi;
      // adjust rhs for this column
      // note: that diag and rhs are in original order
      //       but, the row subscripts of L are  in mindeg order


      degi = node->adj_list.size();
      for (j=0;j<degi;j++) {
	rowj = porder[prows[j]];
	prhs[rowj] -= pa[j]*soli;
      }
    }
  }
  // backward substitution
  // loop for getting solution dim..1
  for (i=dim;i>=1;i--) {
    coli = porder[i];
    diagi = pdiag[coli];
    node  = &nodes[coli];
    prows = node->adj_list.begin();
    pa    = node->a.begin();
    if (diagi == 0) { //lineraly dependent equation
      pfsol[coli] = 0.0;
    }
    else {
      // adjust this rhs for all previous solutions
      degi = node->adj_list.size();  
      soli=psol[coli];
      for (j=0;j<degi;j++) {
	rowj = porder[prows[j]];
	soli -=  pa[j]*pfsol[rowj];
      }
      pfsol[coli] = soli/diagi;
    }
  }
  return fsol;
}

int SparseMatrix::solv(double* x, const double *rhs)
{
  if (!factor_done) throw exception("SparseMatrix::factor must be called first");
  unsigned i,j;
  unsigned coli,rowj;
  unsigned degi;
  double diagi,soli;
  Vector<double> sol,temp_rhs;
  unsigned *prows;
  double *pa;
  NodeStruct *nodes, *node;

  sol.reserve(dim);
  temp_rhs.reserve(dim);
  memcpy(temp_rhs.begin(),rhs,sizeof(double)*dim);

  double   *pdiag  = diag.begin();     pdiag--;
  double   *psol   = sol.begin();       psol--;
  unsigned *porder = order.begin();    porder--;
  double   *prhs   = temp_rhs.begin();  prhs--;
  double   *pfsol; 
  pfsol = x; pfsol--;
 

  
  // forward elimination
  // loop for getting solution 1..dim
  nodes = node_list.begin(); nodes--;
  for (i=1;i<=dim;i++) {
    coli  = porder[i];
    diagi = pdiag[coli];
    node  = &nodes[coli];
    prows = node->adj_list.begin();
    pa    = node->a.begin();

    if (diagi == 0) { //lineraly dependent equation
      psol[coli] = 0.0;
    }
    else {
      soli = psol[coli] = prhs[coli]/diagi;
      // adjust rhs for this column
      // note: that diag and rhs are in original order
      //       but, the row subscripts of L are  in mindeg order


      degi = node->adj_list.size();
      for (j=0;j<degi;j++) {
	rowj = porder[prows[j]];
	prhs[rowj] = prhs[rowj] - pa[j]*soli;
      }
    }
  }
  // backward substitution
  // loop for getting solution dim..1
  for (i=dim;i>=1;i--) {
    coli = porder[i];
    diagi = pdiag[coli];
    node  = &nodes[coli];
    prows = node->adj_list.begin();
    pa    = node->a.begin();
    if (diagi == 0) { //lineraly dependent equation
      pfsol[coli] = 0.0;
    }
    else {
      // adjust this rhs for all previous solutions
      degi = node->adj_list.size();
      for (j=0;j<degi;j++) {
	rowj = porder[prows[j]];
	psol[coli] = psol[coli] - pa[j]*pfsol[rowj];
      }
      pfsol[coli] = psol[coli]/diagi;
    }
  }
  return 1;
}

void SparseMatrix::display_node_list(void) {
  // displays the nodes in the adjacent list
  // using peeling order numbers
  unsigned i,j;
  unsigned count;

  for (i=1;i<=dim;i++) {
    int orderi = order(i);
    std::cout <<"\n" <<"node: " << i <<": ";
    
    count = node_list(orderi).adj_list.size();
    for (j=1;j<=count;j++){
      std::cout << node_list(orderi).adj_list(j)<<" ";
    }
  }
  std::cout << std::endl;

//   unsigned i,j;
//   unsigned count;
//   for (i=1;i<=dim;i++) {
//     std::cout <<"\n" <<"node: "<<i<<" ";
//     count = node_list(i).adj_list.size();
//     for (j=1;j<=count;j++){
//       std::cout << node_list(i).adj_list(j) <<" ";
//     }
//   }
//   std::cout << std::endl;



}
void SparseMatrix::display_lij(void) {
  // displays the lij's
  unsigned i,j;
  unsigned count;

  for (i=1;i<=dim;i++) {
    std::cout <<"\n" <<"node: "<<i<<" ";
    count = node_list(i).adj_list.size();
    for (j=1;j<=count;j++){
      std::cout << node_list(i).a(j) <<" ";
    }
  }
  std::cout << std::endl;
}
void SparseMatrix::display_order(void) {
  // displays the order of diagonal nodes
  unsigned i;
  std::cout << "\n" << "order" << "\n";
  for (i=1;i<=dim;i++) {
    std::cout << order(i) << std::endl;
  }
  std::cout << std::endl;
}

void SparseMatrix::display_inv_order(void) {
  // displays the order of diagonal nodes
  unsigned i;
  std::cout << "\n" << "inv_order" << "\n";
  for (i=1;i<=dim;i++) {
    std::cout << inv_order(i) << std::endl;
  }
  std::cout << std::endl;
}




int SparseMatrix::solvrow(double* x, const double *rhs,unsigned row,Vector<double> &sol,Vector<double> &temp_rhs,unsigned *row_pt,unsigned *col_pt) {
  unsigned i,j,allzero;
  unsigned coli,rowj,corowj;
  unsigned degi;
  double diagi,soli;
  //Vector<double> sol,temp_rhs;
  unsigned *prows;
  double *pa;
  NodeStruct *nodes, *node;
  
  
  if (!factor_done) {
    throw exception("SparseMatrix::factor must be called first");
    return 0;
  }
  //  sol.reserve(dim);
  //  temp_rhs.reserve(dim);
  // unsigned *row_pt=new unsigned [dim];
  unsigned *row_mask=row_pt;
  //  unsigned *col_pt=new unsigned [dim];
  unsigned *col_mask=col_pt;
  row_mask--;
  memset(row_pt,'\0',dim*sizeof(unsigned));
  col_mask--;
  memset(col_pt,'\0',dim*sizeof(unsigned));
  allzero=1;
  double   *pdiag  = diag.begin();     pdiag--;
  double   *psol   = sol.begin();       psol--;
  unsigned *porder = order.begin();    porder--;
  unsigned *pinv_order = inv_order.begin();    pinv_order--;
  double   *prhs   = temp_rhs.begin();  prhs--;
  double   *pfsol; 
  pfsol = x; pfsol--;
  memset(temp_rhs.begin(),'\0',dim*sizeof(double));
  //  prhs[porder[row]]=1;
  prhs[row]=1;
  //  memcpy(temp_rhs.begin(),rhs,dim*sizeof(double));
  memset(x,'\0',dim*sizeof(double));

  
  // forward elimination
  // loop for getting solution 1..dim
  nodes = node_list.begin(); nodes--;
  row_mask[row]=1;
  col_mask[row]=1;
  for (i=row;i<=dim;i++) {
    if(row_mask[i]) {
      coli  = porder[i];
      diagi = pdiag[coli];
      node  = &nodes[coli];
      prows = node->adj_list.begin();
      pa    = node->a.begin(); 
      //row_mask[i]=1;
      if (diagi != 0)  {
	//	soli = psol[coli] = prhs[coli]/diagi;
	soli = psol[i] = prhs[i]/diagi;
	row_mask[i]=1;
	// adjust rhs for this column
	// note: that diag and rhs are in original order
	//       but, the row subscripts of L are  in mindeg order
	
	
	degi = node->adj_list.size();
	for (j=0;j<degi;j++) {
	  row_mask[prows[j]]=1;
	  //	  if(prows[j] == row) col_mask[i]=1;
	  //	  rowj = porder[prows[j]];
	  //	  prhs[rowj] -=  pa[j]*soli;
	  prhs[prows[j]] -=  pa[j]*soli;
	}
      }
    }
  }
  // backward substitution
  // loop for getting solution dim..1
  for (i=dim;i>=row;i--) {
    if(row_mask[i]) {
      coli = porder[i];
      diagi = pdiag[coli];
      if (diagi != 0 )  {
	node  = &nodes[coli];
	prows = node->adj_list.begin();
	pa    = node->a.begin();
	// adjust this rhs for all previous solutions
	degi = node->adj_list.size();
	//	soli=psol[coli];
#pragma ivdep
	for (j=0;j<degi;j++) {
	  //	  rowj = porder[prows[j]];
	  //	  psol[coli] -= pa[j]*pfsol[rowj];
	  psol[i] -= pa[j]*pfsol[prows[j]];
	}
	//	pfsol[coli] = psol[coli]/diagi;
	pfsol[i] = psol[i]/diagi;
      }
    }
  }
  //  for(i=1;i<=dim;i++){
  //    psol[porder[i]]=pfsol[i];
  //  }
  //  delete [] row_pt;
  //  delete [] col_pt;
  return 1;
}

#ifdef NEW_ELIMNODE

int SparseMatrix::elim_node(unsigned i, Vector<unsigned> Aorder)
{
  
  // updates elimination graph for eliminating elimnode
  
  unsigned j,k;
  unsigned naborj,nabork,nabor_nabor;
  unsigned degnode,degnabor;
  unsigned count,*adjpt,*nadjpt,*aadjpt;
  unsigned *maskptr,nmask;
  //  Vector<unsigned> masklist(dim);
  NodeStruct *elimnode, *nabornode;
  elimnode = &node_list(Aorder(i));
  maskptr=mask.begin();
  // loop to update adjacent lists for neigbors of elimnode
  degnode = elimnode->adj_list.size();
  //  memset(mask.begin(),'\0',dim*sizeof(unsigned)); // SDK


 
  adjpt=elimnode->adj_list.begin();
  for(k=1;k<=degnode;k++) {
    naborj = *adjpt++;
    mask[naborj-1]=1;
    masklist[k-1]=naborj-1;
  }
 
  adjpt=elimnode->adj_list.begin();
  for (j=1;j<=degnode;j++) { 
 
    naborj = *adjpt++;
    

    nabornode = &node_list(naborj);
    nmask=degnode;
    degnabor = nabornode->adj_list.size();
    nadjpt = nabornode->adj_list.begin();
    for (k=0;k<degnabor;k++){
      nabor_nabor = *nadjpt++;
      if (nabor_nabor != Aorder(i)) { //elim elimnode
	if(!mask[nabor_nabor-1]) {  //SDK
	  masklist[nmask++]=nabor_nabor-1; //SDK
	}                      //SDK
      }
    }
    count=nmask-1; //SDK
    if (nabornode->adj_list.size() != count) {
      nabornode->adj_list.reserve(count);
      if (nabornode->adj_list.size() < count) return 0;
    } 
    int adj_count=0;
    nadjpt=nabornode->adj_list.begin();
    for(count=0;count<nmask;count++){                  //SDK
      if(masklist[count] != (naborj-1)) {
	adj_count++;
	*nadjpt++ = masklist[count]+1;  //SDK
      }
    }
    if(count != (adj_count+1)) {
      std::cout << count << " " << adj_count << "\n" << std::flush;
      throw exception("Elim node: not possible");
    }//SDK
  }
  
  adjpt=elimnode->adj_list.begin();
  for(k=1;k<=degnode;k++) {
    naborj = *adjpt++;
    mask[naborj-1]=0;
  }
  return 1;
}
#endif






}






