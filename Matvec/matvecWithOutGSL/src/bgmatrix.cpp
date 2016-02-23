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

#include <iostream>
#include "session.h"
#include "util.h"
#include "bgmatrix.h"

namespace matvec {

extern void ludcmp(BG **a,int n,int *indx,int& d, double tol);
extern void lubksb(BG **a,int n,int *indx,BG *b);
extern bool psdefinite(const BG **mat,const unsigned m,const unsigned n,const double tol);


BGMatrix::BGMatrix(const doubleMatrix &a)
{
   resize(a.num_rows(),a.num_cols());
   for (int j,i=0; i<a.num_rows(); ++i) {
     for (j=0; j<a.num_cols(); ++j) {
       Matrix<BG>::me[i][j] = a[i][j];
     }
   }
}

BGMatrix&  BGMatrix::operator = (const Matrix<bool> &a)
{
   resize(a.num_rows(),a.num_cols());
   for (int j,i=0; i<a.num_rows(); ++i) for (j=0; j<a.num_cols(); ++j) Matrix<BG>::me[i][j] = (BG)a[i][j];
   return *this;
}


bool BGMatrix::psd(void) const
{
   if (!this->symmetric()) return false;
   return psdefinite(begin(),Matrix<BG>::nrow,Matrix<BG>::ncol,SESSION.epsilon);
}



BGMatrix& BGMatrix::inv(void)
{
   //////////////////////////////////////////////////////////
   // A itself turns into its inverse
   // returns non_singular indicator
   // REQUIREMENTS:  A must be nonsingular
   /////////////////////////////////////////////////////////
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;

   int i,j,d;
   if ( nrow != ncol ) throw exception("BGMatrix::inv(): matrix must be square");
   Vector<int>    indx;  indx.reserve(nrow);
   Vector<BG> tempcol; tempcol.reserve(nrow);
   Matrix<BG> temp; temp.reserve(nrow,nrow);

   ludcmp(me,nrow,indx.begin(),d, SESSION.epsilon);
   for (j=0; j<nrow; j++) {
      tempcol.assign(0.0);
      tempcol[j] = 1.0;
      lubksb(me,nrow,indx.begin(),tempcol.begin());
      for (i=0; i<nrow;i++) temp[i][j] = tempcol[i];
   }
   //for (i=0;i<nrow; i++) memcpy(me[i],temp[i],sizeof(BG)*nrow);
   for (j=0; j<nrow; j++) { 
     for (i=0; i<nrow;i++) me[i][j] = temp[i][j];
   }
   return *this;
}

BGMatrix& BGMatrix::ginv1(unsigned *irank)
{         // it only works for symmetric-positive-semidefinite matrix
          //  Matrix itself becomes its inverse;
   if (!this->symmetric()) {
     throw exception("BGMatrix::ginv1(): matrix must be symmetric");
   }
   BG lgdet;
   unsigned myrank = ginverse1(begin(),num_rows(),lgdet,1,SESSION.epsilon);
   if (irank) *irank = myrank;
   return *this;
}


BGMatrix BGMatrix::lu_solve(const BGMatrix& rhs)
{
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;

   if (nrow  != ncol) throw exception("BGMatrix::lu_solve(): matrix is not square ");
   int m = rhs.num_rows();
   int n = rhs.num_cols();

   if (ncol != m ) throw exception("BGMatrix::lu_solve(): bad arg");
   BGMatrix solmat(m,n);
   int i,j,d;
   Vector<BG> solvec; solvec.reserve(m);
   Vector<int> indx;  indx.reserve(m);
   ludcmp(me,nrow,indx.begin(),d,SESSION.epsilon);
   warning("A.solve(b), A destroyed");
   for (j=0; j<n; j++) {
      for (i=0; i<m; i++) solvec[i] = rhs[i][j];
      lubksb(me,m,indx.begin(),solvec.begin());
      for (i=0; i<m; i++) solmat[i][j] = solvec[i];
   }
   return solmat;
}

Vector<BG> BGMatrix::lu_solve(const Vector<BG>& rhs)
{
   int n = rhs.size();
   BGMatrix bmat(n,1);
   for (int i=0; i<n; i++) bmat[i][0] = rhs[i];
   bmat = lu_solve(bmat);
   return bmat.vec();
}


void BGMatrix::gs_solve(const BGMatrix& rhs,BGMatrix& solmat,const double relax,
                     const double stopval,const int mxiter)
{
  //   relax    relaxation coefficient, relax=1.2 seems good for some case
  //   stopval  the accuracy at which iteration stops.
  //   mxiter   maximum number of iterations allowed

   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;

   int m = rhs.num_rows();
   int n = rhs.num_cols();
   if (nrow != ncol) throw exception("BGMatrix::gs_solve(): matrix must be square ");
   if (ncol != m ) throw exception("BGMatrix::gs_solve(): matrix and rhs are not conformable");
   if (solmat.num_rows() != m || solmat.num_cols() != n) {
      solmat.resize(m,n);
   }
   int i,j,k;
   BG diag;
   int niter = 0;
   BG a,cmax;
   BG tol = SESSION.epsilon;

   Vector<BG> oldsol(n);
   Vector<BG> newsol(n);
   Vector<BG> cval(n);
   Vector<BG> local(n);
   do {                             // now iteration begins
      cmax = 0.0;
      for (i=0; i<m; i++) {
         for (k=0; k<n; k++) {
            oldsol[k] = solmat[i][k];
            local[k] = 0.0;
         }
         for (j=0; j<m; j++) {
            a = me[i][j];
            if (j != i) for (k=0; k<n; k++) local[k] += a * solmat[j][k];
         }
         diag = me[i][i];
         if (diag > tol) {
            for (k=0; k<n; k++) {
               newsol[k] = (rhs[i][k] - local[k])/diag;
               cval[k] = (newsol[k] - oldsol[k])*relax;
               if (fabs(cval[k]) > cmax) cmax = fabs(cval[k]);
               newsol[k] = solmat[i][k] = oldsol[k] + cval[k];
            }
         }
         else if (diag > -tol) {
            throw exception("BGMatrix::gs_solve(): zero-diagonal found in Matrix::gs_solve()");
         }
         else {
            throw exception("BGMatrix::gs_solve(): matrix is not psd");
         }
      }
      niter += 1;
      if ( (niter % 10) ==0) {
         std::cout << " GS: # of iter = " << niter << ", max_change = "
              << cmax << std::endl;
      }
   } while (niter <= mxiter && cmax > stopval);
   if (cmax > stopval) {
      warning("BGMatrix::gs_solve(): not converge: %20.10f",cmax.f);
   }
}

void BGMatrix::gs_solve(const Vector<BG>& rhs, Vector<BG>& sol, const double relax,
                     const double stopval, const int mxiter)
{
   int n = rhs.size();
   BGMatrix bmat(n,1),solmat(n,1);
   if (sol.size() != n) {
      sol.resize(n);
   }
   int i;
   for (i=0; i<n; i++) {
      bmat[i][0] = rhs[i];
      solmat[i][0] = sol[i];
   }
   gs_solve(bmat,solmat,relax,stopval,mxiter);
   for (i=0; i<n; i++) sol[i] = solmat[i][0];
}


BG BGMatrix::det(void) const
{
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;

   BG detval = 0.0;
   if ( nrow != ncol ) throw exception("BGMatrix::det(): matrix must be square");
   int i,d;
   Vector<int> indx(nrow);
   BGMatrix temp(*this);
   ludcmp(temp.begin(),nrow,indx.begin(),d, SESSION.epsilon);
   detval = (BG)d;
   for (i=0; i<nrow; i++) detval *= temp[i][i];
   return detval;
}


BG BGMatrix::norm(const std::string &s) const
{
   BG retval = 0.0;
   BGMatrix A = this->transpose();
   if (s == "fro") {
      retval = sqrt((A*(*this)).diag(0).sum());
   }
   else {
      throw exception("BGMatrix::norm((): bad arg, it must be either \"inf\" or \"fro\" ");
   }
   return retval;
}

BG BGMatrix::logdet(void)
{
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;

   BG retval = 0.0;
   if (!this->symmetric()) throw exception("BGMatrix::logdet(): matrix must be symmetric");
   Matrix<BG> temp(*this);

   ginverse1(temp.begin(),nrow,retval,3,SESSION.epsilon);
   return retval;
}

BG BGMatrix::quadratic(const Vector<BG> &x, const Vector<BG> &y) const
{
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;
   int m = x.size();
   int n = y.size();

   BG s = 0.0;
   if (m != nrow  || n != ncol)  throw exception("BGMatrix::quadratic(): size not conformable");
   for (int j,i=0;i<m; i++) for (j=0; j<n; j++) s += x[i]*me[i][j]*y[j];
   return s;
}
Vector<BG> BGMatrix::variance(orientation orien) const
{
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;
   int i,j;
   BG ss,s;
   Vector<BG> temp;
   if (orien == COLUMN) {
      assert(nrow > 1);
      temp.resize(ncol);
      for (j=0; j<ncol; j++) {
         for (ss=0.0,s=0.0,i=0; i<ncol; ++i) {
             ss += me[i][j]*me[i][j];
             s  += me[i][j];
         }
         temp[j] = (ss-s*s/nrow)/(nrow - 1);
      }
   } else if (orien == ROW) {
      assert(ncol > 1);
      temp.resize(nrow);
      for (i=0; i<nrow; i++) {
         for (ss=0.0,s=0.0,j=0; j<ncol; ++j) {
             ss += me[i][j]*me[i][j];
             s  += me[i][j];
         }
         temp[i] = (ss-s*s/ncol)/(ncol - 1);
      }
   } else {
      warning("unknown orientation");
   }
   return temp;
}

BGMatrix BGMatrix::covariance(const BGMatrix &B) const
{
   int nrow =  Matrix<BG>::nrow;
   int ncol =  Matrix<BG>::ncol;
   BG** me =  Matrix<BG>::me;

   assert (nrow == B.num_rows() && ncol == B.num_cols());
   BGMatrix temp(ncol,ncol);
   int i,j,k;
   BG ss,sx,sy;
   for (i=0; i<ncol; ++i) {
      for (j=0; j<=i; ++j) {
         for (ss=0.0,sx=0.0,sy=0.0,k=0; k<nrow; ++k) {
            ss += me[k][i]*B[k][j];
            sx += me[k][i];
            sy += B[k][j];
         }
         temp[i][j] = temp[j][i] = (ss - sx*sy/ncol)/(ncol - 1);
      }
   }
   return temp;
}

BGMatrix& BGMatrix::sqrtm(void)
{
   if (!this->symmetric()) throw exception("BGMatrix::sqrtm(): matrix must be symmetric");
   BG lgdet;
   ginverse1(begin(),num_rows(),lgdet,2,SESSION.epsilon);
   return *this;
}

BGMatrix& BGMatrix::identity(const int m,const int n)
{
   resize(m,n);
   BG **me = begin();
   for (int i=0; i<m; ++i) if (i < n) me[i][i] = 1.0;
   return *this;
}

// BGMatrix& BGMatrix::sweep(const int i0,const int i1)
// {
//    int nrow =  Matrix<BG>::nrow;
//    int ncol =  Matrix<BG>::ncol;

//    int n = std::min(nrow,ncol);
//    if (i0<0 || i0 >=n || i0 > i1) throw exception("BGMatrix::sweep(): range error");
//    matvec::sweep(nrow,ncol,begin(),i0,i1,SESSION.epsilon);
//    return *this;
// }


}  //////////// end of namespace

