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

#ifndef MATVEC_DOUBLE_MATRIX_H
#define MATVEC_DOUBLE_MATRIX_H
#include "session.h"
#include "matrix.h"

namespace matvec {
class doubleMatrix : public Matrix<double> {
  friend class Model;
public:
   doubleMatrix(void) : Matrix<double>() {}
   explicit doubleMatrix(const size_type m,const size_type n) : Matrix<double>(m,n)  {}
   doubleMatrix(const size_type m,const size_type n,const double** a) : Matrix<double>(m,n,a) {}
   doubleMatrix(const doubleMatrix& a) : Matrix<double>(a) { }
   doubleMatrix(const Matrix<double>& a) : Matrix<double>(a) { }

   doubleMatrix&   operator = (const Matrix<double> &a) { copy(a); return *this; }
   doubleMatrix&   operator = (const Matrix<bool> &a);

   bool            psd(void) const;
   void            svd(doubleMatrix& u,Vector<double>& s,doubleMatrix& v);
   int             rank(void);
   doubleMatrix&   inv();
   doubleMatrix    ginv0(void) const;
   doubleMatrix    splines(const doubleMatrix knots,const unsigned type=0) const;
   doubleMatrix&   ginv1(unsigned *irank=0);
//   doubleMatrix&   MPGinv();
  doubleMatrix   mat_log(double tol=0) const; //Returns Matrix log
  doubleMatrix   mat_exp(double tol=0) const;  //Returns Matrix exp
  void    mat_exp_der(Vector<doubleMatrix> &der,double tol=0) const; // Returns  partial derivatives
   doubleMatrix    covariance(const doubleMatrix &B) const;
   Vector<double>  variance(orientation orien = COLUMN) const;
   doubleMatrix    lu_solve(const doubleMatrix& rhs);
   Vector<double>  lu_solve(const Vector<double>& rhs);
   void            gs_solve(const doubleMatrix& v,doubleMatrix& solmat,const double relax=1.0,
                      const double stopval=0.001,const int mxiter=1000);
   void            gs_solve(const Vector<double>& v,Vector<double>& solvec,const double relax=1.0,
                      const double stopval=0.001,const int mxiter=1000);

   Vector<double>  eigen(const int job=1);
   double          cond(void);
   double          det(void) const;
   double          norm(const int p) const;
   double          norm(const std::string &s) const;
   double          logdet(void);
   double          quadratic(const Vector<double> &a,const Vector<double> &b) const;
   doubleMatrix&   sqrtm(void);
   doubleMatrix&   identity(const int m,const int n);
   doubleMatrix&   identity(const int m) {return identity(m,m);}
   doubleMatrix&   identity(void) {return identity(num_rows(),num_cols());}
   doubleMatrix&   sweep(const size_type i0,const size_type i1);
   doubleMatrix&   sweep(void) {return sweep(num_rows() - 1, num_cols() - 1);}
  bool            symmetric(void) const {if (!me) return false;
  for (int i=1; i<nrow; ++i) for (int j=0; j<i; ++j) if (fabs(me[i][j]-me[j][i])> (fabs(me[i][j]+me[j][i])*SESSION.epsilon) && fabs(me[i][j]-me[j][i]) > SESSION.epsilon) return false;
   return true;
  } 
//   double MPLogDet();
protected:
  //   int nrow,ncol;
  //   double** me;
};

} ///////// end of namespace matvec
#endif

