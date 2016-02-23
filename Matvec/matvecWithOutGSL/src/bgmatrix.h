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

#ifndef MATVEC_BG_MATRIX_H 
#define MATVEC_BG_MATRIX_H

#include <string>
#include <fstream>
#include "session.h"
#include "matrix.h"
#include "vector.h"
#include "bg.h"

namespace matvec {
class BGMatrix : public Matrix<BG> {
  friend class Model;
public:
   BGMatrix(void) : Matrix<BG>() {}
   BGMatrix(const doubleMatrix &a);
   explicit BGMatrix(const size_type m,const size_type n)     : Matrix<BG>(m,n){}
   BGMatrix(const size_type m,const size_type n,const BG** a) : Matrix<BG>(m,n,a){}
   BGMatrix(const BGMatrix& a) : Matrix<BG>(a) { }
   BGMatrix(const Matrix<BG>& a) : Matrix<BG>(a) { }

   BGMatrix&   operator = (const Matrix<BG> &a) { copy(a); return *this; }
   BGMatrix&   operator = (const Matrix<bool> &a);

   bool            psd(void) const;
   //void            svd(BGMatrix& u,Vector<BG>& s,BGMatrix& v);
   //int             rank(void);
   BGMatrix&   inv();
   BGMatrix    ginv0(void) const;
   BGMatrix    splines(const BGMatrix knots,const unsigned type=0) const;
   BGMatrix&   ginv1(unsigned *irank=0);
   BGMatrix   mat_log(double tol=0.0) const;  //Returns Matrix log
   BGMatrix   mat_exp(double tol=0.0) const;  //Returns Matrix exp
  void    mat_exp_der(Vector<BGMatrix> &der,double tol=0.0) const; // Returns  partial derivatives
   BGMatrix    covariance(const BGMatrix &B) const;
   Vector<BG>  variance(orientation orien = COLUMN) const;
   BGMatrix    lu_solve(const BGMatrix& rhs);
   Vector<BG>  lu_solve(const Vector<BG>& rhs);
   void            gs_solve(const BGMatrix& v,BGMatrix& solmat,const double relax=1.0,
                      const double stopval=0.001,const int mxiter=1000);
   void            gs_solve(const Vector<BG>& v,Vector<BG>& solvec,const double relax=1.0,
                      const double stopval=0.001,const int mxiter=1000);

   //Vector<BG>  eigen(const int job=1);
   BG          cond(void);
   BG          det(void) const;
   BG          norm(const int p) const;
   BG          norm(const std::string &s) const;
   BG          logdet(void);
   BG          quadratic(const Vector<BG> &a,const Vector<BG> &b) const;
   BGMatrix&   sqrtm(void);
   BGMatrix&   identity(const int m,const int n);
   BGMatrix&   identity(const int m) {return identity(m,m);}
   BGMatrix&   identity(void) {return identity(num_rows(),num_cols());}
/*    BGMatrix&   sweep(const size_type i0,const size_type i1); */
/*    BGMatrix&   sweep(void) {return sweep(num_rows() - 1, num_cols() - 1);} */
   bool        symmetric(void) const {
     if (!me) return false;
     for (int i=1; i<nrow; ++i) {
       for (int j=0; j<i; ++j) {
	 if (fabs(me[i][j]-me[j][i])> (fabs(me[i][j]+me[j][i])*SESSION.epsilon) && fabs(me[i][j]-me[j][i]) > SESSION.epsilon) return false;
       }
     }
   return true;
  } ;

protected:
  //   int nrow,ncol;
  //   BG** me;
};

} ///////// end of namespace matvec
#endif

