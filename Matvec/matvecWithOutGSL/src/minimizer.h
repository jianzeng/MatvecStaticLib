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

#ifndef Minimizer_H
#define Minimizer_H

#include "matrix.h"

/*!
   \class Minimizer Minimizer.h
   \brief a minimizer to minimize a function

*/
namespace matvec {

class Minimizer {
   private:
      int prtlevel;
      void    min_first_dir(const int j,const int nits,double* d2,double* x1,
                            double f1, const int fk,Vector<double> &x,const int n,
                            const Matrix<double> &v, Vector<double> &w,Matrix<double> &q01);
      void    quadratic(Vector<double> &x,const int n,const Matrix<double> &v,
                        Vector<double> &w,Matrix<double> &q01);
   protected:
      int minfun_indx;

   public:
      Minimizer(void)  {prtlevel = 2; minfun_indx = 0;}
      virtual ~Minimizer(void) {;}

      virtual double minfun(const Vector<double> &x, const int n) = 0;
      void    min_prtlevel(const int pl) {prtlevel = pl;}
      double  praxis(Vector<double> &x,const int n, int& maxfun,
                     const double tol=1.0e-16,
                     const double epsilon=1.0e-8,const int pl=2);
};
}
#endif
