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

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "minimizer.h"
#include "matrix.h"
#include "session.h"


namespace matvec {
extern double   ranf(void);

static double qd[2];
static double fx,qf1;      // global variables accessable only with this file
  static double macheps, t,m2 ,m4;
  static double SmaLL,vsmall,vlarge,large;
static double dmin,ldt,h;

static int nf,nl;

static void svd(const int n,double eps,const double tol,double** ab,double* q);

static void sort(const int n,double* d, double** v) // d,v in descending order
{
   int k, i, j;
   double s;
   for (i=0; i<n-1; i++) {
       k = i; s = d[i];
       for (j=i+1; j<n; j++) if (d[j] > s) { k = j; s = d[j]; }
       if (k > i) {
	  d[k] = d[i]; d[i] = s;
	  for (j=0; j<n; j++) { s = v[j][i]; v[j][i] = v[j][k]; v[j][k] = s; }
       }
   }
}

static void print(const int n, const Vector<double> &x)
{
   int W = 6 + 6;
   std::cout << "\n";
   std::cout << "-- after " << nf << " function calls\n";
   std::cout << "-- including " << nl << " linear searches\n";
   std::cout << "-- your function f(x) has been reduced to "
        << std::setprecision(W)<< fx <<"\n";
   std::cout.precision(6);
   std::cout << "-- current values of x:";
   x.print();
   std::cout << std::endl;
}

static void flin(const double l, const int j,const Vector<double> &x,const int n,
            const Matrix<double> &v,Vector<double> &w, const Matrix<double> &q01)
{
   int i;
   double qa,qb,qc;
   if (j != -1) {		// linear search
      for (i=0; i<n; i++) w[i] = x[i] + l *v[i][j];
   }
   else {			// search along parabolic space curve
      qa = l*(l-qd[1])/(qd[0]*(qd[0]+qd[1]));
      qb = (l+qd[0])*(qd[1]-l)/(qd[0]*qd[1]);
      qc = l*(l+qd[0])/(qd[1]*(qd[0]+qd[1]));
      for (i=0; i<n; i++) w[i] = qa*q01[0][i]+qb*x[i]+qc*q01[1][i];
   }
}

double Minimizer::praxis(Vector<double> &x, const int n,int& maxfun,const double tol,
                         const double epsilon, const int pl)
{
   ////////////////////////////////////////////////////////////////////
   // references:
   //  - powell, m.j.d., 1964. an efficient method for finding
   //    the minimum of a function in several variables without
   //    calculating derivatives, computer journal, 7, 155-162
   //  - brent, r.p., 1973. algorithms for minimization without
   //    derivatives, prentice hall, englewood cliffs.
   //
   //  problems, suggestions or improvements are always wellcome
   //                    karl gegenfurtner   07/08/87
   //                                        c - version
   //  On Thu Oct 21 19:56:35 CDT 1993, Tianlin made it into C++ version
   //  note that lots of globel variables have been removed.
   //
   //
   // usage: min = praxis(fun, x, n);
   //
   // fun       the function to be minimized. fun is called from
   //           praxis with x and n as arguments
   // x         a double array containing the initial guesses for
   //           the minimum, which will contain the solution on
   //           return
   // n         an integer specifying the number of unknown
   //           parameters
   // min       praxis returns the least calculated value of fun
   //
   // some additional global variables control some more aspects of
   // the inner workings of praxis. setting them is optional, they
   // are all set to some reasonable default values given below.
   //
   //  prin      controls the printed output from the routine.
   //            0 -> no output
   //            1 -> print only starting and final values
   //            2 -> detailed map of the minimization process
   //            3 -> print also eigenvalues and vectors of the
   //                 search directions
   //            the default value is 1
   // tol        is the tolerance allowed for the precision of the
   //            solution. praxis returns if the criterion
   //            2 * ||x[k]-x[k-1]|| <= sqrt(macheps) * ||x[k]|| + tol
   //            is fulfilled more than ktm times.
   //            the default value 1.0e-16
   //            tol is usually equal to  EPISILON*EPISILON
   //            where
   // ktm        see just above. default is 1, and a value of 4 leads
   //            to a very(!) cautious stopping criterion.
   // step       is a steplength parameter and should be set equal
   //            to the expected distance from the solution.
   //            exceptionally small or large values of step lead to
   //            slower convergence on the first few iterations
   //            the default value for step is 1.0
   // scbd       is a scaling parameter. 1.0 is the default and
   //            indicates no scaling. if the scales for the different
   //            parameters are very different, scbd should be set to
   //            a value of about 10.0.
   // illc       should be set to true (1) if the problem is known to
   //            be ill-conditioned. the default is false (0). this
   //            variable is automatically set, when praxis finds
   //            the problem to be ill-conditioned during iterations.
   // maxfun     is the maximum number of calls to fun allowed. praxis
   //            will return after maxfun calls to fun even when the
   //            minimum is not yet found. the default value of 0
   //            indicates no limit on the number of calls.
   //            this return condition is only checked every n
   //            iterations.
   //
   //////////////////////////////////////////////////////////////////////
   int i,j,k,k2, W = 6 + 6;

   int  ktm = 1;
   int illc = 0;
   double scbd = 1.0;
   double step = 10.0;         // after the first n*10 loop, step reset to 1
   Vector<double> y(n);
   Vector<double> z(n);
   Vector<double> d(n);
//   d[0] = 0.0;

   Vector<double> w(n);              // working space for flin
   Matrix<double> v(n,n);
   Matrix<double> q01(2,n);

   qd[0] = 0.0;
   prtlevel = pl;
   macheps = epsilon;
   SmaLL = macheps*macheps;
   vsmall = SmaLL*SmaLL;
   large = 1.0/SmaLL;
   vlarge = 1.0/vsmall;
   m2 = std::sqrt(macheps); m4 = std::sqrt(m2);

   nl = 0;
   nf = 1;
   double t2,sl,dn,df,sf,lds,f1,s;
   h = step;

   dmin = SmaLL;
   t2 = SmaLL + std::fabs(tol);
   t = t2;
   double ldfac = (illc ? 0.1 : 0.01);
   int kt = 0;
   int kl;

   fx = minfun(x, n);
   qf1 = fx;

   if (h < 100.0*t) h = 100.0*t;
   ldt = h;
   for (i=0; i<n; i++) for (j=0; j<n; j++) v[i][j] = (i == j ? 1.0 : 0.0);
   for (i=0; i<n; i++) q01[1][i] = x[i];
   if (prtlevel > 1) {
      std::cout << "\n------------- enter function praxis -----------\n";
      std::cout << "-- current parameter settings\n";
      std::cout << "-- scaling: " << std::setw(W) << scbd   << "\n";
      std::cout << "-- tol    : " << std::setw(W) << t      << "\n";
      std::cout << "-- maxstep: " << std::setw(W) << h      << "\n";
      std::cout << "-- illc   : " << std::setw(W) << illc   << "\n";
      std::cout << "-- ktm    : " << std::setw(W) << ktm    << "\n";
      std::cout << "-- maxite : " << std::setw(W) << maxfun << "\n";
   }
   if (prtlevel) print(n,x);

mloop:
   sf = d[0];
   s = d[0] = 0.0;

   /* minimize along first direction */
   min_first_dir(0,2,&d[0],&s,fx,0,x,n,v,w,q01);
   if (s <= 0.0) for (i=0; i < n; i++) v[i][0] = -v[i][0];
   if ((sf <= (0.9 * d[0])) || ((0.9 * sf) >= d[0])) {
      for (i=1; i<n; i++) d[i] = 0.0;
   }
   for (k=1; k<n; k++) {
      for (i=0; i<n; i++) y[i] = x[i];
      sf = fx;
      illc = illc || (kt > 0);
next:
      kl = k;
      df = 0.0;
      if (illc) {        /* random step to get off resolution valley */
         for (i=0; i<n; i++) {
            z[i] = (0.1*ldt + t2*pow(10.0,(double)kt))*(ranf()-0.5);
            s = z[i];
            for (j=0; j < n; j++) x[j] += s * v[j][i];
  	 }
         fx = minfun(x, n);
         nf++;
      }
       /* minimize along non-conjugate directions */
      for (k2=k; k2<n; k2++) {
         sl = fx;
         s = 0.0;
         min_first_dir(k2,2,&d[k2],&s,fx,0,x,n,v,w,q01);
         if (illc) {
            s = d[k2] * (s + z[k2]) * (s + z[k2]);
         }
         else
            s = sl - fx;
         if (df < s) {
            df = s;
            kl = k2;
         }
      }

      if (!illc && (df < std::fabs(100.0 * macheps * fx))) {
         illc = 1;
         goto next;
      }
      if (k == 1 && prtlevel > 1) {
         std::cout << "\n-- New Direction:\n";
         d.print();
      }
       // minimize along conjugate directions
      for (k2=0; k2<=k-1; k2++) {
         s = 0.0;
         min_first_dir(k2,2,&d[k2],&s,fx,0,x,n,v,w,q01);
      }
      f1 = fx;
      fx = sf;
      lds = 0.0;
      for (i=0; i<n; i++) {
         sl = x[i];
         x[i] = y[i];
         y[i] = sl - y[i];
         sl = y[i];
         lds = lds + sl*sl;
      }
      lds = std::sqrt(lds);
      if (lds > SmaLL) {
         for (i=kl-1; i>=k; i--) {
            for (j=0; j < n; j++) v[j][i+1] = v[j][i];
            d[i+1] = d[i];
         }
         d[k] = 0.0;
         for (i=0; i < n; i++) v[i][k] = y[i] / lds;
         min_first_dir(k,4,&d[k],&lds,f1,1,x,n,v,w,q01);
         if (lds <= 0.0) {
            lds = -lds;
            for (i=0; i<n; i++) v[i][k] = -v[i][k];
         }
      }
      ldt = ldfac * ldt;
      if (ldt < lds) ldt = lds;
      if (prtlevel > 1) print(n,x);
      t2 = 0.0;
      for (i=0; i<n; i++) t2 += x[i]*x[i];
      t2 = m2 * std::sqrt(t2) + t;
      if (ldt > (0.5 * t2))
         kt = 0;
      else
         kt++;
      if (kt > ktm) goto fret;
   }
   /////////////////////////////////////////////
   //  try quadratic extrapolation in case    //
   //  we are stuck in a curved valley        //
   /////////////////////////////////////////////
   quadratic(x,n,v,w,q01);
   dn = 0.0;
   for (i=0; i<n; i++) {
      d[i] = 1.0 / std::sqrt(d[i]);
      if (dn < d[i]) dn = d[i];
   }
   if (prtlevel > 2) {
      std::cout << "\n-- New Matrix of Directions:\n";
      v.print();
   }

   for (j=0; j<n; j++) {
      s = d[j] / dn;
      for (i=0; i < n; i++) v[i][j] *= s;
   }
   if (scbd > 1.0) {               // scale axis to reduce condition number
      s = vlarge;
      for (i=0; i<n; i++) {
         sl = 0.0;
         for (j=0; j < n; j++) sl += v[i][j]*v[i][j];
         z[i] = std::sqrt(sl);
         if (z[i] < m4) z[i] = m4;
         if (s > z[i]) s = z[i];
      }
      for (i=0; i<n; i++) {
         sl = s / z[i];
         z[i] = 1.0 / sl;
         if (z[i] > scbd) {
            sl = 1.0 / scbd;
            z[i] = scbd;
         }
      }
   }
   for (i=1; i<n; i++) {
      for (j=0; j<=i-1; j++) {
         s = v[i][j];
         v[i][j] = v[j][i];
         v[j][i] = s;
      }
   }
   svd(n, macheps, vsmall, v.begin(), d.begin());
   if (scbd > 1.0) {
      for (i=0; i<n; i++) {
         s = z[i];
         for (j=0; j<n; j++) v[i][j] *= s;
      }
      for (i=0; i<n; i++) {
         s = 0.0;
         for (j=0; j<n; j++) s += v[j][i]*v[j][i];
         s = std::sqrt(s);
         d[i] *= s;
         s = 1.0 / s;
         for (j=0; j<n; j++) v[j][i] *= s;
      }
   }
   for (i=0; i<n; i++) {
      if ((dn * d[i]) > large)
         d[i] = vsmall;
      else if ((dn * d[i]) < SmaLL)
         d[i] = vlarge;
      else
         d[i] = pow(dn * d[i],-2.0);
   }
   sort(n,d.begin(),v.begin());               // the new eigenvalues and eigenvectors
   dmin = d[n-1];
   if (dmin < SmaLL) dmin = SmaLL;
   illc = (m2 * d[0]) > dmin;
   if (prtlevel > 2 && scbd > 1.0) {
      std::cout << "\n-- Scale Factors:\n";
      z.print();
   }
   if (prtlevel > 2) {
      std::cout << "\n-- Eigenvalues of A:\n";
      d.print();
      std::cout << "\n-- Eigenvectors of A:\n";
      v.print();
   }

   if (maxfun > 0 && nl > maxfun) {
      if (prtlevel) std::cout << "\n-- maxiimum # of function calls reached \n";
      goto fret;
   }
   if (nf > n*10) step = 1.0;
   goto mloop; 	 // back to main loop

fret:
   if (prtlevel > 0) {
      std::cout << "\n--after " << nf << " function calls\n";
      std::cout << "-- your function f(x) has finally been minimized to "
           << std::setprecision(W) << fx << "\n";
      std::cout.precision(6);
      std::cout << "-- with final x:\n";
      x.print();
   }
   maxfun = nf;
   return fx;
}

void Minimizer::min_first_dir(const int j,const int nits,double *d2,
                              double *x1,
                              double f1, const int fk,Vector<double> &x,const int n,
                              const Matrix<double> &v, Vector<double> &w,Matrix<double> &q01)
{
   int k, i, dz;
   double x2, xm, f0, f2, fm, d1, t2, s, sf1, sx1;

   sf1 = f1; sx1 = *x1;
   k = 0; xm = 0.0; fm = f0 = fx; dz = *d2 < macheps;
   ///////////////////////////////////
   //        find step size        //
   //////////////////////////////////
   s = 0;
   for (i=0; i<n; i++) s += x[i]*x[i];
   s = std::sqrt(s);
   if (dz) { t2 = m4*std::sqrt(std::fabs(fx)/dmin + s*ldt) + m2*ldt; }
   else    { t2 = m4*std::sqrt(std::fabs(fx)/(*d2) + s*ldt) + m2*ldt; }
   s = s*m4 + t;
   if (dz && t2 > s) t2 = s;
   if (t2 < SmaLL) t2 = SmaLL;
   if (t2 > 0.01*h) t2 = 0.01 * h;
   if (fk && f1 <= fm) {
      xm = *x1;
      fm = f1;
   }
   if (!fk || std::fabs(*x1) < t2) {
      *x1 = (*x1 > 0 ? t2 : -t2);
      flin(*x1, j,x,n,v,w,q01);
      nf++;
      f1 = minfun(w,n);
   }
   if (f1 <= fm) {
      xm = *x1;
      fm = f1;
   }
next:
   if (dz) {
      x2 = (f0 < f1 ? -(*x1) : 2*(*x1));
      flin(x2, j,x,n,v,w,q01);
      nf++;
      f2 = minfun(w,n);
      if (f2 <= fm) {
         xm = x2;
	 fm = f2;
      }
      *d2 = (x2*(f1-f0) - (*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
   }
   d1 = (f1-f0)/(*x1) - *x1**d2; dz = 1;
   if (*d2 <= SmaLL) { x2 = (d1 < 0 ? h : -h); }
   else              { x2 = - 0.5*d1/(*d2); }
   if (std::fabs(x2) > h) x2 = (x2 > 0 ? h : -h);
tryit:
   flin(x2, j,x,n,v,w,q01);
   nf++;
   f2 = minfun(w,n);
   if ((k < nits) && (f2 > f0)) {
      k++;
      if ((f0 < f1) && (*x1*x2 > 0.0)) goto next;
      x2 *= 0.5;
      goto tryit;
   }
   nl++;
   if (f2 > fm) x2 = xm; else fm = f2;
   if (std::fabs(x2*(x2-*x1)) > SmaLL) {
      *d2 = (x2*(f1-f0) - *x1*(fm-f0))/(*x1*x2*(*x1-x2));
   }
   else {
      if (k > 0) *d2 = 0;
   }
   if (*d2 <= SmaLL) *d2 = SmaLL;
   *x1 = x2; fx = fm;
   if (sf1 < fx) {
      fx = sf1;
      *x1 = sx1;
   }
   if (j != -1) for (i=0; i<n; i++) x[i] += (*x1)*v[i][j];
}

void Minimizer::quadratic(Vector<double> &x,const int n,const Matrix<double> &v,
                          Vector<double> &w,Matrix<double> &q01)
{                               // look for a minimum along the curve q0, q1 //
   int i;
   double l, s,qa,qb,qc;
   s = fx; fx = qf1; qf1 = s;   // swap fx and qf1
   qd[1] = 0.0;
   for (i=0; i<n; i++) {
       s = x[i]; l = q01[1][i]; x[i] = l; q01[1][i] = s;
       qd[1] += (s-l)*(s-l);
   }
   s = 0.0; qd[1] = std::sqrt(qd[1]); l = qd[1];
   if ( qd[0]>0.0 && qd[1]>0.0 && nl>=3*n*n ) {
      min_first_dir(-1,2, &s,&l,qf1,1,x, n,v,w,q01);
      qa = l*(l-qd[1])/(qd[0]*(qd[0]+qd[1]));
      qb = (l+qd[0])*(qd[1]-l)/(qd[0]*qd[1]);
      qc = l*(l+qd[0])/(qd[1]*(qd[0]+qd[1]));
   }
   else { qa = qb = 0.0; qc = 1.0; }
   qd[0] = qd[1];
   for (i=0; i<n; i++) {
      s = q01[0][i]; q01[0][i] = x[i];
      x[i] = qa*s + qb*x[i] + qc*q01[1][i];
   }
}

/////////////////////////////////////////////
//  singular value decomposition   A = U*S*V'
//  s is returned with values ascendingly
///////////////////////////////////////////////////

static void svd(const int n,double eps,const double tol,double** ab,double* q)
{
   int l, kt, l2, i, j, k;
   double c, f, g, h, s, x, y, z;
   Vector<double> e(n);
   // householder's reduction to bidiagonal form
   l = 0;        // added by  T wang to avoid warning from g++
   x = g = 0.0;
   for (i=0; i<n; i++) {
      e[i] = g; s = 0.0;
      l = i+1;
      for (j=i; j<n; j++) s += ab[j][i] * ab[j][i];
      if (s < tol) {
         g = 0.0;
      }
      else {
	 f = ab[i][i];
         if (f < 0.0) {g = std::sqrt(s); }
	 else         { g = - std::sqrt(s); }
	 h = f*g - s; ab[i][i] = f - g;
         for (j=l; j<n; j++) {
	    f = 0.0;
	    for (k=i; k<n; k++) f += ab[k][i] * ab[k][j];
	    f /= h;
	    for (k=i; k<n; k++) ab[k][j] += f * ab[k][i];
	 }
      }
      q[i] = g; s = 0.0;
      if (i < n) for (j=l; j<n; j++) s += ab[i][j] * ab[i][j];
      if (s < tol) {
	 g = 0.0;
      }
      else {
         f = ab[i][i+1];
	 if (f < 0.0) { g = std::sqrt(s); }
	 else         { g = - std::sqrt(s); }
	 h = f*g - s; ab[i][i+1] = f - g;
	 for (j=l; j<n; j++) e[j] = ab[i][j]/h;
	 for (j=l; j<n; j++) {
            s = 0;
	    for (k=l; k<n; k++) s += ab[j][k]*ab[i][k];
	    for (k=l; k<n; k++) ab[j][k] += s * e[k];
	 }
      }
      y = std::fabs(q[i]) + std::fabs(e[i]);
      if (y > x) x = y;
   }
   /* accumulation of right hand transformations */
   for (i=n-1; i >= 0; i--) {
      if (g != 0.0) {
         h = ab[i][i+1]*g;
	 for (j=l; j<n; j++) ab[j][i] = ab[i][j] / h;
         for (j=l; j<n; j++) {
            s = 0.0;
	    for (k=l; k<n; k++) s += ab[i][k] * ab[k][j];
	    for (k=l; k<n; k++) ab[k][j] += s * ab[k][i];
	 }
      }
      for (j=l; j<n; j++) ab[i][j] = ab[j][i] = 0.0;
      ab[i][i] = 1.0; g = e[i]; l = i;
   }
    // diagonalization to bidiagonal form 
   eps *= x;
   for (k=n-1; k>= 0; k--) {
      kt = 0;
TestFsplitting:
      if (++kt > 30) {
         e[k] = 0.0;
         std::cerr << "\n+++ qr failed\n";
      }
      for (l2=k; l2>=0; l2--) {
         l = l2;
	 if (std::fabs(e[l]) <= eps) goto TestFconvergence;
         if (std::fabs(q[l-1]) <= eps) break;
      }
      c = 0.0; s = 1.0;
      for (i=l; i<=k; i++) {
         f = s * e[i]; e[i] *= c;
	 if (std::fabs(f) <= eps) goto TestFconvergence;
	 g = q[i];
   	 if (std::fabs(f) < std::fabs(g)) {
	    double fg = f/g;
	    h = std::fabs(g)*std::sqrt(1.0+fg*fg);
	 }
	 else {
	    double gf = g/f;
	    h = (f!=0.0 ? std::fabs(f)*std::sqrt(1.0+gf*gf) : 0.0);
         }
         q[i] = h;
	 if (h == 0.0) { h = 1.0; g = 1.0; }
	 c = g/h; s = -f/h;
      }
TestFconvergence:
      z = q[k];
      if (l == k) goto Convergence;
       // shift from bottom 2x2 minor 
      x = q[l]; y = q[k-l]; g = e[k-1]; h = e[k];
      //BRS added to avoid division by y=zero
      if (y) {
	f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h*y);
      }
      else
	f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h);
      //BRS
      g = std::sqrt(f*f+1.0);
      if (f <= 0.0)
         f = ((x-z)*(x+z) + h*(y/(f-g)-h))/x;
      else
         f = ((x-z)*(x+z) + h*(y/(f+g)-h))/x;
       // next qr transformation 
      s = c = 1.0;
      for (i=l+1; i<=k; i++) {
         g = e[i]; y = q[i]; h = s*g; g *= c;
         if (std::fabs(f) < std::fabs(h)) {
	    double fh = f/h;
	    z = std::fabs(h) * std::sqrt(1.0 + fh*fh);
	 }
	 else {
	    double hf = h/f;
	    z = (f!=0.0 ? std::fabs(f)*std::sqrt(1.0+hf*hf) : 0.0);
	 }
	 e[i-1] = z;
	 if (z == 0.0) f = z = 1.0;
	 c = f/z; s = h/z;
	 f = x*c + g*s; g = - x*s + g*c; h = y*s;
	 y *= c;
	 for (j=0; j<n; j++) {
	    x = ab[j][i-1]; z = ab[j][i];
	    ab[j][i-1] = x*c + z*s;
	    ab[j][i] = - x*s + z*c;
	 }
	 if (std::fabs(f) < std::fabs(h)) {
	    double fh = f/h;
	    z = std::fabs(h) * std::sqrt(1.0 + fh*fh);
	 }
	 else {
	    double hf = h/f;
	    z = (f!=0.0 ? std::fabs(f)*std::sqrt(1.0+hf*hf) : 0.0);
	 }
         q[i-1] = z;
	 if (z == 0.0) z = f = 1.0;
	 c = f/z; s = h/z;
	 f = c*g + s*y; x = - s*g + c*y;
      }
      e[l] = 0.0; e[k] = f; q[k] = x;
      goto TestFsplitting;
Convergence:
      if (z < 0.0) {
         q[k] = - z;
         for (j=0; j<n; j++) ab[j][k] = - ab[j][k];
      }
   }
}

} ///// end of namespace matvec

