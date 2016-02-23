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

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "exception.h"
#include "session.h"
#ifdef WIN32
extern "C"{
  extern double erf(double x);
}
#endif

namespace matvec {

extern void warning(const char format[],...);

#define ITMAX 400

// random number generator for uniform [0..1] inclusively.

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(void)
{
   static long ix1,ix2,ix3;
   static double r[98];
   double temp = 0.0;
   static int iff=0;
   int j;

   if (SESSION.rand_seed < 0 || iff == 0) {
      iff=1;
      ix1=(IC1 - SESSION.rand_seed) % M1;
      ix1=(IA1*ix1+IC1) % M1;
      ix2=ix1 % M2;
      ix1=(IA1*ix1+IC1) % M1;
      ix3=ix1 % M3;
      for (j=1;j<=97;j++) {
         ix1=(IA1*ix1+IC1) % M1;
         ix2=(IA2*ix2+IC2) % M2;
         r[j]=(ix1+ix2*RM2)*RM1;
      }
      SESSION.rand_seed=1;
   }
   ix1=(IA1*ix1+IC1) % M1;
   ix2=(IA2*ix2+IC2) % M2;
   ix3=(IA3*ix3+IC3) % M3;
   j = static_cast<int>(1 + ((97*ix3)/M3));
   if (j <= 97 && j >= 1) {
      temp=r[j];
      r[j]=(ix1+ix2*RM2)*RM1;
   }
   else {
      throw exception("ran1(): it can't happen");
   }
   return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

// random number generator for standard normal density

double gasdev(void)
{
   static int iset=0;
   static double gset;
   double fac,r,v1,v2;
   if  (iset == 0) {
      do {
         v1=2.0*ran1()-1.0;
         v2=2.0*ran1()-1.0;
         r=v1*v1+v2*v2;
      } while (r >= 1.0 || r == 0.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else {
      iset=0;
      return gset;
   }
}

double gamdev(const double A,const double R)
{
   /////////////////////////////////////////////////////////////////////
   //   Generates random deviates from the standard gamma distribution
   //    whose density is
   //              1
   //       --------------- * X^(A-1) * Exp(-X/R)
   //        (R^A)*Gamma(A)
   //
   //    Chi-Square can be generated as gamdev(df/2,2.0)
   //
   //    Tianlin Wang at UIUC, May 1,1992
   /////////////////////////////////////////////////////////////////////
   int K,j;
   double am,e,s,v1,v2,x,y;
   K = static_cast<int>(A);
   if (K < 1) throw exception(" gamdev(): bad arg1, value >= 1 expected ");
   if (K < 6) {
      x=1.0;
      for (j=1;j<=K;j++) x *= ran1();
      x = -log(x);
   } else {
      do {
         do {
            do {
               v1=2.0*ran1()-1.0;
               v2=2.0*ran1()-1.0;
            } while (v1*v1+v2*v2 > 1.0);
            y=v2/v1;
            am=K-1;
            s=sqrt(2.0*am+1.0);
            x=s*y+am;
         } while (x <= 0.0);
         e=(1.0+y*y)*exp(am*log(x/am)-s*y);
      } while (ran1() > e);
   }
   return (x*R);
}

double Normal_cdf(const double x)
{
   // ******************************************************************
   // return the cumulative distribution value for a standard normal
   // *****************************************************************
   double pr = 0.5;
   if (x == 0.0) return pr;
   register double c = 1.0/sqrt(2.0);      // c = 1/sqrt(2)
   if (x <0.0) { pr = 0.5 - 0.5*erf(-x*c);}
   else        { pr = 0.5 + 0.5*erf(x*c);}
   return pr;
}

double betacf(const double a,const double b,const double x)
{
   double qap,qam,qab,em,tem,d;
   double bz,bm=1.0,bp,bpp;
   double az=1.0,am=1.0,ap,app,aold;
   int m;

   qab=a+b;
   qap=a+1.0;
   qam=a-1.0;
   bz=1.0-qab*x/qap;
   for (m=1; m<=ITMAX; m++) {
      em = static_cast<double>(m);
      tem = em+em;
      d = em*(b-em)*x/((qam+tem)*(a+tem));
      ap = az+d*am;
      bp = bz+d*bm;
      d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
      app = ap+d*az;
      bpp = bp+d*bz;
      aold = az;
      am = ap/bpp;
      bm = bp/bpp;
      az = app/bpp;
      bz = 1.0;
      if (fabs(az-aold) < (10.0*SESSION.epsilon*fabs(az))) return az;
   }
   warning("betacf(), not converged");
   abort();
   return 0.0;
}

double gammln(const double xx)
{
   if (xx <= 0.0) throw exception(" gammln(): bad arg, value > 0 expected");
   double x,tmp,ser;
   static double cof[6]={76.18009173,-86.50532033,24.01409822,
                         -1.231739516,0.120858003e-2,-0.536382e-5};
   x = xx-1.0;
   tmp = x + 5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.0;
   for (int j=0; j<=5; j++) {
      x += 1.0;
      ser += cof[j]/x;
   }
   return -tmp+log(ser*2.50662827465);
}

double gammin(const double x,const double p)
{
   //////////////////////////////////////////////////////////
   // ALGORITHM AS239 Appl. Statist. (1988) 37(3), GAMMAD()
   // return the incomplete Gamma function for arguments:
   //       x between 0 and +infinite
   //       p positive
   // REQUIREMENTS: Normal_cdf(), gammln()
   ////////////////////////////////////////////////////////
   double zero,one,two,three,nine,tol,oflo,xbig,plimit,elimit;
   double pn1,pn2,pn3,pn4,pn5,pn6,rn,arg,a,b,c,an;
   zero = 0.0; one = 1.0; two = 2.0; three = 3.0; nine = 9.0;
   tol = 1.0e-14; oflo = 1.0e+37; xbig = 1.0e+8;
   plimit = 1000.0; elimit = -88.0;
   double retval = zero;
   if (x < zero || p <= zero) throw exception(" gammin(): bad args");
   if (x == zero) return retval;
   if (p > plimit) {      // use a normal approximation
      pn1 = three*sqrt(p)*(pow(x/p,one/three) + one/(nine*p) - one);
      return Normal_cdf(pn1);
   }
   if (x > xbig) return one;   // x is extremely large compared to p
   if (x <= one || x < p) {    //  use Pearson's series expansion
      arg = p*log(x) - x - gammln(p + one);
      c = one;
      retval = one;
      a = p;
      do {
         a++;
         c *= x/a;
         retval += c;
      } while (c > tol);
      arg += log(retval);
      retval = zero;
      if (arg >= elimit) retval = exp(arg);
   }
   else {
      arg = p*log(x) - x - gammln(p);
      a = one - p;
      b = a + x + one;
      c = zero;
      pn1 = one;   pn2 = x;   pn3 = x + one;   pn4 = x*b;
      retval = pn3/pn4;
      for(;;) {
         a += one; b += two; c += one;
         an = a*c;
         pn5 = b*pn3 - an*pn1;
         pn6 = b*pn4 - an*pn2;
         if (fabs(pn6) > zero) {
            rn = pn5/pn6;
            if (fabs(retval - rn) <= tol*(rn < one ? rn: one)) break;
            retval = rn;
         }
         pn1 = pn3;   pn2 = pn4;   pn3 = pn5;   pn4 = pn6;
         if (fabs(pn5) >= oflo) { // re-scale terms, if terms are large
            pn1 /= oflo;  pn2 /= oflo;  pn3 /= oflo;  pn4 /= oflo;
         }
      }
      arg += log(retval);
      retval = one;
      if (arg >= elimit) retval = one - exp(arg);
   }
   return retval;
}

double betain(const double x,const double a,const double b, const double beta)
{
   //////////////////////////////////////////////////////////
   // return  incomplete beta function for arguments:
   //            x between 0 and 1,
   //            a and b positive,
   //            beta, log of complete beta function value
   ///////////////////////////////////////////////////////////
   if (x < 0.0 || x > 1.0) throw exception(" betain(%f): bad arg1, value between (0,1) expected");
   if (x == 0.0) return 0.0;
   if (x == 1.0) return 1.0;
   double fv,bt;
   if (x < SESSION.epsilon || fabs(x - 1.0) < SESSION.epsilon) {
      bt=0.0;
   }
   else {
      fv = a*log(x) + b*log(1.0-x) - beta;
      if (fv < -300) { bt = 0.0;}
      else           { bt = exp(fv);}
   }
   if (x < (a+1.0)/(a+b+1.0)) { fv =  bt*betacf(a,b,x)/a; }
   else                       { fv = 1.0-bt*betacf(b,a,1.0-x)/b; }
   return fv;
}

double factrl(const long n)
{
   static int ntop = 8;
   static double a[33]={1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0};
   int j;
   if (n > 32) return exp(gammln(static_cast<double>(n+1)));
   while (ntop < n) {
      j = ntop++;
      a[ntop] = a[j]*ntop;
   }
   return a[n];
}

double BinCoef(const long n,const long k)
{

   if (n<0L || k<0L || k > n) throw exception("BinCoef(): bad args");
   double a,b,c;
   a = gammln(static_cast<double>(n+1));
   b = gammln(static_cast<double>(k+1));
   c = gammln(static_cast<double>(n-k));
   return std::floor(0.5 + std::exp(a - b - c));
}

double Gamma_cdf(const double g,const double alfa,const double theta)
{
   double retval = 0.0;
   if (g <= 0.0 ) return retval;
   if (alfa > 0.0 ) {
      if (theta > 0.0 ) {
         retval = gammin(g/theta,alfa);
      }
      else {
         throw exception(" Gamma_cdf(): bad arg3, value > 0 expected");
      }
   }
   else {
      throw exception(" Gamma_cdf(): bad arg2: %f, value > 0 expected");
   }
   return retval;
}

#define ITRMAX 100000
#define ERRMAX 1.0e-10

double ChiSquare_cdf(const double x,const double df,const double theta,
                     int& errcode)
{
   ///////////////////////////////////////////////////////////////////////
   //  ALGORITHM AS 275 APPL.STATIST. (1992), VOL.41, NO.2
   //
   //  Computes the noncentral chi-square cumulative distribution function
   //  with positive real degrees of freedom df and nonnegative
   //  noncentrality parameter theta
   //
   //  errcode   = 0  no error
   //            = 1  parameter range error
   //            = 2  fail to converge
   //  Tianlin Wang, UIUC, Mon Jun  6 14:02:11 CDT 1994
   ////////////////////////////////////////////////////////////////////////
   if (df <= 0.0) throw exception(" ChiSquare_cdf(): bad arg2: %f, value > 0 expected");
   if (theta < 0.0) throw exception(" ChiSquare_cdf(): bad arg3, nonnegative value expected");
   errcode = 0;
   if (x <= 0.0) return 0.0;
   double chi2nc = x;
   if (theta == 0.0) return gammin(chi2nc*0.5, df*0.5);
   double lam,u,v,x2,f2,t,errbd,n2;
   lam = 0.5*theta;
   unsigned n = 1;         // Evaluate the first term
   u = exp(-lam);
   v = u;
   x2 = 0.5*x;
   f2 = 0.5*df;
   t = f2*log(x2) - x2 - gammln(f2+1.0);
   if (t <= -600.0) return 1.0;         // this is added in case exp(t) == 0
   t = exp(t);
   chi2nc = v*t;
   n2 = static_cast<double>(n + n);
   while ((df + n2) <= x) {
      u *= lam/static_cast<double>(n);       // Evaluate the next term of the expansion
      v += u;
      t *= x/(df + n2);
      chi2nc += v*t;
      n++;
      n2 = static_cast<double>(n + n);
   }
   do {
      n2 = static_cast<double>(n + n);
      errbd = t*x/(df + n2 -x);
      if (n>ITRMAX) throw exception ("ChiSquare_cdf(): fail to converge");
      u *= lam/static_cast<double>(n);       // Evaluate the next term of the expansion
      v += u;
      t *= x/(df + n2);
      chi2nc += v*t;
      n++;
   } while (errbd > ERRMAX);
   return chi2nc;
}

double t_cdf(const double t,const double df,const double delta)
{
   ///////////////////////////////////////////////////////////////////
   // ALGORITHM AS 243  APPL. STATIST. (1989), VOL.38, NO. 1
   // Cumulative probability at T of the non-central t-distribution
   // with df (positive real degrees of freedom,may be fractional) and
   // non-centrality parameter DELTA.
   //
   // REQUIREMENTS: gammln(x),betain(x,a,b,beta),errf(x);
   // Tianlin Wang, UIUC, Sun Jun  5 20:12:34 CDT 1994
   ////////////////////////////////////////////////////////////////////
   if (df <= 0.0) throw exception("t_cdf(): bad arg2, value > 0 expected");
   double tnc = 0.5;
   double lnbeta;
   if (delta == 0.0) {
      if (t==0.0) return tnc;
      lnbeta = 0.5*std::log(2.0*std::asin(1.0)) + gammln(0.5*df) - gammln(0.5*(1.0+df));
      tnc = betain(t*t/(t*t+df),0.5,0.5*df,lnbeta);
      if (t < 0.0) {
         tnc = 0.5 - 0.5*tnc;
      }
      else {
         tnc = 0.5 + 0.5*tnc;
      }
      return tnc;
   }
   double r2pi = std::sqrt(1.0/std::asin(1.0));     // 1.0/(Gamma(1.5)*sqrt(2))   // pi = 2.0*std::asin(1.0)
   double lnrpi = std::log(std::sqrt(2.0*std::asin(1.0)));
   double half = 0.5;
   double tt,del,en,x,lambda,p,q,s,a,b,rxb,xodd,godd,xeven,geven,errbd;
   tnc = 0.0;
   if (t<0.0) {
      tt = -t;
      del = -delta;
   }
   else {
      tt = t;
      del = delta;
   }
   // Initialize twin series(Guenther, J. Statis. Computn. Simuln. 6:199-, 1978
   en = 1.0;
   x = tt*tt/(tt*tt+df);
   if (x > 0.0) {
      lambda = del*del;
      p = half*exp(-half*lambda);
      q = r2pi*p*del;
      s = half - p;
      a = half;
      b = half*df;
      rxb = pow(1.0-x,b);
      lnbeta = lnrpi + gammln(b) - gammln(a+b);
      xodd = betain(x,a,b,lnbeta);
      godd = 2.0*rxb*exp(a*log(x)-lnbeta);
      xeven = 1.0 - rxb;
      geven = b*x*rxb;
      tnc = p*xodd + q*xeven;
      do {
         a += 1.0;
         xodd -= godd;
         xeven -= geven;
         godd *= x*(a+b-1.0)/a;
         geven *= x*(a+b-half)/(a+half);
         p *= lambda/(2.0*en);
         q *= lambda/(2.0*en+1.0);
         s -= p;
         en += 1.0;
         tnc += p*xodd + q*xeven;
         errbd = 2.0*s*(xodd-godd);
      } while(errbd > ERRMAX && en <= ITRMAX);
   }
   if (en <= ITRMAX) {
      tnc += 1.0 - Normal_cdf(del);
      if (t < 0.0) tnc = 1.0 - tnc;
   }
   return tnc;
}

double Beta_cdf(const double x,const double alfa,const double beta,
                   const double lambda,int& errcode)
{
   ////////////////////////////////////////////////////////////////////
   //  ALGORITHM AS226 APPL. STATIST. (1987) VOL. 36, NO. 2
   //  Incorporates modification AS R84 from AS vol. 39, pp311-2, 1990
   //
   //  Returns the cumulative probability of X for the non-central beta
   //  distribution with parameters A, B and non-centrality LAMBDA
   //
   //  Auxiliary routines required: ALOGAM - log-gamma function (ACM
   //  291 or AS 245), and BETAIN - incomplete-beta function (AS 63)
   //
   //  errcode   = 0  no error
   //            = 1  parameter range error
   //            = 2  fail to converge
   /////////////////////////////////////////////////////////////////////
   if (alfa <= 0.0 || beta <= 0.0) throw exception(" Beta_cdf(): bad args, value > 0 expected");
   errcode = 0;
   if (x <= 0.0) return 0.0;
   if (x >= 1.0) return 1.0;
   double lnbeta, betanc = 0.0;
   if (lambda == 0.0) {
      lnbeta = gammln(alfa) + gammln(beta) - gammln(alfa+beta);
      betanc = betain(x,alfa,beta,lnbeta);
      return betanc;
   }
   double c,x0,gx,q,xj,ax,sumq,temp,alfa0,errbd;
   c = lambda*0.5;
   x0 = 0.0;              // Initialize the series ...
   q = c - 5.0*sqrt(c);
   if (q > 0.0) x0 = static_cast<int>(q);
   alfa0 = alfa + x0;
   lnbeta = gammln(alfa0) + gammln(beta) - gammln(alfa0 + beta);
   temp = betain(x,alfa0,beta,lnbeta);
   gx = exp(alfa0*log(x) + beta*log(1.0 - x) - lnbeta - log(alfa0));
   if (alfa0 > alfa) {
      q = exp(-c + x0*log(c) - gammln(x0 + 1.0));
   }
   else {
      q = exp(-c);
   }
   xj = 0.0;
   ax = q*temp;
   sumq = 1.0 - q;
   betanc = ax;
   do {   // recur over subsequent terms until convergence is achieved...
      xj++;
      temp -= gx;
      gx *= x*(alfa+ beta + xj - 1.0)/(alfa + xj);
      q *= c/xj;
      sumq -= q;
      ax = temp*q;
      betanc += ax;
      errbd = (temp - gx)*sumq;
   } while(static_cast<int>(xj) < ITMAX && errbd > ERRMAX);
   if (errbd > ERRMAX) {
      warning("Beta_cdf(): fail to converge");
      errcode = 2;
   }
   if (betanc > 1.0) betanc = 1;
   return betanc;
}
#undef ITRMAX
#undef ERRMAX

double F_cdf(const double f,const double df1,const double df2,const double nc,
             int& errcode)
{
   if (df1 <= 0.0 || df2 <= 0.0) throw exception(" F_cdf(): bad args: value > 0 expected");
   errcode = 0;
   double retval;
   if (f <= 0.0) {
      retval = 0.0;
   } else {
      retval = Beta_cdf(df1*f/(df1*f+df2),df1*0.5,df2*0.5,nc,errcode);
   }
   return retval;
}

double Normal_inv(const double p)
{
   //////////////////////////////////////////////////////////////////
   // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
   //
   // Produces the standard normal deviate Z corresponding to a given
   // lower tail area of P; Z is accurate to about 1 part in 10**16.
   //////////////////////////////////////////////////////////////////
   if (p <= 0.0 || p >= 1.0 ) throw exception(" Normal_inv(): bad arg, value in (0,1) expected");
   if (p == 0.5) return 0.0;
   double *a = new double [8];
   double *b = new double [8];
   double *c = new double [8];
   double *d = new double [8];
   double *e = new double [8];
   double *f = new double [8];
   a[0] = 3.3871328727963666080e+0;   a[1] = 1.3314166789178437745e+2;
   a[2] = 1.9715909503065514427e+3;   a[3] = 1.3731693765509461125e+4;
   a[4] = 4.5921953931549871457e+4;   a[5] = 6.7265770927008700853e+4;
   a[6] = 3.3430575583588128105e+4;   a[7] = 2.5090809287301226727e+3;

   b[0] = 0.0;                        b[1] = 4.2313330701600911252e+1;
   b[2] = 6.8718700749205790830e+2;   b[3] = 5.3941960214247511077e+3;
   b[4] = 2.1213794301586595867e+4;   b[5] = 3.9307895800092710610e+4;
   b[6] =  2.8729085735721942674e+4;  b[7] = 5.2264952788528545610e+3;

   c[0] = 1.42343711074968357734e+0;  c[1] = 4.63033784615654529590e+0;
   c[2] = 5.76949722146069140550e+0;  c[3] = 3.64784832476320460504e+0;
   c[4] = 1.27045825245236838258e+0;  c[5] = 2.41780725177450611770e-1;
   c[6] = 2.27238449892691845833e-2;  c[7] = 7.74545014278341407640e-4;

   d[0] = 0.0;                        d[1] = 2.05319162663775882187e+0;
   d[2] = 1.67638483018380384940e+0;  d[3] = 6.89767334985100004550e-1;
   d[4] = 1.48103976427480074590e-1;  d[5] = 1.51986665636164571966e-2;
   d[6] = 5.47593808499534494600e-4;  d[7] = 1.05075007164441684324e-9;

   e[0] = 6.65790464350110377720e+0;  e[1] = 5.46378491116411436990e+0;
   e[2] = 1.78482653991729133580e+0;  e[3] = 2.96560571828504891230e-1;
   e[4] = 2.65321895265761230930e-2;  e[5] = 1.24266094738807843860e-3;
   e[6] = 2.71155556874348757815e-5;  e[7] = 2.01033439929228813265e-7;

   f[0] = 0.0;                        f[1] = 5.99832206555887937690e-1;
   f[2] = 1.36929880922735805310e-1;  f[3] = 1.48753612908506148525e-2;
   f[4] = 7.86869131145613259100e-4;  f[5] = 1.84631831751005468180e-5;
   f[6] = 1.42151175831644588870e-7;  f[7] = 2.04426310338993978564e-15;

   double split1 = 0.425;
   double split2 = 5.0;
   double const1 = 0.180625;
   double const2 = 1.6;

   double q,r,w1,w2,ppnd16;
   q = p - 0.5;
   if (fabs(q) <= split1) {
      r = const1 - q*q;
      w1 = ((((((a[7]*r+a[6])*r+a[5])*r+a[4])*r+a[3])*r+a[2])*r+a[1])*r+a[0];
      w2 = ((((((b[7]*r+b[6])*r+b[5])*r+b[4])*r+b[3])*r+b[2])*r+b[1])*r+1.0;
      ppnd16 = q*w1/w2;
   }
   else {
      if (q < 0.0) {r = p;}
      else         {r = 1.0-p;}
      if (r <= 0.0) {
	if(a){
	  delete [] a;  
	  a=0;
	}
	if(b){
	  delete [] b;  
	  b=0;
	}
	if(c){
	  delete [] c;
	  c=0;
	}
	if(d){
	  delete [] d;  
	  d=0;
	}
	if(e){
	  delete [] e;  
	  e=0;
	}
	if(f){
	  delete [] f;
	  f=0;
	}
         return 0.0;
      }
      r = sqrt(-log(r));
      if (r <= split2) {
	 r -= const2;
         w1= ((((((c[7]*r+c[6])*r+c[5])*r+c[4])*r+c[3])*r+c[2])*r+c[1])*r+c[0];
         w2= ((((((d[7]*r+d[6])*r+d[5])*r+d[4])*r+d[3])*r+d[2])*r+d[1])*r+1.0;
      }
      else {
         r -= split2;
         w1= ((((((e[7]*r+e[6])*r+e[5])*r+e[4])*r+e[3])*r+e[2])*r+e[1])*r+e[0];
         w2= ((((((f[7]*r+f[6])*r+f[5])*r+f[4])*r+f[3])*r+f[2])*r+f[1])*r+1.0;
      }
      ppnd16 = w1/w2;
      if (q< 0.0) ppnd16 = -ppnd16;
   }
   if(a){
     delete [] a;  
     a=0;
   }
   if(b){
     delete [] b;  
     b=0;
   }
   if(c){
     delete [] c;
     c=0;
   }
   if(d){
     delete [] d;  
     d=0;
   }
   if(e){
     delete [] e;  
     e=0;
   }
   if(f){
     delete [] f;
     f=0;
   }
   return ppnd16;
}

double ChiSquare_inv(const double p,const double df,const double nc)
{
   //////////////////////////////////////////////////////////////////
   //
   //        Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
   //
   //        To evaluate the percentage points of the chi-squared
   //        probability distribution function.
   //
   //     Incorporates the suggested changes in AS R85 (vol.40(1), 
   //     pp.233-5, 1991)
   //
   //     Auxiliary routines required: PPND = AS 111 (or AS 241) and
   //     GAMMAD = AS 239.
   // NOTE:
   //  ChiSquare_inv(1.0e-16,100,0.0) returns 52.4813, it seems not accurate
   ///////////////////////////////////////////////////////////////////////
   if (p < 0.0 || p >= 1.0 ) throw exception("ChiSquare_inv(): bad arg1, value in (0,1) expected");
   if (df <= 0.0) throw exception(" ChiSquare_inv(): bad arg2, value > 0 expected");

   if (p == 0.0) return 0.0;
   int errcode = 0;
   double ppchi2 = 0.0;   
   if (nc != 0.0) {
      double eps = 1.0e-7;
      double bot = 0.0;
      double top = 26000.0;
      double p2;
      while (top-bot > eps) {
         ppchi2 = (top + bot)*0.5;
         p2 = ChiSquare_cdf(ppchi2,df,nc,errcode);
         if (errcode != 0) break;
         if (p2 < p ) {
            bot = ppchi2;
         }
         else if (p2 > p) {
            top = ppchi2;
         }
         else {
            break;
         }
      }
      return ppchi2;
   }
   double one = 1.0;
   double aa,e,a,b,d,ch,p1,p2,q,s1,s2,s3,s4,s5,s6,t,x,xx;
   aa = 0.6931471806e+0;
   e = 0.5e-6;
   double* c = new double [38];
    c[0]= 0.01;    c[1]= 0.222222;  c[2]= 0.32;    c[3]= 0.4;     c[4]= 1.24;
    c[5]= 2.2;     c[6]= 4.67;      c[7]= 6.66;    c[8]= 6.73;    c[9]= 13.32;
   c[10]= 60.0;   c[11]= 70.0;     c[12]= 84.0;   c[13]= 105.0;  c[14]= 120.0;
   c[15]= 127.0;  c[16]= 140.0;    c[17]= 1175.0; c[18]= 210.0;  c[19]= 252.0;
   c[20]= 2264.0; c[21]= 294.0;    c[22]= 346.0;  c[23]= 420.0;  c[24]= 462.0;
   c[25]= 606.0;  c[26]= 672.0;    c[27]= 707.0;  c[28]= 735.0;  c[29]= 889.0;
   c[30]= 932.0;  c[31]= 966.0;    c[32]= 1141.0; c[33]= 1182.0; c[34]= 1278.0;
   c[35]= 1740.0; c[36]= 2520.0;   c[37]= 5040.0;

   xx = df*0.5;
   double g = gammln(xx);
   d = xx - one;
       //    starting approximation for small chi-squared
   if (df < -c[4]*log(p)) {
      ch = pow(p*xx*exp(g + xx*aa),one/xx);
      if (ch <= e) { 
	if(c){
	  delete [] c; 
	  c=0;
	}
	return ch;}
   }
   else {
      if (df < c[2]) {       // starting approximation for df <= 0.32
         ch = c[3];
         a = log(one-p);
         do {
            q = ch;
            p1 = one + ch*(c[6] + ch);
            p2 = ch*(c[8] + ch*(c[7] + ch));
            t = (c[6] + 2.0*ch)/p1 - (c[8] + ch*(c[9] + 3.0*ch))/p2 - 0.5;
            ch -= (1.0 - exp(a + g + 0.5*ch + d*aa)*p2/p1)/t;
         } while (fabs(q/ch - one) > c[0]);
      }
      else {
            //  call to algorithm AS 111 - note that p has been tested above.
            // AS 241 could be used as an alternative.
         x = Normal_inv(p);
           // starting approximation using Wilson and Hilferty estimate
         p1 = c[1]/df;
         ch = df*pow(x*sqrt(p1) + one - p1,3.0);
           // starting approximation for p tending to 1
         if (ch > c[5]*df + 6.0) ch = -2.0*(log(one-p) - d*log(ch*0.5) + g);
      }
   }
   int i = 0;
   while (i++ < 20) {
       // call to algorithm AS 239 and calculation of seven term Taylor series
      q = ch;
      p1 = 0.5*ch;
      p2 = p - gammin(p1,xx);
      t = p2*exp(xx*aa + g + p1 - d*log(ch));
      b = t/ch;
      a = 0.5*t - b*d;
      s1 = (c[18]+ a*(c[16]+ a*(c[13]+ a*(c[12]+ a*(c[11]+ c[10]*a)))))/c[23];
      s2 = (c[23]+ a*(c[28]+ a*(c[31]+ a*(c[32] + c[34]*a))))/c[36];
      s3 = (c[18]+ a*(c[24]+ a*(c[27]+ c[30]*a)))/c[36];
      s4 = (c[19]+ a*(c[26]+ c[33]*a)+ d*(c[21]+ a*(c[29]+ c[35]*a)))/c[37];
      s5 = (c[12]+ c[20]*a+ d*(c[17]+ c[25]*a))/c[36];
      s6 = (c[14]+ d*(c[22]+ c[15]*d))/c[37];
      ch += t*(one + 0.5*t*s1- b*d*(s1- b*(s2- b*(s3- b*(s4- b*(s5- b*s6))))));
      if (fabs(q/ch - one) > e) break;
   }
   ppchi2 = ch;
   if(c){
     delete [] c;
     c=0;
   }
   return ppchi2;
}

double t_expected_value(const double df, const double nc)
{
   double retval = 0.0;
   if (nc != 0.0) {
      double v = 0.5*df;
      retval = nc*exp(0.5*log(v)+gammln(v-0.5)-gammln(v));
   }
   return retval;
}

double t_inv(const double p,const double df, const double nc)
{
   ////////////////////////////////////////////////////////////
   // compute percentage point given lower tail area p.
   // binary iteration algoritm, a better algorithm is needed
   // Tianlin Wang, UIUC, Mon Jun  6 14:58:36 CDT 1994
   /////////////////////////////////////////////////////////
   if (p <= 0.0 || p >= 1.0) throw exception(" t_inv(): bad ar1, value in (0,1) expected");
   if (df <= 0.0) throw exception(" t_inv(): bad arg2, value > 0 expected");
   if (p == 0.5 && nc == 0.0) return 0.0;
   double eps = 1.0e-7;
   double p2,bot,top,ppt = t_expected_value(df,nc);
   bot = ppt - 700.0;
   top = ppt + 700.0;

   while (top-bot > eps) {
      ppt = (top + bot)*0.5;
      p2 = t_cdf(ppt,df,nc);
      if (p2 < p ) {
         bot = ppt;
      }
      else if (p2 > p) {
         top = ppt;
      }
      else {
         break;
      }
   }
   return ppt;
}

double F_inv(const double p,const double df1, const double df2,const double nc)
{
   //////////////////////////////////////////////////////////////
   // compute percentage point given lower tail area p.
   // binary iteration algoritm, a better algorithm is needed
   //
   // Tianlin Wang, UIUC, Mon Jun  6 14:58:36 CDT 1994
   ////////////////////////////////////////////////////////////
   if (p < 0.0 || p >= 1.0 ) throw exception(" F_inv(): bad arg1, value in (0,1) expected");
   if (df1 <= 0.0 || df2 <= 0.0 ) throw exception(" F_inv(): bad arg2 or 3, > 0 expected");
   if (p==0.0) return 0.0;
   double p2,ppf;
   if (df2 >= 3) {
      ppf = df2*(df1 + nc)/(df1*(df2 - 2.0));     // mean of F(df1,df2,nc)
   }
   else {
      ppf = 0;
   }
   double eps = 1.0e-7;
   double bot = 0.0;
   double top = ppf + 26000.0;
   int errcode = 0;
   while (top-bot > eps) {
      ppf = (top + bot)*0.5;
      p2 = F_cdf(ppf,df1,df2,nc,errcode);
      if (errcode != 0) break;
      if (p2 < p ) {
         bot = ppf;
      }
      else if (p2 > p) {
         top = ppf;
      }
      else {
         break;
      }
   }
   return ppf;
}
//BRS
/*************************************************************** 
The following are the random number generators taken from
   PEDSIM 2.0 from Mattis Schelling (matthias.schelling@inw.agrl.ethz.ch)
   via Welsey Petersen (wpp@scsc.ethz.ch) and netlib.
   zufall is long-lagged Fibonacci series generator uses lags of 273 and 607.
   Pedsim uses the following call:
double ran1()
{
  double a[1];
  zufall_(1, a);
  return a[0];
}
*/
// Start of the original code from W.P. Petersen, IPS, ETH Zurich
/* Common Block Declarations */
struct klotz0_1_ {
    double buff[607];
    int ptr;
    };
#define klotz0_1 (*(struct klotz0_1_ *) &klotz0_)
#define min(a,b) (a<b)?a:b
    
struct klotz1_1_ {
    double xbuff[1024];
    int first, xptr;
    };


#define klotz1_1 (*(struct klotz1_1_ *) &klotz1_)

/* Initialized data */
struct {
    int fill_1[1214];
    int e_2;
    } klotz0_ = { {0}, 0 };
struct {
    double fill_1[1024];
    int e_2[2];
    double e_3;
    } klotz1_ = { {0}, 0, 0, 0. };


/* Table of constant values */
int zufall_(int n, double *a)
{
    int buffsz = 607;
    int left, aptr, bptr, aptr0, i, k, q;
    double t;
    int nn, vl, qq, k273, k607, kptr=0;

/* portable lagged Fibonacci series uniform random number */
/* generator with "lags" -273 und -607: */
/* W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92 */

    aptr = 0;
    nn = n;
L1:
    if (nn <= 0) {
	return 0;
    }
/* factor nn = q*607 + r */
    q = (nn - 1) / 607;
    left = buffsz - klotz0_1.ptr;
    if (q <= 1) {
/* only one or fewer full segments */
	if (nn < left) {
            kptr = klotz0_1.ptr;
	    for (i = 0; i < nn; ++i) {
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    }
	    klotz0_1.ptr += nn;
	    return 0;
	} else {
            kptr = klotz0_1.ptr;
	    for (i = 0; i < left; ++i) {
		a[i + aptr] = klotz0_1.buff[kptr + i];
	    }
	    klotz0_1.ptr = 0;
	    aptr += left;
	    nn -= left;
/*  buff -> buff case */
	    vl = 273;
	    k273 = 334;
	    k607 = 0;
	    for (k = 0; k < 3; ++k) {
		for (i = 0; i < vl; ++i) {
		   t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		   klotz0_1.buff[k607+i] = t - static_cast<double>(static_cast<int>(t));
		}
		k607 += vl;
		k273 += vl;
		vl = 167;
		if (k == 0) {
		    k273 = 0;
		}
	    }
	    goto L1;
	}
    } else {
/* more than 1 full segment */
        kptr = klotz0_1.ptr;
	for (i = 0; i < left; ++i) {
	    a[i + aptr] = klotz0_1.buff[kptr + i];
	}
	nn -= left;
	klotz0_1.ptr = 0;
	aptr += left;
/* buff -> a(aptr0) */
	vl = 273;
	k273 = 334;
	k607 = 0;
	for (k = 0; k < 3; ++k) {
	    if (k == 0) {
		for (i = 0; i < vl; ++i) {
		    t = klotz0_1.buff[k273+i]+klotz0_1.buff[k607+i];
		    a[aptr + i] = t - static_cast<double>(static_cast<int>(t));
		}
		k273 = aptr;
		k607 += vl;
		aptr += vl;
		vl = 167;
	    } else {
		for (i = 0; i < vl; ++i) {
		    t = a[k273 + i] + klotz0_1.buff[k607 + i];
		    a[aptr + i] = t - static_cast<double>(static_cast<int>(t));
		}
		k607 += vl;
		k273 += vl;
		aptr += vl;
	    }
	}
	nn += -607;
/* a(aptr-607) -> a(aptr) for last of the q-1 segments */
	aptr0 = aptr - 607;
	vl = 607;
	for (qq = 0; qq < q-2; ++qq) {
	    k273 = aptr0 + 334;
	    for (i = 0; i < vl; ++i) {
		t = a[k273 + i] + a[aptr0 + i];
		a[aptr + i] = t - static_cast<double>(static_cast<int>(t));
	    }
	    nn += -607;
	    aptr += vl;
	    aptr0 += vl;
	}
/* a(aptr0) -> buff, last segment before residual */
	vl = 273;
	k273 = aptr0 + 334;
	k607 = aptr0;
	bptr = 0;
	for (k = 0; k < 3; ++k) {
	    if (k == 0) {
		for (i = 0; i < vl; ++i) {
		    t = a[k273 + i] + a[k607 + i];
		    klotz0_1.buff[bptr + i] = t - static_cast<double>(static_cast<int>(t));
		}
		k273 = 0;
		k607 += vl;
		bptr += vl;
		vl = 167;
	    } else {
		for (i = 0; i < vl; ++i) {
		    t = klotz0_1.buff[k273 + i] + a[k607 + i];
		    klotz0_1.buff[bptr + i] = t - static_cast<double>(static_cast<int>(t));
		}
		k607 += vl;
		k273 += vl;
		bptr += vl;
	    }
	}
	goto L1;
    }
} /* zufall_ */


int zufalli_(int seed) {
    /* Initialized data */
    int kl = 9373;
    int ij = 1802;

    /* Local variables */
    int i, j, k, l, m;
    double s, t;
    int ii, jj;

/*  generates initial seed buffer by linear congruential */
/*  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50 */
/*  variable seed should be 0 < seed <31328 */

    if (seed != 0) {
	ij = seed;
    }
    i = ij / 177 % 177 + 2;
    j = ij % 177 + 2;
    k = kl / 169 % 178 + 1;
    l = kl % 169;
    for (ii = 0; ii < 607; ++ii) {
	s = 0.;
	t = .5;
	for (jj = 1; jj <= 24; ++jj) {
	    m = i * j % 179 * k % 179;
	    i = j;
	    j = k;
	    k = m;
	    l = (l * 53 + 1) % 169;
	    if (l * m % 64 >= 32) {
		s += t;
	    }
	    t *= 0.5;
	}
	klotz0_1.buff[ii] = s;
    }
    return 0;
} /* zufalli_ */
//BRS
} /////// end of namespace

