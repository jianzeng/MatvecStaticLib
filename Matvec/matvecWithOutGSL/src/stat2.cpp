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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "exception.h"
#include "session.h"


namespace matvec {

extern void warning(const char format[],...);
extern double ran1(void);

long RAND_STATE_ARRAY[64] = {
	3,
	0x9a319039, 0x32d9c024, 0x9b663182, 0x5da1f342,
	0x7449e56b, 0xbeb1dbb0, 0xab5c5918, 0x946554fd,
	0x8c2e680f, 0xeb3d799f, 0xb11ee0b7, 0x2d436b86,
	0xda672e2a, 0x1588ca88, 0xe369735d, 0x904f35f7,
	0xd7158fd6, 0x6fa6f051, 0x616e6b96, 0xac94efdc,
	0xde3b81e0, 0xdf0a6fb5, 0xf103bc02, 0x48f340fb,
	0x36413f93, 0xc622c298, 0xf5a42ab8, 0x8a88d77b,
	0xf5ad9d0e, 0x8999220b, 0x27fb47b9,
	9,
	0x1b3a5678, 0xc54afd63, 0x8374b4db, 0x4345fd75,
	0x36da4534, 0x27dcdd45, 0xbc630efc, 0x488776ed,
	0xfdd76b12, 0x3452ecf6, 0x9876fecd, 0xa3d644e5,
	0x4567debf, 0x7462dd89, 0xdfe5321b, 0x7ed2ad4a,
	0x63adcdc0, 0x36abd051, 0x6ec56d96, 0x836eeefc,
	0xd83b9de0, 0x3fe88cb5, 0x5bef43f8, 0x5d37d5e3,
	0x35dc66aa, 0x36a74f4c, 0xed5c5548, 0x6eae3922,
	0x53b3ae21, 0x88e3b48b, 0x6c20f1a9};

void set_seed(const unsigned seed)      //  set seed for rand()
{
   srand(seed);
}

double fsign( double num, double sign )
               // Transfers sign of argument sign to argument num
{
   if ( ( sign>0.0 && num<0.0 ) || ( sign<0.0 && num>0.0 ) ) return ( -num );
   else return ( num );
}

double ranf(void)
{
  // wrapper for current uniform random number generator (0,1) 
  // RLF 9/20/02

  double u = SESSION.mtr.rand();
  while (u == 0 || u==1) {
    u = SESSION.mtr.rand();
  }
   return u;
}

double snorm(void)
// ***************************************************************
//
//     (STANDARD-)  N O R M A L  DISTRIBUTION
//
//     FOR DETAILS SEE:
//               AHRENS, J.H. AND DIETER, U.
//               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
//               SAMPLING FROM THE NORMAL DISTRIBUTION.
//               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.
//
//     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
//     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)
//
//     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
//     SUNIF.  The argument IR thus goes away.
//
//     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
//     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
//*****************************************************************

{
static float a[32] = {
    0.0,3.917609e-2,7.841241e-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static float d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static float t[31] = {
    7.673828e-4,2.30687e-3,3.860618e-3,5.438454e-3,7.0507e-3,8.708396e-3,
    1.042357e-2,1.220953e-2,1.408125e-2,1.605579e-2,1.81529e-2,2.039573e-2,
    2.281177e-2,2.543407e-2,2.830296e-2,3.146822e-2,3.499233e-2,3.895483e-2,
    4.345878e-2,4.864035e-2,5.468334e-2,6.184222e-2,7.047983e-2,8.113195e-2,
    9.462444e-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static float h[31] = {
    3.920617e-2,3.932705e-2,3.951e-2,3.975703e-2,4.007093e-2,4.045533e-2,
    4.091481e-2,4.145507e-2,4.208311e-2,4.280748e-2,4.363863e-2,4.458932e-2,
    4.567523e-2,4.691571e-2,4.833487e-2,4.996298e-2,5.183859e-2,5.401138e-2,
    5.654656e-2,5.95313e-2,6.308489e-2,6.737503e-2,7.264544e-2,7.926471e-2,
    8.781922e-2,9.930398e-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
    static long i;
    static double u,s,ustar,aa,w,y,tt;
    double  haha;

    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
        //////////////////  START CENTER  //////////////////////
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-(*(t+i-1)))*(*(h+i-1));
S50:
       //////////////////  EXIT   (BOTH CASES) ///////////////////
    y = aa+w;
    haha = y;
    if(s == 1.0) haha = -y;
    return haha;
S60:
      //////////////////   CENTER CONTINUED   ///////////////////
    u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = ranf();
S80:
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
    //////////////////    START TAIL    ///////////////////
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u*(*(d+i-1));
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

double sexpo(void)
//*****************************************************************
//
//     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION
//
//     FOR DETAILS SEE:
//
//               AHRENS, J.H. AND DIETER, U.
//               COMPUTER METHODS FOR SAMPLING FROM THE
//               EXPONENTIAL AND NORMAL DISTRIBUTIONS.
//               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.
//
//     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM
//     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)
//
//     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
//     SUNIF.  The argument IR thus goes away.
//
//     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
//     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
// ******************************************************************
{
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,1.0
};
	static long i;
	static double *q1 = q;
	static double a,u,ustar,umin;
	double haha;

    a = 0.0;
    u = ranf();
    goto S30;
S20:
    a += *q1;
S30:
    u += u;
    if(u <= 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    haha = a+u;
    return haha;
S60:
    i = 1;
    ustar = ranf();
    umin = ustar;
S70:
    ustar = ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    haha = a+umin*(*q1);
    return haha;
}

double sgamma(const double a)
//*******************************************************************
//
//     (STANDARD-)  G A M M A  DISTRIBUTION
//               PARAMETER  A >= 1.0  !
//
//     FOR DETAILS SEE:
//               AHRENS, J.H. AND DIETER, U.
//               GENERATING GAMMA VARIATES BY A
//               MODIFIED REJECTION TECHNIQUE.
//               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.
//
//     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER
//                                 (STRAIGHTFORWARD IMPLEMENTATION)
//
//     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
//     SUNIF.  The argument IR thus goes away.
//
//               PARAMETER  0.0 < A < 1.0 !
//
//
//     FOR DETAILS SEE:
//
//               AHRENS, J.H. AND DIETER, U.
//               COMPUTER METHODS FOR SAMPLING FROM GAMMA,
//               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.
//               COMPUTING, 12 (1974), 223 - 246.
//
//     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)
//
//     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
//     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
//     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
//     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
//     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
//     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
//     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
//
//  T Wang made lots of changes to avoid all labels statements
//*******************************************************************

{
   static double q1 = 4.166669E-2;
   static double q2 = 2.083148E-2;
   static double q3 = 8.01191E-3;
   static double q4 = 1.44121E-3;
   static double q5 = -7.388E-5;
   static double q6 = 2.4511E-4;
   static double q7 = 2.424E-4;
   static double a1 = 0.3333333;
   static double a2 = -0.250003;
   static double a3 = 0.2000062;
   static double a4 = -0.1662921;
   static double a5 = 0.1423657;
   static double a6 = -0.1367177;
   static double a7 = 0.1233795;
   static double e1 = 1.0;
   static double e2 = 0.4999897;
   static double e3 = 0.166829;
   static double e4 = 4.07753E-2;
   static double e5 = 1.0293E-2;
   static double aa = 0.0;
   static double aaa = 0.0;
   static double sqrt32 = 5.65685424949;
   static double s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;
   double haha =0.0;
	
   if (a != aa) {
      if (a < 1.0) {
           // ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
         aa = 0.0;
         b = 1.0+0.3678794*a;
         for (;;) {
            p = b*ranf();
            if (p >= 1.0) {
               haha = -log((b-p)/ a);
               if(sexpo() >= (1.0-a)*log(haha)) return haha;
            }
            else {
               haha = exp(log(p)/ a);
               if (sexpo() >= haha) return haha;
            }
         }
      }
      else {
         // STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED  //
         aa = a;
         s2 = a-0.5;
         s = sqrt(s2);
         d = sqrt32-12.0*s;
      }
   }
   /////////////////////////////////////////////////////
   //     STEP  2:  T=STANDARD NORMAL DEVIATE,
   //               X=(S,1/2)-NORMAL DEVIATE.
   //               IMMEDIATE ACCEPTANCE (I)
   ////////////////////////////////////////////////////
   t = snorm();
   x = s+0.5*t;
   haha = x*x;
   if (t >= 0.0) return haha;
      //  STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
   u = ranf();
   if (d*u <= t*t*t) return haha;
      //  STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
   if (a != aaa) {
      aaa = a;
      r = 1.0/a;
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
      /////////////////////////////////////////////////////////
      //   APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
      //   THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
      //   C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
      //////////////////////////////////////////////////////////
      if (a > 3.686) {
         if (a > 13.022) {  //  CASE 3:  A .GT. 13.022
            b = 1.77;
            si = 0.75;
            c = 0.1515/s;
         }
         else {            //  CASE 2:  3.686 .LT. A .LE. 13.022
            b = 1.654+0.0076*s2;
            si = 1.68/s+0.275;
            c = 0.062/s+0.024;
         }
      }
      else {               //   CASE 1:  A .LE. 3.686
         b = 0.463+s-0.178*s2;
         si = 1.235;
         c = 0.195/s-0.079+0.016*s;
      }
   }
   //  STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
   if (x > 0.0) {
      // STEP  6:  CALCULATION OF V AND QUOTIENT Q
      v = t/(s+s);
      if (fabs(v) > 0.25) {
         q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
      }
      else {
         q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
      }
      //  STEP  7:  QUOTIENT ACCEPTANCE (Q)
      if(log(1.0-u) <= q) return haha;
   }
   /////////////////////////////////////////////////////////////
   //  STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
   //            U= 0,1 -UNIFORM DEVIATE
   //            T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
   ///////////////////////////////////////////////////////////////
   for (;;) {
      e = sexpo();
      u = ranf();
      u += (u-1.0);
      t = b+fsign(si*e,u);
      //  STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
      if (t >= -0.7187449) {
         //  STEP 10:  CALCULATION OF V AND QUOTIENT Q   //
         v = t/(s+s);
         if (fabs(v) > 0.25) {
            q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
         }
         else {
            q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
         }
         // STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8) //
         if (q > 0.0) {
            if (q > 0.5) { w = exp(q)-1.0;}
            else         { w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q; }
            // IF T IS REJECTED, SAMPLE AGAIN AT STEP 8  //
            if (c*fabs(u) <= w*exp(e-0.5*t*t)) {
               x = s+0.5*t;
               return x*x;
            }
         }
      }
   }
}

double gennor(const double av,const double sd)
//*********************************************************************
//     double gennor(double av,double sd)
//         GENerate random deviate from a NORmal distribution
//                              Function
//     Generates a single random deviate from a normal distribution
//     with mean, AV, and standard deviation, SD.
//                              Arguments
//     av --> Mean of the normal distribution.
//     sd --> Standard deviation of the normal distribution.
//                              Method
//     Renames SNORM from TOMS as slightly modified by BWB to use RANF
//     instead of SUNIF.
//     For details see:
//               Ahrens, J.H. and Dieter, U.
//               Extensions of Forsythe's Method for Random
//              Sampling from the Normal Distribution.
//               Math. Comput., 27,124 (Oct. 1973), 927 - 937.
//********************************************************************
{
    return ( sd*snorm()+av );
}

double gengam(const double a, const double r)

//********************************************************************
//     double gengam(double a,double r)
//           GENerates random deviates from GAMma distribution
//                              Function
//     Generates random deviates from the gamma distribution whose
//     density is
//                 1
//          ----------------* X^(A-1) * Exp(-X/R)
//          (R^A)*Gamma(A)
//                              Arguments
//     R --> Location parameter of Gamma distribution, R>0.0
//     A --> Shape parameter of Gamma distribution,A>=0.0
//                              Method
//     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
//     instead of SUNIF.
//     For details see:
//               (Case A >= 1.0)
//               Ahrens, J.H. and Dieter, U.
//               Generating Gamma Variates by a
//               Modified Rejection Technique.
//               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
//     Algorithm GD
//               (Case 0.0 <= A <= 1.0)
//               Ahrens, J.H. and Dieter, U.
//               Computer Methods for Sampling from Gamma,
//               Beta, Poisson and Binomial Distributions.
//               Computing, 12 (1974), 223-246/
//     Adapted algorithm GS.
//*****************************************************************

{
    return ( sgamma(a)*r );
}

double genchi(const double df)

//****************************************************************
//     double genchi(double df)
//                Generate random value of CHIsquare variable
//                              Function
//     Generates random deviate from the distribution of a chisquare
//     with DF degrees of freedom random variable.
//                              Arguments
//     df --> Degrees of freedom of the chisquare
//            (Must be positive)
//
//                              Method
//     Uses relation between chisquare and gamma.
//****************************************************************
{

   if (df <= 0) throw exception(" genchi(): invalid arg");
   return (2.0*sgamma(df/2.0));
}

double genexp(const double av)

//*********************************************************
//     double genexp(double av)
//                    GENerate EXPonential random deviate
//                              Function
//     Generates a single random deviate from an exponential
//     distribution with mean AV.
//                              Arguments
//     av --> The mean of the exponential distribution from which
//            a random deviate is to be generated.
//                              Method
//     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
//     instead of SUNIF.
//     For details see:
//               Ahrens, J.H. and Dieter, U.
//               Computer Methods for Sampling From the
//               Exponential and Normal Distributions.
//               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
//************************************************************

{

    return ( sexpo()*av );
}

double genf(const int dfn,const int dfd)

//***********************************************************
//     double genf(double dfn,double dfd)
//                GENerate random deviate from the F distribution
//                              Function
//     Generates a random deviate from the F (variance ratio)
//     distribution with DFN degrees of freedom in the numerator
//     and DFD degrees of freedom in the denominator.
//                              Arguments
//     dfn --> Numerator degrees of freedom
//             (Must be positive)
//     dfd --> Denominator degrees of freedom
//             (Must be positive)
//                              Method
//     Directly generates ratio of chisquare variates
//************************************************************

{
   if ( dfn <= 0.0 || dfd <= 0.0) throw exception(" genf(): invalid args");
   static double genfval = 0.0;
   static double xden,xnum;
   xnum = genchi(dfn)/dfn;
          //////  GENF = ( GENCHI(DFN)/DFN )/( GENCHI(DFD)/DFD ) ///////
   xden = genchi(dfd)/dfd;
   if ( xden <= 9.999999999998e-39*xnum) {
      warning("GENF - generated numbers would cause overflow,"
                 " df1 = %16.6E, df2 = %16.6E",xnum,xden);
      genfval = 1.0E38;
   }
   else {
      genfval = xnum/xden;
   }
   return genfval;
}

double gennch(const int df, const double xnonc)

//************************************************************************
//     double gennch(int df,double xnonc)
//           Generate random value of Noncentral CHIsquare variable
//                              Function
//     Generates random deviate  from the  distribution  of a  noncentral
//     chisquare with DF degrees  of freedom and noncentrality  parameter
//     xnonc.
//                              Arguments
//     df --> Degrees of freedom of the chisquare
//           (Must be > 1.0)
//     xnonc --> Noncentrality parameter of the chisquare
//               (Must be >= 0.0)
//                              Method
//     Uses fact that  noncentral chisquare  is  the  sum of a  chisquare
//     deviate with DF-1  degrees of freedom plus the  square of a normal
//     deviate with mean XNONC and standard deviation 1.
//************************************************************************
{
   if ( df > 1 && xnonc >= 0.0) {
      return ( genchi(df-1)+pow(gennor(sqrt(xnonc),1.0),2.0));
   }
   else if (df == 1 && xnonc >= 0.0) {
      double x = snorm() + xnonc;
      return x*x;
   }
   else {
      throw exception(" gennch(): invalid arg");
   }
}

double gennf(const int dfn,const int dfd,const double xnonc)
//*********************************************************************
//     double gennf(int dfn,int dfd,double xnonc)
//           GENerate random deviate from the Noncentral F distribution
//                              Function
//     Generates a random deviate from the  noncentral F (variance ratio)
//     distribution with DFN degrees of freedom in the numerator, and DFD
//     degrees of freedom in the denominator, and noncentrality parameter
//     XNONC.
//                              Arguments
//     dfn --> Numerator degrees of freedom
//             (Must be >= 1)
//     dfd --> Denominator degrees of freedom
//             (Must be positive)
//     xnonc --> Noncentrality parameter
//               (Must be nonnegative)
//                              Method
//    Directly generates ratio of noncentral numerator chisquare variate
//     to central denominator chisquare variate.
//***********************************************************************
{
   if (dfn <= 1 || dfd <= 0 || xnonc < 0.0) throw exception(" gennf(): invalid args");
   static double gennf,xden,xnum;
   xnum = gennch(dfn,xnonc)/dfn;
       /////   GENNF = ( GENNCH(DFN,XNONC)/DFN )/( GENCHI(DFD)/DFD ) /////
   xden = genchi(dfd)/dfd;
   if ( xden <= 9.999999999998e-39*xnum ) {
      warning("gennf(): overflow, df1= %16.6E, df2= %16.6E",xnum,xden);
      gennf = 1.0E38;
   }
   else {
      gennf = xnum/xden;
   }
   return gennf;
}

double genunf(const double low,const double high)
//***********************************************************************
//     double genunf(double low,double high)
//               GeNerate Uniform Real between LOW and HIGH
//                              Function
//     Generates a real uniformly distributed between LOW and HIGH.
//                              Arguments
//     low --> Low bound (exclusive) on real value to be generated
//     high --> High bound (exclusive) on real value to be generated
//***********************************************************************
{
   if (low > high ) throw exception("genunf(): invalid arg");
   return ( low+(high-low)*ranf() );
}

long ignbin(const long n,const double pp)
//****************************************************************
//     long ignbin(long n,double pp)
//                    GENerate BINomial random deviate
//                              Function
//     Generates a single random deviate from a binomial
//     distribution whose number of trials is N and whose
//     probability of an event in each trial is P.
//                              Arguments
//     n  --> The number of trials in the binomial distribution
//            from which a random deviate is to be generated.
//     p  --> The probability of an event in each trial of the
//            binomial distribution from which a random deviate
//            is to be generated.
//     ignbin <-- A random deviate yielding the number of events
//                from N independent trials, each of which has
//                a probability of event P.
//                              Method
//     This is algorithm BTPE from:
//         Kachitvichyanukul, V. and Schmeiser, B. W.
//         Binomial Random Variate Generation.
//         Communications of the ACM, 31, 2
//         (February, 1988) 216.
//
//     SUBROUTINE BTPEC(N,PP,ISEED,JX)
//     BINOMIAL RANDOM VARIATE GENERATOR
//     MEAN .LT. 30 -- INVERSE CDF
//       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
//       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
//       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
//       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
//     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
//     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
//       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
//       USABLE ALGORITHM.
//     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
//       "BINOMIAL RANDOM VARIATE GENERATION,"
//       COMMUNICATIONS OF THE ACM, FORTHCOMING
//     WRITTEN:  SEPTEMBER 1980.
//       LAST REVISED:  MAY 1985, JULY 1987
//     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
//                           GENERATOR
//     ARGUMENTS
//       N : NUMBER OF BERNOULLI TRIALS            (INPUT)
//       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
//       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
//       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
//     VARIABLES
//      PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
//       NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
//       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
//       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
//       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
//       M:  INTEGER VALUE OF THE CURRENT MODE
//       FM:  FLOATING POINT VALUE OF THE CURRENT MODE
//       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
//       P1:  AREA OF THE TRIANGLE
//       C:  HEIGHT OF THE PARALLELOGRAMS
//       XM:  CENTER OF THE TRIANGLE
//       XL:  LEFT END OF THE TRIANGLE
//       XR:  RIGHT END OF THE TRIANGLE
//       AL:  TEMPORARY VARIABLE
//       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
//       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
//       P2:  AREA OF THE PARALLELOGRAMS
//       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
//       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
//       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
//           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
//           FROM THE REGION
//       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
//           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
//           REJECT THE CANDIDATE VALUE
//       IX:  INTEGER CANDIDATE VALUE
//       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
//           AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
//       K:  ABSOLUTE VALUE OF (IX-M)
//       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
//           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
//           ALSO USED IN THE INVERSE TRANSFORMATION
//       R: THE RATIO P/Q
//      G: CONSTANT USED IN CALCULATION OF PROBABILITY
//       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
//            OF F WHEN IX IS GREATER THAN M
//       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
//             CALCULATION OF F WHEN IX IS LESS THAN M
//       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
//       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
//       YNORM: LOGARITHM OF NORMAL BOUND
//       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
//       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
//       USED IN THE FINAL ACCEPT/REJECT TEST
//       QN: PROBABILITY OF NO SUCCESS IN N TRIALS
//     REMARK
//       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
//       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
//       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
//       ARE NOT INVOLVED.
//     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
//     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
//     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
//    *DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY*
//************************************************************************
{
  if (pp < 0.0 || pp>1.0) throw exception("ignbin(): bad arg2, value between (0,1) expected");
   static double psave = -1.0;
   static long nsave = -1;
   static long ignbin,i,ix,ix1,k,m,mp,T1;
   static double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r;
   static double u,v,w,w2,x,x1,x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

   if(pp != psave) goto S10;
   if(n != nsave) goto S20;
   if(xnp < 30.0) goto S150;
   goto S30;
S10:  /////// SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE //////

   psave = pp;
   p = std::min(psave,1.0-psave);
   q = 1.0-p;
S20:
   xnp = n*p;
   nsave = n;
   if(xnp < 30.0) goto S140;
   ffm = xnp+p;
   m = long(ffm);
   fm = m;
   xnpq = xnp*q;
   p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
   xm = fm+0.5;
   xl = xm-p1;
   xr = xm+p1;
   c = 0.134+20.5/(15.3+fm);
   al = (ffm-xl)/(ffm-xl*p);
   xll = al*(1.0+0.5*al);
   al = (xr-ffm)/(xr*q);
   xlr = al*(1.0+0.5*al);
   p2 = p1*(1.0+c+c);
   p3 = p2+c/xll;
   p4 = p3+c/xlr;
S30:   ////////  GENERATE VARIATE  ////////

   u = ranf()*p4;
   v = ranf();

   ///////  TRIANGULAR REGION  ///////

   if(u > p1) goto S40;
   ix = long(xm-p1*v+u);
   goto S170;
S40:     /////  PARALLELOGRAM REGION  /////
   if(u > p2) goto S50;
   x = xl+(u-p1)/c;
   v = v*c+1.0-fabs(xm-x)/p1;
   if(v > 1.0 || v <= 0.0) goto S30;
   ix = long(x);
   goto S70;
S50:       //////  LEFT TAIL  ///////
   if(u > p3) goto S60;
   ix = long(xl+log(v)/xll);
   if(ix < 0) goto S30;
   v *= ((u-p2)*xll);
   goto S70;
S60:     ////  RIGHT TAIL  ////
   ix = long(xr-log(v)/xlr);
   if(ix > n) goto S30;
   v *= ((u-p3)*xlr);
S70:   //// DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST  ////
   k = abs((int)(ix-m));
   if(k > 20 && k < xnpq/2-1) goto S130;

     ///////  EXPLICIT EVALUATION   //////

   f = 1.0;
   r = p/q;
   g = (n+1)*r;
   T1 = m-ix;
   if(T1 < 0) {
      goto S80;
   }
   else if(T1 == 0) {
      goto S120;
   }
   else {
      goto S100;
   }
S80:
   mp = m+1;
   for(i=mp; i<=ix; i++) f *= (g/i-r);
   goto S120;
S100:
   ix1 = ix+1;
   for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
   if(v <= f) goto S170;
   goto S30;
S130:   /////   SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))  ////
   amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
   ynorm = -(k*k/(2.0*xnpq));
   alv = log(v);
   if(alv < ynorm-amaxp) goto S170;
   if(alv > ynorm+amaxp) goto S30;
    ///////////////////////////////////////////////////
    //  STIRLING'S FORMULA TO MACHINE ACCURACY FOR
    // THE FINAL ACCEPTANCE/REJECTION TEST
    ////////////////////////////////////////////////////
   x1 = ix+1.0;
   f1 = fm+1.0;
   z = n+1.0-fm;
   w = n-ix+1.0;
   z2 = z*z;
   x2 = x1*x1;
   f2 = f1*f1;
   w2 = w*w;
   if (alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
     (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
   goto S30;
S140:   ////  INVERSE CDF LOGIC FOR MEAN LESS THAN 30 ////
   qn = pow(q,(double)n);
   r = p/q;
   g = r*(n+1);
S150:
   ix = 0;
   f = qn;
   u = ranf();
S160:
   if(u < f) goto S170;
   if(ix > 110) goto S150;
   u -= f;
   ix += 1;
   f *= (g/ix-r);
   goto S160;
S170:
   if(psave > 0.5) ix = n-ix;
   ignbin = ix;
   return ignbin;
}

long ignpoi(const double mu)
{
   /////////////////////////////////////////////////////////////////////////
   // all labels stattements have been replaced with non-labeled statements
   // T Wang, 4/27/94
   /////////////////////////////////////////////////////////////////////////
   static double a0 = -0.5;
   static double a1 = 0.3333333;
   static double a2 = -0.2500068;
   static double a3 = 0.2000118;
   static double a4 = -0.1661269;
   static double a5 = 0.1421878;
   static double a6 = -0.1384794;
   static double a7 = 0.125006;
   static double muold = 0.0;
   static double muprev = 0.0;
   static double fact[10] = {
    1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
   };
   static long ignpoi,j,k,kflag,l,m;
   static double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g;
   static double omega,p,p0,px,py,q,s,t,u,v,x,xx,pp[35];
   int wflag;    // added by T. Wang

   if (mu != muprev) {
      if (mu < 10.0) {
         //   C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
         muprev = 0.0;
         if (mu != muold) {
            muold = mu;
            m = std::max(1L,(long) (mu));
            l = 0;
            p = exp(-mu);
            q = p0 = p;
         }
         //  STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
         for (;;) {
            u = ranf();
            ignpoi = 0;
            if (u <= p0) return ignpoi;
            //////////////////////////////////////////////////////////
            // STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
            //         PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
            //         (0.458=PP(9) FOR MU=10)
            //////////////////////////////////////////////////////////
            if (l == 0) {
               ////////////////////////////////////////////////////////
               //  STEP C. CREATION OF NEW POISSON PROBABILITIES P
               //          AND THEIR CUMULATIVES Q=PP(K)
               /////////////////////////////////////////////////////////
               l += 1;
               for (k=l; k<=35; k++) {
                  p = p*mu/(double)k;
                  q += p;
                  *(pp+k-1) = q;
                  if(u <= q) {l = ignpoi = k; return ignpoi;}
               }
               l = 35;
            }
            else {
               j = 1;
               if (u > 0.458) j = std::min(l,m);
               for (k=j; k<=l; k++) {
                  if (u <= *(pp+k-1)) {ignpoi = k; return ignpoi;}
               }
               if (l != 35) {
                  l += 1;
                  for (k=l; k<=35; k++) {
                     p = p*mu/(double)k;
                     q += p;
                     *(pp+k-1) = q;
                     if(u <= q) {l = ignpoi = k; return ignpoi;}
                  }
                  l = 35;
               }
            }
         }
      }
      else {
         ////  C A S E  A. (RECALCULATION OF S,D,L IF MU HAS CHANGED) ////
         muprev = mu;
         s = sqrt(mu);
         d = 6.0*mu*mu;
          /////////////////////////////////////////////////////////////
          //   THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
         //      PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
         //      IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
         ///////////////////////////////////////////////////////////////
         l = (long) (mu-1.1484);
      }
   }
   ////   STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE ///
   g = mu+s*snorm();
   if (g >= 0.0) {
      ignpoi = (long) (g);
          //  STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH  //
      if (ignpoi >= l) return ignpoi;
          // STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U  //
      fk = (double)ignpoi;
      difmuk = mu-fk;
      u = ranf();
      if (d*u >= difmuk*difmuk*difmuk) return ignpoi;
   }
   /////////////////////////////////////////////////////////////////////
   //  STEP P. PREPARATIONS FOR STEPS Q AND H.
   //          (RECALCULATIONS OF PARAMETERS IF NECESSARY)
   //          .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
   //          THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
   //          APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
   //          C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
   ///////////////////////////////////////////////////////////////////////
   if (mu != muold) {
      muold = mu;
      omega = 0.3989423/s;
      b1 = 4.166667E-2/mu;
      b2 = 0.3*b1*b1;
      c3 = 0.1428571*b1*b2;
      c2 = b2-15.0*c3;
      c1 = b1-6.0*b2+45.0*c3;
      c0 = 1.0-b1+3.0*b2-15.0*c3;
      c = 0.1069/mu;
   }
   wflag = -1;
   if (g >= 0.0) {  // 'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
      kflag = 0;
      wflag = 0;
   }
   else {
      wflag = 1;
   }

   for (;;) {          //  STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
      if ( wflag > 0) {
         do {
            ////////////////////////////////////////////////////////////////
            //  STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
            //          DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
            //          (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
            ////////////////////////////////////////////////////////////////
            e = sexpo();
            u = ranf();
            u += (u-1.0);
            t = 1.8+fsign(e,u);
         } while (t <= -0.6744);
         ignpoi = (long) (mu+s*t);
         fk = (double)ignpoi;
         difmuk = mu-fk;
         //  'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
         kflag = 1;
      }
      else {
         wflag = 1;
      }
      for (;;) {
         if (ignpoi < 10) {
            px = -mu;
            py = pow(mu,(double)ignpoi)/ *(fact+ignpoi);
            x = (0.5-difmuk)/s;
            xx = x*x;
            fx = -0.5*xx;
            fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
            if (kflag <= 0) {
               if(fy-u*fy <= py*exp(px-fx)) return ignpoi;
               break;
            }
         }
         else {
            ////////////////////////////////////////////////////////////////
            //        CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
            //        A0-A7 FOR ACCURACY WHEN ADVISABLE
            //        .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
            ////////////////////////////////////////////////////////////////
            del = 8.333333E-2/fk;
            del -= (4.8*del*del*del);
            v = difmuk/fk;
            if (fabs(v) > 0.25) {
               px = fk*log(1.0+v)-difmuk-del;
            }
            else {
               px =
                fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
            }
            py = 0.3989423/sqrt(fk);
             x = (0.5-difmuk)/s;
            xx = x*x;
            fx = -0.5*xx;
            fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
            if (kflag <= 0) {
               if(fy-u*fy <= py*exp(px-fx)) return ignpoi;
               break;
            }
         }
               // STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
         if  (c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) {break;}
         else                                        {return ignpoi;}
      }
   }
}

long ignuin(const long low,const long high)
//********************************************************************
//     long ignuin(long low,long high)
//               GeNerate Uniform INteger
//                              Function
//     Generates an integer uniformly distributed between LOW and HIGH.
//                              Arguments
//     low --> Low bound (inclusive) on integer value to be generated
//     high --> High bound (inclusive) on integer value to be generated
//                              Note
//     If (HIGH-LOW) > 2,147,483,561 prints error message on * unit and
//     stops the program.
//
//     IGNLGI generates integers between 1 and 2^31(=2147483648)
//     MAXNUM is 1 less than maximum generable value
//********************************************************************
{
#define maxnum 2147483647L
   if (low > high) throw exception(" ignu(low,high): low > hig");
//    long random();
   static long ignuin,ign,maxnow,range,ranp1;
   range = high-low;
   if (range > maxnum) throw exception(" ignuin(low,high): high too large");
   if (low == high) {
      ignuin = low;
      return ignuin;
   }
   ////////////////////////////////////////////////////////////////////
   //  Number to be generated should be in range 0..RANGE
   //  Set MAXNOW so that the number of integers in 0..MAXNOW is an
   //  integral multiple of the number in 0..RANGE
   ////////////////////////////////////////////////////////////////////
   ranp1 = range+1;
   maxnow = maxnum/ranp1*ranp1;
S40:   //   ign = ignlgi()-1;
   ign = rand();
   if(ign > maxnow) goto S40;
   ignuin = low+ign%ranp1;
   return ignuin;
#undef maxnum
}

double genbet(const double aa,const double bb)
//***********************************************************************
//     double genbet(double aa,double bb)
//               GeNerate BETa random deviate
//                              Function
//     Returns a single random deviate from the beta distribution with
//     parameters A and B.  The density of the beta is
//               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
//                              Arguments
//     aa --> First parameter of the beta distribution, aa > 0.0
//
//     bb --> Second parameter of the beta distribution bb > 0.0
//
//                              Method
//     R. C. H. Cheng
//     Generating Beta Variatew with Nonintegral Shape Parameters
//     Communications of the ACM, 21:317-322  (1978)
//     (Algorithms BB and BC)
// ***********************************************************************
{
#define expmax 89.0
#define infnty 1.0E38
   static double olda = -1.0;
   static double oldb = -1.0;
   static double genbet,a,alpha,b,beta,delta,gamma,k1,k2,r,s,t,u1,u2,v,w,y,z;
   static long qsame;

   qsame = olda == aa && oldb == bb;
   if (qsame) goto S20;
   if (aa <= 0.0 || bb <= 0.0) throw exception("genbet(): bad args");
   olda = aa;
   oldb = bb;
S20:
   if(!(std::min(aa,bb) > 1.0)) goto S100;
                     //  Alborithm BB Initialize
   if(qsame) goto S30;
   a = std::min(aa,bb);
   b = std::max(aa,bb);
   alpha = a+b;
   beta = sqrt((alpha-2.0)/(2.0*a*b-alpha));
   gamma = a+1.0/beta;
S30:
S40:
   u1 = ranf();   //  Step 1
   u2 = ranf();
   v = beta*log(u1/(1.0-u1));
   if(!(v > expmax)) goto S50;
   w = infnty;
   goto S60;
S50:
   w = a*exp(v);
S60:
   z = pow(u1,2.0)*u2;
   r = gamma*v-1.3862944;
   s = a+r-w;
       //  Step 2
   if(s+2.609438 >= 5.0*z) goto S70;
       //  Step 3
   t = log(z);
   if(s > t) goto S70;
         // Step 4
   if(r+alpha*log(alpha/(b+w)) < t) goto S40;
S70:    //  Step 5
   if(!(aa == a)) goto S80;
   genbet = w/(b+w);
   goto S90;
S80:
   genbet = b/(b+w);
S90:
   goto S230;
S100:   //  Algorithm BC Initialize
   if(qsame) goto S110;
   a = std::max(aa,bb);
   b = std::min(aa,bb);
   alpha = a+b;
   beta = 1.0/b;
   delta = 1.0+a-b;
   k1 = delta*(1.38889e-2+4.16667e-2*b)/(a*beta-0.777778);
   k2 = 0.25+(0.5+0.25/delta)*b;
S110:
S120:
   u1 = ranf();
       // Step 1
   u2 = ranf();
   if(u1 >= 0.5) goto S130;
        //Step 2
   y = u1*u2;
   z = u1*y;
   if(0.25*u2+z-y >= k1) goto S120;
   goto S170;
S130:
      // Step 3
   z = pow(u1,2.0)*u2;
   if(!(z <= 0.25)) goto S160;
   v = beta*log(u1/(1.0-u1));
   if(!(v > expmax)) goto S140;
   w = infnty;
   goto S150;
S140:
   w = a*exp(v);
S150:
   goto S200;
S160:
   if(z >= k2) goto S120;
S170:
     // Step 4  Step 5
   v = beta*log(u1/(1.0-u1));
   if(!(v > expmax)) goto S180;
   w = infnty;
   goto S190;
S180:
   w = a*exp(v);
S190:
   if(alpha*(log(alpha/(b+w))+v)-1.3862944 < log(z)) goto S120;
S200:
       // Step 6
   if(!(a == aa)) goto S210;
   genbet = w/(b+w);
   goto S220;
S210:
   genbet = b/(b+w);
S230:
S220:
   return genbet;
#undef expmax
#undef infnty
}
}
