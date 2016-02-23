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

#ifndef MATVEC_STAT_H
#define MATVEC_STAT_H

// random number generators
/*!
   \file stat.h
   \brief all low level statistical founctions

   \sa StatDist
*/
namespace matvec {
extern double ranf(void);
extern double snorm(void);
extern double sexpo(void);
extern double sgamma(const double a);
extern double genbet(const double aa,const double bb);
extern double genchi(const double df);
extern double genexp(const double av);
extern double genf(const int dfn,const int dfd);
extern double gengam(const double a, const double r);
extern double gennch(const int df,const double xnonc);
extern double gennf(const int dfn,const int dfd, const double xnonc);
extern double gennor(const double av,const double sd);
extern double genunf(const double low,const double high);

extern long ignbin(const long n,const double pp);
extern long ignpoi(const double mu);
extern long ignuin(const long low,const long high);

extern void set_seed(const unsigned iseed);
extern double ran1(void);
extern double gasdev(void);
extern double gamdev(const double A,const double R);

extern double factrl(const long n);
extern double BinCoef(const long n,const long k);
extern double gammln(const double xx);
extern double gammin(const double x,const double a);
extern double betain(const double x,const double a,const double b,
                     const double beta);

extern double Normal_cdf(const double x);
extern double Gamma_cdf(const double g,const double alfa,const double theta);
extern double Beta_cdf(const double x,const double alfa,const double beta,
                          const double lambda,int& errcode);
extern double F_cdf(const double f,const double df1,const double df2,
                     const double nc,int& errcode);
extern double ChiSquare_cdf(const double x,const double df,const double theta,
                            int& errcode);
extern double t_cdf(const double t,const double df,const double delta);
extern double t_expected_value(const double df, const double nc);

extern double Normal_inv(const double p);
extern double ChiSquare_inv(const double p,const double df,const double nc);
extern double t_inv(const double p,const double df,const double nc);
extern double F_inv(const double p,const double df1,const double df2,
                    const double nc);

#define chidev(df)  gamdev(df/2.0,2.0)
}
#endif

