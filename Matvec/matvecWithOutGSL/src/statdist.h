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

#ifndef MATVEC_STAT_DIST_H
#define MATVEC_STAT_DIST_H

#include <cstring>
#include <cmath>
#include <iostream>
#include "doublematrix.h"
#include "stat.h"
#include "statdistbase.h"

namespace matvec {

/**
 * uniform statistical distribution.
 *
 *  @see DiscreteUniformDist
 */
class UniformDist: public StatDistBase {
   public:
      UniformDist(const double a=0.0,const double b=1.0)
                                  {distname = "UniformDist"; reset(a,b);}
      UniformDist(const UniformDist& u):StatDistBase()
                    {distname = "UniformDist";reset(u.a_value,u.b_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << a_value
                                        << "," << b_value << ")\n"; return;}
      void    reset(const double a, const double b);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return genunf(a_value,b_value);}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return (a_value + b_value)/2.0;}
      double  variance(void) const {double c=b_value-a_value;return  c*c/12.0;}
      double  pdf(const double x) const {return 1.0/(b_value - a_value);}
      double  cdf(const double x) const{return (x-a_value)/(b_value- a_value);}
      double  mgf(const double t) const;
      double  inv(const double p) const;
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
   protected:
      double   a_value, b_value;
};

/**
 * normal statistical distribution.
 *  \class NormalDist
 *
 */

class NormalDist: public StatDistBase {
   protected:
      double mu_value;
      double sigma2_value;
      double sigma_value;
   public:
      NormalDist(const double mu=0.0, const double sigma2=1.0)
                            {distname="NormalDist"; reset(mu,sigma2);}
      NormalDist(const NormalDist& u):StatDistBase()
               {distname="NormalDist"; reset(u.mu_value,u.sigma2_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << mu_value
                                     << "," << sigma2_value << ")\n"; return;}
      void    reset(const double mu, const double sigma2);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return mu_value + sigma_value*snorm();}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return mu_value;}
      double  variance(void) const {return sigma2_value;}
      double  pdf(const double x) const;
      double  cdf(const double x) const
                              {return Normal_cdf((x-mu_value)/sigma_value);}
      double  mgf(const double t) const
                               {return std::exp(mu_value*t + sigma2_value*t*t/2.0);}
      double  inv(const double p) const
                               {return Normal_inv(p)*sigma_value + mu_value;}
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * log normal statistical distribution.
 *
 *  For more information, see the following references:
 *  <OL>
 *     <LI> Norman L. Johnson, Samuel Kotz (1970) Distributions (volume 2)
 *     <LI> Evans, et al. (1993) Statistical Distributions
 *           (2nd Ed.) Wily pp102-105
 *  </OL>
 *
 * @see NormalDist
 */
class LogNormalDist: public StatDistBase {
   protected:
      double mu_value;
      double sigma2_value;
      double sigma_value;
      double theta_value;
   public:
      LogNormalDist(const double mu=0.0, const double sigma2=1.0,
                    const double theta=1.0)
                   {distname="LogNormalDist"; reset(mu,sigma2,theta);}
      LogNormalDist(const LogNormalDist& u):StatDistBase()
                  {distname="LogNormalDist";
                  reset(u.mu_value,u.sigma2_value,u.theta_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << mu_value
                              << "," << sigma2_value << "," << theta_value
                              << ")\n"; return;}
      void    reset(const double mu,const double sigma2,const double theta);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const
                   {return std::exp(mu_value + sigma_value*snorm()) + theta_value;}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {
                         return std::exp(mu_value + 0.5*sigma2_value) + theta_value;}
      double  variance(void) const;
      double  pdf(const double x) const;
      double  cdf(const double x) const;
      double  mgf(const double t) const
                 {std::cerr << "LogNormalDist::mgf(): not available yet"; return 0.0;}
      double  inv(const double p) const
                 {return std::exp(Normal_inv(p)*sigma_value+ mu_value)+ theta_value;}
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * Chi-Square statistical distribution.
 *
 */
class ChiSquareDist: public StatDistBase {
   protected:
      double   df_value, nc_value;
   public:
      ChiSquareDist(const double df,const double nc=0.0)
                             {distname="ChiSquareDist"; reset(df,nc);}
      ChiSquareDist(const ChiSquareDist& u):StatDistBase()
              {distname="ChiSquareDist";reset(u.df_value,u.nc_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << df_value
                                   << "," << nc_value<< ")\n"; return;}
      void    reset(const double df, const double nc=0.0);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return genchi(df_value);}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return (df_value + nc_value);}
      double  variance(void) const {return 2.0*(df_value + 2.0*nc_value);}
      double  nonc(void) const {return nc_value;}
      double  pdf(const double x) const;
      double  cdf(const double x) const {int errcode=0;
                          return ChiSquare_cdf(x,df_value,nc_value,errcode);}
      double  mgf(const double t) const;
      double  inv(const double p) const {
                              return ChiSquare_inv(p,df_value,nc_value);}
      double  nonct(const double cv,const double p) const;
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * t statistical distribution.
 *
 */
class tDist: public StatDistBase {
   protected:
      double  df_value, nc_value;
   public:
      tDist(const double df, const double nc=0.0) {distname="tDist"; reset(df,nc);}
      tDist(const tDist& u):StatDistBase()
	{distname="tDist"; reset(u.df_value,u.nc_value);}

      void    display(void) const {std::cout << "\t" << distname << "("  << df_value
                                        << "," << nc_value<< ")\n"; return;}
      void    reset(const double df, const double nc=0.0);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const
              {return std::sqrt(df_value)*snorm()/std::sqrt(genchi(int(df_value)));}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return t_expected_value(df_value,nc_value);}
      double  variance(void) const;
      double  nonc(void) const {return nc_value;}
      double  pdf(const double x) const;
      double  cdf(const double x) const
                               {return t_cdf(x,df_value,nc_value);}
      double  mgf(const double t) const;
      double  inv(const double p) const {return t_inv(p,df_value,nc_value);}
      double  nonct(const double cv,const double p) const;
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};


/**
 * F statistical distribution.
 *
 */
class FDist: public StatDistBase {
   protected:
      double  df1_value,df2_value, nc_value;
   public:
      FDist(const double df1,const double df2, const double nc=0.0)
                                {distname="FDist"; reset(df1,df2,nc);}
      FDist(const FDist& u):StatDistBase()
          {distname="FDist";reset(u.df1_value,u.df2_value,u.nc_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << df1_value
                      << "," << df2_value<< ","<< nc_value << ")\n"; return;}
      void    reset(const double df1,const double df2, const double nc=0.0);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return genf(int(df1_value),int(df2_value));}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const;
      double  variance(void) const;
      double  nonc(void) const {return nc_value;}
      double  pdf(const double x) const;
      double  cdf(const double x) const {int errcode=0;
                       return F_cdf(x,df1_value,df2_value,nc_value,errcode);}
      double  mgf(const double t) const;
      double  inv(const double p) const
                          {return F_inv(p,df1_value,df2_value,nc_value);}
      double  nonct(const double cv,const double p) const;
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * gamma statistical distribution.
 *
 */
class GammaDist: public StatDistBase {
   protected:
      double    alpha_value, theta_value;
   public:
    GammaDist(const double a,const double t)
                                   {distname="GammaDist";reset(a,t);}
      GammaDist(const GammaDist& u):StatDistBase()
             {distname="GammaDist";reset(u.alpha_value,u.theta_value);}

      void    display(void) const {std::cout << "\t" << distname << "("
                      << alpha_value << "," << theta_value << ")\n"; return;}
      void    reset(const double a, const double t);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return theta_value*sgamma(alpha_value);}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return alpha_value*theta_value;}
      double  variance(void) const {return alpha_value*theta_value*theta_value;}
      double  pdf(const double x) const;
      double  cdf(const double x) const;
      double  mgf(const double t) const;
      double  inv(const double p) const;
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * exponential statistical distribution.
 *
 */
class ExponentialDist: public GammaDist {

   public:
      ExponentialDist(const double t=1.0): GammaDist(1.0,t) {distname="ExponentialDist";}
      ExponentialDist(const ExponentialDist& u): GammaDist(1.0,u.theta_value) {distname="ExponentialDist";}

      void    display(void) const {std::cout << "\t" << distname << "("
                                        << theta_value << ")\n"; return;}
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * beta statistical distribution.
 *
 */
class BetaDist: public StatDistBase {
   protected:
      double    alpha_value, beta_value, nc_value;
   public:
      BetaDist(const double a,const double b,const double nc=0.0)
                                  {distname="BetaDist"; reset(a,b,nc);}
      BetaDist(const BetaDist& u):StatDistBase() 
	{distname="BetaDist"; reset(u.alpha_value,u.beta_value,u.nc_value);}

      void    display(void) const {std::cout << "\t" << distname << "("
                                << alpha_value << "," << beta_value << ","
                                << nc_value << ")\n"; return;}
      void    reset(const double a, const double b, const double nc=0.0);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return genbet(alpha_value,beta_value);}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const;
      double  variance(void) const;
      double  pdf(const double x) const;
      double  cdf(const double x) const {int errcode=0;
                   return Beta_cdf(x,alpha_value,beta_value,nc_value,errcode);}
      double  mgf(const double t) const;
      double  inv(const double p) const;
      double  nonct(const double cv,const double p) const
                      {std::cerr << "BetaDist::nonct(): not available\n"; return p;}
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * discreteuniform statistical distribution.
 *
 */
class DiscreteUniformDist: public StatDistBase {
   protected:
      long   a_value, b_value;
   public:
      DiscreteUniformDist(const long a,const long b)
                        {distname="DiscreteUniformDist"; reset(a,b);}
      DiscreteUniformDist(const DiscreteUniformDist& u):StatDistBase()
                  {distname="DiscreteUniformDist";
                   reset(u.a_value,u.b_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << a_value
                                        << "," << b_value << ")\n"; return;}
      void    reset(const long a, const long b);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return double(ignuin(a_value,b_value));}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return 0.5*(a_value + b_value);}
      double  variance(void) const {double c=b_value-a_value;return  c*c/12.0;}
      double  pdf(const double x) const {return pdf(long(x));}
      double  cdf(const double x) const {return cdf(long(x));}
      double  pdf(const long x) const {return 1.0/(b_value - a_value);}
      double  cdf(const long x) const
                          {return double(x-a_value)/double(b_value - a_value);}
      double  mgf(const double t) const;
      double  inv(const double p) const
                     {std::cerr << " DiscreteUniformDist::inv(): not available\n"; return p;}
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * binomial statistical distribution.
 *
 */
class BinomialDist: public StatDistBase {
   protected:
      long      n_value;
      double    p_value;
   public:
      BinomialDist(const long n,const double p)
                                  {distname="BinomialDist"; reset(n,p);}
      BinomialDist(const BinomialDist& u):StatDistBase()
                {distname="BinomialDist";reset(u.n_value,u.p_value);}

      void    display(void) const {std::cout << "\t" << distname << "(" << n_value
                                        << "," << p_value << ")\n"; return;}
      void    reset(const long n, const double p);
      void    sample(Vector<double>& x) const;
      void    sample(doubleMatrix& x) const;
      double  sample(void) const {return double(ignbin(n_value,p_value));}
      Vector<double>  sample(unsigned n) const;
      doubleMatrix  sample(unsigned m,unsigned n) const;
      double  mean(void) const {return p_value*n_value;}
      double  variance(void) const {return p_value*(1.0-p_value)*n_value;}
      double  pdf(const double x) const {return pdf(long(x));}
      double  cdf(const double x)  const {return cdf(long(x));}
      double  pdf(const long k) const;
      double  cdf(const long k ) const;
      double  mgf(const double t) const
                            {return std::pow(1.0-p_value + p_value*std::exp(t),(double)n_value);}
      double  inv(const double p) const
                       {std::cerr << " BinomialDist::inv(): not available\n"; return p;}
      double  parameter(const int k) const;
      void    parameter(const int k,const double x);
};

/**
 * Poisson statistical distribution.
 *
 */
class PoissonDist: public StatDistBase {
   protected:
      double    lambda_value;
   public:
      PoissonDist(const double l) {distname="PoissonDist"; reset(l);}
      PoissonDist(const PoissonDist& u):StatDistBase()
                    {distname="PoissonDist";reset(u.lambda_value);}

      void     display(void) const {std::cout << "\t" << distname << "("
                                         << lambda_value << ")\n"; return;}
      void     reset(const double l);
      void     sample(Vector<double>& x) const;
      void     sample(doubleMatrix& x) const;
      double   sample(void) const {return double(ignpoi(lambda_value));}
      Vector<double>   sample(unsigned n) const;
      doubleMatrix   sample(unsigned m,unsigned n) const;
      double   mean(void) const {return lambda_value;}
      double   variance(void) const {return lambda_value;}
      double   pdf(const double x) const {return pdf(long(x));}
      double   cdf(const double x) const {return cdf(long(x));}
      double   pdf(const long k) const;
      double   cdf(const long k) const;
      double   mgf(const double t) const
                                {return std::exp(lambda_value*(std::exp(t)-1.0));}
      double   inv(const double p) const
                      {std::cerr << " PoissonDist::inv(): not available\n"; return p;}
      double   parameter(const int k) const;
      void     parameter(const int k,const double x);
};

/**
 * Geometric statistical distribution.
 *
 */
class GeometricDist: public StatDistBase {
   protected:
      double    p_value;
   public:
      GeometricDist(const double p) {distname="GeometricDist";reset(p);}
      GeometricDist(const GeometricDist& u):StatDistBase()
                          {distname="GeometricDist"; reset(u.p_value);}

      void     display(void) const {std::cout << "\t" << distname << "(" << p_value
                                         << ")\n"; return;}
      void     reset(const double p);
      void     sample(Vector<double>& x) const;
      void     sample(doubleMatrix& x) const;
      double   sample(void) const {return std::floor(std::log(ranf())/std::log(1.0-p_value));}
      Vector<double>   sample(unsigned n) const;
      doubleMatrix   sample(unsigned m,unsigned n) const;
      double   mean(void) const {return (1.0-p_value)/p_value;}
      double   variance(void) const {return (1.0-p_value)/(p_value*p_value);}
      double   pdf(const double x) const {return pdf(long(x));}
      double   cdf(const double x) const {return cdf(long(x));}
      double   pdf(const long k) const;
      double   cdf(const long k) const;
      double   mgf(const double t) const;
      double   inv(const double p) const
                      {std::cerr << " GeometricDist::inv(): not available\n"; return p;}
      double   parameter(const int k) const;
      void     parameter(const int k,const double x);
};

/**
 * NegativeBinomial statistical distribution.
 *  For more information see the reference:
 *   Robert V. Hogg and Elliot A. Tanis (1983) Probability and Statistical
 *    Inference (2nd Ed) pp 80-84
 *
 * @see BinomialDist
 */
class NegativeBinomialDist: public StatDistBase {
   protected:
      long    n_value;
      double  p_value;
   public:
      NegativeBinomialDist(const long n,const double p)
                 {distname="NegativeBinomialDist";reset(n,p);}
      NegativeBinomialDist(const NegativeBinomialDist& u):StatDistBase()
                      {distname="NegativeBinomialDist";
                       reset(u.n_value,u.p_value);}

      void     display(void) const {std::cout << "\t" << distname << "(" << n_value
                                         << "," << p_value << ")\n"; return;}
      void     reset(const long n,const double p);
      void     sample(Vector<double>& x) const;
      void     sample(doubleMatrix& x) const;
      double   sample(void) const;
      Vector<double>   sample(unsigned n) const;
      doubleMatrix   sample(unsigned m,unsigned n) const;
      double   mean(void) const {return (1.0-p_value)*double(n_value)/p_value;}
      double   variance(void) const
                    {return (1.0-p_value)*double(n_value)/(p_value*p_value);}
      double   pdf(const double x) const {return pdf(long(x));}
      double   cdf(const double x) const {return cdf(long(x));}
      double   pdf(const long k) const;
      double   cdf(const long k) const;
      double   mgf(const double t) const;
      double   inv(const double p) const
               {std::cerr << " NegativeBinomialDist::inv(): not available\n"; return p;}
      double   parameter(const int k) const;
      void     parameter(const int k,const double x);
};
}
#endif
