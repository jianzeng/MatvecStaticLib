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

#include "statdist.h"

namespace matvec {

/////////////////// UniformDist class /////////////////////////////

void UniformDist::reset(const double a, const double b)
{
   if (a >= b) throw exception("UniformDist::reset(a,b): a >= b");
      a_value = a;
      b_value = b;
}

void UniformDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   if (a_value==0.0 && b_value==1.0) {
      while (bot < top ) *bot++ = ranf();
   }
   else {
      double c = b_value-a_value;
      while (bot < top ) *bot++ = a_value + c*ranf();
   }
}
void UniformDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   unsigned i;
   if (a_value==0.0 && b_value==1.0) {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = ranf();
      }
   }
   else {
      double c = b_value-a_value;
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = a_value + c*ranf();
      }
   }
}

Vector<double> UniformDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix UniformDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double UniformDist::inv(const double p) const
{
   double retval = 0.0;
   if (p<0.0 || p>1.0) throw exception("UniformDist::inv(): bad arg, value in [0,1] expected");
   retval = a_value + p*(b_value- a_value);
   return retval;
}

double UniformDist::mgf(const double t) const
{
   double mgfval = 0.0;
   if (t > 0.0) {
      mgfval = (std::exp(t*b_value) - std::exp(t*a_value))/(t*(b_value - a_value));
   }
   return mgfval;
}

double UniformDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = a_value;
   }
   else if (k==2) {
      par = b_value;
   }
   else {
      throw exception("UniformDist::parameter(): bad arg value, 1 or 2 is expected");
   }
   return par;
}

void UniformDist::parameter(const int k,const double x)
{
   double par = 0.0;
   if (k==1) {
      a_value = x;
   }
   else if (k==2) {
      b_value = x;
   }
   else {
      throw exception("UniformDist::parameter(): bad arg 1");
   }
   return;
}

///////////////////  NormalDist  class /////////////////////////////

void NormalDist::reset(const double mu, const double sigma2)
{
   if (sigma2 <= 0.0) throw exception("NormalDist::reset(): 2nd arg must be positive");
   mu_value = mu;  sigma2_value = sigma2;
   sigma_value = std::sqrt(sigma2_value);
}

/**
 * The probability density function (pdf) is defined by
 * \f[
 *  f(x) = \frac{1}{\sqrt{2\pi}\sigma} \exp (-\frac{(x-\mu)^2}{2\sigma^2}),
 *   \quad -\infty < x < \infty
 * \f]
 */
double NormalDist::pdf(const double x) const
{
   const double c = 1.0/std::sqrt(4.0*std::asin(1.0));        //    1.0/sqrt(2.0*pi)   pi = 2.0*std::asin(1.0)
   double xx = x - mu_value;
   return  c/sigma_value * std::exp(-0.5*xx*xx/sigma2_value);
}

void NormalDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   if (mu_value==0.0 && sigma_value==1.0) {
      while (bot < top ) *bot++ = snorm();
   }
   else {
      while (bot < top ) *bot++ = mu_value + sigma_value*snorm();
   }
}
void NormalDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   unsigned i;
   if (mu_value==0.0 && sigma2_value==1.0) {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = snorm();
      }
   }
   else {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = mu_value + sigma_value*snorm();
      }
   }
}

Vector<double> NormalDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix NormalDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double NormalDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = mu_value;
   }
   else if (k==2) {
      par = sigma2_value;
   }
   else {
      throw exception("NormalDist::parameter(): bad arg, 1 or 2 is expected");
   }
   return par;
}

void NormalDist::parameter(const int k,const double x)
{
   if (k==1) {
      mu_value = x;
   }
   else if (k==2) {
      if (x > 0.0) {
         sigma2_value = x;
      }
      else {
         throw exception("NormalDist::parameter(): 2nd arg must be positive");
      }
      sigma_value = std::sqrt(sigma2_value);
   }
   else {
      throw exception("NormalDist::parameter(): bad arg");
   }
}

/////////////////// LogNormalDist  class /////////////////////////////
void LogNormalDist::reset(const double mu, const double sigma2,
                          const double theta)
{
   if (sigma2 <= 0.0)  throw exception("LogNormalDist::reset(): 2nd arg must be positive");
   mu_value = mu;  sigma2_value = sigma2; theta_value = theta;
   sigma_value = std::sqrt(sigma2_value);
}

/*!
The probability density function (pdf) is defined by
   \f[
        f(x) = [(x-\theta)\sqrt(2\pi)\sigma]^{-1}
            \exp(-\frac{\log(x-\theta)-\mu)^2}{2\sigma^2}),
                \quad x > \theta
   \f]
*/
double LogNormalDist::pdf(const double x) const
{
   if (x <= theta_value) throw exception("LogNormalDist::pdf(arg): bad arg");
   double tmp = (std::log(x-theta_value)-mu_value)/sigma_value;
   double retval = std::exp(-0.5*tmp*tmp);
   retval /= (x - theta_value)*std::sqrt(4.0*std::asin(1.0))*sigma_value;   // pi = 2.0*std::asin(1.0)
   return retval;
}

double LogNormalDist::cdf(const double x) const
{
   if (x <= theta_value) throw exception("LogNormalDist::pdf(arg): bad arg");
   double retval = (std::log(x - theta_value) - mu_value)/sigma_value;
   retval = Normal_cdf(retval);
   return retval;
}

double LogNormalDist::variance(void) const
{
   double mean_value = std::exp(mu_value + 0.5*sigma2_value) + theta_value;
   return std::exp(2*(mu_value + sigma2_value)) + 2*theta_value*mean_value
             + theta_value*theta_value - mean_value*mean_value;
}

void LogNormalDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   if (mu_value==0.0 && sigma_value==1.0) {
      while (bot < top ) *bot++ = std::exp(snorm()) + theta_value;
   }
   else {
      while (bot<top) *bot++ = std::exp(mu_value+ sigma_value*snorm())+ theta_value;
   }
}
void LogNormalDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   unsigned i;
   if (mu_value==0.0 && sigma2_value==1.0) {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = std::exp(snorm()) + theta_value;
      }
   }
   else {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) {
            *bot++ = std::exp(mu_value + sigma_value*snorm()) + theta_value;
         }
      }
   }
}

Vector<double> LogNormalDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix LogNormalDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double LogNormalDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = mu_value;
   }
   else if (k==2) {
      par = sigma2_value;
   }
   else if (k==3) {
      par = theta_value;
   }
   else {
      throw exception("LogNormalDist::parameter(): bad arg");
   }
   return par;
}

void LogNormalDist::parameter(const int k,const double x)
{
   if (k==1) {
      mu_value = x;
   }
   else if (k==2) {
      sigma2_value = x;
   }
   else if (k==3) {
      theta_value = x;
   }
   else {
      throw exception("LogNormalDist::parameter(): bad arg 1");
   }
}

///////////////////  ChiSquare class /////////////////////////////

void ChiSquareDist::reset(const double df, const double nc)
{
   if (df <= 0.0 || nc < 0.0) throw exception("ChiSquareDist::reset(): bad arg");
   df_value = df;
   nc_value = nc;
}

/*!
The random variable <EM>X</EM> has a non-central \f$\chi^2\f$ distribution if its probability density function (pdf) is defined by
\f[
  f(x) = \frac{\exp(-(x+\lambda)/2)}{2\frac{1}{2}r}
 \sum_{j=0}^\infty \frac{x^{r/2+j-1}\lambda^j}{\Gamma(r/2+j)2^{2j}j!},
 \quad 0 \leq x < \infty
\f]
*/
double ChiSquareDist::pdf(const double x) const
{
   if (nc_value != 0.0) throw exception("ChiSquareDist::pdf(): not available yet: noncentrality");
   if ( x < 0.0) throw exception("ChiSquareDist::pdf(): bad arg, nonnegative value expected");
   double pr;
   double r2 = 0.5*df_value;
   pr = (r2-1.0)*std::log(x) - x/2.0 - gammln(r2) - r2*std::log(2.0);
   return std::exp(pr);
}

double ChiSquareDist::mgf(const double t) const
{
   if (nc_value != 0.0) throw exception("ChiSquareDist::mgf(): not available yet: noncentrality");
   if (t >= 0.5) throw exception("ChiSquareDist::mgf(t):  t must be < 1/2");
   return 1.0/pow(1.0-2.0*t,df_value/2.0);
}

double ChiSquareDist::nonct(const double cv,const double p) const
{
   if (p<0.0 || p>1.0) throw exception("ChiSquareDist::nonct(): bad arg, value in (0,1)  expected");
   int errcode = 0;
   double eps = 1.0e-6;
   double bot = 0.0;
   double top = df_value + 26000.0;
   double p2,ppt;
   while (top-bot > eps) {
      ppt = (top + bot)*0.5;
      p2 = ChiSquare_cdf(cv,df_value,ppt,errcode);
      if (errcode != 0) break;
      if (p2 < p ) {
         top = ppt;
      }
      else if (p2 > p) {
         bot = ppt;
      }
      else {
         break;
      }
   }
   return ppt;
}

void ChiSquareDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   //int df = static_cast<int>(df_value);
   if (nc_value == 0.0) {
      while (bot < top ) *bot++ = genchi(df_value);
   }
   else {
	  int df = static_cast<int>(df_value);
      while (bot < top ) *bot++ = gennch(df,nc_value);
   }
}

void ChiSquareDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   unsigned i;
   //int df = static_cast<int>(df_value);
   if (nc_value == 0.0) {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = genchi(df_value);
      }
   }
   else {
      int df = static_cast<int>(df_value);
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = gennch(df,nc_value);
      }
   }
}

Vector<double> ChiSquareDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix ChiSquareDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double ChiSquareDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = df_value;
   }
   else if (k==2) {
      par = nc_value;
   }
   else {
      throw exception("ChiSquareDist::parameter(): bad arg value, 1 or 2 expected");
   }
   return par;
}

void ChiSquareDist::parameter(const int k,const double x)
{
   if (k==1) {
      df_value = x;
   }
   else if (k==2) {
      if (x >= 0.0) {
         nc_value = x;
      }
      else {
         throw exception("ChiSquareDist::parameter(): bad arg2");
      }
   }
   else {
      throw exception("ChiSquareDist::parameter(): bad arg1");
   }
}


///////////////////  tDist  class /////////////////////////////

void tDist::reset(const double df, const double nc)
{
   if (df <= 0.0) throw exception("tDist::reset(): bad arg1");
   df_value = df;
   nc_value = nc;
}

double tDist::variance(void) const
{
   if (df_value < 3)  throw exception("tDist::variance(): df < 3");
   double u = this->mean();
   u = - u*u + (1.0 + nc_value*nc_value)*df_value/(df_value - 2.0);
   return u;
}

double tDist::pdf(const double x) const
{
   if (nc_value == 0.0) throw exception("tDist::pdf(): noncentrality: not available");
   double pr = 0.0;
   double r = df_value;
   double rr = (df_value + 1.0)/2.0;
   double pi = 2.0*std::asin(1.0);               // pi =  2.0*std::asin(1.0)
   pr = gammln(rr)- 0.5*std::log(pi*r) - gammln(r/2.0) - rr*std::log(1.0+x*x/r);
   return std::exp(pr);
}

double tDist::mgf(const double t) const
{
   throw exception("tDist::mgf(): not available yet");
   return 0.0;
}

double tDist::nonct(const double cv,const double p) const
{
   if (p<0.0 || p>1.0) throw exception("tDist::nonct(): bad arg, value between (0,1) expected");
   double eps = 1.0e-6;
   double p2,ppt;
   double bot = - 700.0;
   double top =  700.0;
   while (top-bot > eps) {
      ppt = (top + bot)*0.5;
      p2 = t_cdf(cv,df_value,ppt);
      if (p2 < p ) {
         top = ppt;
      }
      else if (p2 > p) {
         bot = ppt;
      }
      else {
         break;
      }
   }
   return ppt;
}

void tDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   double c = std::sqrt(df_value);
   while (bot< top ) *bot++ = c*(snorm()+nc_value)/std::sqrt(genchi(static_cast<int>(df_value)));
}

void tDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   double c = std::sqrt(df_value);
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) {
         *bot++ = c*(snorm()+nc_value)/std::sqrt(genchi(static_cast<int>(df_value)));
      }
   }
}
Vector<double> tDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix tDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double tDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = df_value;
   }
   else if (k==2) {
      par = nc_value;
   }
   else {
      throw exception("tDist::parameter(): bad arg value, 1 or 2 is expected");
   }
   return par;
}

void tDist::parameter(const int k,const double x)
{
   if (k==1) {
      if (x > 0.0) {
         df_value = x;
      }
      else {
         throw exception("tDist::parameter(): bad arg2");
      }
   }
   else if (k==2) {
      nc_value = x;
   }
   else {
      throw exception("tDist::parameter(): bad arg1");
   }
}

///////////////////  FDist  class /////////////////////////////

void FDist::reset(const double df1, const double df2,const double nc)
{
   if (nc < 0.0) throw exception("FDist::reset(): bad arg3: nonnegative value expected");
   if (df1 <= 0.0 || df2 <= 0.0) throw exception("FDist::reset(): bad arg1 or 2, positives expected");

   df1_value = df1;
   df2_value = df2;
   nc_value = nc;
}

/*!
If U<SUP>1</SUP> is distributed as \f$ \chi^2(r_1,\delta)\f$ and U<SUP>2</SUP> is ditributed as \f$\chi^2(r_2)\f$ are independent then the random variable
\f[
   X = \frac{U_1/r_1}{U_2/r_2}
\f]
is called the non-central F distribution with r1 (integer) and r2 (integer) degrees of freedom and non-centrality parameter \f$\delta\f$ (real).
*/
double FDist::pdf(const double x) const
{
   if (x < 0.0) throw exception("FDist::pdf(): bad arg, nonnegative value expected");
   if (nc_value != 0.0) throw exception("FDist::pdf():  not available yet: noncentrality");
   double r1 = df1_value;
   double r2 = df2_value;
   double r12 = 0.5*(r1+r2);
   double r3 = 0.5*r1;
   double t1 = gammln(r12) + r3*(std::log(r1)-std::log(r2)) + (r3-1.0)*std::log(x);
   double t2 = gammln(r3) + gammln(0.5*r2) + r12*std::log(1.0+r1*x/r2);
   double pr = t1 - t2;
   return std::exp(pr);
}

double FDist::mgf(const double t) const
{
   throw exception("FDist:mgf(): not available yet");
   return 0.0;
}

double FDist::mean(void) const
{
   if (df2_value < 3 ) throw exception(" FDist::mean(): not exist:  df2 < 3");
   return df2_value*(df1_value + nc_value)/(df1_value*(df2_value - 2.0));
}

double FDist::variance(void) const
{
   if (df2_value < 5) throw exception("FDist::variance(): not exist: df2 < 5");
   double x,v1,v2,t1,t2;
   v1 = df1_value;
   v2 = df2_value;
   x = v2 - 2.0;
   t1 = v2*v2*((v1+nc_value)*(v1+nc_value) + (v1+nc_value*2.0)*x);
   t2 = v1*v1*x*x*(v2 - 4.0);
   return t1/t2*2.0;
}

double FDist::nonct(const double cv,const double p) const
{
   if (p<0.0 || p>1.0) throw exception("FDist::nonct(): bad arg, value between (0,1) expected");
   double p2,ppf;
   if (df2_value >= 3) {
      ppf = df2_value/(df2_value -2.0);
   }
   else {
      ppf = 0.0;
   }
   double eps = 1.0e-6;
   double bot = 0.0;
   double top = ppf + 6000.0;
   int errcode = 0;
   while (top-bot > eps) {
      ppf = (top + bot)*0.5;
      p2 = F_cdf(cv,df1_value,df2_value,ppf,errcode);
      if (errcode != 0) break;
      if (p2 < p ) {
         top = ppf;
      }
      else if (p2 > p) {
         bot = ppf;
      }
      else {
         break;
      }
   }
   return ppf;
}

void FDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   int df1 = static_cast<int>(df1_value);
   int df2 = static_cast<int>(df2_value);
   if (nc_value == 0.0) {
      while (bot < top ) *bot++ = genf(df1,df2);
   }
   else {
      while (bot < top ) *bot++ = gennf(df1,df2,nc_value);
   }
}

void FDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   unsigned i;
   int df1 = static_cast<int>(df1_value);
   int df2 = static_cast<int>(df2_value);
   if (nc_value == 0.0) {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = genf(df1,df2);
      }
   }
   else {
      for (i=0; i<nr; i++) {
         bot = x[i];  top = &bot[nc];
         while (bot < top ) *bot++ = gennf(df1,df2,nc_value);
      }
   }
}

Vector<double> FDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix FDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double FDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = df1_value;
   }
   else if (k==2) {
      par = df2_value;
   }
   else if (k==3) {
      par = nc_value;
   }
   else {
      throw exception("FDist::parameter(): bad arg value, 1(2,3) is expected");
   }
   return par;
}

void FDist::parameter(const int k,const double x)
{
   if (k==1) {
      if (x <= 0.0) {
         throw exception("FDist::parameter(): bad arg2, positives expected");
      }
      else {
         df1_value = x;
      }
   }
   else if (k==2) {
      if (x <= 0.0) {
         throw exception("FDist::parameter(): bad arg2, positives expected");
      }
      else {
         df2_value = x;
      }
   }
   else if (k==3) {
      if (x < 0.0) {
         throw exception("FDist::parameter(): bad arg2, nonnegative expected");
      }
      else {
         nc_value = x;
      }
   }
   else {
      throw exception("FDist::parameter(): bad arg1");
   }
}

///////////////////  GammaDist class /////////////////////////////

void GammaDist::reset(const double a, const double t)
{
   if (a <= 0.0 || t <= 0.0) throw exception("GammaDist::reset(): bad args, they mustbe > 0");
   alpha_value = a;
   theta_value = t;
}

/*!
The random variable <EM>X</EM> has a gamma distribution if its probability density function is defined by
\f[
f(x) = \frac{1}{\Gamma(\alpha)\theta^\alpha} x^{\alpha-1} e^{-x/\theta}, \quad 0 \leq x < \infty.
\f]
*/
double GammaDist::pdf(const double x) const
{
   if (x < 0.0) throw exception("GammaDist::pdf(): bad arg");
   double pr = 0.0;
   double a = gammln(alpha_value) + alpha_value*std::log(theta_value);
   if (x == 0.0) {
      pr = -a;
   }
   else {
      pr = (alpha_value-1.0)*std::log(x) - x/theta_value - a;
   }
   return std::exp(pr);
}

double GammaDist::cdf(const double x) const
{
   if (x < 0.0) throw exception("GammaDist::cdf(): bad arg");
   return Gamma_cdf(x,alpha_value,theta_value);
}

double GammaDist::mgf(const double t) const
{
   if (t >= 1.0/theta_value) throw exception("GammaDist:mgf(t): t must be < 1/theta");
   return 1.0/pow(1.0-theta_value*t,alpha_value);
}

double GammaDist::inv(const double p) const
{
   if (p<0.0 || p>1.0) throw exception("GammaDist::inv(): bad arg, value between [0,1] expected");
   double eps = 1.0e-6;
   double bot = 0.0;
   double top = alpha_value*theta_value + 1200.0;
   double p2,ppt;
   while (top-bot > eps) {
      ppt = (top + bot)*0.5;
      p2 = Gamma_cdf(ppt,alpha_value,theta_value);
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

void GammaDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = theta_value*sgamma(alpha_value);
}

void GammaDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = theta_value*sgamma(alpha_value);
   }
}

Vector<double> GammaDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix GammaDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double GammaDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = alpha_value;
   }
   else if (k==2) {
      par = theta_value;
   }
   else {
      throw exception("GammaDist::parameter(): bad arg value, 1 or 2 is expected");
   }
   return par;
}

void GammaDist::parameter(const int k,const double x)
{
   if (x <= 0.0) throw exception("GammaDist::parameter(): bad arg2, they must > 0");
   if (k == 1) {
      alpha_value = x;
   }
   else if (k == 2) {
      theta_value = x;
   }
   else {
      throw exception("GammaDist::parameter(%d): bad arg1");
   }
}

/////////////////// ExponentialDist class /////////////////////////////

double ExponentialDist::parameter(const int k) const
{
   if (k != 1) throw exception("ExponentialDist::parameter(): bad arg, 1 is expected");
   return theta_value;
}

void ExponentialDist::parameter(const int k,const double x)
{
   if (k != 1) throw exception("ExponentialDist::parameter(): bad arg1");
   theta_value = x;
}


///////////////////  BetaDist class /////////////////////////////

void BetaDist::reset(const double a, const double b,const double nc)
{
   if (a <= 0.0 || b <= 0.0) throw exception("BetaDist::reset(): invalid args");
   alpha_value = a;
   beta_value = b;
   nc_value = nc;
}

double BetaDist::pdf(const double x) const
{
   if (x <=0.0 || x>=1.0) throw exception("BetaDist::pdf(): bad arg, value in (0,1) is expected");
   double pr = 0.0;;
   if (nc_value == 0.0) {
      double a = alpha_value;
      double b = beta_value;
      pr = gammln(a+b) - gammln(a) - gammln(b);
      pr += (a-1.0)*std::log(x) + (b-1.0)*std::log(1.0-x);
   }
   else {
      throw exception("BetaDist::pdf(): not available yet : noncentrality");
   }
   return std::exp(pr);
}

double BetaDist::mgf(const double t) const
{
   throw exception("BetaDist:mgf(t): not available yet");
   return 0.0;
}

double BetaDist::mean(void) const
{
   if (nc_value != 0.0) throw exception("BetaDist::mean(): not availabe: noncentrality");
   return alpha_value/(alpha_value+beta_value);
}

double BetaDist::variance(void) const
{
   if (nc_value != 0.0) throw exception("BetaDist::variance(): not availabe: noncentrality");
   double a = alpha_value*beta_value;
   double b = alpha_value+beta_value;
   return a/((b+1)*b*b);
}

void BetaDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = genbet(alpha_value,beta_value);
}

void BetaDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = genbet(alpha_value,beta_value) ;
   }
}

Vector<double> BetaDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix BetaDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double BetaDist::inv(const double p) const
{
   if (p<0.0 || p>1.0) throw exception("BetaDist::inv(): bad arg, value between [0,1] expected");
   int errcode = 0;
   double eps = 1.0e-6;
   double bot = 0.0;
   double top = 1.0;
   double p2,ppt;
   while (top-bot > eps) {
      ppt = (top + bot)*0.5;
      p2 = Beta_cdf(ppt,alpha_value,beta_value,nc_value,errcode);
      if (errcode != 0) break;
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

double BetaDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = alpha_value;
   }
   else if (k==2) {
      par = beta_value;
   }
   else if (k==3) {
      par = nc_value;
   }
   else {
      throw exception("BetaDist::parameter(): bad arg, value 1(2,3) is expected");
   }
   return par;
}

void BetaDist::parameter(const int k,const double x)
{
   if (x <= 0.0) throw exception("BetaDist::parameter(): bad arg2, it must > 0");

   if (k==1) {
      alpha_value = x;
   }
   else if (k==2) {
      beta_value = x;
   }
   else if (k==3) {
      nc_value = x;
   }
   else {
      throw exception("BetaDist::parameter(): bad arg1");
   }
}

///////////////////   DiscreteUniformDist class /////////////////////////////

void DiscreteUniformDist::reset(const long a, const long b)
{
   if (a > b) throw exception("DiscreUniformDist::reset(a,b): b must be larger than a");
   a_value = a;
   b_value = b;
}

double DiscreteUniformDist::mgf(const double t) const
{
   throw exception("DiscreteUniformDist:mgf(): not available yet");
   return 0.0;
}

void DiscreteUniformDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = static_cast<double>(ignuin(a_value,b_value));
}

void DiscreteUniformDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = static_cast<double>(ignuin(a_value,b_value));
   }
}

Vector<double> DiscreteUniformDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix DiscreteUniformDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double DiscreteUniformDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = static_cast<double>(a_value);
   }
   else if (k==2) {
      par = static_cast<double>(b_value);
   }
   else {
      throw exception("DiscreteUniformDist::parameter(): bad arg, 1 or 2 is expected");
   }
   return par;
}

void DiscreteUniformDist::parameter(const int k, const double x)
{
   if (k==1) {
      if (x <= b_value) {
         a_value = static_cast<long int>(x);
      }
      else {
         throw exception("DiscreUniformDist::parameter(): bad arg2");
      }
   }
   if (k==2) {
      if (a_value <= x) {
         b_value = static_cast<long int>(x);
      }
      else {
         throw exception("DiscreUniformDist::parameter(): bad arg2");
      }
   }
   else {
      throw exception("DiscreteUniformDist::parameter(): bad arg1");
   }
}

///////////////////  BinomialDist class /////////////////////////////

void BinomialDist::reset(const long n, const double p)
{
   if (p < 0.0 || p > 1.0) throw exception("BinomialDist::reset(): bad arg2");
   n_value = n;
   p_value = p;
}

/*!
   The random variable <EM>X</EM> has a binomial distribution if its probability density function (pdf) is defined by
\f[
  f(x) = \frac{n!}{(n-x)!x!}p^x(1-p)^{n-x}, \quad x = 0,1,2,\ldots,n
\f]
*/
double BinomialDist::pdf(const long k) const
{
   if (k > n_value) throw exception("BinomialDist::pdf(): bad arg");
   double pr = std::log(p_value)*k + std::log(1.0-p_value)*(n_value - k);
   pr += gammln(static_cast<double>(n_value + 1))
        - gammln(static_cast<double>(k + 1)) - gammln(static_cast<double>(n_value - k + 1));
   if (pr > 0.0) return 1.0;
   if (pr < -600.00) return 0.0;
   return std::exp(pr);
}

double BinomialDist::cdf(const long k) const
{
   // cdf(k) = Pr(X <= k) = sum{Bin(n,p)), over x=0,1,...k}
   double retval;
   if (k >= n_value) return 1.0;
   if (k > 0L) {
      double a = static_cast<double>(n_value - k);
      double b = static_cast<double>(k+1);
      double lb = gammln(a) + gammln(b) - gammln(a+b);
      retval = betain(1.0 - p_value, a, b, lb);
   }
   else if (k == 0L) {
     retval = pow(1.0 - p_value, static_cast<double>(n_value));
   }
   else {
      retval = 0.0;
   }
   return retval;
}

void BinomialDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = static_cast<double>(ignbin(n_value,p_value));
}

void BinomialDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = static_cast<double>(ignbin(n_value,p_value));
   }
}

Vector<double> BinomialDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix BinomialDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double BinomialDist::parameter(const int k) const
{
   double par =  0.0;
   if (k==1) {
      par = static_cast<double>(n_value);
   }
   else if (k==2) {
      par = p_value;
   }
   else {
      throw exception("BinomialDist::parameter(): bad arg, 1 or 2 is expected");
   }
   return par;
}

void BinomialDist::parameter(const int k, const double x)
{
   if (k==1) {
      n_value = static_cast<long int>(x);
   }
   else if (k==2) {
      if (x >= 0.0 && x <= 1.0) {
         p_value = x;
      }
      else {
         throw exception("BinomialDist::parameter(): bad arg2");
      }
   }
   else {
      throw exception("BinomialDist::parameter(): bad arg1");
   }
}

///////////////////  PoissonDist class /////////////////////////////

void PoissonDist::reset(const double l)
{
   if (l < 0.0) throw exception("PoissonDist::reset(): bad arg, nonnegative expected");
   lambda_value = l;
}

/*!
   The random variable <EM>X</EM> has a Poisson distribution if its probability density function (pdf) is defined by
\f[
  f(x) = \frac{\lambda^xe^{-\lambda}}{x!}, \quad x = 0,1,2,\ldots,
\f]
*/
double PoissonDist::pdf(const long k) const
{
   if (k < 0L) throw exception ("PoissonDist::pdf(): out of range");
   return std::pow(lambda_value,static_cast<double>(k))/(std::exp(lambda_value)*factrl(k));
}

double PoissonDist::cdf(const long k) const
{
   // ********************************************************
   // Pr(X <= k) = sum{Poisson(lambda)), over x=0,1,...k}
   // *******************************************************
   double retval;
   if (k > 0L)  {
      retval = 1.0 - gammin(lambda_value,k+1);
   }
   else if (k == 0L) {
      retval = std::exp(-lambda_value);
   }
   else {
      retval = 0.0;
   }
   return retval;
}

void PoissonDist::sample(Vector<double>& x) const
{
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = static_cast<double>(ignpoi(lambda_value));
}

void PoissonDist::sample(doubleMatrix& x) const
{
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = static_cast<double>(ignpoi(lambda_value));
   }
}

Vector<double> PoissonDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix PoissonDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double PoissonDist::parameter(const int k) const
{
   if (k != 1) throw exception("PoissonDist::parameter(): bad arg, 1 is expected");
   return lambda_value;
}

void PoissonDist::parameter(const int k,const double x)
{
   if (k==1) {
      if (x >= 0.0) {
         lambda_value = x;
      }
      else {
         throw exception("PoissonDist::parameter(): bad arg2");
      }
   }
   else {
      throw exception("PoissonDist::parameter(): bad arg1");
   }
}

///////////////////  GeometricDist class /////////////////////////////

void GeometricDist::reset(const double p)
{
   if (p <= 0.0 && p >= 1.0) throw exception ("GeometricDist::reset(): bad arg");
   p_value = p;
}

/*!
The random variable <EM>X</EM> has a geometric distribution if its probability density function (pdf) is defined by
\f[
  f(x) = (1-p)^x p, \quad x = 0,1,2,\ldots,
\f]
*/
double GeometricDist::pdf(const long k) const
{
   if (k < 0L) throw exception("GeometricDist::pdf(): arg must be non-negative");
   return pow(1.0 - p_value,static_cast<double>(k))*p_value;
}

double GeometricDist::cdf(const long k) const
{
   // ********************************************************
   // Pr(X <= k) = sum{Geometric(p)), over x=0,1,...k}
   // *******************************************************
   if (k < 0L) throw exception("GeometricDist::cd): arg must be non-negative");
   return pow(1.0-p_value,static_cast<double>(k+1));
}

double GeometricDist::mgf(const double t) const
{
   if (t >= -std::log(1.0-p_value)) throw exception ("GeometricDist::mgf(): bad arg");
   return p_value/(1.0-(1.0-p_value)*std::exp(t));
}

void GeometricDist::sample(Vector<double>& x) const
{
   double tmp = std::log(1.0-p_value);
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = std::floor(std::log(ranf())/tmp);
}

void GeometricDist::sample(doubleMatrix& x) const
{
   double tmp = std::log(1.0-p_value);
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = std::floor(std::log(ranf())/tmp);
   }
}

Vector<double> GeometricDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix GeometricDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double GeometricDist::parameter(const int k) const
{
   if (k != 1) throw exception("GeometricDist::parameter(): bad arg, 1 is expected");
   return p_value;
}

void GeometricDist::parameter(const int k,const double x)
{
   if (k != 1) throw exception("GeometricDist::parameter(): bad arg1");
   if (x <= 0.0 || x >= 1.0) throw exception("GeometricDist::parameter(): bad arg2");
   p_value = x;
}

///////////////////  NegativeBinomialDist class //////////////////////
void NegativeBinomialDist::reset(const long n,const double p)
{
   if (n < 1L || p <= 0.0 || p >= 1.0) throw exception("NegativeBinomialDist::reset(): bad args");
   n_value = n;
   p_value = p;
}

/*!
   The random variable <EM>X</EM> has a negative binomial distribution if its probability density function (pdf) is defined by
\f[
  f(x) = \frac{x+n-1)!}{x!(n-1)!}p^n(1-p)^x, \quad x = 0,1,2,\ldots,
\f]
*/
double NegativeBinomialDist::pdf(const long k) const
{
   if (k < 0L) throw exception("NegativeBinomialDist::pdf(): arg must be non-negative");
   double double_n = static_cast<double>(n_value);
   double double_k = static_cast<double>(k);
   double retval = gammln(double_n + double_k) - gammln(double_n)
               - gammln(double_k + 1.0);
   retval += double_n*std::log(p_value) + double_k*std::log(1.0-p_value);
   if (retval < -700.0) return 0.0;
   else return std::exp(retval);
}

/*!
 *  Pr(X <= k) = sum{NegativeBinomial(n,p)), over x=0,1,...k}
 */
double NegativeBinomialDist::cdf(const long k) const
{
   double retval = 0.0;
   for (long i=0; i<=k; i++) retval += pdf(i);
   return retval;
}

double NegativeBinomialDist::mgf(const double t) const
{
   if (t >= -std::log(1.0-p_value)) throw exception("NegativeBinomialDist::mgf(): bad arg");
   double retval = pow(p_value,static_cast<double>(n_value));
   retval /= std::pow(1.0-(1.0-p_value)*std::exp(t),static_cast<double>(n_value));
   return retval;
}

double NegativeBinomialDist::sample(void) const
{
   double x = 0.0;
   double tmp = std::log(1.0-p_value);
   for (long i=0; i<n_value; i++) {
      x += std::floor(std::log(ranf())/tmp);   // x = sum{(GD~Geometric(p)) over k=1..r}
   }
   return x;
}

void NegativeBinomialDist::sample(Vector<double>& x) const
{
   double tmp = std::log(1.0-p_value);
   double *bot = x.begin();
   double *top = x.end();
   while (bot < top ) *bot++ = sample();
}

void NegativeBinomialDist::sample(doubleMatrix& x) const
{
   double tmp = std::log(1.0-p_value);
   int nr = x.num_rows();
   int nc = x.num_cols();
   double *bot,*top;
   for (unsigned i=0; i<nr; i++) {
      bot = x[i];  top = &bot[nc];
      while (bot < top ) *bot++ = sample();
   }
}

Vector<double> NegativeBinomialDist::sample(unsigned n) const
{
   Vector<double> x(n);
   sample(x);
   return x;
}

doubleMatrix NegativeBinomialDist::sample(unsigned m,unsigned n) const
{
   doubleMatrix x(m,n);
   sample(x);
   return x;
}

double NegativeBinomialDist::parameter(const int k) const
{
   double retval =  0.0;
   if (k==1) {
      retval = static_cast<double>(n_value);
   }
   else if (k==2) {
      retval = p_value;
   }
   else {
      throw exception("NegativeBinomialDist::parameter(): bad arg, 1 or 2 expected");
   }
   return retval;
}

void NegativeBinomialDist::parameter(const int k,const double x)
{
   if (k==1) {
      if (x >= 1.0) {
         n_value = static_cast<long int>(x);
      }
      else {
         throw exception("NegativeBinomialDist::parameter(): bad arg2");
      }
   }
   else if (k==2) {
      if (x > 0.0 && x < 1.0) {
         p_value = x;
      }
      else {
         throw exception("NegativeBinomialDist::parameter(): bad arg2");
      }
   }
   else {
      throw exception("NegativeBinomialDist::parameter(): bad arg1");
   }
}
} ////// end of namespace matvec

