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
#include <cmath>
#include <algorithm>
#include "exception.h"
#include "bg.h"

/******************************************************************
*      ginv for positive semi-definite(psd) symmetric matrix
*                                a(0..n-1,0..n-1)
*      mode = 1 inverse, a returns as ginv1(a)
*             2 cholesky only, return a=(L\0)
*             3 log determinent only, a  is destroyed
*     NOTE: user is totally resposible to check if a is symmetric
*******************************************************************/

namespace matvec {
unsigned ginverse1(double **a,const unsigned n,double& lgdet,int mode,
                  const double tol)
{
   unsigned i,j,k,irank = n;
   double sum;
   ///////////////////////////////////////////////////
   // column-by-column, we are doing
   // L (lower triangular) in place such that L*L'=A
   ///////////////////////////////////////////////////
   lgdet = 0.0;
   for (j=0; j<n; ++j) {  // for column j, L[j][j] is calculated first
      sum = a[j][j];
      for (k=0; k<j; ++k) sum -= a[j][k]*a[j][k];
                        //  for other elements L[i,j] in jth column, i=j+1 to n
      if (sum <= -tol) {
	throw exception("ginverse1(): matrix is non-psd");
      }
      if (std::fabs(sum) < tol) {
         irank--;
         for (i=j; i<n; i++) a[i][j]=0.0;
      }
      else {
         lgdet += std::log(sum);
         a[j][j]=std::sqrt(sum);
         for (i=j+1; i<n; i++) {
            sum = a[i][j];
            for (k=0; k<j; ++k) sum -= a[i][k]*a[j][k];
            a[i][j]=sum/a[j][j];
         }
      }
   }
   if (mode != 1) {
      if (mode ==2) {
         for(i=0;i<n-1;i++) memset(&(a[i][i+1]),'\0',sizeof(double)*(n-i-1));
       }
       return irank;
   }
   //////////////////////////////////////////
   //  column-by-column,  we are doing 
   //   inv(L) in place (lower triangl)
   //////////////////////////////////////////
   for (j=0; j<n; j++) { 
      if (std::fabs(a[j][j]) < tol) {
         for (i=j; i<n; i++) a[i][j]=0.0;            
      }
      else {
         a[j][j] = 1./a[j][j];
         for (i=j+1; i<n; i++) {
            if (a[i][i] <= -tol)  throw exception("ginverse1(): matrix is non-psd");
            if (std::fabs(a[i][i]) < tol) {
               a[i][j]=0.0;
            }
            else {
               sum=0.0;
               for (k=j; k<i; k++) sum -= a[i][k]*a[k][j];                      
               a[i][j] = sum/a[i][i];
            }      // end of if-else block
         }         // end of for-loop
      }            // end of if-else block
   }               // end of for-loop

   //////////////////////////////////////
   //  column-by-column, we are doing 
   //  inv(L')*inv(L)= inv(a) in place
   //////////////////////////////////////
   for (j=0; j<n; j++) {
      for (i=j; i<n; i++) {
         for (sum=0.0,k=i; k<n; k++) sum += a[k][i]*a[k][j];    
         a[i][j]=sum;    a[j][i]=sum;
      }
   }
   return irank;
}

unsigned ginverse1(BG **a,const unsigned n,BG& lgdet,int mode,
                  const double tol)
{
   unsigned i,j,k,irank = n;
   BG sum;
   ///////////////////////////////////////////////////
   // column-by-column, we are doing
   // L (lower triangular) in place such that L*L'=A
   ///////////////////////////////////////////////////
   lgdet = 0.0;
   for (j=0; j<n; ++j) {  // for column j, L[j][j] is calculated first
      sum = a[j][j];
      for (k=0; k<j; ++k) sum -= a[j][k]*a[j][k];
                        //  for other elements L[i,j] in jth column, i=j+1 to n
      if (sum <= -tol) {
	throw exception("ginverse1(): matrix is non-psd");
      }
      if (fabs(sum) < tol) {
         irank--;
         for (i=j; i<n; i++) a[i][j]=0.0;
      }
      else {
         lgdet += log(sum);
         a[j][j]=sqrt(sum);
         for (i=j+1; i<n; i++) {
            sum = a[i][j];
            for (k=0; k<j; ++k) sum -= a[i][k]*a[j][k];
            a[i][j]=sum/a[j][j];
         }
      }
   }
   if (mode != 1) {
      if (mode ==2) {
	for(i=0;i<n-1;i++) { 
	  //memset(&(a[i][i+1]),'\0',sizeof(BG)*(n-i-1));
	  for(j=i+1;j<n;j++){
	    a[i][j] = 0.0;
	  }
	}
       }
       return irank;
   }
   //////////////////////////////////////////
   //  column-by-column,  we are doing 
   //   inv(L) in place (lower triangl)
   //////////////////////////////////////////
   for (j=0; j<n; j++) { 
      if (fabs(a[j][j]) < tol) {
         for (i=j; i<n; i++) a[i][j]=0.0;            
      }
      else {
         a[j][j] = 1./a[j][j];
         for (i=j+1; i<n; i++) {
            if (a[i][i] <= -tol)  throw exception("ginverse1(): matrix is non-psd");
            if (fabs(a[i][i]) < tol) {
               a[i][j]=0.0;
            }
            else {
               sum=0.0;
               for (k=j; k<i; k++) sum -= a[i][k]*a[k][j];                      
               a[i][j] = sum/a[i][i];
            }      // end of if-else block
         }         // end of for-loop
      }            // end of if-else block
   }               // end of for-loop

   //////////////////////////////////////
   //  column-by-column, we are doing 
   //  inv(L')*inv(L)= inv(a) in place
   //////////////////////////////////////
   for (j=0; j<n; j++) {
      for (i=j; i<n; i++) {
         for (sum=0.0,k=i; k<n; k++) sum += a[k][i]*a[k][j];    
         a[i][j]=sum;    a[j][i]=sum;
      }
   }
   return irank;
}


bool psdefinite(const double **mat,const unsigned m,const unsigned n,const double tol)
{
   /////// first check to see if it is symmetric or not
   if (m != n) return false;
   int i,j,t;
   bool k = true;
   for (i=1; i<m; i++) {
      for (j=0; j<i; j++) if (std::fabs(mat[i][j] - mat[j][i]) > tol) return false;
   }

   double **a = new double *[m];
   for (i=0; i<m; ++i) a[i] = new double [n];

   for (i=0; i<m; i++) memcpy(a[i],mat[i],sizeof(double)*n);
   double sum;
   for (j=0; j<n; j++) {
      sum = a[j][j];
      for (t=0; t<j; ++t) sum -= a[j][t]*a[j][t];
      if (sum <= -tol) {
         k = false;
         break;
      }
      else if (std::fabs(sum) < tol) {
         for (i=j; i<n; i++) a[i][j]=0.0;
      }
      else {
         a[j][j]=std::sqrt(sum);
         for (i=j+1; i<n; i++) {
            sum = a[i][j];
            for (t=0; t<j; ++t) sum -= a[i][t]*a[j][t];
            a[i][j]=sum/a[j][j];
         }
      }
   }
   for (i=0; i<m; ++i) {
     if(a[i]){
       delete [] a[i];
       a[i]=0;
     }
   }
   if(a){
     delete [] a;
     a=0;
   }
   return k;
}

bool psdefinite(const BG **mat,const unsigned m,const unsigned n,const double tol)
{
   /////// first check to see if it is symmetric or not
   if (m != n) return false;
   int i,j,t;
   bool k = true;
   for (i=1; i<m; i++) {
      for (j=0; j<i; j++) if (fabs(mat[i][j] - mat[j][i]) > tol) return false;
   }

   BG **a = new BG *[m];
   for (i=0; i<m; ++i) a[i] = new BG [n];

   //for (i=0; i<m; i++) memcpy(a[i],mat[i],sizeof(BG)*n);
   for (i=0; i<m; i++) {
     for (j=0; j<m; j++){
       a[i][j] = mat[i][j];
     }
   }
   BG sum;
   for (j=0; j<n; j++) {
      sum = a[j][j];
      for (t=0; t<j; ++t) sum -= a[j][t]*a[j][t];
      if (sum <= -tol) {
         k = false;
         break;
      }
      else if (fabs(sum) < tol) {
         for (i=j; i<n; i++) a[i][j]=0.0;
      }
      else {
         a[j][j]=sqrt(sum);
         for (i=j+1; i<n; i++) {
            sum = a[i][j];
            for (t=0; t<j; ++t) sum -= a[i][t]*a[j][t];
            a[i][j]=sum/a[j][j];
         }
      }
   }
   for (i=0; i<m; ++i) {
     if(a[i]){
       delete [] a[i];
       a[i]=0;
     }
   }
   if(a){
     delete [] a;
     a=0;
   }
   return k;
}

void ludcmp(double **a,int n,int *indx,int& d,double tol)
{
   int i,j,k,imax;
   double big,dum,sum,temp;
   double *vv = new double [n];
   d=1;
   for (i=0;i<n;i++) {
      big = 0.0;
      for (j=0;j<n;j++) if ((temp=std::fabs(a[i][j])) > big) big=temp;
      if (std::fabs(big) < tol) throw exception("ludcmp(): Singular matrix");
      vv[i] = 1.0/big;
   }
   for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
         sum = a[i][j];
         for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      big=0.0;
      imax = j;
      for (i=j;i<n;i++) {
         sum = a[i][j];
         for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
         if ( (dum = vv[i]*std::fabs(sum)) >= big) {
            big = dum;
            imax = i;
         }
      }
      if (j != imax) {
         for (k = 0; k<n; k++) {
            dum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = dum;
         }
         d = -d;
         vv[imax] = vv[j];
      }
      indx[j] = imax;
      if (std::fabs(a[j][j]) < tol) throw exception("ludcmp(): Singular matrix");
      dum = 1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
   }
   if(vv){
     delete []  vv;
     vv=0;
   }
}

void ludcmp(BG **a,int n,int *indx,int& d, double tol)
{
   int i,j,k,imax;
   BG big,dum,sum,temp;
   BG *vv = new BG [n];
   d=1;
   for (i=0;i<n;i++) {
      big = 0.0;
      for (j=0;j<n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
      if (fabs(big) < tol) throw exception("ludcmp(): Singular matrix");
      vv[i] = 1.0/big;
   }
   for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
         sum = a[i][j];
         for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      big=0.0;
      imax = j;
      for (i=j;i<n;i++) {
         sum = a[i][j];
         for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
         if ( (dum = vv[i]*fabs(sum)) >= big) {
            big = dum;
            imax = i;
         }
      }
      if (j != imax) {
         for (k = 0; k<n; k++) {
            dum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = dum;
         }
         d = -d;
         vv[imax] = vv[j];
      }
      indx[j] = imax;
      if (fabs(a[j][j]) < tol) throw exception("ludcmp(): Singular matrix");
      dum = 1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
   }
   if(vv){
     delete []  vv;
     vv=0;
   }
}

void lubksb(double **a,int n,int *indx,double *b)
{
   int i,j;
   int ip, ii = -1;
   double sum;

   for (i=0;i<n;i++) {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii>=0) {
         for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      }
      else if (sum) {
         ii=i;
      }
      b[i] = sum;
   }
   for (i=n-1; i>=0; i--) {
      sum = b[i];
      for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      b[i] = sum/a[i][i];
   }
}

void lubksb(BG **a,int n,int *indx,BG *b)
{
   int i,j;
   int ip, ii = -1;
   BG sum;

   for (i=0;i<n;i++) {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii>=0) {
         for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      }
      else if (sum!=0.0) {
         ii=i;
      }
      b[i] = sum;
   }
   for (i=n-1; i>=0; i--) {
      sum = b[i];
      for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
      b[i] = sum/a[i][i];
   }
}


// PYTHAG computes std::sqrt(a^{2} + b^{2}) without destructive overflow or underflow.

static double at, bt, ct;
#define PYTHAG(a, b) ((at = std::fabs(a)) > (bt = std::fabs(b)) ? \
(ct = bt/at, at*std::sqrt(1.0+ct*ct)): (bt ? (ct = at/bt, bt*std::sqrt(1.0+ct*ct)): 0.0))

static double maxarg1, maxarg2;
#define MAX(a, b) (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? \
(maxarg1) : (maxarg2))

#define SIGN(a, b) ((b) < 0.0 ? -std::fabs(a): std::fabs(a))

void svdcmp(double* a[], int m, int n, double w[], double* v[])
{
   /////////////////////////////////////////////////////////////////////
   // Given a matrix a[m][n], this routine computes its singular value
   // decomposition, A = U*W*V'.  The matrix U replaces a on output.
   // The diagonal matrix of singular values W is output as a vector w[n].
   // The matrix V  is output as v[n][n].
   // m must be greater or equal to n;  if it is smaller, then a should be
   // filled up to square with zero rows.
   /////////////////////////////////////////////////////////////////////////

   int flag, i, its, j, jj, k, l, nm;
   double c, f, h, s, x, y, z;
   double anorm = 0.0, g = 0.0, scale = 0.0;

   if (m < n) throw exception("svdcmp(): Matrix A  was not augmented with extra rows of zeros");

   double* rv1 = new double [n];

    // Householder reduction to bidiagonal form.
   l = 0;    // added by T. Wang to avoid warning in g++
   nm = 0;   // added by T. Wang to avoid warning in g++
   for (i = 0; i < n; i++) {
      l = i + 1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
      if (i < m) {
         for (k = i; k < m; k++) scale += std::fabs(a[k][i]);
            if (scale) {
               for (k = i; k < m; k++) {
                  a[k][i] /= scale;
                  s += a[k][i]*a[k][i];
               }
               f = a[i][i];
               g = -SIGN(std::sqrt(s), f);
               h = f*g - s;
               a[i][i] = f - g;
               if (i != n - 1) {
                  for (j = l; j < n; j++) {
                     for (s  = 0.0, k = i; k < m; k++) s += a[k][i]*a[k][j];
                     f = s/h;
                     for ( k = i; k < m; k++) a[k][j] += f*a[k][i];
                  }
               }
               for (k = i; k < m; k++) a[k][i] *= scale;
            }
         }
         w[i] = scale*g;
         g = s= scale = 0.0;
         if (i < m && i != n - 1) {
            for (k = l; k < n; k++)  scale += std::fabs(a[i][k]);
            if (scale) {
               for (k = l; k < n; k++) {
                  a[i][k] /= scale;
                  s += a[i][k]*a[i][k];
               }
               f = a[i][l];
               g = -SIGN(std::sqrt(s), f);
               h = f*g - s;
               a[i][l] = f - g;
	       //	       if(l==0){
	       //		 cout << "l eq 0" << endl;
	       //	       }
               for (k = l; k < n; k++)  rv1[k] = a[i][k]/h;
               if (i != m - 1) {
                  for (j = l; j < m; j++) {
                     for (s = 0.0, k = l; k < n; k++) s += a[j][k]*a[i][k];
                     for (k = l; k < n; k++) a[j][k] += s*rv1[k];
                  }
               }
               for (k = l; k < n; k++) a[i][k] *= scale;
            }
         }
         anorm = MAX(anorm, (std::fabs(w[i]) + std::fabs(rv1[i])));
      }
        /* Accumulation of right-hand transformations.	*/
      for (i = n - 1; 0 <= i; i--) {
         if (i < n - 1) {
	    if (g) {
               for (j = l; j < n; j++)  v[j][i] = (a[i][j]/a[i][l])/g;
		  /* Double division to avoid possible underflow: */
               for (j = l; j < n; j++) {
                  for (s = 0.0, k = l; k < n; k++) s += a[i][k]*v[k][j];
                  for (k = l; k < n; k++)  v[k][j] += s*v[k][i];
               }
	    }
	    for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
         }
         v[i][i] = 1.0;
         g = rv1[i];
         l = i;
      }
          /* Accumulation of left-hand transformations.	  */
      for (i = n - 1; 0 <= i; i--) {
         l = i + 1;
         g = w[i];
         if (i < n - 1) for (j = l; j < n; j++) a[i][j] = 0.0;
         if (g) {
            g = 1.0/g;
            if (i != n - 1) {
               for (j = l; j < n; j++) {
		  for (s = 0.0, k = l; k < m; k++) s += a[k][i]*a[k][j];
                  f = (s/a[i][i])*g;
                  for (k = i; k < m; k++) a[k][j] += f*a[k][i];
               }
	    }
            for (j = i; j < m; j++)  a[j][i] *= g;
         }
         else {
            for (j = i; j < m; j++) a[j][i] = 0.0;
         }
         a[i][i] += 1.0;   // ++a[i][i]
      }
       /* Diagonalization of the bidiagonal form.  */
      for (k = n - 1; 0 <= k; k--) {        /* Loop over singular values. */
         for (its = 0; its < 30; its++) {    /* Loop over allowed iterations.*/
            flag = 1;
            for (l = k; 0 <= l; l--) {     // Test for splitting:
               nm = l - 1;                 // Note that rv1[1] is always zero
               if ( (std::fabs(rv1[l]) + anorm) == anorm) {
                  flag = 0;
                  break;
               }
               if ( (std::fabs(w[nm]) + anorm) == anorm) break;
            }
            if (flag) {
	      if(l==0) {
		std::cout << " l=0 " <<flag << "\n";
	      }
            c = 0.0;	                   /* Cancellation of rv1[l], if l>0:*/
            s = 1.0;
            for (i = l; i <= k; i++) {
               f = s*rv1[i];
               if ((std::fabs(f) + anorm) != anorm) {
                  g = w[i];
                  h = PYTHAG(f, g);
                  w[i] = h;
                  h = 1.0/h;
                  c = g*h;
                  s = (-f*h);
                  for (j = 0; j < m; j++) {
		    if(nm < 0 || nm >= m){
		      std::cout << nm << " " << m <<" " << l << "\n";
		    }
                     y = a[j][nm];
                     z = a[j][i];
                     a[j][nm] = y*c + z*s;
                     a[j][i]  = z*c - y*s;
                  }
               }
            }
         }
         z = w[k];
         if (l == k) {	     /* Convergence.  */
            if (z < 0.0) {        /* Singular value is made non-negative. */
               w[k] = -z;
               for (j = 0; j < n; j++) v[j][k] = (-v[j][k]);
            }
            break;
         }
         if (its == 29) {
            delete [] rv1;
            exception("svdcmp(): Not convergence in 30 SVDCMP iterations");
         }
         x = w[l];               /* Shift from bottom 2-by-2 minor. */
         nm = k - 1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y);
         g = PYTHAG(f, 1.0);
         f = ((x - z)*(x + z) + h*((y/(f + SIGN(g, f))) - h))/x;
	    /* Next QR transformation:	  */
         c = s = 1.0;
         for (j = l; j <= nm; j++) {
            i = j + 1;
            g = rv1[i];
            y = w[i];
            h = s*g;
            g = c*g;
            z = PYTHAG(f, h);
            rv1[j] = z;
            c = f/z;
            s = h/z;
            f = x*c + g*s;
            g = g*c - x*s;
            h = y*s;
            y = y*c;
            for (jj = 0; jj < n;  jj++) {
               x = v[jj][j];
               z = v[jj][i];
               v[jj][j] = x*c + z*s;
               v[jj][i] = z*c - x*s;
            }
            z = PYTHAG(f, h);
            w[j] = z;        /* Rotation can be arbitrary if z = 0.*/
            if (z) {
               z = 1.0/z;
               c = f*z;
               s = h*z;
            }
            f = (c*g) + (s*y);
            x = (c*y) - (s*g);
            for (jj = 0; jj < m; jj++) {
               y = a[jj][j];
               z = a[jj][i];
               a[jj][j] = y*c + z*s;
               a[jj][i] = z*c - y*s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
   }
      if(rv1){
	delete []  rv1;
	rv1=0;
      }
}

// // PYTHAG computes std::sqrt(a^{2} + b^{2}) without destructive overflow or underflow.

// static BG at, bt, ct;
// #define BG_PYTHAG(a, b) ((at = fabs(a)) > (bt = fabs(b)) ? \
// (ct = bt/at, at*sqrt(1.0+ct*ct)): (bt ? (ct = at/bt, bt*sqrt(1.0+ct*ct)): 0.0))

// static double maxarg1, maxarg2;
// #define BG_MAX(a, b) (maxarg1 = (a), maxarg2 = (b), (maxarg1) > (maxarg2) ? \
// (maxarg1) : (maxarg2))

// #define BG_SIGN(a, b) ((b) < 0.0 ? -fabs(a): fabs(a))

// void svdcmp(double* a[], int m, int n, double w[], double* v[])
// {
//    /////////////////////////////////////////////////////////////////////
//    // Given a matrix a[m][n], this routine computes its singular value
//    // decomposition, A = U*W*V'.  The matrix U replaces a on output.
//    // The diagonal matrix of singular values W is output as a vector w[n].
//    // The matrix V  is output as v[n][n].
//    // m must be greater or equal to n;  if it is smaller, then a should be
//    // filled up to square with zero rows.
//    /////////////////////////////////////////////////////////////////////////

//    int flag, i, its, j, jj, k, l, nm;
//    double c, f, h, s, x, y, z;
//    double anorm = 0.0, g = 0.0, scale = 0.0;

//    if (m < n) throw exception("svdcmp(): Matrix A  was not augmented with extra rows of zeros");

//    double* rv1 = new double [n];

//     // Householder reduction to bidiagonal form.
//    l = 0;    // added by T. Wang to avoid warning in g++
//    nm = 0;   // added by T. Wang to avoid warning in g++
//    for (i = 0; i < n; i++) {
//       l = i + 1;
//       rv1[i] = scale*g;
//       g = s = scale = 0.0;
//       if (i < m) {
//          for (k = i; k < m; k++) scale += fabs(a[k][i]);
//             if (scale) {
//                for (k = i; k < m; k++) {
//                   a[k][i] /= scale;
//                   s += a[k][i]*a[k][i];
//                }
//                f = a[i][i];
//                g = -BG_SIGN(std::sqrt(s), f);
//                h = f*g - s;
//                a[i][i] = f - g;
//                if (i != n - 1) {
//                   for (j = l; j < n; j++) {
//                      for (s  = 0.0, k = i; k < m; k++) s += a[k][i]*a[k][j];
//                      f = s/h;
//                      for ( k = i; k < m; k++) a[k][j] += f*a[k][i];
//                   }
//                }
//                for (k = i; k < m; k++) a[k][i] *= scale;
//             }
//          }
//          w[i] = scale*g;
//          g = s= scale = 0.0;
//          if (i < m && i != n - 1) {
//             for (k = l; k < n; k++)  scale += std::fabs(a[i][k]);
//             if (scale) {
//                for (k = l; k < n; k++) {
//                   a[i][k] /= scale;
//                   s += a[i][k]*a[i][k];
//                }
//                f = a[i][l];
//                g = -BG_SIGN(std::sqrt(s), f);
//                h = f*g - s;
//                a[i][l] = f - g;
// 	       //	       if(l==0){
// 	       //		 cout << "l eq 0" << endl;
// 	       //	       }
//                for (k = l; k < n; k++)  rv1[k] = a[i][k]/h;
//                if (i != m - 1) {
//                   for (j = l; j < m; j++) {
//                      for (s = 0.0, k = l; k < n; k++) s += a[j][k]*a[i][k];
//                      for (k = l; k < n; k++) a[j][k] += s*rv1[k];
//                   }
//                }
//                for (k = l; k < n; k++) a[i][k] *= scale;
//             }
//          }
//          anorm = BG_MAX(anorm, (std::fabs(w[i]) + std::fabs(rv1[i])));
//       }
//         /* Accumulation of right-hand transformations.	*/
//       for (i = n - 1; 0 <= i; i--) {
//          if (i < n - 1) {
// 	    if (g) {
//                for (j = l; j < n; j++)  v[j][i] = (a[i][j]/a[i][l])/g;
// 		  /* Double division to avoid possible underflow: */
//                for (j = l; j < n; j++) {
//                   for (s = 0.0, k = l; k < n; k++) s += a[i][k]*v[k][j];
//                   for (k = l; k < n; k++)  v[k][j] += s*v[k][i];
//                }
// 	    }
// 	    for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
//          }
//          v[i][i] = 1.0;
//          g = rv1[i];
//          l = i;
//       }
//           /* Accumulation of left-hand transformations.	  */
//       for (i = n - 1; 0 <= i; i--) {
//          l = i + 1;
//          g = w[i];
//          if (i < n - 1) for (j = l; j < n; j++) a[i][j] = 0.0;
//          if (g) {
//             g = 1.0/g;
//             if (i != n - 1) {
//                for (j = l; j < n; j++) {
// 		  for (s = 0.0, k = l; k < m; k++) s += a[k][i]*a[k][j];
//                   f = (s/a[i][i])*g;
//                   for (k = i; k < m; k++) a[k][j] += f*a[k][i];
//                }
// 	    }
//             for (j = i; j < m; j++)  a[j][i] *= g;
//          }
//          else {
//             for (j = i; j < m; j++) a[j][i] = 0.0;
//          }
//          a[i][i] += 1.0;   // ++a[i][i]
//       }
//        /* Diagonalization of the bidiagonal form.  */
//       for (k = n - 1; 0 <= k; k--) {        /* Loop over singular values. */
//          for (its = 0; its < 30; its++) {    /* Loop over allowed iterations.*/
//             flag = 1;
//             for (l = k; 0 <= l; l--) {     // Test for splitting:
//                nm = l - 1;                 // Note that rv1[1] is always zero
//                if ( (std::fabs(rv1[l]) + anorm) == anorm) {
//                   flag = 0;
//                   break;
//                }
//                if ( (std::fabs(w[nm]) + anorm) == anorm) break;
//             }
//             if (flag) {
// 	      if(l==0) {
// 		std::cout << " l=0 " <<flag << "\n";
// 	      }
//             c = 0.0;	                   /* Cancellation of rv1[l], if l>0:*/
//             s = 1.0;
//             for (i = l; i <= k; i++) {
//                f = s*rv1[i];
//                if ((std::fabs(f) + anorm) != anorm) {
//                   g = w[i];
//                   h = BG_PYTHAG(f, g);
//                   w[i] = h;
//                   h = 1.0/h;
//                   c = g*h;
//                   s = (-f*h);
//                   for (j = 0; j < m; j++) {
// 		    if(nm < 0 || nm >= m){
// 		      std::cout << nm << " " << m <<" " << l << "\n";
// 		    }
//                      y = a[j][nm];
//                      z = a[j][i];
//                      a[j][nm] = y*c + z*s;
//                      a[j][i]  = z*c - y*s;
//                   }
//                }
//             }
//          }
//          z = w[k];
//          if (l == k) {	     /* Convergence.  */
//             if (z < 0.0) {        /* Singular value is made non-negative. */
//                w[k] = -z;
//                for (j = 0; j < n; j++) v[j][k] = (-v[j][k]);
//             }
//             break;
//          }
//          if (its == 29) {
//             delete [] rv1;
//             exception("svdcmp(): Not convergence in 30 SVDCMP iterations");
//          }
//          x = w[l];               /* Shift from bottom 2-by-2 minor. */
//          nm = k - 1;
//          y = w[nm];
//          g = rv1[nm];
//          h = rv1[k];
//          f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y);
//          g = BG_PYTHAG(f, 1.0);
//          f = ((x - z)*(x + z) + h*((y/(f + BG_SIGN(g, f))) - h))/x;
// 	    /* Next QR transformation:	  */
//          c = s = 1.0;
//          for (j = l; j <= nm; j++) {
//             i = j + 1;
//             g = rv1[i];
//             y = w[i];
//             h = s*g;
//             g = c*g;
//             z = BG_PYTHAG(f, h);
//             rv1[j] = z;
//             c = f/z;
//             s = h/z;
//             f = x*c + g*s;
//             g = g*c - x*s;
//             h = y*s;
//             y = y*c;
//             for (jj = 0; jj < n;  jj++) {
//                x = v[jj][j];
//                z = v[jj][i];
//                v[jj][j] = x*c + z*s;
//                v[jj][i] = z*c - x*s;
//             }
//             z = BG_PYTHAG(f, h);
//             w[j] = z;        /* Rotation can be arbitrary if z = 0.*/
//             if (z) {
//                z = 1.0/z;
//                c = f*z;
//                s = h*z;
//             }
//             f = (c*g) + (s*y);
//             x = (c*y) - (s*g);
//             for (jj = 0; jj < m; jj++) {
//                y = a[jj][j];
//                z = a[jj][i];
//                a[jj][j] = y*c + z*s;
//                a[jj][i] = z*c - y*s;
//             }
//          }
//          rv1[l] = 0.0;
//          rv1[k] = f;
//          w[k] = x;
//       }
//    }
//       if(rv1){
// 	delete []  rv1;
// 	rv1=0;
//       }
// }



void sweep(const int m, const int n, double **a, const int i0,
                  const int i1,const double tol)
/******************************************************************
    sweeps matrix A on the pivots indicated by i0 through i1
    to produce a new matrix.
       for example, suppose that A is partitioned into
          [ R S
            T U ]
      such that R is q by q and nonsingular, U is m-q by n-q.
      then sweep(m,n,A,1,q,tol) returns

        /    -1        -1      \
        |   R         R  S     |
        |      -1          -1  |
        \  -T R       U-T R  S /

  Jan. 20,1992, University of Illinois
  Copyright (C) 1992-1995 Tianlin Wang
********************************************************************/
{
   if (i0>i1 || i0 > std::min(m,n)) exception("sweep(): inappropriate values of arguments");

   int     i,j,k;
   double  d,b, *temp,*ak, *aktemp;

   for (k=i0; k<=i1; k++) {
      ak = a[k];
      d = ak[k];
      if (std::fabs(d) < tol) {
         for (i=0;i<m; i++) a[i][k] = 0.0;
         memset(ak,'\0',sizeof(double)*n);
      }
      else {
         aktemp = ak;
         for (j=0; j<n; j++) *aktemp++ /= d;
         for (i=0; i<m; i++) {
            aktemp = ak;
            if (i != k) {
               b = a[i][k];
               temp = a[i];
               for (j=0; j<n; j++) *temp++ -= *aktemp++ *b;
               a[i][k] = -b/d;
            }
         }
         ak[k] = 1.0/d;
      }
   }
}

void sweep(const int m, const int n, BG **a, const int i0,
                  const int i1,const BG tol)
/******************************************************************
    sweeps matrix A on the pivots indicated by i0 through i1
    to produce a new matrix.
       for example, suppose that A is partitioned into
          [ R S
            T U ]
      such that R is q by q and nonsingular, U is m-q by n-q.
      then sweep(m,n,A,1,q,tol) returns

        /    -1        -1      \
        |   R         R  S     |
        |      -1          -1  |
        \  -T R       U-T R  S /

  Jan. 20,1992, University of Illinois
  Copyright (C) 1992-1995 Tianlin Wang
********************************************************************/
{
   if (i0>i1 || i0 > std::min(m,n)) exception("sweep(): inappropriate values of arguments");

   int     i,j,k;
   BG  d, b, *temp,*ak, *aktemp;

   for (k=i0; k<=i1; k++) {
      ak = a[k];
      d = ak[k];
      if (fabs(d) < tol) {
         for (i=0;i<m; i++) a[i][k] = 0.0;
	 //         memset(ak,'\0',sizeof(BG)*n);
     	 initial_BGarray(ak,n,0.0);
      }
      else {
         aktemp = ak;
         for (j=0; j<n; j++) *aktemp++ /= d;
         for (i=0; i<m; i++) {
            aktemp = ak;
            if (i != k) {
               b = a[i][k];
               temp = a[i];
               for (j=0; j<n; j++) *temp++ -= *aktemp++ *b;
               a[i][k] = -b/d;
            }
         }
         ak[k] = 1.0/d;
      }
   }
}
}  /////////// end of namespace

