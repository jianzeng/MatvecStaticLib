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
#include <complex>
#include "util.h"

namespace matvec {

void balance(double **mat, const int n, int *low, int *hi, double *scale,
             const double base)
{
   ////////////////////////////////////////////////////////////////////
   // Balance a matrix for calculation of eigenvalues and eigenvectors
   //
   ////////////////////////////////////////////////////////////////////

   double c,f,g,r,s;
   int i,j,k,l,done;
        /* search for rows isolating an eigenvalue and push them down */
   for (k = n - 1; k >= 0; k--) {
      for (j = k; j >= 0; j--) {
         for (i = 0; i <= k; i++) {
            if (i != j && fabs(mat[j][i]) != 0) break;
         }

	 if (i > k) {
            scale[k] = j;
            if (j != k) {
               for (i = 0; i <= k; i++) {
                  c = mat[i][j];
                  mat[i][j] = mat[i][k];
                  mat[i][k] = c;
               }

               for (i = 0; i < n; i++) {
                  c = mat[j][i];
                  mat[j][i] = mat[k][i];
                  mat[k][i] = c;
               }
            }
            break;
         }
      }
      if (j < 0) break;
   }

    /* search for columns isolating an eigenvalue and push them left */

   for (l = 0; l <= k; l++) {
      for (j = l; j <= k; j++) {
         for (i = l; i <= k; i++) {
            if (i != j && fabs(mat[i][j]) != 0) break;
	 }

         if (i > k) {
            scale[l] = j;

            if (j != l) {
               for (i = 0; i <= k; i++) {
                  c = mat[i][j];
                  mat[i][j] = mat[i][l];
                  mat[i][l] = c;
               }
               for (i = l; i < n; i++) {
                  c = mat[j][i];
                  mat[j][i] = mat[l][i];
                  mat[l][i] = c;
               }
            }
            break;
         }
      }
      if (j > k) break;
   }

   *hi = k;
   *low = l;

    /* balance the submatrix in rows l through k */

   for (i = l; i <= k; i++) scale[i] = 1;
   do {
      for (done = 1,i = l; i <= k; i++) {
         for (c = 0,r = 0,j = l; j <= k; j++) {
            if (j != i) {
               c += fabs(mat[j][i]);
               r += fabs(mat[i][j]);
            }
         }

         if (c != 0 && r != 0) {
            g = r / base;
            f = 1;
            s = c + r;

            while (c < g) {
               f *= base;
               c *= base * base;
            }
            g = r * base;
            while (c >= g) {
               f /= base;
               c /= base * base;
            }

            if ((c + r) / f < 0.95 * s) {
               done = 0;
               g = 1 / f;
               scale[i] *= f;

               for (j = l; j < n; j++) mat[i][j] *= g;
               for (j = 0; j <= k; j++) mat[j][i] *= f;
            }
         }
      }
   } while (!done);
}


void elemhess(int job,double **mat,const int n,const int low,const int hi,
              double **vr, double **vi)
{
   ////////////////////////////////////////////////////////////////
   // Reduce the submatrix in rows and columns low through hi of
   // real matrix mat to Hessenberg form by elementary similarity
   // transformations
   ////////////////////////////////////////////////////////////////
   int i,j,m;
   double x,y;
   int *array = 0;
   if (hi - low > 1) {
     if(hi>0){
       array = new int [hi];
     }
     else {
       array = 0;
     }
      check_ptr(array);
   }
   for (m = low + 1; m < hi; m++) {
      for (x = 0,i = m,j = m; j <= hi; j++) {
         if (fabs(mat[j][m-1]) > fabs(x)) {
            x = mat[j][m-1];
            i = j;
         }
      }

      if ((array[m] = i) != m) {
         for (j = m - 1; j < n; j++) {
            y = mat[i][j];
            mat[i][j] = mat[m][j];
            mat[m][j] = y;
	 }

	 for (j = 0; j <= hi; j++) {
            y = mat[j][i];
            mat[j][i] = mat[j][m];
            mat[j][m] = y;
	 }
      }

      if (x != 0) {
	 for (i = m + 1; i <= hi; i++) {
            if ((y = mat[i][m-1]) != 0) {
               y = mat[i][m-1] = y / x;
               for (j = m; j < n; j++) mat[i][j] -= y * mat[m][j];
               for (j = 0; j <= hi; j++) mat[j][m] += y * mat[j][i];
            }
         }
      }
   }
   if (job) {
      for (i=0; i<n; i++) {
         for (j=0; j<n; j++) vr[i][j] = 0.0; vi[i][j] = 0.0;
         vr[i][i] = 1.0;
      }

      for (m = hi - 1; m > low; m--) {
         for (i = m + 1; i <= hi; i++) vr[i][m] = mat[i][m-1];
         if ((i = array[m]) != m) {
            for (j = m; j <= hi; j++) {
               vr[m][j] = vr[i][j];
               vr[i][j] = 0.0;
            }
            vr[i][m] = 1.0;
         }
      }
   }
   if (array) {
     delete [] array;
     array=0;
   }
}


void unbalance(const int n,double **vr,double **vi, const int low,
               const int hi, double *scale)
{
   /////////////////////////////////////////////////////
   //  Transform back eigenvectors of a balanced matrix
   //  into the eigenvectors of the original matrix
   /////////////////////////////////////////////////////////

    int i,j,k;
     double tmp;
    for (i = low; i <= hi; i++) {
	for (j = 0; j < n; j++) {
	    vr[i][j] *= scale[i];
	    vi[i][j] *= scale[i];
	}
    }

    for (i = low - 1; i >= 0; i--) {
	if ((k = (int)scale[i]) != i) {
	    for (j = 0; j < n; j++) {
                tmp = vr[i][j];
                vr[i][j] = vr[k][j];
                vr[k][j] = tmp;

                tmp = vi[i][j];
                vi[i][j] = vi[k][j];
                vi[k][j] = tmp;	
	    }
	}
    }

    for (i = hi + 1; i < n; i++) {
	if ((k = (int)scale[i]) != i) {
	    for (j = 0; j < n; j++) {
                tmp = vr[i][j];
                vr[i][j] = vr[k][j];
                vr[k][j] = tmp;

                tmp = vi[i][j];
                vi[i][j] = vi[k][j];
                vi[k][j] = tmp;	
	    }
	}
    }
}


int nonsymm_eigen_core(const int job, double **mat, const int n, const int low,
                      const int hi, double *valr, double *vali, double **vr,
                      double **vi, const double eps,const int maxiter)
{
   ////////////////////////////////////////////////////////////////
   // Calculate eigenvalues and eigenvectors of a real upper
   // Hessenberg matrix
   // Return 1 if converges successfully and 0 otherwise
   /////////////////////////////////////////////////////////////////
   std::complex<double> v;
   double p,q,r,s,t,w,x,y,z,ra,sa,norm;
   int niter,en,i,j,k,l,m;

   for (i=0; i<n; i++) {
      valr[i]=0.0;
      vali[i]=0.0;
   }
      /* store isolated roots and calculate norm */
   for (norm = 0,i = 0; i < n; i++) {
      for (j = std::max(0,i-1); j < n; j++) {
         norm += fabs(mat[i][j]);
      }
      if (i < low || i > hi) valr[i] = mat[i][i];
   }
   t = 0;
   en = hi;

   while (en >= low) {
      niter = 0;
      for (;;) {

       /* look for single small subdiagonal element */

         for (l = en; l > low; l--) {
            s = fabs(mat[l-1][l-1]) + fabs(mat[l][l]);
            if (s == 0) s = norm;
            if (fabs(mat[l][l-1]) <= eps * s) break;
         }

         /* form shift */

         x = mat[en][en];

         if (l == en) {             /* one root found */
            valr[en] = x + t;
            if (job) mat[en][en] = x + t;
            en--;
            break;
         }

         y = mat[en-1][en-1];
         w = mat[en][en-1] * mat[en-1][en];

         if (l == en - 1) {                /* two roots found */
            p = (y - x) / 2;
            q = p * p + w;
            z = sqrt(fabs(q));
            x += t;
            if (job) {
               mat[en][en] = x;
               mat[en-1][en-1] = y + t;
            }
            if (q < 0) {                /* complex pair */
               valr[en-1] = x+p;
               vali[en-1] = z;
               valr[en] = x+p;
               vali[en] = -z;
            }
            else {                      /* real pair */
               z = (p < 0) ? p - z : p + z;
               valr[en-1] = x + z;
               valr[en] = (z == 0) ? x + z : x - w / z;
               if (job) {
                  x = mat[en][en-1];
                  s = fabs(x) + fabs(z);
                  p = x / s;
                  q = z / s;
                  r = sqrt(p*p+q*q);
                  p /= r;
                  q /= r;
                  for (j = en - 1; j < n; j++) {
                     z = mat[en-1][j];
                     mat[en-1][j] = q * z + p *
                     mat[en][j];
                     mat[en][j] = q * mat[en][j] - p*z;
                  }
                  for (i = 0; i <= en; i++) {
                     z = mat[i][en-1];
                     mat[i][en-1] = q * z + p * mat[i][en];
                     mat[i][en] = q * mat[i][en] - p*z;
                  }
                  for (i = low; i <= hi; i++) {
                     z = vr[i][en-1];
                     vr[i][en-1] = q*z + p*vr[i][en];
                     vr[i][en] = q*vr[i][en] - p*z;
                  }
               }
            }
            en -= 2;
            break;
         }
         if (niter == maxiter) return(-1);
         if (niter != 0 && niter % 10 == 0) {
            t += x;
            for (i = low; i <= en; i++) mat[i][i] -= x;
            s = fabs(mat[en][en-1]) + fabs(mat[en-1][en-2]);
            x = y = 0.75 * s;
            w = -0.4375 * s * s;
         }
         niter++;
           /* look for two consecutive small subdiagonal elements */
         for (m = en - 2; m >= l; m--) {
            z = mat[m][m];
            r = x - z;
            s = y - z;
            p = (r * s - w) / mat[m+1][m] + mat[m][m+1];
            q = mat[m+1][m+1] - z - r - s;
            r = mat[m+2][m+1];
            s = fabs(p) + fabs(q) + fabs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m == l || fabs(mat[m][m-1]) * (fabs(q)+fabs(r)) <=
                eps * (fabs(mat[m-1][m-1]) + fabs(z) +
                fabs(mat[m+1][m+1])) * fabs(p)) break;
         }
         for (i = m + 2; i <= en; i++) mat[i][i-2] = 0;
         for (i = m + 3; i <= en; i++) mat[i][i-3] = 0;
             /* double QR step involving rows l to en and columns m to en */
         for (k = m; k < en; k++) {
            if (k != m) {
               p = mat[k][k-1];
               q = mat[k+1][k-1];
               r = (k == en - 1) ? 0 : mat[k+2][k-1];
               if ((x = fabs(p) + fabs(q) + fabs(r)) == 0) continue;
               p /= x;
               q /= x;
               r /= x;
            }
            s = sqrt(p*p+q*q+r*r);
            if (p < 0) s = -s;
            if (k != m) {
               mat[k][k-1] = -s * x;
            }
            else if (l != m) {
               mat[k][k-1] = -mat[k][k-1];
            }
            p += s;
            x = p / s;
            y = q / s;
            z = r / s;
            q /= p;
            r /= p;
                /* row modification */
            for (j = k; j <= (!job ? en : n-1); j++){
               p = mat[k][j] + q * mat[k+1][j];
               if (k != en - 1) {
                  p += r * mat[k+2][j];
                  mat[k+2][j] -= p * z;
               }
               mat[k+1][j] -= p * y;
               mat[k][j] -= p * x;
            }
            j = std::min(en,k+3);
              /* column modification */
            for (i = (!job ? l : 0); i <= j; i++) {
               p = x * mat[i][k] + y * mat[i][k+1];
               if (k != en - 1) {
                  p += z * mat[i][k+2];
                  mat[i][k+2] -= p*r;
               }
               mat[i][k+1] -= p*q;
               mat[i][k] -= p;
            }
            if (job) {             /* accumulate transformations */
               for (i = low; i <= hi; i++) {
                  p = x * vr[i][k] + y * vr[i][k+1];
                  if (k != en - 1) {
                     p += z * vr[i][k+2];
                     vr[i][k+2] -= p*r;
                  }
                  vr[i][k+1] -= p*q;
                  vr[i][k] -= p;
               }
            }
         }
      }
   }

   if (!job) return(0);
   if (norm != 0) {
       /* back substitute to find vectors of upper triangular form */
      for (en = n-1; en >= 0; en--) {
         p = valr[en];
         if ((q = vali[en]) < 0) {            /* complex vector */
            m = en - 1;
            if (fabs(mat[en][en-1]) > fabs(mat[en-1][en])) {
               mat[en-1][en-1] = q / mat[en][en-1];
               mat[en-1][en] = (p - mat[en][en]) /
                     mat[en][en-1];
            }
            else {
               v = std::polar(0.0,-mat[en-1][en])/std::polar(mat[en-1][en-1]-p,q);
               mat[en-1][en-1] = v.real();
               mat[en-1][en] = v.imag();
            }
            mat[en][en-1] = 0;
            mat[en][en] = 1;
            for (i = en - 2; i >= 0; i--) {
               w = mat[i][i] - p;
               ra = 0;
               sa = mat[i][en];
               for (j = m; j < en; j++) {
                  ra += mat[i][j] * mat[j][en-1];
                  sa += mat[i][j] * mat[j][en];
               }
               if (vali[i] < 0) {
                  z = w;
                  r = ra;
                  s = sa;
               }
               else {
                  m = i;
                  if (vali[i] == 0) {
                     v = std::polar(-ra,-sa)/std::polar(w,q);
                     mat[i][en-1] = v.real();
                     mat[i][en] = v.imag();
                  }
                  else {                      /* solve complex equations */
                     x = mat[i][i+1];
                     y = mat[i+1][i];
                     v = std::complex<double>((valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q,
                                               (valr[i] - p)*2*q );

                     if ((fabs(v.real()) + fabs(v.imag())) == 0) {
                        v = std::complex<double>(eps * norm * (fabs(w) +
                                fabs(q) + fabs(x) + fabs(y) + fabs(z)), v.imag());
                     }
                     v = std::polar(x*r-z*ra+q*sa,x*s-z*sa-q*ra) / v;
                     mat[i][en-1] = v.real();
                     mat[i][en] = v.imag();
                     if (fabs(x) > fabs(z) + fabs(q)) {
                        mat[i+1][en-1] =
                             (-ra - w * mat[i][en-1] +
                             q * mat[i][en]) / x;
                        mat[i+1][en] = (-sa - w * mat[i][en] -
                             q * mat[i][en-1]) / x;
                     }
                     else {
                        v = std::polar(-r-y*mat[i][en-1],
                             -s-y*mat[i][en]) / std::polar(z,q);
                        mat[i+1][en-1] = v.real();
                        mat[i+1][en] = v.imag();
                     }
                  }
               }
            }
         }
         else if (q == 0) {                             /* real vector */
            m = en;
            mat[en][en] = 1;
            for (i = en - 1; i >= 0; i--) {
               w = mat[i][i] - p;
               r = mat[i][en];
               for (j = m; j < en; j++) {
                  r += mat[i][j] * mat[j][en];
               }
               if (vali[i] < 0) {
                  z = w;
                  s = r;
               }
               else {
                  m = i;
                  if (vali[i] == 0) {
                     if ((t = w) == 0) t = eps * norm;
                     mat[i][en] = -r / t;
                  }
                  else {            /* solve real equations */
                     x = mat[i][i+1];
                     y = mat[i+1][i];
                     q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
                     t = (x * s - z * r) / q;
                     mat[i][en] = t;
                     if (fabs(x) <= fabs(z)) {
                        mat[i+1][en] = (-s - y * t) / z;
                     }
                     else {
                        mat[i+1][en] = (-r - w * t) / x;
                     }
                  }
               }
            }
         }
      }
             /* vectors of isolated roots */
      for (i = 0; i < n; i++) {
         if (i < low || i > hi) {
            for (j = i; j < n; j++) {
               vr[i][j] = mat[i][j];
            }
         }
      }
       /* multiply by transformation matrix */

      for (j = n-1; j >= low; j--) {
         m = std::min(j,hi);
         for (i = low; i <= hi; i++) {
            for (z = 0,k = low; k <= m; k++) {
               z += vr[i][k] * mat[k][j];
            }
            vr[i][j] = z;
         }
      }
   }
    /* rearrange complex eigenvectors */
   for (j = 0; j < n; j++) {
      if (vali[j] != 0) {
         for (i = 0; i < n; i++) {
            vi[i][j] = vr[i][j+1];
            vr[i][j+1] = vr[i][j];
            vi[i][j+1] = -vi[i][j];
         }
         j++;
      }
   }
   return(0);
}

int nonsymm_eigen(const int job,double **A, const int n,double *rr, double *ri,
                 double **vr, double **vi, const int maxiter)
{    
   //////////////////////////////////////////////////////////////////////////
   //  This nonsym_eigen() works for eigenvalue/vector analysis
   //         for real general square matrix A
   //         A will be destroyed
   //         rr,ri are vectors containing eigenvalues
   //         vr,vi are  matrice containing eigenvectors
   //         maxiter is  max. no. of iterations to converge, usually 30
   //         note that A*v(:,j) = v(:,j)*r(j) where v(:,j) means jth column
   //   job = 0, no eigenvectors are needed
   //         1, eigenvectors are also needed
   //
   //  Algorithm: Handbook for Automatic Computation, vol 2
   //             by Wilkinson and Reinsch, 1971
   //             most of source codes were taken from a public domain
   //             solftware called MATCALC.
   //  Credits:   goes to the authors of MATCALC
   //
   //  return     -2 out of memory
   //             -1 not converged
   //              0 real    eigenvalues/vectors
   //              1 complex eigenvalues/vectors
   //
   //  Tianlin Wang at University of Illinois
   //  Thu May  6 15:22:31 CDT 1993
   //////////////////////////////////////////////////////////////////////////

   double base	= 2.0;           // base of floating point arithmetic
   double digits = 53.0;          // no. of digits to the base in the fraction
   double eps  = pow(base,1.0-digits);
   double tiny =  sqrt(eps);
   int low,hi,i,k=0;
   double *w;
   if(n>0){
     w = new double [n];
   }
   else {
     w = 0;
   }
   if (!w) return -2;
   balance(A,n,&low,&hi,w,base);
   elemhess(job,A,n,low,hi,vr,vi);
   if (-1 == nonsymm_eigen_core(job,A,n,low,hi,rr,ri,vr,vi,eps,maxiter)) {
      k = -1;
   }
   else {
      if (job) unbalance(n,vr,vi,low,hi,w);
      for (i=0; i<n; i++) if (fabs(ri[i]) > tiny) {k = 1; break;}
   }
   if(w){
     delete [] w;
     w=0;
   }
   return k;
}
}  //////// end of namespace matvec

