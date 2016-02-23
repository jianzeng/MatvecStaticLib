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

#include "model.h"
#include "data.h"


namespace matvec {
static void linmin(double *p,double *xi,int n,double& fret);
static double f1dim(double x);
static double brent(double ax,double bx,double cx,double (*f)(double),
                    double tol, double& xmin);
static void mnbrak(double& ax,double& bx,double& cx,double& fa,
                   double& fb,double& fc, double (*func)(double));

static Data *myD;
static Model *myM;
static int nfunk;

static double funk(Vector<double> &x)
{
   int i,k;
   unsigned nterms = myM->nterm();
   myM->vec2var(x);
   if (!myM->residual_var.psd()) return 1.0e30;
   for (i=0; i<nterms; i++) {
      if (myM->term[i].classi() == 'R' || myM->term[i].classi() == 'P') {
         k = myM->term[i].prior->var_matrix()->psd();
         if (!k) return (1.0e30 + 1.0e10*(i+1));
      }
   }
   nfunk++;
   double log_likelihood = myM->restricted_log_likelihood();
   if ( (nfunk % 3) == 0) {
      std::cout << " powell...# of log_likelihood evaluations = " << nfunk << ", "
           << "and its current value = " << log_likelihood << std::endl;
      myM->info(std::cout);
   }
   return log_likelihood*(-1.0);
}

static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

double nr_powell(Data* D,Model* M,Vector<double> &p,Matrix<double> &xi,int n,int& iter,
                double ftol)
{
   myD = D;
   myM = M;
   nfunk = 0;

   int i,ibig,j;
   double t,fptt,fp,del,fret;
   int maxiter = iter;

   Vector<double> pt(n);
   Vector<double> ptt(n);
   Vector<double> xit(n);

   fret = funk(p);
   for (j=0; j<n; j++) pt[j]=p[j];
   for (iter=1; ;iter++) {
      fp = fret;
      ibig = 0;
      del = 0.0;
      for (i=0; i<n; i++) {
         for (j=0; j<n; j++) xit[j] = xi[j][i];
         fptt = fret;
         linmin(p.begin(),xit.begin(),n,fret);
         if (fabs(fptt-fret) > del) {
            del = fabs(fptt-fret);
            ibig = i;
         }
      }
      if (2.0*fabs(fp-fret) <= ftol*(fabs(fp)+fabs(fret))) {
         iter = nfunk;
         return fret;
      }
      if (iter == maxiter) {
         warning(" Too many iterations in routine POWELL\n"
               "%d out of %d",iter,maxiter);
         return fret;
      }
      for (j=0; j<n; j++) {
         ptt[j] = 2.0*p[j] - pt[j];
         xit[j] = p[j] - pt[j];
         pt[j] = p[j];
      }
      fptt = funk(ptt);
      if (fptt < fp) {
         t = 2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
         if (t < 0.0) {
            linmin(p.begin(),xit.begin(),n,fret);
            for (j=0; j<n; j++) xi[j][ibig] = xit[j];
         }
      }
   }
}
#undef SQR

#define TOL 2.0e-8
int ncom=0;   /* defining declarations */
double *pcom = 0;
double *xicom = 0;

static void linmin(double *p,double *xi,int n,double& fret)
{
   int j;
   double xx,xmin,fx,fb,fa,bx,ax;

   ncom = n;
   
   if(n>0){
     pcom = new double [n];
   }
   else {
     pcom = 0;
   }
   if(n>0){
     xicom = new double [n];
   }
   else {
     xicom = 0;
   }
   for (j=0; j<n; j++) {
      pcom[j] = p[j];
      xicom[j] = xi[j];
   }
   ax=0.0;
   xx=1.0;
   bx=0.0; fa=0.0; fx=0.0; fb=0.0; xmin=0.0;
   mnbrak(ax,xx,bx,fa,fx,fb,f1dim);
   fret = brent(ax,xx,bx,f1dim,TOL,xmin);
   for (j=0; j<n; j++) {
      xi[j] *= xmin;
      p[j] += xi[j];
   }
   if(xicom){
     delete [] xicom;
     xicom=0;
   }
   if(pcom){
     delete [] pcom;
     pcom=0;
   }
}
#undef TOL

static double f1dim(double x)
{
   Vector<double> xt(ncom);
   for (int j=0; j<ncom; j++) xt[j] = pcom[j] + x*xicom[j];
   double f = funk(xt);
   return f;
}

#define GOLD 1.618034
#define GLIMIT 200.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static void mnbrak(double& ax,double& bx,double& cx,double& fa,
        double& fb,double& fc,double (*func)(double))
{
   double ulim,u,r,q,fu,dum;

   fa=(*func)(ax);
   fb=(*func)(bx);
   if (fb > fa) {
      SHFT(dum,ax,bx,dum)
      SHFT(dum,fb,fa,dum)
   }
   cx = bx+GOLD*(bx-ax);
   fc = (*func)(cx);
   while (fb > fc) {
      r=(bx-ax)*(fb-fc);
      q=(bx-cx)*(fb-fa);
      u=(bx)-((bx-cx)*q-(bx-ax)*r)/
         (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
      ulim = bx+GLIMIT*(cx-bx);
      if ((bx-u)*(u-cx) > 0.0) {
         fu=(*func)(u);
         if (fu < fc) {
            ax = bx;
            bx = u;
            fa = fb;
            fb = fu;
            return;
         } else if (fu > fb) {
            cx = u;
            fc = fu;
            return;
         }
         u= cx +GOLD*(cx-bx);
         fu=(*func)(u);
      } else if ((cx-u)*(u-ulim) > 0.0) {
         fu=(*func)(u);
         if (fu < fc) {
            SHFT(bx,cx,u,cx+GOLD*(cx-bx))
            SHFT(fb,fc,fu,(*func)(u))
         }
      } else if ((u-ulim)*(ulim-cx) >= 0.0) {
         u=ulim;
         fu=(*func)(u);
      } else {
         u=(cx)+GOLD*(cx-bx);
         fu=(*func)(u);
      }
      SHFT(ax,bx,cx,u)
      SHFT(fa,fb,fc,fu)
   }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT


#define ITMAX 200
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static double brent(double ax,double bx,double cx,double (*f)(double),double tol,double& xmin)
{
   int iter;
   double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;
   double d=0.0;

   a=((ax < cx) ? ax : cx);
   b=((ax > cx) ? ax : cx);
   x=w=v=bx;
   fw=fv=fx=(*f)(x);
   for (iter=1;iter<=ITMAX;iter++) {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
         xmin=x;
         return fx;
      }
      if (fabs(e) > tol1) {
         r=(x-w)*(fx-fv);
         q=(x-v)*(fx-fw);
         p=(x-v)*q-(x-w)*r;
         q=2.0*(q-r);
         if (q > 0.0) p = -p;
         q=fabs(q);
         etemp=e;
         e=d;
         if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
         else {
            d=p/q;
            u=x+d;
            if (u-a < tol2 || b-u < tol2)
               d=SIGN(tol1,xm-x);
         }
      } else {
         d=CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu=(*f)(u);
      if (fu <= fx) {
         if (u >= x) {
            a=x;
         } else {
            b=x;
         }
         SHFT(v,w,x,u)
         SHFT(fv,fw,fx,fu)
      } else {
         if (u < x) {
            a=u;
         } else {
            b=u;
         }
         if (fu <= fw || w == x) {
            v=w;
            w=u;
            fv=fw;
            fw=fu;
         } else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
         }
      }
   }
   warning(" Too many iterations in BRENT");
   xmin=x;
   return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN
} ///////// end of namespace matvec

