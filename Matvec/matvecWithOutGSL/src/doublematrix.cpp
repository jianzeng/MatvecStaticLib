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
#include "session.h"
#include "util.h"
#include "doublematrix.h"


namespace matvec {

extern int nonsymm_eigen(const int job,double **A, const int n,double *rr,
                double *ri, double **vr, double **vi,const int maxiter=50);
extern int symm_eigen(double *vals, double **A, const int n);
extern void svdcmp(double* a[], int m, int n, double w[], double* v[]);
extern void ludcmp(double **a,int n,int *indx,int& d, double tol);
extern void lubksb(double **a,int n,int *indx,double *b);
extern void sweep(const int m, const int n, double **a, const int i0,
                         const int i1,const double tol);
extern bool psdefinite(const double **mat,const unsigned m,const unsigned n,const double tol);

doubleMatrix&  doubleMatrix::operator = (const Matrix<bool> &a)
{
   resize(a.num_rows(),a.num_cols());
   for (int j,i=0; i<a.num_rows(); ++i) for (j=0; j<a.num_cols(); ++j) Matrix<double>::me[i][j] = (double)a[i][j];
   return *this;
}


bool doubleMatrix::psd(void) const
{
   if (!this->symmetric()) return false;
   return psdefinite(begin(),Matrix<double>::nrow,Matrix<double>::ncol,SESSION.epsilon);
}

void doubleMatrix::svd(doubleMatrix& u, Vector<double>& s, doubleMatrix& v)
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   if (nrow < ncol) throw exception("doubleMatrix::svd(), bad args");

   if (u.nrow != nrow || u.ncol != ncol) u.resize(nrow,ncol);
   if (s.size() != ncol ) s.resize(ncol);
   if (v.nrow != ncol || v.ncol != ncol) v.resize(ncol,ncol);

   for (int i=0; i<nrow; i++) memcpy(u[i],me[i],sizeof(double)*ncol);
   svdcmp(u.me, nrow, ncol, s.begin(), v.me);      // Singular value decomposition.
}

int doubleMatrix::rank(void)                      // rank(a) based on svdcmp.cc
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   Vector<double> w(ncol);
   Matrix<double> v(ncol,ncol);
   int m = nrow;
   if (nrow < ncol) m = ncol;                   //u(m,col) where m>=col
   Matrix<double> u(m,ncol);
   int i;
   for (i=0; i<nrow; i++) memcpy(u[i],me[i],sizeof(double)*ncol);
   svdcmp(u.begin(), m, ncol, w.begin(), v.begin());                 // Singular value decomposition.

   int  r = 0;
   for (i=0; i<ncol; i++) if (w[i] > SESSION.epsilon) r++;
   return (r);
}

// Build a design Matrix for cubic splines
// Based on Green and Silverman (1994);Nonparametric Regression and 
//          Generalized Linear Models, Chapman & Hall, Chap. 2
//
doubleMatrix doubleMatrix::splines(const doubleMatrix knots,const unsigned type )
const
{
  int  nk=knots.nrow;
  int  n=nrow;
  double k;
  doubleMatrix R,Q,Delt,H,X;
  R.resize(nk,nk);
  Q.resize(nk,nk);
  H.resize(nk-1,1);
  k=(knots[nk-1][0]-knots[0][0])/((double)(nk-1));
  for(int i=0;i<nk-1;i++) H[i][0]=(knots[i+1][0]-knots[i][0]);
  for(int i=2;i<nk;i++){
    R(i,i)=(H(i,1)+H(i-1,1))/3;
    if (i<(nk-1)) {
      R(i+1,i)=H(i,1)/6;
      R(i,i+1)=H(i,1)/6;
    }
    Q(i-1,i)=1/H(i-1,1);
    Q(i,i)=-(1/H(i-1,1))-(1/H(i,1));
    Q(i+1,i)=1/H(i,1);
  }
  doubleMatrix Xg,Xgamma;
  Xg.resize(n,nk);
  Xgamma.resize(n,nk);
  for(int i=1;i<=n;i++){
    double t=me[i-1][0];
    int j;
    for( j=1;t>knots.me[j][0];j++);
    double th=knots.me[j][0];
    double tl=knots.me[j-1][0];
    Xg(i,j+1)=(t-tl)/H(j,1);
    Xg(i,j)=(th-t)/H(j,1);
    Xgamma(i,j+1)=-(t-tl)*(th-t)*(1+(t-tl)/H(j,1))/6;
    Xgamma(i,j)=-(t-tl)*(th-t)*(1+(th-t)/H(j,1))/6;
  }
  X=Xg+Xgamma*R.ginv1()*Q.transpose();
  if(type==1) {
    doubleMatrix Xs,Delt,G;
    G.resize(nk-2,nk-2);
    Xs.resize(nk,nk);
    Delt.resize(nk,nk-2);
    for(int i=0;i<nk-2;i++){
      Delt[i][i]=1./H[i][0];
      Delt[i+1][i]=-(1./H[i][0]+1./H[i+1][0]);
      Delt[i+2][i]=1./H[i+1][0];
      G[i][i]=(H[i][0]+H[i+1][0])/3.;
      if(i) {
	G[i-1][i]=H[i][0]/6.;
	G[i][i-1]=H[i][0]/6.;
      }
    }
    G.sqrtm();
    doubleMatrix Zs,dpd;
    dpd=Delt.transpose()*Delt;
    Zs=Delt*dpd.inv()*G;
    Zs/=pow(k,1.5);
    for(int i=0;i<nk;i++){
      Xs[i][0]=1.;
      Xs[i][1]=knots[i][0];
      for(int j=0;j<nk-2;j++){
	Xs[i][j+2]=Zs[i][j];
      }

    }
    X=X*Xs;
  }
  //cout << X;

  return(X);
}


doubleMatrix& doubleMatrix::inv(void)
{
   //////////////////////////////////////////////////////////
   // A itself turns into its inverse
   // returns non_singular indicator
   // REQUIREMENTS:  A must be nonsingular
   /////////////////////////////////////////////////////////
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   int i,j,d;
   if ( nrow != ncol ) throw exception("doubleMatrix::inv(): matrix must be square");
   Vector<int>    indx;  indx.reserve(nrow);
   Vector<double> tempcol; tempcol.reserve(nrow);
   Matrix<double> temp; temp.reserve(nrow,nrow);

   ludcmp(me,nrow,indx.begin(),d, SESSION.epsilon);
   for (j=0; j<nrow; j++) {
      tempcol.assign(0.0);
      tempcol[j] = 1.0;
      lubksb(me,nrow,indx.begin(),tempcol.begin());
      for (i=0; i<nrow;i++) temp[i][j] = tempcol[i];
   }
   for (i=0;i<nrow; i++) memcpy(me[i],temp[i],sizeof(double)*nrow);
   return *this;
}




doubleMatrix doubleMatrix::mat_log(double tol) const
{                              // Matrix Log
   if(tol==0) tol=SESSION.epsilon;
   doubleMatrix temp,P;
   if (!this->symmetric()) throw exception("doubleMatrix.mat_exp() must be symmetric\n");
   if (nrow==0) return temp;
   temp.resize(nrow,nrow);
   Vector<double> D;

   P=*this;
   D=P.eigen();
   double dmin=D.min();
   double dmax=D.max();
   if(dmin < dmax*2.*tol) {
     dmin=dmax*2.*tol;
   }
   for(int k=0;k<nrow;k++) {
     if(D[k]<dmin) D[k]=dmin;
     D[k]=std::log(D[k]);
     for(int i=0;i<nrow;i++){
       for(int j=0;j<nrow;j++){
	 temp[i][j]+=P[i][k]*P[j][k]*D[k];
       }
     }
   }
   return temp;
}



doubleMatrix doubleMatrix::mat_exp(double tol) const
{                              // Matrix Exponential
   if(tol==0) tol=SESSION.epsilon;

   doubleMatrix temp,P;
   if (!this->symmetric()) throw exception("doubleMatrix.mat_exp() must be symmetric\n");
   if (nrow==0) return temp;
   temp.resize(nrow,nrow);
   Vector<double> D;

   P=*this;
   D=P.eigen();
   double dmin=D.min();
   double dmax=D.max();
   if(dmin < dmax+2.*std::log(tol)) {
     dmin=dmax+2.*std::log(tol);
   }
   for(int k=0;k<nrow;k++) {
     if(D[k]<dmin) D[k]=dmin;
     D[k]=std::exp(D[k]);
     for(int i=0;i<nrow;i++){
       for(int j=0;j<nrow;j++){
	 temp[i][j]+=P[i][k]*P[j][k]*D[k];
       }
     }
   }
   return temp;
}




void doubleMatrix::mat_exp_der(Vector<doubleMatrix> &temp,double tol) const
{                              // Matrix exp partial derivatives
   if(tol==0) tol=SESSION.epsilon;

   doubleMatrix P;
   if (!this->symmetric()) throw exception("doubleMatrix.mat_exp() must be symmetric\n");
   if (nrow==0) return;
   int nvc=(nrow*(nrow+1))/2;
   temp.resize(nvc);
   for(int k=0;k<nvc;k++) temp[k].resize(nrow,nrow);

   Vector<double> D;

   P=*this;
   D=P.eigen();
   //      cout << "\nP der\n" << P;
   double dmin=D.min();
   double dmax=D.max();
   if(dmin < dmax*2.*tol) {
     dmin=dmax*2.*tol;
   }
   for(int i=0;i<nrow;i++) {
     if(D[i]<dmin) D[i]=dmin;
   }
   doubleMatrix Partkl;
   double g;
   for(int k=0;k<nrow;k++){
     double dk=D[k];
     for(int l=0;l<nrow;l++) {
       double dl=D[l];
       if(std::fabs(dk-dl)< tol) g=std::sqrt(dk*dl);
       else g=(dk-dl)/std::log(dk/dl);
       Partkl.resize(nrow,nrow);

       for(int i=0;i<nrow;i++)
	 for(int j=0;j<nrow;j++)
	   Partkl[i][j]=g*P[i][k]*P[j][l];
       
       
       int vc=0;
       for(int i=0;i<nrow;i++) {
	 for(int j=i;j<nrow;j++,vc++) {
	   temp[vc]+=P[i][k]*P[j][l]*Partkl;
	   if(i!=j) temp[vc]+=P[j][k]*P[i][l]*Partkl;
	   //	    cout << i << ":" << j << " "<< vc << " " << temp[vc] << endl;
	 }
       }

     }
   }
   //      for(int i=0;i<nvc;i++) cout << i << "->" << temp[i] << endl;

   return;
}






doubleMatrix doubleMatrix::ginv0(void) const
{                              // Penrose-Moore inverse of any matrix
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   doubleMatrix temp(ncol,nrow);
   if (ncol==0 || nrow==0) return temp;
   Vector<double>  w; w.reserve(ncol);
   Matrix<double> v; v.reserve(ncol,ncol);
   int m = nrow;
   if (nrow < ncol) m=ncol;                   //u(m,col) where m>=col
   Matrix<double> u; u.reserve(m,ncol);

   int i,j,k;
   double x;
   for (i=0; i<nrow; i++) memcpy(u[i],me[i],sizeof(double)*ncol);
   for (i=0; i<m-nrow; i++) memset(u[nrow+i],'\0',sizeof(double)*ncol);
   svdcmp(u.begin(), m, ncol, w.begin(), v.begin());                  // Singular value decomposition.
   double wmax = w.max();             // Maximum singular value.
   double wmin = wmax*SESSION.epsilon;
   for (k = 0; k < ncol; k++) {
      if (w[k] < wmin) {
         w[k] = 0.0;
      }
      else {
         w[k] = 1.0/w[k];
      }
   }
   for (i = 0; i < ncol; i++) {
      for (k = 0; k < ncol; k++) v[i][k] *= w[k];
      for (j = 0; j < nrow; j++) {
         for (x=0.0,k=0; k<ncol; ++k) x += v[i][k]*u[j][k];
         temp[i][j] = x;
      }
   }
   return temp;
}

doubleMatrix& doubleMatrix::ginv1(unsigned *irank)
{         // it only works for symmetric-positive-semidefinite matrix
          //  Matrix itself becomes its inverse;
   if (!this->symmetric()) {
     throw exception("doubleMatrix::ginv1(): matrix must be symmetric");
   }
   double lgdet;
   unsigned myrank = ginverse1(begin(),num_rows(),lgdet,1,SESSION.epsilon);
   if (irank) *irank = myrank;
   return *this;
}

//doubleMatrix& doubleMatrix::MPGinv(){
//    // Authors: Rohan L. Fernando 
//    // (June, 2006) 
//    // Contributors: 
//
//   bool trnsps = false;
//   int nrow =  Matrix<double>::nrow;
//   int ncol =  Matrix<double>::ncol;
//   double** me =  Matrix<double>::me;
//   if (nrow==1 and ncol==1){
//		me[0][0] = 1.0/me[0][0];
//		return *this;
//   }
//   if (nrow<ncol){
//	  *this = transpose();
//	  nrow = Matrix<double>::nrow;
//	  ncol = Matrix<double>::ncol;
//	   trnsps = true;
//   }
//
//      
//	gsl_matrix *V = gsl_matrix_alloc(ncol, ncol);
//	gsl_vector *S = gsl_vector_alloc(ncol); 
//	gsl_vector *work = gsl_vector_alloc(ncol);
//	gsl_matrix *A = gsl_matrix_alloc(nrow, ncol);
//	for (unsigned i=0;i<nrow;i++){
//		for (unsigned j=0;j<ncol;j++){
//			gsl_matrix_set(A,i,j,me[i][j]);
//		}
//	}
//	
//	gsl_linalg_SV_decomp(A, V, S, work);
//	// gives A = UDV'; U replaces A; 
//	// S contains the diagonals of D
//	
//	
//	double logdet=0.0;
//	for(unsigned i=0; i<ncol; i++){
//		double eigenVal = gsl_vector_get(S, i);
//		if(eigenVal > 1.0e-10){
//			logdet += logf(eigenVal);
//		}
//	}
//
//	resize(ncol,nrow);
//	double HiVal;
//	for(unsigned i=0; i<ncol; i++){
//		for(unsigned j=0; j<nrow; j++){
//			HiVal = 0.0;
//			for(unsigned k=0; k<ncol; k++){
//				double eigenVal  = gsl_vector_get(S, k);
//				double uVal      = gsl_matrix_get(A, j, k);
//				double vtransVal = gsl_matrix_get(V, i, k);
//				if(eigenVal > 1.0e-10){
//					HiVal += uVal * vtransVal / eigenVal;
//
//				}
//			}
//			me[i][j] = HiVal;
//		}
//	}
//	if (trnsps) *this = transpose();
//	return *this;
//}
//
//double doubleMatrix::MPLogDet(){
//    // Authors: Rohan L. Fernando 
//    // (October, 2006) 
//    // Contributors: 
//
//   // returns the log determinant, the matrix becomes the g-inverse
//   bool trnsps = false;
//   int nrow =  Matrix<double>::nrow;
//   int ncol =  Matrix<double>::ncol;
//   double** me =  Matrix<double>::me;
//   if (nrow==1 and ncol==1){
//        double logDet = logf(me[0][0]);
//		me[0][0] = 1.0/me[0][0];
//		return logDet;
//   }
//   if (nrow<ncol){
//	  *this = transpose();
//	  nrow = Matrix<double>::nrow;
//	  ncol = Matrix<double>::ncol;
//	   trnsps = true;
//   }
//
//      
//	gsl_matrix *V     = gsl_matrix_alloc(ncol, ncol);
//	gsl_vector *S     = gsl_vector_alloc(ncol); 
//	gsl_vector *work  = gsl_vector_alloc(ncol);
//	gsl_matrix *A     = gsl_matrix_alloc(nrow, ncol);
//	for (unsigned i=0;i<nrow;i++){
//		for (unsigned j=0;j<ncol;j++){
//			gsl_matrix_set(A,i,j,me[i][j]);
//		}
//	}
//	
//	gsl_linalg_SV_decomp(A, V, S, work);
//	// gives A = UDV'; U replaces A; 
//	// S contains the diagonals of D
//	
//	
//	double logdet=0.0;
//	for(unsigned i=0; i<ncol; i++){
//		double eigenVal = gsl_vector_get(S, i);
//		if(eigenVal > 1.0e-10){
//			logdet += logf(eigenVal);
//		}
//	}
//
//	resize(ncol,nrow);
//	double HiVal;
//	for(unsigned i=0; i<ncol; i++){
//		for(unsigned j=0; j<nrow; j++){
//			HiVal = 0.0;
//			for(unsigned k=0; k<ncol; k++){
//				double eigenVal  = gsl_vector_get(S, k);
//				double uVal      = gsl_matrix_get(A, j, k);
//				double vtransVal = gsl_matrix_get(V, i, k);
//				if(eigenVal > 1.0e-10){
//					HiVal += uVal * vtransVal / eigenVal;
//
//				}
//			}
//			me[i][j] = HiVal;
//		}
//	}
//	if (trnsps) *this = transpose();
//	gsl_vector_free(S);
//	gsl_vector_free(work);
//	gsl_matrix_free(V);
//	gsl_matrix_free(A);
//	
//	return logdet;
//}
//
//


double doubleMatrix::cond(void)
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   if (nrow < ncol) throw exception("doubleMatrix::cond(), M.nrow() >= M.ncol()");

   Vector<double>  w(ncol);
   Matrix<double> v(ncol,ncol);
   Matrix<double> u(nrow,ncol);

   for (int i=0; i<nrow; i++) memcpy(u[i],me[i],sizeof(double)*ncol);
   svdcmp(u.begin(), nrow, ncol, w.begin(), v.begin());         // Singular value decomposition.

   double smax = w.max();
   double smin = w.min();
   if (smin == 0.0) warning("Matrix::cond(): condition is infinite");
   return (smax/smin);
}

doubleMatrix doubleMatrix::lu_solve(const doubleMatrix& rhs)
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   if (nrow  != ncol) throw exception("doubleMatrix::lu_solve(): matrix is not square ");
   int m = rhs.num_rows();
   int n = rhs.num_cols();

   if (ncol != m ) throw exception("doubleMatrix::lu_solve(): bad arg");
   doubleMatrix solmat(m,n);
   int i,j,d;
   Vector<double> solvec; solvec.reserve(m);
   Vector<int> indx;  indx.reserve(m);
   ludcmp(me,nrow,indx.begin(),d,SESSION.epsilon);
   warning("A.solve(b), A destroyed");
   for (j=0; j<n; j++) {
      for (i=0; i<m; i++) solvec[i] = rhs[i][j];
      lubksb(me,m,indx.begin(),solvec.begin());
      for (i=0; i<m; i++) solmat[i][j] = solvec[i];
   }
   return solmat;
}

Vector<double> doubleMatrix::lu_solve(const Vector<double>& rhs)
{
   int n = rhs.size();
   doubleMatrix bmat(n,1);
   for (int i=0; i<n; i++) bmat[i][0] = rhs[i];
   bmat = lu_solve(bmat);
   return bmat.vec();
}


void doubleMatrix::gs_solve(const doubleMatrix& rhs,doubleMatrix& solmat,const double relax,
                     const double stopval,const int mxiter)
{
  //   relax    relaxation coefficient, relax=1.2 seems good for some case
  //   stopval  the accuracy at which iteration stops.
  //   mxiter   maximum number of iterations allowed

   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   int m = rhs.num_rows();
   int n = rhs.num_cols();
   if (nrow != ncol) throw exception("doubleMatrix::gs_solve(): matrix must be square ");
   if (ncol != m ) throw exception("doubleMatrix::gs_solve(): matrix and rhs are not conformable");
   if (solmat.num_rows() != m || solmat.num_cols() != n) {
      solmat.resize(m,n);
   }
   int i,j,k;
   double diag;
   int niter = 0;
   double a,cmax;
   double tol = SESSION.epsilon;

   Vector<double> oldsol(n);
   Vector<double> newsol(n);
   Vector<double> cval(n);
   Vector<double> local(n);
   do {                             // now iteration begins
      cmax = 0.0;
      for (i=0; i<m; i++) {
         for (k=0; k<n; k++) {
            oldsol[k] = solmat[i][k];
            local[k] = 0.0;
         }
         for (j=0; j<m; j++) {
            a = me[i][j];
            if (j != i) for (k=0; k<n; k++) local[k] += a * solmat[j][k];
         }
         diag = me[i][i];
         if (diag > tol) {
            for (k=0; k<n; k++) {
               newsol[k] = (rhs[i][k] - local[k])/diag;
               cval[k] = (newsol[k] - oldsol[k])*relax;
               if (fabs(cval[k]) > cmax) cmax = fabs(cval[k]);
               newsol[k] = solmat[i][k] = oldsol[k] + cval[k];
            }
         }
         else if (diag > -tol) {
            throw exception("doubleMatrix::gs_solve(): zero-diagonal found in Matrix::gs_solve()");
         }
         else {
            throw exception("doubleMatrix::gs_solve(): matrix is not psd");
         }
      }
      niter += 1;
      if ( (niter % 10) ==0) {
         std::cout << " GS: # of iter = " << niter << ", max_change = "
              << cmax << std::endl;
      }
   } while (niter <= mxiter && cmax > stopval);
   if (cmax > stopval) {
      warning("doubleMatrix::gs_solve(): not converge: %20.10f",cmax);
   }
}

void doubleMatrix::gs_solve(const Vector<double>& rhs, Vector<double>& sol, const double relax,
                     const double stopval, const int mxiter)
{
   int n = rhs.size();
   doubleMatrix bmat(n,1),solmat(n,1);
   if (sol.size() != n) {
      sol.resize(n);
   }
   int i;
   for (i=0; i<n; i++) {
      bmat[i][0] = rhs[i];
      solmat[i][0] = sol[i];
   }
   gs_solve(bmat,solmat,relax,stopval,mxiter);
   for (i=0; i<n; i++) sol[i] = solmat[i][0];
}

Vector<double> doubleMatrix::eigen(const int job)
{
   ///////////////////////////////////////////////////////////////////
   // job = 0, no eigenvectors are needed
   //       1, eigenvectors are also needed.
   // return eigenvalues
   // Note that  matrix itself become eigenvectors.
   //   all imaginary numbers in eigenvalus/eigenvectors are ignored
   //   if there are some for non-symmetrix real matrix.
   ////////////////////////////////////////////////////////////////////

  //   int nrow =  Matrix<double>::nrow;
  //   int ncol =  Matrix<double>::ncol;
  //   double** me =  Matrix<double>::me;

   if (nrow != ncol) throw exception( "doubleMatrix::eigen(): matrix must be square" );
   Vector<double> vals(nrow);
   if (this->symmetric()) {
      symm_eigen(vals.begin(),me,nrow);
   }
   else {
     std::cout << *this;
      Vector<double> ri(nrow); 
      Matrix<double> vr;
      Matrix<double> vi;
      if (job == 1) {
          vr.reserve(nrow,nrow);
          vi.reserve(nrow,nrow);
      }
      int i = nonsymm_eigen(job,me,nrow,vals.begin(),ri.begin(),vr.begin(),vi.begin());
      if (i == -2) {
         warning(" Matrix::eigen(): out of memory");
      }
      else if (i == -1) {
         warning(" Matrix::eigen(): not converged");
      }
      else if (i == 1) {
         warning(" Matrix::eigen(): eigenvalues/vectors are complex,\n "
                 " and the imaginary numbers are all ignored");
      }
      if (job == 1) {
         for (i=0; i<nrow; i++) memcpy(me[i],vr[i],sizeof(double)*ncol);
      }
   }
   return vals;
}

double doubleMatrix::det(void) const
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;

   double detval = 0.0;
   if ( nrow != ncol ) throw exception("doubleMatrix::det(): matrix must be square");
   int i,d;
   Vector<int> indx(nrow);
   doubleMatrix temp(*this);
   ludcmp(temp.begin(),nrow,indx.begin(),d, SESSION.epsilon);
   detval = (double)d;
   for (i=0; i<nrow; i++) detval *= temp[i][i];
   return detval;
}

double doubleMatrix::norm(const int p) const
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   double retval = 0.0;
   int i;
   if (p == 1) {
      doubleMatrix A(*this);
      retval = abs(A).sum().max();
   }
   else if (p == 2) {
      Vector<double> w;  w.reserve(ncol);
      Matrix<double> v; v.reserve(ncol,ncol);
      int m = nrow;
      if (nrow < ncol) m = ncol;            //u(m,ccol) where m>=ccol
      Matrix<double> u(m,ncol);

      for (i=0; i<nrow; i++) memcpy(u[i],me[i],sizeof(double)*ncol);
      svdcmp(u.begin(), m, ncol, w.begin(), v.begin());                 // Singular value decomposition.
      retval = w.max();
   }
   else {
      throw exception("doubleMatrix::norm(): bad arg, 1 or 2 expected");
   }
   return retval;
}

double doubleMatrix::norm(const std::string &s) const
{
   double retval = 0.0;
   doubleMatrix A = this->transpose();
   if (s == "inf") {
      retval = abs(A).sum().max();
   }
   else if (s == "fro") {
      retval = std::sqrt((A*(*this)).diag(0).sum());
   }
   else {
      throw exception("doubleMatrix::norm((): bad arg, it must be either \"inf\" or \"fro\" ");
   }
   return retval;
}

double doubleMatrix::logdet(void)
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   double retval = 0.0;
   if (!this->symmetric()) throw exception("doubleMatrix::logdet(): matrix must be symmetric");
   Matrix<double> temp(*this);

   ginverse1(temp.begin(),nrow,retval,3,SESSION.epsilon);
   return retval;
}

double doubleMatrix::quadratic(const Vector<double> &x, const Vector<double> &y) const
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;
   int m = x.size();
   int n = y.size();

   double s = 0.0;
   if (m != nrow  || n != ncol)  throw exception("doubleMatrix::quadratic(): size not conformable");
   for (int j,i=0;i<m; i++) for (j=0; j<n; j++) s += x[i]*me[i][j]*y[j];
   return s;
}
Vector<double> doubleMatrix::variance(orientation orien) const
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;
   int i,j;
   double ss,s;
   Vector<double> temp;
   if (orien == COLUMN) {
      assert(nrow > 1);
      temp.resize(ncol);
      for (j=0; j<ncol; j++) {
         for (ss=0.0,s=0.0,i=0; i<nrow; ++i) {
             ss += me[i][j]*me[i][j];
             s  += me[i][j];
         }
         temp[j] = (ss-s*s/nrow)/(nrow - 1);
      }
   } else if (orien == ROW) {
      assert(ncol > 1);
      temp.resize(nrow);
      for (i=0; i<nrow; i++) {
         for (ss=0.0,s=0.0,j=0; j<ncol; ++j) {
             ss += me[i][j]*me[i][j];
             s  += me[i][j];
         }
         temp[i] = (ss-s*s/ncol)/(ncol - 1);
      }
   } else {
      warning("unknown orientation");
   }
   return temp;
}

doubleMatrix doubleMatrix::covariance(const doubleMatrix &B) const
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;
   double** me =  Matrix<double>::me;

   assert (nrow == B.num_rows() && ncol == B.num_cols());
   doubleMatrix temp(ncol,ncol);
   int i,j,k;
   double ss,sx,sy;
   for (i=0; i<ncol; ++i) {
      for (j=0; j<=i; ++j) {
         for (ss=0.0,sx=0.0,sy=0.0,k=0; k<nrow; ++k) {
            ss += me[k][i]*B[k][j];
            sx += me[k][i];
            sy += B[k][j];
         }
         temp[i][j] = temp[j][i] = (ss - sx*sy/ncol)/(ncol - 1);
      }
   }
   return temp;
}

doubleMatrix& doubleMatrix::sqrtm(void)
{
   if (!this->symmetric()) throw exception("doubleMatrix::sqrtm(): matrix must be symmetric");
   double lgdet;
   ginverse1(begin(),num_rows(),lgdet,2,SESSION.epsilon);
   return *this;
}

doubleMatrix& doubleMatrix::identity(const int m,const int n)
{
   resize(m,n);
   double **me = begin();
   for (int i=0; i<m; ++i) if (i < n) me[i][i] = 1.0;
   return *this;
}

doubleMatrix& doubleMatrix::sweep(const int i0,const int i1)
{
   int nrow =  Matrix<double>::nrow;
   int ncol =  Matrix<double>::ncol;

   int n = std::min(nrow,ncol);
   if (i0<0 || i0 >=n || i0 > i1) throw exception("doubleMatrix::sweep(): range error");
   matvec::sweep(nrow,ncol,begin(),i0,i1,SESSION.epsilon);
   return *this;
}
}  //////////// end of namespace

