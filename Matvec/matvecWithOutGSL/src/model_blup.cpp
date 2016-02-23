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

#include <sstream>
#include "session.h"
#include "model.h"
#include "population.h"
#include "individual.h"

namespace matvec {

extern void getlambda(double **lambda, const int n);

static unsigned lidx(const unsigned i,
                     const unsigned j,
                     const unsigned n)
{
   return ((n+n-j-1)*j+i+i)/2;
}

static void ginverse2(double *a,
                     const unsigned n,
                     const double tol)
{
   unsigned i,j,k,t,ii,jj;
   double sum,aii,ajj;

   ///////////////////////////////////////////////////
   // column-by-column, we are doing
   // L (lower triangular) in place such that L*L'=A
   ///////////////////////////////////////////////////
   for (j=0; j<n; j++) {  // for column j, L[j][j] is calculated first
      jj = lidx(j,j,n);
      sum = ajj = a[jj];
      for (k=0; k<j; ++k) {
        aii = a[lidx(j,k,n)];
	sum -= aii*aii;
      }
                        //  for other elements L[i,j] in jth column, i=j+1 to n
      if (sum <= -tol) {
         throw exception("ginverse2(): matrix is non-psd");
      }
      if (fabs(sum) < tol) {
         for (i=j; i<n; i++) { a[jj++]=0.0; }
      }
      else {
         a[jj]= ajj = std::sqrt(sum);
         for (i=j+1; i<n; i++) {
            sum = a[lidx(i,j,n)];
            for (k=0; k<j; ++k) {
               sum -= a[lidx(i,k,n)]*a[lidx(j,k,n)];
            }
            a[lidx(i,j,n)]=sum/ajj;
         }
      }
   }

   //////////////////////////////////////////
   //  column-by-column,  we are doing 
   //   inv(L) in place (lower triangl)
   //////////////////////////////////////////
   for (j=0; j<n; j++) {
      jj = lidx(j,j,n);
      if (a[jj] < tol) {
         for (i=j; i<n; i++) {
            a[jj++]=0.0;
         }         
      }
      else {
         a[jj] = 1./a[jj];
         for (i=j+1; i<n; i++) {
	    aii = a[lidx(i,i,n)];
	    jj++;
            if (aii <= -tol) {
               throw exception("ginverse2(): matrix is non-psd");
            }
            else if (fabs(aii) < tol) {
               a[jj]=0.0;
            }
            else {
               for (sum=0.0,k=j; k<i; ++k) sum -= a[lidx(i,k,n)]*a[lidx(k,j,n)];
               a[jj] = sum/aii;
            }      // end of if-else block
         }         // end of for-loop
      }            // end of if-else block
   }               // end of for-loop

   //////////////////////////////////////
   //  column-by-column, we are doing 
   //  inv(L')*inv(L)= inv(a) in place
   //////////////////////////////////////
   for (j=0; j<n; j++) {
     jj = lidx(j,j,n);
      for (i=j; i<n; i++) {
         ii = lidx(i,i,n);
         k = lidx(i,j,n);
         for (sum=0.0,t=i; t<n; ++t) sum += a[ii++]*a[k++];    
         a[jj++]=sum;
      }
   }
   return;
}

static void preconditioning(double *result,
                            const double **M,
                            const double *r,
                            const unsigned n)
{
   unsigned row,i,k,t1,t2,nt;
   const double *pc,*dp;
   double *dpt;
   for (row=0,i=0; i<n; ++i) {
      nt = (unsigned)M[i][0];
      pc = &(M[i][1]);
      k = row;
      dpt = result + row;
      for (t1=0; t1<nt; ++t1) {
        *dpt = 0.0;
        dp = r + k;
        for (t2=0; t2<t1; ++t2) {
	  *dpt += pc[lidx(t1,t2,nt)]*(*dp++);
        }
        for (t2=t1; t2<nt; ++t2) {
	  *dpt += pc[lidx(t2,t1,nt)]*(*dp++);
        }
        dpt++;
        row++;
      }
   }
   return;
 }

/**
 * Setup mixed model equation and return it. The
 * \a rhs is also built-up.
 */
SparseMatrix* Model::setup_mme(Vector<double>* rhs)
{
   /////////////////////////////////////////////////
   // build up MME: hmmec, the coefficient matrix,
   //               rellrhs, right-hand-side
   ///////////////////////////////////////////////////
   if (type == bad_model) throw exception("Model::setup_mme(): bad model");
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return 0;
      }
   }
   if (hmmesize==0) return &hmmec;

   if (type == mixed_model) {
      if (residual_var.num_rows() != numtrait) {
         type = bad_model;
         throw exception("Model::setup_mme(): residual variance has not yet assigned");
      }
   }
   else if (type == fixed_model) {
      if (residual_var.num_rows() != numtrait) {
         residual_var.identity(numtrait,numtrait);
      }
   }
   else {
      throw exception("Model::setup_mme(): inappropriate model");
   }
   hmmec.resize(hmmesize,max_nz);

   Vector<double> localrhs;
   Vector<double> *rhsvec = &localrhs;
   if (rhs) rhsvec = rhs;

   rhsvec->resize(hmmesize);
   inverse_residual_var();
   setup_ww(rhsvec);
   add_G_1();
   non_zero = hmmec.nz();
   return &hmmec;
}

/**
 * Return the BLUP solution.
 * \param method choices are ioc,iod,direct
 * \param stopval stop criteria for solver ioc and iod
 * \param relax relaxation factor for solver ioc
 * \param mxiter maximum number of iterations allowed for solver ioc,ios
 */
Vector<double>* Model::blup(const std::string &method,double stopval,
                    double relax, unsigned mxiter)
{
   Vector<double> *retval = 0;
   if (type == bad_model) throw exception("Model::blup(): bad model");
   if (method == "iod") {
     blup_pccg(stopval,mxiter);
     return &blupsol;
   }
   if (hmmec.nz() == 0) {
      try {
         setup_mme(&rellrhs);
      } catch (exception &er) {
         throw;
      }
   }
   if (hmmec.solve(blupsol,rellrhs,method,relax,stopval,mxiter) == 1) {
      retval = &blupsol;
      //hmmec.close();
   }
   // note that after BLUP, HMME remains intack
   //   hmmec.resize(0,0);
   return retval;
}

/**
 * Return BLUE solution for a fixed model.
 * \param method choices are ioc,iod,direct
 * \param stopval stop criteria for solver ioc and iod
 * \param relax relaxation factor for solver ioc
 * \param mxiter maximum number of iterations allowed for solver ioc,ios
 */
Vector<double>* Model::blue(const std::string &method,double stopval,
                    double relax, unsigned mxiter)
{
   Vector<double> *retval = 0;
   if (type != fixed_model) throw exception("Model::blue(): it must be a fixed model");
   retval = blup(method,stopval,relax,mxiter);
   return retval;
}

/**
 * Return the BLUP solution
 *
 *  The difference between blup() and get_blupsol() is that
 *  blup() will always recompute blup again, while
 *  get_blupsol() gets whatsoever in blupsol, Of course, if it's empty
 *  then it will compute blup solution
 *
 * \sa blup,blue
 */
Vector<double>* Model::get_blupsol(const std::string &method,double stopval,
                    double relax, unsigned mxiter)
{
   Vector<double> *retval = 0;
   if (type == bad_model) throw exception("Model::blup(): bad model");
   if (hmmec.nz() == 0) {
      setup_mme(&rellrhs);
      if (type == bad_model) throw exception("Model::blup(): bad model");
      if (hmmec.solve(blupsol,rellrhs,method,relax,stopval,mxiter)==1) {
         retval = &blupsol;
	 //hmmec.close();
      }
   }
   else {
      retval = &blupsol;
   }
   return retval;
}

/**
 * Add G-inverse
 */
void Model::add_G_1(void)
{
  int done_ped=0;
   for (int t=0; t<numterm; t++) {
      if (term[t].classi() == 'P') {
         add_Ag(t);
	 //hmmec.display(12,24,12,24);
      }
      else if (term[t].classi() == 'R') {
         add_Ig(t);
      }
   }
   //hmmec.display(12,24,12,24);
}

/**
 * add diagonals of G-inverse for pccg
 */
void Model::add_G_1_diag(double **M)
{
  int done_ped=0;
   for (int t=0; t<numterm; t++) {
      if (term[t].classi() == 'P') {
         add_Ag_diag(t,M);
	 //hmmec.display(12,24,12,24);
      }
      else if (term[t].classi() == 'R') {
         add_Ig_diag(t,M);
      }
   }
}


/**
 * Setup mixed model equation without adding G-inverse
 */
void Model::setup_ww(Vector<double>* rhsvec)
{
   unsigned i,j,ii,jj,t1,t2,rec;
   double vy,xval,val;
   //memset(rhsvec->begin(),'\0',sizeof(double)*hmmesize);
   rhsvec->resize(hmmesize,0.0);
   double *rhstmp = rhsvec->begin() -1;
   // for gcc-3.2 replaced 
   // std::istrstream modelfile(modelstr, modelpcount);
   // with 
   std::istringstream modelfile(modelstringstr);
   if (!modelfile) throw exception(" Model::setup_ww(): cannot open binary file");
   double **vep;
   Matrix<double> ve(numtrait,numtrait);
   for (yry=0, rec=0; rec<numrec; rec++) {
     input_pos_val(modelfile);
     if (ntermGdist) {
       ii = rec_indid[rec];
       if (ii > 0) trait_vec[0] -= pop->popmember[ii-1]->xbzu();
     }
     vep = rve[pos_term[numterm]].begin();
     val = xval_term[numterm];        // val = weight variable
     for (t1=0; t1<numtrait; t1++) {
       for (t2=0; t2<numtrait; t2++) ve[t1][t2] = val*vep[t1][t2];
     }
     for (t1=0; t1<numtrait; t1++) {
       vy = 0.0;
       for(t2=0; t2<numtrait; t2++) vy += ve[t1][t2]*trait_vec[t2];
       yry += trait_vec[t1]*vy;
       for (i=0; i<numterm; i++) {
	  ii = term[i].addr[t1];
	  if (ii == 0) continue;
	  xval = xval_term[i];
          for (t2=0; t2<numtrait; t2++) {
              val = xval * ve.at(t1,t2);
              for (j=i; j<numterm; j++) {
                 jj = term[j].addr[t2];
                 if (jj >= ii) hmmec.insert(ii,jj,val*xval_term[j]);
              }
          }
	  rhstmp[ii] += vy*xval_term[i];
       }
     }
   }
   //BRS modelfile.close();
}

/**
 * Input position and value for each term in the model
 *
 *  Note that term[k].addr[t] starts from 1 rather than from 0
 */
void Model::input_pos_val(std::istream& modelfile)
{
   int t,k,nt;
   unsigned startaddr;

   modelfile.read((char *)pos_term, (numterm+1)*sizeof(unsigned));
   modelfile.read((char *)xval_term, (numterm+1)*sizeof(double));
   modelfile.read((char *)trait_vec, numtrait*sizeof(double));

   for (k=0; k<numterm; k++) {
     nt=nt_vec[base_effect[k]];
      startaddr = term[k].start + 1 +(pos_term[k]-1)*nt ;
      for (t=0; t<numtrait; t++) {
	if (term[k].trait[t]) term[k].addr[t] = startaddr + t;
	else                  term[k].addr[t] = 0;
      }
   }
}
void Model::input_pos_val_iod(std::vector<pos_val_node>::iterator it)
{
   int t,k,nt;
   unsigned startaddr;
   for (k=0; k<numterm; k++) {
     nt=nt_vec[base_effect[k]];
      startaddr = term[k].start + 1 + (it->pos_term[k]-1)*nt ;
      for (t=0; t<numtrait; t++) {
	if (term[k].trait[t]) term[k].addr[t] = startaddr + t;
	else                  term[k].addr[t] = 0;
      }
   }
}

void Model::add_Ig(int t, Vector<double>*x_p)
{
   unsigned k,i,ii,jj,nt;
   int t1,t2;
   double *x,*r,value;
   r = mme_times_res.begin() - 1;
   if (x_p) x = x_p->begin() - 1;

   double **v = 0;
   doubleMatrix *tmp = term[var_link[t]].prior->var_matrix();
   if (!tmp) throw exception("Model::add_Ig():  you probably forgot to use M.variance(...)\n"
            "  to set variance for each random effect\n");
   nt=nt_vec[var_link[t]];
   v = tmp->begin();
   Matrix<double> rv(nt,nt);
   for (t1=0; t1<nt; t1++) {
      for(t2=0; t2<nt; t2++) rv[t1][t2]=v[t1][t2];
   }
   ginverse1(rv.begin(),nt,lng0vec[t],1,SESSION.epsilon);

   unsigned na = term[t].nlevel();
   unsigned startaddr = term[t].start + 1;
   for (k=0; k<na; k++) {
      i = startaddr + k*nt;
      for (t1=0; t1<nt; t1++) {
         ii = i + t1;
         for (t2=0; t2<nt; t2++) {
            jj = i + t2;
	    if (jj>=ii) {
	       value = rv[t1][t2];
	       if(x_p) {// iteration on data
  		  r[ii] += value*x[jj];
                  if (jj > ii) r[jj] += value*x[ii];                  
               } else {
                  hmmec.insert(ii,jj,value);
               }
	    }
         }
      }
   }
}

void Model::add_Ig_diag(int t,double **M)
{
  unsigned rec,k,i,ii,jj,kk,nt;
  int t1,t2;
  double **v = 0;
  doubleMatrix *tmp = term[var_link[t]].prior->var_matrix();
  if (!tmp) throw exception("Model::add_Ig():  you probably forgot to use M.variance(...)\n"
			    "  to set variance for each random effect\n");
  nt=nt_vec[var_link[t]];
  v = tmp->begin();
  Matrix<double> rv(nt,nt);
  for (t1=0; t1<nt; t1++) {
    for(t2=0; t2<nt; t2++) rv[t1][t2]=v[t1][t2];
  }
  ginverse1(rv.begin(),nt,lng0vec[t],1,SESSION.epsilon);

  unsigned na = term[t].nlevel();
  for (k=0,i=0; i<t; ++i) k += term[i].nlevel();
  for (rec=0; rec<na; rec++) {
    kk = k + rec;
    for (t1=0; t1<nt; ++t1) {
       for (t2=t1; t2<nt; ++t2) {
          M[kk][lidx(t2,t1,nt)+1] += rv[t1][t2];
       }
    }
  }
}


/**
 * Add the inverse of additive relationship matrix to mixed
 * model quation.
 */
void Model::add_Ag(int t, Vector<double>* x_p)
{
   unsigned t1, t2, k, i, j, ii, jj, iii, jjj;
   double dii,val,value;
   double *x,*r;
   r = mme_times_res.begin() - 1;
   if (x_p) x = x_p->begin() - 1;

   double **var = 0;
   doubleMatrix *tmp = term[var_link[t]].prior->var_matrix();
   if (!tmp) throw exception("Model::add_Ag():  you probably forgot to use M.variance(...)\n"
            "  to set variance for each random effect\n");
   var = tmp->begin();
   int nt=nt_vec[var_link[t]];
   Matrix<double> rvarg(nt,nt);
   for (i=0; i<nt;i++) for (j=0; j<nt; j++) rvarg[i][j]= var[i][j];
   ginverse1(rvarg.begin(),nt,lng0vec[t],1,SESSION.epsilon);

   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   unsigned startaddr = term[t].start + 1;
   unsigned asd[3];
   Individual *I;
   unsigned nanim=popsize-numgroup;
   double f1,f2;
   for (k=0; k<nanim; k++) {
      I = pop->member(k);
      //  cout << popsize << " " << " "<<nanim << " " << k << I << "\n";
      //cout << I->id() << " " << I->father_id() << " " << I->mother_id() << "\n";
      asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
      f1= I->father_inbcoef();
      f2= I->mother_inbcoef();
      
      if(asd[1]> nanim) f1=-1.;
      if(asd[2]> nanim) f2=-1.;
      dii = 4.0/(2.0 -f1 -f2);

      for (i=0; i<3; i++) {
         if (asd[i] == 0) continue;
         ii = startaddr + (asd[i]-1)*nt;
         for (j=0; j<3; j++) {
	    if (asd[j] == 0) continue;
            jj = startaddr + (asd[j]-1)*nt;
            val = lambda[i][j]* dii;
            for (t1=0; t1<nt; t1++) {
               iii = ii + t1;
               for (t2=0; t2<nt; t2++) {
                  jjj = jj + t2;
		  if (jjj>=iii) {
		     value = val*rvarg[t1][t2];
		     if (x_p) {  // iteration on data
  		        r[iii] += value*x[jjj];
                        if (jjj > iii) r[jjj] += value*x[iii];
                     } else {
		        hmmec.insert(iii,jjj,value);
                     }
		  }
               }
            }
         }
      }
   }
}

/**
 * add the diagonals of inverse of additive relationship for pccg 
 * 
 */
void Model::add_Ag_diag(int t,double **M)
{
   unsigned t1, t2,i,j,k,kk,rec;
   double dii,val,value;
   double **var = 0;
   doubleMatrix *tmp = term[var_link[t]].prior->var_matrix();
   if (!tmp) throw exception("Model::add_Ag():  you probably forgot to use M.variance(...)\n"
            "  to set variance for each random effect\n");
   var = tmp->begin();
   int nt=nt_vec[var_link[t]];
   Matrix<double> rvarg(nt,nt);
   for (i=0; i<nt;i++) for (j=0; j<nt; j++) rvarg[i][j]= var[i][j];
   ginverse1(rvarg.begin(),nt,lng0vec[t],1,SESSION.epsilon);

   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   unsigned startaddr = term[t].start + 1;
   unsigned asd[3];
   Individual *I;
   unsigned nanim=popsize-numgroup;
   double f1,f2;

   for (k=0,i=0; i<t; ++i) k += term[i].nlevel();
   for (rec=0; rec<nanim; rec++) {
     I = pop->member(rec);
     asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
     f1= I->father_inbcoef();
     f2= I->mother_inbcoef();
      
     if(asd[1]> nanim) f1=-1.;
     if(asd[2]> nanim) f2=-1.;
     dii = 4.0/(2.0 -f1 -f2);

     for (i=0; i<3; i++) {
       if (asd[i] == 0) continue;
       kk = k + asd[i] - 1;
       val = lambda[i][i]* dii;
       for (t1=0; t1<nt; ++t1) {
           for (t2=t1; t2<nt; ++t2) {
               M[kk][lidx(t2,t1,nt)+1] += val*rvarg[t1][t2];
           }
       }
     }
   }
}




/**
 * Compute the inverse of residual variance matrix
 */
void Model::inverse_residual_var(void)
{
   unsigned i,t,t1,t2;
   double **ve;
   for (i=0; i<npattern; i++) {
      ve = rve[i].begin();
      for (t1=0; t1<numtrait; t1++) for (t2=0; t2<numtrait; t2++) {
         ve[t1][t2] = residual_var[t1][t2];
      }
      for (t=0; t<numtrait; t++) {
         if (pattern[i][t]=='0') {
            for (t1=0; t1<numtrait; t1++) ve[t1][t] = 0.0;  // cross col t
            for (t2=0; t2<numtrait; t2++) ve[t][t2] = 0.0;  // cross row t
         }
      }
      ginverse1(ve,numtrait,lnr0vec[i],1,SESSION.epsilon);
   }
}

/**
 * Compute blup by iteration on data using preconditioned conjugate gradient
 * see J.Anim.Sci (2001) 79:1166-1172
 */

void Model::blup_pccg(double tol, unsigned maxiter){

  if (type == bad_model) throw exception("Model::blup_pccg: bad model");
  if (!data_prepared) {
    if (!prepare_data("iod")) {
      type = bad_model;
      throw exception("Model::blup_pccg: bad model");
    }
  }
  if (type == mixed_model) {
    if (residual_var.num_rows() != numtrait) {
      type = bad_model;
      throw exception("Model::blup_pccg: residual variance has not yet assigned");
    }
  }
  else if (type == fixed_model) {
    if (residual_var.num_rows() != numtrait) {
      residual_var.identity(numtrait,numtrait);
    }
  }
  else {
    throw exception("Model::setup_mme(): inappropriate model");
  }

  inverse_residual_var();
  blupsol.resize(hmmesize,0.0);
  rellrhs.reserve(hmmesize);
  unsigned i,j,k;
  unsigned n;

  unsigned nt = (numtrait+1)*numtrait/2;
  for (n=0,i=0; i<numterm; ++i) n += term[i].nlevel();

  double **M;
  M = new(nothrow) double *[n];
  for (i=0; i<n; ++i) {
    M[i] = new(nothrow) double[nt+1];     //M[i][0] is reserved to store the size
    memset(&(M[i][1]),'\0',sizeof(double)*nt);
  }

  k = 0;
  for (i=0; i<numterm; ++i) {
    nt = nt_vec[i];
    for (j=0; j<term[i].nlevel(); ++j)  M[k++][0] = nt;
  }
  compute_rhs_diag_mme(M);

  for (i=0; i<n; ++i) {
    //        nt = (unsigned)M[i][0];
    //        nt = nt*(nt+1)/2;
    //       for (j=1; j<=nt; ++j) cout << " " <<  M[i][j];
    //        cout << "\n";
     ginverse2(&(M[i][1]),(unsigned)M[i][0],1.0e-15);
  }
  mme_times_res.reserve(hmmesize);
  blupsol.resize(hmmesize);
  Vector<double> r;
  r = rellrhs;

  Vector<double> d;
  d.reserve(hmmesize);
  /////////  d = M*r ///////////////
  preconditioning(d.begin(),(const double **)M,r.begin(),n);

  double delta_new, delta_old, convcr = 1.0;
  unsigned iter = 0;
  double alpha,beta;
  delta_new = r.inner_product(d);
  double delta_0 = delta_new*tol*tol;
  bool done = false;
  while (!done){
    iter++;
    std::cout <<"iteration: " << iter <<std::endl;
     /////// q = Ad /////////
    mme_times(d);// the result is stored in mme_times_res
    alpha = delta_new/(d.inner_product(mme_times_res));
    blupsol += alpha*d;
    if(n%50){
       //////// r = r - alpha*q //////////
      r -= alpha*mme_times_res; 
    }
    else{
	 ///// r = b - Ax /////////
      mme_times(blupsol);
      r = rellrhs - mme_times_res;
    }
    /////////  q = M*r ///////////////
    preconditioning(mme_times_res.begin(),(const double **)M,r.begin(),n);

    delta_old = delta_new;
    delta_new = r.inner_product(mme_times_res);
    beta = delta_new/delta_old;
    d = mme_times_res + beta*d;

    if (delta_new > delta_0 && iter < maxiter){
      done = false;
    }
    else if (delta_new <= delta_0){
      ///// r = b - Ax /////////
      mme_times(blupsol);
      r = rellrhs - mme_times_res;

      /////////  q = M*r ///////////////
      preconditioning(mme_times_res.begin(),(const double **)M,r.begin(),n);

      delta_new = r.inner_product(mme_times_res);
      if(delta_new > delta_0 && iter < maxiter){
	done = false;
      }
      else {
	done = true;
      }
    }
    else {
      done = true;
    }
  }
  //  std::cout <<"r = "<<r;
   mme_times_res.clear();
   r.clear();
}
/**
 * Compute rhs and diagonals of mme
 */
void Model::compute_rhs_diag_mme(double **M) {
  unsigned i,j,k,kk,ii,jj,t1,t2,tt,rec;
  std::vector<pos_val_node>::iterator it;
  double **vep,val,vy,xval;
  memset(rellrhs.begin(),'\0',sizeof(double)*hmmesize);
  double *rhstmp = rellrhs.begin()-1;
  Matrix<double> ve(numtrait,numtrait);

  for (it=pos_val_vector.begin(); it!=pos_val_vector.end(); it++) {
    if(!it->pos_term.empty())  {
      input_pos_val_iod(it);
    }
    else {
      continue;
    }

    vep = rve[it->pos_term[numterm]].begin();
    val = it->xval_term[numterm];        // val = weight variable
    for (t1=0; t1<numtrait; t1++) {
      for (t2=0; t2<numtrait; t2++) ve[t1][t2] = val*vep[t1][t2];
    }

    for (t1=0; t1<numtrait; t1++) {
      vy = 0.0;
      for(t2=0; t2<numtrait; t2++) vy += ve[t1][t2]*it->trait_vec[t2];
      for (i=0; i<numterm; i++) {
	ii = term[i].addr[t1];
	if (ii == 0) continue;
        k = (ii-1)/numtrait;
	xval = it->xval_term[i];
        for (t2=t1; t2<numtrait; t2++) {
           jj = term[i].addr[t2];
           if (jj) {
              val = xval* ve[t1][t2];
	      M[k][lidx(t2,t1,nt_vec[i])+1] += val*xval;
           }
        }
        rhstmp[ii] += vy*xval;
      }
    }
  }

  add_G_1_diag(M);
}

/**
 * Compute mme*x without storing mme
 * The result is stored in mme_times_res
 */
void Model::mme_times(Vector<double> &x)
{
   std::vector<pos_val_node>::iterator it;
   unsigned i,j,ii,jj,t1,t2;
   double vy,xval,val,value;

   double *r = mme_times_res.begin();
   memset(r,'\0',sizeof(double)*hmmesize);
   r--;
   double *xp = x.begin() - 1;

   double **vep;
   Matrix<double> ve(numtrait,numtrait);
   for (it=pos_val_vector.begin(); it!=pos_val_vector.end(); it++) {
     if(!it->pos_term.empty()){
       input_pos_val_iod(it);
     }
     else{
       continue;
     }
     if (ntermGdist) {
       throw exception("Model::mme_times: something does not seem to be compatible with iod");
     }
     vep = rve[it->pos_term[numterm]].begin();
     val = it->xval_term[numterm];        // val = weight variable
     for (t1=0; t1<numtrait; t1++) {
       for (t2=0; t2<numtrait; t2++) ve[t1][t2] = val*vep[t1][t2];
     }

     for (t1=0; t1<numtrait; t1++) {
       for (i=0; i<numterm; i++) {
	 ii = term[i].addr[t1];
         if (ii == 0) continue;
	 xval = it->xval_term[i];
         for (t2=0; t2<numtrait; t2++) {
	   val = xval * ve[t1][t2];
           for (j=i; j<numterm; j++) {
               jj = term[j].addr[t2];
               if (jj >= ii) {
	          value = val*it->xval_term[j];
                  r[ii] += value*xp[jj];
	          if (jj > ii) r[jj] += value*xp[ii];
               }
           }
	 }
       }
     }
   }

   for (int t=0; t<numterm; t++) {
     if (term[t].classi() == 'P') {
       add_Ag(t,&x);
     }
     else if (term[t].classi() == 'R') {
       add_Ig(t,&x);
     }
   }

}


} //////// end of namespace matvec

