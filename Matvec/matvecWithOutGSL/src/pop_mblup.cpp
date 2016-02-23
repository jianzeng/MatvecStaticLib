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

#include <list>
#include <sstream>
#include "session.h"
#include "util.h"
#include "population.h"

namespace matvec {
int INPUT_LINE_WIDTH = 1024;

extern void getlambda(double **lambda, const int n);

static double MY_NEG_HUGE = -1.0e+20;

int pdm_val(const unsigned *goff,const unsigned *gsire,const unsigned *gdam,
            Matrix<double> &pdm, const int ori);
void compute_pdq(const double er, const Matrix<double> &pdm,Matrix<double> &pdq);
double pr_osd(const unsigned goff[],const unsigned gsire[],
             const unsigned gdam[],const int ori);

static int inva22(Matrix<double> &a22);

unsigned Population::input_data(const char fname[],GeneticDist* G)
{
   size_t linewidth = INPUT_LINE_WIDTH;
   char *line;
   if(linewidth>0){
     line = new char [linewidth];
   }
   else {
     line = 0;
   }
   std::ifstream in(fname);
   if (!in) {
     if(line){
       delete [] line;
       line=0;
     }
      throw exception("Population::input_data(): cannot open file");
   }
   unsigned numrec = 0;
   while (in.getline(line,linewidth)) if (validline(line)) numrec++;
   in.close();
   resize(numrec,G);
   popsize = numrec;
   nmarker = G->chrom()[0].locus[0].nallele();
   stdid = 1;
   // for gcc-3.2 replaced 
   // std::istrstream str_in(line,linewidth);
   // with 
   std::istringstream str_in(line);
   

   /////////////////// user attention begins ////////////////////
   // the example of data (pedigree) file
   // the requirements:
   // (1) animal id's must be standard (starting from 1, parents precede
   //                                   their offspring, 0 means missing)
   // (2) marker genotype is represented by two columns: marker_allele ids
   //      allele_id must be standard (starting from 1 and 0 means missing)
   // (3) parental origin of marker genotypes is unknown.
   //
   // 1 0 0  80 1 1
   // 2 0 0 120 0 0
   // 3 1 2  90 1 2
   // 4 0 2 110 1 2
   // 5 3 4 115 1 2
   //
   // note that . is not allowed to represent the missing value
   //
   unsigned ind,sire,dam;
   int a1,a2;
   double y;
   Individual* I;
   numrec = 0;
   in.open(fname);
   while (in.getline(line,linewidth)) {
      if (validline(line) ) {
         str_in >> ind >> sire >> dam >> y >> a1 >> a2;
         if (a1<0 || a1>nmarker || a2<0 || a2>nmarker) throw exception(" Pop.input_data(): invalid allele_id");
         I = popmember[numrec++];
         I->genome0.chromosome[0].locus[0].allele = a1-1;
         I->genome1.chromosome[0].locus[0].allele = a2-1;    // locus[0] is ML
         I->myid = ind;
         if (sire == 0) { I->myfather = 0; }
         else           { I->myfather = popmember[sire-1]; }
         if (dam == 0)  { I->mymother = 0; }
         else           { I->mymother = popmember[dam-1]; }
         I->myrecord[0].double_val(y);
         I->p_origin = 0;
         str_in.seekg(0,std::ios::beg);
      }
   }
   in.close();
   if(line){
     delete [] line;
     line=0;
   }
   return numrec;
   /////////////////// user attention end ////////////////////
}

doubleMatrix Population::relv(const double er)
//****************************************************************************
// relationship matrix for v_effects (additive effect of marked QTL allele)
// (0) all animals must numbered sequentially such that elder has smaller id
// (1) both missing parent(s) and ungenotyped Individuals are allowed
// (2) inbreeding is allowed, too.
//****************************************************************************
{
   if (!stdid) renum();
   std::list<int> ungt_inds;
   std::list<int>::iterator it;
   Individual *I;
   unsigned t,i,j, dim = 2*popsize;
   int nm = 0;
   for (t=0; t<popsize; t++) {
      I = popmember[t];
      if (I->genome0.chromosome[0].locus[0].allele == -1 &&
          I->genome1.chromosome[0].locus[0].allele == -1) {   //marker unknown
         nm++;
         ungt_inds.push_back(t);
      }
   }
   if (nm == 0) return relv_nomissing(er);
          // total possible combinations
   unsigned tpc = static_cast<unsigned>(pow(static_cast<double>(nmarker*nmarker),static_cast<double>(nm)));
   Vector<unsigned> vect(nm);
   Vector<unsigned> ngtvec(nm);
   for (i=0; i<nm; i++) ngtvec[i] = nmarker*nmarker;
   Vector<unsigned> binvec(2), bingtvec(2);
   bingtvec[0] = nmarker;
   bingtvec[1] = nmarker;
   doubleMatrix temp(dim,dim);

   doubleMatrix Gt(dim,dim);
   double tmp,sumpr;
   for (sumpr=0.0,t=0; t<tpc; t++) {
      gtindex(t,nm,ngtvec,vect);
      it = ungt_inds.begin();
      for (i=0; i<nm; i++) {
         I = popmember[*it];
         gtindex(vect[i],2,bingtvec,binvec);
         I->genome0.chromosome[0].locus[0].allele = binvec[0];
         I->genome1.chromosome[0].locus[0].allele = binvec[1];
         I->p_origin = 1;
         it++;
      }
      tmp = llhood_genotype(0,0);           // locus[0] is ML
      if (tmp != MY_NEG_HUGE) {
         sumpr += tmp = std::exp(tmp);
         Gt = relv_nomissing(er);
         for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
               temp[i][j] += Gt[i][j]*tmp;
            }
         }
      }
   }
   if (sumpr > 0.0) {
      for (i=0; i<dim; i++) for (j=0; j<dim; j++) temp[i][j] /= sumpr;
   }
   else {
      throw exception("Population::relv(): parent-offspring genotypes incompatible");
   }
   return temp;
}

doubleMatrix Population::relv_nomissing(const double er)
//****************************************************************************
// this is the same as relv(er), except no ungenotyped Individuals are allowed
//****************************************************************************
{
   if (!stdid) renum();
   unsigned j,t,ii,jj,ii1,jj1,dim=2*popsize;
   unsigned s,d, ss,dd;
   Individual *I;
   double x,w1,w2,pt11,pt12,pt21,pt22,sd11,sd12,sd21,sd22;
   double *row1,*row2,*s_row1,*s_row2,*d_row1,*d_row2;
   Matrix<double> pdm(2,4);
   Matrix<double> pdq(2,4);
   doubleMatrix temp(dim,dim);

   for (t=0; t<dim; t++) memset(temp[t],'\0', sizeof(double)*dim);
   for (t=0; t<popsize; t++) {
      I = popmember[t]; s = I->father_id(); d = I->mother_id();
      ii = 2*t; ii1 = ii+1;
      row1 = temp[ii];   row2 = temp[ii1];
      temp[ii][ii] = 1.0; temp[ii1][ii1] = 1.0;

      if (s && d) {
         compute_pdm(I,pdm);
         compute_pdq(er,pdm,pdq);
         ss = 2*(s-1); dd = 2*(d-1);
         s_row1 = temp[ss]; s_row2 = temp[ss+1];
         d_row1 = temp[dd]; d_row2 = temp[dd+1];

         for (j=0; j<t; j++) {
            jj = 2*j; jj1 = jj+1;
            x = pdq[0][0]*s_row1[jj] + pdq[0][1]*s_row2[jj] +
                pdq[0][2]*d_row1[jj] + pdq[0][3]*d_row2[jj];
            temp[jj][ii] = row1[jj] = x;

            x = pdq[0][0]*s_row1[jj1] + pdq[0][1]*s_row2[jj1] +
                pdq[0][2]*d_row1[jj1] + pdq[0][3]*d_row2[jj1];
            temp[jj1][ii] = row1[jj1] = x;

            x = pdq[1][0]*s_row1[jj] + pdq[1][1]*s_row2[jj] +
                pdq[1][2]*d_row1[jj] + pdq[1][3]*d_row2[jj];
            temp[jj][ii1] = row2[jj] = x;

            x = pdq[1][0]*s_row1[jj1] + pdq[1][1]*s_row2[jj1] +
                pdq[1][2]*d_row1[jj1] + pdq[1][3]*d_row2[jj1];
            temp[jj1][ii1] = row2[jj1] = x;
         }
         sd11 = s_row1[dd]; sd12 = s_row1[dd+1];
         sd21 = s_row2[dd]; sd22 = s_row2[dd+1];
         if (sd11 || sd12 || sd21 || sd22) {
            pt11 = pt12 = pt21 = pt22 = 0.0;
            w1 = pdq[0][0] + pdq[0][1];
            if (w1 > 0.1e-9) {
               pt11 = (pdq[0][0]*pdq[1][2])/w1;
               pt12 = (pdq[0][0]*pdq[1][3])/w1;
               pt21 = (pdq[0][1]*pdq[1][2])/w1;
               pt22 = (pdq[0][1]*pdq[1][3])/w1;
            }
            w2 = pdq[0][2] + pdq[0][3];
            if (w2 > 0.1e-9) {
               pt11 += (pdq[0][2]*pdq[1][0])/w2;
               pt12 += (pdq[0][3]*pdq[1][0])/w2;
               pt21 += (pdq[0][2]*pdq[1][1])/w2;
               pt22 += (pdq[0][3]*pdq[1][1])/w2;
            }
            x = pt11*sd11 + pt12*sd12 + pt21*sd21 + pt22*sd22;
            temp[ii1][ii] = temp[ii][ii1] = x;  // x here is the inbcoef
         }
      }
      else if (!s && d ) {             // sire is missing
         compute_pdm(I,pdm);
         compute_pdq(er,pdm,pdq);
         dd = 2*(d - 1);
         d_row1 = temp[dd]; d_row2 = temp[dd+1];
         for (j=0; j<t; j++) {
            jj = 2*j; jj1 = jj+1;
            x = pdq[0][2]*d_row1[jj] + pdq[0][3]*d_row2[jj];
            temp[jj][ii] = row1[jj] = x;

            x = pdq[0][2]*d_row1[jj1] + pdq[0][3]*d_row2[jj1];
            temp[jj1][ii] = row1[jj1] = x;

            x = pdq[1][2]*d_row1[jj] + pdq[1][3]*d_row2[jj];
            temp[jj][ii1] = row2[jj] = x;

            x = pdq[1][2]*d_row1[jj1] + pdq[1][3]*d_row2[jj1];
            temp[jj1][ii1] = row2[jj1] = x;
         }
      }
      else if (s && !d) {         // dam is missing
         compute_pdm(I,pdm);
         compute_pdq(er,pdm,pdq);
         ss = 2*(s-1);
         s_row1 = temp[ss]; s_row2 = temp[ss+1];
         for (j=0; j<t; j++) {
            jj = 2*j; jj1 = jj+1;
            x = pdq[0][0]*s_row1[jj] + pdq[0][1]*s_row2[jj];
            temp[jj][ii] = row1[jj] = x;

            x = pdq[0][0]*s_row1[jj1] + pdq[0][1]*s_row2[jj1];
            temp[jj1][ii] = row1[jj1] = x;

            x = pdq[1][0]*s_row1[jj] + pdq[1][1]*s_row2[jj];
            temp[jj][ii1] = row2[jj] = x;

            x = pdq[1][0]*s_row1[jj1] + pdq[1][1]*s_row2[jj1];
            temp[jj1][ii1] = row2[jj1] = x;
         }
      }
   }
   return temp;
}

doubleMatrix Population::relv_inv(const double er)
//************************************************
// missing genotype is not allowed
// missing parent is allowed;
// inbreeding is ignored
//************************************************
{
   if (!stdid) renum();
   unsigned i,j,k,t,ii,iii,jj,jjj,s,d, dim = 2*popsize;
   Individual *I;
   double v1,v2;
   unsigned sdi[3];
   Matrix<double> pdm(2,4);
   Matrix<double> pdq(2,6);
   Matrix<double> blk(6,6);
   Matrix<double> dinv(2,2);
   doubleMatrix temp(dim,dim);
   unsigned invert = 0;
   pdq[0][4]=1.0; pdq[0][5]=0.0;
   pdq[1][4]=0.0; pdq[1][5]=1.0;

   for (t=0; t<popsize; t++) {
      I = popmember[t]; s = I->father_id(); d = I->mother_id();
      if (s || d) {
         compute_pdm(I,pdm);
         compute_pdq(er,pdm,pdq);

         for (j=0; j<4; j++) {
            pdq[0][j] *= -1.0;
            pdq[1][j] *= -1.0;
         }
         dinv[0][0] = 1.0;  dinv[0][1] = 0.0;
         dinv[1][0] = 0.0;  dinv[1][1] = 1.0;
         for (i=0; i<2; i++) for (j=0; j<2; j++) {
            for (k=0; k<4; k++) dinv[i][j] -= pdq[i][k]*pdq[j][k];
         }
         invert = inva22(dinv);
         if (!invert) throw exception(" dinv is not invertable");
         for (i=0; i<6; i++) {
            v1 = pdq[0][i]*dinv[0][0] + pdq[1][i]*dinv[1][0];
            v2 = pdq[0][i]*dinv[0][1] + pdq[1][i]*dinv[1][1];
            for (j=0; j<6; j++) blk[i][j] = v1*pdq[0][j] + v2*pdq[1][j];
         }
         sdi[0] = s; sdi[1] = d; sdi[2] = I->id();
         for (i=0; i<3; i++) {
            if (sdi[i] != 0) {
               ii = (sdi[i]-1)*2;
               iii = i*2;
               for (j=0; j<3; j++) {
                  if (sdi[j] != 0) {
                     jj = (sdi[j]-1)*2;
                     jjj = j*2;
                     temp[ii][jj] += blk[iii][jjj];
                     temp[ii][jj+1] += blk[iii][jjj+1];
                     temp[ii+1][jj] += blk[iii+1][jjj];
                     temp[ii+1][jj+1] += blk[iii+1][jjj+1];
                  }
               }
            }
         }
      }
      else {                          // both parents are unknown
         ii = 2*t;
         temp[ii][ii] += 1.0;  temp[ii+1][ii+1] += 1.0;
      }
   }
   return temp;
}

void Population::setup_m_ww(SparseMatrix& mme,Vector<double> &rhs,Vector<int> &start_addr,const double ev_e)
{
   if (!stdid) renum();
   unsigned dim = 2*popsize;
   int n_fixed_effects = 1;    // intercept is the only and default fixed effect
   int fixed_effect_level = 1;

   unsigned n = fixed_effect_level + 3*popsize;
   int nterms = n_fixed_effects + 3;
   unsigned maxnze = est_nze(n,1,1);
   mme.resize(n,maxnze);
   rhs.resize(n,0.0);
   Vector<double> term_val(nterms);
   Vector<int>    term_pos(nterms);
   Vector<int>    term_size(nterms);
           term_size[0] = fixed_effect_level;
           term_size[1] = dim;
           term_size[2] = dim;
           term_size[3] = popsize;

   start_addr.reserve(nterms);
   start_addr[0] = 0;
   start_addr[1] = start_addr[0] + term_size[0];
   start_addr[2] = start_addr[1];      // term1 and 2 start from the same point
   start_addr[3] = start_addr[2] + term_size[2];

   int i,j,k,ii,jj;
   double inve,val,y;
   if (ev_e == 0.0) throw exception(" Population::setup_m_ww(): ev_e must be nonzero");
   inve = 1.0/ev_e;
   const DataNode *trait;

   for (k=0; k<popsize; k++) {
      term_pos[0] = 1;      term_val[0] = 1.0;      // intercept only
      term_pos[1] = 2*k+1;  term_val[1] = 1.0;
      term_pos[2] = 2*k+2;  term_val[2] = 1.0;
      term_pos[3] = k+1;    term_val[3] = 1.0;
      trait = popmember[k]->record();
      y = inve*trait[0].double_val();          // no misiing trait is allowed
      for (i=0; i<nterms; i++) {
         ii = start_addr[i] + term_pos[i];
         val = inve*term_val[i];
         for (j=i; j<nterms; j++) {
            jj = start_addr[j] + term_pos[j];
            mme.insert(ii,jj,val*term_val[j]);
         }
         rhs(ii) += term_val[i]*y;
      }
   }
   return;
}

SparseMatrix* Population::setup_m_mme1(Vector<double>& rhs,const double ev_v,
              const double ev_u,const double ev_e,const double er)
//************************************************
// ungenotype Individual is not allowed
// missing parent is allowed;
// inbreeding is ignored
// intercept is the only and default fixed effect
//************************************************
{
   if (popsize == 0) {
      SparseMatrix *retval = new SparseMatrix;
      return retval;
   }
   if (!stdid) renum();
   int n_fixed_effects = 1;                      // intercept only
   Vector<int> start_addr;
   setup_m_ww(hmmec,rhs,start_addr,ev_e);
   add_Gv_inv(hmmec,start_addr[n_fixed_effects], ev_v, er);
   add_Ga_inv(hmmec,start_addr[n_fixed_effects+2], ev_u);
   hmmec.close();
   return &hmmec;
}

SparseMatrix* Population::setup_m_mme(Vector<double>& rhs,const double ev_v,
              const double ev_u,const double ev_e,const double er)
//************************************************
// ungenotype Individual is not allowed
// missing parent is allowed;
// inbreeding is accounted for
// intercept is the only and default fixed effect
//************************************************
{
   if (popsize == 0) {
      SparseMatrix *retval = new SparseMatrix;
      return retval;
   }
   if (!stdid) renum();
   int n_fixed_effects = 1;                     // intercept only
   Vector<int> start_addr;
   setup_m_ww(hmmec,rhs,start_addr,ev_e);
   add_Ga_inv(hmmec,start_addr[n_fixed_effects+2],ev_u);
   double inv_v = 1.0/ev_v;
   doubleMatrix rel_v = relv(er);
   rel_v.ginv1();
   unsigned i,j,ii,jj,k,kk;
   kk = start_addr[n_fixed_effects] + 1;
   k = 2*popsize;
   for (i=0; i<k; i++) {
      ii = kk + i;
      for (j=i; j<k; j++) {
         jj = kk + j;
         hmmec.insert(ii,jj,inv_v*rel_v[i][j]);
      }
   }
   hmmec.close();
   return &hmmec;
}

Vector<double>* Population::mblup1(const double ev_v,const double ev_u,
                           const double ev_e,const double er)
//************************************************
// ungenotype Individual is not allowed
// missing parent is allowed;
// inbreeding is ignored
// intercept is the only and default fixed effect
//************************************************
{
   if (popsize == 0) {
      Vector<double> *retval = new Vector<double>;
      return retval;
   }
   if (!stdid) renum();
   int fixed_effect_level = 1;       // so far only intercept is allowed
   unsigned dim = 3*popsize + fixed_effect_level;
   Vector<double> rhs(dim);
   setup_m_mme1(rhs,ev_v,ev_u,ev_e,er);

   blupsol.resize(dim);
   SESSION.warning = 0;
   hmmec.solve(blupsol,rhs);
   SESSION.warning  = 1;
   hmmec.resize(0,0);

   return &blupsol;
}

Vector<double>* Population::mblup(const double ev_v,const double ev_u,
                           const double ev_e,const double er)
//************************************************
// ungenotype Individual is allowed
// missing parent is allowed;
// inbreeding is accounted for
// intercept is the only and default fixed effect
//************************************************
{
   if (popsize == 0) {
      Vector<double> *retval = new Vector<double>;
      return retval;
   }
   if (!stdid) renum();
   int fixed_effect_level = 1;       // so far only intercept is allowed
   unsigned dim = 3*popsize + fixed_effect_level;
   Vector<double> rhs(dim);
   inbcoef();
   setup_m_mme(rhs,ev_v,ev_u,ev_e,er);

   blupsol.resize(dim);
   SESSION.warning = 0;
   hmmec.solve(blupsol,rhs);
   SESSION.warning = 1;
   hmmec.resize(0,0);

   return &blupsol;
}

void Population::add_Ga_inv(SparseMatrix& mme, const int start_addr,
                            const double ev_a)
{
   int i,j,k,ii,jj;
   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   for (i=0; i<3; i++) for (j=0; j<3; j++) lambda[i][j] /= ev_a;
   Individual *I;
   double dval;
   unsigned asd[3];

   for (k=0; k<popsize; k++) {
      I = popmember[k];
      asd[0] = I->id(); asd[1] = I->father_id(); asd[2] = I->mother_id();
      dval = 4.0/(2.0 - I->father_inbcoef() - I->mother_inbcoef());
      for (i=0; i<3; i++) {
         if (asd[i] != 0) {
            ii = start_addr + asd[i];
            for (j=0; j<3; j++) {
               if (asd[j] != 0) {
                  jj = start_addr + asd[j];
                  if (jj >= ii) mme.insert(ii,jj,lambda[i][j]*dval);
               }
            }
         }
      }
   }
}

void Population::add_Gv_inv(SparseMatrix& mme,const int start_addr,
                            const double ev_v,const double er)
//************************************************
// ungenotype Individual is not allowed
// missing parent is allowed;
// inbreeding is ignored
//*************************************************
{
   int invert,i,j,k,kk,t,ii,jj,iii,jjj;
   double v1,v2;
   unsigned sdi[3];
   Matrix<double> pdm(2,4);
   Matrix<double> pdq(2,6);
   Matrix<double> blk(6,6);
   Matrix<double> dinv(2,2);
   double inv_v = 1.0/ev_v;                    // ev_u must not be 0.0
   invert = 0;
   pdq[0][4]=1.0; pdq[0][5]=0.0;
   pdq[1][4]=0.0; pdq[1][5]=1.0;

   Individual  *I;
   for (t=0; t<popsize; t++) {
      I = popmember[t];
      if (I->genome0.chromosome[0].locus[0].allele == 0 ||
          I->genome1.chromosome[0].locus[0].allele == 0) {
         throw exception(" Population::add_Gv_inv(...): ind is ungenotyped");
      }
      compute_pdm(I,pdm);
      compute_pdq(er,pdm,pdq);

      dinv[0][0] = 1.0;  dinv[0][1] = 0.0;
      dinv[1][0] = 0.0;  dinv[1][1] = 1.0;
      for (i=0; i<2; i++) {
         for (j=0; j<2; j++) {
            for (k=0; k<4; k++) dinv[i][j] -= pdq[i][k]*pdq[j][k];
         }
      }
      invert = inva22(dinv);
      if (!invert) throw exception(" dinv is not invertable");
      for (i=0; i<6; i++) {
         v1 = pdq[0][i]*dinv[0][0] + pdq[1][i]*dinv[1][0];
         v2 = pdq[0][i]*dinv[0][1] + pdq[1][i]*dinv[1][1];
         for (j=0; j<6; j++) {
            blk[i][j] = (v1*pdq[0][j] + v2*pdq[1][j])*inv_v;
         }
      }
      sdi[0] = I->father_id(); sdi[1] = I->mother_id(); sdi[2] = I->id();
      kk = start_addr + 1;
      for (i=0; i<3; i++) {
         if (sdi[i] != 0) {
            ii = kk + (sdi[i]-1)*2;
            iii = i*2;
            for (j=0; j<3; j++) {
               if (sdi[j] != 0) {
                  jj = kk + (sdi[j]-1)*2 ;
                  if (jj >= ii) {
                     jjj = j*2;
                     mme.insert(ii,jj,blk[iii][jjj]);
                     mme.insert(ii,jj+1,blk[iii][jjj+1]);
                     mme.insert(ii+1,jj+1,blk[iii+1][jjj+1]);
                     if (jj > ii) {
                        mme.insert(ii+1,jj,blk[iii+1][jjj]);
                     }
                  }
               }
            }
         }
      }
   }
}

static int inva22(Matrix<double> &a22)
{
   double tem,det;
   det = (a22[0][0]*a22[1][1] - a22[0][1]*a22[1][0]);
   if (fabs(det) < SESSION.epsilon) return 0;
   det = 1.0/det;
   tem=a22[0][0];
   a22[0][0] = det*a22[1][1];
   a22[1][1] = det*tem;
   a22[0][1] *= -det;
   a22[1][0] *= -det;
   return 1;
}

int pdm_val(const unsigned *goff,const unsigned *gsire,const unsigned *gdam,
            Matrix<double> &pdm, const int ori)
    // compute S, PDM matrix of 2 by 4
{
   int i,j,k;

   for (i=0; i<4; i++) {
      pdm[0][i] = 0.0;
      pdm[1][i] = 0.0;
   }

   if (ori==1) {
      k=0;
      if (goff[0]==gsire[0]) {
         pdm[0][0] = 1.0;
         k++;
      }
      if (goff[0]==gsire[1]) {
         pdm[0][1] = 1.0;
         k++;
      }
      if (k == 0) return 1;
      pdm[0][0] /= static_cast<double>(k);
      pdm[0][1] /= static_cast<double>(k);
      k=0;
      if (goff[1] == gdam[0]) {
         pdm[1][2] = 1.0;
         k++;
      }
      if (goff[1] == gdam[1]) {
         pdm[1][3] = 1.0;
         k++;
      }
      if (k == 0) return 1;
      pdm[1][2] /= static_cast<double>(k);
      pdm[1][3] /= static_cast<double>(k);

      return 0;
   }
   int pattern[2][8], nn[2][4];   // working space
   int gts[2], gtd[2];
   gts[0] = gsire[0] + 1000;      // all these magic numbers are coupled with
   gts[1] = gsire[1] + 2000;      // line starting with
   gtd[0] = gdam[0] + 3000;       // if ((pattern[0][i]%1000)==gtoff[0]
   gtd[1] = gdam[1] + 4000;       // the maximun alleles allowed is 1000

   for (i=0; i<2; i++) {
      k = i*2;
      for (j=0; j<2; j++) {
         k += j;
         pattern[0][k]   = gts[i];
         pattern[1][k]   = gtd[j];
         pattern[0][4+k] = gtd[j];
         pattern[1][4+k] = gts[i];
      }
   }
   int n = 0;
   for (i=0; i<8; i++){
      if ((pattern[0][i]%1000)==goff[0] && (pattern[1][i]%1000)==goff[1]) {
         n++;
      }
      else {
         pattern[0][i]=0; pattern[1][i]=0;
      }
   }

   if (n == 0) return 1;
   for (j=0; j<2; j++) {
      nn[j][0] = 0; nn[j][1] = 0; nn[j][2] = 0; nn[j][3] = 0;
      for (i=0; i<8; i++) {
         k= pattern[j][i];
         if (k) {
            if (k == gts[0]) nn[j][0] += 1;
            if (k == gts[1]) nn[j][1] += 1;
            if (k == gtd[0]) nn[j][2] += 1;
            if (k == gtd[1]) nn[j][3] += 1;
         }
      }
   }
   for (i=0; i<4; i++) {
      pdm[0][i] = static_cast<double>(nn[0][i])/n;
      pdm[1][i] = static_cast<double>(nn[1][i])/n;
   }
   return 0;
}

void Population::compute_pdm(const Individual* I,Matrix<double> &pdm)
{
   int j,k1,k2;
   unsigned gtsire[2],gtdam[2],gtoff[2];
   double pr_missing,p0,pros;
   Matrix<double> pdm_i(2,4);
   Individual *father, *mother;
   LocusStruct *ML = &(prior->chrom()[0].locus[0]);

   father = I->father();
   mother = I->mother();
// int ori = I->p_origin;
// parental origin is unknow: allele1 and allele2
   I->genotype(0,0,gtoff);
   for (j=0; j<4; j++) {
      pdm[0][j] = 0.0;  pdm[1][j] = 0.0;
   }
   if (father && mother) {
      father->genotype(0,0,gtsire);
      mother->genotype(0,0,gtdam);
      if (pdm_val(gtoff,gtsire,gtdam,pdm,0))
         throw exception(" incompatible marker genotypes related to animal");
      }
   else if (father && !mother) {
      //**********************************************************
      // pdm = sum (pdm_d * Pr(d|s,o)) over d
      // Pr(d|s,o) = Pr(o|s,d)*Pr(d) / sum(Pr(o|s,d)*Pr(d)) over d
      //***********************************************************
      father->genotype(0,0,gtsire);
      for (pros=0.0, k1=0; k1<nmarker; k1++) {
         gtdam[0] = k1;
         p0 = ML->allele_freq[k1];
         for (k2=0; k2<nmarker; k2++) {
            gtdam[1] = k2;
            pr_missing = pr_osd(gtoff,gtsire,gtdam,0)*p0*ML->allele_freq[k2];
            if (pr_missing != 0.0) {
               pros += pr_missing;
               pdm_val(gtoff,gtsire,gtdam,pdm_i,0);
               for (j=0; j<4; j++) {
                  pdm[0][j] += pdm_i[0][j] * pr_missing;
                  pdm[1][j] += pdm_i[1][j] * pr_missing;
               }
            }
         }
      }
      if (pros > 0.0) {
         for (j=0; j<4; j++) {pdm[0][j] /= pros;   pdm[1][j] /= pros; }
         for (j=2; j<4; j++) {       // the last two columns are set to zero
            pdm[0][j] = 0;           // because mother is unknown
            pdm[1][j] = 0;
         }
      }
      else {
         throw exception("incompatible marker genotypes related to animal");
      }
   }
   else if (!father && mother) {
      //***********************************************************
      // pdm = sum (pdm_s * Pr(s|d,o)) over s
      // Pr(s|d,o) = Pr(o|s,d)*Pr(s) / sum(Pr(o|s,d)*Pr(s)) over s
      //************************************************************
      mother->genotype(0,0,gtdam);
      for (pros=0.0, k1=0; k1<nmarker; k1++) {
         gtsire[0] = k1;
         p0 = ML->allele_freq[k1];
         for (k2=0; k2<nmarker; k2++) {
            gtsire[1] = k2;
            pr_missing = pr_osd(gtoff,gtsire,gtdam,0)*p0*ML->allele_freq[k2];
            if (pr_missing != 0.0) {
               pros += pr_missing;
               pdm_val(gtoff,gtsire,gtdam,pdm_i,0);
               for (j=0; j<4; j++) {
                  pdm[0][j] += pdm_i[0][j] * pr_missing;
                  pdm[1][j] += pdm_i[1][j] * pr_missing;
               }
            }
         }
      }
      if (pros > 0.0) {
         for (j=0; j<4; j++) {
            pdm[0][j] /= pros;   pdm[1][j] /= pros;
         }
         for (j=0; j<2; j++) {      // the first two columns are set to zero
            pdm[0][j] = 0;          // because father is unknown
            pdm[1][j] = 0;
         }
      }
      else {
         throw exception("incompatible marker genotypes related to animal");
      }
   }
}

void compute_pdq(const double er, const Matrix<double> &pdm,Matrix<double> &pdq )
{
   // C = S*R, where C is the PDQ matrix of 2 by 4
   for (int i=0; i<2; i++) {
      pdq[i][0] = (1-er)*pdm[i][0] +    er *pdm[i][1];
      pdq[i][1] =    er *pdm[i][0] + (1-er)*pdm[i][1];
      pdq[i][2] = (1-er)*pdm[i][2] +    er *pdm[i][3];
      pdq[i][3] =    er *pdm[i][2] + (1-er)*pdm[i][3];
   }
}

double pr_osd(const unsigned goff[],const unsigned gsire[],
              const unsigned gdam[], const int ori)
// if ori==0, then goff[0] is paternal allele, and goff[1] is maternal allele
{
   double pr_s,pr_d,retval;
   int k=0;
   if (goff[0] == gsire[0]) k++;
   if (goff[0] == gsire[1]) k++;
   pr_s = static_cast<double>(k)/2.0;
   k = 0;
   if (goff[1] == gdam[0]) k++;
   if (goff[1] == gdam[1]) k++;
   pr_d = static_cast<double>(k)/2.0;
   retval = pr_s*pr_d;
   if (ori == 0) {
      k = 0;
      if (goff[0] == gdam[0]) k++;
      if (goff[0] == gdam[1]) k++;
      pr_s = static_cast<double>(k)/2.0;
      k = 0;
      if (goff[1] == gsire[0]) k++;
      if (goff[1] == gsire[1]) k++;
      pr_d = static_cast<double>(k)/2.0;
      retval += pr_s*pr_d;
   }
   return retval;
}

double Population::cprob_op(const unsigned goff[],const unsigned gparent[],
                         const int ori)
///////////////////////////////////////
//    cprob_op     = pr(goff/gparent)
//////////////////////////////////////
{
   LocusStruct *ML = &(prior->chrom()[0].locus[0]);
   unsigned  gtunknown[2];
   double pros,p0;
   int k1,k2;
   for (pros=0.0, k1=0; k1<nmarker; k1++) {
      gtunknown[0] = k1;
      p0 = ML->allele_freq[k1];
      for (k2=0; k2<nmarker; k2++) {
         gtunknown[1] = k2;
         pros += pr_osd(goff,gparent,gtunknown,ori)*p0*ML->allele_freq[k2];
      }
   }
   return pros;
}

double Population::llhood_genotype(const unsigned c_id,const unsigned l_id)
/**********************************************************************
* log likelihood of genotypes of Individuals in an arbitrary pedigree
* genotypes can be parental origin  known or not (directioned or not)
* or mixed.
* this is generalized than usual likelihood.
***********************************************************************/
{
   int ori;
   LocusStruct *L = &(prior->chrom()[c_id].locus[l_id]);

   Individual *I, *father, *mother;
   unsigned i,igtype[2], sgtype[2], dgtype[2];
   double tmp,prob;

   for (prob=0.0, i=0; i<popsize; i++) {
      I = popmember[i];
      father = I->father();
      mother = I->mother();
      ori = I->p_origin;
      I->genotype(c_id,l_id,igtype);
      if (father && mother ) {
         father->genotype(c_id,l_id,sgtype);
         mother->genotype(c_id,l_id,dgtype);
         tmp = pr_osd(igtype,sgtype,dgtype,ori);
         if (tmp == 0.0) return MY_NEG_HUGE;
         prob += std::log(tmp);
      }
      else if (father && !mother) {
         father->genotype(c_id,l_id,sgtype);
         tmp = cprob_op(igtype,sgtype,ori);
         if (tmp == 0.0) return MY_NEG_HUGE;
         prob += std::log(tmp);
      }
      else if ( !father && mother) {
         mother->genotype(c_id,l_id,dgtype);
         tmp = cprob_op(igtype,dgtype,ori);
         if (tmp == 0.0) return MY_NEG_HUGE;
         prob += std::log(tmp);
      }
      else {
         tmp = L->allele_freq[igtype[0]] * L->allele_freq[igtype[1]];
         if (ori == 0) {
            if (igtype[0] != igtype[1]) tmp *= 2.0;
         }
         if (tmp == 0.0) return MY_NEG_HUGE;
         prob += std::log(tmp);
      }
   }
   return prob;
}
} ///////// end of namespace matvec

