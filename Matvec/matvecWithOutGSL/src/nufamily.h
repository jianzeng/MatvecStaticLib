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

#ifndef NuFamily_H
#define NuFamily_H

#include "individual.h"

   ///////////////////////////////////////////////////////////////////
   //  a nuclei-family is defined as a family where there must be
   // father, mother and at least one offspring.
   // note that an isolated individul is not a nuclei family
   ///////////////////////////////////////////////////////////////////

namespace matvec {
extern NuFamily* new_NuFamily_vec(const unsigned m);
extern double MIXED_TOL;
class  Population;

//RLF
class Sym2x2
{
 public:
  double a11,a12,a22,k;

  Sym2x2(void){initialize();}
  void initialize(void) {a11=a12=a22=k=0;}
};
class Sym3x3
{
 public:
  double a11,a12,a13,a22,a23,a33,k;

  Sym3x3(void){initialize();}
  void initialize(void) {a11=a12=a13=a22=a23=a33=k=0;}
};
class Sym4x4
{
 public:
  double a11,a12,a13,a14,a22,a23,a24,a33,a34,a44,k;

  Sym4x4(void){initialize();}
  void initialize(void)
  {a11=a12=a13=a14=
       a22=a23=a24=
           a33=a34=
               a44=k=0;}
};
//RLF

/*!
   \class NuFamily  NuFamily.h
   \brief a nuclear family in a pedigree

  \sa Population Pedigree
*/

class NuFamily {
   friend class Population;
   protected:
      unsigned    numoffs, kernal, terminal, numcut;
      Individual  *myfather, *mymother,  **myoffspring;
      std::list<Individual *>      connectors;
      // father_indx is the index of the father in the mother's spouse list
      // mother_indx is the index of the mother in the father's spouse list
      unsigned father_indx, mother_indx;

   public:
      double      workvec[2];
      unsigned    tn_genotype,seq_id, tn_qtl;         // tn_genotype must be assigned
      Individual  *in_connector, *out_connector;
      Population  *pop;

      NuFamily(void);
      ~NuFamily(void) {release();}

      void     display(void) const;
      void     father(Individual *sire) { myfather = sire;}
      void     mother(Individual *dam) { mymother = dam;}
      void     offspring(Individual** offs, unsigned noffs);
      void     noffs(const unsigned n) {numoffs = n;}
      void     pretend_missing(int on, const Vector<double>& genotype_freq);
      void     cutting(Individual *unwelcome);
      void     terminal_peeling(doubleMatrix* tr,doubleMatrix& wspace);
      void     terminal_peeling(doubleMatrix* tr,const int iw,doubleMatrix& wspace);
//      void     peeling_initial(const double* genotype_freq); // commented out in new code
      void     iterative_peeling(doubleMatrix& wspace);
      void     anterior(Individual* I,doubleMatrix* tr,doubleMatrix& wspace);
      void     posterior(Individual* I,Individual* J,doubleMatrix* tr,doubleMatrix& wspace);
      void     release(void);
      double   get_posterior(Individual* I, int excJ,Vector<double> &vec);
      double   fullsibs_prob(Individual* excludeI,doubleMatrix* tr,doubleMatrix& wspace);
      double   log_likelihood(doubleMatrix* tr,doubleMatrix& wspace);

      unsigned     build_connectors(void);
      unsigned     nconnector(void);
      unsigned     father_id(void) const;
      unsigned     mother_id(void) const;
      unsigned     father_gid(void) const;
      unsigned     mother_gid(void) const;
      unsigned     noffs(void) {return numoffs;}
      Individual*  father(void) const {return myfather;}
      Individual*  mother(void) const {return mymother;}
      Individual** offspring(void) const {return myoffspring;}
	  void eliminateGenotypes(unsigned locus);

//BRS
      double multi_llh(Dblock& post_mat);
      void multi_posterior(Individual* I, Individual* J, Dblock& post_mat, int i_peel=0);
      double get_m_posterior(Individual* I, int exclJ ,Dblock& post_mat);
      double multi_fullsibs_prob(Individual* excludeI,Dblock& post_mat, int i_peel=0);
      void multi_anterior(Individual* I, Dblock&  post_mat_m, Dblock&  post_mat_f, int i_peel=0);
      void multi_terminal_peeling(Dblock& post_mat_m, Dblock& post_mat_f,const int iw);
      void multi_ant_post(Dblock& post_mat_m, Dblock& post_mat_f);
      void pretend_multi_missing(int on);
      void multi_initialize(doubleMatrix &pen);

// Add classes for multipoint mixture stuff
  void multi_get_tr(Individual* J, int m_q, int f_q, int Mswitch, int Fswitch);
  void multi_sumint_offspring(unsigned i,int i_peel=0);
  void multi_m_posterior(Individual* I, Individual* J,int i_peel=0);
  void multi_m_anterior(Individual* I,int i_peel=0);
  double multi_m_log_likelihood(void);
  void multi_m_terminal_peeling(const int iw);
  void multi_wksp_resize(int tn_qtl, int max_switches);
  void pretend_multi_m_missing(int on);
  void multi_m_ant_post(void);
  void multi_m_initialize(int tn_qtl);
  void display_scale(void) const;
  static Vector<Vector<Sym2x2 > > child_matrix;
  static Vector<Vector<UNormal > > pm;
  static Vector<Vector<UNormal > > wksp_for_gen;
  static Vector<Vector<Sym3x3 >  > spouse_matrix;
  static Vector<Vector<Sym3x3 > > myfather_matrix;
  static Vector<Vector<Sym4x4 > > mymother_matrix;
  static doubleMatrix tr;
  static Vector<double> weight;
  static doubleMatrix m_weight;
  static doubleMatrix penetrance;
  static doubleMatrix wspace;
//BRS		                        

};
}
#endif
