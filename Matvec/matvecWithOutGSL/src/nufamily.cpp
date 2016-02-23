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

#include "nufamily.h"
#include "population.h"
#include "individual.h"
#include "dblock.h"


namespace matvec {

//RLF  tolerance used in mixed approximation

double MIXED_TOL = 0.00000000001;

//RLF
//BRS: Attempt to define static variables for NuFamily here.

Vector<Vector<Sym2x2> >  NuFamily::child_matrix;
Vector<Vector<UNormal> > NuFamily::pm;
Vector<Vector<UNormal> > NuFamily::wksp_for_gen;
Vector<Vector<Sym3x3> >  NuFamily::spouse_matrix;
Vector<Vector<Sym3x3> >  NuFamily::myfather_matrix;
Vector<Vector<Sym4x4> >  NuFamily::mymother_matrix;
Vector<double> NuFamily::weight;
doubleMatrix NuFamily::m_weight;
doubleMatrix NuFamily::tr;
doubleMatrix NuFamily::wspace;
doubleMatrix NuFamily::penetrance;

static int compare_ind_ptr(const void *a, const void *b)
{
   Individual **x,**y;
   x = (Individual **)a; y = (Individual **)b;
   if (x[0] == y[0]) {
      return 0;
   }
   else {
      return 1;
   }
}

NuFamily::NuFamily(void)
{
   seq_id      = 0;
   kernal      = 0;
   terminal    = 0;
   numoffs     = 0;
   numcut      = 0;
   myfather    = 0;
   mymother    = 0;
   myoffspring = 0;
   tn_genotype = 0;
   tn_qtl      = 4;
   pop         = 0;
   in_connector = 0;
   out_connector = 0;
}

unsigned NuFamily::father_id(void) const
{
   if (myfather) {
      return myfather->id();
   }
   else {
      return 0;
   }
}

unsigned NuFamily::mother_id(void) const
{
   if (mymother) {
      return mymother->id();
   }
   else {
      return 0;
   }
}
unsigned NuFamily::father_gid(void) const
{
   if (myfather) {
      return myfather->group_id;
   }
   else {
      return 0;
   }
}

unsigned NuFamily::mother_gid(void) const
{
   if (mymother) {
      return mymother->group_id;
   }
   else {
      return 0;
   }
}

void NuFamily::release(void)
{
   if (pop) {
      unsigned popsize = pop->size();
      if (mymother->id() > popsize) {
	if(mymother){
	  delete mymother;
	  mymother=0;
	}
      }
      if (myfather->id() > popsize) {
	if(myfather){
	  delete myfather;
	  myfather=0;
	}
      }
      if (myoffspring) {
         for (unsigned i=0; i<numoffs; i++) {
            if (myoffspring[i]->id() > popsize) {
	      if(myoffspring[i]){
		delete myoffspring[i];
		myoffspring[i]=0;
	      }
	    }
         }
	 if(myoffspring){
	   delete [] myoffspring; 
	   myoffspring=0;
	 }
      }
   }
   else {
      if (myoffspring) {
         delete [] myoffspring; myoffspring = 0;
      }
   }
}

void NuFamily::offspring(Individual** offs, unsigned noffs)
{
   if (myoffspring) {
     delete [] myoffspring;
     myoffspring=0;
   }
   numoffs = noffs;
   if(numoffs>0){
     myoffspring = new Individual* [numoffs];
   }
   else {
     myoffspring = 0;
   }
   for (unsigned i=0; i<numoffs; i++) myoffspring[i] = offs[i];
}

double NuFamily::get_posterior(Individual* I, int excJ, Vector<double> &vec)
{
   /////////////////////////////////////////////////////////
   //  if excJ < 0, say -1, means taking all posteriors
   /////////////////////////////////////////////////////////
   int i,j,nsp = I->nspouse();
   double retval,*ve;
   for (i=0; i<tn_genotype; i++) vec[i] = 1.0;
   for (retval=0.0,i=0; i<nsp; i++) {
      if (i != excJ) {
         ve = I->posterior[i].begin();
         for (j=0; j<tn_genotype; j++) vec[j] *= *ve++;
         retval += *ve;
      }
   }
   //   std::cout << "retval from get posterior " << retval << std::endl;
   return retval;
}

double NuFamily::fullsibs_prob(Individual* excludeI,doubleMatrix* trans,doubleMatrix& WSpace)
{
   if (!trans) throw exception("NuFamily::fullsibs_prob(arg1,arg2,arg3): null arg2");
   WSpace.assign(1.0);
   unsigned j,um,uf,uj,nsp;
   double val, scale;
   Vector<double> pen_vec(tn_genotype);
   Vector<double> post_vec(tn_genotype);
   Individual *J;

   for (scale = 0.0, j=0; j<numoffs; j++) {
      J = myoffspring[j];
      if (J != excludeI) {
         nsp = J->nspouse();
         if (nsp) scale += get_posterior(J,-1,post_vec);
         scale += J->get_penetrance(pen_vec);
         for (um=0; um<tn_genotype; um++) {
            for (uf=0; uf<tn_genotype; uf++) {
               val = 0.0;
               if (nsp) {
                  for (uj=0; uj<tn_genotype; uj++) {
                     val += trans[uj][um][uf]*pen_vec[uj]*post_vec[uj];
                  }
               }
               else {
                  for (uj=0; uj<tn_genotype; uj++) {
                     val += trans[uj][um][uf]*pen_vec[uj];
                  }
               }
               WSpace[um][uf] *=  val;
            }
         }
      }
   }
   //  std::cout << "scale from fullsibs " << scale << std::endl;
   return scale;
}

void NuFamily::pretend_missing(int on, const Vector<double>& genotype_freq)
{
   unsigned tng = (unsigned) genotype_freq.size();
   const double* gfreq = genotype_freq.begin();
   myfather->pretend_missing(on,gfreq,tng );
   mymother->pretend_missing(on,gfreq,tng);
   for (unsigned i=0; i<numoffs; i++) {
      myoffspring[i]->pretend_missing(on,gfreq,tng);
   }
}

unsigned NuFamily::build_connectors(void)
{
   connectors.clear();
   if (father_gid() > 1) connectors.push_back(myfather);
   if (mother_gid() > 1) connectors.push_back(mymother);
   for (unsigned j=0; j<numoffs; j++) {
      if (myoffspring[j]->gid() > 1) connectors.push_back(myoffspring[j]);
   }
   return connectors.size();
}

unsigned NuFamily::nconnector(void)
{
   std::list<Individual *>::iterator pos;
   pos=connectors.begin(); //BRS flag
   while (pos != connectors.end()) {
      if ((*pos)->group_id <= 1) { pos = connectors.erase(pos);}
      else { pos++; } //BRS
   }
   return connectors.size();
}

void  NuFamily::cutting(Individual *unwelcome){
unsigned damid = mymother->id();
unsigned sireid = myfather->id();
if (damid > pop->size()) damid -= pop->size();
if (sireid > pop->size()) sireid -= pop->size();
//std::cout << " " << pop->ind_name(unwelcome->id())
//     << " " << pop->ind_name(damid) << " " << pop->ind_name(sireid) << std::endl;
   //Vector<double> genotype_freq = pop->get_genotype_freq();
   Individual *I = new Individual(*unwelcome);
   std::string    analysis_type=pop->analysis_type; //BRS
   unwelcome->group_id -= 1;
   unwelcome->connect_keeper = 1;
   I->reset_id(pop->size()+unwelcome->id());
   I->group_id = 1;
   I->loop_connector = 0;
   I->connect_keeper = 0;
   I->related_family.clear();
   I->offs_tree.clear();

   std::list<Individual *>::iterator pos;
   pos=connectors.begin(); //BRS flag
   while (pos != connectors.end()) {
      if ((*pos) == unwelcome) { pos = connectors.erase(pos);}
      else { pos++; } //BRS

   }
   if (I->numoffs_spouse) {
      delete [] I->numoffs_spouse;
      I->numoffs_spouse = 0;
   }
   if (I->spouselist) {
      delete [] I->spouselist;
      I->spouselist = 0;
      I->numspouse = 0;
   }
   if (I->myoffspring) {
      delete [] I->myoffspring;
      I->myoffspring = 0;
//      I->numoffs = numoffs;    // however, I->myoffspring is not needed
   }
   I->anterior_iw = 0;
   if (analysis_type == "multipoint") { //BRS added this if statement
     doubleMatrix pen(pop->P_ndim,4);
     I->initial_multi_anterior(pen);
     if (I->posterior) {
       delete [] I->posterior;
       I->posterior = 0;
       I->initial_multi_posterior(-1);
       //need to delete m_posterior here also but not sure how.
     }
   }
   else if (analysis_type == "multipoint_m") { //BRS added this if statement
     doubleMatrix pen(pop->P_ndim,4);
     I->initial_multi_m_anterior(tn_qtl);
     if (I->posterior) {
       delete [] I->posterior;
       I->posterior = 0;
	I->initial_multi_m_posterior(tn_qtl,-1);
       //need to delete m_posterior here also but not sure how.
     }
   }
   else {
     Vector<double> genotype_freq = pop->get_genotype_freq(); // BRS moved this here
     I->initial_anterior(genotype_freq.begin(),tn_genotype);
     if (I->posterior) {
       delete [] I->posterior;
       I->posterior = 0;
     }
   }
   if (myfather == unwelcome) {
      I->myfather = 0;
      I->mymother = 0;
      I->numspouse = 1;
      I->spouselist = new Individual*[1];
      I->spouselist[0] = mymother;
      I->p_origin = 0;                 // phenotype is extended (not original)
      if (analysis_type == "multipoint") { //BRS added this if statement
	I->initial_multi_posterior(-1);
      }
      if (analysis_type == "multipoint_m") { //BRS added this if statement
	I->initial_multi_m_posterior(tn_qtl,-1);
      }
      else {
	I->initial_posterior(tn_genotype);
      }
      myfather = I;
      unwelcome->posterior_iw[mother_indx] = 3; // unwelcome has extended ped
      mother_indx = 0;
      return;
   }
   if (mymother == unwelcome) {
      I->myfather = 0;
      I->mymother = 0;
      I->numspouse = 1;
      I->spouselist = new Individual*[1];
      I->spouselist[0] = myfather;
      I->p_origin = 0;                 // phenotype is extended (not original)
      if (analysis_type == "multipoint") { //BRS added this if statement
	I->initial_multi_posterior(-1);
      }
      if (analysis_type == "multipoint_m") { //BRS added this if statement
	I->initial_multi_m_posterior(tn_qtl,-1);
      }
      else {
	I->initial_posterior(tn_genotype);
      }
      mymother = I;
      unwelcome->posterior_iw[father_indx] = 3;  // unwelcome has extended ped
      father_indx = 0;
      return;
   }
   for (unsigned i=0; i<numoffs; i++) {
      if (myoffspring[i] == unwelcome) {
         numcut++;
         I->numoffs = 0;
         I->p_origin = 1;             // new offspring keeps original phenotype
         unwelcome->p_origin = 0;     // the connector's phenotype => extended
         myoffspring[i] = I;
         unwelcome->anterior_iw = 3;   // unwelcome has extended ped
         return;
      }
   }
   throw exception("NuFamily::cutting(unwelcome): you have probably found a bug!");
}

void NuFamily::anterior(Individual* I,doubleMatrix* trans,doubleMatrix& WSpace)
{
   /////////////////////////////////////////////////
   // I must be one of offspring in the NuFamily
   //////////////////////////////////////////////////
   if (!trans) throw exception("NuFamily::anterior(arg1,arg2,arg3): null arg2");
   Vector<double> post_vec(tn_genotype);
   Vector<double> pen_vec(tn_genotype);
   Vector<double> apm(tn_genotype);
   Vector<double> apf(tn_genotype);

   unsigned um,uf,ui;
   double ai,val,**trme,scale = 0.0;
   if ( numoffs > 1 ) {
     scale = fullsibs_prob(I,trans,WSpace);
   }
   scale += mymother->anterior[tn_genotype];
   for (um=0; um<tn_genotype; um++) apm[um] = mymother->anterior[um];
   if (mymother->nspouse() > 1) {
     //new version
      scale += get_posterior(mymother,father_indx,post_vec);
      //old version
      //    scale += get_posterior(mymother,mates_indx[1],post_vec);
      for (um=0; um<tn_genotype; um++) apm[um] *= post_vec[um];
   }

   scale += myfather->anterior[tn_genotype];
   for (uf=0; uf<tn_genotype; uf++) apf[uf] = myfather->anterior[uf];
   if (myfather->nspouse() > 1) {
     scale += get_posterior(myfather,mother_indx,post_vec); //new version

     //	 scale += get_posterior(myfather,mates_indx[0],post_vec); //old version
      for (uf=0; uf<tn_genotype; uf++) apf[uf] *= post_vec[uf];
   }

   scale += I->get_penetrance(pen_vec);
   if (numoffs > 1) {
      for (val=0.0,ui=0; ui<tn_genotype; ui++) {
         trme = trans[ui].begin();
         for (ai=0.0,um=0; um<tn_genotype; um++) {
            for(uf=0; uf<tn_genotype; uf++) {
               ai += apm[um]*trme[um][uf]*WSpace[um][uf]*apf[uf];
            }
         }
         val += I->anterior[ui] = ai*pen_vec[ui];
      }
   }
   else {
      for (val=0.0,ui=0; ui<tn_genotype; ui++) {
         trme = trans[ui].begin();
         for (ai=0.0,um=0; um<tn_genotype; um++) {
            for(uf=0; uf<tn_genotype; uf++) ai += apm[um]*trme[um][uf]*apf[uf];
         }
         val += I->anterior[ui] = ai*pen_vec[ui];
      }
   }
   for (ui=0; ui<tn_genotype; ui++) I->anterior[ui] /= val;
   scale += std::log(val);
   I->anterior[tn_genotype] = scale;
}

void NuFamily::posterior(Individual* I,Individual* J,doubleMatrix* trans,
                           doubleMatrix& WSpace)
{
   /////////////////////////////////////////////////////
   // I and J must be parents of this nuclei family
   // posterior should have enouph space
   /////////////////////////////////////////////////////
   if (!trans) throw exception("NuFamily::posterior(arg1,arg2,arg3,arg4): null arg3");
   unsigned jj,ui,uj,excI;
   double val, scale, *ws, *post_j;
   Vector<double> post_vec(tn_genotype);
   Vector<double> apj(tn_genotype);

   if (I==mymother) {jj = father_indx; excI = mother_indx; }
   //  if (I==myfather) {jj = mates_indx[0]; excI = mates_indx[1]; } //old version
   else             {jj = mother_indx; excI = father_indx; }

   scale = fullsibs_prob(0,trans,WSpace);
   for (uj=0; uj<tn_genotype; uj++) apj[uj] = J->anterior[uj];
   scale += J->anterior[tn_genotype];
   if (J->nspouse() > 1) {
      scale += get_posterior(J,excI,post_vec);
      for (uj=0; uj<tn_genotype; uj++) apj[uj] *= post_vec[uj];
   }
   post_j = I->posterior[jj].begin();
   for (ui=0; ui<tn_genotype; ui++) {
      ws = WSpace[ui];
      for (val=0.0, uj=0; uj<tn_genotype; uj++) val += *ws++ * apj[uj];
      post_j[ui] = val;
   }
   for (val=0.0, uj=0; uj<tn_genotype; uj++) val += post_j[uj];
   for (uj=0; uj<tn_genotype; uj++) post_j[uj] /= val;
   scale += std::log(val);
   post_j[tn_genotype] = scale;
}

double NuFamily::log_likelihood(doubleMatrix* trans,doubleMatrix& WSpace)
{
   if (!kernal) throw exception(" NuFamily::log_likelihood(): must be kernal nufamily");
   Vector<double> post_vec(tn_genotype);
   //Note the change in order!!!!!!!
   double llh = mymother->anterior[tn_genotype];
   posterior(mymother,myfather,trans,WSpace);
   llh += get_posterior(mymother,-1,post_vec);
   double scale = mymother->anterior.inner_product(post_vec);
   llh += std::log(scale);
   return llh;
}

void NuFamily::terminal_peeling(doubleMatrix* trans,doubleMatrix& WSpace)
{
   ////////////////////////////////////////////////////////////
   // REQUIREMENTS: initilization for anterior and posterior
   ////////////////////////////////////////////////////////////
   if (!trans) throw exception("NuFamily::terminal_peeling(trans,...): null trans");
   std::list<Individual *>::iterator pos;
   pos = connectors.begin();
   if (pos == connectors.end()) throw exception(" NuFamily::terminal_peeling(): no connector");
   if (*pos == mymother) {
      posterior(*pos,myfather,trans,WSpace);
   }
   else if (*pos == myfather) {
      posterior(*pos,mymother,trans,WSpace);
   }
   else {
      anterior(*pos,trans,WSpace);
   }
}

void NuFamily::terminal_peeling(doubleMatrix* trans,const int iw,doubleMatrix& WSpace)
{
   ///////////////////////////////////////////////////////////
   // REQUIREMENTS: initilization for anterior and posterior
   ///////////////////////////////////////////////////////////
   if (!trans) throw exception("NuFamily::terminal_peeling(trans,...): null trans");
   std::list<Individual *>::iterator pos;
   pos = connectors.begin();
   if (pos == connectors.end()) throw exception(" NuFamily::terminal_peeling(): no connector");
   if (*pos == mymother) {
      posterior(*pos,myfather,trans,WSpace);
      mymother->posterior_iw[father_indx] = iw;
   }
   else if (*pos == myfather) {
      posterior(*pos,mymother,trans,WSpace);
      myfather->posterior_iw[mother_indx] = iw;
   }
   else {
      anterior(*pos,trans,WSpace);
      (*pos)->anterior_iw = iw;
   }
}

void NuFamily::iterative_peeling(doubleMatrix& WSpace)
{
   //////////////////////////////////////////////////////////////
   // REQUIREMENTS: initilization for anterior and posterior
   //////////////////////////////////////////////////////////////
   pop->anterior_posterior(mymother,WSpace);
   pop->anterior_posterior(myfather,WSpace);
   for (unsigned j=0; j<numoffs; j++) {
      pop->anterior_posterior(myoffspring[j],WSpace);
   }
}

void NuFamily::display(void) const {

  int ident;
  int popsize=pop->size();

  std::cout << " Father " << pop->ind_name(father_id());
  std::cout << " Mother " << pop->ind_name(mother_id());
  /*
  for (unsigned i=0; i<numoffs; i++) {
    ident=myoffspring[i]->id();
    std::cout << " child" << i+1;
    if (ident > popsize) {
      // this is a cut individual so need to get the original one
      ident -= popsize;  // correct id;
      std::cout << " (cut) " << pop->ind_name(ident);
    }
  else
    std::cout << "  " << pop->ind_name(ident);
    // std::cout << " or " <<  pop->ind_name(myoffspring[i]->id());
  }
  */
  std::cout << std::endl;
}

//RLF



//RLF

//BRS

double NuFamily::multi_llh(Dblock& post_mat) {
  unsigned i,j,k;
  double scale;

  // std::cout << "entering multi_llh " << std::endl;
  if (!kernal) throw exception(" NuFamily::log_likelihood(): must be kernal nufamily");
  double llh = myfather->m_anterior_scale;
  //  std::cout << "LLh " << llh << " " << (myfather->m_anterior) << std::endl;
  //  std::cout << "multi_llh m_anterior_scale " << llh << std::endl;
  multi_posterior(myfather,mymother,post_mat);
  post_mat.init(1.0);
  llh += get_m_posterior(myfather,-1,post_mat);
  //  std::cout << "llh from get m posterior " << llh << std::endl;
  scale = 0.0;
  for (k=0; k<pop->P_ndim; k++) {
    for (i=0;i< myfather->n_switches(); i++) {
      for (j=0;j<4; j++) {
	scale += ((post_mat[k][i][j])*(myfather->m_anterior[k][i][j]));
      }
    }
  }
  //   std::cout << "Scale " << scale << std::endl;
  llh += std::log(scale);
  // std::cout << "exiting multi_llh " << std::endl;
  return llh;
}

void NuFamily::multi_terminal_peeling(Dblock& post_mat_m, Dblock& post_mat_f,const int iw) {

  //  std::cout << "entering multi_terminal_peeling " << std::endl;
   std::list<Individual *>::iterator pos;
   pos = connectors.begin();
   if (pos == connectors.end()) throw exception(" NuFamily::multi_terminal_peeling(): no connector");
   if (*pos == myfather) {
     multi_posterior(myfather,mymother,post_mat_f);
     myfather->posterior_iw[mother_indx] = iw;
   }
   else if (*pos == mymother) {
     multi_posterior(mymother,myfather,post_mat_m);
      mymother->posterior_iw[father_indx] = iw;
   }
   else {
     multi_anterior(*pos, post_mat_m, post_mat_f);
      (*pos)->anterior_iw = iw;
   }
   // std::cout << "exiting multi_terminal_peeling " << std::endl;
}

void NuFamily::multi_posterior(Individual* I, Individual* J,
                               Dblock&  post_mat, int i_peel)
{
   /////////////////////////////////////////////////////
   // I and J must be parents of this nuclei family
   // posterior should have enouph space
   /////////////////////////////////////////////////////

   unsigned jj,ui,uj,excI,i,j,k,i_sw,i_p,i_q,j_sw,j_p,j_q,i_w,j_w,tag_id,popsize,ndim, drow,dcol;
   unsigned n_switches = I->n_switches();
   unsigned j_switches = J->n_switches();
   double val, scale=0, *ws, *post_j, i_sum=0, temp;

   popsize=pop->size();
   ndim = pop->P_ndim;

   if (I==myfather) {
     jj = mother_indx;
     excI = father_indx;
     scale = (multi_fullsibs_prob(0,post_mat,i_peel));
     //   std::cout << "scale in mp " << setprecision(14) << scale << std::endl;
   }
   else {
     jj = father_indx;
     excI = mother_indx;
     scale = multi_fullsibs_prob(0,post_mat,i_peel);
   }
   //  std::cout << wspace << std::endl;
   if (i_peel) {
     //Iterative peeling so have to ensure that we are dealing with the
     //original offspring not the cut one.
     // check I
     tag_id=I->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       I=pop->popmember[tag_id-1];
     }
     // check J
     tag_id=J->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       J=pop->popmember[tag_id-1];
     }
   }
   for (k=0;k<ndim;k++){
     for (i=0;i<j_switches; i++){
       for (j=0;j<4;j++){
	 post_mat[k][i][j] = J->m_anterior[k][i][j];
       }
     }
   }
   scale += J->m_anterior_scale;
   //   std::cout << "ant scale " << J->m_anterior_scale << " post mat " <<  post_mat << " Workspace " << wspace << std::endl;
   //  std::cout << "scale after anterior added " << scale << std::endl;
   if (J->nspouse() > 1) {
      temp= get_m_posterior(J,excI,post_mat);
      scale += temp;
   }

   //  std::cout << "scale after get m posterior " << scale << std::endl;
   i_w = 0;
   for (i_p=0; i_p<ndim; i_p++){
     for (i_q=0; i_q<4; i_q++){
       for (i_sw=0;i_sw<n_switches; i_sw++) {
	 j_w = 0;
	 val = 0.0;
	 for (j_p=0; j_p<ndim; j_p++){
	   for (j_q=0; j_q<4; j_q++){
	     for (j_sw=0;j_sw<j_switches; j_sw++) {
	       if (I == myfather) { // mother is the row, father is the column
		 drow=j_w;
		 dcol=i_w;
	       }
	       else {
		 drow=i_w;
		 dcol=j_w;
	       }

	       val += (post_mat[j_p][j_sw][j_q]*wspace[drow][dcol]);
	       //   std::cout << j_p << " " << j_q << " " << j_sw << " " << post_mat[j_p][j_sw][j_q] << " " << wspace[i_w][j_w] << std::endl;
	       j_w++;
	     }
	   }
	 }
	 i_w++;
	 I->m_posterior[jj][i_p][i_sw][i_q] = val;
	 i_sum += val;
       //  std::cout << " jj " << jj << " i_sw " << i_sw << " i_q " << i_q << " m_post " << I->m_posterior[jj][i_sw][i_q] << std::endl;
       }
     }
   }

   for (i_p=0; i_p<ndim; i_p++){
     for (i_q=0; i_q<4; i_q++){
       for (i_sw=0;i_sw<n_switches; i_sw++) {
	 I->m_posterior[jj][i_p][i_sw][i_q] /= i_sum;
       }
     }
   }
   scale += std::log(i_sum);
   I->m_posterior_scale[jj] = scale;
   //    std::cout << " Posterior scale " << scale << std::endl;
}

void NuFamily::multi_anterior(Individual* I,
                               Dblock&  post_mat_m, Dblock&  post_mat_f, int i_peel) {

  unsigned i_q,i_sw,f_q,f_sw,m_q,m_sw,w_row,w_col,ndim,i_p,m_p,f_p;
  unsigned m_switches = mymother->n_switches();
  unsigned f_switches = myfather->n_switches();
  double val, scale=0, *ws, i_sum=0,prob_beta, prob_epsilon , t_val, Prob_PGN;
  int EQTable[4][4]={2,3,2,3,0,1,0,1,1,0,1,0,3,2,3,2};
  int BQTable[4][4]={2,2,3,3,0,0,1,1,1,1,0,0,3,3,2,2};
  int Findex_sw=myfather->index_sw;
  int Mindex_sw=mymother->index_sw;
  int Iindex_sw=I->index_sw;
  int eps=I->eps_sw;
  int bet=I->bet_sw;
  int Bqtl, beta_gamete,  Eqtl,epsilon_gamete,Iswitch,Mswitch,Fswitch,i,j,tag_id,popsize;
  Individual *temp_mother, *temp_father, *temp_I;
  temp_father=myfather;
  temp_mother=mymother;
  temp_I=I;
  popsize=pop->size();
  ndim=pop->P_ndim;

   if (i_peel) {
     //Iterative peeling so have to ensure that we are dealing with the
     //original offspring not the cut one.
     // check myfather
     tag_id=myfather->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       temp_father=pop->popmember[tag_id-1];
     }
     // check mymother
     tag_id=mymother->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       temp_mother=pop->popmember[tag_id-1];
     }
     // check I
     tag_id=I->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       temp_I=pop->popmember[tag_id-1];
     }
   }
    unsigned i_switches = temp_I->n_switches(); // ensure correct entry to switch table

   // we are computing fullsibs_prob matrix with rows for mother and columns for father
   if ( numoffs > 1 )
     scale = multi_fullsibs_prob(I, post_mat_m, i_peel);

   //   std::cout << "scale after multi_fullsibs_prob " << scale << std::endl;

   // don't send temp_I because it is not a child of this nu_family
   for (m_p=0; m_p < ndim; m_p++){
     for (m_sw=0;m_sw<m_switches; m_sw++){
       for (m_q=0;m_q<4;m_q++){
	 post_mat_m[m_p][m_sw][m_q] = temp_mother->m_anterior[m_p][m_sw][m_q];
       }
     }
   }
   scale += temp_mother->m_anterior_scale;

   if (temp_mother->nspouse() > 1) {
      scale += get_m_posterior(temp_mother,father_indx,post_mat_m);
   }

   for (f_p=0; f_p < ndim; f_p++){
     for (f_sw=0;f_sw<f_switches; f_sw++){
       for (f_q=0;f_q<4;f_q++){
	 post_mat_f[f_p][f_sw][f_q] = temp_father->m_anterior[f_p][f_sw][f_q];
       }
     }
   }
   scale += temp_father->m_anterior_scale;

   if (temp_father->nspouse() > 1) {
      scale += get_m_posterior(temp_father,mother_indx,post_mat_f);
   }

   scale += I->get_m_penetrance(penetrance);

   if ( numoffs > 1 ) {
     for (i_p=0; i_p < ndim; i_p++){ // Individuals PGN
     for (i_q=0; i_q<4; i_q++) {
       for (i_sw=0;i_sw<i_switches; i_sw++) {   // loop for genotypes of ind
	 w_row = 0;
	 Prob_PGN = 0.0;
	 Iswitch=pop->switch_table[Iindex_sw][i_sw+1];  // switch for ind (i_sw)
	   for (m_p=0; m_p < ndim; m_p++){ // loop for PGN of mother
	     for (m_q=0;m_q<4;m_q++){
	       Bqtl=BQTable[m_q][i_q];                     // Find which type QTL allele from dam
	       //	 if (Bqtl !=3) {                       // need to be compatible
	       for (m_sw=0;m_sw<m_switches; m_sw++){        // loop for genotypes of mother
		 Mswitch=pop->switch_table[Mindex_sw][m_sw+1];              // switch for mother (m_sw)
		 beta_gamete=pop->switch_table_gmt[bet][(Mswitch^Iswitch)]; // get coded mat marker gamete
		 prob_beta=pop->switch_table_prob[beta_gamete][Bqtl];// Get Pr of mat gamete with qtl
		 w_col = 0;
		 for (f_p=0; f_p < ndim; f_p++){  // loop for PGN of father
		   val = 0.0;
		   for (f_q=0;f_q<4;f_q++){
		     Eqtl=EQTable[f_q][i_q];         // Find which type QTL allele from sire
		     //  if (Eqtl != 3) {                          // need to be compatible
		     for (f_sw=0;f_sw<f_switches; f_sw++){             // loop for genotypes of father
		       Fswitch=pop->switch_table[Findex_sw][f_sw+1];     // switch for father (f_sw)
		       epsilon_gamete=pop->switch_table_gmt[eps][(Fswitch^Iswitch)]; // get coded paternal marker gamete
		       prob_epsilon=pop->switch_table_prob[epsilon_gamete][Eqtl];// Get Pr of paternal gamete
		       if ( (Eqtl != 3) &&  (Bqtl !=3)) {
			 val += post_mat_m[m_p][m_sw][m_q]*wspace[w_row][w_col]*post_mat_f[f_p][f_sw][f_q] *prob_beta*prob_epsilon;
		       }
		       w_col++;              // increment w_col only
		     }
		   }
		   Prob_PGN += val*pop->F->getpr(f_p,m_p,i_p);
		 }
		 w_row++;                     // increment w_row only
	       }
	     }
	   }
	   t_val=Prob_PGN*penetrance[i_p][i_q];
	   temp_I->m_anterior[i_p][i_sw][i_q] = t_val; // store value
	   i_sum += t_val;                   // accumulate sum of values for scaling
       }
     }
     }
   }
   else {
     for (i_p=0; i_p < ndim; i_p++){ // Individuals PGN
       for (i_q=0; i_q<4; i_q++) {
	 for (i_sw=0;i_sw<i_switches; i_sw++) {                             // loop for genotypes of ind
	   Prob_PGN=0.0;
	   Iswitch=pop->switch_table[Iindex_sw][i_sw+1];                    // switch for ind (i_sw)
	   for (m_p=0; m_p < ndim; m_p++){ // loop for PGN of mother
	     for (m_q=0;m_q<4;m_q++){
	       Bqtl=BQTable[m_q][i_q];   // Find which type QTL allele from dam
	       if (Bqtl !=3) {                        // need to be compatible
		 for (m_sw=0;m_sw<m_switches; m_sw++){      // loop for genotypes of mother
		   Mswitch=pop->switch_table[Mindex_sw][m_sw+1];    // switch for mother (m_sw)
		   beta_gamete=pop->switch_table_gmt[bet][(Mswitch^Iswitch)]; // get coded maternal marker gamete
		   prob_beta=pop->switch_table_prob[beta_gamete][Bqtl];       // Get prob of maternal gamete with qtl
		   for (f_p=0; f_p < ndim; f_p++){  // loop for PGN of father
		     val = 0.0;
		     for (f_q=0;f_q<4;f_q++){
		       Eqtl=EQTable[f_q][i_q];    // Find which type QTL allele from sire
		       if (Eqtl != 3) {                   // need to be compatible
			 for (f_sw=0;f_sw<f_switches; f_sw++){      // loop for genotypes of father
			   Fswitch=pop->switch_table[Findex_sw][f_sw+1];   // switch for father (f_sw)
			   epsilon_gamete=pop->switch_table_gmt[eps][(Fswitch^Iswitch)]; // get coded paternal marker gamete
			   prob_epsilon=pop->switch_table_prob[epsilon_gamete][Eqtl];  // Get prob of paternal gamete
			   val += post_mat_m[m_p][m_sw][m_q]*post_mat_f[f_p][f_sw][f_q] *prob_beta*prob_epsilon;
			 }
		       }
		     }
		     Prob_PGN += val*pop->F->getpr(f_p,m_p,i_p);
		   }
		 }
	       }
	     }
	   }
	   t_val=Prob_PGN*penetrance[i_p][i_q];
	   temp_I->m_anterior[i_p][i_sw][i_q] = t_val; // store value
	   i_sum += t_val;                   // accumulate sum of values for scaling
	 }
       }
     }
   }
   for (i_p=0; i_p < ndim; i_p++){ // Individuals PGN
     for (i_q=0; i_q<4; i_q++){
       for (i_sw=0;i_sw<i_switches; i_sw++) {
	 temp_I->m_anterior[i_p][i_sw][i_q] /= i_sum;
       }
     }
   }
   scale += std::log(i_sum);
   temp_I->m_anterior_scale = scale;
   //    std::cout << " id " << temp_I->id() << " Anterior scale " << scale << std::endl;
}

double NuFamily::get_m_posterior(Individual* I, int excJ, Dblock& post_mat)
{

   /////////////////////////////////////////////////////////
   //  if excJ < 0, say -1, means taking all posteriors
   /////////////////////////////////////////////////////////

  // make sure post_mat is initialized properly !!!!!!!!!!!!!!!!

  unsigned n_switches = I->n_switches();
  int spouse,i,j,k,nsp = I->nspouse();
  double retval=0.0;

  for (spouse=0; spouse<nsp; spouse++) {
    if (spouse != excJ) {
      for (k=0;k<pop->P_ndim;k++){
	for (i=0;i<n_switches; i++){
	  for (j=0;j<4;j++){
	    post_mat[k][i][j] *= (I->m_posterior[spouse])[k][i][j];
	  }
	}
      }
      //  std::cout << "spouse " << spouse << " scale " << (I->m_posterior_scale[spouse]) << std::endl;
      retval += I->m_posterior_scale[spouse];
    }
    //    std::cout << I->id() << " return value for spouse " << spouse << " is " << retval << std::endl;

  }
  return retval;
}

double NuFamily::multi_fullsibs_prob(Individual* excludeI, Dblock& post_mat, int i_peel) {
// fullsib probs with rows for mother's genotypes and columns for father's

   wspace.assign(1.0);
   unsigned nsp;
   int j, eps, bet, Jswitch, j_switches, w_row, w_col, m_q, m_sw, Fswitch, Mswitch,
   f_q, f_sw, Bqtl, Eqtl, j_q, j_sw, Jindex_sw, beta_gamete, epsilon_gamete, m_p, f_p, j_p;
   double scale=0.0, ProbIP,  ProbIQ, ProbIS, prob_beta, prob_epsilon;
   Individual *J;
   int EQTable[4][4]={2,3,2,3,0,1,0,1,1,0,1,0,3,2,3,2};
   int BQTable[4][4]={2,2,3,3,0,0,1,1,1,1,0,0,3,3,2,2};
   int tag_id;
   int popsize=pop->size();
   int m_switches = mymother->n_switches();
   int f_switches = myfather->n_switches();
   int Findex_sw=myfather->index_sw;
   int Mindex_sw=mymother->index_sw;
   int ndim = pop->P_ndim;

   for (j=0; j<numoffs; j++) {
     J = myoffspring[j];

     if (i_peel) {
       //Iterative peeling so have to ensure that we are dealing with the
       //original offspring not the cut one.
       tag_id=J->id();
       if (tag_id > popsize) {
	 // this is a cut individual so need to get the original one
	 tag_id -= popsize;  // correct id;
	 J=pop->popmember[tag_id-1];
       }
     }
     if (J != excludeI) {
       nsp = J->nspouse();
       if (nsp) {
	 post_mat.init(1.0);
	 scale += get_m_posterior(J,-1,post_mat);
       }

       scale += J->get_m_penetrance(penetrance);
       eps=J->eps_sw;
       bet=J->bet_sw;
       Jindex_sw=J->index_sw;
       j_switches = J->n_switches();
       w_row=0;  // wspace row index
	 for (m_p=0; m_p < ndim; m_p++) { // for each mother  PGN
	   for (m_q=0; m_q < 4; m_q++) { // for each dam qtl
	     for (m_sw=0; m_sw < m_switches; m_sw++){    // loop for genotypes of mother
	       Mswitch=pop->switch_table[Mindex_sw][m_sw+1]; // get current dam switch
	       w_col=0;  // wspace column index
	       for (f_p=0; f_p < ndim; f_p++) { //For each sire PGN
		 for (f_q=0; f_q < 4; f_q++) { //For each sire QTL
		   for (f_sw=0; f_sw < f_switches; f_sw++){   // loop for genotypes of father
		     Fswitch=pop->switch_table[Findex_sw][f_sw+1];  // switch for father (f_sw)
		     ProbIP=0.0;
		     for (j_p=0; j_p < ndim; j_p++) { // for each child  PGN
		       ProbIQ=0.0;
		       for (j_q=0; j_q < 4; j_q++) { // for each child QTL
			 Bqtl=BQTable[m_q][j_q]; // Find which type QTL allele from dam
			 Eqtl=EQTable[f_q][j_q]; // Find which type QTL allele from sire
			 if ((Eqtl != 3) && (Bqtl !=3)) { // need to be compatible
			   ProbIS=0.0;
			   for (j_sw=0; j_sw < j_switches; j_sw++) { // for each child switch
			     Jswitch=pop->switch_table[Jindex_sw][j_sw+1];
			     beta_gamete=pop->switch_table_gmt[bet][(Mswitch^Jswitch)];  // get coded maternal marker gamete
			     prob_beta=pop->switch_table_prob[beta_gamete][Bqtl];    // Get prob of maternal gamete
			     epsilon_gamete=pop->switch_table_gmt[eps][(Fswitch^Jswitch)]; // get coded paternal marker gamete
			     prob_epsilon=pop->switch_table_prob[epsilon_gamete][Eqtl];  // Get prob of paternal gamete
			     if (nsp) {
			       ProbIS += (prob_epsilon*prob_beta*post_mat[j_p][j_sw][j_q]);
			     }
			     else {
			       ProbIS += (prob_epsilon*prob_beta);
			     }
			   }
			   ProbIQ += (ProbIS*penetrance[j_p][j_q]);
			 }
		       }
		       ProbIP += ProbIQ*pop->F->getpr(f_p,m_p,j_p);
		     }
		     wspace[w_row][w_col] *= ProbIP;
		     w_col++;
		   }
		 }
	       }
	       w_row++;
	     }
	   }
	 }
      }
   }
   //   std::cout << "scale from multi_fullsib " << scale << std::endl;
   return scale;
}

void NuFamily::multi_ant_post(Dblock& post_mat_m, Dblock& post_mat_f) {

// mother_indx or indx[0] is the index of the mother in the father's spouse list
// father_indx or indx[1] is the index of the father in the mother's spouse list

  int spouse_index,i,i_peel=1, tag;
  int popsize=pop->size();
  Individual* child;

// calculate posterior for father through mother

  spouse_index = mother_indx;
  if ((!(myfather->base()) && (myfather->loop_connector)) || (myfather->id() > popsize)) {
    if (!(myfather->posterior_iw[spouse_index] == 2 ||
	  myfather->posterior_iw[spouse_index] == 3)) {
      //    std::cout << "Extending myfather " << myfather->id() << std::endl;
      multi_posterior(myfather ,mymother, post_mat_f,i_peel);
    }
  }

// calculate posterior for mother through father

  spouse_index = father_indx;
  if ((!(mymother->base()) && (mymother->loop_connector)) || (mymother->id() > popsize)) {
    if (!(mymother->posterior_iw[spouse_index] == 2 ||
	  mymother->posterior_iw[spouse_index] == 3)) {
      multi_posterior(mymother ,myfather, post_mat_f,i_peel);
    }
  }

// calculate anteriors for children

  for (i=0;i<numoffs; i++){
    child = myoffspring[i];
    tag=child->id();
    //   std::cout << "Extending child original " << tag << " ";
    if (tag > popsize) {
      child=pop->popmember[tag-popsize-1];
      tag -= popsize;
    }
    //   std::cout << "Extending child " << tag << std::endl;
     if (child->loop_connector) {
      if (!(child->anterior_iw == 2 || child->anterior_iw == 3)) {
	multi_anterior(child, post_mat_m, post_mat_f, i_peel);
      }
     }
  }
}

void NuFamily::pretend_multi_missing(int on) {
  //  unsigned tng = (unsigned) genotype_freq.size();
  // const double* gfreq = genotype_freq.ve;
  myfather->pretend_multi_missing(on,penetrance);
  mymother->pretend_multi_missing(on,penetrance);
  for (unsigned i=0; i<numoffs; i++) {
    myoffspring[i]->pretend_multi_missing(on,penetrance);
  }
}

// Approximation to the Mixture Model (based on the Elston-Stewart alg
// for continous genotypes (Elston et al, 1992; Hum. Her. 42:16-27)
// (Fernando, 1996; Proc. Biometrics Conf. in Germany)

void NuFamily::multi_get_tr(Individual* J, int m_q, int f_q, int Mswitch, int Fswitch){

  int j_q,j_sw,j_switches,Jswitch,beta_gamete,epsilon_gamete,bet,eps,Eqtl,Bqtl;
  double prob_beta,prob_epsilon;
  int EQTable[4][4]={2,3,2,3,0,1,0,1,1,0,1,0,3,2,3,2};
  int BQTable[4][4]={2,2,3,3,0,0,1,1,1,1,0,0,3,3,2,2};

  eps=J->eps_sw;
  bet=J->bet_sw;
 int  Jindex_sw=J->index_sw;
  j_switches = J->n_switches();
  for (j_q=0; j_q < 4; j_q++) { // for each child QTL
    Bqtl=BQTable[m_q][j_q]; // Find which type QTL allele from dam
    Eqtl=EQTable[f_q][j_q]; // Find which type QTL allele from sire
    if ((Eqtl != 3) && (Bqtl !=3)) { // need to be compatible
      for (j_sw=0; j_sw < j_switches; j_sw++) { // for each child switch
	Jswitch=pop->switch_table[Jindex_sw][j_sw+1];
	beta_gamete=pop->switch_table_gmt[bet][(Mswitch^Jswitch)];  // get coded maternal marker gamete
	prob_beta=pop->switch_table_prob[beta_gamete][Bqtl];        // Get prob of maternal gamete
	epsilon_gamete=pop->switch_table_gmt[eps][(Fswitch^Jswitch)]; // get coded paternal marker gamete
	prob_epsilon=pop->switch_table_prob[epsilon_gamete][Eqtl];  // Get prob of paternal gamete
	tr[j_q][j_sw] = prob_beta*prob_epsilon;
      }
    }
    else {
      for (j_sw=0; j_sw < j_switches; j_sw++) {
	tr[j_q][j_sw] = 0.0;
      }
    }
  }
}

void NuFamily::multi_sumint_offspring(unsigned i, int i_peel){
//////////////////////////////////////////////////////////////////////////
// Summation over discrete genotypes and integration over continuous 	//
// genotypes for offspring i. Integration is done first. Then, for each	//
// discrete genotype of mother and father, the mixture that results 	//
// from summing over the discrete genotype of i is approximated by a	//
// univariate normal with the same mean and variance as the mixture.	//
// The contributions from this approximation to the posterior or	//
// anterior are accumulated in "unormal_aprox_for_gen"			//
//////////////////////////////////////////////////////////////////////////

// We are adding the summataion over switches to the previous code

  unsigned j,independent;
  Individual *child = myoffspring[i];
  double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));  // pi = 2.0*std::asin(1.0)
  double u11,u12,u22,v,nu, ve;
  double a11,a12,a13,a22,a23,a33,k;
  double P_var = pop->F->var; // polygenic varaince
  int i_q=0,m_q=0,f_q=0;
  int  Iindex_sw=child->index_sw;
  int   i_switches = child->n_switches();
  double sum1=0.0, sum2=0.0;
  double a33_t, k_t, k_q, a11_q, a13_q, a33_q;
  int i_sw, tag_id;
  int popsize=pop->size();
  ve=pop->residual_var[0][0];
  if (i_peel) {
    //Iterative peeling so have to ensure that we are dealing with the
    //original offspring not the cut one.
    tag_id=child->id();
    if (tag_id > popsize) {
      // this is a cut individual so need to get the original one
      tag_id -= popsize;  // correct id;
      child=pop->popmember[tag_id];
    }
  }

////// contributions from transition function of offspring

      v    =  1.0/(0.5*P_var);
      a22  =  0.25*v;
      a23  = -0.5 *v;
      a33_t  =  v;
      k_t    = std::log(inverse_sqrt2pi*std::sqrt(v));

  for (i_q=0;i_q<tn_qtl;i_q++) {
     ///// contributions from penetrance function of offspring
    ///// skip if phenotype of offspring is missing
	
      if (!(child->record()[0].missing)) {
	v    = 1.0/ve;
	nu   = pop->mean_for_genotype[i_q] - child->record()[0].double_val();
	k_q   = std::log(inverse_sqrt2pi*std::sqrt(v));
	a11_q = nu*nu*v;
	a13_q = nu*v;
	a33_q = v;
      }
      else {
	k_q   = 0.0;
	a11_q = 0.0;
	a13_q = 0.0;
	a33_q = 0.0;
      }

    for (i_sw=0; i_sw < i_switches; i_sw++) {
      a11=a12=a13=a33=k=0.0;

    ///// contributions from posteriors of offspring.
    ///// Skip any posteriors that are not yet calculated.
    ///// Skipping will happen only in iterative peeling.
    ///// In termainal and recursive peeling, posteriors
    ///// are always available when required.

      unsigned numspouse = child->nspouse();

      for (j=0; j<numspouse; j++){
	if (child->mix_posterior[j].done) {
	  v    = child->mix_posterior[j].postvec[i_q][i_sw].tsq;
	  nu   = child->mix_posterior[j].postvec[i_q][i_sw].nu;
	  k   += child->mix_posterior[j].postvec[i_q][i_sw].k;
	  a11 += nu*nu*v;
	  a13 += -nu*v;
	  a33 += v;
	}
      }

//  Accumulate a's where needed
      a11 += (a11_q);
      a13 += (a13_q);
      a33 += (a33_q+a33_t);
      k   += (k_q+k_t);
///// Now calculate U matrix and K

      v   = 1.0/(a33);
      u11 = a11 - a13*a13*v;
      u12 = a12 - a13*a23*v;
      u22 = a22 - a23*a23*v;
      k  += std::log(std::sqrt(4.0*std::asin(1.0)*v));   // pi = 2.0*std::asin(1.0)
      independent = 0;
      if (u22*u22 < MIXED_TOL) {
         wksp_for_gen[i_q][i_sw].k = k - 0.5*u11;
	 sum1 +=  wksp_for_gen[i_q][i_sw].k;
         independent =1;
      }
      else {
//// Calculate Normal mean, variance, and k for genotype g
	v   = 1.0/u22;
	nu  = -u12*v;
	wksp_for_gen[i_q][i_sw].tsq = v;
	wksp_for_gen[i_q][i_sw].nu  = nu;
	wksp_for_gen[i_q][i_sw].k = k - 0.5*(u11 - nu*nu*u22) + std::log(std::sqrt(4.0*std::asin(1.0)*v)); // pi = 2.0*std::asin(1.0)
	sum1 +=  wksp_for_gen[i_q][i_sw].k;
      }
    }
  }

//// For each genotype of mother (gm) and of father (f_q):
//// Calculate mean and variance of mixture. Also calculate k for mixture.
//// First calculate weights. sum1 is for scaling; scaling does not
//// affect weights



  sum1 /= (tn_qtl*i_switches);

  // now wksp_for_gen has exp(wksp_for_gen.k - sum1) !!!

  for (i_q=0; i_q<tn_qtl;i_q++) {
    for (i_sw=0; i_sw < i_switches; i_sw++) {
      wksp_for_gen[i_q][i_sw].k = std::exp(wksp_for_gen[i_q][i_sw].k - sum1);
    }
  }

  int m_switches = mymother->n_switches();
  int f_switches = myfather->n_switches();
  int Findex_sw=myfather->index_sw;
  int Mindex_sw=mymother->index_sw;
  int m_sw, Mswitch, f_sw, Fswitch;

//   Loops over the switches for parents go here

  int w_row=0, w_col;
  for (m_q=0; m_q<tn_qtl;m_q++) {
    for (m_sw=0; m_sw < m_switches; m_sw++){    // loop for genotypes of mother
      w_col = 0;
      Mswitch=pop->switch_table[Mindex_sw][m_sw+1]; // get current dam switch
      for (f_q=0; f_q<tn_qtl;f_q++) {
	for (f_sw=0; f_sw < f_switches; f_sw++){   // loop for genotypes of father
	  Fswitch=pop->switch_table[Findex_sw][f_sw+1];  // switch for father (f_sw)
	  multi_get_tr(child,m_q,f_q,Mswitch,Fswitch); // results are stored in tr
	  //	  std::cout << tr;
	  sum2=0.0;
	  for (i_q=0; i_q<tn_qtl;i_q++) {
	    for (i_sw=0; i_sw < i_switches; i_sw++) {
	      m_weight[i_q][i_sw]= tr[i_q][i_sw] * wksp_for_gen[i_q][i_sw].k;
	      //    std::cout << i_q << " i_sw " << i_sw << " m_ weight  = " << m_weight[i_q][i_sw] << " tr " <<  tr[i_q][i_sw] << " wksp_for_gen " <<  wksp_for_gen[i_q][i_sw].k    << std::endl;
	      sum2 += m_weight[i_q][i_sw];
	    }
	  }
	  //	  std::cout << " sum1= " << sum1 << " sum2= " << sum2 << std::endl;
//// Now the mean and variance for mixture are calculated
	  if (!independent) {
	    nu=0.0;
	    v=0.0;
	    for (i_q=0; i_q<tn_qtl;i_q++) {
	      for (i_sw=0; i_sw < i_switches; i_sw++) {
		m_weight[i_q][i_sw] /= sum2;
		nu += wksp_for_gen[i_q][i_sw].nu *m_weight[i_q][i_sw];
	      }
	    }
	    for (i_q=0; i_q<tn_qtl;i_q++) {
	      for (i_sw=0; i_sw < i_switches; i_sw++) {
		v  += wksp_for_gen[i_q][i_sw].tsq*m_weight[i_q][i_sw]
	          + m_weight[i_q][i_sw]*(wksp_for_gen[i_q][i_sw].nu - nu)
		  *(wksp_for_gen[i_q][i_sw].nu - nu);
	      }
	    }
//// k is the constant for a univariate normal (includes 1/sqrt(2*pi*sigma^2))

	    k = std::log(sum2) + sum1 - std::log(std::sqrt(4.0*std::asin(1.0)*v));  // pi = 2.0*std::asin(1.0)

//// contributions to anterior or posterior are accumulated in child_matrix

	    child_matrix[w_row][w_col].a11 +=  nu*nu/v;
	    child_matrix[w_row][w_col].a12 += -nu/v;
	    child_matrix[w_row][w_col].a22 +=  1.0/v;
	    child_matrix[w_row][w_col].k   +=  k;
	  }
	  else {
	    child_matrix[w_row][w_col].k   += std::log(sum2) + sum1;
	  }
	  w_col++;
	}
      }
      w_row++;
    }
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void NuFamily::multi_m_posterior(Individual* I, Individual* J, int i_peel){
 ///////////////////////////////////////////////////
 // I and J must be parents of this nuclei family //
 ///////////////////////////////////////////////////
  // adding the switch stuff


// A separate posterior has to be calculated for each discrete
// genotype of I.
  unsigned j, excI, jj, independent=0;
  double a11, a12, a13, a22, a23, a33,k,ve;
  double u11, u12, u22, sum1, sum2, v, nu,y, v_y, k_y, nu_y, temp;
  double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));            // pi = 2.0*std::asin(1.0)
  int m_q, f_q, j_q, i_q, offspring;
  int popsize=pop->size();
// indx[0] is the index of the mother in the father's spouse list
// indx[1] is the index of the father in the mother's spouse list

  if (I==myfather) {jj = mother_indx; excI = father_indx; }
  else             {jj = father_indx; excI = mother_indx; }

// switch information for parents

  int i_switches = I->n_switches();
  int j_switches = J->n_switches();
  int Iindex_sw=I->index_sw;
  int Jindex_sw=J->index_sw;
  int wrow, wcol, drow, dcol, i_sw, j_sw, tag_id;
ve=pop->residual_var[0][0];

   if (i_peel) {
     //Iterative peeling so have to ensure that we are dealing with the
     //original individual not the cut one.
     // check I
     tag_id=I->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       I=pop->popmember[tag_id];
     }
     // check J
     tag_id=J->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       J=pop->popmember[tag_id];
     }
   }
   //   std::cout << "Post for " << pop->ind_name(I->id()) << " thru " << pop->ind_name(J->id()) << std::endl;

// First, integrate and sum over each offspring
// Approximate results are stored for each i_q and j_q
// in child_matrix

// Initialize child_matrix for this posterior
  wrow = 0;
  for (i_q=0;i_q<tn_qtl;i_q++){
    for (i_sw=0; i_sw < i_switches; i_sw++){    // loop for genotypes of I
      wcol = 0;
      for (j_q=0;j_q<tn_qtl; j_q++){
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
 // First make sure get right access to child_matrix
	  // MOTHER is the ROW
	  // FATHER is the COLUMN
	  if (I == myfather) {
	    drow=wcol;
	    dcol=wrow;
	  }
	  else {
	    drow=wrow;
	    dcol=wcol;
	  }
	  child_matrix[drow][dcol].initialize();
	  wcol++;
	}
      }
      wrow++;
    }
  }

  for (offspring=0; offspring<numoffs; offspring++) { // sum and int. each offspring
    multi_sumint_offspring(offspring);
  }

// Collect results to integrate uj for each (j_q and sj) in spouse

  for (j_q=0;j_q<tn_qtl;j_q++){

    a11=a12=a13=a22=a23=a33=k=0.0;

// contributions from penetrance funtion of J
// skip if phenotype of J is missing

    if (!(J->record()[0].missing)) {
      y    = J->record()[0].double_val();
      v_y    = 1.0/ve;
      k_y   = std::log(inverse_sqrt2pi*std::sqrt(v_y));
      nu_y   = pop->mean_for_genotype[j_q] - y;

    }
    else{
      nu_y = 0.0;
      v_y = 0.0;
      k_y = 0.0;
    }

    for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J

      a11 = nu_y*nu_y*v_y;
      a13 = nu_y*v_y;
      a33 = v_y;
      k   = k_y;

// contributions from J's anterior

      v    =  J->mix_anterior[j_q][j_sw].tsq;
      nu   =  J->mix_anterior[j_q][j_sw].nu;
      k    += J->mix_anterior[j_q][j_sw].k;
      a11  += nu*nu*v;
      a13  += -nu*v;
      a33  += v;


// Contributions from posteriors of J except I.
// Skip over any posteriors not yet available

      unsigned numspouse = J->nspouse();
      for (j=0; j<numspouse; j++){
	if ((J->mix_posterior[j].done) && (j!=excI) ) {
	  v     = J->mix_posterior[j].postvec[j_q][j_sw].tsq;
	  nu    = J->mix_posterior[j].postvec[j_q][j_sw].nu;
	  k    += J->mix_posterior[j].postvec[j_q][j_sw].k;
	  a11  += nu*nu*v;
	  a13  += -nu*v;
	  a33  += v;
	}
      }

// Put results into spouse_matrix

      spouse_matrix[j_q][j_sw].a11 = a11;
      spouse_matrix[j_q][j_sw].a13 = a13;
      spouse_matrix[j_q][j_sw].a33 = a33;
      spouse_matrix[j_q][j_sw].k   = k;
      //     std::cout << "Spouse matrix j_q " << j_q << " j_sw " << j_sw << "  a11 " << a11 << " a13 " << a13 << " a33 " << a33 << " K " << k << std::endl;
    }
  }

// For each (i_q,i_sw) and (j_q,j_sw) complete integration of uj and store results
// in pm
  wrow = 0;
  for (i_q=0;i_q<tn_qtl;i_q++){
    for (i_sw=0; i_sw < i_switches; i_sw++){    // loop for genotypes of I
      wcol = 0;
      for (j_q=0;j_q<tn_qtl; j_q++){
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J

	  a11=a12=a13=a22=a23=a33=k=0.0;

//    Get contributions from spouse_matrix for genotype j_q

	  a11  = spouse_matrix[j_q][j_sw].a11;
	  a13  = spouse_matrix[j_q][j_sw].a13;
	  a33  = spouse_matrix[j_q][j_sw].a33;
	  k    = spouse_matrix[j_q][j_sw].k;

//   Get contributions from child_matrix for genotypes i_q and j_q
	  // First make sure get right access to child_matrix
	  // MOTHER is the ROW
	  // FATHER is the COLUMN
	  if (I == myfather) {
	    drow=wcol;
	    dcol=wrow;
	  }
	  else {
	    drow=wrow;
	    dcol=wcol;
	  }
	  a11 += child_matrix[drow][dcol].a11;
	  a12 += child_matrix[drow][dcol].a12;
	  a13 += child_matrix[drow][dcol].a12;
	  a22 += child_matrix[drow][dcol].a22;
	  a23 += child_matrix[drow][dcol].a22;
	  a33 += child_matrix[drow][dcol].a22;
	  k   += child_matrix[drow][dcol].k;

//   complete integration of uj

	  v    = 1.0/a33;
	  u11  = a11 - a13*a13*v;
	  u12  = a12 - a13*a23*v;
	  u22  = a22 - a23*a23*v;
	  k   += std::log(std::sqrt(4.0*std::asin(1.0)*v));   // pi = 2.0*std::asin(1.0)

	  if (u22 < MIXED_TOL) {
	    pm[wrow][wcol].k = k - 0.5*u11;
	    independent = 1;
	  }
	  else {

// Calculate Normal mean, variance, and k for i_q and j_q

	    v   = 1.0/u22;
	    nu  = -u12*v;
	    pm[wrow][wcol].tsq = v;
	    pm[wrow][wcol].nu  = nu;
	    pm[wrow][wcol].k = k - 0.5*(u11 - nu*nu*u22) + std::log(std::sqrt(4.0*std::asin(1.0)*v));  // pi = 2.0*std::asin(1.0)
	  }
	wcol++;
	}
      }
      wrow++;
    }
  }
  if (independent) {

// The posterior and anterior likelihoods are independent
// This happens due to missing data
    wrow=0;
    for (i_q=0; i_q<tn_qtl; i_q++){
      for (i_sw=0; i_sw < i_switches; i_sw++){    // loop for genotypes of I
	wcol=0;
	sum1=0.0;
	for(j_q=0;j_q<tn_qtl; j_q++){
	  for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	    sum1 += pm[wrow][wcol].k;
	    wcol++;
	  }
	}
	sum1 /= (tn_qtl*j_switches);
        wcol=0;
	sum2=0.0;
	for(j_q=0;j_q<tn_qtl; j_q++){
	  for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	    sum2 += std::exp(pm[wrow][wcol].k - sum1);
	    wcol++;
	  }
	}

	I->mix_posterior[jj].postvec[i_q][i_sw].nu  = 0.0;
	I->mix_posterior[jj].postvec[i_q][i_sw].tsq = 0.0;
	I->mix_posterior[jj].postvec[i_q][i_sw].k   = std::log(sum2) + sum1;
	wrow++;
      }
    }
    I->mix_posterior[jj].done=1;
    return;
  }


// Now for each i_q sum over j_q. The mixture normal that results from summing over
// j_q is approximated by an univariate normal with the same mean and variance as
// the mixture


  wrow=0;
  for (i_q=0; i_q<tn_qtl; i_q++){
    for (i_sw=0; i_sw < i_switches; i_sw++){    // loop for genotypes of I
// Calculate weights of mixture after scaling by sum1
      wcol=0;
      sum1=0.0;
      for(j_q=0;j_q<tn_qtl; j_q++){
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	  sum1 += pm[wrow][wcol].k;
	  wcol++;
	}
      }
      sum1 /= (tn_qtl*j_switches);
      sum2=0.0;
      wcol=0;
      for(j_q=0;j_q<tn_qtl; j_q++){
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	  temp= std::exp(pm[wrow][wcol].k - sum1);
	  sum2 +=temp;
	  m_weight[j_q][j_sw]=temp;
	  wcol++;
	}
      }
      for (j_q=0; j_q<tn_qtl;j_q++){
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	  m_weight[j_q][j_sw] /= sum2;
	}
      }

// Now the mean and variance for mixture are calculated

      nu = v = k = 0;
      wcol=0;
      for (j_q=0; j_q<tn_qtl;j_q++) {
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	  nu += pm[wrow][wcol].nu*m_weight[j_q][j_sw];
	  wcol++;
	}
      }
      wcol=0;
      for (j_q=0; j_q<tn_qtl;j_q++){
	for (j_sw=0; j_sw < j_switches; j_sw++){   // loop for genotypes of J
	  v  += pm[wrow][wcol].tsq*m_weight[j_q][j_sw]
	        + m_weight[j_q][j_sw]*(pm[wrow][wcol].nu - nu)*(pm[wrow][wcol].nu - nu);
	  wcol++;
	}
      }
// k is the constant for a univariate normal (includes 1/sqrt(2*pi*sigma^2))

      k = std::log(sum2) + sum1 - std::log(std::sqrt(4.0*std::asin(1.0)*v));   // pi = 2.0*std::asin(1.0)
      I->mix_posterior[jj].postvec[i_q][i_sw].nu  = nu;
      I->mix_posterior[jj].postvec[i_q][i_sw].tsq = 1.0/v;
      I->mix_posterior[jj].postvec[i_q][i_sw].k   = k;
      wrow++;
      //     std::cout <<  i_q << " " << i_sw << " sum2 " <<  log(sum2)   << " sum1 " << sum1 << " k " << k << std::endl;
    }
  }
  I->mix_posterior[jj].done=1;
}

void NuFamily::multi_m_anterior(Individual* I,int i_peel){

  double a11,a12,a13,a14,a22,a23,a24,a33,a34,a44,k;
  double u11, u12, u22, sum1, sum2, v, nu, y, v_y, k_y, nu_y;
  double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));     // pi = 2.0*std::asin(1.0)
  unsigned i_q,m_q,f_q,j,independent=0;
  int m_sw, f_sw, i_sw, offspring;
  int m_switches = mymother->n_switches();
  int f_switches = myfather->n_switches();
  int Mindex_sw = mymother->index_sw;
  int Findex_sw = myfather->index_sw;
  int EQTable[4][4]={2,3,2,3,0,1,0,1,1,0,1,0,3,2,3,2};
  int BQTable[4][4]={2,2,3,3,0,0,1,1,1,1,0,0,3,3,2,2};
  int eps=I->eps_sw;
  int bet=I->bet_sw;
  int  Iindex_sw=I->index_sw;
  int i_switches = I->n_switches();
  int wrow=0, wcol=0, Iswitch, Bqtl, Mswitch, beta_gamete, Eqtl, Fswitch, epsilon_gamete;
  double prob_beta, prob_epsilon, temp;
  Individual *temp_mother, *temp_father, *temp_I;
  temp_father=myfather;
  temp_mother=mymother;
  temp_I=I;
  int popsize=pop->size();
  int tag_id;
  double P_var = pop->F->var; // polygenic varaince
// indx[0] is the index of the mother in the father's spouse list
// indx[1] is the index of the father in the mother's spouse list

  unsigned excM = mother_indx;
  unsigned excF = father_indx;
double ve=pop->residual_var[0][0];

  if (i_peel) {
     //Iterative peeling so have to ensure that we are dealing with the
     //original individual not the cut one.
     // check myfather
     tag_id=myfather->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       temp_father=pop->popmember[tag_id];
     }
     // check mymother
     tag_id=mymother->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       temp_mother=pop->popmember[tag_id];
     }
     // check I
     tag_id=I->id();
     if (tag_id > popsize) {
       // this is a cut individual so need to get the original one
       tag_id -= popsize;  // correct id;
       temp_I=pop->popmember[tag_id];
     }
   }
  //   std::cout << "Anterior for " <<  pop->ind_name(temp_I->id()) << std::endl;
// First, integrate and sum over each offspring
// Approximate results are stored for each gm and f_q
// in child_matrix

// Initialize child_matrix for this anterior
  for (m_q=0;m_q<tn_qtl;m_q++){
    for (m_sw=0; m_sw < m_switches; m_sw++){    // loop for genotypes of I
      wcol = 0;
      for (f_q=0;f_q<tn_qtl; f_q++){
	for (f_sw=0; f_sw < f_switches; f_sw++){   // loop for genotypes of J
	  child_matrix[wrow][wcol].initialize();
	  wcol++;
	}
      }
      wrow++;
    }
  }

// sum and int. each offspring except I
  for (offspring=0; offspring<numoffs; offspring++) { // sum and int. each offspring
    if ( temp_I != myoffspring[offspring]){
      multi_sumint_offspring(offspring);
    }
  }

// Collect results to integrate uf for each f_q

  for (f_q=0;f_q<tn_qtl;f_q++){
// contributions from penetrance function of father
// skip if phenotype of father is missing

    if (!(temp_father->record()[0].missing)) {
      y    = temp_father->record()[0].double_val();
      v_y  = 1.0/ve;
      k_y  = std::log(inverse_sqrt2pi*std::sqrt(v_y));
      nu_y = pop->mean_for_genotype[f_q] - y;
    }
    else{
      nu_y = 0.0;
      v_y = 0.0;
      k_y = 0.0;
    }

    for (f_sw=0; f_sw < f_switches; f_sw++) {   // loop for genotypes of J

      a11 = nu_y*nu_y*v_y;
      a13 = nu_y*v_y;
      a33 = v_y;
      k   = k_y;
// contributions from Father's anterior

      v   = temp_father->mix_anterior[f_q][f_sw].tsq;
      nu  = temp_father->mix_anterior[f_q][f_sw].nu;
      k   += temp_father->mix_anterior[f_q][f_sw].k;
      a11 += nu*nu*v;
      a13 += -nu*v;
      a33 += v;

// Contributions from posteriors of father except mother.
// Skip over any posteriors not yet available

      unsigned numspouse = temp_father->nspouse();

      for (j=0; j<numspouse; j++){
	if ((temp_father->mix_posterior[j].done) &&  (j!=excM) ) {
	  v     = temp_father->mix_posterior[j].postvec[f_q][f_sw].tsq;
	  nu    = temp_father->mix_posterior[j].postvec[f_q][f_sw].nu;
	  k    += temp_father->mix_posterior[j].postvec[f_q][f_sw].k;
	  a11  += nu*nu*v;
	  a13  += -nu*v;
	  a33  += v;
	}
      }

// Put results into temp_father matrix

      myfather_matrix[f_q][f_sw].a11 = a11;
      myfather_matrix[f_q][f_sw].a13 = a13;
      myfather_matrix[f_q][f_sw].a33 = a33;
      myfather_matrix[f_q][f_sw].k   = k;
    }
  }

// Collect results to integrate um for each m_q

  for (m_q=0;m_q<tn_qtl;m_q++){

// contributions from penetrance funtion of mother
// skip if phenotype of mother is missing

    if (!(temp_mother->record()[0].missing)) {
      y    = temp_mother->record()[0].double_val();
      v_y    = 1.0/ve;
      k_y   = std::log(inverse_sqrt2pi*std::sqrt(v_y));
      nu_y   = pop->mean_for_genotype[m_q] - y;
    }
    else{
      nu_y = 0.0;
      v_y = 0.0;
      k_y = 0.0;
    }
    for (m_sw=0; m_sw < m_switches; m_sw++){   // loop for switches of mother

      a11 = nu_y*nu_y*v_y;
      a14 = nu_y*v_y;
      a44 = v_y;
      k   = k_y;
// contributions from Mother's anterior

      v   = temp_mother->mix_anterior[m_q][m_sw].tsq;
      nu  = temp_mother->mix_anterior[m_q][m_sw].nu;
      k   += temp_mother->mix_anterior[m_q][m_sw].k;
      a11 += nu*nu*v;
      a14 += -nu*v;
      a44 += v;

// Contributions from posteriors of mother except father.
// Skip over any posteriors not yet available

      unsigned numspouse = temp_mother->nspouse();
      for (j=0; j<numspouse; j++){
	if (( temp_mother->mix_posterior[j].done) &&  (j!=excF) ) {
	  v     = temp_mother->mix_posterior[j].postvec[m_q][m_sw].tsq;
	  nu    = temp_mother->mix_posterior[j].postvec[m_q][m_sw].nu;
	  k    += temp_mother->mix_posterior[j].postvec[m_q][m_sw].k;
	  a11  += nu*nu*v;
	  a14  += -nu*v;
	  a44  += v;
	}
      }

// Put results into mother_matrix

      mymother_matrix[m_q][m_sw].a11 = a11;
      mymother_matrix[m_q][m_sw].a14 = a14;
      mymother_matrix[m_q][m_sw].a44 = a44;
      mymother_matrix[m_q][m_sw].k   = k;
    }
  }

// For each i_q and j_q complete integration of um and uf
// and store results in pm
  sum1=0.0;
  wrow=0;
  for (m_q=0; m_q<tn_qtl; m_q++) {
    for (m_sw=0; m_sw < m_switches; m_sw++){   // loop for switches of mother
      wcol=0;
      for (f_q=0; f_q<tn_qtl; f_q++){
	for (f_sw=0; f_sw < f_switches; f_sw++){   // loop for switches of father

	  a11=a12=a13=a14=a22=a23=a24=a33=a34=a44=k=0.0;

//  Get contributions from mother_vector for genotype m_q

	  a11 = mymother_matrix[m_q][m_sw].a11;
	  a14 = mymother_matrix[m_q][m_sw].a14;
	  a44 = mymother_matrix[m_q][m_sw].a44;
	  k   = mymother_matrix[m_q][m_sw].k  ;

//   Get contributions from child_matrix for genotypes f_q and m_q

	  a11 += child_matrix[wrow][wcol].a11;
	  a13 += child_matrix[wrow][wcol].a12;
	  a14 += child_matrix[wrow][wcol].a12;
	  a33 += child_matrix[wrow][wcol].a22;
	  a34 += child_matrix[wrow][wcol].a22;
	  a44 += child_matrix[wrow][wcol].a22;
	  k   += child_matrix[wrow][wcol].k;

// Contributions from transition function of I to integrate um

	  v    = 1.0/(0.5*P_var);
	  a22 += v;
	  a23 += -0.5*v;
	  a24 += -0.5*v;
	  a33 +=  0.25*v;
	  a34 +=  0.25*v;
	  a44 +=  0.25*v;
	  k   +=  std::log(inverse_sqrt2pi*std::sqrt(v));

// Integrate um

	  v   =  1.0/a44;
	  a11 =	a11 - a14*a14*v;
	  a12 =	a12 - a14*a24*v;
	  a13 =	a13 - a14*a34*v;
	  a22 =	a22 - a24*a24*v;
	  a23 =	a23 - a24*a34*v;
	  a33 =	a33 - a34*a34*v;
	  k  += std::log(std::sqrt(4.0*std::asin(1.0)*v));      // pi = 2.0*std::asin(1.0)

// Get contributions from father

	  a11 += myfather_matrix[f_q][f_sw].a11;
	  a13 += myfather_matrix[f_q][f_sw].a13;
	  a33 += myfather_matrix[f_q][f_sw].a33;
	  k   += myfather_matrix[f_q][f_sw].k  ;

 // Integrate uf

	  v   = 1.0/a33;
	  u11 = a11 - a13*a13*v;
	  u12 = a12 - a13*a23*v;
	  u22 = a22 - a23*a23*v;
	  k  += std::log(std::sqrt(4.0*std::asin(1.0)*v));            // pi = 2.0*std::asin(1.0)

	  if(u22 < MIXED_TOL) {
	    pm[wrow][wcol].k   = k - 0.5*u11;
	    independent = 1;
	  }
	  else {

// Calculate Normal mean, variance and k for f_q and m_q

	    v  = 1.0 /u22;
	    nu  = -u12*v;
	    pm[wrow][wcol].tsq = v;
	    pm[wrow][wcol].nu  = nu;
	    pm[wrow][wcol].k   = k - 0.5*(u11 - nu*nu*u22)
                                  + std::log(std::sqrt(4.0*std::asin(1.0)*v));    // pi = 2.0*std::asin(1.0)
	    sum1 +=  pm[wrow][wcol].k;
	  }
	  wcol++;
	}
      }
      wrow++;
    }
  }
  sum1 /= (tn_qtl*f_switches*tn_qtl*m_switches);

  if (independent) {
// The posterior and anterior likelihoods are independent
// This happens due to missing data

    for (i_q=0; i_q<tn_qtl; i_q++) {
      for (i_sw=0; i_sw < i_switches; i_sw++){    // loop for genotypes of I
	Iswitch=pop->switch_table[Iindex_sw][i_sw+1];  // switch for ind (i_sw)
	wrow=0;
	sum2=0.0;
	for(m_q=0;m_q<tn_qtl; m_q++){
	  Bqtl=BQTable[m_q][i_q];        // Find which type QTL allele from dam
	  for (m_sw=0; m_sw < m_switches; m_sw++){   // loop for genotypes of J
	    Mswitch=pop->switch_table[Mindex_sw][m_sw+1]; // switch for mother (m_sw)
	    beta_gamete=pop->switch_table_gmt[bet][(Mswitch^Iswitch)];  // get coded maternal marker gamete
	    prob_beta=pop->switch_table_prob[beta_gamete][Bqtl];        // Get prob of maternal gamete
	    wcol=0;
	    for (f_q=0; f_q<tn_qtl; f_q++) {
	      Eqtl=EQTable[f_q][i_q];         // Find which type QTL allele from sire
	      for (f_sw=0; f_sw < f_switches; f_sw++){    // loop for genotypes of father
		Fswitch=pop->switch_table[Findex_sw][f_sw+1];   // switch for father (f_sw)
		epsilon_gamete=pop->switch_table_gmt[eps][(Fswitch^Iswitch)]; // get coded paternal marker gamete
		prob_epsilon=pop->switch_table_prob[epsilon_gamete][Eqtl];  // Get prob of paternal gamete
		if ( (Eqtl != 3) &&  (Bqtl !=3)) {
		  sum2 += prob_beta*prob_epsilon*std::exp(pm[wrow][wcol].k - sum1);
		}
		wcol++;
	      }
	    }
	    wrow++;
	  }
	}

	temp_I->mix_anterior[i_q][i_sw].nu  = 0.0;
	temp_I->mix_anterior[i_q][i_sw].tsq = 0.0;
	temp_I->mix_anterior[i_q][i_sw].k = std::log(sum2) + sum1;
	//	std::cout <<  i_q << " " << i_sw << " k " << temp_I->mix_anterior[i_q][i_sw].k << std::endl;
      }
    }
    return;
  }

// For each i_q, approximate the mixture that results from summing over
// f_q and m_q by an univariate normal with the same mean and variance
// as the mixture

  for (i_q=0; i_q<tn_qtl; i_q++) {
    for (i_sw=0; i_sw < i_switches; i_sw++){    // loop for genotypes of I
      Iswitch=pop->switch_table[Iindex_sw][i_sw+1];  // switch for ind (i_sw)
      wrow=0;

// Weights of mixture are calculated here.
      sum2=0.0;
      for(m_q=0;m_q<tn_qtl; m_q++){
	Bqtl=BQTable[m_q][i_q];        // Find which type QTL allele from dam
	for (m_sw=0; m_sw < m_switches; m_sw++){   // loop for genotypes of J
	  Mswitch=pop->switch_table[Mindex_sw][m_sw+1]; // switch for mother (m_sw)
	  beta_gamete=pop->switch_table_gmt[bet][(Mswitch^Iswitch)];  // get coded maternal marker gamete
	  prob_beta=pop->switch_table_prob[beta_gamete][Bqtl];        // Get prob of maternal gamete
	  wcol=0;
	  for (f_q=0; f_q<tn_qtl; f_q++) {
	    Eqtl=EQTable[f_q][i_q];         // Find which type QTL allele from sire
	    for (f_sw=0; f_sw < f_switches; f_sw++){    // loop for genotypes of father
	      Fswitch=pop->switch_table[Findex_sw][f_sw+1];   // switch for father (f_sw)
	      epsilon_gamete=pop->switch_table_gmt[eps][(Fswitch^Iswitch)]; // get coded paternal marker gamete
	      prob_epsilon=pop->switch_table_prob[epsilon_gamete][Eqtl];  // Get prob of paternal gamete
	      if ( (Eqtl != 3) &&  (Bqtl !=3)) {
		   temp = prob_beta*prob_epsilon*std::exp(pm[wrow][wcol].k - sum1);
		   sum2 += temp;
		   m_weight[wrow][wcol]=temp;
		}
	      else {
		m_weight[wrow][wcol]=0.0; // because the means are only due the penetrance function
	      }
		wcol++;
	      }
	    }
	    wrow++;
	  }
	}

      nu = v = 0.0;
      wrow=0;
      for (m_q=0; m_q<tn_qtl; m_q++) {
	for (m_sw=0; m_sw < m_switches; m_sw++){   // loop for genotypes of mother
	  wcol=0;
	  for (f_q=0; f_q<tn_qtl; f_q++){
	    for (f_sw=0; f_sw < f_switches; f_sw++){    // loop for genotypes of father
	      m_weight[wrow][wcol] /= sum2;
	      nu += pm[wrow][wcol].nu*m_weight[wrow][wcol];
	      wcol++;
	    }
	  }
	  wrow++;
	}
      }

// Now the variance for mixture is calculated
      wrow=0;
      for (m_q=0; m_q<tn_qtl; m_q++){
	for (m_sw=0; m_sw < m_switches; m_sw++){   // loop for genotypes of mother
	  wcol=0;
	  for (f_q=0; f_q<tn_qtl; f_q++){
	    for (f_sw=0; f_sw < f_switches; f_sw++){    // loop for genotypes of father
	      v += pm[wrow][wcol].tsq*m_weight[wrow][wcol]+m_weight[wrow][wcol]*
		(pm[wrow][wcol].nu-nu)*(pm[wrow][wcol].nu-nu);
	      wcol++;
	    }
	  }
	  wrow++;
	}
      }

// k is the constant for a univariate normal (includes 1/sqrt(2*pi*sigma^2))

      k = std::log(sum2) + sum1 - std::log(std::sqrt(4.0*std::asin(1.0)*v));      // pi = 2.0*std::asin(1.0)

      temp_I->mix_anterior[i_q][i_sw].nu  = nu;
      temp_I->mix_anterior[i_q][i_sw].tsq = 1.0/v;
      temp_I->mix_anterior[i_q][i_sw].k   = k;
      //       std::cout << i_q << " " << i_sw << " sum2 " <<  log(sum2)   << " sum1 " << sum1 << " k " << k << std::endl;
    }
  }
}


void NuFamily::multi_m_terminal_peeling(const int iw){
   ////////////////////////////////////////////////////////////
   // REQUIREMENTS: initilization for anterior and posterior
   ////////////////////////////////////////////////////////////
   std::list<Individual *>::iterator pos;
   pos = connectors.begin();
   if (pos == connectors.end()) throw exception(" NuFamily::mutlti_m_terminal_peeling(): no connector");
   if (*pos == myfather) {
      multi_m_posterior(*pos, mymother);
      myfather->posterior_iw[mother_indx] = iw;
   }
   else if (*pos == mymother) {
      multi_m_posterior(*pos,myfather);
      mymother->posterior_iw[father_indx] = iw;
   }
   else {
      multi_m_anterior(*pos);
      (*pos)->anterior_iw = iw;
   }
}

double NuFamily::multi_m_log_likelihood(){

  unsigned j,gi;
  int f_q, f_sw, f_switches;
  double v, nu, y, v_y, nu_y, k_y;
  double a11,a12,a22,k,sum1=0.0,sum2=0.0;
  double inverse_sqrt2pi = 1.0/std::sqrt(4.0*std::asin(1.0));       // pi = 2.0*std::asin(1.0)
  f_switches=myfather->n_switches();
  if (!kernal) throw exception(" NuFamily::log_likelihood(): must be kernal nufamily");
  multi_m_posterior(myfather,mymother);

// Collect results to integrate uf for each f_q

  for (f_q=0; f_q < tn_qtl; f_q++){
// contributions from penetrance function of father
// skip if phenotype of father is missing

    if (!(myfather->record()[0].missing)) {
      y    = myfather->record()[0].double_val();
      v_y    = 1.0/(pop->residual_var[0][0]);
      k_y    = std::log(inverse_sqrt2pi*std::sqrt(v_y));
      nu_y   = pop->mean_for_genotype[f_q] - y;
    }
    else{
      nu_y = 0.0;
      v_y = 0.0;
      k_y = 0.0;
    }

    for (f_sw=0; f_sw < f_switches; f_sw++){   // loop for genotypes of J

      a11 = nu_y*nu_y*v_y;
      a12 = nu_y*v_y;
      a22 = v_y;
      k   = k_y;
// contributions from Father's anterior

      v   = myfather->mix_anterior[f_q][f_sw].tsq;
      nu  = myfather->mix_anterior[f_q][f_sw].nu;
      k   += myfather->mix_anterior[f_q][f_sw].k;
      a11 += nu*nu*v;
      a12 += -nu*v;
      a22 += v;

// Contributions from posteriors of father except mother.
// Skip over any posteriors not yet available
      unsigned numspouse = myfather->nspouse();
      for (j=0; j<numspouse; j++){
	if ((myfather->mix_posterior[j].done)) {
	  v     = myfather->mix_posterior[j].postvec[f_q][f_sw].tsq;
	  nu    = myfather->mix_posterior[j].postvec[f_q][f_sw].nu;
	  k    += myfather->mix_posterior[j].postvec[f_q][f_sw].k;
	  a11  += nu*nu*v;
	  a12  += -nu*v;
	  a22  += v;
	}
      }

// compute "log likelihood for gi"

    v          = 1.0/a22;
    m_weight[f_q][f_sw]   = -0.5*(a11 - a12*a12*v)
                          + k + std::log(std::sqrt(4.0*std::asin(1.0)*v));         // pi = 2.0*std::asin(1.0)
    sum1 += m_weight[f_q][f_sw];
    }
  }

// compute sum of the "likelihoods"

  sum1 /= (tn_qtl*f_switches);
  for (f_q=0;f_q<tn_qtl;f_q++){
     for (f_sw=0; f_sw < f_switches; f_sw++){   // loop for genotypes of J
       sum2 += (std::exp(m_weight[f_q][f_sw] - sum1));
     }
  }

  return (std::log(sum2) + sum1);
}

void NuFamily::multi_wksp_resize(int t_qtl, int max_switches){
  int i,j,k;

  NuFamily::child_matrix.resize(t_qtl*max_switches);
  NuFamily::pm.resize(t_qtl*max_switches);

  for (i=0; i< (t_qtl*max_switches); i++){
    NuFamily::child_matrix[i].resize(t_qtl*max_switches);
    NuFamily::pm[i].resize(t_qtl*max_switches);
  }

  NuFamily::wksp_for_gen.resize(t_qtl);
  NuFamily::spouse_matrix.resize(t_qtl);
  NuFamily::myfather_matrix.resize(t_qtl);
  NuFamily::mymother_matrix.resize(t_qtl);

  for (i=0; i < t_qtl; i++ ) {
    NuFamily::wksp_for_gen[i].resize(max_switches);
    NuFamily::spouse_matrix[i].resize(max_switches);
    NuFamily::myfather_matrix[i].resize(max_switches);
    NuFamily::mymother_matrix[i].resize(max_switches);
  }

  NuFamily::m_weight.resize(t_qtl*max_switches,t_qtl*max_switches);
  NuFamily::tr.resize(t_qtl,max_switches,0.0);
}

void NuFamily::pretend_multi_m_missing(int on) {
  //  unsigned tng = (unsigned) genotype_freq.size();
  // const double* gfreq = genotype_freq.ve;
  myfather->pretend_multi_m_missing(on,tn_qtl);
  mymother->pretend_multi_m_missing(on,tn_qtl);
  for (unsigned i=0; i<numoffs; i++) {
    myoffspring[i]->pretend_multi_m_missing(on,tn_qtl);
  }
}

void NuFamily::multi_m_ant_post() {

// indx[0] is the index of the mother in the father's spouse list
// indx[1] is the index of the father in the mother's spouse list

int spouse_index,i,i_peel=1, tag;
Individual* child;
int popsize=pop->size();

// calculate posterior for father through mother

  spouse_index = mother_indx;
  if ((!(myfather->base()) && (myfather->loop_connector)) || (myfather->id() > popsize)) {
    if (!(myfather->posterior_iw[spouse_index] == 2 ||
	  myfather->posterior_iw[spouse_index] == 3)) {
      multi_m_posterior(myfather ,mymother, i_peel);
    }
  }

// calculate posterior for mother through father

  spouse_index = father_indx;
  if ((!(mymother->base()) && (mymother->loop_connector)) || (mymother->id() > popsize)) {
    if (!(mymother->posterior_iw[spouse_index] == 2 ||
	  mymother->posterior_iw[spouse_index] == 3)) {
      multi_m_posterior(mymother ,myfather, i_peel);
    }
  }

// calculate anteriors for children

  for (i=0;i<numoffs; i++){
    child = myoffspring[i];
    tag=child->id();
    //   std::cout << "Extending child original " << tag << " ";
    if (tag > popsize) {
      child=pop->popmember[tag-popsize-1];
      tag -= popsize;
    }
    //   std::cout << "Extending child " << tag << std::endl;
    if (child->loop_connector) {
      if (!(child->anterior_iw == 2 || child->anterior_iw == 3)) {
	multi_m_anterior(child, i_peel);
      }
    }
  }
}

void NuFamily::multi_initialize(doubleMatrix &pen){
  // Initialize ALL members of the Nuclear family.
  // Else values are wrong for profile and maximization steps!
  // Would like to do just founders but cut individuals do not have this flag set.
  // WARNING I have reset flags may not be correct
  int nchild;
  Individual *I;

  //   if (myfather->base()) {

    myfather->initial_multi_anterior(pen);
    myfather->initial_multi_posterior(-1);
      myfather->posterior_iw[mother_indx] = 0;
  //      }

 //      if (mymother->base()) {

    mymother->initial_multi_anterior(pen);
    mymother->initial_multi_posterior(-1);
       mymother->posterior_iw[father_indx] = 0;
 //    }

  for (nchild=0; nchild < numoffs; nchild++) {
    I = myoffspring[nchild];
    I->initial_multi_anterior(pen);
    // I->initial_multi_posterior(-1); //BRS Child doesn't need posterior unless parent?
       I->anterior_iw = 0;
    // std::cout << I->id() << " initial ant scale " << I->m_anterior_scale << std::endl;
  }
}

void NuFamily::multi_m_initialize(int t_qtl){
  // Initialize ALL members of the Nuclear family.
  // Else values are wrong for profile and maximization steps!
  // Would like to do just founders but cut individuals do not have this flag set.

  int nchild;
  Individual *I;

  //   if (myfather->base()) {
  myfather->initial_multi_m_posterior(t_qtl,-1);
  myfather->initial_multi_m_anterior(t_qtl);
  myfather->posterior_iw[mother_indx] = 0;
  //      }

 //      if (mymother->base()) {
  mymother->initial_multi_m_posterior(t_qtl,-1);
  mymother->initial_multi_m_anterior(t_qtl);
  mymother->posterior_iw[father_indx] = 0;
 //    }

  for (nchild=0; nchild < numoffs; nchild++) {
    I = myoffspring[nchild];
    I->initial_multi_m_posterior(t_qtl,-1);
    I->initial_multi_m_anterior(t_qtl);
    I->anterior_iw = 0;
  }
}

//BRS

void NuFamily::eliminateGenotypes(unsigned locus){
	set<MaternalPaternalAlleles> momWorkSet, momSaveSet, dadWorkSet, dadSaveSet;
	set<MaternalPaternalAlleles> loopSet;
	SafeSTLVector<set<MaternalPaternalAlleles> > offspringWorkVector, offspringSaveVector;
	SafeSTLVector<set<MaternalPaternalAlleles> > intersectionVector;
	bool sampled        = mymother->genotNodeVector[locus].sampled;
	unsigned genotypeState = mymother->genotNodeVector[locus].genotypeState;
	if(!sampled){
		copy(mymother->genotNodeVector[locus].genotypeVector.begin(),
			 mymother->genotNodeVector[locus].genotypeVector.end(),
			 inserter(momWorkSet,momWorkSet.begin()) );
	} 
	else {
		momWorkSet.insert(mymother->genotNodeVector[locus].genotypeVector[genotypeState]); 
		mymother->genotNodeVector[locus].genotypeState = 0;
	}
	sampled        = myfather->genotNodeVector[locus].sampled;
	genotypeState = myfather->genotNodeVector[locus].genotypeState;
	if(!sampled){
		copy(myfather->genotNodeVector[locus].genotypeVector.begin(),
			 myfather->genotNodeVector[locus].genotypeVector.end(),
			 inserter(dadWorkSet,dadWorkSet.begin()) );
	}
	else{
		dadWorkSet.insert(myfather->genotNodeVector[locus].genotypeVector[genotypeState]); 
		myfather->genotNodeVector[locus].genotypeState = 0;
	}	
	offspringWorkVector.resize(numoffs);
	offspringSaveVector.resize(numoffs);
	intersectionVector.resize(numoffs);
	for (unsigned i=0;i<numoffs;i++){
		Individual *offsPointer = myoffspring[i];
		sampled        = offsPointer->genotNodeVector[locus].sampled;
		genotypeState = offsPointer->genotNodeVector[locus].genotypeState;
		if(!sampled){
			copy(offsPointer->genotNodeVector[locus].genotypeVector.begin(),
				 offsPointer->genotNodeVector[locus].genotypeVector.end(),
				 inserter(offspringWorkVector[i], offspringWorkVector[i].begin()) );
		}
		else {
			offspringWorkVector[i].insert(offsPointer->genotNodeVector[locus].genotypeVector[genotypeState]);
			offsPointer->genotNodeVector[locus].genotypeState = 0;
		}
	}
	set<MaternalPaternalAlleles>::iterator itMom, itDad;
	for (itMom = momWorkSet.begin(); itMom!=momWorkSet.end(); itMom++){
		for (itDad = dadWorkSet.begin(); itDad!=dadWorkSet.end(); itDad++){
			loopSet.clear();
			MaternalPaternalAlleles matPat;
			matPat.maternal = itMom->maternal;
			matPat.paternal = itDad->maternal;
			loopSet.insert(matPat);
			matPat.maternal = itMom->maternal;
			matPat.paternal = itDad->paternal;
			loopSet.insert(matPat);
			matPat.maternal = itMom->paternal;
			matPat.paternal = itDad->maternal;
			loopSet.insert(matPat);
			matPat.maternal = itMom->paternal;
			matPat.paternal = itDad->paternal;
			loopSet.insert(matPat);
			bool saveNot = false;
			for (unsigned i=0;i<numoffs;i++){
				intersectionVector[i].clear();
				set_intersection(loopSet.begin(),loopSet.end(),
								 offspringWorkVector[i].begin(),
								 offspringWorkVector[i].end(),
								 inserter(intersectionVector[i],intersectionVector[i].begin())  );
				if(intersectionVector[i].size()==0) {
					saveNot = true;
					break;
				}
			}
			if (saveNot) continue;
			momSaveSet.insert(*itMom);
			dadSaveSet.insert(*itDad);
			for (unsigned i=0;i<numoffs;i++){
				copy(intersectionVector[i].begin(), intersectionVector[i].end(), 
					 inserter(offspringSaveVector[i],offspringSaveVector[i].begin()) );
			}
		}
	}
	mymother->genotNodeVector[locus].genotypeVector.clear();
	if (momSaveSet.size()==0){
		std::cout <<"Incompatible genotypes in nuclear family\n";
		std::cout <<"Locus: " << locus << endl;
		std::cout <<"Father: " << myfather->myid << endl;
		std::cout <<"Mother: " << mymother->myid << endl;
		throw exception ("NuFamily::eliminateGenotypes: Incompatible genotypes in pedigree");
	}
	copy(momSaveSet.begin(),momSaveSet.end(), inserter(mymother->genotNodeVector[locus].genotypeVector,
													   mymother->genotNodeVector[locus].genotypeVector.begin()) );	
	myfather->genotNodeVector[locus].genotypeVector.clear();
	if (dadSaveSet.size()==0){
		std::cout <<"Incompatible genotypes in nuclear family\n";
		std::cout <<"Locus: " << locus << endl;
		std::cout <<"Father: " << myfather->myid << endl;
		std::cout <<"Mother: " << mymother->myid << endl;
		throw exception ("NuFamily::eliminateGenotypes: Incompatible genotypes in pedigree");
	}
	copy(dadSaveSet.begin(),dadSaveSet.end(), inserter(myfather->genotNodeVector[locus].genotypeVector,
													   myfather->genotNodeVector[locus].genotypeVector.begin()) );
	for (unsigned i=0;i<numoffs;i++){
		Individual *offsPointer = myoffspring[i];	
		offsPointer->genotNodeVector[locus].genotypeVector.clear();
		if (offspringSaveVector[i].size()==0){
			std::cout <<"Incompatible genotypes in nuclear family\n";
			std::cout <<"Locus: " << locus << endl;
			std::cout <<"Father: "    << myfather->myid    << endl;
			std::cout <<"Mother: "    << mymother->myid    << endl;
			std::cout <<"Offspring: " << offsPointer->myid << endl;
			throw exception ("NuFamily::eliminateGenotypes: Incompatible genotypes in pedigree");
		}
		copy(offspringSaveVector[i].begin(),offspringSaveVector[i].end(), 
			 inserter(offsPointer->genotNodeVector[locus].genotypeVector,
					  offsPointer->genotNodeVector[locus].genotypeVector.begin()) );	
	}
}



} ////////// end of namespace matvec
