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

#include <iomanip>

#include "model.h"
#include "population.h"

namespace matvec {

void Model::genotype_dist_peeling(const int prtl,const int estifreq)
{
   //////////////////////////////////////////////////////////////////////////
   // this works for arbitrary pedigrees
   // based on iterative peeling algorithm: Arendonk(1989), Fernando(1993)
   //  prtl     = 0, print nothing
   //             1, default, print genotype distribution and genotypic values
   //  estifreq = 0, not estimate allele frequencies
   //             1, default, estimate allele frequencies
   //////////////////////////////////////////////////////////////////////////
   if (type == bad_model) return;
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return;
      }
   }

   pop->residual_var = residual_var;
   pop->input_data(data);

   blupsol.resize(hmmesize);

   if (numterm) {
      setup_mme(&rellrhs);
      hmmec.solve(blupsol,rellrhs);
   }
   if (ntermGdist) {
      compute_xbzu(blupsol);                   // to adjust trait
      pop->genotype_dist_peeling(prtl,estifreq);
   }

   hmmec.resize(0,0);
   ntermGdist = 0;
}

Vector<double>* Model::genotypic_value_peeling(void)
{
   if (type == bad_model) return 0;
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return 0;
      }
   }

   blupsol.resize(hmmesize + ntermGdist*popsize);
   unsigned i;

   if (ntermGdist ) {
      pop->residual_var = residual_var;
      pop->input_data(data);
      if (numterm) {
         get_blupsol();
         compute_xbzu(blupsol);                   // to adjust trait
      }
      else {
         pop->genotype_dist_peeling(0);
         for (i=0; i<popsize; i++) blupsol[i] = pop->popmember[i]->est_GV;
         return  &(blupsol);
      }
   }
   else {
      setup_mme(&rellrhs);
      hmmec.solve(blupsol,rellrhs);
      return &(blupsol);
   }

   ///////////////////////////////////////////
   //   iterate begins now
   ////////////////////////////////////////////

   double * ve = &(blupsol[hmmesize]);
   unsigned niterate = 0;
   while (niterate++ < 4) {
      pop->genotype_dist_peeling(0);
      for (i=0; i<popsize; i++) ve[i] = pop->popmember[i]->est_GV;
      blupsol.print();
      setup_mme(&rellrhs);
      hmmec.solve(blupsol,rellrhs);
   }
   return &(blupsol);
}

double Model::estimate_peeling(void)
  // double Model::ml_estimate_peeling(void) //old version
{
   if (type == bad_model) throw exception("Model::estimate_peeling(): bad model");
   double max_llhood = 0.0;
   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return max_llhood;
      }
   }

   minfun_indx = 2;
   pop->residual_var = residual_var;
   pop->input_data(data);
   pop->build_pop_gamete();
   unsigned i,n,m,k,j;
   m = pop->tn_gamete - 1;          // because p(A1)+ p(A2) +...+ p(At) = 1
   n = 1 + m;                       // n = residual_var + allele_freq
//   n += pop->tn_genotype;         //     + genotypic_value // uncommented in old version
   Vector<double> x(n);

   x[0] = pop->residual_var[0][0];           // residual_variance

   ChromStruct *chrom = pop->prior->chrom();
   double* ve = chrom[0].locus[0].allele_freq.begin();
   for (i=0; i<m; i++) x[i+1] = ve[i];         // 1 is for residual_variance

   // uncommented in old version
//   k = 1 + m;
//   const double** me = pop->prior->genotypic_val(0,0);
//   for (i=0; i<tn_genotype; i++) for (j=0; j<=i; j++) x[k++] = me[i][j];
   int iter = 0;
   max_llhood = praxis(x,n,iter);
   minfun_indx = 0;
   return max_llhood;
}

double Model::minfun_peeling(const Vector<double> &x, const int n)
{
   // requirement: build_pop_gamete() must have been called
   int i,m = pop->tn_gamete - 1;      // because p(A1)+ p(A2) +...+ p(At) = 1
   double p,pm;

   //////////////////////////////////////////////////////////////
   // constraints for minimization
   //       (0)  residual_variance >= 0;
   //       (1)  0.0 <= p(Ai) <= 1
   //////////////////////////////////////////////////////////////
   if (x[0] < 0.0) return 1.0e+25;
   for (pm=1.0, i=0; i<m; i++) {
      pm -= p = x[1+i];
      if (p>1.0 || p<0.0) return 1.0e+25;
   }
   if (pm>1.0 || pm < 0.0 ) return 1.0e+25;

   //////////////////////////////////////////////////////////////
   // now put back new values in x into appropriate places
   //////////////////////////////////////////////////////////////

//RLF
   /*

   pop->residual_var.me[0][0] = x[0];
   for (i=0; i<m; i++) pop->pop_gamete[0].chromosome[i].freq(x[i+1]);
   pop->pop_gamete[0].chromosome[m].freq(pm);
//   pop->prior->genotypic_val(1,1,&x[1+m]);  // get genotypic_values from x
   */
   pop->residual_var[0][0] = x[0];
   pop->pop_gamete[0].chromosome[0].freq(x[1]);
   pop->pop_gamete[0].chromosome[1].freq(x[2]);
   pop->pop_gamete[0].chromosome[2].freq(pm);
   pop->mean_for_genotype[0]  = x[3];
   pop->mean_for_genotype[1]  = x[4];
   pop->mean_for_genotype[2]  = x[5];
    pop->tn_qtl=4;
//RLF
   double log_llhood = -1.0*pop->log_likelihood_peeling();
// std::cout << "log likelihood value (its negative) = " << log_llhood << std::endl;
   return log_llhood;
}

double Model::log_likelihood_peeling(const unsigned maxiter)
{
   if (type == bad_model) throw exception("Model::log_likelihood_peeling(): bad model");
   if (ntermGdist == 0) {
      return restricted_log_likelihood();
   }

   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return 0.0;
      }
   }

   int pt,i;
   Individual *I;
   double sigma_mu,retval = 0.0;
   double *sol, rr;
   unsigned nl,rank_mme,startaddr,iter = 0;
   pop->analysis_type = "not-multi";
  pop->tn_qtl=4;
   if (numterm) {
      retval += numobs * std::log(residual_var[0][0]);
      hmmec.resize(hmmesize,max_nz);
      rellrhs.resize(hmmesize);
      blupsol.resize(hmmesize);

      pop->genotype_dist_peeling(0);
      for (i=0; i<popsize; i++) {
         I = pop->popmember[i];
         I->xbzu_val = I->est_GV;
      }
      unsigned nvc = totalnvc() - 1;              // number of ratio's
      Vector<double> zz(popsize);
      Vector<double> ratio(nvc);
      Vector<double> vc(nvc+1);

      var2vec(vc);
      for (i=0; i<nvc; i++) ratio[i] = vc[nvc]/vc[i];
      for (pt=-1,i=0; i<numterm; i++) {
         if (term[i].classi() == 'P') {
            pt = i;
            startaddr =  term[pt].start + 1;
            nl = term[pt].nlevel();
            break;
         }
      }

      setup_ww_single_trait(ratio,pt,startaddr,&zz);
      hmmec.reorder();
      rank_mme = hmmec.factorization(5);

      hmmec.solve(blupsol,rellrhs,"ysmp1");
      retval += hmmec.info_vec[0]- rank_mme * std::log(residual_var[0][0]);

      if (pt >= 0) {
         sol = &(blupsol[startaddr-1]);
         for (rr=0.0,i=0; i<nl; i++) {
            rr += (*sol * zz[i] * *sol);  sol++;
         }
         sigma_mu = *(term[pt].prior->var_matrix())[0][0];
         retval += nl * std::log(sigma_mu);
         sigma_mu = hmmec.q(blupsol.begin(),blupsol.begin(),startaddr,startaddr+nl-1)
                       - rr;
         retval += std::log(sigma_mu);
      }
      compute_xbzu(blupsol);                   // to adjust trait
   }
   retval *= -0.5;
   pop->residual_var = residual_var;
   pop->input_data(data);
   retval += pop->log_likelihood_peeling(maxiter);
   return retval;
}

//BRS

void Model::multipoint_setup()
{
  // This computes the multipoint likelihood, first we set up all the necessary parameters and space required except for FPMM

  Q_pos=0;  // Set initial QTL parameter to be in the first position
  new_distance=0;
  N_Parm_List=0;
  pop->tn_genotype=3;
  pop->tn_qtl=4;
  pop->residual_var = residual_var;
  pop->input_markerData(data);
  //  std::cout << "switches \n" << flush;
  pop->set_switches();
  std::cout << "tables \n";
  multipoint_tables();

  RecoTable.clear(); // need to put somewhere!

  multipoint_Rec_table();
  multipoint_QTL_table();
  pop->F=0;
  std::cout << " Multipoint setup done " << std::endl;
  //  std::cout << RecoTable << std::endl;
}

void Model::multipoint_setup(Fpmm& F, char map_fun, int map_par)
{
  // This computes the multipoint likelihood, first we set up all the necessary parameters and space required

  Q_pos=0;  // Set initial QTL parameter to be in the first position
  new_distance=0;
  N_Parm_List=0;
  pop->tn_genotype=3;
  pop->tn_qtl=4;
  pop->residual_var = residual_var;
  pop->input_markerData(data);
  //  std::cout << "switches \n" << flush;
  pop->set_switches();
  std::cout << "tables \n";
  multipoint_tables();
  RecoTable.clear(); // need to put somewhere!
  map_function=map_fun;
  map_parm=map_par;
  multipoint_Rec_table();
  multipoint_QTL_table();
  pop->multipoint_init_parm(F);
  std::cout << " Multipoint setup done " << std::endl;
  //  std::cout << RecoTable << std::endl;
}

double Model::multipoint(int maxiter) {
  double llh;

  pop->analysis_type = "multipoint";
  // This is a bit excessive as means we have two versions of Farg when doing FPMM
  Fpmm Farg(0,0.5,1.0); // setup for 0 polygenic loci.
  if (pop->F == 0) {
    // Need to setup F as not using FPMM stuff
    pop->multipoint_init_parm(Farg);
  }

  llh=pop->multipoint_likelihood(maxiter);
  //  std::cout << "QTL pos " << Q_pos << std::endl;
  std::cout << "Natural log llh= " << std::setprecision(10) << llh << " Log10 llh = " << std::setprecision(10) << (std::log10(std::exp(1.0))*llh) << std::endl;
  /*
    multipoint_Rec_recalc(6);
      multipoint_QTL_table();
 llh=pop->multipoint_likelihood(maxit);
 //  std::cout << "QTL pos " << Q_pos << std::endl;
  std::cout << "Natural log llh= " << std::setprecision(10) << llh << " Log10 llh = " << std::setprecision(10) << (log10(exp(1.0))*llh) << std::endl;
  */
return llh;
}
void Model::multipoint_genot_probs()
{
  int i;
  Individual *I;
  std::cout << "computing genot_probs" << std::endl;
  for (i=0;i<pop->popsize; i++){
    I = pop->popmember[i];
    //I->initial_anterior(genotype_freq.ve,tn_genotype);
    I->initial_posterior(3);
  }
  pop->analysis_type = "multipoint";
  pop->build_nufamily();
  pop->multi_geno_dist_peeling(0);
}

void Model::multipoint_tables()
{
/* Creates the table for all possible switch sets for a given number of loci:
     Nloci is the number of loci involved
     ST is a 'matrix' whose row dimension is the number of switch sets,
          the first column in each row is the number of switches contained in that set
          the next columns are the actual switches possible in that set.

This function first determines which set is being repeated and copies it to the new set as well as
adding the new items which are the addition of the old switch + 2 raised to the power of the locus-1.

*/
  int current, row, EL, TE,i,j,Nel,pow2,pow3;
  double diff;

  int eindex,sindex,etemp,stemp,eval,sval,locus;
  double epow,spow,count=0.0;

  pow2= (int) (1 << pop->n_markerLoci); // pow(2,pop->n_markerLoci));
  pow3= (int) (pow((double) 3,(int) pop->n_markerLoci));

  Matrix<int> ST(pow2,pow2+1);
  Matrix<int> MT(pow3,pow2);

  ST[0][0]=1;
  ST[1][0]=2;
  ST[1][2]=1;

  current=1;
  for (i=2; i <= pop->n_markerLoci; i++)
    {
      diff=(1 << i)-(1 << i-1); // // pow(2,i)-pow(2,i-1);
      for (j=1; j <= diff; j++)
	{
	  current++;
	  row=(int) (current-diff);
	  TE=ST[row][0];
	  Nel=0;
	  for ( EL=1 ; EL <= TE; EL++)
	    {
	      ST[current][EL]=ST[row][EL];
	      ST[current][(TE+EL)]=ST[row][EL]+ ((int) (1 << i-1)); // pow(2,i-1));
	      Nel += 2;
	    }
	  ST[current][0]=Nel;
        }
    }

  pop->switch_table = ST;

  //  std::cout << "ST \n" << ST;
  for(eindex=0; eindex < pow3; eindex++)
    {
      for(sindex=0; sindex < pow2; sindex++)
	{
	  etemp=eindex;
	  stemp=sindex;
	  count =0.0;
	  for (locus=pop->n_markerLoci-1; locus > -1; locus--)
	    {
	      epow=pow((double) 3,(int) locus);
	      eval=(int) (etemp/epow);
	      etemp -= ((int) (epow*eval));
	      spow=(1 << locus); // pow(2,locus);
	      sval=(int) (stemp/spow);
	      stemp -= ((int) (spow*sval));
	      if (eval == 2 )
		count += (2*epow);
	      else if (eval != sval)
		count += epow;
	    }
	  MT[eindex][sindex]= (int) count;
	}
    }
  //  std::cout << "MT \n" << MT;
  pop->switch_table_gmt = MT;
}

void Model::multipoint_QTL_table() {
/* Try to create  the marker and qtl gametic event */

 int current, row, EL, TE,i,j,Nel,pow2,pow3;
  int eindex,sindex,etemp,stemp,eval,sval,locus;
  double epow,spow,count=0.0;
  double diff;

  pow2= (int) (1 << pop->n_markerLoci); // pow(2,pop->n_markerLoci));
  pow3= (int) (pow((double) 3,(int)pop->n_markerLoci));

  int qindex;
  doubleMatrix QT(pow3,3);
  Vector<int> gamete((pop->n_markerLoci)+1);

  for(eindex=0; eindex < pow3; eindex++) {
    etemp=eindex;
    for (locus=(pop->n_markerLoci-1); locus > -1; locus--){
      epow=pow((double) 3,(int) locus);
      eval=(int) (etemp/epow);
      etemp -= ((int) (epow*eval));
      if (((pop->n_markerLoci) - 1 - locus) < (Q_pos))
 	gamete[(pop->n_markerLoci)-1-locus] = eval;
      else
	gamete[(pop->n_markerLoci)-locus]=eval;
      //  std::cout  << eindex << " " << eval << " ";
    }
    // std::cout << std::endl;
    for (qindex=0; qindex < 3; qindex++){
      gamete[Q_pos]=qindex;
      /*
      for (i=0; i < (pop->n_markerLoci)+1 ; i++) {
	std::cout << " " << gamete[i];
      }
      std::cout << std::endl;
      */
      QT[eindex][qindex]=ProbGamete(gamete,((pop->n_markerLoci)+1));
      if (QT[eindex][qindex] < 0) {
	std::cout << "Pg " << ProbGamete(gamete,((pop->n_markerLoci)+1)) << std::endl;
	std::cout << RecoTable << std::endl;
	exit(1);
      }
    }
  }
  //    std::cout << " QT \n" << QT << std::endl;
  pop->switch_table_prob= QT;

}

double Model::ProbGamete (Vector<int> value, int nloci)
{
  int Flocus=0,Slocus;
  double Prob=1.0;


/* Find first element in the gamete that is not 2 */
  while ((Flocus < nloci) && (value[Flocus] == 2)) {
    Flocus++;
  }
/* If all loci are 2 then at this point Flocus is nloci+1 and we end here else we do the IF statement */

  if (Flocus <  nloci) {
    Prob *= 0.5;
    /* Now find the next locus (Slocus)  those gamete is not 2. */
    for (Slocus=Flocus+1; Slocus <  nloci; Slocus++) {
      if (value[Slocus] != 2) {
	if (value[Flocus] == value[Slocus])
	  Prob *= (RecoTable[Slocus][Flocus]);  /* Non recombinant */
	else
	  Prob *= (RecoTable[Flocus][Slocus]);  /* Recombinant */
	Flocus=Slocus;
      }
    }
  }
  return Prob;
}

/*
double Model::MapD(double rate) {
  if (rate < 0.5)
    return (-0.5 * std::log(1-2*rate));
  else
    return 5.0;
}

double Model::MapR(double dist) {
  return (0.5*(1.0 - std::exp(-2.0*dist)));
}
*/
double Model::MapF(double theta, int qr) {
    // Computes different map functions that are multilocus feasible:
    // Model::map_method M is Morgans
    // Model::map_method H (default) is Haldanes
    // Model::map_method B Binomial with extra_parm
    // qr is a flag
    //    0 means theta is a recombination rate
    //    1 means theta is a distance in Morgans no centiMorgrans

   double temp;
   if ((map_function=='M') || (map_function=='m')) {
     //Morgan's function  recombination rate equals distance
     return theta;
    }
   else if ((map_function=='B') || (map_function=='b')) {
     // Binomial Method
        if (qr) {
            temp=std::pow((1.0-2.0*theta),(1.0/map_parm));
            return (0.5*map_parm*(1.0-temp));
            }
        else {
            if (theta < map_parm/2.0) {
               temp=(1.0-2.0*theta/map_parm);
               return 0.5*(1.0-(std::pow(temp,map_parm)));
               }
            else
               return 0.5;
            }
        }
   else  {
     // Haldane's function if method not set
        if (qr) {
           if (theta < 0.5) {
              temp= (-0.5 * std::log(1.0-2.0*theta));
              return temp;
              }
           else
              return 5.0;
            }
        else {
            temp= (0.5*(1.0 - std::exp(-2.0*theta)));
            return temp;
           }
     }
}

void Model::multipoint_Rec_table() {
  int i,j;
  double reco=0.1;
  int nloci=pop->n_markerLoci+1;

  //  RecoTable.me=NULL;  // There is an error here as RecoTable should be NULL but is not
  // Hopefully fixed it by including it Model::initalize
  RecoTable.resize(nloci, nloci);
  RecoTable.assign(reco);

  for (i=0; i< nloci; i++) {
    // we are adding 2.0 to distance so we can place QTL
    // to the left of the leftmost marker
     RecoTable[i][i]=pop->prior->get_distance(1,i+1) + 2.0;
  }
  new_distance=pop->prior->get_distance(1,1) + 2.0; // Store current qtl distance until changed.
  multipoint_compute_Rec_table();
  //   std::cout << RecoTable << std::endl;

}

void Model::multipoint_compute_Rec_table() {
  int nloci=pop->n_markerLoci+1, i,j;
  double temp;
  int Q_start = Q_pos;

  for(i=Q_start; i<nloci-1; i++) {
    if(RecoTable[i][i] > RecoTable[i+1][i+1]) {
      temp = RecoTable[i][i];
      RecoTable[i][i] = RecoTable[i+1][i+1];
      RecoTable[i+1][i+1] = temp;
      Q_pos++;
    }
    else {
      break;
    }
  }

  for(i=Q_start; i>0; i--) {
    if(RecoTable[i][i] < RecoTable[i-1][i-1]) {
      temp = RecoTable[i][i];
      RecoTable[i][i] = RecoTable[i-1][i-1];
      RecoTable[i-1][i-1] = temp;
      Q_pos--;
    }
    else {
      break;
    }
  }

  for (i=0; i < nloci-1; i++){
    for (j=i+1; j < nloci; j++){
      RecoTable[i][j]=MapF(RecoTable[j][j]-RecoTable[i][i]);
      RecoTable[j][i]=1.0-RecoTable[i][j];
      //  std::cout << "i " << i << " j " << j << " Reco " << RecoTable[i][j] << " RT[j] " << RecoTable[j][j] << " RT[i] " << RecoTable[i][i] << std::endl;
    }
  }
  //  std::cout << RecoTable;
}



void Model::multipoint_Rec_recalc(double distance) {
  RecoTable[Q_pos][Q_pos]=distance;
  multipoint_compute_Rec_table();
  //    std::cout << RecoTable << std::endl;
}

void Model::multipoint_profile_distance(int method, int maxiter, double Min, double Max, double step_size, int Print, const std::string &fname) {

  int i, nsteps;
  double distance=0.0, T_dist, temp=0.0;
  doubleMatrix Results;

  std::ofstream Outfile (fname.c_str(), std::ios::app);
    if (!Outfile) throw exception("Cannot open requested filename");
    
  Fpmm Farg(0,0.5,1.0); // setup for 0 polygenic loci.
  if (pop->F == 0) {
    // Need to setup F as not using FPMM stuff
    pop->multipoint_init_parm(Farg);
  }
  std::cout << "Nloci " << pop->n_markerLoci << std::endl;
  T_dist = Max-Min;
  // T_dist = 2.0+RecoTable[pop->n_markerLoci][pop->n_markerLoci];
   temp= (T_dist/step_size);
  nsteps= (int) (std::ceil(temp+1.0));
  std::cout  << "Total dist " << T_dist << " With steps " << nsteps << " Of size " << step_size << std::endl;
  Results.resize(nsteps,2);
  distance = Min -step_size + 2.0;
  //  std::cout << distance << std::endl;

  if (method==0) {
   pop->analysis_type = "multipoint";
       if (pop->P_ndim > 1) { 
       	    std::cout << "Method used : multipoint FPMM" << std::endl;
        }
       else {
   	    std::cout << "Method used : Simple multipoint" << std::endl;
       }
       for (i=0; i < nsteps; i++) { 
       	    distance += (step_size);
            Results[i][0] = (distance-2.0);
       // Find the new location and position of  the QTL and recompute QTL table
            multipoint_Rec_recalc(distance);
            multipoint_QTL_table();
            temp=(pop->multipoint_likelihood(maxiter));
            Results[i][1]=temp;
        } 

       if (Print==1) {
          for (i=0; i < nsteps; i++) { 
  	       Outfile <<  std::setw(10) << std::setprecision(6) << Results[i][0];
      	       Outfile << " " << std::setw(20) << std::setprecision(15);
       	       Outfile <<  Results[i][1] << std::endl; 
          } 
     } 
  } 
  else if (method==1) { 
 	pop->analysis_type = "multipoint_m";
        std::cout << "Method used : multipoint Mixture approx" << std::endl;
        for (i=0; i < nsteps; i++) { 
            distance += (step_size);
	    Results[i][0] = (distance-2.0);
// Find the new location and position of  the QTL and recompute QTL table
            multipoint_Rec_recalc(distance);
	    multipoint_QTL_table();
	    temp = pop->multi_m_log_likelihood_peeling(maxiter);
            Results[i][1]=temp; 
	} 
   if (Print) { 
       for (i=0; i < nsteps; i++) {
 	       Outfile <<  std::setw(10) << std::setprecision(6) << Results[i][0];
       	       Outfile << " " << std::setw(20) << std::setprecision(15);
       	       Outfile <<  Results[i][1] << std::endl; 
       	    } 
   } 
  }
  else {
 	std::cout << "Unknown method for multipoint_profile_distance" << std::endl;
  }
				 

  if (Print==0) {
    std::cout  << Results << std::endl;
  }
  Outfile.close();
}

double Model::multipoint_ml_estimate(const int method, int niter, double tol, double epsilon, int printlevel)
{
      /* 
       * tol  is the tolerance allowed for the precision of the
       *      solution. praxis returns if the criterion
       *      2 * ||x[k]-x[k-1]|| <= sqrt(macheps) * ||x[k]|| + tol
       *      is fulfilled more than ktm times.
       *      the default value 1.0e-16
       *      tol is usually equal to  EPISILON*EPISILON
       *
       * printlevel controls praxis printing output:
       *      0 -> no output
       *      1 -> print only starting and final values
       *      2 -> detailed map of the minimization process
       *      3 -> print also eigenvalues and vectors of the
       *           search directions
       *      the default value is 0
       *                                                  */
   if (type == bad_model) throw exception("Model::ml_estimate_peeling(): bad model");
  // Setting up to run praxis to maximise the llh
  if (method == 0)
      pop->analysis_type = "multipoint";
 else
    pop->analysis_type = "multipoint_m";
  maxit=niter;
  unsigned n;
   double max_llhood = 0.0;

   if (!data_prepared) {
      if (!prepare_data()) {
         type = bad_model;
         return max_llhood;
      }
   }

   minfun_indx = 3;
   //  std::cout << "get x " << N_Parm_List << std::endl;
   Vector<double> x(N_Parm_List);
   for (n=0; n < N_Parm_List; n++)
     x[n]= *(Parm_List[n].Value);
   int iter = 0;
   //  std::cout << "Go to praxis " << std::endl;
   max_llhood = praxis(x,n,iter, tol, epsilon, printlevel);
   minfun_indx = 0;

   return max_llhood;
}

double Model::minfun_multipoint(const Vector<double> &x, const int n) {
  int i;
  // std::cout << "entered minfun_multipoint " << std::endl;
  // setup values from praxis into correct places
  for (i=0; i < N_Parm_List; i++) {
    *(Parm_List[i].Value) = x[i];
    //    std::cout << " Parameter " << i << " is  " << x[i] << std::endl;
    if ((x[i] > Parm_List[i].Max) || (x[i] < Parm_List[i].Min)) {
      std::cout << " Parameter " << i << " is out of bounds  " << x[i] << " Max is " << Parm_List[i].Max;
      std::cout << " Min is " << Parm_List[i].Min << std::endl;
      return 1.0e+25;
    }
  }
    // First need to find the new location and position of  the QTL and recompute QTL table
  // std::cout << "New dist " << new_distance; // << "RecoTable Before change " << RecoTable;
    multipoint_Rec_recalc(new_distance);
    //  std::cout << "After change " << RecoTable;
    //   std::cout << " POP Q FREQS " << pop->Q_freq;
    // Recompute QTL genotype frequencies
    double QTL_freq= pop->prior->chrom()[0].locus[0].allele_freq[0];
    pop->Q_freq(1) = QTL_freq*QTL_freq;
    pop->Q_freq(2) = QTL_freq*(1-QTL_freq);
    pop->Q_freq(3) = QTL_freq*(1-QTL_freq);
    pop->Q_freq(4) = (1-QTL_freq)*(1-QTL_freq);
  // End of old Stuff

  multipoint_QTL_table();
int N_means=5;
      // Find which means was used
 for (i=1; i <= N_Parm_List; i++) {
      // std::cout << "Start:" <<Parm_List[i].Name << ":End" << endl;
     if (Parm_List[i].Name == "Mean_Qq   ") {
        N_means= N_means-2;
       }
     if (Parm_List[i].Name == "Dom_Qq    ") {
       N_means--;
      }
      if (Parm_List[i].Name == "Mean_QQ   ") {
         N_means--;
        }
}
     //Dominance Means
      if (N_means == 3) {
          if (pop->mean_for_genotype[0] > pop->mean_for_genotype[3]){
	     if (pop->mean_for_genotype[1] > pop->mean_for_genotype[0]){
	         pop->mean_for_genotype[1] = pop->mean_for_genotype[0];
	         }
	     else if (pop->mean_for_genotype[3] > pop->mean_for_genotype[1]){
	         pop->mean_for_genotype[1] = pop->mean_for_genotype[3];
	         }
	   }
	   else {
	       if (pop->mean_for_genotype[1] > pop->mean_for_genotype[3]){
	           pop->mean_for_genotype[1] = pop->mean_for_genotype[3];
	           }
	        else if (pop->mean_for_genotype[0] > pop->mean_for_genotype[1]){
	           pop->mean_for_genotype[1] = pop->mean_for_genotype[0];
	           }
	    }
	}
	//Additive means
      else if (N_means == 4) {
	     (pop->mean_for_genotype[1]) = 0.5*((pop->mean_for_genotype[0])+ (pop->mean_for_genotype[3]));
	    }
	//Same means 
      else if (N_means == 5) {
             (pop->mean_for_genotype[1]) = (pop->mean_for_genotype[0]);
             (pop->mean_for_genotype[3]) = (pop->mean_for_genotype[0]);
      }
// General model don't do anything!

      //Force heterozygotes to have same mean
        (pop->mean_for_genotype[2]) = (pop->mean_for_genotype[1]);


      //            cout << "mean 0 = " << pop->mean_for_genotype[0] << endl;
      //            cout << "mean 1 = " << pop->mean_for_genotype[1] << endl;
      //            cout << "mean 2 = " << pop->mean_for_genotype[2] << endl;
      //            cout << "mean 3 = " << pop->mean_for_genotype[3] << endl;
	 
  //  std::cout << "parameters set " << std::endl;
  double log_llhood;
  if (pop->analysis_type == "multipoint") {
//	std::cout << "doing " << pop->analysis_type << " multipoint" << std::endl;
 log_llhood= -1.0*pop->multipoint_likelihood(maxit);
}
   else {
//	std::cout << "doing " << pop->analysis_type << " multipoint" << std::endl;
 log_llhood= -1.0*pop->multi_m_log_likelihood_peeling(maxit);
}
//  std::cout << "log likelihood value (its negative) = " << std::setprecision(14);
//  std::cout << log_llhood;
//  for (i=0; i < N_Parm_List; i++)
//    std::cout << "(" <<  i  << ")" << " = " << x[i] << "  ";
//  exit(1);
  return log_llhood;
}


void Model::multipoint_estimate_Ve(double Initial, double Max, double Min) {
  double tmp;
   if (N_Parm_List == 0) {
    Parm_List.resize(7);
  }
  if (Min < 0.0)
    Min=0.0;
  if (Min > Max) {
    tmp=Min;
    Min=Max;
    Max=tmp;
  }
  if (Initial > Max)
    Initial=Max;
  if (Initial < Min)
    Initial=Min;
  pop->residual_var.assign(Initial);
  N_Parm_List++;

  Parm_List(N_Parm_List).Value = &(pop->residual_var(1,1));
  Parm_List(N_Parm_List).Max = Max;
  Parm_List(N_Parm_List).Min = Min;
  Parm_List(N_Parm_List).Name = "Res_Var   ";
  //  std::cout << " value " <<  *Parm_List(N_Parm_List).Value << " Min " <<  Parm_List(N_Parm_List).Min << " Max " <<  Parm_List(N_Parm_List).Max << " name " <<  Parm_List(N_Parm_List).Name << std::endl;
}

void Model::multipoint_estimate_Distance(double Initial, double Max, double Min) {
  double T_dist,tmp;
  if (N_Parm_List == 0) {
    Parm_List.resize(7);
  }
  if (Min > Max) {
    tmp=Min;
    Min=Max;
    Max=tmp;
  }
  Initial +=2.0;
  Max +=2.0;
  Min +=2.0;
  T_dist = 2.0+RecoTable[pop->n_markerLoci][pop->n_markerLoci];
  // CHECK parameter bounds;
  if (Max > T_dist)
    Max=T_dist;
  if (Min < 0.0)
    Min=0.0;
  if (Initial > Max)
    Initial=Max;
  if (Initial < Min)
    Initial=Min;
  new_distance=Initial;
  N_Parm_List++;
  Parm_List(N_Parm_List).Value = &new_distance;
  Parm_List(N_Parm_List).Max = Max;
  Parm_List(N_Parm_List).Min = Min;
  Parm_List(N_Parm_List).Name = "Distance  ";
  // Need to recompute tables for the new starting location
  multipoint_Rec_recalc(new_distance);
  multipoint_QTL_table();
  // std::cout << T_dist << " min " << Min << " Initial " << Initial << " Max " << Max << std::endl;
  //  std::cout << RecoTable ;
  // exit(1);
  }

void Model::multipoint_estimate_Qmeans(Vector<double> Initial, Vector<double> Max, Vector<double> Min, char mtype) {
     // This maximizes the genotypic means for the genotypes qq, Qq and QQ
     // only two alleles at QTL supported
     // Currently force same heterozygotes mean_Qq is set the same as mean_qQ
     // Future flag maybe added to allow different heterozygote means
     // Current types
     // S is same mean : mean qq == mean Qq == mean QQ
     // A is additive model : mean Qq == average of (mean qq + mean QQ)
     // D is dominance model: mean qq <= mean Qq <= mean QQ
     // G or other value has unrestricted means so that overdominance is possible

double tmp;
int i, j=0;
if (N_Parm_List == 0) {
	    Parm_List.resize(7); // Parm_list is empty therefore resize it
	        }
// Need to check if any other type of QTL mean maximization have been defined
 for (i=1; i <= N_Parm_List; i++) {
 //      cout << "Name " << Parm_List(i).Name << endl;
     if (Parm_List(i).Name == "Mean_qq   ") {
        j=1;
        }
  }
  if (j) {
	  std::cout << "Another estimate Qmeans has been previously defined" << std::endl;
     throw exception("Model::estimate_Qmeans: has been previously defined");
     }
  else {
     N_Parm_List++;
     // For qq genotype
     pop->mean_for_genotype[0] =  Initial[0];
     Parm_List(N_Parm_List).Max = Max[0];
     Parm_List(N_Parm_List).Min = Min[0];
     Parm_List(N_Parm_List).Value = &(pop->mean_for_genotype[0]);
     Parm_List(N_Parm_List).Name = "Mean_qq   ";
     if ((mtype == 'S') || (mtype=='s')){ // Same means
         (pop->mean_for_genotype[1]) = (pop->mean_for_genotype[0]);
         (pop->mean_for_genotype[2]) = (pop->mean_for_genotype[0]);
	 (pop->mean_for_genotype[3]) = (pop->mean_for_genotype[0]);
	 }
     else {
     // For QQ genotype
          N_Parm_List++;
          pop->mean_for_genotype[3] =  Initial[2];
	  Parm_List(N_Parm_List).Max = Max[2];
	  Parm_List(N_Parm_List).Min = Min[2];
	  Parm_List(N_Parm_List).Value = &(pop->mean_for_genotype[3]);
	  Parm_List(N_Parm_List).Name = "Mean_QQ   ";
	  if ((mtype=='A') || (mtype=='a') ) { // Additive
             (pop->mean_for_genotype[1]) = 0.5*((pop->mean_for_genotype[0])+(pop->mean_for_genotype[3]));
	     (pop->mean_for_genotype[2]) = (pop->mean_for_genotype[1]);
	     }
	  else { // set heterozygotes to the same value depending on model
	     N_Parm_List++;
	     if ((mtype == 'D')|| (mtype == 'd')) {// Dominance model
	        Parm_List(N_Parm_List).Name = "Dom_Qq    ";
                // Now check that the Qq means fall within the bounds of QQ and qq means
                tmp= (double) std::max(Max[0],Max[2]);
		if (Max[1] > tmp)
		    Max[1]=tmp;
		tmp=std::min(Min[0],Min[2]);
		if (Min[1] < tmp)
		   Min[1]=tmp;
		tmp= (double) std::max(Initial[0],Initial[2]);
		if (Initial[1] > tmp)
		    Initial[1]=tmp;
		tmp=std::min(Initial[0],Initial[2]);
		if (Initial[1] < tmp)
		    Initial[1]=tmp;
		pop->mean_for_genotype[1] =  Initial[1];
		Parm_List(N_Parm_List).Max = Max[1];
		Parm_List(N_Parm_List).Min = Min[1];
		Parm_List(N_Parm_List).Value = &(pop->mean_for_genotype[1]);
		(pop->mean_for_genotype[2]) = (pop->mean_for_genotype[1]);
		}
	   else {  // General model
              Parm_List(N_Parm_List).Name = "Mean_Qq   ";
              pop->mean_for_genotype[1] =  Initial[1];
	      Parm_List(N_Parm_List).Max = Max[1];
	      Parm_List(N_Parm_List).Min = Min[1];
	      Parm_List(N_Parm_List).Value = &(pop->mean_for_genotype[1]);
	      (pop->mean_for_genotype[2]) = (pop->mean_for_genotype[1]);
	      }
	   pop->mean_for_genotype[2] =  Initial[1];
	  (pop->mean_for_genotype[2]) = (pop->mean_for_genotype[1]);
	  }

      }
 }
  /*
   * for (i=1; i <= N_Parm_List; i++) {
   *          cout << "Name " << Parm_List(i).Name <<  " Value= " << &Parm_List(i).Value
   *          << std::endl;
   *          }
   *          */
}


void Model::multipoint_estimate_Qfreq(double Initial, double Max, double Min)
{
  double QTL_allele_freq =Initial,tmp;
  if (N_Parm_List == 0) {
    Parm_List.resize(7);
  }
  if (Max > 1.0)
    Max=1.0;
  if (Min < 0.0)
    Min=0.0;
  if (Min > Max) {
    tmp=Min;
    Min=Max;
    Max=tmp;
  }
  if (Initial > Max)
    Initial=Max;
  if (Initial < Min)
    Initial=Min;
   N_Parm_List++;
   Parm_List(N_Parm_List).Value = &(pop->prior->chrom()[0].locus[0].allele_freq[0]);
   Parm_List(N_Parm_List).Max = Max;
   Parm_List(N_Parm_List).Min = Min;
   Parm_List(N_Parm_List).Name = "QTL_Freq  ";
  }

double Model::multipoint_mixture_approx(int maxiter, double P_var){

  double llh;
  pop->analysis_type = "multipoint_m";
  std::cout << "lets do multipoint with the mixture approximation!" << std::endl;
  Fpmm Farg(0,0.5,P_var); // setup for 0 polygenic loci.
  // Need to setup F as not using FPMM stuff
  pop->multipoint_init_parm(Farg);

  llh = pop->multi_m_log_likelihood_peeling(maxiter);
    std::cout << "QTL pos " << Q_pos << std::endl;
  std::cout << "Natural log llh= " << std::setprecision(10) << llh << " Log10 llh = "
       << std::setprecision(10) << (std::log10(std::exp(1.0))*llh) << std::endl;
  /*
      multipoint_Rec_recalc(14);
      multipoint_QTL_table();
  llh = pop->multi_m_log_likelihood_peeling(maxiter);
    std::cout << "QTL pos " << Q_pos << std::endl;
  std::cout << "Natural log llh= " << std::setprecision(10) << llh << " Log10 llh = " << std::setprecision(10) << (log10(exp(1.0))*llh) << std::endl;
  */
  return llh;
}

void Model::multipoint_display_parms(void) {
  int nl=0;
  if (pop->analysis_type == "multipoint") {
    if (pop->P_ndim > 1) {
      nl= static_cast<int>(0.5*(pop->P_ndim-1));
      std::cout << "P_dim" << pop->P_ndim << " nl = " << nl << std::endl;
      std::cout << "Multipoint with FPMM with ";
      if (nl >1) {
	std::cout << nl << " loci." << std::endl;
      }
      else {
	std::cout << "1 locus." << std::endl;
      }
    }
    else
      std::cout << "Simple Multipoint" << std::endl;
  }
  else if (pop->analysis_type== "multipoint_m") {
    std::cout << "Multipoint with the PAP style mixture approximation" << std::endl;
  }
  else
    std::cout << "what did you run?" << std::endl;
  
  //Now pretty print the final values
  if (N_Parm_List > 0) {
    std::cout << " Parameter    Final Value " << std::endl;
    int i;
    for (i=1; i <= N_Parm_List; i++) {
      if (Parm_List(i).Name != "Distance  ") 
	std::cout <<  Parm_List(i).Name << "       " <<  *Parm_List(i).Value << std::endl;
      else 
	std::cout <<  Parm_List(i).Name << "       " <<  *Parm_List(i).Value-2.0 << std::endl;
    }
  }
}

void Model::multipoint_estimate_PolyV(double Initial, double Max, double Min) {
    // Maximize the polygenic variance under the mixture model approximation
    // not used by any other method!
	    
double tmp;
if (N_Parm_List == 0) {
	     Parm_List.resize(7); // Parm_list is empty therefore resize it
	         }
if (Min < 0.0)
	    Min=0.0;
if (Min > Max) {
	    tmp=Min;
	        Min=Max;
		    Max=tmp;
		        }
if (Initial > Max)
	    Initial=Max;
if (Initial < Min)
	    Initial=Min;
pop->residual_var[0][0]=Initial;
N_Parm_List++;

Parm_List(N_Parm_List).Value = &(pop->F->var);
Parm_List(N_Parm_List).Max = Max;
Parm_List(N_Parm_List).Min = Min;
Parm_List(N_Parm_List).Name = "Poly_Var  ";
}


//BRS

//MJS

// ***************************************************************
// *                                                             *
// *  iterative peeling adapted to accept graph data as input    *
// *                                                             *
// ***************************************************************


double Model::graph_log_likelihood_peeling(const unsigned maxiter)
{
  if (ntermGdist == 0) {
    throw exception("Model::graph_log_likelihood_peeling(): ntermGdist == 0");
    return 0.0;
  }

  if (type == bad_model) {
    throw exception("Model::graph_log_likelihood_peeling(): bad model");
    return 0.0;
  }
  if (!data_prepared) {
    if (!prepare_data()) {
      type = bad_model;
      return 0.0;
    }
  }

  int pt,i;
  Individual *I;
  double sigma_mu,retval = 0.0;
  double *sol, rr;
  unsigned nl,rank_mme,startaddr,iter = 0;
  pop->analysis_type = "not-multi";
  pop->tn_qtl=4;
  if (numterm) {
    retval += numobs*std::log(residual_var[0][0]);
    hmmec.resize(hmmesize,max_nz);
    rellrhs.resize(hmmesize);
    blupsol.resize(hmmesize);

    pop->genotype_dist_peeling(0);
    for (i=0; i<popsize; i++) {
      I = pop->popmember[i];
      I->xbzu_val = I->est_GV;  
    }
    unsigned nvc = totalnvc() - 1;              // number of ratio's
    Vector<double> zz(popsize);
    Vector<double> ratio(nvc);
    Vector<double> vc(nvc+1);
    
    var2vec(vc);
    for (i=0; i<nvc; i++) ratio[i] = vc[nvc]/vc[i];
    for (pt=-1,i=0; i<numterm; i++) {
      if (term[i].classi() == 'P') {
	pt = i;
	startaddr =  term[pt].start + 1;
	nl = term[pt].nlevel();
	break;
      }
    }
    
    setup_ww_single_trait(ratio,pt,startaddr,&zz);
    hmmec.reorder();
    rank_mme = hmmec.factorization(5);
    
    hmmec.solve(blupsol,rellrhs,"ysmp1");
    retval += hmmec.info_vec[0]- rank_mme*std::log(residual_var[0][0]);

    if (pt >= 0) {
      sol = &(blupsol[startaddr-1]);
      for (rr=0.0,i=0; i<nl; i++) {
	rr += (*sol * zz[i] * *sol);  sol++;
      }
      sigma_mu = *(term[pt].prior->var_matrix())[0][0];
      retval += nl*std::log(sigma_mu);
      sigma_mu = hmmec.q(blupsol.begin(),blupsol.begin(),startaddr,startaddr+nl-1)
	- rr;
      retval += std::log(sigma_mu);
    }
    compute_xbzu(blupsol);                   // to adjust trait
  }
  retval *= -0.5;
  pop->residual_var = residual_var;
  pop->input_data(data);
  //  retval += pop->graph_log_likelihood_peeling(maxiter);
  return retval;
}


} //////// end of namespace matvec

