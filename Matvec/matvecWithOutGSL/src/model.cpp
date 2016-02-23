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
#include <sstream>
#include <iomanip>
#include <cstdarg>
#include "histogram.h"
#include "session.h"
#include  "stat.h"
#include "model.h"
#include "population.h"
#include "parmMap.h"

namespace matvec {


void Model::initialize(void)
{
   num_ped        = 0;
   max_nz         = 0;
   numrec         = 0;
   numobs         = 0;
   numcol         = 0;
   numterm        = 0;
   popsize        = 0;
   numfact        = 0;
   maxorder       = 0;
   numtrait       = 0;
   npattern       = 0;
   hmmesize       = 0;
   numgroup       = 0;
   non_zero       = 0;
   type           = bad_model;    // it becomes good mode after having equations - was model_type
   ntermGdist     = 0;
   pop_created    = 0;
   data_prepared  = 0;
   weightinverse  = 0;
   modelpcount    =0;
   modelstr       = 0;

   pop            = 0;
   rve            = 0;
   data           = 0;
   kvec           = 0;
   rawpos         = 0;
   idlist         = 0;
   pos_term       = 0;
   xval_term      = 0;
   traitname      = 0;
   trait_vec      = 0;
   rec_indid      = 0;
   xact_htable    = 0;
   factor_struct  = 0;
   trait_struct   = 0;
   unknown_prior  = 0;

   // **********************************************
   //  information variables for Model::info(stream)
   // **********************************************
   nfunk_in_dfreml = 0;
   reml_value  = 0.0;
   dfreml_called = 0;
}

void Model::copyfrom(const Model& M)
{
   initialize();
   num_ped        = M.num_ped;
   yry            = M.yry;
   data           = M.data;
   term           = M.term;
   max_nz         = M.max_nz;
   numrec         = M.numrec;
   numobs         = M.numobs;
   numcol         = M.numcol;
   max_nz         = M.max_nz;
   popsize        = M.popsize;
   numterm        = M.numterm;
   numfact        = M.numfact;
   lng0vec        = M.lng0vec;
   lnr0vec        = M.lnr0vec;
   blupsol        = M.blupsol;
   rellrhs        = M.rellrhs;
   popsize        = M.popsize;
   hmmesize       = M.hmmesize;
   numtrait       = M.numtrait;
   npattern       = M.npattern;
   hmmesize       = M.hmmesize;
   numgroup       = M.numgroup;
   maxorder       = M.maxorder;
   non_zero       = M.non_zero;
   rec_indid      = M.rec_indid;
   type     = M.type;
   modelfname     = M.modelfname;
   weightname     = M.weightname;
   ntermGdist     = M.ntermGdist;
   modelstring    = M.modelstring;
   pattern_tree   = M.pattern_tree;
   residual_var   = M.residual_var;
   data_prepared  = M.data_prepared;
   weightinverse  = M.weightinverse;
   nfunk_in_dfreml= M.nfunk_in_dfreml;
   reml_value     = M.reml_value;
   dfreml_called  = M.dfreml_called;
   hmmec.resize(hmmesize,max_nz);

   int i,k;
   // pos_term, xval_term, trait_vec are working spaces
   if (pos_term) {delete [] pos_term; pos_term = 0;}
   if (M.pos_term) pos_term = new unsigned [numterm+1];
   if (xval_term) {delete [] xval_term; xval_term = 0;}
   if (M.xval_term) xval_term = new double [numterm+1];
   if (trait_vec) {delete [] trait_vec;  trait_vec= 0;}
   if (M.trait_vec) trait_vec = new double [numterm+1];

   if (unknown_prior) {
      delete [] unknown_prior;unknown_prior = 0;
   }
   if (M.unknown_prior) {
     if(numterm>0){
       unknown_prior = new UnknownDist[numterm];
     }
     else {
       unknown_prior = 0;
     }
      for (i=0; i<numterm; i++) {
         unknown_prior[i].copyfrom(M.unknown_prior[i]);
         if (strcmp(M.term.termlist[i].prior->name(),"UnknownDist") == 0) {
            term.termlist[i].prior = &(unknown_prior[i]);
         }
         else {
            term.termlist[i].prior = M.term.termlist[i].prior;
         }
      }
   }
   if (traitname) {delete [] traitname; traitname = 0;}
   if (M.traitname) {
     if(numtrait>0){
       traitname = new std::string [numtrait];
     }
     else {
       traitname = 0;
     }
      for (i=0; i<numtrait; i++) traitname[i] = M.traitname[i];
   }
   if (kvec) {delete [] kvec; kvec = 0;}
   if (M.kvec) {
     if(numtrait>0){
       kvec = new unsigned [npattern];
     }
     else {
       kvec = 0;
     }
      for (i=0; i<npattern; i++) kvec[i] = M.kvec[i];
   }
   if (rawpos) {delete [] rawpos; rawpos = 0;}
   if (M.rawpos) {
     if(maxorder>0){
       rawpos = new unsigned [maxorder];  // this is just a working space
     }
     else {
       rawpos = 0;
     }
   }
   pattern.resize(npattern);
   std::set<std::string>::iterator pos;
   for (i=0, pos=pattern_tree.begin(); pos != pattern_tree.end(); ++pos,++i) {
      pattern[i] = *pos;
   }
   pop_created = M.pop_created;
   if (M.pop_created) {
      pop = new Population(*(M.pop));
   }
   else {
      pop = M.pop;
   }
   if (rve) {
     delete [] rve;
     rve = 0;
   }
   if (M.rve) {
     if(npattern>0){
       rve = new Matrix<double> [npattern];
     }
     else {
       rve = 0;
     }
      check_ptr(rve);
      for (i=0; i<npattern; i++) rve[i].reserve(numtrait,numtrait);
   }
   if (factor_struct) {
      delete [] factor_struct; factor_struct = 0;
   }
   if (M.factor_struct) {
     if(numfact>0){
       factor_struct = new FieldStruct [numfact];
     }
     else {
       factor_struct = 0;
     }
      for (i=0; i<numfact; i++) factor_struct[i] = M.factor_struct[i];
   }
   if (trait_struct) {
      delete [] trait_struct; trait_struct = 0;
   }
   if (M.trait_struct) {
     if(numtrait>0){
       trait_struct = new FieldStruct* [numtrait];
     }
     else {
       trait_struct = 0;
     }
      for (i=0; i<numtrait; i++) {
         k = fieldindex(traitname[i],factor_struct,numfact);
         if (k >= 0 && k < numfact) {
            trait_struct[i] = &(factor_struct[k]);
         }
         else {
            throw exception("Model::copyfrom(): you have probably found a bug!");
         }
      }
   }
   if (xact_htable) {delete [] xact_htable; xact_htable = 0;}
   if (M.xact_htable) {
     if(numterm>0){
       xact_htable = new HashTable [numterm];
     }
     else {
       xact_htable = 0;
     }
      for (i=0; i<numterm; i++) xact_htable[i] = M.xact_htable[i];
   }
   if (idlist) {delete [] idlist; idlist = 0;}
   if (M.idlist) {
     if(numcol>0){
       idlist = new HashTable* [numcol];
     }
     else {
       idlist = 0;
     }
      check_ptr(idlist);
      if (data) for (i=0; i<numcol; i++) idlist[i] = data->hashtable[i];
   }
}

int Model::equation(const std::string &modelspecs)
{
  //******************************************************************
  //  take the folowing model as example to show how I process model
  //
  //     model M("y1 =  A + B + A * B,"
  //             "y2 =  B         X(B)"
  //
  // *****************************************************************
   int retval = 1;
   num_ped=0;
   if (modelspecs == "") throw exception("Model::equation(): equation is empty");
   modelstring = modelspecs;
   if (build_model_term(modelstring)) {
      maxorder = term.maxorder();
      if (rawpos) {
	delete [] rawpos;
	rawpos=0;
      }
      if(maxorder>0){
	rawpos = new unsigned [maxorder];  // this is just a working space
      }
      else {
	rawpos = 0;
      }
      lng0vec.resize(numterm);
   }
   else {
      modelstring = "";
      retval = 0;
   }
   if (retval == 1) {
      type = fixed_model;
   }
   else {
      type = bad_model;
   }
   return retval;
}

void Model::fitdata(Data& D)
{
   if (D.num_rows() > 0) {
      data = &D;
   }
   else {
      data = 0;
      type = bad_model;
      warning("Model::fitdata(D): Data is empty");
   }
   max_nz     = 0;
   numrec     = 0;
   numobs     = 0;
   numcol     = 0;
   npattern   = 0;
   hmmesize   = 0;
   data_prepared = 0;
   hmmec.resize(0,0);
}

int Model::display(void) const
{
   std::cout << " model:  " << modelstring << "\n";
   return 1;
}

void Model::release(void)
{
   if (idlist)        {
     delete [] idlist; 
     idlist= 0;
   }
   if (xact_htable)   {
     delete [] xact_htable; 
     xact_htable= 0;
   }
   if (kvec)          {
     delete [] kvec; 
     kvec = 0;
   }
   if (rawpos)        {
     delete [] rawpos; 
     rawpos = 0;
   }
   if (pos_term)      {
     delete [] pos_term; 
     pos_term = 0;
   }
   if (xval_term)     {
     delete [] xval_term; 
     xval_term = 0;
   }
   if (trait_vec)     {
     delete [] trait_vec; 
     trait_vec = 0;
   }
   if (traitname)     {
     delete [] traitname; 
     traitname = 0;
   }
   if (pop_created)   {
     delete pop; pop = 0; 
     pop_created=0;
   }
   if (modelstr)      {
     delete [] modelstr; 
     modelstr=0; modelpcount=0;
   }
   if (trait_struct) {
      delete [] trait_struct; 
      trait_struct=0;
   }
   if (factor_struct) {
      delete [] factor_struct; 
      numfact=0; factor_struct = 0;
   }
   if (unknown_prior) {
      delete [] unknown_prior; 
      unknown_prior=0;
   }
   if(rve){
     delete [] rve;
     rve = 0;
   }
   remove(modelfname.c_str());
}

int Model::prepare_data(const std::string &solver)
{
   if (data_prepared) return 1;
   if (modelstring == "" || !data) throw exception("Model::prepare_data(): empty model or no data");

   // **************************************************************
   // it's necessary to make a hardcopy on disk, because data
   // may be changed temporary due to re-hash, etc.
   // **************************************************************
   if (!data->in_disk()) data->save_datasheet();
   if (!data->in_memory()) data->input_datasheet();
//rlf   modelfname = SESSION.mktemp();

   int k,t,ii,ii1;
   unsigned i,nrow_in_data;
   char CL;
   act_numtrait=numtrait;
   //if(num_ped) act_numtrait *= num_ped;
   // **************************************************************
   //  if there is  B*B in a model, then  B must be a covariate
   // **************************************************************
   ModelTerm *T;
   std::string fn;
   unsigned ped_cnt=0,first_ped;
   for (t=0; t<numterm; t++) {
      T = &(term[t]);
      for (i=0; i<T->order()-1; i++) {
         ii = T->factorindx[i]; ii1 = T->factorindx[i+1];
         if (ii == ii1) {
            if (factor_struct[ii].classi() != 'C') {
               fn =  factor_struct[ii].name();
               throw exception("Model::prepare_data(): Is there a covariate?");
//               break;
            }
         }
      }
   }
   nrow_in_data = numrec  = data->num_rows();
   numcol  = data->num_cols();
   // *********************************************************************
   //  (1) If weight exists, then set associated data-column to be 'W'
   //  (2) reset data record_struct[i].classi with factor_struct[j].classi
   // ********************************************************************
   for (i=0; i<numcol; i++) data->datasheet[i].classi('N');
   if (weightname.size() > 0) {
      k = data->field_index(weightname.c_str());
      if (k >= 0) {
         if (data->datasheet[k].type() != 'S') {
            data->datasheet[k].classi('W');
            data->datasheet[k].value_for_missing(0.0);
         }
         else {
            throw exception("Model::prepare_data(): str column can't be a weight column");
         }
      }
      else {
         throw exception("Model::prepare_data(): weight variable not in data");
      }
   }

   for (i=0; i<numfact; i++)  {

      k = data->field_index(factor_struct[i].name());
      if (k >= 0) {
         CL = factor_struct[i].classi();
         if (CL == 'F') {
            if (data->datasheet[k].type() == 'S') {
               CL = 'U';
            }
            else if (data->datasheet[k].type() == 'I') {
               CL = 'I';
            }
         }
         else if (CL == 'C') {
            if (data->datasheet[k].type() != 'F') {
               throw exception("Model::prepare_data(): a covariate is specified as string");
            }
         }
         else {
            if (data->datasheet[k].type() == 'I') {
               throw exception("Model::prepare_data(): 'intercept' can only be a fixed effect!");
            }
         }
         if (data->datasheet[k].classi() == 'W') {
            throw exception("Model::prepare_data(): a weight variable cann't be a part of model euqation");
         }
         data->datasheet[k].classi(CL);
         factor_struct[i].index(data->datasheet[k].index());
         factor_struct[i].nmiss(data->datasheet[k].nmiss());
      }
      else {
	std::cout << "\n" << factor_struct[i].name() << " " << k <<"\n";
         throw exception("Model::prepare_data(): a variable in model, but not in data");
      }
   }

   // ******************************************************************
   // determine the legality for each record: 0 = perfect, 1 = bad
   // ******************************************************************
   DataNode *column;
   Vector<bool> record_legality(nrow_in_data);

   for (t=1; t<numcol; t++) {
      if (data->datasheet[t].nmiss() == 0) continue;
      CL = data->datasheet[t].classi();
      if (!(CL=='F' || CL=='U' || CL=='P' || CL=='C')) continue;
      column = data->rawcol(t);
      for (i=0; i<nrow_in_data; i++) {
         if (record_legality[i]) continue;
         if (column[i].missing) {
            record_legality[i] = true;
            numrec--;
         }
      }
   }

   re_hash_data(record_legality);
   if (maxorder>1) assign_id_xact(record_legality); //assign ids for interation

   // ***********************************************************
   //  determine how many levels for each main term in the model
   //  it has already done for interaction term.
   // ************************************************************
   for (t=0; t<numfact; t++) {
      i = factor_struct[t].index();
      factor_struct[t].nlevel(data->datasheet[i].nlevel());
   }
   for (t=0; t<numterm; t++) {
      if (term[t].order() == 1) {
         term[t].numlevel = factor_struct[term[t].factorindx[0]].nlevel();
      }
   }

   ////////////////// allocate memory for working spaces //////////////////
   if (pos_term) {
     delete [] pos_term;
     pos_term = 0;
   }
   pos_term = new unsigned [numterm+1];
   if (xval_term) {
     delete [] xval_term;
     xval_term = 0;
   }
   xval_term = new double [numterm+1];
   if (trait_vec) {
     delete [] trait_vec;
     trait_vec = 0;
   }
   if(numtrait>0){
     trait_vec = new double [numtrait];
   }
   else {
     trait_vec = 0;
   }

   // ***********************************************************
   // for each term and trait, say, A B X y1 y2, we'll prepare
   // it's position and value in henderson mixed model equations
   //  1  100  209  1.0 1.0 0.23 y1 y2
   // then we dump them on disk
   // Finally, we release data sheet from the memory
   // ***********************************************************
   save_pos_val(record_legality,solver);
   data->release_datasheet();

   record_legality.clear();

   lnr0vec.resize(npattern);
   if (rve){
     delete [] rve;
     rve = 0;
   }
   if(npattern>0){
     rve = new Matrix<double> [npattern];
   }
   else {
     rve = 0;
   }
   check_ptr(rve);
   for (i=0; i<npattern; i++) rve[i].reserve(numtrait,numtrait);

   // *********************************************
   //   hmmec starting position of each model term
   // **********************************************
   hmmesize=0;
   unsigned end=nt_vec[0]*term[0].numlevel;
   /*   if(term[0].classi() == 'P') {
     ped_cnt=1;
     first_ped=0;
     }*/
   for (term[0].start = 0, i=1; i<numterm; i++) {
     /*if(term[i].classi() == 'P') {
       ped_cnt++;
       if(ped_cnt == 1)	 first_ped=i;
       if(first_ped == 0) {
	 term[i].start=ped_cnt-1;
       }
       else{
       term[i].start = term[first_ped-1].start + numtrait*(term[first_ped-1].numlevel+(ped_cnt-1));
       }
     if(hmmesize < (term[i].start+ act_numtrait*term[i].numlevel))  hmmesize= term[i].start+ act_numtrait*term[i].numlevel;
     }make
     else {*/
     //if(term[i-1].classi() == 'P' ) term[i].start = term[first_ped].start + act_numtrait*term[first_ped].numlevel;
     // else 
     if(pos_vec[i]==0){
       term[i].start=end;
       end+=nt_vec[i]*term[i].numlevel;
     }
     else{
       term[i].start = term[base_effect[i]].start + numtrait*pos_vec[i];
     }
       //  if(hmmesize < (term[i].start+ numtrait*term[i].numlevel))  hmmesize= term[i].start+ numtrait*term[i].numlevel;
       /* }*/
   }
   // *****************************************************
   //  get hmmesize and allocate memory for hmmec and rhs
   // *****************************************************
   for (hmmesize=0,t=0; t<numterm; t++) {
     hmmesize += nt_vec[t]*term[t].nlevel();
     if(nt_vec[t]>act_numtrait) act_numtrait=nt_vec[t];
   }
   if (max_nz == 0)  max_nz = est_nze(hmmesize,act_numtrait,popsize);
   data_prepared = 1;
   return 1;
}

void Model::re_hash_data(Vector<bool> &record_legality)
{
   // ********************************************************************
   // * hash all levels of main factors in data for sequantial integer ids
   // ********************************************************************
   if (!data) return;
   if (!(data->in_memory())) data->input_datasheet();

   unsigned i,id,t,nrow_in_data = data->num_rows();
   char CL,TY;
   const char *strid;

   if (idlist) {
     delete [] idlist;
     idlist = 0;
   }
   if(numcol>0){
     idlist = new HashTable* [numcol];
   }
   else {
     idlist = 0;
   }
   check_ptr(idlist);
   HashTable *dhtable,tempidlist;
   int *tempval = new int;
   DataNode *column;

   for (t=1; t<numcol; t++) {    // first column is reserved for intercept
      dhtable = idlist[t] = data->hashtable[t];
      CL = data->datasheet[t].classi();
      TY = data->datasheet[t].type();
      column = data->rawcol(t);
      switch (CL) {
         case 'F':                  // need hash to get its sequential id
            if (data->hashtable[t]->size() == 0) {
               tempidlist.resize(nrow_in_data,sizeof(int));
               for (i=0; i<nrow_in_data; i++) {
                  if (record_legality[i]) continue;
                  *tempval = static_cast<int>(column[i].double_val());
                  tempidlist.insert(tempval);
               }
               id = tempidlist.size();
               data->datasheet[t].nlevel(id);
               dhtable->resize(id,sizeof(int));
               for (i=0; i<nrow_in_data; i++) {
                  if (record_legality[i]) continue;
                  *tempval = static_cast<int>(column[i].double_val());
                  dhtable->insert(tempval);
               }
            }
            break;
         case 'P':
         case 'G':
            // ***********************************************************
            // the integer ids in data will be changed here according to
            // pedigree, but they will not be saved on disk, ie, we'll
            // release datasheet at the end of save_pos_val();
            // data->datasheet[t].nlevel has changed, it
            // must be changed back at the end of save_pos_val();
            // ***********************************************************
	   if(nrow_in_data>0){
	     rec_indid = new unsigned [nrow_in_data];   // Gibbs
	   }
	   else {
	     rec_indid = 0;
	   }
            data->datasheet[t].nlevel(pop->size());
            for (i=0; i<nrow_in_data; i++) {
               if (record_legality[i]) continue;
               if (TY == 'S') {
                  id = column[i].unsigned_val();
                  strid = (const char *)dhtable->find(id);
                  if (pop->maxnamelen > 0) {
                     id = pop->hashtable.get_id(strid);
                  }
                  else {
                     id = 0;
                  }
               }
               else {
                  id = (unsigned) column[i].double_val();
                  if (pop->maxnamelen > 0) {
		    std::ostringstream ostr;
		    ostr << id;
		    std::string t = ostr.str();
		    strid = t.c_str();
		    id = pop->hashtable.get_id(strid);
                  }
               }
               rec_indid[i] = id;                        // Gibbs
               if (id) {
                  column[i].unsigned_val(id);
               }
               else {
                  column[i].missing = 1;
                  record_legality[i] = true;
                  numrec--;
                  warning("Model::re_hash_data(): %s is in the data, "
                        "but not in the pedigree, so it's ignored",strid);
               }
            }
            break;
         case 'C':
         case 'I':
            data->datasheet[t].nlevel(1);
            break;
         // other possible types: U,T,I,N
      }
   }
   if (tempval){
     delete tempval;
     tempval = 0;
   }
}

void Model::save_pos_val(Vector<bool> &record_legality, const std::string &solver)
{
   // ***************************************************************
   // prepare model specific binary data, each line contains
   //   pos_term[t]   xval_term[t]   trait_vec[i]
   // Note:
   // (1) pos_term[t] contains position ranging from 1 to its nlevel
   // (2) the last column in pos_term[] is the missing patterm index
   // (3) the last column in xval_term[] is the weight value
   // ***************************************************************
   if (!data) return;
   if (!(data->in_memory())) data->input_datasheet();
   unsigned nrow_in_data = data->num_rows();
   double rawxval;
   Vector<unsigned> pos_datacol(numcol);  // working space
   Vector<double> xval_datacol(numcol);
   int *tempval = new int;
   unsigned i,j,k,t,nord;
   unsigned ntm1 = numterm + 1;
   char CL;

   // *****************************************************************
   // (1) find missing pattern,
   // (2) kick out record with all traits being missing,
   // (3) compute the number of full record (ie, no trait is missing
   // ****************************************************************
   char *patstring = new char [numtrait+1];
   patstring[numtrait]='\0';
   pattern_tree.clear();
   DataNode *dn;
   Vector<double> y_mean(numtrait);
   Vector<double> y_std(numtrait);
   Vector<unsigned> y_nlevel(numtrait);
   for (numobs=0,i=0; i<nrow_in_data; i++) {
      if (record_legality[i]) continue;
      for (k=0, t=0; t<numtrait; t++) {
         j = trait_struct[t]->index();
         dn = data->cell(i,j);
         if (dn->missing) {
            patstring[t] = '0';
            k++;
         }
         else {
            patstring[t] = '1';
            y_mean[t] += dn->double_val();
            y_std[t] += dn->double_val()*dn->double_val();
            y_nlevel[t] += 1;
         }
      }
      if (k == numtrait) {
         record_legality[i] = true;
         numrec--;
         for (k=1; k<numcol; k++) {
            CL = data->datasheet[k].classi();
            if (CL != 'C') continue;
            data->cell(i,k)->missing = 1;
            data->datasheet[k].count_miss(1);
         }
      }
      else {
         numobs++;
         pattern_tree.insert(patstring);
      }
   }
   double nr;
   for (t=0; t<numtrait; t++) {
      nr = static_cast<double>(y_nlevel[t]);
      if (nr == 0.0) continue;
      y_mean[t] /= nr;
      if (nr > 1.0) {
         y_std[t] = std::sqrt((y_std[t] - y_mean[t]*y_mean[t]*nr)/(nr-1.0));
      }
      trait_struct[t]->mean(y_mean[t]);
      trait_struct[t]->std(y_std[t]);
      trait_struct[t]->nlevel(y_nlevel[t]);
   }

   // ******************************************************************
   // compute means for covariates
   // ******************************************************************
   for (t=1; t<numcol; t++) {
      if (data->datasheet[t].classi() == 'C') data->datasheet[t].mean();
   }

   npattern = pattern_tree.size();
   pattern.resize(npattern);
   std::set<std::string>::iterator pos;
   for (i=0,pos = pattern_tree.begin(); pos != pattern_tree.end(); ++pos,++i) {
      pattern[i] = *pos;
   }
   // note pattern is in ascending order
   if (kvec) {
     delete [] kvec;
     kvec = 0;
   }
   if(npattern>0){
     kvec = new unsigned [npattern];
   }
   else {
     kvec = 0;
   }
   memset(kvec,'\0',sizeof(unsigned)*npattern);
   DataNode* recd;
   if(numcol>0){
     recd = new DataNode [numcol];
   }
   else {
     recd = 0;
   }

   if (modelstr) {delete [] modelstr; modelstr=0; modelpcount=0;}
   std::stringstream bmodfile;
   if (!bmodfile) throw exception(std::string("Model::save_pos_val(): cannot open file:") + modelfname);
   if(solver=="iod"){
     pos_val_vector.resize(nrow_in_data);
   }
   for (i=0; i<nrow_in_data; i++) {
      if (record_legality[i]) continue;
      xval_term[numterm] = 1.0;             // default weight value
      data->row(i,recd);
      for (j=0; j<numcol; j++) {
         CL = data->datasheet[j].classi();
         if (CL =='F') {                     // non string column
            xval_datacol[j] = 1.0;
            *tempval = static_cast<int>(recd[j].double_val());
            pos_datacol[j] = idlist[j]->get_id(tempval);
         }
         else if (CL=='U' || CL=='P') {    // string column
            xval_datacol[j] = 1.0;
            pos_datacol[j] = recd[j].unsigned_val();
         }
         else if (CL=='C') {
            xval_datacol[j] = recd[j].double_val();
            pos_datacol[j] = 1;
         }
         else if (CL=='W') {
            rawxval = recd[j].double_val();
            if (rawxval >= 0) {           // negative weight is not allowed
               if (weightinverse == 1) {  //  1/weight is demonded by user
                 if (rawxval > 0.0) rawxval = 1.0/rawxval;
               }
               xval_term[numterm] = rawxval;
            }
         }
         else if (CL=='I') {                     // Intercept
            xval_datacol[j] = 1.0;
            pos_datacol[j] = 1;
         }
      }

      for (t=0; t<numterm; t++) {
         rawxval = 1.0;
         nord = term[t].order();
         for (j=0; j<nord; j++) {
            k = factor_struct[term[t].factorindx[j]].index();
            rawpos[j] = pos_datacol[k];
            rawxval *= xval_datacol[k];
         }
         if (nord == 1) {
            pos_term[t] = rawpos[0];
         }
         else {
            pos_term[t] = xact_htable[t].get_id(rawpos);
         }
         xval_term[t] = rawxval;
      }
      for (t=0; t<numtrait; t++) {
         k = trait_struct[t]->index();
         if (recd[k].missing) {
            patstring[t] = '0';
            trait_vec[t] = 0.0;
         }
         else {
            patstring[t] = '1';
            trait_vec[t] = recd[k].double_val();
         }
      }
      for (k=0; k<npattern; ++k) {
         if (pattern[k] == patstring) {break;}
      }
      if (k == npattern) throw exception("Model::save_pos_val(): you have probably found a bug!");
      kvec[k] += 1;
      pos_term[numterm] = k;
      if(solver=="iod"){
	pos_val_vector[i].pos_term.resize(ntm1);
	pos_val_vector[i].xval_term.resize(ntm1);
        pos_val_vector[i].trait_vec.resize(numtrait);
	for (unsigned j=0;j<ntm1;j++){
	  pos_val_vector[i].pos_term[j]  = pos_term[j];
          pos_val_vector[i].xval_term[j] = xval_term[j];
	}
	for (unsigned j=0;j<numtrait;j++){
          pos_val_vector[i].trait_vec[j] = trait_vec[j];
	}
      }
      else{
	bmodfile.write((char *)pos_term, ntm1*sizeof(unsigned));
	bmodfile.write((char *)xval_term, ntm1*sizeof(double));
	bmodfile.write((char *)trait_vec, numtrait*sizeof(double));
      }
   }
//   bmodfile.close();
//   bmodfile.freeze(1); What is this?
   // for gcc-3.2 replaced 
   //    modelstr=bmodfile.str();
   //    modelpcount=bmodfile.pcount();
   // with 
   modelstringstr=bmodfile.str();
   if (recd) {
     delete [] recd;
     recd = 0;
   }
   if (patstring) {
     delete [] patstring;
     patstring = 0;
   }
   if(tempval){
     delete tempval;
     tempval = 0;
   }

   // ***************************************************************
   //  now we change back data->datasheet[t].nlevel, if necessary
   //  then we release datasheet so that the change we made does not
   //  affect the original data
   // ***************************************************************
   for (t=0; t<numcol; t++) {
      if (data->datasheet[t].classi() =='P') {
         data->datasheet[t].nlevel(data->hashtable[t]->size());
      }
   }
   data->release_datasheet();    // datasheet has been released
}

void Model::assign_id_xact(const Vector<bool> &record_legality)
{
   // **********************************************************
   // get sequantial integer ids for interation term in model
   // ***********************************************************
   if (!data) return;
   unsigned nord, nrow_in_data = data->num_rows();

   if (xact_htable) {
     delete [] xact_htable;
     xact_htable = 0;
   }
   if(numterm>0){
     xact_htable = new HashTable [numterm];
   }
   else {
     xact_htable = 0;
   }
   int t;
   for (t=0; t<numterm; t++) {
      if (term[t].order() > 1) {
         xact_htable[t].resize(nrow_in_data,term[t].order()*sizeof(unsigned));
      }
   }
   hashxact(record_legality);
   for (t=0;t<numterm;t++) {
      nord = term[t].order();
      if (nord>1) {
         term[t].numlevel = xact_htable[t].size();
         xact_htable[t].resize(term[t].nlevel(),nord*sizeof(unsigned));
      }
   }
   hashxact(record_legality);
}

void Model::hashxact(const Vector<bool> &record_legality)
{
   if (!data) return;
   char CL;
   int i,j,t,k,nord;
   int *tempval = new int;
   Vector<unsigned> pos_datacol(numcol);
   DataNode* recd; 
   if(numcol>0){
     recd = new DataNode [numcol];
   }
   else {
     recd = 0;
   }
   unsigned nrow_in_data = data->num_rows();

   for (i=0; i<nrow_in_data; i++) {
      if (record_legality[i]) continue;
      data->row(i,recd);
      pos_datacol.initialize(numcol,0);
      for (t=0; t<numterm; t++) {
         nord = term[t].order();
         if (nord > 1) {
            for (j=0; j<nord; j++) {
               k = factor_struct[term[t].factorindx[j]].index();
               if (pos_datacol[k] == 0) {
                  CL = data->datasheet[k].classi();
                  if (CL == 'F') {
                     *tempval = static_cast<int>(recd[k].double_val());
                     pos_datacol[k] = idlist[k]->get_id(tempval);
                  }
                  else if (CL == 'U' ||CL == 'P') {
                     pos_datacol[k] = recd[k].unsigned_val();
                  }
                  else if (CL =='C') {
                     pos_datacol[k] = 1;
                  }
               }
               rawpos[j] = pos_datacol[k];
            }
            xact_htable[t].insert(rawpos);
         }
      }
   }
   if (recd) {
     delete [] recd;
     recd = 0;
   }
   if(tempval){
     delete tempval;
     tempval = 0;
   }
}

void Model::prior_dist(const std::string &termname, GeneticDist *D)
  // ***************************************************************
  //  prior is the method to specify the prior distributions for
  //  each terms in the linear model.
  // ***************************************************************
{
   int k = term.index(termname,factor_struct);
   if (k < 0) throw exception("Model::prior_dist(): term not in the model equation");
   term[k].prior = D;
   term[k].classi('R');
   if (strcmp(D->name(),"GeneticDist") == 0) {
      numterm--;
      ntermGdist++;
      if (k<numterm) term.swap(k,numterm); // move non-Gdist T before Gdist T
   }
}

void Model::prior_dist(const std::string &termname,Pedigree& P, GeneticDist *D)
  // ***************************************************************
  //  prior is the method to specify the prior distributions for
  //  each terms in the linear model.
  // ***************************************************************
{
   int k = term.index(termname,factor_struct);
   if (k<0) throw exception("Model::prior_dist(): term not in the model equation");
   if (D->ntrait() != 1) throw exception("Model::prior_dist(): only one trait allowed");
   term[k].prior = D;
   if (term[k].order() > 1) throw exception("Model::prior_dist(): pedigree cannot be for this term");

   term[k].classi('G');                   // factor classi will override
   factor_struct[term[k].factorindx[0]].classi('G'); // data column classi later
   if (pop_created) delete pop;
   pop = new Population;
   check_ptr(pop);
   pop_created = 1;
   try {
     pop->input_ped(P,D);
   }
   catch (exception &ex) {
      type = bad_model;
      throw;
   }
   numgroup = P.ngroup();
   popsize = pop->size();
   if (strcmp(D->name(),"GeneticDist") == 0) {
      term[k].numlevel = popsize;
      numterm--;
      ntermGdist++;
      if (k<numterm) term.swap(k,numterm); // move non-Gdist T before Gdist T
   }
}

void Model::RSamplerPrior_dist(const std::string &termname,Pedigree& P, Data *D,GeneticDist *G){
  // Authors: L. Radu Totir and Rohan L. Fernando 
  // (June, 2003) 
  // Contributors: 
  int k = term.index(termname,factor_struct);
  if (k<0) throw exception("Model::prior_dist(): term not in the model equation");
  if (G->ntrait() != 1) throw exception("Model::prior_dist(): only one trait allowed");
  term[k].prior = G;
  if (term[k].order() > 1) throw exception("Model::prior_dist(): pedigree cannot be for this term");
  
  term[k].classi('G');                   // factor classi will override
  factor_struct[term[k].factorindx[0]].classi('G'); // data column classi later
  if (pop_created) delete pop;
  pop = new Population;
  check_ptr(pop);
  pop_created = 1;
  pop->model = this;
  try {
    pop->setupRSampler(P,D,G);
    // newHG
    pop->ListAlleleFounders();                         
    pop->SetPossibleHaplotypes(); 
    // newHG
  }
  catch (exception &ex) {
      type = bad_model;
	  cerr << ex.what() << endl;  
      throw;
  }
  numgroup = P.ngroup();
  popsize = pop->size();
  if (strcmp(G->name(),"GeneticDist") == 0) {
    term[k].numlevel = popsize;
    numterm--;
    ntermGdist++;
    if (k<numterm) term.swap(k,numterm); // move non-Gdist T before Gdist T
  }
}
/*! \fn void Model::RSamplerPrior_dist(const std::string
    &termname,Pedigree& P, Data *D,GeneticDist *G)
 *  \brief method to specify the prior distributions for
 each terms in the model for which the R-sampler is invoked; it also
 sets up the structures needed by the peeling process.
*/

void Model::DGSamplerSetup(unsigned numLoci, Data *D){
	// Authors: L. Radu Totir, Chris Stricker
	// (August, 2004) 
	// Contributors:
	pop->descent_graph_setup(D);
	pop->allele_vector1.resize(numLoci); //need to put this somewhere before building
	pop->allele_vector2.resize(numLoci); //allele vectors. Need to do this only once.
	pop->previousAlleleVector1.resize(numLoci);
	pop->previousAlleleVector2.resize(numLoci);
	pop->connected_groups.resize(numLoci);
	pop->previousConnectedGroups.resize(numLoci);
	pop->connect_counter.resize(numLoci);
	pop->previousConnectCounter.resize(numLoci);
	pop->founder_allele_neighbors.resize(numLoci);
	pop->founder_allele_neighbors_all.resize(numLoci);
	pop->previousFounderAlleleNeighbors.resize(numLoci);
	
	pop->founder_allele_counter = 0;
	int popsize = pop->size();
	int zero = 0;  
	for (int i=0;i<popsize;i++) {  
		pop->member(i)->set_founder_alleles(0);
	}
	for(unsigned lcs=1; lcs<=numLoci; lcs++){
		pop->previousAlleleVector1(lcs).resize(pop->founder_allele_counter,0);  // Here we can allocate space for the second dimension of the vector 
		pop->previousAlleleVector2(lcs).resize(pop->founder_allele_counter,0);  // to store the previously sampled allele vector, we want to do this 
																	  // once and not redo it for every DG, as the vector shall not be 
																	  // resized and thus initialized for every new DG as it would loose 
																	  // the values stored in it that we need to restore the previous values.
		pop->previousConnectedGroups(lcs).resize(pop->founder_allele_counter,0);
		pop->previousFounderAlleleNeighbors[lcs-1].resize(pop->founder_allele_counter);
		
  	}
	
	unsigned lastLocus = numLoci;
	for (int i=0; i<popsize; i++) {
		if (pop->member(i)->m_counter_map.empty()){
			int zero = 0;    
			pop->member(i)->m_counter_map.resize(numLoci,zero);
			pop->member(i)->p_counter_map.resize(numLoci,zero);      
		}
//		pop->member(i)->search_heteroz(lastLocus);  
		//cout << pop->member(i)-> id() << " " << pop->member(i)-> ord_heter << endl;
	}       
}

void Model::DGSampler(unsigned PQTL, bool &initSamplerToOrderFounders, unsigned numOfSamples,unsigned numOfSL,unsigned numOfHaplo,unsigned numOfSLCas, char* fname){
	// Authors: L. Radu Totir 
	// (August, 2004) 
	// Contributors: Chris Stricker
	double l_hood = 0.0;
	unsigned k = 0, j = 0, option = 0;
	unsigned k1 = numOfSL;
	unsigned k2 = numOfHaplo;
	unsigned k3 = numOfSLCas;
	cout<< endl << "Descent graph iterations = " << numOfSamples << endl;
	for(int sample=1; sample <= numOfSamples; sample++) {
		if(sample%5000 == 0) {
			cout<<sample<<".....";
			cout.flush();
		}
		//	cout<<endl<< "SAMPLE.................................. " << sample <<endl;  
		if (k >= k1 + k2 + k3)                {k = 0;}
		if (k < k1)                           {option = 1;}  
		else if ((k >= k1) && (k < k1 + k2))  {option = 2;}
		else if (k <= k1 + k2 + k3)           {option = 3;}
		k++;    
		j++;    
		
		l_hood = pop->MH_ibd_sample_map(PQTL,l_hood,option,initSamplerToOrderFounders);    
		pop->sum_descentState_map();  // update pdq's for all marker loci 
									  //    if (j >= 100000) {     
									  //      j = 0;    
									  //      cout << "SI_PDQ sample = " << sample << endl;  
									  //    }  
		pop->accumulateFounderHaplotypeOrigin(PQTL);
		pop->accumulateHaplotypes();  // this will reconstruct the haplotypes for each individual at all markers, use pop->displayHaplotypes(nsamples) to see them
	}
	cout << endl << "Haplotype origin probabilities for adjacent marker intervals, i.e. locus 1-2, 2-3, ..., nLoci-1-nloci for all individuals:"<<endl;
	cout << "Two locus origin m-m, m-p, p-m, p-p, each for maternal and paternal haplotype (the first 4 probabilities are for the maternal, the next 4 are for the paternal haplotype of an individual. "<<endl;
	cout << "Repeated for each interval (within line) and individual (across lines). First number on each line is individual ID." <<endl;

	std::ofstream outfile(fname);
	if (!outfile) throw exception("DGSampler(): cannot open file or it already exists");

	for (int i=0; i<pop->size(); i++) {
		pop->member(i)->display_freq_haplotype(numOfSamples);
		pop->member(i)->pdq_grid(numOfSamples,outfile);
	}
	cout << endl; 
}

void Model::RSamplerInitialDG(const std::string &samplerType,const std::string &whatToCompute){
  // Authors: L. Radu Totir 
  // (August, 2004) 
  // Contributors: 
  unsigned maxsize = myRSamplerParms.maxCutsetSize;  
  pop->getInitialGNodeListSample(maxsize,0,pop->popsize,pop->popsize,samplerType);
  pop->copyGNodeStatesToCandidateStates(samplerType);
  pop->copyCandidateToAccepted(samplerType);
  pop->sampleSegregationIndicators();
  pop->copyCandidateGamete();
}
/*! \fn void Model::RSamplerInitialDG(const std::string &samplerType,const 
std::string &whatToCompute)
*   \brief method to generate the initial descent graph to be used by pop_graph
*/

void Model::RSampler(std::string inputFileName){
	// Authors: Rohan L. Fernando
	// (November, 2004) 
	// Contributors:
	if (type == bad_model){
		return;
	}
	matvec::ParmMap parmMap;
	// These are the default values 
	parmMap["pedBlockSize"]  = "0";
	parmMap["numLoci"]       = "1";
	parmMap["numOfSamples"]  = "1000";
	parmMap["numOfBurnIn"]   = "0";
	parmMap["maxCutsetSize"] = "16384";
	parmMap["resultsFile"]   = "";
	parmMap["samplerType"]   = "genotypic";
	parmMap["samplerUsed"]   = "MH";
	parmMap["whatToCompute"] = "genotypeFreq";
	parmMap["howToSample"]   = "single";
	parmMap["printFlag"]     = "0";
	parmMap["startLocus"]    = "0";
	parmMap["startLocusType"] = "fixed";
	// Now we read in the specific values for this job
	parmMap.inputParms(inputFileName);
	
	matvec::RSamplerParms rSamplerParms;
	rSamplerParms.pedBlockSize  = parmMap.getUnsignedValue("pedBlockSize");
	rSamplerParms.numLoci       = parmMap.getUnsignedValue("numLoci");
	rSamplerParms.numOfSamples  = parmMap.getUnsignedValue("numOfSamples");
	rSamplerParms.numOfBurnIn   = parmMap.getUnsignedValue("numOfBurnIn");
	rSamplerParms.maxCutsetSize = parmMap.getUnsignedValue("maxCutsetSize");
	rSamplerParms.printFlag     = parmMap.getUnsignedValue("printFlag");
	rSamplerParms.startLocus    = parmMap.getUnsignedValue("startLocus");
	rSamplerParms.resultsFile   = parmMap["resultsFile"];
	rSamplerParms.samplerType   = parmMap["samplerType"];
	rSamplerParms.samplerUsed   = parmMap["samplerUsed"];
	rSamplerParms.whatToCompute = parmMap["whatToCompute"];
	rSamplerParms.howToSample   = parmMap["howToSample"];
	rSamplerParms.startLocusType = parmMap["startLocusType"];
	rSamplerParms.display();
	RSampler(rSamplerParms);	
}
	

void Model::RSampler(RSamplerParms rsParms){
	// Authors: Rohan L. Fernando
	// (October, 2004) 
	// Contributors:
	
	myRSamplerParms = rsParms;
	if (myRSamplerParms.samplerType == "genotypic"){
		myRSamplerParms.pedBlockSize = rsParms.pedBlockSize ? rsParms.pedBlockSize : pop->size();
	}
	else if(myRSamplerParms.samplerType == "allelic"){
		myRSamplerParms.pedBlockSize = rsParms.pedBlockSize ? rsParms.pedBlockSize : 2*pop->size();		
	}
	RSampler(myRSamplerParms.pedBlockSize, myRSamplerParms.numLoci,     myRSamplerParms.numOfSamples, myRSamplerParms.samplerType,
			 myRSamplerParms.howToSample,  myRSamplerParms.samplerUsed, myRSamplerParms.whatToCompute);
	
}

void Model::RSampler(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType,
					 const std::string &howToSample, const std::string &samplerUsed, const std::string &whatToCompute){
  // Authors: L. Radu Totir 
  // (September, 2004) 
  // Contributors:
  if(howToSample=="single"){
    if(samplerUsed=="gibbs"){
      RSamplerGibbs(pedBlockSize,numLoci,numOfSamples,samplerType,whatToCompute);
    }
    else if(samplerUsed=="MH"){
      RSamplerMH(pedBlockSize,numLoci,numOfSamples,samplerType,whatToCompute);
    }
    else if(samplerUsed=="gibbsMH"){
      unsigned stepMH = 3;
      RSamplerGibbsMH(pedBlockSize,numLoci,numOfSamples,stepMH,samplerType,whatToCompute);
    }
    else {
      cerr << "samplerUsed must be one of: gibbs, MH or gibbsMH;" << endl;
    }  
  }
  else if(howToSample=="joint"){
    if (samplerType!="allelic"){
      cerr << "For joint sampling, samplerType must be set to allelic;" << endl;
      exit(1);
    } 
    if(samplerUsed=="gibbs"){
      cerr << "For joint sampling, please set the samplerUsed to MH for now;" << endl;
      exit(1);
    }
    else if(samplerUsed=="MH"){
      RSamplerMH(numLoci,numOfSamples,samplerType,whatToCompute);
    }
    else if(samplerUsed=="gibbsMH"){
      cerr << "Not implemented yet for joint sampling, please set samplerUsed to MH for now;" << endl;
      exit(1);
    }
    else {
      cerr << "samplerUsed must be one of: gibbs, MH or gibbsMH;" << endl;
    } 
  }
}      

void Model::RSamplerGibbs(unsigned pedBlockSize,unsigned numLoci,unsigned numOfSamples,const std::string &samplerType,const std::string &whatToCompute){
	// Authors: L. Radu Totir 
	// (June, 2004) 
	// Contributors: 
	// if (pedBlockSize!=pop->popsize) throw exception("Model::RSamplerGibbs(...): size of pedigree block must be equal to the population size");
	double maxNumGNodeElements = 4.0;
	int    maxDimensionCutSet  = 7;
	unsigned maxsize = (unsigned)(std::pow((maxNumGNodeElements),(maxDimensionCutSet))+.5);
	matvec::Matrix<double> IBDCovMatrix;
	IBDCovMatrix.resize(1,1,0.0);
	// cout << maxsize << endl; 
	// next get ready and do the Gibbs sampling in blocks
	unsigned stopBlock=0,startBlock=0, count=1;
	for (int i=0;i<numOfSamples;i++){
		if (i==0){
			// get initial sample by sequential inputation
			// cout << "Sample " << 1 << endl;
			pop->getInitialGNodeListSample(maxsize,0,pop->popsize,pedBlockSize,samplerType);
			if (whatToCompute=="haplotypeFreq"){
				pop->countHaplotypes(samplerType);
			}
			else if(whatToCompute=="probDescHaplo"){
				pop->SetFreqHaploFounders();
			}
			else if(whatToCompute=="ibdCovMatrix"){
				IBDCovMatrix.resize(pop->popsize,pop->popsize,0.0);
			}
		}
		else{
			// ask it to use both left and right loci 
			pop->lookToYourLeft = true;
			pop->lookToYourRight= true;
			// cout << "Sample " << setw(5) << i+1 << endl; 
			do{ 
				if(stopBlock < pop->popsize){
					stopBlock += pedBlockSize;
					startBlock = stopBlock - pedBlockSize;
				}
				else {
					stopBlock  = pedBlockSize;
					startBlock = stopBlock - pedBlockSize;
				}
				//cout << "startBlock = " << startBlock << endl;
				//cout << "stopBlock  = " << stopBlock  << endl;
				pop->getGNodeListSample(maxsize,startBlock,stopBlock,pedBlockSize,samplerType);
				if(whatToCompute=="probDescHaplo"){
					pop->UpdateFreqHaploFounders();
				}
			} while(stopBlock < pop->popsize);
		}
		if (whatToCompute=="haplotypeFreq"){
			pop->countHaplotypes(samplerType);
		}
		else if(whatToCompute=="probDescHaplo"){
			pop->sampleSegregationIndicators(); 
		} 
		else if(whatToCompute=="ibdCovMatrix"){
			IBDCovMatrix += pop->getIBDMatrix();
		}
	}
	if (whatToCompute=="haplotypeFreq"){
		pop->displayHaplotypeFrequencies(numOfSamples);
	}
	else if(whatToCompute=="probDescHaplo"){
		pop->CalcFreqHaploFounders(numOfSamples);
		pop->DisplayFreqHaploFounders();
	}
	else if(whatToCompute=="ibdCovMatrix"){
		char *s = "./tryRES/IBDCovMatrix"; 
		ofstream outfile;
		outfile.open(s);
		outfile << IBDCovMatrix/numOfSamples;
	}
}
/*! \fn void Model::RSamplerGibbs(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType="genotypic",const std::string &whatToCompute="genotypeProb")
*  \brief method to sample ordered genotypes using the Gibbs sampler to sample across loci ( the first sample is obtained by sequential imputation)
*/

void Model::RSamplerMH(unsigned pedBlockSize,unsigned numLoci,unsigned numOfSamples,const std::string &samplerType,const std::string &whatToCompute){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// Contributors: 
	histogram h1;
	unsigned maxsize = myRSamplerParms.maxCutsetSize;
	std::string startLocusType = myRSamplerParms.startLocusType;
	int    badCount  = 0;
	int    goodCount = 1; // 1 because first sample is always accepted
	double alpha      = 1.0;
	double ranNumber  = 0.0;
	pop->gNodeList.logProposal = 0.0;
	pop->gNodeList.logTarget   = 0.0;
	matvec::Matrix<double> IBDCovMatrix, meansVector;
	IBDCovMatrix.resize(1,1,0.0);
	pop->identifyLociToSamplePosFor();
	// next get ready and do the MH sampling in blocks
	unsigned stopBlock=0,startBlock=0;
	if(samplerType=="genotypic"){
		stopBlock = pop->popsize;
	}
	else if(samplerType=="allelic"){
		stopBlock = 2*(pop->popsize);
	} 
	unsigned k=1,k1=1,k2=1;
	for (int i=0;i<numOfSamples;i++){
		if(startLocusType=="random") myRSamplerParms.startLocus = matvec::SESSION.mtr.randInt(numLoci-1);
		if(i==0){
			// get initial sample
		cout << "Sample " << 1 << endl;
			pop->getInitialGNodeListSample(maxsize,startBlock,stopBlock,pedBlockSize,samplerType);
			stopBlock=0;
			logOldTarget = pop->gNodeList.logTarget;
			pop->gNodeList.logOldProposal = pop->gNodeList.logProposal;
			
			pop->copyGNodeStatesToCandidateStates(samplerType);
			pop->storeSampledGametes();
			pop->copyCandidateToAccepted(samplerType);
			if (pop->samplePosForMe!=ULONG_MAX){
				acceptedPosition = pop->prior->chrom()[0].locus[pop->samplePosForMe].distance;
			}
			if(whatToCompute=="probDescHaplo"){
				pop->SetFreqHaploFounders();
			}
			else if(whatToCompute=="ibdCovMatrix"){
				IBDCovMatrix.resize(pop->popsize,pop->popsize,0.0);
				meansVector.resize(pop->popsize ,1           ,0.0);
			}
			else if(whatToCompute=="genotypeFreq"){
				pop->initGenotypeFreq();
			}
			else if(whatToCompute=="pdq"){
				pop->initPDQs();
			}
		}
		else{
			pop->samplePositionOnChromosome();
			if(k<=k1){
				TransitionSet::transmissionType = 1;
				k++;
			}
			else if(k<=k2){
				TransitionSet::transmissionType = 2;
				k++;
			}
			if(k==k2){
				k = 1;
			}
			cout << "Sample " << setw(5) << i+1 << endl; 
			do{ 
				if(stopBlock < pop->popsize){
					stopBlock += pedBlockSize;
					startBlock = stopBlock - pedBlockSize;
				}
				else {
					stopBlock  = pedBlockSize;
					startBlock = stopBlock - pedBlockSize;
				}
				// cout << "startBlock = " << startBlock << endl;
				// cout << "stopBlock  = " << stopBlock  << endl;
				// q(NEW) and pi(NEW)
				pop->gNodeList.logProposal = 0.0;
				pop->gNodeList.logTarget   = 0.0;
				try{
				pop->getGNodeListSample(maxsize,startBlock,stopBlock,pedBlockSize,samplerType);
				pop->copyGNodeStatesToCandidateStates(samplerType);
				pop->storeSampledGametes();
				logNewProposal  = pop->gNodeList.logProposal + pop->logNewProposalForPosition;
				logNewTarget    = pop->gNodeList.logTarget;
				logOldProposal  = pop->gNodeList.logOldProposal + pop->logOldProposalForPosition;
				if(TransitionSet::transmissionType!=1){
					//cout << "calculate q(OLD) for sample " << i+1 << endl;	
					pop->gNodeList.logOldProposal = 0.0;
					pop->getOldGNodeListProbability(maxsize,startBlock,stopBlock,pedBlockSize,samplerType);
					logOldProposal = pop->gNodeList.logOldProposal + pop->logOldProposalForPosition;
				}
				//cout << "logNewTarget   = " << logNewTarget   << endl;
				//cout << "logOldTarget   = " << logOldTarget   << endl;
				//cout << "logNewProposal = " << logNewProposal << endl;
				//cout << "logOldProposal = " << logOldProposal << endl;
				// accept or reject the new sample
				alpha=std::exp(logNewTarget + logOldProposal - logNewProposal - logOldTarget);
				}
				catch(matvec::InvalidSample){
					pop->gNodeList.clearGNodeListForNextLocus();
					alpha = 0.0;
				}
				//cout << " alpha = " << alpha << endl;
				//cout << (*pop->member(4));
				//cout << endl;
				ranNumber=ranf();
				if (ranNumber <= alpha) {
					logOldTarget = logNewTarget;
					if(TransitionSet::transmissionType==1){
						pop->gNodeList.logOldProposal = pop->gNodeList.logProposal;
					}
					pop->copyCandidateToAccepted(samplerType);
					if (pop->samplePosForMe!=ULONG_MAX){ 
						acceptedPosition = pop->prior->chrom()[0].locus[pop->samplePosForMe].distance;
						cout << "Sample no: " << i << endl;
						cout << "Accepted position = " << acceptedPosition << endl;
					}
					goodCount++;
				}
				else {
					//pop->storeSampledGametes();
					badCount++;
					//cout << "badCount = " << badCount << endl;
					if (pop->samplePosForMe!=ULONG_MAX){
						pop->prior->chrom()[0].locus[pop->samplePosForMe].distance=acceptedPosition;
					}
				}
				if (pop->samplePosForMe!=ULONG_MAX){
				   h1.push_back(acceptedPosition);
				}
				//cout << "goodCount = " << goodCount << endl;
				if(whatToCompute=="probDescHaplo"){
					pop->UpdateFreqHaploFounders();
				}
				if (!(i%10000)){
					cout << "No of samples accepted = " << goodCount << " out of a total of " << goodCount+badCount << endl;
					double ratio = double(badCount)/(goodCount+badCount);
					// cout << "Rejection rate = " << ratio << endl;
				}
			} while(stopBlock < pop->popsize);
		}
		if (whatToCompute=="haplotypeFreq"){
			pop->countHaplotypes(samplerType);
		}
		else if (whatToCompute=="genotypeFreq"){
			pop->countGenotypes(samplerType);
		}		
		else if(whatToCompute=="probDescHaplo"){
			pop->sampleSegregationIndicators();
		} 
		else if(whatToCompute=="ibdCovMatrix"){
			IBDCovMatrix += pop->getIBDMatrix();
			meansVector  += pop->getMeans();
		}
		else if(whatToCompute=="initialDG"){
			pop->sampleSegregationIndicators();
		}
		else if(whatToCompute=="pdq"){
			pop->sampleSegregationIndicators();
			pop->sum_descentState(2);
		}
	}
	cout << "No of samples accepted = " << goodCount << endl;
	double ratio = double(badCount)/(goodCount+badCount);
	cout << "Rejection rate = " << ratio << endl;
	if (whatToCompute=="haplotypeFreq"){
		if(myRSamplerParms.resultsFile == ""){
			pop->displayHaplotypeFrequencies(numOfSamples);
		}
		else {
			std::ofstream outfile(myRSamplerParms.resultsFile.c_str());
			pop->displayHaplotypeFrequencies(numOfSamples, outfile);
		}
	}
	else if (whatToCompute=="genotypeFreq"){
		if(myRSamplerParms.resultsFile == ""){
			pop->displayGenotypeFrequencies(numOfSamples);
		}
		else {
			std::ofstream outfile(myRSamplerParms.resultsFile.c_str());
			pop->displayGenotypeFrequencies(numOfSamples, outfile);
		}
	}
	else if (whatToCompute=="initialDG"){
		if(myRSamplerParms.resultsFile == ""){
			pop->displaySegregationIndicators();
		}
		else {
			std::ofstream outfile(myRSamplerParms.resultsFile.c_str());
			pop->displaySegregationIndicators(outfile);
		}
	}
	else if(whatToCompute=="probDescHaplo"){
		pop->CalcFreqHaploFounders(numOfSamples);
		pop->DisplayFreqHaploFounders();
	}
	else if(whatToCompute=="ibdCovMatrix"){
		cout << IBDCovMatrix/numOfSamples;
		cout << meansVector/numOfSamples;
		char *s = "./tryRES/IBDCovMatrix"; 
		ofstream outfile;
		outfile.open(s);
		outfile << (IBDCovMatrix/numOfSamples - meansVector*meansVector.transpose()/numOfSamples/numOfSamples);
		
	}
	else if(whatToCompute=="pdq"){
		pop->output_pdq(numOfSamples,"pdq.out");
	}
	if (pop->samplePosForMe!=ULONG_MAX){
	   h1.plot(20,"plot");
	   h1.P.save("estPosition.ps");
	}
}
/*! \fn void Model::RSamplerMH(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType="genotypic",const std::string &whatToCompute="genotypeProb")
*  \brief method to sample ordered genotypes using the Metropolis Hastings algorithm to accept samples proposed by sequential imputation.
*/

void Model::RSamplerGibbsMH(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,unsigned stepGibbsMH, const std::string &samplerType,const std::string &whatToCompute){
	// Authors: L. Radu Totir 
	// (July, 2004) 
	// Contributors: 
	// if (pedBlockSize!=pop->popsize) throw exception("Model::RSamplerGibbs(...): size of pedigree block must be equal to the population size");
	double maxNumGNodeElements = 4.0;
	int    maxDimensionCutSet  = 7;
	unsigned maxsize = (unsigned)(std::pow((maxNumGNodeElements),(maxDimensionCutSet))+.5);
	matvec::Matrix<double> IBDCovMatrix, meansVector;
	IBDCovMatrix.resize(1,1,0.0);
	int    badCount  = 0;
	int    goodCount = 1; // 1 because first sample is always accepted
	double alpha      = 1.0;
	double ranNumber  = 0.0;
	double logOldProposal, logNewProposal, logOldTarget, logNewTarget;
	pop->gNodeList.logProposal = 0.0;
	pop->gNodeList.logTarget   = 0.0;
	unsigned stopBlock=0,startBlock=0, count=1;
	for (int i=0;i<numOfSamples;i++){
		cout << "Sample " << setw(5) << i+1 << endl; 
		if (i==0){
			// get initial sample
			pop->getInitialGNodeListSample(maxsize,0,pop->popsize,pedBlockSize,samplerType);
			logOldTarget = pop->gNodeList.logTarget;
			logOldProposal = pop->gNodeList.logOldProposal;
			pop->copyGNodeStatesToCandidateStates(samplerType);
			pop->storeSampledGametes();
			pop->copyCandidateToAccepted(samplerType);
			if(whatToCompute=="haplotypeFreq"){
				pop->countHaplotypes(samplerType);
			}
			else if(whatToCompute=="probDescHaplo"){
				pop->SetFreqHaploFounders();
			}
			else if(whatToCompute=="ibdCovMatrix"){
				IBDCovMatrix.resize(pop->popsize,pop->popsize,0.0);
				meansVector.resize(pop->popsize ,1           ,0.0);
			}
		}
		// repeat initial MH sampling procedure after N samples 
		else if (!(i%stepGibbsMH)){
			pop->copyGNodeStatesToCandidateStates(samplerType);
			pop->storeSampledGametes();
			// q(NEW) and pi(NEW)
			pop->gNodeList.logProposal = 0.0;
			pop->gNodeList.logTarget   = 0.0;
			pop->getInitialGNodeListSample(maxsize,0,pop->popsize,pedBlockSize,samplerType);
			pop->copyGNodeStatesToCandidateStates(samplerType);
			pop->storeSampledGametes();
			logNewProposal = pop->gNodeList.logProposal;
			logNewTarget   = pop->gNodeList.logTarget;
			pop->gNodeList.logOldProposal = 0.0;
			pop->getOldGNodeListProbability(maxsize,0,pop->popsize,pedBlockSize,samplerType);
			logOldProposal = pop->gNodeList.logOldProposal;
			cout << "logNewTarget   = " << logNewTarget   << endl;
			cout << "logOldProposal = " << logOldProposal << endl;
			cout << "logOldTarget   = " << logOldTarget   << endl;
			cout << "logNewProposal = " << logNewProposal << endl;
			alpha=std::exp(logNewTarget + logOldProposal - logNewProposal - logOldTarget);
			cout << " alpha = " << alpha << endl;
			ranNumber=ranf();
			if (ranNumber <= alpha) {
				logOldTarget = logNewTarget;
				pop->copyCandidateToAccepted(samplerType);
				goodCount++;
			}
			else {
				pop->retreiveSampledGametes();
				badCount++;
			}
			if(whatToCompute=="probDescHaplo"){
				pop->UpdateFreqHaploFounders();
				pop->sampleSegregationIndicators();
				// switch for 1 column in all the declarations
				// pop->resizeSegregationIndicators(); 
			}
			else if(whatToCompute=="ibdCovMatrix"){
				IBDCovMatrix += pop->getIBDMatrix();
			    meansVector  += pop->getMeans();
			}
			cout << "Got " <<  goodCount << " valid MH samples of " << count++ << " tries." << endl;
		}
		// do Gibbs across loci
		else{
			pop->gNodeList.logProposal = 0.0;
			pop->gNodeList.logTarget   = 0.0;
			do{ 
				if(stopBlock < pop->popsize){
					stopBlock += pedBlockSize;
					startBlock = stopBlock - pedBlockSize;
				}
				else {
					stopBlock  = pedBlockSize;
					startBlock = stopBlock - pedBlockSize;
				}
				// cout << "startBlock = " << startBlock << endl;
				// cout << "stopBlock  = " << stopBlock  << endl;
				pop->getGNodeListSample(maxsize,startBlock,stopBlock,pedBlockSize,samplerType);
				// keep track of the target probability for potential MH step
				logOldTarget   = pop->gNodeList.logTarget;
				pop->copyGNodeStatesToCandidateStates(samplerType);
				pop->storeSampledGametes();
				pop->copyCandidateToAccepted(samplerType);
				if(whatToCompute=="probDescHaplo"){
					pop->UpdateFreqHaploFounders();
				}
			} while(stopBlock < pop->popsize);
			if(whatToCompute=="probDescHaplo"){
				pop->sampleSegregationIndicators();
				// switch for 1 column in all the declarations
				// pop->resizeSegregationIndicators(); 
			} 
			else if(whatToCompute=="ibdCovMatrix"){
				IBDCovMatrix += pop->getIBDMatrix();
				meansVector  += pop->getMeans();
			}
		}
	}
	if(whatToCompute=="probDescHaplo"){
		pop->CalcFreqHaploFounders(numOfSamples);
		pop->DisplayFreqHaploFounders();
	}
	else if(whatToCompute=="ibdCovMatrix"){
		char *s = "./tryRES/IBDCovMatrix"; 
		ofstream outfile;
		outfile.open(s);
		outfile << (IBDCovMatrix/numOfSamples - meansVector*meansVector.transpose()/numOfSamples/numOfSamples);
	}
}

void Model::variance(const std::string &termname,Pedigree& P,const double v00,...)
  // *************************************************************************
  // * variance is the method to get variance components for each random
  // * effect.
  // * the example usage is below:
  // *
  // *    M.variance("animal",P,5.0)
  // *********************************************************************
{

   if (numtrait<1) throw exception("Model::variance(args): no trait(dependent-variable) in the model");
   int t1,t2;
   doubleMatrix var(numtrait,numtrait);
   va_list param_pt;                      // param_pt, an object of va_list
   va_start(param_pt,v00);                // call the setup macro
   var[0][0] = v00;
   for (t2=1; t2<numtrait; t2++) var[0][t2] = va_arg(param_pt,double);
   for (t1=1; t1<numtrait; t1++) for (t2=0; t2<numtrait; t2++) {
       var[t1][t2] = va_arg(param_pt,double);
   }
   va_end(param_pt);
   variance(termname,P,var);
}

void Model::variance(const std::string &termname,Pedigree& P, const doubleMatrix& v)
  // **********************************************************************
  // * variance is the method to get variance components for each random
  // * effect.
  // *********************************************************************
{
   num_ped++;
   if (!factor_struct) {
      std::cerr << "Model::variance(args): no equation(s) in the model\n";
      return;
   }

   int k = term.index(termname,factor_struct);



   if ( k < 0) {
      std::cerr << "Model::variance(): term not in the model\n";
      return;
   }
   //   if (term[k].order() > 1) {
   //      std::cerr << "Model::variance(): Pedigree can only link to single effect\n";
   //      return;
   //   }
   if (term[k].classi() != 'F')  {     // term[k] classi is being changing
     data = 0;
     //   for (int i=0;i<term[k].order();i++) {
     //     factor_struct[term[k].factorindx[i]].classi('F');
     //   }
   }
   term[k].classi('P');                    // factor classi will override
   type = mixed_model;    factor_struct[term[k].factorindx[0]].classi('P'); // data column classi later
   if (pop_created) delete pop;
   pop = new Population;
   check_ptr(pop);
   pop_created = 1;
   term[k].prior->resize(numtrait);     // prior must be UnknownDist
   if (pop->input_ped(P,term[k].prior)) {
      type = bad_model;
      return;
   }
   numgroup = P.ngroup();
   popsize = pop->size();
   if (v.num_rows() == v.num_cols()/* && v.num_rows() == numtrait*/) {
      if (v.psd()) {
         term[k].prior->var_matrix()->copy(v);
      }
      else {
         std::cerr << "Model::variance(): not positive(semi) def: " << termname << "\n";
      }
   }
   else {
      std::cerr << "Model::variance(): invalid variance matrix size\n";
   }
}

void Model::variance(const std::string &termname,const doubleMatrix& v)
  // *********************************************************************
  // * variance is the method to get variance components for each random
  // * term in a model.
  // *********************************************************************
{
   if (!factor_struct) {
      std::cerr << "Model::variance(args): no equation(s) in the model\n";
      return;
   }
   if (v.num_rows() != v.num_cols() || v.num_rows() != numtrait) {
      std::cerr << "Model::variance(): invalid variance matrix size\n";
      return;
   }
   if (!v.psd()) {
      std::cerr << "Model::variance:: not positive(semi) def: " << termname << "\n";
      return;
   }
   if (termname == "residual") {
      residual_var = v;
      return;
   }
   int i,k;
   k = term.index(termname,factor_struct);
   if ( k >= 0) {
      if (term[k].classi() != 'F') {     // term[k] classi is being changing
         data = 0;
	 //         for (i=0;i<term[k].order();i++) {
	 // factor_struct[term[k].factorindx[i]].classi('F');
	 //         }
      }
      term[k].classi('R');                    // factor classi will override
      type = mixed_model;
      term[k].prior->resize(numtrait);
      term[k].prior->var_matrix()->copy(v);
   }
   else {
      std::cerr << "Model::variance(..): not in model: " << termname << "\n";
   }
}

void Model::variance(const std::string &termname,const double v00,...)
  // *************************************************************************
  // * variance is the method to get variance components for each random
  // * term in a model.
  // * the example usage is below:
  // *
  // *    M.variance("animal",1.0,...)
  // *********************************************************************
{
   if (numtrait<1) {
      std::cerr << "Model::variance(args): no trait(dependent-variable) in the model\n";
      return;
   }
   doubleMatrix var(numtrait,numtrait);
   int t1,t2;
   va_list param_pt;                       // an object param_pt
   va_start(param_pt,v00);                // call the setup macro
   var[0][0] = v00;
   for (t2=1; t2<numtrait; t2++) var[0][t2] = va_arg(param_pt,double);

   for (t1=1; t1<numtrait; t1++) for (t2=0; t2<numtrait; t2++) {
      var[t1][t2] = va_arg(param_pt,double);
   }
   va_end(param_pt);
   variance(termname,var);
}




////SDK


void Model::VarLink(const std::string &termname,const std::string &termname2)

  // *********************************************************************
  // * VarLink is the method to get covariancelink a set of variance
  // * components  in a model.
  // *********************************************************************
{
   if (!factor_struct) {
      std::cerr << "Model::VarLink(args): no equation(s) in the model\n";
      return;
   }

   int i,k,k2;
   k = term.index(termname,factor_struct);

   k2 = term.index(termname2,factor_struct);

   if( k < 0 ){
      std::cerr << "Model::VarLink(..): not in model: " << termname << "\n";
      return;
   }
   if( k2 < 0 ){
      std::cerr << "Model::VarLink(..): not in model: " << termname2 << "\n";
      return;
   }
   int bek=base_effect[k];
   if(term[bek].classi() != 'R' && term[bek].classi() != 'P'){
      std::cerr << "Model::VarLink(..): not a random effect: " << termname << "\n";
      return;
   }
   int bek2=base_effect[k2];
   if(term[bek2].classi() != 'R' && term[bek2].classi() != 'P'){
      std::cerr << "Model::VarLink(..): not a random effect: " << termname2 << "\n";
      return;
   }
   if(term[bek].classi() != term[bek2].classi() ){
     std::cerr << "Model::VarLink(..): different types of random effects: \n";
      return;
   }



   int be=base_effect[k];
   int be2=base_effect[k2];
   if(var_link[be]!=var_link[be2]){
     int vl=var_link[be];
     int vl2=var_link[be2];
     if(vl2 <vl) {
       int tmp;
       tmp=vl;
       vl=vl2;
       vl2=tmp;
     }
     var_link[vl2]=vl;  
   }
}


void Model::covariance(const std::string &termname,const std::string &termname2,const doubleMatrix& v)
  // *********************************************************************
  // * covariance is the method to get covariance components for a set
  // * of correlated  random terms in a model.
  // *********************************************************************
{
   if (!factor_struct) {
      std::cerr << "Model::covariance(args): no equation(s) in the model\n";
      return;
   }
   if (v.num_rows() != v.num_cols() || v.num_rows() != numtrait) {
      std::cerr << "Model::covariance(): invalid variance matrix size\n";
      return;
   }
   int i,k,k2;
   k = term.index(termname,factor_struct);

   k2 = term.index(termname2,factor_struct);

   if( k < 0 ){
      std::cerr << "Model::covariance(..): not in model: " << termname << "\n";
      return;
   }
   if( k2 < 0 ){
      std::cerr << "Model::covariance(..): not in model: " << termname2 << "\n";
      return;
   }
   int bek=base_effect[k];
   if(term[bek].classi() != 'R' && term[bek].classi() != 'P'){
      std::cerr << "Model::covariance(..): not a random effect: " << termname << "\n";
      return;
   }
   int bek2=base_effect[k2];
   if(term[bek2].classi() != 'R' && term[bek2].classi() != 'P'){
      std::cerr << "Model::covariance(..): not a random effect: " << termname2 << "\n";
      return;
   }
   if(term[bek].classi() != term[bek2].classi() ){
     std::cerr << "Model::covariance(..): different types of random effects: \n";
      return;
   }


   int be;
   if(base_effect[k]==base_effect[k2]){  
     be=base_effect[k];
   }
   else {
     be=base_effect[k];
     int be2=base_effect[k2];
     if(be2 <be) {
       int tmp;
       tmp=k;
       k=k2;
       k2=tmp;
       tmp=be;
       be=be2;
       be2=tmp;
     }
     base_effect[k2]=be;
     doubleMatrix orig1,orig2;
     orig1.copy(*term[be].prior->var_matrix());
     orig2.copy(*term[be2].prior->var_matrix());
     int newsize=nt_vec[be]+nt_vec[be2];
     Vector<int> orig=trait_map[be];
     trait_map[be].reserve(newsize);
     int ii=0;
     for(int i=0;i<nt_vec[be];i++,ii++) trait_map[be][ii]=orig[i];
     for(int i=0;i<nt_vec[be2];i++,ii++) trait_map[be][ii]=trait_map[be2][i];
     term[be].prior->var_matrix()->resize(newsize,newsize);

     for(int i=0;i<nt_vec[be];i++){
       for(int j=0;j<nt_vec[be];j++){
	 term[be].prior->var_matrix()->me[i][j]=orig1[i][j];
       }
     }

     for(int i=0;i<nt_vec[be2];i++){
       for(int j=0;j<nt_vec[be2];j++){
	 term[be].prior->var_matrix()->me[nt_vec[be]+i][nt_vec[be]+j]=orig2[i][j];
       }
     }
     nt_vec[be]=newsize;
     term[be2].prior->var_matrix()->clear();
     nt_vec[be2]=0;
     term[be2].classi('F');
     int last_eff=be;
     for(;corr_link[last_eff] != -1;)last_eff=corr_link[last_eff];
     corr_link[last_eff]=be2;
     int last_pos=pos_vec[last_eff];
     for(;corr_link[last_eff] != -1;){
       last_eff=corr_link[last_eff];
       pos_vec[last_eff]=pos_vec[last_eff]+last_pos+1;
     }
   }

   int srow,scol;
   srow=pos_vec[k]*numtrait;
   scol=pos_vec[k2]*numtrait;
   for(int irow=0;irow<numtrait;irow++){
     for(int jcol=0;jcol<numtrait;jcol++){
       term[be].prior->var_matrix()->me[srow+irow][scol+jcol]=v[irow][jcol];
       term[be].prior->var_matrix()->me[scol+jcol][srow+irow]=v[irow][jcol];
     }
   }      


   //   std::cout << "Cov " << *term[be].prior->var_matrix();
   

}

void Model::covariance(const std::string &termname,const std::string &termname2,const double v00,...)
  // *************************************************************************
  // * variance is the method to get variance components for each random
  // * term in a model.
  // * the example usage is below:
  // *
  // *    M.variance("animal",1.0,...)
  // *********************************************************************
{
   if (numtrait<1) {
      std::cerr << "Model::variance(args): no trait(dependent-variable) in the model\n";
      return;
   }
   doubleMatrix var(numtrait,numtrait);
   int t1,t2;
   va_list param_pt;                       // an object param_pt
   va_start(param_pt,v00);                // call the setup macro
   var[0][0] = v00;
   for (t2=1; t2<numtrait; t2++) var[0][t2] = va_arg(param_pt,double);

   for (t1=1; t1<numtrait; t1++) for (t2=0; t2<numtrait; t2++) {
      var[t1][t2] = va_arg(param_pt,double);
   }
   va_end(param_pt);
   covariance(termname,termname2,var);
}









////SDK


void Model::covariate(const std::string &covariate_names)
{
  // *********************************************************************
  //  it specifies a list of effect (factor) as covariates.
  //  note that an effect is not a model-term, which consists of effects,
  //  instead an effect is associated with a data column. It is required
  //  that the name of an effect must be the same as the associated
  //  data-column.
  // **********************************************************************
   if (!factor_struct) {
      std::cerr << "Model::covariate(): no equation(s) in the model\n";
      return;
   }

   int j,k,t;
   unsigned i,nc;
   std::vector<std::string> efflist;
   nc = split(covariate_names," ,",&efflist);  // efflist = covariate_names.split(nc,", ");

   for (i=0; i<nc; i++) {
      k = fieldindex(efflist[i],factor_struct,numfact);
      if (k < 0) throw exception("Model::covariate(): not in the model, thus it's ignored");
       // need to check if classi is changed
      if (factor_struct[k].classi() != 'F') data= 0;
      factor_struct[k].classi('C');

      ////////////////////////////////////////////////////////////////////
      // if a term contains at least one covariate, then its classi = 'C'
      // else leave it alone
      //////////////////////////////////////////////////////////////////
      for (t=0; t<numterm; t++) {
         j = term[t].partial_match(efflist[i].c_str(),factor_struct);
         if (j >=0 ) {
            term[t].classi('C');
            break;
         }
      }
   }
}

static void nest_xaction(std::string &s)
{
   int i,j;
   int n = s.size();
   for (i=0; i<n; i++) {
      if ( s[i] == '(' ) s[i] = '*';
      if ( s[i] == ')' ) s[i] = ' ';
   }

      // replace space around * with *
   i = 0;
   while (i<n) {
      if (s[i] == '*' ) {
         j = i;
         while (s[--j] == ' ') s[j]='*';
         while (s[++i] == ' ') s[i]='*';
      }
      i++;
   }
}

int Model::build_model_term(const std::string &modelspec)
{
   int retval = 1;
   unsigned i,j,k;

   std::string mspec = modelspec;  // modelspec:  y1 = A  B + A (B), y2 = B  + X*B
   nest_xaction(mspec);           // now mspec:  y1 = A  B + A**B , y2 = B  + X*B

   std::string sep("=+,* ");
   std::string TMP = mspec;
   ////////   TMP.unique("=+,* ");                // TMP ->y1 y2 A B X
   std::vector<std::string> tmpvec;
   std::vector<std::string>::iterator it;
   std::string::size_type begidx;
   split(TMP,sep,&tmpvec);

   sort(tmpvec.begin(),tmpvec.end());
   it = unique(tmpvec.begin(),tmpvec.end());
   tmpvec.erase(it,tmpvec.end());   
//    for(it = tmpvec.begin(); it != tmpvec.end(); it++){
//      std::cout << *it << endl;
//    }

   categorical = 0;
   numfact =  tmpvec.size();
   if(factor_struct){
     delete [] factor_struct;
     factor_struct = 0;
   }
   if(numfact>0){
     factor_struct = new FieldStruct [numfact];
   }
   else {
     factor_struct = 0;
   }
   check_ptr(factor_struct);
   for (i=0; i<numfact; i++) factor_struct[i].name(tmpvec[i]);

   TMP = mspec;
   numtrait = 0;
   sep = "=";
   begidx = TMP.find(sep,0);
   while(begidx != std::string::npos) {
      numtrait++;
      begidx = TMP.find(sep,begidx+1);
   }
   if (numtrait == 0) {
      std::cerr << "Model::build_model_term(): there is no trait in the model\n";
      return 0;
   }
   if(trait_struct){
     delete [] trait_struct;
     trait_struct = 0;
   }
   if(numtrait>0){
     trait_struct = new FieldStruct* [numtrait];
   }
   else {
     trait_struct = 0;
   }
   check_ptr(trait_struct);
   if(traitname){
     delete [] traitname;
     traitname = 0;
   }
   if(numtrait>0){
     traitname = new std::string [numtrait];
   }
   else {
     traitname = 0;
   }
   check_ptr(traitname);
   std::vector<std::string> model_trait;
   std::vector<std::string> trtnm;
   sep = ",";
   split(TMP,sep,&model_trait);
   if (numtrait !=  model_trait.size()) throw exception("Model::build_model_term(): model for each trait is not separated properly");

   sep = " +";
   unsigned n = 0;
   for (i=0; i<numtrait; i++) {
     // trait name should be before "="
     begidx = model_trait[i].find("=");
     // bbbtrtnmbb is the trait name with possible surrounding blanks
     std::string bbbtrtnmbbb = model_trait[i].substr(0,begidx);
     // getting rid of the surrounding blanks
     split(bbbtrtnmbbb,sep,&trtnm);
     traitname[i]  = trtnm[i];
     model_trait[i]  = model_trait[i].substr(begidx+1);
     n += split(model_trait[i],sep);  /// n += model_trait[i].ntoken("+ ");
     k = fieldindex(traitname[i],factor_struct,numfact);
     if (k >= 0 && k < numfact) {
       factor_struct[k].classi('T');
       trait_struct[i] = &(factor_struct[k]);
     }
     else {
       throw exception("Model::build_model_term(): you have probably found a bug!");
     }
   }
   TermList T(n,numtrait);
   tmpvec.clear();
   sep = " +";
   for (numterm=0,i=0; i<numtrait; i++) {
      tmpvec.clear();
      n = split(model_trait[i],sep,&tmpvec);  //   tmpvec = model_trait[i].split(n,"+ ");
      for (k=0; k<n; k++) {
         for (j=0; j<numterm; j++) {
            if (T[j].match(tmpvec[k],"* ",factor_struct)) break;
         }
         if (j==numterm) {
            T[numterm++].input(tmpvec[k].c_str(),factor_struct,numfact);
         }
         T[j].trait[i]= 1;
      }
   }

   term.resize(numterm,numtrait);
   if (unknown_prior) {
     delete [] unknown_prior;                    // note that unknown_prior[i]
     unknown_prior = 0;
   }
   if(numterm>0){
     unknown_prior = new UnknownDist[numterm];     // needs to resize(numtrait)
   }
   else {
     unknown_prior = 0;
   }
   check_ptr(unknown_prior);
   nt_vec.resize(numterm,numtrait);
   base_effect.reserve(numterm);
   corr_link.resize(numterm,-1);
   pos_vec.resize(numterm);
   trait_map.reserve(numterm);
   var_link.reserve(numterm);
   

   for (i=0; i<numterm; i++) {
     base_effect[i]=i;
     var_link[i]=i;
      term[i] = T[i];
      term[i].prior = &(unknown_prior[i]);
      trait_map[i].resize(numtrait,i);
   }
   return retval;
}

void Model::trait_effect_level(const unsigned mmeaddr, std::string& retval,
                               int info_flag)
{
   // ******************************************************************
   // mmeaddr = [0, hmmesize-1]
   // info_flag = 1 information: (traitname):termname:levelname(default)
   //             2                          termname:levelname
   //             3                                   levelname
   // ******************************************************************
   if (mmeaddr >= hmmesize) throw exception("Model::trait_effect_level(): arg out of range");
   retval = "";
   char CL,strid[11];
   std::string tempstr;
   int nord,*tempval;
   unsigned addr,i,j,k,t,kk,ii,jj,nlevel;
   unsigned *tempuns;
   Vector<unsigned> unsvec(maxorder);

   for (addr=0, i=0; i<numterm; i++) {
      nlevel = term[i].nlevel();
      addr += nlevel*numtrait;
      if (addr > mmeaddr) break;
   }
   if (i >= numterm) throw exception("Model::trait_effect_level(): arg out of range");
   addr = term[i].startaddr();
   nord = term[i].order();
   for (k=addr,j=0; j<nlevel; j++) {
      k += numtrait;
      if (k > mmeaddr) break;
   }
   for (addr= k-numtrait, t=0; t<numtrait; t++) {
      if (mmeaddr == addr++) break;
   }
   if (nord == 1) {
      unsvec[0] = j+1;
   }
   else {
      tempuns = (unsigned *)(xact_htable[i].find(j+1));
      for (ii=0; ii<nord; ii++) unsvec[ii] = tempuns[ii];
   }


   if (term[i].addr[t] == 0) throw exception("Model::trait_effect_level(): not defined address");
   if (numtrait>1 && info_flag == 1) {
      retval = trait_struct[t]->name();
     retval.append(":");
   }

   if (info_flag == 1 || info_flag == 2) {
      for (ii=0; ii<nord; ii++) {
         retval.append(factor_struct[term[i].factorindx[ii]].name());
         if (ii+1 < nord) retval.append("*");
      }
      retval.append(":");
   }
   for (ii=0; ii<nord; ii++) {
      kk = factor_struct[term[i].factorindx[ii]].index();
      CL = data->datasheet[kk].classi();
      jj = unsvec[ii];

      if (CL == 'F') {
         tempval = (int *)(idlist[kk]->find(jj));
         sprintf(strid,"%d",*tempval);
         tempstr = strid;
      }
      else if (CL =='U') {
         tempstr = (const char *)idlist[kk]->find(jj);
      }
      else if ( CL =='C') {
         tempstr = data->datasheet[kk].name();
      }
      else if (CL =='P') {
         tempstr = (const char *)pop->ind_name(jj);
      }
      else if (CL =='I') {
         tempstr = "1";
      }
      else {
         warning("Model::save(): unknown column type: %c",CL);
         break;
      }
      retval.append(tempstr);
      if (ii+1 < nord) retval.append("*");
   }
}

Vector<double> Model::get_solutions(const std::string &termname)
{
  // Authors: Liviu R. Totir and Rohan L. Fernando (November, 2002) 
  // Contributors: 
  Vector<double> retval;
  if (!get_blupsol()) throw exception("Model::get_solutions(): get_blupsol() failed");
  
  int k = term.index(termname,factor_struct);
  if (k < 0) throw exception("Model::get_solutions(): no such term in the model");
  int nlevels = term[k].nlevel()*numtrait;
  unsigned i,startaddr;
  startaddr = term[k].startaddr();
  retval.resize(nlevels);
  for (i=0; i<nlevels; i++) {
    retval[i] = blupsol[startaddr+i];
  }
  return retval;
}

void Model::save(const std::string &fname, const int io_mode)
{
   if (type == bad_model) {
      warning("Model::save(): model is too bad, nothing is saved");
      return;
   }
   if((numterm+ntermGdist) == 0) return;
   if (blupsol.size() != hmmesize + ntermGdist*popsize) return;
   if (!data) {
      warning("Model::save(): no data has been fitted");
      return;
   }
   std::ofstream ofs;
   ofs.open(fname.c_str(),(OpenModeType)io_mode);
   if (!ofs) {
      warning("Model::save: %s cann't open or already exit",fname.c_str());
      return;
   }
   ofs.setf(std::ios::unitbuf | std::ios::showpoint);
   ofs.precision(SESSION.output_precision);
   int W = SESSION.output_precision + 6;        // 6 = +.e+00
   char S = ' ';                                // one blank space 
   std::string SS;
   SS = "          ";

   if (!data->in_memory()) data->input_datasheet();
   this->info(ofs);

   double *ve = blupsol.begin();

   char* tempstr;
   int nord, *tempval;
   unsigned *tempuns;
   unsigned i,j,t,k,ii,kk,jj,nlevel,total_numterm;
   char CL;
                   
   ofs << "\n            BLUP (BLUE, Posterior_Mean) \n";
   ofs << "    ----------------------------------------------- \n\n";
   total_numterm = numterm+ntermGdist;
   for (i=0; i<total_numterm; i++) {   
      ofs << S << std::setw(W)<< "Term";
      for (t=0; t<numtrait; t++) {
         ofs << S << std::setw(W-2) << "Trait" << S << t+1;
      }
      ofs <<"\n";
      
      ofs << S << std::setw(W) << factor_struct[term[i].factorindx[0]].name().c_str();
      for (t=1; t<term[i].order(); t++) {
         ofs << "*" << factor_struct[term[i].factorindx[t]].name().c_str();
      }
      for (t=0; t<numtrait; t++) ofs <<S <<std::setw(W)<<trait_struct[t]->name().c_str();
      ofs << "\n";
      nlevel = term[i].nlevel();
      k = term[i].startaddr();
      nord = term[i].order();
      if (nord == 1) {      // single_factor_term
         kk = factor_struct[term[i].factorindx[0]].index();
         CL = data->datasheet[kk].classi();
         switch(CL) {
            case 'F':
               for (j=0; j<nlevel; j++) {
                  tempval = (int *)(idlist[kk]->find(j+1));
                  ofs << S << std::setw(W) << *tempval;
                  for (t=0; t<numtrait; t++,k++) {
                     ofs << S << std::setw(W);
                     if (term[i].trait[t]) {
                        ofs << ve[k];
                     }
                     else {
                        ofs << SS;
                     }
                  }
		  k+=(nt_vec[base_effect[i]]-numtrait);
                  ofs << "\n";
               }
               break;
            case 'U':
               for (j=0; j<nlevel; j++) {
                  tempstr = (char *)(idlist[kk]->find(j+1));
                  ofs << S << std::setw(W) << tempstr;
                  for (t=0; t<numtrait; t++,k++) {
                     ofs << S << std::setw(W);
                     if (term[i].trait[t]) {
                        ofs << ve[k];
                     }
                     else {
                        ofs << SS;
                     }
                  }
                  ofs << "\n";
               }
               break;
            case 'P':
               for (j=0; j<nlevel; j++) {
                  ofs << S;
                  if (pop->maxnamelen > 0) {
                     ofs << std::setw(W) << (char *)pop->ind_name(j+1);
                  }
                  else {
                     ofs << std::setw(W) << j+1;
                  }
                  for (t=0; t<numtrait; t++,k++) {
                     ofs << S;
                     if (term[i].trait[t]) {
                        ofs << std::setw(W) << ve[k];
                     }
                     else {
                        ofs << std::setw(W) << SS;
                     }
                  }
		  k+=(nt_vec[base_effect[i]]-numtrait);
                  ofs << "\n";
               }
               break;
            case 'C':
               ofs << S << std::setw(W) << data->datasheet[kk].name();
               for (t=0; t<numtrait; t++,k++) {
                  ofs << S << std::setw(W);
                  if (term[i].trait[t]) {
                     ofs << ve[k];
                  }
                  else {
                     ofs << SS;
                  }
               }
               ofs << "\n";
               break;
            case 'I':
               ofs << S << std::setw(W) << "1";
               for (t=0; t<numtrait; t++,k++) {
                  ofs << S << std::setw(W);
                  if (term[i].trait[t]) {
                     ofs << ve[k];
                  }
                  else {
                     ofs << SS;
                  }
               }
               ofs << "\n";
               break;
            default:
               warning("Model.save(): unknown column type: %c",CL);
               break;
         }
         if (i+1 < total_numterm) ofs << "\n";
      }
      else {     // there are interactions
         for (j=0; j<nlevel; j++) {
            tempuns = (unsigned *)(xact_htable[i].find(j+1));
            kk = factor_struct[term[i].factorindx[0]].index();
            CL = data->datasheet[kk].classi();
            jj = tempuns[0];

            if (CL == 'F') {
               tempval = (int *)(idlist[kk]->find(jj));
               ofs << S << std::setw(W) << *tempval;
            }
            else if (CL =='U') {
               tempstr = (char *)(idlist[kk]->find(jj));
               ofs << S << std::setw(W) << tempstr;
            }
            else if ( CL =='C') {
               ofs << S << std::setw(W) << data->datasheet[kk].name();
            }
            else if (CL =='P') {
               ofs << S << std::setw(W) << (char *)pop->ind_name(jj);
            }
            else if (CL =='I') {
               ofs << S << std::setw(W) << "1";
            }
            else {
               warning("Model::save(): unknown column type: %c",CL);
               break;
            }

            for (ii=1; ii<nord; ii++) {
               kk = factor_struct[term[i].factorindx[ii]].index();
               CL = data->datasheet[kk].classi();
               jj = tempuns[ii];

               if (CL == 'F') {
                  tempval = (int *)(idlist[kk]->find(jj));
                  ofs << "*" << *tempval;
               }
               else if (CL =='U') {
                  tempstr = (char *)(idlist[kk]->find(jj));
                  ofs << "*" << tempstr;
               }
               else if (CL =='C' || CL =='I') {
                  ofs << "*" << data->datasheet[kk].name();
               }
               else if (CL =='P') {
                  ofs << "*" << (char *)pop->ind_name(jj);
               }
               else {
                 warning("Model::save(): unknown column type: %c",CL);
               }
            }
            for (t=0; t<numtrait; t++,k++) {
               ofs << S << std::setw(W);
               if (term[i].trait[t]) ofs << ve[k];
               else ofs << SS;
            }
	    k+=(nt_vec[base_effect[i]]-numtrait);
            ofs << "\n";
         }
         if (i+1 < total_numterm) ofs << "\n";
      }
   }
   ofs.close();
   data->release_datasheet();
}

double Model::minfun(const Vector<double> &x, const int n)
{
   double llhood = 0.0;
   switch (minfun_indx) {
      case 1:                //  dfreml vce,   Model_vce.cc
         llhood = minfun_vce(x,n);
         break;
      case 2:                // MLE peeling,  Model_peeling
         llhood = minfun_peeling(x,n);
         break;
// BRS
      case 3:   // Multipoint 
         llhood = minfun_multipoint(x,n);
         break;
// BRS
      default:
         warning("Model::Model::minfun(): unknown minfun_indx: %d",minfun_indx);
   }
   return llhood;
}

void Model::RSamplerMH(unsigned numLoci,unsigned numOfSamples,const std::string &samplerType,const std::string &whatToCompute){
	// Authors: L. Radu Totir and Rohan L. Fernando 
	// (October, 2003) 
	// Contributors: 
	
	unsigned maxsize = myRSamplerParms.maxCutsetSize;
	// cout << maxsize << endl; 
	int    badCount  = 0;
	int    goodCount = 1; // 1 because first sample is always accepted
	double alpha      = 1.0;
	double ranNumber  = 0.0;
	double logOldProposal, logNewProposal, logOldTarget, logNewTarget;
	pop->gNodeList.logProposal = 0.0;
	pop->gNodeList.logTarget   = 0.0;
	matvec::Matrix<double> IBDCovMatrix, meansVector;
	IBDCovMatrix.resize(1,1,0.0);
	for (int i=0;i<numOfSamples;i++){
		if(i==0){
			// get initial sample
			//cout << "Sample " << 1 << endl;
			pop->getInitialGNodeListSample(maxsize,numLoci,samplerType);
			logOldTarget = pop->gNodeList.logTarget;
			logOldProposal = pop->gNodeList.logOldProposal;
			pop->copyGNodeStatesToCandidateStates(samplerType);
			pop->storeSampledGametes();
			pop->copyCandidateToAccepted(samplerType);
			if(whatToCompute=="probDescHaplo"){
				pop->SetFreqHaploFounders();
			}
			else if(whatToCompute=="ibdCovMatrix"){
				IBDCovMatrix.resize(pop->popsize,pop->popsize,0.0);
				meansVector.resize(pop->popsize ,1           ,0.0);
			}
			else if(whatToCompute=="genotypeFreq"){
				pop->initGenotypeFreq();
			}
		}
		else{
			// q(NEW) and pi(NEW)
			pop->gNodeList.logProposal = 0.0;
			pop->gNodeList.logTarget   = 0.0;
			pop->getGNodeListSample(maxsize,numLoci,samplerType);
			pop->copyGNodeStatesToCandidateStates(samplerType);
			pop->storeSampledGametes();
			logNewProposal = pop->gNodeList.logProposal;
			logNewTarget   = pop->gNodeList.logTarget;
			if (!(i%100)){
				cout << "logNewTarget   = " << logNewTarget   << endl;
				cout << "logOldTarget   = " << logOldTarget   << endl;
				cout << "logNewProposal = " << logNewProposal << endl;
				cout << "logOldProposal = " << logOldProposal << endl;
			}
			// accept or reject the new sample
			alpha=std::exp(logNewTarget + logOldProposal - logNewProposal - logOldTarget);
			if (!(i%100)){
				cout << " alpha = " << alpha << endl;
			}
			ranNumber=ranf();
			if (ranNumber <= alpha) {
				logOldTarget = logNewTarget;
				logOldProposal = logNewProposal;
				pop->copyCandidateToAccepted(samplerType);
				goodCount++;
			}
			else {
				pop->storeSampledGametes();
				badCount++;
				//cout << "badCount = " << badCount << endl;
			}
			//cout << "goodCount = " << goodCount << endl;
			if(whatToCompute=="probDescHaplo"){
				pop->UpdateFreqHaploFounders();
			}
			if (!(i%100)){
				cout << "No of samples accepted = " << goodCount << " out of a total of " << goodCount+badCount << endl;
				double ratio = double(badCount)/(goodCount+badCount);
				// cout << "Rejection rate = " << ratio << endl;
			}
		}
		if (whatToCompute=="haplotypeFreq"){
			pop->countHaplotypes(samplerType);
		}
		else if (whatToCompute=="genotypeFreq"){
			pop->countGenotypes(samplerType);
		}
		else if(whatToCompute=="probDescHaplo"){
			pop->sampleSegregationIndicators();
		} 
		else if(whatToCompute=="ibdCovMatrix"){
			IBDCovMatrix += pop->getIBDMatrix();
			meansVector  += pop->getMeans();
		}
	}
	cout << "No of samples accepted = " << goodCount << endl;
	double ratio = double(badCount)/(goodCount+badCount);
	cout << "Rejection rate = " << ratio << endl;
	if (whatToCompute=="haplotypeFreq"){
		pop->displayHaplotypeFrequencies(numOfSamples);
	}
	else if(whatToCompute=="probDescHaplo"){
		pop->CalcFreqHaploFounders(numOfSamples);
		pop->DisplayFreqHaploFounders();
	}
	else if(whatToCompute=="ibdCovMatrix"){
		char *s = "./tryRES/IBDCovMatrix"; 
		ofstream outfile;
		outfile.open(s);
		std::cout << IBDCovMatrix/numOfSamples;
		std::cout << meansVector /numOfSamples;
		outfile << (IBDCovMatrix/numOfSamples - meansVector*meansVector.transpose()/numOfSamples/numOfSamples);
	}
}
/*! \fn void Model::RSamplerMH(unsigned pedBlockSize, unsigned numLoci,unsigned numOfSamples,const std::string &samplerType="genotypic",const std::string &whatToCompute="genotypeProb")
*  \brief method to sample ordered genotypes using the Metropolis Hastings algorithm to accept samples proposed by sequential imputation.
*/
} ////// end of namespace matvec


