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
#include <iomanip>
#include <cstdio>
#include <vector>
#include <string>

#include "session.h"
#include "util.h"
#include "doublematrix.h"
#include "data.h"

namespace matvec {

Data::Data(void)
{
   numcol     = 0;
   maxnumcol  = 0;
   numrec     = 0;
   new_col    = 1;
   datasheet  = 0;
   hashtable  = 0;
   resize(0,0,20);
   tdfname = SESSION.mktemp();
}

Data::Data(Data& D)
{
   numcol     = 0;
   maxnumcol  = 0;
   numrec     = 0;
   datasheet  = 0;
   hashtable  = 0;
   copyfrom(D);
}

const Data& Data::operator=(Data& A)
{
   copyfrom(A);
   return *this;
}

const Data& Data::operator=(const Field& V)
{
   if (new_col == 1) {  // first column intercept can't be overwritten
      resize(V.len(),numcol,20);
   }
   else {
      if (V.len() != numrec) {
         warning("Data = Col: size incompatible");
         return *this;
      }
   }
   if (!data_in_memory) input_datasheet();
   std::string cname = datasheet[new_col].name();
   int indx = datasheet[new_col].index();
   datasheet[new_col] = V;
   datasheet[new_col].name(cname);
   datasheet[new_col].index(indx);
   if (V.hashtable) hashtable[new_col]->copyfrom( *(V.hashtable));
   save_datasheet(0);     // save the changes on disk, but keep them in memory
   return *this;
}

void Data::copyfrom(Data& A)
{
   if (this == &A) return;
   if (A.data_on_disk == 0) A.save_datasheet();
   resize(A.numrec,A.numcol,A.maxnumcol);
   new_col = A.new_col;
   for (unsigned i=0; i<numcol; i++) hashtable[i]->copyfrom(*(A.hashtable[i]));
   tdfname = A.tdfname;
   data_on_disk   = 1;
   data_in_memory = 0;
   input_datasheet();
   save_datasheet(0);              // save changes
}

Data& Data::resize(const unsigned nr,const unsigned nc,const unsigned mc)
{
   if (numrec == nr && numcol == nc && maxnumcol== mc) return *this;
   release();
   numrec    = nr;
   numcol    = nc;
   if (mc < nc) {
      maxnumcol = nc + 10;          // 10 is the buffer columns
   }
   else {
      maxnumcol = mc + 1;           // first column is reserved for intercept
   }
   if (numrec==0) numcol = 0;
   if (numcol==0) numrec = 0;

   data_in_memory = 1;
   data_on_disk = 0;
   hashtable = new HashTable *[maxnumcol];
   check_ptr(hashtable);
   unsigned i;
   for (i=0; i<maxnumcol; i++) {
      hashtable[i] = new HashTable;
      check_ptr(hashtable[i]);
   }
   datasheet = new Field [maxnumcol];
   check_ptr(datasheet);
   datasheet[0].resize(0);         // first column is reserved for intercept
   for (i=1; i<numcol; i++) datasheet[i].resize(numrec);
   return *this;
}

int Data::field_index(const std::string &colname) const
{
   for (unsigned i=0; i<numcol; i++) {
      if (datasheet[i].name() == colname) {
         return datasheet[i].index();
      }
   }
   return -1;
}

void Data::field_index_vec(Vector<int> &ivec,const std::string &fdname)
{
   if (numcol<1) {   // first column is a reserved: intercept
      return;
   }
   unsigned i,nc;
   if (fdname == "") {
      nc = numcol-1;
      ivec.reserve(nc);
      for (i=0; i<nc; i++) ivec[i] = i+1;   // don't print intercept
   }
   else {
      int k,nskip,j;
      std::string sep(" ,");
      std::string fmt(fdname);
      std::vector<std::string> tmpvec;
      nc = split(fmt,sep,&tmpvec);  //////// fmt.split(n,sep);  field name can be any length
      Vector<int> tmpivec(nc);
      for (nskip=0,i=0; i<nc; i++) {
         k = field_index(tmpvec[i]);
         if (k<0) {
            warning("Data::field_index_vec(): %s: unknown, it's skipped",tmpvec[i].c_str());
            nskip++;
         }
         tmpivec[i] = k;
      }

      ivec.reserve(nc - nskip);
      for (j=0,i=0; i<nc; i++) {
         if (tmpivec[i] >= 0) ivec[j++] = tmpivec[i];
      }
   }
   return;
}

void Data::input(const std::string &fname,const std::string &recfmt)
{
   size_t linewidth = 1024;
   char *line = new char [linewidth];
   int k;
   unsigned i,j,nc,nr,id;
   if (recfmt ==  "") {
      warning("Data::input(): no column-name specified");
      return;
   }
   std::string tmpstr;
   tmpstr = recfmt;
   i = 0;
   while (tmpstr[i] == ' ') {i++;}  // find first nonspace char
   if (tmpstr[i] == '$') throw exception("Data::input(): $ is misplaced");
   i = 0;
   while (tmpstr[i]) {             // move $ to the end of each token
      if (tmpstr[i] == '$' ) {
         tmpstr[i] = ' ';
         j = i;
         while (tmpstr[--j] == ' ');
         tmpstr[++j] = '$';
      }
      i++;
   }
   std::string fmt = "intercept ";      // first column is reserved for intercept
   fmt.append(tmpstr);

   std::string sep(" ,");
   std::vector<std::string> tmpvec;
   unsigned tncol = split(fmt,sep,&tmpvec); ///// split(tncol,sep);  tncol >= 1 is required
   nc = tncol;
   for (i=0; i<tncol; ++i) if (tmpvec[i] == "_skip") nc--;
   std::ifstream  in(fname.c_str(),std::ios::in);
   if (!in) {
     if(line){
      delete [] line;
      line=0;
     }
      throw exception("Data::input(): cannot open file");
   }
   if (!in.getline(line,linewidth)) {
      warning("Data::input(): empty datafile: %s",fname.c_str());
      if(line){
	delete [] line;
	line=0;
      }
      return;
   }
   while (!validline(line)) {
      if (!in.getline(line,linewidth)) {
         warning("Data::input(): no real data in datafile: %s",fname.c_str());
	 if(line){
	   delete [] line;
	   line=0;}
         return;
      }
   }
   std::string T(line);
   i = split(T," ");
   if (i < tncol-1) {
     if(line){
       delete [] line;
       line=0;
     }
      throw exception("Data::input(): the # of columns in data < the expected");
      return;
   }
   in.clear();
   in.seekg(0,std::ios::beg);
   nr = 0;
   while (in.getline(line,linewidth)) if (validline(line)) nr++;
   resize(nr,nc);
   int ThereareStrcol = 0;
   Vector<int> intvec(tncol);
   std::string tstr;
   for (i=0; i<tncol; i++) {
      tstr = tmpvec[i];
      if (tstr.find("_skip") >= 0) {
         for (k=i+1; k<tncol; k++) {
            if (tstr == tmpvec[k]) {
	      if(line){
		delete [] line;
		line=0;
	      }
               throw exception("Data::input(): duplicated column names");
            }
         }
      }
   }
   std::string::size_type begidx;
   for (k=0,i=0; i<tncol; i++) {
      if (tmpvec[i] == "_skip") {
         intvec[i] = -1;
      }
      else {
         intvec[i] = k;
         begidx = tmpvec[i].find("$");
         if (begidx != std::string::npos) {
            tmpvec[i].replace(begidx,1,"");
            datasheet[k].type('S');      // string column
            ThereareStrcol = 1;
            hashtable[k]->resize(numrec);
         }
         datasheet[k].name(tmpvec[i]);
         datasheet[k].index(k);
         k++;
      }
   }      // k == numcol-1
   char *token;
   std::fstream tdatfile(tdfname.c_str(),std::ios::out);

   if (!tdatfile) {
     if(line){
      delete [] line;
      line=0;
     }
      throw exception("Data::input(): cannot open file");
   }
   DataNode* dat_cell;
   double x;
   char *endpt;
   j = 0;

   in.clear();
   in.seekg(0L,std::ios::beg);                // rewind data file
   while (in.getline(line,linewidth)) {
      if (validline(line)) {
         token = strtok(line,", ");
         i = 1;
         while (token) {
            if (i >= tncol) break;
            k = intvec[i++];
            if (k > 0) {
               dat_cell = &datasheet[k][j];
               if (strcmp(token,".")) {
                  dat_cell->missing = 0;
                  if (datasheet[k].type() == 'S') {
                     hashtable[k]->insert(token);
                     id = strlen(token)+1;
                     tdatfile.write((char *)&id,sizeof(unsigned));
                     tdatfile.write(token,id);
                  }
                  else {
                     x = strtod(token,&endpt);   // sscanf(token,"%lf",&x);
                     if (*endpt == '\0') {
                        dat_cell->double_val(x);
                     }
                     else {
                        warning("Data::input(): numeric column has non-numerics "
                              "at the corner of row %d and column %d.\n"
                              "  SUGGESTION: claim it as string column in"
                              " D.input() with $ sign",
                              j+1,i-1);
                        resize(0,0);
                        in.close();
                        tdatfile.close();
			if(line){
                        delete [] line;
			line=0;
			}
                        return;
                     }
                  }
               }
               else {
                  dat_cell->missing = 1;
                  datasheet[k].count_miss(1);
               }
            }
            token = strtok('\0',", ");
         }
         j++;
      }  // end of validline(line)
   }
   in.close();
   tdatfile.close();
   datasheet[0].type('I');    // I = type for intercept
   datasheet[0].nlevel(1);    // I = type for intercept
   datasheet[0].nmiss(0);

   /////////////////////////////////////////////////////
   //   now re-hash for each string field, if necessary
   ////////////////////////////////////////////////////
   if (ThereareStrcol) {
      for (i=1; i<numcol; i++) {
         if (datasheet[i].type() != 'S') continue;
         id = hashtable[i]->size();
         hashtable[i]->resize(id);
         datasheet[i].nlevel(id);
      }
      tdatfile.open(tdfname.c_str(),std::ios::in);
      for (i=0; i<numrec; i++) {
         for (j=1; j<numcol; j++) {
            dat_cell = &datasheet[j][i];
            if (datasheet[j].type() == 'S' && !(dat_cell->missing)) {
               tdatfile.read((char *)&id,sizeof(unsigned));
               tdatfile.read(line,id);
               id = hashtable[j]->insert(line);
               dat_cell->unsigned_val(id);
            }
         }
      }
      tdatfile.close();
   }
   if(line){
     delete [] line;
     line=0;
   }
   ////////////////////////////////////////////////////////////////////
   // save a copy of data is a must. Because data could be changed
   // temporarily for some special purposes, the change can be droped
   // by release datasheet
   ////////////////////////////////////////////////////////////////////
   save_datasheet(0);    // save a copy to hard-disk
}

/*
void Data::drop(const char *fdname)
{
   int nc=0;
   int *intvec = field_index_vec(nc,fdname);
   for (int i=0; i<nc; i++) {

   }
   if (intvec) delete [] intvec;

}

void Data::keep(const char *fdname)
{
   int nc=0;
   int *intvec = field_index_vec(nc,fdname);
   for (int i=0; i<nc; i++) {

   }
   if (intvec) delete [] intvec;
}
*/

void Data::value_for_missing(const double vm)
{
   if (!data_in_memory) input_datasheet();
   for (unsigned i=0; i<numcol; i++) datasheet[i].value_for_missing(vm);
}

void Data::save_datasheet(const int relse)
{
   if (!datasheet) {
      warning("Data::save_datasheet(): no data to save");
      return;
   }
   std::ofstream df(tdfname.c_str(),std::ios::out);
   if (!df) throw exception("Data::save_datasheet(): cannot open file");
   for (unsigned i=1; i<numcol; i++) {   // first column is an intercept
      df.write((char *)datasheet[i].dat_vec,numrec*sizeof(DataNode));
   }
   df.close();
   data_on_disk = 1;
   if (relse) release_datasheet();
}

void  Data::input_datasheet(void)
{
   if (data_in_memory) return;
   if (data_on_disk) {
      std::ifstream df(tdfname.c_str());
      if (!df) throw exception("Data::input_datasheet(): cannot open file");
      for (unsigned i=1; i<numcol; i++) {  // first column is an intercept
         datasheet[i].resize(numrec);
         df.read((char *)datasheet[i].dat_vec,numrec*sizeof(DataNode));
      }
      df.close();
      data_in_memory = 1;                // data now is in memorry
   }
   else {
      warning("Data::input_datasheet(): data is not on disk");
   }
}

void Data::release_datasheet(void)
{
   if (datasheet) {
      // any data must have a hard copy in disk
      if (!data_on_disk) save_datasheet();
      for (unsigned i=1; i<numcol; i++) datasheet[i].resize(0);
      data_in_memory = 0;     // data is not in memory, but should  on disk;
   }
}

void Data::row(const unsigned i,DataNode* recd)
{
   if (!data_in_memory) input_datasheet();
   if (!recd) {
     if(numcol>0) {
       recd = new DataNode [numcol];
     }
     else {
       recd = 0;
     }
   }
   for (unsigned j=1; j<numcol; j++) recd[j] = datasheet[j][i];
}

Field Data::col(const std::string &cname)
{
   int k =  field_index(cname);
   if (k<=0) {        // first column intercept is not accessible
      warning("Data::col(%s): no such column",cname.c_str());
      return Field();
   }
   if (!data_in_memory) input_datasheet();

   HashTable *tmp_hashtable = 0;
   DataNode *retval;
   if (numrec>0){
     retval = new DataNode [numrec];
   }
   else {
     retval = 0;
   }
   unsigned i;
   DataNode *colk = datasheet[k].dat_vec;
   if (datasheet[k].type()=='S') {
      tmp_hashtable = new HashTable;
      *tmp_hashtable = *(hashtable[k]);
      for (i=0; i<numrec; i++) {
         if (colk[i].missing) {retval[i].missing = 1;}
         else { retval[i].unsigned_val(colk[i].unsigned_val()); }
      }
   }
   else {
      for (i=0; i<numrec; i++) {
         if (colk[i].missing) { retval[i].missing = 1; }
         else { retval[i].double_val(colk[i].double_val()); }
      }
   }
   return Field(numrec,retval,datasheet[k].col_struct,tmp_hashtable);
}

DataNode* Data::rawcol(const std::string &cname)
{
   int k =  field_index(cname);
   if (k > 0) {              // first column intercept is not accessible
      if (!data_in_memory) input_datasheet();
      return datasheet[k].dat_vec;
   }
   else {
      warning("Data::rawcol(%s): no such column",cname.c_str());
      return 0;
   }
}

DataNode* Data::rawcol(unsigned c)
{
   if (c <= 0 || c >= numcol) throw exception("Data::rawcol(): out of range");
   if (!data_in_memory) input_datasheet();
   return datasheet[c].dat_vec;
}

Data& Data::newcol(const std::string &cname)
{
   if (cname == "") {
      warning("Data::newcol(cname), cname is empty");
      return *this;
   }
   unsigned i;
   int k =  field_index(cname);
   if (k > 0) {         // first column intercept cannot be overwritten
      if (datasheet[k].type() == 'S') {
         warning("Data::newcol(): %s exits, can't overwrite string column",cname.c_str());
         return *this;
      }
      warning("Data.newcol(): %s exits, it's been overwritten",cname.c_str());
      new_col = k;
   }
   else {
      if (!data_in_memory) input_datasheet();   // data must be in memory
      if (numcol == maxnumcol) {               // data sheet is full
         Field *tmp_datasheet = new Field [maxnumcol];
         check_ptr(tmp_datasheet);
         HashTable **tmp_hashtable = new HashTable *[maxnumcol];
         check_ptr(tmp_hashtable);
         for (i=0; i<maxnumcol; i++) {
            tmp_datasheet[i] = datasheet[i];
            tmp_hashtable[i] = hashtable[i];
         }
	 if(datasheet){
	   delete [] datasheet;            // note I do not delete datasheet[i]
	 datasheet=0;
	 }
	 if(hashtable){
	   delete [] hashtable;            // note I do not delete hashtable[i]
	   hashtable=0;
	 }
         maxnumcol += 10;
         datasheet = new Field [maxnumcol];
         check_ptr(datasheet);
         hashtable  = new HashTable *[maxnumcol];
         check_ptr(hashtable);
         for (i=0; i<numcol; i++) {
            datasheet[i] = tmp_datasheet[i];
            hashtable[i] = tmp_hashtable[i];
         }
         for (i=numcol; i<maxnumcol; i++)  {
            hashtable[i] = new HashTable;
            check_ptr(hashtable[i]);
            datasheet[i] = 0;
         }
	 if(tmp_datasheet){
	   delete [] tmp_datasheet;
	   tmp_datasheet=0;
	 }
	 if(tmp_hashtable){
	   delete [] tmp_hashtable;
	   tmp_hashtable=0;
	 }
      }
      new_col = numcol++;
      datasheet[new_col].name(cname);
      datasheet[new_col].type('F');   // floating point number for the colum
      datasheet[new_col].index(new_col);
      datasheet[new_col].resize(numrec);
   }
   return  *this;
}


void Data::newcol(const std::string &cname,const Field& col)
{
   unsigned i,n = col.size();
   if (n != numrec)
      warning("Data::newcol():%d,%d: size not conformable",numrec,n);
   this->newcol(cname);
   datasheet[new_col].col_struct = col.col_struct;
   datasheet[new_col].name(cname);   // cname override Field.name()
   datasheet[new_col].index(new_col);
   if (numrec < n) n = numrec;
   DataNode *tc = datasheet[new_col].dat_vec;
   for (i=0; i<n; i++) tc[i] = col.elem(i);
   for (i=n; i<numrec; i++) tc[i].missing = 1;
   datasheet[new_col].count_miss(numrec-n);
}

Data& Data::adjoin(Data& b)
{
   unsigned n = b.num_rows();
   if (n != numrec) {
      warning("Data::adjoin(b):%d,%d: size unconformable: truncated",numrec,n);
   }
   if (!b.in_memory()) b.input_datasheet();
   unsigned i,nc=b.num_cols();
   for (i=0; i<nc; i++) {
      this->newcol("junk");
      datasheet[new_col] = b.datasheet[i];
      datasheet[new_col].index(new_col);
   }
   save_datasheet(0);
   return *this;
}

Data& Data::stack(Data& b)
{
   warning("Data::stack(b): not yet available");
   return *this;
}

DataNode* Data::cell(const unsigned r,const unsigned c)
{
   if (!data_in_memory) input_datasheet();   // data must be in memory
   if (r>=numrec || (c>=numcol&& c==0)) {
      warning("Data::cell(%d,%d): out of range",c,r);
      return 0;
   }
   else {
      return &datasheet[c][r];
   }
}

void Data::release(void)
{
   if (datasheet) {
     delete [] datasheet;
     datasheet = 0;
   }
   if (hashtable) {
     for (int i=maxnumcol-1; i>=0; i--){
       if(hashtable[i]){
	 delete hashtable[i];
	 hashtable[i]=0;
       }
     }
     if(hashtable){
       delete [] hashtable; 
       hashtable=0;
     }
   }
}


Field Data::max(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].max();
   return xcol;
}

Field Data::min(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].min();
   return xcol;
}

Field Data::mean(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].mean();
   return xcol;
}

Field Data::variance(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].covariance();
   return xcol;
}

Field Data::sum(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].sum();
   return xcol;
}

Field Data::sumsq(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].sumsq();
   return xcol;
}

Field Data::product(const std::string &cname)
{
   if (!data_in_memory) input_datasheet();
   int i,nc=0;
   Vector<int> ivec;
   field_index_vec(ivec,cname);
   nc = ivec.size();
   Field xcol(nc);
   for (i=0; i<nc; i++) xcol[i] = datasheet[ivec[i]].product();
   return xcol;
}

void Data::stat(void)
{
   unsigned W = SESSION.output_precision+6;
   if (!data_in_memory) input_datasheet();
   unsigned i;

   std::cout << "\n  Name";
   for (i=1; i<numcol; i++) {   // first column intercept should be ignored
      if (datasheet[i].type()=='S') continue;
      std::cout << " " << std::setw(W) << datasheet[i].name();
   }
   std::cout << "\n";

   std::cout << "  Nobs";
   for (i=1; i<numcol; i++) {
      if (datasheet[i].type()=='S') continue;
      std::cout << " ";
      if (datasheet[i].type() == 'F') {
         std::cout << std::setw(W) <<  numrec-datasheet[i].nmiss();
      }
      else {
         std::cout << std::setw(W) << ".";
      }
   }
   std::cout << "\n";

   std::cout << "  Min ";
   for (i=1; i<numcol; i++) {
      if (datasheet[i].type()=='S') continue;
      std::cout << datasheet[i].min();
   }
   std::cout << "\n";

   std::cout << "  Max ";
   for (i=1; i<numcol; i++) {
      if (datasheet[i].type()=='S') continue;
      std::cout << datasheet[i].max();
   }
   std::cout << "\n";

   std::cout << "  Mean";
   for (i=1; i<numcol; i++) {
      if (datasheet[i].type()=='S') continue;
      std::cout << datasheet[i].mean();
   }
   std::cout << "\n";

   std::cout << "  S.D.";
   DataNode var;
   for (i=1; i<numcol; i++) {
      if (datasheet[i].type()=='S') continue;
         var = datasheet[i].covariance();
      if (!var.missing) var.double_val(std::sqrt(var.double_val()));
      std::cout << var;
   }
   std::cout << "\n\n";

   return;
}

doubleMatrix Data::mat(void)
{
   doubleMatrix retval(numrec,numcol);
   if (numrec==0) throw exception("Data::mat(): empty data object");
   if (!data_in_memory) input_datasheet();
   DataNode *dat = 0;
   unsigned i,j;
   double *dpt;
   for (i=0; i<numrec; i++) {
      dpt = retval[i];
      for (j=1; j<numcol; j++) {   // first column intercept should be ignored
         if (datasheet[j].type() != 'S') {
            dat = &(datasheet[j].dat_vec[i]);
            if (dat->missing == 0) {
               dpt[j] = dat->double_val();
            }
            else {
               dpt[j] = 0.0;
            }
         }
      }
   }
   release_datasheet();
   return retval;
}

void Data::print(std::ostream& stream,const Vector<int> intvec,const int ic)
{
   if (numrec==0) {
      std::cout << "\t empty data object\n" << std::flush;
      return;
   }
   if (!data_in_memory) input_datasheet();
   int nc = intvec.size();
   int kk;
   unsigned i,j,k,id;
   unsigned W = SESSION.output_precision+6;
   const char *str;
   char ch;
   stream.precision(SESSION.output_precision);
   DataNode *dat = 0;
   for (k=23,i=0; i<numrec; i++) {
      if (ic && i>=k) {
         k += 23;
         stream << "  more ... [q for quit] ";
         std::cin.get(ch);
         std::cin.seekg(0L,std::ios::beg);
         if (ch == 'q') break;
      }
      for (j=0; j<nc; j++) {
         kk = intvec[j];
         if (kk < 1) continue;   // first coloumn is reserved for intercept
         dat = &(datasheet[kk].dat_vec[i]);
         if (dat->missing) {
            stream << " " << std::setw(W) << ".";
         }
         else {
            if (datasheet[kk].type()=='S') {
               id = dat->unsigned_val();
               str = (const char*)(hashtable[kk]->find(id));
               stream << " " << std::setw(W) << str;
            }
            else {
               stream << " " << std::setw(W) << dat->double_val();
            }
         }
      }
      stream << "\n";
   }
   stream << std::flush;
   release_datasheet();
}

void Data::display(const std::string &fdname,const int ic)
{
   Vector<int> intvec;
   field_index_vec(intvec,fdname);
   print(std::cout,intvec,ic);
}

void Data::save(const std::string &fname,
                   const int io_mode )
{
   std::ofstream ofs;
   ofs.open(fname.c_str(),(OpenModeType)io_mode);
   if (!ofs) throw exception("Data::save(): cannot open file");
    Vector<int> intvec;
    field_index_vec(intvec);
    print(ofs,intvec,0);
    ofs.close();
}

std::ostream& operator<<(std::ostream& stream, Data& A)
{
   unsigned nc = A.num_cols() - 1;      // first column is reserved for intercept
   if (nc == 0) return stream;
   Vector<int> intvec(nc);
   for (int i=0; i<nc; i++) intvec[i] = i+1;
   A.print(stream,intvec,1);
   return stream;
}

} ////////// end of namespace matvec

