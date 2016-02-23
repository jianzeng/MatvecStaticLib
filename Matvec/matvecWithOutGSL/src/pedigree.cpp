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

#include <fstream>
#include <set>
#include <string>
#include <iomanip>

#include "session.h"
#include "util.h"
#include "stat.h"
#include "hashtable.h"
#include "pedigree.h"

static int INPUT_LINE_WIDTH = 1024;
namespace matvec {
	Pedigree::ped_type Pedigree::type; 

void getlambda(double **lambda, const int n)
{
   if (lambda) {
      if (n == 3 ) {
         lambda[0][0] =  1.0; lambda[0][1] = -0.5; lambda[0][2] = -0.5;
         lambda[1][0] = -0.5; lambda[1][1] =  0.25; lambda[1][2] = 0.25;
         lambda[2][0] = -0.5; lambda[2][1] =  0.25; lambda[2][2] = 0.25;
      }
      else {
         throw exception("getlambda(a,n),n must be 3");
      }
   }
   else {
      throw exception("getlambda(l,n): l is NULL");
   }
}

int pedcp1(const void* a, const void* b)
{
   PedNode **x,**y;
   x = (PedNode **)a; y = (PedNode **)b;
   return (strcmp((*x)->myname,(*y)->myname));
}

int pedcp2(const void* a, const void* b)
{
   PedNode **x,**y;
   x = (PedNode **)a; y = (PedNode **)b;
   return ((*x)->myid - (*y)->myid);
}

Pedigree::Pedigree(char pedtype[])
{
   /////////////////////////////////////////////////
   // type  = "group"          ->    group
   //       = "standard"
   //       = anything else    ->    nongrouop
   //
   /////////////////////////////////////////////////
   nang          = 0;
   fdone         = 0;
   totalna       = 0;
   numgroup      = 0;
   basesize      = 0;
   maxnamelen    = 0;
   pedmember     = 0;
   ped_storage   = 0;
   sex_confirmed = 0;
   groupingped   = 0;
   type          = Pedigree::raw;
   if (strcmp(pedtype,"group") == 0) groupingped = 1;
   if (strcmp(pedtype,"standard") == 0) type = Pedigree::standard;

   ped_on_disk  = 0;
   ped_in_memory = 1;
   binfname = SESSION.mktemp();
}

Pedigree::Pedigree(Pedigree& P)
{
   totalna     = 0;
   pedmember   = 0;
   ped_storage = 0;
   groupingped = 0;
   copyfrom(P);
}

void Pedigree::copyfrom(Pedigree& P)
{
   if (this == &P) return;
   if (P.ped_in_memory) P.release_pedsheet();
   resize(P.nang);
   fvec          = P.fvec;
   fdone         = P.fdone;
   numgroup      = P.numgroup;
   basesize      = P.basesize;
   totalna       = P.totalna;       // totalna is necessary for save_pedsheet()
   maxnamelen    = P.maxnamelen;
   groupingped   = P.groupingped;
   type          = P.type;
   sex_confirmed = P.sex_confirmed;

   if (P.ped_on_disk) {
      binfname = P.binfname;
      if (P.ped_in_memory) P.release_pedsheet();
      ped_in_memory = 0;
      ped_on_disk = 1;
      input_pedsheet();
   }
   else {
      if (totalna != 0) throw exception("Pedigree::copyfrom(): you have probably found a bug!");
      ped_on_disk = 0;
   }
   binfname = SESSION.mktemp(); // each ped object must have its own binfname
   if (totalna != 0) save_pedsheet(0);  // save changes without memory-release
   ped_in_memory = 1;
}

void Pedigree::resize(const unsigned na)
{
   if (nang == na) return;
   release();
   nang = na;
   if(!(ped_storage = new PedNode[nang])) throw exception("Pedigree.resize(): out of memory");
   if (!(pedmember = new PedNode*[nang])) throw exception("Pedigree.resize(): out of memory");
   PedNode* I;
   for (unsigned i=0; i<nang; i++) {
      I = &(ped_storage[i]);
      I->mysex = '.';
      I->myname = 0;
      I->namelen = 0;
      pedmember[i] = I;
   }
}

unsigned Pedigree::input_standard(std::istream& pedfile,unsigned individual,
                                unsigned mother,unsigned father)
{
   size_t linewidth = INPUT_LINE_WIDTH;
   char *line;
   if(linewidth>0){
     line = new char [linewidth];
   }
   else {
     line = 0;
   }
   size_t isdwidth = sizeof(unsigned)*4;

   pedfile.clear();
   pedfile.seekg(0L,std::ios::beg);
   unsigned i,numrec = 0;
   while (pedfile.getline(line,linewidth)) if (validline(line)) numrec++;
   if (numrec == 0) {
     if(line){
       delete [] line;
       line=0;
     }
      nang  = totalna  = numrec;
      throw exception("Pedigree::input_standard(): pedigree is empty");
   }
   pedfile.clear();
   pedfile.seekg(0L,std::ios::beg);

   resize(numrec);
   PedNode *I;
   char *temp;
   numrec = 0;
   basesize = 0;
   unsigned numsameid = 0;
   while (pedfile.getline(line,linewidth)) {
      if (!validline(line)) continue;
      I = pedmember[numrec++];
      temp = strtok(line,", ");
      i = 0;
      while (temp) {
         if (i==individual) {
            if (!sscanf(temp,"%d",&(I->myid))) {
               numrec--;
               throw exception("Pedigree::input_standard(): non-numerical character"
                     " in your standard pedigree");
            }
         }
         else if (i==father) {
            if (strcmp(temp,".") == 0) {
               I->father_id = 0;
            }
            else {
               if (!sscanf(temp,"%d",&(I->father_id))) {
                  numrec--;
                  throw exception("Pedigree::input_standard(): non-numerical character"
                        " in your standard pedigree");
               }
            }
         }
         else if (i==mother) {
            if (strcmp(temp,".") == 0) {
               I->mother_id = 0;
            }
            else {
               if (!sscanf(temp,"%d",&(I->mother_id))) {
                  numrec--;
                  throw exception("Pedigree::input_standard(): non-numerical character"
                        " in your standard pedigree");
               }
            }
         }
         i++;
         temp = strtok('\0',", ");
      }
      if (I->father_id==0 && I->mother_id==0) basesize++;
   }
   nang  = totalna  = numrec;
   if(line){
     delete [] line;
     line=0;
   }
   fvec.resize(totalna,0.0);
   fdone = 0;
   save_pedsheet(0);
   return totalna;
}

unsigned Pedigree::input(const std::string &pfname,const std::string &recfmt,
                         const std::string &pedtype)
{
   std::string fmt(recfmt);
   std::string::size_type begidx,endidx;
   // replace any $ signs with blanks
   begidx = fmt.find("$");
   while (begidx != std::string::npos) {
      fmt[begidx] = ' ';
      begidx = fmt.find("$",begidx+1);
   }

   std::string tempstr,sep(" ,");
   unsigned i,j,numrec,k;
   PedNode *I;
   int ncol,ind=-1,mother=-1,father=-1,sex=-1;
   ncol = 0;
   begidx = fmt.find_first_not_of(sep);
   while(begidx != std::string::npos) {
      endidx = fmt.find_first_of(sep,begidx);
      if (endidx == std::string::npos) endidx = fmt.length();
      tempstr = fmt.substr(begidx,endidx - begidx);
      if (tempstr == "individual") {
         ind = ncol;
      } else if (tempstr == "father") {
         father = ncol;
      } else if (tempstr == "mother") {
         mother = ncol;
      } else if (tempstr == "sex") {
         sex = ncol;
      }
      ncol++;
      begidx = fmt.find_first_not_of(sep,endidx);
   }
   if (ncol < 3) throw exception("Pedigree::input(fmt): fmt must contain at least three columns");
   if (ind < 0) throw exception("Pedigree::input(fmt): no individual column in fmt");
   if (mother < 0) throw exception("Pedigree::input(fmt): no mother column in fmt");
   if (father < 0) throw exception("Pedigree::input(fmt): no father column in fmt");
   if (sex >= 0) sex_confirmed = 1;

//    std::string tmpfname;
//    tmpfname = SESSION.mktemp();

   std::ifstream pedfile(pfname.c_str(),std::ios::in);
   if (!pedfile) throw exception("Pedigree::input(): can't open pedigree file");
   if (pedtype != "") {
      if (pedtype == "raw") {
         type = Pedigree::raw;
      }
      else if (pedtype == "standard") {
         type = Pedigree::standard;
         numrec = input_standard(pedfile,ind,mother,father);
         pedfile.close();
         return numrec;
      }
      else if (pedtype == "group") {
         groupingped = 1;
      }
      else {
         warning("Pedigree::input(): unknown pedigree type: arg3 is ignored");
      }
   }

   size_t linewidth = INPUT_LINE_WIDTH;
   char *line;
   if(linewidth>0){
     line = new char [linewidth];
   }
   else {
     line = 0;
   }

   numrec = 0;
   while (pedfile.getline(line,linewidth)) if (validline(line)) numrec++;
   if (numrec == 0) {
     if(line){
       delete [] line;
       line=0;
     }
      pedfile.close();
      nang  = totalna  = numrec;
      throw exception("Pedigree::input(): pedigree is empty");
   }
   pedfile.clear();
   pedfile.seekg(0L,std::ios::beg);

   size_t isdwidth = sizeof(unsigned)*3;
   char *wd;
   if(linewidth>0){
     wd = new char [linewidth];
   }
   else {
     wd = 0;
   }
   char *sireid;
   if(linewidth>0){
     sireid= new char [linewidth];
   }
   else {
     sireid = 0;
   }
   char *damid;
   if(linewidth>0){
     damid = new char [linewidth];
   }
   else {
     damid = 0;
   }
   char sexid;
   sexid = '.';
   Vector<unsigned> stdrec(3);

   std::string tmpfname;
   tmpfname = SESSION.mktemp();
   
   std::fstream tpedfile(tmpfname.c_str(),std::ios::out);
   if (!tpedfile) throw exception("can't open internal pedigree file " +  tmpfname);

   resize(numrec);
   HashTable pedigreeTable(numrec,"str");
   char *temp;
   maxnamelen = 0;
   numrec = 0;
   unsigned lineNum = 0;
   while (pedfile.getline(line,linewidth)) {
      lineNum++;
      if (!validline(line)) continue;
      temp = strtok(line,", ");
      i = 0;
      while (temp) {
         if (i == ind) {
            strcpy(wd, temp);
         }
         else if (i == mother) {
            strcpy(damid,temp);
         }
         else if (i == father) {
            strcpy(sireid,temp);
         }
         else if (i == sex) {
            sexid = temp[0];
         }
         i++;
         temp = strtok('\0',", ");
      }
      if (i < ncol) {
         std::cerr << "pedigree " << pfname << ", line " << lineNum << ": invalid\n";
         throw exception("invalid line");
      }
      if (!groupingped && strcmp(damid,sireid) == 0) {
         if (strcmp(sireid,".") && strcmp(sireid,"0")) {
            warning("pedigree %s, line %d, identical parents: %s. It's ignored",
                     pfname.c_str(),lineNum,sireid);
            continue;
         }
      }
      if (pedigreeTable.get_id(wd)) {
         warning("pedigree %s, line %d, duplicated id: %s. It's ignored",
                 pfname.c_str(),lineNum,wd);
         continue;
      }
      else {
         pedigreeTable.insert(wd);
      }
      stdrec[0] = strlen(wd)+1;
      stdrec[1] = strlen(damid)+1;
      stdrec[2] = strlen(sireid)+1;
      I = pedmember[numrec];
      I->myname = new char [stdrec[0]];
      strcpy(I->myname,wd);
      I->namelen = stdrec[0];
      if (stdrec[0] > maxnamelen) maxnamelen = stdrec[0];
      numrec++;
      tpedfile.write((char*)stdrec.begin(),isdwidth);
      tpedfile.write(wd,(size_t)stdrec[0]);
      tpedfile.write(damid,(size_t)stdrec[1]);
      tpedfile.write(sireid,(size_t)stdrec[2]);
      tpedfile.write((char*)&sexid,1);
   }
   if (maxnamelen > 0) maxnamelen--;
   pedigreeTable.release();
   if(line){
     delete [] line;
     line=0;
   }
   pedfile.close();
   tpedfile.close();
   qsort((char *)pedmember,numrec,sizeof(PedNode*),pedcp1);
/*
   for (i=1; i<numrec; i++) {
      if (strcmp(pedmember[i]->myname,pedmember[i-1]->myname) == 0)  {
         cout << " in pedigree " << pfname <<"\n";
         cout << " duplicated id: " << pedmember[i]->myname << "\n";
         delete [] wd;
         delete [] damid;
         delete [] sireid;
         delete [] stdrec;
         release();
         return 0;
      }
   }
*/
   tpedfile.open(tmpfname.c_str(),std::ios::in);

   std::set<std::string> table;             // binary tree
   PedNode **d,**s;
   d = new PedNode*[1];
   s = new PedNode*[1];
   d[0] = new PedNode[1];
   s[0] = new PedNode[1];
   d[0]->myname = damid;
   s[0]->myname = sireid;

   char *c;
   for (i=0; i<numrec; i++) {
      tpedfile.read((char*)stdrec.begin(),isdwidth);
      tpedfile.read(wd,(size_t)stdrec[0]);
      tpedfile.read(damid,(size_t)stdrec[1]);
      tpedfile.read(sireid,(size_t)stdrec[2]);
      tpedfile.read((char*)&sexid,1);
      if (groupingped) {
         c =(char *)bsearch((char *)d,(char *)pedmember,numrec,sizeof(PedNode*),
                             pedcp1);
         if (!c) table.insert(damid);
         c =(char *)bsearch((char *)s,(char *)pedmember,numrec,sizeof(PedNode*),
                           pedcp1);
         if (!c) table.insert(sireid);
      }
      else {
         if (strcmp(damid,"0") && strcmp(damid,".")) {
            c = (char *)bsearch((char *)d,(char *)pedmember,numrec,
                        sizeof(PedNode*),pedcp1);
            if (!c) table.insert(damid);
         }
         if (strcmp(sireid,"0") && strcmp(sireid,".")) {
            c =(char *)bsearch((char *)s,(char *)pedmember,numrec,
                                sizeof(PedNode*),pedcp1);
            if (!c) table.insert(sireid);
         }
      }
   }

   basesize = table.size();
   ///////////////////////////////////////////////////////////////////
   // it is not really # of base animals, it is # of animals which
   // are not listed in the individual column. However, it is # of
   // unknown parents groups for genetic group model.
   // Real basesize will be calculated in renum().
   ///////////////////////////////////////////////////////////////////

   if (groupingped) {
      totalna  = numrec;
      numgroup = basesize;
   }
   else {
      totalna  = basesize + numrec;
      numgroup = 0;
   }
   fvec.resize(totalna,0.0);
   fdone = 0;

   if (basesize ) {      // now hash sire and dam which are not in the id list
      std::fstream bpedfile(binfname.c_str(),std::ios::out);
      if (!bpedfile) {
         release();
         throw exception("Pedigree.input(): can't open internal file");
      }
      for (i=0; i<numrec; i++) {
         I = pedmember[i];
         j = I->namelen;
         bpedfile.write((char*)&j,sizeof(unsigned));
         bpedfile.write(I->myname,j);
      }
      bpedfile.close();
      bpedfile.open(binfname.c_str(),std::ios::in);

      for (i=0; i<numrec; i++) {
	if(ped_storage[i].myname){
	  delete [] ped_storage[i].myname;
	  ped_storage[i].myname=0;
	}
      }
      if(ped_storage){
	delete [] ped_storage;
	ped_storage=0;
      }
      if(pedmember){
	delete [] pedmember;
	pedmember=0;
      }
      nang = basesize + numrec;
      if(nang>0){
	ped_storage = new PedNode[nang];
      }
      else {
	ped_storage = 0;
      }
      check_ptr(ped_storage);
      pedmember = new PedNode*[nang];
      check_ptr(pedmember);
      for (i=0; i<nang; i++) {
         I = &(ped_storage[i]);
         I->mysex = '.';
         I->myname = 0;
         I->namelen = 0;
         pedmember[i] = I;
      }

      for (i=0; i<numrec; i++) {
         I = pedmember[i];
         bpedfile.read((char*)&j,sizeof(unsigned));
         I->namelen = j;
	 if(j>0){
	   I->myname = new char [j];
	 }
	 else {
	   I->myname = 0;
	 }
         check_ptr(I->myname);
         bpedfile.read(I->myname,(size_t)j);
      }
      bpedfile.close();

      qsort((char *)pedmember,numrec,sizeof(PedNode*),pedcp1);
      std::set<std::string>::iterator pos;
      for (i=0,pos=table.begin(); pos != table.end(); ++pos,++i) {
         j = numrec+i;
         I = pedmember[j];
         k = I->namelen = pos->length()+1;
	 if(k>0){
	   I->myname = new char [k];
	 }
	 else {
	   I->myname = 0;
	 }
         strcpy(I->myname, pos->c_str());
      }
   }

   for (i=0; i<nang; i++) {
      I = pedmember[i];
      I->myid = i+1; I->father_id = 0; I->mother_id = 0;
   }
   tpedfile.close();
   tpedfile.open(tmpfname.c_str(),std::ios::in);

   PedNode **w,**pm;
   w = new PedNode*[1];
   w[0] = new PedNode[1];
   w[0]->myname = wd;
   for (i=0; i<numrec; i++) {
      tpedfile.read((char*)stdrec.begin(),isdwidth);
      tpedfile.read(wd,(size_t)stdrec[0]);
      tpedfile.read(damid,(size_t)stdrec[1]);
      tpedfile.read(sireid,(size_t)stdrec[2]);
      tpedfile.read((char*)&sexid,1);
      pm = (PedNode **)bsearch((char *)w,(char *)pedmember,numrec,
                              sizeof(PedNode*),pedcp1);
      I = pedmember[pm[0]->myid - 1];
      if (groupingped) {
         pm = (PedNode **)bsearch((char *)d,(char *)pedmember,numrec,
                                 sizeof(PedNode*),pedcp1);
         if (!pm) {
            pm = (PedNode **)bsearch((char *)d,(char *)&(pedmember[numrec]),
                                 basesize, sizeof(PedNode*),pedcp1);
            if (!pm) throw exception("Pedigree::input(): you have probably found a bug!");
         }
         I->mother_id = pm[0]->myid;

         pm = (PedNode **)bsearch((char *)s,(char *)pedmember,numrec,
                                  sizeof(PedNode*),pedcp1);
         if (!pm) {
            pm = (PedNode **)bsearch((char *)s,(char *)&(pedmember[numrec]),
                                 basesize, sizeof(PedNode*),pedcp1);
           if (!pm) throw exception("Pedigree::input(): you have probably found a bug!");
         }
         I->father_id = pm[0]->myid;
      }
      else {
         if (strcmp(damid,"0") && strcmp(damid,".")) {
            pm = (PedNode **)bsearch((char *)d,(char *)pedmember,numrec,
                                    sizeof(PedNode*),pedcp1);
            if (!pm) {
               pm =(PedNode **)bsearch((char *)d,(char *)&(pedmember[numrec]),
                                    basesize, sizeof(PedNode*),pedcp1);
               if (!pm) throw exception("Pedigree::input(): you have probably found a bug!");
            }
            I->mother_id = pm[0]->myid;
         }
         if (strcmp(sireid,"0") && strcmp(sireid,".") ) {
            pm = (PedNode **)bsearch((char *)s,(char *)pedmember,numrec,
                                     sizeof(PedNode*),pedcp1);
            if (!pm) {
               pm =(PedNode **)bsearch((char *)s,(char *)&(pedmember[numrec]),
                                    basesize, sizeof(PedNode*),pedcp1);
               if (!pm) throw exception("Pedigree::input(): you have probably found a bug!");
            }
            I->father_id = pm[0]->myid;
         }
      }
      I->mysex = sexid;
   }
   tpedfile.close();
   if(wd){
     delete [] wd;
     wd=0;
   }
   if(damid){
     delete [] damid;
     damid=0;
   }
   if(sireid){
     delete [] sireid;
     sireid=0;
   }
   if(w[0]){
     delete [] w[0];  
     w[0]=0;
   }
   if(w){
     delete [] w;
     w=0;
   }
   if(d[0]){
     delete [] d[0];  
     d[0]=0;
   }
   if(d){
     delete [] d;
     d=0;
   }
   if(s[0]){
     delete [] s[0];  
     s[0]=0;
   }
   if(s){
     delete [] s;
     s=0;
   }
   renum();
   save_pedsheet();
   remove(tmpfname.c_str());
   return numrec;
}

Pedigree& Pedigree::sample(const unsigned maxsize,const unsigned maxg,
                          const unsigned dam0,const unsigned sire0,
                          const double imrate,const unsigned parent,
                          const int no_po, const int no_fsib,
                          const double sexratio)
{
   if (maxg < 1) throw exception("Pedigree::sample(maxsize,maxg,...): maxg must >= 1");
   size_t isdsexwidth = 4*sizeof(unsigned);
   unsigned tgsire, tsire = sire0;
   unsigned tgdam, tdam = dam0;
   if ((tdam+tsire) > maxsize) throw exception("Pedigree::sample(): out of memory");
   unsigned s,d,i,j,max_male,max_female;
   d = maxsize-tsire-tdam;
   s = static_cast<unsigned>(d*sexratio);
   max_female = tdam+d-s+200;
   max_male = tsire+s+200;           // 200 is the buffer for randomness
   Vector<unsigned> female(max_female--);
   Vector<unsigned> male(max_male--);

   for (i=0; i<=tdam; i++) female[i]=i;
   for (male[0]=0,i=1; i<=tsire; i++) male[i] = tdam + i;
   totalna = tdam + tsire;
   unsigned npg = (maxsize-tdam-tsire)/maxg +1;
   Matrix<unsigned> pedholder(3,maxsize+1);

  // sexcode:  female = 1, male = 2
   for (i=1; i<=tdam; i++) pedholder[2][i] = 1;
   for (i=tdam+1; i<=totalna; i++) pedholder[2][i] = 2;

   long selecteddam,selectedsire;
   unsigned ds,dd,ss,sd;
   int tryagain,ntry;
   for (tgdam=tdam, tgsire=tsire, i=0; i<maxg; i++) {
      if (totalna >= maxsize) break;
      for (j=0; j<npg; j++) {
         if (totalna >= maxsize) break;
         if (ran1() > imrate) {
            ntry = 0;
            do {
               selecteddam = ignbin(long(tdam),0.6);
               selectedsire = ignbin(long(tsire),0.6);
               if (selecteddam == 0) selecteddam = parent;
               if (selectedsire == 0) selectedsire = parent;
               tryagain = 0;
               d = female[selecteddam];
               s = male[selectedsire];
               dd = pedholder[0][d];
               ds = pedholder[1][d];
               sd = pedholder[0][s];
               ss = pedholder[1][s];
               if (no_po) {
                  if ( d>0 && sd == d) {
                     tryagain = 1;
                  }
                  else if (s > 0 && ds == s) {
                     tryagain = 1;
                  }
               }
               if (no_fsib) {
                  if (ss > 0 && sd >0) if (ss == ds && sd == dd) tryagain = 1;
               }
               if (tryagain) ntry++;
               if (ntry >100) {
                  tryagain = 0;
                 selecteddam = 0;
                 selectedsire = 0;
               }
            } while (tryagain);
         }
         else {
            selecteddam = 0;
            selectedsire = 0;
         }
         totalna++;
         pedholder[0][totalna] = female[selecteddam];
         pedholder[1][totalna] = male[selectedsire];

         if (ran1()<sexratio && tgdam < max_female) {
            tgdam++;
            female[tgdam] = totalna;
            pedholder[2][totalna] = 1;
         }
         else {
            if (tgsire < max_male) {
               tgsire++;
               male[tgdam] = totalna;
               pedholder[2][totalna] = 2;
            }
            else {
               tgdam++;
               female[tgdam] = totalna;
               pedholder[2][totalna] = 1;
            }
         }
      }
      tdam = tgdam;
      tsire = tgsire;
   }

   std::ofstream outfile(binfname.c_str(),std::ios::out);
   if (!outfile) throw exception("can't open internal file");
   Vector<unsigned> stdrec(3);
   basesize = 0;
   for (j=1; j<=totalna; j++) {
      stdrec[0] = j;
      stdrec[1] = pedholder[0][j];
      stdrec[2] = pedholder[1][j];
      stdrec[3] = pedholder[2][j];
      if (stdrec[1] == 0 && stdrec[2] == 0) basesize++;
      outfile.write((char *)stdrec.begin(),isdsexwidth);
   }
   outfile.close();
   ped_in_memory = 0;
   ped_on_disk = 1;
   nang = totalna;
   numgroup = 0;
   fdone       = 0;
   fvec.resize(totalna,0.0);
   maxnamelen = 0;
   sex_confirmed = 1;
   return *this;
}

void Pedigree::save_pedsheet(const int relse)
{
  //cout << "Save pedsheet\n";
   std::ofstream pedfile(binfname.c_str(),std::ios::out);
   if (!pedfile) throw exception("can't open internal file");
   PedNode *I;
   unsigned idssex[4];
   size_t idssexwidth = 4*sizeof(unsigned);
   unsigned i;
   for (i=0; i<totalna; i++) {
      I = pedmember[i];
      idssex[0] = I->myid; idssex[1] = I->mother_id; idssex[2] = I->father_id;
      idssex[3] = I->mysex;
      pedfile.write((char *)idssex,idssexwidth);
   }
   if (maxnamelen > 0) {
      for (i=0; i<nang; i++) {
         I = pedmember[i];
         pedfile.write((char *)&(I->namelen),sizeof(unsigned));
         pedfile.write(I->myname,I->namelen);
      }
   }
   pedfile.close();
   ped_on_disk = 1;
   if (relse) release_pedsheet();
}

void Pedigree::input_pedsheet(void)
{
   if (ped_in_memory) return;
   if (!ped_on_disk) throw exception("Pedigree::input_pedsheet(): Pedigre is not on disk");
   std::ifstream pedfile(binfname.c_str(),std::ios::in);
   if (!pedfile) throw exception("(1) most likely, pedigree has not been inputed\n"
            " (2) or probably, it can't open internal file");
   unsigned i,len;
   PedNode *I;
   if(nang>0){
     ped_storage = new PedNode[nang];
   }
   else {
     ped_storage = 0;
   }
   check_ptr(ped_storage);
   pedmember = new PedNode*[nang];
   check_ptr(pedmember);

   for (i=0; i<nang; i++) {
      I = &(ped_storage[i]);
      I->mysex = '.';
      I->myname = 0;
      I->namelen = 0;
      pedmember[i] = I;
   }
   unsigned idssex[4];
   size_t idswidth = 4*sizeof(unsigned);
   for (i=0; i<totalna; i++) {
      pedfile.read((char *)idssex,idswidth);
      I = pedmember[i];
      I->myid = idssex[0]; I->mother_id = idssex[1]; I->father_id = idssex[2];
      I->mysex = (char)idssex[3];
      //cout << "input_pedsheet "<< I << I->myid << " " << I->father_id << " " << I->mother_id <<"\n";
   }
   if (maxnamelen > 0) {
      for (i=0; i<nang; i++) {
         I = pedmember[i];
         pedfile.read((char *)&len,sizeof(unsigned));
	 if(len>0){
	   I->myname = new char [len];
	 }
	 else {
	   I->myname = 0;
	 }
         pedfile.read(I->myname,(size_t)len);
      }
   }
   pedfile.close();
   ped_in_memory = 1;                // pedigree data now is in memorry
}

void Pedigree::release_pedsheet(void)
{
   if (ped_storage) {
      for (unsigned i=0; i<nang; i++) {
         if (ped_storage[i].myname) {
	   if(ped_storage[i].myname){
	     delete [] ped_storage[i].myname;
	     ped_storage[i].myname=0;
	   }
	 }
      }
      if(ped_storage){
	delete [] ped_storage;
	ped_storage=0;
      }
      ped_storage = 0;
   }
   if (pedmember) {delete [] pedmember; pedmember = 0;}
   ped_in_memory = 0;   // pedigree is not in memory, but should  on disk;
}

void Pedigree::renum(void)
{
   if (totalna == 0) return;
   unsigned i,j,k,d,s;
   PedNode *I;
   Vector<unsigned> agen; agen.reserve(nang+1);
   agen[0] = 0;
   for (i=1; i<=totalna; i++) agen[i] = 1;
   for (i=(totalna+1); i<=nang; i++) agen[i] = 0;    // group is assigned 0

   int okk = 1;
   int g = 0;
   while (okk) {
      okk = 0;
      for (i=0; i<totalna; i++) {
         I = pedmember[i];
         d = agen[I->mother_id];
         s = agen[I->father_id];
         k = std::max(d,s);
         if (agen[i+1] <= k) {
            agen[i+1]= k+1;
            okk = I->myid;
         }
      }
      if (g++ >100) {
         std::cout << " A B C" << std::endl;
         std::cout << " B A D" << std::endl;
         std::cout << "Is it strange ?  but it's in your pedigree\n";
         std::cout << "check the families with " << pedmember[okk-1]->myname << std::endl;
         exit(1);
      }
   }
   unsigned  mg = agen.max();
   k = 0;
   for (i=1; i<=mg; i++) {
      for (j=0; j<totalna; j++) {
         if (agen[j+1] == i) pedmember[j]->myid = ++k;
      }
   }
   basesize = 0;
   if (!groupingped) {
      for (i=1; i<=totalna; i++) if (agen[i]==1) basesize++;
   }
   for (i=0; i<totalna; i++) {
      I = pedmember[i];
      d = I->mother_id; s = I->father_id;
      if (d>0 && d<=totalna) I->mother_id = pedmember[d-1]->myid;
      if (s>0 && s<=totalna) I->father_id = pedmember[s-1]->myid;
      //      std::cout << "renum " << i << " " <<I->father_id <<" " <<I->mother_id<< "\n";
   }
   qsort((char *)pedmember,totalna,sizeof(PedNode*),pedcp2);
}

void Pedigree::release(void)
{
   release_pedsheet();
   nang = 0;
   fdone = 0;
   totalna = 0;
   basesize = 0;
   maxnamelen = 0;
   ped_on_disk = 0;
   ped_in_memory = 0;
   remove(binfname.c_str());
}

Vector<double>* Pedigree::inbcoef_quaas(int keepinmem)
{
   if (totalna == 0) return &fvec;
   if (fdone) return &fvec;
   if (!ped_in_memory) input_pedsheet();

   unsigned i,j,s,d;
   PedNode *I,*J;
   double t;
   Vector<double> work(totalna);     // working space
   double *inbc = fvec.begin();
   for (j=0; j<totalna; j++) {
      t = 0.0;
      J = pedmember[j];
       d = J->mother_id; s = J->father_id;
      if (d> 0 && d <= totalna) t += inbc[d-1];
      if (s> 0 && s <= totalna) t += inbc[s-1];
      t = 1.0 - 0.25*t;
      inbc[j] += t;
      work[j] = std::sqrt(t);
      for (i=j+1; i<totalna; i++) {
         t = 0.0;
         I = pedmember[i];
         d = I->mother_id; s = I->father_id;
         if (d > j && j<=totalna) t += .5*work[d-1];
         if (s > j && j<=totalna) t += .5*work[s-1];
         work[i] = t;
         inbc[i] += t*t;
      }
   }
   for (i=0; i<totalna; i++) *inbc++ -= 1.0;
   fdone = 1;
   if (!keepinmem) release_pedsheet();
   return &fvec;
}

Vector<double>* Pedigree::inbcoef_meuwissen(int keepinmem)
{
   if (totalna == 0) return &fvec;
   if (fdone) return &fvec;
   if (!ped_in_memory) input_pedsheet();
   double *L;
   if(totalna>0){
     L = new double [totalna];
   }
   else {
     L = 0;
   }
   double *D;
   if(totalna>0){
     D = new double [totalna];
   }
   else {
     D = 0;
   }
   unsigned *point;
   if(totalna>0){
     point = new unsigned [totalna];
   }
   else {
     point = 0;
   }

   memset(L,'\0',sizeof(double)*totalna);
   memset(D,'\0',sizeof(double)*totalna);
   memset(point,'\0',sizeof(unsigned)*totalna);

   L--; D--; point--;
   double *F = fvec.begin()-1;      // now F starting from 1, no F[0] at all

   unsigned i,j,k,d,s,kd,ks,kk,previous_d,previous_s;
   double r,fi;
   PedNode *I;

   previous_s = 0;  previous_d = 0;
   for (i=1; i<=totalna; i++) {
      I = pedmember[i-1];
      d = I->mother_id;  if (d > totalna) d = 0;
      s = I->father_id;  if (s > totalna) s = 0;

      if (s == 0 || d == 0) {
         if (d+s == 0) {
            D[i] = 1.0;
         }
         else {
            D[i] = 0.75;
         }
         F[i] = 0.0;
      }
      else if ((d == previous_d) && (s == previous_s)) {
         F[i] = F[i-1];
         D[i] = D[i-1];
      }
      else {
         D[i] = 0.5 - 0.25*(F[d] + F[s]);
         fi = -1.0;
         L[i] = 1.0;
         j = i;
         while (j != 0) {
            k = j;
            r = 0.5*L[k];
            I = pedmember[k-1];
            kd = I->mother_id;
            ks = I->father_id;
            if (kd > ks) {kk = ks; ks = kd; kd = kk;}
            if (ks > 0  && ks <= totalna) {
               while (point[k] > ks) k = point[k];
               L[ks] += r;
               if (ks != point[k]) {
                  point[ks] = point[k];
                  point[k] = ks;
               }
               if (kd > 0 && kd <=totalna) {
                  while (point[k] > kd) k = point[k];
                  L[kd] += r;
                  if (kd != point[k]) {
                     point[kd] = point[k];
                     point[k] = kd;
                  }
               }
            }
            fi += L[j]*L[j]*D[j];
            L[j] = 0.0;
            k = j;
            j = point[j];
            point[k] = 0;
         }
         F[i] = fi;
      }
      previous_d = d;  previous_s = s;
   }
   L++;     
   if(L){
     delete [] L;
     L=0;
   }
   D++;     
   if(D){
     delete [] D;
     D=0;
   }
   point++; 
   if(point){
     delete [] point;
     point=0;
   }

   fdone = 1;
   if (!keepinmem) release_pedsheet();
   return &fvec;
}

double Pedigree::logdet(int keepinmem)
{
   ////////////////////////////////////////////////////////
   // ln of determinant of additive relationship matrix A
   ////////////////////////////////////////////////////////
   double retval = 0.0;
   if (totalna==0) throw exception("Pedigree::logdet(): empty pedigree");
   if (!fdone) inbcoef_meuwissen(1);
   double *F = fvec.begin()-1;      // now F starting from 1, no F[0] at all
   PedNode *I;
   unsigned father,mother;

   for (unsigned i=0; i<totalna; i++) {
      I = pedmember[i];
      mother = I->mother_id;
      father = I->father_id;
      if (mother && father) {
         retval += std::log((2.0-F[mother]-F[father])/4.0);
      }
      else if (mother && !father) {
         retval += std::log((3.0-F[mother])/4.0);
      }
      else if (!mother && father) {
         retval += std::log((3.0-F[father])/4.0);
      }
   }
   if (!keepinmem) release_pedsheet();
   return retval;
}

doubleMatrix Pedigree::ainv(int keepinmem)
{
   doubleMatrix ainvmat(totalna,totalna);
   if (totalna == 0) return ainvmat;
   if (!fdone) inbcoef(1);
   if (!ped_in_memory) input_pedsheet();

   Matrix<double> lambda(3,3);
   getlambda(lambda.begin(),3);
   double f,x;
   Vector<unsigned> ids(3);
   unsigned i,j,k,ii,jj;
   PedNode *I;

   for (k=0; k<totalna; k++) memset(ainvmat[k],'\0',sizeof(double)*totalna);
   for (k=0; k<totalna; k++) {
      I = pedmember[k];
      ids[0] = I->myid;
      ids[1] = I->mother_id;
      ids[2] = I->father_id;
      for (i=1; i<3; i++) if (ids[i]>totalna) ids[i] = 0;
      f = 2.0;
      if (ids[1] == 0) {f -= -1.0; }
      else             {f -= fvec[ids[1]-1]; }

      if (ids[2] == 0) {f -= -1.0; }
      else             {f -= fvec[ids[2]-1]; }

      f = 4.0/f;
      for (i=0; i<3; i++) {
         ii = ids[i];
         if (ii) {
            ii--;
            for (j=0; j<3; j++) {
               jj = ids[j];
               if (jj) {
                  jj--;
                  x = lambda[i][j]*f;
                  ainvmat[ii][jj] += x;
               }
            }
         }
      }
   }
   if (!keepinmem) release_pedsheet();
   return ainvmat;
}

doubleMatrix Pedigree::rela(int keepinmem)
{
   doubleMatrix amat(totalna,totalna);
   if (totalna == 0) return amat;
   if (!ped_in_memory) input_pedsheet();

   PedNode *I;
   unsigned i,j,a,d,s;
   for (i=0; i<totalna; i++) {
      memset(amat[i],'\0',sizeof(double)*totalna);
      I = pedmember[i];
      a = I->myid-1;
      d = I->mother_id;
      s = I->father_id;
      amat[a][a] = 1.0;
      if ( d>totalna ) d = 0;
      if ( s>totalna ) s = 0;
      if (d && s) {
         d--; s--;
         for (j=0; j<a; j++) {
            amat[j][a] = amat[a][j] = (amat[d][j]+amat[s][j])/2.0;
         }
         amat[a][a] += amat[d][s]/2.0;
      }
      else if (s && !d) {
         s--;
         for (j=0; j<a; j++) amat[j][a] = amat[a][j] = amat[s][j]/2.0;
      }
      else if (!s && d) {
         d--;
         for (j=0; j<a; j++) amat[j][a] = amat[a][j] = amat[d][j]/2.0;
      }
   }
   if (!keepinmem) release_pedsheet();
   return amat;
}

doubleMatrix Pedigree::reld(int keepinmem)
{
   doubleMatrix temp(totalna,totalna);

   if (totalna == 0) return temp;
   warning("P.reld() is appropriate only if pedigree is non-inbred");
   doubleMatrix A(totalna,totalna);
   A = rela(1);
   double **me = A.begin();
   PedNode *I,*J;
   unsigned i,j,di,si,dj,sj;

   for (i=0; i<totalna; i++) memset(temp[i],'\0',sizeof(double)*totalna);
   for (i=0; i<totalna-1; i++) {
      temp[i][i] = 1.0;
      I = pedmember[i];
      di = I->mother_id;
      si = I->father_id;
      if (di>0 && di<=totalna && si>0 && si<=totalna) {
         di--; si--;
         for (j=i+1; j<totalna; j++) {
            J = pedmember[j];
            dj = J->mother_id;
            sj = J->father_id;
            if ( dj>0 && dj<=totalna && sj>0 && sj<=totalna) {
               dj--; sj--;
               temp[j][i] = temp[i][j] = .25*(me[di][dj]*me[si][sj]+
                                              me[di][sj]*me[si][dj]);
            }
         }
      }
   }
   temp[totalna-1][totalna-1] = 1.0;
   if (!keepinmem) release_pedsheet();
   return temp;
}

doubleMatrix Pedigree::mat(void)
{
   doubleMatrix retval(totalna,3);
   if (totalna == 0) throw exception("Pedigree::mat(): empty pedigree objec");

   if (!ped_in_memory) input_pedsheet();
   double *dpt;
   PedNode *I;
   for (unsigned i=0; i<totalna; i++) {
      I = pedmember[i];
      dpt = retval[i];
      *dpt++ = static_cast<double>(I->myid);
      *dpt++ = static_cast<double>(I->mother_id);
      *dpt++ = static_cast<double>(I->father_id);
   }
   release_pedsheet();
   return retval;
}

void Pedigree::out_to_stream(std::ostream& os,const int style,const int ic)
{
   if (totalna == 0) {
      os << "\t empty pedigree object\n" << std::flush;
      return;
   }
   char line[10];
   unsigned i,k;
   int W = SESSION.output_precision+6;              // 6 = +.e+00
   char S = ' ';                              // one blank space
   PedNode *I;

   if (!ped_in_memory) input_pedsheet();
   if (maxnamelen > 0)  {
     if (style==0) {
         os << "   string_ID   renumed_ID    mother_ID    father_ID"
            << " sex";
         if (fdone) os << "      inbcoef";
         os << "\n";
         for (k=23,i=0; i<totalna; i++) {
            if (ic && i>=k) {
               k += 23;
               os << "  more ... [q for quit] ";
               std::cin.getline(line,sizeof(line));
               if (line[0]=='q') break;
            }
            I = pedmember[i];
            os << std::setw(maxnamelen) << I->myname
               << S << std::setw(W)<< I->myid
               << S << std::setw(W) << I->mother_id
               << S << std::setw(W) << I->father_id
               << S << std::setw(1) << I->mysex;
            if (fdone) os << S << std::setw(W) << fvec[i];
            os << "\n";
         }
      }
      else if (style==1) {
         for (k=23,i=0; i<totalna; i++) {
            if (ic && i>=k) {
               k += 23;
               os << "  more ... [q for quit] ";
               std::cin.getline(line,sizeof(line));
               if (line[0]=='q') break;
            }
            I = pedmember[i];
            os << std::setw(maxnamelen) << I->myname;
            if (I->mother_id) {
               os << S << std::setw(maxnamelen) << pedmember[I->mother_id-1]->myname;
            }
            else {
               os << S << std::setw(maxnamelen) << "0";
            }
            if (I->father_id) {
               os << S << std::setw(maxnamelen) << pedmember[I->father_id-1]->myname;
            }
            else {
               os << S << std::setw(maxnamelen) << "0";
            }
            os << S << std::setw(1) << I->mysex;
            if (fdone) os << S << std::setw(W) << fvec[i];
            os << "\n";
         }
      }
      else {
         warning("Pedigree::out_to_stream(,style): invalid style");
         return;
      }
   }
   else {
      if (style==0) {
         os << "  renumed_ID    mother_ID    father_ID sex";
         if (fdone) os << "      inbcoef";
         os << "\n";
      }
      for (k=23,i=0; i<totalna; i++) {
         if (ic && i>=k) {
            k += 23;
            os << "  more ... [q for quit] ";
            std::cin.getline(line,sizeof(line));
            if (line[0]=='q') break;
         }
         I = pedmember[i];
         os << std::setw(W)<< I->myid
            << S << std::setw(W) << I->mother_id << S << std::setw(W) << I->father_id
            << S << std::setw(1) << I->mysex;
         if (fdone) os << S << std::setw(W) << fvec[i];
         os << "\n";
      }
   }
   os << std::flush;
   release_pedsheet();
}

void Pedigree::save(const std::string &fname,const int style, const int io_mode)
{
   std::ofstream opedf;
   opedf.open(fname.c_str(),(OpenModeType)io_mode);
   if (!opedf) throw exception("Pedigree::save(): cannot open file");
   out_to_stream(opedf,style);
   opedf.close();
}

std::ostream& operator<<(std::ostream& stream, Pedigree& A)
{
   A.out_to_stream(stream);
   return stream;
}

void Pedigree::display(const std::string msg,const int style,const int ic)
{
   if (msg != "") std::cout << "\n" << msg << "\n" << std::endl;
   out_to_stream(std::cout,style,ic);
}

} /////// end of namespace matvec

