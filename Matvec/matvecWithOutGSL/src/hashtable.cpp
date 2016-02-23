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

#define HASH_CONST   1482907	      // prime close to 2^{20.5}
#include <fstream>
#include "exception.h"
#include "util.h"
#include "hashtable.h"

namespace matvec {

HashNode::HashNode(void)
{
   data = 0;
   datasize = 0;
   id = 0;
}

void HashNode::copyfrom(const HashNode& A)
{
   if (this == &A) return;
   release();
   if (A.data) {
      datasize = A.datasize;
      if(datasize>0){
	data = new char [datasize];
      }
      else {
	data = 0;
      }
      memcpy(data,A.data,datasize);
      id = A.id;
   }
}

void HashNode::resize(const size_t s)
{
  if (datasize == s) return;
  datasize = s;
  if (data) {
    delete [] data;
    data=0;
  }
  if(datasize>0){
    data = new char [datasize];
  }
  else {
    data = 0;
  }
  id = 0;
}

const HashNode& HashNode::operator=(const HashNode& A)
{
   copyfrom(A);
   return *this;
}

int HashNode::equal(const char *v)
{
   for (unsigned i=0; i<datasize; i++) if (data[i] != v[i]) return 0;
   return 1;
}

void HashNode::insert(const char *v,const size_t ds,const unsigned idno)
   // if nodesize is constant, then memory should has been allocated already
{
   if (datasize != ds) resize(ds);
   memcpy(data,v,datasize);
   id = idno;
}

void HashNode::release(void)
{
   if (data) { delete [] data; data = 0;}
   datasize=0;
}

//////////////////////////////////////////////////////////
//   this routine was taken from bibindex.c by
//    Nelson H. F. Beebe <beebe@math.utah.edu>
//  long next_prime(long n)
//  Return the next prime number after n.
/////////////////////////////////////////////////////////
long next_prime(long n)
{
   long prime;                  // tentative prime
   long factor;                 // prime factor
   int is_prime;		// 'prime' is a prime number

   n = (n < 0L) ? -n : n;       // be safe -- force n positive
   prime = 2L*(n/2L) + 1L;      // next odd number
   is_prime = (prime <= 3L);
   while (!is_prime) {
      factor = 5L;
      is_prime = (prime % 2L) && (prime % 3L);
      while (is_prime && (factor*factor <= prime)) {
	 if ((prime % factor) == 0L)
            is_prime = 0;
	 else if ((prime % (factor + 2L)) == 0L)
            is_prime = 0;
         else               // factor+4 is divisible by 3 (every 3rd odd is)
            factor += 6L;
      }
      if (!is_prime) prime += 2L;
   }
   return (prime);
}

HashTable::HashTable(void)
{
   datasize      = 0;
   tablesize     = 0;
   id_changed    = 0;
   ext_tablesize = 0;
   act_tablesize = 0;

   hash_table    = 0;
   hash_storage  = 0;
}

HashTable::HashTable(const unsigned n, const size_t s)
{
   datasize      = 0;
   tablesize     = 0;
   id_changed    = 0;
   ext_tablesize = 0;

   hash_table    = 0;
   hash_storage  = 0;
   resize(n,s);
}

HashTable::HashTable(const unsigned n, const char str[])
{

   // str could be any string
   datasize      = 0;
   tablesize     = 0;
   id_changed    = 0;
   ext_tablesize = 0;

   hash_table    = 0;
   hash_storage  = 0;
   resize(n,0);
}

void HashTable::copyfrom(const HashTable& A)
{
   if (this == &A) return;
   resize(A.tablesize,A.datasize);
   maxsize(A.ext_tablesize);
   id_changed = A.id_changed;
   act_tablesize = A.act_tablesize;
   HashNode *node;
   for (unsigned i=0; i<ext_tablesize; i++) {
      node = &A.hash_storage[i];
      if (node->data) {
         hash_storage[i] = *node;
         hash_table[node->id_no()-1] = &hash_storage[i];
      }
   }
}

const HashTable& HashTable::operator=(const HashTable& A)
{
   copyfrom(A);
   return *this;
}

unsigned HashTable::hash(const char* v, HashTable::hashaction action)
{
   register int i;
   unsigned long h,new_h,skip;
   unsigned strsize;
   HashNode *node;

   if (datasize==0) {
      strsize = strlen(v)+1;
      if (strsize==1) {
         warning("HashTable.hash(v): empty v or zero size of v");
         return 0;
      }
   }
   else strsize = datasize;
   for (h=0, skip=1, i=0; i<strsize; i++) {
      h = (h*HASH_CONST + v[i]) % ext_tablesize;
      skip += 2*h;
   }
   node = &hash_storage[h];
   for (int ntry=0; ntry<500; ntry++) {
      if ( node->id == 0) {                // node is empty
         if (action == HashTable::INSERT) {
            hash_table[act_tablesize] = node;
            act_tablesize++;
            node->insert(v,strsize,act_tablesize);
            return act_tablesize;
         }
         else if (action == HashTable::GETIDNO) {
            return 0;
         }
         else {
            warning("HashTable.hash(), unknown hash action type");
         }
      }
      else {   // node is not empty
         if (node->equal(v)) {
            return node->id_no();     // for both INSERT and GETIDNO
         }
         else {
            new_h = (h + skip) % ext_tablesize;
            h = (new_h == h) ? (h+1)% ext_tablesize : new_h;
            node = &hash_storage[h];
         }
      }
   }
   return 0;
}

unsigned HashTable::insert(const void* vv)
{
   if (ext_tablesize==0) {
      warning("HashTable.insert(): HashTable is null, do resize() first");
      return 0;
   }
   const char *v = (const char *)vv;
   unsigned id;
   id = hash(v,HashTable::INSERT);   // return its id of vv in HashTable
   if (id == 0) {
      unsigned k = act_tablesize;
      HashTable TMP;
      TMP.copyfrom(*this);
      maxsize(static_cast<unsigned>(next_prime(static_cast<long>(ext_tablesize+51))));
      for (unsigned i=0; i<k; i++) {
         id = hash(TMP.hash_table[i]->data,HashTable::INSERT);
         if (id == 0) throw exception("hashtable size too small");
      }
      id = hash(v,HashTable::INSERT);
      if (id == 0) throw exception("hashtable size too small");
   }
   return id;
}

void HashTable::change_id(const unsigned oldid,const unsigned newid)
{
   ///////////////////////////////////////////////////////////////////////
   //  once calling this routine, you must call it for each non-empty
   //  HashNode (id starts from 1 to act_tablesize). newif must starts
   //  1 to act_tablesize, too.
   //  After calling for each non-empty HashNode, then must call reorder()
   //  to re-order hash_table so that hash_table[i] points HashNode
   //  with id i+1 where i starts from 0 to act_tablesize.
   ///////////////////////////////////////////////////////////////////////
   if (oldid >act_tablesize || newid>act_tablesize) throw exception("HashTable.change_id(oldid,newid): invalid oldid or newid");
   hash_table[oldid-1]->id = newid;
   id_changed = 1;
}

void HashTable::reorder(void)
{
   if (id_changed == 0) return;
   unsigned i,id;
   for (i=0; i<ext_tablesize; i++) {
      id = hash_storage[i].id_no();
      if (id > 0) hash_table[id-1] = &(hash_storage[i]);
   }
}

unsigned HashTable::get_id(const void *vv)
{
   const char *v = (const char *)vv;
   return hash(v,HashTable::GETIDNO); // if 0, means it cannot found
}

const void* HashTable::find(const unsigned id) const
{
   if (id == 0 || id > act_tablesize) {
     throw exception("HashTable.find(id): range error");
   }
   return (const void *)(hash_table[id-1]->data); // index starts from 0
}

HashTable& HashTable::resize(const unsigned n, const size_t s)
{
   if (n==0) {
      release();
      return *this;
   }
   unsigned i;
   unsigned tmpsize = static_cast<unsigned>(next_prime(static_cast<long>(1.30*n + 50.0)));
                                                    //  30% of tablesize
   if (tmpsize != ext_tablesize || datasize != s) {   // definitely delete them
      ext_tablesize = tmpsize;
      datasize = s;
      if(hash_storage){
	delete [] hash_storage;
	hash_storage=0;
      }
      if (ext_tablesize>0){
	hash_storage = new HashNode [ext_tablesize];
      }
      else{
	hash_storage = 0;
      }
      check_ptr(hash_storage);
      for (int i=0; i<ext_tablesize; ++i) hash_storage[i].resize(datasize);
      if(hash_table){
	delete [] hash_table;
	hash_table=0;
      }
      if(ext_tablesize>0){
	hash_table = new HashNode*[ext_tablesize];
      }
      else {
	hash_table = 0;
      }
      check_ptr(hash_table);
   }
   else {
      for (i=0; i<act_tablesize; i++) hash_table[i]->id = 0;
   }
   tablesize = n;
   act_tablesize = 0;
   return *this;
}

void  HashTable::maxsize(const unsigned ms)
{
   if (ext_tablesize == ms) return;
   ext_tablesize = ms;
   if(hash_storage){
     delete [] hash_storage;
     hash_storage=0;
   }
   if (ext_tablesize>0){
     hash_storage = new HashNode [ext_tablesize];
   }
   else{
     hash_storage = 0;
   }
   check_ptr(hash_storage);
   for (int i=0; i<ext_tablesize; ++i) hash_storage[i].resize(datasize);

   if(hash_table){
     delete [] hash_table;
     hash_table=0;
   }
   if(ext_tablesize>0){
     hash_table = new HashNode*[ext_tablesize];
   }
   else {
     hash_table = 0;
   }
   check_ptr(hash_table);
   act_tablesize = 0;
}

void HashTable::release(void)
{
   if (hash_storage) {delete [] hash_storage; hash_storage = 0;}
   if (hash_table)   {delete [] hash_table; hash_table = 0;}
   tablesize     = 0;
   datasize      = 0;
   ext_tablesize = 0;
}

void HashTable::save(const char fname[],const int io_mode)
{
   std::ofstream hfile;
   hfile.open(fname,(OpenModeType)io_mode);
   if (!hfile) throw exception(" HashTable::save(): cannot open file");
   save_to_disk(hfile,0);
   hfile.close();
}

void HashTable::input(const char fname[])
{
   std::ifstream hfile(fname,std::ios::in);
   if (!hfile) throw exception(" HashTable::input(): cannot open file");
   input_from_disk(hfile);
   hfile.close();
}

/**********************************************************************
* it save the hash_storage and other necessary information into file hfile
* so that they can be loaded into other HashTable objects.
*  hash_storage and hash_table have been released from memory
***********************************************************************/
void HashTable::save_to_disk(std::ostream& stream,const int relse)
{
   unsigned ts,act_ts,ext_ts,tpos;
   size_t   ds;
   ts = tablesize;
   act_ts = act_tablesize;
   ext_ts = ext_tablesize;
   ds = datasize;

   stream.write((char*)&ts,sizeof(unsigned));
   stream.write((char*)&act_ts,sizeof(unsigned));
   stream.write((char*)&ext_ts,sizeof(unsigned));
   stream.write((char*)&ds,sizeof(size_t));

   HashNode *node;
   for (tpos=0; tpos<ext_tablesize; tpos++) {
      node = &hash_storage[tpos];
      if (node->data) {
         stream.write((char*) &tpos,sizeof(unsigned));
         stream.write((char*) &(node->id),sizeof(unsigned));
         stream.write((char*) &(node->datasize),sizeof(size_t));
         stream.write((char*)node->data,node->datasize);
      }
   }
   if (relse) release();
}

void HashTable::input_from_disk(std::istream& stream)
{
   unsigned ts,act_ts,ext_ts,tpos,id;
   size_t ds = 0;
   ts = act_ts = ext_ts = tpos = id = 0;
   stream.read((char*) &ts,sizeof(unsigned));
   stream.read((char*) &act_ts,sizeof(unsigned));
   stream.read((char*) &ext_ts,sizeof(unsigned));
   stream.read((char*) &ds,sizeof(size_t));

   resize(ts,ds);       // allocate momery for hash_storage, hash_table
   maxsize(ext_ts);
   act_tablesize = act_ts;

   HashNode *node;
   for (ts=0; ts<act_tablesize; ts++) {
      stream.read((char*) &tpos,sizeof(unsigned));
      stream.read((char*) &(id),sizeof(unsigned));
      stream.read((char*) &(ds),sizeof(size_t));
      node = &hash_storage[tpos];
      node->id = id;
      if (datasize == 0) {
         node->datasize = ds;
	 if(ds>0){
	   node->data = new char [ds];
	 }
	 else {
	   node->data = 0;
	 }
      }
      stream.read((char*)node->data,ds);
      hash_table[id-1] = node;
   }
}

void HashTable::display(displaytype dply)
{
   for (unsigned i=0; i<act_tablesize; i++) dply(hash_table[i]->data);
}
}
#undef HASH_CONST
