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

#ifndef HashTable_H
#define HashTable_H

#include <cstdlib>
#include <iostream>

namespace matvec {
class HashTable;
/*!
   \class HashNode  hashtable.h
   \brief a node for the hashtable

   \sa HashTable
*/

class HashNode {
   friend  class HashTable;
   private:
      size_t   datasize;
      char     *data;
      unsigned id;
      void     copyfrom(const HashNode& A);

   public:
      HashNode(void);
      HashNode(const HashNode& A){data = 0; datasize = 0; copyfrom(A);}
      ~HashNode(void){release();}

      const    HashNode& operator=(const HashNode& A);
      int      equal(const char *v);
      void     insert(const char *v,const size_t ds,const unsigned idno);
      void     resize(const size_t s);
      unsigned id_no(void) {return id;}
      void     release(void);
};

HashNode* new_HashNode_vec(const unsigned n,const size_t s);

typedef void (*displaytype)(const void*);

/*!
   \class HashTable  HashTable.h
   \brief a hash table

*/

class HashTable {
    public:
       enum hashaction {INSERT,GETIDNO};
       HashTable(void);
       HashTable(const unsigned n, const size_t s);
       HashTable(const unsigned n, const char str[]);
       HashTable(const HashTable& A){copyfrom(A);}
       ~HashTable(void){release();}

       const HashTable& operator=(const HashTable& A);

       HashTable&  resize(const unsigned n, const size_t s);
       HashTable&  resize(const unsigned n) {return resize(n,0);}

       unsigned    insert(const void *v); // return its id of v in HashTable
       size_t      data_len(void){return datasize;}
       unsigned    get_id(const void *v);
       unsigned    size(void) const {return act_tablesize;}
       void        change_id(const unsigned oldid,const unsigned newid);
       void        reorder(void);
       void        maxsize(const unsigned ms);
       void        release(void);
       void        save(const char *fname,
                        const int io_mode=std::ios::out); //SDK ios::noreplace is missing 
       void        input(const char fname[]);
       void        save_to_disk(std::ostream& stream,const int relse=1);
       void        input_from_disk(std::istream& stream);
       void        display(displaytype dply);
       void        copyfrom(const HashTable& A);

       const void* find(const unsigned id) const;

   private:
      unsigned  tablesize,ext_tablesize,act_tablesize,datasize, id_changed;
      HashNode  **hash_table;
      HashNode  *hash_storage;

   protected:
      unsigned hash(const char* v, HashTable::hashaction action);
};

}
#endif
