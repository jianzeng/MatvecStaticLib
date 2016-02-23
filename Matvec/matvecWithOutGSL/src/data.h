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

#ifndef MATVEC_DATA_H
#define MATVEC_DATA_H

#include "datanode.h"
#include "field.h"

namespace matvec {
class doubleMatrix;
/*!
 *  \class Data  Data.h
 *  \brief A set of rows of data records
 *
 */

class Data {
   protected:
      int         data_on_disk, data_in_memory;
      unsigned    numcol,new_col,maxnumcol, numrec;
      std::string      tdfname;

      void       copyfrom(Data& A);

   public:
      HashTable    **hashtable;
      Field       *datasheet;

      Data(void);
      Data(Data& D);
      ~Data(void){release();}

      const Data&   operator=(Data& A);
      const Data&   operator=(const Field& V);

      Data&  resize(const unsigned nr,const unsigned nc,const unsigned mc=0);
      friend std::ostream&  operator<<(std::ostream& stream, Data& A);

      int        in_memory(void) const {return data_in_memory;}
      int        in_disk(void) const {return data_on_disk;}
      int        field_index(const std::string &colname) const;
      void       field_index_vec(Vector<int> &intvec,const std::string &fdname="");
      void       value_for_missing(const double vm);
//      void     symbol_for_missing(const char sm[]);
      void       input(const std::string & fname,const std::string & recfmt);
      void       save_datasheet(const int relse=1);
      void       input_datasheet(void);
      void       release_datasheet(void);
      void       release(void);
      void       print(std::ostream& stream,const Vector<int> intvec, const int ic=0);
      void       save(const std::string &fname,
                      const int io_mode = std::ios::out);//SDK|ios::noreplace);
      void       display(const std::string &fieldnames="",const int ic=0);
      void       newcol(const std::string &cname,const Field& col);
      void       row(const unsigned i,DataNode* recd);
      Field      col(const std::string &cname);
      Data&      newcol(const std::string & cname);
      Data&      stack(Data& b);
      Data&      adjoin(Data& b);

      DataNode*   cell(const unsigned r,const unsigned c);
      DataNode*   rawcol(unsigned c);
      DataNode*   rawcol(const std::string &cname);
      unsigned    size(void) const {return numrec;}
      unsigned    num_cols(void) const {return numcol;}
      unsigned    num_rows(void) const {return numrec;}
      Field       mean(const std::string &cname = "");
      Field       variance(const std::string &cname = "");
      Field       sum(const std::string &cname = "");
      Field       sumsq(const std::string &cname = "");
      Field       product(const std::string &cname = "");
      Field       max(const std::string &cname = "");
      Field       min(const std::string &cname = "");
      void        stat(void);
      doubleMatrix      mat(void);
};
}

#endif
