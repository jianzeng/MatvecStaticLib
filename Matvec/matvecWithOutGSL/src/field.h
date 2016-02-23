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

#ifndef MATVEC_FIELD_H
#define MATVEC_FIELD_H

#include <string>
#include "util.h"
#include "vector.h"
#include "datanode.h"
#include "fieldstruct.h"
#include "hashtable.h"

namespace matvec {
/*!
   \class Field Field.h
   \brief a column for a data set

   \sa Data DataNode
*/

class Field {
   friend class Data;
   protected:
      unsigned     ne;
      DataNode     mean_value;
      DataNode*    dat_vec;
      HashTable*   hashtable;

      void copyfrom(const Field& a);

   public:
      FieldStruct  col_struct;

      Field(void);                                    // Constructor 1
      Field(const unsigned n);                        // Constructor 2
      Field(const Field& a);                          // Constructor 3
      Field(const unsigned n,DataNode *a,FieldStruct st,
             HashTable* ht = 0);
      ~Field(void) {release();}                               // Destructor

      Field&          resize(const unsigned n);
      Field&          resize(const unsigned n, DataNode *a);  //dat_vec=0

      const Field&    operator=(const Field& a);

      Field&           operator+=(const Field& a);
      Field&           operator-=(const Field& a);
      Field&           operator*=(const Field& a);
      Field&           operator/=(const Field& a);

      Field&           operator+=(const double s);
      Field&           operator-=(const double s);
      Field&           operator*=(const double s);
      Field&           operator/=(const double s);

      DataNode&        operator()(const unsigned i);
      DataNode&        operator[](const unsigned i) {return dat_vec[i];}

      Field            operator!(void) const;
      Field            operator-(void) const;               //unary minus
      Field&           operator+(void) {return *this;};     //unary plus

      Vector<bool>           operator==(const Field &a) const;
      Vector<bool>           operator==(const double x) const;
      Vector<bool>           operator<(const Field &a) const;
      Vector<bool>           operator>(const Field &a) const;
      Vector<bool>           operator!=(const Field &a) const;
      Vector<bool>           operator<=(const Field &a) const;
      Vector<bool>           operator>=(const Field &a) const;

      friend Vector<bool>    operator==(const double a, const Field &b);
      friend Vector<bool>    operator<(const Field &a,const double b);
      friend Vector<bool>    operator<(const double a,const Field &b);
      friend Vector<bool>    operator!=(const Field &a,const double b);
      friend Vector<bool>    operator!=(const double a,const Field &b);
      friend Vector<bool>    operator>(const Field & a,const double b);
      friend Vector<bool>    operator>(const double a,const Field & b);
      friend Vector<bool>    operator<=(const Field &a,const double b);
      friend Vector<bool>    operator<=(const double a,const Field b);
      friend Vector<bool>    operator>=(const Field &a,const double b);
      friend Vector<bool>    operator>=(const double a,const Field &b);

      friend Field     operator+(const Field& a, const Field& b);
      friend Field     operator+(const Field& a, const double b);
      friend Field     operator+(const double b, const Field& a);

      friend Field     operator-(const Field& a, const Field& b);
      friend Field     operator-(const Field& a, const double b);
      friend Field     operator-(const double b, const Field& a);

      friend Field     operator / (const Field& a, const Field& b);
      friend Field     operator / (const Field& a, const double s);
      friend Field     operator / (const double s, const Field& b);

      friend Field     operator * (const Field& a, const Field& b);
      friend Field     operator * (const Field& a, const double b);
      friend Field     operator * (const double a, const Field& b);

      friend std::ostream&  operator<<(std::ostream& stream, const Field& a);

      friend Field     sin(Field& a);
      friend Field     asin(Field& a);
      friend Field     cos(Field& a);
      friend Field     acos(Field& a);
      friend Field     tan(Field& a);
      friend Field     atan(Field& a);
      friend Field     ceil(Field& a);
      friend Field     floor(Field& a);
      friend Field     log(Field& a);
      friend Field     log10(Field& a);
      friend Field     exp(Field& a);
      friend Field     sqrt(Field& a);
      friend Field     abs(Field& a);
      friend Field     erf(Field& a);
      friend Field     erfc(Field& a);

      void             value_for_missing(const double vm);
      void             out_to_stream(std::ostream& stream, const int ic) const;
      void             save(const std::string &ffname,
			    const int io_mode = std::ios::out) const;//SDK | ios::noreplace) const;
      void             display(const std::string &meg = "", const int ic=0) const;
      void             set_missing(const unsigned k);
      void             pretend_missing(const unsigned k);
      void             recover_missing(const unsigned k);
      Field            map(double (*f)(double)) const;
      Field            sub(const int i1,const int i2) const;
      Field&           zeros(void);
      Field&           ones(void);
      Field&           sort(void);
      DataNode         max(void) const;
      DataNode         min(void) const;
      DataNode         sum(void) const;
      DataNode         sumsq(void) const;
      DataNode         product(void) const;
      DataNode         mean(const int flag=1);
      DataNode         covariance(const Field *B=0) const;
      DataNode         elem(const unsigned i) const;

      unsigned         len(void) const {return ne;}
      unsigned         size(void) const {return ne;}
      void             release(void);

      void             index(const unsigned k)   {col_struct.index(k);}
      void             nlevel(const unsigned k)  {col_struct.nlevel(k);}
      void             nmiss(const unsigned k)   {col_struct.nmiss(k);}
      void             count_miss(const int k)   {col_struct.count_miss(k);}
      void             classi(const char c)      {col_struct.classi(c);}
      void             type(const char c)        {col_struct.type(c);}
      void             name(const std::string &n)      {col_struct.name(n);}
      unsigned         index(void)  const        {return col_struct.index();}
      unsigned         nlevel(void) const        {return col_struct.nlevel();}
      unsigned         nmiss(void)  const        {return col_struct.nmiss();}
      char             classi(void) const        {return col_struct.classi();}
      char             type(void)   const        {return col_struct.type();}
      const std::string      name(void)   const  {return col_struct.name();}

};

//public:

inline Field::Field(void)                            //Constructor 1
{ ne = 0; dat_vec = 0; hashtable = 0; }

inline Field::Field(const unsigned n)                //Constructor 2
{ 
  ne = n;  
  if(n>0){
    dat_vec = new DataNode [n]; 
  }
  else {
    dat_vec = 0;
  }
  hashtable = 0; 
}

inline DataNode& Field::operator()(const unsigned i)
{
   if (i - 1 >= ne ) throw exception("Field(): subscript out of range");
   return dat_vec[i - 1];
}

inline DataNode Field::elem(const unsigned i) const
{
   if (i - 1 >= ne) throw exception("Field.elem(): subscript out of range");
   return dat_vec[i - 1];
}

inline Field& Field::zeros(void)
{
   for (unsigned i=0; i<ne; i++) dat_vec[i].double_val(0.0);
   return *this;
}

inline Field& Field::ones(void)
{
   for (unsigned i=0; i<ne; i++) dat_vec[i].double_val(1.0);
   return *this;
}

} //////// end of namespace matvec

#endif
