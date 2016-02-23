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

#ifndef DATANODE_H
#define DATANODE_H

namespace matvec {
/*!
  \class DataNode  DataNode.h
  \brief a node in a data structure

  \sa Data
*/
class DataNode {
   friend class Individual;
   friend class Field;
   protected:
      union {
         int       int_value;
         double    double_value;
         unsigned  unsigned_value;
      } data;

      void copyfrom(const DataNode& A);

   public:
      int  missing;

      DataNode(void) {missing=1;}              // missing is the default
      DataNode(const DataNode& A);             // copy constructor

      friend std::ostream&  operator<<(std::ostream& stream, const DataNode& A);

      const DataNode&  operator=(const DataNode& A);
      const DataNode&  operator=(const double x) {data.double_value = x;
                                                  missing = 0; return *this;}

      const DataNode&  operator=(const int x) {data.int_value = x; missing=0;
                                               return *this;}

      const DataNode&  operator=(const unsigned x) {data.unsigned_value=x;
                                                    missing=0; return *this;}
      const DataNode&  operator+=(const double x) {data.double_value += x;
                                                   return *this;}
      const DataNode&  operator-=(const double x) {data.double_value -= x;
                                                   return *this;}
      const DataNode&  operator*=(const double x) {data.double_value *= x;
                                                   return *this;}
      const DataNode&  operator/=(const double x) {data.double_value /= x;
                                                   return *this;}

      const DataNode&  operator+=(const int x) {data.int_value += x;
                                                return *this;}
      const DataNode&  operator-=(const int x) {data.int_value -= x;
                                                return *this;}
      const DataNode&  operator*=(const int x) {data.int_value *= x;
                                                return *this;}
      const DataNode&  operator/=(const int x) {data.int_value /= x;
                                                return *this;}

      const DataNode&  operator+=(const unsigned x) {data.unsigned_value += x;
                                                     return *this;}
      const DataNode&  operator-=(const unsigned x) {data.unsigned_value -= x;
                                                     return *this;}
      const DataNode&  operator*=(const unsigned x) {data.unsigned_value *= x;
                                                     return *this;}
      const DataNode&  operator/=(const unsigned x) {data.unsigned_value /= x;
                                                     return *this;}

      void     double_val(const double x) {data.double_value=x; missing=0;}
      void     int_val(const int x) {data.int_value = x; missing=0;}
      void     unsigned_val(const unsigned x){data.unsigned_value=x; missing=0;}
      double   double_val(void) const {return data.double_value;}
      int      int_val(void) const {return data.int_value;}
      unsigned unsigned_val(void) const {return data.unsigned_value;}
};

}
#endif
