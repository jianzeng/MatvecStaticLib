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

#include <cmath>
#include <iomanip>
#include "field.h"

extern "C"{
extern double erf(double x);
extern double erfc(double x);
}

namespace matvec {

Field::Field(const unsigned n, DataNode *a,FieldStruct st,HashTable* ht)      {
   ne = n;
   dat_vec = a;
   col_struct = st;
   hashtable = ht;
}

Field::Field(const Field& a)     //Constructor 3
{
   ne = 0; dat_vec = 0; hashtable = 0;
   copyfrom(a);
}

void Field::release(void)
{
   if (dat_vec)  { delete [] dat_vec; dat_vec = 0; }
   if (hashtable) { delete hashtable; hashtable = 0; }
   ne = 0;
}

void Field::copyfrom(const Field& A)
{
   if (this == &A) return;
   resize(A.ne);
   memcpy(dat_vec,A.dat_vec,sizeof(DataNode)*ne);
   col_struct = A.col_struct;
   if (A.hashtable) {
      if (hashtable) {
	delete hashtable;
	hashtable=0;
      }
      hashtable = new HashTable;
      *hashtable = *(A.hashtable);
   }
}

const Field& Field::operator=(const Field& a)
{
   copyfrom(a);
   return *this;
}

Field& Field::resize(const unsigned n)
{
   if (ne == n) return *this;
   if (dat_vec) {
     delete [] dat_vec;
     dat_vec=0;
   }
   ne = n;
   if (ne>0){
     dat_vec = new DataNode [ne];
   }
   else {
     dat_vec = 0;
   }
     
   if (hashtable) {delete hashtable; hashtable = 0;}
   return *this;
}

Field& Field::resize(const unsigned n, DataNode *A)
{
   if (dat_vec) {
     delete [] dat_vec;
     dat_vec=0;
   }
   ne = n;
   dat_vec = A;
   return *this;
}


Field& Field::operator+=(const Field& a)
{
   if (ne == a.ne ) {
      if (col_struct.type() == 'S') {
         warning("Field += Field: string column");
      }
      else {
         for (unsigned i=0; i<ne; i++) {
            if (!a.dat_vec[i].missing) {
               dat_vec[i].data.double_value += a.dat_vec[i].data.double_value;
            }
         }
      }
   }
   else {
      warning("operator +=, size not conformable");
   }
   return *this;
}

Field& Field::operator-=(const Field& a)
{
   if (ne == a.ne ) {
      if (col_struct.type() == 'S') {
         warning("Field -= Field: string column");
      }
      else {
         for (unsigned i=0; i<ne; i++) {
            if (!a.dat_vec[i].missing) {
               dat_vec[i].data.double_value -= a.dat_vec[i].data.double_value;
            }
         }
      }
   }
   else {
      warning("operator -=, Fields not conformable");
   }
   return *this;
}

Field& Field::operator*=(const Field& a)
{
   if (ne == a.ne ) {
      if (col_struct.type() == 'S') {
         warning("Field *= Field: string column");
      }
      else {
         for (unsigned i=0; i<ne; i++) {
            if (!a.dat_vec[i].missing) {
               dat_vec[i].data.double_value *= a.dat_vec[i].data.double_value;
            }
         }
      }
   }
   else {
      warning("operator *=, Fields not conformable");
   }
   return *this;
}

Field& Field::operator/=(const Field& a)
{
   if (ne != a.ne ) {
      if (col_struct.type() == 'S') {
         warning("Field /= Field: string column");
      }
      else {
         for (unsigned i=0; i<ne; i++) {
            if (!a.dat_vec[i].missing) {
               dat_vec[i].data.double_value /= a.dat_vec[i].data.double_value;
            }
         }
      }
   }
   else {
      warning("operator /=, Fields not conformable");
   }
   return *this;
}

Field& Field::operator+=(const double s)
{
   if (col_struct.type() == 'S') {
      warning("Field += s: failed for string column");
   }
   else {
      for (unsigned i=0; i<ne; i++) dat_vec[i].data.double_value += s;
   }
   return *this;
}
Field& Field::operator-=(const double s)
{
   if (col_struct.type() == 'S') {
      warning("Field -= s: failed for string column");
   }
   else {
      for (unsigned i=0; i<ne; i++) dat_vec[i].data.double_value -= s;
   }
   return *this;
}

Field& Field::operator*=(const double s)
{
   if (col_struct.type() == 'S') {
      warning("Field *= s: failed for string column");
   }
   else {
      for (unsigned i=0; i<ne; i++) {
         dat_vec[i].data.double_value *= s;
      }
   }
   return *this;
}

Field& Field::operator/=(const double s)
{
   if (col_struct.type() == 'S') {
      warning("Field /= s: failed for string column");
   }
   else {
      for (unsigned i=0; i<ne; i++) dat_vec[i].data.double_value /= s;
   }
   return *this;
}

Vector<bool> Field::operator==(const Field& a) const
{
   Vector<bool> temp(ne);
   if (ne != a.len() ) {
      warning("Field::operator==(), size not conformable");
      return temp;
   }

   if (col_struct.type() != a.col_struct.type() ) return temp;
   if (col_struct.nlevel() != a.col_struct.nlevel()) return temp;

   if (col_struct.type() == 'S') {
      char *str, *a_str;
      for (unsigned i=0; i<ne; i++) {
         if (!dat_vec[i].missing && !a.dat_vec[i].missing) {
            str = (char *)hashtable->find(dat_vec[i].unsigned_val());
            a_str = (char *)a.hashtable->find(dat_vec[i].unsigned_val());
            if (strcmp(str,a_str)==0) temp[i] = true;
         }
      }
   }
   else if (col_struct.type() == 'F') {
      for (unsigned i=0; i<ne; i++) {
         if (dat_vec[i].double_val()==a.dat_vec[i].double_val()) temp[i]= true;
      }
   }
   return temp;
}

Vector<bool> Field::operator==(const double x) const
{
   Vector<bool> temp(ne);

   if (col_struct.type() == 'S') return temp;
   for (unsigned i=0; i<ne; i++) {
      if (dat_vec[i].double_val()==x) temp[i] = true;
   }
   return temp;
}

Vector<bool> Field::operator!=(const Field& a) const {return !(*this == a);}
Vector<bool> Field::operator<=(const Field& a) const {return !(*this > a);}
Vector<bool> Field::operator>=(const Field& a) const {return !(*this < a);}

Vector<bool> Field::operator<(const Field& a) const
{
   Vector<bool> temp(ne);
   if (ne != a.len() ) {
      warning("A==B, two Fields are not conformable");
      return temp;
   }

   if (col_struct.type() != a.col_struct.type() ) return temp;
   if (col_struct.nlevel() != a.col_struct.nlevel()) return temp;

   if (col_struct.type() == 'S') {
      char *str, *a_str;
      for (unsigned i=0; i<ne; i++) {
         if (!dat_vec[i].missing && !a.dat_vec[i].missing) {
            str = (char *)hashtable->find(dat_vec[i].unsigned_val());
            a_str = (char *)a.hashtable->find(dat_vec[i].unsigned_val());
            if (strcmp(str,a_str) < 0) temp[i] = true;
         }
      }
   }
   else if (col_struct.type() == 'F') {
      for (unsigned i=0; i<ne; i++) {
         if (dat_vec[i].double_val()<a.dat_vec[i].double_val()) temp[i]= true;
      }
   }
   return temp;

}

Vector<bool> Field::operator > (const Field& a) const
{
   Vector<bool> temp(ne);
   if (ne != a.len() ) {
      warning("A==B, two Fields are not conformable");
      return temp;
   }

   if (col_struct.type() != a.col_struct.type() ) return temp;
   if (col_struct.nlevel() != a.col_struct.nlevel()) return temp;

   if (col_struct.type() == 'S') {
      char *str, *a_str;
      for (unsigned i=0; i<ne; i++) {
         if (!dat_vec[i].missing && !a.dat_vec[i].missing) {
            str = (char *)hashtable->find(dat_vec[i].unsigned_val());
            a_str = (char *)a.hashtable->find(dat_vec[i].unsigned_val());
            if (strcmp(str,a_str) > 0) temp[i] = true;
         }
      }
   }
   else if (col_struct.type() == 'F') {
      for (unsigned i=0; i<ne; i++) {
         if (dat_vec[i].double_val()>a.dat_vec[i].double_val()) temp[i]= true;
      }
   }
   return temp;
}

Vector<bool> operator==(const double a, const Field& b) {return (b==a);}
Vector<bool> operator!=(const Field& a,const double b) {return !(a==b);}
Vector<bool> operator!=(const double a,const Field& b) {return !(b==a); }
Vector<bool> operator>(const Field& a,const double b) {return (b<a);}
Vector<bool> operator>(const double a,const Field& b) {return (b<a);}
Vector<bool> operator<=(const Field& a,const double b) {return !(b<a);}
Vector<bool> operator<=(const double a,const Field& b) {return !(b<a);}
Vector<bool> operator>=(const Field& a,const double b) {return !(a<b);}
Vector<bool> operator>=(const double a,const Field& b) {return !(a<b);}

Vector<bool> operator<(const Field& A,const double x)
{
   unsigned n = A.ne;
   Vector<bool> temp(n);

   if (A.col_struct.type() == 'S') return temp;
   for (unsigned i=0; i<n; i++) {
      if (A.dat_vec[i].double_val() < x) temp[i] = true;
   }
   return temp;
}

Vector<bool> operator<(const double x,const Field& A)
{
   unsigned n = A.ne;
   Vector<bool> temp(n);

   if (A.col_struct.type() == 'S') return temp;
   for (unsigned i=0; i<n; i++) {
      if (x < A.dat_vec[i].double_val()) temp[i] = true;
   }
   return temp;
}

Field operator+(const Field& A,const Field& B)
{
   unsigned n = A.ne;
   if (n != B.ne) {
      warning("Col + Col: size not conformable");
      return Field();
   }
   if (A.col_struct.type()=='S' || B.col_struct.type()=='S') {
      warning("Col + Col: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(n>0){
     temp = new DataNode [n];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<n; i++) {
      if (A.dat_vec[i].missing || B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val() +
                            B.dat_vec[i].double_val());
      }
   }
   return Field(n,temp,A.col_struct);
}

Field operator-(const Field& A,const Field& B)
{
   unsigned n = A.ne;
   if (n != B.ne) {
      warning("Col - Col: size not conformable");
      return Field();
   }
   if (A.col_struct.type()=='S' || B.col_struct.type()=='S') {
      warning("Col - Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(n>0){
     temp = new DataNode [n];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<n; i++) {
      if (A.dat_vec[i].missing || B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val()-
                            B.dat_vec[i].double_val());
      }
   }
   return Field(n,temp,A.col_struct);
}

Field operator*(const Field& A,const Field& B)
{
   unsigned n = A.ne;
   if (n != B.ne) {
      warning("Col * Col: size not conformable");
      return Field();
   }
   if (A.col_struct.type()=='S' || B.col_struct.type()=='S') {
      warning("Col * Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(n>0){
     temp = new DataNode [n];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<n; i++) {
      if (A.dat_vec[i].missing || B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val()*
                            B.dat_vec[i].double_val());
      }
   }
   return Field(n,temp,A.col_struct);
}

Field operator/(const Field& A,const Field& B)
{
   unsigned n = A.ne;
   if (n != B.ne) {
      warning("Col / Col: size not conformable");
      return Field();
   }
   if (A.col_struct.type()=='S' || B.col_struct.type()=='S') {
      warning("Field / Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(n>0){
     temp = new DataNode [n];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<n; i++) {
      if (A.dat_vec[i].missing || B.dat_vec[i].missing){
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val()/
                            B.dat_vec[i].double_val());
      }
   }
   return Field(n,temp,A.col_struct);
}

Field operator+(const Field& A,const double s)
{
   if (A.col_struct.type()=='S') {
      warning("Field + s: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(A.ne>0){
     temp = new DataNode [A.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<A.ne; i++) {
      if (A.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val() + s);
      }
   }
   return Field(A.ne,temp,A.col_struct);
}

Field operator-(const Field& A,const double s)
{
   if (A.col_struct.type()=='S') {
      warning("Field - s: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(A.ne>0){
     temp = new DataNode [A.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<A.ne; i++) {
      if (A.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val() - s);
      }
   }
   return Field(A.ne,temp,A.col_struct);
}

Field operator*(const Field& A,const double s)
{
   if (A.col_struct.type()=='S') {
      warning("Field * s: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(A.ne>0){
     temp = new DataNode [A.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<A.ne; i++) {
      if (A.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val() * s);
      }
   }
   return Field(A.ne,temp,A.col_struct);
}

Field operator/(const Field& A,const double s)
{
   if (A.col_struct.type()=='S') {
      warning("Field / s: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(A.ne>0){
     temp = new DataNode [A.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<A.ne; i++) {
      if (A.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(A.dat_vec[i].double_val() / s);
      }
   }
   return Field(A.ne,temp,A.col_struct);
}

Field operator+(const double s,const Field& B)
{
   if (B.col_struct.type()=='S') {
      warning("s + Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(B.ne>0){
     temp = new DataNode [B.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<B.ne; i++) {
      if (B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(s + B.dat_vec[i].double_val());
      }
   }
   return Field(B.ne,temp,B.col_struct);
}

Field operator-(const double s,const Field& B)
{
   if (B.col_struct.type()=='S') {
      warning("s - Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(B.ne>0){
     temp = new DataNode [B.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<B.ne; i++) {
      if (B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(s - B.dat_vec[i].double_val());
      }
   }
   return Field(B.ne,temp,B.col_struct);
}

Field operator*(const double s,const Field& B)
{
   if (B.col_struct.type()=='S') {
      warning("s * Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(B.ne>0){
     temp = new DataNode [B.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<B.ne; i++) {
      if (B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val( s * B.dat_vec[i].double_val());
      }
   }
   return Field(B.ne,temp,B.col_struct);
}
Field operator/(const double s,const Field& B)
{
   if (B.col_struct.type()=='S') {
      warning("s / Col: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(B.ne>0){
     temp = new DataNode [B.ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<B.ne; i++) {
      if (B.dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(s / B.dat_vec[i].double_val());
      }
   }
   return Field(B.ne,temp,B.col_struct);
}

Field Field::operator!(void) const  //unary minus
{
   if (col_struct.type()=='S') {
      warning("!Field: failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(ne>0){
     temp = new DataNode [ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<ne; i++) {
      if (dat_vec[i].missing) {
         temp[i].missing = 1;
      }
      else {
         temp[i].double_val(-(dat_vec[i].double_val()));
      }
   }
   return Field(ne,temp,col_struct);
}

Field Field::map(double (*f)(double)) const
{
   if (col_struct.type()=='S') {
      warning("Field:map(): failed for string column");
      return Field();
   }
   DataNode *temp; 
   if(ne>0){
     temp = new DataNode [ne];
   }
   else {
     temp = 0;
   }
   for (unsigned i=0; i<ne; i++) {
      if (dat_vec[i].missing)  {temp[i].missing = 1; }
      else { temp[i].double_val((*f)(dat_vec[i].double_val())); }
   }
   return Field(ne,temp,col_struct);
}

Field Field::sub(const int i1,const int i2) const
{
   if (ne == 0) return Field();
   int newne = i2 - i1 + 1;
   if ( newne >ne || i2 > ne) {
      warning("Field::s(), subscript out of range");
      return Field();
   }
   unsigned id,i,k;
   HashTable *tmp_hashtable = 0;
   FieldStruct tmp_struct = col_struct;

   DataNode *temp; 
   if(newne>0){
     temp = new DataNode [newne];
   }
   else {
     temp = 0;
   }
   if (col_struct.type()=='S') {
      tmp_hashtable = new HashTable;
      tmp_hashtable->resize(col_struct.nlevel());
      for (i=i1-1; i<i2; i++) {
         tmp_hashtable->insert((char *)hashtable->find(i+1));
      }
      id = tmp_hashtable->size();
      tmp_struct.nlevel(id);
      tmp_hashtable->resize(id);
      for (k=0,i=i1-1; i<i2; i++) {
         id = tmp_hashtable->insert((char *)hashtable->find(i+1));
         temp[k++].unsigned_val(id);
      }
   }
   else {
      memcpy(temp,&dat_vec[i1-1],sizeof(DataNode)*newne);
   }
   return Field(newne,temp,tmp_struct,tmp_hashtable);
}

Field sin(Field& a){return a.map(std::sin);}
Field asin(Field& a){return a.map(std::asin);}
Field cos(Field& a){return a.map(std::cos);}
Field acos(Field& a){return a.map(std::acos);}
Field tan(Field& a){return a.map(std::tan);}
Field atan(Field& a){return a.map(std::atan);}
Field ceil(Field& a){return a.map(std::ceil);}
Field floor(Field& a){return a.map(std::floor);}
Field log(Field& a){return a.map(std::log);}
Field log10(Field& a){return a.map(std::log10);}
Field exp(Field& a){return a.map(std::exp);}
Field sqrt(Field& a){return a.map(std::sqrt);}
Field abs(Field& a){return a.map(std::fabs);}
Field erf(Field& a){return a.map(::erf);}
Field erfc(Field& a){return a.map(::erfc);}

DataNode Field::max(void) const
{
   DataNode retval;
   retval.missing = 1;

   if (col_struct.type() != 'S') {
      unsigned i,k=0;
      while (dat_vec[k].missing) k++;
      double y,mx = dat_vec[k].double_val();
      for (i=k; i<ne; i++) {
         if (!dat_vec[i].missing) {
            y = dat_vec[i].double_val();
            if ( y > mx) mx = y;
         }
      }
      if (k < ne) {
         retval.missing = 0;
         retval.double_val(mx);
      }
   }
   return retval;
}

DataNode Field::min(void) const
{
   DataNode retval;
   retval.missing = 1;

   if (col_struct.type() != 'S') {
      unsigned i, k=0;
      while (dat_vec[k].missing) k++;
      double y,mx = dat_vec[k].double_val();
      for (i=k; i<ne; i++) {
         if (!dat_vec[i].missing) {
            y = dat_vec[i].double_val();
            if ( y < mx) mx = y;
         }
      }
      if (k < ne) {
         retval.missing = 0;
         retval.double_val(mx);
      }
   }
   return retval;
}

DataNode Field::sum(void) const
{
   DataNode retval;
   retval.missing = 1;

   if (col_struct.type() != 'S') {
      unsigned i,k;
      double x;
      for (x=0.0, k=0,i=0; i<ne; i++) {
         if (!dat_vec[i].missing) {
            k++;
            x += dat_vec[i].double_val();
         }
      }
      if (k) {
         retval.missing = 0;
         retval.double_val(x);
      }
   }
   return retval;
}

DataNode Field::sumsq(void) const
{
   DataNode retval;
   retval.missing = 1;
   if (col_struct.type() !='S') {
      unsigned i,k;
      double y,x;
      for (x=0.0, k=0,i=0; i<ne; i++) {
         if (!dat_vec[i].missing) {
            k++;
            y = dat_vec[i].double_val();
            x += y*y;
         }
      }
      if (k) {
         retval.missing = 0;
         retval.double_val(x);
      }
   }
   return retval;
}

DataNode Field::product(void) const
{
   DataNode retval;
   retval.missing = 1;
   if (col_struct.type() != 'S') {
      unsigned i,k=0;
      double x;
      for (x=1.0,k=0,i=0; i<ne; i++) {
         if (!dat_vec[i].missing) {
            k++;
            x *= dat_vec[i].double_val();
         }
      }
      if (k) {
         retval.missing = 0;
         retval.double_val(x);
      }
   }
   return retval;
}

DataNode Field::mean(const int flag)
{
   if (flag == 0 || col_struct.type()=='S') return mean_value;

   unsigned i,k;
   double x;
   for (x=0.0,k=0,i=0; i<ne; i++) {
      if (!dat_vec[i].missing) {
         k++;
         x += dat_vec[i].double_val();
      }
   }
   if (k >=1 ) {
      x /= k;
      mean_value.missing = 0;
      mean_value.double_val(x);
   }
   return mean_value;
}


DataNode Field::covariance(const Field *B) const
{
   if (B) throw exception("Field::covariance(B): not yet implemented");
   DataNode retval;
   retval.missing = 1;

   if (col_struct.type() != 'S') {
      unsigned i,k;
      for (k=0,i=0; i<ne; i++) if (!dat_vec[i].missing) k++;
      if (k >= 2) {
         Vector<double> tmp(k);
         for (k=0,i=0; i<ne; i++) {
            if (!dat_vec[i].missing) tmp[k++] = dat_vec[i].double_val();
         }
         retval.double_val(tmp.variance());
         retval.missing = 0;
      }
   }
   return retval;
}

void Field::value_for_missing(const double vm)
{
   if (col_struct.type() == 'F') {
      for (unsigned i=0; i<ne; i++) {
          if (dat_vec[i].missing) dat_vec[i].double_val(vm);
      }
   }
   else {
       warning("Field::value_for_missing(): won't work for string column");
   }
}

void Field::set_missing(const unsigned k)
{
   if (k >= ne) throw exception("Field::set_missing(): out of range");
   if (!dat_vec[k].missing) col_struct.count_miss(1);
   dat_vec[k].missing = 1;
}

void Field::pretend_missing(const unsigned k)
{
   if (k >= ne) throw exception("Field::pretend_missing(): out of range");
   dat_vec[k].missing++;
   col_struct.count_miss(1);
}

void Field::recover_missing(const unsigned k)
{
   if (k >= ne) throw exception("Field::recover_missing(): out of range");
   if (dat_vec[k].missing <= 0) {
      warning("recover::recover_missing(): can't recover from missing");
   }
   else {
      dat_vec[k].missing--;
      col_struct.count_miss(-1);
   }
}

int compare_DataNode(const void* A,const void* B)
{
   DataNode *x = (DataNode *)A;
   DataNode *y = (DataNode *)B;
   if (x->missing && y->missing) {
       return 0;
   }
   else if (x->missing && !y->missing) {
       return 1;
   }
   else if (!x->missing && y->missing) {
      return -1;
   }
   else {
      if (x->double_val() < y->double_val()) {
         return -1;
      }
      else if (x->double_val() == y->double_val()) {
         return 0;
      }
      else {
         return 1;
      }
   }
}

Field& Field::sort(void)
{
   qsort((DataNode*)dat_vec,ne,sizeof(DataNode),compare_DataNode);
   return *this;
}

void Field::out_to_stream(std::ostream& stream,const int ic) const
{
   stream << "\n";
   unsigned i,k;
   stream.precision(6);
   unsigned W = 6+6;
   char ch;

   if (col_struct.type()=='S') {
      for (k=23,i=0; i<ne; i++) {
         if (ic && i>=k) {
            k += 23;
            stream << "  more ... [q for quit] ";
            std::cin.get(ch);
            std::cin.seekg(0L,std::ios::beg);
            if (ch == 'q') break;
         }
         if (dat_vec[i].missing) {
            stream << " " << std::setw(W) << ".";
         }
         else {
            stream << " " << std::setw(W)
                          << (char *)hashtable->find(dat_vec[i].unsigned_val());
         }
         stream << "\n";
      }
   }
   else {
      for (k=23,i=0; i<ne; i++) {
         if (ic && i>=k) {
            k += 23;
            stream << "  more ... [q for quit] ";
            std::cin.get(ch);
            std::cin.seekg(0L,std::ios::beg);
            if (ch == 'q') break;
         }
         if (dat_vec[i].missing) {
            stream << " " << std::setw(W) << ".";
         }
         else {
            stream << " " << std::setw(W) << dat_vec[i].double_val();
         }
         stream << "\n";
      }
   }
}

void Field::save(const std::string &fname,const int io_mode) const
{
   std::ofstream ofs;
   ofs.open(fname.c_str(),(OpenModeType)io_mode);
   if (!ofs) throw exception("Field::save(): cannot open file");
   this->out_to_stream(ofs,0);
   ofs.close();
}

void Field::display(const std::string &msg,const int ic) const
{
   if (msg != "") std::cout <<"\n     " << msg << "\n";
   this->out_to_stream(std::cout,ic);
}

std::ostream& operator<<(std::ostream& stream, const Field& V)
{
   V.out_to_stream(stream,1);
   return stream;
}
}

