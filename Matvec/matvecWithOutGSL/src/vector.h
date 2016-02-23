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

#ifndef MATVEC_VECTOR_H
#define MATVEC_VECTOR_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include "session.h"
#include "exception.h"
#include "util.h"
/**
 * Namespace matvec.
 */
namespace matvec {
/**
 *  A class template for one-dimensional array.
 *  Vector is also a container. For example:
 *
 *  matvec::Vector<int> a(5);
 *
 * @see Matrix
 */ 
template<class T> class Vector {

public:
   typedef         int size_type;
   typedef         T   value_type;
   typedef         T   element_type;
   typedef         T*  pointer;
   typedef         T*  iterator;
   typedef         T&  reference;
   typedef const   T*  const_iterator;
   typedef const   T*  const_pointer;
   typedef const   T&  const_reference;
   typedef std::reverse_iterator<iterator>reverse_iterator;
   typedef std::reverse_iterator<const_iterator>const_reverse_iterator;

   Vector(void)                         { initialize(0,0); }        ///< Constructs an empty vector of type T
   explicit Vector(const size_type n)   { initialize(n,0); }        ///< Constructs a vector of size n with each element being T()
   Vector(const size_type n,const T *a) { initialize(n,a); }        ///< Constructs a vector of size n initialized by *a.
   Vector(const Vector &a)              { initialize(a.ne,a.ve); }  ///< Copy constructor.

   iterator               begin(void)        { return ve; }         ///< returns a random access iterator for the first element.
   iterator               end()              { if (ve) return ve + ne; else return 0; }    ///< returns a random access iterator for the position after the last element.
  /*const*/ iterator         begin(void)  const { return ve; }
  /* const*/ iterator         end(void)    const { if (ve) return ve + ne; else return 0; }
   reverse_iterator       rbegin(void)       { return reverse_iterator(end());}         ///< returns a reverse iterator for the first element of a reverse iteration.
   reverse_iterator       rend(void)         { return reverse_iterator(begin()); }      ///< returns a reverse iterator for the position after the last element of a reverse iteration.
   const reverse_iterator rbegin(void) const { return const_reverse_iterator(end()); }
   const reverse_iterator rend(void)   const { return const_reverse_iterator(begin());}

   virtual ~Vector(void){clear();}     ///< Destructs a vector.

   void             copy(const Vector &a);
   Vector&          assign(const Vector &a);
   Vector&          assign(const T &x);

   Vector&          resize(const size_type n, const T &val = T());
   Vector&          resize(const size_type n, const T *a);
   Vector&          resize(const Vector& a){resize(a.ne); return *this;}
   Vector&          reserve(const size_type n);

   Vector&          operator =  (const Vector &a) { return assign(a); }
   Vector&          operator =  (const T &x)      { return assign(x); }

   Vector           operator +  (const Vector &a) const;
   Vector           operator -  (const Vector &a) const;
   Vector           operator *  (const Vector &a) const;
   Vector           operator /  (const Vector &a) const;

   Vector           operator +  (const T &x) const;
   Vector           operator -  (const T &x) const;
   Vector           operator *  (const T &x) const;
   Vector           operator /  (const T &x) const;

   Vector&          operator += (const Vector &a);
   Vector&          operator -= (const Vector &a);
   Vector&          operator *= (const Vector &a);
   Vector&          operator /= (const Vector &a);
   Vector&          operator += (const T &x);
   Vector&          operator -= (const T &x);
   Vector&          operator *= (const T &x);
   Vector&          operator /= (const T &x);

   T&               operator [] (const size_type i) {return ve[i];}
   const T&         operator [] (const size_type i) const {return ve[i];}
   T&               operator () (const size_type i);
   Vector           operator - (void) const;
   Vector&          operator + (void) {return *this;}

   Vector<bool>     operator !  (void) const;

   Vector<bool>     operator == (const Vector &a) const;
   Vector<bool>     operator <  (const Vector &a) const;
   Vector<bool>     operator >  (const Vector &a) const;
   Vector<bool>     operator != (const Vector &a) const;
   Vector<bool>     operator <= (const Vector &a) const;
   Vector<bool>     operator >= (const Vector &a) const;

   Vector<bool>     operator == (const T& x) const;
   Vector<bool>     operator <  (const T& x) const;
   Vector<bool>     operator >  (const T& x) const;
   Vector<bool>     operator != (const T& x) const;
   Vector<bool>     operator <= (const T& x) const;
   Vector<bool>     operator >= (const T& x) const;

   bool             all(void) const;
   bool             any(void) const;
   bool             empty(void) const { return ne == 0;}
   Vector           apply(T (*f)(const T&)) const;
   Vector           apply(T (*f)(const T&,const T&),const Vector &b) const;
   Vector           apply(T (*f)(const T&,const T&),const T &b) const;
   Vector           apply(T (*f)(T)) const;
   Vector           apply(T (*f)(T,T),const Vector &b) const;
   Vector           apply(T (*f)(T,T),const T &b) const;
   Vector&          append(const T &x);
   Vector&          append(const Vector &a);
   Vector&          sort(void) {std::sort(begin(),end()); return *this;}
   T                sum(T init = T()) const {return std::accumulate(begin(),end(),init);}
   T                inner_product(const T  *b, T init = T()) const;
   T                inner_product(const Vector &b,T init = T()) const {return inner_product(b.begin(),init);}
   T                inner_product(T init = T()) const { return inner_product(this->begin(),init); }
   T                sumsq(T init = T()) const { return inner_product(this->begin(),init); }
   T                variance(T init = T()) const;
   Vector           subvec(const size_type idx = 0,const unsigned len = 0) const;
   Vector           subvec(const Vector<bool> &v) const;
   Vector           find(const Vector<bool> &v) const;
   int              find(const T &x) const;
   T                max(void) const { return *(std::max_element(begin(),end())); }
   T                min(void) const { return *(std::min_element(begin(),end())); }
   std::ostream&         print(std::ostream &os = std::cout) const;
   T&               at(const size_type i);
   T&               front(void) { return *begin(); }
   const T&         front(void) const { return *begin(); }
   T&               back() { return *(end() - 1); }
   const T&         back() const { return *(end() - 1); }

   int              size(void) const {return ne;}
   void             clear(void) { ///< Releases memory.
                      if (ne > 0) {
			delete [] ve; 
			ne = 0; 
			ve = 0;
		      } 
                    } 
   void             input(const std::string &fname,const size_type n);
  bool             save(const std::string &fname,const int io_mode = std::ios::out) const;//SDK |std::ios::noreplace) const;
   void initialize(size_type n,const T* a);
protected:
   T *ve;           ///< it store the pointer to the memory allocated for @c ne element.
   int ne;          ///< it stores the number of elements.
};

/**
 * A private member used by constructors.
 */
template<class T> inline void Vector<T>::initialize(const size_type n,const T *a)
{
   if (n == 0) {
      ne = 0; ve = 0;
   } else {
      ne = n;
      ve = new T [ne];
      check_ptr(ve);
      T *bot = ve, *top = ve + ne;
      if (a) {
         while (bot < top) *bot++ = *a++;
      } else {
         while (bot < top) *bot++ = T();
      }
   }
}

/**
 * It returns a reference to element i with range-cheching.
 * Note that element index starting from 0.
 */
template<class T> inline T& Vector<T>::at(const size_type i)
{
   if (i < 0 || i >= ne) throw exception("Vector<T>::at(i): out of range");
   return ve[i];
}

/**
 * It returns a reference to element i with range-cheching.
 * Note that element index starting from 1.
 */
template<class T> inline T& Vector<T>::operator () (const size_type i)
{
   if (i <= 0 || i > ne) throw exception("Vector<T>::operator(): out of range");
   return ve[i-1];
}

/**
 * Input data from an external text file.
 * @param fname Filename
 * @param n number of elements you want to read
 * @see save()
 */
template<class T> inline void Vector<T>::input(const std::string &fname,const size_type n)
{
   reserve(n);
   std::ifstream infile(fname.c_str(),std::ios::in);//SDK | std::ios::nocreate);
   if (!infile) throw exception("Vector<T>::input(): cannot open file");
   int i = 0;
   while (!infile.eof()) {
      if (i >= n) break;
      infile >> ve[i++];
   }
   infile.close();
   return;
}

/**
 * Save the content into an text file.
 * @param fname Filename
 * @param io_mode determines how the file is opened. It must be one (or more by ORing) of these values:
 * <ul>
 *   <li> ios::in        Opens for reading;
 *   <li> ios::out       Opens for writing;
 *   <li> ios::ate       Seeks to EOS after file is create;
 *   <li> ios::app       All writes added to end of file;
 *   <li> ios::truct     If file already exists, truncates;
 *   <li> ios::nocreate  Won't open if file does not exit;
 *   <li> ios::noreplace Won't  open if file does exit;
 *   <li> ios::binary    Opens file in binary mode (default text)
 * </ul>
 * For instance ios::app|ios::nocreate  append all writes to an existing file, it won't open if file does not exit.
 *
 * @see input()
 */
template<class T> inline bool Vector<T>::save(const std::string &fname,const int io_mode) const
{
   std::ofstream ofs;
   ofs.open(fname.c_str(),(OpenModeType)io_mode);
   if (!ofs) throw exception("Vector<T>::save(): cannot open file");
   print(ofs);
   ofs.close();
   return true;
}

template<class T> inline void Vector<T>::copy(const Vector<T>& a)
{
   if (this == &a) return;
   reserve(a.ne);
   // memcpy did not work correctly for the BG class (RLF)
   //memcpy(ve,a.ve,sizeof(T)*ne);
   for (unsigned i=0;i<ne;i++){
     ve[i] = a.ve[i];
   }
   return;
}

template<class T> inline Vector<T>& Vector<T>::assign(const Vector<T>& a)
{
   copy(a);
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::assign(const T& a)
{
   if (ve) {
      resize(ne,a);
   }
   else {
      resize(1,a);
   }
   return *this;
}

/**
 * It resizes vector without initialization.
 * @see resize
 */
template<class T> inline Vector<T>& Vector<T>::reserve(const size_type n)
{
   if (ne != n) {
      clear();
      ne = n;
      if(ne>0){
	ve = new T [ne];
      }
      else {
	ve = 0;
      }
      check_ptr(ve);
   }
   return *this;
}

/**
 * It resizes vector with each element being val.
 * @see resize
 */
template<class T> inline Vector<T>& Vector<T>::resize(const size_type n,const T &val)
{
   if (n == 0) {
      clear();
   } else {
      if (ne != n) {
         clear();
         ne = n;
         ve = new T [ne];
         check_ptr(ve);
      }
      iterator it=begin(), endit = end();
      while  (it != endit) *it++ = val;
   }
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::resize(const size_type n, const T *a)
{
   if (n == 0) {
      clear();
   } else {
      if (ne != n) {
         clear();
         ne = n;
         ve = new T [ne];
         assert( ve != 0 );
      }
      iterator it = begin(), endit = end();
      while (it != endit) *it++ = *a++;
   }
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator += (const Vector<T> &a)
{
   assert( ne == a.ne );
   iterator it = begin(), endit = end(), a_it = a.begin();
   while (it != endit) *it++ += *a_it++;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator -= (const Vector<T> &a)
{
   assert( ne == a.ne );
   iterator it = begin(), endit = end(), a_it = a.begin();
   while (it != endit) *it++ -= *a_it++;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator *= (const Vector<T> &a)
{
   assert( ne == a.ne );
   iterator it = begin(), endit = end(), a_it = a.begin();
   while (it != endit) *it++ *= *a_it++;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator /= (const Vector<T> &a)
{
   assert( ne == a.ne );
   iterator it = begin(), endit = end(), a_it = a.begin();
   while (it != endit) *it++ /= *a_it++;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator += (const T &s)
{
   iterator it = begin(), endit = end();
   while (it != endit) *it++ += s;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator -= (const T &s)
{
   iterator it = begin(), endit = end();
   while (it != endit) *it++ -= s;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator *= (const T &s)
{
   iterator it = begin(), endit = end();
   while (it != endit) *it++ *= s;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::operator /= (const T &s)
{
   iterator it = begin(), endit = end();
   while (it != endit) *it++ /= s;
   return *this;
}

template<class T> inline Vector<bool> Vector<T>::operator == (const Vector<T> &a) const
{
   assert(ne == a.size());
   Vector<bool> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin();
   bool *t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ == *a_it++;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator == (const T &x) const
{
   Vector<bool> temp(ne);
   iterator it = begin(), endit = end();
   bool *t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ == x;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator < (const Vector<T> &a) const
{
   assert (ne == a.ne );
   Vector<bool> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin();
   bool *t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ < *a_it++;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator < (const T &x) const
{
   Vector<bool> temp(ne);
   iterator it = begin(), endit = end();
   bool *t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ < x;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator > (const Vector<T> &a) const
{
   assert (ne != a.ne );
   Vector<bool> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin();
   bool *t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ > *a_it++;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator > (const T &x) const
{
   Vector<bool> temp(ne);
   iterator it = begin(), endit = end();
   bool *t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ > x;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator != (const Vector<T> &a) const
{ return !(*this == a); }

template<class T> inline Vector<bool> Vector<T>::operator != (const T &x) const
{ return !(*this == x); }

template<class T> inline Vector<bool> Vector<T>::operator <= (const Vector<T> &a) const
{ return !(*this > a); }

template<class T> inline Vector<bool> Vector<T>::operator <= (const T &x) const
{ return !(*this > x); }


template<class T> inline Vector<bool> Vector<T>::operator >= (const Vector<T> &a) const
{ return !(*this < a); }

template<class T> inline Vector<bool> Vector<T>::operator >= (const T &x) const
{ return !(*this < x); }

template<class T> inline Vector<T> Vector<T>::operator + (const Vector<T> &a) const
{
   assert ( ne == a.ne );
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ + *a_it++;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator - (const Vector<T> &a) const
{
   assert ( ne == a.ne );
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ - *a_it++;
   return temp;
}


template<class T> inline Vector<T> Vector<T>::operator * (const Vector<T> &a) const
{
   assert ( ne == a.ne );
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ * *a_it++;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator / (const Vector<T> &a) const
{
   assert ( ne == a.ne );
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(), a_it = a.begin(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ / *a_it++;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator + (const T &a) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ + a;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator - (const T &a) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ - a;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator * (const T &a) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ * a;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator / (const T &a) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = *it++ / a;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::operator - (void) const  //unary minus
{
   if (ne == 0) return *this;
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = - *it++;
   return temp;
}

template<class T> inline Vector<bool> Vector<T>::operator ! (void) const  //unary negative
{
   Vector<bool> temp(ne);
   if (ne == 0) return temp;
   for (int i=0; i<ne; ++i) temp[i] = !ve[i];
   return temp;
}

template<class T> inline Vector<T> Vector<T>::apply(T (*f)(const T&)) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = (*f)(*it++);
   return temp;
}

template<class T> inline Vector<T> Vector<T>::apply(T (*f)(const T&,const T&),const Vector &a) const
{
   assert(ne == a.size());
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin(),a_it = a.begin();
   while (it != endit) *t_it++ = (*f)(*it++,*a_it++);
   return temp;
}

template<class T> inline Vector<T> Vector<T>::apply(T (*f)(const T&,const T&),const T &a) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin(),a_it = a.begin();
   while (it != endit) *t_it++ = (*f)(*it++,a);
   return temp;
}

template<class T> inline Vector<T> Vector<T>::apply(T (*f)(T)) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = (*f)(*it++);
   return temp;
}

template<class T> inline Vector<T> Vector<T>::apply(T (*f)(T,T),const Vector &a) const
{
   assert(ne == a.size());
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin(),a_it = a.begin();
   while (it != endit) *t_it++ = (*f)(*it++,*a_it++);
   return temp;
}

template<class T> inline Vector<T> Vector<T>::apply(T (*f)(T,T),const T &a) const
{
   Vector<T> temp(ne);
   iterator it = begin(), endit = end(),t_it = temp.begin();
   while (it != endit) *t_it++ = (*f)(*it++,a);
   return temp;
}

template<class T> inline Vector<T>& Vector<T>::append(const T &x)
{
   Vector<T> temp(ne,ve);
	 int n(ne);
   reserve(ne+1);
   for (int i=0; i<n; ++i) ve[i] = temp[i];
   ve[ne] = x;
   return *this;
}

template<class T> inline Vector<T>& Vector<T>::append(const Vector &a)
{
   Vector<T> temp(ne,ve);
   int i,n(ne);
   reserve(n + a.ne);
   for (i=0; i<n; ++i) ve[i] = temp[i];
   for (i=0; i<a.ne; ++i) ve[n+i] = a[i];
   return *this;
}

template<class T> inline bool Vector<T>::all(void) const
{
   if (!ve) return false;
   iterator it = begin(), endit = end();
   while (it != endit) if (*it++ == false) return false;
   return true;
}

template<class T> inline bool Vector<T>::any(void) const
{
   if (!ve) return false;
   iterator it = begin(), endit = end();
   while (it != endit) if (*it++ == true) return true;
   return false;
}

template<class T> inline T Vector<T>::inner_product(const T  *b, T init) const
{
   iterator it = begin(), endit = end();
   while (it != endit) init = init + (*it++ * *b++);
   return init;
}

template<class T> inline T Vector<T>::variance(T init) const
{
   assert( ne > 1);
   T x(sum(init));
   init = (inner_product(init) - x*x/ne)/(ne - 1);
   return init;
}

template<class T> inline Vector<T> Vector<T>::subvec(const size_type idx,const unsigned len) const
{
   assert (idx >= 0 && idx < ne && len <= (ne - idx));
   int newne = len;
   if (newne == 0) {newne = ne - idx;}
   Vector<T> temp(newne);
   T *bot = &ve[idx], *top = bot + newne, *t_pt = temp.begin();
   while (bot < top) *t_pt++ = *bot++;
   return temp;
}

template<class T> inline Vector<T> Vector<T>::subvec(const Vector<bool> &v) const
{
   int n = 0;
   int m = v.size();
   bool *ve_pt = v.ve;
   int i;
   for (i=0; i<m; i++) if (*ve_pt++) n++;

   Vector<T> temp(n);
   T *t_pt = temp.begin();
   ve_pt = v.ve;
   for (i=0; i<m; i++) if (*ve_pt++) *t_pt++ = ve[i];
   return temp;
}

template<class T> inline int Vector<T>::find(const T& x) const
{
   int i = -1;
   for (i=0; i<ne; ++i) {
      if (ve[i] == x) break;
   }
   if (i == ne) i = -1;
   return i;
}

template<class T> inline Vector<T> Vector<T>::find(const Vector<bool> &a) const
{
   assert (ne == a.size() );
   int n = 0;
   bool *bot = a.begin(), *top = a.end();
   while (bot < top) if (*bot++) n++;
   Vector<T> temp(n);
   n = 0;
   bot = a.begin();
   T *pt = ve;
   while (bot < top) {
      if (*bot++) temp[n++] = *pt;
      pt++;
   }
   return temp;
}

template<class T> inline std::ostream & Vector<T>::print(std::ostream &os) const
{
  for (int i=0; i<ne; ++i) {
    os.width(8);
    os.precision(SESSION.output_precision); 
    os <<ve[i] << "\n";
  }
   return os;
}

/////////// the following are  functions/operator template ///////////////

/**
 * @relates Vector
 * returns a vector of sin(a[i]).
 */
template<class T> Vector<T> sin  (const Vector<T> &a) {return a.apply(std::sin);}

/**
 * @relates Vector
 * returns a vector of asin(a[i]).
 */
template<class T> Vector<T> asin (const Vector<T> &a) {return a.apply(std::asin);}

/**
 * @relates Vector
 * returns a vector of cos(a[i]).
 */

template<class T> Vector<T> cos  (const Vector<T> &a) {return a.apply(std::cos);}

/**
 * @relates Vector
 * returns a vector of acos(a[i]).
 */
template<class T> Vector<T> acos (const Vector<T> &a) {return a.apply(std::acos);}

/**
 * @relates Vector
 * returns a vector of tan(a[i]).
 */
template<class T> Vector<T> tan  (const Vector<T> &a) {return a.apply(std::tan);}

/**
 * @relates Vector
 * returns a vector of atan(a[i]).
 */
template<class T> Vector<T> atan (const Vector<T> &a) {return a.apply(std::atan);}

/**
 * @relates Vector
 * returns a vector of ceil(a[i]).
 */
template<class T> Vector<T> ceil (const Vector<T> &a) {return a.apply(std::ceil);}

/**
 * @relates Vector
 * returns a vector of floor(a[i]).
 */
template<class T> Vector<T> floor(const Vector<T> &a) {return a.apply(std::floor);}

/**
 * @relates Vector
 * returns a vector of log(a[i]).
 */
template<class T> Vector<T> log  (const Vector<T> &a) {return a.apply(std::log);}

/**
 * @relates Vector
 * returns a vector of log10(a[i]).
 */
template<class T> Vector<T> log10(const Vector<T> &a) {return a.apply(std::log10);}

/**
 * @relates Vector
 * returns a vector of exp(a[i]).
 */
template<class T> Vector<T> exp  (const Vector<T> &a) {return a.apply(std::exp);}

/**
 * @relates Vector
 * returns a vector of sqrt(a[i]).
 */
template<class T> Vector<T> sqrt (const Vector<T> &a) {return a.apply(std::sqrt);}

/**
 * @relates Vector
 * returns a vector of abs(a[i]).
 */
template<class T> Vector<T> abs  (const Vector<T> &a) {return a.apply(std::abs);}

/**
 * @relates Vector
 * returns a vector of erf(a[i]).
 */
template<class T> Vector<T> erf  (const Vector<T> &a) {return a.apply(erf);}

/**
 * @relates Vector
 * returns a vector of erfc(a[i]).
 */
template<class T> Vector<T> erfc (const Vector<T> &a) {return a.apply(erfc);}


/**
 * @relates Vector
 * If each element is true return true, otherwise return false.
 */
template<class T> bool all (const Vector<T> &a) {return a.all();}


/**
 * @relates Vector
 * If each element is false return false, otherwise return true.
 */
template<class T> bool any (const Vector<T> &a) {return a.any();}


/**
 * @relates Vector
 * returns a + b.
 */
template<class T> Vector<T> operator + (const T &a,const Vector<T> &b) {return b+a;}


/**
 * @relates Vector
 * returns a - b.
 */
template<class T> Vector<T> operator - (const T &a,const Vector<T> &b) {return -b+a;}


/**
 * @relates Vector
 * returns a * b..
 */
template<class T> Vector<T> operator * (const T &a,const Vector<T> &b) {return b*a;}


/**
 * @relates Vector
 * returns a  / b
 */
template<class T> Vector<T> operator / (const T &a,const Vector<T> &b)
{
   Vector<T> temp(b.size());
   for (int i=0; i<b.size(); ++i) temp[i] = a / b[i];
   return temp;
}

/**
 * @relates Vector
 * It print Vector a to an output stream.
 */
template<class T> std::ostream&  operator << (std::ostream &os, const Vector<T> &a) {return a.print(os);}
} ///////////// end of namespace matvec

#endif
