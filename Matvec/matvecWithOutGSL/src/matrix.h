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

#ifndef MATVEC_MATRIX_H
#define MATVEC_MATRIX_H

#include <string>
#include <fstream>
#include <cassert>
#include "session.h"
#include "vector.h"


namespace matvec {
/*!
   \class Matrix  vector.h
   \brief  A vector is a on-dimensional array with double precision

  \sa Matrix
*/


template<class T> class Matrix {

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

   enum orientation {ROW,COLUMN};

   Matrix(void) { 
     //std::cout <<"constructor 1\n";
     initialize(0,0,0);  //Constructor 1
   } 
   Matrix(const size_type m,const size_type n) { 
     //std::cout <<"constructor 2\n";
     initialize(m,n,0);   //Constructor 2
   }
   Matrix(const size_type m,const size_type n,const T** a) { 
     //std::cout <<"constructor 3\n";
     initialize(m,n,a);   //Constructor 3
   }
   Matrix(const Matrix& a){ 
     //std::cout <<"constructor 4\n";
     initialize(a.nrow,a.ncol,(const T**)a.me);   //Constructor 4
   }

   T**       begin(void)       { return me; }
   const T** begin(void) const { return (const T**)me; }

   virtual ~Matrix(void){clear();}                                  // Destructor

   void              copy(const Matrix &a);
   Matrix&           assign(const Matrix &a) {copy(a); return *this;}
   Matrix&           assign(const T &x);

   Matrix&           resize(const size_type n,const size_type m, const T &val = T());
   Matrix&           resize(const size_type n,const size_type m, const T *a);
   Matrix&           resize(const Matrix& a){resize(a.ne); return *this;}
   Matrix&           reserve(const size_type n,const size_type m);

   Matrix&           operator =  (const Matrix &a) { return assign(a); }
   Matrix&           operator =  (const T &x)      { return assign(x); }

   Matrix            operator +  (const Matrix &a) const;
   Matrix            operator -  (const Matrix &a) const;
   Matrix            operator *  (const Matrix &a) const;
   Matrix            operator /  (const Matrix &a) const;

   Vector<T>         operator *  (const Vector<T> &a) const;

   Matrix            operator +  (const T& x) const;
   Matrix            operator -  (const T& x) const;
   Matrix            operator *  (const T& x) const;
   Matrix            operator /  (const T& x) const;

   Matrix&           operator += (const Matrix& a);
   Matrix&           operator -= (const Matrix& a);
   Matrix&           operator *= (const Matrix& a);
   Matrix&           operator /= (const Matrix& a);
   Matrix&           operator += (const T& x);
   Matrix&           operator -= (const T& x);
   Matrix&           operator *= (const T& x);
   Matrix&           operator /= (const T& x);

   T*                operator [] (const size_type i) { return me[i]; }
   const T*          operator [] (const size_type i) const { return me[i]; }
   T&                operator () (const size_type i,const size_type j);
   Matrix            operator -  (void) const;               //unary minus
   Matrix&           operator +  (void) { return *this; }      //unary plus

   Matrix<bool>      operator !  (void) const;

   Matrix<bool>      operator == (const Matrix &a) const;
   Matrix<bool>      operator <  (const Matrix &a) const;
   Matrix<bool>      operator >  (const Matrix &a) const;
   Matrix<bool>      operator != (const Matrix &a) const;
   Matrix<bool>      operator <= (const Matrix &a) const;
   Matrix<bool>      operator >= (const Matrix &a) const;

   Matrix<bool>      operator == (const T &x) const;
   Matrix<bool>      operator <  (const T &x) const;
   Matrix<bool>      operator >  (const T &x) const;
   Matrix<bool>      operator != (const T &x) const;
   Matrix<bool>      operator <= (const T &x) const;
   Matrix<bool>      operator >= (const T &x) const;

   bool              all(void) const;
   bool              any(void) const;
   bool              symmetric(void) const;
   bool              empty(void) const { return (nrow == 0 || ncol == 0); }
  bool              save(const std::string &fname,const int io_mode = std::ios::out) const;//SDK|std::ios::noreplace) const;

   Matrix            apply(T (*f)(T)) const;
   Matrix            apply(T (*f)(T,T),const Matrix &b) const;
   Matrix            apply(T (*f)(T,T),const T &b) const;
   Matrix            submat(const size_type idx1,const size_type idx2, const int len1 = 0, const int len2 = 0) const;
   Matrix            submat(const Matrix<bool>& a) const;
   Matrix            transpose(void) const;
   Matrix            kron(const Matrix &b) const;
   Matrix            reshape(const size_type m,const size_type n,orientation orient = COLUMN) const;
   Matrix            diag(void) const;
   Vector<T>         diag(const size_type k) const;
   Vector<T>         vec(orientation orient = COLUMN) const;
   Vector<T>         sum(orientation orient = COLUMN) const;
   Vector<T>         sumsq(orientation orient = COLUMN) const;
   Vector<T>         max(orientation orient = COLUMN) const;
   Vector<T>         min(orientation orient = COLUMN) const;
   virtual std::ostream&  print(std::ostream &os = std::cout) const;
   T                 trace(T init = T()) const;
   T&                at(const size_type i,const size_type j);

   Matrix&           sortby(orientation orient,const size_type k);
   int               size(void)     const { return nrow*ncol; }
   int               num_rows(void) const { return nrow; }
   int               num_cols(void) const { return ncol; }
   void              clear(void);
   void              input(const std::string &fname,
                           const size_type m,
                           const size_type n);

   //protected:
   T **me;
   int nrow,ncol;

   void initialize(const size_type m,const size_type n,const T **a);
};

template<class T> inline void Matrix<T>::initialize(const size_type m,const size_type n,const T **a)
{
   if (m == 0 || n == 0) {
      nrow = 0; ncol = 0; me = 0;
   } else {
      nrow = m; ncol = n;
      me = new T* [nrow];
      assert( me != 0 );
      for (int j,i=0; i<nrow; ++i) {
         me[i] = new T [ncol];
         assert( me[i] != 0 );
         if (a) {
            for (j=0; j<ncol; ++j) me[i][j] = a[i][j];
         } else {
            for (j=0; j<ncol; ++j) me[i][j] = T();
         }
      }
   }
   return;
}

template<class T> inline T& Matrix<T>::at(const size_type i,const size_type j)
{
   if (i < 0 || i >= nrow || j < 0 || j >= ncol) throw exception("Matrix<T>::at(): out of range");
   return me[i][j];
}

template<class T> inline void Matrix<T>::clear(void)
{
   if (me) {
      for (int i=0; i<nrow; ++i) delete [] me[i];
      delete [] me;
      me = 0;
   }
   nrow = 0;
   ncol = 0;
}

template<class T> inline void Matrix<T>::input(const std::string &fname,
                                               const size_type m,
                                               const size_type n)
{
   resize(m,n);
   std::ifstream infile(fname.c_str(),std::ios::in);//SDK | std::ios::nocreate);
   if (!infile) {
      std::cerr << fname << ": cannot open\n";
   } else {
      for (int j,i=0; i<m; ++i) for (j=0; j<n; ++j) infile >> me[i][j];
      infile.close();
   }
   return;
}

template<class T> inline void Matrix<T>::copy(const Matrix<T> & a)
{
   if (this == &a) return;
   reserve(a.nrow,a.ncol);
/*    for (int i=0; i<nrow; ++i) memcpy(me[i],a.me[i],sizeof(T)*ncol);  */
/*    memcopy does not work correctly for complex objects RLF */
   for (int i=0; i<nrow; ++i) {
     for (int j=0; j<ncol; j++){
       me[i][j] = a.me[i][j];
     }
   }
   return;
}

template<class T> inline Matrix<T>& Matrix<T>::assign(const T& a)
{
   if (me) {
      resize(nrow,ncol,a);
   } else {
     resize(1,1,a);
   }
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::reserve(const size_type m,const size_type n)
{
   if (m != nrow || n != ncol) {
      clear();
      nrow = m;
      ncol = n;
      if(nrow>0){
	me = new T* [nrow];
      }
      else {
	me = 0;
	return *this;
      }
      assert( me != 0 );
      for (int i=0; i<nrow; ++i) {
	if(ncol>0){
	  me[i] = new T [ncol];
	}
	else {
	  me[i] = 0;
	}
         assert( me[i] != 0 );
      }
   }
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::resize(const size_type m,const size_type n,const T& val)
{
   if (m == 0 || n == 0) {
      clear();
   } else {
      int i,j;
      if (m != nrow || n != ncol) {
         clear();
         nrow = m;
         ncol = n;
         me = new T* [nrow];
         assert( me != 0 );
         for (i=0; i<nrow; ++i) {
            me[i] = new T [ncol];
            assert( me[i] != 0 );
         }
      }
      for (i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] = val;
   }
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::resize(const size_type m,const size_type n,const T *a)
{
   if (m == 0 || n == 0) {
      clear();
   } else {
      int i,j;
      if (m != nrow || n != ncol) {
         clear();
         nrow = m;
         ncol = n;
         me = new T* [nrow];
         assert( me != 0 );
         for (i=0; i<nrow; ++i) {
            me[i] = new T [ncol];
            assert( me[i] != 0 );
         }
      }
      for (i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] = *a++;
   }
   return *this;
}

template<class T> inline T& Matrix<T>::operator () (const size_type i,size_type j)
{
   if (i < 1 || i > nrow || j < 1 || j > ncol) throw exception("Matrix<T>::operator(): out of range");
   return me[i-1][j-1];
}

template<class T> inline Matrix<T>& Matrix<T>::operator += (const Matrix<T> &a)
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator+=: size not conformable");
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] += a[i][j];
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator -= (const Matrix<T> &a)
{
   if ( nrow != a.nrow || ncol != a.ncol ) throw exception("Matrix<T>::operator-=: size not conformable");
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] -= a[i][j];
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator *= (const Matrix<T> &a)
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator*=: size not conformable");
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] *= a[i][j];
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator /= (const Matrix<T> &a)
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator/=: size not conformable");
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] /= a[i][j];
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator += (const T &s)
{
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] += s;
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator -= (const T &s)
{
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] -= s;
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator *= (const T &s)
{
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] *= s;
   return *this;
}

template<class T> inline Matrix<T>& Matrix<T>::operator /= (const T &s)
{
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) me[i][j] /= s;
   return *this;
}

template<class T> inline Matrix<bool> Matrix<T>::operator == (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol)  throw exception("Matrix<T>::operator==: size not conformable");
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] == a[i][j]);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator == (const T &x) const
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] == x);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator < (const Matrix<T>& a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator<: size not conformable");
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] < a[i][j]);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator < (const T &x) const
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] < x);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator > (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator>: size not conformable");
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] > a[i][j]);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator > (const T &x) const
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] > x);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator != (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator!=: size not conformable");
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] != a[i][j]);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator != (const T &x) const
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] != x);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator <= (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator<=: size not conformable");
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] <= a[i][j]);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator <= (const T &x) const
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] <= x);
   return temp;
}


template<class T> inline Matrix<bool> Matrix<T>::operator >= (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator>=: size not conformable");
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] >= a[i][j]);
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator >= (const T &x) const
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (me[i][j] >= x);
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator + (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator+: size not conformable");
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] + a.me[i][j];
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator - (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator-: size not conformable");
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] - a.me[i][j];
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator * (const Matrix<T> &a) const
{
   if (ncol != a.nrow) throw exception("Matrix<T>::operator*: size not conformable");
   Matrix<T> temp(nrow,a.ncol);
   T x;
   for (int t,j,i=0; i<nrow; ++i) {
      for (j=0; j<a.ncol; ++j) {
         x = T();
         for (t=0; t<ncol; ++t) x += me[i][t]*a[t][j];
         temp[i][j] = x;
      }
   }
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::operator * (const Vector<T> &a) const
{
   if (ncol != a.size()) throw exception("Matrix<T>::operator*: size not conformable");
   Vector<T> temp(nrow);
   T x;
   for (int j,i=0; i<nrow; ++i) {
      x = T();
      for (j=0; j<ncol; ++j)  x += me[i][j]*a[j];
      temp[i] = x;
   }
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator / (const Matrix<T> &a) const
{
   if (nrow != a.nrow || ncol != a.ncol) throw exception("Matrix<T>::operator/: size not conformable");
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] / a.me[i][j];
   return temp;
}


template<class T> inline Matrix<T> Matrix<T>::operator + (const T &a) const
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] + a;
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator - (const T &a) const
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] - a;
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator * (const T &a) const
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] * a;
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator / (const T &a) const
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = me[i][j] / a;
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::operator - (void) const  //unary minus
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = - me[i][j];
   return temp;
}

template<class T> inline Matrix<bool> Matrix<T>::operator ! (void) const  //unary negative
{
   Matrix<bool> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = ! me[i][j];
   return temp;
}

template<class T> inline bool Matrix<T>::all(void) const
{
   if (!me) return false;
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) if (me[i][j] == false) return false;
   return true;
}

template<class T> inline bool Matrix<T>::any(void) const
{
   if (!me) return false;
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) if (me[i][j] == true) return true;
   return false;
}

template<class T> inline bool Matrix<T>::symmetric(void) const
{
   if (!me) return false;
   for (int j,i=1; i<nrow; ++i) for (j=0; j<i; ++j) if (me[i][j] != me[j][i]) return false;
   return true;
}

template<class T> inline bool Matrix<T>::save(const std::string &fname,const int io_mode) const
{
   std::ofstream ofs(fname.c_str(),(OpenModeType)io_mode);
   if (ofs) {
      print(ofs);
      ofs.close();
   } else {
      std::cerr << fname << ": cannot open\n";
   }
   return true;
}

template<class T> inline T Matrix<T>::trace(T init) const
{
   int n = std::min(nrow,ncol);
   for (int i=0; i<n; ++i) init += me[i][i];
   return init;
}

template<class T> inline Matrix<T> Matrix<T>::apply(T (*f)(T)) const
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (*f)(me[i][j]);
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::apply(T (*f)(T,T),const Matrix &b) const
{
   if (nrow != b.nrow || ncol != b.ncol) throw exception("Matrix<T>::apply(): size not conformable");
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (*f)(me[i][j],b[i][j]);
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::apply(T (*f)(T,T),const T &b) const
{
   Matrix<T> temp(nrow,ncol);
   for (int j,i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[i][j] = (*f)(me[i][j],b);
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::submat(const size_type idx1,const size_type idx2, const int len1, const int len2) const
{
   int ir,ic,i,j,t,k;
   if (len1 == 0) {
      ir = nrow - 1;
   } else {
      ir = idx1 + len1 - 1;
   }
   if (len2 == 0) {
      ic = ncol - 1;
   } else {
      ic = idx2 + len2 - 1;
   }
   Matrix<T> temp(len1,len2);
   for (i=idx1; i<= ir; ++i) {
      t = i - idx1;
      for (j=idx2; j<=ic; ++j) {
         k = j - idx2;
         temp[t][k] = me[i][j];
      }
   }
   return temp;
}
  /*
template<class T> inline Vector<T> Matrix<T>::submat(const Matrix<bool> &v) const {
    int n = 0;
    int nr = v.nrow;
    int nc = v.ncol;
    int i,j;
    for (i=0; i<nr; i++) for (j=0; j<nc; j++) if (v.me[i][j]) n++;

    Vector<T> temp(n);
    T *t_pt = temp.begin();
    for (i=0; i<nr; i++) for (j=0; j<nc; ++j) if (v.me[i][j]) *t_pt++ = me[i][j];
    return temp;
}
  */
template<class T> inline Matrix<T> Matrix<T>::transpose(void) const
{
   Matrix<T> temp(ncol,nrow);
   for (int j,i=0; i<ncol; ++i) for (j=0; j<nrow; ++j) temp[i][j] = me[j][i];
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::kron(const Matrix<T> &b) const
{
   Matrix<T> temp(nrow*b.nrow,ncol*b.ncol);
   int i,j,k,t,ii,iii,jj,jjj;
   for (i=0; i<nrow; ++i) {
      ii = i*b.nrow;
      for (j=0; j<ncol; ++j) {
         jj = j*b.ncol;
         iii = ii;
         for (k=0; k<b.nrow; ++k) {
            jjj = jj;
            for (t=0; t<b.ncol; ++t) temp[iii][jjj++] = me[i][j]*b[k][t];
            iii++;
         }
      }
   }
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::reshape(const size_type m,const size_type n,orientation orient) const
{
   Matrix<T> temp(m,n);
   int i,j,t,k;
   if (orient == ROW) {
      for (t=0,k=0,i=0; i<m; ++i) for (j=0; j<n; ++j) {
         if (k >= ncol) {
            t++;
            if (t >= nrow) return temp;
            k = 0;
         }
         temp[i][j] = me[t][k++];
      }
   } else if (orient == COLUMN) {
      for (t=0,k=0,i=0; i<m; ++i) for (j=0; j<n; ++j) {
         if (t >= nrow) {
            k++;
            if (k >= ncol) return temp;
            t = 0;
         }
         temp[i][j] = me[t++][k];
      }
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return temp;
}

template<class T> inline Matrix<T> Matrix<T>::diag(void) const
{
   int n = std::min(nrow,ncol);
   Matrix<T> temp(n,n);
   for (int i=0; i<n; ++i) temp[i][i] = me[i][i];
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::diag(const size_type k) const
{
   int i,n = std::min(nrow,ncol) - std::abs(k);
   Vector<T> temp(n);
   if (k >= 0) {
      for (i=0; i<n; i++) temp[i] = me[i][i+k];
   } else {
      for (i=0; i<n; i++) temp[i] = me[i-k][i];
   }
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::vec(orientation orient) const
{
   Vector<T> temp;
   if (!me) return temp;
   temp.reserve(nrow*ncol);
   int i,j,t=0;
   if (orient == ROW) {
      for (i=0; i<nrow; ++i) for (j=0; j<ncol; ++j) temp[t++] = me[i][j];
   } else if (orient == COLUMN)  {
      for (j=0; j<ncol; ++j) for (i=0; i<nrow; ++i) temp[t++] = me[i][j];
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::sum(orientation orient) const
{
   Vector<T> temp;
   if (!me) return temp;
   int i,j;
   T init;
   if (orient == ROW) {
      temp.reserve(nrow);
      for (i=0; i<nrow; ++i) {
         init = T();
         for (j=0; j<ncol; ++j) init += me[i][j];
         temp[i] = init;
      }
   } else if (orient == COLUMN)  {
      temp.reserve(ncol);
      for (j=0; j<ncol; ++j) {
         init = T();
         for (i=0; i<nrow; ++i) init += me[i][j];
         temp[j] = init;
      }
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::sumsq(orientation orient) const
{
   Vector<T> temp;
   if (!me) return temp;
   int i,j;
   T init;
   if (orient == ROW) {
      temp.reserve(nrow);
      for (i=0; i<nrow; ++i) {
         init = T();
         for (j=0; j<ncol; ++j) init += me[i][j]*me[i][j];
         temp[i] = init;
      }
   } else if (orient == COLUMN)  {
      temp.reserve(ncol);
      for (j=0; j<ncol; ++j) {
         init = T();
         for (i=0; i<nrow; ++i) init += me[i][j]*me[i][j];
         temp[j] = init;
      }
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::max(orientation orient) const
{
   Vector<T> temp;
   if (!me) return temp;
   int i,j;
   if (orient == ROW) {
      temp.reserve(nrow);
      for (i=0; i<nrow; ++i) {
         temp[i] = me[i][0];
         for (j=0; j<ncol; ++j) if (me[i][j] > temp[i]) temp[i] = me[i][j];
      }
   } else if (orient == COLUMN)  {
      temp.reserve(ncol);
      for (j=0; j<ncol; ++j) {
         temp[j] = me[0][j];
         for (i=0; i<nrow; ++i) if (me[i][j] > temp[j]) temp[j] = me[i][j];
      }
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return temp;
}

template<class T> inline Vector<T> Matrix<T>::min(orientation orient) const
{
   Vector<T> temp;
   if (!me) return temp;
   int i,j;
   if (orient == ROW) {
      temp.reserve(nrow);
      for (i=0; i<nrow; ++i) {
         temp[i] = me[i][0];
         for (j=1; j<ncol; ++j) if (me[i][j] < temp[i]) temp[i] = me[i][j];
      }
   } else if (orient == COLUMN)  {
      temp.reserve(ncol);
      for (j=0; j<ncol; ++j) {
         temp[j] = me[0][j];
         for (i=1; i<nrow; ++i) if (me[i][j] < temp[j]) temp[j] = me[i][j];
      }
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return temp;
}

template<class T> inline std::ostream& Matrix<T>::print(std::ostream &os) const
{
   int i,j;
   for (i=0; i<nrow; ++i) {
     for (j=0; j<ncol; ++j)  {
        os << ' ' ;
	os.width(8) ;
	os.precision(SESSION.output_precision) ;
	os << me[i][j];
     }
      os <<  "\n";
   }
   return os;
}

template<class T> inline Matrix<T>& Matrix<T>::sortby(orientation orient,const size_type t)
{
    int n = 1;
    Vector<int> by(n);
    by[0] = t;

    int nstack = 100;
    int l,r,i,j,k,j1;
    T a1,a2;
    Vector<T> x(n);                         // allocating for x
    Matrix<int> stack(2,nstack);              // allocating for stack
    int s = 0;
    stack[s][s] = 0;

   if (orient == COLUMN) {
      stack[1][0] = nrow-1;
      do {
         l = stack[0][s]; r = stack[1][s]; s -= 1;
         do {
            i = l; j = r; j1 = (l+r)/2;
            for (k=0; k<n; k++)  x[k] = me[j1][by[k]];
            do {
               k = 0;
               while (k < n) {
                  a1 = me[i][by[k]]; a2 = x[k];
                  if (a1 < a2) {
                     i += 1; k = 0;
                  }
                  else if (a1 == a2) {
                     k += 1;
                  }
                  else if (a1 > a2) break;
               }
               k = 0;
               while (k < n) {
                  a1 = me[j][by[k]]; a2 = x[k];
                  if (a1 > a2) {
                     j -= 1; k = 0;
                  }
                  else if (a1 == a2) {
                     k += 1;
                  }
                  else if (a1 < a2) break;
               }

               if (i <= j) {
                  for (k=0; k<ncol; k++) {
                     a1 = me[i][k];
                     me[i][k] = me[j][k];
                     me[j][k] = a1;
                  }
                  i += 1;  j -= 1;
               }
            } while (i <= j);
            if (i < r) {
               s += 1;
               if (s >= nstack-1) {
                  std::cerr << " NSTACK " << nstack << " too small in SORTQR\n";
                  s = -1;
                  break;
               }
               stack[0][s] = i;  stack[1][s] = r;
            }
            r = j;
         } while (l < r);
      } while (s >= 0);
   } else if (orient == ROW) {
      stack[1][0] = ncol - 1;
      do {
         l = stack[0][s]; r = stack[1][s]; s -= 1;
         do {
            i = l; j = r; j1 = (l+r)/2;
            for (k=0; k<n; k++)  x[k] = me[by[k]][j1];
            do {
               k = 0;
               while (k < n) {
                  a1 = me[by[k]][i]; a2 = x[k];
                  if (a1 < a2) {
                     i += 1; k = 0;
                  }
                  else if (a1 == a2) {
                     k += 1;
                  }
                  else if (a1 > a2) break;
               }
               k = 0;
               while (k < n) {
                  a1 = me[by[k]][j]; a2 = x[k];
                  if (a1 > a2) {
                     j -= 1; k = 0;
                  }
                  else if (a1 == a2) {
                     k += 1;
                  }
                  else if (a1 < a2) break;
               }
               if (i <= j) {
                  for (k=0; k<nrow; k++) {
                     a1 = me[k][i];
                     me[k][i] = me[k][j];
                     me[k][j] = a1;
                  }
                  i += 1;  j -= 1;
               }
            } while (i <= j);
            if (i < r) {
               s += 1;
               if (s >= nstack-1) {
                  std::cerr << " NSTACK " << nstack << " too small in SORTQR\n";
                  s = -1;
                  break;
               }
               stack[0][s] = i;  stack[1][s] = r;
            }
            r = j;
         } while (l < r);
      } while (s >= 0);
   } else {
      std::cerr << "unknow orientation\n";
      exit(1);
   }
   return *this;
}

/////////// the following are friend functions/operator /////////

/////////// the following are template functions/operator ///////////////
template<class T> inline Matrix<T> operator / (const T &a,const Matrix<T> &b)
{
   Matrix<T> temp(b.num_rows(),b.num_cols());
   for (int j,i=0; i<b.num_rows(); ++i) for (j=0; j<b.num_cols(); ++j) temp[i][j] = a/b[i][j];
   return temp;
}

template<class T> inline Vector<T> operator * ( Vector<T> &a, Matrix<T> &b)
{
   if (a.size() != b.num_rows()) throw exception("Matrix<T>::operator*: size not conformable");
   Vector<T> temp(b.num_cols());
   T** me = b.begin();
   T x;
   for (int i,j=0; j<b.num_cols(); ++j) {
      x = T();
      for (i=0; i<b.num_rows(); ++i)  x += a[i] * me[i][j];
      temp[j] = x;
   }
   return temp;
}

template<class T> inline Matrix<T> hadjoin (const Matrix<T> &a,const Matrix<T> &b)
{
   if (a.num_rows() != b.num_rows()) throw exception("Matrix<T>::hadjoin(): size not conformable");
   Matrix<T> temp(a.num_rows(),a.num_cols() + b.num_cols());
   for (int i=0; i<a.num_rows(); ++i) {
      memcpy(temp[i],a[i],sizeof(T)*a.num_cols());
      memcpy(&(temp[i][a.num_cols()]),b[i],sizeof(T)*b.num_cols());
   }
   return temp;
}

template<class T> inline Matrix<T> hadjoin (const Matrix<T> &a,const Vector<T> &b)
{
   if (a.num_rows() != b.size()) throw exception("Matrix<T>::hadjoin(): size not conformable");
   Matrix<T> temp(a.num_rows(),a.num_cols() + 1);
   for (int i=0; i<a.num_rows(); ++i) {
      memcpy(temp[i],a[i],sizeof(T)*a.num_cols());
      temp[i][a.num_cols()] = b[i];
   }
   return temp;
}

template<class T> inline Matrix<T> hadjoin (const Vector<T> &a,const Matrix<T> &b)
{
   if (a.size() != b.num_rows()) throw exception("Matrix<T>::hadjoin(): size not conformable");
   Matrix<T> temp(a.size(),1 + b.num_cols());
   for (int i=0; i<a.size(); ++i) {
      temp[i][0] = a[i];
      memcpy(&(temp[i][1]),b[i],sizeof(T)*b.num_cols());
   }
   return temp;
}

template<class T> inline Matrix<T> vadjoin (const Matrix<T> &a,const Matrix<T> &b)
{
   if (a.num_cols() != b.num_cols()) throw exception("Matrix<T>::vadjoin(): size not conformable");
   Matrix<T> temp(a.num_rows() + b.num_rows(),a.num_cols());
   for (int i,j=0; j<a.num_cols(); ++j) {
      for (i=0; i<a.num_rows(); ++i) temp[i][j] = a[i][j];
      for (i=0; i<b.num_rows(); ++i) temp[a.num_rows() + i][j] = b[i][j];
   }
   return temp;
}

template<class T> inline Matrix<T> vadjoin (const Matrix<T> &a,const Vector<T> &b)
{
   if (a.num_cols() != b.size())  throw exception("Matrix<T>::vadjoin(): size not conformable");
   Matrix<T> temp(a.num_rows() + 1,a.num_cols());
   for (int i,j=0; j<a.num_cols(); ++j) {
      for (i=0; i<a.num_rows(); ++i) temp[i][j] = a[i][j];
      temp[a.num_rows()][j] = b[j];
   }
   return temp;
}

template<class T> inline Matrix<T> vadjoin (const Vector<T> &a,const Matrix<T> &b)
{
   if (a.size() != b.num_cols())  throw exception("Matrix<T>::vadjoin(): size not conformable");
   Matrix<T> temp(1 + b.num_rows(),b.num_cols());
   for (int i,j=0; j<b.num_cols(); ++j) {
      temp[0][j] = a[j];
      for (i=0; i<b.num_rows(); ++i) temp[1 + i][j] = b[i][j];
   }
   return temp;
}

template<class T> inline Matrix<T> diag (const Vector<T> &a)
{
   Matrix<T> temp(a.size(),a.size());
   for (int i=0; i<a.size(); ++i) temp[i][i] = a[i];
   return temp;
}

template<class T> inline Matrix<T> diag (const Matrix<T> &a) { return a.diag(); }


template<class T> Matrix<T> sin  (const Matrix<T> &a) { return a.apply(std::sin); }
template<class T> Matrix<T> asin (const Matrix<T> &a) { return a.apply(std::asin); }
template<class T> Matrix<T> cos  (const Matrix<T> &a) { return a.apply(std::cos); }
template<class T> Matrix<T> acos (const Matrix<T> &a) { return a.apply(std::acos); }
template<class T> Matrix<T> tan  (const Matrix<T> &a) { return a.apply(std::tan); }
template<class T> Matrix<T> atan (const Matrix<T> &a) { return a.apply(std::atan); }
template<class T> Matrix<T> ceil (const Matrix<T> &a) { return a.apply(std::ceil); }
template<class T> Matrix<T> floor(const Matrix<T> &a) { return a.apply(std::floor); }
template<class T> Matrix<T> log  (const Matrix<T> &a) { return a.apply(std::log); }
template<class T> Matrix<T> log10(const Matrix<T> &a) { return a.apply(std::log10); }
template<class T> Matrix<T> exp  (const Matrix<T> &a) { return a.apply(std::exp); }
template<class T> Matrix<T> sqrt (const Matrix<T> &a) { return a.apply(std::sqrt); }
template<class T> Matrix<T> abs  (const Matrix<T> &a) { return a.apply(std::abs); }
template<class T> Matrix<T> erf  (const Matrix<T> &a) { return a.apply(erf); }
template<class T> Matrix<T> erfc (const Matrix<T> &a) { return a.apply(erfc); }

template<class T> bool all (const Matrix<T> &a)       { return a.all();}
template<class T> bool any (const Matrix<T> &a)       { return a.any();}

template<class T> Matrix<T> operator + (const T &a,const Matrix<T> &b) { return b + a; }
template<class T> Matrix<T> operator - (const T &a,const Matrix<T> &b) { return -b + a; }
template<class T> Matrix<T> operator * (const T &a,const Matrix<T> &b) { return b * a; }
template<class T> Matrix<T> operator / (const T &a,const Matrix<T> &b);
template<class T> std::ostream&  operator << (std::ostream &os, const Matrix<T> &a) { return a.print(os); }

} ////// end of namespace matvec


#endif
