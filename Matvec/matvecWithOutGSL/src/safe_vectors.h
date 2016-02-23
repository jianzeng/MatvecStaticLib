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



#ifndef hsafe_vectors_h   
#define hsafe_vectors_h
#include <vector>
#include <string>
#include <algorithm>
#include "session.h"
#include "vector.h"
#include "exception.h"
using namespace std;



template <class T>

class SafeSTLVector:public vector<T> {
public:
	vector<T> *base;
	string name;
	typedef size_t size_type;
	/**
	 * It returns a reference to element i with range-cheching.
	 * Note that element index starting from 0.
	 */

	SafeSTLVector(void){name = "unknown";}
	T&  operator [] (const size_type i){
		if (i >= this->size()) {
			cout << "out of range error in SafeSTLVector: " << name << endl;
			cout << "i = " << i << endl;
			cout << "size of vector = " << this->size() << endl;
			throw matvec::exception("SafeSTLVector: out of range error");
			exit(1);
		}
		return *(this->begin() + i) ;
	}
};



template <class T>

class SafeMVVector:public matvec::Vector<T> {
  public:
  matvec::Vector<T> *base;
  string name;
  typedef size_t size_type;
/**
 * It returns a reference to element i with range-cheching.
 * Note that element index starting from 0.
 */
  T&  operator [] (const size_t i){
    if (i >= this->ne) {
      cout << "out of range error in SafeMVVector: " << name << endl;
      exit(1);
    }
    base = this;
    return *(this->begin() + i) ;
  }
  SafeMVVector &  assign(const SafeMVVector &a);
  SafeMVVector &  assign(const matvec::Vector<T> &a);
  SafeMVVector& operator =  (const SafeMVVector &a)   { return assign(a); }
  SafeMVVector& operator =  (const matvec::Vector<T> &a) { return assign(a); }


};



template<class T> inline SafeMVVector<T>& SafeMVVector<T>::assign(const matvec::Vector<T> &a)
{
   matvec::Vector<T>::copy(a);
   return *this;
}



template<class T> inline SafeMVVector<T>& SafeMVVector<T>::assign(const SafeMVVector<T>& a)

{
   matvec::Vector<T>::copy(a);
   return *this;
}


#endif

