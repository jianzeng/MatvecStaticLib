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

#ifndef MIMUTIL_H
#define MIMUTIL_H
#include <map>
#include <string>
#include "safe_vectors.h"
#include <iostream>
#include <iomanip>

using namespace std;

namespace matvec{
	template <class T>
    class Recoder : public map<T,unsigned> {
		unsigned count;
public:
		Recoder(void){count=0;}
		unsigned code(T s){
			typename map<T,unsigned>::iterator mapit = this->find(s);
			if(mapit == this->end()){
				(*this)[s] = ++count;
				return count;
			}
			else {
				return (*mapit).second;
			}
		}
		void display_codes(ostream & os = std::cout){
			typename Recoder::iterator it;
			for (it=this->begin(); it!=this->end();it++){
				os << (*it).first << " " << (*it).second << std::endl;
			}
		}
	};
  
  class Tokenizer:public SafeSTLVector<string> {
  public:
    void getTokens(const string &str, const string &sep);
    int  getIndex(string str);
  };
  
  
  class idx {
  public:
    unsigned i,j;
    bool operator<(const idx y) const
    {
      if (i < y.i) {
        return 1;
      }
      if (i > y.i) {
        return 0;
      }
      else {
        if (j < y.j) {
          return 1;
        }
        else {
          return 0;
        }
      }
    }
  };

  class SparseCij {
    map<const idx , double> C;
  public:
    double retrieve_cij(const unsigned i, const unsigned j);
    void   put_cij(const unsigned i, const unsigned j, double cij);
    void   clear(void){C.clear();}
  };
}
#endif  
