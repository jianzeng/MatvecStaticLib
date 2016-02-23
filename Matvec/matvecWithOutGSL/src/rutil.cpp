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

#include "rutil.h"
using namespace std;

namespace matvec{
  void Tokenizer::getTokens(const string &str, const string &sep){
    // Authors: Rohan L. Fernando
    // this is essentially the split funciton from matvec
    // I guess split was written by Tianlin Wang
    // (2005)
    // Contributors:
    clear();
    string::size_type begidx,endidx;
    begidx = str.find_first_not_of(sep);
    while (begidx != string::npos) {
      endidx = str.find_first_of(sep,begidx);
      if (endidx == string::npos) endidx = str.length();
      push_back(str.substr(begidx,endidx - begidx));
      begidx = str.find_first_not_of(sep,endidx);
    }
  }

  int Tokenizer::getIndex(string str){
    // Authors: Rohan L. Fernando
    // (2005) 
    // Contributors:
    for (unsigned i=0;i<size();i++){
      if((*this)[i]==str){
	return i;
      }
    }
    return -1;
  }

  double SparseCij::retrieve_cij(const unsigned i, const unsigned j){
    // Authors: Rohan L. Fernando
    // (2005) 
    // Contributors:
    idx ind;
    unsigned ii,jj;
    ii = i;
    jj = j;
    if (i>j) {
      // symmetry is assumed, and only upper triangular elements saved
      ii = j;
      jj = i;
    } 
    ind.i = ii;
    ind.j = jj;
    map<const idx , double>::iterator cit = C.find(ind);
    if (cit != C.end()) {
      return cit->second;
    }
    else {
      return -1.0;
    }
  }

  void SparseCij::put_cij(const unsigned i, const unsigned j, double cij){
    // Authors: Rohan L. Fernando
    // (2005) 
    // Contributors:
    idx ind;
    unsigned ii,jj;
    ii = i;
    jj = j;
    if (i>j) {
      // symmetry is assumed, and only upper triangular elements saved
      ii = j;
      jj = i;
    }
    ind.i = ii;
    ind.j = jj;
    C[ind] = cij;
  }

} ////// end of namespace matvec

