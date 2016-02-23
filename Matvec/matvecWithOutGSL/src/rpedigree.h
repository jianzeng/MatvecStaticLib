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

#ifndef PMap_H
#define PMap_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cctype>
#include <sstream>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include "rutil.h"
#include "doublematrix.h"

using namespace std; 

namespace matvec{
  class PNode {
  public:
    int    ind, sire, dam;
    double f;
    string ind_str, sire_str, dam_str;
    
    PNode(string indstr, string sirestr, string damstr){
      
      ind  = -1; 
      sire = -1; 
      dam  = -1;
      f    = -1.0;
      
      ind_str  = indstr;
      sire_str = sirestr;
      dam_str  = damstr;
    }
  };
  
  struct pcomp:public binary_function<double, double, bool> {
    bool operator()(PNode *x, PNode *y) { return x->ind < y->ind; }
  };
  
  class RPedigree : public map<string,PNode*> {
  public:
    RPedigree(void){
		orderedPed = false;
	}
    unsigned COUNT;
    SparseCij SpCij;
    SafeSTLVector <PNode*> pedVector;
    Recoder<string> coder; 
	Tokenizer colName;
	bool orderedPed;
	void putColNames(string str);
    void inputPed(string fileString);
    void displayPed(void);
    void generateEntriesforParents(void);
	void seqnPed(void);
    void codePed();
    void code(PNode *ptr);
    void calc_inbreeding(void);
    void makePedVector(void);
    void fillCoder(void);
    double get_rij(int i, int j);
    void output(char* ped);
    void addAinv(matvec::doubleMatrix& lhs, unsigned startRow, unsigned statCol, double ratio);
  };
} 
#endif
