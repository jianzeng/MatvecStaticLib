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

#ifndef MATVEC_DBLOCK_H
#define MATVEC_DBLOCK_H

#include <iostream>
#include <vector>
#include "doublematrix.h"

namespace matvec {
 /*!
   \class Dblock
   \brief This is a array of Matrix

   \sa Array Matrix
*/

class Dblock
{
protected:
  int ndim,nrow,ncol;
  Vector < doubleMatrix > block;

public:
  int iw; // just a flag
  Dblock(void) {ndim=0,nrow=0,ncol=0;}          //Constructor 1
  Dblock(int nd, int nr, int nc);               //Constructor 2

  void resize(int nd, int nr, int nc);
  void resize (int nd, int nr, int nc, double init);
  doubleMatrix & operator[](int i){return block[i];};
  double& operator()(int i,int j, int k){return block(i)(j,k);};
  int get_ndim(void) {return ndim;}
  int get_nrow(void) {return nrow;}
  int get_ncol(void) {return ncol;}
  void init(double val);
  friend std::ostream& operator << (std::ostream &stream, Dblock& D);
}; 
}

#endif
