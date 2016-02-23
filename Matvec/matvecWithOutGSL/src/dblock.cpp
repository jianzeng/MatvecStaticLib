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

#include "dblock.h"

namespace matvec {

Dblock::Dblock (int nd, int nr, int nc){ 
  if (nd<=0) {
    std::cerr <<"Depth of block is not positive \n" << std::flush;
  }
  else if (nr<=0) {
    std::cerr <<"Number of rows Block is not positive \n" << std::flush;
  }
  else if (nc<=0) {
    std::cerr <<"Number of columns Block is not positive \n" << std::flush;
  }
  else {
    ndim=nd;
    nrow=nr; 
    ncol=nc;
    block.resize(nd);
    for (int i=1;i<=nd;i++){
      block(i).resize(nr,nc);
    }
  }
}

void Dblock::resize (int nd, int nr, int nc){
  if (nd<=0) {
    std::cerr <<"Depth of block is not positive \n" << std::flush;
  }
  else if (nr<=0) {
    std::cerr <<"Number of rows Block is not positive \n" << std::flush;
  }
  else if (nc<=0) {
    std::cerr <<"Number of columns Block is not positive \n" << std::flush;
  }
  else {
    ndim=nd;
    nrow=nr; 
    ncol=nc;
    block.resize(nd);
    for (int i=1;i<=nd;i++){
      block(i).resize(nr,nc);
    }
  }
}

void Dblock::resize (int nd, int nr, int nc, double init){
  if (nd<=0) {
    std::cerr <<"Depth of block is not positive \n" << std::flush;
  }
  else if (nr<=0) {
    std::cerr <<"Number of rows Block is not positive \n" << std::flush;
  }
  else if (nc<=0) {
    std::cerr <<"Number of columns Block is not positive \n" << std::flush;
  }
  else {
    ndim=nd;
    nrow=nr; 
    ncol=nc;
    block.resize(nd);
    for (int i=1;i<=nd;i++){
      block(i).resize(nr,nc,init);
    }
  }
}

void Dblock::init(double val){
  int i,j,k;
  for (i=0;i<ndim;i++){
       for (j=0;j<nrow;j++){
	 for (k=0;k<ncol;k++){
             block[i][j][k]=val;
	 }
       }
  }
}

std::ostream& operator << (std::ostream &stream, Dblock& D)
{
   stream << "\n";
   int k,i,j,jstart,jend,ncols_remaining;

   for (k=1;k<=D.ndim; k++){
    ncols_remaining = D.ncol;
    jend = 0;
    while (ncols_remaining > 0) {
      jstart = jend + 1;
      jend = jstart + ncols_remaining - 1;
      ncols_remaining -= 10;
      if (ncols_remaining < 0) {
        ncols_remaining = 0;
      }
      else {
        jend = jstart + 9;
      }	
      stream << "\n Block " << k <<"\n" << "          ";
      stream.width(2);
      for (j=jstart;j<=jend;j++){
        stream << " column " << j;
      }
      stream << "\n";
      for (i=1;i<=D.nrow; i++) {
        stream << "\n" << "Row ";
        stream.width(4);
        stream << i << "  ";
        for (j=jstart;j<=jend;j++){
          stream.width(9);
	  stream << D(k,i,j);
        } 
      }
    }
   stream << "\n" ;     
   }
   stream << "\n" << "\n" ;     
   return stream;
}
}
