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

#ifndef MATVEC_MATRIXWBG_H
#define MATVEC_MATRIXWBG_H

#include <string>
#include <fstream>
#include "session.h"
#include "vector.h"
#include "bg.h"
#include "matrix.h"

namespace matvec {
  /*!
    \class Matrix  vector.h
    \brief  A vector is a on-dimensional array with double precision

    \sa Matrix
  */

  template<class T> class Matrixwbg : public Matrix <T> {
  
  public:
    typedef size_t size_type;
   Matrixwbg(void)                                            { initialize(0,0,0); }  //Constructor 1
   Matrixwbg(const size_type m,const size_type n)             { initialize(m,n,0); }  //Constructor 2
   Matrixwbg(const size_type m,const size_type n,const T** a) { initialize(m,n,a); }  //Constructor 3
   Matrixwbg(const Matrixwbg<T>& a)                           { initialize(a.nrow,a.ncol,(const T**)a.me); }  //Constructor 4
   Matrixwbg(const Matrix<T>& a)                              { initialize(a.nrow,a.ncol,(const T**)a.me); }  //Constructor 5
  };


} ////// end of namespace matvec


#endif
