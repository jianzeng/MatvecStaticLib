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

#include "fieldstruct.h"

namespace matvec {

int fieldindex(const std::string &fw,const FieldStruct* field_vec,const unsigned n)
{
   int k = -1;
   for (unsigned i=0; i<n; i++) {
      if (field_vec[i].name() == fw) {
         k = i; break;
      }
   }
   return k;
}

FieldStruct::FieldStruct(const FieldStruct& A) {copyfrom(A);}

void FieldStruct::copyfrom(const FieldStruct& A)
{
   if (this == &A) return;
   mytype   = A.mytype;
   myclassi = A.myclassi;
   mynlevel = A.mynlevel;
   mynmiss  = A.mynmiss;
   myindex  = A.myindex;
   myname   = A.myname;
}

const FieldStruct& FieldStruct::operator=(const FieldStruct& A)
{
   copyfrom(A);
   return *this;
}

} ///////// end of namespace matvec

