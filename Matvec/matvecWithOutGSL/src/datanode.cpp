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

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#include "datanode.h"

namespace matvec {

DataNode::DataNode(const DataNode& A)
{
   copyfrom(A);
}

void DataNode::copyfrom(const DataNode& A)
{
   if (this == &A) return;
   data = A.data;
   missing = A.missing;
}

const DataNode& DataNode::operator=(const DataNode& A)
{
   copyfrom(A);
   return *this;
}

std::ostream&  operator<<(std::ostream& stream, const DataNode& A)
{
   std::cout.precision(6);
   if (A.missing) {
      stream << " " << std::setw(12) << ".";
   }
   else {
      stream << " " << std::setw(12) << A.data.double_value;
   }
   return stream;
}
}  //////////// end of namespace matvec
