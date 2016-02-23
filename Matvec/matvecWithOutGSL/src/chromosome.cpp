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

#include "util.h"
#include "chromosome.h"

namespace matvec {

/////////////////////////// Locus class /////////////////////////////

const Locus&  Locus::operator=(const Locus& A)
{
   if (this != &A) {
      allele  = A.allele;
      effect  = A.effect;
   }
   return *this;
}

/////////////////////////// Chromosome class /////////////////////////////

const Chromosome&  Chromosome::operator=(const Chromosome& A)
{
   copyfrom(A);
   return *this;
}

void  Chromosome::copyfrom(const Chromosome& A)
{
   if (this == &A) return;
   chrom_id   = A.chrom_id;
   chrom_freq = A.chrom_freq;
   resize(A.numloci);
   for (unsigned i=0; i<numloci; i++) locus[i] = A.locus[i];
}

void Chromosome::resize(const unsigned nl)
{
   if (numloci == nl) return;
   numloci = nl;
   if (locus) {
     delete locus;
     locus = 0;
   }
   if (numloci>0){
     locus = new Locus [numloci];
   }
   else {
     locus = 0;
   }
   assert(locus);
}
} //////// end of namespace matvec

