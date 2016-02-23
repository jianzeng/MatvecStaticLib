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

#ifndef MATVEC_CHROMOSOME_H
#define MATVEC_CHROMOSOME_H

#include <cassert>
#include <cstdlib>
#include <cstdio>

namespace matvec {

/////////////////////// Locus class /////////////////////////////////
/*!
   \class Locus Chromosome.h
   \brief locus

   \sa Chromosome
*/
class Locus {
   public:
      double    effect;
      int       allele;       // can be QTL allele or ML allele

      Locus(void){effect=0.0; allele=0;}

      const Locus& operator=(const Locus& A);
};

/////////////////////// Chromosome class /////////////////////////////////
/*!
   \class Chromosome
   \brief A chromosome

   \sa Genome
*/
class Chromosome {
   protected:
      unsigned  chrom_id, numloci;
      double    chrom_freq;
      void      copyfrom(const Chromosome& A);

   public:
      Locus   *locus;

      Chromosome(void){numloci=0; locus = 0;}
      Chromosome(const Chromosome& A){copyfrom(A);}
      ~Chromosome(void){ if(locus) {delete [] locus; locus = 0;}}

      const Chromosome& operator=(const Chromosome& A);

      void      id(const unsigned newid) {chrom_id = newid;}
      void      freq(const double new_f) {chrom_freq = new_f;}
      void      resize(const unsigned nl);

      double    freq(void) const {return chrom_freq;}
      unsigned  nloci(void) {return numloci;}
      unsigned  id(void) const {return chrom_id;}
};

extern Chromosome* new_Chromosome_vec(const unsigned n);
}

#endif
