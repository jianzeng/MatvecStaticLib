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

#ifndef Genome_H
#define Genome_H

#include "geneticdist.h"
#include "chromosome.h"

namespace matvec {
class Chromosome;
/*!
   \class Genome  Genome.h
   \brief a genome consists of chromosomes

   \sa Chromosome
*/

class Genome {
   protected:
      unsigned  numchrom;
      void      copyfrom(const Genome& A);
   public:
      Chromosome    *chromosome;

      Genome(void) {numchrom=0; chromosome = 0;}
      Genome(GeneticDist *D);
      Genome(const Genome& A);             // copy constructor
      ~Genome(void){release();};

      const Genome&  operator=(const Genome& A);

      unsigned    nchrom(void) {return numchrom;}
      unsigned    size(void) {return numchrom;}
      int         chrom_id(const Chromosome& chrm) const;
      void        release(void);
      void        remodel(GeneticDist *D);
      void        resize(const unsigned nc, const unsigned nl);
};
extern Genome* new_Genome_vec(const unsigned m,
                              GeneticDist *D = 0);
}
#endif
