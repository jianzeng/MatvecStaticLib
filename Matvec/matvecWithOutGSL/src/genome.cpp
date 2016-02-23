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

#include "genome.h"

namespace matvec {

Genome* new_Genome_vec(const unsigned m, GeneticDist *D)
{
   Genome *T = 0;
   if (m) {
      T = new Genome[m];
      check_ptr(T);
   }
   if (D) {
      for (unsigned i=0; i<m; i++) T[i].remodel(D);
   }
   return T;
}

Genome::Genome(GeneticDist *D)
{
   numchrom = 0;
   chromosome = 0;
   remodel(D);
}

Genome::Genome(const Genome& A)
{
   numchrom = 0;
   chromosome = 0;
   copyfrom(A);
}

void Genome::copyfrom(const Genome& A)
{
   if (this == &A) return;
   numchrom = A.numchrom;
   if (chromosome) {
     delete [] chromosome;
     chromosome=0;
   }
   if(numchrom>0){
     chromosome = new Chromosome [numchrom];  // each chrom should resize
   }
   else{
     chromosome = 0;
   }
   for (unsigned i=0; i<numchrom; i++) chromosome[i] = A.chromosome[i];
}

const Genome& Genome::operator=(const Genome& A)
{
   copyfrom(A);
   return *this;
}

void Genome::remodel(GeneticDist *D)
{
   numchrom = D->nchrom();
   ChromStruct *Chrom = D->chrom();
   if (chromosome) {
     delete [] chromosome;
     chromosome=0;
   }
   if (numchrom >0) {
     chromosome = new Chromosome [numchrom];  // each chrom should resize
   }
   else {
     chromosome = 0;
   }
   for (unsigned i=0; i<numchrom; i++) chromosome[i].resize(Chrom[i].nloci());
}

void Genome::resize(const unsigned nc, const unsigned nl)
{
   numchrom = nc;
   if (chromosome) {
     delete [] chromosome;
     chromosome=0;
   }
   if(numchrom>0){
     chromosome = new Chromosome [numchrom];  // each chrom should resize
   }
   else{
     chromosome = 0;
   }
   for (unsigned i=0; i<numchrom; i++) chromosome[i].resize(nl);
}

void Genome::release(void)
{
   if (chromosome) { delete [] chromosome; chromosome = 0; }
}

int Genome::chrom_id(const Chromosome& chrm) const
{
   int retval = -1;
   Chromosome *C;
   unsigned i,j,nl;
   for (i=0; i<numchrom; i++) {
      C = &(chromosome[i]);
      nl = C->nloci();
      for (j=0; j<nl; j++) {
         if (chrm.locus[j].allele != C->locus[j].allele) break; 
      }
      if (j==nl) {
         retval = C->id();
         break;
      }
   }
   return retval;
}

} /////////// end of namespace matvec
