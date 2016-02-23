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

#ifndef MATVEC_PEDIGREE_H
#define MATVEC_PEDIGREE_H

#include <iostream>
#include <string>

#include "doublematrix.h"

//////////////////////////////////////////////////////////////////
// nang is the length of ped_storage, nang = totalna + numgroup
//     if it is a not grouping pedogree, then nang = totalna
///////////////////////////////////////////////////////////////////

/*!
 * \struct PedNode mvpedigree.h
 * \brief a node structure used in Pedigree object
 *
 * \sa Pedigree
 */
namespace matvec {
struct PedNode {
   char mysex, *myname;
   unsigned myid,mother_id,father_id,namelen;
};

/*!
 *  \class Pedigree  mvpedigree.h
 *  \brief a pedigree
 *
 *  \sa Population
 */
class Pedigree { 
   protected:
      std::string binfname;
      unsigned   totalna, basesize, numgroup,nang;
      Vector<double>     fvec;
      int        fdone, ped_on_disk, ped_in_memory,groupingped;
      PedNode    *ped_storage;
      void       renum(void);

   public:
      enum ped_type {raw, standard};   //default = raw pedigree
      static Pedigree::ped_type type;


      PedNode    **pedmember;
      unsigned   maxnamelen,sex_confirmed;

      Pedigree(char pedtype[]="");
      Pedigree(Pedigree& P);
      ~Pedigree(void){release();}

      const Pedigree&  operator=(Pedigree& P) {copyfrom(P); return *this;}
      friend std::ostream&  operator<<(std::ostream& stream, Pedigree& A);

      int        inbcoef_done(void) const {return fdone;}
      int        in_memory(void) const {return ped_in_memory;}
      void       copyfrom(Pedigree& P);
      void       save_pedsheet(const int relse=1);
      void       input_pedsheet(void);
      void       release_pedsheet(void);
      void       resize(const unsigned na);
      unsigned   input(const std::string &pfname,
                       const std::string &fmt = "individual mother father",
                       const std::string &pedtype = "raw");
      unsigned   input_standard(std::istream& pedfile,unsigned ind,
                                unsigned mother,unsigned father);
      unsigned   ngroup(void) const {return numgroup;}
      unsigned   size(void) const {return totalna;}
      unsigned   nbase(void) const {return basesize;}
      Vector<double>*    inbcoef_quaas(int keepinmem=0);
      Vector<double>*    inbcoef_meuwissen(int keepinmem=0);
      Vector<double>*    inbcoef(int keepinmem=0) {return inbcoef_meuwissen(keepinmem);}
      doubleMatrix     mat(void);
      double     logdet(int keepinmem=0);
      doubleMatrix     rela(int keepinmem=0);
      doubleMatrix     ainv(int keepinmem=0);
      doubleMatrix     reld(int keepinmem=0);
      Pedigree&  sample(const unsigned maxsize,const unsigned maxg,
                        const unsigned sire0=1, const unsigned dam0=1,
                        const double imrate=0.1,const unsigned parent=1,
                        const int no_po=1, const int no_fsib=1,
                        const double sexratio=0.5);
      void       out_to_stream(std::ostream& stream,const int style=0,
                               const int ic=0);
      void       save(const std::string &fname,const int style=0,
                      const int io_mode=std::ios::out);//SDK | std::ios::noreplace);
      void       display(const std::string msg= "",const int style=0,const int ic=0);
      void       release(void);
      const std::string diskfname(void) const {return binfname; }
};
}
#endif
