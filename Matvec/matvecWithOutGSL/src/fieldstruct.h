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

#ifndef MATVEC_FIELDSTRUCT_H
#define MATVEC_FIELDSTRUCT_H

#include <string>

namespace matvec {
/*!
   \class FieldStruct  fieldstruct.h
   \brief a structure for a data column

   \sa Field Data
*/

class FieldStruct {
   protected:
      char     mytype,myclassi; // myclassi='I','T','C','U','P','N','F'(default)
      unsigned mynlevel, myindex,mynmiss;  // mytype = 'S','I', 'F'(default)
      double   mymean, mystd;

      std::string   myname;                     // myindex = field_id - 1

      void copyfrom(const FieldStruct& A);

   public:

      FieldStruct(void);
      FieldStruct(const FieldStruct& A);             // copy constructor
      ~FieldStruct(void) {;}

      const FieldStruct&  operator=(const FieldStruct& A);

      void        index(const unsigned k)   {myindex = k;}
      void        nlevel(const unsigned k)  {mynlevel = k;}
      void        nmiss(const unsigned k)   {mynmiss = k;}
      void        mean(const double x)      {mymean = x;}
      void        std(const double x)       {mystd = x;}
      void        count_miss(const int k)   {mynmiss += k;}
      void        classi(const char c)      {myclassi = c;}
      void        type(const char c)        {mytype = c;}
      void        name(const std::string &n)      {myname = n;}

      unsigned    index(void)  const {return myindex;}
      unsigned    nlevel(void) const {return mynlevel;}
      unsigned    nmiss(void)  const {return mynmiss;}
      double      mean(void)   const {return mymean;}
      double      std(void)    const {return mystd;}
      char        classi(void) const {return myclassi;}
      char        type(void)   const {return mytype;}
      const std::string name(void)   const {return myname;}

};

inline FieldStruct::FieldStruct(void)
{ mytype='F'; myclassi='F'; mynlevel=0; myindex=0; mynmiss=0; }

extern int fieldindex(const std::string &fw,const FieldStruct* field_vec,
                      const unsigned n);
}
#endif
