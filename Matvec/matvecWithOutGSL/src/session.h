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

#ifndef MATVEC_SESSION_H
#define MATVEC_SESSION_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <dirent.h>
#include "MersenneTwister.h"
 
namespace matvec {
/**
 * A session class.
 */
class Session {
public:
   MTRand mtr;
   int       rand_seed,warning,output_precision,output_line_width;
   double    epsilon;

   Session(void) : tmpdir("") { }
   ~Session(void) { clear(); }

   Session&   operator = (const Session &s) { tmpdir = s.tmpdir; tmpfile=s.tmpfile;return(*this);}

   void        initialize(const std::string &td);
   void        clear(void);
   std::string mktemp(void);

private:
   std::string tmpdir;
   std::vector<std::string> tmpfile;
};

extern Session SESSION;   ///// one and only one matvec global variable
} //////// end of namespace  matvec
#endif
