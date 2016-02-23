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

#ifdef WIN32
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#include "exception.h"
#include "session.h"
#include <unistd.h>
#include <stdlib.h>


namespace matvec{
Session SESSION;

void Session::initialize(const std::string &d)
{

   tmpdir = d;
#ifdef WIN32
   _fmode=O_BINARY;
#endif
   DIR *dir = opendir(tmpdir.c_str());
   if (!dir) throw exception("Session::initialize(): " + d + ": No such directory");
   closedir(dir);


   warning = 1;
   rand_seed = -987654321;
   output_precision = 6;
   output_line_width = 80;
   epsilon = 1.0e-14; // was -14
   return;
}

std::string Session::mktemp(void)
{
   if (tmpdir == "") throw exception("Session::mktemp():  session is not initialized");
   char *fname;
#ifndef WIN32
   fname= new char [tmpdir.length() + 11];
   strcpy(fname,tmpdir.c_str());
   strcat(fname,"/mv.XXXXXX");
#endif

   int tfile;
   if (
#ifdef WIN32
       !(fname=tempnam(tmpdir.c_str(),"mv."))
#else
       mkstemp(fname)==-1
#endif
       )   throw exception("Session::mktemp(): cannot make temporary filename");
   tmpfile.push_back(fname);
#ifndef WIN32
   if(fname){
     delete [] fname;
     fname=0;
   }
#endif
   return tmpfile.back();
}

void Session::clear(void)
{
   if (tmpfile.size() == 0) return;
   std::vector<std::string>::iterator it = tmpfile.begin();
   for (it = tmpfile.begin(); it != tmpfile.end(); ++it) std::remove(it->c_str());
   tmpfile.clear();
   return;
}
} //////// end of namespace matvec
