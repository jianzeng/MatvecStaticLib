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

#include <string>
#include <cstring>
#include <cstdarg>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "util.h"

namespace matvec {

static const int MATVEC_MESSAGE_BUF_LEN = 1024;

void check_ptr(void *a)
{
   ;
}

/**
 * split string into pieces of tokens into a vector.
 * @param str a string to be splitted
 * @param sep the separator
 * @param strvec a pointer to vector,which contains tokens.
 * Note that strvec is cleared first, so that the first token is always referred by strvec->begin().
 */
int split(const std::string &str,const std::string &sep, std::vector<std::string> *strvec)
{
   std::string::size_type begidx,endidx;
   begidx = str.find_first_not_of(sep);
   int ntoken = 0;
   while (begidx != std::string::npos) {
      ntoken++;
      endidx = str.find_first_of(sep,begidx);
      if (endidx == std::string::npos) endidx = str.length();
      if (strvec) strvec->push_back(str.substr(begidx,endidx - begidx));
      begidx = str.find_first_not_of(sep,endidx);
   }
   return ntoken;
}

void warning(const char format[],...)
{
//   if (!WARNING ) return;
   char *buf = new char [MATVEC_MESSAGE_BUF_LEN];
   va_list args;
   va_start(args,format);
   std::cerr << "\n ***WARNING***\n";
#ifdef SYSV
   int len;
   len = vsprintf(buf,format,args);
   assert(len < MATVEC_MESSAGE_BUF_LEN);
#else
   //    std::vsprintf(buf,format,(void *)args);
   std::vsprintf(buf,format,args); //changed to compile (LRT 11/30/00)
#endif
   va_end(args);
   std::cerr << " " << buf << std::endl;
   if(buf){
     delete [] buf;
     buf=0;
   }
}

bool validline(char *line)
{
   char *s;
   if (s = strstr(line,"//")) s[0] = '\0';
   if(s=strpbrk(line,"\r\n")) s[0] = '\0';
   int i=0;
   while (line[i]==' ') i++;
   if (line[i]=='\0') return false;
   else return true;
}

unsigned est_nze(const unsigned dim, const unsigned ntrait,const int pedsize)
{
   // ***********************************************
   // *   get a suggested number of max_nz elements
   // *  (1) calculate max_nz for single trait model
   // *  (2) then multiply it by ntrait*ntrait
   // ************************************************

   if (dim==0) return 0;
   double onetraitdim = static_cast<double>(dim)/static_cast<double>(ntrait);
   double nze, halfsize = onetraitdim*(onetraitdim - 1.0)/2.0;
   if(halfsize < 1000.0) {
      nze = halfsize;
   }
   else if(halfsize < 5000.0) {
      nze = 0.7 * halfsize;
   }
   else if(halfsize < 1.0e+4) {
      nze = 0.4 * halfsize;
   }
   else if(halfsize < 1.0e+5) {
      nze = 5e-1 * halfsize;
   }
   else if(halfsize < 1.0e+6) {
      nze = 5e-2 * halfsize;
   }
   else if(halfsize < 5.0e+6) {
      nze = 1e-2 * halfsize;
   }
   else if(halfsize < 1.0e+7) {
      nze = 8e-3 * halfsize;
   }
   else if(halfsize < 5.0e+7) {
      nze = 6e-3 * halfsize;
   }
   else if(halfsize < 1.0e+8) {
      nze = 4e-3 * halfsize;
   }
   else if(halfsize < 5.0e+8) {
      nze = 2e-3 * halfsize;
   }
   else if(halfsize < 1.0e+9) {
      nze = 9e-4 * halfsize;
   }
   else if(halfsize < 3.0e+9) {
      nze = 7e-4 * halfsize;
   }
   else if(halfsize < 5.0e+9) {
      nze = 5e-4 * halfsize;
   }
   else if(halfsize < 7.0e+9) {
      nze = 3e-4 * halfsize;
   }
   else if(halfsize < 1.0e+10) {
      nze = 9e-5 * halfsize;
   }
   else if(halfsize < 3.0e+10) {
      nze = 7e-5 * halfsize;
   }
   else if(halfsize < 5.0e+10) {
      nze = 5e-5 * halfsize;
   }
   else if(halfsize < 7.0e+10) {
      nze = 3e-5 * halfsize;
   }
   else if(halfsize < 9.0e+10) {
      nze = 1e-5 * halfsize;
   }
   else {
      nze = 1e-6 * halfsize;
   }

   nze += static_cast<double>(pedsize*3);
   nze = nze*ntrait*(ntrait+1.0)/2.0 + dim;
   return static_cast<unsigned>(nze);
}

//end namespace matvec
}

