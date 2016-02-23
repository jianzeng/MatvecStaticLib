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

#ifndef MATVEC_UTIL_H
#define MATVEC_UTIL_H

#include <vector>
#include <iostream>
#include <string>
#include <cctype>
#include <sstream>

namespace matvec {
/*!
   \file mvutil.h
   \brief all utility routines are defined here
*/

class BG;

extern void check_ptr(void *a);
extern int split(const std::string &str,const std::string &sep, std::vector<std::string> *strvec=0);
extern void warning(const char format[],...);

extern bool     validline(char *line);
extern unsigned ginverse1(double **a,const unsigned n,double& lgdet,int mode,const double tol);
extern unsigned ginverse1(BG **a,const unsigned n,BG& lgdet,int mode,const double tol);
extern long     next_prime(long n);
extern unsigned est_nze(const unsigned dim,const unsigned ntrait,const int ped);

template <class T> std::string toString(T t, std::ios_base & (*f)(std::ios_base&)){
	std::ostringstream oss;
	oss<< f << t;
	return oss.str();
}

} /////////// end of namespace

#endif
