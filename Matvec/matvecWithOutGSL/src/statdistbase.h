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

#ifndef MATVEC_STAT_DIST_BASE_H
#define MATVEC_STAT_DIST_BASE_H

#include "vector.h"

namespace matvec {
class doubleMatrix;


/**
 * A base class for statistical distributions.
 *
 * @see GeneticDist
 */

class StatDistBase {
   protected:
      std::string  distname;

   public:
      StatDistBase(void) {distname="NullDist";}
      virtual ~StatDistBase(void) {;}

      virtual const std::string name(void) const {return distname;}
      virtual void              display(void) const  = 0;
      virtual void              sample(Vector<double>& x) const = 0;;
      virtual void              sample(doubleMatrix& x) const = 0;;
      virtual double            sample(void) const  = 0;
      virtual Vector<double>    sample(unsigned n) const  = 0;
      virtual doubleMatrix      sample(unsigned m,unsigned n) const  = 0;
      virtual double            mean(void) const  = 0;
      virtual double            variance(void) const  = 0;
      virtual double            pdf(const double x) const  = 0;
      virtual double            cdf(const double x) const  = 0;
      virtual double            mgf(const double x) const  = 0;
      virtual double            inv(const double p) const  = 0;
      virtual double            nonct(const double cv,const double p) {std::cerr << "error\n"; return 0.0;}
      virtual double            parameter(const int k) const  = 0;
      virtual void              parameter(const int k,const double x)  = 0;
};
}

#endif

