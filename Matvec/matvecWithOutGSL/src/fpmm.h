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

#ifndef FPMM_H
#define FPMM_H

#include "dblock.h"
#include "vector.h"

/*!
   \class Fpmm
   \brief Calculate transmission matrix for FPMM model

   We do this "recursively".
   We first compute the joint probabliities for the polygenic numbers
   for mother, father, and offspring for just one locus. Then we assuming the
   second locus is independent, we can compute the polygenic numbers for mother, father,
   and offspring for two loci.
*/

namespace matvec {
class Fpmm {
   protected:
      Dblock trm;   // transmission matrix 
      int ndim;     // 2 times nloci + 1 (number of polygenic classes)
                    // goes from 0,1,...,2xnloci
      Vector<double> genval;
      Vector<double> genfreq;
      double freq;


   public:

      Fpmm(void) {ndim=0;freq=0.5;var=1.0;}
      Fpmm(Fpmm& F);
      Fpmm(int nl, double frq, double vr) {
	setup_trm(nl, frq);
	setup_genval(vr);
      }
      double var;
      void initial(int nl, double frq, double vr) {
//        cout << "inside Fpmm.initial" << endl;
	setup_trm(nl, frq);
	setup_genval(vr);
      //  cout << "exiting Fpmm.initial" << endl;
      }
      void setup_trm(int nl, double frq);
      void setup_genval(double vr);
      double getpr(int fpgn, int mpgn, int opgn){return trm[fpgn][mpgn][opgn];};
      Vector<double> getgv(void){return genval;}
      Vector<double> getgenfreq(void){return genfreq;}
      void print_trm(void){std::cout << trm;}
      void print_genval(void){std::cout << genval;}
      void print_genfreq(void){std::cout << genfreq;}
      int get_dim(void) {return ndim;}
      int display(void);
};
}

//BRS
#endif


