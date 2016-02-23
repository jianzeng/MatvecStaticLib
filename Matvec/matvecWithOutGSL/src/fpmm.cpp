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

#include "fpmm.h"

namespace matvec {

Fpmm::Fpmm(Fpmm& F){
 trm  = F.trm;
 ndim = F.ndim;
 genval = F.genval;
 genfreq = F.genfreq;
 freq = F.freq;
}

void Fpmm::setup_trm(int nl, double fq){

  int locus,m, f, o,m1, f1, o1;
  int ndimt;
  Vector<double> gf(3); // genotypic freqs for 1 locus
  double sum;
  Dblock trmt;  // temporary transmission matrix; used to update trm
  Dblock trm1;  // transmission matrix for one locus
 
  freq = 1 - fq; // to be consistent with SALP freq is the frequency of the bad allele
  ndim=2*nl+1;

  trm.resize(ndim,ndim,ndim);
  genfreq.resize(ndim);
  if (ndim == 1) {
    // We have no polygenic component so ndim=ncol=nrow=1 so that
    // trm is 1 and genfreq=1
    trm[0][0][0]=1.0;
    genfreq[0]=1.0;
  }
  else {
  trmt.resize(ndim,ndim,ndim);
  trm1.resize(3,3,3);
 
  gf(1) = (1-freq)*(1-freq);
  gf(2) = 2*freq*(1-freq);
  gf(3) = freq*freq;

  // transition probs; father, mother, child
  double A[3][3][3] = 
{1.0, 0.0, 0.0,   0.5,  0.5, 0.0,    0.0, 1.0, 0.0,
 0.5, 0.5, 0.0,   0.25, 0.5, 0.25,   0.0, 0.5, 0.5,
 0.0, 1.0, 0.0,   0.0,  0.5, 0.5,    0.0,0.0,1.0};

// calculating trm1

 for (f=0;f<3;f++){
    for (m=0;m<3;m++){
       for (o=0;o<3;o++){
	 trm1[f][m][o] = A[f][m][o]*gf[f]*gf[m];
         trmt[f][m][o] = trm1[f][m][o];
         trm[f][m][o] = trm1[f][m][o];
       }
    }
 }

 for (locus=2; locus<=nl; locus++){
   // inittialize trm to zero
   ndimt = 2*locus + 1;
    for (f=0;f<ndimt;f++){
      for (m=0;m<ndimt;m++){
       for (o=0;o<ndimt;o++){
	 trm[f][m][o] = 0.0;
       }
      }
    }
   // go through all combinations of trmt
    ndimt = 2*(locus-1) + 1;
    for (f=0;f<ndimt;f++){
      for (m=0;m<ndimt;m++){
       for (o=0;o<ndimt;o++){
       // make appropriate contributions for each comb. of trm1
	 for (f1=0;f1<3;f1++){
	   for (m1=0;m1<3;m1++){
	     for (o1=0;o1<3;o1++){
	       trm[f+f1][m+m1][o+o1] += trmt[f][m][o]*trm1[f1][m1][o1];
	     }
	   }
	 }
       }
      }
    }
    ndimt = 2*locus + 1;
    for (f=0;f<ndimt;f++){
      for (m=0;m<ndimt;m++){
       for (o=0;o<ndimt;o++){
	 trmt[f][m][o] = trm[f][m][o];
       }
      }
    }
 }
    for (f=0;f<ndim;f++){
      genfreq[f] = 0.0;
      for (m=0;m<ndim;m++){
        sum = 0.0;
	for (o=0;o<ndim;o++){ 
	  sum += trm[f][m][o];
	}
	genfreq[f] += sum;
	for (o=0;o<ndim;o++){ 
	  trm[f][m][o] /= sum;
	}
      }
    }
  }
  //  std::cout << genfreq;

}

void Fpmm::setup_genval(double vr){
  var = vr;  
  int gen;
  double mean=0.0, vart=0.0,scale;

  genval.resize(ndim);
  if (ndim == 1) {
    // no polygenic component
    genval[0]=0.0;
  }
  else {
    for (gen=0; gen<ndim; gen++){
      mean += gen*genfreq[gen];
    }
    for (gen=0; gen<ndim; gen++){
      vart += (gen - mean)*(gen - mean)*genfreq[gen];
    }
    scale = std::sqrt(var/vart);
    for (gen=0; gen<ndim; gen++){
      genval[gen] = (gen - mean)*scale;
    }
  }
}

int Fpmm::display(void){
  std::cout << "Fpmm Object: Polygenic Values" << std::endl;
  std::cout << genval << std::endl;
  return 1;
}

} ///////// end of namespace matvec

