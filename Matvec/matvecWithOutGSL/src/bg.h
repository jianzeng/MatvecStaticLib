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

#include "doublematrix.h"
#include <cmath>
#include <iostream>
#ifndef BG_H
#define BG_H

namespace matvec {


  class BG{
  public:
    static unsigned       dimen;
    static doubleMatrix   x;
    double                f;
    doubleMatrix          d;
    doubleMatrix          D;

    BG();
    BG(unsigned x_i, double b);
    BG(double a);
    BG(const BG &a);
    
    static void set_dimen(unsigned i){
      dimen = i;
      x.resize(dimen,1,0.0);
    }
    BG& initialize(unsigned i, double b);

    friend BG operator+(const BG &a, const BG &b);
    friend BG operator-(const BG &a, const BG &b);
    friend BG operator*(const BG &a, const BG &b);
    friend BG operator/(const BG &a, const BG &b); 
    friend BG operator-(const BG &a); 
    friend BG operator/(double a, const BG &b);
    friend BG operator+(double a, const BG &b);
    friend BG operator-(double a, const BG &b);
    friend BG operator*(double a, const BG &b);

    friend BG   operator+ (const BG &a, double b);
    friend BG   operator- (const BG &a, double b);
    friend BG   operator* (const BG &a, double b);
    friend BG   operator/ (const BG &a, double b);
    
    friend bool operator<=(const BG &a, double b){return a.f<=b;};
    friend bool operator< (const BG &a, double b){return a.f< b;};
    friend bool operator>=(const BG &a, double b){return a.f>=b;};
    friend bool operator> (const BG &a, double b){return a.f> b;};
    friend bool operator==(const BG &a, double b){return a.f==b;};
    friend bool operator!=(const BG &a, double b){return a.f!=b;};
				   
    friend bool operator<=(const BG &a, const BG &b){return a.f<=b.f;};
    friend bool operator< (const BG &a, const BG &b){return a.f< b.f;};
    friend bool operator>=(const BG &a, const BG &b){return a.f>=b.f;};
    friend bool operator> (const BG &a, const BG &b){return a.f> b.f;};
    friend bool operator==(const BG &a, const BG &b){return a.f==b.f;};
    friend bool operator!=(const BG &a, const BG &b){return a.f!=b.f;};

    friend bool operator<=(double b, const BG &a){return b<=a.f;};
    friend bool operator< (double b, const BG &a){return b< a.f;};
    friend bool operator>=(double b, const BG &a){return b>=a.f;};
    friend bool operator> (double b, const BG &a){return b> a.f;};
    friend bool operator==(double b, const BG &a){return b==a.f;};
    friend bool operator!=(double b, const BG &a){return b!=a.f;};


    friend BG& operator+=(BG &a, const BG &b);
    friend BG& operator-=(BG &a, const BG &b);
    friend BG& operator*=(BG &a, const BG &b);
    friend BG& operator/=(BG &a, const BG &b);
    friend BG& operator+=(BG &a, double b);
    friend BG& operator-=(BG &a, double b);
    friend BG& operator*=(BG &a, double b);
    friend BG& operator/=(BG &a, double b);
    

    BG& operator=(double b);
    //BG& operator=(BG& a);
    //BG& operator=(BG a);

    friend BG sqrt(const BG &a);
    friend BG log (const BG &a);
    friend BG exp (const BG &a);
    BG power(double y) const {BG r = exp(y*log(*this)); return r;}

    friend double fabs(const BG &a){return std::fabs(a.f);};

    void NRupdate(void);

    void display(void);

    friend std::ostream& operator<<(std::ostream &os, const BG &a){
      os << a.f ; 
      return os;
    }


  };

 inline void initial_BGarray( BG *v, unsigned n, double a) {
  for (unsigned i = 0; i<n; i++) {
     v[i] = a;
  }
}

}  ///////// end of namespace matvec
#endif  
