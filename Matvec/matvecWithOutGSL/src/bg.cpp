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

#include <iostream>
#include <fstream>
#include "bg.h"

namespace matvec{

  unsigned BG::dimen=3;
  matvec::doubleMatrix BG::x;

  BG::BG(void){
    if(dimen==0){
      std::cout <<"BG has not been initialized yet 1\n";
      exit (1);
    } 
    else {
      f = 0.0;
      d.resize(dimen,1,0.0);
      D.resize(dimen,dimen,0.0);
    }
  }

  BG::BG(unsigned i, double b){
    if(dimen==0){
      std::cout <<"BG has not been initialized yet 2\n";
      exit (1);
    }
    else {
      f = x[i][0]*b;
      d.resize(dimen,1,0.0);
      d[i][0] = b;
      D.resize(dimen,dimen,0.0);
    }
  }

  BG::BG(double a){
    if(dimen==0){
      std::cout <<"BG has not been initialized yet 3\n";
      exit (1);
    }
    else {
      f = a;
      d.resize(dimen,1,0.0);
      D.resize(dimen,dimen,0.0);
    }
  }

  BG::BG(const BG &a){
    if(dimen==0){
      std::cout <<"BG has not been initialized yet 4\n";
      exit (1);
    }
    else {
      f = a.f;
      d = a.d;
      D = a.D;
    }
  }
    
    

  BG operator+(const BG &a, const BG &b){
    BG r;
    r.f = a.f + b.f;
    r.d = a.d + b.d;
    r.D = a.D + b.D;
    return r;
  }

  BG operator-(const BG &a, const BG &b){
    BG r;
    r.f = a.f - b.f;
    r.d = a.d - b.d;
    r.D = a.D - b.D;
    return r;
  }

  BG operator*(const BG &a, const BG &b){
    BG r;
    r.f = a.f * b.f;
    r.d = a.f*b.d + a.d*b.f;
    r.D = a.f*b.D + a.d*b.d.transpose() + b.d*a.d.transpose() + a.D*b.f;
    return r;
  }

  BG operator/(double a, const BG &b){
    BG r;
    r.f = a/b.f;
    r.d = -a/(b.f*b.f) * b.d;
    r.D = -a*( 1.0/(b.f*b.f) * b.D - 2.0/(b.f*b.f*b.f)*b.d*b.d.transpose() );
    return r;
  }

  BG operator-(const BG &a){
    BG r;
    r.f = -a.f;
    r.d = -a.d;
    r.D = -a.D;
    return r;
  }

  BG operator/(const BG &a, const BG &b){
    BG r = a * (1.0/b);
    return r;
  }

  BG sqrt(const BG &a){
    BG r;
    r.f = std::sqrt(a.f);
    r.d = 0.5/std::sqrt(a.f)*a.d;
    r.D = 0.5*( 1.0/std::sqrt(a.f)*a.D - 0.5/std::sqrt(a.f*a.f*a.f) * a.d*a.d.transpose() );
    return r;
  }

  BG log(const BG &a){
    BG r; 
    r.f = std::log(a.f);
    r.d = 1.0/a.f*a.d;
    r.D = 1.0/a.f*a.D - 1.0/(a.f*a.f) * a.d*a.d.transpose();
    return r;
  }

  BG exp(const BG &a){
    BG r; 
    r.f = std::exp(a.f);
    r.d = r.f*a.d;
    r.D = r.f*a.D + r.d*a.d.transpose();
    return r;
  }

  BG operator+(double a, const BG &b){
    BG r;
    r.f = a + b.f;
    r.d = b.d;
    r.D = b.D;
    return r;
  }

  BG operator-(double a, const BG &b){
    BG r;
    r.f = a - b.f;
    r.d = -b.d;
    r.D = -b.D;
    return r;
  }

  BG operator*(double a, const BG &b){
    BG r;
    r.f = a*b.f;
    r.d = a*b.d;
    r.D = a*b.D;
    return r;
  }

  BG operator+(const BG &a, double b){
    BG r = b + a;
    return r;
  }

  BG operator-(const BG &a, double b){
    BG r;
    r.f = a.f - b;
    r.d = a.d;
    r.D = a.D;
    return r;
  }

  BG operator*(const BG &a, double b){
    BG r = b * a;
    return r;
  }

  BG operator/(const BG &a, double b){
    BG r = a * (1.0/b);
    return r;
  }  

  BG& operator+=(BG &a, const BG &b){
    a = a + b;
    return a;
  }

  BG& operator-=(BG &a, const BG &b){
    a = a - b;
    return a;
  }

  BG& operator*=(BG &a, const BG &b){
    a = a*b;
    return a;
  }

  BG& operator/=(BG &a, const BG &b){
    a = a/b;
    return a;
  }

  BG& operator+=(BG &a, double b){
    a = a + b;
    return a;
  }

  BG& operator-=(BG &a, double b){
    a = a - b;
    return a;
  }

  BG& operator*=(BG &a, double b){
    a = a * b;
    return a;
  }

  BG& operator/=(BG &a, double b){
    a = a / b;
    return a;
  }

  BG& BG::operator=(double b){
    f = b;
    d.resize(BG::dimen,1,0.0);
    D.resize(BG::dimen,BG::dimen,0.0);
    return (*this);
  }

//   BG& BG::operator=(BG& a){
//     f = a.f;
//     d = a.d;
//     D = a.D;
//     return (*this);
//   }

//   BG& BG::operator=(BG a){
//     f = a.f;
//     d = a.d;
//     D = a.D;
//     return (*this);
//   }

  BG& BG::initialize(unsigned i, double b){
    f = b*x[i][0];
    d.resize(BG::dimen,1,0.0);
    d[i][0] = b;
    D.resize(BG::dimen,BG::dimen,0.0);
    return (*this);
  }


  void BG::display(void){
    std::cout <<"Function value is: " << f << std::endl << std::endl;
    std::cout <<"First Derivatives are: " << std::endl;
    std::cout << d << std::endl;
    std::cout <<"Second Derivatives are: " << std::endl;
    std::cout << D << std::endl;
  }
  void BG::NRupdate(void){
    x -= D.ginv0()*d;
  }


} /// end of namespace matvec
