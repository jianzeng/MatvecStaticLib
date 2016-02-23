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
#include <iomanip>
#include <fstream>
#include "histogram.h"
#include "safe_vectors.h"
#include "plotter.h"
#include "doublematrix.h"

namespace matvec {
	
	void histogram::makeHistogram(unsigned n){
		
		// Authors: Rohan L. Fernando 
		// (September, 2003) 
		// Contributors: 
		// find min and max values
		
		histogram::iterator hiter;
		float max_val = -1e200;
		float min_val = +1e200;
		for (hiter=begin();hiter!=end();hiter++){
			if(*hiter >= max_val) max_val = *hiter;
			if(*hiter <= min_val) min_val = *hiter;
		}
		float diff = max_val - min_val;
		Vector<int> y;
		y.resize(n,0);
		Vector<double>   x;
		x.resize(n,0.0);
		for (hiter=begin();hiter!=end();hiter++){
			unsigned i = (unsigned)( (*hiter - min_val)/diff * n);
			i = (i>n-1) ? n-1 : i;
			y[i]++;
		}
		data.resize(n,2);
		float incr = diff/n;
		float start = min_val + incr/2.0;
		for (unsigned i=0;i<n;i++){
			data[i][0] = start + incr*i;
			data[i][1] = y[i]/double(size());
		}
	}
	
	void histogram::plot(unsigned n,string flag){
		// Authors: Rohan L. Fernando 
		// (September, 2003) 
		// Contributors: 
		
		makeHistogram(n);
		if(flag=="plot"){
			P.plot(data,matvec::Plotter::plt,"with impulses");
		}
		else {
			P.plot(data,matvec::Plotter::replt,"with impulses");
		}
	}
    
} //end of namespace matvec

