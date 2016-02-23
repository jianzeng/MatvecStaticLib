/*
 *  RSamplerParms.cpp
 *  Used to input parameters for the RSampler
 *  Created by Rohan Fernando on 10/24/04.
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 
 *   You should have received a copy of the GNU Library General Public
 *   License along with this library; if not, write to the Free
 *   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 *   MA 02111-1307, USA 
 *
 */

#include "rsamplerparms.h"

void matvec::RSamplerParms::display(std::ostream &os){
	os << "sampler type:        " << samplerType    << std::endl;
	os << "sampler used:        " << samplerUsed    << std::endl;
	os << "computed quantity:   " << whatToCompute  << std::endl;
	os << "how sampled:         " << howToSample    << std::endl;
	os << "number of samples:   " << numOfSamples   << std::endl;
	os << "number for burnin:   " << numOfBurnIn    << std::endl;
	os << "maximum cutset size: " << maxCutsetSize  << std::endl;
	os << "number of loci:      " << numLoci        << std::endl;
	os << "start from locus:    " << startLocus     << std::endl;
	os << "start Locus Type:    " << startLocusType << std::endl;
}
