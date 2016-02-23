/*
 *  RSamplerParms.h
 *  Used to input and store the parameters for RSampler
 *
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

#ifndef MATVEC_RSAMPLER_PARMS_H
#define MATVEC_RSAMPLER_PARMS_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


namespace matvec {
	class RSamplerParms {
		public:
		
		unsigned pedBlockSize;
		unsigned numLoci;
		unsigned numOfSamples;
		unsigned numOfBurnIn;
		unsigned maxCutsetSize;
		unsigned printFlag;
		unsigned startLocus;
		
		std::string resultsFile;
		std::string samplerType;
		std::string whatToCompute;
		std::string samplerUsed;
		std::string howToSample;
		std::string startLocusType;
		
		RSamplerParms(void){
			pedBlockSize = 0;
			numLoci = 1;
			numOfSamples = 1000;
			numOfBurnIn  = 0;
			maxCutsetSize = 16384;
			printFlag = 0;
			startLocus = 0;
			resultsFile = "";
			samplerType = "genotypic";
			samplerUsed = "MH";
			whatToCompute = "genotypeFreq";
			howToSample = "single";
			startLocusType = "fixed";
		}
		void display(std::ostream& os = std::cout);
	};				
}  ///////// end of namespace matvec
#endif 		
	
	
