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

#include "parmMap.h"
#include "exception.h"
#include "rutil.h"

namespace matvec {
	
//	void ParmMap::inputParms(std::string fileName){
//		std::ifstream inFile(fileName.c_str());
//		if(!inFile) {
//			std::cerr << "Couldn't open parmFile: " <<fileName << std::endl;
//			exit (-1);
//		}
//
//		std:: string parmName, parmValue;
//		ParmMap::iterator mapit;
//		while (inFile >> parmName >> parmValue){
//			if (parmName == "//") continue; // comment
//			mapit = find(parmName);
//			if(mapit== end()){
//				(*this)[parmName]=parmValue;
//			}
//			else {
//				mapit->second = parmValue;
//			}
//		} 
//	}
//	
//	void ParmMap::inputParms(std::string fileName){
//		std::ifstream inFile(fileName.c_str());
//		if(!inFile) {
//			std::cerr << "Couldn't open parmFile: " <<fileName << std::endl;
//			exit (-1);
//		}
//        SafeSTLVector<std::string> colData;
//		std:: string parmName, parmValue, inputStr,sep;
//		unsigned commentCount = 0;
//		sep = " ,\t";
//		ParmMap::iterator mapit; 
//		while (getline(inFile,inputStr)){
//			string::size_type begidx,endidx;
//			begidx = inputStr.find_first_not_of(sep);
//			colData.clear();
//			while (begidx != string::npos) {
//				endidx = inputStr.find_first_of(sep,begidx);
//				if (endidx == string::npos) endidx = inputStr.length();
//				colData.push_back(inputStr.substr(begidx,endidx - begidx));
//				begidx = inputStr.find_first_not_of(sep,endidx);
//			}
//			
//			if (colData.size()<2) continue;
//			parmName = colData[0];
//			if (parmName == "//") {
//			   std::ostringstream S;
//			   commentCount++;
//				S << commentCount;
//				std::string countStr;
//				countStr = S.str();
//			    parmName = "H"+countStr+"  ";
//				parmValue = inputStr;
//			}
//			else {
//				parmValue = colData[1];
//			}
//			(*this)[parmName] = parmValue;
//		}
//	}
//
	void ParmMap::inputParms(std::string fileName){
		std::ifstream inFile(fileName.c_str());
		if(!inFile) {
			std::cerr << "Couldn't open parmFile: " <<fileName << std::endl;
			exit (-1);
		}
//        SafeSTLVector<std::string> colData;
		std:: string parmName, parmValue, inputStr,sep;
		unsigned commentCount = 0;
		sep = " ,\t";
		Tokenizer colData;
		ParmMap::iterator mapit; 
		while (getline(inFile,inputStr)){
			colData.getTokens(inputStr,sep);
			if (colData.size()<2) continue;
			parmName = colData[0];
			if (parmName == "//") {
			   std::ostringstream S;
			   commentCount++;
				S << commentCount;
				std::string countStr;
				countStr = S.str();
			    parmName = "H"+countStr+"  ";
				parmValue = inputStr;
			}
			else {
				parmValue = colData[1];
			}
			(*this)[parmName] = parmValue;
		}
	}


	
	void ParmMap::display(void){
		ParmMap::iterator mapit;
		for(mapit=this->begin(); mapit!= this->end(); mapit++){
			std::cout<<mapit->first<<"\t"<<mapit->second<<std::endl;
		}
	}
	
	void ParmMap::display(std::string fileName){
		std::ofstream outFile;
		outFile.open(fileName.c_str());
		if(!outFile) {
			std::cerr << "In matvec::parmMap::display: Couldn't open out file: " << fileName << std::endl;
			exit (-1);
		}
		ParmMap::iterator mapit;
		for(mapit=this->begin(); mapit!= this->end(); mapit++){
			outFile << setw(20) << setiosflags (ios::right | ios::fixed) <<mapit->first<<" ";
			outFile << setw(10) << setiosflags (ios::right | ios::fixed) <<mapit->second<<std::endl;
		}
	}

	
	double ParmMap::getDoubleValue(std::string parmName){
		ParmMap::iterator mapit = find(parmName);
		if(mapit== end()){
			std::cerr << "Unrecognized parameter: " << parmName << std::endl;
			throw exception("ParmMap::getDoubleValue: Error\n");
		}
		std::string parmValue = mapit->second;
		std::istringstream istr(parmValue);
		double doubleValue;
		istr >> doubleValue;
		return doubleValue;
	}
	
	unsigned ParmMap::getUnsignedValue(std::string parmName){
		ParmMap::iterator mapit = find(parmName);
		if(mapit== end()){
			std::cerr << "Unrecognized parameter: " << parmName << std::endl;
			throw exception("ParmMap::getUnsignedValue: Error\n");;
		}
		std::string parmValue = mapit->second;
		std::istringstream istr(parmValue);
		unsigned unsignedValue;
		istr >> unsignedValue;
		return unsignedValue;
	}
	
	std::string ParmMap::getStringValue(std::string parmName){
		ParmMap::iterator mapit = find(parmName);
		if(mapit== end()){
			std::cerr << "Unrecognized parameter: " << parmName << std::endl;
			throw exception("ParmMap::getStringValue: Error\n");;
		}
		std::string parmValue = mapit->second;
		return parmValue;
	}
	
	
	char* ParmMap::getCharPtr(std::string parmName){
		ParmMap::iterator mapit = find(parmName);
		if(mapit== end()){
			std::cerr << "Unrecognized parameter: " << parmName << std::endl;
			throw exception("ParmMap::getCharPtr: Error\n");;
		}
		char* parmValue = (char*)(mapit->second).c_str();
		return parmValue;
	}


}  ///////// end of namespace matvec
