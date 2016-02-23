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

#ifndef MATVEC_PARM_MAP_H
#define MATVEC_PARM_MAP_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <exception>


namespace matvec {
	class ParmMap: public std::map<std::string,std::string> {
public:
		void inputParms(std::string fileName);
		double getDoubleValue(std::string parmName);
		unsigned getUnsignedValue(std::string parmName);
		std::string getStringValue(std::string parmName);
		char* getCharPtr(std::string parmName);
		void display(void);
		void display(std::string fileName);		
	};
			
}  ///////// end of namespace matvec
#endif 				
