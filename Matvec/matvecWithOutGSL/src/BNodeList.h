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
 
#ifndef BNL_H
#define BNL_H

#include <iostream>
#include "safe_vectors.h"

namespace matvec {
	
	
	class BlockNode {
	public:
		GNode *myNode;
		CutSet *generatedSet;
		CutSet *neighborSet;
		SafeSTLVector<GNodeSet*> vectorGNSets;
		void peel();
		void reverseSampleGNode(void);
	};
	
	class BlockNodeList:public SafeSTLVector<BlockNode> {
	public:
		void initNodes(GNodeList &blockList);
		void displayGNodeSets(void);
		void peelNoCutNoOrder(void);
		void sampleSegregationIndicators(void);
		void reverseSample(void);
		void setGNodeSampleFlags(bool state);
		void saveGNodeStates(void);
		void updateCounters(void);
	};
	
	
}////// end of namespace matvec


#endif
