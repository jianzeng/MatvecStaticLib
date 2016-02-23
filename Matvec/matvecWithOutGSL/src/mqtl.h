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

#ifndef MQTL_H
#define MQTL_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "rpedigree.h"
#include "rutil.h"

using namespace std;

namespace matvec{
  class QNode {
  public:
    unsigned mLevel, pLevel;
    double   mPDQ,   pPDQ, f;
	double prQ11, prQ12, prQ21, prQ22, varvi, varvj, covvivj;
	double matQ2Prob, patQ2Prob;
    
    QNode(){ 
      mLevel = 0;
      pLevel = 0;
	  varvi = 0.0;
	  varvj = 0.0;
      f=-1;
    }
  };
  
  class MQTL:public SafeSTLVector<QNode*> {
  public:
	string   mLevelName,  pLevelName;
	string   probName;
    unsigned count;	
    MQTL(RPedigree& P){
		myPedPtr = &P; 
		MQTLStart = -1;
		QTLProbStart = -1;
		pseudoLevel = 0;
	}
    static RPedigree* myPedPtr;
    SparseCij SpCij;
    string MQTLNames;
    matvec::Vector<double> regCoeff;
    Recoder<string> myRecoder;
    double variance;
    int MQTLStart;
	int QTLProbStart;
	unsigned pseudoLevel;
    
    void inputPDQ(char *filename);
	void inputQTLProbs(char *filename);
    void generateMQTLLevels(void);
    void codeMQTL(unsigned i);
	void generateQTLLevels(void);
    void codeQTL(unsigned i);
    void calcMQTLInbreeding(void);
	void calcQTLVars(double mu);
	void calcQTLcov(void);

    double getMQTLrij(unsigned i, unsigned j, unsigned pi, unsigned pj);
	double getQTLcovvivj(unsigned i, unsigned j, unsigned pi, unsigned pj);
    void display(void);
    string getString (unsigned i);
    void codeLevels(void);
    void addGinv(matvec::doubleMatrix& lhs, unsigned startRow, unsigned statCol, double ratio);
	void addGinvDis(matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol);
    void addtoGinv(unsigned pos[], double q[], double v,
		   matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol, 
		   double ratio);
	void addtoGinv(unsigned pos[], double q[], double v,
		   matvec::doubleMatrix& lhs, unsigned startRow, unsigned startCol);

  };
}
#endif
