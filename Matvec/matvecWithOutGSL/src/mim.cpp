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

#include "mim.h"
#include "parmMap.h"
#include "rsamplerparms.h"
#include "population.h"
#include "stat.h"
#include "BNodeList.h"

using namespace std;

namespace matvec{
	
	MIM* MIModelTerm::myMIMPtr;
	string MIM::solMethod;
	Vector<double> *MIM::vec, MIM::diag, MIM::res;
	ParmMap MIM::parmMap;
	doubleMatrix MIM::IBDCovMatrix, MIM::meansVector;
	IBDMatrix MIM::ibdMatrix;
	
	
	void MIM::initialize(void)
	{
		pop = 0;
		numTerms = 0;
	}
	
	void MIM::release(void)
	{
		if (pop) {
			delete pop; 
			pop = 0; 
		}
		for (unsigned i=0;i<numTerms;i++){
			if (modelTrmVec[i].myRecoderPtr){
				if (!modelTrmVec[i].myMQTLPtr  && !modelTrmVec[i].pedTerm) { // Recoder for MQTL objects don't have to be deleted
					if(modelTrmVec[i].secondMQTLEffect == false) { // the Recoder Object for the second MQTLallele is pointing to the same recoder as the first MQTLallele, this is avoiding the error of releasing it twice!
						delete modelTrmVec[i].myRecoderPtr; 
					}
					modelTrmVec[i].myRecoderPtr = 0;
				}
			}
		}
	}
	
	void MIModelTerm::putFactors(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		Tokenizer tokens;
		string sep("*");
		tokens.getTokens(str,sep);
		factors.clear();
		for (unsigned i=0;i<tokens.size();i++){
			unsigned factorIndex = myMIMPtr->colName.getIndex(tokens[i]);
			if (factorIndex == -1){
				cerr <<"Independent Variable " 
				<< tokens[i]
				<< " not in list of column names \n";
				exit(-1);
			}
			else {
				factors.push_back(factorIndex);
			}
		}
	}
	
	string MIModelTerm::getTermString(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		unsigned numFactors = factors.size();
		string trmStr;
		unsigned factorIndex = factors[0];
		if(myMIMPtr->colType[factorIndex]=="COV"){
			trmStr = myMIMPtr->colName[factorIndex];
		}
		else {
			trmStr = myMIMPtr->colData[factorIndex];
		}
		for (unsigned i=1;i<numFactors;i++){
			factorIndex = factors[i];
			if(myMIMPtr->colType[factorIndex]=="COV"){
				trmStr += "*" + myMIMPtr->colName[factorIndex];
			}
			else{
				trmStr += "*" + myMIMPtr->colData[factorIndex];
			}
		}
		return trmStr;
	}
	
	double MIModelTerm::getTermValue(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		unsigned numFactors = factors.size();
		double value = 1.0;  
		for (unsigned i=0;i<numFactors;i++){
			unsigned factorIndex = factors[i];
			if (myMIMPtr->colType[factorIndex]=="COV"){
				string covStr = myMIMPtr->colData[factorIndex];
				value *= MIM::getDouble(covStr);
			}
		}
		return value;
	}
	
	void MIM::putColNames(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" \t");
		str = "intercept " + str;
		colName.getTokens(str,sep);
		numCols = colName.size();
	}
	
	void MIM::putColTypes(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" \t");
		str = "CLASS " + str; 
		colType.getTokens(str,sep);
		if (numCols!=colType.size()){
			cerr <<"number of column names and column types do not match\n";
			exit (-1);
		}
		unsigned n = 0;
		for (unsigned i=0;i<numCols;i++){
			if (colType[i] == "DEP") {
				depVar.push_back(colName[i]);
				n++;
			}
		}
		numTraits = n;
	}
	
	void MIM::putMarkerColNames(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" \t");
		markerColName.getTokens(str,sep);
		markerNumCols = markerColName.size();
	}
	
	void MIM::putModels(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		MIModelTerm::myMIMPtr = this;
		Tokenizer models;
		string sep =";";
		models.getTokens(str,sep);
		for (unsigned i=0; i<models.size();i++){
			putModel(models[i]);
		}
		for(unsigned i=0;i<MQTLVec.size();i++){
			putMQTLStuffInModelTerms(*MQTLVec[i]);
		}
	}
	
	void MIM::putModel(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" =+"); 
		Tokenizer modelTokens;
		modelTokens.getTokens(str,sep);
		unsigned nTokens = modelTokens.size();
		int depVarIndex = colName.getIndex(modelTokens[0]);
		if (depVarIndex == -1){
			cerr << "Dependent Variable " 
			<< modelTokens[0] 
			<< " not in list of column names \n";
			exit (-1);
		}
		MIModelTerm modelTrm;
		modelTrm.secondMQTLEffect = false;
		modelTrm.pedTerm = false;
		modelTrm.myMQTLPtr = 0;
		modelTrm.depVarName = modelTokens[0];
		modelTrm.trait = depVar.getIndex(modelTokens[0]);
		for (unsigned i=1;i<nTokens;i++){
			modelTrm.myRecoderPtr = new Recoder<string>;
			modelTrm.name = modelTokens[i];
			modelTrm.putFactors(modelTrm.name);
			modelTrmVec.push_back(modelTrm);
		}
	}
	
	void MIM::isMIMReady(void){
		if (colName.size()==0){
			throw exception("Did you forget to use MIM::putColNames?\n");
		}
		if (colType.size()==0){
			throw exception("Did you forget to use MIM::putColTypes?\n");
		}
		if (modelTrmVec.size()==0){
			throw exception("Did you forget to use MIM::putModels?\n");
		}
	}
	
	void MIM::inputData(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		if(dataVec.size()) return;
		SafeSTLVector<unsigned> multVector, missVector; // for generating missing value index (missIndex) 
		multVector.resize(numTraits);
		missVector.resize(numTraits);
		multVector[0]=1;
		for (unsigned i=1;i<numTraits;i++){
			multVector[i] = multVector[i-1]*2;
		}
		
		unsigned missIndex;
		cout <<"\n processing data file \n";
		MimDataNode dataNode;
		numTerms = modelTrmVec.size();
		dataNode.trmVec.resize(numTerms);
		dataNode.depVec.resize(numTraits);
		dataNode.misVec.resize(numTraits);
		ifstream datafile;
		datafile.open(phenDataFileName.c_str());
		if(!datafile) {
			cerr << "Couldn't open data file: " <<phenDataFileName << endl;
			exit (-1);
		}
		std::string sep(" \t,");
		std::string inputStr;
		unsigned lineNumber = 0;
		double rec = 0, rec1 = 0;
		while (getline(datafile,inputStr)){
			rec++;
			if(rec==1000){
				cout<<rec+rec1<<"\r";
				cout.flush();
				rec1 += rec;
				rec = 0;
			}
		    lineNumber++;
			inputStr = "--- " + inputStr;
			colData.getTokens(inputStr,sep);
			unsigned n = colData.size();
			if (n != numCols){
				cerr << "Line " << lineNumber << " of data file "
				<< phenDataFileName<< "  has "   << n <<" columns; " 
				<< numCols <<" expected " << endl;
				throw exception("Error in MIM::inputData()\n");
			}
			unsigned j=0;
			for (unsigned i=0;i<numCols;i++){
				if (colType[i]=="DEP") {
					if(colData[i]=="."){
						missVector[j] = 0;
						dataNode.misVec[j]   = true;
						dataNode.depVec[j++] = 0.0;
					}
					else {
						missVector[j] = 1;
						dataNode.misVec[j]   = false;
						dataNode.depVec[j++] = getDouble(colData[i]);
					}
				}
				else if (colType[i]=="CLASS" || colType[i]=="COV"){
					if (colData[i]=="."){
						continue;
					}
				}
			}
			missIndex = 0;
			for (unsigned i=0;i<numTraits;i++){
				missIndex += missVector[i]*multVector[i];
			}
			if (missIndex==0) continue;
			map<unsigned,doubleMatrix>::iterator RMapit = RMap.find(missIndex);
			if (RMapit==RMap.end()){
				RMap[missIndex] = getMissR(missVector);
			}
			dataNode.RiPtr = &(RMap[missIndex]);
			
			for (unsigned i=0;i<numTerms;i++){ 
				dataNode.trmVec[i].level = modelTrmVec[i].getTermLevel();
				dataNode.trmVec[i].value = modelTrmVec[i].getTermValue();
				if (modelTrmVec[i].myMQTLPtr){
					unsigned tr = modelTrmVec[i].trait;
					dataNode.trmVec[i].value *= modelTrmVec[i].myMQTLPtr->regCoeff[tr];
				}
			}
			dataVec.push_back(dataNode);
		}
		datafile.close();
	}
	
	doubleMatrix MIM::getMissR(SafeSTLVector<unsigned> missVector){
		doubleMatrix Ri = R;
		for (unsigned k=0;k<numTraits;k++){
			if (missVector[k]==0){
				for (unsigned j=0;j<numTraits;j++){
					Ri[k][j] = Ri[j][k] = 0;
				}
			}
		}
		Ri.ginv1();
		return (Ri);
	}
	
	double MIM::getDouble(string& Str) {
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		istringstream inputStrStream(Str.c_str());
		double val;
		inputStrStream >> val; 
		return val;
	}
	int MIM::getInteger(string& Str) {
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		istringstream inputStrStream(Str.c_str());
		int val;
		inputStrStream >> val; 
		return val;
	}
	
	void MIM::calcStarts(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		unsigned prevI =0,  maxI = 0;
		modelTrmVec[0].start = 0;
		for (unsigned i=1;i<numTerms;i++){
			if (modelTrmVec[i].myMQTLPtr) {
				if(modelTrmVec[i].myMQTLPtr->probName != modelTrmVec[i].name ){
					if(modelTrmVec[i].myMQTLPtr->MQTLStart == -1){
						modelTrmVec[i].start = modelTrmVec[prevI].start + modelTrmVec[prevI].nLevels();
						prevI = i;
						modelTrmVec[i].myMQTLPtr->MQTLStart = modelTrmVec[i].start;
					}
					else {
						modelTrmVec[i].start = modelTrmVec[i].myMQTLPtr->MQTLStart;
					}
				}
				else {
					if(modelTrmVec[i].myMQTLPtr->MQTLStart == -1){			   
						modelTrmVec[i].start = modelTrmVec[prevI].start + modelTrmVec[prevI].nLevels();
						prevI = i;
						modelTrmVec[i].myMQTLPtr->QTLProbStart = modelTrmVec[i].start;
					}
					else {
						modelTrmVec[i].start = modelTrmVec[i].myMQTLPtr->QTLProbStart;
					}
				}
			}
			else {
				modelTrmVec[i].start = modelTrmVec[prevI].start + modelTrmVec[prevI].nLevels();
				prevI = i;
			}
			if (modelTrmVec[i].start > modelTrmVec[maxI].start){
				maxI = i;
			}
		}
		mmeSize = modelTrmVec[maxI].start + modelTrmVec[maxI].nLevels();
	}
	
	void MIM::calcWPW(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		doubleMatrix* Ri;
		unsigned ii,jj,ti,tj;
		double vi,vj;
		bool direct = (solMethod=="direct") ? true : false;
		if(!direct){
			diag.resize(mmeSize,0.0);
			res.resize(mmeSize,0.0);
		}
		rhs.resize(mmeSize,0.0);
		if(printResVar){
			cout <<"Residual Cov Matrix = \n" << R << endl;
			printResVar = false;
		}
		double rec = 0, rec1 = 0;
		if (solMethod=="direct"){
			for (unsigned i=0;i<dataVec.size();i++){
				rec++;
				if(rec==1000){
					cout<< rec1+rec <<"\r";
					cout.flush();
					rec1 += rec;
					rec = 0;
				}		
				Ri = dataVec[i].RiPtr;
				for (unsigned mi=0;mi<numTerms;mi++){
					ii = modelTrmVec[mi].start + dataVec[i].trmVec[mi].level - 1;
					ti =  modelTrmVec[mi].trait;
					vi = dataVec[i].trmVec[mi].value;
					for (unsigned k=0;k<numTraits;k++){
						rhs[ii] += vi*(*Ri)[ti][k]*dataVec[i].depVec[k];
					}
					for (unsigned mj=0;mj<numTerms;mj++){
						jj = modelTrmVec[mj].start + dataVec[i].trmVec[mj].level - 1; 
						tj =  modelTrmVec[mj].trait;
						vj = dataVec[i].trmVec[mj].value;
						lhs[ii][jj] += vi*(*Ri)[ti][tj]*vj; 
					}
				}
			}
		}
		else {
			vector<unsigned> xxVec; 
			vector<double>   vxVec; 
			xxVec.resize(numTerms);
			vxVec.resize(numTerms);
			MimDataNode dataNode;
			double cij;
			for (unsigned i=0;i<dataVec.size();i++){
				dataNode = dataVec[i];
				rec++;
				if(rec==1000){
					cout<<rec1+rec <<"\r";					
					cout.flush();
					rec1 += rec;
					rec = 0;
				}
				for (unsigned mx=0;mx<numTerms;mx++){
					xxVec[mx] = modelTrmVec[mx].start + dataNode.trmVec[mx].level - 1;
					vxVec[mx] = dataNode.trmVec[mx].value; 
				}
				Ri = dataVec[i].RiPtr;
				for (unsigned mi=0;mi<numTerms;mi++){ 
					ii  = xxVec[mi];
					ti  =  modelTrmVec[mi].trait;
					vi  = vxVec[mi];
					diag[ii] += vi*(*Ri)[ti][ti]*vi;
					res[ii] += vi*(*Ri)[ti][ti]*vi*(*vec)[ii];
					for (unsigned k=0;k<numTraits;k++){
						rhs[ii] += vi*(*Ri)[ti][k]*dataNode.depVec[k];
					}
					for (unsigned mj=mi+1;mj<numTerms;mj++){
						jj = xxVec[mj]; 
						tj =  modelTrmVec[mj].trait;
						vj = vxVec[mj];
						cij = vi*(*Ri)[ti][tj]*vj;
						res[ii] += cij*(*vec)[jj];
						res[jj] += cij*(*vec)[ii];
					}
				}
			}
		}
		cout << endl;
	}
	void MIM::addGinv(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		for (unsigned i=0;i<covBlockVec.size();i++){
			covBlockVec[i].addGinv();
		}
	}
	
	void MIM::initSetup(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		inputData();
		calcStarts();
		if(solMethod=="direct"){
			lhs.resize(mmeSize,mmeSize,0.0);
		}
		else {
			sol.resize(mmeSize,0.0);
		}
	}
	
	void MIM::getDirectSolution(){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		solMethod = "direct";
		initSetup();	
		calcWPW();
		addGinv();
		cout << "solving \n";	
		sol = lhs.ginv0()*rhs;
	}
	
	void MIM::display(std::string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:  
		if (str=="pop"){
			pop->display("marker");
		}
		if (sol.size()==0) return;
		if (solMethod=="direct"){
			cout << "Diagonals of LHS " << endl;
			for (unsigned i = 0;i<mmeSize;i++){
				for (unsigned j = 0;j<mmeSize;j++){
					cout << setw(8) << setprecision (3) 
					<< setiosflags (ios::right | ios::fixed) << lhs[i][j] <<" ";
				}
				cout << endl;
			}
			cout << "RHS " << endl;
			cout << rhs << endl;
			for (unsigned i=0;i<modelTrmVec.size();i++){
				cout << "Solutions for " << modelTrmVec[i].name 
				<< ", Trait: " << modelTrmVec[i].depVarName << endl;
				Recoder<string>::iterator it;
				for (it=modelTrmVec[i].myRecoderPtr->begin(); it!=modelTrmVec[i].myRecoderPtr->end();it++){
					unsigned ii = modelTrmVec[i].start + it->second - 1;
					cout << setw(10) << it->first << " " << sol[ii] << endl;
				}
			}
		}
		else {
			for (unsigned i=0;i<modelTrmVec.size();i++){
				cout << "Solutions for " << modelTrmVec[i].name 
				<< ", Trait: " << modelTrmVec[i].depVarName << endl;
				cout << setw(20) << " Name "     << " ";
				cout << setw(10) << " Solution " << " ";  
				cout << setw(10) << " N "        << " "; 
				cout << setw(10) << " e/rhs \n"; 
				Recoder<string>::iterator it;
				for (it=modelTrmVec[i].myRecoderPtr->begin(); it!=modelTrmVec[i].myRecoderPtr->end();it++){
					unsigned ii = modelTrmVec[i].start + it->second - 1;
					double standardResid = (rhs[ii]==0) ? resid[ii] : resid[ii]/rhs[ii]; 
					cout << setw(20)                      << setiosflags (ios::right | ios::fixed) << it->first << " "; 
					cout << setw(10) << setprecision (3)  << setiosflags (ios::right | ios::fixed) << sol[ii]   << " "; 
					cout << setw(10) << setprecision (3)  << setiosflags (ios::right | ios::fixed) << diag[ii]   << " "; 
					cout << setw(10) << setprecision (3)  << setiosflags (ios::right | ios::fixed) << standardResid << endl;
				}
			}
		}
	}
	
	void MIM::getJacobiSolution(double p){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		solMethod = "jacobi";	
		initSetup();
		mmeTimes(sol);                   // result goes into res, also creates rhs and diag
		resid = (rhs - res);
		tempSol = resid/diag + sol;
		double diff = resid.sumsq();
		unsigned iter = 0;
		while(diff/mmeSize > .0000000001 && ++iter<2000){
			sol = p*tempSol + (1-p)*sol;
			cout <<"Iteration : " << iter <<" diff = "<<diff/mmeSize << endl << endl;
			mmeTimes(sol);
			resid = (rhs - res);
			tempSol = resid/diag + sol;
			diff = resid.sumsq();
		}
		cout << resid << endl;
	}
	
	void MIM::getCGSolution(double eps){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		solMethod = "cg";	
		initSetup();
		mmeTimes(sol);                            // result goes into res, also creates rhs and diag
		resid = (res-rhs);
		matvec::Vector<double> d = resid;
		double oldDiffSq;
		double newDiffSq = resid.sumsq();
		unsigned iter = 0;
		while(newDiffSq/mmeSize > eps && ++iter<2*mmeSize){
			cout <<"Iteration : " << iter <<" diff = "<<newDiffSq/mmeSize << endl << endl;
			mmeTimes(d);
			double alpha = -newDiffSq/(res*d).sum();
			sol += alpha*d;
			if (iter%10){
				resid += alpha*res;
			}
			else {
				mmeTimes(sol);
				resid = (res-rhs);
			}
			oldDiffSq = newDiffSq;
			newDiffSq = resid.sumsq();
			double beta = -newDiffSq/oldDiffSq;
			d = resid - beta*d;
		}
		cout << resid << endl;
	}
	
	void MIM::getPCCGSolution(double eps){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		cout << "\n PCCG...\n";
		solMethod = "pccg";	
		initSetup();
		cout <<"\n start of iteration \n";
		mmeTimes(sol); // result goes into res, also creates rhs and diag
		resid = (res - rhs);
		matvec::Vector<double> d = resid/diag;
		double oldDiffSq;
		double newDiffSq = (resid*d).sum();
		unsigned iter = 0;
		while(newDiffSq/mmeSize > eps && ++iter<2*mmeSize){
			clock_t startTime = clock();
			cout <<"Iteration : " << iter <<" diff = "<<newDiffSq/mmeSize << endl << endl;
			mmeTimes(d);
			double alpha = -newDiffSq/(res*d).sum();
			sol += alpha*d;
			if (iter % 25){
				resid += alpha*res;
			}
			else {
				mmeTimes(sol);
				resid = (res - rhs);
			}
			matvec::Vector<double> s = resid/diag;
			oldDiffSq = newDiffSq;
			newDiffSq = (resid*s).sum();
			double beta = -newDiffSq/oldDiffSq;
			d = s - beta*d;
			clock_t endTime = clock();
			clock_t compTime = endTime - startTime; 
			cout << "\n computing time in seconds = " << compTime/CLOCKS_PER_SEC << endl <<endl;
		}
	}
	
	void MIM::mmeTimes(matvec::Vector<double>& x){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		vec = &x;
		calcWPW();
		addGinv();
		for (unsigned i=0;i<mmeSize;i++){
			diag[i] = (diag[i]==0) ? 1.0 : diag[i];
		}
	}
	
	void MIM::mergeMQTLLevelsBy(string mergeStr){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:	
		unsigned numColsInFile = numCols - 1;
		int indexMergeStr = colName.getIndex(mergeStr) - 1;
		if(indexMergeStr == -1){
			cerr << mergeStr <<" not in Column Names \n";
			exit(-1);
		}
		ifstream datafile;
		ofstream outfile;
		string outFileName = phenDataFileName + ".mrgMQTLQ2prob";
		datafile.open(phenDataFileName.c_str());
		if(!datafile) {
			cerr << "Couldn't open data file: " << phenDataFileName
			<< endl;
			exit (-1);
		}
		outfile.open(outFileName.c_str());
		if(!outfile) {
			
			cerr << "Couldn't open file: " << outFileName << endl;
			exit (-1);
		}
		size_t linewidth = 1024;
		char *line = new char [linewidth];
		string sep(" \t,");
		unsigned lineNumber=0;
		while (datafile.getline(line,linewidth)){
			lineNumber++;
			string inputStr(line);
			colData.getTokens(inputStr,sep);
			unsigned n = colData.size();
			if (n != numColsInFile){
				cerr << "Line " << lineNumber << " of data file has " << n <<" columns" << endl;
				cerr << numColsInFile <<" expected " << endl;
				exit(-1);
			}
			string mergeVar = colData[indexMergeStr];
			int mergeID = MQTLVec[0]->myPedPtr->coder.code(mergeVar);
			if (mergeID>MQTLVec[0]->myPedPtr->size()){
				cerr << "Line " << lineNumber << " of data file has merge string "<<mergeVar <<endl;
				cerr << "which is not in Pedigree file\n"; 
			}
			outfile << inputStr << " ";
			for (unsigned i=0;i<MQTLVec.size();i++){
				outfile << (*MQTLVec[i])[mergeID-1]->matQ2Prob + (*MQTLVec[i])[mergeID-1]->patQ2Prob  << " "
				<< (*MQTLVec[i])[mergeID-1]->pLevel << " "  
				<< (*MQTLVec[i])[mergeID-1]->mLevel << " ";
			}
			outfile << endl;
		}
		datafile.close();
		outfile.close();
		phenDataFileName = outFileName;
	}
//	
//	This is the version for equilibrium
//	
//	
//	// Contributors:
//	unsigned numColsInFile = numCols - 1;
//	int indexMergeStr = colName.getIndex(mergeStr);
//	if(indexMergeStr == -1){
//		cerr << mergeStr <<" not in Column Names \n";
//		exit(-1);
//	}
//	indexMergeStr -= 1; // data file does not have a column for the inercept
//	ifstream datafile;
//	ofstream outfile;
//	string outFileName = phenDataFileName + ".mrgMQTL";
//	datafile.open(phenDataFileName.c_str());
//	if(!datafile) {
//		cerr << "Couldn't open data file: " << phenDataFileName << endl;
//		exit (-1);
//	}
//	outfile.open(outFileName.c_str());
//	if(!outfile) {
//		cerr << "Couldn't open file: " << outFileName << endl;
//		exit (-1);
//	}
//	size_t linewidth = 1024;
//	std::string sep(" \t,");
//	string inputStr;
//	unsigned lineNumber=0;
//	while (getline(datafile,inputStr)){
//		lineNumber++;
//		colData.getTokens(inputStr,sep);
//		unsigned n = colData.size();
//		if (n != numColsInFile){
//			cerr << "Line " << lineNumber << " of data file " << phenDataFileName << " has " << n <<" columns" << endl;
//			cerr << numColsInFile <<" expected " << endl;
//			exit(-1);
//		}
//		string mergeVar = colData[indexMergeStr];
//		int mergeID = MQTLVec[0]->myPedPtr->coder.code(mergeVar);
//		if (mergeID>MQTLVec[0]->myPedPtr->size()){
//			cerr << "Line " << lineNumber << " of data file has merge string "<<mergeVar <<endl;
//			cerr << "which is not in Pedigree file\n"; 
//		}
//		outfile << inputStr << " ";
//		for (unsigned i=0;i<MQTLVec.size();i++){
//			outfile << (*MQTLVec[i])[mergeID-1]->pLevel << " "  
//			<< (*MQTLVec[i])[mergeID-1]->mLevel << " ";
//		}
//		outfile << endl;
//	}
//	datafile.close();
//	outfile.close();
//	phenDataFileName = outFileName;
//}
	
	void MIM::addColNames(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" \t");
		Tokenizer tokens;
		tokens.getTokens(str,sep);
		Tokenizer::iterator it;
		for(it=tokens.begin();it!=tokens.end();it++){
			colName.push_back(*it);
		}
		numCols = colName.size();
	}
	
	void MIM::addColTypes(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" \t");
		Tokenizer tokens;
		tokens.getTokens(str,sep);
		Tokenizer::iterator it;
		for(it=tokens.begin();it!=tokens.end();it++){
			colType.push_back(*it);
		}
		if (numCols!=colType.size()){
			cerr <<"number of column names and column types do not match\n";
			exit (-1);
		}
	}
	
	void MIM::putMQTLStuffInModelTerms(MQTL& mQTL){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		Tokenizer names;
		string sep = " ,";
		names.getTokens(mQTL.MQTLNames,sep);
		mQTL.probName    = names[0];
		mQTL.pLevelName  = names[1];
		mQTL.mLevelName  = names[2];
		for (unsigned i=0;i<modelTrmVec.size();i++){
			if(modelTrmVec[i].name==mQTL.pLevelName || modelTrmVec[i].name==mQTL.mLevelName){
				modelTrmVec[i].myMQTLPtr = &mQTL;
				delete modelTrmVec[i].myRecoderPtr; 
				modelTrmVec[i].myRecoderPtr = &mQTL.myRecoder;
			}
			else if(modelTrmVec[i].name==mQTL.probName){
				modelTrmVec[i].myMQTLPtr = &mQTL;
			}
		}
	}	
	
	void CovBlock::buildModelTrmVec(string str){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		string sep(" "); 
		Tokenizer modelTokens;
		modelTokens.getTokens(str,sep);
		unsigned nTokens = modelTokens.size();
		unsigned numModelTrms = MIModelTerm::myMIMPtr->modelTrmVec.size();
		for (unsigned i=0;i<nTokens;i++){
			for (unsigned j=0;j<numModelTrms;j++){
				if (modelTokens[i]==MIModelTerm::myMIMPtr->modelTrmVec[j].name){
					modelTrmPtrVec.push_back(&(MIModelTerm::myMIMPtr->modelTrmVec[j]));
					if (pedPtr){
						delete MIModelTerm::myMIMPtr->modelTrmVec[j].myRecoderPtr;
						MIModelTerm::myMIMPtr->modelTrmVec[j].myRecoderPtr = &pedPtr->coder;
						MIModelTerm::myMIMPtr->modelTrmVec[j].pedTerm = true;
					}
				}
			}
		}
	}
	
	void CovBlock::addGinv(void){
		// Authors: Rohan L. Fernando
		// (2005) 
		// Contributors:
		bool direct = (MIM::solMethod=="direct") ? true : false;
		if(myMQTLPtr){
			myMQTLPtr->addGinvDis(MIModelTerm::myMIMPtr->lhs,
							   myMQTLPtr->MQTLStart,
							   myMQTLPtr->MQTLStart);
			return;
		}
		Vari = Var;
		Vari.inv();
		unsigned n = modelTrmPtrVec.size();
		if (pedPtr && MIM::solMethod!="direct"){
			double q[3], vjj, vii;
			double d,fs,fd;
			unsigned pos[3];
			SafeSTLVector<PNode*>::iterator it;
			if (pedPtr->coder.size() > pedPtr->pedVector.size()){
				cout << "RPedigree is not complete \n";
				exit(-1);
			}
			for (it=pedPtr->pedVector.begin();it!=pedPtr->pedVector.end();it++){
				pos[0] = (*it)->sire;
				pos[1] = (*it)->dam;
				pos[2] = (*it)->ind;
				if((*it)->sire && (*it)->dam){
					q[0] = -0.5;
					q[1] = -0.5;
					q[2] =  1.0;
					fs = pedPtr->pedVector[pos[0]-1]->f;
					fd = pedPtr->pedVector[pos[1]-1]->f;
					d = 4.0/(2 - fs - fd);
				}
				else if((*it)->sire){
					q[0] = -0.5;
					q[1] =  0.0;
					q[2] =  1.0;
					fs = pedPtr->pedVector[pos[0]-1]->f;
					d = 4.0/(3-fs);
				}
				else if((*it)->dam){
					q[0] =  0.0;
					q[1] = -0.5;
					q[2] =  1.0;
					fd = pedPtr->pedVector[pos[1]-1]->f;
					d = 4.0/(3-fd);
				}
				else{
					q[0] =  0.0;
					q[1] =  0.0;
					q[2] =  1.0;
					d    =  1.0;
				}
				for (unsigned c=0;c<n;c++){
					MIModelTerm* mtermcPtr = modelTrmPtrVec[c];
					unsigned startCol = mtermcPtr->start;
					for (unsigned j=0;j<3;j++){
						if(pos[j]) { 
							unsigned jj = startCol + pos[j] - 1;
							MIM::diag[jj] += Vari[c][c]*q[j]*d*q[j];
							vjj = d*q[j]*(*MIM::vec)[jj];
							for (unsigned r=0;r<n;r++){
								double ratio = Vari[r][c];
								MIModelTerm* mtermrPtr = modelTrmPtrVec[r];
								unsigned startRow = mtermrPtr->start;
								for (unsigned i=0;i<3;i++){
									if(pos[i]){
										vii = ratio*q[i];
										unsigned ii = startRow + pos[i] - 1;
										MIM::res[ii] += vii*vjj;
										
									}
								}
							}
						}
					}				
				}
			}
		}
		else if (pedPtr && MIM::solMethod=="direct"){
			double q[3];
			double d,fs,fd;
			unsigned pos[3];
			SafeSTLVector<PNode*>::iterator it;
			if (pedPtr->coder.size() > pedPtr->pedVector.size()){
				cout << "RPedigree is not complete \n";
				exit(-1);
			}
			for (it=pedPtr->pedVector.begin();it!=pedPtr->pedVector.end();it++){
				pos[0] = (*it)->sire;
				pos[1] = (*it)->dam;
				pos[2] = (*it)->ind;
				if((*it)->sire && (*it)->dam){
					q[0] = -0.5;
					q[1] = -0.5;
					q[2] =  1.0;
					fs = pedPtr->pedVector[pos[0]-1]->f;
					fd = pedPtr->pedVector[pos[1]-1]->f;
					d = 4.0/(2 - fs - fd);
				}
				else if((*it)->sire){
					q[0] = -0.5;
					q[1] =  0.0;
					q[2] =  1.0;
					fs = pedPtr->pedVector[pos[0]-1]->f;
					d = 4.0/(3-fs);
				}
				else if((*it)->dam){
					q[0] =  0.0;
					q[1] = -0.5;
					q[2] =  1.0;
					fd = pedPtr->pedVector[pos[1]-1]->f;
					d = 4.0/(3-fd);
				}
				else{
					q[0] =  0.0;
					q[1] =  0.0;
					q[2] =  1.0;
					d    =  1.0;
				}
				for (unsigned r=0;r<n;r++){
					MIModelTerm* mtermrPtr = modelTrmPtrVec[r];
					unsigned startRow = mtermrPtr->start;
					for (unsigned c=0;c<n;c++){
						MIModelTerm* mtermcPtr = modelTrmPtrVec[c];
						unsigned startCol = mtermcPtr->start;
						double ratio = Vari[r][c];
						for (unsigned i=0;i<3;i++){
							if(pos[i]){
								unsigned ii = startRow + pos[i] - 1;
								for (unsigned j=0;j<3;j++){
									if(pos[j]) { 
										unsigned jj = startCol + pos[j] - 1;
										MIModelTerm::myMIMPtr->lhs[ii][jj] += ratio*q[i]*d*q[j];
									}
								}
							}
						}
					}
				} 
			}
		}
		else{
			for (unsigned i=0;i<n;i++){
				MIModelTerm* mtermiPtr = modelTrmPtrVec[i];
				unsigned starti = mtermiPtr->start;
				for (unsigned j=0;j<n;j++){
					MIModelTerm* mtermjPtr = modelTrmPtrVec[j];
					unsigned startj = mtermjPtr->start;
					unsigned numLevels = mtermiPtr->nLevels();
					for (unsigned k=0;k<numLevels;k++){
						unsigned ii = starti + k;
						unsigned jj = startj + k;
						if (direct) {
							MIModelTerm::myMIMPtr->lhs[ii][jj] += Vari[i][j];
						}
						else {
							MIM::res[ii] += Vari[i][j] * (*MIM::vec)[jj];
							if (ii==jj)  MIM::diag[ii] += Vari[i][j];
						}	
					}
				}
			}
		}
	}

	
//	void CovBlock::addGinv(void){
//		// Authors: Rohan L. Fernando
//		// (2005) 
//		// Contributors:
//		bool direct = (MIM::solMethod=="direct") ? true : false;
//		if(myMQTLPtr){
//			myMQTLPtr->addGinv(MIModelTerm::myMIMPtr->lhs,
//							   myMQTLPtr->myStart,
//							   myMQTLPtr->myStart,
//							   1.0/myMQTLPtr->variance);
//			return;
//		}
//		Vari = Var;
//		Vari.inv();
//		unsigned n = modelTrmPtrVec.size();
//		for (unsigned i=0;i<n;i++){
//			MIModelTerm* mtermiPtr = modelTrmPtrVec[i];
//			unsigned starti = mtermiPtr->start;
//			for (unsigned j=0;j<n;j++){
//				MIModelTerm* mtermjPtr = modelTrmPtrVec[j];
//				unsigned startj = mtermjPtr->start;
//				if (pedPtr){
//					pedPtr->addAinv(MIModelTerm::myMIMPtr->lhs,starti,startj,Vari[i][j]);
//				} 
//				else{
//					unsigned numLevels = mtermiPtr->nLevels();
//					for (unsigned k=0;k<numLevels;k++){
//						unsigned ii = starti + k;
//						unsigned jj = startj + k;
//						if (direct) {
//							MIModelTerm::myMIMPtr->lhs[ii][jj] += Vari[i][j];
//						}
//						else {
//							MIM::res[ii] += Vari[i][j] * (*MIM::vec)[jj];
//							if (ii==jj)  MIM::diag[ii] += Vari[i][j];
//						}	
//					}
//				}
//			}	
//		}
//	}
	
	///GNodeSampler methods:
		
	void MIM::runAnalysis(std::string inputFileName,RPedigree& P,GeneticDist& G){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (August, 2005) 
		// Contributors: 
		
		putGNodeSamplerParms(inputFileName);
		setupPopulation(P,G);
		if(monoMendel(G)){       
			transferTraitData();
		}
		if (parmMap["GNodeSampler"]=="GS1"){
			GNodeSamplerGS1();
		}
		else if	(parmMap["GNodeSampler"]=="OBGSL"){
			GNodeSamplerOBGSL();
		}
		else if	(parmMap["GNodeSampler"]=="OBGSL1"){
			GNodeSamplerOBGSL1();
		}
		else if	(parmMap["GNodeSampler"]=="EEGSL"){
			GNodeSamplerEEGSL();
		}
		else if (parmMap["GNodeSampler"]=="JNC"){
			GNodeSamplerJNC();
		}
		else if (parmMap["GNodeSampler"]=="JNCIDG"){
			GNodeSamplerJNCIDG();
		}
		else if (parmMap["GNodeSampler"]=="LG"){
			GNodeSamplerLG();
		}
		else if (parmMap["GNodeSampler"]=="SLNC"){
			GNodeSamplerSLNC();
		}
		else if(parmMap["GNodeSampler"]=="JointNC"){
			GNodeSamplerJointNC();
		}
		else if(parmMap["GNodeSampler"]=="OBGML"){
			GNodeSamplerOBGML();
		}
		unsigned nSamples = parmMap.getUnsignedValue("chainLength");
		pop->outputCounters(nSamples);
	}
		
	void MIM::putGNodeSamplerParms(std::string inputFileName){
	// Authors: Rohan L. Fernando
	// (November, 2004) 
	// Contributors:
	
	// These are the default values 
	parmMap["locusBlockSize"]    = "1";
	parmMap["pedigreeBlockSize"] = "0";
	parmMap["typeOfGNodes"]      = "allelic";
	parmMap["chainLength"]       = "1000";
	parmMap["burnInLength"]      = "0";
	parmMap["maxCutSetSize"]     = "16384";
	parmMap["resultsFile"]       = "";
	parmMap["GNodeSampler"]      = "GS1";
	parmMap["printFlag"]         = "0";
	parmMap["startLocus"]        = "1";
	parmMap["startLocusType"]    = "fixed";
	parmMap["penetranceModel"]   = "";
	parmMap["whatToCompute"]     = "nothing";
	parmMap["numGibbsSteps"]     = "10";
	parmMap["numMHTrials"]       = "1";
	parmMap["blockWidth"]        = "3"; 
	parmMap["prunePedigree"]     = "no";
	parmMap["leftLocus"]         = "1";
	parmMap["rightLocus"]        = "2";
	parmMap["pivotInd"]          = "1";
	parmMap["pedigreeBlockWidth"]   = "3";	
	parmMap["pedigreeBlockOverlap"] = "1";
	parmMap["locusBlockWidth"]   = "3";	
	parmMap["locusBlockOverlap"] = "1";

	
	// Now we read in the specific values for this job
	parmMap.inputParms(inputFileName);
	}
	
	bool MIM::monoMendel(GeneticDist& G){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (August, 2005) 
		// Contributors: 
		unsigned count = 0;
		for(unsigned i=0;i<G.nloci_chrom(1);i++){	
			if(G.chrom()[0].locus[i].qtl_ml=='r') count++;
		}
		if(count>1){
			throw exception("Only one recessive locus allowed in MIM");
		}
		if(count==1 && modelTrmVec.size()>1){
			throw exception("More than one term in MIM::monoMendel()");
		}
		if(count==1){
			individualModelTerm = 0;
			return true;
		}
		else {
			return false;
		}
	}
		
	void MIM::setupPopulation(RPedigree& P,GeneticDist& G){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (August, 2005) 
		// Contributors: 
		isMIMReady();
		G.numtrait = numTraits;	
		if (pop) delete pop;
		pop = new Population;
		check_ptr(pop);
		try {
			pop->setupGNodeSampler(P,*this,G);
			//pop->ListAlleleFounders();        // don't need this for sampling genotypes RLF                        
			//pop->SetPossibleHaplotypes(); 
		}
		catch (exception &ex) {
			cerr << " Failed to setupGNodeSampler! " << endl;  
			throw;
		}
	}
		
	void MIM::transferTraitData(){
		// Authors: L. Radu Totir and Rohan L. Fernando 
		// (August, 2005) 
		// Contributors: 
		individualModelTerm = colName.getIndex("individual");
		inputData();
		DataNode missData, tempData;
		for(unsigned i=0;i<dataVec.size();i++){
			unsigned ind = dataVec[i].trmVec[individualModelTerm].level;
			for(unsigned j=0;j<numTraits;j++){
				if (dataVec[i].misVec[j]){
					pop->popmember[ind-1]->myrecord[j] = missData;
				}
				else{
					tempData = dataVec[i].depVec[j];
					pop->popmember[ind-1]->myrecord[j] = tempData;
				}
			}
		}
	}
	
void MIM::GNodeSamplerOBGSL(void){
    // Authors: Rohan L. Fernando and  L. Radu Totir
    // (October, 2005) 
    // Contributors: 
	cout << "\n Entered Overlapping Blocking Gibbs Sampler for Single Locus (GNodeSamplerOBGSL)\n\n"; 
	if(pop->prior->nloci_chrom(1) > 1){
		throw exception("MIM::GNodeSamplerOBGSL: cannot have more than one locus \n");
	}
	string typeOfGNode = "genotypic";
	unsigned nAlleles = pop->prior->chrom()[0].locus[0].nallele();
	if (nAlleles > 2){
		pop->getInitialGNodeListSampleMIM();
	}
	else {
		pop->setMissingToHet(0);
		pop->prior->chrom()[0].locus[0].peelOrder.resize(pop->popsize);
	}
	pop->initGenotypeNodeList(0);
	pop->gNodeList.saveGNodeStates();
	//pop->gNodeList.displayGNodeSets();
	GNodeSet pivotSet = pop->getSirePivots(0);
	std::cout << "Number of blocks:" << pivotSet.size() << std::endl;
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned printFlag     = parmMap.getUnsignedValue("printFlag");
	unsigned blockWidth    = parmMap.getUnsignedValue("blockWidth");
	if (printFlag>4) pop->gNodeList.displayGNodeSets();
	string prune = parmMap["prunePedigree"];
	if(prune=="yes"){
		pop->gNodeList.pruneUninformativeTerminalGenotypeNodes();
		if (printFlag>4){
			cout << "after pruning\n";
			pop->gNodeList.displayGNodeSets();
		}
	}
	pop->initCounters();
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	GNodeList blockList, bigList;
	for (unsigned i=1;i<=chainLength;i++){
		bool gotSample = false;
		unsigned nTry  = 0;
		while (!gotSample) {
			try {
				cout << "\n OBGSL Round = " << i << endl;
				GNodeSet::iterator pivotIt;
				unsigned blockj = 0;
				for (pivotIt=pivotSet.begin();pivotIt!=pivotSet.end();pivotIt++){
					unsigned tryWidth = blockWidth;
					bool peeled = false;
					do {
					    blockList.clear();
						blockList = pop->gNodeList.findNeighbors(*pivotIt,tryWidth);
						if (printFlag>0){
							cout<< " Block: " << ++blockj
						        << " GNode: " << (*pivotIt)->id 							    
								<< " Block size " << blockList.size() 
							    << " width " << tryWidth << "\r";
							cout.flush();		
                        }
						bigList.saveGNodeSets();
						blockList.setGNodeSampleFlags(false); 
						peeled = blockList.peelNoCut(maxCutSetSize);
						if (peeled){
							for (unsigned j=1;j<=1;j++){
								blockList.setGNodeSampleFlags(false);
								blockList.reverseSample();
								blockList.sampleSegregationIndicators();
								blockList.saveGNodeStates();
								blockList.updateCounters();
							}
						}
						else {
							blockList.recoverSavedGNodeStates();
							blockList.setGNodeSampleFlags(true);
							tryWidth--;
						}
						bigList.resetGNodeList();
						GNode::gNodeListPtr->releaseCutSets();
					} while (!peeled && tryWidth);
					if (tryWidth==0){
						throw exception("\n MIM::GNodeSamplerOBGSL: failed in peeling \n");
					}
				}
				GNode::gNodeListPtr->setUninformativeGNodesNotSampled();
				GNode::gNodeListPtr->sampleUpdateUninformativeTerminalGNodes();
				gotSample = true;
			}
			catch (matvec::NormalizeFailed){
				if (nTry++ > 5) {
					cerr << "\n Normalize has failed more than five times in this sample \n";
					cerr << "Suspect something is wrong \n";
					cerr << "Going to quit here :-(.\n";
					exit (0);
				}
				else {
					cerr << "\n Normalize has failed    \n";
					cerr << "Will try sampling again \n";
					blockList.recoverSavedGNodeStates();
					blockList.setGNodeSampleFlags(true);
					bigList.resetGNodeList();
					GNode::gNodeListPtr->releaseCutSets();
				}
			}
		}
	}
}

void MIM::GNodeSamplerEEGSL(void){
    // Authors: Rohan L. Fernando and  L. Radu Totir
    // (October, 2005) 
    // Contributors: 
	cout << "\n Entered EE Blocking Gibbs Sampler for Single Locus (GNodeSamplerOBGSL)\n\n"; 
	if(pop->prior->nloci_chrom(1) > 1){
		throw exception("MIM::GNodeSamplerEEGSL: cannot have more than one locus \n");
	}
	string typeOfGNode = "genotypic";
	unsigned nAlleles = pop->prior->chrom()[0].locus[0].nallele();
	if (nAlleles > 2){
		pop->getInitialGNodeListSampleMIM();
	}
	else {
		pop->setMissingToHet(0);
		pop->prior->chrom()[0].locus[0].peelOrder.resize(pop->popsize);
	}
	pop->initGenotypeNodeList(0);
	pop->gNodeList.saveGNodeStates();
	//pop->gNodeList.displayGNodeSets();
	GNodeSet pivotSet = pop->getSirePivots(0);
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned printFlag     = parmMap.getUnsignedValue("printFlag");
	unsigned blockWidth    = parmMap.getUnsignedValue("blockWidth");
	string prune = parmMap["prunePedigree"];
	if(prune=="yes") pop->gNodeList.pruneUninformativeTerminalGenotypeNodes();

	pop->initCounters();
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	for (unsigned i=1;i<=chainLength;i++){
	    cout << "Round = " << i << endl;
		GNodeSet::iterator pivotIt;
		for (pivotIt=pivotSet.begin();pivotIt!=pivotSet.end();pivotIt++){
			unsigned tryWidth = blockWidth;
			bool peeled = false;
			do {
				GNodeList blockList = pop->gNodeList.findNeighbors(*pivotIt,tryWidth);
				cout << "Blocking with pivot GNode: " << (*pivotIt)->id 
					 << " Block size " << blockList.size() 
					 << " width " << tryWidth << endl; 
			    if (printFlag>1) blockList.displayGNodeSets(); 
			    GNodeList bigList = blockList.getBigList();
			    bigList.saveGNodeSets();
			    blockList.setGNodeSampleFlags(false); 
			    peeled = blockList.peelNoCut(maxCutSetSize);
				if (peeled){
					blockList.EEUpdateCounters();
					blockList.reverseSample();
				}
				else {
					blockList.recoverSavedGNodeStates();
				    blockList.setGNodeSampleFlags(true);
					tryWidth--;
				}
				bigList.resetGNodeList();
				GNode::gNodeListPtr->releaseCutSets();
				blockList.resetBackSets();
			} while (!peeled && tryWidth);
			if (tryWidth==0){
				throw exception("MIM::GNodeSamplerOBGSL: failed in peeling \n");
			}
		}
		GNode::gNodeListPtr->setUninformativeGNodesNotSampled();
		GNode::gNodeListPtr->sampleUpdateUninformativeTerminalGNodes();
	}
}

//void MIM::GNodeSamplerEEGSL(void){
//    // Authors: Rohan L. Fernando
//    // (October, 2005) 
//    // Contributors: 
//	cout << "\n Entered  EEGS for single locus\n\n"; 
//	if(pop->prior->nloci_chrom(1) > 1){
//		throw exception("MIM::GNodeSamplerEEGSL: cannot have more than one locus \n");
//	}
//	string typeOfGNode = "genotypic";
//	unsigned nAlleles = pop->prior->chrom()[0].locus[0].nallele();
//	if (nAlleles > 2){
//		pop->getInitialGNodeListSampleMIM();
//	}
//	else {
//		pop->setMissingToHet(0);
//		pop->prior->chrom()[0].locus[0].peelOrder.resize(pop->popsize);
//	}
//	pop->initGenotypeNodeList(0);
//	GNodeSet pivotSet = pop->getSirePivots(0);
//	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
//	unsigned printFlag     = parmMap.getUnsignedValue("printFlag");
//	unsigned blockWidth    = parmMap.getUnsignedValue("blockWidth");
//	if (printFlag>1) pop->gNodeList.displayGNodeSets();
//	pop->initCounters();
//	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
//	for (unsigned i=1;i<=chainLength;i++){
//	    cout << "Round = " << i << endl;
//		GNodeSet::iterator pivotIt;
//		for (pivotIt=pivotSet.begin();pivotIt!=pivotSet.end();pivotIt++){
//		    bool peeled = false;
//			GNodeList blockList = pop->gNodeList.findNeighbors(*pivotIt,blockWidth);
//			cout << "Blocking with pivot GNode " << (*pivotIt)->id 
//				<< " Block size " << blockList.size() << endl; 
//			if (printFlag>1) blockList.displayGNodeSets();
//			GNodeList bigList = blockList.getBigList();
//			bigList.saveGNodeSets();
//			blockList.setGNodeSampleFlags(false); 
//			peeled = blockList.peelNoCut(maxCutSetSize);			
//			//blockList.displayGNodeSets();
//			if (!peeled){
//				throw exception("GNodeList::orderAndPeelMIM: failed to peel without cutting\n");
//			}
//			blockList.EEUpdateCounters();
//			blockList.reverseSample();
//			bigList.resetGNodeList();
//			GNode::gNodeListPtr->releaseCutSets();
//			blockList.resetBackSets();
//		}
//	}
//}

void MIM::GNodeSamplerGS1(void){
    // Authors: L. Radu Totir and Rohan L. Fernando 
    // (September, 2005) 
    // Contributors: 
    int    badCount  = 0;
    int    goodCount = 0;
    double alpha      = 1.0;
    double ranNumber  = 0.0;
    double logOldProposal, logNewProposal, logOldTarget, logNewTarget;
    pop->gNodeList.logProposal = 0.0;
    pop->gNodeList.logTarget   = 0.0;
    string typeOfGNode = parmMap["typeOfGNodes"];
    // get an initial sample
    pop->getInitialGNodeListSampleMIM();
    logOldTarget   = pop->gNodeList.logTarget;
    logOldProposal = pop->gNodeList.logProposal;
    pop->copyGNodeStatesToCandidateStates(typeOfGNode);
    pop->storeSampledGametes();
    pop->copyCandidateToAccepted(typeOfGNode);
    // initialize structures for capturing results
    pop->initCounters();
    // continue the sampling process
    unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
    unsigned numGibbsSteps = parmMap.getUnsignedValue("numGibbsSteps");
    unsigned numMHTrials   = parmMap.getUnsignedValue("numMHTrials");	
    unsigned count=0;
    while (count<chainLength){
        // Gibbs step 
        for (unsigned j =0;j<numGibbsSteps;j++){
            pop->gNodeList.logProposal = 0.0;
            pop->gNodeList.logTarget   = 0.0;
            pop->getGNodeListSampleGibbsMIM();
            if(j==(numGibbsSteps-1)){
                // keep track of the target probability for MH step
                logOldTarget   = pop->gNodeList.logTarget;
            }
			pop->copyGNodeStatesToCandidateStates(typeOfGNode);
			pop->storeSampledGametes();
			pop->copyCandidateToAccepted(typeOfGNode);
			pop->updateCounters();
            count++;
			cout << "Sample " << setw(5) << count << endl; 
        }	
        // MH step	
        for (unsigned j =0;j<numMHTrials;j++){
            pop->gNodeList.logProposal = 0.0;
            pop->gNodeList.logTarget   = 0.0;
            pop->getGNodeListSampleSeqInMIM();
            pop->copyGNodeStatesToCandidateStates(typeOfGNode);
            pop->storeSampledGametes();
            logNewProposal = pop->gNodeList.logProposal;
            logNewTarget   = pop->gNodeList.logTarget;
            pop->gNodeList.logOldProposal = 0.0;
            pop->getOldGNodeListProbabilityMIM();
            logOldProposal = pop->gNodeList.logOldProposal;
            cout << "logNewTarget   = " << logNewTarget   << endl;
            cout << "logOldProposal = " << logOldProposal << endl;
            cout << "logOldTarget   = " << logOldTarget   << endl;
            cout << "logNewProposal = " << logNewProposal << endl;
            alpha=std::exp(logNewTarget + logOldProposal - logNewProposal - logOldTarget);
            cout << " alpha = " << alpha << endl;
            ranNumber=ranf();
            if (ranNumber <= alpha) {
                logOldTarget = logNewTarget;
                logOldProposal = logNewProposal;
                pop->copyCandidateToAccepted(typeOfGNode);
                goodCount++;
                count++;
                pop->updateCounters();	
                cout << "Got a valid MH sample after " << j+1 << " tries." << endl;
                break;
            }
            else {
                pop->retreiveSampledGametes();
                badCount++;
                count++;
                pop->updateCounters();	
            }
			cout << "Sample " << setw(5) << count << endl; 
        }
    }
    cout << "No. of additional MH samples besides the initial one = " << goodCount << endl;
    cout << "Each new MH sample is followed by: " << numGibbsSteps << " Gibbs samples " << endl; 
}

void MIM::GNodeSamplerJNC(void){
    // Authors: Rohan L. Fernando
    // (January, 2006) 
    // Contributors: 
	cout << "\n Entered Joint Sampler without cutting (GNodeSamplerJNC)\n\n";
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	pop->initJointAlleleNodeList(numLoci);
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	pop->gNodeList.setGNodeSampleFlags(false);
	bool peeled = pop->gNodeList.peelNoCut2(maxCutSetSize);
	pop->initCounters();
	if (!peeled) {
		throw exception("MIM::GNodesamplerJNC: failed to peel \n");
	}
	for (unsigned i=1;i<=chainLength;i++){
		cout << "Round = " << i << endl;
		pop->gNodeList.setGNodeSampleFlags(false);
		pop->gNodeList.reverseSample();
		pop->gNodeList.saveGNodeStates();
		pop->updateCounters();
	}
}

void MIM::GNodeSamplerJNCIDG(void){
    // Authors: Rohan L. Fernando
    // (January, 2006) 
    // Contributors: 
	cout << "\n Entered Joint Sampler without cutting (GNodeSamplerJNCIDG)\n\n";
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	pop->initJointAlleleNodeList(numLoci);
	unsigned maxCutSetSize   = parmMap.getUnsignedValue("maxCutSetSize");
	pop->gNodeList.setGNodeSampleFlags(false);
	bool peeled = pop->gNodeList.peelNoCut(maxCutSetSize);
	if (!peeled) {
		throw exception("MIM::GNodesamplerJNC: failed to peel \n");
	}
	pop->gNodeList.setGNodeSampleFlags(false);
	pop->gNodeList.reverseSample();
	string outFile = parmMap["resultsFile"];
	pop->gNodeList.saveGNodeStates();
	pop->output_descentGraphJP(outFile);
}

void MIM::displayGenotypicValues(RPedigree &ped, std::ostream &outfile){

    for (unsigned mtrm=0;mtrm<modelTrmVec.size();mtrm++){
		if (modelTrmVec[mtrm].name != "individual") continue;
		outfile << "Genotypic values for trait: " << modelTrmVec[mtrm].depVarName << endl;
		outfile << setw(25) << "Individual " << setw(25) << "Genotypic Value" << endl;
		for (unsigned i=0;i<ped.pedVector.size();i++){
			unsigned ii = modelTrmVec[mtrm].start + i;
			unsigned ti = modelTrmVec[mtrm].trait;
			double   gv = sol[ii];
			for (unsigned mqtl=0;mqtl<MQTLVec.size();mqtl++){
				MQTL *ptr = MQTLVec[mqtl];
				ii = ptr->MQTLStart + (*ptr)[i]->mLevel - 1;
				gv += sol[ii];
				ii = ptr->MQTLStart + (*ptr)[i]->pLevel - 1;	
				gv += sol[ii];
				ii = ptr->QTLProbStart;
				gv += sol[ii]*ptr->regCoeff[ti];
			}	
			outfile << setw(25) << ped.pedVector[i]->ind_str 
					<< setw(25) << gv << endl;		
		}
	}	
}

void MIM::GNodeSamplerOBGSL1(void){
    // Authors: Rohan L. Fernando and  L. Radu Totir
    // (October, 2005) 
    // Contributors: 
	cout << "\n Entered Overlapping Blocking Gibbs Sampler for Single Locus (GNodeSamplerOBGSL)\n\n"; 
	if(pop->prior->nloci_chrom(1) > 1){
		throw exception("MIM::GNodeSamplerOBGSL: cannot have more than one locus \n");
	}
	string typeOfGNode = "genotypic";
	unsigned nAlleles = pop->prior->chrom()[0].locus[0].nallele();
	if (nAlleles > 2){
		pop->getInitialGNodeListSampleMIM();
	}
	else {
		pop->setMissingToHet(0);
		pop->prior->chrom()[0].locus[0].peelOrder.resize(pop->popsize);
	}
	pop->initGenotypeNodeList(0);
	pop->gNodeList.saveGNodeStates();
	//pop->gNodeList.displayGNodeSets();
	GNodeSet pivotSet = pop->getSirePivots(0);
	SafeSTLVector<GNodeList> vectorBlockList;
	GNodeSet::iterator pivotIt;
	std::cout << "Number of blocks:" << pivotSet.size() << std::endl;
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned printFlag     = parmMap.getUnsignedValue("printFlag");
	unsigned blockWidth    = parmMap.getUnsignedValue("blockWidth");
	string prune = parmMap["prunePedigree"];
	cout<<"GnodeList before pruning:"<<endl;
	if (printFlag>4) pop->gNodeList.displayGNodeSets();
	if(prune=="yes") pop->gNodeList.pruneUninformativeTerminalGenotypeNodes();
	cout<<"\nGnodeList after pruning:"<<endl;
	if (printFlag>4) pop->gNodeList.displayGNodeSets();
	for (pivotIt=pivotSet.begin();pivotIt!=pivotSet.end();pivotIt++){
		vectorBlockList.push_back(pop->gNodeList.findNeighbors(*pivotIt,blockWidth));
	}
	
	pop->initCounters();
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	GNodeList bigList;
	for (unsigned i=1;i<=chainLength;i++){
		bool gotSample = false;
		while (!gotSample) {
			try {
				cout << "Round = " << i << endl;
				pivotIt = pivotSet.begin();
				for (unsigned blocki=0;blocki<vectorBlockList.size();blocki++){
					if (printFlag>0) {
						cout << "Blocking with pivot GNode: " << (*pivotIt)->id 
						<< " Block size " << vectorBlockList[blocki].size() 
						<< " width " << blockWidth << endl;
					} 
					pivotIt++;
					bool peeled = false;
					bigList.clear();
					bigList = vectorBlockList[blocki].getBigList();
					bigList.saveGNodeSets();
					vectorBlockList[blocki].setGNodeSampleFlags(false); 
					if (printFlag>1) vectorBlockList[blocki].displayGNodeSets();
					if(i == 1){
						peeled = vectorBlockList[blocki].peelNoCut(maxCutSetSize);
					}
					else{
						vectorBlockList[blocki].peelNoOrder(maxCutSetSize);
						peeled = true;					 	
					}
					if (peeled){
						for (unsigned j=1;j<=1;j++){
							vectorBlockList[blocki].setGNodeSampleFlags(false);
							vectorBlockList[blocki].reverseSample();
							vectorBlockList[blocki].sampleSegregationIndicators();
							vectorBlockList[blocki].saveGNodeStates();
							vectorBlockList[blocki].updateCounters();
						}
					}
					else {
						throw exception("MIM::GNodeSamplerOBGSLnew: failed in peeling \n");
					}
					bigList.resetGNodeList();
					GNode::gNodeListPtr->releaseCutSets();
				}
				GNode::gNodeListPtr->setUninformativeGNodesNotSampled();
				GNode::gNodeListPtr->sampleUpdateUninformativeTerminalGNodes();
				gotSample = true;
			}
			catch (matvec::NormalizeFailed){
				cerr << "Normalize has failed more than five times in this sample \n";
				cerr << "Suspect something is wrong \n";
				cerr << "Continuing for now :-(.\n";
				continue;
			}
		}
	}
}

void MIM::GNodeSamplerLG(void){
    // Authors: Rohan L. Fernando and Radu Totir
    // (May, 2006) 
    // Contributors: 
	cout << "\n Entered LG Sampler (GNodeSamplerLG)\n\n";
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	pop->initJointAlleleNodeList(numLoci);
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	bool peeled = pop->gNodeList.peelNoOrder(maxCutSetSize);
	pop->initCounters();
	if (!peeled) {
		throw exception("MIM::GNodesamplerLG: failed to peel \n");
	}
	for (unsigned i=1;i<=chainLength;i++){
		cout << "Round = " << i << endl;
		pop->gNodeList.setGNodeSampleFlags(false);
		pop->gNodeList.reverseSample();
		pop->gNodeList.saveGNodeStates();
		pop->updateCounters();
	}
}

void MIM::GNodeSamplerSLNC(void){
    // Authors: Rohan L. Fernando and David Habier
    // (July, 2006) 
    // Contributors: 
	cout << "\n Entered Single Locus NoCut Sampler (GNodeSamplerSLNC)\n\n";
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	string typeOfGNode = parmMap["typeOfGNodes"];
	pop->initGenotypeNodeList(0);
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	bool peeled = pop->gNodeList.peelNoCut(maxCutSetSize);
	pop->initCounters();
	if (!peeled) {
		throw exception("MIM::GNodesamplerLG: failed to peel \n");
	}
	for (unsigned i=1;i<=chainLength;i++){
		cout << "Round = " << i << endl;
		pop->gNodeList.setGNodeSampleFlags(false);
		pop->gNodeList.reverseSample();
		pop->copyGNodeStatesToAcceptedStates(typeOfGNode);
		pop->updateCounters();
	}
}

void MIM::GNodeSamplerJointNC(void){
    // Authors: Rohan L. Fernando
    // (July, 2006) 
    // Contributors: 
	cout << "\n Entered Joint Sampler without cutting (GNodeSamplerJointNC)\n\n";
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	pop->initJointNodeList(numLoci);
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	pop->gNodeList.setGNodeSampleFlags(false);	
	bool peeled = pop->gNodeList.peelNoCut2(maxCutSetSize);
	pop->initCounters();
	if (!peeled) {
		throw exception("MIM::GNodesamplerJointNC: failed to peel \n");
	}
	for (unsigned i=1;i<=chainLength;i++){
		cout << "Round = " << i << endl;
		pop->gNodeList.setGNodeSampleFlags(false);
		pop->gNodeList.reverseSample();
		pop->updateCountersSimple();
	}
}

void MIM::GNodeProbsAPESJointNC(std::string inputFileName,RPedigree& P,GeneticDist& G){
    // Authors: Rohan L. Fernando and Radu Totir
    // (July, 2006) 
    // Contributors: 
	cout << "\n Entered GNode probability calculations MIM::GNodeProbsAPESJointNC(void)\n\n";
	putGNodeSamplerParms(inputFileName);
	setupPopulation(P,G);
	if(monoMendel(G)){       
		transferTraitData();
	}
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	pop->initJointNodeList(numLoci);
	unsigned maxCutSetSize   = MIM::parmMap.getUnsignedValue("maxCutSetSize");
	pop->gNodeList.setGNodeSampleFlags(false);	
	bool peeled = pop->gNodeList.peelNoCut(maxCutSetSize);
	if (!peeled) {
		throw exception(" MIM::GNodeProbsAPESJointNC(void) failed to peel \n");
	}
	string fname = parmMap["resultsFile"];
	if (fname==""){
		pop->displayGNodeProbs(std::cout);
	}
	else {
		std::ofstream outFile(fname.c_str());
		pop->displayGNodeProbs(outFile);
	}	
}

//SafeMVVector<double> MIM::getYHat(void){
//	SafeMVVector<double> yhat;
//	yhat.resize(dataVec.size(),0.0);
//	for (unsigned i=0;i<dataVec.size();i++){
//		for (unsigned mi=0;mi<numTerms;mi++){
//			ii = modelTrmVec[mi].start + dataVec[i].trmVec[mi].level - 1;
//			ti =  modelTrmVec[mi].trait;
//			vi = dataVec[i].trmVec[mi].value;
//			yhat[i] += vi*sol[ii];
//		}
//	}
//}

void MIM::DGSampleIBDMatrix(unsigned numOfSamples, unsigned numSL, unsigned numHaplo, unsigned numSLCas, unsigned leftLocus,
                            unsigned rightLocus, string fileName){
	// Authors: L. Rohan L. Fernando
	// (August, 2006) 
	// Contributors: 
	double l_hood = 0.0;
	unsigned k = 0, j = 0, option = 0;
	unsigned k1 = numSL;
	unsigned k2 = numHaplo;
	unsigned k3 = numSLCas;
	bool order = false;
	ibdMatrix.initialize(leftLocus,rightLocus);
	cout<< endl << "Descent graph iterations = " << numOfSamples << endl;
	for(int sample=1; sample <= numOfSamples; sample++) {
		if(sample%100 == 0) {
			cout<<sample<<".....";
			cout.flush();
		}
		//	cout<<endl<< "SAMPLE.................................. " << sample <<endl;  
		if (k >= k1 + k2 + k3)                {k = 0;}
		if (k < k1)                           {option = 1;}  
		else if ((k >= k1) && (k < k1 + k2))  {option = 2;}
		else if (k <= k1 + k2 + k3)           {option = 3;}
		k++;    
		j++;    
        l_hood = pop->M_H_sample_ext(2,l_hood, option, order);  
		ibdMatrix.sampleSimple();
	}
	ibdMatrix.output(fileName);


}

void MIM::GNodeSamplerOBGML(void){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
	cout << "\n Entered Multilocus Overlapping Blocking Gibbs Sampler (GNodeSamplerOBGML)\n\n";
	unsigned printFlag = parmMap.getUnsignedValue("printFlag");
	unsigned chainLength   = parmMap.getUnsignedValue("chainLength");
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	for(unsigned locus=0;locus<numLoci;locus++){
		unsigned nAlleles = pop->prior->chrom()[0].locus[locus].nallele();
		if (nAlleles > 2){
			pop->getInitialSampleForLocusMIM(locus);
		}
		else {
			pop->setMissingToHet(locus);
		}
		pop->setSegregationIndex(locus,"genotypic");
	}
	pop->sampleSegregationIndicatorsSimple();
	cout <<"making pedigree blocks " << endl;
	makePedigreeBlocks();
	if (printFlag>1) {
		for (unsigned ind=0;ind<pop->popsize;ind++){
			if (pop->popmember[ind]->pedBlock.size()){
				cout <<"pedBlock for dist:  " << pop->popmember[ind]->distanceToPivot0 << " pivot ind: " << pop->popmember[ind]->id() << endl;
				for (unsigned i=0;i<pop->popmember[ind]->pedBlock.size();i++){
					cout << pop->popmember[ind]->pedBlock[i]->id() << endl;
				}
			}
		}	
	}
	pop->initJointNodeList(numLoci);
	pop->gNodeList.setGNodeSampleFlags(true);
	makeBlockGNodeListVector();
	pop->initCounters();
	for (unsigned i=1; i<=chainLength;i++){
		for (unsigned blocki=0; blocki<blockNodeListVector.size(); blocki++){
			blockNodeListVector[blocki].setGNodeSampleFlags(false); 
			if (printFlag>1) blockNodeListVector[blocki].displayGNodeSets();
			blockNodeListVector[blocki].peelNoCutNoOrder();
			blockNodeListVector[blocki].setGNodeSampleFlags(false);
			blockNodeListVector[blocki].reverseSample();
			//blockNodeListVector[blocki].saveGNodeStates();
		}
		pop->updateCountersSimple(); // should consider doing this by blockNodeList
	}
}

void MIM::makePedigreeBlocks(void){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
	unsigned maxDist=0;
	string pivotIndStr = parmMap["pivotInd"];
	pop->calcDistanceToIndividual(pivotIndStr);
	for (unsigned ind=0;ind<pop->popsize;ind++){
		unsigned distToPivot0 = pop->popmember[ind]->distanceToPivot0;
		cout << pedPtr->pedVector[ind]->ind_str <<" " << distToPivot0 << endl;
		if(distToPivot0>maxDist) maxDist = distToPivot0;
	}
	unsigned blockWidth = parmMap.getUnsignedValue("pedigreeBlockWidth");
	unsigned overlap    = parmMap.getUnsignedValue("pedigreeBlockOverlap");
	unsigned fromDist   = 0;
	unsigned toDist     = blockWidth;
	toDist = (toDist>maxDist) ? maxDist : toDist;
	SafeSTLVector<Individual*> pedBlockPivots;
	while(toDist <= maxDist){
		pedBlockPivots.clear();
		for (unsigned ind=0;ind<pop->popsize;ind++){
			if(pop->popmember[ind]->distanceToPivot0 == fromDist){
				pop->popmember[ind]->getPedBlockUpto(toDist);
				if (pop->popmember[ind]->pedBlock.size()){
					Individual* pivot = pop->popmember[ind];
					pedBlockPivots.push_back(pivot);
				}
			}
		}
		mergeAdjacentPedBlocks(pedBlockPivots);
		fromDist += blockWidth + 1 - overlap;
		unsigned width = ((maxDist - fromDist) < blockWidth) ?  maxDist - fromDist : blockWidth;
		if (width==0) break;
		toDist    = fromDist + width;
	}
	removeRedundantBlocks();
}

void MIM::mergeAdjacentPedBlocks(SafeSTLVector<Individual*>& pedBlockPivots){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 
    SafeSTLVector<Individual*>::iterator first, last;
	for (unsigned i=0;i<pedBlockPivots.size()-1;i++){
		Individual* iPivot = pedBlockPivots[i];
		for (unsigned j=i+1;j<pedBlockPivots.size();j++){
			Individual* jPivot = pedBlockPivots[j];
			if (iPivot->pedBlockAdjacent(jPivot->pedBlock)){
				first = iPivot->pedBlock.begin();
				last  = iPivot->pedBlock.end();
				jPivot->pedBlock.insert(jPivot->pedBlock.begin(),first,last);
				jPivot->pedBlock.insert(jPivot->pedBlock.begin(),iPivot);
				iPivot->pedBlock.clear();
				break;
			}
		}
	}
}

void MIM::removeRedundantBlocks(void){
    // Authors: Rohan L. Fernando 
    // (October, 2006) 
    // Contributors: 

	// in addition to removing redundant blocks, this method
	// also includes pivot in its block
	
	unsigned overlap = parmMap.getUnsignedValue("pedigreeBlockOverlap");
	for (unsigned i=0;i<pop->popsize;i++){
		Individual* ind = pop->popmember[i];
		if (ind->pedBlock.size()==0) continue;
		ind->pedBlock.insert(ind->pedBlock.begin(),ind);
		unsigned pivotDistance = ind->distanceToPivot0;
		unsigned maxDist = 0;
		for (unsigned j=0;j<ind->pedBlock.size();j++){
			if (maxDist<ind->pedBlock[j]->distanceToPivot0) maxDist = ind->pedBlock[j]->distanceToPivot0;
		}
		if (maxDist < pivotDistance+overlap){
			ind->pedBlock.clear();
		}
	}	
}

void MIM::makeBlockGNodeListVector(void){
    // Authors: Rohan L. Fernando
    // (October, 2006) 
    // Contributors: 
	
	blockNodeListVector.name = "blockNodeListVector";	
	unsigned numLoci = Population::prior->chromosome[0].nloci();
	unsigned blockWidth = parmMap.getUnsignedValue("locusBlockWidth");
	unsigned overlap    = parmMap.getUnsignedValue("locusBlockOverlap");
	for (unsigned startLocus=0;startLocus<numLoci;startLocus+=(blockWidth-overlap)){
		unsigned endLocus = (startLocus+blockWidth<numLoci) ? startLocus+blockWidth-1 : numLoci-1;
		addBlocksFor(startLocus,endLocus);
	}
} 

void MIM::addBlocksFor(unsigned startLocus, unsigned endLocus){
    // Authors: Rohan L. Fernando
    // (October, 2006) 
    // Contributors: 
	cout <<"addBlocksFor locus "<< startLocus <<" to " << endLocus << endl;  
	unsigned maxCutSetSize = parmMap.getUnsignedValue("maxCutSetSize");
	GNodeList blockList, bigList;
	BlockNodeList blockNodeList;
	for (unsigned ind=0;ind<pop->popsize;ind++){
		if (pop->popmember[ind]->pedBlock.size()){
			blockList.clear();
			for (unsigned j=0;j<pop->popmember[ind]->pedBlock.size();j++){
				Individual* indjPtr = pop->popmember[ind]->pedBlock[j];
				for (unsigned locus=startLocus;locus<=endLocus;locus++){
					if (pop->prior->chrom()[0].locus[locus].gnodeType=="allelic"){
						AlleleStateNode* mAlleleState = &indjPtr->malleleStateNodeVector[locus];
						AlleleStateNode* pAlleleState = &indjPtr->palleleStateNodeVector[locus];
						blockList.push_back(mAlleleState);
						blockList.push_back(pAlleleState);
					}
					else if (pop->prior->chrom()[0].locus[locus].gnodeType=="genotypic"){
						GenotypeNode* indGenotype = &indjPtr->genotNodeVector[locus];
                        blockList.push_back(indGenotype);
				    }
					if (indjPtr->mymother){
						AlleleOriginNode *mAlleleOrigin = &indjPtr->malleleOriginNodeVector[locus];
						AlleleOriginNode *pAlleleOrigin = &indjPtr->palleleOriginNodeVector[locus];
						blockList.push_back(mAlleleOrigin);
						blockList.push_back(pAlleleOrigin);
					}
				}
			}
			bigList.clear();
			bigList = blockList.getBigList();
			bigList.saveGNodeSets();
			blockList.setGNodeSampleFlags(false); 
			bool peeled = blockList.peelNoCut2(maxCutSetSize);
			if (peeled){
				blockNodeList.initNodes(blockList);
				blockNodeListVector.push_back(blockNodeList);
				blockList.setGNodeSampleFlags(false);
				blockList.reverseSample();
				blockList.saveGNodeStates();
				//blockList.updateCounters();
			}
			else {
				throw exception("MIM::GNodeSamplerOBGML: failed in peeling \n");
			}
			bigList.resetGNodeList();
		}	
		
	}
}
	

	
} ////// end of namespace matvec

