#include "rootIncludes.inc"
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TROOT.h"

#include "upsilon_MCMassFit.C"

using namespace onia;

//====================================
int main(int argc, char** argv){

	int nState=999;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
	    }

	Double_t nSigma = 2.; //sigma_mass cut for preparation of figures
	Char_t *fileNameIn = "tmpFiles/selEvents_data_Ups.root";

	char ProcessStateMacro[200];
	sprintf(ProcessStateMacro,".L upsilon_MCMassFit.C+");

	  gROOT->ProcessLine(ProcessStateMacro);

  for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 0; iPT < max; iPT++){

    	upsilon_MCMassFit(iRap, iPT, nSigma, fileNameIn, nState);

    }
  }

  return 0;
}
