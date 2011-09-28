#include <iostream>
#include <string>
#include <sstream>
using namespace std;


#include "rootIncludes.inc"

#include "TROOT.h"

#include "TrimEventContent.C"

using namespace onia;

//====================================
int main(int argc, char** argv){

	Double_t fracL = 0.5;
	Double_t nSigma = 2.;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("FracLSB") != std::string::npos) {char* fracLchar = argv[i]; char* fracLchar2 = strtok (fracLchar, "p"); fracL = atof(fracLchar2); fracL=fracL/100; cout<<"fracLSB = "<<fracL<<endl;}
		    if(std::string(argv[i]).find("nSigma") != std::string::npos) {char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "p"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl;}
	    }

  gROOT->ProcessLine(".L TrimEventContent.C+");

  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
	TrimEventContent(iRap, iPT, fracL, nSigma, iState);
      }
    }
  }
}
