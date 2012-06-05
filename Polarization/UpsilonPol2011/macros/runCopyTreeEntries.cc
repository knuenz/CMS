#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"

#include "TROOT.h"

#include "CopyTreeEntries.C"

using namespace onia;

//====================================
int main(int argc, char** argv){

	bool UpsMC=false;
	bool DoCPUconsumingPlots=false;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("UpsMC=1") != std::string::npos) {UpsMC=true;}
		    if(std::string(argv[i]).find("DoCPUconsumingPlots=1") != std::string::npos) {DoCPUconsumingPlots=true;}
	    }

  Char_t *fileNameIn = "tmpFiles/selEvents_data_Ups.root";

  gROOT->ProcessLine(".L CopyTreeEntries.C+");

  for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 0; iPT < max; iPT++){
      CopyTreeEntries(iRap, iPT, fileNameIn, UpsMC, DoCPUconsumingPlots);
    }
  }
  return 0;
}
