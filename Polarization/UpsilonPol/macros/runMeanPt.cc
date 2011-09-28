#include <iostream>
#include <string>
#include <sstream>
using namespace std;


#include "rootIncludes.inc"
#include "commonVar.h"
#include "calcMeanPt.C"

#include "TROOT.h"

//void calcMeanPt(Int_t iRapBin, Int_t iPTBin, Double_t nSigma, Int_t nUpsState, Int_t output);

//====================================
int main(int argc, char** argv){

	Double_t nSigma = 2.;
	Double_t fracL = 0.5;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("FracLSB") != std::string::npos) {char* fracLchar = argv[i]; char* fracLchar2 = strtok (fracLchar, "p"); fracL = atof(fracLchar2); fracL=fracL/100; cout<<"fracLSB = "<<fracL<<endl;}
		    if(std::string(argv[i]).find("nSigma") != std::string::npos) {char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "p"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl;}
	    }

	gROOT->ProcessLine(".L calcMeanPt.C+");

	char filename[200];
	sprintf(filename,"AllStates_%1.2fSigma_FracLSB%dPercent/meanPtlll.txt", nSigma, int(fracL*100));
	gSystem->Unlink(filename);


	  int numEvents=0;

		  for(int iState = 0; iState < 3; iState++){

			  cout<<"Upsilon("<<iState+1<<"):"<<endl;

  cout<<endl;
  cout<<endl;

  cout<<"double ptCentre[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,0, fracL);
    	  if(iPT<max-1) {cout<<", ";  }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{";}
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;

  cout<<"int numEvents[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,1, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};"; }
      }
    }
//  }
  cout<<endl;
  cout<<endl;

  cout<<"double fracBackground[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,2, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;

		  }

  return 0;

}

