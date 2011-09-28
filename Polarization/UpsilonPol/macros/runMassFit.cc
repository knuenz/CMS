#include "rootIncludes.inc"

#include "TROOT.h"

#include "upsilon_2StepFit.C"

using namespace onia;

//====================================
int main(){

	Double_t nSigma = 2.; //sigma_mass cut for preparation of figures
	Char_t *fileNameIn = "tmpFiles/selEvents_data_Ups.root";

  gROOT->ProcessLine(".L upsilon_2StepFit.C+");

  for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 0; iPT < max; iPT++){
      upsilon_2StepFit(iRap, iPT, nSigma, fileNameIn);
    }
  }

  return 0;
}
