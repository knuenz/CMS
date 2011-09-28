#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"

#include "TROOT.h"

#include "CopyTreeEntries.C"

using namespace onia;

//====================================
int main(){

  Char_t *fileNameIn = "tmpFiles/selEvents_data_Ups.root";

  gROOT->ProcessLine(".L CopyTreeEntries.C+");

  for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 1; iPT < max; iPT++){
      CopyTreeEntries(iRap, iPT, fileNameIn);
    }
  }
  return 0;
}
