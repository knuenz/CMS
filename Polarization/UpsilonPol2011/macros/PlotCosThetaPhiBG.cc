#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"

#include "TH2.h"

enum{L,R};
Char_t *bgLabel[2] = {"L", "R"};
TH2F *hCosThetaPhi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2];
void LoadHistos(Int_t iRapBin, Int_t iPTBin);
void PlotHistos(Int_t iRapBin, Int_t iPTBin, Int_t iFrame, Int_t iWindow);
//===========================
int main(){

	gStyle->SetPalette(1);

  for(int iRap = 0; iRap <= 2; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 0; iPT < max; iPT++){
      //for(int iPT = 1; iPT < 2; iPT++){
      LoadHistos(iRap, iPT);
      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
	PlotHistos(iRap, iPT, iFrame, L);
	PlotHistos(iRap, iPT, iFrame, R);
      }
    }
  }
  return 0;
}
//===========================
void PlotHistos(Int_t iRapBin, Int_t iPTBin, Int_t iFrame, Int_t iWindow){


  gStyle->SetPadRightMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetFrameBorderMode(0);

  Char_t name[100];
  sprintf(name, "c1_%s_rap%d_pT%d_%s", onia::frameLabel[iFrame], iRapBin, iPTBin, bgLabel[iWindow]);
  TCanvas *c1 = new TCanvas(name, "", 500, 500);
  c1->SetFillColor(kWhite);

  double yOffset=1.4;
  hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow]->GetYaxis()->SetTitleOffset(yOffset);
  gPad->SetLeftMargin(0.125);

  hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow]->Draw("colz");
  sprintf(name, "Figures/cosThetaPhi_%s_rap%d_pT%d_%s.pdf", onia::frameLabel[iFrame], iRapBin, iPTBin, bgLabel[iWindow]);
  c1->Print(name);
}

//===========================
void LoadHistos(Int_t iRapBin, Int_t iPTBin){


  Char_t name[100];
  sprintf(name, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);

  TFile *fIn = new TFile(name);

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
	  cout<<endl;
    sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
    hCosThetaPhi[iRapBin][iPTBin][iFrame][L] = (TH2F *) gDirectory->Get(name);

    hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->Print();
    sprintf(name, "hCosThetaPhi_%s_pT%d_rap%d_L", onia::frameLabel[iFrame], iRapBin, iPTBin);
    hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->SetName(name);
    //
    sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
    hCosThetaPhi[iRapBin][iPTBin][iFrame][R] = (TH2F *) gDirectory->Get(name);
    sprintf(name, "hCosThetaPhi_%s_pT%d_rap%d_R", onia::frameLabel[iFrame], iRapBin, iPTBin);
    hCosThetaPhi[iRapBin][iPTBin][iFrame][R]->SetName(name);
  }
}
