#include "rootIncludes.inc"

#include "TLorentzVector.h"
#include "TSystem.h"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};

void calcMeanPt(Int_t iRapBin = 1,
		      Int_t iPTBin = 1,
		      Double_t nSigma = 2.,
		      Int_t nUpsState=0,//[0]... 1S, [1]... 2S, [2]... 3S
		      Int_t output=999,
		      Double_t fracL = 0.5){


  Char_t fileNameIn[100];
  sprintf(fileNameIn, "AllStates_%1.2fSigma_FracLSB%dPercent/data_%dSUps_rap%d_pT%d.root", nSigma, int(fracL*100), nUpsState+1, iRapBin, iPTBin);
  //==============================
  //read inputs from input file:
  TFile *fIn = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");
  //==============================

  //==========================================================
  //reading fit parameters to establish signal mass window
  //as well as the L and R sideband window for the 3D BG histo
  //==========================================================

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();

  TH1D*  onia_pThisto   = new TH1D( "onia_pThisto","", 100, 0, 100);

  TH1D* backfrac_ = (TH1D*)fIn->Get("background_fraction");
  double backfrac = backfrac_->GetBinContent(1);

  int Ncut=0;
  double pTcut=3;
  bool Cut(false);

  Double_t onia_Pt;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);

    *onia = *(lepP) + *(lepN);
    onia_Pt = onia->Pt();

    Cut=false;
    	if(lepP->Pt() < pTcut || lepN->Pt() < pTcut) {Ncut++;continue;}
//   	if(lepP->Eta() < 0.3 && lepP->Eta() > 0.2 || lepN->Eta() < 0.3 && lepN->Eta() > 0.2) {Ncut++;continue;}

//    if(lepP->Eta() < 0.3 && lepP->Pt() < 6. || lepP->Eta() > 0.3 && lepP->Pt() < 4.6) Cut=true;
//    	else if(lepN->Eta() < 0.3 && lepN->Pt() < 6. || lepN->Eta() > 0.3 && lepN->Pt() < 4.6) Cut=true;
//    	if(Cut){Ncut++;continue;}

//    	if(lepP->Pt() < pTcut || lepN->Pt() < pTcut) {Ncut++;continue;}
   	onia_pThisto->Fill( onia_Pt );

  }

	char filename[200];
	sprintf(filename, "AllStates_%1.2fSigma_FracLSB%dPercent/meanPtlll.txt", nSigma, int(fracL*100));


  if(output==0) {cout<<onia_pThisto->GetMean();}
  if(output==1) {cout<<int((treeIn->GetEntries()-Ncut)*(1-backfrac));}
  if(output==2) {cout<<backfrac;}


  fIn->Close();

}
