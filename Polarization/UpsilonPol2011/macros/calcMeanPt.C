#include "rootIncludes.inc"

#include "TLorentzVector.h"
#include "TSystem.h"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};
int numEvents_[4][100][100];
double meanPt[4][100][100];
double meanRap[4][100][100];
double BackgroundFrac[4][100][100];


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
  if(gDirectory->Get("selectedData")==NULL){
    return;
  }

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
  TH1D*  onia_raphisto   = new TH1D( "onia_raphisto","", 100, 0, 1.2);

  TH1D* backfrac_ = (TH1D*)fIn->Get("background_fraction");
  double backfrac = backfrac_->GetBinContent(1);

  int Ncut=0;
  double pTcut=3;
  bool Cut(false);

  Double_t onia_Pt;
  Double_t onia_Rap;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);

    *onia = *(lepP) + *(lepN);
    onia_Pt = onia->Pt();
    onia_Rap = onia->Rapidity();

    Cut=false;
    	if(lepP->Pt() < pTcut || lepN->Pt() < pTcut) {Ncut++;continue;}
//   	if(lepP->Eta() < 0.3 && lepP->Eta() > 0.2 || lepN->Eta() < 0.3 && lepN->Eta() > 0.2) {Ncut++;continue;}

//    if(lepP->Eta() < 0.3 && lepP->Pt() < 6. || lepP->Eta() > 0.3 && lepP->Pt() < 4.6) Cut=true;
//    	else if(lepN->Eta() < 0.3 && lepN->Pt() < 6. || lepN->Eta() > 0.3 && lepN->Pt() < 4.6) Cut=true;
//    	if(Cut){Ncut++;continue;}

//    	if(lepP->Pt() < pTcut || lepN->Pt() < pTcut) {Ncut++;continue;}
       	onia_pThisto->Fill( onia_Pt );
       	onia_raphisto->Fill( TMath::Abs(onia_Rap) );

  }

	char filename[200];
	sprintf(filename, "AllStates_%1.2fSigma_FracLSB%dPercent/meanPtlll.txt", nSigma, int(fracL*100));


/*	  sprintf(fileNameIn2, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
	  TFile *fIn2 = new TFile(fileNameIn2);
	  fIn2->cd();
	  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
	  TF1 *fUps[kNbSpecies], *fBG = 0;
	  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
	  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
	  treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
	  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
	  treeFitPar->SetBranchAddress("fBG", &fBG);
	  treeFitPar->LoadTree(0);
	  treeFitPar->GetEntry(0);

	  Double_t mass[kNbSpecies], sigma[kNbSpecies],massMin[kNbSpecies],massMax[kNbSpecies];
	  for(int iState = 0; iState < kNbSpecies; iState++){
	    mass[iState] = fUps[iState]->GetParameter(1);
	    sigma[iState] = fUps[iState]->GetParameter(2);
	    massMin[iState] = mass[iState]-nSigma*sigma[iState];
	    massMax[iState] = mass[iState]+nSigma*sigma[iState];
	  }

	  Char_t name[100];
	  sprintf(name, "Reco_Onia_mass_rap%d_pT%d", iRapBin, iPTBin);
	  TH1F* hMass = (TH1F*) gDirectory->Get(name);
	  hMass->Rebin(2);
	  double binWidth = hMass->GetBinWidth(1); //valid only for an equal bin histogram!
	  printf("binwidth = %1.2e\n", binWidth);

	  nUps[0] = fUps[0]->Integral(massMin[0], massMax[0])/binWidth;
	  nUps[1] = fUps[1]->Integral(massMin[1], massMax[1])/binWidth;
	  nUps[2] = fUps[2]->Integral(massMin[2], massMax[2])/binWidth;
*/
  if(output==0) {cout<<onia_pThisto->GetMean();meanPt[nUpsState][iRapBin][iPTBin]=onia_pThisto->GetMean();}
  if(output==1) {cout<<int((treeIn->GetEntries()-Ncut)*(1-backfrac)); numEvents_[nUpsState][iRapBin][iPTBin]=int((treeIn->GetEntries()-Ncut)*(1-backfrac));}
  if(output==2) {cout<<backfrac;BackgroundFrac[nUpsState][iRapBin][iPTBin]=100*backfrac;}
  if(output==5) {cout<<onia_raphisto->GetMean();meanRap[nUpsState][iRapBin][iPTBin]=onia_raphisto->GetMean();}

  TH1D* hMassScanInfo = (TH1D*)fIn->Get("hMassScanInfo");
  if(output==3) {cout<<hMassScanInfo->GetBinContent(1);}
  if(output==4) {cout<<hMassScanInfo->GetBinContent(2);}



  fIn->Close();

}
