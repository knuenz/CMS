#include "rootIncludes.inc"
#include "TLorentzVector.h"
#include "calcPol.C"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};
//=====================================================
void CopyTreeEntries(Int_t iRapBin = 1, 
		     Int_t iPTBin = 1,
		     Char_t *fileNameIn = "RootFiles/selEvents_data_Ups_noCowboys_3Aug2011.root"){

  Char_t name[100], title[100];  
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
  printf("updating file %s\n", fileNameOut);

  //==============================
  TFile *fIn = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");
  //==============================

  //==============================
  //definition of output variables 
  TFile *fOut = new TFile(fileNameOut, "UPDATE");
  gStyle->SetPadRightMargin(0.2);
  TTree *treeOut = treeIn->CloneTree(0);
  TH2D *hBG_cosThetaPhi[onia::kNbFrames][2];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    //book the 2D (cosTheta, phi) histos for the L and R mass sideband
    sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    hBG_cosThetaPhi[iFrame][L] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
			       onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhi[iFrame][L]->Sumw2();
    //
    sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    hBG_cosThetaPhi[iFrame][R] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
			       onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhi[iFrame][R]->Sumw2();
  }

  //==============================
  //reading info from input file:
  //===============================
  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
  TF1 *fUps[kNbSpecies], *fBG = 0;
  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
  //treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
  treeFitPar->LoadTree(0);
  treeFitPar->GetEntry(0);
  fUps[0]->Print();

  Double_t mass1S = fUps[UPS1S]->GetParameter(1);
  Double_t sigma1S = fUps[UPS1S]->GetParameter(2);
  Double_t mass3S = fUps[UPS3S]->GetParameter(1);
  Double_t sigma3S = fUps[UPS3S]->GetParameter(2);
  printf("1S: mass = %1.3f, sigma = %1.3f\n", mass1S, sigma1S);
  printf("3S: mass = %1.3f, sigma = %1.3f\n", mass3S, sigma3S);
  Double_t massMin[2], massMax[2];
  massMin[L] = onia::massMinL;
  massMax[L] = mass1S - onia::nSigmaL*sigma1S;
  massMin[R] = mass3S + onia::nSigmaR*sigma3S;
  massMax[R] = onia::massMaxR;
  printf("--> L mass window: %1.3f < M < %1.3f GeV\n", massMin[L], massMax[L]);
  printf("--> R mass window: %1.3f < M < %1.3f GeV\n", massMin[R], massMax[R]);

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();

  Double_t onia_mass;
  Int_t index;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 100000 == 0)
      cout << "entry " << iEntry << " out of " << treeIn->GetEntries() << endl;

    *onia = *(lepP) + *(lepN);
    if(onia->Pt() > onia::pTRange[iRapBin][iPTBin-1] && onia->Pt() < onia::pTRange[iRapBin][iPTBin] &&
       TMath::Abs(onia->Rapidity()) > onia::rapForPTRange[iRapBin-1] && TMath::Abs(onia->Rapidity()) < onia::rapForPTRange[iRapBin]){

      treeOut->Fill(); //stores TLorenzVectors of the two muons in the given pT and rap cell

      //store now the cosTheta and phi distributions of the BG:
      onia_mass = onia->M();
      if(onia_mass > massMin[L] && onia_mass < massMax[L])
	index = L;
      else if(onia_mass > massMin[R] && onia_mass < massMax[R])
	index = R;
      else 
	continue;

      calcPol(*lepP, *lepN);

      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
	hBG_cosThetaPhi[iFrame][index]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
    }
  }
  fOut->cd();

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
    hBG_cosThetaPhi[iFrame][L]->Write();
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
    hBG_cosThetaPhi[iFrame][R]->Write();
  treeOut->Write();  
  fOut->Close();
}
