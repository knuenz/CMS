#include "rootIncludes.inc"
#include "calcPol.C"
#include "TH3D.h"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};

void TrimEventContent(Int_t iRapBin = 1, 
		      Int_t iPTBin = 1,
		      Double_t fracL = 0.5, Double_t nSigma = 2., 
		      Int_t nUpsState=0//[0]... 1S, [1]... 2S, [2]... 3S
		      ){


  printf("\n\n\nfracL = %1.1f, nSigma = %1.1f, iState = %d, rap %d, pT %d\n", fracL, nSigma, nUpsState, iRapBin, iPTBin);

  Char_t name[100], title[100];
  Char_t fileNameIn[100];
  sprintf(fileNameIn, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
  //==============================
  //read inputs from input file:
  TFile *fIn = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");
  TH2D *hBG_cosThetaPhiLR[onia::kNbFrames][2];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
    hBG_cosThetaPhiLR[iFrame][L] = (TH2D *) gDirectory->Get(name);
    sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
    hBG_cosThetaPhiLR[iFrame][R] = (TH2D *) gDirectory->Get(name);
  }
  //==============================

  //definition of output variables 
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "AllStates_%1.2fSigma_FracLSB%dPercent/data_%dSUps_rap%d_pT%d.root", nSigma, int(fracL*100), nUpsState+1, iRapBin, iPTBin);
  TFile *fOut = new TFile(fileNameOut, "RECREATE");
  gStyle->SetPadRightMargin(0.2);
  TTree *treeOut = treeIn->CloneTree(0);
  // treeOut->SetName("data");
  TH2D *hBG_cosThetaPhi[onia::kNbFrames];
  // TH2D *hBG_cosThetaPhiSignal[onia::kNbFrames];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    // //book the histo for the signal
    // sprintf(name, "total_%s", onia::frameLabel[iFrame]);
    // sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    // hBG_cosThetaPhiSignal[iFrame] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
    // 					  onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    // hBG_cosThetaPhiSignal[iFrame]->Sumw2();
    //copy the L and R sideband histos into one output BG histogram
    hBG_cosThetaPhiLR[iFrame][L]->Scale(fracL/hBG_cosThetaPhiLR[iFrame][L]->Integral());
    hBG_cosThetaPhiLR[iFrame][R]->Scale((1.-fracL)/hBG_cosThetaPhiLR[iFrame][R]->Integral());
    sprintf(name, "background_costhphi%s", onia::frameLabel[iFrame]);
    hBG_cosThetaPhi[iFrame] = (TH2D *) hBG_cosThetaPhiLR[iFrame][L]->Clone(name);
    hBG_cosThetaPhi[iFrame]->Add(hBG_cosThetaPhiLR[iFrame][R]);
  }

  //==========================================================
  //reading fit parameters to establish signal mass window
  //as well as the L and R sideband window for the 3D BG histo
  //==========================================================
  fIn->cd();
  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
  TF1 *fUps[kNbSpecies], *fBG = 0;
  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
  treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
  treeFitPar->SetBranchAddress("fBG", &fBG);
  treeFitPar->LoadTree(0);
  treeFitPar->GetEntry(0);

  Double_t mass[kNbSpecies], sigma[kNbSpecies];
  for(int iState = 0; iState < kNbSpecies; iState++){
    mass[iState] = fUps[iState]->GetParameter(1);
    sigma[iState] = fUps[iState]->GetParameter(2);
  }
  printf("1S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS1S], sigma[UPS1S]);
  printf("2S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS2S], sigma[UPS2S]);
  printf("3S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS3S], sigma[UPS3S]);
  Double_t poleMass = mass[nUpsState], massMin, massMax;
  massMin = poleMass - nSigma*sigma[nUpsState];
  massMax = poleMass + nSigma*sigma[nUpsState];

  printf("--> signal mass window: %1.3f < M < %1.3f GeV\n", massMin, massMax);

  //calculate the fraction of BG under the signal mass window
  Double_t nBG = fBG->Integral(massMin, massMax);
  Double_t nSignal = fUps[nUpsState]->Integral(massMin, massMax);
  Double_t fracBG = nBG / (nBG + nSignal);
  sprintf(name, ";;fraction of BG in %1.1f sigma window", nSigma);
  TH1D *hFracBG = new TH1D("background_fraction", name, 1, 0., 1.);
  hFracBG->SetBinContent(1, fracBG);

  //calculate the L and R mass windows:
  Double_t massMinBG[2], massMaxBG[2];
  massMinBG[L] = onia::massMinL;
  massMaxBG[L] = mass[UPS1S] - onia::nSigmaL*sigma[UPS1S];
  massMinBG[R] = mass[UPS3S] + onia::nSigmaR*sigma[UPS3S];
  massMaxBG[R] = onia::massMaxR;
  printf("--> L mass window: %1.3f < M < %1.3f GeV\n", massMinBG[L], massMaxBG[L]);
  printf("--> R mass window: %1.3f < M < %1.3f GeV\n", massMinBG[R], massMaxBG[R]);

  if(massMaxBG[L] > massMin){
    printf("the right sideband window is LARGER than the left signal window!!!!\n\n\n\n");
    exit(0);
  }
  if(massMinBG[R] < massMax){
    printf("the left sideband window is SMALLER than the right signal window!!!!\n\n\n\n");
    exit(0);
  }

  //build the 3D (pT, |y|, M) histos for the L and R mass sideband 
  TH3D *hBG_pTRapMass[2];
  hBG_pTRapMass[L] = new TH3D("hBG_pTRapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]", 
			      7, onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin],
			      2, onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin],
			      7, massMin, massMax);//*signal* mass window!
  hBG_pTRapMass[L]->Sumw2();
  //
  hBG_pTRapMass[R] = new TH3D("hBG_pTRapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]", 
			      7, onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin],
			      2, onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin],
			      7, massMin, massMax);//*signal* mass window!
  hBG_pTRapMass[R]->Sumw2();

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();

  Double_t onia_mass;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 10000 == 0)
      cout << "entry " << iEntry << " out of " << treeIn->GetEntries() << endl;

    *onia = *(lepP) + *(lepN);
    onia_mass = onia->M();

    if(onia_mass  > massMinBG[L] && onia_mass < massMaxBG[L])
      hBG_pTRapMass[L]->Fill(onia->Pt(), TMath::Abs(onia->Rapidity()), fBG->GetRandom(massMin, massMax));
    else if(onia_mass > massMinBG[R] && onia_mass < massMaxBG[R])
      hBG_pTRapMass[R]->Fill(onia->Pt(), TMath::Abs(onia->Rapidity()), fBG->GetRandom(massMin, massMax));

    if(onia_mass > massMin && onia_mass < massMax){
      treeOut->Fill(); //stores TLorenzVectors of the two muons

      // //store now the cosTheta and phi distributions of the signal window:
      // calcPol(*lepP, *lepN);

      // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
      // 	hCosThetaPhiSignal[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
    }
  }

  //now, add the L and R 3D (pT, |y|, M) histos
  TH3D *hBG_pTRapMassSum;
  hBG_pTRapMass[L]->Scale(fracL/hBG_pTRapMass[L]->Integral());
  hBG_pTRapMass[R]->Scale((1.-fracL)/hBG_pTRapMass[R]->Integral());
  sprintf(name, "background_pTrapMass");
  hBG_pTRapMassSum = (TH3D*) hBG_pTRapMass[L]->Clone(name);
  hBG_pTRapMassSum->Add(hBG_pTRapMass[R]);

  //write the output
  fOut->cd();
  treeOut->Write();

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    hBG_cosThetaPhi[iFrame]->Write();
    // hCosThetaPhiSignal[iFrame]->Write();
  }
  hBG_pTRapMassSum->Write();
  hFracBG->Write();

  fOut->Close();
  fIn->Close();
}
