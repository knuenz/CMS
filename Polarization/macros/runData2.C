#ifndef __CINT__
// #include "RooGlobalFunc.h"
// #include "RooDataSet.h"
#endif

#include "TChain.h"
#include "../interface/rootIncludes.inc"
#include "PolData2.C"

void BookHistosReco(Char_t *oniaLabel);
void WriteHistosReco(Char_t *fNameOut);
//==========================================
//usage: root 'runData.C+' or
//       root 'runData.C+("pol_MC_HLT_Mu0Track0Jpsi.root", kTRUE)' (e.g.)
//==========================================
void runData2(Char_t *fNameOut = "pol_data_HLT_Mu0TkMu0Jpsi.root",
	     Bool_t newOutputFile = kTRUE, //allows to create a new file or to append info
	     Char_t *nameDataSet = "data", //"data" or "recoData"
	     Int_t selDimuType = 4, //0...only GG, 1... only GT, 2... only TT, 3...GG+GT, 4...GG+GT+TT
	     Bool_t writeOutEvents = kFALSE, //writes out Run, LS, Ev.Nb for any J/psi candidate
	     Char_t *oniaLabel = "J/#psi"){//"Ups(1S)"

  TChain *chain = new TChain(nameDataSet);
  chain->Add("/home/hermine/CMS/Work/Polarization/Florian/7Dec2010/TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root");
  chain->Add("/home/hermine/CMS/Work/Polarization/Florian/7Dec2010/TTree_pol_noTriggerFilter_Run2010B-Nov4ReReco_v1_06Dec2010.root");
  TTree *tree = chain;

  TFile *fOut;
  if(newOutputFile) fOut = new TFile(fNameOut, "RECREATE");
  else fOut = new TFile(fNameOut, "UPDATE");

  PolData2 treeReco(tree);
  BookHistosReco(oniaLabel);
  printf("after booking of histo\n");
  treeReco.Loop(selDimuType, writeOutEvents);
  WriteHistosReco(fNameOut);

  fOut->Close();
}
//==========================================
void BookHistosReco(Char_t *oniaLabel){

  //mass
  Int_t nBinsMass = 80;
  Double_t massMin = 8.0, massMax = 12.0;
  if(strncmp(oniaLabel, "J/#psi", 6) == 0){
    nBinsMass = 160;
    massMin = 2.5;
    massMax = 4.1;
  }

  //pt
  Int_t nBinsPt = 100, nBinsPtGamma = 30;
  Double_t pTMin = 0., pTMaxOnia = 30., pTMaxGamma = 3.0;
  //energy
  Int_t nBinsE = 50, nBinsEOnia = 200;
  Double_t enMin = 0., enMax = 5., enMaxOnia = 50.;
  //phi
  Int_t nBinsPhi = 157; //314
  Double_t phiMin = -3.14, phiMax = 3.14;
  //rap
  Int_t nBinsRap = 100;
  Double_t rapMin = -2.5, rapMax = 2.5;
  //pseudo-rap
  Int_t nBinsEtaGamma = 120;
  Double_t etaMinGamma = -3.0, etaMaxGamma = 3.0;

  Char_t name[100], title[300];
  //statistics
  Reco_StatEv = new TH1F("Reco_StatEv", "", 12, 0., 12.);

  //reconstruction variables for the Onia
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
//       printf("mass and phi: rap %d, pT %d\n", iRapBin, iPTBin);
      //Mass:
      sprintf(name, "Reco_Onia_mass_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";M [GeV/c^{2}]");
      Reco_Onia_mass[iPTBin][iRapBin] = new TH1F(name, title, nBinsMass,massMin,massMax);
      Reco_Onia_mass[iPTBin][iRapBin]->Sumw2();
      //phi-lab:
      sprintf(name, "Reco_Onia_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, ";%s #phi (rad)", oniaLabel);
      Reco_Onia_phi[iPTBin][iRapBin] = new TH1F(name, title, nBinsPhi,phiMin,phiMax);
      Reco_Onia_phi[iPTBin][iRapBin]->Sumw2();
    }
  }
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    //pT
    sprintf(name, "Reco_Onia_pt_rap%d", iRapBin);
    sprintf(title, ";%s p_{T} [GeV/c]", oniaLabel);
    Reco_Onia_pt[iRapBin]  = new TH1F(name, title, nBinsPt,pTMin,pTMaxOnia);
    Reco_Onia_pt[iRapBin]->Sumw2();
  }

  //2D Onia histos:
  sprintf(name, "Reco_Onia_rap_pt");
  sprintf(title, ";%s y;p_{T} [GeV/c]", oniaLabel);
  Reco_Onia_rap_pT = new TH2F(name, title, nBinsRap,rapMin,rapMax, nBinsPt,pTMin,pTMaxOnia);
  Reco_Onia_rap_pT->Sumw2();

  //debugging histos (single Muons):
  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
      sprintf(name, "Reco_mupl_pt_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;p_{T}(#mu^{+})[GeV/c]",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      Reco_mupl_pt[iPTBin][iRapBin]  = new TH1F(name, title,nBinsPt,pTMin,pTMaxOnia);
      sprintf(name,"Reco_mupl_eta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#eta(#mu^{+})",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      Reco_mupl_eta[iPTBin][iRapBin] = new TH1F(name,title,nBinsRap,rapMin,rapMax);
      sprintf(name,"Reco_mupl_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#phi(#mu^{+})",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      Reco_mupl_phi[iPTBin][iRapBin] = new TH1F(name,title, nBinsPhi,phiMin,phiMax);

      sprintf(name,"Reco_mumi_pt_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;p_{T}(#mu^{-}) [GeV/c]",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      Reco_mumi_pt[iPTBin][iRapBin]  = new TH1F(name, title,nBinsPt,pTMin,pTMaxOnia);
      sprintf(name,"Reco_mumi_eta_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#eta(#mu^{-})",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      Reco_mumi_eta[iPTBin][iRapBin] = new TH1F(name,title,nBinsRap,rapMin,rapMax);
      sprintf(name,"Reco_mumi_phi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#phi(#mu^{-})",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      Reco_mumi_phi[iPTBin][iRapBin] = new TH1F(name,title, nBinsPhi,phiMin,phiMax);

      Reco_mupl_pt[iPTBin][iRapBin] ->Sumw2();
      Reco_mupl_eta[iPTBin][iRapBin]->Sumw2();
      Reco_mupl_phi[iPTBin][iRapBin]->Sumw2();

      Reco_mumi_pt[iPTBin][iRapBin] ->Sumw2();
      Reco_mumi_eta[iPTBin][iRapBin]->Sumw2();
      Reco_mumi_phi[iPTBin][iRapBin]->Sumw2();

      sprintf(name, "Reco_PhiPos_PhiNeg_pT%d_rap%d" , iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < |p_{T}(%s)| < %1.1f GeV/c;#phi(#mu^{-});#phi(#mu^{+})",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      hPhiPos_PhiNeg[iPTBin][iRapBin] = new TH2F(name, title, 60,-180.,180., 60,-180.,180.);
      sprintf(name, "Reco_PtPos_PtNeg_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < |p_{T}(%s)| < %1.1f GeV/c;p_{T}(#mu^{-}) [GeV/c];p_{T}(#mu^{+}) [GeV/c]",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      hPtPos_PtNeg[iPTBin][iRapBin] = new TH2F(name, title, 80, 0., 20., 80, 0., 20.);
      sprintf(name, "Reco_EtaPos_EtaNeg_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < |y(%s)| < %1.1f, %1.1f < |p_{T}(%s)| < %1.1f GeV/c;#eta(#mu^{-});#eta(#mu^{+})",
	      jpsi::rapForPTRange[iRapBin-1], oniaLabel, jpsi::rapForPTRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      hEtaPos_EtaNeg[iPTBin][iRapBin] = new TH2F(name, title, 48, -2.4, 2.4, 48, -2.4, 2.4);
    }
  }
  // for(int iRapBin = 1; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++)
  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){

      sprintf(name, "Reco_DeltaPhi_pT%d_rap%d", iPTBin, iRapBin);
      sprintf(title, "%1.1f < y(%s) < %1.1f, %1.1f < p_{T}(%s) < %1.1f GeV/c;#phi(#mu^{+}) - #phi(#mu^{-})",
	      jpsi::rapRange[iRapBin-1], oniaLabel, jpsi::rapRange[iRapBin],
	      jpsi::pTRange[iRapBin][iPTBin-1], oniaLabel, jpsi::pTRange[iRapBin][iPTBin]);
      hDeltaPhi[iPTBin][iRapBin] = new TH1F(name, title, 96, -1.6, 1.6);
      hDeltaPhi[iPTBin][iRapBin]->Sumw2();
    }
  }

  //polarization histos:
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){

    //pT and y double differential pol histos:
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s}", jpsi::frameLabel[iFrame]);
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol] = new TH1F(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax);
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Sumw2();
	//
	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol] = new TH1F(name, title, jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Sumw2();
	//
	sprintf(name, "Reco_Onia_cos2Phi_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos(2#phi_{%s})", jpsi::frameLabel[iFrame]);
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cos2PhiPol] = new TH1F(name, title, jpsi::kNbBinsCos2Phi, jpsi::cos2PhiMin, jpsi::cos2PhiMax);
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cos2PhiPol]->Sumw2();
	//2D histo:
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2F(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								   jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
    //pT and y double differential pol histos: FWD and BWD rapidities separately
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t matchRapBin = fabs(jpsi::kNbRapForPTBins - iRapBin);
      if(iRapBin >= jpsi::kNbRapForPTBins) matchRapBin += 1;
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[matchRapBin]+1; iPTBin++){
	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s}", jpsi::frameLabel[iFrame]);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::cosThPol] = new TH1F(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Sumw2();
	//
	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";#phi_{%s} [deg]", jpsi::frameLabel[iFrame]);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::phiPol] = new TH1F(name, title, jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Sumw2();
	//
	sprintf(name, "Reco_Onia_cos2Phi_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos(2#phi_{%s})", jpsi::frameLabel[iFrame]);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::cos2PhiPol] = new TH1F(name, title, jpsi::kNbBinsCos2Phi, jpsi::cos2PhiMin, jpsi::cos2PhiMax);
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::cos2PhiPol]->Sumw2();
	//2D histo:
	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d_FWDBWD", jpsi::frameLabel[iFrame], iPTBin, iRapBin);
	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", jpsi::frameLabel[iFrame], jpsi::frameLabel[iFrame]);
	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin] = new TH2F(name, title, jpsi::kNbBinsCosT, jpsi::cosTMin, jpsi::cosTMax,
								   jpsi::kNbBinsPhiPol, jpsi::phiPolMin, jpsi::phiPolMax);
	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin]->Sumw2();
      }
    }
  }

  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
      //checking the rotation angle between HX and CS:
      sprintf(name, "Reco_hDelta_pT%d_rap%d", iPTBin, iRapBin);
      hDelta[iPTBin][iRapBin] = new TH1F(name, ";#delta(HX --> CS) [#circ]", 180, 0., 180.);
      hDelta[iPTBin][iRapBin]->Sumw2();
      sprintf(name, "Reco_hSin2Delta_pT%d_rap%d", iPTBin, iRapBin);
      hSin2Delta[iPTBin][iRapBin] = new TH1F(name, ";sin^{2}#delta(HX --> CS)", 100, 0., 1.);
      hSin2Delta[iPTBin][iRapBin]->Sumw2();
    }
  }
}

//==========================================
void WriteHistosReco(Char_t *fNameOut){

  Reco_StatEv->Write();

  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
      Reco_Onia_mass[iPTBin][iRapBin]->Write();
    }
  }
  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++)
    for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++)
      Reco_Onia_phi[iPTBin][iRapBin]->Write();

  for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++)
    Reco_Onia_pt[iRapBin]->Write();

  Reco_Onia_rap_pT->Write();

  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
       //debugging histos (single muons):
       Reco_mupl_pt[iPTBin][iRapBin] ->Write();
       Reco_mupl_eta[iPTBin][iRapBin]->Write();
       Reco_mupl_phi[iPTBin][iRapBin]->Write();

       Reco_mumi_pt[iPTBin][iRapBin] ->Write();
       Reco_mumi_eta[iPTBin][iRapBin]->Write();
       Reco_mumi_phi[iPTBin][iRapBin]->Write();

       hPhiPos_PhiNeg[iPTBin][iRapBin]->Write();
       hPtPos_PtNeg[iPTBin][iRapBin]->Write();
       hEtaPos_EtaNeg[iPTBin][iRapBin]->Write();
     }
  }
  // for(int iRapBin = 1; iRapBin < 2*jpsi::kNbRapBins+1; iRapBin++)
  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++)
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++)
       hDeltaPhi[iPTBin][iRapBin]->Write();

  //polarization histos: Reco
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Write();
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Write();
	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][jpsi::cos2PhiPol]->Write();
 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }
  for(int iFrame = 0; iFrame < jpsi::kNbFrames; iFrame++){
    for(int iRapBin = 0; iRapBin < 2*jpsi::kNbRapBins; iRapBin++){
      Int_t matchRapBin = fabs(jpsi::kNbRapForPTBins - iRapBin);
      if(iRapBin >= jpsi::kNbRapForPTBins) matchRapBin += 1;
      for(int iPTBin = 0; iPTBin < jpsi::kNbPTBins[matchRapBin]+1; iPTBin++){
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::cosThPol]->Write();
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::phiPol]->Write();
	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][jpsi::cos2PhiPol]->Write();
 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin]->Write();
      }
    }
  }

  for(int iRapBin = 1; iRapBin < jpsi::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 1; iPTBin < jpsi::kNbPTBins[iRapBin]+1; iPTBin++){
      hDelta[iPTBin][iRapBin]->Write();
      hSin2Delta[iPTBin][iRapBin]->Write();
    }
  }
}
