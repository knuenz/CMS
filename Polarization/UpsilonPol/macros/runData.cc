#ifndef __CINT__
#endif

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TChain.h"
#include "rootIncludes.inc"
#include "PolData.C"

void BookHistosReco();
void WriteHistosReco(Char_t *fNameOut);
//==========================================
int main(int argc, char** argv){
  
	Char_t *fNameOut = "tmpFiles/selEvents_data_Ups.root";
	Bool_t rejectCowboys = kTRUE;

	Char_t *inputTree1 = "Default";
	Char_t *inputTree2 = "Default";
	Char_t *inputTree3 = "Default";
	Char_t *inputTree4 = "Default";
	int FidCuts=999;

	int inputTrees=0;




	  for( int i=0;i < argc; ++i ) {
	    if(std::string(argv[i]).find("rejectCowboys=false") != std::string::npos) {rejectCowboys=kFALSE; cout<<"reject Cowboys"<<endl;}
	    if(std::string(argv[i]).find("inputTree1") != std::string::npos) {inputTrees++; char* inputTree1char = argv[i]; char* inputTree1char2 = strtok (inputTree1char, "="); inputTree1 = inputTree1char2; cout<<"inputTree1 = "<<inputTree1<<endl;}
	    if(std::string(argv[i]).find("inputTree2") != std::string::npos) {inputTrees++; char* inputTree2char = argv[i]; char* inputTree2char2 = strtok (inputTree2char, "="); inputTree2 = inputTree2char2; cout<<"inputTree2 = "<<inputTree2<<endl;}
	    if(std::string(argv[i]).find("inputTree3") != std::string::npos) {inputTrees++; char* inputTree3char = argv[i]; char* inputTree3char2 = strtok (inputTree3char, "="); inputTree3 = inputTree3char2; cout<<"inputTree3 = "<<inputTree3<<endl;}
	    if(std::string(argv[i]).find("inputTree4") != std::string::npos) {inputTrees++; char* inputTree4char = argv[i]; char* inputTree4char2 = strtok (inputTree4char, "="); inputTree4 = inputTree4char2; cout<<"inputTree4 = "<<inputTree4<<endl;}
	    if(std::string(argv[i]).find("FidCuts") != std::string::npos) {char* FidCutschar = argv[i]; char* FidCutschar2 = strtok (FidCutschar, "p"); FidCuts = atof(FidCutschar2); cout<<"FidCuts = "<<FidCuts<<endl;}
	    }

	  cout<<"Number of Input Trees = "<<inputTrees<<endl;

  TChain *chain = new TChain("data");
  chain->Add(inputTree1);
  if(inputTrees>1) chain->Add(inputTree2);
  if(inputTrees>2) chain->Add(inputTree3);
  if(inputTrees>3) chain->Add(inputTree4);

  TTree *tree = chain;
  
  TFile *fOut = fOut = new TFile(fNameOut, "RECREATE");

  PolData treeReco(tree);
  BookHistosReco();
  printf("after booking of histo\n");
  Int_t selDimuType = 4; //0...only GG, 1... only GT, 2... only TT, 3...GG+GT, 4...GG+GT+TT
  treeReco.Loop(selDimuType, rejectCowboys, FidCuts);
  printf("writing out the histograms\n");
  WriteHistosReco(fNameOut);

  fOut->Close();

  return 0;
}
//==========================================
void BookHistosReco(){

  //mass
  Int_t nBinsMass = 320;
  Double_t massMin = 8.4, massMax = 11.6;
  //pt
  Int_t nBinsPt = 1000;
  Double_t pTMin = 0., pTMaxOnia = 100.;
  //rap
  Int_t nBinsRap = 100;
  Double_t rapMin = -2.5, rapMax = 2.5;

  Char_t name[100], title[300];
  //statistics
  Reco_StatEv = new TH1F("Reco_StatEv", "", 12, 0., 12.);

  //reconstruction variables for the Onia
  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
      //Mass:
      sprintf(name, "Reco_Onia_mass_rap%d_pT%d", iRapBin, iPTBin);
      sprintf(title, ";M [GeV/c^{2}]");
      Reco_Onia_mass[iPTBin][iRapBin] = new TH1F(name, title, nBinsMass, massMin, massMax);
      Reco_Onia_mass[iPTBin][iRapBin]->Sumw2();
    }
  }

  sprintf(name, "Reco_Onia_rap_pt");
  sprintf(title, ";y(#mu#mu);p_{T}^{#mu#mu} [GeV/c]");
  Reco_Onia_rap_pT = new TH2F(name, title, nBinsRap,rapMin,rapMax, nBinsPt,pTMin,pTMaxOnia);
  Reco_Onia_rap_pT->Sumw2();

  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
    //pT
    sprintf(name, "Reco_Onia_pt_rap%d", iRapBin);
    sprintf(title, ";p_{T}^{#mu#mu} [GeV/c]");
    Reco_Onia_pt[iRapBin]  = new TH1F(name, title, nBinsPt,pTMin,pTMaxOnia);
    Reco_Onia_pt[iRapBin]->Sumw2();
  }
  for(int iPTBin = 0; iPTBin < onia::kNbPTMaxBins+1; iPTBin++){
    //rap
    sprintf(name, "Reco_Onia_rap_pT%d", iPTBin);
    sprintf(title, ";y(#mu#mu)");
    Reco_Onia_rap[iPTBin]  = new TH1F(name, title, nBinsRap, rapMin,rapMax);
    Reco_Onia_rap[iPTBin]->Sumw2();
  }

  //prepare the branches for the output tree
  treeOut = new TTree ("selectedData", "selected events");
  lepP = new TLorentzVector();
  lepN = new TLorentzVector();
  treeOut->Branch("lepP", "TLorentzVector", &lepP);
  treeOut->Branch("lepN", "TLorentzVector", &lepN);

  // //polarization histos:
  // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

  //   //pT and y double differential pol histos:
  //   for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  // 	//2D histo:
  // 	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
  // 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin] = new TH2F(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
  // 								   onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
  // 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Sumw2();
  //     }
  //   }
  //   //pT and y double differential pol histos: FWD and BWD rapidities separately
  //   for(int iRapBin = 0; iRapBin < 2*onia::kNbRapBins; iRapBin++){
  //     Int_t matchRapBin = fabs(onia::kNbRapForPTBins - iRapBin);
  //     if(iRapBin >= onia::kNbRapForPTBins) matchRapBin += 1;
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[matchRapBin]+1; iPTBin++){
  // 	sprintf(name, "Reco_Onia_cosTh_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos#theta_{%s}", onia::frameLabel[iFrame]);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cosThPol] = new TH1F(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cosThPol]->Sumw2();
  // 	//
  // 	sprintf(name, "Reco_Onia_phi_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";#phi_{%s} [deg]", onia::frameLabel[iFrame]);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::phiPol] = new TH1F(name, title, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::phiPol]->Sumw2();
  // 	//
  // 	sprintf(name, "Reco_Onia_cos2Phi_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos(2#phi_{%s})", onia::frameLabel[iFrame]);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cos2PhiPol] = new TH1F(name, title, onia::kNbBinsCos2Phi, onia::cos2PhiMin, onia::cos2PhiMax);
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cos2PhiPol]->Sumw2();
  // 	//2D histo:
  // 	sprintf(name, "Reco2D_Onia_%s_pT%d_rap%d_FWDBWD", onia::frameLabel[iFrame], iPTBin, iRapBin);
  // 	sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
  // 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin] = new TH2F(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
  // 								   onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
  // 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin]->Sumw2();
  //     }
  //   }
  // }

  // for(int iRapBin = 1; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 1; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  //     //checking the rotation angle between HX and CS:
  //     sprintf(name, "Reco_hDelta_pT%d_rap%d", iPTBin, iRapBin);
  //     hDelta[iPTBin][iRapBin] = new TH1F(name, ";#delta(HX --> CS) [#circ]", 180, 0., 180.);
  //     hDelta[iPTBin][iRapBin]->Sumw2();
  //     sprintf(name, "Reco_hSin2Delta_pT%d_rap%d", iPTBin, iRapBin);
  //     hSin2Delta[iPTBin][iRapBin] = new TH1F(name, ";sin^{2}#delta(HX --> CS)", 100, 0., 1.);
  //     hSin2Delta[iPTBin][iRapBin]->Sumw2();
  //   }
  // }
}

//==========================================
void WriteHistosReco(Char_t *fNameOut){

  treeOut->Write();

  Reco_StatEv->Write();

  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
    for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){  
      Reco_Onia_mass[iPTBin][iRapBin]->Write();
    }
  }

  Reco_Onia_rap_pT->Write();
  for(int iPTBin = 0; iPTBin < onia::kNbPTMaxBins+1; iPTBin++)
    Reco_Onia_rap[iPTBin]->Write();
  for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++)
    Reco_Onia_pt[iRapBin]->Write();


  // //polarization histos: Reco
  // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
  //   for(int iRapBin = 0; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  // 	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][onia::cosThPol]->Write();
  // 	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][onia::phiPol]->Write();
  // 	Reco_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin][onia::cos2PhiPol]->Write();
  // 	Reco2D_Onia_pol_pT_rap[iFrame][iPTBin][iRapBin]->Write();
  //     }
  //   }
  // }
  // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
  //   for(int iRapBin = 0; iRapBin < 2*onia::kNbRapBins; iRapBin++){
  //     Int_t matchRapBin = fabs(onia::kNbRapForPTBins - iRapBin);
  //     if(iRapBin >= onia::kNbRapForPTBins) matchRapBin += 1;
  //     for(int iPTBin = 0; iPTBin < onia::kNbPTBins[matchRapBin]+1; iPTBin++){
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cosThPol]->Write();
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::phiPol]->Write();
  // 	Reco_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin][onia::cos2PhiPol]->Write();
  // 	Reco2D_Onia_pol_pT_rap_FWDBWD[iFrame][iPTBin][iRapBin]->Write();
  //     }
  //   }
  // }

  // for(int iRapBin = 1; iRapBin < onia::kNbRapForPTBins+1; iRapBin++){
  //   for(int iPTBin = 1; iPTBin < onia::kNbPTBins[iRapBin]+1; iPTBin++){
  //     hDelta[iPTBin][iRapBin]->Write();
  //     hSin2Delta[iPTBin][iRapBin]->Write();
  //   }
  // }

  // for(int iCut = 0; iCut < 4; iCut++){
  //   Reco_muHLT_pT_eta[iCut]->Write();
  //   Reco_muHLT_p_eta[iCut]->Write();
  //   Reco_muTM_pT_eta[iCut]->Write();
  //   Reco_muTM_p_eta[iCut]->Write();
  // }
}
