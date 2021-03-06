#include <iostream>
#include <sstream>

//J/Psi common vars
#include "calcPol.C"

//#include "MakeClassTry.h"
//Fitting routine
//#include "CompositeModelBuilder.h"

// RooFit Includes
//#include "RooAddPdf.h"
//#include "RooExtendPdf.h"
#include "RooRealVar.h"
//#include "RooPlot.h"
//#include "RooArgSet.h"
#include "RooDataSet.h"
//#include "RooFitResult.h"

//ROOT includes
//#include "TH2F.h"
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <TLorentzVector.h>

//#include "TRoot.h"

//#include "TDirectory.h"


int main(int argc, char** argv) {
  using namespace RooFit;
//	TFile* fIn = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root");//,"READ");
	TFile* fIn = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/TTree_MCprompt_Jpsi_Fall10.root");//,"READ");


    TTree* fChain = (TTree*)fIn->Get("data");



    TLorentzVector * JpsiP;
    TLorentzVector * muPosP;
    TLorentzVector * muNegP;



    Double_t phiFolded[jpsi::kNbFrames];
    Double_t thetaAdjusted[jpsi::kNbFrames];

    //TTree * fChain=tree;

     Long64_t nentries = fChain->GetEntries();
     cout<<nentries<<endl;
     Long64_t nbytes = 0;
     Int_t nb = 0;

     Double_t JpsiMass;
     Double_t JpsiPt;
     Double_t JpsiRap;

     Double_t muPosPt;
     Double_t muPosEta;
     Double_t muPosPhi;

     Double_t muNegPt;
     Double_t muNegEta;
     Double_t muNegPhi;

     Int_t runNb;
     Double_t Jpsict;
     Double_t JpsictErr;

     Int_t           HLT_Mu0_Track0_Jpsi;
     Int_t           HLT_Mu3_Track0_Jpsi;
     Int_t           HLT_Mu5_Track0_Jpsi;
     Int_t           HLT_Mu0_TkMu0_Jpsi;
     Int_t           HLT_Mu3_TkMu0_Jpsi;
     Int_t           HLT_Mu5_TkMu0_Jpsi;
     Int_t           HLT_Mu0_TkMu0_OST_Jpsi;
     Int_t           HLT_Mu3_TkMu0_OST_Jpsi;
     Int_t           HLT_Mu5_TkMu0_OST_Jpsi;
     Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;
     Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;
     Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;
     Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
     Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;
     Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;
     Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;
     Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;
     Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;
     Int_t           HLT_Mu3_Track3_Jpsi;
     Int_t           HLT_Mu3_Track3_Jpsi_v2;
     Int_t           HLT_Mu3_Track5_Jpsi_v1;
     Int_t           HLT_Mu3_Track5_Jpsi_v2;
     Int_t           HLT_DoubleMu0;
     Int_t           HLT_DoubleMu0_Quarkonium_v1;
     Int_t           HLT_DoubleMu0_Quarkonium_LS_v1;
     Int_t           HLT_L1DoubleMuOpen;
     Int_t           HLT_L1DoubleMuOpen_Tight;
     Int_t           HLT_DoubleMu3;
     Int_t           HLT_Mu3;
     Int_t           HLT_Mu5;
     Int_t           HLT_Mu7;
     Int_t           HLT_Mu9;
     Int_t           HLT_Mu11;
 	cout << "Entry 100 :" << fChain->GetEntry(100) << endl;

     TBranch        *b_Jpsict;
     TBranch        *b_JpsictErr;
     TBranch        *b_runNb;
     TBranch        *b_HLT_Mu0_Track0_Jpsi;   //!
     TBranch        *b_HLT_Mu3_Track0_Jpsi;   //!
     TBranch        *b_HLT_Mu5_Track0_Jpsi;   //!
     TBranch        *b_HLT_Mu0_TkMu0_Jpsi;   //!
     TBranch        *b_HLT_Mu3_TkMu0_Jpsi;   //!
     TBranch        *b_HLT_Mu5_TkMu0_Jpsi;   //!
     TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi;   //!
     TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi;   //!
     TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi;   //!
     TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;   //!
     TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;   //!
     TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;   //!
     TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;   //!
     TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;   //!
     TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;   //!
     TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;   //!
     TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;   //!
     TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;   //!
     TBranch        *b_HLT_Mu3_Track3_Jpsi;   //!
     TBranch        *b_HLT_Mu3_Track3_Jpsi_v2;   //!
     TBranch        *b_HLT_Mu3_Track5_Jpsi_v1;   //!
     TBranch        *b_HLT_Mu3_Track5_Jpsi_v2;   //!
     TBranch        *b_HLT_DoubleMu0;   //!
     TBranch        *b_HLT_DoubleMu0_Quarkonium_v1;   //!
     TBranch        *b_HLT_DoubleMu0_Quarkonium_LS_v1;   //!
     TBranch        *b_HLT_L1DoubleMuOpen;   //!
     TBranch        *b_HLT_L1DoubleMuOpen_Tight;   //!
     TBranch        *b_HLT_DoubleMu3;   //!
     TBranch        *b_HLT_Mu3;   //!
     TBranch        *b_HLT_Mu5;   //!
     TBranch        *b_HLT_Mu7;   //!
     TBranch        *b_HLT_Mu9;   //!
     TBranch        *b_HLT_Mu11;   //!
     TBranch        *b_JpsiP;   //!
     TBranch        *b_muPosP;   //!
     TBranch        *b_muNegP;   //!

 	cout << "Entry 100 :" << fChain->GetEntry(100) << endl;

   	cout << "Entry 100 :" << fChain->GetEntry(100) << endl;

  	 fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
  	 fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
   	 fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
     fChain->SetBranchAddress("HLT_Mu0_Track0_Jpsi", &HLT_Mu0_Track0_Jpsi, &b_HLT_Mu0_Track0_Jpsi);
     fChain->SetBranchAddress("HLT_Mu3_Track0_Jpsi", &HLT_Mu3_Track0_Jpsi, &b_HLT_Mu3_Track0_Jpsi);
     fChain->SetBranchAddress("HLT_Mu5_Track0_Jpsi", &HLT_Mu5_Track0_Jpsi, &b_HLT_Mu5_Track0_Jpsi);
     fChain->SetBranchAddress("HLT_Mu0_TkMu0_Jpsi", &HLT_Mu0_TkMu0_Jpsi, &b_HLT_Mu0_TkMu0_Jpsi);
     fChain->SetBranchAddress("HLT_Mu3_TkMu0_Jpsi", &HLT_Mu3_TkMu0_Jpsi, &b_HLT_Mu3_TkMu0_Jpsi);
     fChain->SetBranchAddress("HLT_Mu5_TkMu0_Jpsi", &HLT_Mu5_TkMu0_Jpsi, &b_HLT_Mu5_TkMu0_Jpsi);
     fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi", &HLT_Mu0_TkMu0_OST_Jpsi, &b_HLT_Mu0_TkMu0_OST_Jpsi);
     fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi", &HLT_Mu3_TkMu0_OST_Jpsi, &b_HLT_Mu3_TkMu0_OST_Jpsi);
     fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi", &HLT_Mu5_TkMu0_OST_Jpsi, &b_HLT_Mu5_TkMu0_OST_Jpsi);
     fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1);
     fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1);
     fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1);
     fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2);
     fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2);
     fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2);
     fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3);
     fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3);
     fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3);
     fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi", &HLT_Mu3_Track3_Jpsi, &b_HLT_Mu3_Track3_Jpsi);
     fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v2", &HLT_Mu3_Track3_Jpsi_v2, &b_HLT_Mu3_Track3_Jpsi_v2);
     fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v1", &HLT_Mu3_Track5_Jpsi_v1, &b_HLT_Mu3_Track5_Jpsi_v1);
     fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v2", &HLT_Mu3_Track5_Jpsi_v2, &b_HLT_Mu3_Track5_Jpsi_v2);
     fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
     fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &HLT_DoubleMu0_Quarkonium_v1, &b_HLT_DoubleMu0_Quarkonium_v1);
     fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_LS_v1", &HLT_DoubleMu0_Quarkonium_LS_v1, &b_HLT_DoubleMu0_Quarkonium_LS_v1);
     fChain->SetBranchAddress("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, &b_HLT_L1DoubleMuOpen);
     fChain->SetBranchAddress("HLT_L1DoubleMuOpen_Tight", &HLT_L1DoubleMuOpen_Tight, &b_HLT_L1DoubleMuOpen_Tight);
     fChain->SetBranchAddress("HLT_DoubleMu3", &HLT_DoubleMu3, &b_HLT_DoubleMu3);
     fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
     fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
     fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
     fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
     fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);

 	cout << "Entry 100 :" << fChain->GetEntry(100) << endl;

    RooRealVar* Jpsi_Mass = new RooRealVar("JpsiMass","J/psi mass",0,50,"GeV/c^{2}");
    RooRealVar* Jpsi_Pt = new RooRealVar("JpsiPt","J/psi pt",0,100,"GeV/c");
    RooRealVar* Jpsi_Rap = new RooRealVar("JpsiRap","J/psi eta",-10,10);
    RooRealVar* Jpsi_ct = new RooRealVar("Jpsict","J/psi ctau",-10,10,"mm");
    RooRealVar* Jpsi_ctErr = new RooRealVar("JpsictErr","J/psi ctau error",-10,10,"mm");

    RooRealVar* costh_CS = new RooRealVar("costh_CS","costh_CS",-1.,1.);
    RooRealVar* phi_CS = new RooRealVar("phi_CS","phi_CS",-180.,180.,"deg");

    RooRealVar* costh_HX = new RooRealVar("costh_HX","costh_HX",-1.,1.);
    RooRealVar* phi_HX = new RooRealVar("phi_HX","phi_HX",-180.,180.,"deg");

    RooRealVar* costh_PHX = new RooRealVar("costh_PHX","costh_PHX",-1.,1.);
    RooRealVar* phi_PHX = new RooRealVar("phi_PHX","phi_PHX",-180.,180.,"deg");

    RooRealVar* costh_sGJ = new RooRealVar("costh_sGJ","costh_sGJ",-1.,1.);
    RooRealVar* phi_sGJ = new RooRealVar("phi_sGJ","phi_sGJ",-180.,180.,"deg");

    RooRealVar* costh_GJ1 = new RooRealVar("costh_GJ1","costh_GJ1",-1.,1.);
    RooRealVar* phi_GJ1 = new RooRealVar("phi_GJ1","phi_GJ1",-180.,180.,"deg");

    RooRealVar* costh_GJ2 = new RooRealVar("costh_GJ2","costh_GJ2",-1.,1.);
    RooRealVar* phi_GJ2 = new RooRealVar("phi_GJ2","phi_GJ2",-180.,180.,"deg");

    RooRealVar* muPos_Pt = new RooRealVar("muPos_Pt","muPos_Pt",-1000.,1000.);
    RooRealVar* muNeg_Pt = new RooRealVar("muNeg_Pt","muNeg_Pt",-1000.,1000.);
    RooRealVar* muPos_Eta = new RooRealVar("muPos_Eta","muPos_Eta",-10.,10.);
    RooRealVar* muNeg_Eta = new RooRealVar("muNeg_Eta","muNeg_Eta",-10.,10.);
    RooRealVar* muPos_Phi = new RooRealVar("muPos_Phi","muPos_Phi",-4.,4.);
    RooRealVar* muNeg_Phi = new RooRealVar("muNeg_Phi","muNeg_Phi",-4.,4.);





    RooArgList varlist(*Jpsi_Mass,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap);
    varlist.add(*Jpsi_ctErr);
    varlist.add(*costh_CS); varlist.add(*phi_CS); //varlist.add(*phi_CS_prime);//varlist.add(*costh_CS_prime);
    varlist.add(*costh_HX); varlist.add(*phi_HX);// varlist.add(*costh_HX_prime);  varlist.add(*phi_HX_prime);
    varlist.add(*costh_PHX); varlist.add(*phi_PHX);// varlist.add(*costh_PHX_prime); varlist.add(*phi_PHX_prime);
    varlist.add(*costh_sGJ); varlist.add(*phi_sGJ);// varlist.add(*costh_sGJ_prime); varlist.add(*phi_sGJ_prime);
    varlist.add(*costh_GJ1); varlist.add(*phi_GJ1);// varlist.add(*costh_GJ1_prime); varlist.add(*phi_GJ1_prime);
    varlist.add(*costh_GJ2); varlist.add(*phi_GJ2);// varlist.add(*costh_GJ2_prime); varlist.add(*phi_GJ2_prime);
    varlist.add(*muPos_Pt); varlist.add(*muPos_Eta); varlist.add(*muPos_Phi); varlist.add(*muNeg_Pt); varlist.add(*muNeg_Eta); varlist.add(*muNeg_Phi);

    RooDataSet* data = new RooDataSet("data","A sample",varlist);
    RooDataSet* data_unf = new RooDataSet("data","A sample",varlist);

	char outputfilename[200];
	sprintf(outputfilename,"PrepareFitTreeCuts.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");

    int FiducialCutCountPos = 0;
    int FiducialCutCountNeg = 0;
    int TriggerLossCount = 0;
    int RunNbLossCount = 0;

    cout << "TChain " << fChain << endl;
    //cout << "TChain->GetEntry(0)" << fChain->Show() << endl;
    cout << "TChain->GetEntry(0)" << fChain->GetEntry(0) << endl;

    fChain->SetBranchAddress("JpsiP", &JpsiP,&b_JpsiP);
    fChain->SetBranchAddress("muPosP", &muPosP,&b_muPosP);
    fChain->SetBranchAddress("muNegP", &muNegP,&b_muNegP);

    for (Long64_t jentry=1000; jentry<10000;jentry++) {
      		    	 if(jentry % 1000 == 0) printf("event %d\n", (Int_t) jentry);
      		 //  	cout << "Entry 100 :" << fChain->GetEntry(100) << endl;



      		    	 cout<<" -- 1 where are u, you shitty bug?"<<endl;
      		   // 	cout<<" -- 1a where are u, you shitty bug?"<< fChain->GetEntry(jentry) << endl;
      		    	 nb = fChain->GetEntry(jentry);
      		    	 cout<<"-- 2 where are u, you shitty bug?"<<endl;
      		    	 nbytes += nb;
      		    	 cout<<"-- 3 where are u, you shitty bug?"<<endl;

//////////// kinematical fiducial cuts ///////////////////////////////////

/*      	  Double_t etaMuPos = muPosP->Eta();
      	  Double_t etaMuNeg = muNegP->Eta();
      	  Double_t pTMuPos = muPosP->Pt();
      	  Double_t pTMuNeg = muNegP->Pt();
      	  Double_t pMuPos = muPosP->P();
      	  Double_t pMuNeg = muNegP->P();

      	  if((fabs(etaMuPos) < jpsi::etaPS[0] && pTMuPos < jpsi::pTMuMin[0]) || //mid-rapidity cut
      	     (fabs(etaMuPos) > jpsi::etaPS[0] && fabs(etaMuPos) < jpsi::etaPS[1] && pMuPos < jpsi::pMuMin[1]) ||
      	     (fabs(etaMuPos) > jpsi::etaPS[1] && fabs(etaMuPos) < jpsi::etaPS[2] && pTMuPos < jpsi::pTMuMin[1]))
      	    {FiducialCutCountPos++; continue;}
      	  //(b) on the negative muon
      	  if((fabs(etaMuNeg) < jpsi::etaPS[0] && pTMuNeg < jpsi::pTMuMin[0]) || //mid-rapidity cut
      	     (fabs(etaMuNeg) > jpsi::etaPS[0] && fabs(etaMuNeg) < jpsi::etaPS[1] && pMuNeg < jpsi::pMuMin[1]) ||
      	     (fabs(etaMuNeg) > jpsi::etaPS[1] && fabs(etaMuNeg) < jpsi::etaPS[2] && pTMuNeg < jpsi::pTMuMin[1]))
      	    {FiducialCutCountNeg++; continue;}
*/
////////////// TRIGGER requirements //////////////////////////////////////////////////////////////////
          cout<<"where are u?"<<endl;


      	  if(runNb >= 140116 && runNb <= 144114){
      	    if(HLT_Mu0_TkMu0_Jpsi != 1) {TriggerLossCount++;continue;}
      	  }
      	  else if(runNb >= 146428 && runNb <= 148058){
      	    if(HLT_Mu0_TkMu0_OST_Jpsi != 1) {TriggerLossCount++;continue;}
      	  }
      	  else if(runNb >= 148819 && runNb <= 149182){
      	    if(HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2 != 1) {TriggerLossCount++;continue;}
      	  }
      	  else if(runNb >= 149291 && runNb <= 149442){
      	    if(HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3 != 1) {TriggerLossCount++;continue;}
      	  }
      	  else{
      	    printf("rejecting events in run %d\n", runNb);
      	  RunNbLossCount++;
      	  continue;
      	  }


////////////// calculate angular distributions //////////////////////////////////////////////////////////
          cout<<"where are u?"<<endl;

      	   calcPol(*muPosP,*muNegP);

///// Folding:

      	    for (int iFrame=0; iFrame<jpsi::kNbFrames;iFrame++) {

       	phiFolded [iFrame] = thisPhi[iFrame];
      	thetaAdjusted [iFrame] = thisCosTh[iFrame];
      	  if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.)
      	    phiFolded[iFrame]*= -1;
      	  else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
      	    phiFolded[iFrame] = 180. - thisPhi[iFrame];
      	    thetaAdjusted [iFrame]*= -1;
      	  }
      	  else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
      	    phiFolded [iFrame]= 180. + thisPhi[iFrame];
      	    thetaAdjusted [iFrame]*= -1;
      	  }

      	    }


///////////////// get values from tree and store in dataset ///////////////////////////////////////////////////
      	    cout<<"where are u?"<<endl;

/*      	   JpsiMass=JpsiP->M(); JpsiPt=JpsiP->Pt(); JpsiRap=JpsiP->Eta();

      	   muPosPt=muPosP->Pt(); muNegPt=muNegP->Pt(); muPosEta=muPosP->Eta(); muNegEta=muNegP->Eta(); muPosPhi=muPosP->Phi(); muNegPhi=muNegP->Phi();

      	   muPos_Pt->setVal(muPosPt); muPos_Eta->setVal(muPosEta); muPos_Phi->setVal(muPosPhi); muNeg_Pt->setVal(muNegPt); muNeg_Eta->setVal(muNegEta); muNeg_Phi->setVal(muNegPhi);
*/
      	   Jpsi_Mass->setVal(JpsiMass); Jpsi_Pt->setVal(JpsiPt); Jpsi_Rap->setVal(JpsiRap); Jpsi_ct->setVal(Jpsict); Jpsi_ctErr->setVal(JpsictErr);

           costh_CS->setVal(thisCosTh[jpsi::CS]); phi_CS->setVal(thisPhi[jpsi::CS]);

           costh_HX->setVal(thisCosTh[jpsi::HX]); phi_HX->setVal(thisPhi[jpsi::HX]);

           costh_PHX->setVal(thisCosTh[jpsi::PHX]); phi_PHX->setVal(thisPhi[jpsi::PHX]);

           costh_sGJ->setVal(thisCosTh[jpsi::sGJ]); phi_sGJ->setVal(thisPhi[jpsi::sGJ]);

           costh_GJ1->setVal(thisCosTh[jpsi::GJ1]); phi_GJ1->setVal(thisPhi[jpsi::GJ1]);

           costh_GJ2->setVal(thisCosTh[jpsi::GJ2]); phi_GJ2->setVal(thisPhi[jpsi::GJ2]);

           data_unf->add(varlist);




          costh_CS->setVal(thetaAdjusted [0]); phi_CS->setVal(phiFolded [0]);

          costh_HX->setVal(thetaAdjusted [1]); phi_HX->setVal(phiFolded [1]);

          costh_PHX->setVal(thetaAdjusted [2]); phi_PHX->setVal(phiFolded [2]);

          costh_sGJ->setVal(thetaAdjusted [3]); phi_sGJ->setVal(phiFolded [3]);

          costh_GJ1->setVal(thetaAdjusted [4]); phi_GJ1->setVal(phiFolded [4]);

          costh_GJ2->setVal(thetaAdjusted [5]); phi_GJ2->setVal(phiFolded [5]);

          data->add(varlist);

       }

//    RooDataSet * reddata=data->reduce("Jpsi_Mass > 2.7 && Jpsi_Mass < 3.5");


	fprintf(outputFile, "Pos muons cut: %d\n",FiducialCutCountPos);
	fprintf(outputFile, "Neg muons cut: %d\n",FiducialCutCountNeg);
	fprintf(outputFile, "Trigger cuts: %d\n",TriggerLossCount);
	fprintf(outputFile, "RunNb cuts: %d\n",RunNbLossCount);

    TTree* datatree = (TTree*)data->tree();
    TTree* datatree_unf = (TTree*)data_unf->tree();

//    data->SaveAs("/scratch/knuenz/Polarization/RootInput/RooDataSet_pol_final_MCprompt_Jpsi_Fall10_folded.root");
//    data_unf->SaveAs("/scratch/knuenz/Polarization/RootInput/RooDataSet_pol_final_MCprompt_Jpsi_Fall10.root");

    datatree->SaveAs("/scratch/knuenz/Polarization/RootInput/TTree_pol_final_MCprompt_Jpsi_Fall10_folded.root");
    datatree_unf->SaveAs("/scratch/knuenz/Polarization/RootInput/TTree_pol_final_MCprompt_Jpsi_Fall10.root");

    datatree->Print();



return 0;

}
