#include <iostream>
#include <sstream>
#include <cstring>

//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"

// RooFit Includes
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPolarizationPdf.h"
#include "RooGaussian.h"
#include "RooPolarizationConstraint.h"
#include "RooProdPdf.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooNumIntConfig.h"
#include "RooCmdArg.h"
#include "RooLinkedList.h"
#include "RooAbsArg.h"
#include "RooHistFunc.h"
#include "RooFormulaVar.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TText.h"


int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

//  bool doprompt(true), dononprompt(true), dobkg(true),dopol(true),pereverr(false),newdata(false),domass(false),dolifetime(false);

 /* if(argc < 2) {
    std::cout << "Usage: extractMassShape /path/to/sample1.root /path/to/sample2.root ..." << std::endl;
    std::cout << "You must ensure that the MC samples are properly weighted so they can be combined." << std::endl;
    return 1;
  }  

  for( int i=0;i < argc; ++i ) { 
    if(std::string(argv[i]).find("--noNonPrompt") != std::string::npos) dononprompt = false;
    if(std::string(argv[i]).find("--noPrompt") != std::string::npos) doprompt = false;
    if(std::string(argv[i]).find("--noBackground") != std::string::npos) dobkg = false;
    if(std::string(argv[i]).find("--noPol") != std::string::npos) dopol = false;
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
    if(std::string(argv[i]).find("--noMass") != std::string::npos) domass = false;
    if(std::string(argv[i]).find("--noLifetime") != std::string::npos) dolifetime = false;
  }
*/



  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar JpsiVprob("JpsiVprob","JpsiVprob",0.01,100000);
  RooRealVar HLT_DoubleMu0("HLT_DoubleMu0","HLT_DoubleMu0",0,3);
  RooRealVar HLT_Mu0_TkMu0_OST_Jpsi("HLT_Mu0_TkMu0_OST_Jpsi","HLT_Mu0_TkMu0_OST_Jpsi",0,3);

  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(JpsictErr);
  varlist.add(JpsiVprob);
  varlist.add(HLT_DoubleMu0);
  varlist.add(HLT_Mu0_TkMu0_OST_Jpsi);

  bool eff(true);

  char rationame[200];
  if(eff) sprintf(rationame,"SysRatioHistsEff.root");
  else sprintf(rationame,"SysRatioHists.root");

  TFile *HistFile = new TFile("SysChecksHists.root","UPDATE");
	TFile *RatioFile = new TFile(rationame,"RECREATE");
	  TFile *EffHistFileDouble = new TFile("trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
	  TFile *EffHistFileTrk = new TFile("trigEffHistos_Mu0TkMu0_OST_7Feb2011_phiFolded_reBinned_zeroBinsCorrected.root","UPDATE");

//  TFile* fInMC = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_prep_PR_.root");
  //TFile* fInMC = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_HXth_plus1.root");

//  TTree* dataTreesMC = (TTree*)fInMC->Get("data");

//  TFile* fInData = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_prep_RunB.root");
  //TFile* fInData = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_HXth_plus1.root");

//  TTree* dataTreesData = (TTree*)fInData->Get("data");

//  RooDataSet *dataMC = new RooDataSet("dataMC","dataMC",varlist,Import(*dataTreesMC));
//  RooDataSet *dataData = new RooDataSet("dataData","dataData",varlist,Import(*dataTreesData));

//  cout<<"dataMC"<<endl;
//  dataMC->Print();
//  cout<<"dataData"<<endl;
//  dataData->Print();

  for(int yBin = 0; yBin < 2 /*jpsi::kNbRapForPTBins*/; ++yBin) {
    for(int ptBin = 5; ptBin < 8 /*jpsi::kNbPTBins[yBin+1]*/; ++ptBin) {

    	if (yBin==0 && ptBin==7) continue;

/*
    	      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1]
		<< " && JpsiMass > " << jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass 
		<< " && JpsiMass < " << jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass;


      JpsiMass.setRange("lowBand",2.7,
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow);
      JpsiMass.setRange("signalRegion",
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass,
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass);
      JpsiMass.setRange("highBand",
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh,3.5);

      std::cout << cutString.str() << std::endl;
      

	RooAbsData *thisBinMC = dataMC->reduce((cutString.str()).c_str());
	RooAbsData *thisBinData = dataData->reduce((cutString.str()).c_str());

	  cout<<"thisBinMC"<<endl;
	  thisBinMC->Print();
	  cout<<"thisBinData"<<endl;
	  thisBinData->Print();

	RooAbsData *thisBinMC_Double = thisBinMC->reduce("0.5<HLT_DoubleMu0 && HLT_DoubleMu0 < 1.5");
	RooAbsData *thisBinMC_TkMu0 = thisBinMC->reduce("0.5<HLT_Mu0_TkMu0_OST_Jpsi && HLT_Mu0_TkMu0_OST_Jpsi < 1.5");

	RooAbsData *thisBinData_Double = thisBinData->reduce("0.5<HLT_DoubleMu0 && HLT_DoubleMu0 < 1.5");
	RooAbsData *thisBinData_TkMu0 = thisBinData->reduce("0.5<HLT_Mu0_TkMu0_OST_Jpsi && HLT_Mu0_TkMu0_OST_Jpsi < 1.5");

	  cout<<"thisBinMC_Double"<<endl;
	  thisBinMC_Double->Print();
	  cout<<"thisBinMC_TkMu0"<<endl;
	  thisBinMC_TkMu0->Print();

	  cout<<"thisBinData_Double"<<endl;
	  thisBinData_Double->Print();
	  cout<<"thisBinData_TkMu0"<<endl;
	  thisBinData_TkMu0->Print();

	

	  char histName[200];
		  sprintf(histName,"thisBinMC_Double_hist_CS_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinMC_Double_hist_CS = (TH2F*)thisBinMC_Double->createHistogram(histName,costh_CS,Binning(40),YVar(phi_CS,Binning(36)));
		  sprintf(histName,"thisBinMC_Double_hist_HX_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinMC_Double_hist_HX = (TH2F*)thisBinMC_Double->createHistogram(histName,costh_HX,Binning(40),YVar(phi_HX,Binning(36)));

		  sprintf(histName,"thisBinMC_TkMu0_hist_CS_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinMC_TkMu0_hist_CS = (TH2F*)thisBinMC_TkMu0->createHistogram(histName,costh_CS,Binning(40),YVar(phi_CS,Binning(36)));
		  sprintf(histName,"thisBinMC_TkMu0_hist_HX_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinMC_TkMu0_hist_HX = (TH2F*)thisBinMC_TkMu0->createHistogram(histName,costh_HX,Binning(40),YVar(phi_HX,Binning(36)));

		  sprintf(histName,"thisBinData_Double_hist_CS_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinData_Double_hist_CS = (TH2F*)thisBinData_Double->createHistogram(histName,costh_CS,Binning(40),YVar(phi_CS,Binning(36)));
		  sprintf(histName,"thisBinData_Double_hist_HX_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinData_Double_hist_HX = (TH2F*)thisBinData_Double->createHistogram(histName,costh_HX,Binning(40),YVar(phi_HX,Binning(36)));

		  sprintf(histName,"thisBinData_TkMu0_hist_CS_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinData_TkMu0_hist_CS = (TH2F*)thisBinData_TkMu0->createHistogram(histName,costh_CS,Binning(40),YVar(phi_CS,Binning(36)));
		  sprintf(histName,"thisBinData_TkMu0_hist_HX_rap%d_pT%d",yBin+1,ptBin+1);
		TH2F * thisBinData_TkMu0_hist_HX = (TH2F*)thisBinData_TkMu0->createHistogram(histName,costh_HX,Binning(40),YVar(phi_HX,Binning(36)));



		HistFile->Add(thisBinMC_Double_hist_CS);
		HistFile->Add(thisBinMC_Double_hist_HX);
		HistFile->Add(thisBinMC_TkMu0_hist_CS);
		HistFile->Add(thisBinMC_TkMu0_hist_HX);
		HistFile->Add(thisBinData_Double_hist_CS);
		HistFile->Add(thisBinData_Double_hist_HX);
		HistFile->Add(thisBinData_TkMu0_hist_CS);
		HistFile->Add(thisBinData_TkMu0_hist_HX);
	
*/

  	  char histName[200];
		  sprintf(histName,"thisBinMC_Double_hist_CS_rap%d_pT%d__costh_CS_phi_CS",yBin+1,ptBin+1);
		TH2F * thisBinMC_Double_hist_CS = (TH2F*)HistFile->Get(histName);
		sprintf(histName,"thisBinMC_Double_hist_HX_rap%d_pT%d__costh_HX_phi_HX",yBin+1,ptBin+1);
		TH2F * thisBinMC_Double_hist_HX = (TH2F*)HistFile->Get(histName);
		  sprintf(histName,"thisBinMC_TkMu0_hist_CS_rap%d_pT%d__costh_CS_phi_CS",yBin+1,ptBin+1);
		TH2F * thisBinMC_TkMu0_hist_CS = (TH2F*)HistFile->Get(histName);
		  sprintf(histName,"thisBinMC_TkMu0_hist_HX_rap%d_pT%d__costh_HX_phi_HX",yBin+1,ptBin+1);
		TH2F * thisBinMC_TkMu0_hist_HX = (TH2F*)HistFile->Get(histName);

		sprintf(histName,"thisBinData_Double_hist_CS_rap%d_pT%d__costh_CS_phi_CS",yBin+1,ptBin+1);
		TH2F * thisBinData_Double_hist_CS = (TH2F*)HistFile->Get(histName);
		  sprintf(histName,"thisBinData_Double_hist_HX_rap%d_pT%d__costh_HX_phi_HX",yBin+1,ptBin+1);
		TH2F * thisBinData_Double_hist_HX = (TH2F*)HistFile->Get(histName);
		  sprintf(histName,"thisBinData_TkMu0_hist_CS_rap%d_pT%d__costh_CS_phi_CS",yBin+1,ptBin+1);
		TH2F * thisBinData_TkMu0_hist_CS = (TH2F*)HistFile->Get(histName);
		  sprintf(histName,"thisBinData_TkMu0_hist_HX_rap%d_pT%d__costh_HX_phi_HX",yBin+1,ptBin+1);
		TH2F * thisBinData_TkMu0_hist_HX = (TH2F*)HistFile->Get(histName);

	      TH2F *DoubleRatioCS = (TH2F*) thisBinData_Double_hist_CS->Clone("DoubleRatioCS");
	      TH2F *DoubleRatioHX = (TH2F*) thisBinData_Double_hist_HX->Clone("DoubleRatioHX");
	      TH2F *TkMuRatioCS = (TH2F*) thisBinData_TkMu0_hist_CS->Clone("TkMuRatioCS");
	      TH2F *TkMuRatioHX = (TH2F*) thisBinData_TkMu0_hist_HX->Clone("TkMuRatioHX");


/*		thisBinMC_Double_hist_CS->Sumw2();
		thisBinMC_Double_hist_HX->Sumw2();
		thisBinMC_TkMu0_hist_CS->Sumw2();
		thisBinMC_TkMu0_hist_HX->Sumw2();

		thisBinData_Double_hist_CS->Sumw2();
		thisBinData_Double_hist_HX->Sumw2();
		thisBinData_TkMu0_hist_CS->Sumw2();
		thisBinData_TkMu0_hist_HX->Sumw2();
*/

//     	1) real data costh.vs.phi with DoubleMu0 trigger
//     	2) real data costh.vs.phi with TkMu0 trigger
//
//     	3) MC costh.vs.phi with DoubleMu0 trigger
//     	4) MC costh.vs.phi with TkMu0 trigger
//
//     	The double ratio is
//
//     	 1 / 3
//     	-------
//     	 2 / 4
//

	      if (eff){

		char orignameCS[200];
		sprintf(orignameCS,"hAcc2D_CS_pT%d_rap%d",ptBin+1,yBin+1);
		char orignameHX[200];
		sprintf(orignameHX,"hAcc2D_HX_pT%d_rap%d",ptBin+1,yBin+1);

		thisBinMC_Double_hist_CS = (TH2F*)EffHistFileDouble->Get(orignameCS);
		thisBinMC_Double_hist_HX = (TH2F*)EffHistFileDouble->Get(orignameHX);
		thisBinMC_TkMu0_hist_CS = (TH2F*)EffHistFileTrk->Get(orignameCS);
		thisBinMC_TkMu0_hist_HX = (TH2F*)EffHistFileTrk->Get(orignameHX);

	      for(int costhbin=1;costhbin<41+1;costhbin++){
	          for(int phibin=1;phibin<37+1;phibin++){
	        	  if(thisBinMC_Double_hist_CS->GetBinContent(costhbin,phibin)<0.05) thisBinMC_Double_hist_CS->SetBinContent(costhbin,phibin,0);
	        	  if(thisBinMC_Double_hist_HX->GetBinContent(costhbin,phibin)<0.05) thisBinMC_Double_hist_HX->SetBinContent(costhbin,phibin,0);
	        	  if(thisBinMC_TkMu0_hist_CS->GetBinContent(costhbin,phibin)<0.05) thisBinMC_TkMu0_hist_CS->SetBinContent(costhbin,phibin,0);
	        	  if(thisBinMC_TkMu0_hist_HX->GetBinContent(costhbin,phibin)<0.05) thisBinMC_TkMu0_hist_HX->SetBinContent(costhbin,phibin,0);
	      }}


	      double effDouble_CS[thisBinMC_Double_hist_CS->GetNbinsX()+1][thisBinMC_Double_hist_CS->GetNbinsY()+1];
	      double effDouble_HX[thisBinMC_Double_hist_CS->GetNbinsX()+1][thisBinMC_Double_hist_CS->GetNbinsY()+1];
	      double effTrk_CS[thisBinMC_Double_hist_CS->GetNbinsX()+1][thisBinMC_Double_hist_CS->GetNbinsY()+1];
	      double effTrk_HX[thisBinMC_Double_hist_CS->GetNbinsX()+1][thisBinMC_Double_hist_CS->GetNbinsY()+1];

	      for(int costhbin=1;costhbin<thisBinMC_TkMu0_hist_HX->GetNbinsX()+1;costhbin++){
	          for(int phibin=1;phibin<thisBinMC_TkMu0_hist_HX->GetNbinsY()+1;phibin++){
	        	  effDouble_CS[costhbin][phibin]=thisBinMC_Double_hist_CS->GetBinContent(costhbin,phibin);
	        	  effDouble_HX[costhbin][phibin]=thisBinMC_Double_hist_HX->GetBinContent(costhbin,phibin);
	        	  effTrk_CS[costhbin][phibin]=thisBinMC_Double_hist_CS->GetBinContent(costhbin,phibin);
	        	  effTrk_HX[costhbin][phibin]=thisBinMC_TkMu0_hist_HX->GetBinContent(costhbin,phibin);

	      }
	      }

	      int binfactor_costh;
	      int binfactor_phi;
	      int binscosth;
	      int binsphi;

	      cout<<"REBINNING"<<endl;
	      binfactor_costh=2;
	      binfactor_phi=2;
	      binscosth=40;
	      binsphi=36;


	      TH2F *effDouble_CS_ = new TH2F("effDouble_CS_","effDouble_CS_",binscosth,-1,1,binsphi,0,90);
	      TH2F *effDouble_HX_ = new TH2F("effDouble_HX_","effDouble_HX_",binscosth,-1,1,binsphi,0,90);
	      TH2F *effTrk_CS_ = new TH2F("effTrk_CS_","effTrk_CS_",binscosth,-1,1,binsphi,0,90);
	      TH2F *effTrk_HX_ = new TH2F("effTrk_HX_","effTrk_HX_",binscosth,-1,1,binsphi,0,90);

	      for(int costhbin=1;costhbin<thisBinMC_TkMu0_hist_HX->GetNbinsX()+1;costhbin++){
	          for(int phibin=1;phibin<thisBinMC_TkMu0_hist_HX->GetNbinsY()+1;phibin++){
	        	  for(int i=1;i<binfactor_costh+1;i++){
	            	  for(int j=1;j<binfactor_phi+1;j++){

	            		  effDouble_CS_->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,effDouble_CS[costhbin][phibin]);
	            		  effDouble_HX_->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,effDouble_HX[costhbin][phibin]);
	            		  effTrk_CS_->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,effTrk_CS[costhbin][phibin]);
	            		  effTrk_HX_->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,effTrk_HX[costhbin][phibin]);

	            	  }
	        	  }
	          }
	 	  }

		    DoubleRatioCS->Divide(effDouble_CS_);
		    DoubleRatioHX->Divide(effDouble_HX_);

		    TkMuRatioCS->Divide(effTrk_CS_);
		    TkMuRatioHX->Divide(effTrk_HX_);




	      }


	      else{
	    DoubleRatioCS->Divide(thisBinMC_Double_hist_CS);
	    DoubleRatioHX->Divide(thisBinMC_Double_hist_HX);

	    TkMuRatioCS->Divide(thisBinMC_TkMu0_hist_CS);
	    TkMuRatioHX->Divide(thisBinMC_TkMu0_hist_HX);
	      }

	    TH2F *CSratio = (TH2F*) DoubleRatioCS->Clone("CSratio");
	    TH2F *HXratio = (TH2F*) DoubleRatioHX->Clone("HXratio");



	    CSratio->SetNormFactor(1/CSratio->GetSumOfWeights());
	    HXratio->SetNormFactor(1/CSratio->GetSumOfWeights());
	    TkMuRatioCS->SetNormFactor(1/TkMuRatioCS->GetSumOfWeights());
	    TkMuRatioHX->SetNormFactor(1/TkMuRatioHX->GetSumOfWeights());


		CSratio->Divide(TkMuRatioCS);
		HXratio->Divide(TkMuRatioHX);



		char CSRATname[200];
		sprintf(CSRATname,"CSratio_rap%d_pT%d",yBin+1,ptBin+1);
	    TH2F *CSRAT = (TH2F*) CSratio->Clone(CSRATname);
		char HXRATname[200];
		sprintf(HXRATname,"HXratio_rap%d_pT%d",yBin+1,ptBin+1);
	    TH2F *HXRAT = (TH2F*) HXratio->Clone(HXRATname);


		RatioFile->Add(CSRAT);
		RatioFile->Add(HXRAT);



//		CSratio->Print("all");

		double errorFromZero=10;

	      for(int costhbin=1;costhbin<41+1;costhbin++){
	          for(int phibin=1;phibin<37+1;phibin++){
//	        	  if(CSratio->GetBinContent(costhbin,phibin)==0) CSratio->SetBinError(costhbin,phibin,errorFromZero);
//	        	  else CSratio->SetBinError(costhbin,phibin,0);
//	        	  if(HXratio->GetBinContent(costhbin,phibin)==0) HXratio->SetBinError(costhbin,phibin,errorFromZero);
//	        	  else HXratio->SetBinError(costhbin,phibin,0);
	      }}

//	    CSratio->Print("all");

		gStyle->SetPalette(1);
		char FilenameAcc[200];


		TCanvas* DoubleCanvasCS = new TCanvas("DoubleCanvasCS","DoubleCanvasCS",3200,3200);
		DoubleCanvasCS->Divide(3,3);  DoubleCanvasCS->SetFillColor(kWhite);

		DoubleCanvasCS->cd(1) ; gPad->SetFillColor(kWhite); thisBinData_Double_hist_CS->Draw("colz");
		DoubleCanvasCS->cd(2) ; gPad->SetFillColor(kWhite); thisBinMC_Double_hist_CS->Draw("colz");
		DoubleCanvasCS->cd(3) ; gPad->SetFillColor(kWhite); DoubleRatioCS->Draw("colz");

		DoubleCanvasCS->cd(4) ; gPad->SetFillColor(kWhite); thisBinData_TkMu0_hist_CS->Draw("colz");
		DoubleCanvasCS->cd(5) ; gPad->SetFillColor(kWhite); thisBinMC_TkMu0_hist_CS->Draw("colz");
		DoubleCanvasCS->cd(6) ; gPad->SetFillColor(kWhite); TkMuRatioCS->Draw("colz");

		DoubleCanvasCS->cd(7) ; gPad->SetFillColor(kWhite); DoubleRatioCS->Draw("colz");
		DoubleCanvasCS->cd(8) ; gPad->SetFillColor(kWhite); TkMuRatioCS->Draw("colz");
		DoubleCanvasCS->cd(9) ; gPad->SetFillColor(kWhite); CSratio->Draw("colz");

	    if(eff) sprintf(FilenameAcc,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/SysChecks/SysCheckHistos_eff_CS_rap%d_pt%d.png",yBin+1,ptBin+1);
	    else sprintf(FilenameAcc,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/SysChecks/SysCheckHistos_CS_rap%d_pt%d.png",yBin+1,ptBin+1);

	    DoubleCanvasCS->SaveAs(FilenameAcc);
	    DoubleCanvasCS->Close();


		TCanvas* DoubleCanvasHX = new TCanvas("DoubleCanvasHX","DoubleCanvasHX",3200,3200);
		DoubleCanvasHX->Divide(3,3);  DoubleCanvasHX->SetFillColor(kWhite);

		DoubleCanvasHX->cd(1) ; gPad->SetFillColor(kWhite); thisBinData_Double_hist_HX->Draw("colz");
		DoubleCanvasHX->cd(2) ; gPad->SetFillColor(kWhite); thisBinMC_Double_hist_HX->Draw("colz");
		DoubleCanvasHX->cd(3) ; gPad->SetFillColor(kWhite); DoubleRatioHX->Draw("colz");

		DoubleCanvasHX->cd(4) ; gPad->SetFillColor(kWhite); thisBinData_TkMu0_hist_HX->Draw("colz");
		DoubleCanvasHX->cd(5) ; gPad->SetFillColor(kWhite); thisBinMC_TkMu0_hist_HX->Draw("colz");
		DoubleCanvasHX->cd(6) ; gPad->SetFillColor(kWhite); TkMuRatioHX->Draw("colz");

		DoubleCanvasHX->cd(7) ; gPad->SetFillColor(kWhite); DoubleRatioHX->Draw("colz");
		DoubleCanvasHX->cd(8) ; gPad->SetFillColor(kWhite); TkMuRatioHX->Draw("colz");
		DoubleCanvasHX->cd(9) ; gPad->SetFillColor(kWhite); HXratio->Draw("colz");

	    if(eff) sprintf(FilenameAcc,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/SysChecks/SysCheckHistos_eff_HX_rap%d_pt%d.png",yBin+1,ptBin+1);
	    else sprintf(FilenameAcc,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/SysChecks/SysCheckHistos_HX_rap%d_pt%d.png",yBin+1,ptBin+1);

	    DoubleCanvasHX->SaveAs(FilenameAcc);
	    DoubleCanvasHX->Close();


	    RooDataHist* CSratioHist = new RooDataHist("CSratioHist","CSratioHist",RooArgList(costh_CS,phi_CS),CSratio);
	    RooDataHist* HXratioHist = new RooDataHist("HXratioHist","HXratioHist",RooArgList(costh_HX,phi_HX),HXratio);

///////////// FIT ////////////

	    double bound=100;

	    RooRealVar lambdathetaCS("lambdathetaCS","lambdathetaCS",0,-bound,bound);
	    RooRealVar lambdaphiCS("lambdaphiCS","lambdaphiCS",0,-bound,bound);
	    RooRealVar lambdathetaphiCS("lambdathetaphiCS","lambdathetaphiCS",0,-bound,bound);

	    RooPolarizationPdf* polpdfCS = new RooPolarizationPdf("polpdfCS","polpdfCS",costh_CS,phi_CS,lambdathetaCS,lambdaphiCS,lambdathetaphiCS);

	    RooRealVar lambdathetaHX("lambdathetaHX","lambdathetaHX",0,-bound,bound);
	    RooRealVar lambdaphiHX("lambdaphiHX","lambdaphiHX",0,-bound,bound);
	    RooRealVar lambdathetaphiHX("lambdathetaphiHX","lambdathetaphiHX",0,-bound,bound);

	    RooPolarizationPdf* polpdfHX = new RooPolarizationPdf("polpdfHX","polpdfHX",costh_HX,phi_HX,lambdathetaHX,lambdaphiHX,lambdathetaphiHX);


//		RooFitResult* fitResCS=polpdfCS->chi2FitTo(*CSratioHist,DataError(0));

//		RooFitResult* fitResHX=polpdfHX->chi2FitTo(*HXratioHist,DataError(0));


  double lamthCS [2][8]={
		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.486461,-0.0352042},
		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.520544,-0.0808785,-0.0319447}
  };

  double lamphCS [2][8]={
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,0.367306,-0.166666},
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,0.402154,0.112628,-0.0383048}
   };

  double lamthphCS [2][8]={
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.00370167,0.0407368},
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.10504,-0.0405559,-0.0052447}
   };


  double lamthHX [2][8]={
		  {0.000000,0.000000,0.000000,0.000000,0.000000,0.852886,0.132954},
		  {0.000000,0.000000,0.000000,0.000000,0.000000,1.69006,0.0520585,0.0212357}
  };

  double lamphHX [2][8]={
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.0763282, -0.0762038},
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.118597,-0.129362,-0.0416266}
   };

  double lamthphHX [2][8]={
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.0383324,0.0206791},
 		  {0.000000,0.000000,0.000000,0.000000,0.000000,-0.179886,-0.0433619,-0.198789}
   };



  lambdathetaCS.setVal(lamthCS[yBin][ptBin]);
  lambdaphiCS.setVal(lamphCS[yBin][ptBin]);
  lambdathetaphiCS.setVal(lamthphCS[yBin][ptBin]);

  lambdathetaHX.setVal(lamthHX[yBin][ptBin]);
  lambdaphiHX.setVal(lamphHX[yBin][ptBin]);
  lambdathetaphiHX.setVal(lamthphHX[yBin][ptBin]);

		double norm=1;

	    RooPlot* Frame_costh_CS = new RooPlot;
	    Frame_costh_CS = costh_CS.frame() ;
	    CSratioHist->plotOn(Frame_costh_CS,DataError(RooAbsData::SumW2),MarkerSize(0.4));
//	    polpdfCS->plotOn(Frame_costh_CS,LineWidth(2),Normalization(norm));
//	    polpdfCS->paramOn(Frame_costh_CS,Format("NE",AutoPrecision(2)),Layout(0.10,0.9,1));
	    Frame_costh_CS->SetTitle(0);


	    RooPlot* Frame_phi_CS = new RooPlot;
	    Frame_phi_CS = phi_CS.frame() ;
	    CSratioHist->plotOn(Frame_phi_CS,DataError(RooAbsData::SumW2),MarkerSize(0.4));
//	    polpdfCS->plotOn(Frame_phi_CS,LineWidth(2),Normalization(norm));
//	    polpdfCS->paramOn(Frame_phi_CS,Format("NE",AutoPrecision(2)),Layout(0.10,0.9,1));
	    Frame_phi_CS->SetTitle(0);

	    RooPlot* Frame_costh_HX = new RooPlot;
	    Frame_costh_HX = costh_HX.frame() ;
	    HXratioHist->plotOn(Frame_costh_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
//	    polpdfHX->plotOn(Frame_costh_HX,LineWidth(2),Normalization(norm));
//	    polpdfHX->paramOn(Frame_costh_HX,Format("NE",AutoPrecision(2)),Layout(0.10,0.9,1));
	    Frame_costh_HX->SetTitle(0);

	    RooPlot* Frame_phi_HX = new RooPlot;
	    Frame_phi_HX = phi_HX.frame() ;
	    HXratioHist->plotOn(Frame_phi_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
//	    polpdfHX->plotOn(Frame_phi_HX,LineWidth(2),Normalization(norm));
//	    polpdfHX->paramOn(Frame_phi_HX,Format("NE",AutoPrecision(2)),Layout(0.10,0.9,1));
	    Frame_phi_HX->SetTitle(0);

		TCanvas* FitCanvas = new TCanvas("FitCanvas","FitCanvas",3200,2600);
		FitCanvas->Divide(2,2);  FitCanvas->SetFillColor(kWhite);
		FitCanvas->cd(1) ; gPad->SetFillColor(kWhite); gPad->SetTopMargin(0.15);Frame_costh_CS->Draw();
		FitCanvas->cd(2) ; gPad->SetFillColor(kWhite); gPad->SetTopMargin(0.15);Frame_phi_CS->Draw();
		FitCanvas->cd(3) ; gPad->SetFillColor(kWhite); gPad->SetTopMargin(0.15);Frame_costh_HX->Draw();
		FitCanvas->cd(4) ; gPad->SetFillColor(kWhite); gPad->SetTopMargin(0.15);Frame_phi_HX->Draw();

	    if(eff) sprintf(FilenameAcc,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/SysChecks/SysCheckHistos_eff_1Dnomodel_rap%d_pt%d.png",yBin+1,ptBin+1);
	    else sprintf(FilenameAcc,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/SysChecks/SysCheckHistos_1Dnomodel_rap%d_pt%d.png",yBin+1,ptBin+1);

	    FitCanvas->SaveAs(FilenameAcc);
	    FitCanvas->Close();




	delete CSratioHist;
	delete HXratioHist;

      }
    }



//	HistFile->Write();

	RatioFile->Write();
	RatioFile->Close();
	delete RatioFile;

	EffHistFileDouble->Close();
	delete EffHistFileDouble;
	EffHistFileTrk->Close();
	delete EffHistFileTrk;

	HistFile->Close();
	delete HistFile;

//  delete fInMC;
//  delete dataMC;
//  delete fInData;
//  delete dataData;


  return 0;
}
