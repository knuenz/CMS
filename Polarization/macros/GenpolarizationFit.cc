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

  bool grid(false);
  if(grid) cout<<"Hello Grid :)"<<endl;

  bool doprompt(true), dononprompt(false), dobkg(false),dopol(true)/*,pereverr(true),newdata(false)*/,fitHX(true),fitCS(true);

  gStyle->SetPalette(1);
  gStyle->SetTitleFillColor(kWhite);

  bool ploterrorband(false);

  bool initialhesse(false);

  bool migrad(true);
  bool improve(false);
  bool simplex(false);
  bool minos(false);

  bool hesse(true);

  bool scan(false);
  bool fluc(false);

//  bool fast(true);

  bool useAccMaps(false);

  bool rebin(true);

  bool rebuildmodel(false);

  bool Prompt(false);

  bool JobIDinput(false);

  char minimizerType [200];
  char minimizerAlgo [200];
  sprintf(minimizerType,"Minuit2");
  sprintf(minimizerAlgo,"Scan");

  bool minimizeBFGS(false);

/// define generated statistics
  bool genLevel(false);
  bool accLevel(false);
  bool effLevel(false);
  bool fixNum(false); int fixEvents = 100000;

  bool RealDataPRNPstats(false);
  bool DoubleMu0stage(true);
  bool moreBins(true);

  int n_errsigma = 3;

///////////////////////////////
  double pTin;
  double yin;
  char* JobIDinput_;

  for( int i=0;i < argc; ++i ) {
    if(std::string(argv[i]).find("--noHX") != std::string::npos) fitHX = false;
    if(std::string(argv[i]).find("--noCS") != std::string::npos) fitCS = false;
    if(std::string(argv[i]).find("--ploterrorband") != std::string::npos) ploterrorband = true;
    if(std::string(argv[i]).find("pT") != std::string::npos) {char* pTchar = argv[i]; char* pTchar2 = strtok (pTchar, "p"); pTin = atof(pTchar2); cout<<"pTin = "<<pTin<<endl;}
    if(std::string(argv[i]).find("rap") != std::string::npos) {char* ychar = argv[i]; char* ychar2 = strtok (ychar, "r"); yin = atof(ychar2); cout<<"yin = "<<yin<<endl;}
    if(std::string(argv[i]).find("--grid") != std::string::npos) grid = false;
    if(std::string(argv[i]).find("serial") != std::string::npos) {JobIDinput = true; JobIDinput_ = argv[i];}

  }

  if(grid) moreBins=false;

  if(ploterrorband) {fitHX = false; fitCS = false;}

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.4,2.4);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);



//  TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/jPsiFit_PR_FSR_ga0701re0702te0702_ZBCsmsmsm.root","UPDATE");
//  TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/jPsiFit_Bins_PR_FSR_geomAcc07Jan.root","UPDATE");
//  TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/jPsiFit_geomAcc07Jan_ZBC_sm.root","UPDATE");
//    TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/jPsiFit_PR_FSR_GA07_ZBC.root","UPDATE");
//  TFile *fitInput = new TFile("jPsiFit_PR_FSR_ga0701re0702te0702_ZBCsmsmsm.root","UPDATE");
//  TFile *fitInputmodel = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/jPsiFit_GA08RETE_smsmsm_pythia_less.root","UPDATE");
//  TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/jPsiFit_GA08RETE_smsmsm_pythia_less.root","UPDATE");
//  TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/jPsiFit_GA08RETE_smsmsm_cascade_less.root","UPDATE");
//    TFile *fitInput = new TFile("jPsiFit_PR_FSR_ga0701re0702te0702_ZBCsmsmsm.root","UPDATE");
//    TFile *fitInputmodel = new TFile("jPsiFit_PR_FSR_ga0701re0702te0702_ZBCsmsmsm.root","UPDATE");


//  fitInput->Print();


//  RooDataSet *data = NULL;
//  CompositeModelBuilder *hx, *cs;


  int iterations = 200;
  int generate_gt=0;

//  char output_cs_[200];
//	  sprintf(output_cs_,"jPsiFitFinal_cs_ToyMC.root");
//	  char output_hx_[200];
//	  sprintf(output_hx_,"jPsiFitFinal_hx_ToyMC.root");
//	  TFile *output_cs = new TFile(output_cs_,"RECREATE");
//	  TFile *output_hx = new TFile(output_hx_,"RECREATE");

//	  TDirectory *current_in,*current_in_model,*current_out_hx,*current_out_cs;

		  TFile *resfile = new TFile("RooFitResult.root","RECREATE");

		  TFile *DataGenFile = new TFile("Model.root","UPDATE");

		  TFile *ErrBandPlots = new TFile("ErrBandPlots_GEN_minus05.root","RECREATE");



//		  TFile *accFile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
//		  TFile *recoEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
//		  TFile *trigEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");

//		  TFile *accFile = new TFile("/scratch/knuenz/Polarization/RootInput/geomAccHistos_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded.root","UPDATE");
//		  TFile *recoEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/recoEffHistos_ATLASPT_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
//		  TFile *trigEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");

		  TFile *accFile = new TFile("/scratch/knuenz/Polarization/RootInput/geomAccHistos_NP_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root","UPDATE");
		  TFile *recoEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/recoEffHistos_NP_ATLASPT_19March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
		  TFile *trigEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root","UPDATE");

//		  TFile *accFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/geomAccsmoothSingle.root","UPDATE");
//		  TFile *recoEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/recoEffsmoothSingle.root","UPDATE");
//		  TFile *trigEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/trigEffsmoothSingle.root","UPDATE");

		  TFile *accFileExtract = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
//		  TFile *accFileExtract = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/geomAccCascadesmooth11032011.root","UPDATE");

		  char resultname[200];


  for(int yBin =yin-1; yBin < jpsi::kNbRapForPTBins+yin-5; ++yBin) {
 	  for(int ptBin =pTin-1; ptBin < jpsi::kNbPTBins[yBin+1]+pTin-6-yin; ++ptBin) {

//  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-4; ++yBin) {
// 	  for(int ptBin = 5; ptBin < jpsi::kNbPTBins[yBin+1]-1; ++ptBin) {

//		  if(ptBin==4 && yBin==0)continue;
// 		  if(ptBin==7 && yBin==1)continue;

// 		  if(yBin==0)continue;
// 		  if(ptBin>4)continue;

 		  char inputfilenamepseudo[200];
 		  sprintf(inputfilenamepseudo,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/saveTrees/TTree_GEN200_UNSMOOTHnlb200_scen3_rap%d_pT%d.root",yBin+1,ptBin+1);
		  TFile *PseudoInput = new TFile(inputfilenamepseudo,"UPDATE");


 			char JobID[200];

 			sprintf(JobID,"multi_Fast_plus05_1_2_DoubleMu0_effstat");
//  			sprintf(JobID,"multi_try");

			char dirStruct[200];
 			sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/CRAB/%s",JobID);
 			if(JobIDinput) {sprintf(JobID,"%s",JobIDinput_); cout<<JobID<<endl; sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/SERIAL/%s",JobID); gSystem->mkdir(dirStruct); cout<<dirStruct<<endl; sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/SERIAL/%s_res",JobID);}
 			gSystem->mkdir(dirStruct);

 			char inputfilename[200];
 			sprintf(inputfilename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/src/JPsiPolarizationSave3/CRAB/%s/RooFitResult_rap%d_pt%d.root",JobID,yBin+1,ptBin+1);
 			if(JobIDinput) sprintf(inputfilename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/SERIAL/%s_res/RooFitResult_rap%d_pt%d.root",JobID,yBin+1,ptBin+1);

 			TFile *resfile_bins = new TFile(inputfilename,"UPDATE");


		      TH2F *AccMap_CS,*AccMap_HX;
		      TH2F *AccMapExtract_CS,*AccMapExtract_HX;
		      TH2F *RecoEffMap_CS,*RecoEffMap_HX;
		      TH2F *TrigEffMap_CS,*TrigEffMap_HX;

		      std::stringstream nameHXp,nameCSp;
		      nameHXp << "hAcc2D_HX_pT" << ptBin+1 << "_rap" << yBin+1;
		      nameCSp << "hAcc2D_CS_pT" << ptBin+1 << "_rap" << yBin+1;
		      AccMap_CS = (TH2F*)accFile->Get(nameCSp.str().c_str());
		      AccMap_HX = (TH2F*)accFile->Get(nameHXp.str().c_str());
		      AccMapExtract_CS = (TH2F*)accFileExtract->Get(nameCSp.str().c_str());
		      AccMapExtract_HX = (TH2F*)accFileExtract->Get(nameHXp.str().c_str());
		      RecoEffMap_CS = (TH2F*)recoEffFile->Get(nameCSp.str().c_str());
		      RecoEffMap_HX = (TH2F*)recoEffFile->Get(nameHXp.str().c_str());
		      TrigEffMap_CS = (TH2F*)trigEffFile->Get(nameCSp.str().c_str());
		      TrigEffMap_HX = (TH2F*)trigEffFile->Get(nameHXp.str().c_str());

		      TH2F *Corr_CS = (TH2F*) AccMap_CS->Clone("Corr_CS");
		      TH2F *Corr_HX = (TH2F*) AccMap_HX->Clone("Corr_HX");
		      TH2F *CorrExtract_CS = (TH2F*) AccMapExtract_CS->Clone("CorrExtract_CS");
		      TH2F *CorrExtract_HX = (TH2F*) AccMapExtract_HX->Clone("CorrExtract_HX");

	if(!rebin){
		      RecoEffMap_CS->Multiply(TrigEffMap_CS);
		      Corr_CS->Multiply(RecoEffMap_CS);
		      RecoEffMap_HX->Multiply(TrigEffMap_HX);
		      Corr_HX->Multiply(RecoEffMap_HX);

		      CorrExtract_CS->Multiply(RecoEffMap_CS);
		      CorrExtract_HX->Multiply(RecoEffMap_HX);

	}



	if(rebin){
		      double recoEffPrompt_CS[RecoEffMap_CS->GetNbinsX()+1][RecoEffMap_CS->GetNbinsY()+1];
		      double recoEffPrompt_HX[RecoEffMap_HX->GetNbinsX()+1][RecoEffMap_HX->GetNbinsY()+1];
		      double trigEffPrompt_CS[TrigEffMap_CS->GetNbinsX()+1][TrigEffMap_CS->GetNbinsY()+1];
		      double trigEffPrompt_HX[TrigEffMap_HX->GetNbinsX()+1][TrigEffMap_HX->GetNbinsY()+1];

		       for(int costhbin=1;costhbin<RecoEffMap_CS->GetNbinsX()+1;costhbin++){
		           for(int phibin=1;phibin<RecoEffMap_CS->GetNbinsY()+1;phibin++){
		        	   recoEffPrompt_CS[costhbin][phibin]=RecoEffMap_CS->GetBinContent(costhbin,phibin);
		        	   recoEffPrompt_HX[costhbin][phibin]=RecoEffMap_HX->GetBinContent(costhbin,phibin);
		        	   trigEffPrompt_CS[costhbin][phibin]=TrigEffMap_CS->GetBinContent(costhbin,phibin);
		        	   trigEffPrompt_HX[costhbin][phibin]=TrigEffMap_HX->GetBinContent(costhbin,phibin);
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

		       TH2F *promptRecoEffRebin_CS = new TH2F("promptRecoEffRebin_CS","promptRecoEffRebin_CS",binscosth,-1,1,binsphi,0,90);
		       TH2F *promptRecoEffRebin_HX = new TH2F("promptRecoEffRebin_HX","promptRecoEffRebin_HX",binscosth,-1,1,binsphi,0,90);
		       TH2F *promptTrigEffRebin_CS = new TH2F("promptTrigEffRebin_CS","promptTrigEffRebin_CS",binscosth,-1,1,binsphi,0,90);
		       TH2F *promptTrigEffRebin_HX = new TH2F("promptTrigEffRebin_HX","promptTrigEffRebin_HX",binscosth,-1,1,binsphi,0,90);

		       for(int costhbin=1;costhbin<RecoEffMap_CS->GetNbinsX()+1;costhbin++){
		           for(int phibin=1;phibin<RecoEffMap_CS->GetNbinsY()+1;phibin++){
		         	  for(int i=1;i<binfactor_costh+1;i++){
		             	  for(int j=1;j<binfactor_phi+1;j++){

		             		promptRecoEffRebin_CS->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,recoEffPrompt_CS[costhbin][phibin]);
		             		promptRecoEffRebin_HX->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,recoEffPrompt_HX[costhbin][phibin]);
		             		promptTrigEffRebin_CS->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,trigEffPrompt_CS[costhbin][phibin]);
		             		promptTrigEffRebin_HX->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,trigEffPrompt_HX[costhbin][phibin]);

		             	  }
		         	  }
		           }
		  	  }


		      promptRecoEffRebin_CS->Multiply(promptTrigEffRebin_CS);
		      Corr_CS->Multiply(promptRecoEffRebin_CS);
		      promptRecoEffRebin_HX->Multiply(promptTrigEffRebin_HX);
		      Corr_HX->Multiply(promptRecoEffRebin_HX);
		      Corr_CS->SaveAs("Corr_CS.root");
                      
	}


		      char Corr_CS_hist_char[200];
		      char Corr_CS_histfunc_char[200];
		      char map_uniform_CS_char[200];
		      char map_theta_CS_char[200];
		      char map_phi_CS_char[200];
		      char map_thetaphi_CS_char[200];
		      char CorrExtract_CS_hist_char[200];
		      char CorrExtract_CS_histfunc_char[200];
		      char map_uniformExtract_CS_char[200];
		      char map_thetaExtract_CS_char[200];
		      char map_phiExtract_CS_char[200];
		      char map_thetaphiExtract_CS_char[200];
		      char Corr_HX_hist_char[200];
		      char Corr_HX_histfunc_char[200];
		      char map_uniform_HX_char[200];
		      char map_theta_HX_char[200];
		      char map_phi_HX_char[200];
		      char map_thetaphi_HX_char[200];
		      char CorrExtract_HX_hist_char[200];
		      char CorrExtract_HX_histfunc_char[200];
		      char map_uniformExtract_HX_char[200];
		      char map_thetaExtract_HX_char[200];
		      char map_phiExtract_HX_char[200];
		      char map_thetaphiExtract_HX_char[200];

		      sprintf(Corr_CS_hist_char,"Corr_CS_hist_%d_%d",yBin+1,ptBin+1);
		      sprintf(Corr_CS_histfunc_char,"Corr_CS_histfunc_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_uniform_CS_char,"map_uniform_CS_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_theta_CS_char,"map_theta_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_phi_CS_char,"map_phi_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaphi_CS_char,"map_thetaphi_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(CorrExtract_CS_hist_char,"CorrExtract_CS_hist_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(CorrExtract_CS_histfunc_char,"CorrExtract_CS_histfunc_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_uniformExtract_CS_char,"map_uniformExtract_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaExtract_CS_char,"map_thetaExtract_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_phiExtract_CS_char,"map_phiExtract_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaphiExtract_CS_char,"map_thetaphiExtract_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(Corr_HX_hist_char,"Corr_HX_hist_%d_%d",yBin+1,ptBin+1);
		      sprintf(Corr_HX_histfunc_char,"Corr_HX_histfunc_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_uniform_HX_char,"map_uniform_HX_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_theta_HX_char,"map_theta_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_phi_HX_char,"map_phi_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaphi_HX_char,"map_thetaphi_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(CorrExtract_HX_hist_char,"CorrExtract_HX_hist_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(CorrExtract_HX_histfunc_char,"CorrExtract_HX_histfunc_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_uniformExtract_HX_char,"map_uniformExtract_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaExtract_HX_char,"map_thetaExtract_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_phiExtract_HX_char,"map_phiExtract_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaphiExtract_HX_char,"map_thetaphiExtract_HX_char%d_pt%d",yBin+1,ptBin+1);

		      int linsmooth=0;

		  	RooDataHist Corr_CS_hist(Corr_CS_hist_char,Corr_CS_hist_char,RooArgSet(costh_CS,phi_CS),Corr_CS);
		  	RooDataHist Corr_HX_hist(Corr_HX_hist_char,Corr_HX_hist_char,RooArgSet(costh_HX,phi_HX),Corr_HX);
		  	RooHistFunc Corr_CS_histfunc = RooHistFunc(Corr_CS_histfunc_char,Corr_CS_histfunc_char, RooArgSet(costh_CS,phi_CS),Corr_CS_hist,linsmooth);
		  	RooHistFunc Corr_HX_histfunc = RooHistFunc(Corr_HX_histfunc_char,Corr_HX_histfunc_char, RooArgSet(costh_HX,phi_HX),Corr_HX_hist,linsmooth);
		  	RooHistFunc map_uniform_CS = RooHistFunc(map_uniform_CS_char,map_uniform_CS_char,RooArgSet(costh_CS,phi_CS),Corr_CS_hist,linsmooth);
		  	RooFormulaVar map_theta_CS = RooFormulaVar(map_theta_CS_char,map_theta_CS_char,"@2(@0,@1)*@0*@0",RooArgList(costh_CS, phi_CS, Corr_CS_histfunc));
		  	RooFormulaVar map_phi_CS = RooFormulaVar(map_phi_CS_char,map_phi_CS_char,"@2(@0,@1)*(1-@0*@0)*cos(@1*2*pi/180)",RooArgList(costh_CS, phi_CS, Corr_CS_histfunc));
			RooFormulaVar map_thetaphi_CS = RooFormulaVar(map_thetaphi_CS_char,map_thetaphi_CS_char,"@2(@0,@1)*2*sqrt(1-@0*@0)*@0*cos(@1*pi/180)",RooArgList(costh_CS, phi_CS, Corr_CS_histfunc));
		  	RooHistFunc map_uniform_HX = RooHistFunc(map_uniform_HX_char,map_uniform_HX_char,RooArgSet(costh_HX,phi_HX),Corr_HX_hist,linsmooth);
		  	RooFormulaVar map_theta_HX = RooFormulaVar(map_theta_HX_char,map_theta_HX_char,"@2(@0,@1)*@0*@0",RooArgList(costh_HX, phi_HX, Corr_HX_histfunc));
		  	RooFormulaVar map_phi_HX = RooFormulaVar(map_phi_HX_char,map_phi_HX_char,"@2(@0,@1)*(1-@0*@0)*cos(@1*2*pi/180)",RooArgList(costh_HX, phi_HX, Corr_HX_histfunc));
			RooFormulaVar map_thetaphi_HX = RooFormulaVar(map_thetaphi_HX_char,map_thetaphi_HX_char,"@2(@0,@1)*2*sqrt(1-@0*@0)*@0*cos(@1*pi/180)",RooArgList(costh_HX, phi_HX, Corr_HX_histfunc));

		  	RooDataHist CorrExtract_CS_hist(CorrExtract_CS_hist_char,CorrExtract_CS_hist_char,RooArgSet(costh_CS,phi_CS),CorrExtract_CS);
		  	RooDataHist CorrExtract_HX_hist(CorrExtract_HX_hist_char,CorrExtract_HX_hist_char,RooArgSet(costh_HX,phi_HX),CorrExtract_HX);
		  	RooHistFunc CorrExtract_CS_histfunc = RooHistFunc(CorrExtract_CS_histfunc_char,CorrExtract_CS_histfunc_char, RooArgSet(costh_CS,phi_CS),CorrExtract_CS_hist,linsmooth);
		  	RooHistFunc CorrExtract_HX_histfunc = RooHistFunc(CorrExtract_CS_histfunc_char,CorrExtract_CS_histfunc_char, RooArgSet(costh_HX,phi_HX),CorrExtract_HX_hist,linsmooth);
		  	RooHistFunc map_uniformExtract_CS = RooHistFunc(map_uniformExtract_CS_char,map_uniformExtract_CS_char,RooArgSet(costh_CS,phi_CS),CorrExtract_CS_hist,linsmooth);
		  	RooFormulaVar map_thetaExtract_CS = RooFormulaVar(map_thetaExtract_CS_char,map_thetaExtract_CS_char,"@2(@0,@1)*@0*@0",RooArgList(costh_CS, phi_CS, CorrExtract_CS_histfunc));
		  	RooFormulaVar map_phiExtract_CS = RooFormulaVar(map_phiExtract_CS_char,map_phiExtract_CS_char,"@2(@0,@1)*(1-@0*@0)*cos(@1*2*pi/180)",RooArgList(costh_CS, phi_CS, CorrExtract_CS_histfunc));
			RooFormulaVar map_thetaphiExtract_CS = RooFormulaVar(map_thetaphiExtract_CS_char,map_thetaphiExtract_CS_char,"@2(@0,@1)*2*sqrt(1-@0*@0)*@0*cos(@1*pi/180)",RooArgList(costh_CS, phi_CS, CorrExtract_CS_histfunc));
		  	RooHistFunc map_uniformExtract_HX = RooHistFunc(map_uniformExtract_HX_char,map_uniformExtract_HX_char,RooArgSet(costh_HX,phi_HX),CorrExtract_HX_hist,linsmooth);
			RooFormulaVar map_thetaExtract_HX = RooFormulaVar(map_thetaExtract_HX_char,map_thetaExtract_HX_char,"@2(@0,@1)*@0*@0",RooArgList(costh_HX, phi_HX, CorrExtract_HX_histfunc));
			RooFormulaVar map_phiExtract_HX = RooFormulaVar(map_phiExtract_HX_char,map_phiExtract_HX_char,"@2(@0,@1)*(1-@0*@0)*cos(@1*2*pi/180)",RooArgList(costh_HX, phi_HX, CorrExtract_HX_histfunc));
			RooFormulaVar map_thetaphiExtract_HX = RooFormulaVar(map_thetaphiExtract_HX_char,map_thetaphiExtract_HX_char,"@2(@0,@1)*2*sqrt(1-@0*@0)*@0*cos(@1*pi/180)",RooArgList(costh_HX, phi_HX, CorrExtract_HX_histfunc));












 			  for (int iteration = 1; iteration < iterations+1; iteration++){
 				  cout<<"GENERATION "<<iteration<<endl;

      std::stringstream binName,cutString, binNameOut;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      binNameOut << "pt" << ptBin+1 << "_rapidity" << yBin+1 << "_gen" << iteration;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];
      std::cout << cutString.str() << std::endl;


/*
      JpsiMass.setRange("lowBand",2.7,
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow);
      JpsiMass.setRange("signalRegionCS",
			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass,
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass);
      JpsiMass.setRange("signalRegionHX",
      			jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass,
      			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigMass);
      JpsiMass.setRange("highBand",
			jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh,3.5);

*/
      bool writeProgress(true);
      if(writeProgress){
  	char outputfilename[200];
  	sprintf(outputfilename,"FitProgressToyMC_%s.txt",JobID);
  	if(JobIDinput) sprintf(outputfilename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/SERIAL/FitProgressToyMC_%s.txt",JobIDinput_);
//  	if(JobIDinput) sprintf(outputfilename,"/afs/hephy.at/user/k/knuenz/www/Polarization/FitProgress/FitProgressToyMC_%s.txt",JobIDinput_);
  	FILE *outputFile = fopen(outputfilename,"w");
	fprintf(outputFile, "rapidity%d_pt%d, iteration %d\n",yBin+1,ptBin+1,iteration);
	fclose(outputFile);
      }



//      if(current_in = fitInput->GetDirectory(binName.str().c_str())) {

    int numEvents_;
    if(effLevel) numEvents_=jpsi::numEventsBin_HLT_Mu0_TkMu0_OST_Jpsi[yBin][ptBin];
    if(accLevel || genLevel) numEvents_=jpsi::numEventsBin[yBin][ptBin];
    if(fixNum) numEvents_=fixEvents;
    if(RealDataPRNPstats){
    	if(Prompt)numEvents_=jpsi::numEventsPR[yBin][ptBin];
    	if(!Prompt)numEvents_=jpsi::numEventsNP[yBin][ptBin];
    }
    if(DoubleMu0stage){
    	if(Prompt)numEvents_=jpsi::numEventsBin_DoubleMu0[yBin][ptBin]*(1-jpsi::Bfrac[yBin][ptBin]);
    	
if(!Prompt)numEvents_=jpsi::numEventsBin_DoubleMu0[yBin][ptBin]*(1-jpsi::fracPRregion[yBin][ptBin])*jpsi::fNPinNP_actual[yBin][ptBin];//jpsi::numEventsBin_DoubleMu0[yBin][ptBin]*jpsi::Bfrac[yBin][ptBin];

    }
    cout<<"number of events to be generated: "<<numEvents_<<endl;

/*	hx = new CompositeModelBuilder("HX","");
	cs = new CompositeModelBuilder("CS","");

	hx->setUsePrompt(doprompt);
	hx->setUseNonPrompt(dononprompt);
	hx->setUseBkg(dobkg);
	hx->setUsePol(dopol);
	hx->setUseAcceptanceMaps(useAccMaps);
	hx->setUseMass(false);
	hx->setUseLifetime(false);


	cs->setUsePrompt(doprompt);
	cs->setUseNonPrompt(dononprompt);
	cs->setUseBkg(dobkg);
	cs->setUsePol(dopol);
	cs->setUseAcceptanceMaps(useAccMaps);
	cs->setUseMass(false);
	cs->setUseLifetime(false);
*/

//	cs->loadParameters(*current_in);
//	hx->loadParameters(*current_in);

//	cs->getMassModel()->fix("CBm");
//	cs->getMassModel()->fix("CBs");


//	hx->getMassModel()->fix("CBm");
//	hx->getMassModel()->fix("CBs");


//	  cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);//,map_uniform_CS,map_theta_CS,map_phi_CS,map_thetaphi_CS);
//	  hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX);//,map_uniform_HX,map_theta_HX,map_phi_HX,map_thetaphi_HX);

//	  cs->setVal("nPrompt",numEvents_); cs->fix("nPrompt");
//	  hx->setVal("nPrompt",numEvents_); hx->fix("nPrompt");

//	  cs->setVal("nNonPrompt",numEvents_); cs->fix("nNonPrompt");

	double injected_CS_theta = 0;
	double injected_CS_phi = 0;
	double injected_CS_thetaphi = 0;
	double injected_HX_theta = 0;
	double injected_HX_phi = 0;
	double injected_HX_thetaphi = 0;

    if(RealDataPRNPstats){
    	if(Prompt){
    		injected_CS_theta=jpsi::inj_PR_CS_th[yBin][ptBin];
    		injected_CS_phi=jpsi::inj_PR_CS_ph[yBin][ptBin];
    		injected_CS_thetaphi=jpsi::inj_PR_CS_thph[yBin][ptBin];
    		injected_HX_theta=jpsi::inj_PR_HX_th[yBin][ptBin];
    		injected_HX_phi=jpsi::inj_PR_HX_ph[yBin][ptBin];
    		injected_HX_thetaphi=jpsi::inj_PR_HX_thph[yBin][ptBin];
    	}
    	if(!Prompt){
    		injected_CS_theta=jpsi::inj_NP_CS_th[yBin][ptBin];
    		injected_CS_phi=jpsi::inj_NP_CS_ph[yBin][ptBin];
    		injected_CS_thetaphi=jpsi::inj_NP_CS_thph[yBin][ptBin];
    		injected_HX_theta=jpsi::inj_NP_HX_th[yBin][ptBin];
    		injected_HX_phi=jpsi::inj_NP_HX_ph[yBin][ptBin];
    		injected_HX_thetaphi=jpsi::inj_NP_HX_thph[yBin][ptBin];
    	}
    }

    cout<<"thetaCS "<<injected_CS_theta<<endl;
    cout<<"phiCS "<<injected_CS_phi<<endl;
    cout<<"thetaphiCS "<<injected_CS_thetaphi<<endl;
    cout<<"thetaHX "<<injected_HX_theta<<endl;
    cout<<"phiHX "<<injected_HX_phi<<endl;
    cout<<"thetaphiHX "<<injected_HX_thetaphi<<endl;

/*	cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",injected_CS_phi);
	cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",injected_CS_theta);
	cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",injected_CS_thetaphi);
	hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",injected_HX_phi);
	hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",injected_HX_theta);
	hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",injected_HX_thetaphi);
*/
/*	cs->getPromptPolarizationModel()->setVal("nonpromptlambda_phi_CS",-injected_CS_phi);
	cs->getPromptPolarizationModel()->setVal("nonpromptlambda_theta_CS",-injected_CS_theta);
	cs->getPromptPolarizationModel()->setVal("nonpromptlambda_thetaphi_CS",-injected_CS_thetaphi);
	hx->getPromptPolarizationModel()->setVal("nonpromptlambda_phi_HX",-injected_HX_phi);
	hx->getPromptPolarizationModel()->setVal("nonpromptlambda_theta_HX",-injected_HX_theta);
	hx->getPromptPolarizationModel()->setVal("nonpromptlambda_thetaphi_HX",-injected_HX_thetaphi);
*/
//	costh_CS.setBins(40);
//	phi_CS.setBins(36);

	RooArgList* varlist2 = new RooArgList(costh_CS,phi_CS);

	RooRealVar coef("coef","coef",1);
	RooRealVar coef2("coef2","coef2",1);
	RooArgList* varlist3 = new RooArgList(coef,coef2);



	  RooRealVar lambda_th_CS("lambda_th_CS","lambda_th_CS",0,-2,2);
	  RooRealVar lambda_ph_CS("lambda_ph_CS","lambda_ph_CS",0,-2,2);
	  RooRealVar lambda_thph_CS("lambda_thph_CS","lambda_thph_CS",0,-2,2);

	  RooPolarizationPdf CSpolPDF("CSpolPDF","CSpolPDF",costh_CS,phi_CS,lambda_th_CS,lambda_ph_CS,lambda_thph_CS,map_uniform_CS,map_theta_CS,map_phi_CS,map_thetaphi_CS);

	  RooRealVar lambda_th_HX("lambda_th_HX","lambda_th_HX",0,-2,2);
	  RooRealVar lambda_ph_HX("lambda_ph_HX","lambda_ph_HX",0,-2,2);
	  RooRealVar lambda_thph_HX("lambda_thph_HX","lambda_thph_HX",0,-2,2);

	  RooPolarizationPdf HXpolPDF("HXpolPDF","HXpolPDF",costh_HX,phi_HX,lambda_th_HX,lambda_ph_HX,lambda_thph_HX,map_uniform_HX,map_theta_HX,map_phi_HX,map_thetaphi_HX);






//	RooAbsData *thisBinGenCS = CSpolPDF.generate(RooArgSet(costh_CS,phi_CS),numEvents_);//thisBin->tree()->GetEntries());
//	RooAbsData *thisBinGenHX = HXpolPDF.generate(RooArgSet(costh_HX,phi_HX),numEvents_);//thisBin->tree()->GetEntries());




	  RooRealVar generation("generation","generation",iteration-0.5,iteration+0.5);
	  RooRealVar MCType("MCType","MCType",-0.5,4);
	  if(Prompt) RooRealVar MCType("MCType","MCType",-0.5,+0.5);
	  if(!Prompt) {
	RooRealVar MCType("MCType","MCType",+0.5,+1.5);
	Jpsict.setMin(0.1);
	}
	  RooArgSet PseudoVarlistCS(costh_CS,phi_CS,generation,MCType,Jpsict);
	  RooArgSet PseudoVarlistHX(costh_HX,phi_HX,generation,MCType,Jpsict);

	  TTree* thisBinGentree = (TTree*)PseudoInput->Get("data");

	  RooDataSet *thisBinGenCS = new RooDataSet("thisBinGenCS","thisBinGenCS",PseudoVarlistCS,Import(*thisBinGentree));
	  RooDataSet *thisBinGenHX = new RooDataSet("thisBinGenHX","thisBinGenHX",PseudoVarlistHX,Import(*thisBinGentree));



	  lambda_th_CS.setVal(0);
	  lambda_ph_CS.setVal(0);
	  lambda_thph_CS.setVal(0);
	  lambda_th_HX.setVal(0);
	  lambda_ph_HX.setVal(0);
	  lambda_thph_HX.setVal(0);





	costh_HX.setBins(200);
	phi_HX.setBins(180);

/*	char TH2name[100];
	sprintf(TH2name,"GenData_CS_rap%d_pt%d_generation%d",yBin+1,ptBin+1,iteration);
	TH2F * thisBinGenCSTH = (TH2F*)thisBinGenCS->createHistogram(TH2name,costh_CS,Binning(200),YVar(phi_CS,Binning(180)));
	sprintf(TH2name,"GenData_HX_rap%d_pt%d_generation%d",yBin+1,ptBin+1,iteration);
	TH2F * thisBinGenHXTH = (TH2F*)thisBinGenHX->createHistogram(TH2name,costh_HX,Binning(200),YVar(phi_HX,Binning(180)));


	DataGenFile->Add(thisBinGenCSTH);
	DataGenFile->Add(thisBinGenHXTH);
	if(iteration==1) DataGenFile->Add(thisBinGenHX);
	if(iteration==1) DataGenFile->Add(thisBinGenCS);


	continue;
*/


//	RooDataHist *thisBinGenCS_hist = new RooDataHist("thisBinGenCS_hist","thisBinGenCS_hist",RooArgSet(costh_CS,phi_CS),*thisBinGenCS);

	cout<<"Generation Character CS: "<<thisBinGenCS->mean(costh_CS)<<" "<<thisBinGenCS->mean(phi_CS)<<endl;
	cout<<"Generation Character HX: "<<thisBinGenHX->mean(costh_HX)<<" "<<thisBinGenHX->mean(phi_HX)<<endl;



/*	  RooRealVar mean("mean","mean",.5);
	  RooRealVar sigma("sigma","sigma",2.5);

	  RooGaussian* gausstry2 = new RooGaussian("gausstry2","gausstry2",JpsiMass,mean,sigma);
	  TDirectory *direct;
	  direct->Add(gausstry2);

	  delete gausstry2;
*/


/*

		RooPlot* AccFramePR_costh_HX = new RooPlot;
		AccFramePR_costh_HX = costh_HX.frame() ;
		thisBinGenHX->plotOn(AccFramePR_costh_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
		hx->model()->plotOn(AccFramePR_costh_HX,LineWidth(2),Normalization(1.0));

		RooPlot* AccFramePR_phi_HX = new RooPlot;
		AccFramePR_phi_HX = phi_HX.frame() ;
		thisBinGenHX->plotOn(AccFramePR_phi_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
		hx->model()->plotOn(AccFramePR_phi_HX,LineWidth(2),Normalization(1.0));



		TCanvas* Acc1DCanvas = new TCanvas("Acc1DCanvas","Acc1DCanvas",2000,1000);
		Acc1DCanvas->SetFillColor(kWhite); Acc1DCanvas->Divide(2);
		Acc1DCanvas->cd(1) ; gPad->SetFillColor(kWhite); AccFramePR_costh_HX->Draw();
		Acc1DCanvas->cd(2) ; gPad->SetFillColor(kWhite); AccFramePR_phi_HX->Draw();

		char FilenameAcc[200];
		sprintf(FilenameAcc,"Plots/AcceptanceMaps/AccCompTry.png");

		Acc1DCanvas->SaveAs(FilenameAcc);
		Acc1DCanvas->Close();

*/


////////////// FIT MODEL ////////////////////////////////////////////////////////////////////////////////////

	if(rebuildmodel){

/*		delete cs;
		delete hx;

		hx = new CompositeModelBuilder("HX","");
		cs = new CompositeModelBuilder("CS","");

		hx->setUsePrompt(doprompt);
		hx->setUseNonPrompt(dononprompt);
		hx->setUseBkg(dobkg);
		hx->setUsePol(dopol);
		hx->setUseAcceptanceMaps(useAccMaps);
		hx->setUseMass(false);
		hx->setUseLifetime(false);

		cs->setUsePrompt(doprompt);
		cs->setUseNonPrompt(dononprompt);
		cs->setUseBkg(dobkg);
		cs->setUsePol(dopol);
		cs->setUseAcceptanceMaps(useAccMaps);
		cs->setUseMass(false);
		cs->setUseLifetime(false);
*/

//	current_in_model = fitInputmodel->GetDirectory(binName.str().c_str());

//	cs->loadParameters(*current_in_model);
//	hx->loadParameters(*current_in_model);

//	cs->getMassModel()->fix("CBm");
//	cs->getMassModel()->fix("CBs");
//	cs->setVal("nPrompt",numEvents_); cs->fix("nPrompt");


//	hx->getMassModel()->fix("CBm");
//	hx->getMassModel()->fix("CBs");
//	hx->setVal("nPrompt",numEvents_); hx->fix("nPrompt");

//	  cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS,map_uniformExtract_CS,map_thetaExtract_CS,map_phiExtract_CS,map_thetaphiExtract_CS);
//	  hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX,map_uniformExtract_HX,map_thetaExtract_HX,map_phiExtract_HX,map_thetaphiExtract_HX);


	}

//	cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",0);
//	cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",0);
//	cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",0);
//	hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",0);
//	hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",0);
//	hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",0);




/*
	RooPlot* AccFramePR_costh_HX = new RooPlot;
	AccFramePR_costh_HX = costh_HX.frame() ;
	thisBinGenHX->plotOn(AccFramePR_costh_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	hx->model()->plotOn(AccFramePR_costh_HX,LineWidth(2),Normalization(1.0));

	RooPlot* AccFramePR_phi_HX = new RooPlot;
	AccFramePR_phi_HX = phi_HX.frame() ;
	thisBinGenHX->plotOn(AccFramePR_phi_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	hx->model()->plotOn(AccFramePR_phi_HX,LineWidth(2),Normalization(1.0));



	TCanvas* Acc1DCanvas = new TCanvas("Acc1DCanvas","Acc1DCanvas",2000,1000);
	Acc1DCanvas->SetFillColor(kWhite); Acc1DCanvas->Divide(2);
	Acc1DCanvas->cd(1) ; gPad->SetFillColor(kWhite); AccFramePR_costh_HX->Draw();
	Acc1DCanvas->cd(2) ; gPad->SetFillColor(kWhite); AccFramePR_phi_HX->Draw();

	char FilenameAcc[200];
	sprintf(FilenameAcc,"Plots/AcceptanceMaps/AccCompTry.png");

	Acc1DCanvas->SaveAs(FilenameAcc);
	Acc1DCanvas->Close();

*/



//////////////////////////////// Set Integration method /////////////////////////////////////////




		RooNumIntConfig* conf = RooAbsReal::defaultIntegratorConfig();
		conf->setEpsAbs(1e-13);
		conf->setEpsRel(1e-13);



/////////////////////////////////////////////////////////////////////////////////////////////////




	if(fitCS){

	std::cout << "Collins-Soper model before fit:" << std::endl;
//	cs->Print();

	std::cout << "Fitting CS..." << std::endl;
//	  if(pereverr) {

		  if(iteration>generate_gt){

//			 RooRealVar* lambda_th=cs->getPromptPolarizationModel()->lambdathetaCS();
//		 	 RooRealVar* lambda_ph=cs->getPromptPolarizationModel()->lambdaphiCS();
//		 	 RooRealVar* lambda_thph=cs->getPromptPolarizationModel()->lambdathetaphiCS();


//		 	 RooPolarizationConstraint* constr = new RooPolarizationConstraint("constr","constr",*lambda_th,*lambda_ph,*lambda_thph);

//		 	 RooArgSet constrSet;
//		 	 constrSet.add(*constr);

	 		 RooAbsReal* nll=CSpolPDF.createNLL(*thisBinGenCS/*,RooFit::Range("signalRegionCS"),RooFit::ConditionalObservables(JpsictErr)/*,RooFit::ExternalConstraints(constrSet)*/);

//	 		 if(fluc) RooAbsReal* nll=final->createNLL(*thisBinGenCS,RooFit::Range("signalRegionCS"),RooFit::ConditionalObservables(JpsictErr),RooFit::ExternalConstraints(constrSet));

//	 		 RooNLLVar *nll_CS = new RooNLLVar("nll_CS","nll_CS",*cs->model(),*thisBinGenCS);

	 		 RooMinuit* minuit = new RooMinuit(*nll);

//	 		 RooMinimizer* minuit = new RooMinimizer(*nll);

	 		 minuit->setStrategy(2);
	 		 minuit->setPrintEvalErrors(-1);
	 		 minuit->setPrintLevel(3);
//	 		 minuit->setVerbose(true);
//	 		 minuit->setEps(1);


if(initialhesse){
	 		 cout<<"InitialHesse_"<<endl;
	 		 minuit->hesse();
	 		 sprintf(resultname,"cs_initialhesse_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
	 		 RooFitResult* initialhesseresult=minuit->save(resultname,resultname);
	 		 cout<<"INITIAL HESSE RESULT:"<<endl;
	 		 cout<<initialhesseresult->edm()<<endl;
	 		 cout<<initialhesseresult->status()<<endl;
	 		 cout<<"InitialHesse-covQual"<<endl;
	 		 cout<<initialhesseresult->covQual()<<endl;
	 		 if(initialhesseresult->covQual()==3)cout<<"INITIAL HESSE CONV_ERGED"<<endl;
	 		 else cout<<"INITIAL HESSE FAIL_ED"<<endl;
	 		 initialhesseresult->Print();

	 		 resfile->Add(initialhesseresult);
	 		if(moreBins)resfile_bins->Add(initialhesseresult);

}

if(simplex){
		 		 cout<<"Simplex_"<<endl;
		 		 minuit->simplex();
		 		 sprintf(resultname,"cs_simplex_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
		 		 RooFitResult* simplexresult=minuit->save(resultname,resultname);
		 		 cout<<"SIMPLEX RESULT:"<<endl;
		 		 cout<<simplexresult->edm()<<endl;
		 		 cout<<simplexresult->status()<<endl;
		 		 cout<<"Simplex-covQual"<<endl;
		 		 cout<<simplexresult->covQual()<<endl;
		 		 if(simplexresult->covQual()==3)cout<<"SIMPLEX CONV_ERGED"<<endl;
		 		 else cout<<"SIMPLEX FAIL_ED"<<endl;
		 		 simplexresult->Print();

		 		 resfile->Add(simplexresult);
		 		if(moreBins)resfile_bins->Add(simplexresult);

}

if(migrad){
	 		 cout<<"Migrad_"<<endl;
	 		 minuit->migrad();
	 		 sprintf(resultname,"cs_migrad_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
	 		 RooFitResult* migradresult=minuit->save(resultname,resultname);
	 		 cout<<"MIGRAD RESULT:"<<endl;
	 		 cout<<migradresult->edm()<<endl;
	 		 cout<<migradresult->status()<<endl;
	 		 cout<<"Migrad-covQual"<<endl;
	 		 cout<<migradresult->covQual()<<endl;
	 		 if(migradresult->covQual()==3)cout<<"MIGRAD CONV_ERGED"<<endl;
	 		 else cout<<"MIGRAD FAIL_ED"<<endl;
	 		 migradresult->Print();

	 		 resfile->Add(migradresult);
	 		if(moreBins)resfile_bins->Add(migradresult);

}
if(improve){
      		 cout<<"Improve_"<<endl;
      		 minuit->improve();
      		 sprintf(resultname,"cs_improve_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
      		 RooFitResult* improveresult=minuit->save(resultname,resultname);
      		 cout<<"IMPROVE RESULT:"<<endl;
      		 cout<<improveresult->edm()<<endl;
      		 cout<<improveresult->status()<<endl;
      		 cout<<"Improve-covQual"<<endl;
      		 cout<<improveresult->covQual()<<endl;
      		 if(improveresult->covQual()==3)cout<<"IMPROVE CONV_ERGED"<<endl;
	 		 else cout<<"IMPROVE FAIL_ED"<<endl;
      		 improveresult->Print();

             resfile->Add(improveresult);
             if(moreBins)resfile_bins->Add(improveresult);

}

/*
if(minimizeBFGS){

		 cout<<"Bfgs_"<<endl;
 		 minuit->minimize(minimizerType,minimizerAlgo);
 		 sprintf(resultname,"cs_bfgs_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
 		 RooFitResult* bfgsresult=minuit->save(resultname,resultname);
 		 cout<<"BFGS RESULT:"<<endl;
 		 cout<<bfgsresult->edm()<<endl;
 		 cout<<bfgsresult->status()<<endl;
 		 cout<<"Bfgs-covQual"<<endl;
 		 cout<<bfgsresult->covQual()<<endl;
 		 if(bfgsresult->covQual()==3)cout<<"BFGS CONV_ERGED"<<endl;
		 else cout<<"BFGS FAIL_ED"<<endl;
 		bfgsresult->Print();

        resfile->Add(bfgsresult);

}
*/
RooFitResult* hesseresult;
if(hesse){
	 		 cout<<"Hesse_"<<endl;
	 		 minuit->hesse();
	 		 sprintf(resultname,"cs_hesse_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
	 		 hesseresult=minuit->save(resultname,resultname);
	 		 cout<<"HESSE RESULT:"<<endl;
	 		 cout<<hesseresult->edm()<<endl;
	 		 cout<<hesseresult->status()<<endl;
	 		 cout<<"Hesse-covQual"<<endl;
	 		 cout<<hesseresult->covQual()<<endl;
	 		 if(hesseresult->covQual()==3)cout<<"HESSE CONV_ERGED"<<endl;
	 		 else cout<<"HESSE FAIL_ED"<<endl;
	 		 hesseresult->Print();

	 		 resfile->Add(hesseresult);
	 		if(moreBins)resfile_bins->Add(hesseresult);

}

if(minos){
	 		 cout<<"Minos_"<<endl;
	 		 minuit->minos();
	 		 sprintf(resultname,"cs_minos_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
	 		 RooFitResult* minosresult=minuit->save(resultname,resultname);
	 		 cout<<"MINOS RESULT:"<<endl;
	 		 cout<<minosresult->edm()<<endl;
	 		 cout<<minosresult->status()<<endl;
	 		 cout<<"Minos-covQual"<<endl;
	 		 cout<<minosresult->covQual()<<endl;
	 		 if(minosresult->covQual()==3)cout<<"MINOS CONV_ERGED"<<endl;
	 		 else cout<<"MINOS FAIL_ED"<<endl;
	 		 minosresult->Print();

	 		 resfile->Add(minosresult);
	 		if(moreBins)resfile_bins->Add(minosresult);

}

	 					delete minuit;
//	 					delete constr;



















		  }
//	  }

	    else{
	    }


	std::cout << "Collins-Soper model after fit:" << std::endl;
//	cs->Print();


	}

	if(fitHX){

	std::cout << std::endl << "Helicity Frame model before fit:" << std::endl;
//	hx->Print();

	std::cout << "Fitting HX..." << std::endl;
//	  if(pereverr){
		  if(iteration>generate_gt){


//			     RooRealVar* lambda_th=hx->getPromptPolarizationModel()->lambdathetaHX();
//			 	 RooRealVar* lambda_ph=hx->getPromptPolarizationModel()->lambdaphiHX();
//			 	 RooRealVar* lambda_thph=hx->getPromptPolarizationModel()->lambdathetaphiHX();



//			 	 RooPolarizationConstraint* constr = new RooPolarizationConstraint("constr","constr",*lambda_th,*lambda_ph,*lambda_thph);

//			 	 RooArgSet constrSet;

//			 	 constrSet.add(*constr);


		 		 RooAbsReal* nll=HXpolPDF.createNLL(*thisBinGenHX/*,RooFit::Range("signalRegionHX"),RooFit::ConditionalObservables(JpsictErr)/*,RooFit::ExternalConstraints(constrSet)*/);
//		 		 RooAbsReal* nll=hx->model()->createNLL(*compositePDFData,RooFit::Range("signalRegionHX"),RooFit::ConditionalObservables(JpsictErr)/*,RooFit::ExternalConstraints(constrSet)*/);

//		 		 RooNLLVar *nll_HX = new RooNLLVar("nll_HX","nll_HX",*hx->model(),*thisBinGenHX);

		 		 RooMinuit* minuit = new RooMinuit(*nll);
		 		 minuit->setStrategy(2);
		 		 minuit->setPrintEvalErrors(-1);
		 		 minuit->setPrintLevel(3);
//		 		 minuit->setVerbose(true);


if(initialhesse){
	 		     cout<<"InitialHesse_"<<endl;
	 		     minuit->hesse();
	 		     sprintf(resultname,"hx_initialhesse_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
	 		     RooFitResult* initialhesseresult=minuit->save(resultname,resultname);
	 		     cout<<"INITIAL HESSE RESULT:"<<endl;
	 		     cout<<initialhesseresult->edm()<<endl;
	 		     cout<<initialhesseresult->status()<<endl;
	 		     cout<<"InitialHesse-covQual"<<endl;
	 		     cout<<initialhesseresult->covQual()<<endl;
	 		     if(initialhesseresult->covQual()==3)cout<<"INITIAL HESSE CONV_ERGED"<<endl;
		 		 else cout<<"INITIAL HESSE FAIL_ED"<<endl;
	 		     initialhesseresult->Print();

	 		     resfile->Add(initialhesseresult);
	 		    if(moreBins)resfile_bins->Add(initialhesseresult);

}

if(simplex){
		 		 cout<<"Simplex_"<<endl;
		 		 minuit->simplex();
		 		 sprintf(resultname,"hx_simplex_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
		 		 RooFitResult* simplexresult=minuit->save(resultname,resultname);
		 		 cout<<"SIMPLEX RESULT:"<<endl;
		 		 cout<<simplexresult->edm()<<endl;
		 		 cout<<simplexresult->status()<<endl;
		 		 cout<<"Simplex-covQual"<<endl;
		 		 cout<<simplexresult->covQual()<<endl;
		 		 if(simplexresult->covQual()==3)cout<<"SIMPLEX CONV_ERGED"<<endl;
		 		 else cout<<"SIMPLEX FAIL_ED"<<endl;
		 		 simplexresult->Print();

		 		 resfile->Add(simplexresult);
		 		if(moreBins)resfile_bins->Add(simplexresult);

}

if(migrad){
		 		 cout<<"Migrad_"<<endl;
		 		 minuit->migrad();
		 		 sprintf(resultname,"hx_migrad_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
		 		 RooFitResult* migradresult=minuit->save(resultname,resultname);
		 		 cout<<"MIGRAD RESULT:"<<endl;
		 		 cout<<migradresult->edm()<<endl;
		 		 cout<<migradresult->status()<<endl;
		 		 cout<<"Migrad-covQual"<<endl;
		 		 cout<<migradresult->covQual()<<endl;
		 		 if(migradresult->covQual()==3)cout<<"MIGRAD CONV_ERGED"<<endl;
		 		 else cout<<"MIGRAD FAIL_ED"<<endl;
		 		 migradresult->Print();

		 		 resfile->Add(migradresult);
			 		if(moreBins)resfile_bins->Add(migradresult);

}
if(improve){
				 cout<<"Improve_"<<endl;
				 minuit->improve();
				 sprintf(resultname,"hx_improve_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
				 RooFitResult* improveresult=minuit->save(resultname,resultname);
				 cout<<"IMPROVE RESULT:"<<endl;
				 cout<<improveresult->edm()<<endl;
				 cout<<improveresult->status()<<endl;
				 cout<<"Improve-covQual"<<endl;
				 cout<<improveresult->covQual()<<endl;
				 if(improveresult->covQual()==3)cout<<"IMPROVE CONV_ERGED"<<endl;
		 		 else cout<<"IMPROVE FAIL_ED"<<endl;
				 improveresult->Print();

				 resfile->Add(improveresult);
				 if(moreBins)resfile_bins->Add(improveresult);

}

if(hesse){
		 		 cout<<"Hesse_"<<endl;
		 		 minuit->hesse();
		 		 sprintf(resultname,"hx_hesse_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
		 		 RooFitResult* hesseresult=minuit->save(resultname,resultname);
		 		 cout<<"HESSE RESULT:"<<endl;
		 		 cout<<hesseresult->edm()<<endl;
		 		 cout<<hesseresult->status()<<endl;
		 		 cout<<"Hesse-covQual"<<endl;
		 		 cout<<hesseresult->covQual()<<endl;
		 		 if(hesseresult->covQual()==3)cout<<"HESSE CONV_ERGED"<<endl;
		 		 else cout<<"HESSE FAIL_ED"<<endl;
		 		 hesseresult->Print();

		 		 resfile->Add(hesseresult);
		 		if(moreBins)resfile_bins->Add(hesseresult);

}

if(minos){
		 		 cout<<"Minos_"<<endl;
		 		 minuit->minos();
		 		 sprintf(resultname,"hx_minos_rap%d_pt%d_gen%d",yBin+1,ptBin+1,iteration);
		 		 RooFitResult* minosresult=minuit->save(resultname,resultname);
		 		 cout<<"MINOS RESULT:"<<endl;
		 		 cout<<minosresult->edm()<<endl;
		 		 cout<<minosresult->status()<<endl;
		 		 cout<<"Minos-covQual"<<endl;
		 		 cout<<minosresult->covQual()<<endl;
		 		 if(minosresult->covQual()==3)cout<<"MINOS CONV_ERGED"<<endl;
		 		 else cout<<"MINOS FAIL_ED"<<endl;
		 		 minosresult->Print();

		 		 resfile->Add(minosresult);
		 		if(moreBins)resfile_bins->Add(minosresult);

}

		 		 delete minuit;
//				 delete constr;






		  }

//	  }

	  else{
	  }



	std::cout << std::endl << "Helicity Frame model after fit:" << std::endl;
//	hx->Print();

	}

//	current_out_hx = output_hx->mkdir(binNameOut.str().c_str());
//	current_out_cs = output_cs->mkdir(binNameOut.str().c_str());


//	cs->saveParameter("nPrompt",*current_out_cs);
//	cs->saveParameter("nNonPrompt",*current_out_cs);
//	cs->saveParameter("nBackground",*current_out_cs);

//	hx->saveParameter("nPrompt",*current_out_hx);
//	hx->saveParameter("nNonPrompt",*current_out_hx);
//	hx->saveParameter("nBackground",*current_out_hx);

//	cs->saveParameters(*current_out_cs);
//	hx->saveParameters(*current_out_hx);














/*

if(ploterrorband && iteration>generate_gt){


		char JobID[200];
		sprintf(JobID,"multi_MI_GARETE_trigstat");
//		sprintf(JobID,"multi_MI_GEN_plus05");
//		sprintf(JobID,"multi_MI_GA_minus05");

		char miniCrit[200];
		sprintf(miniCrit,"migrad");
		char miniCrit2[200];
		sprintf(miniCrit2,"hesse");

		char inputfilename[200];
		sprintf(inputfilename,"CRAB/%s/RooFitResult_rap%d_pt%d.root",JobID,yBin+1,ptBin+1);

		TFile *resfile2 = new TFile(inputfilename,"UPDATE");
		resfile2->Print();

		char dirStruct[200];
		sprintf(dirStruct,"Plots/ToyMC/%s/ErrorPlots_%dsigma",JobID,n_errsigma);

		gSystem->mkdir(dirStruct);


		RooFitResult* roofitres;
		RooFitResult* roofitres2;

		char canvasname[200];
		sprintf(canvasname,"canvas_rap%d_pt%d_iter%d",yBin+1,ptBin+1,iteration);

		TCanvas* costh_errorband_canvas = new TCanvas(canvasname,canvasname,3600,3200);
		costh_errorband_canvas->Divide(2,2);  costh_errorband_canvas->SetFillColor(kWhite);

		RooPlot* costh_errorband_cs = new RooPlot;
		RooPlot* costh_errorband_hx = new RooPlot;
		RooPlot* phi_errorband_cs = new RooPlot;
		RooPlot* phi_errorband_hx = new RooPlot;
//		RooPlot* param_CS = new RooPlot;
//		RooPlot* param_HX = new RooPlot;




	    for(int frameDex = 0; frameDex<2; ++frameDex){

	  char resname[200];
	  char res2name[200];

	  if(frameDex==0){sprintf(resname,"cs_%s_rap%d_pt%d_gen%d",miniCrit,yBin+1,ptBin+1,iteration); sprintf(res2name,"cs_%s_rap%d_pt%d_gen%d",miniCrit2,yBin+1,ptBin+1,iteration);}
	  if(frameDex==1){sprintf(resname,"hx_%s_rap%d_pt%d_gen%d",miniCrit,yBin+1,ptBin+1,iteration); sprintf(res2name,"hx_%s_rap%d_pt%d_gen%d",miniCrit2,yBin+1,ptBin+1,iteration);}

	  roofitres = (RooFitResult*) resfile2->Get(resname);
	  roofitres2 = (RooFitResult*) resfile2->Get(res2name);

	  bool onlyCrit1converged(false);
	  bool onlyCrit2converged(false);
	  bool bothconverged(false);

	  bool criterium(false);

	  if(roofitres->covQual() == 3 || roofitres2->covQual() == 3) {criterium=true; if(roofitres2->covQual() != 3) onlyCrit1converged=true;if(roofitres->covQual() != 3) onlyCrit2converged=true;}
	  if(roofitres->covQual() == 3 && roofitres2->covQual() == 3) bothconverged=true;
	  RooArgList params = roofitres2->floatParsFinal();

	  char resultCrit[200];
	  sprintf(resultCrit,"%s",miniCrit2);
	  char convergenceString[200];

	  if(!criterium) {sprintf(convergenceString," (not converged)"); continue;}
	  if(onlyCrit1converged) {sprintf(convergenceString," (only %s converged)",miniCrit);}
	  if(onlyCrit2converged) {sprintf(convergenceString," (only %s converged)",miniCrit2);}
	  if(bothconverged) {sprintf(convergenceString," (%s and %s converged)",miniCrit,miniCrit2);}


	  if(onlyCrit1converged){
		  RooArgList params = roofitres->floatParsFinal();
		  sprintf(resultCrit,"%s",miniCrit);
	  }
	  RooRealVar* lambdaphi_;
	  RooRealVar* lambdatheta_;
	  RooRealVar* lambdathetaphi_;

	  lambdaphi_ =  (RooRealVar*)params.at(0);
	  lambdatheta_ =  (RooRealVar*)params.at(1);
	  lambdathetaphi_ =  (RooRealVar*)params.at(2);

	  if(frameDex==0){
		cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",lambdaphi_->getVal());
		cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",lambdatheta_->getVal());
		cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",lambdathetaphi_->getVal());
		cs->getPromptPolarizationModel()->setErr("promptlambda_phi_CS",lambdaphi_->getError());
		cs->getPromptPolarizationModel()->setErr("promptlambda_theta_CS",lambdatheta_->getError());
		cs->getPromptPolarizationModel()->setErr("promptlambda_thetaphi_CS",lambdathetaphi_->getError());
	  }
	  if(frameDex==1){
		hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",lambdaphi_->getVal());
		hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",lambdatheta_->getVal());
		hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",lambdathetaphi_->getVal());
		hx->getPromptPolarizationModel()->setErr("promptlambda_phi_HX",lambdaphi_->getError());
		hx->getPromptPolarizationModel()->setErr("promptlambda_theta_HX",lambdatheta_->getError());
		hx->getPromptPolarizationModel()->setErr("promptlambda_thetaphi_HX",lambdathetaphi_->getError());
	  }


//	 RooRealVar* lamthCS_=cs->getPromptPolarizationModel()->lambdathetaCS();
//	 RooRealVar* lamphCS_=cs->getPromptPolarizationModel()->lambdaphiCS();
//	 RooRealVar* lamthphCS_=cs->getPromptPolarizationModel()->lambdathetaphiCS();

//	 lamthCS_->SetName("#lambda_{#theta CS}");
//	 lamphCS_->SetName("#lambda_{#phi CS}");
//	 lamthphCS_->SetName("#lambda_{#theta #phi CS}");




//	  cs->getPromptPolarizationModel()->lambdathetaCS()->SetName("#lambda_{#theta CS}");
//	  cs->getPromptPolarizationModel()->lambdaphiCS()->SetName("#lambda_{#phi CS}");
//	  cs->getPromptPolarizationModel()->lambdathetaphiCS()->SetName("#lambda_{#theta #phi CS}");

		  char frametitle[200];
		  double titlesize = 0.0125;

			sprintf(frametitle,"%s fit result with %d  #sigma error band%s",resultCrit,n_errsigma,convergenceString);

			if(frameDex==0) {
			costh_errorband_cs = costh_CS.frame() ;
			costh_errorband_cs->SetTitle(0);
			thisBinGenCS->plotOn(costh_errorband_cs,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			cs->model()->plotOn(costh_errorband_cs,VisualizeError(*roofitres,n_errsigma,false),RooFit::FillColor(kGreen-6));
			thisBinGenCS->plotOn(costh_errorband_cs,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			cs->model()->plotOn(costh_errorband_cs,RooFit::LineColor(kGreen+3));

			phi_errorband_cs = phi_CS.frame() ;
			phi_errorband_cs->SetTitle(0);
			thisBinGenCS->plotOn(phi_errorband_cs,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			cs->model()->plotOn(phi_errorband_cs,VisualizeError(*roofitres,n_errsigma,false),RooFit::FillColor(kGreen-6));
			thisBinGenCS->plotOn(phi_errorband_cs,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			cs->model()->plotOn(phi_errorband_cs,RooFit::LineColor(kGreen+3));

			cs->model()->paramOn(costh_errorband_cs,Format("NE",AutoPrecision(2)),Label(frametitle),Layout(0.10,0.9,1));
			costh_errorband_cs->getAttText()->SetTextSize(0.025);

			}
			if(frameDex==1) {
			costh_errorband_hx = costh_HX.frame() ;
			costh_errorband_hx->SetTitle(0);
			thisBinGenHX->plotOn(costh_errorband_hx,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			hx->model()->plotOn(costh_errorband_hx,VisualizeError(*roofitres,n_errsigma,false),RooFit::FillColor(kGreen-6));
			thisBinGenHX->plotOn(costh_errorband_hx,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			hx->model()->plotOn(costh_errorband_hx,RooFit::LineColor(kGreen+3));

			phi_errorband_hx = phi_HX.frame() ;
			phi_errorband_hx->SetTitle(0);
			thisBinGenHX->plotOn(phi_errorband_hx,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			hx->model()->plotOn(phi_errorband_hx,VisualizeError(*roofitres,n_errsigma,false),RooFit::FillColor(kGreen-6));
			thisBinGenHX->plotOn(phi_errorband_hx,DataError(RooAbsData::SumW2),MarkerSize(0.4));
			hx->model()->plotOn(phi_errorband_hx,RooFit::LineColor(kGreen+3));

			hx->model()->paramOn(costh_errorband_hx,Format("NE",AutoPrecision(2)),Label(frametitle),Layout(0.1,0.9,1));
			costh_errorband_hx->getAttText()->SetTextSize(0.025);

			}








	    }





//		 TCanvas *PromptPullCanvas = new TCanvas("PromptCanvas","",1000,1000);

		  TPad* toppad = new TPad("toppad","",0,0.5,1,1);
		  toppad->SetFillColor(0);

		  TPad* bottompad = new TPad("bottompad","",0,0,1,0.5);
		  bottompad->SetFillColor(0);

		  toppad->Draw();
		  toppad->cd();
	//	  toppad->SetBottomMargin(0.65);
		  gPad->SetFillColor(kWhite);
		  costh_errorband_cs->Draw();

		  bottompad->Draw();
		  bottompad->cd();
	//	  bottompad->SetTopMargin(0.15);
		  gPad->SetFillColor(kWhite);
		  phi_errorband_cs->Draw();










		costh_errorband_canvas->cd(1) ; gPad->SetFillColor(kWhite); gPad->SetTopMargin(0.3); costh_errorband_cs->GetYaxis()->SetTitleOffset(1.5) ; costh_errorband_cs->Draw();
	 	costh_errorband_canvas->cd(3) ; gPad->SetFillColor(kWhite); gPad->SetBottomMargin(0.3); phi_errorband_cs->GetYaxis()->SetTitleOffset(1.5) ; phi_errorband_cs->Draw();
	 	costh_errorband_canvas->cd(2) ; gPad->SetFillColor(kWhite); gPad->SetTopMargin(0.3); costh_errorband_hx->GetYaxis()->SetTitleOffset(1.5) ; costh_errorband_hx->Draw();
	 	costh_errorband_canvas->cd(4) ; gPad->SetFillColor(kWhite); gPad->SetBottomMargin(0.3); phi_errorband_hx->GetYaxis()->SetTitleOffset(1.5) ;phi_errorband_hx->Draw();

//		costh_errorband_canvas->cd(5) ; gPad->SetFillColor(kWhite); gPad->SetBottomMargin(0.4); gPad->SetTopMargin(2); gPad->SetLeftMargin(0.4); gPad->SetRightMargin(0.4); param_CS->Draw();
//		costh_errorband_canvas->cd(6) ; gPad->SetFillColor(kWhite); gPad->SetBottomMargin(0.4); gPad->SetTopMargin(2); gPad->SetLeftMargin(0.4); gPad->SetRightMargin(0.4); param_HX->Draw();

		char costh_errorband_name[200];
		sprintf(costh_errorband_name,"%s/costh_errorband_rapidity%d_pt%d_iter%d.png",dirStruct,yBin+1,ptBin+1,iteration);

		costh_errorband_canvas->SaveAs(costh_errorband_name);
		costh_errorband_canvas->Close();




//	   delete costh_errorband_canvas;
	   delete phi_errorband_cs;
	   delete costh_errorband_cs;
	   delete phi_errorband_hx;
	   delete costh_errorband_hx;
//	   delete param_CS;
//	   delete param_HX;
//	   delete toppad;
//	   delete bottompad;

	   resfile2->Close();


}




*/












//	delete cs;
//	delete hx;
	delete thisBinGenCS;
	delete thisBinGenHX;

//      }
//       else cout<<"skippedFit"<<endl;
}

 			 resfile_bins->Write();
 			resfile_bins->Close();
 			  PseudoInput->Close();

 	  }








  }

	DataGenFile->Write();
	DataGenFile->Close();

	ErrBandPlots->Write();
	ErrBandPlots->Close();


  resfile->Write();
//  cout<<"RooFitResults written"<<endl;
 resfile->Close();
//  delete resfile;
//  output_hx->Write();
//   output_cs->Write();

//  output_hx->Close();
//  output_cs->Close();
//  delete output_hx;
//   delete output_cs;
//  fitInput->Close();

//  delete fitInput;

  return 0;
}
