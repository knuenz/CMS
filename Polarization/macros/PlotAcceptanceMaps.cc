#include <iostream>
#include <sstream>

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
#include "RooHistPdf.h"
#include "RooDataHist.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TClass.h"
#include "TLegend.h"
#include "TSystem.h"


int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;






  gSystem->mkdir("Plots/AcceptanceMaps");

   gStyle->SetPalette(1);
   gStyle->SetPadLeftMargin(0.20);
//   gStyle->SetPadRightMargin(0.02);
   gStyle->SetTitleFillColor(kWhite);


  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5); 
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);
  RooRealVar MCweight("MCweight","MCweight",1);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar JpsiVprob("JpsiVprob","JpsiVprob",0.01,100000);
  RooRealVar HLT_DoubleMu0("HLT_DoubleMu0","HLT_DoubleMu0",0.5,1.5);

  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCweight);
  varlist.add(JpsictErr);
  varlist.add(JpsiVprob);
  varlist.add(HLT_DoubleMu0);

  
  char ID[200];
  sprintf(ID,"NP_garete_DM0_ATLAS_sm");

  char dirstruct[200];
  sprintf(dirstruct,"Plots/AcceptanceMaps/%s",ID);

  gSystem->mkdir(dirstruct);

//  TFile *fInput = new TFile("jPsiFit_gareteDM0_PR.root","UPDATE");

  Char_t *fileNameInPR = "/scratch/knuenz/Polarization/RootInput/TTree_prep_NP_.root";
// Char_t *fileNameInPR = "/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo.root";
//  Char_t *fileNameInPR = "/afs/hephy.at/scratch/k/knuenz/TTree_final_notrigger_MCprompt_Jpsi_Fall10_folded_.root";

  TFile* fInPR = new TFile(fileNameInPR);
  TTree* dataTreesPR = (TTree*)fInPR->Get("data");
  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",varlist,Import(*dataTreesPR),WeightVar(MCweight));


  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-3; ++yBin) {
	  for(int ptBin = 1; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {


      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];

  	std::cout << cutString.str() << std::endl;


	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);// && JpsictErr > 0.0004 && JpsictErr < 0.999999

	  RooDataSet* dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);

	  TH2F * promptMap_CS_MC = (TH2F*)dataPRbin->createHistogram("promptMap_CS_MC",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol)));
	  TH2F * promptMap_HX_MC = (TH2F*)dataPRbin->createHistogram("promptMap_HX_MC",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol)));

      char AccCSPRTitle[200];
      sprintf(AccCSPRTitle,"Prompt CS Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char AccHXPRTitle[200];
      sprintf(AccHXPRTitle,"Prompt HX Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char AccCSNPTitle[200];
      sprintf(AccCSNPTitle,"Non Prompt CS Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char AccHXNPTitle[200];
      sprintf(AccHXNPTitle,"Non Prompt HX Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char MCCSPRTitle[200];
      sprintf(MCCSPRTitle,"MC Prompt CS %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char MCHXPRTitle[200];
      sprintf(MCHXPRTitle,"MC Prompt HX %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);


      char FilenameAcc[200];









   bool rebin(false);

//		  TFile *accFile = new TFile("/scratch/knuenz/Polarization/RootInput/geomAccHistos_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded.root","UPDATE");
//		  TFile *recoEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/recoEffHistos_ATLASPT_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
//		  TFile *trigEffFile = new TFile("/scratch/knuenz/Polarization/RootInput/trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");

//		  TFile *accFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/geomAccHistos_NP_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root","UPDATE");
//		  TFile *recoEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/recoEffHistos_NP_ATLASPT_19March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
//		  TFile *trigEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root","UPDATE");

//		  TFile *accFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/PRgeomAccATLASpT24Mar.root","UPDATE");
//		  TFile *recoEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/PRrecoEffATLASpT24Mar.root","UPDATE");
//		  TFile *trigEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/PRtrigEffATLASpT24Mar.root","UPDATE");

		  TFile *accFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/NPgeomAccATLASpT24Mar.root","UPDATE");
		  TFile *recoEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/NPrecoEffATLASpT24Mar.root","UPDATE");
		  TFile *trigEffFile = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/NPtrigEffATLASpT24Mar.root","UPDATE");

		      TH2F *AccMap_CS,*AccMap_HX;
		      TH2F *AccMapExtract_CS,*AccMapExtract_HX;
		      TH2F *RecoEffMap_CS,*RecoEffMap_HX;
		      TH2F *TrigEffMap_CS,*TrigEffMap_HX;

		      std::stringstream nameHXp,nameCSp;
		      nameHXp << "hAcc2D_HX_pT" << ptBin+1 << "_rap" << yBin+1;
		      nameCSp << "hAcc2D_CS_pT" << ptBin+1 << "_rap" << yBin+1;
		      AccMap_CS = (TH2F*)accFile->Get(nameCSp.str().c_str());
		      AccMap_HX = (TH2F*)accFile->Get(nameHXp.str().c_str());
		      RecoEffMap_CS = (TH2F*)recoEffFile->Get(nameCSp.str().c_str());
		      RecoEffMap_HX = (TH2F*)recoEffFile->Get(nameHXp.str().c_str());
		      TrigEffMap_CS = (TH2F*)trigEffFile->Get(nameCSp.str().c_str());
		      TrigEffMap_HX = (TH2F*)trigEffFile->Get(nameHXp.str().c_str());

		      TH2F *Corr_CS = (TH2F*) AccMap_CS->Clone("Corr_CS");
		      TH2F *Corr_HX = (TH2F*) AccMap_HX->Clone("Corr_HX");

	if(!rebin){


	//	Corr_CS->Rebin2D(2,2);
	//	Corr_HX->Rebin2D(2,2);

			  RecoEffMap_CS->Multiply(TrigEffMap_CS);
		      Corr_CS->Multiply(RecoEffMap_CS);
		      RecoEffMap_HX->Multiply(TrigEffMap_HX);
		      Corr_HX->Multiply(RecoEffMap_HX);



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

	}


		      char Corr_CS_hist_char[200];
		      char Corr_CS_histfunc_char[200];
		      char map_uniform_CS_char[200];
		      char map_theta_CS_char[200];
		      char map_phi_CS_char[200];
		      char map_thetaphi_CS_char[200];

		      char Corr_HX_hist_char[200];
		      char Corr_HX_histfunc_char[200];
		      char map_uniform_HX_char[200];
		      char map_theta_HX_char[200];
		      char map_phi_HX_char[200];
		      char map_thetaphi_HX_char[200];


		      sprintf(Corr_CS_hist_char,"Corr_CS_hist_%d_%d",yBin+1,ptBin+1);
		      sprintf(Corr_CS_histfunc_char,"Corr_CS_histfunc_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_uniform_CS_char,"map_uniform_CS_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_theta_CS_char,"map_theta_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_phi_CS_char,"map_phi_CS_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaphi_CS_char,"map_thetaphi_CS_char%d_pt%d",yBin+1,ptBin+1);

		      sprintf(Corr_HX_hist_char,"Corr_HX_hist_%d_%d",yBin+1,ptBin+1);
		      sprintf(Corr_HX_histfunc_char,"Corr_HX_histfunc_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_uniform_HX_char,"map_uniform_HX_char_%d_%d",yBin+1,ptBin+1);
		      sprintf(map_theta_HX_char,"map_theta_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_phi_HX_char,"map_phi_HX_char%d_pt%d",yBin+1,ptBin+1);
		      sprintf(map_thetaphi_HX_char,"map_thetaphi_HX_char%d_pt%d",yBin+1,ptBin+1);

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



			CompositeModelBuilder* hx = new CompositeModelBuilder("HX","");
			CompositeModelBuilder* cs = new CompositeModelBuilder("CS","");

			hx->setUsePrompt(true);
			hx->setUseNonPrompt(false);
			hx->setUseBkg(false);
			hx->setUsePol(true);
			hx->setUseAcceptanceMaps(false);
			hx->setUseMass(false);
			hx->setUseLifetime(false);


			cs->setUsePrompt(true);
			cs->setUseNonPrompt(false);
			cs->setUseBkg(false);
			cs->setUsePol(true);
			cs->setUseAcceptanceMaps(false);
			cs->setUseMass(false);
			cs->setUseLifetime(false);



			  cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS,map_uniform_CS,map_theta_CS,map_phi_CS,map_thetaphi_CS);
			  hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX,map_uniform_HX,map_theta_HX,map_phi_HX,map_thetaphi_HX);









/*		         RooAbsPdf* promptCShistPdf_ = cs->model();
		         RooAbsPdf* promptHXhistPdf_ = hx->model();
		//         RooHistPdf* nonpromptCShistPdf_ = modelCS->getAcceptanceMaps()->nonPrompt();
		//         RooHistPdf* nonpromptHXhistPdf_ = modelHX->getAcceptanceMaps()->nonPrompt();


		    	TH1 * hist_CS_PR_model = promptCShistPdf_->createHistogram("hist_CS_PR_model",costh_CS,Binning(250),YVar(phi_CS,Binning(250))) ;

					RooRealVar norm1("norm1","norm1",1);
		     		RooDataHist accmodelCShist("accmodelCShist","accmodelCShist", RooArgList(costh_CS,phi_CS),hist_CS_PR_model);
		     		RooHistPdf accmodelCShistpdf("accmodelCShistpdf","accmodelCShistpdf",RooArgList(costh_CS,phi_CS),accmodelCShist,1);
		     		RooExtendPdf accmodelCS("accmodelCS","accmodelCS",accmodelCShistpdf,norm1);

		        TH1 * hist_HX_PR_model = promptHXhistPdf_->createHistogram("hist_HX_PR_model",costh_HX,Binning(250),YVar(phi_HX,Binning(250))) ;

		         	RooDataHist accmodelHXhist("accmodelHXhist","accmodelHXhist", RooArgList(costh_HX,phi_HX),hist_HX_PR_model);
		         	RooHistPdf accmodelHXhistpdf("accmodelHXhistpdf","accmodelHXhistpdf",RooArgList(costh_HX,phi_HX),accmodelHXhist,1);
		         	RooExtendPdf accmodelHX("accmodelHX","accmodelHX",accmodelHXhistpdf,norm1);

*/







   	  TH1 * promptCShistPdf = cs->model()->createHistogram("promptCShistPdf",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;
   	  TH1 * promptHXhistPdf = hx->model()->createHistogram("promptHXhistPdf",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol))) ;
//   	  TH1 * nonpromptCShistPdf = nonpromptCShistPdf_->createHistogram("nonpromptCShistPdf",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;
//   	  TH1 * nonpromptHXhistPdf = nonpromptHXhistPdf_->createHistogram("nonpromptHXhistPdf",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol))) ;


   	TCanvas* AcceptanceHistCanvas = new TCanvas("AcceptanceHistCanvas","AcceptanceHistCanvas",1600, 700);
   	AcceptanceHistCanvas->Divide(2);  AcceptanceHistCanvas->SetFillColor(kWhite);
   	AcceptanceHistCanvas->cd(1) ; promptCShistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptCShistPdf->SetStats(0); promptCShistPdf->GetXaxis()->SetTitle("cos #theta_{CS}"); promptCShistPdf->GetYaxis()->SetTitle("#phi_{CS} [deg]"); promptCShistPdf->SetTitle(AccCSPRTitle); promptCShistPdf->Draw("colz");
   	AcceptanceHistCanvas->cd(2) ; promptHXhistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptHXhistPdf->SetStats(0); promptHXhistPdf->GetXaxis()->SetTitle("cos #theta_{HX}"); promptHXhistPdf->GetYaxis()->SetTitle("#phi_{HX} [deg]"); promptHXhistPdf->SetTitle(AccHXPRTitle); promptHXhistPdf->Draw("colz");

   	      sprintf(FilenameAcc,"%s/AccHists_rapidity%d_pt%d.png",dirstruct,yBin+1,ptBin+1);

   	   AcceptanceHistCanvas->SaveAs(FilenameAcc);
   	AcceptanceHistCanvas->Close();


    TCanvas* AcceptanceCanvas2 = new TCanvas("AcceptanceCanvas2","AcceptanceCanvas2",1600, 1400);
    AcceptanceCanvas2->Divide(2,2);  AcceptanceCanvas2->SetFillColor(kWhite);
    AcceptanceCanvas2->cd(1) ; promptCShistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptCShistPdf->SetStats(0); promptCShistPdf->GetXaxis()->SetTitle("cos #theta_{CS}"); promptCShistPdf->GetYaxis()->SetTitle("#phi_{CS} [deg]"); promptCShistPdf->SetTitle(AccCSPRTitle); promptCShistPdf->Draw("colz");
    AcceptanceCanvas2->cd(2) ; promptHXhistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptHXhistPdf->SetStats(0); promptHXhistPdf->GetXaxis()->SetTitle("cos #theta_{HX}"); promptHXhistPdf->GetYaxis()->SetTitle("#phi_{HX} [deg]"); promptHXhistPdf->SetTitle(AccHXPRTitle); promptHXhistPdf->Draw("colz");
    AcceptanceCanvas2->cd(3) ; promptMap_CS_MC->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptMap_CS_MC->SetStats(0); promptMap_CS_MC->GetXaxis()->SetTitle("cos #theta_{CS}"); promptMap_CS_MC->GetYaxis()->SetTitle("#phi_{CS} [deg]"); promptMap_CS_MC->SetTitle(MCCSPRTitle); promptMap_CS_MC->Draw("colz");
    AcceptanceCanvas2->cd(4) ; promptMap_HX_MC->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptMap_HX_MC->SetStats(0); promptMap_HX_MC->GetXaxis()->SetTitle("cos #theta_{HX}"); promptMap_HX_MC->GetYaxis()->SetTitle("#phi_{HX} [deg]"); promptMap_HX_MC->SetTitle(MCHXPRTitle); promptMap_HX_MC->Draw("colz");

    sprintf(FilenameAcc,"%s/AccMaps_rapidity%d_pt%d.png",dirstruct,yBin+1,ptBin+1);

    AcceptanceCanvas2->SaveAs(FilenameAcc);

    AcceptanceCanvas2->Close();



    RooPlot* AccFramePR_costh_CS = new RooPlot;
    AccFramePR_costh_CS = costh_CS.frame() ;
	dataPRbin->plotOn(AccFramePR_costh_CS,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	cs->model()->plotOn(AccFramePR_costh_CS,LineWidth(2)/*,Normalization(1.0)*/);

	RooPlot* AccFramePR_phi_CS = new RooPlot;
    AccFramePR_phi_CS = phi_CS.frame() ;
	dataPRbin->plotOn(AccFramePR_phi_CS,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	cs->model()->plotOn(AccFramePR_phi_CS,LineWidth(2),Normalization(1.0));

    RooPlot* AccFramePR_costh_HX = new RooPlot;
    AccFramePR_costh_HX = costh_HX.frame() ;
 	dataPRbin->plotOn(AccFramePR_costh_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 	hx->model()->plotOn(AccFramePR_costh_HX,LineWidth(2),Normalization(1.0));

 	RooPlot* AccFramePR_phi_HX = new RooPlot;
    AccFramePR_phi_HX = phi_HX.frame() ;
 	dataPRbin->plotOn(AccFramePR_phi_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 	hx->model()->plotOn(AccFramePR_phi_HX,LineWidth(2),Normalization(1.0));

    TCanvas* Acc1DCanvas = new TCanvas("Acc1DCanvas","Acc1DCanvas",1600,1000);
    Acc1DCanvas->Divide(2,2);  Acc1DCanvas->SetFillColor(kWhite);
    Acc1DCanvas->cd(1) ; gPad->SetFillColor(kWhite); AccFramePR_costh_CS->Draw();
    Acc1DCanvas->cd(2) ; gPad->SetFillColor(kWhite); AccFramePR_phi_CS->Draw();
    Acc1DCanvas->cd(3) ; gPad->SetFillColor(kWhite); AccFramePR_costh_HX->Draw();
    Acc1DCanvas->cd(4) ; gPad->SetFillColor(kWhite); AccFramePR_phi_HX->Draw();

    sprintf(FilenameAcc,"%s/AccHistsCS1D_rapidity%d_pt%d.png",dirstruct,yBin+1,ptBin+1);

    Acc1DCanvas->SaveAs(FilenameAcc);
    Acc1DCanvas->Close();

    delete AcceptanceCanvas2;
    delete AcceptanceHistCanvas;
    delete Acc1DCanvas;
    delete AccFramePR_costh_CS;
    delete AccFramePR_phi_CS;
    delete AccFramePR_costh_HX;
    delete AccFramePR_phi_HX;

    accFile->Close();
    recoEffFile->Close();
    trigEffFile->Close();

    delete accFile;
    delete recoEffFile;
    delete trigEffFile;

    }
  }
  


//  delete fInput;

  return 0;
}
