#include <iostream>
//#include <stdio.h>
//#include <direct.h>
//#include <sys/stat.h>
//#include <string>

//#include <sstream>
//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"

// RooFit Includes
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooExponential.h"
#include "RooAbsPdf.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooGaussian.h"


//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"



int main(int argc, char**argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

  bool nofrac(false);
  bool nonorms(false);

  for( int i=0;i < argc; ++i ) {
    if(std::string(argv[i]).find("--noFrac") != std::string::npos) nofrac = true;
    if(std::string(argv[i]).find("--noNorms") != std::string::npos) nonorms = true;

  }

    cout<<"Options"<<endl;
	cout<<"        --noFrac       -> Fraction integrals are not calculated and plotted"<<endl;
	cout<<"        --noNorms      -> PR, NP and BK normalizations and B fractions are not plotted"<<endl;

   gSystem->mkdir("Plots/ParameterPlots");

   gStyle->SetPalette(1,0);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadRightMargin(0.15);

   gStyle->SetTickLength(-0.02, "xyz");
   gStyle->SetLabelOffset(0.02, "x");
   gStyle->SetLabelOffset(0.02, "y");
   gStyle->SetTitleOffset(1.3, "x");
   gStyle->SetTitleOffset(1.4, "y");
   gStyle->SetTitleFillColor(kWhite);

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,500);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
  varlist.add(JpsictErr);

	char outputfilename[200];
	sprintf(outputfilename,"Plots/ParameterPlots/PolarizationResults.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");

   Char_t *fileNameInputCS = "jPsiFitFinal_cs_GA07RE17.root";
   TFile* fInputCS = new TFile(fileNameInputCS);

   Char_t *fileNameInputHX = "jPsiFitFinal_hx_GA07RE17.root";
   TFile* fInputHX = new TFile(fileNameInputHX);

   TFile* PromptFile = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root");

   TTree *PromptTree = (TTree*)PromptFile->Get("data");
   RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",PromptTree,varlist);//,0,"MCweight");


   bool pereverr(false);

     for(int i=0;i<argc; ++i) {
       if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
     }

//// bin numbers to be plotted /////
   const int ptBinMin=2;
   const int ptBinMax=0;//number of pT Bins not to be plotted starting from the highest ptBin
   const int yBinMin=1;
   const int yBinMax=2;

   double offsetValue=-10000;

	RooRealVar *nNonPromptCS, *nPromptCS, *nBackground, *nNonPromptHX, *nPromptHX;
	RooRealVar *lambda_theta_CS_PR, *lambda_phi_CS_PR, *lambda_thetaphi_CS_PR, *lambda_theta_HX_PR, *lambda_phi_HX_PR, *lambda_thetaphi_HX_PR;

	double nNonPromptCS_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errnNonPromptCS_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nPromptCS_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errnPromptCS_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nBackground_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errnBackground_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double BfractionCS_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errBfractionCS_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nNonPromptHX_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errnNonPromptHX_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nPromptHX_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errnPromptHX_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double BfractionHX_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errBfractionHX_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double sigOvbkg_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errsigOvbkg_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double lambda_theta_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errlambda_theta_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double lambda_phi_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errlambda_phi_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double lambda_theta_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errlambda_theta_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double lambda_phi_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errlambda_phi_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double invariant_lambda_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errinvariant_lambda_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double invariant_lambda_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errinvariant_lambda_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double invariant_F_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errinvariant_F_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double invariant_F_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errinvariant_F_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double prompt_lifetime_cut[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nonprompt_lifetime_frac[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double prompt_lifetime_frac[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nonprompt_lifetime_cut[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double background_lifetime_frac[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double nonprompt_lifetime_absfrac[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double lambda_thetaphi_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errlambda_thetaphi_CS_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double lambda_thetaphi_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];
	double errlambda_thetaphi_HX_PR_[jpsi::kNbPTMaxBins][jpsi::kNbRapForPTBins];

	double ptMean[jpsi::kNbPTMaxBins]={-9999,-9999,-9999,-9999,-9999,-9999};
	double errptMean[jpsi::kNbPTMaxBins];
	double BfractionCS_rap[jpsi::kNbPTMaxBins];
	double errBfractionCS_rap[jpsi::kNbPTMaxBins];
	double BfractionHX_rap[jpsi::kNbPTMaxBins];
	double errBfractionHX_rap[jpsi::kNbPTMaxBins];
	double nBackground_rap[jpsi::kNbPTMaxBins];
	double errnBackground_rap[jpsi::kNbPTMaxBins];
	double sigOvbkg_rap[jpsi::kNbPTMaxBins];
	double errsigOvbkg_rap[jpsi::kNbPTMaxBins];
	double lambda_theta_CS_PR_rap[jpsi::kNbPTMaxBins];
	double errlambda_theta_CS_PR_rap[jpsi::kNbPTMaxBins];
	double lambda_phi_CS_PR_rap[jpsi::kNbPTMaxBins];
	double errlambda_phi_CS_PR_rap[jpsi::kNbPTMaxBins];
	double lambda_theta_HX_PR_rap[jpsi::kNbPTMaxBins];
	double errlambda_theta_HX_PR_rap[jpsi::kNbPTMaxBins];
	double lambda_phi_HX_PR_rap[jpsi::kNbPTMaxBins];
	double errlambda_phi_HX_PR_rap[jpsi::kNbPTMaxBins];
	double invariant_lambda_CS_PR_rap[jpsi::kNbPTMaxBins];
	double errinvariant_lambda_CS_PR_rap[jpsi::kNbPTMaxBins];
	double invariant_lambda_HX_PR_rap[jpsi::kNbPTMaxBins];
	double errinvariant_lambda_HX_PR_rap[jpsi::kNbPTMaxBins];
	double invariant_F_CS_PR_rap[jpsi::kNbPTMaxBins];
	double errinvariant_F_CS_PR_rap[jpsi::kNbPTMaxBins];
	double invariant_F_HX_PR_rap[jpsi::kNbPTMaxBins];
	double errinvariant_F_HX_PR_rap[jpsi::kNbPTMaxBins];
	double prompt_lifetime_cut_rap[jpsi::kNbPTMaxBins];
	double errprompt_lifetime_cut_rap[jpsi::kNbPTMaxBins];
	double nonprompt_lifetime_frac_rap[jpsi::kNbPTMaxBins];
	double errnonprompt_lifetime_frac_rap[jpsi::kNbPTMaxBins];
	double prompt_lifetime_frac_rap[jpsi::kNbPTMaxBins];
	double errprompt_lifetime_frac_rap[jpsi::kNbPTMaxBins];
	double nonprompt_lifetime_cut_rap[jpsi::kNbPTMaxBins];
	double errnonprompt_lifetime_cut_rap[jpsi::kNbPTMaxBins];
	double background_lifetime_frac_rap[jpsi::kNbPTMaxBins];
	double errbackground_lifetime_frac_rap[jpsi::kNbPTMaxBins];
	double nonprompt_lifetime_absfrac_rap[jpsi::kNbPTMaxBins];
	double errnonprompt_lifetime_absfrac_rap[jpsi::kNbPTMaxBins];
	double puffer_rap[jpsi::kNbPTMaxBins];
	double errpuffer_rap[jpsi::kNbPTMaxBins];
	double lambda_thetaphi_CS_PR_rap[jpsi::kNbPTMaxBins];
	double errlambda_thetaphi_CS_PR_rap[jpsi::kNbPTMaxBins];
	double lambda_thetaphi_HX_PR_rap[jpsi::kNbPTMaxBins];
	double errlambda_thetaphi_HX_PR_rap[jpsi::kNbPTMaxBins];

	Double_t XS_Bfrac_mid[3] = {0.162, 0.257, 0.369};
	Double_t errXS_Bfrac_mid[3] = {0.071, 0.033, 0.041};
	Double_t XS_pTMean_mid[3] = {5.59, 7.88, 13.53};
	Double_t XS_Bfrac_fwd[5] = {0.098, 0.112, 0.165, 0.203, 0.331};
	Double_t errXS_Bfrac_fwd[5] = {0.058, 0.024, 0.029, 0.029, 0.057};
	Double_t XS_pTMean_fwd[5] = {1.23, 2.86, 4.89, 7.59, 13.14};

	int nonConvergeptBinRap1CS1 = 2;int nonConvergeptBinRap1CS2 = 5;int nonConvergeptBinRap1CS3 = 6;int nonConvergeptBinRap1CS4 = 999;

	int nonConvergeptBinRap2CS1 = 3;int nonConvergeptBinRap2CS2 = 4;int nonConvergeptBinRap2CS3 = 5;int nonConvergeptBinRap2CS4 = 6;

	int nonConvergeptBinRap3CS1 = 999;int nonConvergeptBinRap3CS2 = 999;int nonConvergeptBinRap3CS3 = 999;int nonConvergeptBinRap3CS4 = 999;

	int nonConvergeptBinRap4CS1 = 999;int nonConvergeptBinRap4CS2 = 999;int nonConvergeptBinRap4CS3 = 999;int nonConvergeptBinRap4CS4 = 999;

	int nonConvergeptBinRap5CS1 = 999;int nonConvergeptBinRap5CS2 = 999;int nonConvergeptBinRap5CS3 = 999;int nonConvergeptBinRap5CS4 = 999;


	int nonConvergeptBinRap1HX1 = 3;int nonConvergeptBinRap1HX2 = 7;int nonConvergeptBinRap1HX3 = 999;int nonConvergeptBinRap1HX4 = 999;

	int nonConvergeptBinRap2HX1 = 2;int nonConvergeptBinRap2HX2 = 7;int nonConvergeptBinRap2HX3 = 999;int nonConvergeptBinRap2HX4 = 999;

	int nonConvergeptBinRap3HX1 = 999;int nonConvergeptBinRap3HX2 = 999;int nonConvergeptBinRap3HX3 = 999;int nonConvergeptBinRap3HX4 = 999;

	int nonConvergeptBinRap4HX1 = 999;int nonConvergeptBinRap4HX2 = 999;int nonConvergeptBinRap4HX3 = 999;int nonConvergeptBinRap4HX4 = 999;

	int nonConvergeptBinRap5HX1 = 999;int nonConvergeptBinRap5HX2 = 999;int nonConvergeptBinRap5HX3 = 999;int nonConvergeptBinRap5HX4 = 999;


	double ptplotMin = 0;



	  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
		  for(int ptBin = 0; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {

      
      
      
      
      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char DirectoryPath[200];
	  sprintf(DirectoryPath,"/pt%d_rapidity%d",ptBin+1,yBin+1);
	  char DirectoryPathPol[200];
	  sprintf(DirectoryPathPol,"/pt%d_rapidity%d/PromptPolarizationModel",ptBin+1,yBin+1);






	  if(fInputCS->GetDirectory(DirectoryPath)!=NULL){
		  TDirectory *InputDirectoryCS = (TDirectory*)fInputCS->GetDirectory(DirectoryPath);

		  if(InputDirectoryCS->Get("nPrompt")!=NULL){
	  nPromptCS = (RooRealVar*)InputDirectoryCS->Get("nPrompt");
	  nPromptCS_[ptBin][yBin]=nPromptCS->getVal();
	  errnPromptCS_[ptBin][yBin]=nPromptCS->getError();
		  }

		  if(InputDirectoryCS->Get("nNonPrompt")!=NULL){
	  nNonPromptCS = (RooRealVar*)InputDirectoryCS->Get("nNonPrompt");
	  nNonPromptCS_[ptBin][yBin]=nNonPromptCS->getVal();
	  errnNonPromptCS_[ptBin][yBin]=nNonPromptCS->getError();
		  }

	  BfractionCS_[ptBin][yBin]=nNonPromptCS_[ptBin][yBin]/(nNonPromptCS_[ptBin][yBin]+nPromptCS_[ptBin][yBin]);
	  errBfractionCS_[ptBin][yBin]=(1/(nPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin])+nNonPromptCS_[ptBin][yBin]/pow(nPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin],2))*errnNonPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin]/pow(nPromptCS_[ptBin][yBin]+nNonPromptCS_[ptBin][yBin],2)*errnPromptCS_[ptBin][yBin];

	  sigOvbkg_[ptBin][yBin]=(nNonPromptCS_[ptBin][yBin]+nPromptCS_[ptBin][yBin])/nBackground_[ptBin][yBin];
	  errsigOvbkg_[ptBin][yBin] = sigOvbkg_[ptBin][yBin] * sqrt(pow((errnNonPromptCS_[ptBin][yBin]+errnPromptCS_[ptBin][yBin])/(nNonPromptCS_[ptBin][yBin]+nPromptCS_[ptBin][yBin]),2) + pow(errnBackground_[ptBin][yBin]/nBackground_[ptBin][yBin],2));

	  if(InputDirectoryCS->Get("nBackground")!=NULL){
	  nBackground = (RooRealVar*)InputDirectoryCS->Get("nBackground");
	  nBackground_[ptBin][yBin]=nBackground->getVal();
	  errnBackground_[ptBin][yBin]=nBackground->getError();
	  }


	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);

	  RooDataSet *dataPRbin;

	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);


	  if(!nofrac){

	  CompositeModelBuilder *PromptModel = new CompositeModelBuilder("CS");
	  PromptModel->setUseLifetime(true);
	  PromptModel->setUseMass(false);
	  PromptModel->setUsePol(false);

	  PromptModel->loadParameters(*InputDirectoryCS);

		  if(pereverr)
			  PromptModel->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
		 else
			  PromptModel->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

			  RooAbsReal* integralPrompt;
			  double RangeMax = -1;
			  double DeltaRangeMax = 0.1;
			  double PromptCondition = 0;

			  for(int runDex = 0; runDex < 10000; runDex++) {

	  RangeMax=RangeMax+DeltaRangeMax;
	  Jpsict.setRange("PromptRange", -100000, RangeMax);
	  if(pereverr){
      integralPrompt =  PromptModel->getPromptModel()->createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("PromptRange"));
	  }
      else
	  integralPrompt =  PromptModel->getPromptModel()->createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("PromptRange"));//,NormSet(Jpsict)
//	  cout<<"int"<<integralPrompt->getVal()<<" upperborder"<<RangeMax<<endl;
//	  cout<<integralPrompt->getVal()<<endl;
	  if(DeltaRangeMax < 0.001 & integralPrompt->getVal() > PromptCondition) break;

	  if(integralPrompt->getVal() > PromptCondition) {
		  RangeMax = RangeMax - DeltaRangeMax;
		  DeltaRangeMax = DeltaRangeMax/10;
	  }


				}


	 prompt_lifetime_cut[ptBin][yBin] = RangeMax;

	 Jpsict.setRange("NonPromptRange", RangeMax, 2.5);

	 nonprompt_lifetime_frac[ptBin][yBin] = PromptModel->getNonPromptModel()->createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("NonPromptRange"))->getVal();
	 background_lifetime_frac[ptBin][yBin] = PromptModel->getBkgModel()->createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("NonPromptRange"))->getVal();

	 cout<<"integral="<<integralPrompt->getVal()<<" UpperBoarder="<<RangeMax<<", Fraction of Non Prompt accepted by cut = "<<nonprompt_lifetime_frac[ptBin][yBin]<<", Fraction of Background accepted by cut = "<<background_lifetime_frac[ptBin][yBin]<<endl;;

	 nonprompt_lifetime_absfrac[ptBin][yBin]=nonprompt_lifetime_frac[ptBin][yBin]*nNonPromptCS_[ptBin][yBin]/(nonprompt_lifetime_frac[ptBin][yBin]*nNonPromptCS_[ptBin][yBin]+(1-PromptCondition)*nPromptCS_[ptBin][yBin]+nBackground_[ptBin][yBin]*background_lifetime_frac[ptBin][yBin]);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	 RooAbsReal* integralNonPrompt;
	 RangeMax = 2.5;
	 DeltaRangeMax = 0.1;

	  for(int runDex = 0; runDex < 10000; runDex++) {

	  RangeMax=RangeMax-DeltaRangeMax;
	  Jpsict.setRange("NonPromptRangeCut", -100, RangeMax);

	  integralNonPrompt =  PromptModel->getNonPromptModel()->createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("NonPromptRangeCut"));//

//	  cout<<"int"<<integralNonPrompt->getVal()<<"delta"<<DeltaRangeMax<<endl;

	  if(DeltaRangeMax < 0.000005 & integralNonPrompt->getVal() < 0.01) break;

	  if(integralNonPrompt->getVal() < 0){//0.01) {
		  RangeMax = RangeMax + DeltaRangeMax;
		  DeltaRangeMax = DeltaRangeMax/10;
	  }


				}

	  nonprompt_lifetime_cut[ptBin][yBin] = RangeMax;

	 	 Jpsict.setRange("PromptRangeFrac", -100000, RangeMax);

	 	 prompt_lifetime_frac[ptBin][yBin] = PromptModel->getPromptModel()->createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("PromptRangeFrac"))->getVal();

	 	 cout<<"integral="<<integralNonPrompt->getVal()<<" UpperBoarder="<<RangeMax<<", Fraction of Prompt accepted by cut = "<<prompt_lifetime_frac[ptBin][yBin]<<endl;


	  delete PromptModel;


	  }
	  }
	  else cout<<"TDirectory "<<DirectoryPath<<" does not exist in File"<<fileNameInputCS<<endl;


	  if(fInputHX->GetDirectory(DirectoryPath)!=NULL){
		  TDirectory *InputDirectoryHX = (TDirectory*)fInputHX->GetDirectory(DirectoryPath);

		  if(InputDirectoryHX->Get("nPrompt")!=NULL){
	  nPromptHX = (RooRealVar*)InputDirectoryHX->Get("nPrompt");
	  nPromptHX_[ptBin][yBin]=nPromptHX->getVal();
	  errnPromptHX_[ptBin][yBin]=nPromptHX->getError();
		  }

		  if(InputDirectoryHX->Get("nPrompt")!=NULL){
	  nNonPromptHX = (RooRealVar*)InputDirectoryHX->Get("nNonPrompt");
	  nNonPromptHX_[ptBin][yBin]=nNonPromptHX->getVal();
	  errnNonPromptHX_[ptBin][yBin]=nNonPromptHX->getError();
		  }

	  BfractionHX_[ptBin][yBin]=nNonPromptHX_[ptBin][yBin]/(nNonPromptHX_[ptBin][yBin]+nPromptHX_[ptBin][yBin]);
	  errBfractionHX_[ptBin][yBin]=(1/(nPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin])+nNonPromptHX_[ptBin][yBin]/pow(nPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin],2))*errnNonPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin]/pow(nPromptHX_[ptBin][yBin]+nNonPromptHX_[ptBin][yBin],2)*errnPromptHX_[ptBin][yBin];
	  }
	  else cout<<"TDirectory "<<DirectoryPath<<" does not exist in File "<<fileNameInputHX<<endl;



	  if(fInputCS->GetDirectory(DirectoryPathPol)!=NULL){
		  TDirectory *InputDirectoryPolCS = (TDirectory*)fInputCS->GetDirectory(DirectoryPathPol);

	  lambda_theta_CS_PR = (RooRealVar*)InputDirectoryPolCS->Get("promptlambda_theta_CS");
	  lambda_theta_CS_PR_[ptBin][yBin]=lambda_theta_CS_PR->getVal();
	  errlambda_theta_CS_PR_[ptBin][yBin]=lambda_theta_CS_PR->getError();

	  lambda_phi_CS_PR = (RooRealVar*)InputDirectoryPolCS->Get("promptlambda_phi_CS");
	  lambda_phi_CS_PR_[ptBin][yBin]=lambda_phi_CS_PR->getVal();
	  errlambda_phi_CS_PR_[ptBin][yBin]=lambda_phi_CS_PR->getError();

	  lambda_thetaphi_CS_PR = (RooRealVar*)InputDirectoryPolCS->Get("promptlambda_thetaphi_CS");
	  lambda_thetaphi_CS_PR_[ptBin][yBin]=lambda_thetaphi_CS_PR->getVal();
	  errlambda_thetaphi_CS_PR_[ptBin][yBin]=lambda_thetaphi_CS_PR->getError();

	  invariant_lambda_CS_PR_[ptBin][yBin]=(lambda_theta_CS_PR_[ptBin][yBin]+3*lambda_phi_CS_PR_[ptBin][yBin])/(1-lambda_phi_CS_PR_[ptBin][yBin]);
	  invariant_F_CS_PR_[ptBin][yBin]=(1+lambda_theta_CS_PR_[ptBin][yBin]+2*lambda_phi_CS_PR_[ptBin][yBin])/(3-lambda_theta_CS_PR_[ptBin][yBin]);
	  }
	  else cout<<"TDirectory "<<DirectoryPathPol<<" does not exist in File "<<fileNameInputCS<<endl;

	  if(fInputHX->GetDirectory(DirectoryPathPol)!=NULL){
		  TDirectory *InputDirectoryPolHX = (TDirectory*)fInputHX->GetDirectory(DirectoryPathPol);

	  lambda_theta_HX_PR = (RooRealVar*)InputDirectoryPolHX->Get("promptlambda_theta_HX");
	  lambda_theta_HX_PR_[ptBin][yBin]=lambda_theta_HX_PR->getVal();
	  errlambda_theta_HX_PR_[ptBin][yBin]=lambda_theta_HX_PR->getError();

	  lambda_phi_HX_PR = (RooRealVar*)InputDirectoryPolHX->Get("promptlambda_phi_HX");
	  lambda_phi_HX_PR_[ptBin][yBin]=lambda_phi_HX_PR->getVal();
	  errlambda_phi_HX_PR_[ptBin][yBin]=lambda_phi_HX_PR->getError();

	  lambda_thetaphi_HX_PR = (RooRealVar*)InputDirectoryPolHX->Get("promptlambda_thetaphi_HX");
	  lambda_thetaphi_HX_PR_[ptBin][yBin]=lambda_thetaphi_HX_PR->getVal();
	  errlambda_thetaphi_HX_PR_[ptBin][yBin]=lambda_thetaphi_HX_PR->getError();

	  invariant_lambda_HX_PR_[ptBin][yBin]=(lambda_theta_HX_PR_[ptBin][yBin]+3*lambda_phi_HX_PR_[ptBin][yBin])/(1-lambda_phi_HX_PR_[ptBin][yBin]);
	  invariant_F_HX_PR_[ptBin][yBin]=(1+lambda_theta_HX_PR_[ptBin][yBin]+2*lambda_phi_HX_PR_[ptBin][yBin])/(3-lambda_theta_HX_PR_[ptBin][yBin]);
	  }
	  else cout<<"TDirectory "<<DirectoryPathPol<<" does not exist in File "<<fileNameInputHX<<endl;

    }



	}

////////////////////////////////// 	B - FRACTION CS ///////////////////////////////////////////////////////////////////////////////////



	TLegend* bfracLegendCS=new TLegend(0.15,0.7,0.65,0.89);
	bfracLegendCS->SetFillColor(kWhite);
	bfracLegendCS->SetTextFont(72);
	bfracLegendCS->SetTextSize(0.02);
	bfracLegendCS->SetBorderSize(0);

	char legendrap[200];

	TCanvas *BfractionCanvasCS = new TCanvas("fraction of J/#psi from B hadrons","fraction of J/#psi from B hadrons",1000,700);
	BfractionCanvasCS->SetFillColor(kWhite);
//	BfractionCanvasCS->SetGrid();
	BfractionCanvasCS->GetFrame()->SetFillColor(kWhite);
	BfractionCanvasCS->GetFrame()->SetBorderSize(0);
//	BfractionCanvasCS->SetLeftMargin(0.15) ;

	TH1F *BbfracHistoCS = BfractionCanvasCS->DrawFrame(0,0,20,1);
	BbfracHistoCS->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	BbfracHistoCS->SetYTitle("fraction of J/#psi from B hadrons");
	BbfracHistoCS->GetYaxis()->SetTitleOffset(1.5);


    for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
    	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	  					   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
	  					   errptMean[ptBin]=0;
	  					 BfractionCS_rap[ptBin]=BfractionCS_[ptBin][yBin];
	  					if (ptBin==0 && yBin==0 | yBin==1) BfractionCS_rap[ptBin] = -1;
	  					if (ptBin==1 && yBin==0) BfractionCS_rap[ptBin] = -1;
	  					 errBfractionCS_rap[ptBin]=errBfractionCS_[ptBin][yBin];

    	}

	TGraphErrors *bfracGraphCS = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,BfractionCS_rap,errptMean,errBfractionCS_rap);
	bfracGraphCS->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	bfracGraphCS->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphCS->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	sprintf(legendrap,"Rapidity %d",yBin+1);
	bfracLegendCS->AddEntry(bfracGraphCS,legendrap,"ple");
	bfracGraphCS->Draw("P");

	  		}

	TGraphErrors *bfracGraphCS = new TGraphErrors(3,XS_pTMean_mid,XS_Bfrac_mid,errptMean,errXS_Bfrac_mid);
	bfracGraphCS->SetMarkerColor(kBlack);
//	bfracGraphCS->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphCS->SetMarkerStyle(13);
	sprintf(legendrap,"XS cross check, |y|<1.4 (100nb^{-1})");
	bfracLegendCS->AddEntry(bfracGraphCS,legendrap,"ple");
	bfracGraphCS->Draw("P");

	TGraphErrors *bfracGraphCS2 = new TGraphErrors(5,XS_pTMean_fwd,XS_Bfrac_fwd,errptMean,errXS_Bfrac_fwd);
	bfracGraphCS2->SetMarkerColor(kBlue);
//	bfracGraphCS2->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphCS2->SetMarkerStyle(24);
	sprintf(legendrap,"XS cross check, 1.4<|y|<2.4 (100nb^{-1})");
	bfracLegendCS->AddEntry(bfracGraphCS2,legendrap,"ple");
	bfracGraphCS2->Draw("P");

	  		bfracLegendCS->Draw();

	  		  char Filename[200];
	  		  sprintf(Filename,"Plots/ParameterPlots/pT_BfractionCS.png");

	  		if(!nonorms){


	  		BfractionCanvasCS->SaveAs(Filename);

	  		BfractionCanvasCS->Close();

////////////////////////////////// 	B - FRACTION HX ///////////////////////////////////////////////////////////////////////////////////


	TLegend* bfracLegendHX=new TLegend(0.15,0.7,0.65,0.89);
	bfracLegendHX->SetFillColor(kWhite);
	bfracLegendHX->SetTextFont(72);
	bfracLegendHX->SetTextSize(0.02);
	bfracLegendHX->SetBorderSize(0);


	TCanvas *BfractionCanvasHX = new TCanvas("fraction of J/#psi from B hadrons","fraction of J/#psi from B hadrons",1000,700);
	BfractionCanvasHX->SetFillColor(kWhite);
//	BfractionCanvasHX->SetGrid();
	BfractionCanvasHX->GetFrame()->SetFillColor(kWhite);
	BfractionCanvasHX->GetFrame()->SetBorderSize(0);
//	BfractionCanvasHX->SetLeftMargin(0.15) ;

	TH1F *BbfracHistoHX = BfractionCanvasHX->DrawFrame(0,0,20,1);
	BbfracHistoHX->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	BbfracHistoHX->SetYTitle("fraction of J/#psi from B hadrons");
	BbfracHistoHX->GetYaxis()->SetTitleOffset(1.5);

	for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	    	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

		    				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
	  					   errptMean[ptBin]=0;
	  					 BfractionHX_rap[ptBin]=BfractionHX_[ptBin][yBin];
	  					if (ptBin==0 && yBin==0 | yBin==1) BfractionHX_rap[ptBin] = -1;
	  					if (ptBin==1 && yBin==0) BfractionHX_rap[ptBin] = -1;

	  					errBfractionHX_rap[ptBin]=errBfractionHX_[ptBin][yBin];
	  					  }

	TGraphErrors *bfracGraphHX = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,BfractionHX_rap,errptMean,errBfractionHX_rap);
	bfracGraphHX->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	bfracGraphHX->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphHX->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	sprintf(legendrap,"Rapidity %d",yBin+1);
	bfracLegendHX->AddEntry(bfracGraphHX,legendrap,"ple");
	bfracGraphHX->Draw("P");

	  		}



	TGraphErrors *bfracGraphHX3 = new TGraphErrors(3,XS_pTMean_mid,XS_Bfrac_mid,errptMean,errXS_Bfrac_mid);
	bfracGraphHX3->SetMarkerColor(kBlack);
//	bfracGraphHX3->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphHX3->SetMarkerStyle(13);
	sprintf(legendrap,"XS cross check, |y|<1.4 (100nb^{-1})");
	bfracLegendHX->AddEntry(bfracGraphHX3,legendrap,"ple");
	bfracGraphHX3->Draw("P");

	TGraphErrors *bfracGraphHX4 = new TGraphErrors(5,XS_pTMean_fwd,XS_Bfrac_fwd,errptMean,errXS_Bfrac_fwd);
	bfracGraphHX4->SetMarkerColor(kBlue);
//	bfracGraphHX4->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	bfracGraphHX4->SetMarkerStyle(24);
	sprintf(legendrap,"XS cross check, 1.4<|y|<2.4 (100nb^{-1})");
	bfracLegendHX->AddEntry(bfracGraphHX4,legendrap,"ple");
	bfracGraphHX4->Draw("P");



	  		bfracLegendHX->Draw();

	  		  sprintf(Filename,"Plots/ParameterPlots/pT_BfractionHX.png");

	  		BfractionCanvasHX->SaveAs(Filename);

	  		BfractionCanvasHX->Close();

////////////////////////////////// 	Number of Background Events ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* BackgroundLegend=new TLegend(0.875,0.6,1,0.89);
	 BackgroundLegend->SetFillColor(kWhite);
	 BackgroundLegend->SetTextFont(72);
	 BackgroundLegend->SetTextSize(0.02);
	 BackgroundLegend->SetBorderSize(0);

	 TCanvas *BackgroundCanvas = new TCanvas("number of Background events","",1000,700);
	 BackgroundCanvas->SetFillColor(kWhite);
//   BackgroundCanvas->SetGrid();
	 BackgroundCanvas->GetFrame()->SetFillColor(kWhite);
	 BackgroundCanvas->GetFrame()->SetBorderSize(0);
	 BackgroundCanvas->SetLeftMargin(0.15) ;

	 TH1F *BackgroundhistoHisto = BackgroundCanvas->DrawFrame(0,0,30,65000);
	 BackgroundhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 BackgroundhistoHisto->SetYTitle("number of Background events");
	 BackgroundhistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	    			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
	 	  					nBackground_rap[ptBin]=nBackground_[ptBin][yBin];
	 	  					errnBackground_rap[ptBin]=errnBackground_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *BackgroundGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,nBackground_rap,errptMean,errnBackground_rap);
	 BackgroundGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	 BackgroundGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 BackgroundGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 BackgroundLegend->AddEntry(BackgroundGraph,legendrap,"ple");
	 BackgroundGraph->Draw("P");

	 	  		}

	 	  		BackgroundLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_nBackground.png");

	 	  		BackgroundCanvas->SaveAs(Filename);

	 	  		BackgroundCanvas->Close();

////////////////////////////////// 	Number of PromptCS Events ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* PromptCSLegend=new TLegend(0.875,0.6,1,0.89);
	 PromptCSLegend->SetFillColor(kWhite);
	 PromptCSLegend->SetTextFont(72);
	 PromptCSLegend->SetTextSize(0.02);
	 PromptCSLegend->SetBorderSize(0);

	 TCanvas *PromptCSCanvas = new TCanvas("number of PromptCS events","",1000,700);
	 PromptCSCanvas->SetFillColor(kWhite);
//   PromptCSCanvas->SetGrid();
	 PromptCSCanvas->GetFrame()->SetFillColor(kWhite);
	 PromptCSCanvas->GetFrame()->SetBorderSize(0);
	 PromptCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptCShistoHisto = PromptCSCanvas->DrawFrame(0,0,30,65000);
	 PromptCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptCShistoHisto->SetYTitle("number of PromptCS events");
	 PromptCShistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	    			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
	 	  					puffer_rap[ptBin]=nPromptCS_[ptBin][yBin];
	 	  					errpuffer_rap[ptBin]=errnPromptCS_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *PromptCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,puffer_rap,errptMean,errpuffer_rap);
	 PromptCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	 PromptCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 PromptCSLegend->AddEntry(PromptCSGraph,legendrap,"ple");
	 PromptCSGraph->Draw("P");

	 	  		}

	 	  		PromptCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_nPromptCS.png");

	 	  		PromptCSCanvas->SaveAs(Filename);

	 	  		PromptCSCanvas->Close();


////////////////////////////////// 	Number of NonPromptCS Events ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* NonPromptCSLegend=new TLegend(0.875,0.6,1,0.89);
	 NonPromptCSLegend->SetFillColor(kWhite);
	 NonPromptCSLegend->SetTextFont(72);
	 NonPromptCSLegend->SetTextSize(0.02);
	 NonPromptCSLegend->SetBorderSize(0);

	 TCanvas *NonPromptCSCanvas = new TCanvas("number of NonPromptCS events","",1000,700);
	 NonPromptCSCanvas->SetFillColor(kWhite);
//   NonPromptCSCanvas->SetGrid();
	 NonPromptCSCanvas->GetFrame()->SetFillColor(kWhite);
	 NonPromptCSCanvas->GetFrame()->SetBorderSize(0);
	 NonPromptCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *NonPromptCShistoHisto = NonPromptCSCanvas->DrawFrame(0,0,30,65000);
	 NonPromptCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 NonPromptCShistoHisto->SetYTitle("number of NonPromptCS events");
	 NonPromptCShistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	    			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
	 	  					puffer_rap[ptBin]=nNonPromptCS_[ptBin][yBin];
	 	  					errpuffer_rap[ptBin]=errnNonPromptCS_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *NonPromptCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,puffer_rap,errptMean,errpuffer_rap);
	 NonPromptCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	 NonPromptCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 NonPromptCSLegend->AddEntry(NonPromptCSGraph,legendrap,"ple");
	 NonPromptCSGraph->Draw("P");

	 	  		}

	 	  		NonPromptCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_nNonPromptCS.png");

	 	  		NonPromptCSCanvas->SaveAs(Filename);

	 	  		NonPromptCSCanvas->Close();



////////////////////////////////// 	SIG/BKG RATIO ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* sigOvbkgLegend=new TLegend(0.875,0.6,1,0.89);
	 sigOvbkgLegend->SetFillColor(kWhite);
	 sigOvbkgLegend->SetTextFont(72);
	 sigOvbkgLegend->SetTextSize(0.02);
	 sigOvbkgLegend->SetBorderSize(0);

	 TCanvas *sigOvbkgCanvas = new TCanvas("Sig/Bkg ratio","",1000,700);
	 sigOvbkgCanvas->SetFillColor(kWhite);
//   sigOvbkgCanvas->SetGrid();
	 sigOvbkgCanvas->GetFrame()->SetFillColor(kWhite);
	 sigOvbkgCanvas->GetFrame()->SetBorderSize(0);
//	 sigOvbkgCanvas->SetLeftMargin(0.15) ;

	 TH1F *sigOvbkghistoHisto = sigOvbkgCanvas->DrawFrame(0,0,30,20);
	 sigOvbkghistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 sigOvbkghistoHisto->SetYTitle("Sig/Bkg ratio");
	 sigOvbkghistoHisto->GetYaxis()->SetTitleOffset(1.5);

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	     			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
	 	  					sigOvbkg_rap[ptBin]=sigOvbkg_[ptBin][yBin];
	 	  					errsigOvbkg_rap[ptBin]=errsigOvbkg_[ptBin][yBin];
	 	  					  }

	 TGraphErrors *sigOvbkgGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,sigOvbkg_rap,errptMean,errsigOvbkg_rap);
	 sigOvbkgGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	 sigOvbkgGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 sigOvbkgGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 sigOvbkgLegend->AddEntry(sigOvbkgGraph,legendrap,"ple");
	 sigOvbkgGraph->Draw("P");

	 	  		}

	 	  		sigOvbkgLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_sigOvbkg.png");

	 	  		sigOvbkgCanvas->SaveAs(Filename);

	 	  		sigOvbkgCanvas->Close();
	 }
////////////////////////////////// LAMBDA THETA PROMPT CS ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaThetaPRCSLegend=new TLegend(0.875,0.6,1,0.89);
	 LambdaThetaPRCSLegend->SetFillColor(kWhite);
	 LambdaThetaPRCSLegend->SetTextFont(72);
	 LambdaThetaPRCSLegend->SetTextSize(0.02);
	 LambdaThetaPRCSLegend->SetBorderSize(0);

	 TCanvas *LambdaThetaPRCSCanvas = new TCanvas("Prompt #lambda_{#theta_{CS}}","",1000,700);
	 LambdaThetaPRCSCanvas->SetFillColor(kWhite);
     LambdaThetaPRCSCanvas->SetGrid(0,1);
	 LambdaThetaPRCSCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaThetaPRCSCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaThetaPRCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaThetaPRCShistoHisto = LambdaThetaPRCSCanvas->DrawFrame(4,-1.3,30,1.3);
	 LambdaThetaPRCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaThetaPRCShistoHisto->SetYTitle("Prompt #lambda_{#theta_{CS}}");
	 LambdaThetaPRCShistoHisto->GetYaxis()->SetTitleOffset(1.5);
		fprintf(outputFile, "\nLambda_Theta_Prompt_CS\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_theta_CS_PR_rap[ptBin]=lambda_theta_CS_PR_[ptBin][yBin];
	 	  					errlambda_theta_CS_PR_rap[ptBin]=errlambda_theta_CS_PR_[ptBin][yBin];
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.6f +- %1.6f\n",ptBin+1,yBin+1,lambda_theta_CS_PR_rap[ptBin],errlambda_theta_CS_PR_rap[ptBin]);

	 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) lambda_theta_CS_PR_rap[ptBin] = offsetValue;

		  					if (ptBin+1==nonConvergeptBinRap1CS1 |ptBin+1==nonConvergeptBinRap1CS2 | ptBin+1==nonConvergeptBinRap1CS3 | ptBin+1==nonConvergeptBinRap1CS4 && yBin==0) lambda_theta_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap2CS1 |ptBin+1==nonConvergeptBinRap2CS2 | ptBin+1==nonConvergeptBinRap2CS3 | ptBin+1==nonConvergeptBinRap2CS4 && yBin==1) lambda_theta_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap3CS1 |ptBin+1==nonConvergeptBinRap3CS2 | ptBin+1==nonConvergeptBinRap3CS3 | ptBin+1==nonConvergeptBinRap3CS4 && yBin==2) lambda_theta_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap4CS1 |ptBin+1==nonConvergeptBinRap4CS2 | ptBin+1==nonConvergeptBinRap4CS3 | ptBin+1==nonConvergeptBinRap4CS4 && yBin==3) lambda_theta_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap5CS1 |ptBin+1==nonConvergeptBinRap5CS2 | ptBin+1==nonConvergeptBinRap5CS3 | ptBin+1==nonConvergeptBinRap5CS4 && yBin==4) lambda_theta_CS_PR_rap[ptBin] = offsetValue;


	     	}

	 TGraphErrors *LambdaThetaPRCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,lambda_theta_CS_PR_rap,errptMean,errlambda_theta_CS_PR_rap);
	 LambdaThetaPRCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaThetaPRCSLegend->AddEntry(LambdaThetaPRCSGraph,legendrap,"ple");
	 LambdaThetaPRCSGraph->Draw("P");

	 	  		}

	 	  		LambdaThetaPRCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_LambdaThetaPRCS.png");

	 	  		LambdaThetaPRCSCanvas->SaveAs(Filename);

	 	  		LambdaThetaPRCSCanvas->Close();

////////////////////////////////// LAMBDA PHI PROMPT CS ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaPhiPRCSLegend=new TLegend(0.875,0.6,1,0.89);
	 LambdaPhiPRCSLegend->SetFillColor(kWhite);
	 LambdaPhiPRCSLegend->SetTextFont(72);
	 LambdaPhiPRCSLegend->SetTextSize(0.02);
	 LambdaPhiPRCSLegend->SetBorderSize(0);

	 TCanvas *LambdaPhiPRCSCanvas = new TCanvas("Prompt #lambda_{#phi_{CS}}","",1000,700);
	 LambdaPhiPRCSCanvas->SetFillColor(kWhite);
     LambdaPhiPRCSCanvas->SetGrid(0,1);
	 LambdaPhiPRCSCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaPhiPRCSCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaPhiPRCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaPhiPRCShistoHisto = LambdaPhiPRCSCanvas->DrawFrame(4,-1.3,30,1.3);
	 LambdaPhiPRCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaPhiPRCShistoHisto->SetYTitle("Prompt #lambda_{#phi_{CS}}");
	 LambdaPhiPRCShistoHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nLambda_Phi_Prompt_CS\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	     			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_phi_CS_PR_rap[ptBin]=lambda_phi_CS_PR_[ptBin][yBin];
	 	  					errlambda_phi_CS_PR_rap[ptBin]=errlambda_phi_CS_PR_[ptBin][yBin];
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.6f +- %1.6f\n",ptBin+1,yBin+1,lambda_phi_CS_PR_rap[ptBin],errlambda_phi_CS_PR_rap[ptBin]);

		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) lambda_phi_CS_PR_rap[ptBin] = offsetValue;


			  					if (ptBin+1==nonConvergeptBinRap1CS1 |ptBin+1==nonConvergeptBinRap1CS2 | ptBin+1==nonConvergeptBinRap1CS3 | ptBin+1==nonConvergeptBinRap1CS4 && yBin==0) lambda_phi_CS_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap2CS1 |ptBin+1==nonConvergeptBinRap2CS2 | ptBin+1==nonConvergeptBinRap2CS3 | ptBin+1==nonConvergeptBinRap2CS4 && yBin==1) lambda_phi_CS_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap3CS1 |ptBin+1==nonConvergeptBinRap3CS2 | ptBin+1==nonConvergeptBinRap3CS3 | ptBin+1==nonConvergeptBinRap3CS4 && yBin==2) lambda_phi_CS_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap4CS1 |ptBin+1==nonConvergeptBinRap4CS2 | ptBin+1==nonConvergeptBinRap4CS3 | ptBin+1==nonConvergeptBinRap4CS4 && yBin==3) lambda_phi_CS_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap5CS1 |ptBin+1==nonConvergeptBinRap5CS2 | ptBin+1==nonConvergeptBinRap5CS3 | ptBin+1==nonConvergeptBinRap5CS4 && yBin==4) lambda_phi_CS_PR_rap[ptBin] = offsetValue;


	     	}

	 TGraphErrors *LambdaPhiPRCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,lambda_phi_CS_PR_rap,errptMean,errlambda_phi_CS_PR_rap);
	 LambdaPhiPRCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaPhiPRCSLegend->AddEntry(LambdaPhiPRCSGraph,legendrap,"ple");
	 LambdaPhiPRCSGraph->Draw("P");

	 	  		}

	 	  		LambdaPhiPRCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_LambdaPhiPRCS.png");

	 	  		LambdaPhiPRCSCanvas->SaveAs(Filename);

	 	  		LambdaPhiPRCSCanvas->Close();


////////////////////////////////// LAMBDA THETA PROMPT CS ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaThetaPhiPRCSLegend=new TLegend(0.875,0.6,1,0.89);
	 LambdaThetaPhiPRCSLegend->SetFillColor(kWhite);
	 LambdaThetaPhiPRCSLegend->SetTextFont(72);
	 LambdaThetaPhiPRCSLegend->SetTextSize(0.02);
	 LambdaThetaPhiPRCSLegend->SetBorderSize(0);

	 TCanvas *LambdaThetaPhiPRCSCanvas = new TCanvas("Prompt #lambda_{#theta #phi_{CS}}","",1000,700);
	 LambdaThetaPhiPRCSCanvas->SetFillColor(kWhite);
     LambdaThetaPhiPRCSCanvas->SetGrid(0,1);
	 LambdaThetaPhiPRCSCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaThetaPhiPRCSCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaThetaPhiPRCSCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaThetaPhiPRCShistoHisto = LambdaThetaPhiPRCSCanvas->DrawFrame(4,-1.3,30,1.3);
	 LambdaThetaPhiPRCShistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaThetaPhiPRCShistoHisto->SetYTitle("Prompt #lambda_{#theta #phi_{CS}}");
	 LambdaThetaPhiPRCShistoHisto->GetYaxis()->SetTitleOffset(1.5);
		fprintf(outputFile, "\nLambda_ThetaPhi_Prompt_CS\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_thetaphi_CS_PR_rap[ptBin]=lambda_thetaphi_CS_PR_[ptBin][yBin];
	 	  					errlambda_thetaphi_CS_PR_rap[ptBin]=errlambda_thetaphi_CS_PR_[ptBin][yBin];
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.6f +- %1.6f\n",ptBin+1,yBin+1,lambda_thetaphi_CS_PR_rap[ptBin],errlambda_thetaphi_CS_PR_rap[ptBin]);

	 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) lambda_thetaphi_CS_PR_rap[ptBin] = offsetValue;

		  					if (ptBin+1==nonConvergeptBinRap1CS1 |ptBin+1==nonConvergeptBinRap1CS2 | ptBin+1==nonConvergeptBinRap1CS3 | ptBin+1==nonConvergeptBinRap1CS4 && yBin==0) lambda_thetaphi_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap2CS1 |ptBin+1==nonConvergeptBinRap2CS2 | ptBin+1==nonConvergeptBinRap2CS3 | ptBin+1==nonConvergeptBinRap2CS4 && yBin==1) lambda_thetaphi_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap3CS1 |ptBin+1==nonConvergeptBinRap3CS2 | ptBin+1==nonConvergeptBinRap3CS3 | ptBin+1==nonConvergeptBinRap3CS4 && yBin==2) lambda_thetaphi_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap4CS1 |ptBin+1==nonConvergeptBinRap4CS2 | ptBin+1==nonConvergeptBinRap4CS3 | ptBin+1==nonConvergeptBinRap4CS4 && yBin==3) lambda_thetaphi_CS_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap5CS1 |ptBin+1==nonConvergeptBinRap5CS2 | ptBin+1==nonConvergeptBinRap5CS3 | ptBin+1==nonConvergeptBinRap5CS4 && yBin==4) lambda_thetaphi_CS_PR_rap[ptBin] = offsetValue;


	     	}

	 TGraphErrors *LambdaThetaPhiPRCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,lambda_thetaphi_CS_PR_rap,errptMean,errlambda_thetaphi_CS_PR_rap);
	 LambdaThetaPhiPRCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPhiPRCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPhiPRCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaThetaPhiPRCSLegend->AddEntry(LambdaThetaPhiPRCSGraph,legendrap,"ple");
	 LambdaThetaPhiPRCSGraph->Draw("P");

	 	  		}

	 	  		LambdaThetaPhiPRCSLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_LambdaThetaPhiPRCS.png");

	 	  		LambdaThetaPhiPRCSCanvas->SaveAs(Filename);

	 	  		LambdaThetaPhiPRCSCanvas->Close();














////////////////////////////////// LAMBDA THETA PROMPT HX ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaThetaPRHXLegend=new TLegend(0.875,0.6,1,0.89);
	 LambdaThetaPRHXLegend->SetFillColor(kWhite);
	 LambdaThetaPRHXLegend->SetTextFont(72);
	 LambdaThetaPRHXLegend->SetTextSize(0.02);
	 LambdaThetaPRHXLegend->SetBorderSize(0);

	 TCanvas *LambdaThetaPRHXCanvas = new TCanvas("Prompt #lambda_{#theta_{HX}}","",1000,700);
	 LambdaThetaPRHXCanvas->SetFillColor(kWhite);
     LambdaThetaPRHXCanvas->SetGrid(0,1);
	 LambdaThetaPRHXCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaThetaPRHXCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaThetaPRHXCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaThetaPRHXhistoHisto = LambdaThetaPRHXCanvas->DrawFrame(4,-1.3,30,1.3);
	 LambdaThetaPRHXhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaThetaPRHXhistoHisto->SetYTitle("Prompt #lambda_{#theta_{HX}}");
	 LambdaThetaPRHXhistoHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nLambda_Theta_Prompt_HX\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_theta_HX_PR_rap[ptBin]=lambda_theta_HX_PR_[ptBin][yBin];
	 	  					errlambda_theta_HX_PR_rap[ptBin]=errlambda_theta_HX_PR_[ptBin][yBin];
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.6f +- %1.6f\n",ptBin+1,yBin+1,lambda_theta_HX_PR_rap[ptBin],errlambda_theta_HX_PR_rap[ptBin]);

		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) lambda_theta_HX_PR_rap[ptBin] = offsetValue;

			  					if (ptBin+1==nonConvergeptBinRap1HX1 |ptBin+1==nonConvergeptBinRap1HX2 | ptBin+1==nonConvergeptBinRap1HX3 | ptBin+1==nonConvergeptBinRap1HX4 && yBin==0) lambda_theta_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap2HX1 |ptBin+1==nonConvergeptBinRap2HX2 | ptBin+1==nonConvergeptBinRap2HX3 | ptBin+1==nonConvergeptBinRap2HX4 && yBin==1) lambda_theta_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap3HX1 |ptBin+1==nonConvergeptBinRap3HX2 | ptBin+1==nonConvergeptBinRap3HX3 | ptBin+1==nonConvergeptBinRap3HX4 && yBin==2) lambda_theta_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap4HX1 |ptBin+1==nonConvergeptBinRap4HX2 | ptBin+1==nonConvergeptBinRap4HX3 | ptBin+1==nonConvergeptBinRap4HX4 && yBin==3) lambda_theta_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap5HX1 |ptBin+1==nonConvergeptBinRap5HX2 | ptBin+1==nonConvergeptBinRap5HX3 | ptBin+1==nonConvergeptBinRap5HX4 && yBin==4) lambda_theta_HX_PR_rap[ptBin] = offsetValue;

	     	}

	 TGraphErrors *LambdaThetaPRHXGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,lambda_theta_HX_PR_rap,errptMean,errlambda_theta_HX_PR_rap);
	 LambdaThetaPRHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPRHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaThetaPRHXLegend->AddEntry(LambdaThetaPRHXGraph,legendrap,"ple");
	 LambdaThetaPRHXGraph->Draw("P");

	 	  		}

	 	  		LambdaThetaPRHXLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_LambdaThetaPRHX.png");

	 	  		LambdaThetaPRHXCanvas->SaveAs(Filename);

	 	  		LambdaThetaPRHXCanvas->Close();

////////////////////////////////// LAMBDA PHI PROMPT HX ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaPhiPRHXLegend=new TLegend(0.875,0.6,1,0.89);
	 LambdaPhiPRHXLegend->SetFillColor(kWhite);
	 LambdaPhiPRHXLegend->SetTextFont(72);
	 LambdaPhiPRHXLegend->SetTextSize(0.02);
	 LambdaPhiPRHXLegend->SetBorderSize(0);

	 TCanvas *LambdaPhiPRHXCanvas = new TCanvas("Prompt #lambda_{#phi_{HX}}","",1000,700);
	 LambdaPhiPRHXCanvas->SetFillColor(kWhite);
     LambdaPhiPRHXCanvas->SetGrid(0,1);
	 LambdaPhiPRHXCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaPhiPRHXCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaPhiPRHXCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaPhiPRHXhistoHisto = LambdaPhiPRHXCanvas->DrawFrame(4,-1.3,30,1.3);
	 LambdaPhiPRHXhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaPhiPRHXhistoHisto->SetYTitle("Prompt #lambda_{#phi_{HX}}");
	 LambdaPhiPRHXhistoHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nLambda_Phi_Prompt_HX\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_phi_HX_PR_rap[ptBin]=lambda_phi_HX_PR_[ptBin][yBin];
	 	  					cout<<"pt,y"<<ptBin<<","<<yBin<<lambda_phi_HX_PR_rap[ptBin]<<endl;
		  					 errlambda_phi_HX_PR_rap[ptBin]=errlambda_phi_HX_PR_[ptBin][yBin];
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.6f +- %1.6f\n",ptBin+1,yBin+1,lambda_phi_HX_PR_rap[ptBin],errlambda_phi_HX_PR_rap[ptBin]);

		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) lambda_phi_HX_PR_rap[ptBin] = offsetValue;

			  					if (ptBin+1==nonConvergeptBinRap1HX1 |ptBin+1==nonConvergeptBinRap1HX2 | ptBin+1==nonConvergeptBinRap1HX3 | ptBin+1==nonConvergeptBinRap1HX4 && yBin==0) lambda_phi_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap2HX1 |ptBin+1==nonConvergeptBinRap2HX2 | ptBin+1==nonConvergeptBinRap2HX3 | ptBin+1==nonConvergeptBinRap2HX4 && yBin==1) lambda_phi_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap3HX1 |ptBin+1==nonConvergeptBinRap3HX2 | ptBin+1==nonConvergeptBinRap3HX3 | ptBin+1==nonConvergeptBinRap3HX4 && yBin==2) lambda_phi_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap4HX1 |ptBin+1==nonConvergeptBinRap4HX2 | ptBin+1==nonConvergeptBinRap4HX3 | ptBin+1==nonConvergeptBinRap4HX4 && yBin==3) lambda_phi_HX_PR_rap[ptBin] = offsetValue;
			  					if (ptBin+1==nonConvergeptBinRap5HX1 |ptBin+1==nonConvergeptBinRap5HX2 | ptBin+1==nonConvergeptBinRap5HX3 | ptBin+1==nonConvergeptBinRap5HX4 && yBin==4) lambda_phi_HX_PR_rap[ptBin] = offsetValue;

		 	  					cout<<"pt,y"<<ptBin<<","<<yBin<<lambda_phi_HX_PR_rap[ptBin]<<endl;
		 	  					cout<<ptMean[ptBin]<<endl;
 	}

	 TGraphErrors *LambdaPhiPRHXGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,lambda_phi_HX_PR_rap,errptMean,errlambda_phi_HX_PR_rap);
	 LambdaPhiPRHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaPhiPRHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaPhiPRHXLegend->AddEntry(LambdaPhiPRHXGraph,legendrap,"ple");
	 LambdaPhiPRHXGraph->Draw("P");

	 	  		}

	 	  		LambdaPhiPRHXLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_LambdaPhiPRHX.png");

	 	  		LambdaPhiPRHXCanvas->SaveAs(Filename);

	 	  		LambdaPhiPRHXCanvas->Close();



////////////////////////////////// LAMBDA THETA PROMPT HX ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* LambdaThetaPhiPRHXLegend=new TLegend(0.875,0.6,1,0.89);
	 LambdaThetaPhiPRHXLegend->SetFillColor(kWhite);
	 LambdaThetaPhiPRHXLegend->SetTextFont(72);
	 LambdaThetaPhiPRHXLegend->SetTextSize(0.02);
	 LambdaThetaPhiPRHXLegend->SetBorderSize(0);

	 TCanvas *LambdaThetaPhiPRHXCanvas = new TCanvas("Prompt #lambda_{#theta #phi_{HX}}","",1000,700);
	 LambdaThetaPhiPRHXCanvas->SetFillColor(kWhite);
     LambdaThetaPhiPRHXCanvas->SetGrid(0,1);
	 LambdaThetaPhiPRHXCanvas->GetFrame()->SetFillColor(kWhite);
	 LambdaThetaPhiPRHXCanvas->GetFrame()->SetBorderSize(0);
//	 LambdaThetaPhiPRHXCanvas->SetLeftMargin(0.15) ;

	 TH1F *LambdaThetaPhiPRHXhistoHisto = LambdaThetaPhiPRHXCanvas->DrawFrame(4,-1.3,30,1.3);
	 LambdaThetaPhiPRHXhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 LambdaThetaPhiPRHXhistoHisto->SetYTitle("Prompt #lambda_{#theta #phi_{HX}}");
	 LambdaThetaPhiPRHXhistoHisto->GetYaxis()->SetTitleOffset(1.5);
		fprintf(outputFile, "\nLambda_ThetaPhi_Prompt_HX\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin];
		  					   errptMean[ptBin]=0;
		  					 lambda_thetaphi_HX_PR_rap[ptBin]=lambda_thetaphi_HX_PR_[ptBin][yBin];
	 	  					errlambda_thetaphi_HX_PR_rap[ptBin]=errlambda_thetaphi_HX_PR_[ptBin][yBin];
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.6f +- %1.6f\n",ptBin+1,yBin+1,lambda_thetaphi_HX_PR_rap[ptBin],errlambda_thetaphi_HX_PR_rap[ptBin]);

	 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) lambda_thetaphi_HX_PR_rap[ptBin] = offsetValue;

		  					if (ptBin+1==nonConvergeptBinRap1HX1 |ptBin+1==nonConvergeptBinRap1HX2 | ptBin+1==nonConvergeptBinRap1HX3 | ptBin+1==nonConvergeptBinRap1HX4 && yBin==0) lambda_thetaphi_HX_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap2HX1 |ptBin+1==nonConvergeptBinRap2HX2 | ptBin+1==nonConvergeptBinRap2HX3 | ptBin+1==nonConvergeptBinRap2HX4 && yBin==1) lambda_thetaphi_HX_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap3HX1 |ptBin+1==nonConvergeptBinRap3HX2 | ptBin+1==nonConvergeptBinRap3HX3 | ptBin+1==nonConvergeptBinRap3HX4 && yBin==2) lambda_thetaphi_HX_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap4HX1 |ptBin+1==nonConvergeptBinRap4HX2 | ptBin+1==nonConvergeptBinRap4HX3 | ptBin+1==nonConvergeptBinRap4HX4 && yBin==3) lambda_thetaphi_HX_PR_rap[ptBin] = offsetValue;
		  					if (ptBin+1==nonConvergeptBinRap5HX1 |ptBin+1==nonConvergeptBinRap5HX2 | ptBin+1==nonConvergeptBinRap5HX3 | ptBin+1==nonConvergeptBinRap5HX4 && yBin==4) lambda_thetaphi_HX_PR_rap[ptBin] = offsetValue;


	     	}

	 TGraphErrors *LambdaThetaPhiPRHXGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,lambda_thetaphi_HX_PR_rap,errptMean,errlambda_thetaphi_HX_PR_rap);
	 LambdaThetaPhiPRHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPhiPRHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 LambdaThetaPhiPRHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 LambdaThetaPhiPRHXLegend->AddEntry(LambdaThetaPhiPRHXGraph,legendrap,"ple");
	 LambdaThetaPhiPRHXGraph->Draw("P");

	 	  		}

	 	  		LambdaThetaPhiPRHXLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_LambdaThetaPhiPRHX.png");

	 	  		LambdaThetaPhiPRHXCanvas->SaveAs(Filename);

	 	  		LambdaThetaPhiPRHXCanvas->Close();














////////////////////////////////// INVARIANT LAMBDA PROMPT ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* PromptInvariantLambdaLegend=new TLegend(0.875,0.6,1,0.89);
	 PromptInvariantLambdaLegend->SetFillColor(kWhite);
	 PromptInvariantLambdaLegend->SetTextFont(72);
	 PromptInvariantLambdaLegend->SetTextSize(0.02);
	 PromptInvariantLambdaLegend->SetBorderSize(0);

	 TCanvas *PromptInvariantLambdaCanvas = new TCanvas("Invariant Prompt #lambda_{CS}","",1000,700);
	 PromptInvariantLambdaCanvas->SetFillColor(kWhite);
     PromptInvariantLambdaCanvas->SetGrid(0,1);
	 PromptInvariantLambdaCanvas->GetFrame()->SetFillColor(kWhite);
	 PromptInvariantLambdaCanvas->GetFrame()->SetBorderSize(0);
//	 PromptInvariantLambdaCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptInvariantLambdahistoHisto = PromptInvariantLambdaCanvas->DrawFrame(0,-1.3,30,1.3);
	 PromptInvariantLambdahistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptInvariantLambdahistoHisto->SetYTitle("Invariant Prompt #lambda_{CS}");
	 PromptInvariantLambdahistoHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nInvariant_Lambda_Prompt\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      			   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					invariant_lambda_CS_PR_rap[ptBin]=invariant_lambda_CS_PR_[ptBin][yBin];
	 	  					errinvariant_lambda_CS_PR_rap[ptBin]=0;
	 	  					invariant_lambda_HX_PR_rap[ptBin]=invariant_lambda_HX_PR_[ptBin][yBin];
	 	  					errinvariant_lambda_HX_PR_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d CS %1.4f +- %1.4f\n",ptBin+1,yBin+1,invariant_lambda_CS_PR_rap[ptBin],errinvariant_lambda_CS_PR_rap[ptBin]);
	 	  					fprintf(outputFile, "pt%d_rapidity%d HX %1.4f +- %1.4f\n",ptBin+1,yBin+1,invariant_lambda_HX_PR_rap[ptBin],errinvariant_lambda_HX_PR_rap[ptBin]);

		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) invariant_lambda_HX_PR_rap[ptBin] = -10000;
		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) invariant_lambda_CS_PR_rap[ptBin] = -10000;

			 	  				if (ptBin==2 |ptBin==3 | ptBin==6 && yBin==0) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
				 	  			if (ptBin==5 && yBin==1)  {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
				 	  			if (ptBin==4 | ptBin==7 && yBin==2) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
				 	  			if (ptBin==2 && yBin==0)  {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
				 	  			if (ptBin==4 && yBin==2) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}

			  					if (ptBin+1==nonConvergeptBinRap1HX1 |ptBin+1==nonConvergeptBinRap1HX2 | ptBin+1==nonConvergeptBinRap1HX3 | ptBin+1==nonConvergeptBinRap1HX4 && yBin==0) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap2HX1 |ptBin+1==nonConvergeptBinRap2HX2 | ptBin+1==nonConvergeptBinRap2HX3 | ptBin+1==nonConvergeptBinRap2HX4 && yBin==1) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap3HX1 |ptBin+1==nonConvergeptBinRap3HX2 | ptBin+1==nonConvergeptBinRap3HX3 | ptBin+1==nonConvergeptBinRap3HX4 && yBin==2) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap4HX1 |ptBin+1==nonConvergeptBinRap4HX2 | ptBin+1==nonConvergeptBinRap4HX3 | ptBin+1==nonConvergeptBinRap4HX4 && yBin==3) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap5HX1 |ptBin+1==nonConvergeptBinRap5HX2 | ptBin+1==nonConvergeptBinRap5HX3 | ptBin+1==nonConvergeptBinRap5HX4 && yBin==4) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}

			  					if (ptBin+1==nonConvergeptBinRap1CS1 |ptBin+1==nonConvergeptBinRap1CS2 | ptBin+1==nonConvergeptBinRap1CS3 | ptBin+1==nonConvergeptBinRap1CS4 && yBin==0) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap2CS1 |ptBin+1==nonConvergeptBinRap2CS2 | ptBin+1==nonConvergeptBinRap2CS3 | ptBin+1==nonConvergeptBinRap2CS4 && yBin==1) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap3CS1 |ptBin+1==nonConvergeptBinRap3CS2 | ptBin+1==nonConvergeptBinRap3CS3 | ptBin+1==nonConvergeptBinRap3CS4 && yBin==2) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap4CS1 |ptBin+1==nonConvergeptBinRap4CS2 | ptBin+1==nonConvergeptBinRap4CS3 | ptBin+1==nonConvergeptBinRap4CS4 && yBin==3) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}
			  					if (ptBin+1==nonConvergeptBinRap5CS1 |ptBin+1==nonConvergeptBinRap5CS2 | ptBin+1==nonConvergeptBinRap5CS3 | ptBin+1==nonConvergeptBinRap5CS4 && yBin==4) {invariant_lambda_HX_PR_rap[ptBin] = -10000;invariant_lambda_CS_PR_rap[ptBin] = -10000;}


	     	}

	 TGraphErrors *PromptInvariantLambdaCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,invariant_lambda_CS_PR_rap,errptMean,errinvariant_lambda_CS_PR_rap);
	 PromptInvariantLambdaCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d, CS",yBin+1);
	 PromptInvariantLambdaLegend->AddEntry(PromptInvariantLambdaCSGraph,legendrap,"ple");
	 PromptInvariantLambdaCSGraph->Draw("P");

	 TGraphErrors *PromptInvariantLambdaHXGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,invariant_lambda_HX_PR_rap,errptMean,errinvariant_lambda_HX_PR_rap);
	 PromptInvariantLambdaHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantLambdaHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]+3);
	 sprintf(legendrap,"Rapidity %d, HX",yBin+1);
	 PromptInvariantLambdaLegend->AddEntry(PromptInvariantLambdaHXGraph,legendrap,"ple");
	 PromptInvariantLambdaHXGraph->Draw("P");

	 	  		}

	 	  		PromptInvariantLambdaLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_PromptInvariantLambda.png");

	 	  		PromptInvariantLambdaCanvas->SaveAs(Filename);

	 	  		PromptInvariantLambdaCanvas->Close();


////////////////////////////////// INVARIANT F PROMPT ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* PromptInvariantFLegend=new TLegend(0.875,0.6,1,0.89);
	 PromptInvariantFLegend->SetFillColor(kWhite);
	 PromptInvariantFLegend->SetTextFont(72);
	 PromptInvariantFLegend->SetTextSize(0.02);
	 PromptInvariantFLegend->SetBorderSize(0);

	 TCanvas *PromptInvariantFCanvas = new TCanvas("Invariant Prompt F_{CS}","",1000,700);
	 PromptInvariantFCanvas->SetFillColor(kWhite);
     PromptInvariantFCanvas->SetGrid(0,1);
	 PromptInvariantFCanvas->GetFrame()->SetFillColor(kWhite);
	 PromptInvariantFCanvas->GetFrame()->SetBorderSize(0);
//	 PromptInvariantFCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptInvariantFhistoHisto = PromptInvariantFCanvas->DrawFrame(0,-1.3,30,1.3);
	 PromptInvariantFhistoHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptInvariantFhistoHisto->SetYTitle("Invariant Prompt F_{CS}");
	 PromptInvariantFhistoHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nInvariant_F_Prompt\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 invariant_F_CS_PR_rap[ptBin]=invariant_F_CS_PR_[ptBin][yBin];
	 	  					errinvariant_F_CS_PR_rap[ptBin]=0;
	 	  					invariant_F_HX_PR_rap[ptBin]=invariant_F_HX_PR_[ptBin][yBin];
	 	  					errinvariant_F_HX_PR_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d CS %1.4f +- %1.4f\n",ptBin+1,yBin+1,invariant_F_CS_PR_rap[ptBin],errinvariant_F_CS_PR_rap[ptBin]);
	 	  					fprintf(outputFile, "pt%d_rapidity%d HX %1.4f +- %1.4f\n",ptBin+1,yBin+1,invariant_F_HX_PR_rap[ptBin],errinvariant_F_HX_PR_rap[ptBin]);

		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) invariant_F_CS_PR_rap[ptBin] = -10000;
		 	  		 		  if(jpsi::pTRange[yBin+1][ptBin]<ptplotMin) invariant_F_HX_PR_rap[ptBin] = -10000;


		  					if (ptBin+1==nonConvergeptBinRap1HX1 |ptBin+1==nonConvergeptBinRap1HX2 | ptBin+1==nonConvergeptBinRap1HX3 | ptBin+1==nonConvergeptBinRap1HX4 && yBin==0) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap2HX1 |ptBin+1==nonConvergeptBinRap2HX2 | ptBin+1==nonConvergeptBinRap2HX3 | ptBin+1==nonConvergeptBinRap2HX4 && yBin==1) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap3HX1 |ptBin+1==nonConvergeptBinRap3HX2 | ptBin+1==nonConvergeptBinRap3HX3 | ptBin+1==nonConvergeptBinRap3HX4 && yBin==2) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap4HX1 |ptBin+1==nonConvergeptBinRap4HX2 | ptBin+1==nonConvergeptBinRap4HX3 | ptBin+1==nonConvergeptBinRap4HX4 && yBin==3) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap5HX1 |ptBin+1==nonConvergeptBinRap5HX2 | ptBin+1==nonConvergeptBinRap5HX3 | ptBin+1==nonConvergeptBinRap5HX4 && yBin==4) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}

		  					if (ptBin+1==nonConvergeptBinRap1CS1 |ptBin+1==nonConvergeptBinRap1CS2 | ptBin+1==nonConvergeptBinRap1CS3 | ptBin+1==nonConvergeptBinRap1CS4 && yBin==0) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap2CS1 |ptBin+1==nonConvergeptBinRap2CS2 | ptBin+1==nonConvergeptBinRap2CS3 | ptBin+1==nonConvergeptBinRap2CS4 && yBin==1) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap3CS1 |ptBin+1==nonConvergeptBinRap3CS2 | ptBin+1==nonConvergeptBinRap3CS3 | ptBin+1==nonConvergeptBinRap3CS4 && yBin==2) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap4CS1 |ptBin+1==nonConvergeptBinRap4CS2 | ptBin+1==nonConvergeptBinRap4CS3 | ptBin+1==nonConvergeptBinRap4CS4 && yBin==3) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}
		  					if (ptBin+1==nonConvergeptBinRap5CS1 |ptBin+1==nonConvergeptBinRap5CS2 | ptBin+1==nonConvergeptBinRap5CS3 | ptBin+1==nonConvergeptBinRap5CS4 && yBin==4) {invariant_F_HX_PR_rap[ptBin] = -10000;invariant_F_CS_PR_rap[ptBin] = -10000;}


	     	}

	 TGraphErrors *PromptInvariantFCSGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,invariant_F_CS_PR_rap,errptMean,errinvariant_F_CS_PR_rap);
	 PromptInvariantFCSGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFCSGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFCSGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d, CS",yBin+1);
	 PromptInvariantFLegend->AddEntry(PromptInvariantFCSGraph,legendrap,"ple");
	 PromptInvariantFCSGraph->Draw("P");

	 TGraphErrors *PromptInvariantFHXGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,invariant_F_HX_PR_rap,errptMean,errinvariant_F_HX_PR_rap);
	 PromptInvariantFHXGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFHXGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptInvariantFHXGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]+3);
	 sprintf(legendrap,"Rapidity %d, HX",yBin+1);
	 PromptInvariantFLegend->AddEntry(PromptInvariantFHXGraph,legendrap,"ple");
	 PromptInvariantFHXGraph->Draw("P");

	 	  		}

	 	  		PromptInvariantFLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_PromptInvariantF.png");

	 	  		PromptInvariantFCanvas->SaveAs(Filename);

	 	  		PromptInvariantFCanvas->Close();




////////////////////////////////// Prompt l_J/psi cut ///////////////////////////////////////////////////////////////////////////////////

	 if(!nofrac){

	 TLegend* PromptCutLegend=new TLegend(0.875,0.6,1,0.89);
	 PromptCutLegend->SetFillColor(kWhite);
	 PromptCutLegend->SetTextFont(72);
	 PromptCutLegend->SetTextSize(0.02);
	 PromptCutLegend->SetBorderSize(0);

	 TCanvas *PromptCutCanvas = new TCanvas("Prompt l_{J/#psi} [mm]","",1000,700);
	 PromptCutCanvas->SetFillColor(kWhite);
	 PromptCutCanvas->SetGrid(0,1);
     PromptCutCanvas->GetFrame()->SetFillColor(kWhite);
     PromptCutCanvas->GetFrame()->SetBorderSize(0);
//	 PromptCutCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptCutHisto = PromptCutCanvas->DrawFrame(0,-1,30,2.5);
	 PromptCutHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptCutHisto->SetYTitle("l_{J/#psi}_{max} [mm]");
	 PromptCutHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nPrompt l_J/psi cut\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 prompt_lifetime_cut_rap[ptBin]=prompt_lifetime_cut[ptBin][yBin];
	 	  					errprompt_lifetime_cut_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.4f\n",ptBin+1,yBin+1,prompt_lifetime_cut_rap[ptBin]);

	     	}

	 TGraphErrors *PromptCutGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,prompt_lifetime_cut_rap,errptMean,errprompt_lifetime_cut_rap);
	 PromptCutGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptCutGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptCutGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 PromptCutLegend->AddEntry(PromptCutGraph,legendrap,"ple");
	 PromptCutGraph->Draw("P");


	 	  		}

	 PromptCutLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_PromptCut.png");

	 	  		PromptCutCanvas->SaveAs(Filename);

	 	  		PromptCutCanvas->Close();


////////////////////////////////// Non Prompt fraction above cut ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* NonPromptFracLegend=new TLegend(0.875,0.6,1,0.89);
	 NonPromptFracLegend->SetFillColor(kWhite);
	 NonPromptFracLegend->SetTextFont(72);
	 NonPromptFracLegend->SetTextSize(0.02);
	 NonPromptFracLegend->SetBorderSize(0);

	 TCanvas *NonPromptFracCanvas = new TCanvas("Prompt l_{J/#psi} [mm]","",1000,700);
	 NonPromptFracCanvas->SetFillColor(kWhite);
	 NonPromptFracCanvas->SetGrid(0,1);
	 NonPromptFracCanvas->GetFrame()->SetFillColor(kWhite);
	 NonPromptFracCanvas->GetFrame()->SetBorderSize(0);
//	 NonPromptFracCanvas->SetLeftMargin(0.15) ;

	 TH1F *NonPromptFracHisto = PromptCutCanvas->DrawFrame(0,0,30,1);
	 NonPromptFracHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 NonPromptFracHisto->SetYTitle("Fraction of Non Prompt J/  #psi accepted by Cut");
	 NonPromptFracHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nNonPrompt Fraction accepted by Cut\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 nonprompt_lifetime_frac_rap[ptBin]=nonprompt_lifetime_frac[ptBin][yBin];
	 	  					errnonprompt_lifetime_frac_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.4f\n",ptBin+1,yBin+1,nonprompt_lifetime_frac_rap[ptBin]);

	     	}

	 TGraphErrors *NonPromptFracGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,nonprompt_lifetime_frac_rap,errptMean,errnonprompt_lifetime_frac_rap);
	 NonPromptFracGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptFracGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptFracGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 NonPromptFracLegend->AddEntry(NonPromptFracGraph,legendrap,"ple");
	 NonPromptFracGraph->Draw("P");


	 	  		}

	 NonPromptFracLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_NonPromptFrac.png");

	 	  		NonPromptFracCanvas->SaveAs(Filename);

	 	  		NonPromptFracCanvas->Close();



////////////////////////////////// Fraction of NP J/#Psi's in region ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* NonPromptAbsFracLegend=new TLegend(0.875,0.6,1,0.89);
	 NonPromptAbsFracLegend->SetFillColor(kWhite);
	 NonPromptAbsFracLegend->SetTextFont(72);
	 NonPromptAbsFracLegend->SetTextSize(0.02);
	 NonPromptAbsFracLegend->SetBorderSize(0);

	 TCanvas *NonPromptAbsFracCanvas = new TCanvas("Prompt l_{J/#psi} [mm]","",1000,700);
	 NonPromptAbsFracCanvas->SetFillColor(kWhite);
	 NonPromptAbsFracCanvas->SetGrid(0,1);
	 NonPromptAbsFracCanvas->GetFrame()->SetFillColor(kWhite);
	 NonPromptAbsFracCanvas->GetFrame()->SetBorderSize(0);
//	 NonPromptAbsFracCanvas->SetLeftMargin(0.15) ;

	 TH1F *NonPromptAbsFracHisto = PromptCutCanvas->DrawFrame(0,0,30,1);
	 NonPromptAbsFracHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 NonPromptAbsFracHisto->SetYTitle("Fraction of NP J/  #Psi's in region");
	 NonPromptAbsFracHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nFraction of NP J/Psi's in region\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 nonprompt_lifetime_absfrac_rap[ptBin]=nonprompt_lifetime_absfrac[ptBin][yBin];
	 	  					errnonprompt_lifetime_absfrac_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.4f\n",ptBin+1,yBin+1,nonprompt_lifetime_absfrac_rap[ptBin]);

	     	}

	 TGraphErrors *NonPromptAbsFracGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,nonprompt_lifetime_absfrac_rap,errptMean,errnonprompt_lifetime_absfrac_rap);
	 NonPromptAbsFracGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptAbsFracGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptAbsFracGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 NonPromptAbsFracLegend->AddEntry(NonPromptAbsFracGraph,legendrap,"ple");
	 NonPromptAbsFracGraph->Draw("P");


	 	  		}

	 NonPromptAbsFracLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_NonPromptAbsFrac.png");

	 	  		NonPromptAbsFracCanvas->SaveAs(Filename);

	 	  		NonPromptAbsFracCanvas->Close();


////////////////////////////////// Background fraction above cut ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* BackgroundFracLegend=new TLegend(0.875,0.6,1,0.89);
	 BackgroundFracLegend->SetFillColor(kWhite);
	 BackgroundFracLegend->SetTextFont(72);
	 BackgroundFracLegend->SetTextSize(0.02);
	 BackgroundFracLegend->SetBorderSize(0);

	 TCanvas *BackgroundFracCanvas = new TCanvas("Prompt l_{J/#psi} [mm]","",1000,700);
	 BackgroundFracCanvas->SetFillColor(kWhite);
	 BackgroundFracCanvas->SetGrid(0,1);
	 BackgroundFracCanvas->GetFrame()->SetFillColor(kWhite);
	 BackgroundFracCanvas->GetFrame()->SetBorderSize(0);
//	 BackgroundFracCanvas->SetLeftMargin(0.15) ;

	 TH1F *BackgroundFracHisto = PromptCutCanvas->DrawFrame(0,0,30,1);
	 BackgroundFracHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 BackgroundFracHisto->SetYTitle("Fraction of Background accepted by Cut");
	 BackgroundFracHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nBackground Fraction accepted by Cut\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 background_lifetime_frac_rap[ptBin]=background_lifetime_frac[ptBin][yBin];
	 	  					errbackground_lifetime_frac_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.4f\n",ptBin+1,yBin+1,background_lifetime_frac_rap[ptBin]);

	     	}

	 TGraphErrors *BackgroundFracGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,background_lifetime_frac_rap,errptMean,errbackground_lifetime_frac_rap);
	 BackgroundFracGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 BackgroundFracGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 BackgroundFracGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 BackgroundFracLegend->AddEntry(BackgroundFracGraph,legendrap,"ple");
	 BackgroundFracGraph->Draw("P");


	 	  		}

	 BackgroundFracLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_BackgroundFrac.png");

	 	  		BackgroundFracCanvas->SaveAs(Filename);

	 	  		BackgroundFracCanvas->Close();





////////////////////////////////// NonPrompt l_J/psi cut ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* NonPromptCutLegend=new TLegend(0.875,0.6,1,0.89);
	 NonPromptCutLegend->SetFillColor(kWhite);
	 NonPromptCutLegend->SetTextFont(72);
	 NonPromptCutLegend->SetTextSize(0.02);
	 NonPromptCutLegend->SetBorderSize(0);

	 TCanvas *NonPromptCutCanvas = new TCanvas("NonPrompt l_{J/#psi} [mm]","",1000,700);
	 NonPromptCutCanvas->SetFillColor(kWhite);
	 NonPromptCutCanvas->SetGrid(0,1);
     NonPromptCutCanvas->GetFrame()->SetFillColor(kWhite);
     NonPromptCutCanvas->GetFrame()->SetBorderSize(0);
//	 NonPromptCutCanvas->SetLeftMargin(0.15) ;

	 TH1F *NonPromptCutHisto = NonPromptCutCanvas->DrawFrame(0,-0.7,30,0);
	 NonPromptCutHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 NonPromptCutHisto->SetYTitle("l_{J/#psi}_min [mm]");
	 NonPromptCutHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nNonPrompt l_J/psi cut\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 nonprompt_lifetime_cut_rap[ptBin]=nonprompt_lifetime_cut[ptBin][yBin];
	 	  					errnonprompt_lifetime_cut_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.4f\n",ptBin+1,yBin+1,nonprompt_lifetime_cut_rap[ptBin]);

	     	}

	 TGraphErrors *NonPromptCutGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,nonprompt_lifetime_cut_rap,errptMean,errnonprompt_lifetime_cut_rap);
	 NonPromptCutGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptCutGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 NonPromptCutGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 NonPromptCutLegend->AddEntry(NonPromptCutGraph,legendrap,"ple");
	 NonPromptCutGraph->Draw("P");


	 	  		}

	 NonPromptCutLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_NonPromptCut.png");

	 	  		NonPromptCutCanvas->SaveAs(Filename);

	 	  		NonPromptCutCanvas->Close();


////////////////////////////////// Prompt fraction below cut ///////////////////////////////////////////////////////////////////////////////////

	 TLegend* PromptFracLegend=new TLegend(0.875,0.6,1,0.89);
	 PromptFracLegend->SetFillColor(kWhite);
	 PromptFracLegend->SetTextFont(72);
	 PromptFracLegend->SetTextSize(0.02);
	 PromptFracLegend->SetBorderSize(0);

	 TCanvas *PromptFracCanvas = new TCanvas("Prompt l_{J/#psi} [mm]","",1000,700);
	 PromptFracCanvas->SetFillColor(kWhite);
	 PromptFracCanvas->SetGrid(0,1);
	 PromptFracCanvas->GetFrame()->SetFillColor(kWhite);
	 PromptFracCanvas->GetFrame()->SetBorderSize(0);
//	 PromptFracCanvas->SetLeftMargin(0.15) ;

	 TH1F *PromptFracHisto = PromptCutCanvas->DrawFrame(0,0,30,1);
	 PromptFracHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	 PromptFracHisto->SetYTitle("Fraction of Prompt J/  #psi accepted by cut");
	 PromptFracHisto->GetYaxis()->SetTitleOffset(1.5);

		fprintf(outputFile, "\nPrompt Fraction below Prompt Cut\n");

	 for(int yBin = yBinMin-1; yBin < yBinMax; ++yBin) {
	     	for(int ptBin = ptBinMin-1; ptBin < jpsi::kNbPTBins[yBin+1]-ptBinMax; ++ptBin) {

	 	      				   ptMean[ptBin]=jpsi::pTWCentre_rap[yBin+1][ptBin]+(yBin)*0.1;
		  					   errptMean[ptBin]=0;
		  					 prompt_lifetime_frac_rap[ptBin]=prompt_lifetime_frac[ptBin][yBin];
	 	  					errprompt_lifetime_frac_rap[ptBin]=0;
	 	  					fprintf(outputFile, "pt%d_rapidity%d %1.4f\n",ptBin+1,yBin+1,prompt_lifetime_frac_rap[ptBin]);

	     	}

	 TGraphErrors *PromptFracGraph = new TGraphErrors(jpsi::kNbPTBins[yBin+1],ptMean,prompt_lifetime_frac_rap,errptMean,errprompt_lifetime_frac_rap);
	 PromptFracGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptFracGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	 PromptFracGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	 sprintf(legendrap,"Rapidity %d",yBin+1);
	 PromptFracLegend->AddEntry(PromptFracGraph,legendrap,"ple");
	 PromptFracGraph->Draw("P");


	 	  		}

	 PromptFracLegend->Draw();

	 	  		  sprintf(Filename,"Plots/ParameterPlots/pT_PromptFrac.png");

	 	  		PromptFracCanvas->SaveAs(Filename);

	 	  		PromptFracCanvas->Close();


	 }















  delete fInputCS;
  delete fInputHX;

  fclose(outputFile);

  return 0;
}
