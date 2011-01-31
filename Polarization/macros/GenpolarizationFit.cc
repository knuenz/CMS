///////
// This program is meant to extract the polarization parameters from the data.
// It takes the input parameters from the previously run "extract*" programs
// and then runs the constrained fits to extract the polarization of the prompt
// and non-prompt components.
// NOTE: For now you need to change the input file names in the code and recompile right now.... 
//       Will write a nice commandline parser later.
// \author Lindsey Gray (UW Madison)
//////

#include <iostream>
#include <sstream>

//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"

// RooFit Includes
#include "RooAddPdf.h"
#include "RooExtendPdf.h"///////
// This program is meant to extract the polarization parameters from the data.
// It takes the input parameters from the previously run "extract*" programs
// and then runs the constrained fits to extract the polarization of the prompt
// and non-prompt components.
// NOTE: For now you need to change the input file names in the code and recompile right now....
//       Will write a nice commandline parser later.
// \author Lindsey Gray (UW Madison)
//////

#include <iostream>
#include <sstream>

//J/Psi common vars
//#include "commonVar.h"

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

int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

  cout<<"Hello Grid :)"<<endl;

  bool doprompt(true), dononprompt(true), dobkg(true),dopol(true),pereverr(false),newdata(false),fitHX(true),fitCS(true);

  gSystem->mkdir("fitRoots");



  gStyle->SetPalette(1);
  gStyle->SetTitleFillColor(kWhite);

  dononprompt=false;
  dobkg=false;
  pereverr=true;

  for( int i=0;i < argc; ++i ) {
    if(std::string(argv[i]).find("--noNonPrompt") != std::string::npos) dononprompt = false;
    if(std::string(argv[i]).find("--noPrompt") != std::string::npos) doprompt = false;
    if(std::string(argv[i]).find("--noBackground") != std::string::npos) dobkg = false;
    if(std::string(argv[i]).find("--noPol") != std::string::npos) dopol = false;
    if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
    if(std::string(argv[i]).find("--newTTree") != std::string::npos) newdata = true;
    if(std::string(argv[i]).find("--noHX") != std::string::npos) fitHX = false;
    if(std::string(argv[i]).find("--noCS") != std::string::npos) fitCS = false;

  }



  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.4,2.4);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",0.5,1.5);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,500);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  if(pereverr)
    varlist.add(JpsictErr);
  if(newdata)
    varlist.add(HLT_Mu0_TkMu0_Jpsi);
  varlist.add(MCweight);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);


  TFile *fitInput = new TFile("$CMSSW_BASE/src/JPsiPolarizationSave3/data/jPsiFit_Bins_PR_FSR_geomAcc07Jan.root","UPDATE");



  RooDataSet *data = NULL;
  CompositeModelBuilder *hx, *cs;

  char minimizerChar[200];
  sprintf(minimizerChar,"simplex");
  char minimizerType[200];
  sprintf(minimizerType,"Minuit2");

  int iterations = 10;
  int generate_gt=0;

  char output_cs_[200];
	  sprintf(output_cs_,"jPsiFitFinal_cs_ToyMC.root");
	  char output_hx_[200];
	  sprintf(output_hx_,"jPsiFitFinal_hx_ToyMC.root");
	  TFile *output_cs = new TFile(output_cs_,"UPDATE");
	  TFile *output_hx = new TFile(output_hx_,"UPDATE");

	  TDirectory *current_in,*current_out_hx,*current_out_cs;

  for (int iteration = 1; iteration < iterations+1; iteration++){

	  cout<<"GENERATION "<<iteration<<endl;




  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-4; ++yBin) {
 	  for(int ptBin = 1; ptBin < jpsi::kNbPTBins[yBin+1]-5; ++ptBin) {

// 		  if(jpsi::pTRange[yBin+1][ptBin]<6.) continue;


      std::stringstream binName,cutString, binNameOut;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      binNameOut << "pt" << ptBin+1 << "_rapidity" << yBin+1 << "_gen" << iteration;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];
      std::cout << cutString.str() << std::endl;



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


      bool writeProgress(true);
      if(writeProgress){
  	char outputfilename[200];
  	sprintf(outputfilename,"FitProgressToyMC_rap%d_pt%d.txt",yBin+1,ptBin+1);
  	FILE *outputFile = fopen(outputfilename,"w");
	fprintf(outputFile, "rapidity%d_pt%d, iteration %d\n",yBin+1,ptBin+1,iteration);
	fclose(outputFile);
      }



      if(current_in = fitInput->GetDirectory(binName.str().c_str())) {


	hx = new CompositeModelBuilder("HX","");
	cs = new CompositeModelBuilder("CS","");

	hx->setUsePrompt(doprompt);
	hx->setUseNonPrompt(dononprompt);
	hx->setUseBkg(dobkg);
	hx->setUsePol(dopol);
	hx->setUseAcceptanceMaps(dopol);
	hx->setUseMass(false);
	hx->setUseLifetime(false);

	cs->setUsePrompt(doprompt);
	cs->setUseNonPrompt(dononprompt);
	cs->setUseBkg(dobkg);
	cs->setUsePol(dopol);
	cs->setUseAcceptanceMaps(dopol);
	cs->setUseMass(false);

	cs->setUseLifetime(false);

	cs->loadParameters(*current_in);
	hx->loadParameters(*current_in);

	cs->getMassModel()->fix("CBm");
	cs->getMassModel()->fix("CBs");


	hx->getMassModel()->fix("CBm");
	hx->getMassModel()->fix("CBs");


	if(pereverr) {
	  cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
	  hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX);
	} else {
	  cs->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
	  hx->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);
	}

	RooAbsData *thisBinGenCS = cs->model()->generate(RooArgSet(costh_CS,phi_CS),jpsi::numEventsBin[yBin][ptBin]);//thisBin->tree()->GetEntries());
	RooAbsData *thisBinGenHX = hx->model()->generate(RooArgSet(costh_HX,phi_HX),jpsi::numEventsBin[yBin][ptBin]);//thisBin->tree()->GetEntries());

	cout<<"Generation Character: "<<thisBinGenCS->mean(costh_CS)<<" "<<thisBinGenCS->mean(phi_CS)<<endl;

//	cs->getPromptPolarizationModel()->fix("promptlambda_phi_CS");
//	cs->getPromptPolarizationModel()->fix("promptlambda_theta_CS");
//	cs->getPromptPolarizationModel()->fix("promptlambda_thetaphi_CS");


//	hx->getPromptPolarizationModel()->fix("promptlambda_phi_HX");
//	hx->getPromptPolarizationModel()->fix("promptlambda_theta_HX");
//	hx->getPromptPolarizationModel()->fix("promptlambda_thetaphi_HX");


	if(fitCS){

	std::cout << "Collins-Soper model before fit:" << std::endl;
	cs->Print();

	std::cout << "Fitting CS..." << std::endl;
	  if(pereverr) {

		  if(iteration>generate_gt){
		  RooFitResult* csResult;

		  for(int iMig=0;iMig<5;iMig++){
		  cout<<"MigradTry#"<<iMig+1<<endl;
		  csResult = cs->model()->fitTo(*thisBinGenCS,Save(true),RooFit::NumCPU(2),RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			    //   RooFit::Minos(1),
			       RooFit::Strategy(2),
			       RooFit::Range("signalRegionCS"));

	    csResult->Print();

		  }
		  }
	  }

	    else{
	    }


	std::cout << "Collins-Soper model after fit:" << std::endl;
	cs->Print();


	}

	if(fitHX){

	std::cout << std::endl << "Helicity Frame model before fit:" << std::endl;
	hx->Print();

	std::cout << "Fitting HX..." << std::endl;
	  if(pereverr){
		  if(iteration>generate_gt){
			  RooFitResult* hxResult;
			  for(int iMig=0;iMig<5;iMig++){
			  cout<<"MigradTry#"<<iMig+1<<endl;
	    hxResult =  hx->model()->fitTo(*thisBinGenHX,Save(true),RooFit::NumCPU(2),RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			    //   RooFit::Minos(1),
			       RooFit::Strategy(2),
			       RooFit::Range("signalRegionHX"));

	    hxResult->Print();
			  }
		  }

	  }

	  else{
	  }



	std::cout << std::endl << "Helicity Frame model after fit:" << std::endl;
	hx->Print();

	}

	current_out_hx = output_hx->mkdir(binNameOut.str().c_str());
	current_out_cs = output_cs->mkdir(binNameOut.str().c_str());


	cs->saveParameter("nPrompt",*current_out_cs);
	cs->saveParameter("nNonPrompt",*current_out_cs);
	cs->saveParameter("nBackground",*current_out_cs);

	hx->saveParameter("nPrompt",*current_out_hx);
	hx->saveParameter("nNonPrompt",*current_out_hx);
	hx->saveParameter("nBackground",*current_out_hx);

	cs->saveParameters(*current_out_cs);
	hx->saveParameters(*current_out_hx);

	delete cs;
	delete hx;
      }
    }
  }



  }

  output_hx->Write();
   output_cs->Write();

  output_hx->Close();
  output_cs->Close();
  delete output_hx;
   delete output_cs;
  fitInput->Close();

  delete fitInput;

  return 0;
}
