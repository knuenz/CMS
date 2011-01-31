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
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"

int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

  bool doprompt(true), dononprompt(true), dobkg(true),dopol(true),pereverr(false),newdata(false),fitHX(true),fitCS(true);

  gSystem->mkdir("fitRoots");

  dononprompt=false;
  dobkg=false;
  pereverr=true;

  if(argc < 2) {
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
  RooRealVar MCweight("MCweight","MCweight",1);
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  if(pereverr) 
    varlist.add(JpsictErr);
  if(newdata)
    varlist.add(HLT_Mu0_TkMu0_Jpsi);
  varlist.add(MCweight);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);  


//  TChain *samples = new TChain("data");
  TFile *fitInput = new TFile("jPsiFit_PR_FSR_GA07RE17smsm0NEW.root","UPDATE");

  RooDataSet *data = NULL;
  CompositeModelBuilder *hx, *cs;
  
  TFile* fIn;
  double entries[2][8];

  for( unsigned int arg = 1; arg < argc; ++arg ) {
    if(std::string(argv[arg]).find(".root") != std::string::npos) {
      std::cout << "Adding: " << argv[arg] << std::endl;

      fIn = new TFile(argv[arg]);

//      samples->Add(argv[arg]);
    }
  }
  
//  Char_t *fileNameIn = "/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root";
//  TFile* fIn = new TFile(fileNameIn);
  TTree* dataTreesPR = (TTree*)fIn->Get("data");

  data = new RooDataSet("data","Supplied Data Prompt",varlist,Import(*dataTreesPR),WeightVar(MCweight));

//  data = new RooDataSet("data","Concatenated Samples",samples,varlist,0,"MCweight");

  char minimizerChar[200];
  sprintf(minimizerChar,"migrad");

  int iterations = 1;

  for (int iteration = 1; iteration < iterations+1; iteration++){

	  cout<<"GENERATION "<<iteration<<endl;

	  char output_cs_[200];
	  sprintf(output_cs_,"jPsiFitFinal_cs_GA07RE17smsm0NEW.root");
	  char output_hx_[200];
	  sprintf(output_hx_,"jPsiFitFinal_hx_GA07RE17smsm0NEW.root");
	  TFile *output_cs = new TFile(output_cs_,"RECREATE");
	  TFile *output_hx = new TFile(output_hx_,"RECREATE");

  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-3; ++yBin) {
 	  for(int ptBin = 1; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {
      
// 		  if(jpsi::pTRange[yBin+1][ptBin]<6.) continue;
// 		  if(ptBin==0 && yBin==0) continue;


      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];
//      << " && Jpsict < " << 0.1;
      //<< " && JpsiMass > " << jpsi::JpsiMassMin[yBin+1] << " && JpsiMass < " << jpsi::JpsiMassMax[yBin+1];
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

//      costh_CS.setRange("signalRegionCS",-1,1);
//      costh_HX.setRange("signalRegionHX",-1,1);
//      phi_CS.setRange("signalRegionCS",0,90);
//      phi_HX.setRange("signalRegionHX",0,90);


      
  	char outputfilename[200];
  	sprintf(outputfilename,"FitProgress_RealMC_smsm0.txt");
  	FILE *outputFile = fopen(outputfilename,"w");
	fprintf(outputFile, "rapidity%d_pt%d, iteration %d\n",yBin+1,ptBin+1,iteration);
	fclose(outputFile);

	TDirectory *current_in,*current_out_hx,*current_out_cs;

      if(current_in = fitInput->GetDirectory(binName.str().c_str())) {


	RooAbsData *thisBin = data->reduce((cutString.str()).c_str());

	entries[yBin][ptBin]=thisBin->numEntries();

	thisBin->Print();

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

//	RooAbsData *thisBinGenCS = cs->model()->generate(RooArgSet(costh_CS,phi_CS),jpsi::numEventsBin[yBin][ptBin]);//thisBin->tree()->GetEntries());
//	RooAbsData *thisBinGenHX = hx->model()->generate(RooArgSet(costh_HX,phi_HX),jpsi::numEventsBin[yBin][ptBin]);//thisBin->tree()->GetEntries());


//	cs->getPromptPolarizationModel()->fix("promptlambda_phi_CS");
//	cs->getPromptPolarizationModel()->fix("promptlambda_theta_CS");
//	cs->getPromptPolarizationModel()->fix("promptlambda_thetaphi_CS");


//	hx->getPromptPolarizationModel()->fix("promptlambda_phi_HX");
//	hx->getPromptPolarizationModel()->fix("promptlambda_theta_HX");
//	hx->getPromptPolarizationModel()->fix("promptlambda_thetaphi_HX");


	/*
	cs->getLifetimeModel()->unfix("MeanPrompt1");
	cs->getLifetimeModel()->unfix("CoefPrompt1");
	cs->getLifetimeModel()->unfix("CoefPrompt2");
	*/
	//cs->getLifetimeModel()->unfix("SigmaPrompt1");
	/*
	cs->getLifetimeModel()->unfix("SigmaPrompt2");

	hx->getLifetimeModel()->unfix("MeanPrompt1");
	hx->getLifetimeModel()->unfix("CoefPrompt1");
	hx->getLifetimeModel()->unfix("CoefPrompt2");
	*/
	//hx->getLifetimeModel()->unfix("SigmaPrompt1");
	/*
	hx->getLifetimeModel()->unfix("SigmaPrompt2");
	*/

	if(fitCS){

	std::cout << "Collins-Soper model before fit:" << std::endl;
	cs->Print();
		
	std::cout << "Fitting CS..." << std::endl;
	if(!(ptBin==0 && (yBin == 0))) {
	  if(pereverr) {
	    RooFitResult* csResult = cs->model()->fitTo(*thisBin,Save(true),RooFit::NumCPU(2),//RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       RooFit::Minos(1),
			       //RooFit::InitialHesse(true),
			       RooFit::Strategy(2),
			       RooFit::Range("signalRegionCS"));
			       //RooFit::Strategy(1));
//			       RooFit::PrintEvalErrors(-1));
				   //RooFit::PrintLevel(3);
	    //RooFit::ConditionalObservables(RooArgSet(JpsictErr)));

	    csResult->Print();

	  }
	    else{
	    	RooFitResult* csResult = cs->model()->fitTo(*thisBin,Save(true),RooFit::Minimizer("Minuit",minimizerChar),RooFit::NumCPU(2),RooFit::Timer(true),
				       RooFit::Extended(true),RooFit::SumW2Error(true),
				       //RooFit::Minos(1),
				       //RooFit::InitialHesse(true),
				       RooFit::Strategy(2),
				       RooFit::Range("signalRegionCS"));
				       //RooFit::Strategy(1));
	//			       RooFit::PrintEvalErrors(-1));
					   //RooFit::PrintLevel(3);
	    csResult->Print();

	    }
	    }
	
	std::cout << "Collins-Soper model after fit:" << std::endl;
	cs->Print();


	}

	if(fitHX){

	std::cout << std::endl << "Helicity Frame model before fit:" << std::endl;
	hx->Print();

	std::cout << "Fitting HX..." << std::endl;
	if(!(ptBin==0 && (yBin == 0))) {
	  if(pereverr){
		  RooFitResult* hxResult =  hx->model()->fitTo(*thisBin,Save(true),RooFit::NumCPU(2),//RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       RooFit::Minos(1),
			       //RooFit::InitialHesse(true),
			       RooFit::Strategy(2),
			       RooFit::Range("signalRegionHX"));
			       //RooFit::Strategy(1));
//			       RooFit::PrintEvalErrors(-1));
				   //RooFit::PrintLevel(3);
	  // RooFit::ConditionalObservables(RooArgSet(JpsictErr)));
	    hxResult->Print();

	  }
	  else{
		  RooFitResult* hxResult = hx->model()->fitTo(*thisBin,Save(true),RooFit::Minimizer("Minuit",minimizerChar),RooFit::NumCPU(2),RooFit::Timer(true),
			       RooFit::Extended(true),RooFit::SumW2Error(true),
			       //RooFit::Minos(1),
			       //RooFit::InitialHesse(true),
			       RooFit::Strategy(2),
			       RooFit::Range("signalRegionHX"));
			       //RooFit::Strategy(1));
//			       RooFit::PrintEvalErrors(-1));
				   //RooFit::PrintLevel(3);

		    hxResult->Print();
	  }
	}
	
	
	std::cout << std::endl << "Helicity Frame model after fit:" << std::endl;
	hx->Print();

	}

	current_out_hx = output_hx->mkdir(binName.str().c_str());
	current_out_cs = output_cs->mkdir(binName.str().c_str());
	

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
	delete thisBin;
      }
    }
  }


  cout<<entries[0][0]<<","<<entries[0][1]<<","<<entries[0][2]<<","<<entries[0][3]<<","<<entries[0][4]<<","<<entries[0][5]<<","<<entries[0][6]<<endl;
  cout<<entries[1][0]<<","<<entries[1][1]<<","<<entries[1][2]<<","<<entries[1][3]<<","<<entries[1][4]<<","<<entries[1][5]<<","<<entries[1][6]<<","<<entries[1][7]<<endl;

  output_hx->Write();
  output_hx->Close();

  output_cs->Write();
  output_cs->Close();
  delete output_hx;
  delete output_cs;
  }
  fitInput->Close();

  delete fitInput;

  delete data;
//  delete samples;
  return 0;
}

