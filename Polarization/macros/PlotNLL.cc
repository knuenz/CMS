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
#include "TColor.h"
#include "TF2.h"
#include "TAxis.h"
#include "TLatex.h"


void useNiceColorPalette( Int_t NCont = 255 ) {
 const Int_t NRGBs = 5;
 Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
 Double_t red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
 Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
 Double_t blue[NRGBs]= { 0.51, 1.00, 0.12, 0.00, 0.00 };
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 gStyle->SetNumberContours(NCont);
}


int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

  bool doprompt(true), dononprompt(true), dobkg(true),dopol(true),pereverr(false),newdata(false),fitHX(true),fitCS(true);

  gSystem->mkdir("fitRoots");
  gSystem->mkdir("Plots/LogLikelihood");
  gSystem->mkdir("Results/LogLikelihood");

//  gStyle->SetPalette(1);
//  gStyle->SetTitleFillColor(kWhite);

//  useNiceColorPalette(104);

  dononprompt=false;
  dobkg=false;
  pereverr=true;

/*  if(argc < 2) {
    std::cout << "Usage: extractMassShape /path/to/sample1.root /path/to/sample2.root ..." << std::endl;
    std::cout << "You must ensure that the MC samples are properly weighted so they can be combined." << std::endl;
    return 1;
  }
*/
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


//  TChain *samples = new TChain("data");
  TFile *fitInput_cs = new TFile("jPsiFitFinal_cs_ToyMC.root","UPDATE");
  TFile *fitInput_hx = new TFile("jPsiFitFinal_hx_ToyMC.root","UPDATE");

  RooDataSet *data = NULL;
  CompositeModelBuilder *hx, *cs;


  bool RealMC(false);



  Char_t *fileNameIn;
	if(RealMC) fileNameIn = "/scratch/knuenz/Polarization/RootInput/TTree_final_notrigger_MCprompt_Jpsi_Fall10_folded_.root";
	else fileNameIn = "/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo.root";
    TFile* fIn = new TFile(fileNameIn);
    TTree* dataTreesPR = (TTree*)fIn->Get("data");

    data = new RooDataSet("data","Supplied Data Prompt",varlist,Import(*dataTreesPR),WeightVar(MCweight));


  int iterations = 10;
  int generate_gt=0;
  bool plotNLL(true);

  char outputfilename2[200];
  sprintf(outputfilename2,"Results/LogLikelihood/parametersNLL.txt");
  FILE *outputFile2 = fopen(outputfilename2,"w");


  for (int iteration = 1; iteration < iterations+1; iteration++){

	  cout<<"GENERATION "<<iteration<<endl;


  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-4; ++yBin) {
 	  for(int ptBin = 1; ptBin < jpsi::kNbPTBins[yBin+1]-5; ++ptBin) {


      std::stringstream binName,cutString;
      if(RealMC) binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      else binName << "pt" << ptBin+1 << "_rapidity" << yBin+1 << "_gen" << iteration;

      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];
      std::cout << cutString.str() << std::endl;

	  cout<<"?"<<endl;


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

      bool writeProgress(true);
      if(writeProgress){
  	char outputfilename[200];
  	sprintf(outputfilename,"PlotNllToyMC_rap%d_pt%d.txt",yBin+1,ptBin+1);
  	FILE *outputFile = fopen(outputfilename,"w");
	fprintf(outputFile, "rapidity%d_pt%d, iteration %d\n",yBin+1,ptBin+1,iteration);
	fclose(outputFile);
      }

	TDirectory *current_in_cs,*current_in_hx;

//      if(current_in_cs = fitInput_cs->GetDirectory(binName.str().c_str()) && current_in_hx = fitInput_hx->GetDirectory(binName.str().c_str())) {


    	  current_in_cs = fitInput_cs->GetDirectory(binName.str().c_str());
    	  current_in_hx = fitInput_hx->GetDirectory(binName.str().c_str());

	RooAbsData *thisBin = data->reduce((cutString.str()).c_str());

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

	if(RealMC){
	cs->loadParameters(*current_in_cs);
	hx->loadParameters(*current_in_hx);
	}
	else{
	cs->loadParameters(*current_in_cs);
	hx->loadParameters(*current_in_hx);
	}

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
	hx->Print();
	cs->Print();

	RooAbsData *thisBinGenCS = cs->model()->generate(RooArgSet(costh_CS,phi_CS),jpsi::numEventsBin[yBin][ptBin]);//thisBin->tree()->GetEntries());
	RooAbsData *thisBinGenHX = hx->model()->generate(RooArgSet(costh_HX,phi_HX),jpsi::numEventsBin[yBin][ptBin]);//thisBin->tree()->GetEntries());

	cout<<"Generation Character: "<<thisBinGenCS->mean(costh_CS)<<" "<<thisBinGenCS->mean(phi_CS)<<endl;

//	cout<<cs->getPromptPolarizationModel()->lambdatheta()->getVal()<<endl;

//	cout<<hx->getPromptPolarizationModel()->lambdatheta()->getVal()<<endl;

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

	    if(plotNLL){



	    RooAddPdf* model;

	    char FilenameNLL[200];


		double nPrompt_fitvalue=cs->promptNorm()->getVal();
		double lambdatheta_fitvalue=cs->getPromptPolarizationModel()->lambdathetaCS()->getVal();
		double lambdaphi_fitvalue=cs->getPromptPolarizationModel()->lambdaphiCS()->getVal();
		double lambdathetaphi_fitvalue=cs->getPromptPolarizationModel()->lambdathetaphiCS()->getVal();

		double lambdatheta_fiterror=cs->getPromptPolarizationModel()->lambdathetaCS()->getError();
		double lambdaphi_fiterror=cs->getPromptPolarizationModel()->lambdaphiCS()->getError();
		double lambdathetaphi_fiterror=cs->getPromptPolarizationModel()->lambdathetaphiCS()->getError();

		double nsigma=3;

		double thetastep=nsigma*lambdatheta_fiterror;
		double phistep=nsigma*lambdaphi_fiterror;
		double thetaphistep=nsigma*lambdathetaphi_fiterror;

		cs->setVal("nPrompt",nPrompt_fitvalue); cs->fix("nPrompt");

		double lambdatheta_plotmean=lambdatheta_fitvalue;
		double lambdaphi_plotmean=lambdaphi_fitvalue;
		double lambdathetaphi_plotmean=lambdathetaphi_fitvalue;

		thetastep=0.1;
	    phistep=0.1;
	    thetaphistep=0.1;

	    lambdatheta_plotmean=0;
	    lambdaphi_plotmean=0;
	    lambdathetaphi_plotmean=0;


		cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",lambdatheta_plotmean);//lambdatheta_fitvalue);
		cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",lambdaphi_plotmean);//lambdaphi_fitvalue);
		cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",lambdathetaphi_plotmean);//lambdathetaphi_fitvalue);

		cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);

		RooRealVar* lambda_th=cs->getPromptPolarizationModel()->lambdathetaCS();
		RooRealVar* lambda_ph=cs->getPromptPolarizationModel()->lambdaphiCS();
		RooRealVar* lambda_thph=cs->getPromptPolarizationModel()->lambdathetaphiCS();
		RooRealVar* nPrompt=cs->promptNorm();

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",0);
	    cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",0);
	    cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",0);

	    model=cs->model();

	    gStyle->SetPadRightMargin(0.25);
		//gStyle->SetPalette(1);
		useNiceColorPalette();

	    TF2 *nllfunc_thph_fix_ = new TF2("nllfunc_thph_fix_","sin(x)+cos(y)-1500",1.6,3.14,2,4);
	    TF2 *nllfunc_ph_fix_ = new TF2("nllfunc_thph_fix_","sin(x)+sin(y)-1500",-2,2,-2,2);
	    TF2 *nllfunc_th_fix_ = new TF2("nllfunc_thph_fix_","sin(x)+cos(y)-1500",-5,5,-5,5);
		TCanvas* nll_canvas = new TCanvas("nll_canvas_thph","nll_canvas_thph",1200,1000);
		TPad* pad1 = new TPad("pad1","Gouraud shading",0,0,1,1,kWhite);
		pad1->Draw();
		pad1->cd();
		nllfunc_thph_fix_->GetHistogram()->SetXTitle("#lambda_{#theta CS}");
		nllfunc_thph_fix_->GetHistogram()->SetYTitle("#lambda_{#phi CS}");
		nllfunc_thph_fix_->SetNpx(1000);
		nllfunc_thph_fix_->Draw("colz");
	    sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_THPHconst_rapidity%d_pt%d_iter%d_try.png",yBin+1,ptBin+1,iteration);
//	    nll_canvas->SaveAs(FilenameNLL);
		nll_canvas->Close();
	    delete nll_canvas;


	    cout<<"?"<<endl;

	    RooAbsReal* nll_thph;
	    if(RealMC) nll_thph=model->createNLL(*thisBin);
	    else nll_thph=model->createNLL(*thisBinGenCS);

	    cout<<"?"<<endl;



	    fprintf(outputFile2, "Rapidity %d, pT %d, generation %d:\n",yBin+1,ptBin+1,iteration);
		fprintf(outputFile2, "NLL value = %f\n",nll_thph->getVal());


		cout<<"Start NLL thph"<<endl;
		TF1* nllfunc_thph_fix = new TF1;// nllfunc_thph_fix->SetNpx(10);
		nllfunc_thph_fix=nll_thph->asTF(RooArgList(*lambda_th,*lambda_ph),RooArgList(*lambda_thph,*nPrompt));
		nllfunc_thph_fix->SetRange(lambdatheta_plotmean-thetastep,lambdaphi_plotmean-phistep,lambdatheta_plotmean+thetastep,lambdaphi_plotmean+phistep);

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;



		cout<<"Start Canvas"<<endl;
	    cout<<"?"<<endl;

		TCanvas* nll_canvas_thph = new TCanvas("nll_canvas_thph","nll_canvas_thph",1200,1000);
	//	nllfunc_thph_fix->GetHistogram()->SetXTitle("#lambda_{#theta CS}");
	//	nllfunc_thph_fix->GetHistogram()->SetYTitle("#lambda_{#phi CS}");
		nllfunc_thph_fix->SetNpx(10); //nllfunc_thph_fix->SetNpy(10);
		nll_canvas_thph->SetFillColor(kWhite); gPad->SetFillColor(kWhite); nllfunc_thph_fix->Draw("cont1z");
		cout<<"Pad1 done"<<endl;
	    sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_THPHconst_rapidity%d_pt%d_iter%d.png",yBin+1,ptBin+1,iteration);
	    nll_canvas_thph->SaveAs(FilenameNLL);
		cout<<"Save1 done"<<endl;
	    nll_canvas_thph->Close();
	    delete nll_canvas_thph;
	    delete nllfunc_thph_fix;

	    cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
	    cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
	    cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",lambdatheta_plotmean);//lambdatheta_fitvalue);
	    cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",lambdaphi_plotmean);//lambdaphi_fitvalue);
	    cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",lambdathetaphi_plotmean);//lambdathetaphi_fitvalue);


	    cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",0);
	    cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",0);
	    cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",0);


	    cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);

	    model=cs->model();

	    RooAbsReal* nll_th;
	    if(RealMC) nll_th=model->createNLL(*thisBin);
	    else nll_th=model->createNLL(*thisBinGenCS);

		cout<<"Start NLL th"<<endl;
		TF1* nllfunc_th_fix = new TF1;
		nllfunc_th_fix=nll_th->asTF(RooArgList(*lambda_ph,*lambda_thph),RooArgList(*lambda_th,*nPrompt));
		nllfunc_th_fix->SetRange(lambdaphi_plotmean-phistep,lambdathetaphi_plotmean-thetaphistep,lambdaphi_plotmean+phistep,lambdathetaphi_plotmean+thetaphistep);

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
				cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
				cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

		TCanvas* nll_canvas_th = new TCanvas("nll_canvas_th","nll_canvas_th",1200,1000);
//		nllfunc_th_fix->GetHistogram()->SetXTitle("#lambda_{#phi CS}");
//		nllfunc_th_fix->GetHistogram()->SetYTitle("#lambda_{#thetaphi CS}");
		nllfunc_th_fix->SetNpx(100);
		nll_canvas_th->SetFillColor(kWhite); gPad->SetFillColor(kWhite);  nllfunc_th_fix->Draw("cont1z");
		cout<<"Pad2 done"<<endl;
		sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_THconst_rapidity%d_pt%d_iter%d.png",yBin+1,ptBin+1,iteration);
		nll_canvas_th->SaveAs(FilenameNLL);
		cout<<"Save2 done"<<endl;
		nll_canvas_th->Close();
		delete nll_canvas_th;
	    delete nllfunc_th_fix;

	    cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
	    		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
	    		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    		cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",lambdatheta_plotmean);//lambdatheta_fitvalue);
	    		cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",lambdaphi_plotmean);//lambdaphi_fitvalue);
	    		cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",lambdathetaphi_plotmean);//lambdathetaphi_fitvalue);

	    	    cs->getPromptPolarizationModel()->setVal("promptlambda_theta_CS",0);
	    	    cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",0);
	    	    cs->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_CS",0);

	    		cs->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);

	    	    model=cs->model();

	    RooAbsReal* nll_ph;
	    if(RealMC) nll_ph=model->createNLL(*thisBin);
	    else nll_ph=model->createNLL(*thisBinGenCS);

		cout<<"Start NLL ph"<<endl;
		TF1* nllfunc_ph_fix = new TF1;
		nllfunc_ph_fix=nll_ph->asTF(RooArgList(*lambda_th,*lambda_thph),RooArgList(*lambda_ph,*nPrompt));
		nllfunc_ph_fix->SetRange(lambdatheta_plotmean-thetastep,lambdathetaphi_plotmean-thetaphistep,lambdatheta_plotmean+thetastep,lambdathetaphi_plotmean+thetaphistep);

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
				cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
				cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

		TCanvas* nll_canvas_ph = new TCanvas("nll_canvas_ph","nll_canvas_ph",1200,1000);
//		nllfunc_ph_fix->GetHistogram()->SetXTitle("#lambda_{#theta CS}");
//		nllfunc_ph_fix->GetHistogram()->SetYTitle("#lambda_{#thetaphi CS}");
		nllfunc_ph_fix->SetNpx(1000);
		nll_canvas_ph->SetFillColor(kWhite); gPad->SetFillColor(kWhite); nllfunc_ph_fix->Draw("cont1z");
		cout<<"Pad3 done"<<endl;
		sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_PHconst_rapidity%d_pt%d_iter%d.png",yBin+1,ptBin+1,iteration);
		nll_canvas_ph->SaveAs(FilenameNLL);
		cout<<"Save3 done"<<endl;
	    nll_canvas_ph->Close();
	    delete nll_canvas_ph;
	    delete nllfunc_ph_fix;

	    cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
	    		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
	    		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    cs->getPromptPolarizationModel()->setVal("promptlambda_phi_CS",0.9);//lambdaphi_fitvalue);









	    cout<<"Canvas done"<<endl;

	  }

//	  }
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
	  }
	  else{
	  }









/*







	    if(plotNLL){

	    	cout<<"?"<<endl;

	    RooAddPdf* model;

	    char FilenameNLL[200];
    	cout<<"?"<<endl;


    	cout<<"?"<<endl;double nPrompt_fitvalue=hx->promptNorm()->getVal();
    	cout<<"?"<<endl;double lambdatheta_fitvalue=hx->getPromptPolarizationModel()->lambdatheta()->getVal();
    	cout<<"?"<<endl;double lambdaphi_fitvalue=hx->getPromptPolarizationModel()->lambdaphi()->getVal();
    	cout<<"?"<<endl;double lambdathetaphi_fitvalue=hx->getPromptPolarizationModel()->lambdathetaphi()->getVal();

    	cout<<"?"<<endl;double lambdatheta_fiterror=hx->getPromptPolarizationModel()->lambdatheta()->getError();
    	cout<<"?"<<endl;double lambdaphi_fiterror=hx->getPromptPolarizationModel()->lambdaphi()->getError();
    	cout<<"?"<<endl;double lambdathetaphi_fiterror=hx->getPromptPolarizationModel()->lambdathetaphi()->getError();
    	cout<<"?"<<endl;

		double nsigma=3;

		double thetastep=nsigma*lambdatheta_fiterror;
		double phistep=nsigma*lambdaphi_fiterror;
		double thetaphistep=nsigma*lambdathetaphi_fiterror;

		hx->setVal("nPrompt",nPrompt_fitvalue); hx->fix("nPrompt");

		double lambdatheta_plotmean=lambdatheta_fitvalue;
		double lambdaphi_plotmean=lambdaphi_fitvalue;
		double lambdathetaphi_plotmean=lambdathetaphi_fitvalue;

		thetastep=0.1;
	    phistep=0.1;
	    thetaphistep=0.1;

	    lambdatheta_plotmean=0;
	    lambdaphi_plotmean=0;
	    lambdathetaphi_plotmean=0;


		hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",lambdatheta_plotmean);//lambdatheta_fitvalue);
		hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",lambdaphi_plotmean);//lambdaphi_fitvalue);
		hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",lambdathetaphi_plotmean);//lambdathetaphi_fitvalue);

		hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX);

		RooRealVar* lambda_th=hx->getPromptPolarizationModel()->lambdatheta();
		RooRealVar* lambda_ph=hx->getPromptPolarizationModel()->lambdaphi();
		RooRealVar* lambda_thph=hx->getPromptPolarizationModel()->lambdathetaphi();
		RooRealVar* nPrompt=hx->promptNorm();

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",0);
	    hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",0);
	    hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",0);
    	cout<<"?"<<endl;

	    model=hx->model();

	    gStyle->SetPadRightMargin(0.25);
		//gStyle->SetPalette(1);
		useNiceColorPalette();

	    TF2 *nllfunc_thph_fix_ = new TF2("nllfunc_thph_fix_","sin(x)+cos(y)-1500",1.6,3.14,2,4);
	    TF2 *nllfunc_ph_fix_ = new TF2("nllfunc_thph_fix_","sin(x)+sin(y)-1500",-2,2,-2,2);
	    TF2 *nllfunc_th_fix_ = new TF2("nllfunc_thph_fix_","sin(x)+cos(y)-1500",-5,5,-5,5);
		TCanvas* nll_canvas = new TCanvas("nll_canvas_thph","nll_canvas_thph",1200,1000);
		TPad* pad1 = new TPad("pad1","Gouraud shading",0,0,1,1,kWhite);
		pad1->Draw();
		pad1->cd();
		nllfunc_thph_fix_->GetHistogram()->SetXTitle("#lambda_{#theta HX}");
		nllfunc_thph_fix_->GetHistogram()->SetYTitle("#lambda_{#phi HX}");
		nllfunc_thph_fix_->SetNpx(1000);
		nllfunc_thph_fix_->Draw("colz");
	    sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_THPHconst_rapidity%d_pt%d_iter%d_try.png",yBin+1,ptBin+1,iteration);
//	    nll_canvas->SaveAs(FilenameNLL);
		nll_canvas->Close();
	    delete nll_canvas;


	    cout<<"?"<<endl;

	    RooAbsReal* nll_thph;
	    if(RealMC) nll_thph=model->createNLL(*thisBin);
	    else nll_thph=model->createNLL(*thisBinGenHX);

	    cout<<"?"<<endl;



	    fprintf(outputFile2, "Rapidity %d, pT %d, generation %d:\n",yBin+1,ptBin+1,iteration);
		fprintf(outputFile2, "NLL value = %f\n",nll_thph->getVal());


		cout<<"Start NLL thph"<<endl;
		TF1* nllfunc_thph_fix = new TF1;// nllfunc_thph_fix->SetNpx(10);
		nllfunc_thph_fix=nll_thph->asTF(RooArgList(*lambda_th,*lambda_ph),RooArgList(*lambda_thph,*nPrompt));
		nllfunc_thph_fix->SetRange(lambdatheta_plotmean-thetastep,lambdaphi_plotmean-phistep,lambdatheta_plotmean+thetastep,lambdaphi_plotmean+phistep);

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;



		cout<<"Start Canvas"<<endl;
	    cout<<"?"<<endl;

		TCanvas* nll_canvas_thph = new TCanvas("nll_canvas_thph","nll_canvas_thph",1200,1000);
	//	nllfunc_thph_fix->GetHistogram()->SetXTitle("#lambda_{#theta HX}");
	//	nllfunc_thph_fix->GetHistogram()->SetYTitle("#lambda_{#phi HX}");
		nllfunc_thph_fix->SetNpx(10); //nllfunc_thph_fix->SetNpy(10);
		nll_canvas_thph->SetFillColor(kWhite); gPad->SetFillColor(kWhite); nllfunc_thph_fix_->Draw("cont1z");
		cout<<"Pad1 done"<<endl;
	    sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_THPHconst_rapidity%d_pt%d_iter%d.png",yBin+1,ptBin+1,iteration);
	    nll_canvas_thph->SaveAs(FilenameNLL);
		cout<<"Save1 done"<<endl;
	    nll_canvas_thph->Close();
	    delete nll_canvas_thph;
	    delete nllfunc_thph_fix;

	    cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
	    cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
	    cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",lambdatheta_plotmean);//lambdatheta_fitvalue);
	    hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",lambdaphi_plotmean);//lambdaphi_fitvalue);
	    hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",lambdathetaphi_plotmean);//lambdathetaphi_fitvalue);


	    hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",0);
	    hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",0);
	    hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",0);


	    hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX);

	    model=hx->model();

	    RooAbsReal* nll_th;
	    if(RealMC) nll_th=model->createNLL(*thisBin);
	    else nll_th=model->createNLL(*thisBinGenHX);

		cout<<"Start NLL th"<<endl;
		TF1* nllfunc_th_fix = new TF1;
		nllfunc_th_fix=nll_th->asTF(RooArgList(*lambda_ph,*lambda_thph),RooArgList(*lambda_th,*nPrompt));
		nllfunc_th_fix->SetRange(lambdaphi_plotmean-phistep,lambdathetaphi_plotmean-thetaphistep,lambdaphi_plotmean+phistep,lambdathetaphi_plotmean+thetaphistep);

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
				cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
				cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

		TCanvas* nll_canvas_th = new TCanvas("nll_canvas_th","nll_canvas_th",1200,1000);
//		nllfunc_th_fix->GetHistogram()->SetXTitle("#lambda_{#phi HX}");
//		nllfunc_th_fix->GetHistogram()->SetYTitle("#lambda_{#thetaphi HX}");
		nllfunc_th_fix->SetNpx(100);
		nll_canvas_th->SetFillColor(kWhite); gPad->SetFillColor(kWhite);  nllfunc_th_fix_->Draw("cont1z");
		cout<<"Pad2 done"<<endl;
		sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_THconst_rapidity%d_pt%d_iter%d.png",yBin+1,ptBin+1,iteration);
		nll_canvas_th->SaveAs(FilenameNLL);
		cout<<"Save2 done"<<endl;
		nll_canvas_th->Close();
		delete nll_canvas_th;
	    delete nllfunc_th_fix;

	    cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
	    		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
	    		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    		hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",lambdatheta_plotmean);//lambdatheta_fitvalue);
	    		hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",lambdaphi_plotmean);//lambdaphi_fitvalue);
	    		hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",lambdathetaphi_plotmean);//lambdathetaphi_fitvalue);

	    	    hx->getPromptPolarizationModel()->setVal("promptlambda_theta_HX",0);
	    	    hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",0);
	    	    hx->getPromptPolarizationModel()->setVal("promptlambda_thetaphi_HX",0);

	    		hx->initModel(JpsiMass,Jpsict,JpsictErr,costh_HX,phi_HX);

	    	    model=hx->model();

	    RooAbsReal* nll_ph;
	    if(RealMC) nll_ph=model->createNLL(*thisBin);
	    else nll_ph=model->createNLL(*thisBinGenHX);

		cout<<"Start NLL ph"<<endl;
		TF1* nllfunc_ph_fix = new TF1;
		nllfunc_ph_fix=nll_ph->asTF(RooArgList(*lambda_th,*lambda_thph),RooArgList(*lambda_ph,*nPrompt));
		nllfunc_ph_fix->SetRange(lambdatheta_plotmean-thetastep,lambdathetaphi_plotmean-thetaphistep,lambdatheta_plotmean+thetastep,lambdathetaphi_plotmean+thetaphistep);

		cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
				cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
				cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

		TCanvas* nll_canvas_ph = new TCanvas("nll_canvas_ph","nll_canvas_ph",1200,1000);
//		nllfunc_ph_fix->GetHistogram()->SetXTitle("#lambda_{#theta HX}");
//		nllfunc_ph_fix->GetHistogram()->SetYTitle("#lambda_{#thetaphi HX}");
		nllfunc_ph_fix->SetNpx(1000);
		nll_canvas_ph->SetFillColor(kWhite); gPad->SetFillColor(kWhite); nllfunc_ph_fix_->Draw("cont1z");
		cout<<"Pad3 done"<<endl;
		sprintf(FilenameNLL,"Plots/LogLikelihood/NLL_RealMC_PHconst_rapidity%d_pt%d_iter%d.png",yBin+1,ptBin+1,iteration);
		nll_canvas_ph->SaveAs(FilenameNLL);
		cout<<"Save3 done"<<endl;
	    nll_canvas_ph->Close();
	    delete nll_canvas_ph;
	    delete nllfunc_ph_fix;

	    cout<<"lambda_theta_val"<<lambda_th->getVal()<<endl;
	    		cout<<"lambda_phi_val"<<lambda_ph->getVal()<<endl;
	    		cout<<"lambda_thetaphi_val"<<lambda_thph->getVal()<<endl;

	    hx->getPromptPolarizationModel()->setVal("promptlambda_phi_HX",0.9);//lambdaphi_fitvalue);









	    cout<<"Canvas done"<<endl;

	  }
*/








	std::cout << std::endl << "Helicity Frame model after fit:" << std::endl;
	hx->Print();

	}

	delete cs;
	delete hx;
//	delete thisBin;
    // if von 220  }
    }
  }

//  output_hx->Write();
//  output_hx->Close();

//  output_cs->Write();
//  output_cs->Close();
//  delete output_hx;
//  delete output_cs;
  }
//  fitInput->Close();
  fclose(outputFile2);
  delete fitInput_cs;
  delete fitInput_hx;

  delete fIn;

  delete data;
//  delete samples;
  return 0;
}

