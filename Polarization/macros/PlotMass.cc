#include <iostream>
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
#include "RooCBShape.h"



//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TFrame.h"



int main(int argc, char**argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

  bool dopull(true);

for( int i=0;i < argc; ++i ) {
     if(std::string(argv[i]).find("--noPull") != std::string::npos) dopull = false;
   }

	cout<<"Options: --noPull -> no Pull distributions are plotted"<<endl;


  gSystem->mkdir("Plots/Mass");
  gStyle->SetTitleFillColor(kWhite);

  gStyle->SetPalette(1,0);
/*  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.15);

  gStyle->SetTickLength(-0.02, "xyz");
  gStyle->SetLabelOffset(0.02, "x");
  gStyle->SetLabelOffset(0.02, "y");
  gStyle->SetTitleOffset(1.3, "x");
  gStyle->SetTitleOffset(1.4, "y");
*/

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1.,2.5);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,1000);
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",0.5,1.5);

  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
  //varlist.add(HLT_Mu0_TkMu0_Jpsi);
  varlist.add(MCweight);


  Char_t *fileNameInput = "jPsiFit_Bins_PR_FSR.root";
  TFile* fInput = new TFile(fileNameInput);

  Char_t *fileNameIn = "/scratch/knuenz/Polarization/RootInput/TTree_weightTry.root";
  TFile* fIn = new TFile(fileNameIn);
  TTree* dataTreesPR = (TTree*)fIn->Get("data");

  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",varlist,Import(*dataTreesPR),WeightVar(MCweight));

  RooDataSet *dataPRbin;

  MassModel* massinp_;

  RooPlot* FSRMassFrame;

  int MassBins = 50;

  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-4; ++yBin) {
	  for(int ptBin = 3; ptBin < jpsi::kNbPTBins[yBin+1]-3; ++ptBin) {


      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char reduceNP[200];
	  sprintf(reduceNP,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reduceBK[200];
	  sprintf(reduceBK,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiMass < %f | JpsiMass > %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow,jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh);

	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);

	  JpsiMass.setBins(50);

	  char DirectoryPathMass[200];
	  sprintf(DirectoryPathMass,"/pt%d_rapidity%d/MassModel",ptBin+1,yBin+1);

	  TDirectory *InputDirectoryMass = (TDirectory*)fInput->Get(DirectoryPathMass);


	  massinp_ = new MassModel();
	  if(fInput->GetDirectory(DirectoryPathMass)!=NULL){
	  massinp_->loadParameters(*InputDirectoryMass);
	  }
	  else cout<<"TDirectory "<<DirectoryPathMass<<" does not exist in File "<<fileNameInput<<endl;
	  massinp_->initModels(JpsiMass);



  FSRMassFrame = new RooPlot;
  FSRMassFrame = JpsiMass.frame(Bins(50)) ;
  dataPRbin->plotOn(FSRMassFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
  massinp_->prompt()->plotOn(FSRMassFrame,Normalization(1.0), LineWidth(2));
  double chi2_FSRMassFrame=FSRMassFrame->chiSquare();
  FSRMassFrame->SetMinimum(1);
  char FSRMassTitle[200];
  sprintf(FSRMassTitle,"MC MassFit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_FSRMassFrame);
  FSRMassFrame->SetTitle(FSRMassTitle);


  TCanvas* FSRMassCanvas = new TCanvas("FSRMassCanvas","FSRMassCanvas",1400, 600);
  FSRMassCanvas->Divide(2);  FSRMassCanvas->SetFillColor(kWhite);
  FSRMassCanvas->cd(1)->SetLogy(1);
  FSRMassCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); FSRMassFrame->GetYaxis()->SetTitleOffset(2) ; FSRMassFrame->Draw();
  FSRMassCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); FSRMassFrame->GetYaxis()->SetTitleOffset(2) ; FSRMassFrame->Draw();

  char FilenamePRMass[200];
  sprintf(FilenamePRMass,"Plots/Mass/massPrompt_rapidity%d_pt%d.png",yBin+1,ptBin+1);


  FSRMassCanvas->SaveAs(FilenamePRMass);

  FSRMassCanvas->Close();


////////////////////////////////// 	"Pull" ///////////////////////////////////////////////////////////////////////////////////

if(dopull){

  TH1 * PromptModelHist =  massinp_->prompt()->createHistogram("JpsiMass",MassBins) ;
  TH1 * PromptHist = dataPRbin->createHistogram("JpsiMass",MassBins);


  double deltaPrompt[MassBins];
  double errdeltaPrompt[MassBins];
  double FitValue;
  double deltaPromptMax=0;
  double massMean[MassBins];
  double errmassMean[MassBins];

  for(int BinRun = 1; BinRun < MassBins+1; ++BinRun) {

	  FitValue=PromptModelHist->GetBinContent(BinRun)*dataPRbin->sumEntries()/PromptModelHist->GetSumOfWeights();

	  if(PromptHist->GetBinError(BinRun)==0) deltaPrompt[BinRun-1]=100000;
	  else deltaPrompt[BinRun-1]=((PromptHist->GetBinContent(BinRun)-FitValue)/PromptHist->GetBinError(BinRun));


	  errdeltaPrompt[BinRun-1]=0;//PromptHist->GetBinError(BinRun);

	  if (fabs(deltaPrompt[BinRun-1])+fabs(errdeltaPrompt[BinRun-1]) > deltaPromptMax && fabs(deltaPrompt[BinRun-1])+fabs(errdeltaPrompt[BinRun-1]) < 99999)
	   	    		  deltaPromptMax=fabs(deltaPrompt[BinRun-1])+fabs(errdeltaPrompt[BinRun-1]);

	  massMean[BinRun-1]=2.7+0.8/MassBins*BinRun-0.4/MassBins;
	   	    errmassMean[BinRun-1]=0;
  }








  	    deltaPromptMax*=1.05;

	 TCanvas *MassCanvas = new TCanvas("MassCanvas","",1000,1000);

	  TPad* MassPad = new TPad("MassPad","",0,0,1,1);
	  MassPad->SetFillColor(0);

	  TPad* MassPullPad = new TPad("MassPullPad","",0,0,1,0.3);
	  MassPullPad->SetFillColor(0);

	  MassPad->Draw();
	  MassPad->cd();
	  MassPad->SetBottomMargin(0.4);
	  gPad->SetLeftMargin(0.15) ; gPad->SetFillColor(kWhite); FSRMassFrame->GetYaxis()->SetTitleOffset(1.5) ;
	  FSRMassFrame->Draw();

	  MassPullPad->Draw();
	  MassPullPad->cd();
	  MassPullPad->SetGrid();
	  gPad->SetLeftMargin(0.15) ; gPad->SetFillColor(kWhite); FSRMassFrame->GetYaxis()->SetTitleOffset(1.5);

	  TH1F *MassPullhistoHisto = MassPullPad->DrawFrame(2.7,-deltaPromptMax,3.5,deltaPromptMax);
	  MassPullhistoHisto->SetYTitle("(Data - Fit)/Sigma");
	  MassPullhistoHisto->GetYaxis()->SetTitleOffset(1.5);

	  TGraphErrors *MassPullGraph = new TGraphErrors(MassBins,massMean,deltaPrompt,errmassMean,errdeltaPrompt);
	  MassPullGraph->SetMarkerColor(kBlack);
	  MassPullGraph->SetMarkerStyle(10);
	  MassPullGraph->Draw("P");

	  char Filename[50];
	  sprintf(Filename,"Plots/Mass/MassPull_rapidity%d_pt%d.png",yBin+1,ptBin+1);

   MassCanvas->SaveAs(Filename);
   MassCanvas->Close();

}

  delete FSRMassFrame;
  delete massinp_;


	  }
	}

  delete dataTreesPR;
  delete dataPR;
  delete fInput;
  delete fIn;



  return 0;
}
