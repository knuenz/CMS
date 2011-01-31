#include <iostream>
//#include <string>
#include "TString.h"

#include <sstream>
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
#include "RooExponential.h"
#include "RooAbsPdf.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"

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

  gSystem->mkdir("Plots/Lifetime");
  gStyle->SetTitleFillColor(kWhite);

  bool pereverr(false);
  bool noPR(false);
  bool noNP(false);
  bool noBK(false);
  bool noPull(false);
  bool newTTree(false);

   for(int i=0;i<argc; ++i) {
     if(std::string(argv[i]).find("--perEventErrors") != std::string::npos) pereverr = true;
     if(std::string(argv[i]).find("--noPR") != std::string::npos) noPR = true;
     if(std::string(argv[i]).find("--noNP") != std::string::npos) noNP = true;
     if(std::string(argv[i]).find("--noBK") != std::string::npos) noBK = true;
     if(std::string(argv[i]).find("--noPull") != std::string::npos) noPull = true;
     if(std::string(argv[i]).find("--newTTree") != std::string::npos) newTTree = true;
   }

// CHANGES: RooHist und RooCurve - dings... Option-output terminal, pulls erweitern auf NP,BK

  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,360);
  RooRealVar MCType_idx("MCType_idx","MCType_idx",0,2);//0=PR,1=NP,2=BK
  RooRealVar JpsiType_idx("JpsiType_idx","JpsiType_idx",0,2.5);//0=GG,1=GT,2=TT
  RooRealVar MCweight("MCweight","MCweight",0,500);
  RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
  RooRealVar HLT_Mu0_TkMu0_Jpsi("HLT_Mu0_TkMu0_Jpsi","Passes HLT_Mu0_TkMu0_Jpsi Trigger",0.5,1.5);

  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX,JpsictErr);
  varlist.add(MCweight);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
if(newTTree)
  varlist.add(HLT_Mu0_TkMu0_Jpsi);

//  JpsictErr.setBins(100);

  int ctauBins=50;


  Char_t *fileNameInPR = "/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root";
  TFile* fInPR = new TFile(fileNameInPR);
  TTree* dataTreesPR = (TTree*)fInPR->Get("data");
  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",varlist,Import(*dataTreesPR),WeightVar(MCweight));

  Char_t *fileNameInNP = "/scratch/knuenz/Polarization/RootInput/TTree_NPerr.root";
  TFile* fInNP = new TFile(fileNameInNP);
  TTree* dataTreesNP = (TTree*)fInNP->Get("data");
  RooDataSet *dataNP = new RooDataSet("dataNP","Supplied Data Non Prompt",varlist,Import(*dataTreesNP),WeightVar(MCweight));

  Char_t *fileNameInBK = "/scratch/knuenz/Polarization/RootInput/TTree_Data29Sept.root";
  TFile* fInBK = new TFile(fileNameInBK);
  TTree* dataTreesBK = (TTree*)fInBK->Get("data");
  RooDataSet *dataBK = new RooDataSet("dataBK","Supplied Real Data",varlist,Import(*dataTreesBK),WeightVar(MCweight));


   Char_t *fileNameInput = "jPsiFit_Bins_PR.root";
   TFile* fInput = new TFile(fileNameInput);

 
  RooDataSet *dataNPbin, *dataPRbin, *dataBKbin, *dataREALbin;

  LifetimeModel* lifetinp_;

  RooPlot* BackgroundFrame;
  RooPlot* NonPromptFrame;
  RooPlot* PromptFrame;

	char outputfilename[200];
	sprintf(outputfilename,"Results/LifetimeParameterResults.txt");
	printf("output filename is: %s\n", outputfilename);
	FILE *outputFile = fopen(outputfilename,"w");

 RooDataHist *dataPRbinHist;
  RooDataHist *dataNPbinHist;
  RooDataHist *dataBKbinHist;




//  TFile* fInData = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root");

  RooRealVar *ctMeanPrompt1,*ctSigmaPrompt1;
  RooRealVar *ssTauNonPrompt;
  RooRealVar *ssTauBkg,*dsTauBkg,*fTauBkg,*ssCoefBkg,*fCoefBkg;





  for(int yBin = 1; yBin < jpsi::kNbRapForPTBins-3; ++yBin) {
	  for(int ptBin = 0; ptBin < jpsi::kNbPTBins[yBin+1]-7; ++ptBin) {

	  fprintf(outputFile, "pt%d_rapidity%d\n",ptBin+1,yBin+1);

      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char reduceNP[200];
	  sprintf(reduceNP,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);// && JpsictErr > 0.0004 && JpsictErr < 0.999999
	  char reduceBK[200];
	  sprintf(reduceBK,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiMass < %f | JpsiMass > %f",
		  jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],
		  jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow,
		  jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh);

	  cout<<reducePR<<endl;

	  dataNPbin = (RooDataSet*)dataNP->reduce(reduceNP);
	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);
	  dataBKbin = (RooDataSet*)dataBK->reduce(reduceBK);

	  dataNPbinHist = new RooDataHist("dataNPbinHist","dataNPbinHist",RooArgSet(Jpsict,JpsictErr),*dataNPbin);
	  dataPRbinHist = new RooDataHist("dataPRbinHist","dataPRbinHist",RooArgSet(Jpsict,JpsictErr),*dataPRbin);
	  dataBKbinHist = new RooDataHist("dataBKbinHist","dataBKbinHist",RooArgSet(Jpsict,JpsictErr),*dataBKbin);

	  char DirectoryPath[200];
	 	  sprintf(DirectoryPath,"/pt%d_rapidity%d/LifetimeModel",ptBin+1,yBin+1);

	 	  TDirectory *InputDirectory = (TDirectory*)fInput->Get(DirectoryPath);


	    	TTree *treeData = (TTree*)dataPRbin->tree();//fInData->Get("data");

	    	TTree *fChain;
	    	fChain=treeData;

	    	Long64_t nentries = fChain->GetEntries();
	    	cout<<nentries<<endl;
	    	Long64_t nb = 0;

	    	Double_t        JpsictErr_;

	    	TBranch        *b_JpsictErr;   //!

	    	fChain->SetBranchAddress("JpsictErr", &JpsictErr_, &b_JpsictErr);

	    	RooRealVar* mean_ = (RooRealVar*)InputDirectory->Get("ctMeanPrompt1");
	    	double mean=mean_->getVal();
	    	RooRealVar* sigma_ = (RooRealVar*)InputDirectory->Get("ctSigmaPrompt1");
	    	double sigma=sigma_->getVal();

//	    	cout<<mean<<endl;
//	    	cout<<sigma<<endl;

	    	RooExtendPdf**egauss=new RooExtendPdf*[nentries];
	    	RooGaussian**gauss=new RooGaussian*[nentries];
	    	RooRealVar**sigmaErr_=new RooRealVar*[nentries];
	    	RooRealVar**gaussWeight=new RooRealVar*[nentries];


	    	RooArgList gaussList;
	    	RooArgList coefList;

//	    	RooRealVar mean__("mean__","mean__",mean);
//	    	RooRealVar sigma__("sigma__","sigma__",sigma);

	        for (Long64_t jentry=0; jentry<nentries;jentry++) {


	       	   nb = fChain->GetEntry(jentry);

				//   cout<<JpsictErr_<<endl;

				   double sigmaErr=JpsictErr_*sigma;
				   double weightErr=1/JpsictErr_;

	       	   gaussWeight[jentry]=new RooRealVar(Form("GaussWeight_%i",jentry),Form("GaussWeight_%i",jentry),weightErr);
	       	   sigmaErr_[jentry]=new RooRealVar(Form("sigmaErr_%i",jentry),Form("sigmaErr_%i",jentry),sigmaErr);
	       	   gauss[jentry]=new RooGaussian(Form("gauss_%i",jentry),Form("gauss_%i",jentry),Jpsict,*mean_,*sigmaErr_[jentry]);//*ErrorArray[jentry]);Form("gauss_%i",jentry),Form("gauss_%i",jentry)
	       	   egauss[jentry]=new RooExtendPdf(Form("egauss_%i",jentry),Form("egauss_%i",jentry),*gauss[jentry],*gaussWeight[jentry]);

	       	   gaussList.add(*egauss[jentry]);
	       	   coefList.add(*gaussWeight[jentry]);

	        }

	        RooAddPdf gaussSum("gaussSum","gaussSum",gaussList);//,coefList);


	  	  lifetinp_ = new LifetimeModel();
	  	  lifetinp_->loadParameters(*InputDirectory);
	  	  lifetinp_->initModels(Jpsict,JpsictErr);

	  CompositeModelBuilder *LifetimeModel = new CompositeModelBuilder("","");
		  LifetimeModel->setUseLifetime(true);
		  LifetimeModel->setUseMass(false);
		  LifetimeModel->setUsePol(false);
		  if(fInput->GetDirectory(DirectoryPath)!=NULL){
		  LifetimeModel->loadParameters(*InputDirectory);
		  }
		  else cout<<"TDirectory "<<DirectoryPath<<" does not exist in File "<<fInput<<endl;

		  if(pereverr)
		    LifetimeModel->initModel(JpsiMass,Jpsict,JpsictErr,costh_CS,phi_CS);
		  else
		    LifetimeModel->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

//		  RooDataSet* generateSet = LifetimeModel->getPromptModel()->generate(Jpsict,10000);


if(!noPR){
	  fprintf(outputFile, "Prompt Parameters:\n");

	  PromptFrame = new RooPlot;
	  PromptFrame = Jpsict.frame(Bins(ctauBins)) ;
	  PromptFrame->setPadFactor(0.5);
	  double chi2_PromptFrame;
	  dataPRbin->plotOn(PromptFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr){
	  lifetinp_->prompt()->plotOn(PromptFrame, RooFit::ProjWData(RooArgSet(JpsictErr),*dataPRbin), Normalization(1.0), LineWidth(3),LineColor(kBlue));//
	  chi2_PromptFrame=PromptFrame->chiSquare();
      }
	  else{
	  lifetinp_->prompt()->plotOn(PromptFrame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
	  chi2_PromptFrame=PromptFrame->chiSquare();
//	  lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss1()), RooFit::LineStyle(kSolid),RooFit::LineColor(kBlack),LineWidth(2));
//	  lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss2()), RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),LineWidth(2));
//	  lifetinp_->prompt()->plotOn(PromptFrame, Components(*lifetinp_->promptGauss3()), RooFit::LineStyle(kSolid),RooFit::LineColor(kGreen),LineWidth(2));
	  }
//	  gaussSum.plotOn(PromptFrame, Normalization(1.0), LineWidth(3),LineColor(kRed));
	  PromptFrame->SetMinimum(1e-3);
	  char PromptTitle[200];
	  sprintf(PromptTitle,"MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_PromptFrame);
	  PromptFrame->SetTitle(PromptTitle);

	 ctMeanPrompt1 = (RooRealVar*)InputDirectory->Get("ctMeanPrompt1");
     fprintf(outputFile, "ctMeanPrompt1 = %1.5f +- %1.5f\n",ctMeanPrompt1->getVal(),ctMeanPrompt1->getError());
     ctSigmaPrompt1 = (RooRealVar*)InputDirectory->Get("ctSigmaPrompt1");
     fprintf(outputFile, "ctSigmaPrompt1 = %1.5f +- %1.5f\n",ctSigmaPrompt1->getVal(),ctSigmaPrompt1->getError());
}

if(!noNP){
	  fprintf(outputFile, "Prompt Parameters:\n");

	  NonPromptFrame = new RooPlot;
	  NonPromptFrame = Jpsict.frame(Bins(ctauBins)) ;
	  dataNPbin->plotOn(NonPromptFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr)
	  lifetinp_->nonPrompt()->plotOn(NonPromptFrame, RooFit::ProjWData(RooArgSet(JpsictErr),*dataNPbin),Normalization(1.0), LineWidth(3));
	  else
	  lifetinp_->nonPrompt()->plotOn(NonPromptFrame,Normalization(1.0), LineWidth(3));
	  double chi2_NonPromptFrame=NonPromptFrame->chiSquare();
	  NonPromptFrame->SetMinimum(1e-1);
	  char NonPromptTitle[200];
	  sprintf(NonPromptTitle,"MC Non Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_NonPromptFrame);
	  NonPromptFrame->SetTitle(NonPromptTitle);

	  ssTauNonPrompt = (RooRealVar*)InputDirectory->Get("ssTauNonPrompt");
	  fprintf(outputFile, "ssTauNonPrompt = %1.5f +- %1.5f\n",ssTauNonPrompt->getVal(),ssTauNonPrompt->getError());

}

if(!noBK){
	  fprintf(outputFile, "Prompt Parameters:\n");

	  BackgroundFrame = new RooPlot;
	  BackgroundFrame = Jpsict.frame(Bins(ctauBins)) ;
	  dataBKbin->plotOn(BackgroundFrame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	  if(pereverr)
	  lifetinp_->background()->plotOn(BackgroundFrame,Normalization(1.0), LineWidth(3));//,RooFit::ProjWData(RooArgSet(JpsictErr),*dataBKbinHist)
	  else
	  lifetinp_->background()->plotOn(BackgroundFrame,Normalization(1.0), LineWidth(3));
	  double chi2_BackgroundFrame=BackgroundFrame->chiSquare();
//	  lifetinp_->background()->plotOn(BackgroundFrame, Components(*lifetinp_->backgroundDecayss()), RooFit::LineStyle(kSolid),RooFit::LineColor(kBlack),LineWidth(2));
//	  lifetinp_->background()->plotOn(BackgroundFrame, Components(*lifetinp_->backgroundDecayf()), RooFit::LineStyle(kSolid),RooFit::LineColor(kRed),LineWidth(2));
//	  lifetinp_->background()->plotOn(BackgroundFrame, Components(*lifetinp_->backgroundDecayds()), RooFit::LineStyle(kSolid),RooFit::LineColor(kGreen),LineWidth(2));
	  BackgroundFrame->SetMinimum(1e-1);
	  char BackgroundTitle[200];
	  sprintf(BackgroundTitle,"Data Background Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_BackgroundFrame);
	  BackgroundFrame->SetTitle(BackgroundTitle);

	  ssTauBkg = (RooRealVar*)InputDirectory->Get("ssTauBkg");
	  fprintf(outputFile, "ssTauBkg = %1.5f +- %1.5f\n",ssTauBkg->getVal(),ssTauBkg->getError());
	  fTauBkg = (RooRealVar*)InputDirectory->Get("fTauBkg");
	  fprintf(outputFile, "fTauBkg = %1.5f +- %1.5f\n",fTauBkg->getVal(),fTauBkg->getError());
	  dsTauBkg = (RooRealVar*)InputDirectory->Get("dsTauBkg");
	  fprintf(outputFile, "dsTauBkg = %1.5f +- %1.5f\n",dsTauBkg->getVal(),dsTauBkg->getError());
	  ssCoefBkg = (RooRealVar*)InputDirectory->Get("ssCoefBkg");
	  fprintf(outputFile, "ssCoefBkg = %1.5f +- %1.5f\n",ssCoefBkg->getVal(),ssCoefBkg->getError());
	  fCoefBkg = (RooRealVar*)InputDirectory->Get("fCoefBkg");
	  fprintf(outputFile, "fCoefBkg = %1.5f +- %1.5f\n",fCoefBkg->getVal(),fCoefBkg->getError());
}

  TH1* legendBlue = dataPR->createHistogram("legendBlue",JpsiMass,Binning(100)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(3) ;
  TH1* legendRed = dataPR->createHistogram("legendRed",JpsiMass,Binning(100)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(kSolid) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = dataPR->createHistogram("legendBlack",JpsiMass,Binning(100)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = dataPR->createHistogram("legendGreen",JpsiMass,Binning(100)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;

  TLegend* PromptLegend=new TLegend(0.6,0.6,0.9,0.9);
  PromptLegend->SetFillColor(kWhite);
  PromptLegend->SetTextFont(72);
  PromptLegend->SetTextSize(0.025);
  PromptLegend->AddEntry(legendBlue,"MC Prompt fit","l");
  if(!pereverr){
  PromptLegend->AddEntry(legendBlack,"Gauss1","l");
  PromptLegend->AddEntry(legendRed,"Gauss2","l");
  PromptLegend->AddEntry(legendGreen,"Gauss3","l");
  }
  TLegend* NonPromptLegend=new TLegend(0.6,0.6,0.9,0.9);
  NonPromptLegend->SetFillColor(kWhite);
  NonPromptLegend->SetTextFont(72);
  NonPromptLegend->SetTextSize(0.025);
  NonPromptLegend->AddEntry(legendBlue,"MC Non Prompt fit","l");
//  NonPromptLegend->AddEntry(legendBlack,"SingleSided Decay","l");
//  NonPromptLegend->AddEntry(legendRed,"Flipped Decay","l");

  TLegend* BackgroundLegend=new TLegend(0.6,0.6,0.9,0.9);
  BackgroundLegend->SetFillColor(kWhite);
  BackgroundLegend->SetTextFont(72);
  BackgroundLegend->SetTextSize(0.025);
  BackgroundLegend->AddEntry(legendBlue,"Data background fit","l");
  BackgroundLegend->AddEntry(legendBlack,"SingleSided Decay","l");
  BackgroundLegend->AddEntry(legendRed,"Flipped Decay","l");
  BackgroundLegend->AddEntry(legendGreen,"DoubleSided Decay","l");

  char FilenameNP[200];
  sprintf(FilenameNP,"Plots/Lifetime/ctNonPrompt_rapidity%d_pt%d.png",yBin+1,ptBin+1);
  char FilenamePR[200];
  sprintf(FilenamePR,"Plots/Lifetime/ctPrompt_rapidity%d_pt%d.png",yBin+1,ptBin+1);
  char FilenameBK[200];
  sprintf(FilenameBK,"Plots/Lifetime/data_ctBackground_rapidity%d_pt%d.png",yBin+1,ptBin+1);

  if(!noNP){
  TCanvas* NonPromptCanvas = new TCanvas("NonPromptCanvas","NonPromptCanvas",1400, 600);
  NonPromptCanvas->Divide(2);  NonPromptCanvas->SetFillColor(kWhite);
  NonPromptCanvas->cd(1)->SetLogy(1);
  NonPromptCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); NonPromptFrame->GetYaxis()->SetTitleOffset(2) ; NonPromptFrame->Draw();
  NonPromptCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); NonPromptFrame->GetYaxis()->SetTitleOffset(2) ; NonPromptFrame->Draw();
  NonPromptCanvas->SaveAs(FilenameNP);
  NonPromptCanvas->Close();
  }
  if(!noPR){
  TCanvas* PromptCanvas = new TCanvas("PromptCanvas","PromptCanvas",1400, 600);
  PromptCanvas->Divide(2);  PromptCanvas->SetFillColor(kWhite);
  PromptCanvas->cd(1)->SetLogy(1);
  PromptCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); PromptFrame->GetYaxis()->SetTitleOffset(2) ; PromptFrame->Draw();
  PromptCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); PromptFrame->GetYaxis()->SetTitleOffset(2) ; PromptFrame->Draw(); if(!pereverr) PromptLegend->Draw();
  PromptCanvas->SaveAs(FilenamePR);
  PromptCanvas->Close();
  }
  if(!noBK){
  TCanvas* BackgroundCanvas = new TCanvas("BackgroundCanvas","BackgroundCanvas",1400, 600);
  BackgroundCanvas->Divide(2);  BackgroundCanvas->SetFillColor(kWhite);
  BackgroundCanvas->cd(1)->SetLogy(1);
  BackgroundCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); BackgroundFrame->GetYaxis()->SetTitleOffset(2) ; BackgroundFrame->Draw();
  BackgroundCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); BackgroundFrame->GetYaxis()->SetTitleOffset(2) ; BackgroundFrame->Draw(); BackgroundLegend->Draw();
  BackgroundCanvas->SaveAs(FilenameBK);
  BackgroundCanvas->Close();
  }


  cout<<"hallo"<<endl;




////////////////////////////////// 	"Pull" ///////////////////////////////////////////////////////////////////////////////////
if(!noPull){

	if(!noPR){


		  TH1 * PromptModelHist =  gaussSum.createHistogram("Jpsict",ctauBins) ;
		  TH1 * PromptHist = dataPRbin->createHistogram("Jpsict",ctauBins);

		  double deltaPrompt[ctauBins];
		  double errdeltaPrompt[ctauBins];
		  double FitValue;
		  double deltaPromptMax=0;
		  double promptMean[ctauBins];
		  double errpromptMean[ctauBins];

		  for(int BinRun = 1; BinRun < ctauBins+1; ++BinRun) {

			  FitValue=PromptModelHist->GetBinContent(BinRun)*dataPRbin->sumEntries()/PromptModelHist->GetSumOfWeights();

			  if(PromptHist->GetBinError(BinRun)==0) deltaPrompt[BinRun-1]=100000;
			  else deltaPrompt[BinRun-1]=((PromptHist->GetBinContent(BinRun)-FitValue)/PromptHist->GetBinError(BinRun));
			  cout<<deltaPrompt[BinRun-1]<<endl;

			  errdeltaPrompt[BinRun-1]=0;//PromptHist->GetBinError(BinRun);

			  if (fabs(deltaPrompt[BinRun-1])+fabs(errdeltaPrompt[BinRun-1]) > deltaPromptMax && fabs(deltaPrompt[BinRun-1])+fabs(errdeltaPrompt[BinRun-1]) < 99999)
			   	    		  deltaPromptMax=fabs(deltaPrompt[BinRun-1])+fabs(errdeltaPrompt[BinRun-1]);

			  promptMean[BinRun-1]=-1+3.5/ctauBins*BinRun-1.75/ctauBins;
			   	    errpromptMean[BinRun-1]=0;
		  }

		  	    deltaPromptMax*=1.05;

			 TCanvas *PromptPullCanvas = new TCanvas("PromptCanvas","",1000,1000);

			  TPad* PromptPad = new TPad("PromptPad","",0,0,1,1);
			  PromptPad->SetFillColor(0);

			  TPad* PromptPullPad = new TPad("PromptPullPad","",0,0,1,0.3);
			  PromptPullPad->SetFillColor(0);

			  PromptPad->Draw();
			  PromptPad->cd();
			  PromptPad->SetBottomMargin(0.4);
			  gPad->SetLeftMargin(0.15) ; gPad->SetFillColor(kWhite); PromptFrame->GetYaxis()->SetTitleOffset(1.5) ;
			  PromptFrame->Draw();

			  PromptPullPad->Draw();
			  PromptPullPad->cd();
			  PromptPullPad->SetGrid();
			  gPad->SetLeftMargin(0.15) ; gPad->SetFillColor(kWhite); PromptFrame->GetYaxis()->SetTitleOffset(1.5);

			  TH1F *PromptPullhistoHisto = PromptPullPad->DrawFrame(-1,-deltaPromptMax,2.5,deltaPromptMax);
			  PromptPullhistoHisto->SetYTitle("(Data - Fit)/Sigma");
			  PromptPullhistoHisto->GetYaxis()->SetTitleOffset(1.5);

			  TGraphErrors *PromptPullGraph = new TGraphErrors(ctauBins,promptMean,deltaPrompt,errpromptMean,errdeltaPrompt);
			  PromptPullGraph->SetMarkerColor(kBlack);
			  PromptPullGraph->SetMarkerStyle(10);
			  PromptPullGraph->Draw("P");

			  char Filename[50];
			  sprintf(Filename,"Plots/Lifetime/PromptPull_rapidity%d_pt%d.png",yBin+1,ptBin+1);

		   PromptPullCanvas->SaveAs(Filename);
		   PromptPullCanvas->Close();

/*

     TCanvas* CSPol2DCanvas = new TCanvas("CSPol2DCanvas","CSPol2DCanvas",1400, 700);
     CSPol2DCanvas->Divide(2);  CSPol2DCanvas->SetFillColor(kWhite);
     CSPol2DCanvas->cd(1) ; PromptModelHist->Draw();
     CSPol2DCanvas->cd(2) ; PromptHist->Draw();

     char FilenameCScos[200];
     sprintf(FilenameCScos,"Plots/Lifetime/PromptHists_rapidity%d_pt%d.png",yBin+1,ptBin+1);


     CSPol2DCanvas->SaveAs(FilenameCScos);
*/

}





}

if(!noBK) delete BackgroundFrame;
if(!noNP) delete NonPromptFrame;
if(!noPR) delete PromptFrame;

  delete lifetinp_;
  delete PromptLegend;
	  }
	}



  delete dataTreesNP;
  delete dataNP;
  delete dataTreesPR;
  delete dataPR;
//  delete dataTreesBK;
//  delete dataBK;
  delete fInput;
  delete dataNPbinHist;



  return 0;
}
