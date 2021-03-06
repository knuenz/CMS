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



//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"



int main(int argc, char**argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;

  gSystem->mkdir("Plots/Polarization");

  gStyle->SetPalette(1);

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
  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCType_idx);
  varlist.add(JpsiType_idx);
  varlist.add(MCweight);

  RooRealVar costh_CS_Gen("costh_CS_Gen","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS_Gen("phi_CS_Gen","#phi_{CS} [deg]",0,360);
  RooRealVar costh_HX_Gen("costh_HX_Gen","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX_Gen("phi_HX_Gen","#phi_{HX} [deg]",0,360);

  varlist.add(costh_CS_Gen);
  varlist.add(phi_CS_Gen);
  varlist.add(costh_HX_Gen);
  varlist.add(phi_HX_Gen);


  TChain *dataTreesPR = new TChain("data");
//  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root");
//  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo_CSth_minus1.root");
//  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo_CSth_plus1.root");
//  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo_HXth_minus1.root");
//  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo_HXth_plus1.root");
//  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_pol_Mu0Track0Jpsi_MCprompt.root");
  dataTreesPR->Add("/scratch/knuenz/Polarization/RootInput/TTree_pol_Mu0Track0Jpsi_MCprompt.root");


   Char_t *fileNameInputCS = "jPsiFitFinal_cs_Gen.root";
   TFile* fInputCS = new TFile(fileNameInputCS);

   Char_t *fileNameInputHX = "jPsiFitFinal_hx_Gen.root";
   TFile* fInputHX = new TFile(fileNameInputHX);


  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",dataTreesPR,varlist);//(,0,"MCweight");
  dataPR->Print();

  RooDataSet *dataNPbin, *dataPRbin, *dataBKbin, *dataREALbin;
//  RooDataHist *dataPRbinhist;

  CompositeModelBuilder* cs;
  CompositeModelBuilder* hx;

  RooPlot* CS_cos_Frame;
  RooPlot* CS_phi_Frame;
  RooPlot* HX_cos_Frame;
  RooPlot* HX_phi_Frame;

  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins; ++yBin) {
	  for(int ptBin = 0; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {





      cout<<"pT"<<ptBin+1<<"rapidity"<<yBin+1<<endl;

	  char reduceNP[200];
	  sprintf(reduceNP,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
	  char reduceBK[200];
	  sprintf(reduceBK,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiMass < %f | JpsiMass > %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],jpsi::polMassJpsi[yBin+1]-jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgLow,jpsi::polMassJpsi[yBin+1]+jpsi::sigmaMassJpsi[yBin+1]*jpsi::nSigBkgHigh);

	  dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);
	  RooDataHist dataPRbinhistCS("dataPRbinhistCS","dataPRbinhistCS", RooArgSet(costh_CS_Gen, phi_CS_Gen));
	  RooDataHist dataPRbinhistHX("dataPRbinhistHX","dataPRbinhistHX", RooArgSet(costh_HX_Gen, phi_HX_Gen));
	  dataPRbin->Print();

	  dataPRbinhistCS.add(*dataPRbin);
	  dataPRbinhistHX.add(*dataPRbin);

//	  TTree *TreedataPRbin = (TTree*)dataPRbin->tree();

	  costh_CS_Gen.setBins(jpsi::kNbBinsCosT);
	  phi_CS_Gen.setBins(jpsi::kNbBinsPhiPol);

	  costh_HX_Gen.setBins(jpsi::kNbBinsCosT);
	  phi_HX_Gen.setBins(jpsi::kNbBinsPhiPol);


	//  TH2F * hist_CS_PR = (TH2F*)dataPRbinhistCS.createHistogram("hist_CS_PR",costh_CS,Binning(jpsi::nBinsCosT),YVar(phi_CS,Binning(jpsi::nBinsPhiPol)));
	//  TH2F * hist_HX_PR = (TH2F*)dataPRbinhistHX.createHistogram("hist_HX_PR",costh_HX,Binning(jpsi::nBinsCosT),YVar(phi_HX,Binning(jpsi::nBinsPhiPol)));

	  TH2F * hist_CS_PR = dataPRbin->createHistogram(costh_CS_Gen,phi_CS_Gen);//"hist_CS_PR",costh_CS,Binning(jpsi::nBinsCosT),YVar(phi_CS,Binning(jpsi::nBinsPhiPol))) ;
	  TH2F * hist_HX_PR = dataPRbin->createHistogram(costh_HX_Gen,phi_HX_Gen);//"hist_CS_PR",costh_CS,Binning(jpsi::nBinsCosT),YVar(phi_CS,Binning(jpsi::nBinsPhiPol))) ;


	  char DirectoryPath[200];
	  sprintf(DirectoryPath,"/pt%d_rapidity%d/",ptBin+1,yBin+1);


	  TDirectory *InputDirectoryCS = (TDirectory*)fInputCS->GetDirectory(DirectoryPath);
	  TDirectory *InputDirectoryHX = (TDirectory*)fInputHX->GetDirectory(DirectoryPath);


	  cs = new CompositeModelBuilder("CS_Gen");
	  cs->setUsePrompt(true);
	  cs->setUseNonPrompt(false);
	  cs->setUseBkg(false);
	  cs->setUseMass(false);
	  cs->setUseLifetime(false);
	  cs->setUsePol(true);
	  cs->setUseAcceptanceMaps(false);

	  if(fInputCS->GetDirectory(DirectoryPath)!=NULL){
	  cs->loadParameters(*InputDirectoryCS);
	  }
	  else cout<<"TDirectory "<<DirectoryPath<<" does not exist in File "<<fileNameInputCS<<endl;
	  cs->initModel(JpsiMass,Jpsict,costh_CS_Gen,phi_CS_Gen);

	  TH1 * hist_CS_PR_model = cs->model()->createHistogram("hist_CS_PR_model",costh_CS_Gen,Binning(jpsi::kNbBinsCosT),YVar(phi_CS_Gen,Binning(jpsi::kNbBinsPhiPol))) ;

      hx = new CompositeModelBuilder("HX_Gen");
 	  hx->setUsePrompt(true);
 	  hx->setUseNonPrompt(false);
	  hx->setUseBkg(false);
 	  hx->setUseMass(false);
 	  hx->setUseLifetime(false);
 	  hx->setUsePol(true);
	  hx->setUseAcceptanceMaps(false);

	  if(fInputHX->GetDirectory(DirectoryPath)!=NULL){
	  hx->loadParameters(*InputDirectoryHX);
	  }
	  else cout<<"TDirectory "<<DirectoryPath<<" does not exist in File "<<fileNameInputHX<<endl;
 	  hx->initModel(JpsiMass,Jpsict,costh_HX_Gen,phi_HX_Gen);

 	  TH1 * hist_HX_PR_model = hx->model()->createHistogram("hist_HX_PR_model",costh_HX_Gen,Binning(jpsi::kNbBinsCosT),YVar(phi_HX_Gen,Binning(jpsi::kNbBinsPhiPol))) ;


  CS_cos_Frame = new RooPlot;
  CS_cos_Frame = costh_CS_Gen.frame(Bins(50)) ;
  dataPRbin->plotOn(CS_cos_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
  cs->model()->plotOn(CS_cos_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
  double chi2_CS_cos_Frame=CS_cos_Frame->chiSquare();
  CS_cos_Frame->SetMinimum(1);
  char CS_cos_Title[200];
  sprintf(CS_cos_Title,"Cos#theta_{CS} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_CS_cos_Frame);
  CS_cos_Frame->SetTitle(CS_cos_Title);

  CS_phi_Frame = new RooPlot;
  CS_phi_Frame = phi_CS_Gen.frame(Bins(50)) ;
  dataPRbin->plotOn(CS_phi_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
  cs->model()->plotOn(CS_phi_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
  double chi2_CS_phi_Frame=CS_phi_Frame->chiSquare();
  CS_phi_Frame->SetMinimum(1);
  char CS_phi_Title[200];
  sprintf(CS_phi_Title,"#phi_{CS} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_CS_phi_Frame);
  CS_phi_Frame->SetTitle(CS_phi_Title);


  char CS2DTitle[200];
  sprintf(CS2DTitle,"MC %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
  char CS2DFitTitle[200];
  sprintf(CS2DFitTitle,"Fit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);

  cout<<"CanvasCS"<<endl;

  TCanvas* CS_cos_Canvas = new TCanvas("CS_cos_Canvas","CS_cos_Canvas",700, 600);
  CS_cos_Canvas->SetFillColor(kWhite);
  CS_cos_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); CS_cos_Frame->GetYaxis()->SetTitleOffset(2) ; CS_cos_Frame->Draw(); //PromptLegend->Draw();

  TCanvas* CS_phi_Canvas = new TCanvas("CS_phi_Canvas","CS_phi_Canvas",700, 600);
  CS_phi_Canvas->SetFillColor(kWhite);
  CS_phi_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); CS_phi_Frame->GetYaxis()->SetTitleOffset(2) ; CS_phi_Frame->Draw(); //PromptLegend->Draw();


  TCanvas* CSPol2DCanvas = new TCanvas("CSPol2DCanvas","CSPol2DCanvas",1400, 700);
  CSPol2DCanvas->Divide(2);  CSPol2DCanvas->SetFillColor(kWhite);
  CSPol2DCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_CS_PR->SetStats(0); hist_CS_PR->GetXaxis()->SetTitle("cos #theta_{CS}"); hist_CS_PR->GetYaxis()->SetTitle("#phi_{CS} [deg]"); hist_CS_PR->GetYaxis()->SetTitleOffset(2); hist_CS_PR->SetTitle(CS2DTitle); hist_CS_PR->Draw("colz");
  CSPol2DCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_CS_PR_model->SetStats(0); hist_CS_PR_model->SetTitle(CS2DFitTitle); hist_CS_PR_model->GetYaxis()->SetTitleOffset(2); hist_CS_PR_model->Draw("colz");

  char FilenameCScos[200];
  sprintf(FilenameCScos,"Plots/Polarization/CS_costh_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenameCSphi[200];
  sprintf(FilenameCSphi,"Plots/Polarization/CS_phi_pt%d_rapidity%d.png",ptBin+1,yBin+1);
  char FilenameCS2D[200];
  sprintf(FilenameCS2D,"Plots/Polarization/CS_2D_pt%d_rapidity%d.png",ptBin+1,yBin+1);

  CS_phi_Canvas->SaveAs(FilenameCSphi);
  CS_cos_Canvas->SaveAs(FilenameCScos);
  CSPol2DCanvas->SaveAs(FilenameCS2D);

  CS_phi_Canvas->Close();
  CS_cos_Canvas->Close();
  CSPol2DCanvas->Close();



  HX_cos_Frame = new RooPlot;
   HX_cos_Frame = costh_HX_Gen.frame(Bins(50)) ;
   dataPRbin->plotOn(HX_cos_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
   hx->model()->plotOn(HX_cos_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
   double chi2_HX_cos_Frame=HX_cos_Frame->chiSquare();
   HX_cos_Frame->SetMinimum(1);
   char HX_cos_Title[200];
   sprintf(HX_cos_Title,"Cos#theta_{HX} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_HX_cos_Frame);
   HX_cos_Frame->SetTitle(HX_cos_Title);

   HX_phi_Frame = new RooPlot;
   HX_phi_Frame = phi_HX_Gen.frame(Bins(50)) ;
   dataPRbin->plotOn(HX_phi_Frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
   hx->model()->plotOn(HX_phi_Frame, Normalization(1.0), LineWidth(3),LineColor(kBlue));
   double chi2_HX_phi_Frame=HX_phi_Frame->chiSquare();
   HX_phi_Frame->SetMinimum(1);
   char HX_phi_Title[200];
   sprintf(HX_phi_Title,"#phi_{HX} MC Prompt Component %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f, chi2 = %1.4f ",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1],chi2_HX_phi_Frame);
   HX_phi_Frame->SetTitle(HX_phi_Title);


   char HX2DTitle[200];
   sprintf(HX2DTitle,"MC %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
   char HX2DFitTitle[200];
   sprintf(HX2DFitTitle,"Fit %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);

   cout<<"CanvasHX"<<endl;


   TCanvas* HX_cos_Canvas = new TCanvas("HX_cos_Canvas","HX_cos_Canvas",700, 600);
   HX_cos_Canvas->SetFillColor(kWhite);
   HX_cos_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); HX_cos_Frame->GetYaxis()->SetTitleOffset(2) ; HX_cos_Frame->Draw(); //PromptLegend->Draw();

   TCanvas* HX_phi_Canvas = new TCanvas("HX_phi_Canvas","HX_phi_Canvas",700, 600);
   HX_phi_Canvas->SetFillColor(kWhite);
   HX_phi_Canvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); HX_phi_Frame->GetYaxis()->SetTitleOffset(2) ; HX_phi_Frame->Draw(); //PromptLegend->Draw();


   TCanvas* HXPol2DCanvas = new TCanvas("HXPol2DCanvas","HXPol2DCanvas",1400, 700);
   HXPol2DCanvas->Divide(2);  HXPol2DCanvas->SetFillColor(kWhite);
   HXPol2DCanvas->cd(1) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_HX_PR->SetStats(0); hist_HX_PR->GetXaxis()->SetTitle("cos #theta_{HX}"); hist_HX_PR->GetYaxis()->SetTitle("#phi_{HX} [deg]"); hist_HX_PR->GetYaxis()->SetTitleOffset(2); hist_HX_PR->SetTitle(HX2DTitle); hist_HX_PR->Draw("colz");
   HXPol2DCanvas->cd(2) ; gPad->SetLeftMargin(0.2) ; gPad->SetFillColor(kWhite); hist_HX_PR_model->SetStats(0); hist_HX_PR_model->SetTitle(HX2DFitTitle); hist_HX_PR_model->GetYaxis()->SetTitleOffset(2); hist_HX_PR_model->Draw("colz");

   char FilenameHXcos[200];
   sprintf(FilenameHXcos,"Plots/Polarization/HX_costh_pt%d_rapidity%d.png",ptBin+1,yBin+1);
   char FilenameHXphi[200];
   sprintf(FilenameHXphi,"Plots/Polarization/HX_phi_pt%d_rapidity%d.png",ptBin+1,yBin+1);
   char FilenameHX2D[200];
   sprintf(FilenameHX2D,"Plots/Polarization/HX_2D_pt%d_rapidity%d.png",ptBin+1,yBin+1);

   HX_phi_Canvas->SaveAs(FilenameHXphi);
   HX_cos_Canvas->SaveAs(FilenameHXcos);
   HXPol2DCanvas->SaveAs(FilenameHX2D);

   HX_phi_Canvas->Close();
   HX_cos_Canvas->Close();
   HXPol2DCanvas->Close();

	  }
	}



  delete dataTreesPR;
  delete dataPR;

  delete fInputCS;
  delete fInputHX;



  return 0;
}
