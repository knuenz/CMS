/*
 * minErrorToys.cc
 *
 *  Created on: Apr 19, 2012
 *      Author: valentinknuenz
 */

#include <iostream>
#include <sstream>
#include <iomanip>


// RooFit Includes
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooArgusBG.h"
#include "RooAbsReal.h"
#include "RooChebychev.h"
#include "RooBinning.h"
#include "RooMinuit.h"
#include "RooAddition.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFormula.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TLine.h"

int main(int argc, char** argv) {


	Char_t *FitID = "Default"; //Storage Directory
	int nToy;
	int niBkg;
	int niSig;
	int nSig0;
	int delta_nSig;
	int nBkg0;
	int delta_nBkg;

	  for( int i=0;i < argc; ++i ) {
		  if(std::string(argv[i]).find("FitID") != std::string::npos) {char* FitIDchar = argv[i]; char* FitIDchar2 = strtok (FitIDchar, "="); FitID = FitIDchar2; cout<<"FitID = "<<FitID<<endl;}
		  if(std::string(argv[i]).find("nToy") != std::string::npos) {char* nToychar = argv[i]; char* nToychar2 = strtok (nToychar, "p"); nToy = atof(nToychar2); cout<<"nToy = "<<nToy<<endl;}
		  if(std::string(argv[i]).find("niBkg") != std::string::npos) {char* niBkgchar = argv[i]; char* niBkgchar2 = strtok (niBkgchar, "p"); niBkg = atof(niBkgchar2); cout<<"niBkg = "<<niBkg<<endl;}
		  if(std::string(argv[i]).find("niSig") != std::string::npos) {char* niSigchar = argv[i]; char* niSigchar2 = strtok (niSigchar, "p"); niSig = atof(niSigchar2); cout<<"niSig = "<<niSig<<endl;}
		  if(std::string(argv[i]).find("nSig0") != std::string::npos) {char* nSig0char = argv[i]; char* nSig0char2 = strtok (nSig0char, "p"); nSig0 = atof(nSig0char2); cout<<"nSig0 = "<<nSig0<<endl;}
		  if(std::string(argv[i]).find("nBkg0") != std::string::npos) {char* nBkg0char = argv[i]; char* nBkg0char2 = strtok (nBkg0char, "p"); nBkg0 = atof(nBkg0char2); cout<<"nBkg0 = "<<nBkg0<<endl;}
		  if(std::string(argv[i]).find("delta_nSig") != std::string::npos) {char* delta_nSigchar = argv[i]; char* delta_nSigchar2 = strtok (delta_nSigchar, "p"); delta_nSig = atof(delta_nSigchar2); cout<<"delta_nSig = "<<delta_nSig<<endl;}
		  if(std::string(argv[i]).find("delta_nBkg") != std::string::npos) {char* delta_nBkgchar = argv[i]; char* delta_nBkgchar2 = strtok (delta_nBkgchar, "p"); delta_nBkg = atof(delta_nBkgchar2); cout<<"delta_nBkg = "<<delta_nBkg<<endl;}
	  }


		  char dirstruct[200];
	  sprintf(dirstruct,"ToyFigures/%s",FitID);
	  gSystem->mkdir(dirstruct);
	  gStyle->SetPalette(1);

	  char filename[200];
	  sprintf(filename,"%s/ToyResults.root",dirstruct);
	  TFile* resultsFile = new TFile(filename, "READ");


	  TTree* Results = (TTree*)resultsFile->Get("Results");

	  TH1D* h_S_scale   = new TH1D( "h_S_scale", "h_S_scale", 100, 0.,1.5);
	  Results->Draw("S/S_in>>h_S_scale");
	  TH1D* h_B_scale   = new TH1D( "h_B_scale", "h_B_scale", 100, 0.,1.5);
	  Results->Draw("B/B_in>>h_B_scale");

	  double S_scale=h_S_scale->GetMean();
	  double S_scaleBins=4;
	  cout<<"S_Scale "<<S_scale<<endl;
	  double B_scale=h_B_scale->GetMean();
	  double B_scaleBins=4;
	  cout<<"B_Scale "<<B_scale<<endl;


	  double BkgOverSigPlotAxis=(nBkg0+delta_nBkg*niBkg)*B_scale/((nSig0+delta_nSig*niSig)*S_scale);

	  char DrawTree[200];
	  sprintf(DrawTree,"S_in*%f:B_in*%f>>Sin_vs_Bin",S_scale,B_scale);
  TH2D* Sin_vs_Bin   = new TH2D( "Sin_vs_Bin", "Sin_vs_Bin", niBkg, nBkg0*B_scale, (nBkg0+delta_nBkg*niBkg)*B_scale, niSig, nSig0*S_scale, (nSig0+delta_nSig*niSig)*S_scale);
  Results->Draw(DrawTree,"Merr*(Merr<0.01)","colz");

//  TH2D* S_vs_B   = new TH2D( "S_vs_B", "S_vs_B", B_scaleBins*niBkg, nBkg0, (nBkg0+delta_nBkg*niBkg)*B_scale, S_scaleBins*niSig, nSig0, (nSig0+delta_nSig*niSig)*S_scale);
//  Results->Draw("S:B>>S_vs_B","Merr*(Merr<0.01)","colz");


  double SBPlotmax=(nSig0+delta_nSig*niSig)*S_scale/nBkg0/B_scale;
  cout<<"SBmax "<<SBPlotmax<<endl;

//  TH2D* S_vs_SB   = new TH2D( "S_vs_SB", "S_vs_SB", S_scaleBins*niSig, 0., SBPlotmax, S_scaleBins*niSig, nSig0, nSig0+delta_nSig*niSig*S_scale);
//  TH2D* S_vs_SB   = new TH2D( "S_vs_SB", "S_vs_SB", S_scaleBins*niSig, 0.,500,S_scaleBins*niSig, nSig0, nSig0+delta_nSig*niSig*S_scale);
//  Results->Draw("S:SB>>S_vs_SB","Merr*(Merr<0.01)","colz");

  TCanvas *c2 = new TCanvas("c2","c2",1200,800);

//  S_scale=1.;
//  B_scale=1.;

  for(int iSin_vs_BinPlots=1;iSin_vs_BinPlots<4;iSin_vs_BinPlots++){

  Sin_vs_Bin->SetYTitle("Number of signal events in 1.5#sigma region");
  Sin_vs_Bin->SetXTitle("Number of background events in 1.5#sigma region");
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.1);
  gPad->SetFillColor(kWhite);
  Sin_vs_Bin->GetXaxis()->SetTitleOffset(1.25);
  Sin_vs_Bin->GetYaxis()->SetTitleOffset(1.4);
  Sin_vs_Bin->SetTitle(0);
  Sin_vs_Bin->SetStats(0);
  Sin_vs_Bin->Draw("COLZ");

//  TLine* SBline;
  char iLineLegend[200];
  char signFormula[200];
  char SBFormula[200];
  const int nLine=10;
  double nSBplot[nLine]={.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1};
  double nSign[nLine]={4,5,6,7,8,9,10,12,14,16};
  TLatex *SBtext[nLine];
  TLatex *Sigtext[nLine];
  TF1 *SigLine[nLine];
  TF1 *SBLine[nLine];

  for(int iLine=0;iLine<nLine;iLine++){

	  double SigMax=nSBplot[iLine]*(nBkg0+delta_nBkg*niBkg)*B_scale;
	  double BkgMax=(nBkg0+delta_nBkg*niBkg)*B_scale;

	  double whereTexteInPlotX=BkgMax*0.8/TMath::Sqrt(iLine+1);
	  double whereTexteInPlotY;
	  double whereTexteInPlotXScale=0.925;

	  if(iSin_vs_BinPlots==2){
	  double k=SigMax/BkgMax;
  sprintf(SBFormula,"%f*x",k);
  SBLine[iLine]=new TF1("SB_formula",SBFormula,0,BkgMax);
  SBLine[iLine]->SetLineColor(kWhite);
  SBLine[iLine]->SetLineWidth( 2 );
  SBLine[iLine]->Draw("same");


  whereTexteInPlotY=SBLine[iLine]->Eval(whereTexteInPlotX)*whereTexteInPlotXScale;
  sprintf (iLineLegend,"S/B = %1.1f",nSBplot[iLine]);
  SBtext[iLine] = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,iLineLegend);
  SBtext[iLine]->SetTextSize(0.015);
  SBtext[iLine]->SetTextColor(kWhite);
  SBtext[iLine]->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, TMath::ATan(SBLine[iLine]->Derivative(whereTexteInPlotX)*BkgOverSigPlotAxis)*180./TMath::Pi(), 0.015, iLineLegend);
  SBtext[iLine]->Draw( "same" );

	  }

	  if(iSin_vs_BinPlots==3){
  double p=nSign[iLine]*nSign[iLine]*(-1);
  sprintf(signFormula,"(-%f/2.+TMath::Sqrt((%f/2.)*(%f/2.)-%f*x))",p,p,p,p);
  SigLine[iLine]=new TF1("sig_formula",signFormula,0,BkgMax);
  SigLine[iLine]->SetLineColor(kWhite);
  SigLine[iLine]->SetLineWidth( 2 );
  SigLine[iLine]->Draw("same");

  whereTexteInPlotY=SigLine[iLine]->Eval(whereTexteInPlotX)*whereTexteInPlotXScale;
  sprintf (iLineLegend,"S/Sqrt(S+B) = %1.1f",nSign[iLine]);
  Sigtext[iLine] = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,iLineLegend);
  Sigtext[iLine]->SetTextSize(0.015);
  Sigtext[iLine]->SetTextColor(kWhite);
  Sigtext[iLine]->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, TMath::ATan(SigLine[iLine]->Derivative(whereTexteInPlotX)*BkgOverSigPlotAxis)*180./TMath::Pi(), 0.015, iLineLegend);
  Sigtext[iLine]->Draw( "same" );

//  cout<<"angle "<<TMath::ATan(SigLine[iLine]->Derivative(whereTexteInPlotX))*360./TMath::Pi()/BkgOverSigPlotAxis<<endl;
}
  }

  if(iSin_vs_BinPlots==1) sprintf(filename,"%s/h_Sin_VS_Bin.pdf",dirstruct);
  if(iSin_vs_BinPlots==2) sprintf(filename,"%s/h_Sin_VS_Bin_SBimposed.pdf",dirstruct);
  if(iSin_vs_BinPlots==3) sprintf(filename,"%s/h_Sin_VS_Bin_Sigimposed.pdf",dirstruct);
  c2->SaveAs(filename);

  }

/*  c2 = new TCanvas("c2","c2",1200,800);
  S_vs_B->SetYTitle("Number of signal events in 1.5#sigma region");
  S_vs_B->SetXTitle("Number of background events in 1.5#sigma region");
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.1);
  gPad->SetFillColor(kWhite);
  S_vs_B->GetXaxis()->SetTitleOffset(1.25);
  S_vs_B->GetYaxis()->SetTitleOffset(1.4);
  S_vs_B->SetTitle(0);
  S_vs_B->SetStats(0);
  S_vs_B->Draw("COLZ");
  sprintf(filename,"%s/h_S_VS_B.pdf",dirstruct);
  c2->SaveAs(filename);

  c2 = new TCanvas("c2","c2",1200,800);
  S_vs_SB->SetYTitle("Number of signal events in 1.5#sigma region");
  S_vs_SB->SetXTitle("S/B ratio in 1.5#sigma region");
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.1);
  gPad->SetFillColor(kWhite);
  S_vs_SB->GetXaxis()->SetTitleOffset(1.25);
  S_vs_SB->GetYaxis()->SetTitleOffset(1.4);
  S_vs_SB->SetTitle(0);
  S_vs_SB->SetStats(0);
  S_vs_SB->Draw("COLZ");
  sprintf(filename,"%s/h_S_VS_SB.pdf",dirstruct);
  c2->SaveAs(filename);
*/
	return 0;
}

