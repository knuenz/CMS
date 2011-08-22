#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TROOT.h"

#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"



// calculation of the frame-invariants lambdatheta^star and lambdaphi^star:

void calcLambdastar( double& lthstar, double& lphstar,
                     double& lth,     double& lph,     double& ltp ) {

  double LamPlus = 0.25 * ( lth - lph + TMath::Sqrt(TMath::Power( lth - lph, 2. ) + 4. * ltp*ltp ) );
  double LamMnus = 0.25 * ( lth - lph - TMath::Sqrt(TMath::Power( lth - lph, 2. ) + 4. * ltp*ltp ) );

  lthstar = ( lth - 3.*LamPlus ) / ( 1. + LamPlus );
  lphstar = ( lph + LamPlus )    / ( 1. + LamPlus );

  double lphstarMnus = ( lph + LamMnus ) / ( 1. + LamMnus );

  if ( TMath::Abs(lphstarMnus) < TMath::Abs(lphstar) ) {
     lthstar = ( lth - 3.*LamMnus ) / ( 1. + LamMnus );
     lphstar = lphstarMnus;
  }
}


	  const int NrapBins_=2;
	  const int NptBins_=2;

double Gauss(double* x, double* par)
{
	double a=2*TMath::Pi();
	double gaussval=par[2]/(par[1]*TMath::Sqrt(a))*TMath::Exp(-pow(x[0]-par[0],2)/(2*pow(par[1],2)));
  return gaussval;
}

void PlotObject(TH1D* object_histo, char Xaxistitle[200], double l_min_object, double l_max_object, double l_step_1D_object, double object_histo_mean, double object_histo_meanerr, double object_histo_sigma, double object_histo_sigmaerr, char filename [200], int nGenerations, bool pull, double injvalue){

	  TCanvas* c1 = new TCanvas("c1", "c1", 10, 28, 588,563);
	  c1->Range(-1.370833,-114.012,1.029167,745.8285);
	  c1->SetFillColor(0);
	  c1->SetBorderMode(0);
	  c1->SetBorderSize(0);
	  c1->SetLeftMargin(0.10215278);
	  c1->SetRightMargin(0.10215278);
	  c1->SetTopMargin(0.01841621);
	  c1->SetBottomMargin(0.1325967);
	  c1->SetFrameBorderMode(0);

	 double plotborder = 1.1;
	 double plotMax = plotborder * object_histo->GetMaximum();
	 object_histo->SetMinimum(0.);
	 object_histo->SetMaximum(plotMax);

	 object_histo->GetXaxis()->SetTitle(Xaxistitle);
	 object_histo->GetXaxis()->SetLabelOffset(0.028);
	 object_histo->GetXaxis()->SetTitleSize(0.05);
	 object_histo->GetXaxis()->SetTickLength(-0.03);
	 object_histo->GetXaxis()->SetTitleOffset(1.20);
//	 object_histo->GetYaxis()->SetTitle("Pull Distribution");
	 object_histo->GetYaxis()->SetLabelOffset(0.032);
	 object_histo->GetYaxis()->SetTitleSize(0.05);
	 object_histo->GetYaxis()->SetTickLength(-0.03);
	 object_histo->GetYaxis()->SetTitleOffset(1.55);
	 object_histo->Draw("");
	 object_histo->Draw( "same" );

	 TF1* gauss = new TF1("gauss",Gauss,l_min_object,l_max_object,3);
	 gauss->SetParameter(0,object_histo_mean);
	 gauss->SetParameter(1,object_histo_sigma);
	 gauss->SetParameter(2,nGenerations*l_step_1D_object);
	 gauss->SetLineColor( kGreen+2 );
	 gauss->SetLineWidth( 2 );
	 gauss->SetLineStyle( 1 );
	 gauss->Draw( "same" );

	 if (pull){
		 TF1* gaussRef = new TF1("gaussRef",Gauss,l_min_object,l_max_object,3);
		 gaussRef->SetParameter(0,0);
		 gaussRef->SetParameter(1,1);
		 gaussRef->SetParameter(2,nGenerations*l_step_1D_object);
		 gaussRef->SetLineColor( kRed );
		 gaussRef->SetLineWidth( 1 );
		 gaussRef->SetLineStyle( 3 );
		 gaussRef->Draw( "same" );
	 }

	 if (!pull){
			TLine* injLine = new TLine( injvalue, 0, injvalue, plotMax );
			injLine->SetLineWidth( 1 );
			injLine->SetLineStyle( 2 );
			injLine->SetLineColor( kRed );
			injLine->Draw( "same" );
	 }
	 char meantext[200], sigmatext[200];
	 if (!pull) sprintf(meantext, "#splitline{#mu_{#lambda} = %1.3f #pm %1.3f}{#sigma_{#lambda} = %1.3f #pm %1.3f}",object_histo_mean,object_histo_meanerr,object_histo_sigma,object_histo_sigmaerr);
	 if (pull) sprintf(meantext, "#splitline{#mu_{z} = %1.3f #pm %1.3f}{#sigma_{z} = %1.3f #pm %1.3f}",object_histo_mean,object_histo_meanerr,object_histo_sigma,object_histo_sigmaerr);

	 TLatex *text = new TLatex(l_min_object+0.05*(l_max_object-l_min_object),plotMax*0.9,meantext);
	 text->SetTextSize(0.035);
	 text->Draw( "same" );

	 if (!pull){
		 sprintf(meantext, "#splitline{#lambda^{Inj.} = %1.3f}{#delta_{#lambda} = %1.1f #sigma_{#lambda}}",injvalue,TMath::Abs(object_histo_mean-injvalue)/object_histo_sigma);
		 text = new TLatex(l_min_object+0.05*(l_max_object-l_min_object),plotMax*0.75,meantext);
		 text->SetTextSize(0.035);
		 text->Draw( "same" );
	 }
	 if (pull){
		 sprintf(meantext, "#splitline{#delta_{#mu_{z}} = %1.1f #sigma_{#mu_{z}}}{#delta_{#sigma_{z}} = %1.1f #sigma_{#sigma_{z}}}",object_histo_mean/object_histo_meanerr,TMath::Abs(object_histo_sigma-1)/object_histo_sigmaerr);
		 text = new TLatex(l_min_object+0.05*(l_max_object-l_min_object),plotMax*0.75,meantext);
		 text->SetTextSize(0.035);
		 text->Draw( "same" );
	 }
	 c1->Print( filename );

	 delete gauss ;
	 delete c1 ;


}


void PlotRapPt(int ptBinMin, int ptBinMax, int rapBin, double yMin, double yMax, char yAxisTitle[200], char filename[200], double lmean[ToyMC::nPtBins], double lmeanerr[ToyMC::nPtBins], bool DrawLine, double lamline){

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

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	TH1F *plotHisto = new TH1F;
	plotHisto = plotCanvas->DrawFrame(ToyMC::ptRange[ptBinMin-1],yMin,ToyMC::ptRange[ptBinMax],yMax);
	plotHisto->SetXTitle("p_{T} [GeV/c]");
	plotHisto->SetYTitle(yAxisTitle);
	plotHisto->GetYaxis()->SetTitleOffset(1.5);


	double yLine=lamline;

	TLine* RefLine = new TLine( ToyMC::ptRange[ptBinMin-1], yLine, ToyMC::ptRange[ptBinMax], yLine );
	RefLine->SetLineWidth( 2 );
	RefLine->SetLineStyle( 2 );
	RefLine->SetLineColor( kGreen+2 );
	if(DrawLine) RefLine->Draw( "same" );

	int nBinspT=ptBinMax-ptBinMin+1;
	double ptCentre_[nBinspT];
	double ptCentreErr_[nBinspT];

	int pt=0;
	for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
		ptCentre_[pt]=ToyMC::ptCentre[ptBin-1];
		ptCentreErr_[pt]=ToyMC::ptCentreErr[ptBin-1];
	pt++;
	}

 	TGraphErrors *plotGraph = new TGraphErrors(nBinspT,ptCentre_,lmean,ptCentreErr_,lmeanerr);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor[rapBin]);
	plotGraph->SetLineColor(ToyMC::MarkerColor[rapBin]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle[rapBin]);
	plotGraph->SetMarkerSize(2);
	plotGraph->Draw("P");


	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

}






void PlotComparisonRapPt(int ptBinMin, int ptBinMax, int rapBin, double yMin, double yMax, char yAxisTitle[200], char filename[200], double lmean1[ToyMC::nPtBins], double lmean1err[ToyMC::nPtBins], double lmean2[ToyMC::nPtBins], double lmean2err[ToyMC::nPtBins], double lmean3[ToyMC::nPtBins], double lmean3err[ToyMC::nPtBins], double lamline){

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

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	TH1F *plotHisto = new TH1F;
	plotHisto = plotCanvas->DrawFrame(ToyMC::ptRange[ptBinMin-1],yMin,ToyMC::ptRange[ptBinMax],yMax);
	plotHisto->SetXTitle("p_{T} [GeV/c]");
	plotHisto->SetYTitle(yAxisTitle);
	plotHisto->GetYaxis()->SetTitleOffset(1.5);

	TLegend* plotLegend=new TLegend(0.825,0.75,0.95,0.9);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.035);
	plotLegend->SetBorderSize(1);

	char legendentry[200];

	int nBinspT=ptBinMax-ptBinMin+1;
	double ptCentre_[nBinspT];
	double ptCentreErr_[nBinspT];

	double yLine=lamline;

	TLine* RefLine = new TLine( ToyMC::ptRange[ptBinMin-1], yLine, ToyMC::ptRange[ptBinMax], yLine );
	RefLine->SetLineWidth( 2 );
	RefLine->SetLineStyle( 2 );
	RefLine->SetLineColor( kGreen+2 );
	RefLine->Draw( "same" );

	int pt=0;
	for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
		ptCentre_[pt]=ToyMC::ptCentre[ptBin-1]+1;
		ptCentreErr_[pt]=ToyMC::ptCentreErr[ptBin-1];
	pt++;
	}

 	TGraphErrors *plotGraph = new TGraphErrors(nBinspT,ptCentre_,lmean1,ptCentreErr_,lmean1err);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor[1]);
	plotGraph->SetLineColor(ToyMC::MarkerColor[1]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle[rapBin]);
	plotGraph->SetMarkerSize(2);
	plotGraph->Draw("P");
	sprintf(legendentry,"CS");
	plotLegend->AddEntry(plotGraph,legendentry,"ple");

	pt=0;
	for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
		ptCentre_[pt]=ToyMC::ptCentre[ptBin-1];
		ptCentreErr_[pt]=ToyMC::ptCentreErr[ptBin-1];
	pt++;
	}

 	plotGraph = new TGraphErrors(nBinspT,ptCentre_,lmean2,ptCentreErr_,lmean2err);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor2[1]);
	plotGraph->SetLineColor(ToyMC::MarkerColor2[1]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle2[rapBin]);
	plotGraph->SetMarkerSize(2);
	plotGraph->Draw("P");
	sprintf(legendentry,"HX");
	plotLegend->AddEntry(plotGraph,legendentry,"ple");

	pt=0;
	for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
		ptCentre_[pt]=ToyMC::ptCentre[ptBin-1]-1;
		ptCentreErr_[pt]=ToyMC::ptCentreErr[ptBin-1];
	pt++;
	}

 	plotGraph = new TGraphErrors(nBinspT,ptCentre_,lmean3,ptCentreErr_,lmean3err);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor3[1]);
	plotGraph->SetLineColor(ToyMC::MarkerColor3[1]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle3[rapBin]);
	plotGraph->SetMarkerSize(2);
	plotGraph->Draw("P");
	sprintf(legendentry,"PX");
	plotLegend->AddEntry(plotGraph,legendentry,"ple");


	plotLegend->Draw();

	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

}











int main(int argc, char** argv) {


	int nGenerations=1;
	int polScen=1;
	int frame=1;
	Char_t *dirstruct = "ToyDirectory_Default";
	int rapBinMin=1;
	int rapBinMax=1;
	int ptBinMin=1;
	int ptBinMax=1;



	  for( int i=0;i < argc; ++i ) {
	    if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
	    if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
	    if(std::string(argv[i]).find("rapBinMin") != std::string::npos) {char* rapBinMinchar = argv[i]; char* rapBinMinchar2 = strtok (rapBinMinchar, "p"); rapBinMin = atof(rapBinMinchar2); cout<<"rapBinMin = "<<rapBinMin<<endl;}
	    if(std::string(argv[i]).find("rapBinMax") != std::string::npos) {char* rapBinMaxchar = argv[i]; char* rapBinMaxchar2 = strtok (rapBinMaxchar, "p"); rapBinMax = atof(rapBinMaxchar2); cout<<"rapBinMax = "<<rapBinMax<<endl;}
	    if(std::string(argv[i]).find("Sigframe") != std::string::npos) {char* framechar = argv[i]; char* framechar2 = strtok (framechar, "p"); frame = atof(framechar2); cout<<"frame = "<<frame<<endl;}
	    if(std::string(argv[i]).find("nGenerations") != std::string::npos) {char* nGenerationschar = argv[i]; char* nGenerationschar2 = strtok (nGenerationschar, "p"); nGenerations = atof(nGenerationschar2); cout<<"nGenerations = "<<nGenerations<<endl;}
	    if(std::string(argv[i]).find("polScen") != std::string::npos) {char* polScenchar = argv[i]; char* polScenchar2 = strtok (polScenchar, "p"); polScen = atof(polScenchar2); cout<<"polScen = "<<polScen<<endl;}
	    if(std::string(argv[i]).find("Sig_") != std::string::npos) {dirstruct = argv[i];cout<<"dirstruct = "<<dirstruct<<endl;}
	  }










	cout<<"";

	  gROOT->Reset();
	  gROOT->SetStyle("Plain");
	  gStyle->SetOptFit(1111);
	  gStyle->SetOptStat(0);

	  char NatFrameName[200];
	  sprintf(NatFrameName,"%s",onia::frameLabel[frame-1]);





	  bool HX_is_natural;
	  bool PX_is_natural;
	  if(frame==2) HX_is_natural=true; if(frame==3) PX_is_natural=true; //else CS is the natural frame by default

	  char filename [500];

	  int NrapBins = rapBinMax-rapBinMin+1;
	  int NptBins = ptBinMax-ptBinMin+1;

	  bool emptyBin[NrapBins][NptBins];

	  double pull_lth_CS_mean[NrapBins][NptBins],	pull_lth_CS_meanerr[NrapBins][NptBins],	pull_lth_CS_sigma[NrapBins][NptBins],	pull_lth_CS_sigmaerr[NrapBins][NptBins];
	  double pull_lph_CS_mean[NrapBins][NptBins],	pull_lph_CS_meanerr[NrapBins][NptBins],	pull_lph_CS_sigma[NrapBins][NptBins],	pull_lph_CS_sigmaerr[NrapBins][NptBins];
	  double pull_ltp_CS_mean[NrapBins][NptBins],	pull_ltp_CS_meanerr[NrapBins][NptBins],	pull_ltp_CS_sigma[NrapBins][NptBins],	pull_ltp_CS_sigmaerr[NrapBins][NptBins];
	  double pull_lthstar_CS_mean[NrapBins][NptBins],	pull_lthstar_CS_meanerr[NrapBins][NptBins],	pull_lthstar_CS_sigma[NrapBins][NptBins],	pull_lthstar_CS_sigmaerr[NrapBins][NptBins];
	  double pull_lphstar_CS_mean[NrapBins][NptBins],	pull_lphstar_CS_meanerr[NrapBins][NptBins],	pull_lphstar_CS_sigma[NrapBins][NptBins],	pull_lphstar_CS_sigmaerr[NrapBins][NptBins];
	  double pull_ltilde_CS_mean[NrapBins][NptBins],	pull_ltilde_CS_meanerr[NrapBins][NptBins],	pull_ltilde_CS_sigma[NrapBins][NptBins],	pull_ltilde_CS_sigmaerr[NrapBins][NptBins];

	  double pull_lth_HX_mean[NrapBins][NptBins],	pull_lth_HX_meanerr[NrapBins][NptBins],	pull_lth_HX_sigma[NrapBins][NptBins],	pull_lth_HX_sigmaerr[NrapBins][NptBins];
	  double pull_lph_HX_mean[NrapBins][NptBins],	pull_lph_HX_meanerr[NrapBins][NptBins],	pull_lph_HX_sigma[NrapBins][NptBins],	pull_lph_HX_sigmaerr[NrapBins][NptBins];
	  double pull_ltp_HX_mean[NrapBins][NptBins],	pull_ltp_HX_meanerr[NrapBins][NptBins],	pull_ltp_HX_sigma[NrapBins][NptBins],	pull_ltp_HX_sigmaerr[NrapBins][NptBins];
	  double pull_lthstar_HX_mean[NrapBins][NptBins],	pull_lthstar_HX_meanerr[NrapBins][NptBins],	pull_lthstar_HX_sigma[NrapBins][NptBins],	pull_lthstar_HX_sigmaerr[NrapBins][NptBins];
	  double pull_lphstar_HX_mean[NrapBins][NptBins],	pull_lphstar_HX_meanerr[NrapBins][NptBins],	pull_lphstar_HX_sigma[NrapBins][NptBins],	pull_lphstar_HX_sigmaerr[NrapBins][NptBins];
	  double pull_ltilde_HX_mean[NrapBins][NptBins],	pull_ltilde_HX_meanerr[NrapBins][NptBins],	pull_ltilde_HX_sigma[NrapBins][NptBins],	pull_ltilde_HX_sigmaerr[NrapBins][NptBins];

	  double pull_lth_PX_mean[NrapBins][NptBins],	pull_lth_PX_meanerr[NrapBins][NptBins],	pull_lth_PX_sigma[NrapBins][NptBins],	pull_lth_PX_sigmaerr[NrapBins][NptBins];
	  double pull_lph_PX_mean[NrapBins][NptBins],	pull_lph_PX_meanerr[NrapBins][NptBins],	pull_lph_PX_sigma[NrapBins][NptBins],	pull_lph_PX_sigmaerr[NrapBins][NptBins];
	  double pull_ltp_PX_mean[NrapBins][NptBins],	pull_ltp_PX_meanerr[NrapBins][NptBins],	pull_ltp_PX_sigma[NrapBins][NptBins],	pull_ltp_PX_sigmaerr[NrapBins][NptBins];
	  double pull_lthstar_PX_mean[NrapBins][NptBins],	pull_lthstar_PX_meanerr[NrapBins][NptBins],	pull_lthstar_PX_sigma[NrapBins][NptBins],	pull_lthstar_PX_sigmaerr[NrapBins][NptBins];
	  double pull_lphstar_PX_mean[NrapBins][NptBins],	pull_lphstar_PX_meanerr[NrapBins][NptBins],	pull_lphstar_PX_sigma[NrapBins][NptBins],	pull_lphstar_PX_sigmaerr[NrapBins][NptBins];
	  double pull_ltilde_PX_mean[NrapBins][NptBins],	pull_ltilde_PX_meanerr[NrapBins][NptBins],	pull_ltilde_PX_sigma[NrapBins][NptBins],	pull_ltilde_PX_sigmaerr[NrapBins][NptBins];

	  double param_lth_CS_mean[NrapBins][NptBins],	param_lth_CS_meanerr[NrapBins][NptBins],	param_lth_CS_sigma[NrapBins][NptBins],	param_lth_CS_sigmaerr[NrapBins][NptBins];
	  double param_lph_CS_mean[NrapBins][NptBins],	param_lph_CS_meanerr[NrapBins][NptBins],	param_lph_CS_sigma[NrapBins][NptBins],	param_lph_CS_sigmaerr[NrapBins][NptBins];
	  double param_ltp_CS_mean[NrapBins][NptBins],	param_ltp_CS_meanerr[NrapBins][NptBins],	param_ltp_CS_sigma[NrapBins][NptBins],	param_ltp_CS_sigmaerr[NrapBins][NptBins];
	  double param_lthstar_CS_mean[NrapBins][NptBins],	param_lthstar_CS_meanerr[NrapBins][NptBins],	param_lthstar_CS_sigma[NrapBins][NptBins],	param_lthstar_CS_sigmaerr[NrapBins][NptBins];
	  double param_lphstar_CS_mean[NrapBins][NptBins],	param_lphstar_CS_meanerr[NrapBins][NptBins],	param_lphstar_CS_sigma[NrapBins][NptBins],	param_lphstar_CS_sigmaerr[NrapBins][NptBins];
	  double param_ltilde_CS_mean[NrapBins][NptBins],	param_ltilde_CS_meanerr[NrapBins][NptBins],	param_ltilde_CS_sigma[NrapBins][NptBins],	param_ltilde_CS_sigmaerr[NrapBins][NptBins];

	  double param_lth_HX_mean[NrapBins][NptBins],	param_lth_HX_meanerr[NrapBins][NptBins],	param_lth_HX_sigma[NrapBins][NptBins],	param_lth_HX_sigmaerr[NrapBins][NptBins];
	  double param_lph_HX_mean[NrapBins][NptBins],	param_lph_HX_meanerr[NrapBins][NptBins],	param_lph_HX_sigma[NrapBins][NptBins],	param_lph_HX_sigmaerr[NrapBins][NptBins];
	  double param_ltp_HX_mean[NrapBins][NptBins],	param_ltp_HX_meanerr[NrapBins][NptBins],	param_ltp_HX_sigma[NrapBins][NptBins],	param_ltp_HX_sigmaerr[NrapBins][NptBins];
	  double param_lthstar_HX_mean[NrapBins][NptBins],	param_lthstar_HX_meanerr[NrapBins][NptBins],	param_lthstar_HX_sigma[NrapBins][NptBins],	param_lthstar_HX_sigmaerr[NrapBins][NptBins];
	  double param_lphstar_HX_mean[NrapBins][NptBins],	param_lphstar_HX_meanerr[NrapBins][NptBins],	param_lphstar_HX_sigma[NrapBins][NptBins],	param_lphstar_HX_sigmaerr[NrapBins][NptBins];
	  double param_ltilde_HX_mean[NrapBins][NptBins],	param_ltilde_HX_meanerr[NrapBins][NptBins],	param_ltilde_HX_sigma[NrapBins][NptBins],	param_ltilde_HX_sigmaerr[NrapBins][NptBins];

	  double param_lth_PX_mean[NrapBins][NptBins],	param_lth_PX_meanerr[NrapBins][NptBins],	param_lth_PX_sigma[NrapBins][NptBins],	param_lth_PX_sigmaerr[NrapBins][NptBins];
	  double param_lph_PX_mean[NrapBins][NptBins],	param_lph_PX_meanerr[NrapBins][NptBins],	param_lph_PX_sigma[NrapBins][NptBins],	param_lph_PX_sigmaerr[NrapBins][NptBins];
	  double param_ltp_PX_mean[NrapBins][NptBins],	param_ltp_PX_meanerr[NrapBins][NptBins],	param_ltp_PX_sigma[NrapBins][NptBins],	param_ltp_PX_sigmaerr[NrapBins][NptBins];
	  double param_lthstar_PX_mean[NrapBins][NptBins],	param_lthstar_PX_meanerr[NrapBins][NptBins],	param_lthstar_PX_sigma[NrapBins][NptBins],	param_lthstar_PX_sigmaerr[NrapBins][NptBins];
	  double param_lphstar_PX_mean[NrapBins][NptBins],	param_lphstar_PX_meanerr[NrapBins][NptBins],	param_lphstar_PX_sigma[NrapBins][NptBins],	param_lphstar_PX_sigmaerr[NrapBins][NptBins];
	  double param_ltilde_PX_mean[NrapBins][NptBins],	param_ltilde_PX_meanerr[NrapBins][NptBins],	param_ltilde_PX_sigma[NrapBins][NptBins],	param_ltilde_PX_sigmaerr[NrapBins][NptBins];


	  double lambda_theta_injected_CS[NrapBins][NptBins];
	  double lambda_phi_injected_CS[NrapBins][NptBins];
	  double lambda_thetaphi_injected_CS[NrapBins][NptBins];
	  double lambda_thetastar_injected_CS[NrapBins][NptBins];
	  double lambda_phistar_injected_CS[NrapBins][NptBins];
	  double lambda_tilde_injected_CS[NrapBins][NptBins];

	  double lambda_theta_injected_HX[NrapBins][NptBins];
	  double lambda_phi_injected_HX[NrapBins][NptBins];
	  double lambda_thetaphi_injected_HX[NrapBins][NptBins];
	  double lambda_thetastar_injected_HX[NrapBins][NptBins];
	  double lambda_phistar_injected_HX[NrapBins][NptBins];
	  double lambda_tilde_injected_HX[NrapBins][NptBins];

	  double lambda_theta_injected_PX[NrapBins][NptBins];
	  double lambda_phi_injected_PX[NrapBins][NptBins];
	  double lambda_thetaphi_injected_PX[NrapBins][NptBins];
	  double lambda_thetastar_injected_PX[NrapBins][NptBins];
	  double lambda_phistar_injected_PX[NrapBins][NptBins];
	  double lambda_tilde_injected_PX[NrapBins][NptBins];


	  double l_min;
	  double l_max;
	  double l_step_1D;

	  double l_min_pull;
	  double l_max_pull;
	  double l_step_1D_pull;

	  double delta_l;

		int rap=0;
	    for(int iRap = rapBinMin; iRap < rapBinMax+1; iRap++){
			rap++;int pt=0;
	        for(int iPt = ptBinMin; iPt < ptBinMax+1; iPt++){
			pt++;



///////////////////// Extraction of mean injected parameters from GenResults.root //////////////////////////////////////////////////

			  l_min = -1.4;
			  l_max =  1.4;
			  l_step_1D = 0.02;

			  l_min_pull = -5;
			  l_max_pull =  5;
			  l_step_1D_pull = 0.5;

			  TH1D* inj_lth_CS = new TH1D( "inj_lth_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			  TH1D* inj_lph_CS = new TH1D( "inj_lph_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			  TH1D* inj_ltp_CS = new TH1D( "inj_ltp_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

			  TH1D* inj_lth_HX = new TH1D( "inj_lth_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			  TH1D* inj_lph_HX = new TH1D( "inj_lph_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			  TH1D* inj_ltp_HX = new TH1D( "inj_ltp_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

			  TH1D* inj_lth_PX = new TH1D( "inj_lth_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			  TH1D* inj_lph_PX = new TH1D( "inj_lph_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			  TH1D* inj_ltp_PX = new TH1D( "inj_ltp_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );


			  for(int iGen=1;iGen<nGenerations+1;iGen++){

			  sprintf(filename,"%s/rap%d_pT%d/Generation%d/GenResults.root",dirstruct,iRap,iPt,iGen);
			  TFile* GenResultFile = new TFile(filename);

			  if(GenResultFile->Get("GenResults") == NULL) {GenResultFile->Close(); continue;}

			  TTree* GenResults = (TTree*) GenResultFile->Get("GenResults");


			  double lamth_CS; double lamph_CS; double lamtp_CS;
			  double lamth_HX; double lamph_HX; double lamtp_HX;
			  double lamth_PX; double lamph_PX; double lamtp_PX;

			  TBranch        *b_inj_lth_CS;				GenResults->SetBranchAddress("lthCS", &lamth_CS, &b_inj_lth_CS);
			  TBranch        *b_inj_lph_CS;				GenResults->SetBranchAddress("lphCS", &lamph_CS, &b_inj_lph_CS);
			  TBranch        *b_inj_ltp_CS;				GenResults->SetBranchAddress("ltpCS", &lamtp_CS, &b_inj_ltp_CS);
			  TBranch        *b_inj_lth_HX;				GenResults->SetBranchAddress("lthHX", &lamth_HX, &b_inj_lth_HX);
			  TBranch        *b_inj_lph_HX;				GenResults->SetBranchAddress("lphHX", &lamph_HX, &b_inj_lph_HX);
			  TBranch        *b_inj_ltp_HX;				GenResults->SetBranchAddress("ltpHX", &lamtp_HX, &b_inj_ltp_HX);
			  TBranch        *b_inj_lth_PX;				GenResults->SetBranchAddress("lthPX", &lamth_PX, &b_inj_lth_PX);
			  TBranch        *b_inj_lph_PX;				GenResults->SetBranchAddress("lphPX", &lamph_PX, &b_inj_lph_PX);
			  TBranch        *b_inj_ltp_PX;				GenResults->SetBranchAddress("ltpPX", &lamtp_PX, &b_inj_ltp_PX);

			  GenResults->GetEvent( 0 );


			  inj_lth_CS->Fill( lamth_CS );
			  inj_lph_CS->Fill( lamph_CS );
			  inj_ltp_CS->Fill( lamtp_CS );

			  inj_lth_HX->Fill( lamth_HX );
			  inj_lph_HX->Fill( lamph_HX );
			  inj_ltp_HX->Fill( lamtp_HX );

			  inj_lth_PX->Fill( lamth_PX );
			  inj_lph_PX->Fill( lamph_PX );
			  inj_ltp_PX->Fill( lamtp_PX );


			  GenResultFile->Close();

			  }

			  lambda_theta_injected_CS[rap-1][pt-1]=inj_lth_CS->GetMean();
			  lambda_phi_injected_CS[rap-1][pt-1]=inj_lph_CS->GetMean();
			  lambda_thetaphi_injected_CS[rap-1][pt-1]=inj_ltp_CS->GetMean();
			  calcLambdastar(lambda_thetastar_injected_CS[rap-1][pt-1],lambda_phistar_injected_CS[rap-1][pt-1],lambda_theta_injected_CS[rap-1][pt-1],lambda_phi_injected_CS[rap-1][pt-1],lambda_thetaphi_injected_CS[rap-1][pt-1]);
			  lambda_tilde_injected_CS[rap-1][pt-1]=(lambda_theta_injected_CS[rap-1][pt-1]+3.*lambda_phi_injected_CS[rap-1][pt-1])/(1-lambda_phi_injected_CS[rap-1][pt-1]);

			  lambda_theta_injected_HX[rap-1][pt-1]=inj_lth_HX->GetMean();
			  lambda_phi_injected_HX[rap-1][pt-1]=inj_lph_HX->GetMean();
			  lambda_thetaphi_injected_HX[rap-1][pt-1]=inj_ltp_HX->GetMean();
			  calcLambdastar(lambda_thetastar_injected_HX[rap-1][pt-1],lambda_phistar_injected_HX[rap-1][pt-1],lambda_theta_injected_HX[rap-1][pt-1],lambda_phi_injected_HX[rap-1][pt-1],lambda_thetaphi_injected_HX[rap-1][pt-1]);
			  lambda_tilde_injected_HX[rap-1][pt-1]=(lambda_theta_injected_HX[rap-1][pt-1]+3.*lambda_phi_injected_HX[rap-1][pt-1])/(1-lambda_phi_injected_HX[rap-1][pt-1]);

			  lambda_theta_injected_PX[rap-1][pt-1]=inj_lth_PX->GetMean();
			  lambda_phi_injected_PX[rap-1][pt-1]=inj_lph_PX->GetMean();
			  lambda_thetaphi_injected_PX[rap-1][pt-1]=inj_ltp_PX->GetMean();
			  calcLambdastar(lambda_thetastar_injected_PX[rap-1][pt-1],lambda_phistar_injected_PX[rap-1][pt-1],lambda_theta_injected_PX[rap-1][pt-1],lambda_phi_injected_PX[rap-1][pt-1],lambda_thetaphi_injected_PX[rap-1][pt-1]);
			  lambda_tilde_injected_PX[rap-1][pt-1]=(lambda_theta_injected_PX[rap-1][pt-1]+3.*lambda_phi_injected_PX[rap-1][pt-1])/(1-lambda_phi_injected_PX[rap-1][pt-1]);

///////////////////// End of Extraction ////////////////////////////////////////////////////////////////////////////////////////////



	  l_min = -1.4;
	  l_max =  1.4;

	  l_min_pull = -10;
	  l_max_pull =  10;
	  l_step_1D_pull = 2*l_max_pull/50.;

	  delta_l = 1.2;
	  l_step_1D = 2*delta_l/50.;

	  double l_min_lth_CS = lambda_theta_injected_CS[rap-1][pt-1] - delta_l;
	  double l_max_lth_CS = lambda_theta_injected_CS[rap-1][pt-1] + delta_l;
	  double l_min_lph_CS = lambda_phi_injected_CS[rap-1][pt-1] - delta_l;
	  double l_max_lph_CS = lambda_phi_injected_CS[rap-1][pt-1] + delta_l;
	  double l_min_ltp_CS = lambda_thetaphi_injected_CS[rap-1][pt-1] - delta_l;
	  double l_max_ltp_CS = lambda_thetaphi_injected_CS[rap-1][pt-1] + delta_l;
	  double l_min_lthstar_CS = lambda_thetastar_injected_CS[rap-1][pt-1] - delta_l;
	  double l_max_lthstar_CS = lambda_thetastar_injected_CS[rap-1][pt-1] + delta_l;
	  double l_min_lphstar_CS = lambda_phistar_injected_CS[rap-1][pt-1] - delta_l;
	  double l_max_lphstar_CS = lambda_phistar_injected_CS[rap-1][pt-1] + delta_l;
	  double l_min_ltilde_CS = lambda_tilde_injected_CS[rap-1][pt-1] - delta_l;
	  double l_max_ltilde_CS = lambda_tilde_injected_CS[rap-1][pt-1] + delta_l;

	  double l_min_lth_HX = lambda_theta_injected_HX[rap-1][pt-1] - delta_l;
	  double l_max_lth_HX = lambda_theta_injected_HX[rap-1][pt-1] + delta_l;
	  double l_min_lph_HX = lambda_phi_injected_HX[rap-1][pt-1] - delta_l;
	  double l_max_lph_HX = lambda_phi_injected_HX[rap-1][pt-1] + delta_l;
	  double l_min_ltp_HX = lambda_thetaphi_injected_HX[rap-1][pt-1] - delta_l;
	  double l_max_ltp_HX = lambda_thetaphi_injected_HX[rap-1][pt-1] + delta_l;
	  double l_min_lthstar_HX = lambda_thetastar_injected_HX[rap-1][pt-1] - delta_l;
	  double l_max_lthstar_HX = lambda_thetastar_injected_HX[rap-1][pt-1] + delta_l;
	  double l_min_lphstar_HX = lambda_phistar_injected_HX[rap-1][pt-1] - delta_l;
	  double l_max_lphstar_HX = lambda_phistar_injected_HX[rap-1][pt-1] + delta_l;
	  double l_min_ltilde_HX = lambda_tilde_injected_HX[rap-1][pt-1] - delta_l;
	  double l_max_ltilde_HX = lambda_tilde_injected_HX[rap-1][pt-1] + delta_l;

	  double l_min_lth_PX = lambda_theta_injected_PX[rap-1][pt-1] - delta_l;
	  double l_max_lth_PX = lambda_theta_injected_PX[rap-1][pt-1] + delta_l;
	  double l_min_lph_PX = lambda_phi_injected_PX[rap-1][pt-1] - delta_l;
	  double l_max_lph_PX = lambda_phi_injected_PX[rap-1][pt-1] + delta_l;
	  double l_min_ltp_PX = lambda_thetaphi_injected_PX[rap-1][pt-1] - delta_l;
	  double l_max_ltp_PX = lambda_thetaphi_injected_PX[rap-1][pt-1] + delta_l;
	  double l_min_lthstar_PX = lambda_thetastar_injected_PX[rap-1][pt-1] - delta_l;
	  double l_max_lthstar_PX = lambda_thetastar_injected_PX[rap-1][pt-1] + delta_l;
	  double l_min_lphstar_PX = lambda_phistar_injected_PX[rap-1][pt-1] - delta_l;
	  double l_max_lphstar_PX = lambda_phistar_injected_PX[rap-1][pt-1] + delta_l;
	  double l_min_ltilde_PX = lambda_tilde_injected_PX[rap-1][pt-1] - delta_l;
	  double l_max_ltilde_PX = lambda_tilde_injected_PX[rap-1][pt-1] + delta_l;

	  TH1D* param_lth_CS = new TH1D( "param_lth_CS", "", int((l_max_lth_CS-l_min_lth_CS)/l_step_1D), l_min_lth_CS, l_max_lth_CS );
	  TH1D* param_lph_CS = new TH1D( "param_lph_CS", "", int((l_max_lph_CS-l_min_lph_CS)/l_step_1D), l_min_lph_CS, l_max_lph_CS );
	  TH1D* param_ltp_CS = new TH1D( "param_ltp_CS", "", int((l_max_ltp_CS-l_min_ltp_CS)/l_step_1D), l_min_ltp_CS, l_max_ltp_CS );
	  TH1D* param_ltilde_CS = new TH1D( "param_ltilde_CS", "", int((l_max_ltilde_CS-l_min_ltilde_CS)/l_step_1D), l_min_ltilde_CS, l_max_ltilde_CS );
	  TH1D* param_lphstar_CS = new TH1D( "param_lphstar_CS", "", int((l_max_lphstar_CS-l_min_lphstar_CS)/l_step_1D), l_min_lphstar_CS, l_max_lphstar_CS );
	  TH1D* param_lthstar_CS = new TH1D( "param_lthstar_CS", "", int((l_max_lthstar_CS-l_min_lthstar_CS)/l_step_1D), l_min_lthstar_CS, l_max_lthstar_CS );

	  TH1D* param_lth_HX = new TH1D( "param_lth_HX", "", int((l_max_lth_HX-l_min_lth_HX)/l_step_1D), l_min_lth_HX, l_max_lth_HX );
	  TH1D* param_lph_HX = new TH1D( "param_lph_HX", "", int((l_max_lph_HX-l_min_lph_HX)/l_step_1D), l_min_lph_HX, l_max_lph_HX );
	  TH1D* param_ltp_HX = new TH1D( "param_ltp_HX", "", int((l_max_ltp_HX-l_min_ltp_HX)/l_step_1D), l_min_ltp_HX, l_max_ltp_HX );
	  TH1D* param_ltilde_HX = new TH1D( "param_ltilde_HX", "", int((l_max_ltilde_HX-l_min_ltilde_HX)/l_step_1D), l_min_ltilde_HX, l_max_ltilde_HX );
	  TH1D* param_lphstar_HX = new TH1D( "param_lphstar_HX", "", int((l_max_lphstar_HX-l_min_lphstar_HX)/l_step_1D), l_min_lphstar_HX, l_max_lphstar_HX );
	  TH1D* param_lthstar_HX = new TH1D( "param_lthstar_HX", "", int((l_max_lthstar_HX-l_min_lthstar_HX)/l_step_1D), l_min_lthstar_HX, l_max_lthstar_HX );

	  TH1D* param_lth_PX = new TH1D( "param_lth_PX", "", int((l_max_lth_PX-l_min_lth_PX)/l_step_1D), l_min_lth_PX, l_max_lth_PX );
	  TH1D* param_lph_PX = new TH1D( "param_lph_PX", "", int((l_max_lph_PX-l_min_lph_PX)/l_step_1D), l_min_lph_PX, l_max_lph_PX );
	  TH1D* param_ltp_PX = new TH1D( "param_ltp_PX", "", int((l_max_ltp_PX-l_min_ltp_PX)/l_step_1D), l_min_ltp_PX, l_max_ltp_PX );
	  TH1D* param_ltilde_PX = new TH1D( "param_ltilde_PX", "", int((l_max_ltilde_PX-l_min_ltilde_PX)/l_step_1D), l_min_ltilde_PX, l_max_ltilde_PX );
	  TH1D* param_lphstar_PX = new TH1D( "param_lphstar_PX", "", int((l_max_lphstar_PX-l_min_lphstar_PX)/l_step_1D), l_min_lphstar_PX, l_max_lphstar_PX );
	  TH1D* param_lthstar_PX = new TH1D( "param_lthstar_PX", "", int((l_max_lthstar_PX-l_min_lthstar_PX)/l_step_1D), l_min_lthstar_PX, l_max_lthstar_PX );

	  TH1D* pull_lth_CS = new TH1D( "pull_lth_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lph_CS = new TH1D( "pull_lph_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_ltp_CS = new TH1D( "pull_ltp_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_ltilde_CS = new TH1D( "pull_ltilde_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lphstar_CS = new TH1D( "pull_lphstar_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lthstar_CS = new TH1D( "pull_lthstar_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );

	  TH1D* pull_lth_HX = new TH1D( "pull_lth_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lph_HX = new TH1D( "pull_lph_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_ltp_HX = new TH1D( "pull_ltp_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_ltilde_HX = new TH1D( "pull_ltilde_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lphstar_HX = new TH1D( "pull_lphstar_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lthstar_HX = new TH1D( "pull_lthstar_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );

	  TH1D* pull_lth_PX = new TH1D( "pull_lth_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lph_PX = new TH1D( "pull_lph_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_ltp_PX = new TH1D( "pull_ltp_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_ltilde_PX = new TH1D( "pull_ltilde_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lphstar_PX = new TH1D( "pull_lphstar_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
	  TH1D* pull_lthstar_PX = new TH1D( "pull_lthstar_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );


	  for(int iGen=1;iGen<nGenerations+1;iGen++){

	  sprintf(filename,"%s/rap%d_pT%d/Generation%d/results.root",dirstruct,iRap,iPt,iGen);
	  TFile* results = new TFile(filename);

	  if(results->Get("Results") == NULL) {results->Close(); continue;}

	  TTree* Results = (TTree*) results->Get("Results");


  	cout<<"";


	  double lth_CS; double lph_CS; double ltp_CS; double lthstar_CS; double lphstar_CS; double ltilde_CS;
	  double lth_HX; double lph_HX; double ltp_HX; double lthstar_HX; double lphstar_HX; double ltilde_HX;
	  double lth_PX; double lph_PX; double ltp_PX; double lthstar_PX; double lphstar_PX; double ltilde_PX;
	  int positivity_CS; int positivity_HX; int positivity_PX;

	  double err_lth_CS; double err_lph_CS; double err_ltp_CS; double err_lthstar_CS; double err_lphstar_CS; double err_ltilde_CS;
	  double err_lth_HX; double err_lph_HX; double err_ltp_HX; double err_lthstar_HX; double err_lphstar_HX; double err_ltilde_HX;
	  double err_lth_PX; double err_lph_PX; double err_ltp_PX; double err_lthstar_PX; double err_lphstar_PX; double err_ltilde_PX;

	  TBranch        *b_lth_CS;				Results->SetBranchAddress("lthCS", &lth_CS, &b_lth_CS);
	  TBranch        *b_lph_CS;				Results->SetBranchAddress("lphCS", &lph_CS, &b_lph_CS);
	  TBranch        *b_ltp_CS;				Results->SetBranchAddress("ltpCS", &ltp_CS, &b_ltp_CS);
	  TBranch        *b_lthstar_CS;			Results->SetBranchAddress("lthstarCS", &lthstar_CS, &b_lthstar_CS);
	  TBranch        *b_lphstar_CS;			Results->SetBranchAddress("lphstarCS", &lphstar_CS, &b_lphstar_CS);
	  TBranch        *b_ltilde_CS;			Results->SetBranchAddress("ltildeCS", &ltilde_CS, &b_ltilde_CS);
	  TBranch        *b_positivity_CS;		Results->SetBranchAddress("positivityCS", &positivity_CS, &b_positivity_CS);
	  TBranch        *b_lth_HX;				Results->SetBranchAddress("lthHX", &lth_HX, &b_lth_HX);
	  TBranch        *b_lph_HX;				Results->SetBranchAddress("lphHX", &lph_HX, &b_lph_HX);
	  TBranch        *b_ltp_HX;				Results->SetBranchAddress("ltpHX", &ltp_HX, &b_ltp_HX);
	  TBranch        *b_lthstar_HX;			Results->SetBranchAddress("lthstarHX", &lthstar_HX, &b_lthstar_HX);
	  TBranch        *b_lphstar_HX;			Results->SetBranchAddress("lphstarHX", &lphstar_HX, &b_lphstar_HX);
	  TBranch        *b_ltilde_HX;			Results->SetBranchAddress("ltildeHX", &ltilde_HX, &b_ltilde_HX);
	  TBranch        *b_positivity_HX;		Results->SetBranchAddress("positivityHX", &positivity_HX, &b_positivity_HX);
	  TBranch        *b_lth_PX;				Results->SetBranchAddress("lthPX", &lth_PX, &b_lth_PX);
	  TBranch        *b_lph_PX;				Results->SetBranchAddress("lphPX", &lph_PX, &b_lph_PX);
	  TBranch        *b_ltp_PX;				Results->SetBranchAddress("ltpPX", &ltp_PX, &b_ltp_PX);
	  TBranch        *b_lthstar_PX;			Results->SetBranchAddress("lthstarPX", &lthstar_PX, &b_lthstar_PX);
	  TBranch        *b_lphstar_PX;			Results->SetBranchAddress("lphstarPX", &lphstar_PX, &b_lphstar_PX);
	  TBranch        *b_ltilde_PX;			Results->SetBranchAddress("ltildePX", &ltilde_PX, &b_ltilde_PX);
	  TBranch        *b_positivity_PX;		Results->SetBranchAddress("positivityPX", &positivity_PX, &b_positivity_PX);

	  TBranch        *b_err_lth_CS;			Results->SetBranchAddress("err_lthCS",		 	&err_lth_CS, 			&b_err_lth_CS);
	  TBranch        *b_err_lph_CS;			Results->SetBranchAddress("err_lphCS",		 	&err_lph_CS, 			&b_err_lph_CS);
	  TBranch        *b_err_ltp_CS;			Results->SetBranchAddress("err_ltpCS",			&err_ltp_CS, 			&b_err_ltp_CS);
	  TBranch        *b_err_lthstar_CS;		Results->SetBranchAddress("err_lthstarCS",		&err_lthstar_CS, 		&b_err_lthstar_CS);
	  TBranch        *b_err_lphstar_CS;		Results->SetBranchAddress("err_lphstarCS",		&err_lphstar_CS, 		&b_err_lphstar_CS);
	  TBranch        *b_err_ltilde_CS;		Results->SetBranchAddress("err_ltildeCS",		&err_ltilde_CS, 		&b_err_ltilde_CS);
	  TBranch        *b_err_lth_HX;			Results->SetBranchAddress("err_lthHX", 			&err_lth_HX, 			&b_err_lth_HX);
	  TBranch        *b_err_lph_HX;			Results->SetBranchAddress("err_lphHX", 			&err_lph_HX, 			&b_err_lph_HX);
	  TBranch        *b_err_ltp_HX;			Results->SetBranchAddress("err_ltpHX", 			&err_ltp_HX, 			&b_err_ltp_HX);
	  TBranch        *b_err_lthstar_HX;		Results->SetBranchAddress("err_lthstarHX",		&err_lthstar_HX, 		&b_err_lthstar_HX);
	  TBranch        *b_err_lphstar_HX;		Results->SetBranchAddress("err_lphstarHX",		&err_lphstar_HX, 		&b_err_lphstar_HX);
	  TBranch        *b_err_ltilde_HX;		Results->SetBranchAddress("err_ltildeHX", 		&err_ltilde_HX, 		&b_err_ltilde_HX);
	  TBranch        *b_err_lth_PX;			Results->SetBranchAddress("err_lthPX", 			&err_lth_PX, 			&b_err_lth_PX);
	  TBranch        *b_err_lph_PX;			Results->SetBranchAddress("err_lphPX", 			&err_lph_PX, 			&b_err_lph_PX);
	  TBranch        *b_err_ltp_PX;			Results->SetBranchAddress("err_ltpPX", 			&err_ltp_PX, 			&b_err_ltp_PX);
	  TBranch        *b_err_lthstar_PX;		Results->SetBranchAddress("err_lthstarPX", 		&err_lthstar_PX, 		&b_err_lthstar_PX);
	  TBranch        *b_err_lphstar_PX;		Results->SetBranchAddress("err_lphstarPX", 		&err_lphstar_PX, 		&b_err_lphstar_PX);
	  TBranch        *b_err_ltilde_PX;		Results->SetBranchAddress("err_ltildePX", 		&err_ltilde_PX, 		&b_err_ltilde_PX);

	  Results->GetEvent( 0 );

	  	cout<<"";


	    param_lth_CS->Fill( lth_CS );
	    param_lph_CS->Fill( lph_CS );
	    param_ltp_CS->Fill( ltp_CS );
	    param_ltilde_CS->Fill( ltilde_CS );
	    param_lphstar_CS->Fill( lphstar_CS );
	    param_lthstar_CS->Fill( lthstar_CS );

	    param_lth_HX->Fill( lth_HX );
	    param_lph_HX->Fill( lph_HX );
	    param_ltp_HX->Fill( ltp_HX );
	    param_ltilde_HX->Fill( ltilde_HX );
	    param_lphstar_HX->Fill( lphstar_HX );
	    param_lthstar_HX->Fill( lthstar_HX );

	    param_lth_PX->Fill( lth_PX );
	    param_lph_PX->Fill( lph_PX );
	    param_ltp_PX->Fill( ltp_PX );
	    param_ltilde_PX->Fill( ltilde_PX );
	    param_lphstar_PX->Fill( lphstar_PX );
	    param_lthstar_PX->Fill( lthstar_PX );


	    pull_lth_CS->Fill( (lth_CS-lambda_theta_injected_CS[rap-1][pt-1])/err_lth_CS );
	    pull_lph_CS->Fill( (lph_CS-lambda_phi_injected_CS[rap-1][pt-1])/err_lph_CS );
	    pull_ltp_CS->Fill( (ltp_CS-lambda_thetaphi_injected_CS[rap-1][pt-1])/err_ltp_CS );
	    pull_lthstar_CS->Fill( (lthstar_CS-lambda_thetastar_injected_CS[rap-1][pt-1])/err_lthstar_CS );
	    pull_lphstar_CS->Fill( (lphstar_CS-lambda_phistar_injected_CS[rap-1][pt-1])/err_lphstar_CS );
	    pull_ltilde_CS->Fill( (ltilde_CS-lambda_tilde_injected_CS[rap-1][pt-1])/err_ltilde_CS );

	    pull_lth_HX->Fill( (lth_HX-lambda_theta_injected_HX[rap-1][pt-1])/err_lth_HX );
	    pull_lph_HX->Fill( (lph_HX-lambda_phi_injected_HX[rap-1][pt-1])/err_lph_HX );
	    pull_ltp_HX->Fill( (ltp_HX-lambda_thetaphi_injected_HX[rap-1][pt-1])/err_ltp_HX );
	    pull_lthstar_HX->Fill( (lthstar_HX-lambda_thetastar_injected_HX[rap-1][pt-1])/err_lthstar_HX );
	    pull_lphstar_HX->Fill( (lphstar_HX-lambda_phistar_injected_HX[rap-1][pt-1])/err_lphstar_HX );
	    pull_ltilde_HX->Fill( (ltilde_HX-lambda_tilde_injected_HX[rap-1][pt-1])/err_ltilde_HX );

	    pull_lth_PX->Fill( (lth_PX-lambda_theta_injected_PX[rap-1][pt-1])/err_lth_PX );
	    pull_lph_PX->Fill( (lph_PX-lambda_phi_injected_PX[rap-1][pt-1])/err_lph_PX );
	    pull_ltp_PX->Fill( (ltp_PX-lambda_thetaphi_injected_PX[rap-1][pt-1])/err_ltp_PX );
	    pull_lthstar_PX->Fill( (lthstar_PX-lambda_thetastar_injected_PX[rap-1][pt-1])/err_lthstar_PX );
	    pull_lphstar_PX->Fill( (lphstar_PX-lambda_phistar_injected_PX[rap-1][pt-1])/err_lphstar_PX );
	    pull_ltilde_PX->Fill( (ltilde_PX-lambda_tilde_injected_PX[rap-1][pt-1])/err_ltilde_PX );


	  results->Close();

	  }

	  emptyBin[rap-1][pt-1]=false;
	  if(param_lth_CS->GetMean()==0 || param_lth_HX->GetMean()==0 || param_lth_PX->GetMean()==0) emptyBin[rap-1][pt-1]=true;

	    pull_lth_CS_mean[rap-1][pt-1]=pull_lth_CS->GetMean(); pull_lth_CS_meanerr[rap-1][pt-1]=pull_lth_CS->GetMeanError(); pull_lth_CS_sigma[rap-1][pt-1]=pull_lth_CS->GetRMS(); pull_lth_CS_sigmaerr[rap-1][pt-1]=pull_lth_CS->GetRMSError();
	    pull_lph_CS_mean[rap-1][pt-1]=pull_lph_CS->GetMean(); pull_lph_CS_meanerr[rap-1][pt-1]=pull_lph_CS->GetMeanError(); pull_lph_CS_sigma[rap-1][pt-1]=pull_lph_CS->GetRMS(); pull_lph_CS_sigmaerr[rap-1][pt-1]=pull_lph_CS->GetRMSError();
	    pull_ltp_CS_mean[rap-1][pt-1]=pull_ltp_CS->GetMean(); pull_ltp_CS_meanerr[rap-1][pt-1]=pull_ltp_CS->GetMeanError(); pull_ltp_CS_sigma[rap-1][pt-1]=pull_ltp_CS->GetRMS(); pull_ltp_CS_sigmaerr[rap-1][pt-1]=pull_ltp_CS->GetRMSError();
	    pull_lthstar_CS_mean[rap-1][pt-1]=pull_lthstar_CS->GetMean(); pull_lthstar_CS_meanerr[rap-1][pt-1]=pull_lthstar_CS->GetMeanError(); pull_lthstar_CS_sigma[rap-1][pt-1]=pull_lthstar_CS->GetRMS(); pull_lthstar_CS_sigmaerr[rap-1][pt-1]=pull_lthstar_CS->GetRMSError();
	    pull_lphstar_CS_mean[rap-1][pt-1]=pull_lphstar_CS->GetMean(); pull_lphstar_CS_meanerr[rap-1][pt-1]=pull_lphstar_CS->GetMeanError(); pull_lphstar_CS_sigma[rap-1][pt-1]=pull_lphstar_CS->GetRMS(); pull_lphstar_CS_sigmaerr[rap-1][pt-1]=pull_lphstar_CS->GetRMSError();
	    pull_ltilde_CS_mean[rap-1][pt-1]=pull_ltilde_CS->GetMean(); pull_ltilde_CS_meanerr[rap-1][pt-1]=pull_ltilde_CS->GetMeanError(); pull_ltilde_CS_sigma[rap-1][pt-1]=pull_ltilde_CS->GetRMS(); pull_ltilde_CS_sigmaerr[rap-1][pt-1]=pull_ltilde_CS->GetRMSError();

	    pull_lth_HX_mean[rap-1][pt-1]=pull_lth_HX->GetMean(); pull_lth_HX_meanerr[rap-1][pt-1]=pull_lth_HX->GetMeanError(); pull_lth_HX_sigma[rap-1][pt-1]=pull_lth_HX->GetRMS(); pull_lth_HX_sigmaerr[rap-1][pt-1]=pull_lth_HX->GetRMSError();
	    pull_lph_HX_mean[rap-1][pt-1]=pull_lph_HX->GetMean(); pull_lph_HX_meanerr[rap-1][pt-1]=pull_lph_HX->GetMeanError(); pull_lph_HX_sigma[rap-1][pt-1]=pull_lph_HX->GetRMS(); pull_lph_HX_sigmaerr[rap-1][pt-1]=pull_lph_HX->GetRMSError();
	    pull_ltp_HX_mean[rap-1][pt-1]=pull_ltp_HX->GetMean(); pull_ltp_HX_meanerr[rap-1][pt-1]=pull_ltp_HX->GetMeanError(); pull_ltp_HX_sigma[rap-1][pt-1]=pull_ltp_HX->GetRMS(); pull_ltp_HX_sigmaerr[rap-1][pt-1]=pull_ltp_HX->GetRMSError();
	    pull_lthstar_HX_mean[rap-1][pt-1]=pull_lthstar_HX->GetMean(); pull_lthstar_HX_meanerr[rap-1][pt-1]=pull_lthstar_HX->GetMeanError(); pull_lthstar_HX_sigma[rap-1][pt-1]=pull_lthstar_HX->GetRMS(); pull_lthstar_HX_sigmaerr[rap-1][pt-1]=pull_lthstar_HX->GetRMSError();
	    pull_lphstar_HX_mean[rap-1][pt-1]=pull_lphstar_HX->GetMean(); pull_lphstar_HX_meanerr[rap-1][pt-1]=pull_lphstar_HX->GetMeanError(); pull_lphstar_HX_sigma[rap-1][pt-1]=pull_lphstar_HX->GetRMS(); pull_lphstar_HX_sigmaerr[rap-1][pt-1]=pull_lphstar_HX->GetRMSError();
	    pull_ltilde_HX_mean[rap-1][pt-1]=pull_ltilde_HX->GetMean(); pull_ltilde_HX_meanerr[rap-1][pt-1]=pull_ltilde_HX->GetMeanError(); pull_ltilde_HX_sigma[rap-1][pt-1]=pull_ltilde_HX->GetRMS(); pull_ltilde_HX_sigmaerr[rap-1][pt-1]=pull_ltilde_HX->GetRMSError();

	    pull_lth_PX_mean[rap-1][pt-1]=pull_lth_PX->GetMean(); pull_lth_PX_meanerr[rap-1][pt-1]=pull_lth_PX->GetMeanError(); pull_lth_PX_sigma[rap-1][pt-1]=pull_lth_PX->GetRMS(); pull_lth_PX_sigmaerr[rap-1][pt-1]=pull_lth_PX->GetRMSError();
	    pull_lph_PX_mean[rap-1][pt-1]=pull_lph_PX->GetMean(); pull_lph_PX_meanerr[rap-1][pt-1]=pull_lph_PX->GetMeanError(); pull_lph_PX_sigma[rap-1][pt-1]=pull_lph_PX->GetRMS(); pull_lph_PX_sigmaerr[rap-1][pt-1]=pull_lph_PX->GetRMSError();
	    pull_ltp_PX_mean[rap-1][pt-1]=pull_ltp_PX->GetMean(); pull_ltp_PX_meanerr[rap-1][pt-1]=pull_ltp_PX->GetMeanError(); pull_ltp_PX_sigma[rap-1][pt-1]=pull_ltp_PX->GetRMS(); pull_ltp_PX_sigmaerr[rap-1][pt-1]=pull_ltp_PX->GetRMSError();
	    pull_lthstar_PX_mean[rap-1][pt-1]=pull_lthstar_PX->GetMean(); pull_lthstar_PX_meanerr[rap-1][pt-1]=pull_lthstar_PX->GetMeanError(); pull_lthstar_PX_sigma[rap-1][pt-1]=pull_lthstar_PX->GetRMS(); pull_lthstar_PX_sigmaerr[rap-1][pt-1]=pull_lthstar_PX->GetRMSError();
	    pull_lphstar_PX_mean[rap-1][pt-1]=pull_lphstar_PX->GetMean(); pull_lphstar_PX_meanerr[rap-1][pt-1]=pull_lphstar_PX->GetMeanError(); pull_lphstar_PX_sigma[rap-1][pt-1]=pull_lphstar_PX->GetRMS(); pull_lphstar_PX_sigmaerr[rap-1][pt-1]=pull_lphstar_PX->GetRMSError();
	    pull_ltilde_PX_mean[rap-1][pt-1]=pull_ltilde_PX->GetMean(); pull_ltilde_PX_meanerr[rap-1][pt-1]=pull_ltilde_PX->GetMeanError(); pull_ltilde_PX_sigma[rap-1][pt-1]=pull_ltilde_PX->GetRMS(); pull_ltilde_PX_sigmaerr[rap-1][pt-1]=pull_ltilde_PX->GetRMSError();

	    param_lth_CS_mean[rap-1][pt-1]=param_lth_CS->GetMean(); param_lth_CS_meanerr[rap-1][pt-1]=param_lth_CS->GetMeanError(); param_lth_CS_sigma[rap-1][pt-1]=param_lth_CS->GetRMS(); param_lth_CS_sigmaerr[rap-1][pt-1]=param_lth_CS->GetRMSError();
	    param_lph_CS_mean[rap-1][pt-1]=param_lph_CS->GetMean(); param_lph_CS_meanerr[rap-1][pt-1]=param_lph_CS->GetMeanError(); param_lph_CS_sigma[rap-1][pt-1]=param_lph_CS->GetRMS(); param_lph_CS_sigmaerr[rap-1][pt-1]=param_lph_CS->GetRMSError();
	    param_ltp_CS_mean[rap-1][pt-1]=param_ltp_CS->GetMean(); param_ltp_CS_meanerr[rap-1][pt-1]=param_ltp_CS->GetMeanError(); param_ltp_CS_sigma[rap-1][pt-1]=param_ltp_CS->GetRMS(); param_ltp_CS_sigmaerr[rap-1][pt-1]=param_ltp_CS->GetRMSError();
	    param_lthstar_CS_mean[rap-1][pt-1]=param_lthstar_CS->GetMean(); param_lthstar_CS_meanerr[rap-1][pt-1]=param_lthstar_CS->GetMeanError(); param_lthstar_CS_sigma[rap-1][pt-1]=param_lthstar_CS->GetRMS(); param_lthstar_CS_sigmaerr[rap-1][pt-1]=param_lthstar_CS->GetRMSError();
	    param_lphstar_CS_mean[rap-1][pt-1]=param_lphstar_CS->GetMean(); param_lphstar_CS_meanerr[rap-1][pt-1]=param_lphstar_CS->GetMeanError(); param_lphstar_CS_sigma[rap-1][pt-1]=param_lphstar_CS->GetRMS(); param_lphstar_CS_sigmaerr[rap-1][pt-1]=param_lphstar_CS->GetRMSError();
	    param_ltilde_CS_mean[rap-1][pt-1]=param_ltilde_CS->GetMean(); param_ltilde_CS_meanerr[rap-1][pt-1]=param_ltilde_CS->GetMeanError(); param_ltilde_CS_sigma[rap-1][pt-1]=param_ltilde_CS->GetRMS(); param_ltilde_CS_sigmaerr[rap-1][pt-1]=param_ltilde_CS->GetRMSError();

	    param_lth_HX_mean[rap-1][pt-1]=param_lth_HX->GetMean(); param_lth_HX_meanerr[rap-1][pt-1]=param_lth_HX->GetMeanError(); param_lth_HX_sigma[rap-1][pt-1]=param_lth_HX->GetRMS(); param_lth_HX_sigmaerr[rap-1][pt-1]=param_lth_HX->GetRMSError();
	    param_lph_HX_mean[rap-1][pt-1]=param_lph_HX->GetMean(); param_lph_HX_meanerr[rap-1][pt-1]=param_lph_HX->GetMeanError(); param_lph_HX_sigma[rap-1][pt-1]=param_lph_HX->GetRMS(); param_lph_HX_sigmaerr[rap-1][pt-1]=param_lph_HX->GetRMSError();
	    param_ltp_HX_mean[rap-1][pt-1]=param_ltp_HX->GetMean(); param_ltp_HX_meanerr[rap-1][pt-1]=param_ltp_HX->GetMeanError(); param_ltp_HX_sigma[rap-1][pt-1]=param_ltp_HX->GetRMS(); param_ltp_HX_sigmaerr[rap-1][pt-1]=param_ltp_HX->GetRMSError();
	    param_lthstar_HX_mean[rap-1][pt-1]=param_lthstar_HX->GetMean(); param_lthstar_HX_meanerr[rap-1][pt-1]=param_lthstar_HX->GetMeanError(); param_lthstar_HX_sigma[rap-1][pt-1]=param_lthstar_HX->GetRMS(); param_lthstar_HX_sigmaerr[rap-1][pt-1]=param_lthstar_HX->GetRMSError();
	    param_lphstar_HX_mean[rap-1][pt-1]=param_lphstar_HX->GetMean(); param_lphstar_HX_meanerr[rap-1][pt-1]=param_lphstar_HX->GetMeanError(); param_lphstar_HX_sigma[rap-1][pt-1]=param_lphstar_HX->GetRMS(); param_lphstar_HX_sigmaerr[rap-1][pt-1]=param_lphstar_HX->GetRMSError();
	    param_ltilde_HX_mean[rap-1][pt-1]=param_ltilde_HX->GetMean(); param_ltilde_HX_meanerr[rap-1][pt-1]=param_ltilde_HX->GetMeanError(); param_ltilde_HX_sigma[rap-1][pt-1]=param_ltilde_HX->GetRMS(); param_ltilde_HX_sigmaerr[rap-1][pt-1]=param_ltilde_HX->GetRMSError();

	    param_lth_PX_mean[rap-1][pt-1]=param_lth_PX->GetMean(); param_lth_PX_meanerr[rap-1][pt-1]=param_lth_PX->GetMeanError(); param_lth_PX_sigma[rap-1][pt-1]=param_lth_PX->GetRMS(); param_lth_PX_sigmaerr[rap-1][pt-1]=param_lth_PX->GetRMSError();
	    param_lph_PX_mean[rap-1][pt-1]=param_lph_PX->GetMean(); param_lph_PX_meanerr[rap-1][pt-1]=param_lph_PX->GetMeanError(); param_lph_PX_sigma[rap-1][pt-1]=param_lph_PX->GetRMS(); param_lph_PX_sigmaerr[rap-1][pt-1]=param_lph_PX->GetRMSError();
	    param_ltp_PX_mean[rap-1][pt-1]=param_ltp_PX->GetMean(); param_ltp_PX_meanerr[rap-1][pt-1]=param_ltp_PX->GetMeanError(); param_ltp_PX_sigma[rap-1][pt-1]=param_ltp_PX->GetRMS(); param_ltp_PX_sigmaerr[rap-1][pt-1]=param_ltp_PX->GetRMSError();
	    param_lthstar_PX_mean[rap-1][pt-1]=param_lthstar_PX->GetMean(); param_lthstar_PX_meanerr[rap-1][pt-1]=param_lthstar_PX->GetMeanError(); param_lthstar_PX_sigma[rap-1][pt-1]=param_lthstar_PX->GetRMS(); param_lthstar_PX_sigmaerr[rap-1][pt-1]=param_lthstar_PX->GetRMSError();
	    param_lphstar_PX_mean[rap-1][pt-1]=param_lphstar_PX->GetMean(); param_lphstar_PX_meanerr[rap-1][pt-1]=param_lphstar_PX->GetMeanError(); param_lphstar_PX_sigma[rap-1][pt-1]=param_lphstar_PX->GetRMS(); param_lphstar_PX_sigmaerr[rap-1][pt-1]=param_lphstar_PX->GetRMSError();
	    param_ltilde_PX_mean[rap-1][pt-1]=param_ltilde_PX->GetMean(); param_ltilde_PX_meanerr[rap-1][pt-1]=param_ltilde_PX->GetMeanError(); param_ltilde_PX_sigma[rap-1][pt-1]=param_ltilde_PX->GetRMS(); param_ltilde_PX_sigmaerr[rap-1][pt-1]=param_ltilde_PX->GetRMSError();


		sprintf(filename,"%s/Figures",dirstruct); gSystem->mkdir(filename);

bool plotDist(true);
if(plotDist){

		 sprintf(filename,"%s/Figures/Pull_lth_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lth_CS, "z[#lambda^{CS}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lth_CS_mean[rap-1][pt-1], pull_lth_CS_meanerr[rap-1][pt-1], pull_lth_CS_sigma[rap-1][pt-1], pull_lth_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lph_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lph_CS, "z[#lambda^{CS}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lph_CS_mean[rap-1][pt-1], pull_lph_CS_meanerr[rap-1][pt-1], pull_lph_CS_sigma[rap-1][pt-1], pull_lph_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_ltp_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_ltp_CS, "z[#lambda^{CS}_{#theta#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltp_CS_mean[rap-1][pt-1], pull_ltp_CS_meanerr[rap-1][pt-1], pull_ltp_CS_sigma[rap-1][pt-1], pull_ltp_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lthstar_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lthstar_CS, "z[#lambda^{*CS}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lthstar_CS_mean[rap-1][pt-1], pull_lthstar_CS_meanerr[rap-1][pt-1], pull_lthstar_CS_sigma[rap-1][pt-1], pull_lthstar_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lphstar_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lphstar_CS, "z[#lambda^{*CS}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lphstar_CS_mean[rap-1][pt-1], pull_lphstar_CS_meanerr[rap-1][pt-1], pull_lphstar_CS_sigma[rap-1][pt-1], pull_lphstar_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_ltilde_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_ltilde_CS, "z[#tilde{#lambda}^{CS}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltilde_CS_mean[rap-1][pt-1], pull_ltilde_CS_meanerr[rap-1][pt-1], pull_ltilde_CS_sigma[rap-1][pt-1], pull_ltilde_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);


		 sprintf(filename,"%s/Figures/Pull_lth_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lth_HX, "z[#lambda^{HX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lth_HX_mean[rap-1][pt-1], pull_lth_HX_meanerr[rap-1][pt-1], pull_lth_HX_sigma[rap-1][pt-1], pull_lth_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lph_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lph_HX, "z[#lambda^{HX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lph_HX_mean[rap-1][pt-1], pull_lph_HX_meanerr[rap-1][pt-1], pull_lph_HX_sigma[rap-1][pt-1], pull_lph_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_ltp_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_ltp_HX, "z[#lambda^{HX}_{#theta#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltp_HX_mean[rap-1][pt-1], pull_ltp_HX_meanerr[rap-1][pt-1], pull_ltp_HX_sigma[rap-1][pt-1], pull_ltp_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lthstar_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lthstar_HX, "z[#lambda^{*HX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lthstar_HX_mean[rap-1][pt-1], pull_lthstar_HX_meanerr[rap-1][pt-1], pull_lthstar_HX_sigma[rap-1][pt-1], pull_lthstar_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lphstar_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lphstar_HX, "z[#lambda^{*HX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lphstar_HX_mean[rap-1][pt-1], pull_lphstar_HX_meanerr[rap-1][pt-1], pull_lphstar_HX_sigma[rap-1][pt-1], pull_lphstar_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_ltilde_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_ltilde_HX, "z[#tilde{#lambda}^{HX}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltilde_HX_mean[rap-1][pt-1], pull_ltilde_HX_meanerr[rap-1][pt-1], pull_ltilde_HX_sigma[rap-1][pt-1], pull_ltilde_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);


		 sprintf(filename,"%s/Figures/Pull_lth_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lth_PX, "z[#lambda^{PX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lth_PX_mean[rap-1][pt-1], pull_lth_PX_meanerr[rap-1][pt-1], pull_lth_PX_sigma[rap-1][pt-1], pull_lth_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lph_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lph_PX, "z[#lambda^{PX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lph_PX_mean[rap-1][pt-1], pull_lph_PX_meanerr[rap-1][pt-1], pull_lph_PX_sigma[rap-1][pt-1], pull_lph_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_ltp_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_ltp_PX, "z[#lambda^{PX}_{#theta#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltp_PX_mean[rap-1][pt-1], pull_ltp_PX_meanerr[rap-1][pt-1], pull_ltp_PX_sigma[rap-1][pt-1], pull_ltp_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lthstar_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lthstar_PX, "z[#lambda^{*PX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lthstar_PX_mean[rap-1][pt-1], pull_lthstar_PX_meanerr[rap-1][pt-1], pull_lthstar_PX_sigma[rap-1][pt-1], pull_lthstar_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_lphstar_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_lphstar_PX, "z[#lambda^{*PX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lphstar_PX_mean[rap-1][pt-1], pull_lphstar_PX_meanerr[rap-1][pt-1], pull_lphstar_PX_sigma[rap-1][pt-1], pull_lphstar_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);

		 sprintf(filename,"%s/Figures/Pull_ltilde_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(pull_ltilde_PX, "z[#tilde{#lambda}^{PX}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltilde_PX_mean[rap-1][pt-1], pull_ltilde_PX_meanerr[rap-1][pt-1], pull_ltilde_PX_sigma[rap-1][pt-1], pull_ltilde_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,true,0);




		 sprintf(filename,"%s/Figures/Param_lth_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lth_CS, "#lambda^{CS}_{#theta}", l_min_lth_CS, l_max_lth_CS, l_step_1D, param_lth_CS_mean[rap-1][pt-1], param_lth_CS_meanerr[rap-1][pt-1], param_lth_CS_sigma[rap-1][pt-1], param_lth_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_theta_injected_CS[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lph_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lph_CS, "#lambda^{CS}_{#phi}", l_min_lph_CS, l_max_lph_CS, l_step_1D, param_lph_CS_mean[rap-1][pt-1], param_lph_CS_meanerr[rap-1][pt-1], param_lph_CS_sigma[rap-1][pt-1], param_lph_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_phi_injected_CS[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_ltp_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_ltp_CS, "#lambda^{CS}_{#theta#phi}", l_min_ltp_CS, l_max_ltp_CS, l_step_1D, param_ltp_CS_mean[rap-1][pt-1], param_ltp_CS_meanerr[rap-1][pt-1], param_ltp_CS_sigma[rap-1][pt-1], param_ltp_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_thetaphi_injected_CS[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lthstar_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lthstar_CS, "#lambda^{*CS}_{#theta}", l_min_lthstar_CS, l_max_lthstar_CS, l_step_1D, param_lthstar_CS_mean[rap-1][pt-1], param_lthstar_CS_meanerr[rap-1][pt-1], param_lthstar_CS_sigma[rap-1][pt-1], param_lthstar_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_thetastar_injected_CS[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lphstar_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lphstar_CS, "#lambda^{*CS}_{#phi}", l_min_lphstar_CS, l_max_lphstar_CS, l_step_1D, param_lphstar_CS_mean[rap-1][pt-1], param_lphstar_CS_meanerr[rap-1][pt-1], param_lphstar_CS_sigma[rap-1][pt-1], param_lphstar_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_phistar_injected_CS[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_ltilde_CS_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_ltilde_CS, "#tilde{#lambda}^{CS}", l_min_ltilde_CS, l_max_ltilde_CS, l_step_1D, param_ltilde_CS_mean[rap-1][pt-1], param_ltilde_CS_meanerr[rap-1][pt-1], param_ltilde_CS_sigma[rap-1][pt-1], param_ltilde_CS_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_tilde_injected_CS[rap-1][pt-1]);


		 sprintf(filename,"%s/Figures/Param_lth_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lth_HX, "#lambda^{HX}_{#theta}", l_min_lth_HX, l_max_lth_HX, l_step_1D, param_lth_HX_mean[rap-1][pt-1], param_lth_HX_meanerr[rap-1][pt-1], param_lth_HX_sigma[rap-1][pt-1], param_lth_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_theta_injected_HX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lph_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lph_HX, "#lambda^{HX}_{#phi}", l_min_lph_HX, l_max_lph_HX, l_step_1D, param_lph_HX_mean[rap-1][pt-1], param_lph_HX_meanerr[rap-1][pt-1], param_lph_HX_sigma[rap-1][pt-1], param_lph_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_phi_injected_HX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_ltp_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_ltp_HX, "#lambda^{HX}_{#theta#phi}", l_min_ltp_HX, l_max_ltp_HX, l_step_1D, param_ltp_HX_mean[rap-1][pt-1], param_ltp_HX_meanerr[rap-1][pt-1], param_ltp_HX_sigma[rap-1][pt-1], param_ltp_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_thetaphi_injected_HX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lthstar_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lthstar_HX, "#lambda^{*HX}_{#theta}", l_min_lthstar_HX, l_max_lthstar_HX, l_step_1D, param_lthstar_HX_mean[rap-1][pt-1], param_lthstar_HX_meanerr[rap-1][pt-1], param_lthstar_HX_sigma[rap-1][pt-1], param_lthstar_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_thetastar_injected_HX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lphstar_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lphstar_HX, "#lambda^{*HX}_{#phi}", l_min_lphstar_HX, l_max_lphstar_HX, l_step_1D, param_lphstar_HX_mean[rap-1][pt-1], param_lphstar_HX_meanerr[rap-1][pt-1], param_lphstar_HX_sigma[rap-1][pt-1], param_lphstar_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_phistar_injected_HX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_ltilde_HX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_ltilde_HX, "#tilde{#lambda}^{HX}", l_min_ltilde_HX, l_max_ltilde_HX, l_step_1D, param_ltilde_HX_mean[rap-1][pt-1], param_ltilde_HX_meanerr[rap-1][pt-1], param_ltilde_HX_sigma[rap-1][pt-1], param_ltilde_HX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_tilde_injected_HX[rap-1][pt-1]);


		 sprintf(filename,"%s/Figures/Param_lth_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lth_PX, "#lambda^{PX}_{#theta}", l_min_lth_PX, l_max_lth_PX, l_step_1D, param_lth_PX_mean[rap-1][pt-1], param_lth_PX_meanerr[rap-1][pt-1], param_lth_PX_sigma[rap-1][pt-1], param_lth_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_theta_injected_PX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lph_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lph_PX, "#lambda^{PX}_{#phi}", l_min_lph_PX, l_max_lph_PX, l_step_1D, param_lph_PX_mean[rap-1][pt-1], param_lph_PX_meanerr[rap-1][pt-1], param_lph_PX_sigma[rap-1][pt-1], param_lph_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_phi_injected_PX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_ltp_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_ltp_PX, "#lambda^{PX}_{#theta#phi}", l_min_ltp_PX, l_max_ltp_PX, l_step_1D, param_ltp_PX_mean[rap-1][pt-1], param_ltp_PX_meanerr[rap-1][pt-1], param_ltp_PX_sigma[rap-1][pt-1], param_ltp_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_thetaphi_injected_PX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lthstar_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lthstar_PX, "#lambda^{*PX}_{#theta}", l_min_lthstar_PX, l_max_lthstar_PX, l_step_1D, param_lthstar_PX_mean[rap-1][pt-1], param_lthstar_PX_meanerr[rap-1][pt-1], param_lthstar_PX_sigma[rap-1][pt-1], param_lthstar_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_thetastar_injected_PX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_lphstar_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_lphstar_PX, "#lambda^{*PX}_{#phi}", l_min_lphstar_PX, l_max_lphstar_PX, l_step_1D, param_lphstar_PX_mean[rap-1][pt-1], param_lphstar_PX_meanerr[rap-1][pt-1], param_lphstar_PX_sigma[rap-1][pt-1], param_lphstar_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_phistar_injected_PX[rap-1][pt-1]);

		 sprintf(filename,"%s/Figures/Param_ltilde_PX_rap%dpt%d.pdf",dirstruct,iRap,iPt);
		 PlotObject(param_ltilde_PX, "#tilde{#lambda}^{PX}", l_min_ltilde_PX, l_max_ltilde_PX, l_step_1D, param_ltilde_PX_mean[rap-1][pt-1], param_ltilde_PX_meanerr[rap-1][pt-1], param_ltilde_PX_sigma[rap-1][pt-1], param_ltilde_PX_sigmaerr[rap-1][pt-1], filename, nGenerations,false,lambda_tilde_injected_PX[rap-1][pt-1]);


}

	delete param_lth_CS ;
	delete param_lph_CS ;
	delete param_ltp_CS ;
	delete param_ltilde_CS ;
	delete param_lphstar_CS ;
	delete param_lthstar_CS  ;
	delete param_lth_HX      ;
	delete param_lph_HX      ;
	delete param_ltp_HX      ;
	delete param_ltilde_HX   ;
	delete param_lphstar_HX  ;
	delete param_lthstar_HX  ;
	delete param_lth_PX      ;
	delete param_lph_PX      ;
	delete param_ltp_PX      ;
	delete param_ltilde_PX   ;
	delete param_lphstar_PX  ;
	delete param_lthstar_PX  ;
	delete pull_lth_CS         ;
	delete pull_lph_CS         ;
	delete pull_ltp_CS         ;
	delete pull_ltilde_CS      ;
	delete pull_lphstar_CS     ;
	delete pull_lthstar_CS     ;
	delete pull_lth_HX         ;
	delete pull_lph_HX         ;
	delete pull_ltp_HX         ;
	delete pull_ltilde_HX      ;
	delete pull_lphstar_HX     ;
	delete pull_lthstar_HX     ;
	delete pull_lth_PX         ;
	delete pull_lph_PX         ;
	delete pull_ltp_PX         ;
	delete pull_ltilde_PX      ;
	delete pull_lphstar_PX     ;
	delete pull_lthstar_PX     ;


	        }}

//////// Rap/Pt Summary Plots ////////////////


	    double lammin;
	    double lammax;

	    double pmumin=-5.1;
	    double pmumax=5.1;
	    double psigmin=0;
	    double psigmax=2;

	    double dmumin=-1.1;
	    double dmumax=1.1;
	    double dsigmin=0;
	    double dsigmax=0.25;
	    bool DrawLine(false);
	    double lamline=999;

	    double ddeltamumin=-delta_l;
	    double ddeltamumax=delta_l;

	    double KinDepRap[NptBins];
	    double KinDepRaperr[NptBins];
	    double KinDepRap2[NptBins];
	    double KinDepRap2err[NptBins];
	    double KinDepRap3[NptBins];
	    double KinDepRap3err[NptBins];

	    int nLamPlots=105;
	    char axislabel[200];

	    for(int iLam=1; iLam<nLamPlots+1;iLam++){

	    DrawLine=false;

	    if(iLam<13)continue;

		int rap=0;
	    for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++) {


	    	if(iLam==13){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lth_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{CS}_{#theta}");}
	    	if(iLam==14){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lth_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{CS}_{#theta}}");}
	    	if(iLam==15){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lph_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{CS}_{#phi}");}
	    	if(iLam==16){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lph_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{CS}_{#phi}}");}
	    	if(iLam==17){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_ltp_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{CS}_{#theta#phi}");}
	    	if(iLam==18){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_ltp_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{CS}_{#theta#phi}}");}
	    	if(iLam==19){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{*CS}_{#theta}");}
	    	if(iLam==20){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{*CS}_{#theta}}");}
	    	if(iLam==21){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{*CS}_{#phi}");}
	    	if(iLam==22){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{*CS}_{#phi}}");}
	    	if(iLam==23){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#tilde{#lambda}^{CS}");}
	    	if(iLam==24){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_ltilde_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#tilde{#lambda}^{CS}}");}

	    	if(iLam==25){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lth_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{HX}_{#theta}");}
	    	if(iLam==26){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lth_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{HX}_{#theta}}");}
	    	if(iLam==27){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lph_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{HX}_{#phi}");}
	    	if(iLam==28){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lph_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{HX}_{#phi}}");}
	    	if(iLam==29){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_ltp_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{HX}_{#theta#phi}");}
	    	if(iLam==30){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_ltp_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{HX}_{#theta#phi}}");}
	    	if(iLam==31){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{*HX}_{#theta}");}
	    	if(iLam==32){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{*HX}_{#theta}}");}
	    	if(iLam==33){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{*HX}_{#phi}");}
	    	if(iLam==34){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{*HX}_{#phi}}");}
	    	if(iLam==35){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#tilde{#lambda}^{HX}");}
	    	if(iLam==36){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_ltilde_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#tilde{#lambda}^{HX}}");}

	    	if(iLam==37){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lth_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{PX}_{#theta}");}
	    	if(iLam==38){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lth_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{PX}_{#theta}}");}
	    	if(iLam==39){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lph_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{PX}_{#phi}");}
	    	if(iLam==40){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lph_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{PX}_{#phi}}");}
	    	if(iLam==41){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_ltp_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{PX}_{#theta#phi}");}
	    	if(iLam==42){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_ltp_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{PX}_{#theta#phi}}");}
	    	if(iLam==43){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{*PX}_{#theta}");}
	    	if(iLam==44){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{*PX}_{#theta}}");}
	    	if(iLam==45){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#lambda^{*PX}_{#phi}");}
	    	if(iLam==46){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#lambda^{*PX}_{#phi}}");}
	    	if(iLam==47){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#tilde{#lambda}^{PX}");}
	    	if(iLam==48){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_ltilde_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{#tilde{#lambda}^{PX}}");}

	    	if(iLam==49){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lth_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{CS}_{#theta}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==50){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lth_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{CS}_{#theta}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==51){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lph_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{CS}_{#phi}]}");					DrawLine=true;lamline=0;}
	    	if(iLam==52){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lph_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{CS}_{#phi}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==53){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_ltp_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{CS}_{#theta#phi}]}");			DrawLine=true;lamline=0;}
	    	if(iLam==54){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_ltp_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{CS}_{#theta#phi}]}");			DrawLine=true;lamline=1;}
	    	if(iLam==55){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*CS}_{#theta}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==56){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*CS}_{#theta}]}");			DrawLine=true;lamline=1;}
	    	if(iLam==57){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*CS}_{#phi}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==58){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*CS}_{#phi}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==59){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#tilde{#lambda}^{CS}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==60){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#tilde{#lambda}^{CS}]}");				DrawLine=true;lamline=1;}

	    	if(iLam==61){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lth_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{HX}_{#theta}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==62){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lth_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{HX}_{#theta}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==63){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lph_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{HX}_{#phi}]}");					DrawLine=true;lamline=0;}
	    	if(iLam==64){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lph_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{HX}_{#phi}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==65){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_ltp_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{HX}_{#theta#phi}]}");			DrawLine=true;lamline=0;}
	    	if(iLam==66){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_ltp_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{HX}_{#theta#phi}]}");			DrawLine=true;lamline=1;}
	    	if(iLam==67){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*HX}_{#theta}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==68){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*HX}_{#theta}]}");			DrawLine=true;lamline=1;}
	    	if(iLam==69){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*HX}_{#phi}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==70){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*HX}_{#phi}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==71){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#tilde{#lambda}^{HX}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==72){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#tilde{#lambda}^{HX}]}");				DrawLine=true;lamline=1;}

	    	if(iLam==73){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lth_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{PX}_{#theta}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==74){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lth_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{PX}_{#theta}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==75){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lph_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{PX}_{#phi}]}");					DrawLine=true;lamline=0;}
	    	if(iLam==76){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lph_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{PX}_{#phi}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==77){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_ltp_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{PX}_{#theta#phi}]}");			DrawLine=true;lamline=0;}
	    	if(iLam==78){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_ltp_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{PX}_{#theta#phi}]}");			DrawLine=true;lamline=1;}
	    	if(iLam==79){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*PX}_{#theta}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==80){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*PX}_{#theta}]}");			DrawLine=true;lamline=1;}
	    	if(iLam==81){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*PX}_{#phi}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==82){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*PX}_{#phi}]}");				DrawLine=true;lamline=1;}
	    	if(iLam==83){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#mu_{[z[#tilde{#lambda}^{PX}]}");				DrawLine=true;lamline=0;}
	    	if(iLam==84){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#sigma_{z[#tilde{#lambda}^{PX}]}");				DrawLine=true;lamline=1;}


	    	if(iLam==85) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lth_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{CS}_{#theta}}");		DrawLine=true;lamline=0;}
	    	if(iLam==86) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lph_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{CS}_{#phi}}");			DrawLine=true;lamline=0;}
	    	if(iLam==87) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_ltp_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{CS}_{#theta#phi}}");	DrawLine=true;lamline=0;}
	    	if(iLam==88) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#lambda^{*CS}_{#theta}}");		DrawLine=true;lamline=0;}
	    	if(iLam==89) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#lambda^{*CS}_{#phi}}");			DrawLine=true;lamline=0;}
	    	if(iLam==90) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#tilde{#lambda}^{CS}}");			DrawLine=true;lamline=0;}

	    	if(iLam==91) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lth_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{HX}_{#theta}}");		DrawLine=true;lamline=0;}
	    	if(iLam==92) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lph_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{HX}_{#phi}}");			DrawLine=true;lamline=0;}
	    	if(iLam==93) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_ltp_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{HX}_{#theta#phi}}");	DrawLine=true;lamline=0;}
	    	if(iLam==94) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lthstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#lambda^{*HX}_{#theta}}");		DrawLine=true;lamline=0;}
	    	if(iLam==95) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#lambda^{*HX}_{#phi}}");			DrawLine=true;lamline=0;}
	    	if(iLam==96) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#tilde{#lambda}^{HX}}");			DrawLine=true;lamline=0;}

	    	if(iLam==97) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lth_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{PX}_{#theta}}");		DrawLine=true;lamline=0;}
	    	if(iLam==98) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lph_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{PX}_{#phi}}");			DrawLine=true;lamline=0;}
	    	if(iLam==99) {lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_ltp_rap%d.pdf",dirstruct,rapBin); 			sprintf(axislabel,"#delta_{#lambda^{PX}_{#theta#phi}}");	DrawLine=true;lamline=0;}
	    	if(iLam==100){lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#delta_{#lambda^{*PX}_{#theta}}");		DrawLine=true;lamline=0;}
	    	if(iLam==101){lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lphstar_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#lambda^{*PX}_{#phi}}");			DrawLine=true;lamline=0;}
	    	if(iLam==102){lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#delta_{#tilde{#lambda}^{PX}}");			DrawLine=true;lamline=0;}

	    	if(iLam==103){lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_mean_Comparison_ltilde_rap%d.pdf",dirstruct,rapBin); 		sprintf(axislabel,"#tilde{#lambda}"); lamline=lambda_tilde_injected_CS[rapBinMin-1][ptBinMin-1];}
	    	if(iLam==104){lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_mean_Comparison_lthstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#lambda*_{#theta}"); lamline=lambda_thetastar_injected_CS[rapBinMin-1][ptBinMin-1];}
	    	if(iLam==105){lammin=ddeltamumin;	lammax=ddeltamumax;	sprintf(filename,"%s/Figures/KinDep_param_mean_Comparison_lphstar_rap%d.pdf",dirstruct,rapBin); 	sprintf(axislabel,"#lambda*_{#phi}"); lamline=lambda_phistar_injected_CS[rapBinMin-1][ptBinMin-1];}

					int pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {



						if(iLam==13){ KinDepRap[pt]= param_lth_CS_mean[rap][pt]; 		KinDepRaperr[pt]=param_lth_CS_meanerr[rap][pt];}
						if(iLam==14){ KinDepRap[pt]= param_lth_CS_sigma[rap][pt]; 		KinDepRaperr[pt]=param_lth_CS_sigmaerr[rap][pt];}
						if(iLam==15){ KinDepRap[pt]= param_lph_CS_mean[rap][pt]; 		KinDepRaperr[pt]=param_lph_CS_meanerr[rap][pt];}
						if(iLam==16){ KinDepRap[pt]= param_lph_CS_sigma[rap][pt]; 		KinDepRaperr[pt]=param_lph_CS_sigmaerr[rap][pt];}
						if(iLam==17){ KinDepRap[pt]= param_ltp_CS_mean[rap][pt]; 		KinDepRaperr[pt]=param_ltp_CS_meanerr[rap][pt];}
						if(iLam==18){ KinDepRap[pt]= param_ltp_CS_sigma[rap][pt]; 		KinDepRaperr[pt]=param_ltp_CS_sigmaerr[rap][pt];}
						if(iLam==19){ KinDepRap[pt]= param_lthstar_CS_mean[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_CS_meanerr[rap][pt];}
						if(iLam==20){ KinDepRap[pt]= param_lthstar_CS_sigma[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_CS_sigmaerr[rap][pt];}
						if(iLam==21){ KinDepRap[pt]= param_lphstar_CS_mean[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_CS_meanerr[rap][pt];}
						if(iLam==22){ KinDepRap[pt]= param_lphstar_CS_sigma[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_CS_sigmaerr[rap][pt];}
						if(iLam==23){ KinDepRap[pt]= param_ltilde_CS_mean[rap][pt]; 	KinDepRaperr[pt]=param_ltilde_CS_meanerr[rap][pt];}
						if(iLam==24){ KinDepRap[pt]= param_ltilde_CS_sigma[rap][pt]; 	KinDepRaperr[pt]=param_ltilde_CS_sigmaerr[rap][pt];}

						if(iLam==25){ KinDepRap[pt]= param_lth_HX_mean[rap][pt]; 		KinDepRaperr[pt]=param_lth_HX_meanerr[rap][pt];}
						if(iLam==26){ KinDepRap[pt]= param_lth_HX_sigma[rap][pt]; 		KinDepRaperr[pt]=param_lth_HX_sigmaerr[rap][pt];}
						if(iLam==27){ KinDepRap[pt]= param_lph_HX_mean[rap][pt]; 		KinDepRaperr[pt]=param_lph_HX_meanerr[rap][pt];}
						if(iLam==28){ KinDepRap[pt]= param_lph_HX_sigma[rap][pt]; 		KinDepRaperr[pt]=param_lph_HX_sigmaerr[rap][pt];}
						if(iLam==29){ KinDepRap[pt]= param_ltp_HX_mean[rap][pt]; 		KinDepRaperr[pt]=param_ltp_HX_meanerr[rap][pt];}
						if(iLam==30){ KinDepRap[pt]= param_ltp_HX_sigma[rap][pt]; 		KinDepRaperr[pt]=param_ltp_HX_sigmaerr[rap][pt];}
						if(iLam==31){ KinDepRap[pt]= param_lthstar_HX_mean[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_HX_meanerr[rap][pt];}
						if(iLam==32){ KinDepRap[pt]= param_lthstar_HX_sigma[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_HX_sigmaerr[rap][pt];}
						if(iLam==33){ KinDepRap[pt]= param_lphstar_HX_mean[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_HX_meanerr[rap][pt];}
						if(iLam==34){ KinDepRap[pt]= param_lphstar_HX_sigma[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_HX_sigmaerr[rap][pt];}
						if(iLam==35){ KinDepRap[pt]= param_ltilde_HX_mean[rap][pt]; 	KinDepRaperr[pt]=param_ltilde_HX_meanerr[rap][pt];}
						if(iLam==36){ KinDepRap[pt]= param_ltilde_HX_sigma[rap][pt]; 	KinDepRaperr[pt]=param_ltilde_HX_sigmaerr[rap][pt];}

						if(iLam==37){ KinDepRap[pt]= param_lth_PX_mean[rap][pt]; 		KinDepRaperr[pt]=param_lth_PX_meanerr[rap][pt];}
						if(iLam==38){ KinDepRap[pt]= param_lth_PX_sigma[rap][pt]; 		KinDepRaperr[pt]=param_lth_PX_sigmaerr[rap][pt];}
						if(iLam==39){ KinDepRap[pt]= param_lph_PX_mean[rap][pt]; 		KinDepRaperr[pt]=param_lph_PX_meanerr[rap][pt];}
						if(iLam==40){ KinDepRap[pt]= param_lph_PX_sigma[rap][pt]; 		KinDepRaperr[pt]=param_lph_PX_sigmaerr[rap][pt];}
						if(iLam==41){ KinDepRap[pt]= param_ltp_PX_mean[rap][pt]; 		KinDepRaperr[pt]=param_ltp_PX_meanerr[rap][pt];}
						if(iLam==42){ KinDepRap[pt]= param_ltp_PX_sigma[rap][pt]; 		KinDepRaperr[pt]=param_ltp_PX_sigmaerr[rap][pt];}
						if(iLam==43){ KinDepRap[pt]= param_lthstar_PX_mean[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_PX_meanerr[rap][pt];}
						if(iLam==44){ KinDepRap[pt]= param_lthstar_PX_sigma[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_PX_sigmaerr[rap][pt];}
						if(iLam==45){ KinDepRap[pt]= param_lphstar_PX_mean[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_PX_meanerr[rap][pt];}
						if(iLam==46){ KinDepRap[pt]= param_lphstar_PX_sigma[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_PX_sigmaerr[rap][pt];}
						if(iLam==47){ KinDepRap[pt]= param_ltilde_PX_mean[rap][pt]; 	KinDepRaperr[pt]=param_ltilde_PX_meanerr[rap][pt];}
						if(iLam==48){ KinDepRap[pt]= param_ltilde_PX_sigma[rap][pt]; 	KinDepRaperr[pt]=param_ltilde_PX_sigmaerr[rap][pt];}

						if(iLam==49){ KinDepRap[pt]= pull_lth_CS_mean[rap][pt]; 		KinDepRaperr[pt]=pull_lth_CS_meanerr[rap][pt];}
						if(iLam==50){ KinDepRap[pt]= pull_lth_CS_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_lth_CS_sigmaerr[rap][pt];}
						if(iLam==51){ KinDepRap[pt]= pull_lph_CS_mean[rap][pt]; 		KinDepRaperr[pt]=pull_lph_CS_meanerr[rap][pt];}
						if(iLam==52){ KinDepRap[pt]= pull_lph_CS_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_lph_CS_sigmaerr[rap][pt];}
						if(iLam==53){ KinDepRap[pt]= pull_ltp_CS_mean[rap][pt]; 		KinDepRaperr[pt]=pull_ltp_CS_meanerr[rap][pt];}
						if(iLam==54){ KinDepRap[pt]= pull_ltp_CS_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_ltp_CS_sigmaerr[rap][pt];}
						if(iLam==55){ KinDepRap[pt]= pull_lthstar_CS_mean[rap][pt]; 	KinDepRaperr[pt]=pull_lthstar_CS_meanerr[rap][pt];}
						if(iLam==56){ KinDepRap[pt]= pull_lthstar_CS_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_lthstar_CS_sigmaerr[rap][pt];}
						if(iLam==57){ KinDepRap[pt]= pull_lphstar_CS_mean[rap][pt]; 	KinDepRaperr[pt]=pull_lphstar_CS_meanerr[rap][pt];}
						if(iLam==58){ KinDepRap[pt]= pull_lphstar_CS_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_lphstar_CS_sigmaerr[rap][pt];}
						if(iLam==59){ KinDepRap[pt]= pull_ltilde_CS_mean[rap][pt]; 		KinDepRaperr[pt]=pull_ltilde_CS_meanerr[rap][pt];}
						if(iLam==60){ KinDepRap[pt]= pull_ltilde_CS_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_ltilde_CS_sigmaerr[rap][pt];}

						if(iLam==61){ KinDepRap[pt]= pull_lth_HX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_lth_HX_meanerr[rap][pt];}
						if(iLam==62){ KinDepRap[pt]= pull_lth_HX_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_lth_HX_sigmaerr[rap][pt];}
						if(iLam==63){ KinDepRap[pt]= pull_lph_HX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_lph_HX_meanerr[rap][pt];}
						if(iLam==64){ KinDepRap[pt]= pull_lph_HX_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_lph_HX_sigmaerr[rap][pt];}
						if(iLam==65){ KinDepRap[pt]= pull_ltp_HX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_ltp_HX_meanerr[rap][pt];}
						if(iLam==66){ KinDepRap[pt]= pull_ltp_HX_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_ltp_HX_sigmaerr[rap][pt];}
						if(iLam==67){ KinDepRap[pt]= pull_lthstar_HX_mean[rap][pt]; 	KinDepRaperr[pt]=pull_lthstar_HX_meanerr[rap][pt];}
						if(iLam==68){ KinDepRap[pt]= pull_lthstar_HX_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_lthstar_HX_sigmaerr[rap][pt];}
						if(iLam==69){ KinDepRap[pt]= pull_lphstar_HX_mean[rap][pt]; 	KinDepRaperr[pt]=pull_lphstar_HX_meanerr[rap][pt];}
						if(iLam==70){ KinDepRap[pt]= pull_lphstar_HX_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_lphstar_HX_sigmaerr[rap][pt];}
						if(iLam==71){ KinDepRap[pt]= pull_ltilde_HX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_ltilde_HX_meanerr[rap][pt];}
						if(iLam==72){ KinDepRap[pt]= pull_ltilde_HX_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_ltilde_HX_sigmaerr[rap][pt];}

						if(iLam==73){ KinDepRap[pt]= pull_lth_PX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_lth_PX_meanerr[rap][pt];}
						if(iLam==74){ KinDepRap[pt]= pull_lth_PX_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_lth_PX_sigmaerr[rap][pt];}
						if(iLam==75){ KinDepRap[pt]= pull_lph_PX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_lph_PX_meanerr[rap][pt];}
						if(iLam==76){ KinDepRap[pt]= pull_lph_PX_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_lph_PX_sigmaerr[rap][pt];}
						if(iLam==77){ KinDepRap[pt]= pull_ltp_PX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_ltp_PX_meanerr[rap][pt];}
						if(iLam==78){ KinDepRap[pt]= pull_ltp_PX_sigma[rap][pt]; 		KinDepRaperr[pt]=pull_ltp_PX_sigmaerr[rap][pt];}
						if(iLam==79){ KinDepRap[pt]= pull_lthstar_PX_mean[rap][pt]; 	KinDepRaperr[pt]=pull_lthstar_PX_meanerr[rap][pt];}
						if(iLam==80){ KinDepRap[pt]= pull_lthstar_PX_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_lthstar_PX_sigmaerr[rap][pt];}
						if(iLam==81){ KinDepRap[pt]= pull_lphstar_PX_mean[rap][pt]; 	KinDepRaperr[pt]=pull_lphstar_PX_meanerr[rap][pt];}
						if(iLam==82){ KinDepRap[pt]= pull_lphstar_PX_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_lphstar_PX_sigmaerr[rap][pt];}
						if(iLam==83){ KinDepRap[pt]= pull_ltilde_PX_mean[rap][pt]; 		KinDepRaperr[pt]=pull_ltilde_PX_meanerr[rap][pt];}
						if(iLam==84){ KinDepRap[pt]= pull_ltilde_PX_sigma[rap][pt]; 	KinDepRaperr[pt]=pull_ltilde_PX_sigmaerr[rap][pt];}


						if(iLam==85){ KinDepRap[pt]= param_lth_CS_mean[rap][pt]-lambda_theta_injected_CS[rap][pt]; 			KinDepRaperr[pt]=param_lth_CS_sigma[rap][pt];}
						if(iLam==86){ KinDepRap[pt]= param_lph_CS_mean[rap][pt]-lambda_phi_injected_CS[rap][pt]; 			KinDepRaperr[pt]=param_lph_CS_sigma[rap][pt];}
						if(iLam==87){ KinDepRap[pt]= param_ltp_CS_mean[rap][pt]-lambda_thetaphi_injected_CS[rap][pt]; 		KinDepRaperr[pt]=param_ltp_CS_sigma[rap][pt];}
						if(iLam==88){ KinDepRap[pt]= param_lthstar_CS_mean[rap][pt]-lambda_thetastar_injected_CS[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_CS_sigma[rap][pt];}
						if(iLam==89){ KinDepRap[pt]= param_lphstar_CS_mean[rap][pt]-lambda_phistar_injected_CS[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_CS_sigma[rap][pt];}
						if(iLam==90){ KinDepRap[pt]= param_ltilde_CS_mean[rap][pt]-lambda_tilde_injected_CS[rap][pt]; 		KinDepRaperr[pt]=param_ltilde_CS_sigma[rap][pt];}

						if(iLam==91){ KinDepRap[pt]= param_lth_HX_mean[rap][pt]-lambda_theta_injected_HX[rap][pt]; 			KinDepRaperr[pt]=param_lth_HX_sigma[rap][pt];}
						if(iLam==92){ KinDepRap[pt]= param_lph_HX_mean[rap][pt]-lambda_phi_injected_HX[rap][pt]; 			KinDepRaperr[pt]=param_lph_HX_sigma[rap][pt];}
						if(iLam==93){ KinDepRap[pt]= param_ltp_HX_mean[rap][pt]-lambda_thetaphi_injected_HX[rap][pt]; 		KinDepRaperr[pt]=param_ltp_HX_sigma[rap][pt];}
						if(iLam==94){ KinDepRap[pt]= param_lthstar_HX_mean[rap][pt]-lambda_thetastar_injected_HX[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_HX_sigma[rap][pt];}
						if(iLam==95){ KinDepRap[pt]= param_lphstar_HX_mean[rap][pt]-lambda_phistar_injected_HX[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_HX_sigma[rap][pt];}
						if(iLam==96){ KinDepRap[pt]= param_ltilde_HX_mean[rap][pt]-lambda_tilde_injected_HX[rap][pt]; 		KinDepRaperr[pt]=param_ltilde_HX_sigma[rap][pt];}

						if(iLam==97){ KinDepRap[pt]= param_lth_PX_mean[rap][pt]-lambda_theta_injected_PX[rap][pt]; 			KinDepRaperr[pt]=param_lth_PX_sigma[rap][pt];}
						if(iLam==98){ KinDepRap[pt]= param_lph_PX_mean[rap][pt]-lambda_phi_injected_PX[rap][pt]; 			KinDepRaperr[pt]=param_lph_PX_sigma[rap][pt];}
						if(iLam==99){ KinDepRap[pt]= param_ltp_PX_mean[rap][pt]-lambda_thetaphi_injected_PX[rap][pt]; 		KinDepRaperr[pt]=param_ltp_PX_sigma[rap][pt];}
						if(iLam==100){ KinDepRap[pt]= param_lthstar_PX_mean[rap][pt]-lambda_thetastar_injected_PX[rap][pt]; 	KinDepRaperr[pt]=param_lthstar_PX_sigma[rap][pt];}
						if(iLam==101){ KinDepRap[pt]= param_lphstar_PX_mean[rap][pt]-lambda_phistar_injected_PX[rap][pt]; 	KinDepRaperr[pt]=param_lphstar_PX_sigma[rap][pt];}
						if(iLam==102){ KinDepRap[pt]= param_ltilde_PX_mean[rap][pt]-lambda_tilde_injected_PX[rap][pt]; 		KinDepRaperr[pt]=param_ltilde_PX_sigma[rap][pt];}

						if(iLam==103){  KinDepRap[pt]= param_ltilde_CS_mean[rap][pt]; 		 KinDepRaperr[pt]=param_ltilde_CS_sigma[rap][pt];	 KinDepRap2[pt]= param_ltilde_HX_mean[rap][pt]; 		 KinDepRap2err[pt]=param_ltilde_HX_sigma[rap][pt];	 KinDepRap3[pt]= param_ltilde_PX_mean[rap][pt]; 		 KinDepRap3err[pt]=param_ltilde_PX_sigma[rap][pt];}
						if(iLam==104){ KinDepRap[pt]= param_lthstar_CS_mean[rap][pt]; 		KinDepRaperr[pt]=param_lthstar_CS_sigma[rap][pt];	KinDepRap2[pt]= param_lthstar_HX_mean[rap][pt]; 		KinDepRap2err[pt]=param_lthstar_HX_sigma[rap][pt];	KinDepRap3[pt]= param_lthstar_PX_mean[rap][pt]; 		KinDepRap3err[pt]=param_lthstar_PX_sigma[rap][pt];}
						if(iLam==105){ KinDepRap[pt]= param_lphstar_CS_mean[rap][pt]; 		KinDepRaperr[pt]=param_lphstar_CS_sigma[rap][pt];	KinDepRap2[pt]= param_lphstar_HX_mean[rap][pt]; 		KinDepRap2err[pt]=param_lphstar_HX_sigma[rap][pt];	KinDepRap3[pt]= param_lphstar_PX_mean[rap][pt]; 		KinDepRap3err[pt]=param_lphstar_PX_sigma[rap][pt];}

						if(emptyBin[rap][pt]){ KinDepRap[pt]=999; KinDepRap2[pt]=999; KinDepRap3[pt]=999;}

						pt++;
					}


					if(iLam<103) PlotRapPt(ptBinMin, ptBinMax, rapBin, lammin, lammax, axislabel, filename, KinDepRap, KinDepRaperr, DrawLine, lamline);
					else PlotComparisonRapPt(ptBinMin, ptBinMax, rapBin, lammin, lammax, axislabel, filename, KinDepRap, KinDepRaperr, KinDepRap2, KinDepRap2err, KinDepRap3, KinDepRap3err, lamline);

					rap++;
	    }



	    }

///////////////// TABLE PRODUCTION /////////////////////////////////////////////////////////////////////////



		char NumFileName[200];
		sprintf(NumFileName,"ToyNumericalResults_%s.tex",dirstruct);
		FILE *NumFile = fopen(NumFileName,"w");

		fprintf(NumFile, "\n");
		fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\geometry{a4paper,left=20mm,right=20mm, top=1.5cm, bottom=1.5cm}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");





	    double lth_tab;
	    double ltherr_tab;
	    double lph_tab;
	    double lpherr_tab;
	    double ltp_tab;
	    double ltperr_tab;
	    double ltilde_tab;
	    double ltildeerr_tab;
	    double lthstar_tab;
	    double lthstarerr_tab;
	    double lphstar_tab;
	    double lphstarerr_tab;

	    char framerap[200];

		int nTables=3;

	    for(int iTab=1; iTab<nTables+1;iTab++){

			fprintf(NumFile, "\n\n\n\n");

			if(iTab==1){
				fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean of the distribution of the standard score $z(\\lambda_{i})$: $\\mu_{z(\\lambda_{i})} \\pm \\sigma_{\\mu_{z(\\lambda_{i})}}$}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
				fprintf(NumFile, "$p_{T}$ [GeV] & $\\mu[z(\\lambda_{\\vartheta})]$ & $\\mu[z(\\lambda_{\\varphi})]$ &  $\\mu[z(\\lambda_{\\vartheta \\varphi})]$ & $\\mu[z(\\tilde{\\lambda})]$ & $\\mu[z(\\lambda^*_{\\vartheta})]$ & $\\mu[z(\\lambda^*_{\\varphi})]$ \\\\\n");
			}
			if(iTab==2){
				fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{R.m.s. of the distribution of the standard score $z(\\lambda_{i})$: $\\sigma_{z(\\lambda_{i})} \\pm \\sigma_{\\sigma_{z(\\lambda_{i})}}$}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
				fprintf(NumFile, "$p_{T}$ [GeV] & $\\mu[z(\\lambda_{\\vartheta})]$ & $\\mu[z(\\lambda_{\\varphi})]$ &  $\\mu[z(\\lambda_{\\vartheta \\varphi})]$ & $\\mu[z(\\tilde{\\lambda})]$ & $\\mu[z(\\lambda^*_{\\vartheta})]$ & $\\mu[z(\\lambda^*_{\\varphi})]$ \\\\\n");
			}
			if(iTab==3){
				fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean Deviation of the distribution of the parameter estimates $\\lambda_{i}$: $\\delta_{\\lambda_{i}}=(\\mu_{\\lambda_{i}}-\\lambda^{Truth}_{i}) \\pm  \\sigma_{\\lambda_{i}}$}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
				fprintf(NumFile, "$p_{T}$ [GeV] & $\\mu[z(\\lambda_{\\vartheta})]$ & $\\mu[z(\\lambda_{\\varphi})]$ &  $\\mu[z(\\lambda_{\\vartheta \\varphi})]$ & $\\mu[z(\\tilde{\\lambda})]$ & $\\mu[z(\\lambda^*_{\\vartheta})]$ & $\\mu[z(\\lambda^*_{\\varphi})]$ \\\\\n");
			}

	    for(int iFrame=1; iFrame<4; iFrame++){


		int rap=0;
	    for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++) {

	    	if(iFrame==1) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{CS frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
	    	if(iFrame==2) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{HX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
	    	if(iFrame==3) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{PX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}

				int pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {


						if(iTab==1){
							if(iFrame==1){
								lth_tab=pull_lth_CS_mean[rap][pt]; ltherr_tab=pull_lth_CS_meanerr[rap][pt];
								lph_tab=pull_lph_CS_mean[rap][pt]; lpherr_tab=pull_lph_CS_meanerr[rap][pt];
								ltp_tab=pull_ltp_CS_mean[rap][pt]; ltperr_tab=pull_ltp_CS_meanerr[rap][pt];
								ltilde_tab=pull_ltilde_CS_mean[rap][pt]; ltildeerr_tab=pull_ltilde_CS_meanerr[rap][pt];
								lthstar_tab=pull_lthstar_CS_mean[rap][pt]; lthstarerr_tab=pull_lthstar_CS_meanerr[rap][pt];
								lphstar_tab=pull_lphstar_CS_mean[rap][pt]; lphstarerr_tab=pull_lphstar_CS_meanerr[rap][pt];
							}
							if(iFrame==2){
								lth_tab=pull_lth_HX_mean[rap][pt]; ltherr_tab=pull_lth_HX_meanerr[rap][pt];
								lph_tab=pull_lph_HX_mean[rap][pt]; lpherr_tab=pull_lph_HX_meanerr[rap][pt];
								ltp_tab=pull_ltp_HX_mean[rap][pt]; ltperr_tab=pull_ltp_HX_meanerr[rap][pt];
								ltilde_tab=pull_ltilde_HX_mean[rap][pt]; ltildeerr_tab=pull_ltilde_HX_meanerr[rap][pt];
								lthstar_tab=pull_lthstar_HX_mean[rap][pt]; lthstarerr_tab=pull_lthstar_HX_meanerr[rap][pt];
								lphstar_tab=pull_lphstar_HX_mean[rap][pt]; lphstarerr_tab=pull_lphstar_HX_meanerr[rap][pt];
							}
							if(iFrame==3){
								lth_tab=pull_lth_PX_mean[rap][pt]; ltherr_tab=pull_lth_PX_meanerr[rap][pt];
								lph_tab=pull_lph_PX_mean[rap][pt]; lpherr_tab=pull_lph_PX_meanerr[rap][pt];
								ltp_tab=pull_ltp_PX_mean[rap][pt]; ltperr_tab=pull_ltp_PX_meanerr[rap][pt];
								ltilde_tab=pull_ltilde_PX_mean[rap][pt]; ltildeerr_tab=pull_ltilde_PX_meanerr[rap][pt];
								lthstar_tab=pull_lthstar_PX_mean[rap][pt]; lthstarerr_tab=pull_lthstar_PX_meanerr[rap][pt];
								lphstar_tab=pull_lphstar_PX_mean[rap][pt]; lphstarerr_tab=pull_lphstar_PX_meanerr[rap][pt];
							}
						}

						if(iTab==2){
							if(iFrame==1){
								lth_tab=pull_lth_CS_sigma[rap][pt]; ltherr_tab=pull_lth_CS_sigmaerr[rap][pt];
								lph_tab=pull_lph_CS_sigma[rap][pt]; lpherr_tab=pull_lph_CS_sigmaerr[rap][pt];
								ltp_tab=pull_ltp_CS_sigma[rap][pt]; ltperr_tab=pull_ltp_CS_sigmaerr[rap][pt];
								ltilde_tab=pull_ltilde_CS_sigma[rap][pt]; ltildeerr_tab=pull_ltilde_CS_sigmaerr[rap][pt];
								lthstar_tab=pull_lthstar_CS_sigma[rap][pt]; lthstarerr_tab=pull_lthstar_CS_sigmaerr[rap][pt];
								lphstar_tab=pull_lphstar_CS_sigma[rap][pt]; lphstarerr_tab=pull_lphstar_CS_sigmaerr[rap][pt];
							}
							if(iFrame==2){
								lth_tab=pull_lth_HX_sigma[rap][pt]; ltherr_tab=pull_lth_HX_sigmaerr[rap][pt];
								lph_tab=pull_lph_HX_sigma[rap][pt]; lpherr_tab=pull_lph_HX_sigmaerr[rap][pt];
								ltp_tab=pull_ltp_HX_sigma[rap][pt]; ltperr_tab=pull_ltp_HX_sigmaerr[rap][pt];
								ltilde_tab=pull_ltilde_HX_sigma[rap][pt]; ltildeerr_tab=pull_ltilde_HX_sigmaerr[rap][pt];
								lthstar_tab=pull_lthstar_HX_sigma[rap][pt]; lthstarerr_tab=pull_lthstar_HX_sigmaerr[rap][pt];
								lphstar_tab=pull_lphstar_HX_sigma[rap][pt]; lphstarerr_tab=pull_lphstar_HX_sigmaerr[rap][pt];
							}
							if(iFrame==3){
								lth_tab=pull_lth_PX_sigma[rap][pt]; ltherr_tab=pull_lth_PX_sigmaerr[rap][pt];
								lph_tab=pull_lph_PX_sigma[rap][pt]; lpherr_tab=pull_lph_PX_sigmaerr[rap][pt];
								ltp_tab=pull_ltp_PX_sigma[rap][pt]; ltperr_tab=pull_ltp_PX_sigmaerr[rap][pt];
								ltilde_tab=pull_ltilde_PX_sigma[rap][pt]; ltildeerr_tab=pull_ltilde_PX_sigmaerr[rap][pt];
								lthstar_tab=pull_lthstar_PX_sigma[rap][pt]; lthstarerr_tab=pull_lthstar_PX_sigmaerr[rap][pt];
								lphstar_tab=pull_lphstar_PX_sigma[rap][pt]; lphstarerr_tab=pull_lphstar_PX_sigmaerr[rap][pt];
							}
						}

						if(iTab==3){
							if(iFrame==1){
								lth_tab=param_lth_CS_mean[rap][pt]			-lambda_theta_injected_CS[rap][pt];  	ltherr_tab=param_lth_CS_sigma[rap][pt];
								lph_tab=param_lph_CS_mean[rap][pt]			-lambda_phi_injected_CS[rap][pt]; 		lpherr_tab=param_lph_CS_sigma[rap][pt];
								ltp_tab=param_ltp_CS_mean[rap][pt]			-lambda_thetaphi_injected_CS[rap][pt]; 	ltperr_tab=param_ltp_CS_sigma[rap][pt];
								ltilde_tab=param_ltilde_CS_mean[rap][pt]	-lambda_tilde_injected_CS[rap][pt]; 	ltildeerr_tab=param_ltilde_CS_sigma[rap][pt];
								lthstar_tab=param_lthstar_CS_mean[rap][pt]	-lambda_thetastar_injected_CS[rap][pt];	lthstarerr_tab=param_lthstar_CS_sigma[rap][pt];
								lphstar_tab=param_lphstar_CS_mean[rap][pt]	-lambda_phistar_injected_CS[rap][pt];	lphstarerr_tab=param_lphstar_CS_sigma[rap][pt];
							}
							if(iFrame==2){
								lth_tab=param_lth_HX_mean[rap][pt]			-lambda_theta_injected_HX[rap][pt];  	ltherr_tab=param_lth_HX_sigma[rap][pt];
								lph_tab=param_lph_HX_mean[rap][pt]			-lambda_phi_injected_HX[rap][pt]; 		lpherr_tab=param_lph_HX_sigma[rap][pt];
								ltp_tab=param_ltp_HX_mean[rap][pt]			-lambda_thetaphi_injected_HX[rap][pt]; 	ltperr_tab=param_ltp_HX_sigma[rap][pt];
								ltilde_tab=param_ltilde_HX_mean[rap][pt]	-lambda_tilde_injected_HX[rap][pt]; 	ltildeerr_tab=param_ltilde_HX_sigma[rap][pt];
								lthstar_tab=param_lthstar_HX_mean[rap][pt]	-lambda_thetastar_injected_HX[rap][pt];	lthstarerr_tab=param_lthstar_HX_sigma[rap][pt];
								lphstar_tab=param_lphstar_HX_mean[rap][pt]	-lambda_phistar_injected_HX[rap][pt];	lphstarerr_tab=param_lphstar_HX_sigma[rap][pt];
							}
							if(iFrame==3){
								lth_tab=param_lth_PX_mean[rap][pt]			-lambda_theta_injected_PX[rap][pt];  	ltherr_tab=param_lth_PX_sigma[rap][pt];
								lph_tab=param_lph_PX_mean[rap][pt]			-lambda_phi_injected_PX[rap][pt]; 		lpherr_tab=param_lph_PX_sigma[rap][pt];
								ltp_tab=param_ltp_PX_mean[rap][pt]			-lambda_thetaphi_injected_PX[rap][pt]; 	ltperr_tab=param_ltp_PX_sigma[rap][pt];
								ltilde_tab=param_ltilde_PX_mean[rap][pt]	-lambda_tilde_injected_PX[rap][pt]; 	ltildeerr_tab=param_ltilde_PX_sigma[rap][pt];
								lthstar_tab=param_lthstar_PX_mean[rap][pt]	-lambda_thetastar_injected_PX[rap][pt];	lthstarerr_tab=param_lthstar_PX_sigma[rap][pt];
								lphstar_tab=param_lphstar_PX_mean[rap][pt]	-lambda_phistar_injected_PX[rap][pt];	lphstarerr_tab=param_lphstar_PX_sigma[rap][pt];
							}
						}


						if(emptyBin[rap][pt])fprintf(NumFile, "%1.0f--%1.0f   &  $nan$  & $nan$  &  $nan$ &  $nan$  & $nan$  &  $nan$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]);

						else fprintf(NumFile, "%1.0f--%1.0f   &  $%1.3f \\pm %1.3f$  & $%1.3f \\pm %1.3f$  &  $%1.3f \\pm %1.3f$ &  $%1.3f \\pm %1.3f$  & $%1.3f \\pm %1.3f$  &  $%1.3f \\pm %1.3f$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],lth_tab,ltherr_tab   ,lph_tab,lpherr_tab   ,ltp_tab,ltperr_tab   ,ltilde_tab,ltildeerr_tab   ,lthstar_tab,lthstarerr_tab   ,lphstar_tab,lphstarerr_tab);


						pt++;
					}



					rap++;

	    }//end rapBin

	    }//end iFrame


		fprintf(NumFile, "\\hline\n");
		fprintf(NumFile, "\\end{tabular}\n");
		fprintf(NumFile, "\\label{tab:syst_acceptance}\n");
		fprintf(NumFile, "\\end{table}\n");
		fprintf(NumFile, "\n");

	    }//end iTab


			fprintf(NumFile, "\\end{document}");

			fclose(NumFile);


















				return 0;
}
