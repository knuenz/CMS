#include <iostream>
#include <string>
#include <sstream>
using namespace std;


#include "rootIncludes.inc"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH1.h"

#include "TH2F.h"
#include "TTree.h"
#include "TDirectory.h"
#include "Riostream.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"

#include "TrimEventContent.C"

using namespace onia;

//====================================
int main(int argc, char** argv){

	Double_t fracL = 0.5;
	Double_t nSigma = 2.;
	bool UpsMC=false;
	bool f_BG_zero=false;
	bool ProjectLSBdata=false;
	bool ProjectRSBdata=false;
	bool CombineSignalPeaks=false;
	bool Y1Sto2S_SB=false;
	bool RightSided=false;
	bool LeftSided=false;
	bool MassScan=false;
	bool adjustOverlapBorders=true;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("FracLSB") != std::string::npos) {char* fracLchar = argv[i]; char* fracLchar2 = strtok (fracLchar, "p"); fracL = atof(fracLchar2); fracL=fracL/100.; cout<<"fracLSB = "<<fracL<<endl;}
		    if(std::string(argv[i]).find("nSigma") != std::string::npos) {char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "p"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl;}
		    if(std::string(argv[i]).find("UpsMC=1") != std::string::npos) {UpsMC=true;}
		    if(std::string(argv[i]).find("f_BG_zero=1") != std::string::npos) {f_BG_zero=true;}
		    if(std::string(argv[i]).find("ProjectLSBdata=1") != std::string::npos) {ProjectLSBdata=true;}
		    if(std::string(argv[i]).find("ProjectRSBdata=1") != std::string::npos) {ProjectRSBdata=true;}
		    if(std::string(argv[i]).find("CombineSignalPeaks=1") != std::string::npos) {CombineSignalPeaks=true;}
		    if(std::string(argv[i]).find("Y1Sto2S_SB=1") != std::string::npos) {Y1Sto2S_SB=true;}
		    if(std::string(argv[i]).find("LeftSided=1") != std::string::npos) {LeftSided=true;}
		    if(std::string(argv[i]).find("RightSided=1") != std::string::npos) {RightSided=true;}
		    if(std::string(argv[i]).find("MassScan=1") != std::string::npos) {MassScan=true;}
		    if(std::string(argv[i]).find("adjustOverlapBorders=0") != std::string::npos) {adjustOverlapBorders=false; cout<<"adjustOverlapBorders=false"<<endl;}
	    }

	  double fracLnS[3][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1];

	  double contamination2Sin1S_[onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double contamination1Sin2S_[onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double contamination3Sin2S_[onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double contamination2Sin3S_[onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double backgroundFrac_[3][onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double err_backgroundFrac_[3][onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double SignalEvents_[3][onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double ActualMassMin_[3][onia::kNbRapForPTBins][onia::kNbPTMaxBins];
	  double ActualMassMax_[3][onia::kNbRapForPTBins][onia::kNbPTMaxBins];

	  double ActualMassMin_Ups1S_rap1[onia::kNbPTMaxBins];
	  double ActualMassMin_Ups2S_rap1[onia::kNbPTMaxBins];
	  double ActualMassMin_Ups3S_rap1[onia::kNbPTMaxBins];
	  double ActualMassMax_Ups1S_rap1[onia::kNbPTMaxBins];
	  double ActualMassMax_Ups2S_rap1[onia::kNbPTMaxBins];
	  double ActualMassMax_Ups3S_rap1[onia::kNbPTMaxBins];
	  double ActualMassMin_Ups1S_rap2[onia::kNbPTMaxBins];
	  double ActualMassMin_Ups2S_rap2[onia::kNbPTMaxBins];
	  double ActualMassMin_Ups3S_rap2[onia::kNbPTMaxBins];
	  double ActualMassMax_Ups1S_rap2[onia::kNbPTMaxBins];
	  double ActualMassMax_Ups2S_rap2[onia::kNbPTMaxBins];
	  double ActualMassMax_Ups3S_rap2[onia::kNbPTMaxBins];

  gROOT->ProcessLine(".L TrimEventContent.C+");


  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 0; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 0; iPT < max; iPT++){
    	  fracLnS[iState][iRap][iPT]=TrimEventContent(iRap, iPT, fracL, nSigma, iState,UpsMC,f_BG_zero,ProjectLSBdata,ProjectRSBdata,CombineSignalPeaks,Y1Sto2S_SB,LeftSided,RightSided,MassScan,adjustOverlapBorders);
     	  if(iRap>0&&iPT>0){
     	  if(iState==0)contamination2Sin1S_[iRap-1][iPT-1]=contamination2Sin1S;
     	  if(iState==1)contamination1Sin2S_[iRap-1][iPT-1]=contamination1Sin2S;
     	  if(iState==1)contamination3Sin2S_[iRap-1][iPT-1]=contamination3Sin2S;
     	  if(iState==2)contamination2Sin3S_[iRap-1][iPT-1]=contamination2Sin3S;
    	  backgroundFrac_[iState][iRap-1][iPT-1]=backgroundFrac;
    	  err_backgroundFrac_[iState][iRap-1][iPT-1]=err_backgroundFrac;
    	  SignalEvents_[iState][iRap-1][iPT-1]=SignalEvents;
    	  ActualMassMin_[iState][iRap-1][iPT-1]=ActualMassMin;
    	  ActualMassMax_[iState][iRap-1][iPT-1]=ActualMassMax;
    	  cout<<"Number of signal events of Y"<<iState+1<<", rap"<<iRap<<", pT"<<iPT<<": "<< SignalEvents_[iState][iRap-1][iPT-1]<<endl;
    	  cout<<"Background fraction under  Y"<<iState+1<<", rap"<<iRap<<", pT"<<iPT<<": "<< backgroundFrac_[iState][iRap-1][iPT-1]<<" +- "<<err_backgroundFrac_[iState][iRap-1][iPT-1]<<endl;
    	  }
      }
    }
  }


  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 6; iPT < 11; iPT++){
    	  cout<<"Background fraction under  Y"<<iState+1<<", rap"<<iRap<<", pT"<<iPT<<": "<< backgroundFrac_[iState][iRap-1][iPT-1]<<" +- "<<err_backgroundFrac_[iState][iRap-1][iPT-1]<<" --> Relative uncertainty errfBG/fBG = "<<err_backgroundFrac_[iState][iRap-1][iPT-1]/backgroundFrac_[iState][iRap-1][iPT-1]*100.<<" %"<<endl;
      }
    }
  }


//  const int nBinspT=5;
//  double ptCentre_graph[nBinspT]={11.,14.,18.,25.,40.};

  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 6; iPT < 11; iPT++){
    	  cout<<"Ups"<<iState+1<<"S, rap"<<iRap<<", pT"<<iPT<<": MassMin = "<< ActualMassMin_[iState][iRap-1][iPT-1]<<", MassMax = "<<ActualMassMax_[iState][iRap-1][iPT-1]<<endl;
      }
    }
  }

  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
    	cout<<endl;
  	  cout<<"MassMin - Ups"<<iState+1<<"S, "<<"rap"<<iRap<<":"<<endl;
     Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 6; iPT < 11; iPT++){
    	  cout<<ActualMassMin_[iState][iRap-1][iPT-1]<<", ";
      }
    }
  }

  for(int iState = 0; iState < 3; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
    	cout<<endl;
  	  cout<<"MassMax - Ups"<<iState+1<<"S, "<<"rap"<<iRap<<":"<<endl;
     Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 6; iPT < 11; iPT++){
    	  cout<<ActualMassMax_[iState][iRap-1][iPT-1]<<", ";
      }
    }
  }

// Get all the values to plot

	int iPTmin=6;
	int iPTmax=onia::kNbPTMaxBins;
	int nBinspT=iPTmax-iPTmin+1;
  	int MarkerStyle[4][3]={{0,0,0},{0,33,27},{0,20,24},{0,21,25}}; // for each state, rapBin (1= closed, 2=open)
  	int MarkerColor[4] = {0,600,632,418}; // for each frame
  	double MarkerSize[4][3]={{0,0,0},{0,2.75,2.75},{0,1.65,1.65},{0,1.65,1.65}};// for each state, rapBin

	double ptCentre_[nBinspT];
	double ptCentreErr_low[nBinspT];
	double ptCentreErr_high[nBinspT];

  double fracL1S_rap1[nBinspT];
  double fracL1S_rap2[nBinspT];
  double fracL2S_rap1[nBinspT];
  double fracL2S_rap2[nBinspT];
  double fracL3S_rap1[nBinspT];
  double fracL3S_rap2[nBinspT];

  double backgroundFrac_1S_rap1[nBinspT];
  double backgroundFrac_1S_rap2[nBinspT];
  double backgroundFrac_2S_rap1[nBinspT];
  double backgroundFrac_2S_rap2[nBinspT];
  double backgroundFrac_3S_rap1[nBinspT];
  double backgroundFrac_3S_rap2[nBinspT];

  double err_backgroundFrac_1S_rap1[nBinspT];
  double err_backgroundFrac_1S_rap2[nBinspT];
  double err_backgroundFrac_2S_rap1[nBinspT];
  double err_backgroundFrac_2S_rap2[nBinspT];
  double err_backgroundFrac_3S_rap1[nBinspT];
  double err_backgroundFrac_3S_rap2[nBinspT];

  double SignalEvents_1S_rap1[nBinspT];
  double SignalEvents_1S_rap2[nBinspT];
  double SignalEvents_2S_rap1[nBinspT];
  double SignalEvents_2S_rap2[nBinspT];
  double SignalEvents_3S_rap1[nBinspT];
  double SignalEvents_3S_rap2[nBinspT];

  double contamination2Sin1S_rap1[nBinspT];
  double contamination2Sin1S_rap2[nBinspT];
  double contamination1Sin2S_rap1[nBinspT];
  double contamination1Sin2S_rap2[nBinspT];
  double contamination3Sin2S_rap1[nBinspT];
  double contamination3Sin2S_rap2[nBinspT];
  double contamination2Sin3S_rap1[nBinspT];
  double contamination2Sin3S_rap2[nBinspT];

	double frac1Smean=0;
	double frac2Smean=0;
	double frac3Smean=0;

	int pt=0;
      for(int iPT = iPTmin; iPT < iPTmax+1; iPT++){
    	  fracL1S_rap1[pt]=fracLnS[0][1][iPT];
    	  fracL1S_rap2[pt]=fracLnS[0][2][iPT];
    	  fracL2S_rap1[pt]=fracLnS[1][1][iPT];
    	  fracL2S_rap2[pt]=fracLnS[1][2][iPT];
    	  fracL3S_rap1[pt]=fracLnS[2][1][iPT];
    	  fracL3S_rap2[pt]=fracLnS[2][2][iPT];

       	  backgroundFrac_1S_rap1[pt]=backgroundFrac_[0][0][iPT-1];
          backgroundFrac_1S_rap2[pt]=backgroundFrac_[0][1][iPT-1];
          backgroundFrac_2S_rap1[pt]=backgroundFrac_[1][0][iPT-1];
          backgroundFrac_2S_rap2[pt]=backgroundFrac_[1][1][iPT-1];
          backgroundFrac_3S_rap1[pt]=backgroundFrac_[2][0][iPT-1];
          backgroundFrac_3S_rap2[pt]=backgroundFrac_[2][1][iPT-1];

       	  err_backgroundFrac_1S_rap1[pt]=err_backgroundFrac_[0][0][iPT-1];
          err_backgroundFrac_1S_rap2[pt]=err_backgroundFrac_[0][1][iPT-1];
          err_backgroundFrac_2S_rap1[pt]=err_backgroundFrac_[1][0][iPT-1];
          err_backgroundFrac_2S_rap2[pt]=err_backgroundFrac_[1][1][iPT-1];
          err_backgroundFrac_3S_rap1[pt]=err_backgroundFrac_[2][0][iPT-1];
          err_backgroundFrac_3S_rap2[pt]=err_backgroundFrac_[2][1][iPT-1];

       	  SignalEvents_1S_rap1[pt]=SignalEvents_[0][0][iPT-1];
          SignalEvents_1S_rap2[pt]=SignalEvents_[0][1][iPT-1];
          SignalEvents_2S_rap1[pt]=SignalEvents_[1][0][iPT-1];
          SignalEvents_2S_rap2[pt]=SignalEvents_[1][1][iPT-1];
          SignalEvents_3S_rap1[pt]=SignalEvents_[2][0][iPT-1];
          SignalEvents_3S_rap2[pt]=SignalEvents_[2][1][iPT-1];

       	  contamination2Sin1S_rap1[pt]=contamination2Sin1S_[0][iPT-1];
          contamination2Sin1S_rap2[pt]=contamination2Sin1S_[1][iPT-1];
          contamination1Sin2S_rap1[pt]=contamination1Sin2S_[0][iPT-1];
          contamination1Sin2S_rap2[pt]=contamination1Sin2S_[1][iPT-1];
          contamination3Sin2S_rap1[pt]=contamination3Sin2S_[0][iPT-1];
          contamination3Sin2S_rap2[pt]=contamination3Sin2S_[1][iPT-1];
          contamination2Sin3S_rap1[pt]=contamination2Sin3S_[0][iPT-1];
          contamination2Sin3S_rap2[pt]=contamination2Sin3S_[1][iPT-1];

        ptCentre_[pt]=(onia::pTRange[1][iPT-1]+onia::pTRange[1][iPT])/2;
  		ptCentreErr_low[pt]=TMath::Abs((onia::pTRange[1][iPT-1]-onia::pTRange[1][iPT])/2);
  		ptCentreErr_high[pt]=TMath::Abs((onia::pTRange[1][iPT-1]-onia::pTRange[1][iPT])/2);
  		frac1Smean+=fracL1S_rap1[pt]+fracL1S_rap2[pt];
  		frac2Smean+=fracL2S_rap1[pt]+fracL2S_rap2[pt];
  		frac3Smean+=fracL3S_rap1[pt]+fracL3S_rap2[pt];
  		pt++;
      }

      frac1Smean=frac1Smean/nBinspT/2.;
      frac2Smean=frac2Smean/nBinspT/2.;
      frac3Smean=frac3Smean/nBinspT/2.;


// Plot FracL

	   gStyle->SetPadLeftMargin(0.175);

  	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

  	plotCanvas->SetFillColor(kWhite);
  	plotCanvas->SetGrid();
  	plotCanvas->GetFrame()->SetFillColor(kWhite);
  	plotCanvas->GetFrame()->SetBorderSize(0);
  	plotCanvas->SetRightMargin(0.05) ;

  	TH1F *plotHisto = new TH1F;
  	plotHisto = plotCanvas->DrawFrame(onia::pTRange[1][iPTmin-1],0,onia::pTRange[1][iPTmax],1);
  	plotHisto->SetXTitle("p_{T} [GeV]");
  	plotHisto->SetYTitle("frac_{LSB}");
  	plotHisto->GetYaxis()->SetTitleOffset(1.6);
//  	plotHisto->SetLeftMargin(0.3);


  	TLine* RefLine1S = new TLine( onia::pTRange[1][iPTmin-1], frac1Smean, onia::pTRange[1][iPTmax], frac1Smean );
  	RefLine1S->SetLineWidth( 2 );
  	RefLine1S->SetLineStyle( 2 );
  	RefLine1S->SetLineColor( MarkerColor[1] );
  	RefLine1S->Draw( "same" );
  	TLine* RefLine2S = new TLine( onia::pTRange[1][iPTmin-1], frac2Smean, onia::pTRange[1][iPTmax], frac2Smean );
  	RefLine2S->SetLineWidth( 2 );
  	RefLine2S->SetLineStyle( 2 );
  	RefLine2S->SetLineColor( MarkerColor[2] );
  	RefLine2S->Draw( "same" );
  	TLine* RefLine3S = new TLine( onia::pTRange[1][iPTmin-1], frac3Smean, onia::pTRange[1][iPTmax], frac3Smean );
  	RefLine3S->SetLineWidth( 2 );
  	RefLine3S->SetLineStyle( 2 );
  	RefLine3S->SetLineColor( MarkerColor[3] );
  	RefLine3S->Draw( "same" );

   	TGraphAsymmErrors *graph_fracL1S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,fracL1S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
   	graph_fracL1S_rap1->SetMarkerColor(MarkerColor[1]);
   	graph_fracL1S_rap1->SetLineColor(MarkerColor[1]);
   	graph_fracL1S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
   	graph_fracL1S_rap1->SetMarkerSize(MarkerSize[1][1]);
   	graph_fracL1S_rap1->Draw("P");
   	TGraphAsymmErrors *graph_fracL1S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,fracL1S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
   	graph_fracL1S_rap2->SetMarkerColor(MarkerColor[1]);
   	graph_fracL1S_rap2->SetLineColor(MarkerColor[1]);
   	graph_fracL1S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
   	graph_fracL1S_rap2->SetMarkerSize(MarkerSize[1][2]);
   	graph_fracL1S_rap2->Draw("P");

   	TGraphAsymmErrors *graph_fracL2S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,fracL2S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
   	graph_fracL2S_rap1->SetMarkerColor(MarkerColor[2]);
   	graph_fracL2S_rap1->SetLineColor(MarkerColor[2]);
   	graph_fracL2S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
   	graph_fracL2S_rap1->SetMarkerSize(MarkerSize[1][1]);
   	graph_fracL2S_rap1->Draw("P");
   	TGraphAsymmErrors *graph_fracL2S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,fracL2S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
   	graph_fracL2S_rap2->SetMarkerColor(MarkerColor[2]);
   	graph_fracL2S_rap2->SetLineColor(MarkerColor[2]);
   	graph_fracL2S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
   	graph_fracL2S_rap2->SetMarkerSize(MarkerSize[1][2]);
   	graph_fracL2S_rap2->Draw("P");

   	TGraphAsymmErrors *graph_fracL3S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,fracL3S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
   	graph_fracL3S_rap1->SetMarkerColor(MarkerColor[3]);
   	graph_fracL3S_rap1->SetLineColor(MarkerColor[3]);
   	graph_fracL3S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
   	graph_fracL3S_rap1->SetMarkerSize(MarkerSize[1][1]);
   	graph_fracL3S_rap1->Draw("P");
   	TGraphAsymmErrors *graph_fracL3S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,fracL3S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
   	graph_fracL3S_rap2->SetMarkerColor(MarkerColor[3]);
   	graph_fracL3S_rap2->SetLineColor(MarkerColor[3]);
   	graph_fracL3S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
   	graph_fracL3S_rap2->SetMarkerSize(MarkerSize[1][2]);
   	graph_fracL3S_rap2->Draw("P");

   	char meantext[200];
   	double xText=30;
   	sprintf(meantext,"mean frac_{LSB} Y(1S) = %1.3f",frac1Smean);
	TLatex *text1S = new TLatex(xText,frac1Smean+0.05,meantext);
	text1S->SetTextSize(0.03);
	text1S->Draw( "same" );
	sprintf(meantext,"mean frac_{LSB} Y(2S) = %1.3f",frac2Smean);
	TLatex *text2S = new TLatex(xText,frac2Smean+0.05,meantext);
	text2S->SetTextSize(0.03);
	text2S->Draw( "same" );
	sprintf(meantext,"mean frac_{LSB} Y(3S) = %1.3f",frac3Smean);
	TLatex *text3S = new TLatex(xText,frac3Smean+0.05,meantext);
	text3S->SetTextSize(0.03);
	text3S->Draw( "same" );

	TLegend* plotLegend=new TLegend(0.175,0.1,0.525,0.3);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.03);
	plotLegend->SetBorderSize(1);
	TLegend* plotLegend2=new TLegend(0.525,0.1,0.95,0.3);
	plotLegend2->SetFillColor(kWhite);
	plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.03);
	plotLegend2->SetBorderSize(1);

	char legendentry[200];
	sprintf(legendentry,"frac_{LSB} Y(1S), |y| < 0.6");
	plotLegend->AddEntry(graph_fracL1S_rap1,legendentry,"pl");
	sprintf(legendentry,"frac_{LSB} Y(1S), 0.6 < |y| < 1.2");
	plotLegend2->AddEntry(graph_fracL1S_rap2,legendentry,"pl");

	sprintf(legendentry,"frac_{LSB} Y(2S), |y| < 0.6");
	plotLegend->AddEntry(graph_fracL2S_rap1,legendentry,"pl");
	sprintf(legendentry,"frac_{LSB} Y(2S), 0.6 < |y| < 1.2");
	plotLegend2->AddEntry(graph_fracL2S_rap2,legendentry,"pl");

	sprintf(legendentry,"frac_{LSB} Y(3S), |y| < 0.6");
	plotLegend->AddEntry(graph_fracL3S_rap1,legendentry,"pl");
	sprintf(legendentry,"frac_{LSB} Y(3S), 0.6 < |y| < 1.2");
	plotLegend2->AddEntry(graph_fracL3S_rap2,legendentry,"pl");


	plotLegend->Draw();
	plotLegend2->Draw();

  	char filename[200];
  	sprintf(filename, "AllStates_%1.2fSigma_FracLSB%dPercent/Figures/fracL.pdf", nSigma, int(fracL*100));

  	plotCanvas->SaveAs(filename);
  	plotCanvas->Close();







// Plot Contaminations

  	for(int iRap=1;iRap<3;iRap++){

  	  	plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

  	  	plotCanvas->SetFillColor(kWhite);
  	  	plotCanvas->SetGrid();
  	  	plotCanvas->GetFrame()->SetFillColor(kWhite);
  	  	plotCanvas->GetFrame()->SetBorderSize(0);
  	  	plotCanvas->SetRightMargin(0.05) ;

  	  	plotHisto = new TH1F;
  	  	plotHisto = plotCanvas->DrawFrame(onia::pTRange[1][iPTmin-1],0,onia::pTRange[1][iPTmax],0.07);
  	  	plotHisto->SetXTitle("p_{T} [GeV]");
  	  	plotHisto->SetYTitle("Signal contamination");
  	  	plotHisto->GetYaxis()->SetTitleOffset(1.6);
//  	  	plotHisto->SetLeftMargin(0.3);



  	   	TGraphAsymmErrors *graph_Cont2Sin1S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination2Sin1S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
  	graph_Cont2Sin1S_rap1->SetMarkerColor(MarkerColor[1]);
  	graph_Cont2Sin1S_rap1->SetLineColor(MarkerColor[1]);
  	graph_Cont2Sin1S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	graph_Cont2Sin1S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	if(iRap==1) graph_Cont2Sin1S_rap1->Draw("P");
  	   	TGraphAsymmErrors *graph_Cont2Sin1S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination2Sin1S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
  	graph_Cont2Sin1S_rap2->SetMarkerColor(MarkerColor[1]);
  	graph_Cont2Sin1S_rap2->SetLineColor(MarkerColor[1]);
  	graph_Cont2Sin1S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	graph_Cont2Sin1S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	if(iRap==2) graph_Cont2Sin1S_rap2->Draw("P");

	   	TGraphAsymmErrors *graph_Cont1Sin2S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination1Sin2S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
	graph_Cont1Sin2S_rap1->SetMarkerColor(MarkerColor[2]);
	graph_Cont1Sin2S_rap1->SetLineColor(MarkerColor[2]);
	graph_Cont1Sin2S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
	graph_Cont1Sin2S_rap1->SetMarkerSize(MarkerSize[1][1]);
	if(iRap==1) graph_Cont1Sin2S_rap1->Draw("P");
	   	TGraphAsymmErrors *graph_Cont1Sin2S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination1Sin2S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
	graph_Cont1Sin2S_rap2->SetMarkerColor(MarkerColor[2]);
	graph_Cont1Sin2S_rap2->SetLineColor(MarkerColor[2]);
	graph_Cont1Sin2S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
	graph_Cont1Sin2S_rap2->SetMarkerSize(MarkerSize[1][2]);
	if(iRap==2) graph_Cont1Sin2S_rap2->Draw("P");

	   	TGraphAsymmErrors *graph_Cont3Sin2S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination3Sin2S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
	graph_Cont3Sin2S_rap1->SetMarkerColor(MarkerColor[3]);
	graph_Cont3Sin2S_rap1->SetLineColor(MarkerColor[3]);
	graph_Cont3Sin2S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
	graph_Cont3Sin2S_rap1->SetMarkerSize(MarkerSize[1][1]);
	if(iRap==1) graph_Cont3Sin2S_rap1->Draw("P");
	   	TGraphAsymmErrors *graph_Cont3Sin2S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination3Sin2S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
	graph_Cont3Sin2S_rap2->SetMarkerColor(MarkerColor[3]);
	graph_Cont3Sin2S_rap2->SetLineColor(MarkerColor[3]);
	graph_Cont3Sin2S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
	graph_Cont3Sin2S_rap2->SetMarkerSize(MarkerSize[1][2]);
	if(iRap==2) graph_Cont3Sin2S_rap2->Draw("P");

	   	TGraphAsymmErrors *graph_Cont2Sin3S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination2Sin3S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
	graph_Cont2Sin3S_rap1->SetMarkerColor(kMagenta);
	graph_Cont2Sin3S_rap1->SetLineColor(kMagenta);
	graph_Cont2Sin3S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
	graph_Cont2Sin3S_rap1->SetMarkerSize(MarkerSize[1][1]);
	if(iRap==1) graph_Cont2Sin3S_rap1->Draw("P");
	   	TGraphAsymmErrors *graph_Cont2Sin3S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,contamination2Sin3S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
	graph_Cont2Sin3S_rap2->SetMarkerColor(kMagenta);
	graph_Cont2Sin3S_rap2->SetLineColor(kMagenta);
	graph_Cont2Sin3S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
	graph_Cont2Sin3S_rap2->SetMarkerSize(MarkerSize[1][2]);
	if(iRap==2) graph_Cont2Sin3S_rap2->Draw("P");

  		plotLegend=new TLegend(0.55,0.65,0.95,0.9);
  		plotLegend->SetFillColor(kWhite);
  		plotLegend->SetTextFont(72);
  		plotLegend->SetTextSize(0.035);
  		plotLegend->SetBorderSize(1);

  		sprintf(legendentry,"Y(2S) in Y(1S) sample");
  		if(iRap==1) plotLegend->AddEntry(graph_Cont2Sin1S_rap1,legendentry,"pl");
  		sprintf(legendentry,"Y(2S) in Y(1S) sample");
  		if(iRap==2) plotLegend->AddEntry(graph_Cont2Sin1S_rap2,legendentry,"pl");

  		sprintf(legendentry,"Y(1S) in Y(2S) sample");
  		if(iRap==1) plotLegend->AddEntry(graph_Cont1Sin2S_rap1,legendentry,"pl");
  		sprintf(legendentry,"Y(1S) in Y(2S) sample");
  		if(iRap==2) plotLegend->AddEntry(graph_Cont1Sin2S_rap2,legendentry,"pl");

  		sprintf(legendentry,"Y(3S) in Y(2S) sample");
  		if(iRap==1) plotLegend->AddEntry(graph_Cont3Sin2S_rap1,legendentry,"pl");
  		sprintf(legendentry,"Y(3S) in Y(2S) sample");
  		if(iRap==2) plotLegend->AddEntry(graph_Cont3Sin2S_rap2,legendentry,"pl");

  		sprintf(legendentry,"Y(2S) in Y(3S) sample");
  		if(iRap==1) plotLegend->AddEntry(graph_Cont2Sin3S_rap1,legendentry,"pl");
  		sprintf(legendentry,"Y(2S) in Y(3S) sample");
  		if(iRap==2) plotLegend->AddEntry(graph_Cont2Sin3S_rap2,legendentry,"pl");

  		plotLegend->Draw();


  		char texTex[200];
  		if(iRap==1) sprintf(texTex,"      |y| < 0.6");
  		if(iRap==2) sprintf(texTex,"0.6 < |y| < 1.2");
		 TLatex *text = new TLatex(15,0.055,texTex);
  		 text->SetTextSize(0.035);
  		 text->Draw( "same" );

  	  	sprintf(filename, "AllStates_%1.2fSigma_FracLSB%dPercent/Figures/SignalContamination_%d.pdf", nSigma, int(fracL*100),iRap);
  	  	plotCanvas->SaveAs(filename);
  	  	plotCanvas->Close();

  	}



// Plot Background fraction


  	    	plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

  	    	plotCanvas->SetFillColor(kWhite);
  	    	plotCanvas->SetGrid();
  	    	plotCanvas->GetFrame()->SetFillColor(kWhite);
  	    	plotCanvas->GetFrame()->SetBorderSize(0);
  	    	plotCanvas->SetRightMargin(0.05) ;

  	    	plotHisto = new TH1F;
  	    	plotHisto = plotCanvas->DrawFrame(onia::pTRange[1][iPTmin-1],0,onia::pTRange[1][iPTmax],0.3);
  	    	plotHisto->SetXTitle("p_{T} [GeV]");
  	    	plotHisto->SetYTitle("Background fraction");
  	    	plotHisto->GetYaxis()->SetTitleOffset(1.6);
//  	    	plotHisto->SetLeftMargin(0.3);


  	     	TGraphAsymmErrors *graph_backgroundFrac_1S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,backgroundFrac_1S_rap1,ptCentreErr_low,ptCentreErr_high,err_backgroundFrac_1S_rap1,err_backgroundFrac_1S_rap1);
  	     	graph_backgroundFrac_1S_rap1->SetMarkerColor(MarkerColor[1]);
  	     	graph_backgroundFrac_1S_rap1->SetLineColor(MarkerColor[1]);
  	     	graph_backgroundFrac_1S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	     	graph_backgroundFrac_1S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	     	graph_backgroundFrac_1S_rap1->Draw("P");
  	     	TGraphAsymmErrors *graph_backgroundFrac_1S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,backgroundFrac_1S_rap2,ptCentreErr_low,ptCentreErr_high,err_backgroundFrac_1S_rap2,err_backgroundFrac_1S_rap2);
  	     	graph_backgroundFrac_1S_rap2->SetMarkerColor(MarkerColor[1]);
  	     	graph_backgroundFrac_1S_rap2->SetLineColor(MarkerColor[1]);
  	     	graph_backgroundFrac_1S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	     	graph_backgroundFrac_1S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	     	graph_backgroundFrac_1S_rap2->Draw("P");

  	     	TGraphAsymmErrors *graph_backgroundFrac_2S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,backgroundFrac_2S_rap1,ptCentreErr_low,ptCentreErr_high,err_backgroundFrac_2S_rap1,err_backgroundFrac_2S_rap1);
  	     	graph_backgroundFrac_2S_rap1->SetMarkerColor(MarkerColor[2]);
  	     	graph_backgroundFrac_2S_rap1->SetLineColor(MarkerColor[2]);
  	     	graph_backgroundFrac_2S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	     	graph_backgroundFrac_2S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	     	graph_backgroundFrac_2S_rap1->Draw("P");
  	     	TGraphAsymmErrors *graph_backgroundFrac_2S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,backgroundFrac_2S_rap2,ptCentreErr_low,ptCentreErr_high,err_backgroundFrac_2S_rap2,err_backgroundFrac_2S_rap2);
  	     	graph_backgroundFrac_2S_rap2->SetMarkerColor(MarkerColor[2]);
  	     	graph_backgroundFrac_2S_rap2->SetLineColor(MarkerColor[2]);
  	     	graph_backgroundFrac_2S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	     	graph_backgroundFrac_2S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	     	graph_backgroundFrac_2S_rap2->Draw("P");

  	     	TGraphAsymmErrors *graph_backgroundFrac_3S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,backgroundFrac_3S_rap1,ptCentreErr_low,ptCentreErr_high,err_backgroundFrac_3S_rap1,err_backgroundFrac_3S_rap1);
  	     	graph_backgroundFrac_3S_rap1->SetMarkerColor(MarkerColor[3]);
  	     	graph_backgroundFrac_3S_rap1->SetLineColor(MarkerColor[3]);
  	     	graph_backgroundFrac_3S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	     	graph_backgroundFrac_3S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	     	graph_backgroundFrac_3S_rap1->Draw("P");
  	     	TGraphAsymmErrors *graph_backgroundFrac_3S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,backgroundFrac_3S_rap2,ptCentreErr_low,ptCentreErr_high,err_backgroundFrac_3S_rap2,err_backgroundFrac_3S_rap2);
  	     	graph_backgroundFrac_3S_rap2->SetMarkerColor(MarkerColor[3]);
  	     	graph_backgroundFrac_3S_rap2->SetLineColor(MarkerColor[3]);
  	     	graph_backgroundFrac_3S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	     	graph_backgroundFrac_3S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	     	graph_backgroundFrac_3S_rap2->Draw("P");

   	  	plotLegend=new TLegend(0.55,0.6,0.95,0.9);
  	  	plotLegend->SetFillColor(kWhite);
  	  	plotLegend->SetTextFont(72);
  	  	plotLegend->SetTextSize(0.025);
  	  	plotLegend->SetBorderSize(1);

  	  	sprintf(legendentry,"Y(1S) sample, |y| < 0.6");
  	  	plotLegend->AddEntry(graph_backgroundFrac_1S_rap1,legendentry,"pl");
  	  	sprintf(legendentry,"Y(1S) sample, 0.6 < |y| < 1.2");
  	  	plotLegend->AddEntry(graph_backgroundFrac_1S_rap2,legendentry,"pl");

  	  	sprintf(legendentry,"Y(2S) sample, |y| < 0.6");
  	  	plotLegend->AddEntry(graph_backgroundFrac_2S_rap1,legendentry,"pl");
  	  	sprintf(legendentry,"Y(2S) sample, 0.6 < |y| < 1.2");
  	  	plotLegend->AddEntry(graph_backgroundFrac_2S_rap2,legendentry,"pl");

  	  	sprintf(legendentry,"Y(3S) sample, |y| < 0.6");
  	  	plotLegend->AddEntry(graph_backgroundFrac_3S_rap1,legendentry,"pl");
  	  	sprintf(legendentry,"Y(3S) sample, 0.6 < |y| < 1.2");
  	  	plotLegend->AddEntry(graph_backgroundFrac_3S_rap2,legendentry,"pl");


  	  	plotLegend->Draw();

  	    	sprintf(filename, "AllStates_%1.2fSigma_FracLSB%dPercent/Figures/fracBG.pdf", nSigma, int(fracL*100));

  	    	plotCanvas->SaveAs(filename);
  	    	plotCanvas->Close();






// Plot SignalEvents


  	    	  	    	plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

  	    	  	    	plotCanvas->SetFillColor(kWhite);
  	    	  	    	plotCanvas->SetGrid();
  	    	  	    	plotCanvas->GetFrame()->SetFillColor(kWhite);
  	    	  	    	plotCanvas->GetFrame()->SetBorderSize(0);
  	    	  	    	plotCanvas->SetRightMargin(0.05) ;

  	    	  	    	plotHisto = new TH1F;
  	    	  	    	plotHisto = plotCanvas->DrawFrame(onia::pTRange[1][iPTmin-1],0,onia::pTRange[1][iPTmax],50000);
  	    	  	    	plotHisto->SetXTitle("p_{T} [GeV]");
  	    	  	    	plotHisto->SetYTitle("Number of signal events");
  	    	  	    	plotHisto->GetYaxis()->SetTitleOffset(1.6);


  	    	  	     	TGraphAsymmErrors *graph_SignalEvents_1S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,SignalEvents_1S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
  	    	  	     	graph_SignalEvents_1S_rap1->SetMarkerColor(MarkerColor[1]);
  	    	  	     	graph_SignalEvents_1S_rap1->SetLineColor(MarkerColor[1]);
  	    	  	     	graph_SignalEvents_1S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	    	  	     	graph_SignalEvents_1S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	    	  	     	graph_SignalEvents_1S_rap1->Draw("P");
  	    	  	     	TGraphAsymmErrors *graph_SignalEvents_1S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,SignalEvents_1S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
  	    	  	     	graph_SignalEvents_1S_rap2->SetMarkerColor(MarkerColor[1]);
  	    	  	     	graph_SignalEvents_1S_rap2->SetLineColor(MarkerColor[1]);
  	    	  	     	graph_SignalEvents_1S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	    	  	     	graph_SignalEvents_1S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	    	  	     	graph_SignalEvents_1S_rap2->Draw("P");

  	    	  	     	TGraphAsymmErrors *graph_SignalEvents_2S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,SignalEvents_2S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
  	    	  	     	graph_SignalEvents_2S_rap1->SetMarkerColor(MarkerColor[2]);
  	    	  	     	graph_SignalEvents_2S_rap1->SetLineColor(MarkerColor[2]);
  	    	  	     	graph_SignalEvents_2S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	    	  	     	graph_SignalEvents_2S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	    	  	     	graph_SignalEvents_2S_rap1->Draw("P");
  	    	  	     	TGraphAsymmErrors *graph_SignalEvents_2S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,SignalEvents_2S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
  	    	  	     	graph_SignalEvents_2S_rap2->SetMarkerColor(MarkerColor[2]);
  	    	  	     	graph_SignalEvents_2S_rap2->SetLineColor(MarkerColor[2]);
  	    	  	     	graph_SignalEvents_2S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	    	  	     	graph_SignalEvents_2S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	    	  	     	graph_SignalEvents_2S_rap2->Draw("P");

  	    	  	     	TGraphAsymmErrors *graph_SignalEvents_3S_rap1 = new TGraphAsymmErrors(nBinspT,ptCentre_,SignalEvents_3S_rap1,ptCentreErr_low,ptCentreErr_high,0,0);
  	    	  	     	graph_SignalEvents_3S_rap1->SetMarkerColor(MarkerColor[3]);
  	    	  	     	graph_SignalEvents_3S_rap1->SetLineColor(MarkerColor[3]);
  	    	  	     	graph_SignalEvents_3S_rap1->SetMarkerStyle(MarkerStyle[1][1]);
  	    	  	     	graph_SignalEvents_3S_rap1->SetMarkerSize(MarkerSize[1][1]);
  	    	  	     	graph_SignalEvents_3S_rap1->Draw("P");
  	    	  	     	TGraphAsymmErrors *graph_SignalEvents_3S_rap2 = new TGraphAsymmErrors(nBinspT,ptCentre_,SignalEvents_3S_rap2,ptCentreErr_low,ptCentreErr_high,0,0);
  	    	  	     	graph_SignalEvents_3S_rap2->SetMarkerColor(MarkerColor[3]);
  	    	  	     	graph_SignalEvents_3S_rap2->SetLineColor(MarkerColor[3]);
  	    	  	     	graph_SignalEvents_3S_rap2->SetMarkerStyle(MarkerStyle[1][2]);
  	    	  	     	graph_SignalEvents_3S_rap2->SetMarkerSize(MarkerSize[1][2]);
  	    	  	     	graph_SignalEvents_3S_rap2->Draw("P");

  	    	     	plotLegend=new TLegend(0.55,0.6,0.95,0.9);
  	    	  	  	plotLegend->SetFillColor(kWhite);
  	    	  	  	plotLegend->SetTextFont(72);
  	    	  	  	plotLegend->SetTextSize(0.025);
  	    	  	  	plotLegend->SetBorderSize(1);

  	    	  	  	sprintf(legendentry,"Y(1S) sample, |y| < 0.6");
  	    	  	  	plotLegend->AddEntry(graph_SignalEvents_1S_rap1,legendentry,"pl");
  	    	  	  	sprintf(legendentry,"Y(1S) sample, 0.6 < |y| < 1.2");
  	    	  	  	plotLegend->AddEntry(graph_SignalEvents_1S_rap2,legendentry,"pl");

  	    	  	  	sprintf(legendentry,"Y(2S) sample, |y| < 0.6");
  	    	  	  	plotLegend->AddEntry(graph_SignalEvents_2S_rap1,legendentry,"pl");
  	    	  	  	sprintf(legendentry,"Y(2S) sample, 0.6 < |y| < 1.2");
  	    	  	  	plotLegend->AddEntry(graph_SignalEvents_2S_rap2,legendentry,"pl");

  	    	  	  	sprintf(legendentry,"Y(3S) sample, |y| < 0.6");
  	    	  	  	plotLegend->AddEntry(graph_SignalEvents_3S_rap1,legendentry,"pl");
  	    	  	  	sprintf(legendentry,"Y(3S) sample, 0.6 < |y| < 1.2");
  	    	  	  	plotLegend->AddEntry(graph_SignalEvents_3S_rap2,legendentry,"pl");


  	    	  	  	plotLegend->Draw();

  	    	  	    	sprintf(filename, "AllStates_%1.2fSigma_FracLSB%dPercent/Figures/SignalEvents.pdf", nSigma, int(fracL*100));

  	    	  	    	plotCanvas->SaveAs(filename);
  	    	  	    	plotCanvas->Close();






  	  	delete plotCanvas;


  	cout<<"fLSB 1S, integrated in rap, pT"<<fracLnS[0][0][0]<<endl;
  	cout<<"fLSB 2S, integrated in rap, pT"<<fracLnS[1][0][0]<<endl;
  	cout<<"fLSB 3S, integrated in rap, pT"<<fracLnS[2][0][0]<<endl;

cout<<"SigTry "<<SignalEvents_[0][0][5]<<endl;


		FILE *NumFile;

		bool SaveTables=true;
	    if(SaveTables){


	    char framerap[200];

		char NumFileName[200];
		sprintf(NumFileName,"AllStates_%1.2fSigma_FracLSB%dPercent/numEvent_fBG.tex", nSigma, int(fracL*100));
		NumFile = fopen(NumFileName,"w");

		fprintf(NumFile, "\n");
		fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

		fprintf(NumFile, "\n\n\n\n");

	    	int nTables=2;
		    for(int iTab=2; iTab<nTables+1;iTab++){

				fprintf(NumFile, "\n\n\n\n");

				if(iTab==2){
					fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Estimated number of signal events and background fraction, as a function of $p_{T}$, in a 1$\\sigma$ window around the pole mass}\n \\begin{tabular}{|c|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$p_{T}$ [GeV] & $N_{\\Upsilon(1S)}$ & $f_{BG}^{\\Upsilon(1S)} [\\%] $ & $N_{\\Upsilon(2S)}$ & $f_{BG}^{\\Upsilon(2S)} [\\%]$& $N_{\\Upsilon(3S)}$ & $f_{BG}^{\\Upsilon(3S)} [\\%]$ \\\\\n");
				}



			int rap=0;
		    for(int rapBin = 1; rapBin < 3; rapBin++) {

		    	sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{$%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);

					int pt=0;
						for(int ptBin = 6; ptBin < 11; ptBin++) {

							if(iTab==2){
								fprintf(NumFile, "%1.0f--%1.0f   &  $%1.0f  $  & $%1.1f $  &  $%1.0f $ &  $%1.1f $ &  $%1.0f$ &  $%1.1f $\\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],SignalEvents_[0][rapBin-1][ptBin-1] ,100*backgroundFrac_[0][rapBin-1][ptBin-1],SignalEvents_[1][rapBin-1][ptBin-1] ,100*backgroundFrac_[1][rapBin-1][ptBin-1],SignalEvents_[2][rapBin-1][ptBin-1] ,100*backgroundFrac_[2][rapBin-1][ptBin-1]);
								printf("%1.0f--%1.0f   &  $%1.0f  $  & $%1.1f $  &  $%1.0f $ &  $%1.1f $ &  $%1.0f$ &  $%1.1f $\\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],SignalEvents_[0][rapBin-1][ptBin-1] ,100*backgroundFrac_[0][rapBin-1][ptBin-1],SignalEvents_[1][rapBin-1][ptBin-1] ,100*backgroundFrac_[1][rapBin-1][ptBin-1],SignalEvents_[2][rapBin-1][ptBin-1] ,100*backgroundFrac_[2][rapBin-1][ptBin-1]);

}


							pt++;
						}



						rap++;

		    }//end rapBin




			fprintf(NumFile, "\\hline\n");
			fprintf(NumFile, "\\end{tabular}\n");
			if(iTab==2) fprintf(NumFile, "\\label{tab:NumEventsFracBG}\n");
			fprintf(NumFile, "\\end{table}\n");
			fprintf(NumFile, "\n");

		    }//end iTab



	    }

			fprintf(NumFile, "\\end{document}");

			fclose(NumFile);






}
