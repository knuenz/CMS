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

	char dirStruct[200];
	sprintf(dirStruct,"Plots/ParameterPlots");

	gSystem->mkdir(dirStruct);

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

	double plotvalues_[jpsi::kNbPTMaxBins];
	double errplotvalues_[jpsi::kNbPTMaxBins];
	double ptCentre[jpsi::kNbRapForPTBins][jpsi::kNbPTMaxBins]={
			{3,6.5,7.5,9,12.5,17.5,25},
			{2,5,6.5,7.5,9,12.5,17.5,25}
	};
	double ptMean[jpsi::kNbPTMaxBins]={-9999,-9999,-9999,-9999,-9999,-9999};
	double errptMean[jpsi::kNbPTMaxBins]={NULL};

	char legendrap[200];
	char Filename[200];
	char yAxis[200];
	double yMin;
	double yMax;

//	double plotvalues[2][12];
//	double errplotvalues[2][12];

	for(int i=1;i<17;i++){


////////////////////////////// FILL INS /////////////////////////////
		double plotvalues[22][2][12]={
				{{-0.07240,0.11677,-0.14533,-0.29752,-0.04619,-9999.00000},{0.09032,0.11540,-0.02027,-0.07337,-0.13592,-0.10115,0.09064}},
				{{0.11646,-0.22619,-0.10316,0.14406,0.04314,-9999.00000},{0.00101,0.04654,0.34408,0.20649,0.24520,0.15611,0.18207}},
				{{-0.01986,0.01301,-0.01714,0.04978,-0.13467,-9999.00000},{0.03394,-0.12039,0.08061,-0.10259,0.11117,0.01947,0.13219}},
				{{0.16535,0.03471,0.09044,-0.19271,-0.31755,-9999.00000},{0.06379,0.14039,0.21751,0.17525,0.13289,-0.01600,0.03353}},
				{{-0.13676,0.11479,-0.16688,0.29657,-0.20683,-9999.00000},{0.05664,-0.18328,-0.05448,-0.00938,0.00961,-0.08935,0.00800}},
				{{-0.06855,0.21371,0.17602,-0.14533,-0.14850,-9999.00000},{0.12266,-0.02246,-0.04664,-0.00511,-0.04197,0.04235,0.08434}},
				{{0.84285,1.12189,0.87579,1.06521,1.06482,-9999.00000},{0.87995,0.91128,1.03967,0.89443,0.91656,1.01205,0.92472}},
				{{0.95640,0.92226,1.03505,0.89825,0.96664,-9999.00000},{1.03002,1.06554,0.92142,0.92426,0.99351,1.08098,1.00305}},
				{{1.05581,0.98838,1.09055,0.69054,0.90197,-9999.00000},{1.04265,1.08439,1.09067,1.09411,1.07211,1.01879,0.85151}},
				{{0.94977,0.99739,0.87443,1.04632,0.99787,-9999.00000},{1.02570,1.11737,1.16007,1.22825,1.06840,1.11417,0.88095}},
				{{0.94789,1.01579,1.00604,0.88485,1.11341,-9999.00000},{0.93295,1.03066,1.04311,0.86823,1.05470,1.04909,1.01993}},
				{{1.16568,0.89397,0.99386,1.20574,0.95077,-9999.00000},{1.01179,0.91177,0.94019,0.88844,0.83012,1.04642,1.09326}},
				{{50,50,50,50,50,0},{50,50,50,50,50,50,50}},
				{{50,50,50,50,50,0},{50,50,50,50,50,50,50}},
				{{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000}},
				{{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000}}

	};

	double errplotvalues[22][2][12]={


			{{0.11993,0.14486,0.13455,0.12094,0.09807,0.11071},{0.16557,0.17942,0.11901,0.09310,0.09926,0.11261,0.10169}},
			 {{0.11919,0.15864,0.12385,0.15063,0.15057,0.00000},{0.12443,0.12886,0.14702,0.12648,0.12961,0.14311,0.13076}},
							{{0.13524,0.13042,0.14636,0.12702,0.13669,0.00000},{0.14565,0.15068,0.13029,0.13070,0.14049,0.15285,0.14184}},
							{{0.14930,0.13976,0.15420,0.09765,0.12755,0.00000},{0.14744,0.15333,0.15422,0.15471,0.15160,0.14406,0.12041}},
							{{0.13431,0.14104,0.12365,0.14795,0.14111,0.00000},{0.14504,0.15800,0.16404,0.17367,0.15108,0.15755,0.12457}},
							{{0.13404,0.14364,0.14226,0.12513,0.15744,0.00000},{0.13193,0.14575,0.14750,0.12278,0.14914,0.14835,0.14423}},
							{{0.16483,0.12642,0.14054,0.17049,0.13445,0.00000},{0.14307,0.12893,0.13295,0.12563,0.11739,0.14797,0.15459}},
							{{0.08424,0.11212,0.08754,0.10643,0.10633,0.00000},{0.08796,0.09104,0.10389,0.08943,0.09162,0.10115,0.09244}},
							{{0.09562,0.09221,0.10341,0.08981,0.09663,0.00000},{0.10296,0.10656,0.09212,0.09243,0.09932,0.10807,0.10028}},
							{{0.11652,0.08937,0.09936,0.12053,0.09501,0.00000},{0.10108,0.09116,0.09399,0.08882,0.08301,0.10462,0.10920}},
							{{0.09493,0.09963,0.08743,0.10459,0.09971,0.00000},{0.10255,0.11171,0.11598,0.12274,0.10682,0.11133,0.08804}},
							{{0.09477,0.10151,0.10056,0.08847,0.11128,0.00000},{0.09329,0.10301,0.10424,0.08676,0.10545,0.10483,0.10196}},
							{{0.10556,0.09882,0.10892,0.06903,0.09020,0.00000},{0.10418,0.10842,0.10904,0.10934,0.10712,0.10186,0.08513}},
							{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
							{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
							{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
							{{0,0,0,0,0,0},{0,0,0,0,0,0,0}}




	};


/////// mean /////
	if(i<6.5) {yMin=-1.5;yMax=1.5;}
////// sigma /////
	if(i>6.5 && i<12.5) {yMin=0.1;yMax=1.9;}
////// Conv counter ////
	if(i>12.5 && i<14.5) {yMin=20;yMax=55;}
	if(i>14.5 && i<16.5) {yMin=0;yMax=12;}
	if(i>16.5) {yMin=-1;yMax=1;}

	if(i==1) sprintf(yAxis,"mean  #lambda_{#phi CS}");
	if(i==2) sprintf(yAxis,"mean  #lambda_{#theta CS}");
	if(i==3) sprintf(yAxis,"mean  #lambda_{#theta #phi CS}");
	if(i==4) sprintf(yAxis,"mean  #lambda_{#phi HX}");
	if(i==5) sprintf(yAxis,"mean  #lambda_{#theta HX}");
	if(i==6) sprintf(yAxis,"mean  #lambda_{#theta #phi HX}");
	if(i==7) sprintf(yAxis,"sigma  #lambda_{#phi CS}");
	if(i==8) sprintf(yAxis,"sigma  #lambda_{#theta CS}");
	if(i==9) sprintf(yAxis,"sigma  #lambda_{#theta #phi CS}");
	if(i==10) sprintf(yAxis,"sigma  #lambda_{#phi HX}");
	if(i==11) sprintf(yAxis,"sigma  #lambda_{#theta HX}");
	if(i==12) sprintf(yAxis,"sigma  #lambda_{#theta #phi HX}");
	if(i==13) sprintf(yAxis,"Converged fits CS");
	if(i==14) sprintf(yAxis,"Converged fits HX");
	if(i==15) sprintf(yAxis,"Overflow CS");
	if(i==16) sprintf(yAxis,"Overflow HX");
	if(i==17) sprintf(yAxis,"#lambda_{#theta CS}");
	if(i==18) sprintf(yAxis,"#lambda_{#theta HX}");
	if(i==19) sprintf(yAxis,"#lambda_{#phi CS}");
	if(i==20) sprintf(yAxis,"#lambda_{#phi HX}");
	if(i==21) sprintf(yAxis,"#lambda_{#theta #phi CS}");
	if(i==22) sprintf(yAxis,"#lambda_{#theta #phi HX}");

//	if(i==15) sprintf(yAxis,"nEvents");
//	if(i==16) sprintf(yAxis,"nEvents");

	if(i==1) sprintf(Filename,"%s/pT_ToyMCmean_phi_CS.png",dirStruct);
	if(i==2) sprintf(Filename,"%s/pT_ToyMCmean_theta_CS.png",dirStruct);
	if(i==3) sprintf(Filename,"%s/pT_ToyMCmean_thetaphi_CS.png",dirStruct);
	if(i==4) sprintf(Filename,"%s/pT_ToyMCmean_phi_HX.png",dirStruct);
	if(i==5) sprintf(Filename,"%s/pT_ToyMCmean_theta_HX.png",dirStruct);
	if(i==6) sprintf(Filename,"%s/pT_ToyMCmean_thetaphi_HX.png",dirStruct);
	if(i==7) sprintf(Filename,"%s/pT_ToyMCsigma_phi_CS.png",dirStruct);
	if(i==8) sprintf(Filename,"%s/pT_ToyMCsigma_theta_CS.png",dirStruct);
	if(i==9) sprintf(Filename,"%s/pT_ToyMCsigma_thetaphi_CS.png",dirStruct);
	if(i==10) sprintf(Filename,"%s/pT_ToyMCsigma_phi_HX.png",dirStruct);
	if(i==11) sprintf(Filename,"%s/pT_ToyMCsigma_theta_HX.png",dirStruct);
	if(i==12) sprintf(Filename,"%s/pT_ToyMCsigma_thetaphi_HX.png",dirStruct);
	if(i==13) sprintf(Filename,"%s/pT_ToyMCconvCS.png",dirStruct);
	if(i==14) sprintf(Filename,"%s/pT_ToyMCconvHX.png",dirStruct);
	if(i==15) sprintf(Filename,"%s/pT_ToyMCoverflowCS.png",dirStruct);
	if(i==16) sprintf(Filename,"%s/pT_ToyMCoverflowHX.png",dirStruct);
	if(i==17) sprintf(Filename,"%s/pT_RealMClambda_theta_CS.png",dirStruct);
	if(i==18) sprintf(Filename,"%s/pT_RealMClambda_theta_HX.png",dirStruct);
	if(i==19) sprintf(Filename,"%s/pT_RealMClambda_phi_CS.png",dirStruct);
	if(i==20) sprintf(Filename,"%s/pT_RealMClambda_phi_HX.png",dirStruct);
	if(i==21) sprintf(Filename,"%s/pT_RealMClambda_thetaphi_CS.png",dirStruct);
	if(i==22) sprintf(Filename,"%s/pT_RealMClambda_thetaphi_HX.png",dirStruct);

/////////////////////////////////////////////////////////////////////////////////////////////////////

	TLegend* plotLegend=new TLegend(0.75,0.75,0.95,0.9);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.035);
	plotLegend->SetBorderSize(1);


	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	TH1F *plotHisto = new TH1F;
//	TH1F *
	plotHisto = plotCanvas->DrawFrame(4,yMin,30,yMax);
	plotHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	plotHisto->SetYTitle(yAxis);
	plotHisto->GetYaxis()->SetTitleOffset(1.5);

	for(int yBin = 0; yBin < 2; ++yBin) {
	    	for(int ptBin = 0; ptBin < 10; ++ptBin) {

		    				   ptMean[ptBin]=ptCentre[yBin][ptBin]+0.25*yBin;
	  					   errptMean[ptBin]=0;
	  					 plotvalues_[ptBin+1]=plotvalues[i-1][yBin][ptBin];
	  					errplotvalues_[ptBin+1]=errplotvalues[i-1][yBin][ptBin];
	  				//	cout<<"i="<<i<<" "<<plotvalues[i-1][yBin][ptBin]<<endl;
	  					  }

	TGraphErrors *plotGraph = new TGraphErrors(jpsi::kNbPTBins[yBin]+1,ptMean,plotvalues_,errptMean,errplotvalues_);
	plotGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
//	plotGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
	plotGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
	sprintf(legendrap,"Rapidity %d",yBin+1);
	plotLegend->AddEntry(plotGraph,legendrap,"ple");
	plotGraph->Draw("P");
	//delete plotGraph;

	  		}






	  		plotLegend->Draw();

	//  		TPad* PaintLines = new TPad;
	//  		PaintLines->PaintLine(0.1,0.4,0.8,0.5);
	//  		PaintLines->Draw();

	  		plotCanvas->SaveAs(Filename);

	  		plotCanvas->Close();

	  		delete plotLegend;
	  		delete plotCanvas;
//	  		delete plotHisto;
	}
	  	  return 0;
	  	}
