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

   gSystem->mkdir("Plots/ParameterPlots");

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

	for(int i=1;i<23;i++){


////////////////////////////// FILL INS /////////////////////////////
		double plotvalues[22][2][12]={
				{{0.04469,0.05873,-0.33852,0.20095,0.07878,0.13100},{0.23269,0.23359,-0.19622,0.01873,-0.08477,-0.10432,-0.02330}},
				{{-0.05447,0.03048,0.26528,-0.15412,0.00332,-0.02896},{-0.17737,-0.10829,-0.05041,-0.05292,-0.19150,-0.04838,-0.01421}},
				{{0.00207,0.03729,0.09803,0.03559,-0.20080,0.07553},{-0.13249,-0.13423,0.20015,0.04578,-0.07222,-0.03566,-0.15591}},
				{{0.14851,0.13650,0.13785,0.07573,0.00464,-0.09050},{0.02267,-0.09895,-0.03712,0.02372,-0.01600,0.00760,0.01024}},
				{{0.10635,-0.10146,0.07211,0.07939,-0.02314,-0.12347},{0.01227,0.17414,-0.00660,-0.14739,0.03183,-0.00872,0.03942}},
				{{0.20941,0.02991,0.07852,-0.10006,0.12297,0.03033},{-0.16939,0.11519,-0.06880,0.02963,0.01416,0.08868,0.09939}},
				{{1.07886,0.88612,1.09233,0.99002,0.97587,1.03858},{1.73674,1.35993,0.95447,0.92642,0.98765,1.02001,0.97442}},
				{{1.16743,0.83323,1.07672,1.04801,0.83224,1.05978},{0.99948,1.01897,0.97160,0.97400,0.96982,0.97708,0.96863}},
				{{1.03142,0.83823,0.82007,0.98876,0.93282,0.97877},{1.86219,1.30765,1.04801,1.07748,1.05945,0.95079,0.91756}},
				{{1.02117,1.03620,1.07911,1.06466,1.02775,0.94515},{0.91427,0.90533,1.14151,1.09185,1.05516,1.02425,1.07052}},
				{{1.27859,1.02200,0.84612,1.08334,0.99245,0.97848},{1.13069,0.90094,1.02395,0.95624,1.05047,0.93191,1.09460}},
				{{1.08111,1.10637,0.87653,0.94183,1.10609,1.00321},{1.13413,1.05203,1.03666,0.95170,0.88479,1.00900,0.94925}},
				{{78,41,70,67,100,89},{72,61,65,100,100,85,95}},
				{{79,79,92,54,100,97},{48,58,95,100,100,96,97}},
				{{0.00000,1.00000,0.00000,0.00000,0.00000,0.00000},{0.00000,0.00000,1.00000,0.00000,0.00000,0.00000,0.00000}},
				{{7.00000,11.00000,1.00000,0.00000,0.00000,0.00000},{8.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000}},

				{{0.0640,9999,9999,9999,9999,9999},{9999,9999,9999,9999,9999,9999,9999}},
				{{0.1771,9999,9999,0.1578,0.0849,9999},{9999,0.2963,0.1697,9999,0.0482,9999,0.2445}},
				{{-0.0654,9999,9999,9999,9999,9999},{9999,9999,9999,9999,9999,9999,9999}},
				{{-0.0002,9999,9999,0.0047,-0.0094,9999},{9999,-0.0083,-0.0193,9999,-0.0001,9999,-0.0287}},
				{{0.0548,9999,9999,9999,9999,9999},{9999,9999,9999,9999,9999,9999,9999}},
				{{0.0214,9999,9999,-0.0156,-0.0124,9999},{9999,0.0419,0.0034,9999,-0.0055,9999,-0.0063}}



	};

	double errplotvalues[22][2][12]={


			{{0.12287,0.14177,0.13141,0.12088,0.09805,0.11059},{0.20432,0.17383,0.11923,0.09308,0.09921,0.11123,0.09993}},
			{{0.13294,0.13332,0.12954,0.12797,0.08362,0.11291},{0.11772,0.13040,0.12138,0.09785,0.09743,0.10656,0.09934}},
			{{0.11747,0.13412,0.09868,0.12072,0.09363,0.10429},{0.21896,0.16728,0.13092,0.10824,0.10643,0.10369,0.09410}},
			{{0.12027,0.12648,0.11363,0.14614,0.10325,0.09643},{0.14443,0.11984,0.11766,0.10968,0.10600,0.10504,0.10920}},
			{{0.15054,0.12477,0.08916,0.14867,0.09970,0.09982},{0.17854,0.11926,0.10556,0.09607,0.10553,0.09558,0.11166}},
			{{0.12732,0.13506,0.09236,0.12928,0.11111,0.10234},{0.17908,0.13924,0.10687,0.09561,0.08890,0.10347,0.09675}},
			{{0.08692,0.10029,0.09292,0.08551,0.06935,0.07823},{0.14474,0.12307,0.08428,0.06582,0.07015,0.07866,0.07069}},
			{{0.09407,0.09432,0.09166,0.09052,0.05913,0.07987},{0.08328,0.09224,0.08585,0.06919,0.06889,0.07534,0.07021}},
			{{0.09001,0.09548,0.06532,0.09145,0.07854,0.07237},{0.12677,0.09852,0.07561,0.06758,0.06287,0.07319,0.06845}},
			{{0.08507,0.08952,0.08042,0.10337,0.07302,0.06820},{0.10218,0.08473,0.08324,0.07755,0.07498,0.07430,0.07722}},
			{{0.10652,0.08826,0.06304,0.10518,0.07052,0.07061},{0.12638,0.08436,0.07467,0.06793,0.07462,0.06759,0.07896}},
			{{0.08305,0.09488,0.06980,0.08535,0.06625,0.07377},{0.15511,0.11837,0.09259,0.07655,0.07527,0.07335,0.06653}},
			{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
			{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
			{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
			{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},


			{{0.0849,0,0,0,0,0},{0,0,0,0,0,0,0}},
			{{0.2180,0,0,0.0492,0.0637,0},{0,0.1090,0.0907,0,0.0493,0,0.1220}},
			{{0.0985,0,0,0,0,0},{0,0,0,0,0,0,0}},
			{{0.0098,0,0,0.0051,0.0111,0},{0,0.0093,0.0090,0,0.0085,0,0.0325}},
			{{0.0772,0,0,0,0,0},{0,0,0,0,0,0,0}},
			{{0.0305,0,0,0.0107,0.0169,0},{0,0.0187,0.0166,0,0.0123,0,0.0404}}



	};


/////// mean /////
	if(i<6.5) {yMin=-0.5;yMax=0.5;}
////// sigma /////
	if(i>6.5 && i<12.5) {yMin=0;yMax=2.5;}
////// Conv counter ////
	if(i>12.5 && i<14.5) {yMin=0;yMax=105;}
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

	if(i==1) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCmean_phi_CS.png");
	if(i==2) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCmean_theta_CS.png");
	if(i==3) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCmean_thetaphi_CS.png");
	if(i==4) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCmean_phi_HX.png");
	if(i==5) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCmean_theta_HX.png");
	if(i==6) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCmean_thetaphi_HX.png");
	if(i==7) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCsigma_phi_CS.png");
	if(i==8) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCsigma_theta_CS.png");
	if(i==9) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCsigma_thetaphi_CS.png");
	if(i==10) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCsigma_phi_HX.png");
	if(i==11) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCsigma_theta_HX.png");
	if(i==12) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCsigma_thetaphi_HX.png");
	if(i==13) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCconvCS.png");
	if(i==14) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCconvHX.png");
	if(i==15) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCoverflowCS.png");
	if(i==16) sprintf(Filename,"Plots/ParameterPlots/pT_ToyMCoverflowHX.png");
	if(i==17) sprintf(Filename,"Plots/ParameterPlots/pT_RealMClambda_theta_CS.png");
	if(i==18) sprintf(Filename,"Plots/ParameterPlots/pT_RealMClambda_theta_HX.png");
	if(i==19) sprintf(Filename,"Plots/ParameterPlots/pT_RealMClambda_phi_CS.png");
	if(i==20) sprintf(Filename,"Plots/ParameterPlots/pT_RealMClambda_phi_HX.png");
	if(i==21) sprintf(Filename,"Plots/ParameterPlots/pT_RealMClambda_thetaphi_CS.png");
	if(i==22) sprintf(Filename,"Plots/ParameterPlots/pT_RealMClambda_thetaphi_HX.png");

/////////////////////////////////////////////////////////////////////////////////////////////////////

	TLegend* plotLegend=new TLegend(0.875,0.7,1,0.89);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.02);
	plotLegend->SetBorderSize(0);


	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,700);
	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.15) ;

	TH1F *plotHisto = new TH1F;
//	TH1F *
	plotHisto = plotCanvas->DrawFrame(4,yMin,30,yMax);
	plotHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
	plotHisto->SetYTitle(yAxis);
	plotHisto->GetYaxis()->SetTitleOffset(1.5);

	for(int yBin = 0; yBin < 2; ++yBin) {
	    	for(int ptBin = 0; ptBin < 10; ++ptBin) {

		    				   ptMean[ptBin]=ptCentre[yBin][ptBin];
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

	  		plotCanvas->SaveAs(Filename);

	  		plotCanvas->Close();

	  		delete plotLegend;
	  		delete plotCanvas;
//	  		delete plotHisto;
	}
	  	  return 0;
	  	}
