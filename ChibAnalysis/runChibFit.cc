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
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
//#include "ProfileLikelihoodCalculator.h"
//#include "HypoTestResult.h"
#include "BkgMixer.C"

double NormalizedIntegral(RooAbsPdf * function, RooRealVar& integrationVar, double lowerLimit, double upperLimit)
{
using namespace RooFit;
integrationVar.setRange("integralRange", lowerLimit, upperLimit);
RooAbsReal* integral = (*function).createIntegral(integrationVar,NormSet(integrationVar), Range("integralRange"));
double normalizedIntegralValue = integral->getVal();
delete integral;
return normalizedIntegralValue;
}
int main(int argc, char** argv) {

//	using namespace RooStats;

	Char_t *FitID = "Default"; //Storage Directory
	Char_t *cutName = "Default"; //Code Directory
  	Char_t *fileName = "Default";
  	int nState=1000;
  	char cutName_[500];
    int nToy=5000000;
    int nMix=5000000;

  	bool nYsigCut=false;
  	bool vtxProbCut=false;
  	bool ctSigCut=false;
  	bool gammaptCut=false;
  	bool ctCut=false;
  	bool RConvUpperCut=false;
  	bool RConvLowerCut=false;
  	bool pTCut=false;
  	bool vtxChi2ProbGammaCut=false;
  	bool vtxChi2ProbGammaLogCut=false;
    bool useSBforBKGmodel=false;
    bool SaveAll=false;
    bool alteredToy=false;
    bool SetSignalToZero=false;
    bool DrawTextOnPlots=false;
    bool useLeftSB=false;
    bool useRightSB=false;
    bool PlotStep1=false;
    bool PlotStep2=false;
    bool finalPlots=false;
    bool FixStuffForStep2=false;
    bool MCsample=false;
    bool useAnalyticalBKG=false;
    bool BkgMixerBool=false;
    bool BkgToy=false;
	bool useExistingMixFile=false;
	bool Pi0Cut=false;
	bool dzSigCut=false;
    bool OneSPlotModel=true;
    bool dzCut=false;
    bool restrictMaximum=false;
    bool NarrowOptimization=false;
    bool gammaptDCut=false;
    bool determine3PparametersDuringStep1=false;
    bool BkgToy1SBkgMixer2S=false;
    bool useNewExistingMixFile=false;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("FitID") != std::string::npos) {char* FitIDchar = argv[i]; char* FitIDchar2 = strtok (FitIDchar, "="); FitID = FitIDchar2; cout<<"FitID = "<<FitID<<endl;}
		    if(std::string(argv[i]).find("cutName") != std::string::npos) {char* cutNamechar = argv[i]; char* cutNamechar2 = strtok (cutNamechar, "="); cutName = cutNamechar2; cout<<"cutName = "<<cutName<<endl;}
		    if(std::string(argv[i]).find("fileName") != std::string::npos) {char* fileNamechar = argv[i]; char* fileNamechar2 = strtok (fileNamechar, "="); fileName = fileNamechar2; cout<<"fileName = "<<fileName<<endl;}
		    if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
		    if(std::string(argv[i]).find("nToy") != std::string::npos) {char* nToychar = argv[i]; char* nToychar2 = strtok (nToychar, "p"); nToy = atof(nToychar2); cout<<"nToy = "<<nToy<<endl;}
		    if(std::string(argv[i]).find("nMix") != std::string::npos) {char* nMixchar = argv[i]; char* nMixchar2 = strtok (nMixchar, "p"); nMix = atof(nMixchar2); cout<<"nMix = "<<nMix<<endl;}

		    if(std::string(argv[i]).find("nYsigCut") != std::string::npos) {nYsigCut=true; cout<<"Optimizing nYsigCut"<<endl;}
		    if(std::string(argv[i]).find("vtxProbCut") != std::string::npos) {vtxProbCut=true; cout<<"Optimizing vtxProbCut"<<endl;}
		    if(std::string(argv[i]).find("ctSigCut") != std::string::npos) {ctSigCut=true; cout<<"Optimizing ctSigCut"<<endl;}
		    if(std::string(argv[i]).find("ctCut") != std::string::npos) {ctCut=true; cout<<"Optimizing ctCut"<<endl;}
		    if(std::string(argv[i]).find("gammaptCut") != std::string::npos) {gammaptCut=true; cout<<"Optimizing gammaptCut"<<endl;}
		    if(std::string(argv[i]).find("gammaptDCut") != std::string::npos) {gammaptDCut=true; cout<<"Optimizing gammaptDCut"<<endl;}
		    if(std::string(argv[i]).find("RConvUpperCut") != std::string::npos) {RConvUpperCut=true; cout<<"Optimizing RConvUpperCut"<<endl;}
		    if(std::string(argv[i]).find("RConvLowerCut") != std::string::npos) {RConvLowerCut=true; cout<<"Optimizing RConvLowerCut"<<endl;}
		    if(std::string(argv[i]).find("pTCut") != std::string::npos) {pTCut=true; cout<<"Optimizing pTCut"<<endl;}
		    if(std::string(argv[i]).find("vtxChi2ProbGammaCut") != std::string::npos) {vtxChi2ProbGammaCut=true; cout<<"Optimizing vtxChi2ProbGammaCut"<<endl;}
		    if(std::string(argv[i]).find("vtxChi2ProbGammaLogCut") != std::string::npos) {vtxChi2ProbGammaLogCut=true; cout<<"Optimizing vtxChi2ProbGammaLogCut"<<endl;}
		    if(std::string(argv[i]).find("Pi0Cut") != std::string::npos) {Pi0Cut=true; cout<<"Optimizing Pi0Cut"<<endl;}
		    if(std::string(argv[i]).find("dzSigCut") != std::string::npos) {dzSigCut=true; cout<<"Optimizing dzSigCut"<<endl;}
		    if(std::string(argv[i]).find("dzCut") != std::string::npos) {dzCut=true; cout<<"Optimizing dzCut"<<endl;}
		    if(std::string(argv[i]).find("restrictMaximum=true") != std::string::npos) {restrictMaximum=true; cout<<"restrictMaximum"<<endl;}
		    if(std::string(argv[i]).find("NarrowOptimization=true") != std::string::npos) {NarrowOptimization=true; cout<<"NarrowOptimization"<<endl;}
		    if(std::string(argv[i]).find("determine3PparametersDuringStep1=true") != std::string::npos) {determine3PparametersDuringStep1=true; cout<<"determine3PparametersDuringStep1"<<endl;}

		    if(std::string(argv[i]).find("useSBforBKGmodel=true") != std::string::npos) {useSBforBKGmodel=true; cout<<"use SB for BKG model"<<endl;}
		    if(std::string(argv[i]).find("SaveAll=true") != std::string::npos) {SaveAll=true; cout<<"SaveAll"<<endl;}
		    if(std::string(argv[i]).find("SetSignalToZero=true") != std::string::npos) {SetSignalToZero=true; cout<<"SetSignalToZero"<<endl;}
		    if(std::string(argv[i]).find("DrawTextOnPlots=true") != std::string::npos) {DrawTextOnPlots=true; cout<<"DrawTextOnPlots"<<endl;}
		    if(std::string(argv[i]).find("useLeftSB=true") != std::string::npos) {useLeftSB=true; cout<<"useLeftSB"<<endl;}
		    if(std::string(argv[i]).find("useRightSB=true") != std::string::npos) {useRightSB=true; cout<<"useRightSB"<<endl;}
		    if(std::string(argv[i]).find("alteredToy=true") != std::string::npos) {alteredToy=true; cout<<"alteredToy"<<endl;}
		    if(std::string(argv[i]).find("PlotStep1=true") != std::string::npos) {PlotStep1=true; cout<<"PlotStep1"<<endl;}
		    if(std::string(argv[i]).find("PlotStep2=true") != std::string::npos) {PlotStep2=true; cout<<"PlotStep2"<<endl;}
		    if(std::string(argv[i]).find("finalPlots=true") != std::string::npos) {finalPlots=true; cout<<"finalPlots"<<endl;}
		    if(std::string(argv[i]).find("FixStuffForStep2=true") != std::string::npos) {FixStuffForStep2=true; cout<<"FixStuffForStep2"<<endl;}
		    if(std::string(argv[i]).find("MCsample=true") != std::string::npos) {MCsample=true; cout<<"MCsample"<<endl;}
		    if(std::string(argv[i]).find("useAnalyticalBKG=true") != std::string::npos) {useAnalyticalBKG=true; cout<<"useAnalyticalBKG"<<endl;}
		    if(std::string(argv[i]).find("BkgMixer=true") != std::string::npos) {BkgMixerBool=true; cout<<"BkgMixerBool"<<endl;}
		    if(std::string(argv[i]).find("BkgToy=true") != std::string::npos) {BkgToy=true; cout<<"BkgToy"<<endl;}
		    if(std::string(argv[i]).find("useExistingMixFile=true") != std::string::npos) {useExistingMixFile=true; cout<<"useExistingMixFile"<<endl;}
		    if(std::string(argv[i]).find("OneSPlotModel=false") != std::string::npos) {OneSPlotModel=false; cout<<"OneSPlotModel"<<endl;}
		    if(std::string(argv[i]).find("BkgToy1SBkgMixer2S=true") != std::string::npos) {BkgToy1SBkgMixer2S=true; cout<<"BkgToy1SBkgMixer2S"<<endl;}
		    if(std::string(argv[i]).find("useNewExistingMixFile=true") != std::string::npos) {useNewExistingMixFile=true; cout<<"useNewExistingMixFile"<<endl;}

	  }



	  bool Optimize=false;
	  if(gammaptDCut || dzCut || dzSigCut || Pi0Cut|| nYsigCut||vtxProbCut||ctSigCut||ctCut||gammaptCut||RConvUpperCut||RConvLowerCut||pTCut||vtxChi2ProbGammaCut||vtxChi2ProbGammaLogCut) {cout<<"Optimize!"<<endl; Optimize=true;}

	  char dirstruct[200];
	  sprintf(dirstruct,"Figures/%s",FitID);
	  gSystem->mkdir(dirstruct);
	  bool PlotBothSteps=false;
	  if(PlotStep1&&PlotStep2) PlotBothSteps=true;

// Plotting settings

	  int MarkerStyle_nS[4]={0,20,25,24}; // for each [nS] state
	  int MarkerColor_nS[4] = {0,600,418,632};// for each [nS] state
	  int LineColor_nS[4] = {0,600,418,632};// for each [nS] state
	  double MarkerSize_nS[4]={0,1.,1.,1.};// for each [nS] state

	  int LineColor_nJ[4] = {616,416,632,1};// for each [nj] state 0,1,2 and 3 for background
	  int LineColor_nP[4] = {0,410,808,869};// for each [nP] state 1,2, 3 and 0 for 2P->2S, reflected in 1S

	  double restrictMaximum_nPlots[8]={400,110,60,400,320,88,48,320};//0=1S, 1=2S, 2=3S, 3=1S2S3S, 4=1SNE, 5=2SNE, 6=3SNE, 7=1S2S3SNE

	  //Non-Equ-Plots:
	  double ScalePlot1S_NE=10.;
	  double ScalePlot2S_NE=10.;
	  double ScalePlot3S_NE=10.;

	  //Normal Plots:
	  double binWidth=12.5;

	  double PlotMinimum=0.1;

//load RooFit library
	using namespace RooFit;
    gSystem->Load("libRooFit");
    gROOT->SetBatch(1);
    gStyle->SetFrameBorderMode(0);

//open file
    char filename_[500];
    sprintf(filename_,"rooDS_%s.root",fileName);
    TFile data_file(filename_);

//grab roodataset
//    RooDataSet* data_=new RooDataSet();
//    data_=(RooDataSet*)data_file.Get("d");

//declare RooDataSet data variables
    char invmName[200];

    double Ymass1S;
    double Ymass2S;
    double Ymass3S;
    Ymass1S=9.4603;
    Ymass2S=10.02326;
    Ymass3S=10.3552;
	double M_PDG_chib1_1P=9.89278;
	double M_PDG_chib2_1P=9.91221;
	double M_PDG_chib1_2P=10.25546;
	double M_PDG_chib2_2P=10.26865;

    double invm_min1S;
    double invm_max1S;
    double invm_min2S;
    double invm_max2S;
    double invm_min3S;
    double invm_max3S;

    invm_min1S=9.5; invm_max1S=11.;
    invm_min2S=10; invm_max2S=11.;
    invm_min3S=10.3; invm_max3S=11.;

    if(MCsample) invm_min1S=9.7;
    if(MCsample) invm_max1S=10.35;


    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(1S)^{PDG}} [GeV]");
    RooRealVar invm1S("invm1S",invmName,invm_min1S,invm_max1S);
    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(2S)^{PDG}} [GeV]");
    RooRealVar invm2S("invm2S",invmName,invm_min2S,invm_max2S);
    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(3S)^{PDG}} [GeV]");
    RooRealVar invm3S("invm3S",invmName,invm_min3S,invm_max3S);

    invm1S.setRange("TotalRange1S",invm_min1S,invm_max1S);
    invm2S.setRange("TotalRange2S",invm_min2S,invm_max2S);
    invm3S.setRange("TotalRange3S",invm_min3S,invm_max3S);

    RooRealVar jpsipt   = RooRealVar("jpsipt", "jpsipt",0,200);
    RooRealVar jpsimass   = RooRealVar("jpsimass", "jpsimass",8,12);
    RooRealVar gammapt   = RooRealVar("gammapt", "gammapt",0,100);
    RooRealVar deltaRChiJpsi  = RooRealVar("deltaRChiJpsi","deltaRChiJpsi",0,10);
    RooRealVar deltaRJpsig    = RooRealVar("deltaRJpsig","deltaRJpsig",0,10);
    RooRealVar jpsieta  = RooRealVar("jpsieta","jpsieta",-5,5);
    RooRealVar ctpv     = RooRealVar("ctpv",   "ctpv",-5,5);
    RooRealVar ctpverr  = RooRealVar("ctpverr","ctpverr",0,10);
    RooRealVar ctbs     = RooRealVar("ctbs",   "ctbs",-5,5);
    RooRealVar ctbserr  = RooRealVar("ctbserr","ctbserr",0,10);
    RooRealVar ctpvsig     = RooRealVar("ctpvsig",   "ctpvsig",0,1000);
    RooRealVar ctbssig     = RooRealVar("ctbssig",   "ctbssig",0,1000);
    RooRealVar weight   = RooRealVar("weight", "weight",0,1);
    RooRealVar Rconv    = RooRealVar("Rconv","Rconv",0,100);
    RooRealVar jpsiVprob    = RooRealVar("jpsiVprob","jpsiVprob",0,1.5);
    RooRealVar Y1Smass_nSigma    = RooRealVar("Y1Smass_nSigma","Y1Smass_nSigma",-100,100);
    RooRealVar Y2Smass_nSigma    = RooRealVar("Y2Smass_nSigma","Y2Smass_nSigma",-100,100);
    RooRealVar Y3Smass_nSigma    = RooRealVar("Y3Smass_nSigma","Y3Smass_nSigma",-100,100);
    RooRealVar vertexChi2ProbGamma    = RooRealVar("vertexChi2ProbGamma","vertexChi2ProbGamma",0,100);

    RooRealVar jpsipx   = RooRealVar("jpsipx", "jpsipx",-100,100);
    RooRealVar jpsipy   = RooRealVar("jpsipy", "jpsipy",-100,100);
    RooRealVar jpsipz   = RooRealVar("jpsipz", "jpsipz",-100,100);
    RooRealVar gammapx   = RooRealVar("gammapx", "gammapx",-100,100);
    RooRealVar gammapy   = RooRealVar("gammapy", "gammapy",-100,100);
    RooRealVar gammapz   = RooRealVar("gammapz", "gammapz",-100,100);
    RooRealVar gammaeta   = RooRealVar("gammaeta", "gammaeta",-100,100);
    RooRealVar Q   = RooRealVar("Q", "Q",0,100);
    RooRealVar Pi0Mass   = RooRealVar("Pi0Mass", "Pi0Mass",0,5);

    RooRealVar vtxNsigmadz   = RooRealVar("vtxNsigmadz", "vtxNsigmadz",-15,15);
    RooRealVar vtxdz   = RooRealVar("vtxdz", "vtxdz",-15,15);
    RooRealVar vtxerrdz   = RooRealVar("vtxerrdz", "vtxerrdz",0,15);

    RooRealVar UpsIndex   	= RooRealVar("UpsIndex", "UpsIndex",0,10000000);
	RooRealVar GammaIndex   	= RooRealVar("GammaIndex", "GammaIndex",0,500);

	RooRealVar RunNb     = RooRealVar("RunNb", "RunNb",0,1e9);
	RooRealVar EventNb     = RooRealVar("EventNb", "EventNb",0,1e15);

    RooArgSet* FullArgSet=new RooArgSet();
    FullArgSet->add(invm1S);
    FullArgSet->add(invm2S);
    FullArgSet->add(invm3S);
    FullArgSet->add(jpsipt);
    FullArgSet->add(jpsimass);
    FullArgSet->add(gammapt);
    FullArgSet->add(deltaRChiJpsi);
    FullArgSet->add(deltaRJpsig);
    FullArgSet->add(jpsieta);
    FullArgSet->add(ctpv);
    FullArgSet->add(ctpverr);
    FullArgSet->add(ctbs);
    FullArgSet->add(ctbserr);
    FullArgSet->add(ctpvsig);
    FullArgSet->add(ctbssig);
    FullArgSet->add(weight);
    FullArgSet->add(Rconv);
    FullArgSet->add(jpsiVprob);
    FullArgSet->add(Y1Smass_nSigma);
    FullArgSet->add(Y2Smass_nSigma);
    FullArgSet->add(Y3Smass_nSigma);
    FullArgSet->add(vertexChi2ProbGamma);
    FullArgSet->add(jpsipx );
    FullArgSet->add(jpsipy );
    FullArgSet->add(jpsipz );
    FullArgSet->add(gammapx);
    FullArgSet->add(gammapy);
    FullArgSet->add(gammapz);
    FullArgSet->add(gammaeta);
    FullArgSet->add(Q);
    FullArgSet->add(Pi0Mass);
    FullArgSet->add(vtxNsigmadz);
    FullArgSet->add(vtxdz);
    FullArgSet->add(vtxerrdz);
    FullArgSet->add(UpsIndex);
    FullArgSet->add(GammaIndex);
    FullArgSet->add(RunNb);
    FullArgSet->add(EventNb);

    int n_Cut=1;
    double n_Cut_0;
    double deltaCut;
    int n_Cut_Default=5;

// nYsig
    if(nYsigCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=1.5;
        deltaCut=0.02;
        if(NarrowOptimization){
            n_Cut_0= 1.5;
            deltaCut=0.5;
        }
    }
// vtxProb
    if(vtxProbCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.001;
        deltaCut=0.0005;
    }
// ctSig
    if(ctSigCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.5;
        deltaCut=0.05;
        if(NarrowOptimization){
            n_Cut_0= 2.;
            deltaCut=0.5;
        }
    }
// ct
    if(ctCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.005;
        deltaCut=0.0005;
    }
// gammapt
    if(gammaptCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.3;
        deltaCut=0.017;
    }
// gammaptD
    if(gammaptDCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.3;
        deltaCut=0.017;
        if(NarrowOptimization){
            n_Cut_0= -0.3;
            deltaCut=0.15;
        }
    }

// RConvUpperCut
    if(RConvUpperCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=5;
        deltaCut=0.6;
    }
// RConvLowerCut
    if(RConvLowerCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=1.5;
        deltaCut=0.035;
        if(NarrowOptimization){
            n_Cut_0= 1.5;
            deltaCut=0.5;
        }
    }
// pT
    if(pTCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0;
        deltaCut=0.15;
    }
// vtxChi2ProbGamma
    if(vtxChi2ProbGammaCut){
        n_Cut=n_Cut_Default;
        n_Cut_0= 0;
        deltaCut=0.000005;
    }
// vtxChi2ProbGammaLog
        if(vtxChi2ProbGammaLogCut){
            n_Cut=n_Cut_Default;
            n_Cut_0= -50;
            deltaCut=0.5;
        }

// Pi0Cut
        if(Pi0Cut){
            n_Cut=n_Cut_Default;
            n_Cut_0= 0;
            deltaCut=0.00018;
            if(NarrowOptimization){
                n_Cut_0= 0.010;
                deltaCut=0.0025;
            }
        }

// dzSigCut
        if(dzSigCut){
            n_Cut=n_Cut_Default;
            n_Cut_0= 1.;
            deltaCut=0.2;
            if(NarrowOptimization){
                n_Cut_0= 2.;
                deltaCut=0.5;
            }
        }

// dzCut
        if(dzCut){
            n_Cut=n_Cut_Default;
            n_Cut_0= 0.05;
            deltaCut=0.025;
            if(NarrowOptimization){
                n_Cut_0= 0.25;
                deltaCut=0.25;
            }
        }

    TH1F  *hCut_1Psig = new TH1F("hCut_1Psig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1PSB = new TH1F("hCut_1PSB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1PS = new TH1F("hCut_1PS","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1PB = new TH1F("hCut_1PB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Psig = new TH1F("hCut_2Psig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Psig_2S = new TH1F("hCut_2Psig_2S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PSB = new TH1F("hCut_2PSB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PS = new TH1F("hCut_2PS","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PB = new TH1F("hCut_2PB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_sig = new TH1F("hCut_sig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_sigLL = new TH1F("hCut_sigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_covQual1 = new TH1F("hCut_covQual1","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_1Pmean = new TH1F("hCut_1Pmean","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1Pwidth = new TH1F("hCut_1Pwidth","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1Pnevt = new TH1F("hCut_1Pnevt","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pmean = new TH1F("hCut_2Pmean","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pwidth = new TH1F("hCut_2Pwidth","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3Pwidth = new TH1F("hCut_3Pwidth","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3Pwidth_1S = new TH1F("hCut_3Pwidth_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pwidth_2S = new TH1F("hCut_2Pwidth_2S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pnevt = new TH1F("hCut_2Pnevt","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pnevt_2S = new TH1F("hCut_2Pnevt_2S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_CBa = new TH1F("hCut_CBa","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_CBn = new TH1F("hCut_CBn","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_BGnevt1S = new TH1F("hCut_BGnevt1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_BGnevt2S = new TH1F("hCut_BGnevt2S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_Ntotal1S = new TH1F("hCut_Ntotal1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_Ntotal2S = new TH1F("hCut_Ntotal2S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_2PSover1PS = new TH1F("hCut_2PSover1PS","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_3Psig = new TH1F("hCut_3Psig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PSB = new TH1F("hCut_3PSB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PS = new TH1F("hCut_3PS","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PB = new TH1F("hCut_3PB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3Pmass = new TH1F("hCut_3Pmass","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PMerr = new TH1F("hCut_3PMerr","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_3Psig_1S = new TH1F("hCut_3Psig_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PSB_1S = new TH1F("hCut_3PSB_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PS_1S = new TH1F("hCut_3PS_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PB_1S = new TH1F("hCut_3PB_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_3PsigALL = new TH1F("hCut_3PsigALL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_PhotonEnergyScale = new TH1F("hCut_PhotonEnergyScale","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_PhotonEnergyScale2P = new TH1F("hCut_PhotonEnergyScale2P","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_PhotonSigmaScale = new TH1F("hCut_PhotonSigmaScale","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_PhotonSigmaScale2P = new TH1F("hCut_PhotonSigmaScale2P","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_1PsigLL = new TH1F("hCut_1PsigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PsigLL = new TH1F("hCut_2PsigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_sigLLhalf = new TH1F("hCut_sigLLhalf","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PsigLLhalf = new TH1F("hCut_3PsigLLhalf","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);



    for(int inCut=0; inCut<n_Cut; inCut++){

///////////////////////////////////////////////
//////////// CUTS /////////////////////////////
///////////////////////////////////////////////
        RooDataSet* data_=new RooDataSet();
        data_=(RooDataSet*)data_file.Get("d");

    	cout<<"Cut number "<<inCut<<endl;
        data_->Print();

// Y cuts

    double cut_rap=1000.;
    double nYsig=2.5;
    double cut_vtxProb=0.01;
    double cut_ctSig=3.;
    double cut_ct=1000.;//0.01;
    double cut_Ypt=0.;

// Gamma cuts

    double cut_gammapt = 0.;//1.1;
    double cut_gammapt1S = cut_gammapt;
    double cut_gammapt_max = 10000;//1.1;
    double gammaptK = 1.111;//1.388//1.111//1.666
    double gammaptD = 0.0111;//0.088//0.0111//-0.4333
    double cut_gammaeta=1000.;

//    AN Apr6
//    gammapt>1.388*Q+0.088 && Rconv>1.5 && ctpv<0.01 && abs(Y1Smass_nSigma)<2.5 && abs(Pi0Mass-0.135)>0.025 && vertexChi2ProbGamma>5e-4 && jpsiVprob>0.01

    double cut_RconvMin = 1.5;
    double cut_RconvMax = 200;
    double cut_vtxChi2ProbGamma = 0.;//5e-4;//5e-4;//000005;
    double cut_vtxChi2ProbGammaLog = -1000.;//1e-5;//000005;

    double cut_Pi0 = 0.015;//0.03
    double Pi0_mean = 0.134;//0.134

// Y - Gamma cuts

    double MINcosalphaCut=-1000;
    double cut_dzSig=15.;
    double cut_dz=0.5;//1.

//1S extreme cuts
/*    nYsig=3;
    Pi0_mean=0.135;
    cut_Pi0=0.1;
    cut_Ypt=10;
    cut_gammaeta=1.1;
    cut_rap=1000;
    gammaptD=-1000;
    gammaptK=0;
    cut_dz=1000;
    cut_ctSig=1000;
    cut_vtxProb=0.01;
*/
/*//2S extreme cuts
     nYsig=2;
    Pi0_mean=0.165;
    cut_Pi0=0.03;
    cut_Ypt=6;
    cut_gammaeta=1.1;
    cut_rap=1000;
    gammaptD=-1000;
    gammaptK=0;
    cut_dz=1000;
    cut_ctSig=1000;
    cut_vtxProb=0.01;
*/
    //Ernest Cuts
    double GammaEtaBorder=0.9;
    double cut_Pi0_Midrap=0.;//0.0154413;
    double cut_Pi0_Forward=0.;//0.0204159;
    double cut_gammapt_Midrap=0.;//0.392851;
    double cut_gammapt_Forward=0.;//0.320923;

		cut_rap=1.25;
        gammaptD=-1000;
        gammaptK=0;

        cut_dzSig=5.;
        cut_dz=1.35;
        cut_ctSig=1000.;
        cut_vtxProb=0.01;

        nYsig=2.5;
        cut_Ypt=9.5;
        cut_gammaeta=1.4;
        Pi0_mean=0.137;
        cut_Pi0=0.018;
        cut_gammapt=0.;//0.4;




///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////

    if(nYsigCut)    	nYsig=n_Cut_0+inCut*deltaCut;
    if(vtxProbCut)    cut_vtxProb=n_Cut_0+inCut*deltaCut;
    if(ctSigCut)    cut_ctSig=n_Cut_0+inCut*deltaCut;
    if(gammaptCut)    cut_gammapt=n_Cut_0+inCut*deltaCut;
    if(ctCut)    	cut_ct=n_Cut_0+inCut*deltaCut;
    if(RConvUpperCut)	cut_RconvMax=n_Cut_0+inCut*deltaCut;
    if(RConvLowerCut)	cut_RconvMin=n_Cut_0+inCut*deltaCut;
    if(pTCut)    	cut_Ypt=n_Cut_0+inCut*deltaCut;
    if(vtxChi2ProbGammaCut)    	cut_vtxChi2ProbGamma=n_Cut_0+inCut*deltaCut;
    if(vtxChi2ProbGammaLogCut)    	cut_vtxChi2ProbGammaLog=n_Cut_0+inCut*deltaCut;
    if(Pi0Cut)    	cut_Pi0=n_Cut_0+inCut*deltaCut;
    if(dzSigCut)    	cut_dzSig=n_Cut_0+inCut*deltaCut;
    if(dzCut)    	cut_dz=n_Cut_0+inCut*deltaCut;
    if(gammaptDCut)    gammaptD=n_Cut_0+inCut*deltaCut;

//Cut RooDataSet


    char DScutChar[2000];

	 RooArgSet CutVars("CutVars");
	 CutVars.add(invm1S);
	 CutVars.add(jpsimass);
	 CutVars.add(invm2S);
	 CutVars.add(invm3S);
	 CutVars.add(Y1Smass_nSigma);
	 CutVars.add(Y2Smass_nSigma);
	 CutVars.add(Y3Smass_nSigma);
	 CutVars.add(gammapx);
	 CutVars.add(gammapy);
	 CutVars.add(gammapz);
	 CutVars.add(gammapt);
	 CutVars.add(jpsipx);
	 CutVars.add(jpsipy);
	 CutVars.add(jpsipz);
	 CutVars.add(jpsieta);
	 CutVars.add(gammaeta);
	 CutVars.add(Rconv);

	    TTree* treeBeforeAllCuts=new TTree();
	    if(SaveAll){
	    	treeBeforeAllCuts=(TTree*)data_->tree();
	    }

//    sprintf(DScutChar,"jpsipt>%f && gammapt>%f && gammapt<%f && Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && (TMath::Sign(-1.,Q-1.65)+1)/2.*0.8+(TMath::Sign(-1.,1.65-Q)+1)/2.*(gammapt>%f*Q+%f) && (gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>%f",cut_Ypt,cut_gammapt,cut_gammapt_max,cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,gammaptK,gammaptD,MINcosalphaCut);

	    char FinalCutChar[2000];

	 cout<<"dataAfterBasicCuts"<<endl;
    sprintf(DScutChar,"jpsipt>%f && gammapt>%f && gammapt<%f && Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && gammapt>%f*Q+%f && (gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>%f && TMath::Abs(vtxNsigmadz) < %f &&TMath::Abs(Pi0Mass-%f)>%f &&TMath::Abs(vtxdz)<%f &&TMath::Abs(ctpvsig)<%f && TMath::Abs(gammaeta)<%f && (TMath::Abs(gammaeta)<%f && TMath::Abs(Pi0Mass-%f)>%f || TMath::Abs(gammaeta)>%f && TMath::Abs(Pi0Mass-%f)>%f) && (TMath::Abs(gammaeta)<%f && gammapt>%f || TMath::Abs(gammaeta)>%f && gammapt>%f )",cut_Ypt,cut_gammapt,cut_gammapt_max,cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,gammaptK,gammaptD,MINcosalphaCut,cut_dzSig, Pi0_mean, cut_Pi0, cut_dz, cut_ctSig, cut_gammaeta, GammaEtaBorder, Pi0_mean, cut_Pi0_Midrap, GammaEtaBorder, Pi0_mean, cut_Pi0_Forward, GammaEtaBorder, cut_gammapt_Midrap, GammaEtaBorder, cut_gammapt_Forward);
//    sprintf(DScutChar,"jpsipt>%f && gammapt>%f && gammapt<%f && Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && gammapt>%f*Q+%f && (gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>%f",cut_Ypt,cut_gammapt,cut_gammapt_max,cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,gammaptK,gammaptD,MINcosalphaCut);
    RooDataSet* dataAfterBasicCuts=(RooDataSet*)data_->reduce(Cut(DScutChar));
    cout<<DScutChar<<endl;cout<<endl;
    dataAfterBasicCuts->Print();

    sprintf(cutName_,"vtxY%d_ct%d_ctSig%d_nYsig%d_gampt%d_gammaptD%d_RMin%d_RMax%d_Ypt%d_vtxgam%d_vtxgamLog%d_Pi0%d_dzSig%d_dz%d",int(1000000*cut_vtxProb),int(1000000*cut_ct),int(1000000*cut_ctSig),int(1000000*nYsig),int(1000000*cut_gammapt),int(1000000*gammaptD),int(1000000*cut_RconvMin),int(1000000*cut_RconvMax),int(1000000*cut_Ypt),int(100000000*cut_vtxChi2ProbGamma),-int(10000*cut_vtxChi2ProbGammaLog),int(100000*cut_Pi0),int(100000*cut_dzSig),int(100000*cut_dz));
    sprintf(FinalCutChar,"%s",DScutChar);

    RooDataSet* dataAfterBasicCutsExceptGammaPtCut;
    if(BkgMixerBool&&!useExistingMixFile){
   	 cout<<"dataAfterBasicCutsExceptGammaPtCut"<<endl;
    	char DScutCharNoGammaPtCut[2000];
//        sprintf(DScutCharNoGammaPtCut,"Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && (gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>%f && TMath::Abs(vtxNsigmadz) < %f &&TMath::Abs(Pi0Mass-%f)>%f &&TMath::Abs(vtxdz)<%f &&TMath::Abs(ctpvsig)<%f",cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,MINcosalphaCut,cut_dzSig, Pi0_mean, cut_Pi0, cut_dz, cut_ctSig);
        sprintf(DScutCharNoGammaPtCut,"jpsipt>%f && gammapt>%f && gammapt<%f && Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && gammapt>%f*Q+%f && (gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>%f && TMath::Abs(vtxNsigmadz) < %f &&TMath::Abs(Pi0Mass-%f)>%f &&TMath::Abs(vtxdz)<%f &&TMath::Abs(ctpvsig)<%f && TMath::Abs(gammaeta)<%f && (TMath::Abs(gammaeta)<%f && TMath::Abs(Pi0Mass-%f)>%f || TMath::Abs(gammaeta)>%f && TMath::Abs(Pi0Mass-%f)>%f) && (TMath::Abs(gammaeta)<%f && gammapt>%f || TMath::Abs(gammaeta)>%f && gammapt>%f )",cut_Ypt,cut_gammapt,cut_gammapt_max,cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,gammaptK,gammaptD,MINcosalphaCut,cut_dzSig, Pi0_mean, cut_Pi0, cut_dz, cut_ctSig, cut_gammaeta, GammaEtaBorder, Pi0_mean, cut_Pi0_Midrap, GammaEtaBorder, Pi0_mean, cut_Pi0_Forward, GammaEtaBorder, cut_gammapt_Midrap, GammaEtaBorder, cut_gammapt_Forward);
        dataAfterBasicCutsExceptGammaPtCut=(RooDataSet*)data_->reduce(Cut(DScutCharNoGammaPtCut));
        cout<<"without gammapt cut:"<<endl;
        dataAfterBasicCutsExceptGammaPtCut->Print();
    }

/////// Producing SB dataset

        double LowerSBmin;
        double nSigSBlow;
        double nSigSBhigh;
        double UpperSBmax;


        double LowerSBmin2S;
        double nSigSBlow2S;
        double nSigSBhigh2S;
        double UpperSBmax2S;

        if(useLeftSB&&useRightSB){
        	LowerSBmin=8.7;
        	nSigSBlow=-3.5;
        	nSigSBhigh=3.5;
        	UpperSBmax=10.85;

        	LowerSBmin2S=8.7;
        	nSigSBlow2S=-3.5;
        	nSigSBhigh2S=3.5;
        	UpperSBmax2S=10.85;
        }

        if(useLeftSB&&!useRightSB){
        	LowerSBmin=8.7;
        	nSigSBlow=-3.5;
        	nSigSBhigh=1000;
        	UpperSBmax=0;

        	LowerSBmin2S=8.7;
        	nSigSBlow2S=-3.5;
        	nSigSBhigh2S=1000;
        	UpperSBmax2S=0;
        }

        if(!useLeftSB&&useRightSB){
        	LowerSBmin=1000;
        	nSigSBlow=-1000;
        	nSigSBhigh=3.5;
        	UpperSBmax=10.85;

        	LowerSBmin2S=1000;
        	nSigSBlow2S=-1000;
        	nSigSBhigh2S=3.5;
        	UpperSBmax2S=10.85;
        }

      	 cout<<"dataSB1S"<<endl;
    sprintf(DScutChar,"jpsimass > %f && Y1Smass_nSigma < %f || jpsimass < %f && Y2Smass_nSigma > %f",LowerSBmin,nSigSBlow,UpperSBmax,nSigSBhigh);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB1S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    dataSB1S->Print();

 	 cout<<"dataSB2S"<<endl;
    sprintf(DScutChar,"jpsimass > %f && Y1Smass_nSigma < %f || jpsimass < %f && Y2Smass_nSigma > %f",LowerSBmin2S,nSigSBlow2S,UpperSBmax2S,nSigSBhigh2S);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB2S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    dataSB2S->Print();


/////// Producing Signal dataset
	 cout<<"data1S"<<endl;
    sprintf(DScutChar,"Y1Smass_nSigma > %f && Y1Smass_nSigma < %f && Y2Smass_nSigma < %f",-nYsig,nYsig,10000.);
    cout<<DScutChar<<endl;
    RooDataSet* data1S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    data1S->Print();

    RooDataSet* data1SExceptGammaPtCut;
    if(BkgMixerBool&&!useExistingMixFile){
   	 cout<<"data1SExceptGammaPtCut"<<endl;
        data1SExceptGammaPtCut=(RooDataSet*)dataAfterBasicCutsExceptGammaPtCut->reduce(Cut(DScutChar));
        cout<<"without gammapt cut:"<<endl;
        data1SExceptGammaPtCut->Print();
    }

  	 cout<<"data2S"<<endl;
    sprintf(DScutChar,"Y2Smass_nSigma > %f && Y2Smass_nSigma < %f ",-nYsig*Ymass2S/Ymass1S,nYsig*Ymass2S/Ymass1S);
    cout<<DScutChar<<endl;
    RooDataSet* data2S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    data2S->Print();

    RooDataSet* data2SExceptGammaPtCut;
    if(BkgMixerBool&&!useExistingMixFile){
     	 cout<<"data2SExceptGammaPtCut"<<endl;
        data2SExceptGammaPtCut=(RooDataSet*)dataAfterBasicCutsExceptGammaPtCut->reduce(Cut(DScutChar));
        cout<<"without gammapt cut:"<<endl;
        data2SExceptGammaPtCut->Print();
    }

	 cout<<"data3S"<<endl;
    sprintf(DScutChar,"Y3Smass_nSigma > %f && Y3Smass_nSigma < %f && Y2Smass_nSigma > %f",-nYsig*Ymass3S/Ymass2S,nYsig*Ymass3S/Ymass2S,-10000.);
    cout<<DScutChar<<endl;
    RooDataSet* data3S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    data3S->Print();

    RooDataSet* data3SExceptGammaPtCut;
    if(BkgMixerBool&&!useExistingMixFile){
     	 cout<<"data3SExceptGammaPtCut"<<endl;
     	data3SExceptGammaPtCut=(RooDataSet*)dataAfterBasicCutsExceptGammaPtCut->reduce(Cut(DScutChar));
        cout<<"without gammapt cut:"<<endl;
        data3SExceptGammaPtCut->Print();
    }

    if(SetSignalToZero&&useSBforBKGmodel){
    	data1S=dataSB1S;
    	data2S=dataSB2S;
    	data3S=dataSB2S;
    }
    TTree* tree1S=new TTree();
    TTree* tree2S=new TTree();
    TTree* tree3S=new TTree();
    TTree* treeAllCuts=new TTree();

    TTree* tree1SExceptGammaPtCut=new TTree();
    TTree* tree2SExceptGammaPtCut=new TTree();
    TTree* tree3SExceptGammaPtCut=new TTree();
    if(BkgMixerBool&&!useExistingMixFile){
    	tree1SExceptGammaPtCut=(TTree*)data1SExceptGammaPtCut->tree();
    	tree2SExceptGammaPtCut=(TTree*)data2SExceptGammaPtCut->tree();
    	tree3SExceptGammaPtCut=(TTree*)data3SExceptGammaPtCut->tree();
    }

    tree1S=(TTree*)data1S->tree();
    tree2S=(TTree*)data2S->tree();
    tree3S=(TTree*)data3S->tree();
    treeAllCuts=(TTree*)dataAfterBasicCuts->tree();

    if(useSBforBKGmodel){
    	tree1S=(TTree*)dataSB1S->tree();
    	tree2S=(TTree*)dataSB2S->tree();
    	tree3S=(TTree*)dataSB2S->tree();
    }

    if(SetSignalToZero){
    	data1S=dataSB1S;
    	data2S=dataSB2S;
    	data3S=dataSB2S;
    }

    if(MCsample){
        tree2S=(TTree*)data1S->tree();
        tree3S=(TTree*)data1S->tree();
    	data2S=data1S;
    	data3S=data1S;
    }

    char treeName[200];
    sprintf(treeName,"Figures/%s/tree_Y1S.root",FitID);
    if(SaveAll) tree1S->SaveAs(treeName);
    sprintf(treeName,"Figures/%s/tree_Y2S.root",FitID);
     if(SaveAll) tree2S->SaveAs(treeName);
     sprintf(treeName,"Figures/%s/tree_Y3S.root",FitID);
      if(SaveAll) tree3S->SaveAs(treeName);
    sprintf(treeName,"Figures/%s/treeAllCuts.root",FitID);
    if(SaveAll) treeAllCuts->SaveAs(treeName);
    sprintf(treeName,"Figures/%s/tree1SExceptGammaPtCut.root",FitID);

    if(SaveAll){

/*    TH1F  *MassHisto1S = new TH1F("MassHisto1S","",100,-20,20);
    treeAllCuts->Draw("Y1Smass_nSigma>>MassHisto1S");
    MassCanvas->SetFillColor(kWhite);
    MassHisto1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); MassHisto1S->SetStats(0);MassHisto1S->GetXaxis()->SetTitle("n_{#sigma} Y(1S)");	MassHisto1S->Draw();
    MassCanvas->Modified();
    sprintf(saveMass,"Figures/%s/nSigmaY1SMass.pdf",FitID);
    MassCanvas->SaveAs(saveMass);

    TH1F  *MassHisto2S = new TH1F("MassHisto2S","",100,-20,20);
    treeAllCuts->Draw("Y2Smass_nSigma>>MassHisto2S");
    MassCanvas->SetFillColor(kWhite);
    MassHisto2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); MassHisto2S->SetStats(0);MassHisto2S->GetXaxis()->SetTitle("n_{#sigma} Y(2S)");	MassHisto2S->Draw();
    MassCanvas->Modified();
    sprintf(saveMass,"Figures/%s/nSigmaY2SMass.pdf",FitID);
    MassCanvas->SaveAs(saveMass);

    TH1F  *MassHisto3S = new TH1F("MassHisto3S","",100,-20,20);
    treeAllCuts->Draw("Y3Smass_nSigma>>MassHisto3S");
    MassCanvas->SetFillColor(kWhite);
    MassHisto3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); MassHisto3S->SetStats(0);MassHisto3S->GetXaxis()->SetTitle("n_{#sigma} Y(3S)");	MassHisto3S->Draw();
    MassCanvas->Modified();
    sprintf(saveMass,"Figures/%s/nSigmaY3SMass.pdf",FitID);
    MassCanvas->SaveAs(saveMass);
*/

        char saveMass[200];
        char InvmCutChar[200];
        sprintf(InvmCutChar,"invm1S<%f",invm_max1S);
        double yOff=1.3;

    TH1F  *MassHisto = new TH1F("MassHisto","",112,8.6,11.4);//25MeV
    TH1F  *GammaPtHisto = new TH1F("GammaPtHisto","",100,0,15);//150MeV
    TH1F  *GammaEtaHisto = new TH1F("GammaEtaHisto","",100,-2.6,2.6);
    TH1F  *UpsPtHisto = new TH1F("UpsPtHisto","",100,0,50);//500MeV
    TH1F  *UpsRapHisto = new TH1F("UpsRapHisto","",100,-1.75,1.75);
    TH1F  *RConvHisto = new TH1F("RConvHisto","",150,0,75);//0.5cm

    TH1F  *MassHistoInvmCut = new TH1F("MassHistoInvmCut","",112,8.6,11.4);//25MeV
    TH1F  *GammaPtHistoInvmCut = new TH1F("GammaPtHistoInvmCut","",100,0,5);//50MeV
    TH1F  *GammaEtaHistoInvmCut = new TH1F("GammaEtaHistoInvmCut","",100,-2.6,2.6);
    TH1F  *UpsPtHistoInvmCut = new TH1F("UpsPtHistoInvmCut","",100,0,50);//500MeV
    TH1F  *UpsRapHistoInvmCut = new TH1F("UpsRapHistoInvmCut","",100,-1.75,1.75);
    TH1F  *RConvHistoInvmCut = new TH1F("RConvHistoInvmCut","",100,0,75);//0.5cm

    treeAllCuts->Draw("jpsimass>>MassHisto");
    treeAllCuts->Draw("sqrt(gammapx*gammapx+gammapy*gammapy)>>GammaPtHisto");
    treeAllCuts->Draw("0.5*TMath::Log((sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)+gammapz)/(sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)-gammapz))>>GammaEtaHisto");
    treeAllCuts->Draw("jpsipt>>UpsPtHisto");
    treeAllCuts->Draw("jpsieta>>UpsRapHisto");
    treeAllCuts->Draw("Rconv>>RConvHisto");

    treeAllCuts->Draw("jpsimass>>MassHistoInvmCut",InvmCutChar);
    treeAllCuts->Draw("sqrt(gammapx*gammapx+gammapy*gammapy)>>GammaPtHistoInvmCut",InvmCutChar);
    treeAllCuts->Draw("0.5*TMath::Log((sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)+gammapz)/(sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)-gammapz))>>GammaEtaHistoInvmCut",InvmCutChar);
    treeAllCuts->Draw("jpsipt>>UpsPtHistoInvmCut",InvmCutChar);
    treeAllCuts->Draw("jpsieta>>UpsRapHistoInvmCut",InvmCutChar);
    treeAllCuts->Draw("Rconv>>RConvHistoInvmCut",InvmCutChar);


    TH1F  *AN_MassHisto = new TH1F("AN_MassHisto","",112,8.6,11.4);//25MeV
    TH1F  *AN_UpsPtHisto = new TH1F("AN_UpsPtHisto","",100,0,50);//500MeV
    TH1F  *AN_UpsRapHisto = new TH1F("AN_UpsRapHisto","",100,-1.75,1.75);
    TH1F  *AN_DimuonVertexProbHisto = new TH1F("AN_DimuonVertexProbHisto","",100,0,100);

    TH1F  *AN_1S_GammaVertexProbHisto = new TH1F("AN_1S_GammaVertexProbHisto","",100,0,100);
    TH1F  *AN_1S_dzsigHisto = new TH1F("AN_1S_dzsigHisto","",100,-5,5);
    TH1F  *AN_1S_dzHisto = new TH1F("AN_1S_dzHisto","",100,-5,5);
    TH1F  *AN_1S_GammaEtaHisto = new TH1F("AN_1S_GammaEtaHisto","",100,-2.5,2.5);
    TH1F  *AN_1S_GammaPtHisto = new TH1F("AN_1S_GammaPtHisto","",67,0,10.05);//150MeV
    TH1F  *AN_1S_Pi0MassHistoBroad = new TH1F("AN_1S_Pi0MassHistoBroad","",100,0,1);//10MeV
    TH1F  *AN_1S_Pi0MassHistoNarrow = new TH1F("AN_1S_Pi0MassHistoNarrow","",100,0.035,0.235);//2MeV
    TH1F  *AN_1S_RConvHisto = new TH1F("AN_1S_RConvHisto","",150,0,75);//0.5cm

    TH1F  *AN_2S_GammaVertexProbHisto = new TH1F("AN_2S_GammaVertexProbHisto","",100,0,100);
    TH1F  *AN_2S_dzsigHisto = new TH1F("AN_2S_dzsigHisto","",100,-5,5);
    TH1F  *AN_2S_dzHisto = new TH1F("AN_2S_dzHisto","",100,-5,5);
    TH1F  *AN_2S_GammaEtaHisto = new TH1F("AN_2S_GammaEtaHisto","",100,-2.5,2.5);
    TH1F  *AN_2S_GammaPtHisto = new TH1F("AN_2S_GammaPtHisto","",67,0,10.05);//150MeV
    TH1F  *AN_2S_Pi0MassHistoBroad = new TH1F("AN_2S_Pi0MassHistoBroad","",100,0,1);//10MeV
    TH1F  *AN_2S_Pi0MassHistoNarrow = new TH1F("AN_2S_Pi0MassHistoNarrow","",100,0.035,0.235);//2MeV
    TH1F  *AN_2S_RConvHisto = new TH1F("AN_2S_RConvHisto","",150,0,75);//0.5cm

    TH1F  *AN_3S_GammaVertexProbHisto = new TH1F("AN_3S_GammaVertexProbHisto","",100,0,100);
    TH1F  *AN_3S_dzsigHisto = new TH1F("AN_3S_dzsigHisto","",100,-5,5);
    TH1F  *AN_3S_dzHisto = new TH1F("AN_3S_dzHisto","",100,-5,5);
    TH1F  *AN_3S_GammaEtaHisto = new TH1F("AN_3S_GammaEtaHisto","",100,-2.5,2.5);
    TH1F  *AN_3S_GammaPtHisto = new TH1F("AN_3S_GammaPtHisto","",67,0,10.05);//150MeV
    TH1F  *AN_3S_Pi0MassHistoBroad = new TH1F("AN_3S_Pi0MassHistoBroad","",100,0,1);//10MeV
    TH1F  *AN_3S_Pi0MassHistoNarrow = new TH1F("AN_3S_Pi0MassHistoNarrow","",100,0.035,0.235);//2MeV
    TH1F  *AN_3S_RConvHisto = new TH1F("AN_3S_RConvHisto","",150,0,75);//0.5cm


    char PlotCutChar[1000];
    char PlotCutChar1S[1000];
    char PlotCutChar2S[1000];
    char PlotCutChar3S[1000];

    sprintf(PlotCutChar,"");
    treeBeforeAllCuts->Draw("jpsimass>>AN_MassHisto",PlotCutChar);
    treeBeforeAllCuts->Draw("jpsipt>>AN_UpsPtHisto",PlotCutChar);
    treeBeforeAllCuts->Draw("jpsieta>>AN_UpsRapHisto",PlotCutChar);
    treeBeforeAllCuts->Draw("100*jpsiVprob>>AN_DimuonVertexProbHisto",PlotCutChar);

    sprintf(PlotCutChar1S,"TMath::Abs(Y1Smass_nSigma)<%f && jpsipt>%f && jpsiVprob>%f && TMath::Abs(jpsieta)<%f", nYsig, cut_Ypt, cut_vtxProb, cut_rap);
    treeBeforeAllCuts->Draw("100*vertexChi2ProbGamma>>AN_1S_GammaVertexProbHisto",PlotCutChar1S);
    treeBeforeAllCuts->Draw("vtxNsigmadz>>AN_1S_dzsigHisto",PlotCutChar1S);
    treeBeforeAllCuts->Draw("vtxdz>>AN_1S_dzHisto",PlotCutChar1S);
    treeBeforeAllCuts->Draw("Rconv>>AN_1S_RConvHisto",PlotCutChar1S);

    sprintf(PlotCutChar1S,"%s && TMath::Abs(vtxNsigmadz)<%f && TMath::Abs(vtxdz)<%f && invm1S<%f", PlotCutChar1S, cut_dzSig, cut_dz,invm_max1S);
    treeBeforeAllCuts->Draw("Pi0Mass>>AN_1S_Pi0MassHistoBroad",PlotCutChar1S);
    treeBeforeAllCuts->Draw("Pi0Mass>>AN_1S_Pi0MassHistoNarrow",PlotCutChar1S);

    sprintf(PlotCutChar1S,"%s && TMath::Abs(Pi0Mass-%f)>%f", PlotCutChar1S, Pi0_mean, cut_Pi0);
    treeBeforeAllCuts->Draw("gammapt>>AN_1S_GammaPtHisto",PlotCutChar1S);

    sprintf(PlotCutChar1S,"%s && gammapt>%f", PlotCutChar1S, cut_gammapt);
    treeBeforeAllCuts->Draw("gammaeta>>AN_1S_GammaEtaHisto",PlotCutChar1S);


    sprintf(PlotCutChar2S,"TMath::Abs(Y2Smass_nSigma)<%f && jpsipt>%f && jpsiVprob>%f && TMath::Abs(jpsieta)<%f", nYsig*Ymass2S/Ymass1S, cut_Ypt, cut_vtxProb, cut_rap);
    treeBeforeAllCuts->Draw("100*vertexChi2ProbGamma>>AN_2S_GammaVertexProbHisto",PlotCutChar2S);
    treeBeforeAllCuts->Draw("vtxNsigmadz>>AN_2S_dzsigHisto",PlotCutChar2S);
    treeBeforeAllCuts->Draw("vtxdz>>AN_2S_dzHisto",PlotCutChar2S);
    treeBeforeAllCuts->Draw("Rconv>>AN_2S_RConvHisto",PlotCutChar2S);

    sprintf(PlotCutChar2S,"%s && TMath::Abs(vtxNsigmadz)<%f && TMath::Abs(vtxdz)<%f && invm2S<%f", PlotCutChar2S, cut_dzSig, cut_dz,invm_max2S);
    treeBeforeAllCuts->Draw("Pi0Mass>>AN_2S_Pi0MassHistoBroad",PlotCutChar2S);
    treeBeforeAllCuts->Draw("Pi0Mass>>AN_2S_Pi0MassHistoNarrow",PlotCutChar2S);

    sprintf(PlotCutChar2S,"%s && TMath::Abs(Pi0Mass-%f)>%f", PlotCutChar2S, Pi0_mean, cut_Pi0);
    treeBeforeAllCuts->Draw("gammapt>>AN_2S_GammaPtHisto",PlotCutChar2S);

    sprintf(PlotCutChar2S,"%s && gammapt>%f", PlotCutChar2S, cut_gammapt);
    treeBeforeAllCuts->Draw("gammaeta>>AN_2S_GammaEtaHisto",PlotCutChar2S);



    sprintf(PlotCutChar3S,"TMath::Abs(Y3Smass_nSigma)<%f && jpsipt>%f && jpsiVprob>%f && TMath::Abs(jpsieta)<%f", nYsig*Ymass3S/Ymass2S, cut_Ypt, cut_vtxProb, cut_rap);
    treeBeforeAllCuts->Draw("100*vertexChi2ProbGamma>>AN_3S_GammaVertexProbHisto",PlotCutChar3S);
    treeBeforeAllCuts->Draw("vtxNsigmadz>>AN_3S_dzsigHisto",PlotCutChar3S);
    treeBeforeAllCuts->Draw("vtxdz>>AN_3S_dzHisto",PlotCutChar3S);
    treeBeforeAllCuts->Draw("Rconv>>AN_3S_RConvHisto",PlotCutChar3S);

    sprintf(PlotCutChar3S,"%s && TMath::Abs(vtxNsigmadz)<%f && TMath::Abs(vtxdz)<%f && invm3S<%f", PlotCutChar3S, cut_dzSig, cut_dz,invm_max3S);
    treeBeforeAllCuts->Draw("Pi0Mass>>AN_3S_Pi0MassHistoBroad",PlotCutChar3S);
    treeBeforeAllCuts->Draw("Pi0Mass>>AN_3S_Pi0MassHistoNarrow",PlotCutChar3S);

    sprintf(PlotCutChar3S,"%s && TMath::Abs(Pi0Mass-%f)>%f", PlotCutChar3S, Pi0_mean, cut_Pi0);
    treeBeforeAllCuts->Draw("gammapt>>AN_3S_GammaPtHisto",PlotCutChar3S);

    sprintf(PlotCutChar3S,"%s && gammapt>%f", PlotCutChar3S, cut_gammapt);
    treeBeforeAllCuts->Draw("gammaeta>>AN_3S_GammaEtaHisto",PlotCutChar3S);




/*    nYsig=2.5;
    cut_Ypt=9.5;
    cut_gammaeta=1.4;
    cut_rap=1000;
    gammaptD=-1000;
    gammaptK=0;
    cut_dzSig=5.;
    cut_dz=1.35;
    cut_ctSig=1000000.;
    cut_vtxProb=0.01;
    Pi0_mean=0.134;
    cut_Pi0=0.017;
    cut_gammapt=0.4;//0.4;
*/

    TCanvas* PlotCanvas = new TCanvas("PlotCanvas","PlotCanvas",1600, 1200);
    PlotCanvas->SetFillColor(kWhite);
    gPad->SetLeftMargin(0.125);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
//    gPad->SetBottomMargin(0.1);
    gPad->SetFillColor(kWhite);



	  TLine* cutLine;
	  cutLine = new TLine( 0,0,0,0 );
	  cutLine->SetLineWidth( 1 );
	  cutLine->SetLineStyle( 2 );
	  cutLine->SetLineColor( kGreen+2 );
	  cutLine->SetY1(0);
	  TLine* cutLine2;
	  cutLine2 = new TLine( 0,0,0,0 );
	  cutLine2->SetLineWidth( 1 );
	  cutLine2->SetLineStyle( 2 );
	  cutLine2->SetLineColor( kGreen+2 );
	  cutLine2->SetY1(0);

	  double MaxFact=1.2;
	  yOff=1.5;
	  double PlotMinDist=1.;

	  AN_MassHisto->SetMinimum(PlotMinDist); AN_MassHisto->SetMaximum(AN_MassHisto->GetMaximum()*MaxFact);
	  AN_UpsPtHisto->SetMinimum(PlotMinDist); AN_UpsPtHisto->SetMaximum(AN_UpsPtHisto->GetMaximum()*MaxFact);
	  AN_UpsRapHisto->SetMinimum(PlotMinDist); AN_UpsRapHisto->SetMaximum(AN_UpsRapHisto->GetMaximum()*MaxFact);
	  AN_DimuonVertexProbHisto->SetMinimum(PlotMinDist); AN_DimuonVertexProbHisto->SetMaximum(AN_DimuonVertexProbHisto->GetMaximum()*MaxFact);

	  AN_1S_GammaVertexProbHisto->SetMinimum(PlotMinDist*100); AN_1S_GammaVertexProbHisto->SetMaximum(AN_1S_GammaVertexProbHisto->GetMaximum()*MaxFact);
	  AN_1S_dzsigHisto->SetMinimum(PlotMinDist*10); AN_1S_dzsigHisto->SetMaximum(AN_1S_dzsigHisto->GetMaximum()*MaxFact);
	  AN_1S_dzHisto->SetMinimum(PlotMinDist); AN_1S_dzHisto->SetMaximum(AN_1S_dzHisto->GetMaximum()*MaxFact);
	  AN_1S_GammaEtaHisto->SetMinimum(PlotMinDist); AN_1S_GammaEtaHisto->SetMaximum(AN_1S_GammaEtaHisto->GetMaximum()*MaxFact);
	  AN_1S_GammaPtHisto->SetMinimum(PlotMinDist); AN_1S_GammaPtHisto->SetMaximum(AN_1S_GammaPtHisto->GetMaximum()*MaxFact);
	  AN_1S_Pi0MassHistoBroad->SetMinimum(PlotMinDist); AN_1S_Pi0MassHistoBroad->SetMaximum(AN_1S_Pi0MassHistoBroad->GetMaximum()*MaxFact);
	  AN_1S_Pi0MassHistoNarrow->SetMinimum(PlotMinDist); AN_1S_Pi0MassHistoNarrow->SetMaximum(AN_1S_Pi0MassHistoNarrow->GetMaximum()*MaxFact);
	  AN_1S_RConvHisto->SetMinimum(PlotMinDist); AN_1S_RConvHisto->SetMaximum(AN_1S_RConvHisto->GetMaximum()*MaxFact);

/*	  AN_2S_GammaVertexProbHisto->SetMinimum(PlotMinDist*100); AN_2S_GammaVertexProbHisto->SetMaximum(AN_2S_GammaVertexProbHisto->GetMaximum()*MaxFact);
	  AN_2S_dzsigHisto->SetMinimum(PlotMinDist*100); AN_2S_dzsigHisto->SetMaximum(AN_2S_dzsigHisto->GetMaximum()*MaxFact);
	  AN_2S_dzHisto->SetMinimum(PlotMinDist); AN_2S_dzHisto->SetMaximum(AN_2S_dzHisto->GetMaximum()*MaxFact);
	  AN_2S_GammaEtaHisto->SetMinimum(PlotMinDist); AN_2S_GammaEtaHisto->SetMaximum(AN_2S_GammaEtaHisto->GetMaximum()*MaxFact);
	  AN_2S_GammaPtHisto->SetMinimum(PlotMinDist); AN_2S_GammaPtHisto->SetMaximum(AN_2S_GammaPtHisto->GetMaximum()*MaxFact);
	  AN_2S_Pi0MassHistoBroad->SetMinimum(PlotMinDist); AN_2S_Pi0MassHistoBroad->SetMaximum(AN_2S_Pi0MassHistoBroad->GetMaximum()*MaxFact);
	  AN_2S_Pi0MassHistoNarrow->SetMinimum(PlotMinDist); AN_2S_Pi0MassHistoNarrow->SetMaximum(AN_2S_Pi0MassHistoNarrow->GetMaximum()*MaxFact);
	  AN_2S_RConvHisto->SetMinimum(PlotMinDist); AN_2S_RConvHisto->SetMaximum(AN_2S_RConvHisto->GetMaximum()*MaxFact);

	  AN_3S_GammaVertexProbHisto->SetMinimum(PlotMinDist*100); AN_3S_GammaVertexProbHisto->SetMaximum(AN_3S_GammaVertexProbHisto->GetMaximum()*MaxFact);
	  AN_3S_dzsigHisto->SetMinimum(PlotMinDist*100); AN_3S_dzsigHisto->SetMaximum(AN_3S_dzsigHisto->GetMaximum()*MaxFact);
	  AN_3S_dzHisto->SetMinimum(PlotMinDist); AN_3S_dzHisto->SetMaximum(AN_3S_dzHisto->GetMaximum()*MaxFact);
	  AN_3S_GammaEtaHisto->SetMinimum(PlotMinDist); AN_3S_GammaEtaHisto->SetMaximum(AN_3S_GammaEtaHisto->GetMaximum()*MaxFact);
	  AN_3S_GammaPtHisto->SetMinimum(PlotMinDist); AN_3S_GammaPtHisto->SetMaximum(AN_3S_GammaPtHisto->GetMaximum()*MaxFact);
	  AN_3S_Pi0MassHistoBroad->SetMinimum(PlotMinDist); AN_3S_Pi0MassHistoBroad->SetMaximum(AN_3S_Pi0MassHistoBroad->GetMaximum()*MaxFact);
	  AN_3S_Pi0MassHistoNarrow->SetMinimum(PlotMinDist); AN_3S_Pi0MassHistoNarrow->SetMaximum(AN_3S_Pi0MassHistoNarrow->GetMaximum()*MaxFact);
	  AN_3S_RConvHisto->SetMinimum(PlotMinDist); AN_3S_RConvHisto->SetMaximum(AN_3S_RConvHisto->GetMaximum()*MaxFact);
*/

	  MaxFact=1.;



    AN_MassHisto->GetYaxis()->SetTitleOffset(yOff); AN_MassHisto->SetStats(0);AN_MassHisto->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");	AN_MassHisto->GetYaxis()->SetTitle("Counts per 25 MeV"); AN_MassHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/AN_MassHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_UpsPtHisto->GetYaxis()->SetTitleOffset(yOff); AN_UpsPtHisto->SetStats(0);AN_UpsPtHisto->GetXaxis()->SetTitle("p_{T}^{#mu#mu} [GeV]");	AN_UpsPtHisto->GetYaxis()->SetTitle("Counts per 500 MeV"); AN_UpsPtHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
    sprintf(saveMass,"Figures/%s/AN_UpsPtHisto.pdf",FitID);
	cutLine->SetX1(cut_Ypt); cutLine->SetX2(cut_Ypt); cutLine->SetY2(AN_UpsPtHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
    PlotCanvas->SaveAs(saveMass);

    AN_UpsRapHisto->GetYaxis()->SetTitleOffset(yOff); AN_UpsRapHisto->SetStats(0);AN_UpsRapHisto->GetXaxis()->SetTitle("y_{#mu#mu}");	AN_UpsRapHisto->GetYaxis()->SetTitle("Counts"); AN_UpsRapHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
	cutLine->SetX1(-cut_rap); cutLine->SetX2(-cut_rap); cutLine->SetY2(AN_UpsRapHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(cut_rap); cutLine2->SetX2(cut_rap); cutLine2->SetY2(AN_UpsRapHisto->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_UpsRapHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_DimuonVertexProbHisto->GetYaxis()->SetTitleOffset(yOff); AN_DimuonVertexProbHisto->SetStats(0);AN_DimuonVertexProbHisto->GetXaxis()->SetTitle("#mu#mu-vertex #chi^{2} prob. [%]");	AN_DimuonVertexProbHisto->GetYaxis()->SetTitle("Counts"); AN_DimuonVertexProbHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/AN_DimuonVertexProbHisto.pdf",FitID);
	cutLine->SetX1(cut_vtxProb*100); cutLine->SetX2(cut_vtxProb*100); cutLine->SetY2(AN_DimuonVertexProbHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
    PlotCanvas->SaveAs(saveMass);





    AN_1S_GammaVertexProbHisto->GetYaxis()->SetTitleOffset(yOff); AN_1S_GammaVertexProbHisto->SetStats(0);AN_1S_GammaVertexProbHisto->GetXaxis()->SetTitle("conv.-vertex #chi^{2} prob. [%]");	AN_1S_GammaVertexProbHisto->GetYaxis()->SetTitle("Counts"); AN_1S_GammaVertexProbHisto->Draw("E");
    AN_1S_GammaVertexProbHisto->SetStats(0);
    AN_1S_GammaVertexProbHisto->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_GammaVertexProbHisto->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_GammaVertexProbHisto->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_GammaVertexProbHisto->Draw("E");
    AN_2S_GammaVertexProbHisto->SetStats(0);
    AN_2S_GammaVertexProbHisto->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_GammaVertexProbHisto->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_GammaVertexProbHisto->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_GammaVertexProbHisto->Draw("same,E");
	AN_3S_GammaVertexProbHisto->SetStats(0);
	AN_3S_GammaVertexProbHisto->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_GammaVertexProbHisto->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_GammaVertexProbHisto->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_GammaVertexProbHisto->Draw("same,E");

	TLegend* plotcompLegend=new TLegend(0.8,0.725,0.875,0.925);
	plotcompLegend->SetFillColor(0);
	plotcompLegend->SetTextFont(72);
	plotcompLegend->SetTextSize(0.04);
	plotcompLegend->SetBorderSize(0);
	char complegendentry[200];
	sprintf(complegendentry,"#Upsilon(1S)");
	plotcompLegend->AddEntry(AN_1S_GammaVertexProbHisto,complegendentry,"elp");
	sprintf(complegendentry,"#Upsilon(2S)");
	if(nState>1.5) plotcompLegend->AddEntry(AN_2S_GammaVertexProbHisto,complegendentry,"elp");
	sprintf(complegendentry,"#Upsilon(3S)");
	if(nState>2.5) plotcompLegend->AddEntry(AN_3S_GammaVertexProbHisto,complegendentry,"elp");

	if(nState>1.5) plotcompLegend->Draw();
	PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
    sprintf(saveMass,"Figures/%s/AN_GammaVertexProbHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_dzsigHisto->GetYaxis()->SetTitleOffset(yOff); AN_1S_dzsigHisto->SetStats(0);AN_1S_dzsigHisto->GetXaxis()->SetTitle("dz/#sigma_{dz}");	AN_1S_dzsigHisto->GetYaxis()->SetTitle("Counts"); AN_1S_dzsigHisto->Draw("E");
    AN_1S_dzsigHisto->SetStats(0);
    AN_1S_dzsigHisto->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_dzsigHisto->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_dzsigHisto->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_dzsigHisto->Draw("E");
    AN_2S_dzsigHisto->SetStats(0);
    AN_2S_dzsigHisto->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_dzsigHisto->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_dzsigHisto->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_dzsigHisto->Draw("same,E");
	AN_3S_dzsigHisto->SetStats(0);
	AN_3S_dzsigHisto->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_dzsigHisto->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_dzsigHisto->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_dzsigHisto->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
	cutLine->SetX1(-cut_dzSig); cutLine->SetX2(-cut_dzSig); cutLine->SetY2(AN_1S_dzsigHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(cut_dzSig); cutLine2->SetX2(cut_dzSig); cutLine2->SetY2(AN_1S_dzsigHisto->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_dzsigHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_dzHisto->GetYaxis()->SetTitleOffset(yOff); AN_1S_dzHisto->SetStats(0);AN_1S_dzHisto->GetXaxis()->SetTitle("dz [cm]");	AN_1S_dzHisto->GetYaxis()->SetTitle("Counts"); AN_1S_dzHisto->Draw("E");
    AN_1S_dzHisto->SetStats(0);
    AN_1S_dzHisto->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_dzHisto->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_dzHisto->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_dzHisto->Draw("E");
    AN_2S_dzHisto->SetStats(0);
    AN_2S_dzHisto->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_dzHisto->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_dzHisto->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_dzHisto->Draw("same,E");
	AN_3S_dzHisto->SetStats(0);
	AN_3S_dzHisto->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_dzHisto->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_dzHisto->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_dzHisto->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
	cutLine->SetX1(-cut_dz); cutLine->SetX2(-cut_dz); cutLine->SetY2(AN_1S_dzHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(cut_dz); cutLine2->SetX2(cut_dz); cutLine2->SetY2(AN_1S_dzHisto->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_dzHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_GammaEtaHisto->GetYaxis()->SetTitleOffset(yOff); AN_1S_GammaEtaHisto->SetStats(0);AN_1S_GammaEtaHisto->GetXaxis()->SetTitle("#eta_{#gamma}");	AN_1S_GammaEtaHisto->GetYaxis()->SetTitle("Counts"); AN_1S_GammaEtaHisto->Draw("E");
    AN_1S_GammaEtaHisto->SetStats(0);
    AN_1S_GammaEtaHisto->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_GammaEtaHisto->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_GammaEtaHisto->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_GammaEtaHisto->Draw("E");
    AN_2S_GammaEtaHisto->SetStats(0);
    AN_2S_GammaEtaHisto->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_GammaEtaHisto->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_GammaEtaHisto->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_GammaEtaHisto->Draw("same,E");
	AN_3S_GammaEtaHisto->SetStats(0);
	AN_3S_GammaEtaHisto->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_GammaEtaHisto->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_GammaEtaHisto->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_GammaEtaHisto->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
	cutLine->SetX1(-cut_gammaeta); cutLine->SetX2(-cut_gammaeta); cutLine->SetY2(AN_1S_GammaEtaHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(cut_gammaeta); cutLine2->SetX2(cut_gammaeta); cutLine2->SetY2(AN_1S_GammaEtaHisto->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_GammaEtaHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_GammaPtHisto->GetYaxis()->SetTitleOffset(yOff); AN_1S_GammaPtHisto->SetStats(0);AN_1S_GammaPtHisto->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");	AN_1S_GammaPtHisto->GetYaxis()->SetTitle("Counts per 150 MeV"); AN_1S_GammaPtHisto->Draw("E");
    AN_1S_GammaPtHisto->SetStats(0);
    AN_1S_GammaPtHisto->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_GammaPtHisto->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_GammaPtHisto->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_GammaPtHisto->Draw("E");
    AN_2S_GammaPtHisto->SetStats(0);
    AN_2S_GammaPtHisto->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_GammaPtHisto->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_GammaPtHisto->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_GammaPtHisto->Draw("same,E");
	AN_3S_GammaPtHisto->SetStats(0);
	AN_3S_GammaPtHisto->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_GammaPtHisto->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_GammaPtHisto->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_GammaPtHisto->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
	cutLine->SetX1(cut_gammapt); cutLine->SetX2(cut_gammapt); cutLine->SetY2(AN_1S_GammaPtHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_GammaPtHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_Pi0MassHistoBroad->GetYaxis()->SetTitleOffset(yOff); AN_1S_Pi0MassHistoBroad->SetStats(0);AN_1S_Pi0MassHistoBroad->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");	AN_1S_Pi0MassHistoBroad->GetYaxis()->SetTitle("Counts per 10 MeV"); AN_1S_Pi0MassHistoBroad->Draw("E");
    AN_1S_Pi0MassHistoBroad->SetStats(0);
    AN_1S_Pi0MassHistoBroad->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_Pi0MassHistoBroad->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_Pi0MassHistoBroad->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_Pi0MassHistoBroad->Draw("E");
    AN_2S_Pi0MassHistoBroad->SetStats(0);
    AN_2S_Pi0MassHistoBroad->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_Pi0MassHistoBroad->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_Pi0MassHistoBroad->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_Pi0MassHistoBroad->Draw("same,E");
	AN_3S_Pi0MassHistoBroad->SetStats(0);
	AN_3S_Pi0MassHistoBroad->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_Pi0MassHistoBroad->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_Pi0MassHistoBroad->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_Pi0MassHistoBroad->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
	cutLine->SetX1(Pi0_mean-cut_Pi0); cutLine->SetX2(Pi0_mean-cut_Pi0); cutLine->SetY2(AN_1S_Pi0MassHistoBroad->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(Pi0_mean+cut_Pi0); cutLine2->SetX2(Pi0_mean+cut_Pi0); cutLine2->SetY2(AN_1S_Pi0MassHistoBroad->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_Pi0MassHistoBroad.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_Pi0MassHistoNarrow->GetYaxis()->SetTitleOffset(yOff); AN_1S_Pi0MassHistoNarrow->SetStats(0);AN_1S_Pi0MassHistoNarrow->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");	AN_1S_Pi0MassHistoNarrow->GetYaxis()->SetTitle("Counts per 2 MeV"); AN_1S_Pi0MassHistoNarrow->Draw("E");
    AN_1S_Pi0MassHistoNarrow->SetStats(0);
    AN_1S_Pi0MassHistoNarrow->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_Pi0MassHistoNarrow->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_Pi0MassHistoNarrow->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_Pi0MassHistoNarrow->Draw("E");
    AN_2S_Pi0MassHistoNarrow->SetStats(0);
    AN_2S_Pi0MassHistoNarrow->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_Pi0MassHistoNarrow->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_Pi0MassHistoNarrow->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_Pi0MassHistoNarrow->Draw("same,E");
	AN_3S_Pi0MassHistoNarrow->SetStats(0);
	AN_3S_Pi0MassHistoNarrow->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_Pi0MassHistoNarrow->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_Pi0MassHistoNarrow->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_Pi0MassHistoNarrow->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
	cutLine->SetX1(Pi0_mean-cut_Pi0); cutLine->SetX2(Pi0_mean-cut_Pi0); cutLine->SetY2(AN_1S_Pi0MassHistoNarrow->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(Pi0_mean+cut_Pi0); cutLine2->SetX2(Pi0_mean+cut_Pi0); cutLine2->SetY2(AN_1S_Pi0MassHistoNarrow->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_Pi0MassHistoNarrow.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);


    TF1* Pi0Gauss = new TF1("Pi0Gauss","gaus",0.134-0.012,0.134+0.012);
    Pi0Gauss->SetParameter(0,AN_1S_Pi0MassHistoNarrow->GetEntries());
    Pi0Gauss->SetParameter(1,Pi0_mean);
    Pi0Gauss->SetParameter(2,cut_Pi0);
    AN_1S_Pi0MassHistoNarrow->Fit("Pi0Gauss","EFNR");
    Pi0Gauss->SetLineWidth(0.4);
    Pi0Gauss->SetLineColor(kRed);
    Pi0Gauss->SetLineStyle(1);

    double Pi0Gauss_p0 = Pi0Gauss->GetParameter(1);
    double err_Pi0Gauss_p0 = Pi0Gauss->GetParError(1);
    double Pi0Gauss_p1 = Pi0Gauss->GetParameter(2);
    double err_Pi0Gauss_p1 = Pi0Gauss->GetParError(2);
    double Pi0Gauss_chi2=Pi0Gauss->GetChisquare();
    double Pi0Gauss_NDF=Pi0Gauss->GetNDF();

    double highest=AN_1S_Pi0MassHistoNarrow->GetMaximum()*0.97;

    char text[500];
    sprintf(text,"#color[2]{Fitting Gauss:} #mu = %1.4f #pm %1.4f, #sigma = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d", Pi0Gauss_p0, err_Pi0Gauss_p0, Pi0Gauss_p1, err_Pi0Gauss_p1,Pi0Gauss_chi2,int(Pi0Gauss_NDF));
    TLatex textPi0Gauss1 = TLatex(0.045,highest,text);
    textPi0Gauss1.SetTextSize(0.03)                                                                                                                                                                                                                                             ;


    AN_1S_Pi0MassHistoNarrow->GetYaxis()->SetTitleOffset(yOff); AN_1S_Pi0MassHistoNarrow->SetStats(0);AN_1S_Pi0MassHistoNarrow->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");	AN_1S_Pi0MassHistoNarrow->GetYaxis()->SetTitle("Counts per 2 MeV"); AN_1S_Pi0MassHistoNarrow->Draw("E");
    AN_1S_Pi0MassHistoNarrow->SetStats(0);
    AN_1S_Pi0MassHistoNarrow->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_Pi0MassHistoNarrow->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_Pi0MassHistoNarrow->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_Pi0MassHistoNarrow->Draw("E");
    AN_2S_Pi0MassHistoNarrow->SetStats(0);
    AN_2S_Pi0MassHistoNarrow->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_Pi0MassHistoNarrow->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_Pi0MassHistoNarrow->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_Pi0MassHistoNarrow->Draw("same,E");
	AN_3S_Pi0MassHistoNarrow->SetStats(0);
	AN_3S_Pi0MassHistoNarrow->SetMarkerColor(MarkerColor_nS[3]);
	AN_3S_Pi0MassHistoNarrow->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_Pi0MassHistoNarrow->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_Pi0MassHistoNarrow->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
	cutLine->SetX1(Pi0_mean-cut_Pi0); cutLine->SetX2(Pi0_mean-cut_Pi0); cutLine->SetY2(AN_1S_Pi0MassHistoNarrow->GetMaximum()*MaxFact); cutLine->Draw( "same" );
	cutLine2->SetX1(Pi0_mean+cut_Pi0); cutLine2->SetX2(Pi0_mean+cut_Pi0); cutLine2->SetY2(AN_1S_Pi0MassHistoNarrow->GetMaximum()*MaxFact); cutLine2->Draw( "same" );
	Pi0Gauss->Draw("same");
    textPi0Gauss1.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(saveMass,"Figures/%s/AN_Pi0MassHistoNarrow_WithFit.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    AN_1S_RConvHisto->GetYaxis()->SetTitleOffset(yOff); AN_1S_RConvHisto->SetStats(0);AN_1S_RConvHisto->GetXaxis()->SetTitle("#rho_{conv.} [cm]");	AN_1S_RConvHisto->GetYaxis()->SetTitle("Counts per 5 mm"); AN_1S_RConvHisto->Draw();
    AN_1S_RConvHisto->SetStats(0);
    AN_1S_RConvHisto->SetMarkerColor(MarkerColor_nS[1]);
    AN_1S_RConvHisto->SetLineColor(MarkerColor_nS[1]);
    AN_1S_RConvHisto->SetMarkerStyle(MarkerStyle_nS[1]);
    AN_1S_RConvHisto->SetMarkerSize(MarkerSize_nS[1]);
    AN_1S_RConvHisto->Draw("E");
    AN_2S_RConvHisto->SetStats(0);
    AN_2S_RConvHisto->SetMarkerColor(MarkerColor_nS[2]);
    AN_2S_RConvHisto->SetLineColor(MarkerColor_nS[2]);
    AN_2S_RConvHisto->SetMarkerStyle(MarkerStyle_nS[2]);
    AN_2S_RConvHisto->SetMarkerSize(MarkerSize_nS[2]);
	if(nState>1.5) AN_2S_RConvHisto->Draw("same,E");
	AN_3S_RConvHisto->SetStats(0);
	AN_3S_RConvHisto->SetMarkerColor(MarkerColor_nS[3]);
    AN_3S_RConvHisto->SetLineColor(MarkerColor_nS[3]);
	AN_3S_RConvHisto->SetMarkerStyle(MarkerStyle_nS[3]);
	AN_3S_RConvHisto->SetMarkerSize(MarkerSize_nS[3]);
	if(nState>2.5) AN_3S_RConvHisto->Draw("same,E");
	if(nState>1.5) plotcompLegend->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
	cutLine->SetX1(cut_RconvMin); cutLine->SetX2(cut_RconvMin); cutLine->SetY2(AN_1S_RConvHisto->GetMaximum()*MaxFact); cutLine->Draw( "same" );
    sprintf(saveMass,"Figures/%s/AN_RConvHisto.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);










    GammaPtHisto->GetYaxis()->SetTitleOffset(yOff); GammaPtHisto->SetStats(0);GammaPtHisto->GetXaxis()->SetTitle("p_{T#gamma} [GeV]");	GammaPtHisto->GetYaxis()->SetTitle("Events per 150 MeV"); GammaPtHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
    sprintf(saveMass,"Figures/%s/NoInvmCut_GammaPt.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    GammaEtaHisto->GetYaxis()->SetTitleOffset(yOff); GammaEtaHisto->SetStats(0);GammaEtaHisto->GetXaxis()->SetTitle("#eta_{#gamma}");	GammaEtaHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/NoInvmCut_GammaEta.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    UpsPtHisto->GetYaxis()->SetTitleOffset(yOff); UpsPtHisto->SetStats(0);UpsPtHisto->GetXaxis()->SetTitle("p_{T#mu#mu} [GeV]"); UpsPtHisto->GetYaxis()->SetTitle("Events per 500 MeV");	UpsPtHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
    sprintf(saveMass,"Figures/%s/NoInvmCut_UpsPt.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    UpsRapHisto->GetYaxis()->SetTitleOffset(yOff); UpsRapHisto->SetStats(0);UpsRapHisto->GetXaxis()->SetTitle("y_{#mu#mu}");	UpsRapHisto->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/NoInvmCut_UpsRap.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    RConvHisto->GetYaxis()->SetTitleOffset(yOff); RConvHisto->SetStats(0);RConvHisto->GetXaxis()->SetTitle("#rho_{conv} [cm]");	RConvHisto->GetYaxis()->SetTitle("Conversions per 0.5 cm");RConvHisto->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/NoInvmCut_RConv.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

     MassHisto->GetYaxis()->SetTitleOffset(yOff); MassHisto->SetStats(0);MassHisto->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");	MassHisto->GetYaxis()->SetTitle("Events per 25 MeV");MassHisto->Draw("E");
     PlotCanvas->Modified();
     PlotCanvas->SetLogy(false);
     sprintf(saveMass,"Figures/%s/NoInvmCut_UpsMass.pdf",FitID);
     PlotCanvas->SaveAs(saveMass);


    GammaPtHistoInvmCut->GetYaxis()->SetTitleOffset(yOff); GammaPtHistoInvmCut->SetStats(0);GammaPtHistoInvmCut->GetXaxis()->SetTitle("p_{T#gamma} [GeV]");		GammaPtHistoInvmCut->GetYaxis()->SetTitle("Events per 50 MeV"); GammaPtHistoInvmCut->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
    sprintf(saveMass,"Figures/%s/InvmCut_GammaPt.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    GammaEtaHistoInvmCut->GetYaxis()->SetTitleOffset(yOff); GammaEtaHistoInvmCut->SetStats(0);GammaEtaHistoInvmCut->GetXaxis()->SetTitle("#eta_{#gamma}");	GammaEtaHistoInvmCut->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/InvmCut_GammaEta.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    UpsPtHistoInvmCut->GetYaxis()->SetTitleOffset(yOff); UpsPtHistoInvmCut->SetStats(0);UpsPtHistoInvmCut->GetXaxis()->SetTitle("p_{T#mu#mu} [GeV]");UpsPtHistoInvmCut->GetYaxis()->SetTitle("Events per 500 MeV");	UpsPtHistoInvmCut->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(true);
    sprintf(saveMass,"Figures/%s/InvmCut_UpsPt.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    UpsRapHistoInvmCut->GetYaxis()->SetTitleOffset(yOff); UpsRapHistoInvmCut->SetStats(0);UpsRapHistoInvmCut->GetXaxis()->SetTitle("y_{#mu#mu}");	UpsRapHistoInvmCut->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/InvmCut_UpsRap.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    RConvHistoInvmCut->GetYaxis()->SetTitleOffset(yOff); RConvHistoInvmCut->SetStats(0);RConvHistoInvmCut->GetXaxis()->SetTitle("#rho_{conv} [cm]");	RConvHistoInvmCut->GetYaxis()->SetTitle("Conversions per 0.5 cm"); RConvHistoInvmCut->Draw();
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/InvmCut_RConv.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    MassHistoInvmCut->GetYaxis()->SetTitleOffset(yOff); MassHistoInvmCut->SetStats(0);MassHistoInvmCut->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");MassHistoInvmCut->GetYaxis()->SetTitle("Events per 25 MeV");	MassHistoInvmCut->Draw("E");
    PlotCanvas->Modified();
    PlotCanvas->SetLogy(false);
    sprintf(saveMass,"Figures/%s/InvmCut_UpsMass.pdf",FitID);
    PlotCanvas->SaveAs(saveMass);

    delete PlotCanvas;
    delete MassHisto;
    delete GammaPtHisto;
    delete GammaEtaHisto;
    delete UpsPtHisto;
    delete UpsRapHisto;
    delete RConvHisto;
    delete MassHistoInvmCut;
    delete GammaPtHistoInvmCut;
    delete GammaEtaHistoInvmCut;
    delete UpsPtHistoInvmCut;
    delete UpsRapHistoInvmCut;
    delete RConvHistoInvmCut;

   }



    /////////////////////////////////////////////
    ///////////// THE Background Saga ///////////
    /////////////////////////////////////////////

    RooAbsPdf* background1S;
    RooAbsPdf* background2S;
    RooAbsPdf* background3S;


    // Ernest's ToyMC for background estimation:

    double bb_thresh=10.58;
    double chib1Pmin=9.85;
    double chib1Pmax=9.95;
    double chib2Pmin=10.15;
    double chib2Pmax=10.325;
    double chib3Pmin=10.4;
    double chib3Pmax=10.55;


        int nHistBins=100;
        float m_gamma = 0.;
        char DrawChar[500];
        char DrawChar1S[500];
        char DrawChar2S[500];
        char DrawChar3S[500];

        sprintf(DrawChar1S,"invm1S<%f | invm1S>%f&&invm1S<%f|invm1S>%f &&invm1S<%f|invm1S>%f",chib1Pmin,chib1Pmax,chib2Pmin,chib2Pmax,chib3Pmin,chib3Pmax);
        sprintf(DrawChar2S,"invm2S<%f|invm2S>%f &&invm2S<%f|invm2S>%f",chib2Pmin,chib2Pmax,chib3Pmin,chib3Pmax);
        sprintf(DrawChar3S,"invm3S<%f|invm3S>%f",chib3Pmin,chib3Pmax);
        if(useSBforBKGmodel){
        	sprintf(DrawChar1S,"invm1S<10000");
        	sprintf(DrawChar2S,"invm2S<10000");
        	sprintf(DrawChar3S,"invm3S<10000");
        }

        TH1F  *hYmass1S = new TH1F("hYmass1S","",nHistBins,8.7,11.4);
        TH1F  *hYmass_check1S = new TH1F("hYmass_check1S","",nHistBins,8.7,11.4);
        TH1F  *hYmass2S = new TH1F("hYmass2S","",nHistBins,8.7,11.4);
        TH1F  *hYmass_check2S = new TH1F("hYmass_check2S","",nHistBins,8.7,11.4);
        TH1F  *hYmass3S = new TH1F("hYmass3S","",nHistBins,8.7,11.4);
        TH1F  *hYmass_check3S = new TH1F("hYmass_check3S","",nHistBins,8.7,11.4);

        tree1S->Draw("jpsimass>>hYmass1S",DrawChar1S);
        tree2S->Draw("jpsimass>>hYmass2S",DrawChar2S);
        tree3S->Draw("jpsimass>>hYmass3S",DrawChar3S);

        TH1F  *hGammaP1S = new TH1F("hGammaP1S","",10*nHistBins,0,150);
        TH1F  *hGammaP_check1S = new TH1F("hGammaP_check1S","",10*nHistBins,0,100);
        TH1F  *hGammaP2S = new TH1F("hGammaP2S","",10*nHistBins,0,100);
        TH1F  *hGammaP_check2S = new TH1F("hGammaP_check2S","",10*nHistBins,0,100);
        TH1F  *hGammaP3S = new TH1F("hGammaP3S","",10*nHistBins,0,100);
        TH1F  *hGammaP_check3S = new TH1F("hGammaP_check3S","",10*nHistBins,0,100);

        tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP1S",DrawChar1S);
        tree2S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP2S",DrawChar2S);
        tree3S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP3S",DrawChar3S);

        TH1F  *hUpsP1S = new TH1F("hUpsP1S","",3*nHistBins,0,100);
        TH1F  *hUpsP_check1S = new TH1F("hUpsP_check1S","",3*nHistBins,0,100);
        TH1F  *hUpsP2S = new TH1F("hUpsP2S","",1.5*nHistBins,0,100);
        TH1F  *hUpsP_check2S = new TH1F("hUpsP_check2S","",1.5*nHistBins,0,100);
        TH1F  *hUpsP3S = new TH1F("hUpsP3S","",1.5*nHistBins,0,100);
        TH1F  *hUpsP_check3S = new TH1F("hUpsP_check3S","",1.5*nHistBins,0,100);

        tree1S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP1S",DrawChar1S);
        tree2S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP2S",DrawChar2S);
        tree3S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP3S",DrawChar3S);

        TH1F  *hCosAlphaP1S = new TH1F("hCosAlphaP1S","",1.5*nHistBins,-1,1);
        TH1F  *hCosAlphaP_check1S = new TH1F("hCosAlphaP_check1S","",1.5*nHistBins,-1,1);
        TH1F  *hCosAlphaP2S = new TH1F("hCosAlphaP2S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check2S = new TH1F("hCosAlphaP_check2S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP3S = new TH1F("hCosAlphaP3S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check3S = new TH1F("hCosAlphaP_check3S","",nHistBins,-1,1);

        tree1S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP1S",DrawChar1S);
        tree2S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP2S",DrawChar2S);
        tree3S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP3S",DrawChar3S);

        hYmass1S->Print();
        hGammaP1S->Print();
        hUpsP1S->Print();
        hCosAlphaP1S->Print();
        hYmass2S->Print();
        hGammaP2S->Print();
        hUpsP2S->Print();
        hCosAlphaP2S->Print();
        hYmass3S->Print();
        hGammaP3S->Print();
        hUpsP3S->Print();
        hCosAlphaP3S->Print();

        bool saveBkgToy=false;

        char bkgToyname[200];
        sprintf(bkgToyname,"%s/bkgToy.root",dirstruct);
        bool existingBKGfile=false;
        TFile *fQ = new TFile(bkgToyname,"RECREATE");
        TTree* tQ = new TTree("tQ","toy MC");


//        bool existingBKGfile=true;
//        TFile *fQ = new TFile("/Users/valentinknuenz/usr/local/workspace/Chib2011_/Chi2OniaGamma/analysis/bkgToy.root","READ");
//        TTree* tQ = (TTree*)fQ->Get("tQ");



        if(existingBKGfile&&BkgToy) {
        	cout<<"using existing bkg-toyMC file"<<endl;
        }

        if(!existingBKGfile&&BkgToy){
        	cout<<"producing new bkg-toyMC file"<<endl;

        hYmass1S->Write();
        hGammaP1S->Write();
        hUpsP1S->Write();
        hCosAlphaP1S->Write();
        hYmass2S->Write();
        hGammaP2S->Write();
        hUpsP2S->Write();
        hCosAlphaP2S->Write();
        hYmass3S->Write();
        hGammaP3S->Write();
        hUpsP3S->Write();
        hCosAlphaP3S->Write();


        float Q1S,Q2S,Q3S;
        float M1S,M2S,M3S;

        tQ->Branch("invm1S",&M1S,"invm1S/F");
        tQ->Branch("invm2S",&M2S,"invm2S/F");
        tQ->Branch("invm3S",&M3S,"invm3S/F");

        TRandom3 *randOM = new TRandom3();

        for(int i=0;i<nToy;i++){

          float Ymass_1S = hYmass1S->GetRandom();
          float p_gamma1S = hGammaP1S->GetRandom();
          float p_Ups1S = hUpsP1S->GetRandom();

          if(!alteredToy) Ymass_1S=Ymass1S;

//          float cosAlpha1S = 2*randOM->Uniform()-1;
          float cosAlpha1S = hCosAlphaP1S->GetRandom();
          float e_gamma1S = sqrt(m_gamma*m_gamma+p_gamma1S*p_gamma1S);
          float e_Ups1S = sqrt(Ymass_1S*Ymass_1S+p_Ups1S*p_Ups1S);

//          if(!alteredToy) Ymass_1S=Ymass1S;

          M1S = sqrt(m_gamma*m_gamma + Ymass_1S*Ymass_1S + 2*e_gamma1S*e_Ups1S - 2*p_gamma1S*p_Ups1S*cosAlpha1S)-Ymass_1S+Ymass1S;

       	  hUpsP_check1S->Fill(p_Ups1S);
          hGammaP_check1S->Fill(p_gamma1S);
          hYmass_check1S->Fill(Ymass_1S);
          hCosAlphaP_check1S->Fill(cosAlpha1S);

          float Ymass_2S = hYmass2S->GetRandom();
          float p_gamma2S = hGammaP2S->GetRandom();
          float p_Ups2S = hUpsP2S->GetRandom();

          if(!alteredToy) Ymass_2S=Ymass2S;

//          float cosAlpha2S = 2*randOM->Uniform()-1;
          float cosAlpha2S = hCosAlphaP2S->GetRandom();
          float e_gamma2S = sqrt(m_gamma*m_gamma+p_gamma2S*p_gamma2S);
          float e_Ups2S = sqrt(Ymass_2S*Ymass_2S+p_Ups2S*p_Ups2S);

          M2S = sqrt(m_gamma*m_gamma + Ymass_2S*Ymass_2S + 2*e_gamma2S*e_Ups2S - 2*p_gamma2S*p_Ups2S*cosAlpha2S)-Ymass_2S+Ymass2S;

       	  hUpsP_check2S->Fill(p_Ups2S);
          hGammaP_check2S->Fill(p_gamma2S);
          hYmass_check2S->Fill(Ymass_2S);
          hCosAlphaP_check2S->Fill(cosAlpha2S);

          float Ymass_3S = hYmass3S->GetRandom();
          float p_gamma3S = hGammaP3S->GetRandom();
          float p_Ups3S = hUpsP3S->GetRandom();

          if(!alteredToy) Ymass_3S=Ymass3S;

//          float cosAlpha3S = 2*randOM->Uniform()-1;
          float cosAlpha3S = hCosAlphaP3S->GetRandom();
          float e_gamma3S = sqrt(m_gamma*m_gamma+p_gamma3S*p_gamma3S);
          float e_Ups3S = sqrt(Ymass_3S*Ymass_3S+p_Ups3S*p_Ups3S);

          M3S = sqrt(m_gamma*m_gamma + Ymass_3S*Ymass_3S + 2*e_gamma3S*e_Ups3S - 2*p_gamma3S*p_Ups3S*cosAlpha3S)-Ymass_3S+Ymass3S;

       	  hUpsP_check3S->Fill(p_Ups3S);
          hGammaP_check3S->Fill(p_gamma3S);
          hYmass_check3S->Fill(Ymass_3S);
          hCosAlphaP_check3S->Fill(cosAlpha3S);

          if (M1S<11.4 && M2S<11.4 && M3S<11.4) {
        	  tQ->Fill();
           	if (i%10000==9999) cout << i+1 << endl;
          }

          else i--;


        }



        hCosAlphaP_check1S->Write();
        hUpsP_check1S->Write();
        hGammaP_check1S->Write();
        hYmass_check1S->Write();
        hCosAlphaP_check2S->Write();
        hUpsP_check2S->Write();
        hGammaP_check2S->Write();
        hYmass_check2S->Write();
        hCosAlphaP_check3S->Write();
        hUpsP_check3S->Write();
        hGammaP_check3S->Write();
        hYmass_check3S->Write();

        tQ->Write();
    	cout<<"finalising bkg-toyMC file"<<endl;
        delete randOM;


        }





        if(SaveAll&&BkgToy){
              char saveBkgToy[200];
                 TCanvas* BkgToyCanvas = new TCanvas("BkgToyCanvas","BkgToyCanvas",1600, 800);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hCosAlphaP1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCosAlphaP1S->SetStats(0);hCosAlphaP1S->GetXaxis()->SetTitle("cos#alpha");	hCosAlphaP1S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hCosAlphaP1S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hCosAlphaP_check1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCosAlphaP_check1S->SetStats(0);hCosAlphaP_check1S->GetXaxis()->SetTitle("cos#alpha");	hCosAlphaP_check1S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hCosAlphaP_check1S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hUpsP1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hUpsP1S->SetStats(0);hUpsP1S->GetXaxis()->SetTitle("p_{Y(1S)}");	hUpsP1S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hUpsP1S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hUpsP_check1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hUpsP_check1S->SetStats(0);hUpsP_check1S->GetXaxis()->SetTitle("p_{Y(1S)}");	hUpsP_check1S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hUpsP_check1S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hGammaP1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hGammaP1S->SetStats(0);hGammaP1S->GetXaxis()->SetTitle("p_{#gamma}");	hGammaP1S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hGammaP1S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hGammaP_check1S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hGammaP_check1S->SetStats(0);hGammaP_check1S->GetXaxis()->SetTitle("p_{#gamma}");	hGammaP_check1S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hGammaP_check1S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);


                 BkgToyCanvas->SetFillColor(kWhite);
                 hCosAlphaP2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCosAlphaP2S->SetStats(0);hCosAlphaP2S->GetXaxis()->SetTitle("cos#alpha");	hCosAlphaP2S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hCosAlphaP2S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hCosAlphaP_check2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCosAlphaP_check2S->SetStats(0);hCosAlphaP_check2S->GetXaxis()->SetTitle("cos#alpha");	hCosAlphaP_check2S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hCosAlphaP_check2S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hUpsP2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hUpsP2S->SetStats(0);hUpsP2S->GetXaxis()->SetTitle("p_{Y(2S)}");	hUpsP2S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hUpsP2S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hUpsP_check2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hUpsP_check2S->SetStats(0);hUpsP_check2S->GetXaxis()->SetTitle("p_{Y(2S)}");	hUpsP_check2S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hUpsP_check2S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hGammaP2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hGammaP2S->SetStats(0);hGammaP2S->GetXaxis()->SetTitle("p_{#gamma}");	hGammaP2S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hGammaP2S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hGammaP_check2S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hGammaP_check2S->SetStats(0);hGammaP_check2S->GetXaxis()->SetTitle("p_{#gamma}");	hGammaP_check2S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hGammaP_check2S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hCosAlphaP3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCosAlphaP3S->SetStats(0);hCosAlphaP3S->GetXaxis()->SetTitle("cos#alpha");	hCosAlphaP3S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hCosAlphaP3S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hCosAlphaP_check3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCosAlphaP_check3S->SetStats(0);hCosAlphaP_check3S->GetXaxis()->SetTitle("cos#alpha");	hCosAlphaP_check3S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hCosAlphaP_check3S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hUpsP3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hUpsP3S->SetStats(0);hUpsP3S->GetXaxis()->SetTitle("p_{Y(3S)}");	hUpsP3S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hUpsP3S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hUpsP_check3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hUpsP_check3S->SetStats(0);hUpsP_check3S->GetXaxis()->SetTitle("p_{Y(3S)}");	hUpsP_check3S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hUpsP_check3S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

                 BkgToyCanvas->SetFillColor(kWhite);
                 hGammaP3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hGammaP3S->SetStats(0);hGammaP3S->GetXaxis()->SetTitle("p_{#gamma}");	hGammaP3S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hGammaP3S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);
                 BkgToyCanvas->SetFillColor(kWhite);
                 hGammaP_check3S->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hGammaP_check3S->SetStats(0);hGammaP_check3S->GetXaxis()->SetTitle("p_{#gamma}");	hGammaP_check3S->Draw();
                 BkgToyCanvas->Modified();
                 sprintf(saveBkgToy,"Figures/%s/hGammaP_check3S.pdf",FitID);
                 BkgToyCanvas->SaveAs(saveBkgToy);

        }





        TH1F  *BkgToyHist1S = new TH1F("BkgToyHist1S","",100,9.5,12);
        tQ->Draw("invm1S>>BkgToyHist1S");
        TH1F  *BkgToyHist2S = new TH1F("BkgToyHist2S","",100,9.5,12);
        tQ->Draw("invm2S>>BkgToyHist2S");
        TH1F  *BkgToyHist3S = new TH1F("BkgToyHist3S","",100,9.5,12);
        tQ->Draw("invm3S>>BkgToyHist3S");


        invm1S.setBins(50);
        invm2S.setBins(50);
        invm3S.setBins(50);

        RooDataSet* RooBkgToySet = new RooDataSet("RooBkgToySet","RooBkgToySet",tQ,RooArgList(invm1S,invm2S,invm3S));
        RooBkgToySet->Print();
        RooDataHist* RooBkgToyHist1S = new RooDataHist("RooBkgToyHist1S","RooBkgToyHist1S",RooArgList(invm1S),*RooBkgToySet);
        RooBkgToyHist1S->Print();
        if(BkgToy) background1S = new RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S,1);
        BkgToyHist1S->Write();
        RooDataHist* RooBkgToyHist2S = new RooDataHist("RooBkgToyHist2S","RooBkgToyHist2S",RooArgList(invm2S),*RooBkgToySet);
        RooBkgToyHist2S->Print();
        if(BkgToy) background2S = new RooHistPdf("background2S","background2S",RooArgSet(invm2S),*RooBkgToyHist2S,1);
        BkgToyHist2S->Write();
        RooDataHist* RooBkgToyHist3S = new RooDataHist("RooBkgToyHist3S","RooBkgToyHist3S",RooArgList(invm3S),*RooBkgToySet);
        RooBkgToyHist3S->Print();
        if(BkgToy) background3S = new RooHistPdf("background3S","background3S",RooArgSet(invm3S),*RooBkgToyHist3S,1);
        BkgToyHist3S->Write();

        fQ->Write();
        fQ->Close();



        delete hYmass1S;
        delete hGammaP1S;
        delete hUpsP1S;
        delete hCosAlphaP1S;
        delete hYmass2S;
        delete hGammaP2S;
        delete hUpsP2S;
        delete hCosAlphaP2S;
        delete hCosAlphaP_check1S;
        delete hUpsP_check1S;
        delete hGammaP_check1S;
        delete hYmass_check1S;
        delete hCosAlphaP_check2S;
        delete hUpsP_check2S;
        delete hGammaP_check2S;
        delete hYmass_check2S;

        delete hYmass3S;
        delete hGammaP3S;
        delete hUpsP3S;
        delete hCosAlphaP3S;

        delete hCosAlphaP_check3S;
        delete hUpsP_check3S;
        delete hGammaP_check3S;
        delete hYmass_check3S;


////////////////// New Bkg mixer /////////////////

        if(BkgMixerBool&&!useNewExistingMixFile){
        	cout<<"Using BkgMixer tool to build background shape"<<endl;
     	    gROOT->ProcessLine(".L BkgMixer.C++");
     	    gROOT->ProcessLine(".L BkgMixer.C++");
     	    gROOT->ProcessLine(".L BkgMixer.C++");
     	    char MixerInputTree[200];
     	    int nStateToMix;
     	    char ExcludeSignal[200];
     	    int nBinsCosAlphaTry;

     	    //useExistingMixFile=true;///
     	    nStateToMix=1;
     	    nBinsCosAlphaTry=1000;
     	    sprintf(ExcludeSignal,"invm1S<%f | invm1S>%f&&invm1S<%f|invm1S>%f &&invm1S<%f|invm1S>%f",chib1Pmin,chib1Pmax,chib2Pmin,chib2Pmax,chib3Pmin,chib3Pmax);
     	    BkgMixer(tree1SExceptGammaPtCut,nStateToMix,invm_min1S,invm_max1S,useExistingMixFile,nMix,dirstruct,ExcludeSignal,cut_gammapt,gammaptK,gammaptD,Ymass1S,nBinsCosAlphaTry,GammaEtaBorder, cut_gammapt_Midrap, cut_gammapt_Forward , cut_Ypt, cut_gammaeta);
     	    TH1F *hInvM1S=(TH1F*)hInvMnS->Clone("hInvM1S");
     	    hInvMnS->Print();
     	    hInvM1S->Print();

          RooDataHist* RooBkgToyHist1S_Mixer = new RooDataHist("RooBkgToyHist1S_Mixer","RooBkgToyHist1S_Mixer",RooArgList(invm1S),hInvM1S);
          RooBkgToyHist1S_Mixer->Print();
          if(!BkgToy1SBkgMixer2S) background1S = new  RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S_Mixer,1);

          char copyFrom[200];
          char copyTo[200];
          if(Optimize){
              cout<<"Saving eventMixer file for 1S background..."<<endl;
        	  sprintf(copyFrom,"%s/eventMixer1S.root",dirstruct);
        	  sprintf(copyTo,"%s/eventMixer1S_Cut%d.root",dirstruct,inCut);
        	  gSystem->CopyFile(copyFrom,copyTo,kTRUE);
          }

          fMix->Close();
          cout<<"Finished Bkg mixing for 1S..."<<endl;

          if(nState<1.5) background2S=background1S;
          if(nState<2.5) background3S=background1S;

          TH1F *hInvM2S;
          if(nState>1.5){
          //useExistingMixFile=true;///
          nStateToMix=2;
   	      nBinsCosAlphaTry=1000;
		  sprintf(ExcludeSignal,"invm2S>%f && invm2S<%f|invm2S>%f",chib2Pmax,chib3Pmin,chib3Pmax);
          BkgMixer(tree2SExceptGammaPtCut,nStateToMix,invm_min2S,invm_max2S,useExistingMixFile,nMix,dirstruct,ExcludeSignal,cut_gammapt,gammaptK,gammaptD,Ymass2S,nBinsCosAlphaTry, GammaEtaBorder, cut_gammapt_Midrap, cut_gammapt_Forward , cut_Ypt, cut_gammaeta);
          hInvM2S=(TH1F*)hInvMnS->Clone("hInvM2S");
          hInvMnS->Print();
          hInvM2S->Print();

          if(Optimize){
              cout<<"Saving eventMixer file for 2S background..."<<endl;
        	  sprintf(copyFrom,"%s/eventMixer2S.root",dirstruct);
        	  sprintf(copyTo,"%s/eventMixer2S_Cut%d.root",dirstruct,inCut);
        	  gSystem->CopyFile(copyFrom,copyTo,kTRUE);
          }
          cout<<"Finished Bkg mixing for 2S..."<<endl;

          hInvM2S=(TH1F*)hInvMnS->Clone("hInvM2S");
          RooDataHist* RooBkgToyHist2S_Mixer = new RooDataHist("RooBkgToyHist2S_Mixer","RooBkgToyHist2S_Mixer",RooArgList(invm2S),hInvM2S);
          RooBkgToyHist2S_Mixer->Print();
          background2S = new  RooHistPdf("background2S","background2S",RooArgSet(invm2S),*RooBkgToyHist2S_Mixer,1);
          fMix->Close();

          }


          TH1F *hInvM3S;
          if(nState>2.5){
          //useExistingMixFile=false;///
          nStateToMix=3;
   	      nBinsCosAlphaTry=1000;
		  sprintf(ExcludeSignal,"invm3S<%f|invm3S>%f",chib3Pmin,chib3Pmax);
          BkgMixer(tree3SExceptGammaPtCut,nStateToMix,invm_min3S,invm_max3S,useExistingMixFile,nMix,dirstruct,ExcludeSignal,cut_gammapt,gammaptK,gammaptD,Ymass3S,nBinsCosAlphaTry, GammaEtaBorder, cut_gammapt_Midrap, cut_gammapt_Forward , cut_Ypt, cut_gammaeta);
          hInvM3S=(TH1F*)hInvMnS->Clone("hInvM3S");
          hInvMnS->Print();
          hInvM3S->Print();

          if(Optimize){
              cout<<"Saving eventMixer file for 2S background..."<<endl;
        	  sprintf(copyFrom,"%s/eventMixer3S.root",dirstruct);
        	  sprintf(copyTo,"%s/eventMixer3S_Cut%d.root",dirstruct,inCut);
        	  gSystem->CopyFile(copyFrom,copyTo,kTRUE);
          }
          cout<<"Finished Bkg mixing for 3S..."<<endl;
          hInvM3S=(TH1F*)hInvMnS->Clone("hInvM3S");
          RooDataHist* RooBkgToyHist3S_Mixer = new RooDataHist("RooBkgToyHist3S_Mixer","RooBkgToyHist3S_Mixer",RooArgList(invm3S),hInvM3S);
          RooBkgToyHist3S_Mixer->Print();
          background3S = new  RooHistPdf("background3S","background3S",RooArgSet(invm3S),*RooBkgToyHist3S_Mixer,1);

          fMix->Close();
          }


        }

	      if(useNewExistingMixFile){
	 	        TFile *NewBkgMixerFile;
	 	        NewBkgMixerFile= new TFile("bkgMixer.root","READ");

		          TH1F *hInvM1S;
		          TH1F *hInvM2S;
		          TH1F *hInvM3S;

	    	hInvM1S=(TH1F*)NewBkgMixerFile->Get("hInvM1S");
	    	hInvM1S->Print();
	        RooDataHist* RooBkgToyHist1S_Mixer = new RooDataHist("RooBkgToyHist1S_Mixer","RooBkgToyHist1S_Mixer",RooArgList(invm1S),hInvM1S);
	        RooBkgToyHist1S_Mixer->Print();
	        background1S = new  RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S_Mixer,1);

	        if(nState<1.5) background2S=background1S;
	        if(nState<2.5) background3S=background1S;

	        if(nState>1.5){
		    	hInvM2S=(TH1F*)NewBkgMixerFile->Get("hInvM2S");
		    	hInvM2S->Print();
		        RooDataHist* RooBkgToyHist2S_Mixer = new RooDataHist("RooBkgToyHist2S_Mixer","RooBkgToyHist2S_Mixer",RooArgList(invm2S),hInvM2S);
		        RooBkgToyHist2S_Mixer->Print();
		        background2S = new  RooHistPdf("background2S","background2S",RooArgSet(invm2S),*RooBkgToyHist2S_Mixer,1);
	        }

	        if(nState>2.5){
		    	hInvM3S=(TH1F*)NewBkgMixerFile->Get("hInvM3S");
		    	hInvM3S->Print();
		        RooDataHist* RooBkgToyHist3S_Mixer = new RooDataHist("RooBkgToyHist3S_Mixer","RooBkgToyHist3S_Mixer",RooArgList(invm3S),hInvM3S);
		        RooBkgToyHist3S_Mixer->Print();
		        background3S = new  RooHistPdf("background3S","background3S",RooArgSet(invm3S),*RooBkgToyHist3S_Mixer,1);
	        }

	        NewBkgMixerFile->Close();
	    }

        cout<<"Fully finished Bkg mixing..."<<endl;

        /*
        if(BkgMixerBool){

      	  TH1F *hInvM1S = new TH1F("hInvM1S",";Q + M(#Upsilon(1S)) [GeV]",100,invm_min1S-0.1,invm_max1S+0.1);
      	TFile *fMix;

        	if(!useExistingMixFile){

        	  int n_evt=nMix;

        	  TFile *f = new TFile("BkgMixTree/tree_Y1S_Pt0.root");
        	  TTree *t = (TTree*)f->Get("d");
              sprintf(bkgToyname,"%s/eventMixer.root",dirstruct);
        	  fMix = new TFile(bkgToyname,"recreate");

        	  char gammaPtCutChar[200];
        	  sprintf(gammaPtCutChar,"gammapt>%f && gammapt>%f *Q+%f",cut_gammapt,gammaptK,gammaptD);
        	  char ExcludeSignal[200];
        	  sprintf(ExcludeSignal,"invm1S<%f | invm1S>%f&&invm1S<%f|invm1S>%f &&invm1S<10.4|invm1S>10.55",chib1Pmin,chib1Pmax,chib2Pmin,chib2Pmax);
        	  TTree *t2 = (TTree*)tree1S->CopyTree(ExcludeSignal);




              const int nBinsCosAlphaDataHist=25;

              TH1F *hWNonEqu = new TH1F("hWNonEqu","; weights for cos#alpha",nBinsCosAlphaDataHist,-1,1);
              TH1F *hCosAlphaNonEqu = new TH1F("hCosAlphaNonEqu",";cos#alpha",nBinsCosAlphaDataHist,-1,1);
              TH1F *hCosAlphaDataNonEqu = new TH1F("hCosAlphaDataNonEqu",";cos#alpha",nBinsCosAlphaDataHist,-1,1);

              cout<<"Begin weighting of cosAlpha histo"<<endl;


              TH1F *hCosAlphaDataBuffer;
              int nBinsCosAlphaTry=1000;
              double BinBordersCosAlpha[nBinsCosAlphaDataHist+1];
              int n_t2 = t2->GetEntries(ExcludeSignal);
              int targetNumEventsPerBin=n_t2/nBinsCosAlphaDataHist;

              BinBordersCosAlpha[0]=-1.;
              int iRealBin=1;
              for(int iBins=1;iBins<nBinsCosAlphaTry+1;iBins++){
              hCosAlphaDataBuffer = new TH1F("hCosAlphaDataBuffer",";cos#alpha",50,-1,1);
              double cosAlphaBufferCut=-1+2.*double(iBins)/double(nBinsCosAlphaTry);
              char AlphaCut[200];
              sprintf(AlphaCut,"(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz) < %f",cosAlphaBufferCut);
              t2->Draw("(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz)>>hCosAlphaDataBuffer",AlphaCut);

              if(hCosAlphaDataBuffer->GetEntries()>iRealBin*targetNumEventsPerBin){
      //        	cout<<"hCosAlphaDataBuffer->GetEntries() = "<<hCosAlphaDataBuffer->GetEntries()<<endl;
              	cout<<"iCosAlphaBin = "<<iRealBin<<" / "<<nBinsCosAlphaDataHist<<endl;
              	BinBordersCosAlpha[iRealBin]=cosAlphaBufferCut-1./double(nBinsCosAlphaTry);
      //        	cout<<"BinBordersCosAlpha[iRealBin] = "<<BinBordersCosAlpha[iRealBin]<<endl;
              	iRealBin++;
              }
              if(iRealBin>nBinsCosAlphaDataHist-1) break;
              delete hCosAlphaDataBuffer;
              }
              BinBordersCosAlpha[nBinsCosAlphaDataHist]=1.;

              double xBinsCosAlpha[nBinsCosAlphaDataHist+1];

              for(int iBinCosAlpha=0;iBinCosAlpha<nBinsCosAlphaDataHist+1;iBinCosAlpha++){
              	xBinsCosAlpha[iBinCosAlpha]=BinBordersCosAlpha[iBinCosAlpha];
      //        	cout<<iBinCosAlpha<<" xBinsCosAlpha[iBinCosAlpha]"<<xBinsCosAlpha[iBinCosAlpha]<<endl;
              }


              hCosAlphaDataNonEqu->GetXaxis()->Set(nBinsCosAlphaDataHist, xBinsCosAlpha);
              hCosAlphaNonEqu->GetXaxis()->Set(nBinsCosAlphaDataHist, xBinsCosAlpha);
              hWNonEqu->GetXaxis()->Set(nBinsCosAlphaDataHist, xBinsCosAlpha);

              t2->Draw("(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz)>>hCosAlphaDataNonEqu");




        	  double jpsimass,jpsipx,jpsipy,jpsipz,gammapt_,gammapx,gammapy,gammapz;

        	  t->SetBranchAddress("jpsimass",&jpsimass);
        	  t->SetBranchAddress("jpsipx",&jpsipx);
        	  t->SetBranchAddress("jpsipy",&jpsipy);
        	  t->SetBranchAddress("jpsipz",&jpsipz);
        	  t->SetBranchAddress("gammapt",&gammapt_);
        	  t->SetBranchAddress("gammapx",&gammapx);
        	  t->SetBranchAddress("gammapy",&gammapy);
        	  t->SetBranchAddress("gammapz",&gammapz);

//              RooRealVar invm1S_write = RooRealVar("invm1S_write", "invm1S_write",0,100);
//              RooRealVar gammapt_write = RooRealVar("gammapt_write", "gammapt_write",0,100);
//              RooRealVar cosalpha_write = RooRealVar("cosalpha_write", "cosalpha_write",-1,1);
//              RooRealVar Q_write = RooRealVar("Q_write", "Q_write",0,100);
//              RooArgSet argSet = RooArgSet();
//              argSet.add(invm1S_write);
//              argSet.add(gammapt_write);
//              argSet.add(cosalpha_write);
//              argSet.add(Q_write);
//
//              RooDataSet rds = RooDataSet("d","d",argSet);


        	  int i_evt=0;
        	  int bin;

        	  TRandom3 *r = new TRandom3();

        	  int n = t->GetEntries();

        	  for(int iIter=1;iIter<3;iIter++){

        		  if(iIter==1) cout<<"Start mixing..."<<endl;

        		  if(iIter==2){
            		for(int iDivide=1;iDivide<nBinsCosAlphaDataHist+1;iDivide++){
//            			cout<<"dataHist "<<hCosAlphaDataNonEqu->GetBinContent(iDivide)<<endl;
//            			cout<<" mixHist "<<hCosAlphaNonEqu->GetBinContent(iDivide)<<endl;
            		hWNonEqu->SetBinContent(iDivide,hCosAlphaDataNonEqu->GetBinContent(iDivide)/hCosAlphaNonEqu->GetBinContent(iDivide));
            		}
            		cout<<"Finish weighting of cosAlpha histo"<<endl;
        		  }

        		  i_evt=0;

        	  while(i_evt<n_evt){
        	    double m_jpsimass,m_jpsipx,m_jpsipy,m_jpsipz;
        	    int i=n*r->Uniform();
        	    if (i==n) i=0;
        	    t->GetEntry(i);
        	    m_jpsimass = jpsimass;
        	    m_jpsipx = jpsipx;
        	    m_jpsipy = jpsipy;
        	    m_jpsipz = jpsipz;
        	    double m_jpsip = sqrt(m_jpsipx*m_jpsipx + m_jpsipy*m_jpsipy + m_jpsipz*m_jpsipz);
        	    double m_jpsie = sqrt(m_jpsimass*m_jpsimass + m_jpsip*m_jpsip);
        	    int j=n*r->Uniform();
        	    if (j==n) j=0;
        	    if (j!=i){
        	      t->GetEntry(j);
        	      double gammae = sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz);
        	      double chibe = gammae + m_jpsie;
        	      double chibpx = gammapx + m_jpsipx;
        	      double chibpy = gammapy + m_jpsipy;
        	      double chibpz = gammapz + m_jpsipz;
        	      double chibmass = sqrt(chibe*chibe - chibpx*chibpx - chibpy*chibpy - chibpz*chibpz);
        	      double Q = chibmass - m_jpsimass;
        	      if (gammapt_>cut_gammapt&&gammapt_>gammaptK*Q+gammaptD){
        		double cosalpha = (m_jpsipx*gammapx + m_jpsipy*gammapy + m_jpsipz*gammapz)/gammae/m_jpsip;
        		double chibmass_corr = Q + Ymass1S;
        		if(iIter==1) hCosAlphaNonEqu->Fill(cosalpha);
        		if(iIter==2) {
         	    	bin = hWNonEqu->FindBin(cosalpha);
              	    hInvM1S->Fill(chibmass_corr,hWNonEqu->GetBinContent(bin));
//              	    cout<<"chibmass_corr "<<chibmass_corr<<endl;
//              	    cout<<"bin "<<bin<<endl;
//              	    cout<<"hWNonEqu->GetBinContent(bin) "<<hWNonEqu->GetBinContent(bin)<<endl;
        		}

       		if (i_evt%10000==9999) cout << i_evt+1 << endl;
        		i_evt ++;
        	      }
        	    }
        	  }
        	  }

          	  hInvM1S->Write();
          	  hCosAlphaNonEqu->Write();
          	  hCosAlphaDataNonEqu->Write();
          	  hWNonEqu->Write();


        	}

        	else {
        		sprintf(bkgToyname,"%s/eventMixer.root",dirstruct);
        		fMix = new TFile(bkgToyname,"update");
              	hInvM1S=(TH1F*)fMix->Get("hInvM1S");
        	}




//        	hInvM1S->Rebin(2);

        RooDataHist* RooBkgToyHist1S_Mixer = new RooDataHist("RooBkgToyHist1S_Mixer","RooBkgToyHist1S_Mixer",RooArgList(invm1S),hInvM1S);
        RooBkgToyHist1S_Mixer->Print();
        background1S = new  RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S_Mixer,1);

      	fMix->Close();



        }
*/

////////////////// End New Bkg mixer /////////////////




        char BGformula[200];
/*        RooRealVar BGalpha_1S= RooRealVar("#alpha_{BG,1S}","",1.,0.,3.);
        RooRealVar BGbeta_1S= RooRealVar("#beta_{BG,1S}","",-1.,-3.,0.);
        sprintf(BGformula,"pow((@0-%f),@1)*exp((@0-%f)*@2)",q01S,q01S);
        RooGenericPdf background1S=RooGenericPdf("background1S","background1S",BGformula,RooArgList(invm1S,BGalpha_1S,BGbeta_1S));
*/

        double q01S_Start=9.55;
        RooRealVar alpha1("alpha1","alpha1",1.5,0.2,3.5);
        RooRealVar beta1("beta1","beta1",-2.5,-7.,0.);
        RooRealVar gamma1("gamma1","gamma1",0.,0.,0.2);
        RooRealVar q01S("q01S","q01S",q01S_Start,9.5,9.7);
        RooFormulaVar a1("a1","TMath::Abs(@0-@1)",RooArgList(invm1S,q01S));
        RooFormulaVar b1("b1","@0*(@1-@2)",RooArgList(beta1,invm1S,q01S));
        RooFormulaVar c1("c1","TMath::Abs(@0*(@1-@2))",RooArgList(gamma1,invm1S,q01S));
        RooFormulaVar signum1("signum1","(TMath::Sign(-1.,@0-@1)+1)/2.",RooArgList(invm1S,q01S));
        RooFormulaVar signum1_2("signum1_2","(TMath::Sign(-1.,@4*pow(@0,@1)*exp(@2)+@3)+1)/2.",RooArgList(a1,alpha1,b1,c1,signum1));

//        RooGenericPdf background1S("background1S","(signum1*pow(a1,alpha1)*exp(b1)+c1)*signum1_2",RooArgSet(signum1,a1,alpha1,b1,c1,signum1_2));

	if(useAnalyticalBKG) background1S = new RooGenericPdf("background1S","signum1*pow(a1,alpha1)*exp(b1)",RooArgSet(signum1,a1,alpha1,b1));

//        q01S.setConstant();
//        gamma1.setConstant();

        double q02S_Start=10.25;
        RooRealVar alpha2("alpha2","alpha2",0.6,0.1,15.);
        RooRealVar beta2("beta2","beta2",-2.5,-15.,5.);
        RooRealVar gamma2("gamma2","gamma2",0.,-0.1,0.1);
        RooRealVar q02S("q02S","q02S",q02S_Start,10.10,10.4);
        RooFormulaVar a2("a2","TMath::Abs(@0-@1)",RooArgList(invm2S,q02S));
        RooFormulaVar b2("b2","@0*(@1-@2)",RooArgList(beta2,invm2S,q02S));
        RooFormulaVar c2("c2","@0*(@1-@2)",RooArgList(gamma2,invm2S,q02S));
        RooFormulaVar signum2("signum2","(TMath::Sign(-1.,@0-@1)+1)/2.",RooArgList(invm2S,q02S));
        RooFormulaVar signum2_2("signum2_2","(TMath::Sign(-1.,@4*pow(@0,@1)*exp(@2)+@3)+1)/2.",RooArgList(a2,alpha2,b2,c2,signum2));

        if(useAnalyticalBKG) background2S = new RooGenericPdf("background2S","signum2*pow(a2,alpha2)*exp(b2)",RooArgSet(signum2,a2,alpha2,b2));

//declare fit variables

    RooRealVar m_chib1P= RooRealVar("Mean #chi_{b}(1P)","Mean #chi_{b}(1P)",9.888,9.85,9.925);
    RooRealVar m_chib2P= RooRealVar("Mean #chi_{b}(2P)","Mean #chi_{b}(2P)",10.25,10.225,10.275);

    RooRealVar sigma  = RooRealVar("#sigma","#sigma",0.01,0.005,0.02);
    RooRealVar alpha  = RooRealVar("#alpha","#alpha",0.654,0.,5.);
    RooRealVar n      = RooRealVar("n","n",2.58,.5,15.);
    RooRealVar sigma1  = RooRealVar("#sigma1","#sigma1",0.01,0.005,0.02);
    RooRealVar sigma2  = RooRealVar("#sigma2","#sigma2",0.01,0.005,0.02);
    RooRealVar sigma3  = RooRealVar("#sigma3","#sigma3",1.2094e-02);



    double m_chib3P_START=10.51;
    double m_chib3P_1S_START=10.51;
    double m_chib3P_3S_START=10.51;

    RooRealVar m_chib3P= RooRealVar("M_{#chi_{b}(3P)}","Mean #chi_{b}(3P)",m_chib3P_START,10.325,10.56);
    RooRealVar m_chib3P_1S= RooRealVar("M_{#chi_{b}(3P),1S}","Mean #chi_{b}(3P) in 1S sample",m_chib3P_1S_START,10.325,10.56);
    RooRealVar m_chib3P_3S= RooRealVar("M_{#chi_{b}(3P),3S}","Mean #chi_{b}(3P) in 3S sample",m_chib3P_3S_START,10.325,10.56);

    m_chib3P.setVal(m_chib3P_START);
    m_chib3P_1S.setVal(m_chib3P_1S_START);
    m_chib3P_3S.setVal(m_chib3P_3S_START);



/*
	RooCBShape chib1P_sig        = RooCBShape ("chib1P_sig","chib1P_sig",invm1S,m_chib1P,sigma,alpha,n);
	RooRealVar chib1P_nevt       = RooRealVar ("N#chi_{b}(1P)","N#chi_{b}(1P)",100,0,500);

    RooCBShape chib2P_sig      = RooCBShape ("chib2P_sig","chib2P_sig",invm1S,m_chib2P,sigma,alpha,n);
    RooRealVar chib2P_nevt     = RooRealVar ("N#chi_{b}(2P)","N#chi_{b}(2P)",50,0,500);

    RooCBShape chib3P_sig      = RooCBShape ("chib3P_sig","chib3P_sig",invm2S,m_chib3P,sigma,alpha,n);
    RooRealVar chib3P_nevt     = RooRealVar ("N#chi_{b}(3P)","N#chi_{b}(3P)",20,0,50);
*/

/////////// TwoCB adventure //////////////////

    char MassScaleFormula[200];
    char SigmaScaleFormula[200];

    double sigmaInd_START=0.010;
    double sigmaInd2P_START=0.010;
    double sigmaInd2P_2S_START=0.010;
    double sigmaInd3P_START=0.010;
    double sigmaInd3P_1S_START=0.010;
    double alpha_3J_START=0.6;
    double n_3J_START=3.;
    double PhotonMassScale_START=0.9840;
    double PhotonMassScale2P_START=0.9835;
    double PhotonMassScale3P_START=0.985;
    double PhotonSigmaScale_START=0.0171;
    double PhotonSigmaScale2P_START=0.0170;
    double PhotonSigmaScale3P_START=0.016561;

//    PhotonMassScale3P_START+=0;
//    PhotonMassScale3P_START+=0.0005;
//    PhotonSigmaScale3P_START+=0.;
//    PhotonSigmaScale3P_START+=-0.000498;

    RooRealVar sigmaInd  = RooRealVar("sigmaInd","sigmaInd",sigmaInd_START,0.003,0.015);
    RooRealVar sigmaInd2P  = RooRealVar("sigmaInd2P","sigmaInd2P",sigmaInd2P_START,0.003,0.025);
    RooRealVar sigmaInd2P_2S  = RooRealVar("sigmaInd2P_2S","sigmaInd2P_2S",sigmaInd2P_2S_START,0.003,0.025);
    RooRealVar sigmaInd3P  = RooRealVar("sigmaInd3P","sigmaInd3P",sigmaInd3P_START,0.003,0.025);
    RooRealVar sigmaInd3P_1S  = RooRealVar("sigmaInd3P_1S","sigmaInd3P_1S",sigmaInd3P_1S_START,0.003,0.04);
    RooRealVar alpha_3J  = RooRealVar("#alpha","#alpha_3J",alpha_3J_START,0.,3.);//0.2,1.);
    RooRealVar alpha_3J_2P  = RooRealVar("#alpha_{1P}","#alpha_3J",alpha_3J_START,0.,3.);
    RooRealVar alpha_3J_3P  = RooRealVar("#alpha_{3P}","#alpha_3J",alpha_3J_START,0.,3.);
    RooRealVar n_3J      = RooRealVar("n","n_3J",n_3J_START,.5,15.);
    RooRealVar n_3J_2P      = RooRealVar("n_3J_2P","n_3J_2P",n_3J_START,.5,15.);
    RooRealVar PhotonMassScale= RooRealVar("PES_{1P#rightarrow 1S+#gamma}","PhotonMassScale",PhotonMassScale_START,0.95,1.05);//0.97,0.99);
    RooRealVar PhotonMassScale2P= RooRealVar("PES_{2P#rightarrow 1S+#gamma}","PhotonMassScale2P",PhotonMassScale2P_START,0.95,1.05);
    RooRealVar PhotonMassScale3P= RooRealVar("PES_{3P#rightarrow nS+#gamma}","PhotonMassScale3P",PhotonMassScale3P_START,0.95,1.05);
    RooRealVar PhotonMassScaleX= RooRealVar("PES_{X}","PhotonMassScaleX",PhotonMassScale3P_START,0.95,1.05);
    RooRealVar PhotonSigmaScale= RooRealVar("#sigma_{Q}/Q_{1P#rightarrow 1S+#gamma}","PhotonSigmaScale",PhotonSigmaScale_START,0.005,0.04);//0.012,0.02);
    RooRealVar PhotonSigmaScale2P= RooRealVar("#sigma_{Q}/Q_{2P#rightarrow 1S+#gamma}","PhotonSigmaScale",PhotonSigmaScale2P_START,0.005,0.04);
    RooRealVar PhotonSigmaScale3P= RooRealVar("#sigma_{Q}/Q_{3P#rightarrow nS+#gamma}","PhotonSigmaScale",PhotonSigmaScale3P_START,0.01,0.025);
    RooRealVar PhotonSigmaScaleX= RooRealVar("#sigma_{Q}/Q_{X}","PhotonSigmaScale",PhotonSigmaScale3P_START,0.01,0.05);
    RooRealVar PhotonMassScale2P_2Sin1S= RooRealVar("PES_{2P#rightarrow 1S+X+#gamma}","PhotonMassScale2P_2Sin1S",PhotonMassScale3P_START,0.95,1.05);
    RooRealVar PhotonMassScale2P_2S= RooRealVar("PES_{2P#rightarrow 2S+#gamma}","PhotonMassScale2P_1S",PhotonMassScale3P_START,0.95,1.05);
    RooRealVar PhotonSigmaScale2P_2Sin1S= RooRealVar("#sigma_{Q}/Q_{2P#rightarrow 1S+X+#gamma}","PhotonSigmaScale",PhotonSigmaScale3P_START,0.01,0.05);
    RooRealVar PhotonSigmaScale2P_2S= RooRealVar("#sigma_{Q}/Q_{2P#rightarrow 2S+#gamma}","PhotonSigmaScale",PhotonSigmaScale3P_START,0.01,0.05);

    sigmaInd.setVal(sigmaInd_START);
    sigmaInd2P.setVal(sigmaInd2P_START);
    sigmaInd2P_2S.setVal(sigmaInd2P_2S_START);
    sigmaInd3P.setVal(sigmaInd3P_START);
    sigmaInd3P_1S.setVal(sigmaInd3P_1S_START);
    alpha_3J.setVal(alpha_3J_START);
    alpha_3J_3P.setVal(alpha_3J_START);
    n_3J.setVal(n_3J_START);
    n_3J_2P.setVal(n_3J_START);
    PhotonMassScale.setVal(PhotonMassScale_START);
    PhotonMassScale2P.setVal(PhotonMassScale2P_START);
    PhotonMassScale3P.setVal(PhotonMassScale3P_START);
    PhotonSigmaScale.setVal(PhotonSigmaScale_START);
    PhotonSigmaScale2P.setVal(PhotonSigmaScale2P_START);
    PhotonSigmaScale3P.setVal(PhotonSigmaScale3P_START);
    PhotonMassScale2P_2Sin1S.setVal(PhotonSigmaScale3P_START);
    PhotonSigmaScale2P_2Sin1S.setVal(PhotonSigmaScale3P_START);
    PhotonMassScale2P_2S.setVal(PhotonSigmaScale3P_START);
    PhotonSigmaScale2P_2S.setVal(PhotonSigmaScale3P_START);

    if(FixStuffForStep2){
        PhotonMassScale3P.setVal(0.9862);
        PhotonSigmaScale3P.setVal(0.015346);
        PhotonMassScale3P.setConstant();
        PhotonSigmaScale3P.setConstant();
    }

    if(!MCsample) n_3J.setConstant();
//    alpha_3J.setConstant();

    double ratio_J2overJ1_START=.5;
    RooRealVar ratio_J2overJ1= RooRealVar("N_{#chi_{b2}}/N_{#chi_{b1}}","ratio_J2overJ1",ratio_J2overJ1_START,0.,5.);
    ratio_J2overJ1.setVal(ratio_J2overJ1_START);
	RooFormulaVar fracJ1("fracJ1", "1./(1.+@0)", RooArgList(ratio_J2overJ1));

    ratio_J2overJ1.setConstant();

    RooRealVar ratio_J2overJ1_1P= RooRealVar("N_{#chi_{b2}}/N_{#chi_{b1}},1P","ratio_J2overJ1",ratio_J2overJ1_START,0.,5.);
    ratio_J2overJ1_1P.setVal(ratio_J2overJ1_START);
	RooFormulaVar fracJ1_1P("fracJ1_1P", "1./(1.+@0)", RooArgList(ratio_J2overJ1_1P));



	RooRealVar fracJ1_2P= RooRealVar("fracJ1_2P","fracJ1_2P",0.666);

    RooRealVar m_chib1P1fix= RooRealVar("m_chib1P1fix","m_chib1P1fix",M_PDG_chib1_1P);
    RooRealVar m_chib1P2fix= RooRealVar("m_chib1P2fix","m_chib1P2fix",M_PDG_chib2_1P);//mean~9.9, Q=0.44

    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass1S,Ymass1S);
    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass1S);

    RooFormulaVar m_chib1P1float("m_chib1P1float", MassScaleFormula, RooArgList(m_chib1P1fix, PhotonMassScale));
	RooFormulaVar m_chib1P2float("m_chib1P2float", MassScaleFormula, RooArgList(m_chib1P2fix, PhotonMassScale));

    RooFormulaVar w_chib1P1float("w_chib1P1float", SigmaScaleFormula, RooArgList(m_chib1P1float, PhotonSigmaScale));
    RooFormulaVar w_chib1P2float("w_chib1P2float", SigmaScaleFormula, RooArgList(m_chib1P2float, PhotonSigmaScale));


    //	RooCBShape chib1P1      = RooCBShape ("chib1P1","chib1P1",invm1S,m_chib1P1float,sigmaInd,alpha_3J,n_3J);
//    RooCBShape chib1P2      = RooCBShape ("chib1P2","chib1P2",invm1S,m_chib1P2float,sigmaInd,alpha_3J,n_3J);
	RooCBShape chib1P1      = RooCBShape ("chib1P1","chib1P1",invm1S,m_chib1P1float,w_chib1P1float,alpha_3J,n_3J);
    RooCBShape chib1P2      = RooCBShape ("chib1P2","chib1P2",invm1S,m_chib1P2float,w_chib1P2float,alpha_3J,n_3J);

    double n1PnJ_Max=6000;
    if(MCsample) n1PnJ_Max=6000;
    double n1PnJ_START=400;
    RooRealVar n1PnJ= RooRealVar("N_{1P#rightarrow 1S+#gamma}","n1PJ",n1PnJ_START,100,n1PnJ_Max);
    n1PnJ.setVal(n1PnJ_START);
    RooAddPdf chib1P_sigInd= RooAddPdf("chib1P_sigInd","chib1P_sigInd",RooArgList(chib1P1,chib1P2),RooArgList(fracJ1_1P));

    RooFormulaVar chib1P_nevt("chib1P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n1PnJ,fracJ1_1P));
    RooAddPdf chib1P_sig= RooAddPdf("chib1P_sig","chib1P_sig",RooArgList(chib1P_sigInd),RooArgList(chib1P_nevt));



    RooRealVar m_chib2P1fix= RooRealVar("m_chib2P1fix","m_chib2P1fix",M_PDG_chib1_2P);
    RooRealVar m_chib2P2fix= RooRealVar("m_chib2P2fix","m_chib2P2fix",M_PDG_chib2_2P);//mean~10.260, Q=0.8

	RooFormulaVar m_chib2P1float("m_chib2P1float", MassScaleFormula, RooArgList(m_chib2P1fix, PhotonMassScale2P));
	RooFormulaVar m_chib2P2float("m_chib2P2float", MassScaleFormula, RooArgList(m_chib2P2fix, PhotonMassScale2P));

    RooFormulaVar w_chib2P1float("w_chib2P1float", SigmaScaleFormula, RooArgList(m_chib2P1float, PhotonSigmaScale2P));
    RooFormulaVar w_chib2P2float("w_chib2P2float", SigmaScaleFormula, RooArgList(m_chib2P2float, PhotonSigmaScale2P));

//    RooCBShape chib2P1      = RooCBShape ("chib2P1","chib2P1",invm1S,m_chib2P1float,sigmaInd2P,alpha_3J,n_3J);
//    RooCBShape chib2P2      = RooCBShape ("chib2P2","chib2P2",invm1S,m_chib2P2float,sigmaInd2P,alpha_3J,n_3J);
    RooCBShape chib2P1      = RooCBShape ("chib2P1","chib2P1",invm1S,m_chib2P1float,w_chib2P1float,alpha_3J,n_3J);
    RooCBShape chib2P2      = RooCBShape ("chib2P2","chib2P2",invm1S,m_chib2P2float,w_chib2P2float,alpha_3J,n_3J);

    double n2PnJ_Max=3000;
    if(MCsample) n2PnJ_Max=15000;
    double n2PnJ_START=260;
    RooRealVar n2PnJ= RooRealVar("N_{2P#rightarrow 1S+#gamma}","n2PJ",n2PnJ_START,100,n2PnJ_Max);
    n2PnJ.setVal(n2PnJ_START);
    RooAddPdf chib2P_sigInd= RooAddPdf("chib2P_sigInd","chib2P_sigInd",RooArgList(chib2P1,chib2P2),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt("chib2P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n2PnJ,fracJ1));
    RooAddPdf chib2P_sig= RooAddPdf("chib2P_sig","chib2P_sig",RooArgList(chib2P_sigInd),RooArgList(chib2P_nevt));


    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass2S,Ymass2S);
    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass2S);

	RooFormulaVar m_chib2P1float_2S("m_chib2P1float_2S", MassScaleFormula, RooArgList(m_chib2P1fix, PhotonMassScale2P_2S));
	RooFormulaVar m_chib2P2float_2S("m_chib2P2float_2S", MassScaleFormula, RooArgList(m_chib2P2fix, PhotonMassScale2P_2S));

    RooFormulaVar w_chib2P1float_2S("w_chib2P1float_2S", SigmaScaleFormula, RooArgList(m_chib2P1float_2S, PhotonSigmaScale2P_2S));
    RooFormulaVar w_chib2P2float_2S("w_chib2P2float_2S", SigmaScaleFormula, RooArgList(m_chib2P2float_2S, PhotonSigmaScale2P_2S));

//	RooCBShape chib2P1_2S      = RooCBShape ("chib2P1_2S","chib2P1_2S",invm2S,m_chib2P1float_2S,sigmaInd2P_2S,alpha_3J,n_3J);
//    RooCBShape chib2P2_2S      = RooCBShape ("chib2P2_2S","chib2P2_2S",invm2S,m_chib2P2float_2S,sigmaInd2P_2S,alpha_3J,n_3J);
	RooCBShape chib2P1_2S      = RooCBShape ("chib2P1_2S","chib2P1_2S",invm2S,m_chib2P1float_2S,w_chib2P1float_2S,alpha_3J,n_3J);
    RooCBShape chib2P2_2S      = RooCBShape ("chib2P2_2S","chib2P2_2S",invm2S,m_chib2P2float_2S,w_chib2P2float_2S,alpha_3J,n_3J);

    double n2PnJ_2S_START=8;
    RooRealVar n2PnJ_2S= RooRealVar("N_{2P#rightarrow 2S+#gamma}","n2PnJ_2S",n2PnJ_2S_START,0,300);
    n2PnJ_2S.setVal(n2PnJ_2S_START);
    RooAddPdf chib2P_sigInd_2S= RooAddPdf("chib2P_sigInd_2S","chib2P_sigInd_2S",RooArgList(chib2P1_2S,chib2P2_2S),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt_2S("chib2P_nevt_2S", "@1*@0+(1-@1)*@0", RooArgList(n2PnJ_2S,fracJ1));
    RooAddPdf chib2P_sig_2S= RooAddPdf("chib2P_sig_2S","chib2P_sig_2S",RooArgList(chib2P_sigInd_2S),RooArgList(chib2P_nevt_2S));


    sprintf(MassScaleFormula,"(@0-%f-0.012*(1-@2))*@1+%f",Ymass2S,Ymass2S);
	RooFormulaVar m_chib3P1float("m_chib3P1float", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale, fracJ1));
    sprintf(MassScaleFormula,"(@0-%f+0.012*@2)*@1+%f",Ymass2S,Ymass2S);
	RooFormulaVar m_chib3P2float("m_chib3P2float", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale, fracJ1));

    RooFormulaVar w_chib3P1float("w_chib3P1float", SigmaScaleFormula, RooArgList(m_chib3P1float, PhotonSigmaScale));
    RooFormulaVar w_chib3P2float("w_chib3P2float", SigmaScaleFormula, RooArgList(m_chib3P2float, PhotonSigmaScale));

//    RooCBShape chib3P1      = RooCBShape ("chib3P1","chib3P1",invm2S,m_chib3P1float,sigmaInd3P,alpha_3J,n_3J);
//    RooCBShape chib3P2      = RooCBShape ("chib3P2","chib3P2",invm2S,m_chib3P2float,sigmaInd3P,alpha_3J,n_3J);
    RooCBShape chib3P1      = RooCBShape ("chib3P1","chib3P1",invm2S,m_chib3P1float,w_chib3P1float,alpha_3J,n_3J);
    RooCBShape chib3P2      = RooCBShape ("chib3P2","chib3P2",invm2S,m_chib3P2float,w_chib3P1float,alpha_3J,n_3J);

    double n3PnJ_START=8;
    RooRealVar n3PnJ= RooRealVar("N_{3P#rightarrow 2S+#gamma},2S","n3PnJ",n3PnJ_START,0,300);
    n3PnJ.setVal(n3PnJ_START);
    RooAddPdf chib3P_sigInd= RooAddPdf("chib3P_sigInd","chib3P_sigInd",RooArgList(chib3P1,chib3P2),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt("chib3P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ,fracJ1));
    RooAddPdf chib3P_sig= RooAddPdf("chib3P_sig","chib3P_sig",RooArgList(chib3P_sigInd),RooArgList(chib3P_nevt));





    sprintf(MassScaleFormula,"(@0-%f-0.012*(1-@2))*@1+%f",Ymass1S,Ymass1S);
	RooFormulaVar m_chib3P1float_1S("m_chib3P1float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale3P, fracJ1));
    sprintf(MassScaleFormula,"(@0-%f+0.012*@2)*@1+%f",Ymass1S,Ymass1S);
	RooFormulaVar m_chib3P2float_1S("m_chib3P2float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale3P, fracJ1));

    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass1S);
    RooFormulaVar w_chib3P1float_1S("w_chib3P1float_1S", SigmaScaleFormula, RooArgList(m_chib3P1float_1S, PhotonSigmaScale3P));
    RooFormulaVar w_chib3P2float_1S("w_chib3P2float_1S", SigmaScaleFormula, RooArgList(m_chib3P2float_1S, PhotonSigmaScale3P));

//    RooCBShape chib3P1_1S      = RooCBShape ("chib3P1_1S","chib3P1_1S",invm1S,m_chib3P1float_1S,sigmaInd3P_1S,alpha_3J,n_3J);
//    RooCBShape chib3P2_1S     = RooCBShape ("chib3P2_1S","chib3P2_1S",invm1S,m_chib3P2float_1S,sigmaInd3P_1S,alpha_3J,n_3J);
    RooCBShape chib3P1_1S      = RooCBShape ("chib3P1_1S","chib3P1_1S",invm1S,m_chib3P1float_1S,w_chib3P1float_1S,alpha_3J,n_3J);
    RooCBShape chib3P2_1S     = RooCBShape ("chib3P2_1S","chib3P2_1S",invm1S,m_chib3P2float_1S,w_chib3P2float_1S,alpha_3J,n_3J);

    double n3PnJ_1S_START=100;
    RooRealVar n3PnJ_1S= RooRealVar("N_{3P#rightarrow 1S+#gamma}","n3PnJ_1S",n3PnJ_1S_START,0,1000);
    n3PnJ_1S.setVal(n3PnJ_1S_START);
    RooAddPdf chib3P_sigInd_1S= RooAddPdf("chib3P_sigInd_1S","chib3P_sigInd_1S",RooArgList(chib3P1_1S,chib3P2_1S),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt_1S("chib3P_nevt_1S", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ_1S,fracJ1));
    RooAddPdf chib3P_sig_1S= RooAddPdf("chib3P_sig_1S","chib3P_sig_1S",RooArgList(chib3P_sigInd_1S),RooArgList(chib3P_nevt_1S));




    sprintf(MassScaleFormula,"(@0-%f-0.012*(1-@2))*@1+%f",Ymass3S,Ymass3S);
	RooFormulaVar m_chib3P1float_3S("m_chib3P1float_3S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale3P, fracJ1));
    sprintf(MassScaleFormula,"(@0-%f+0.012*@2)*@1+%f",Ymass3S,Ymass3S);
	RooFormulaVar m_chib3P2float_3S("m_chib3P2float_3S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale3P, fracJ1));

    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass3S);
    RooFormulaVar w_chib3P1float_3S("w_chib3P1float_3S", SigmaScaleFormula, RooArgList(m_chib3P1float_3S, PhotonSigmaScale3P));
    RooFormulaVar w_chib3P2float_3S("w_chib3P2float_3S", SigmaScaleFormula, RooArgList(m_chib3P2float_3S, PhotonSigmaScale3P));

//    RooCBShape chib3P1_3S      = RooCBShape ("chib3P1_3S","chib3P1_3S",invm3S,m_chib3P1float_3S,sigmaInd3P_3S,alpha_3J,n_3J);
//    RooCBShape chib3P2_3S     = RooCBShape ("chib3P2_3S","chib3P2_3S",invm3S,m_chib3P2float_3S,sigmaInd3P_3S,alpha_3J,n_3J);
    RooCBShape chib3P1_3S      = RooCBShape ("chib3P1_3S","chib3P1_3S",invm3S,m_chib3P1float_3S,w_chib3P1float_3S,alpha_3J,n_3J);
    RooCBShape chib3P2_3S     = RooCBShape ("chib3P2_3S","chib3P2_3S",invm3S,m_chib3P2float_3S,w_chib3P2float_3S,alpha_3J,n_3J);

    double n3PnJ_3S_START=100;
    RooRealVar n3PnJ_3S= RooRealVar("N_{3P#rightarrow 3S+#gamma}","n3PnJ_3S",n3PnJ_3S_START,0,1000);
    n3PnJ_3S.setVal(n3PnJ_3S_START);
    RooAddPdf chib3P_sigInd_3S= RooAddPdf("chib3P_sigInd_3S","chib3P_sigInd_3S",RooArgList(chib3P1_3S,chib3P2_3S),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt_3S("chib3P_nevt_3S", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ_3S,fracJ1));
    RooAddPdf chib3P_sig_3S= RooAddPdf("chib3P_sig_3S","chib3P_sig_3S",RooArgList(chib3P_sigInd_3S),RooArgList(chib3P_nevt_3S));




    RooRealVar mX= RooRealVar("mX","mX",9.69,9.66,9.72);

/*    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass1S,Ymass1S);
    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass1S);
	RooFormulaVar m_Xfloat("m_Xfloat", MassScaleFormula, RooArgList(mX, PhotonMassScale3P));
    RooFormulaVar w_Xfloat("w_Xfloat", SigmaScaleFormula, RooArgList(m_Xfloat, PhotonSigmaScale3P));
    RooRealVar nX= RooRealVar("N_{X}","nX",10,0,50);
    RooCBShape XCB      = RooCBShape ("XCB","XCB",invm1S,m_Xfloat,w_Xfloat,alpha_3J,n_3J);
    RooAddPdf Xsig= RooAddPdf("Xsig","Xsig",RooArgList(XCB),RooArgList(nX));
*/




    double Mass2P1_2Sin1S=m_chib2P1fix.getVal()+Ymass1S-Ymass2S;
    double Mass2P2_2Sin1S=m_chib2P2fix.getVal()+Ymass1S-Ymass2S;

    RooRealVar m_chib2P1fix_2Sin1S= RooRealVar("m_chib2P1fix_2Sin1S","m_chib2P1fix_2Sin1S",Mass2P1_2Sin1S);
    RooRealVar m_chib2P2fix_2Sin1S= RooRealVar("m_chib2P2fix_2Sin1S","m_chib2P2fix_2Sin1S",Mass2P2_2Sin1S);//mean~9.9, Q=0.44

    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass1S,Ymass1S);
    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass1S);

	RooFormulaVar m_chib2P1float_2Sin1S("m_chib2P1float_2Sin1S", MassScaleFormula, RooArgList(m_chib2P1fix_2Sin1S, PhotonMassScale2P_2Sin1S));
	RooFormulaVar m_chib2P2float_2Sin1S("m_chib2P2float_2Sin1S", MassScaleFormula, RooArgList(m_chib2P2fix_2Sin1S, PhotonMassScale2P_2Sin1S));

//	RooFormulaVar m_chib2P1float_2Sin1S("m_chib2P1float_2Sin1S", "@0", RooArgList(mX));
//	RooFormulaVar m_chib2P2float_2Sin1S("m_chib2P2float_2Sin1S", "@0+0.01319", RooArgList(mX));

    RooFormulaVar w_chib2P1float_2Sin1S("w_chib2P1float_2Sin1S", SigmaScaleFormula, RooArgList(m_chib2P1float_2Sin1S, PhotonSigmaScale2P_2Sin1S));
    RooFormulaVar w_chib2P2float_2Sin1S("w_chib2P2float_2Sin1S", SigmaScaleFormula, RooArgList(m_chib2P1float_2Sin1S, PhotonSigmaScale2P_2Sin1S));

	RooCBShape chib2P1_2Sin1S      = RooCBShape ("chib2P1_2Sin1S","chib2P1_2Sin1S",invm1S,m_chib2P1float_2Sin1S,w_chib2P1float_2Sin1S,alpha_3J,n_3J);
    RooCBShape chib2P2_2Sin1S      = RooCBShape ("chib2P2_2Sin1S","chib2P2_2Sin1S",invm1S,m_chib2P2float_2Sin1S,w_chib2P2float_2Sin1S,alpha_3J,n_3J);

    double nX_START=10;
    RooRealVar nX= RooRealVar("N_{2P#rightarrow 1S+X+#gamma}","nX",nX_START,0,200);
    nX.setVal(nX_START);
    RooAddPdf chib2P_sigInd_2Sin1S= RooAddPdf("chib2P_sigInd_2Sin1S","chib2P_sigInd_2Sin1S",RooArgList(chib2P1_2Sin1S,chib2P2_2Sin1S),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt_2Sin1S("chib2P_nevt_2Sin1S", "@1*@0+(1-@1)*@0", RooArgList(nX,fracJ1));
    RooAddPdf Xsig= RooAddPdf("Xsig","Xsig",RooArgList(chib2P_sigInd_2Sin1S),RooArgList(chib2P_nevt_2Sin1S));


////////////////////////////////////////////////////////


    RooRealVar background_nevt1S = RooRealVar("N_{bkg,1S+#gamma}","N_{bkg}1S",1500,0,100000);
    RooRealVar background_nevtSB1S = RooRealVar("background_nevtSB1S","background_nevtSB1S",1500,0,10000);
    RooAddPdf backgroundSB1S= RooAddPdf("backgroundSB1S","backgroundSB1S",RooArgList(*background1S),RooArgList(background_nevtSB1S));
    RooAddPdf backgroundE1S= RooAddPdf("backgroundE1S","backgroundE1S",RooArgList(*background1S),RooArgList(background_nevt1S));
    RooRealVar background_nevt2S = RooRealVar("N_{bkg,2S+#gamma}","N_{bkg}2S",1500,0,40000);
    RooRealVar background_nevtSB2S = RooRealVar("background_nevtSB2S","background_nevtSB2S",1500,0,10000);
    RooAddPdf backgroundSB2S= RooAddPdf("backgroundSB2S","backgroundSB2S",RooArgList(*background2S),RooArgList(background_nevtSB2S));
    RooAddPdf backgroundE2S= RooAddPdf("backgroundE2S","backgroundE2S",RooArgList(*background2S),RooArgList(background_nevt2S));

    RooRealVar background_nevt3S = RooRealVar("N_{bkg,3S+#gamma}","N_{bkg}3S",1500,0,100000);
    RooRealVar background_nevtSB3S = RooRealVar("background_nevtSB3S","background_nevtSB3S",1500,0,10000);
    RooAddPdf backgroundSB3S= RooAddPdf("backgroundSB3S","backgroundSB3S",RooArgList(*background3S),RooArgList(background_nevtSB3S));


    sprintf(DScutChar,"invm1S > %f && invm1S < %f",invm_min1S,invm_max1S);
    RooDataSet* data1S_masswindow=(RooDataSet*)data1S->reduce(Cut(DScutChar));
    sprintf(DScutChar,"invm2S > %f && invm2S < %f",invm_min2S,invm_max2S);
    RooDataSet* data2S_masswindow=(RooDataSet*)data2S->reduce(Cut(DScutChar));
    sprintf(DScutChar,"invm3S > %f && invm3S < %f",invm_min3S,invm_max3S);
    RooDataSet* data3S_masswindow=(RooDataSet*)data3S->reduce(Cut(DScutChar));




////////// Y1S //////////////////////

	 RooArgSet Set_modelPdf1S_1P2P("Set_modelPdf1S_1P2P");
	 Set_modelPdf1S_1P2P.add(chib1P_sig);
	 Set_modelPdf1S_1P2P.add(chib2P_sig);
	 if(!MCsample) Set_modelPdf1S_1P2P.add(*background1S);

	 RooArgSet Set_Coef_modelPdf1S_1P2P("Set_Coef_modelPdf1S_1P2P");
	 Set_Coef_modelPdf1S_1P2P.add(chib1P_nevt);
	 Set_Coef_modelPdf1S_1P2P.add(chib2P_nevt);
	 if(!MCsample) Set_Coef_modelPdf1S_1P2P.add(background_nevt1S);


//    RooAddPdf modelPdf1S= RooAddPdf("modelPdf1S","modelPdf1S",RooArgList(chib1P_sig,chib2P_sig,chib3P_sig_1S,*background1S),RooArgList(chib1P_nevt,chib2P_nevt,chib3P_nevt_1S,background_nevt1S));
//    RooAddPdf modelPdf1S_1P2P= RooAddPdf("modelPdf1S_1P2P","modelPdf1S_1P2P",Set_modelPdf1S_1P2P,Set_Coef_modelPdf1S_1P2P);

// MC:
//    RooAddPdf modelPdf1S_1P2P= RooAddPdf("modelPdf1S_1P2P","modelPdf1S_1P2P",RooArgList(chib1P_sig,chib2P_sig),RooArgList(chib1P_nevt,chib2P_nevt));

// Include X:
    RooAddPdf modelPdf1S= RooAddPdf("modelPdf1S","modelPdf1S",RooArgList(Xsig,chib1P_sig,chib2P_sig,chib3P_sig_1S,*background1S),RooArgList(nX,chib1P_nevt,chib2P_nevt,chib3P_nevt_1S,background_nevt1S));
    RooAddPdf modelPdf1S_1P2P= RooAddPdf("modelPdf1S_1P2P","modelPdf1S_1P2P",RooArgList(Xsig,chib1P_sig,chib2P_sig,*background1S),RooArgList(nX,chib1P_nevt,chib2P_nevt,background_nevt1S));

    RooAddPdf modelPdf1S_3P= RooAddPdf("modelPdf1S_3P","modelPdf1S_3P",RooArgList(chib3P_sig_1S,*background1S),RooArgList(chib3P_nevt_1S,background_nevt1S));

////////// Y2S //////////////////////

    RooAddPdf modelPdf2S= RooAddPdf("modelPdf2S","modelPdf2S",RooArgList(chib3P_sig,chib2P_sig_2S,*background2S),RooArgList(chib3P_nevt,chib2P_nevt_2S,background_nevt2S));
    RooAddPdf modelPdf2S_red= RooAddPdf("modelPdf2S_red","modelPdf2S_red",RooArgList(chib3P_sig,chib2P_sig_2S,*background2S),RooArgList(chib3P_nevt,chib2P_nevt_2S,background_nevt2S));
    RooAddPdf modelPdf2S_2P= RooAddPdf("modelPdf2S_2P","modelPdf2S_2P",RooArgList(chib2P_sig_2S,*background2S),RooArgList(chib2P_nevt_2S,background_nevt2S));

////////// Y3S //////////////////////

    RooAddPdf modelPdf3S= RooAddPdf("modelPdf3S","modelPdf3S",RooArgList(chib3P_sig_3S,*background3S),RooArgList(chib3P_nevt_3S,background_nevt3S));
    RooAddPdf modelPdf3S_red= RooAddPdf("modelPdf3S_red","modelPdf3S_red",RooArgList(chib3P_sig_3S,*background3S),RooArgList(chib3P_nevt_3S,background_nevt3S));


    invm1S.setRange("SBregion2",chib1Pmax,chib2Pmin);
    invm1S.setRange("SBregion31S",bb_thresh,invm_max1S);

    invm1S.setRange("Chib1Pregion",invm_min1S,chib1Pmax);
    invm1S.setRange("Chib2Pregion",chib2Pmin,chib2Pmax);
    invm1S.setRange("Chib3Pregion",chib2Pmax,bb_thresh);
    invm1S.setRange("Above2P",chib2Pmax,invm_max1S);

    invm1S.setRange("ChibSBregion1",invm_min1S,chib1Pmin);
    invm1S.setRange("ChibSBregion2",chib1Pmax,chib2Pmin);
    invm1S.setRange("ChibSBregion3",bb_thresh,invm_max1S);

    invm1S.setRange("TotalBlindedRegion1",invm_min1S,chib2Pmax);
    invm1S.setRange("TotalBlindedRegion2",bb_thresh,invm_max1S);

    invm1S.setRange("Chib1P2Pregion",invm_min1S,chib2Pmax);

    invm2S.setRange("Y2SAll",invm_min2S,invm_max2S);
    invm1S.setRange("Y1SAll",invm_min1S,invm_max1S);
    invm2S.setRange("SBregion32S",bb_thresh,invm_max2S);
    invm2S.setRange("Chib3Pregion2S",invm_min2S,bb_thresh);

    char SBregion2Char[200];
    sprintf(SBregion2Char,"invm1S > %f && invm1S < %f",chib1Pmax,chib2Pmin);
    char SBregion31SChar[200];
    sprintf(SBregion31SChar,"invm1S > %f && invm1S < %f",bb_thresh,invm_max1S);
    char SBregion32SChar[200];
    sprintf(SBregion32SChar,"invm2S > %f && invm2S < %f",bb_thresh,invm_max2S);
    char SBregion3SChar[200];
    sprintf(SBregion3SChar,"invm3S > %f && invm3S < %f",bb_thresh,invm_max3S);

    double Nbkg_1S_SB1=data1S->sumEntries(SBregion2Char);
    double Nbkg_1S_SB2=data1S->sumEntries(SBregion31SChar);
    double Nbkg_2S_SB1=data2S->sumEntries(SBregion32SChar);
    double Nbkg_3S_SB1=data3S->sumEntries(SBregion3SChar);
    cout<<"background_nevt1S sumentriesrange = " << Nbkg_1S_SB1<<endl;
    cout<<"background_nevt1S sumentriesrange2 = " << Nbkg_1S_SB2<<endl;
    cout<<"background_nevt2S sumentriesrange = " << Nbkg_2S_SB1<<endl;
    cout<<"background_nevt3S sumentriesrange = " << Nbkg_3S_SB1<<endl;


    double BkgSBInt1S;
    BkgSBInt1S = NormalizedIntegral(&backgroundSB1S, invm1S, bb_thresh, invm_max1S);
    cout<<BkgSBInt1S<<endl;
    background_nevt1S.setVal((Nbkg_1S_SB2)/BkgSBInt1S);
    background_nevt1S.setConstant();
    cout<<"background normalization 1S = " << background_nevt1S.getVal()<<endl;

    double BkgSBInt2S;
    BkgSBInt2S = NormalizedIntegral(&backgroundSB2S, invm2S, bb_thresh, invm_max2S);
    cout<<BkgSBInt2S<<endl;
    background_nevt2S.setVal(Nbkg_2S_SB1/BkgSBInt2S);
    background_nevt2S.setConstant();
    cout<<"background normalization 2S = " << background_nevt2S.getVal()<<endl;

    double BkgSBInt3S;
    BkgSBInt3S = NormalizedIntegral(&backgroundSB3S, invm3S, bb_thresh, invm_max3S);
    cout<<BkgSBInt3S<<endl;
    background_nevt3S.setVal(Nbkg_3S_SB1/BkgSBInt3S);
    background_nevt3S.setConstant();
    cout<<"background normalization 3S = " << background_nevt3S.getVal()<<endl;

    double covQual_BG2S=0;
    double covQual_1P2P2S=0;
    double covQual_BG1S=0;
    double covQual_1P2P1S=0;

    double Nevt_buff2;
    double Nevt_buff=n3PnJ_1S.getVal();
    n3PnJ_1S.setVal(0);
    n3PnJ_1S.setConstant();
    m_chib3P.setConstant();
    sigmaInd3P_1S.setConstant();

    char SignormRegion[200];
    sprintf(SignormRegion,"Chib1P2Pregion,ChibSBregion3");
    char OnePTwoPDatacut[200];
    sprintf(OnePTwoPDatacut,"invm1S<%f",chib2Pmax);

    char BGnormRegion[200];
     sprintf(BGnormRegion,"ChibSBregion1,ChibSBregion2,ChibSBregion3");

     if(useAnalyticalBKG) background_nevt1S.setConstant(kFALSE);
     if(useAnalyticalBKG) background_nevt2S.setConstant(kFALSE);

     if(MCsample){
     	background_nevt1S.setVal(0.001);
     	background_nevt1S.setConstant();
//     	n2PnJ.setVal(0.001);
//     	n2PnJ.setConstant();
     	n3PnJ_1S.setVal(0.001);
     	n3PnJ_1S.setConstant();
//     	PhotonMassScale2P.setConstant();
//     	PhotonSigmaScale2P.setConstant();
     	alpha1.setConstant();
     	beta1.setConstant();
     	q01S.setConstant();
     	ratio_J2overJ1.setConstant(kFALSE);
        m_chib3P.setConstant();
     }

/*    RooFitResult* BG1SFitResult = background1S.fitTo(*data1S,Save(1),Range(BGnormRegion),PrintEvalErrors(-1),PrintLevel(-1),Warnings(kFALSE));
    BG1SFitResult->Print();

    alpha1.setConstant();
    beta1.setConstant();
    q0.setConstant();

*/

 /*       RooFitResult* Fit1P2P = modelPdf1S_1P2P.fitTo(*data1S,Save(1),Range(SignormRegion),PrintEvalErrors(-1),PrintLevel(1),Warnings(kFALSE));
        Fit1P2P->Print();
        covQual_1P2P1S = Fit1P2P->covQual();
*/

 //   alpha1.setConstant();
  //  beta1.setConstant();


	 //	 RooFitResult* r_sig = modelPdf1S_1P2P.fitTo(*data1S,Save(kTRUE),Range("Y1SAll"),PrintEvalErrors(-1),PrintLevel(-1),Warnings(kFALSE));
	 //	 r_sig->Print();

	 data1S->Print();
	 RooDataSet* data1Sexclude3P=(RooDataSet*)data1S->reduce(Cut("invm1S<1000"));//->reduce(Cut("invm1S<10.325||invm1S>10.58"));
	 data1Sexclude3P->Print();

	 RooAbsReal* nll1S_1P2P_=modelPdf1S_1P2P.createNLL(*data1S,Range(invm_min1S,invm_max1S));




	 RooAbsReal* nll1S_1P2P=modelPdf1S_1P2P.createNLL(*data1S,Range(invm_min1S,chib2Pmax));//chib2Pmax
	 RooAbsReal* nll1S_BB=modelPdf1S_1P2P.createNLL(*data1S,Range(bb_thresh,invm_max1S));//bb_thresh

//	 RooAbsReal* nll1S_1P2P=modelPdf1S_1P2P.createNLL(*data1S,Range(SignormRegion),SplitRange(kTRUE));//chib2Pmax


	 RooArgSet simNLL_Set("simNLL_Set");
	 simNLL_Set.add(*nll1S_1P2P);
	 if(useAnalyticalBKG) simNLL_Set.add(*nll1S_BB);

     RooAddition simNLL_innitial = RooAddition("add_innitial","add_innitial",simNLL_Set);//RooArgSet(*nll1S_1P2P));

	 RooMinuit* minuit1 = new RooMinuit(simNLL_innitial);
	 if(useAnalyticalBKG) minuit1 = new RooMinuit(*nll1S_1P2P_);

	 minuit1->setStrategy(2);
	 minuit1->setPrintEvalErrors(-1);
	 minuit1->setPrintLevel(-1);
	 minuit1->setNoWarn();

	 char result1name[200];

	 cout<<"InitialHesse_"<<endl;
	 minuit1->hesse();
	 sprintf(result1name,"initialhesse");
	 RooFitResult* initialhesseresult1=minuit1->save(result1name,result1name);
	 cout<<"INITIAL HESSE RESULT:"<<endl;
	 cout<<initialhesseresult1->edm()<<endl;
	 cout<<initialhesseresult1->status()<<endl;
	 cout<<"InitialHesse-covQual"<<endl;
	 cout<<initialhesseresult1->covQual()<<endl;
	 initialhesseresult1->Print();
       delete initialhesseresult1;

	 cout<<"Migrad_"<<endl;
	 minuit1->migrad();
	 sprintf(result1name,"migrad");
	 RooFitResult* migradresult1=minuit1->save(result1name,result1name);
	 cout<<"MIGRAD RESULT:"<<endl;
	 cout<<migradresult1->edm()<<endl;
	 cout<<migradresult1->status()<<endl;
	 cout<<"Migrad-covQual"<<endl;
	 cout<<migradresult1->covQual()<<endl;
	 migradresult1->Print();

	 if(migradresult1->edm()>0.1){
	 cout<<"Migrad_2"<<endl;
	 minuit1->migrad();
	 sprintf(result1name,"migrad");
	 migradresult1=minuit1->save(result1name,result1name);
	 cout<<"MIGRAD RESULT 2:"<<endl;
	 cout<<migradresult1->edm()<<endl;
	 cout<<migradresult1->status()<<endl;
	 cout<<"Migrad-covQual 2"<<endl;
	 cout<<migradresult1->covQual()<<endl;
	 migradresult1->Print();
	 }
	 int covQualOptimization=migradresult1->covQual();
     delete migradresult1;


	 cout<<"Hesse_"<<endl;
	 minuit1->hesse();
	 sprintf(result1name,"hesse");
	 RooFitResult* hesseresult1=minuit1->save(result1name,result1name);
	 cout<<"HESSE RESULT:"<<endl;
	 cout<<hesseresult1->edm()<<endl;
	 cout<<hesseresult1->status()<<endl;
	 cout<<"Hesse-covQual"<<endl;
	 cout<<hesseresult1->covQual()<<endl;
	 hesseresult1->Print();
       delete hesseresult1;

/*        RooFitResult* Fit1P2P = modelPdf1S_1P2P.fitTo(*data1S,Save(1),Range("Chib1P2Pregion,ChibSBregion3"),PrintEvalErrors(-1),PrintLevel(1),Warnings(kFALSE));
        Fit1P2P->Print();
        covQual_1P2P1S = Fit1P2P->covQual();
*/
        n3PnJ_1S.setConstant(kFALSE);
        n3PnJ_1S.setVal(Nevt_buff);
        m_chib3P.setConstant(kFALSE);
        sigmaInd3P_1S.setConstant(kFALSE);

        sigma.setConstant();
        alpha.setConstant();
        n.setConstant();

        sigmaInd.setConstant();
        sigmaInd2P.setConstant();
        alpha_3J.setConstant();
        n_3J.setConstant();
        PhotonMassScale2P.setConstant();
        ratio_J2overJ1.setConstant();
        ratio_J2overJ1_1P.setConstant();
        PhotonSigmaScale.setConstant();
        PhotonMassScale.setConstant();
        PhotonMassScale2P_2Sin1S.setConstant();
        PhotonSigmaScale2P_2Sin1S.setConstant();

        alpha1.setConstant();
        beta1.setConstant();
        background_nevt1S.setConstant();
        q01S.setConstant();
        gamma1.setConstant();

        mX.setConstant();
        nX.setConstant();

        if(MCsample){
            n3PnJ_1S.setVal(0.01);
            n3PnJ_1S.setConstant();
            m_chib3P.setConstant();
//            n2PnJ.setVal(0.01);
//            n2PnJ.setConstant();

        }

   	 RooAbsReal* nll2S_2P=modelPdf2S_2P.createNLL(*data2S,Range(invm_min2S,10.55));
   	 RooAbsReal* nll2S_BB=modelPdf2S_2P.createNLL(*data2S,Range(10.55,invm_max2S));

        RooAddition simNLL_innitial2S = RooAddition("simNLL_innitial2S","simNLL_innitial2S",RooArgSet(*nll2S_2P,*nll2S_BB));

/*   	 RooMinuit* minuit2 = new RooMinuit(simNLL_innitial2S);

   	 minuit2->setStrategy(1);
   	 minuit2->setPrintEvalErrors(-1);
   	 minuit2->setPrintLevel(-1);
   	 minuit2->setNoWarn();

   	 char result2name[200];

   	 cout<<"InitialHesse_"<<endl;
   	 minuit2->hesse();
   	 sprintf(result2name,"initialhesse");
   	 RooFitResult* initialhesseresult2=minuit2->save(result2name,result2name);
   	 cout<<"INITIAL HESSE RESULT:"<<endl;
   	 cout<<initialhesseresult2->edm()<<endl;
   	 cout<<initialhesseresult2->status()<<endl;
   	 cout<<"InitialHesse-covQual"<<endl;
   	 cout<<initialhesseresult2->covQual()<<endl;
   	 initialhesseresult2->Print();
          delete initialhesseresult2;

   	 cout<<"Migrad_"<<endl;
   	 minuit2->migrad();
   	 sprintf(result2name,"migrad");
   	 RooFitResult* migradresult2=minuit2->save(result2name,result2name);
   	 cout<<"MIGRAD RESULT:"<<endl;
   	 cout<<migradresult2->edm()<<endl;
   	 cout<<migradresult2->status()<<endl;
   	 cout<<"Migrad-covQual"<<endl;
   	 cout<<migradresult2->covQual()<<endl;
   	 migradresult2->Print();
          delete migradresult2;


   	 cout<<"Hesse_"<<endl;
   	 minuit2->hesse();
   	 sprintf(result2name,"hesse");
   	 RooFitResult* hesseresult2=minuit2->save(result2name,result2name);
   	 cout<<"HESSE RESULT:"<<endl;
   	 cout<<hesseresult2->edm()<<endl;
   	 cout<<hesseresult2->status()<<endl;
   	 cout<<"Hesse-covQual"<<endl;
   	 cout<<hesseresult2->covQual()<<endl;
   	 hesseresult2->Print();
          delete hesseresult2;
*/
        alpha2.setConstant();
        beta2.setConstant();
        gamma2.setConstant();
        q02S.setConstant();

//        PhotonMassScale.setVal(PhotonMassScale.getVal()-PhotonMassScale.getError());


        background_nevt1S.setConstant();
        background_nevt2S.setConstant();
        background_nevt3S.setConstant();

//        delete Fit1P2P;
//        sigmaInd2P.setVal(sigmaInd.getVal());

        double k=(sigmaInd.getVal()-sigmaInd2P.getVal())/(430.-783.);
        double d=sigmaInd.getVal()-k*430.;
//        double k=((sigmaInd.getVal()-sigmaInd.getError())-(sigmaInd2P.getVal()+sigmaInd2P.getError()))/(430.-783.);
//        double d=(sigmaInd.getVal()-sigmaInd.getError())-k*430.;
        sigmaInd2P_2S.setVal(k*220.+d);
        sigmaInd2P_2S.setConstant();
        sigmaInd3P.setVal(k*487+d);
        sigmaInd3P.setConstant();
        sigmaInd3P_1S.setVal(k*1050.+d);
        sigmaInd3P_1S.setConstant();

/*        cout<<"1P width = "<<sigmaInd.getVal()<<endl;
        cout<<"2P width = "<<sigmaInd2P.getVal()<<endl;
        cout<<"k = "<<k<<endl;
        cout<<"d = "<<d<<endl;
        cout<<"3P width (1S) = "<<sigmaInd3P_1S.getVal()<<endl;
        cout<<"2P width (2S) = "<<sigmaInd2P_2S.getVal()<<endl;
        cout<<"3P width (2S) = "<<sigmaInd3P.getVal()<<endl;
*/
        cout<<"1P width = "<<sigmaInd.getVal()<<endl;
        cout<<"2P width = "<<sigmaInd2P.getVal()<<endl;
        cout<<"k = "<<k<<endl;
        cout<<"d = "<<d<<endl;
        cout<<"3P width (1S) = "<<sigmaInd3P_1S.getVal()<<endl;
        cout<<"2P width (2S) = "<<sigmaInd2P_2S.getVal()<<endl;
        cout<<"3P width (2S) = "<<sigmaInd3P.getVal()<<endl;

cout<<"deltaM 1P, J2: "<<m_chib1P2float.getVal()-m_chib1P2fix.getVal()<<endl;
cout<<"deltaM 1P, J1: "<<m_chib1P1float.getVal()-m_chib1P1fix.getVal()<<endl;
cout<<"deltaM 2P, J2: "<<m_chib2P2float.getVal()-m_chib2P2fix.getVal()<<endl;
cout<<"deltaM 2P, J1: "<<m_chib2P1float.getVal()-m_chib2P1fix.getVal()<<endl;
cout<<"deltaM 2P, 2S sample, J2: "<<m_chib2P2float_2S.getVal()-m_chib2P2fix.getVal()<<endl;
cout<<"deltaM 2P, 2S sample, J1: "<<m_chib2P1float_2S.getVal()-m_chib2P1fix.getVal()<<endl;

m_chib1P.setVal(m_chib1P1float.getVal()*fracJ1.getVal()+m_chib1P2float.getVal()*(1-fracJ1.getVal()));
m_chib2P.setVal(m_chib2P1float.getVal()*fracJ1.getVal()+m_chib2P2float.getVal()*(1-fracJ1.getVal()));
cout<<"m_chib1P: "<<m_chib1P.getVal()<<endl;
cout<<"m_chib2P: "<<m_chib2P.getVal()<<endl;


if(Optimize||determine3PparametersDuringStep1){
	cout<<"Set 3P PES and PWS to the mean of the 1P/2P estimates"<<endl;
PhotonMassScale3P.setVal( (PhotonMassScale.getVal()+PhotonMassScale2P.getVal())/2.);
PhotonSigmaScale3P.setVal( (PhotonSigmaScale.getVal()+PhotonSigmaScale2P.getVal())/2.);
}



		 RooAbsReal* nll1S=modelPdf1S_3P.createNLL(*data1S,Range(chib2Pmax,bb_thresh));
		 RooAbsReal* nll2S=modelPdf2S_red.createNLL(*data2S,Range(invm_min2S,bb_thresh));
		 RooAbsReal* nll3S=modelPdf3S_red.createNLL(*data3S,Range(invm_min3S,bb_thresh));

		 RooArgSet simNLL_Set2("simNLL_Set2");
		 simNLL_Set2.add(*nll1S);
		 if(nState>1) simNLL_Set2.add(*nll2S);
		 if(nState>2) simNLL_Set2.add(*nll3S);


		 RooAddition simNLL = RooAddition("add","add",simNLL_Set2);
//         RooAddition simNLL = RooAddition("add","add",RooArgSet(*nll1S));
         double Nllbefore3P=simNLL.getVal();

		 RooMinuit* minuit = new RooMinuit(simNLL);

		 minuit->setStrategy(1);
		 minuit->setPrintEvalErrors(-1);
		 minuit->setPrintLevel(-1);
		 minuit->setNoWarn();

		 char resultname[200];

		 cout<<"InitialHesse_"<<endl;
		 minuit->hesse();
		 sprintf(resultname,"initialhesse");
		 RooFitResult* initialhesseresult=minuit->save(resultname,resultname);
		 cout<<"INITIAL HESSE RESULT:"<<endl;
		 cout<<initialhesseresult->edm()<<endl;
		 cout<<initialhesseresult->status()<<endl;
		 cout<<"InitialHesse-covQual"<<endl;
		 cout<<initialhesseresult->covQual()<<endl;
		 initialhesseresult->Print();
	        delete initialhesseresult;

		 cout<<"Migrad_"<<endl;
		 minuit->migrad();
		 sprintf(resultname,"migrad");
		 RooFitResult* migradresult=minuit->save(resultname,resultname);
		 cout<<"MIGRAD RESULT:"<<endl;
		 cout<<migradresult->edm()<<endl;
		 cout<<migradresult->status()<<endl;
		 cout<<"Migrad-covQual"<<endl;
		 cout<<migradresult->covQual()<<endl;
		 migradresult->Print();
	        delete migradresult;


		 cout<<"Hesse_"<<endl;
		 minuit->hesse();
		 sprintf(resultname,"hesse");
		 RooFitResult* hesseresult=minuit->save(resultname,resultname);
		 cout<<"HESSE RESULT:"<<endl;
		 cout<<hesseresult->edm()<<endl;
		 cout<<hesseresult->status()<<endl;
		 cout<<"Hesse-covQual"<<endl;
		 cout<<hesseresult->covQual()<<endl;
		 hesseresult->Print();
	        delete hesseresult;


			 cout<<"Minos_"<<endl;
			 minuit->minos(RooArgSet(m_chib3P));
			 sprintf(resultname,"minos");
			 RooFitResult* minosresult=minuit->save(resultname,resultname);
			 cout<<"MINOS RESULT:"<<endl;
			 cout<<minosresult->edm()<<endl;
			 cout<<minosresult->status()<<endl;
			 cout<<"Minos-covQual"<<endl;
			 cout<<minosresult->covQual()<<endl;
			 minosresult->Print();
		        delete minosresult;



//////// Calculate significance of 3P state
/*
//		        modelProd.fitTo(data1S,Extended(kTRUE),Save());

//		        mean.setConstant(kTRUE);
//		        width.setConstant(kTRUE);
		        m_chib3P.setConstant(kTRUE);
		        n3PnJ_1S.setConstant(kTRUE);

		        RooArgSet POI(RooArgList(n3PnJ_1S,m_chib3P));

		        ProfileLikelihoodCalculator plc(data1S, modelPdf1S_3P, POI);

		        // Create a copy of the POI parameters to set the values to zero
		        RooArgSet nullparams;
		        nullparams.addClone(POI);
		        ((RooRealVar *) (nullparams.find("N_{#chi_{b}(3P),1S}")))->setVal(0);
		        ((RooRealVar *) (nullparams.find("M_{#chi_{b}(3P),1S}")))->setConstant(kTRUE);
		        plc.SetNullParameters(nullparams);

		        HypoTestResult* htr = plc.GetHypoTest();

		        cout << "significance: " << htr->NullPValue() << " " << htr->Significance() << endl;


*/

		        m_chib3P.setError((m_chib3P.getAsymErrorHi()-m_chib3P.getAsymErrorLo())/2.);
//		        m_chib3P_1S.setError((m_chib3P_1S.getAsymErrorHi()-m_chib3P_1S.getAsymErrorLo())/2.);

		        background_nevt2S.setConstant();


		 double Nevt_buff3;
         double Nllafter3P=simNLL.getVal();
         Nevt_buff=n3PnJ.getVal();
         Nevt_buff2=n3PnJ_1S.getVal();
         Nevt_buff3=n3PnJ_3S.getVal();
 		 n3PnJ.setVal(0.001);
 		 n3PnJ_1S.setVal(0.001);
 		 n3PnJ_3S.setVal(0.001);
 		 Nllbefore3P=simNLL.getVal();
 		 n3PnJ.setVal(Nevt_buff);
 		 n3PnJ_1S.setVal(Nevt_buff2);
 		 n3PnJ_3S.setVal(Nevt_buff3);

         cout<<"minNllNoSigAll "<<Nllbefore3P<<endl;
         cout<<"minNll3PAll    "<<Nllafter3P<<endl;
         double logLikelihoodRatio3PAll=TMath::Log(TMath::Exp(Nllafter3P-Nllbefore3P));
         cout<<"logLikelihoodRatio3PAll "<<logLikelihoodRatio3PAll<<endl;
         double ThreePsigALL=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3PAll));
         cout<<"ThreePsigALL "<<ThreePsigALL<<endl;




        RooAbsReal* model2SNLL=modelPdf2S.createNLL(*data2S);
        Nevt_buff=n3PnJ.getVal();
		double minNll3P=model2SNLL->getVal();
		n3PnJ.setVal(0);
		double minNllNoSig2S=model2SNLL->getVal();
		n3PnJ.setVal(Nevt_buff);

        cout<<"minNllNoSig2S "<<minNllNoSig2S<<endl;
        cout<<"minNll3P    "<<minNll3P<<endl;
        double logLikelihoodRatio3P=TMath::Log(TMath::Exp(minNll3P-minNllNoSig2S));
        cout<<"logLikelihoodRatio3P "<<logLikelihoodRatio3P<<endl;
        double ThreePsig=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3P));
        cout<<"ThreePsig "<<ThreePsig<<endl;


        RooAbsReal* model1SNLL=modelPdf1S.createNLL(*data1S);
        RooAbsReal* model1SNLL_3P=modelPdf1S_3P.createNLL(*data1S,Range(chib2Pmax,bb_thresh));
        Nevt_buff=n3PnJ_1S.getVal();
		double minNll3P_1S=model1SNLL_3P->getVal();
		n3PnJ_1S.setVal(0.001);
		double minNllNoSig1S=model1SNLL_3P->getVal();
		n3PnJ_1S.setVal(Nevt_buff);

		cout<<"minNllNoSig1S "<<minNllNoSig1S<<endl;
        cout<<"minNll3P_1S    "<<minNll3P_1S<<endl;
        double logLikelihoodRatio3P_1S=TMath::Log(TMath::Exp(minNll3P_1S-minNllNoSig1S));
        cout<<"logLikelihoodRatio3P_1S "<<logLikelihoodRatio3P_1S<<endl;
        double ThreePsig_1S=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3P_1S));
        cout<<"ThreePsig_1S "<<ThreePsig_1S<<endl;


        RooAbsReal* model3SNLL=modelPdf3S.createNLL(*data3S);
        Nevt_buff=n3PnJ_3S.getVal();
		double minNll3P_3S=model3SNLL->getVal();
		n3PnJ_3S.setVal(0.001);
		double minNllNoSig3S=model3SNLL->getVal();
		n3PnJ_3S.setVal(Nevt_buff);

		cout<<"minNllNoSig3S "<<minNllNoSig3S<<endl;
        cout<<"minNll3P_3S    "<<minNll3P_3S<<endl;
        double logLikelihoodRatio3P_3S=TMath::Log(TMath::Exp(minNll3P_3S-minNllNoSig3S));
        cout<<"logLikelihoodRatio3P_3S "<<logLikelihoodRatio3P_3S<<endl;
        double ThreePsig_3S=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3P_3S));
        cout<<"ThreePsig_3S "<<ThreePsig_3S<<endl;



        Nevt_buff=n1PnJ.getVal();
		double minNll1P=model1SNLL->getVal();
		n1PnJ.setVal(0);
		double minNllNo1P1S=model1SNLL->getVal();
		n1PnJ.setVal(Nevt_buff);

		cout<<"minNllNo1P1S "<<minNllNo1P1S<<endl;
        cout<<"minNll1P    "<<minNll1P<<endl;
        double logLikelihoodRatio1P=TMath::Log(TMath::Exp(minNll1P-minNllNo1P1S));
        cout<<"logLikelihoodRatio1P "<<logLikelihoodRatio1P<<endl;
        double OnePsig=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio1P));
        cout<<"OnePsig "<<OnePsig<<endl;

        Nevt_buff=n2PnJ.getVal();
		double minNll2P=model1SNLL->getVal();
		n2PnJ.setVal(0);
		double minNllNo2P1S=model1SNLL->getVal();
		n2PnJ.setVal(Nevt_buff);

		cout<<"minNllNo2P1S "<<minNllNo2P1S<<endl;
        cout<<"minNll2P    "<<minNll2P<<endl;
        double logLikelihoodRatio2P=TMath::Log(TMath::Exp(minNll2P-minNllNo2P1S));
        cout<<"logLikelihoodRatio2P "<<logLikelihoodRatio2P<<endl;
        double TwoPsig=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio2P));
        cout<<"TwoPsig "<<TwoPsig<<endl;



        double Nevt_buff_;
        Nevt_buff=n1PnJ.getVal();
        Nevt_buff_=n2PnJ.getVal();
		double minNll1P2Pcomb=model1SNLL->getVal();
		n2PnJ.setVal(0);
		n1PnJ.setVal(0);
		double minNllNo1P2Pcomv=model1SNLL->getVal();
		n1PnJ.setVal(Nevt_buff);
		n2PnJ.setVal(Nevt_buff_);

		cout<<"minNllNo1P2Pcomv "<<minNllNo1P2Pcomv<<endl;
        cout<<"minNll1P2Pcomb    "<<minNll1P2Pcomb<<endl;
        double logLikelihoodRatio1P2Pcomb=TMath::Log(TMath::Exp(minNll1P2Pcomb-minNllNo1P2Pcomv));
        cout<<"logLikelihoodRatio1P2Pcomb "<<logLikelihoodRatio1P2Pcomb<<endl;
        double OnePTwoPsigComb=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio1P2Pcomb));
        cout<<"OnePTwoPsigComb "<<OnePTwoPsigComb<<endl;


        Nevt_buff=n2PnJ_2S.getVal();
		double minNll2P_2S=model2SNLL->getVal();
		n2PnJ_2S.setVal(0);
		double minNllNo2P_2S=model2SNLL->getVal();
		n2PnJ_2S.setVal(Nevt_buff);

		cout<<"minNllNo2P_2S "<<minNllNo2P_2S<<endl;
        cout<<"minNll2P_2S    "<<minNll2P_2S<<endl;
        double logLikelihoodRatio2P_2S=TMath::Log(TMath::Exp(minNll2P_2S-minNllNo2P_2S));
        cout<<"logLikelihoodRatio2P_2S "<<logLikelihoodRatio2P_2S<<endl;
        double TwoPsig_2S=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio2P_2S));
        cout<<"TwoPsig_2S "<<TwoPsig_2S<<endl;





        Nevt_buff=nX.getVal();
		double minNll_X=model1SNLL->getVal();
		nX.setVal(0);
		double minNll_NoX=model1SNLL->getVal();
		nX.setVal(Nevt_buff);

		cout<<"minNll_NoX "<<minNll_NoX<<endl;
        cout<<"minNll_X    "<<minNll_X<<endl;
        double logLikelihoodRatio_X=TMath::Log(TMath::Exp(minNll_X-minNll_NoX));
        cout<<"logLikelihoodRatio_X "<<logLikelihoodRatio_X<<endl;
        double Xsignificance=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio_X));
        cout<<"Xsignificance "<<Xsignificance<<endl;

		delete model2SNLL;
		delete model1SNLL;
		delete model3SNLL;

        cout<<"chib1P_nevt: "<<chib1P_nevt.getVal()<<endl;
        cout<<"chib2P_nevt: "<<chib2P_nevt.getVal()<<endl;
        cout<<"chib3P_nevt_1S: "<<chib3P_nevt_1S.getVal()<<endl;
        cout<<"chib3P_nevt_3S: "<<chib3P_nevt_3S.getVal()<<endl;
        cout<<"chib3P_nevt: "<<chib3P_nevt.getVal()<<endl;
        cout<<"chib2P_nevt_2S: "<<chib2P_nevt_2S.getVal()<<endl;


     double TotalRealEvents1S = data1S_masswindow->sumEntries();
     double totalEventsInFit1S;
     double TotalRealEvents2S = data2S_masswindow->sumEntries();
     double totalEventsInFit2S;
     double TotalRealEvents3S = data3S_masswindow->sumEntries();
     double totalEventsInFit3S;

     totalEventsInFit1S = nX.getVal()+chib1P_nevt.getVal()+chib2P_nevt.getVal()+chib3P_nevt_1S.getVal()+background_nevt1S.getVal();
     totalEventsInFit2S = chib3P_nevt.getVal()+chib2P_nevt_2S.getVal()+background_nevt2S.getVal();
     totalEventsInFit3S = chib3P_nevt_3S.getVal()+background_nevt3S.getVal();

     cout<< "N_tot_Fit1S =                            "<<totalEventsInFit1S <<endl;
     cout<< "N_tot_Samlpe1S =                         "<<TotalRealEvents1S <<endl;
     cout<< "Delta_N_tot1S = N_tot_Sample1S-N_tot_Fit1S = "<<TotalRealEvents1S-totalEventsInFit1S <<endl;

     cout<< "N_tot_Fit2S =                            "<<totalEventsInFit2S <<endl;
     cout<< "N_tot_Samlpe2S =                         "<<TotalRealEvents2S <<endl;
     cout<< "Delta_N_tot2S = N_tot_Sample2S-N_tot_Fit2S = "<<TotalRealEvents2S-totalEventsInFit2S <<endl;

     cout<< "N_tot_Fit3S =                            "<<totalEventsInFit3S <<endl;
     cout<< "N_tot_Samlpe3S =                         "<<TotalRealEvents3S <<endl;
     cout<< "Delta_N_tot3S = N_tot_Sample3S-N_tot_Fit3S = "<<TotalRealEvents3S-totalEventsInFit3S <<endl;

     double onePwidth=(m_chib1P1fix.getVal()-Ymass1S)*PhotonSigmaScale.getVal();
     double onePwidtherr=(m_chib1P1fix.getVal()-Ymass1S)*PhotonSigmaScale.getError();
     cout<< "onePwidth =                            "<<onePwidth <<endl;
     cout<< "onePwidtherr =                            "<<onePwidtherr <<endl;
// likelihood gaussian check

     if(SaveAll) {

    	 char SaveLogs[200];

     Nevt_buff=n3PnJ.getVal();

     int nEvt3P2S=Nevt_buff;
     TH1F  *logL_3P2S = new TH1F("logL_3P2S","",2*nEvt3P2S,0,2*nEvt3P2S);
     double Nllafter3P2Svar;
     for(int i=0;i<2*nEvt3P2S+1;i++){
     n3PnJ.setVal(i);
	 Nllafter3P2Svar=simNLL.getVal();
	 cout<<"Nsig3P_2S = "<<i<<" Nll = "<<Nllafter3P2Svar<<endl;
	 logL_3P2S->SetBinContent(i,-Nllafter3P2Svar);
     }
     sprintf(SaveLogs,"Figures/%s/logL_3P2S.root",FitID);
     logL_3P2S->SaveAs(SaveLogs);

     n3PnJ.setVal(Nevt_buff);


     Nevt_buff=n3PnJ_1S.getVal();

     int nEvt3P1S=Nevt_buff;
     TH1F  *logL_3P1S = new TH1F("logL_3P1S","",2*nEvt3P1S,0,2*nEvt3P1S);
     double Nllafter3P1Svar;
     for(int i=0;i<2*nEvt3P1S+1;i++){
     n3PnJ_1S.setVal(i);
	 Nllafter3P1Svar=simNLL.getVal();
	 cout<<"Nsig3P_1S = "<<i<<" Nll = "<<Nllafter3P1Svar<<endl;
	 logL_3P1S->SetBinContent(i,-Nllafter3P1Svar);
     }
     sprintf(SaveLogs,"Figures/%s/logL_3P1S.root",FitID);
     logL_3P1S->SaveAs(SaveLogs);

     n3PnJ_1S.setVal(Nevt_buff);

     delete logL_3P2S;
     delete logL_3P1S;
     }

     double nSig=1.5;

     double fmchib1P = m_chib1P.getVal();
     double fmchib2P = m_chib2P.getVal();
     double fmchib3P = m_chib3P.getVal();
     double err_fmchib3P = m_chib3P.getError();
     double fmchib3P_1S = m_chib3P.getVal();//m_chib3P_1S.getVal();
     double err_fmchib3P_1S = m_chib3P.getError();//m_chib3P_1S.getError();
     double fmchib3P_3S = m_chib3P.getVal();//m_chib3P_1S.getVal();
     double err_fmchib3P_3S = m_chib3P.getError();//m_chib3P_1S.getError();

     double fmchib1_3P=fmchib3P-0.012*(1-1./(1+ratio_J2overJ1.getVal()));
     double fmchib2_3P=fmchib3P+0.012*1./(1+ratio_J2overJ1.getVal());


     double OneSigmaCL = 0.682689492137;

     double sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib1P_sig, invm1S, fmchib1P-sigmaCalc, fmchib1P+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma1=sigmaCalc;

     sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib2P_sig, invm1S, fmchib2P-sigmaCalc, fmchib2P+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma2=sigmaCalc;

     sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib3P_sig, invm2S, (fmchib3P-Ymass2S)*PhotonMassScale3P.getVal()+Ymass2S-sigmaCalc, (fmchib3P-Ymass2S)*PhotonMassScale3P.getVal()+Ymass2S+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma3=sigmaCalc;

     sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib3P_sig_1S, invm1S, (fmchib3P_1S-Ymass1S)*PhotonMassScale3P.getVal()+Ymass1S-sigmaCalc, (fmchib3P_1S-Ymass1S)*PhotonMassScale3P.getVal()+Ymass1S+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma3_1S=sigmaCalc;

     sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib3P_sig_3S, invm3S, (fmchib3P_3S-Ymass3S)*PhotonMassScale3P.getVal()+Ymass3S-sigmaCalc, (fmchib3P_3S-Ymass3S)*PhotonMassScale3P.getVal()+Ymass3S+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma3_3S=sigmaCalc;



     m_chib2P.setVal(m_chib2P1float.getVal()*fracJ1.getVal()+m_chib2P2float.getVal()*(1-fracJ1.getVal()));


     cout<<"m_chib1P: "<<m_chib1P.getVal()<<endl;
     cout<<"m_chib2P: "<<m_chib2P.getVal()<<endl;
     cout<<"m_chib3P: "<<m_chib3P.getVal()<<endl;
     cout<<"m_chib3P_1S: "<<m_chib3P_1S.getVal()<<endl;
     cout<<"m_chib3P_3S: "<<m_chib3P_3S.getVal()<<endl;
     cout<<"fsigma1: "<<fsigma1<<endl;
     cout<<"fsigma2: "<<fsigma2<<endl;
     cout<<"fsigma3: "<<fsigma3<<endl;
     cout<<"fsigma3_1S: "<<fsigma3_1S<<endl<<endl;
     cout<<"fsigma3_3S: "<<fsigma3_3S<<endl<<endl;

     cout<<"m_chib1_2Pin1S: "<<m_chib2P1float_2Sin1S.getVal()<<endl;
     cout<<"m_chib2_2Pin1S: "<<m_chib2P2float_2Sin1S.getVal()<<endl;

     double m_chib1_2Pin1S_corrected=(m_chib2P1float_2Sin1S.getVal()-Ymass1S)/PhotonMassScale3P.getVal()+Ymass1S;
     cout<<"m_chib1_2Pin1S corrected: "<<m_chib1_2Pin1S_corrected<<endl;
     double expectedMass_m_chib1_2Pin1S_corrected=m_chib2P1fix.getVal()-Ymass2S+Ymass1S;
     cout<<"expectedMass_m_chib1_2Pin1S_corrected: "<<expectedMass_m_chib1_2Pin1S_corrected<<endl;
     cout<<"shift m_chib1_2Pin1S = "<<m_chib1_2Pin1S_corrected-expectedMass_m_chib1_2Pin1S_corrected<<" +- "<<mX.getError()<<endl;

     double m_chib2_2Pin1S_corrected=(m_chib2P2float_2Sin1S.getVal()-Ymass1S)/PhotonMassScale3P.getVal()+Ymass1S;
     cout<<"m_chib2_2Pin1S corrected: "<<m_chib2_2Pin1S_corrected<<endl;
     double expectedMass_m_chib2_2Pin1S_corrected=m_chib2P2fix.getVal()-Ymass2S+Ymass1S;
     cout<<"expectedMass_m_chib2_2Pin1S_corrected: "<<expectedMass_m_chib2_2Pin1S_corrected<<endl;
     cout<<"shift m_chib2_2Pin1S = "<<m_chib2_2Pin1S_corrected-expectedMass_m_chib2_2Pin1S_corrected<<" +- "<<mX.getError()<<endl;

//     cout<<"shift m_chib1_2Pin1S = "<<mX.getVal()-(m_chib2P1fix.getVal()-Ymass2S)*PhotonMassScale3P.getVal()-Ymass1S<<" +- "<<mX.getError()<<endl;
//     cout<<"shift m_chib2_2Pin1S = "<<mX.getVal()-(m_chib2P2fix.getVal()-Ymass2S)*PhotonMassScale3P.getVal()-Ymass1S<<" +- "<<mX.getError()<<endl;

     double fnchib1P = chib1P_nevt.getVal();
     double fnchib2P = chib2P_nevt.getVal();
     double fnchib3P = chib3P_nevt.getVal();
     double fnchib3P_1S = chib3P_nevt_1S.getVal();
     double fnchib3P_3S = chib3P_nevt_3S.getVal();
     double fnbg1S = background_nevt1S.getVal();
     double fnbg2S = background_nevt2S.getVal();
     double fnbg3S = background_nevt3S.getVal();
     double rchib1P,rchib2P,rchib3P,rbgchib1P,rbgchib2P,rbgchib3P,rsigchib1P,rsigchib2P,rsigchib3P,rchib3P_1S,rchib3P_3S,rbgchib3P_1S,rbgchib3P_3S,rsigchib3P_BityukovKrasnikov,rsigchib3P_BityukovKrasnikov_1S,rsigchib3P_ScL,rsigchib3P_ScL_1S;
     double rratchib1P,rratchib2P,rratchib3P,rratchib3P_1S,rsigchib3P_1S,rsigchib3P_3S,rratchib3P_3S;
     rchib1P=fnchib1P*NormalizedIntegral(&chib1P_sig, invm1S, fmchib1P-nSig*fsigma1, fmchib1P+nSig*fsigma1);
     rchib2P=fnchib2P*NormalizedIntegral(&chib2P_sig, invm1S, fmchib2P-nSig*fsigma2, fmchib2P+nSig*fsigma2);
     rchib3P=fnchib3P*NormalizedIntegral(&chib3P_sig, invm2S, fmchib3P-nSig*fsigma3, fmchib3P+nSig*fsigma3);
     rchib3P_1S=fnchib3P_1S*NormalizedIntegral(&chib3P_sig_1S, invm1S, fmchib3P_1S-nSig*fsigma3_1S, fmchib3P_1S+nSig*fsigma3_1S);
     rchib3P_3S=fnchib3P_3S*NormalizedIntegral(&chib3P_sig_3S, invm3S, fmchib3P_3S-nSig*fsigma3_3S, fmchib3P_3S+nSig*fsigma3_3S);
     rbgchib1P=fnbg1S*NormalizedIntegral(&*background1S, invm1S, fmchib1P-nSig*fsigma1, fmchib1P+nSig*fsigma1);
     rbgchib2P=fnbg1S*NormalizedIntegral(&*background1S, invm1S, fmchib2P-nSig*fsigma2, fmchib2P+nSig*fsigma2);
     rbgchib3P=fnbg2S*NormalizedIntegral(&*background2S, invm2S, fmchib3P-nSig*fsigma3, fmchib3P+nSig*fsigma3);
     rbgchib3P_1S=fnbg1S*NormalizedIntegral(&*background1S, invm1S, fmchib3P_1S-nSig*fsigma3_1S, fmchib3P_1S+nSig*fsigma3_1S);
     rbgchib3P_3S=fnbg3S*NormalizedIntegral(&*background3S, invm3S, fmchib3P_3S-nSig*fsigma3_3S, fmchib3P_3S+nSig*fsigma3_3S);
     rratchib1P=rchib1P/rbgchib1P;
     rratchib2P=rchib2P/rbgchib2P;
     rratchib3P=rchib3P/rbgchib3P;
     rratchib3P_1S=rchib3P_1S/rbgchib3P_1S;
     rratchib3P_3S=rchib3P_3S/rbgchib3P_3S;
     rsigchib1P=rchib1P/sqrt(rchib1P+rbgchib1P);
     rsigchib2P=rchib2P/sqrt(rchib2P+rbgchib2P);

     rsigchib3P=rchib3P/sqrt(rchib3P+rbgchib3P);
     rsigchib3P_1S=rchib3P_1S/sqrt(rchib3P_1S+rbgchib3P_1S);

     rsigchib3P_BityukovKrasnikov=2*(sqrt(rchib3P+rbgchib3P)-sqrt(rbgchib3P));
     rsigchib3P_BityukovKrasnikov_1S=2*(sqrt(rchib3P_1S+rbgchib3P_1S)-sqrt(rbgchib3P_1S));

     rsigchib3P_ScL=TMath::Sqrt(2*(rchib3P+rbgchib3P)*TMath::Log(1+rchib3P/rbgchib3P)-rchib3P);
     rsigchib3P_ScL_1S=sqrt(2*(rchib3P_1S+rbgchib3P_1S)*TMath::Log(1+rchib3P_1S/rbgchib3P_1S)-rchib3P_1S);

     cout<<"rchib3P: "<<rchib3P<<endl;
     cout<<"rbgchib3P: "<<rbgchib3P<<endl;
     cout<<"rratchib3P: "<<rratchib3P<<endl<<endl;
     cout<<"rchib3P_1S: "<<rchib3P_1S<<endl;
     cout<<"rbgchib3P_1S: "<<rbgchib3P_1S<<endl;
     cout<<"rratchib3P_1S: "<<rratchib3P_1S<<endl<<endl;
     cout<<"rsigchib3P_1S: "<<rsigchib3P_1S<<endl;
     cout<<"rsigchib3P: "<<rsigchib3P<<endl<<endl;
     cout<<"rsigchib3P_BityukovKrasnikov: "<<rsigchib3P_BityukovKrasnikov<<endl;
     cout<<"rsigchib3P_BityukovKrasnikov_1S: "<<rsigchib3P_BityukovKrasnikov_1S<<endl<<endl;
     cout<<"rsigchib3P_ScL: "<<rsigchib3P_ScL<<endl;
     cout<<"rsigchib3P_ScL_1S: "<<rsigchib3P_ScL_1S<<endl<<endl;

     int nHistBins_Gamma=25;
     char MasswindowChar[200];
     TH1F  *hGammaE_1P = new TH1F("hGammaE_1P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaE_2P = new TH1F("hGammaE_2P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaE_3P = new TH1F("hGammaE_3P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaE_3P_1S = new TH1F("hGammaE_3P_1S","",nHistBins_Gamma,0,5);
     TH1F  *hGammaE_1PSB2P = new TH1F("hGammaE_1PSB2P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaE_2PSB3P = new TH1F("hGammaE_2PSB3P","",nHistBins_Gamma,0,5);

     TH1F  *hGammaPt_1P = new TH1F("hGammaPt_1P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaPt_2P = new TH1F("hGammaPt_2P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaPt_3P = new TH1F("hGammaPt_3P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaPt_3P_1S = new TH1F("hGammaPt_3P_1S","",nHistBins_Gamma,0,5);
     TH1F  *hGammaPt_1PSB2P = new TH1F("hGammaPt_1PSB2P","",nHistBins_Gamma,0,5);
     TH1F  *hGammaPt_2PSB3P = new TH1F("hGammaPt_2PSB3P","",nHistBins_Gamma,0,5);

     TH1F  *hChibPt_1P = new TH1F("hChibPt_1P","",nHistBins_Gamma,0,50);
     TH1F  *hChibPt_2P = new TH1F("hChibPt_2P","",nHistBins_Gamma,0,50);

     hGammaE_1P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaE_2P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaE_3P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaE_3P_1S->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaE_1PSB2P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaE_2PSB3P->GetYaxis()->SetTitle("Events per 200 MeV");

     hGammaPt_1P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaPt_2P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaPt_3P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaPt_3P_1S->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaPt_1PSB2P->GetYaxis()->SetTitle("Events per 200 MeV");
     hGammaPt_2PSB3P->GetYaxis()->SetTitle("Events per 200 MeV");

     hChibPt_1P->GetYaxis()->SetTitle("Events per 2 GeV");
     hChibPt_2P->GetYaxis()->SetTitle("Events per 2 GeV");

     sprintf(MasswindowChar,"invm1S < %f && invm1S > %f",fmchib1P+nSig*fsigma1,fmchib1P-nSig*fsigma1);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_1P",MasswindowChar);
     tree1S->Draw("gammapt>>hGammaPt_1P",MasswindowChar);
     tree1S->Draw("sqrt((gammapx+jpsipx)*(gammapx+jpsipx)+(gammapy+jpsipy)*(gammapy+jpsipy))>>hChibPt_1P",MasswindowChar);
     sprintf(MasswindowChar,"invm1S < %f && invm1S > %f",fmchib2P+nSig*fsigma2,fmchib2P-nSig*fsigma2);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_2P",MasswindowChar);
     tree1S->Draw("gammapt>>hGammaPt_2P",MasswindowChar);
     tree1S->Draw("sqrt((gammapx+jpsipx)*(gammapx+jpsipx)+(gammapy+jpsipy)*(gammapy+jpsipy))>>hChibPt_2P",MasswindowChar);
     sprintf(MasswindowChar,"invm1S < %f && invm1S > %f",fmchib3P+nSig*fsigma3_1S,fmchib3P-nSig*fsigma3_1S);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_3P_1S",MasswindowChar);
     tree1S->Draw("gammapt>>hGammaPt_3P_1S",MasswindowChar);
     sprintf(MasswindowChar,"invm2S < %f && invm2S > %f",fmchib3P+nSig*fsigma3,fmchib3P-nSig*fsigma3);
     tree2S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_3P",MasswindowChar);
     tree2S->Draw("gammapt>>hGammaPt_3P",MasswindowChar);
     sprintf(MasswindowChar,"invm1S > %f && invm1S < %f",fmchib1P+nSig*fsigma1,fmchib2P-nSig*fsigma2);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_1PSB2P",MasswindowChar);
     tree1S->Draw("gammapt>>hGammaPt_1PSB2P",MasswindowChar);
     sprintf(MasswindowChar,"invm1S > %f && invm1S < %f",fmchib2P+nSig*fsigma2,fmchib3P-nSig*fsigma3);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_2PSB3P",MasswindowChar);
     tree1S->Draw("gammapt>>hGammaPt_2PSB3P",MasswindowChar);

     hGammaE_1P->Print();
     hGammaE_2P->Print();
     hGammaE_3P->Print();
     hGammaE_3P_1S->Print();

     char saveName[200];

     TCanvas* PhotonEnergyCanvas = new TCanvas("PhotonEnergyCanvas","PhotonEnergyCanvas",1600, 1200);
     PhotonEnergyCanvas->SetFillColor(kWhite);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_1P->GetXaxis()->SetTitle("E_{#gamma} in 1P peak (1S sample)");	hGammaE_1P->Draw("E");
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_1P_1S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_2P->GetXaxis()->SetTitle("E_{#gamma} in 2P peak (1S sample)");	hGammaE_2P->Draw("E");
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_2P_1S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_3P->GetXaxis()->SetTitle("E_{#gamma} in 3P peak (2S sample)");	hGammaE_3P->Draw("E");
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_3P_2S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_3P_1S->GetXaxis()->SetTitle("E_{#gamma} in 3P peak (1S sample)");	hGammaE_3P_1S->Draw("E");
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy3P_1S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_1PSB2P->GetXaxis()->SetTitle("E_{#gamma} of events in mass region between 1P and 2P peaks");	hGammaE_1PSB2P->Draw("E");
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy1PSB2P.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_2PSB3P->GetXaxis()->SetTitle("E_{#gamma} of events in mass region between 2P and 3P peaks");	hGammaE_2PSB3P->Draw("E");
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy2PSB3P.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);

     TCanvas* PhotonPtCanvas = new TCanvas("PhotonPtCanvas","PhotonPtCanvas",1600, 1200);
     PhotonPtCanvas->SetFillColor(kWhite);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_1P->GetXaxis()->SetTitle("p_{T#gamma} in 1P peak (1S sample)");	hGammaPt_1P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonPt_1P_1S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_2P->GetXaxis()->SetTitle("p_{T#gamma} in 2P peak (1S sample)");	hGammaPt_2P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonPt_2P_1S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_3P->GetXaxis()->SetTitle("p_{T#gamma} in 3P peak (2S sample)");	hGammaPt_3P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonPt_3P_2S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_3P_1S->GetXaxis()->SetTitle("p_{T#gamma} in 3P peak (1S sample)");	hGammaPt_3P_1S->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonPt_3P_1S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_1PSB2P->GetXaxis()->SetTitle("p_{T#gamma} of events in mass region between 1P and 2P peaks");	hGammaPt_1PSB2P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonPt_1PSB2P.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_2PSB3P->GetXaxis()->SetTitle("p_{T#gamma} of events in mass region between 2P and 3P peaks");	hGammaPt_2PSB3P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonPt_2PSB3P.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);

     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hChibPt_1P->GetXaxis()->SetTitle("p_{T#chi_{b}} in 1P peak");	hChibPt_1P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/ChibPt_1P_1S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hChibPt_2P->GetXaxis()->SetTitle("p_{T#chi_{b}} in 2P peak (1S sample)");	hChibPt_2P->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/ChibPt_2P_1S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);

     int nHistBins_Gamma_2D=50;
     TH2F  *hGammaPt_invm1S = new TH2F("hGammaPt_invm1S","",nHistBins_Gamma_2D,invm_min1S,invm_max1S,nHistBins_Gamma_2D,0,5);
     tree1S->Draw("gammapt:invm1S>>hGammaPt_invm1S");
     TH2F  *hGammaPt_invm1S_broad = new TH2F("hGammaPt_invm1S_broad","",100,invm_min1S,14,100,0,7);
     tree1S->Draw("gammapt:invm1S>>hGammaPt_invm1S_broad");
     TH2F  *hCosAlpha_invm1S = new TH2F("hCosAlpha_invm1S","",100,invm_min1S,invm_max1S,100,-1,1);
     tree1S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz):invm1S>>hCosAlpha_invm1S");
     gStyle->SetPalette(1);
     hGammaPt_invm1S->SetStats(0);
     hGammaPt_invm1S_broad->SetStats(0);
     hCosAlpha_invm1S->SetStats(0);

     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_invm1S->GetXaxis()->SetTitle(invm1S.getTitle());hGammaPt_invm1S->GetYaxis()->SetTitle("p_{T#gamma}");	hGammaPt_invm1S->Draw("colz");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/GammaPt_vs_invm1S.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_invm1S->GetXaxis()->SetTitle(invm1S.getTitle());hGammaPt_invm1S->GetYaxis()->SetTitle("p_{T#gamma}");	hGammaPt_invm1S->Draw("box");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/GammaPt_vs_invm1S_Scatter.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaPt_invm1S_broad->GetXaxis()->SetTitle(invm1S.getTitle());hGammaPt_invm1S_broad->GetYaxis()->SetTitle("p_{T#gamma}");	hGammaPt_invm1S_broad->Draw("box");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/GammaPt_vs_invm1S_Scatter_Broad.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCosAlpha_invm1S->GetXaxis()->SetTitle(invm1S.getTitle());hCosAlpha_invm1S->GetYaxis()->SetTitle("cos#alpha");	hCosAlpha_invm1S->Draw("colz");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/CosAlpha_vs_invm1S_colz.root",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);

     sprintf(saveName,"Figures/%s/GammaPt_vs_invm1S_Histo.root",FitID);
     if(SaveAll) hGammaPt_invm1S->SaveAs(saveName);

     int nHistBins_MeanPt=50;
     TH1F  *hGammaMeanPt = new TH1F("hGammaMeanPt","",nHistBins_MeanPt,invm_min1S,invm_max1S);
     double massrange=invm_max1S-invm_min1S;

     int nHistBins_GammaForMeanPt=25;
     for(int i=1;i<nHistBins_MeanPt+1;i++){
     TH1F  *hGammaPtBuffer = new TH1F("hGammaPtBuffer","",nHistBins_GammaForMeanPt,0,5);
     sprintf(MasswindowChar,"invm1S > %f && invm1S < %f",invm_min1S+(i-1)*massrange/nHistBins_MeanPt,invm_min1S+(i)*massrange/nHistBins_MeanPt);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaPtBuffer",MasswindowChar);

     double numEntriesGammaMeanPtBuffer=hGammaPtBuffer->GetSumOfWeights();
     double median;
     double binContentBuffer=0;
     for(int j=1;j<nHistBins_GammaForMeanPt+1;j++){
    	 binContentBuffer+=hGammaPtBuffer->GetBinContent(j);
    	 if(binContentBuffer/numEntriesGammaMeanPtBuffer>0.5){median=hGammaPtBuffer->GetBinCenter(j);break;}
     }
     hGammaMeanPt->SetBinContent(i,hGammaPtBuffer->GetMean());
     hGammaMeanPt->SetBinError(i,hGammaPtBuffer->GetMeanError());
//     hGammaMeanPt->SetBinContent(i,median);
//     hGammaMeanPt->SetBinError(i,0.0001);
     }

     hGammaMeanPt->SetStats(0);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaMeanPt->GetXaxis()->SetTitle(invm1S.getTitle());hGammaMeanPt->GetYaxis()->SetTitle("#hat{p_{T#gamma}}");	hGammaMeanPt->Draw("E");
     PhotonPtCanvas->Modified();
     sprintf(saveName,"Figures/%s/GammaPtMean.pdf",FitID);
     if(SaveAll) PhotonPtCanvas->SaveAs(saveName);




     delete hGammaE_1P;
     delete hGammaE_2P;
     delete hGammaE_3P;
     delete hGammaE_3P_1S;
     delete hGammaE_1PSB2P;
     delete hGammaE_2PSB3P;
     delete hGammaPt_1P;
     delete hGammaPt_2P;
     delete hGammaPt_3P;
     delete hGammaPt_3P_1S;
     delete hGammaPt_1PSB2P;
     delete hGammaPt_2PSB3P;
     delete hChibPt_1P;
     delete hChibPt_2P;





//////////////////////////////////////////////////////
///////// Calculate reflections //////////////////////
//////////////////////////////////////////////////////

// Mass values:

//     Ymass1S=9.4603;
//     Ymass2S=10.02326;
//     Ymass3S=10.3552;
// 	   M_PDG_chib1_1P=9.89278;
// 	   M_PDG_chib2_1P=9.91221;
// 	   M_PDG_chib1_2P=10.25546;
// 	   M_PDG_chib2_2P=10.26865;
//	   fmchib3P_1S
//	   fmchib1_3P
//	   fmchib2_3P

     double M_PDG_chib0_1P=9.85944;
     double M_PDG_chib0_2P=10.2325;

// Q values 'standard peaks':

     double Q_PDG_chib1_1P_1S=M_PDG_chib1_1P-Ymass1S;
     double Q_PDG_chib2_1P_1S=M_PDG_chib2_1P-Ymass1S;
     double Q_PDG_chib1_2P_1S=M_PDG_chib1_2P-Ymass1S;
     double Q_PDG_chib2_2P_1S=M_PDG_chib2_2P-Ymass1S;
     double Q_PDG_chib1_2P_2S=M_PDG_chib1_2P-Ymass2S;
     double Q_PDG_chib2_2P_2S=M_PDG_chib2_2P-Ymass2S;
     double Q_PDG_chib1_3P_1S=fmchib1_3P-Ymass1S;
     double Q_PDG_chib2_3P_1S=fmchib2_3P-Ymass1S;
     double Q_PDG_chib1_3P_2S=fmchib1_3P-Ymass2S;
     double Q_PDG_chib2_3P_2S=fmchib2_3P-Ymass2S;
     double Q_PDG_chib1_3P_3S=fmchib1_3P-Ymass3S;
     double Q_PDG_chib2_3P_3S=fmchib2_3P-Ymass3S;

     double Q_PDG_chib0_1P_1S=M_PDG_chib0_1P-Ymass1S;
     double Q_PDG_chib0_2P_1S=M_PDG_chib0_2P-Ymass1S;

     double Q_MEAS_chib1_1P_1S=Q_PDG_chib1_1P_1S*PhotonMassScale.getVal();
     double Q_MEAS_chib2_1P_1S=Q_PDG_chib2_1P_1S*PhotonMassScale.getVal();
     double Q_MEAS_chib1_2P_1S=Q_PDG_chib1_2P_1S*PhotonMassScale2P.getVal();
     double Q_MEAS_chib2_2P_1S=Q_PDG_chib2_2P_1S*PhotonMassScale2P.getVal();
     double Q_MEAS_chib1_2P_2S=Q_PDG_chib1_2P_2S*PhotonMassScale2P_2S.getVal();
     double Q_MEAS_chib2_2P_2S=Q_PDG_chib2_2P_2S*PhotonMassScale2P_2S.getVal();
     double Q_MEAS_chib1_3P_1S=Q_PDG_chib1_3P_1S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_3P_1S=Q_PDG_chib2_3P_1S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_3P_2S=Q_PDG_chib1_3P_2S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_3P_2S=Q_PDG_chib2_3P_2S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_3P_3S=Q_PDG_chib1_3P_3S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_3P_3S=Q_PDG_chib2_3P_3S*PhotonMassScale3P.getVal();

     double Q_MEAS_chib0_1P_1S=Q_PDG_chib0_1P_1S*PhotonMassScale.getVal();
     double Q_MEAS_chib0_2P_1S=Q_PDG_chib0_2P_1S*PhotonMassScale2P.getVal();

     double M_MEAS_chib1_1P_1S=Q_MEAS_chib1_1P_1S+Ymass1S;
     double M_MEAS_chib2_1P_1S=Q_MEAS_chib2_1P_1S+Ymass1S;
     double M_MEAS_chib1_2P_1S=Q_MEAS_chib1_2P_1S+Ymass1S;
     double M_MEAS_chib2_2P_1S=Q_MEAS_chib2_2P_1S+Ymass1S;
     double M_MEAS_chib1_2P_2S=Q_MEAS_chib1_2P_2S+Ymass2S;
     double M_MEAS_chib2_2P_2S=Q_MEAS_chib2_2P_2S+Ymass2S;
     double M_MEAS_chib1_3P_1S=Q_MEAS_chib1_3P_1S+Ymass1S;
     double M_MEAS_chib2_3P_1S=Q_MEAS_chib2_3P_1S+Ymass1S;
     double M_MEAS_chib1_3P_2S=Q_MEAS_chib1_3P_2S+Ymass2S;
     double M_MEAS_chib2_3P_2S=Q_MEAS_chib2_3P_2S+Ymass2S;
     double M_MEAS_chib1_3P_3S=Q_MEAS_chib1_3P_3S+Ymass3S;
     double M_MEAS_chib2_3P_3S=Q_MEAS_chib2_3P_3S+Ymass3S;

     double M_MEAS_chib0_1P_1S=Q_MEAS_chib0_1P_1S+Ymass1S;
     double M_MEAS_chib0_2P_1S=Q_MEAS_chib0_2P_1S+Ymass1S;

// Q values 'reflection peaks':

     double Q_PDG_chib1_2P_2Sin1S=M_PDG_chib1_2P-Ymass2S;
     double Q_PDG_chib2_2P_2Sin1S=M_PDG_chib2_2P-Ymass2S;
     double Q_PDG_chib1_3P_2Sin1S=fmchib1_3P-Ymass2S;
     double Q_PDG_chib2_3P_2Sin1S=fmchib2_3P-Ymass2S;
     double Q_PDG_chib1_3P_3Sin1S=fmchib1_3P-Ymass3S;
     double Q_PDG_chib2_3P_3Sin1S=fmchib2_3P-Ymass3S;
     double Q_PDG_chib1_3P_3Sin2S=fmchib1_3P-Ymass3S;
     double Q_PDG_chib2_3P_3Sin2S=fmchib2_3P-Ymass3S;

     double Q_MEAS_chib1_2P_2Sin1S=Q_PDG_chib1_2P_2Sin1S*PhotonMassScale2P_2S.getVal();
     double Q_MEAS_chib2_2P_2Sin1S=Q_PDG_chib2_2P_2Sin1S*PhotonMassScale2P_2S.getVal();
     double Q_MEAS_chib1_3P_2Sin1S=Q_PDG_chib1_3P_2Sin1S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_3P_2Sin1S=Q_PDG_chib2_3P_2Sin1S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_3P_3Sin1S=Q_PDG_chib1_3P_3Sin1S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_3P_3Sin1S=Q_PDG_chib2_3P_3Sin1S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_3P_3Sin2S=Q_PDG_chib1_3P_3Sin2S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_3P_3Sin2S=Q_PDG_chib2_3P_3Sin2S*PhotonMassScale3P.getVal();

     double M_MEAS_chib1_2P_2Sin1S=Q_MEAS_chib1_2P_2Sin1S+Ymass1S;
     double M_MEAS_chib2_2P_2Sin1S=Q_MEAS_chib2_2P_2Sin1S+Ymass1S;
     double M_MEAS_chib1_3P_2Sin1S=Q_MEAS_chib1_3P_2Sin1S+Ymass1S;
     double M_MEAS_chib2_3P_2Sin1S=Q_MEAS_chib2_3P_2Sin1S+Ymass1S;
     double M_MEAS_chib1_3P_3Sin1S=Q_MEAS_chib1_3P_3Sin1S+Ymass1S;
     double M_MEAS_chib2_3P_3Sin1S=Q_MEAS_chib2_3P_3Sin1S+Ymass1S;
     double M_MEAS_chib1_3P_3Sin2S=Q_MEAS_chib1_3P_3Sin2S+Ymass2S;
     double M_MEAS_chib2_3P_3Sin2S=Q_MEAS_chib2_3P_3Sin2S+Ymass2S;

// Q values 'mirror peaks':

     double Q_PDG_chib1_1P_1Sfrom3S=Ymass3S-M_PDG_chib1_1P;
     double Q_PDG_chib2_1P_1Sfrom3S=Ymass3S-M_PDG_chib2_1P;
     double Q_PDG_chib1_1P_1Sfrom2S=Ymass2S-M_PDG_chib1_1P;
     double Q_PDG_chib2_1P_1Sfrom2S=Ymass2S-M_PDG_chib2_1P;
     double Q_PDG_chib1_2P_1Sfrom3S=Ymass3S-M_PDG_chib1_2P;
     double Q_PDG_chib2_2P_1Sfrom3S=Ymass3S-M_PDG_chib2_2P;
     double Q_PDG_chib1_2P_2Sfrom3S=Ymass3S-M_PDG_chib1_2P;
     double Q_PDG_chib2_2P_2Sfrom3S=Ymass3S-M_PDG_chib2_2P;

     double Q_MEAS_chib1_1P_1Sfrom3S=Q_PDG_chib1_1P_1Sfrom3S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_1P_1Sfrom3S=Q_PDG_chib2_1P_1Sfrom3S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_1P_1Sfrom2S=Q_PDG_chib1_1P_1Sfrom2S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_1P_1Sfrom2S=Q_PDG_chib2_1P_1Sfrom2S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_2P_1Sfrom3S=Q_PDG_chib1_2P_1Sfrom3S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_2P_1Sfrom3S=Q_PDG_chib2_2P_1Sfrom3S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib1_2P_2Sfrom3S=Q_PDG_chib1_2P_2Sfrom3S*PhotonMassScale3P.getVal();
     double Q_MEAS_chib2_2P_2Sfrom3S=Q_PDG_chib2_2P_2Sfrom3S*PhotonMassScale3P.getVal();

     double M_MEAS_chib1_1P_1Sfrom3S=Q_MEAS_chib1_1P_1Sfrom3S+Ymass1S;
     double M_MEAS_chib2_1P_1Sfrom3S=Q_MEAS_chib2_1P_1Sfrom3S+Ymass1S;
     double M_MEAS_chib1_1P_1Sfrom2S=Q_MEAS_chib1_1P_1Sfrom2S+Ymass1S;
     double M_MEAS_chib2_1P_1Sfrom2S=Q_MEAS_chib2_1P_1Sfrom2S+Ymass1S;
     double M_MEAS_chib1_2P_1Sfrom3S=Q_MEAS_chib1_2P_1Sfrom3S+Ymass1S;
     double M_MEAS_chib2_2P_1Sfrom3S=Q_MEAS_chib2_2P_1Sfrom3S+Ymass1S;
     double M_MEAS_chib1_2P_2Sfrom3S=Q_MEAS_chib1_2P_2Sfrom3S+Ymass2S;
     double M_MEAS_chib2_2P_2Sfrom3S=Q_MEAS_chib2_2P_2Sfrom3S+Ymass2S;

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


     hCut_1Psig->SetBinContent(inCut+1,rsigchib1P);
     if(rratchib1P<100) hCut_1PSB->SetBinContent(inCut+1,rratchib1P);
     hCut_1PS->SetBinContent(inCut+1,rchib1P);
     hCut_1PB->SetBinContent(inCut+1,rbgchib1P);
     hCut_2Psig->SetBinContent(inCut+1,rsigchib2P);
     if(rratchib2P<100) hCut_2PSB->SetBinContent(inCut+1,rratchib2P);
     hCut_2PS->SetBinContent(inCut+1,rchib2P);
     hCut_2PB->SetBinContent(inCut+1,rbgchib2P);
     hCut_sig->SetBinContent(inCut+1,(rsigchib1P+rsigchib2P)/2.);
     hCut_covQual1->SetBinContent(inCut+1,covQualOptimization);
     hCut_1Pmean ->SetBinContent(inCut+1,fmchib1P);
     hCut_1Pwidth->SetBinContent(inCut+1,sigmaInd.getVal());
     hCut_1Pnevt ->SetBinContent(inCut+1,fnchib1P);
     hCut_2Pmean ->SetBinContent(inCut+1,fmchib2P);
     hCut_2Pwidth->SetBinContent(inCut+1,sigmaInd2P.getVal());
     hCut_2Pnevt ->SetBinContent(inCut+1,fnchib2P);
     hCut_BGnevt1S ->SetBinContent(inCut+1,fnbg1S);
     hCut_BGnevt2S ->SetBinContent(inCut+1,fnbg2S);
     hCut_Ntotal1S ->SetBinContent(inCut+1,data1S_masswindow->sumEntries());
     hCut_Ntotal2S ->SetBinContent(inCut+1,data2S_masswindow->sumEntries());
     if(rchib1P>0.9) hCut_2PSover1PS->SetBinContent(inCut+1,rchib2P/rchib1P);

     hCut_3Psig->SetBinContent(inCut+1,ThreePsig);
     if(rratchib3P<100) hCut_3PSB->SetBinContent(inCut+1,rratchib3P);
     hCut_3PS->SetBinContent(inCut+1,rchib3P);
     hCut_3PB->SetBinContent(inCut+1,rbgchib3P);
     hCut_3Pmass->SetBinContent(inCut+1,fmchib3P);
     hCut_3Pmass->SetBinError(inCut+1,err_fmchib3P);
     hCut_3PMerr->SetBinContent(inCut+1,err_fmchib3P);

     hCut_CBa    ->SetBinContent(inCut+1,alpha_3J.getVal());
     hCut_CBa    ->SetBinError(inCut+1,alpha_3J.getError());
     hCut_CBn    ->SetBinContent(inCut+1,n_3J.getVal());
     hCut_CBn    ->SetBinError(inCut+1,n_3J.getError());
     hCut_PhotonEnergyScale->SetBinContent(inCut+1,PhotonMassScale.getVal());
     hCut_PhotonEnergyScale->SetBinError(inCut+1,PhotonMassScale.getError());
     hCut_PhotonEnergyScale2P->SetBinContent(inCut+1,PhotonMassScale2P.getVal());
     hCut_PhotonEnergyScale2P->SetBinError(inCut+1,PhotonMassScale2P.getError());
     hCut_PhotonSigmaScale->SetBinContent(inCut+1,PhotonSigmaScale.getVal());
     hCut_PhotonSigmaScale->SetBinError(inCut+1,PhotonSigmaScale.getError());
     hCut_PhotonSigmaScale2P->SetBinContent(inCut+1,PhotonSigmaScale2P.getVal());
     hCut_PhotonSigmaScale2P->SetBinError(inCut+1,PhotonSigmaScale2P.getError());

     hCut_3Psig_1S->SetBinContent(inCut+1,ThreePsig_1S);
     if(rratchib3P_1S<100) hCut_3PSB_1S->SetBinContent(inCut+1,rratchib3P_1S);
     hCut_3PS_1S->SetBinContent(inCut+1,rchib3P_1S);
     hCut_3PB_1S->SetBinContent(inCut+1,rbgchib3P_1S);

     hCut_3PsigALL->SetBinContent(inCut+1,ThreePsigALL);

     hCut_1PsigLL->SetBinContent(inCut+1,OnePsig);
     hCut_2PsigLL->SetBinContent(inCut+1,TwoPsig);

     hCut_3Pwidth->SetBinContent(inCut+1,sigmaInd3P.getVal());
     hCut_3Pwidth_1S->SetBinContent(inCut+1,sigmaInd3P_1S.getVal());
     hCut_2Pwidth_2S->SetBinContent(inCut+1,sigmaInd2P_2S.getVal());

     hCut_2Pnevt_2S->SetBinContent(inCut+1,n2PnJ_2S.getVal());
     hCut_sigLL->SetBinContent(inCut+1,OnePTwoPsigComb);
     hCut_2Psig_2S->SetBinContent(inCut+1,TwoPsig_2S);

     hCut_sigLLhalf->SetBinContent(inCut+1,(OnePsig+OnePsig)/2.);
     hCut_3PsigLLhalf->SetBinContent(inCut+1,(ThreePsig_1S+ThreePsig)/2.);

     double linewidth = 1;
     double chi21S;
     double chi22S;
     if(MCsample) binWidth=5.;
     int FrameBins1S=(invm_max1S-invm_min1S)*1000/binWidth;
     int FrameBins2Sin1S=(invm_max1S-invm_min1S)*1000/binWidth;
     int FrameBins2S=(invm_max1S-invm_min1S)*1000/binWidth;
     int FrameBins3Sin1S=(invm_max1S-invm_min1S)*1000/binWidth;
     int FrameBins3S=(invm_max1S-invm_min1S)*1000/binWidth;

     if(PlotStep2){
     n3PnJ_1S.setConstant(kFALSE);
     m_chib3P.setConstant(kFALSE);
     }
     if(!PlotStep2){
     n3PnJ_1S.setConstant(kTRUE);
     m_chib3P.setConstant(kTRUE);
     }

     if(PlotStep1){
     alpha_3J.setConstant(kFALSE);
     PhotonMassScale2P.setConstant(kFALSE);
     PhotonSigmaScale.setConstant(kFALSE);
     PhotonMassScale.setConstant(kFALSE);
     PhotonSigmaScale2P.setConstant(kFALSE);
     alpha1.setConstant(kFALSE);
     beta1.setConstant(kFALSE);
     q01S.setConstant(kFALSE);
     background_nevt1S.setConstant(kFALSE);
     n2PnJ.setConstant(kFALSE);
     n1PnJ.setConstant(kFALSE);
     gamma1.setConstant(kFALSE);

     PhotonMassScale2P_2Sin1S.setConstant(kFALSE);
     PhotonSigmaScale2P_2Sin1S.setConstant(kFALSE);
     mX.setConstant(kFALSE);
     nX.setConstant(kFALSE);
    }
     if(!PlotStep1){
     alpha_3J.setConstant(kTRUE);
     PhotonMassScale2P.setConstant(kTRUE);
     PhotonSigmaScale.setConstant(kTRUE);
     PhotonMassScale.setConstant(kTRUE);
     PhotonSigmaScale2P.setConstant(kTRUE);
     alpha1.setConstant(kTRUE);
     beta1.setConstant(kTRUE);
     q01S.setConstant(kTRUE);
     background_nevt1S.setConstant(kTRUE);
     n2PnJ.setConstant(kTRUE);
     n1PnJ.setConstant(kTRUE);
     gamma1.setConstant(kTRUE);
     PhotonMassScale2P_2Sin1S.setConstant(kTRUE);
     PhotonSigmaScale2P_2Sin1S.setConstant(kTRUE);
     }
     if(PlotStep1&&PlotStep2&&MCsample){
     alpha_3J.setConstant(kFALSE);
     PhotonMassScale2P.setConstant(kFALSE);
     PhotonSigmaScale.setConstant(kFALSE);
     PhotonMassScale.setConstant(kFALSE);
     PhotonSigmaScale2P.setConstant(kFALSE);
     alpha1.setConstant();
     beta1.setConstant();
     q01S.setConstant();
     background_nevt1S.setConstant();
     n2PnJ.setConstant(kFALSE);
     n1PnJ.setConstant(kFALSE);
     gamma1.setConstant();
     n3PnJ_1S.setConstant(kTRUE);
     m_chib3P.setConstant(kTRUE);
     mX.setConstant();
     nX.setConstant();
  	 ratio_J2overJ1.setConstant(kFALSE);
     n_3J.setConstant(kFALSE);
     PhotonMassScale2P_2Sin1S.setConstant(kFALSE);
     PhotonSigmaScale2P_2Sin1S.setConstant(kFALSE);
    }

     char plotYtitle[200];

     if(SetSignalToZero){
     n1PnJ.setVal(1e-8);
     n2PnJ.setVal(1e-8);
     n3PnJ.setVal(1e-8);
     n2PnJ_2S.setVal(1e-8);
     n3PnJ_1S.setVal(1e-8);
     }

     int CanvasSize=1900;
//     if(!OneSPlotModel) CanvasSize=1330;
    TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",CanvasSize, 800);
    ChibCanvas1S->Divide(1);
    ChibCanvas1S->SetFillColor(kWhite);
    ChibCanvas1S->cd(1);
    /*if(OneSPlotModel)*/ gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    double FontSize=0.0325;



    //    invm1S.setRange("SBregion31S",bb_thresh,invm_max1S);
    //    invm1S.setRange("Chib1P2Pregion",invm_min1S,chib2Pmax);
    //    invm1S.setRange("Chib3Pregion",chib2Pmax,bb_thresh);
    //    PlotStep1


    char PlotRangeSteps[200];
    double totalEventsInFit1S_plot;
    double PlotInt;

    if(PlotStep1) sprintf(PlotRangeSteps,"Chib1P2Pregion,SBregion31S");
    if(PlotStep2) sprintf(PlotRangeSteps,"Chib3Pregion");
    if(PlotBothSteps) sprintf(PlotRangeSteps,"Y1SAll");

    if(PlotStep1){
        PlotInt = NormalizedIntegral(&modelPdf1S, invm1S, invm_min1S, chib2Pmax)+NormalizedIntegral(&modelPdf1S, invm1S, bb_thresh, invm_max1S);
        totalEventsInFit1S_plot=totalEventsInFit1S*PlotInt/NormalizedIntegral(&modelPdf1S, invm1S, invm_min1S, invm_max1S);
    }

    if(PlotStep2){
        PlotInt = NormalizedIntegral(&modelPdf1S, invm1S, chib2Pmax, bb_thresh);
        totalEventsInFit1S_plot=totalEventsInFit1S*PlotInt/NormalizedIntegral(&modelPdf1S, invm1S, invm_min1S, invm_max1S);
    }

    if(PlotBothSteps) totalEventsInFit1S_plot=totalEventsInFit1S;
    double totalEventsInFit2S_plot=totalEventsInFit2S;
    double totalEventsInFit3S_plot=totalEventsInFit3S;

    cout<<"Plotting in Range "<<PlotRangeSteps<<endl;
    cout<<"Normalizing to "<<totalEventsInFit1S_plot<<" Events , totalEventsInFit1S = "<<totalEventsInFit1S<<endl;

    RooPlot* frame1S= invm1S.frame(FrameBins1S);
    frame1S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    frame1S->GetXaxis()->SetLimits(invm_min1S,invm_max1S);
    data1S->plotOn(frame1S,MarkerStyle(MarkerStyle_nS[1]),MarkerColor(MarkerColor_nS[1]),MarkerSize(MarkerSize_nS[1]));
    sprintf(plotYtitle,"Events per %1.1f MeV",binWidth);
    frame1S->SetYTitle(plotYtitle);
    frame1S->SetTitleOffset(1.15,"X");
    frame1S->SetTitleOffset(0.825,"Y");
    if(OneSPlotModel) modelPdf1S.plotOn(frame1S,Range(PlotRangeSteps),LineWidth(linewidth),Normalization(totalEventsInFit1S_plot,2));;
    chi21S = frame1S->chiSquare();
    /////////// Two CB adventure //////////////////
    if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(frame1S, Components("chib2P1_2Sin1S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(frame1S, Components("chib2P2_2Sin1S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(frame1S, Components("chib1P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(frame1S, Components("chib1P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(frame1S, Components("chib2P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(frame1S, Components("chib2P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep2) modelPdf1S.plotOn(frame1S, Components("chib3P1_1S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel&&PlotStep2) modelPdf1S.plotOn(frame1S, Components("chib3P2_1S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel) modelPdf1S.plotOn(frame1S, Components("background1S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    if(OneSPlotModel) modelPdf1S.plotOn(frame1S,Range(PlotRangeSteps),LineColor(LineColor_nS[1]),LineWidth(linewidth),Normalization(totalEventsInFit1S_plot,2));

    bool plotParamOnPlot=true;
    if(!OneSPlotModel) plotParamOnPlot=false;
    if(PlotStep1&&PlotStep2&&finalPlots) plotParamOnPlot=false;
    if(plotParamOnPlot) modelPdf1S.paramOn(frame1S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

    frame1S->SetTitle(0);
    frame1S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
    if(restrictMaximum) frame1S->SetMaximum(restrictMaximum_nPlots[0]);
    frame1S->Draw();

    double xText=invm_max1S-0.68;
    double highestText=frame1S->GetMaximum();
    double deltaText=0.06;
    if(MCsample) xText=9.55;

    char text[200];

    	cout<<"DRAW LATEX"<<endl;
    sprintf(text,"S/B #chi_{b}(1P) = %1.0f / %1.0f, #chi_{b}(1P)_{sig.} = %1.3f",rchib1P,rbgchib1P,OnePsig);
    TLatex text1 = TLatex(xText,highestText*(0.95-1*deltaText),text);
    text1.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&PlotBothSteps&&!finalPlots&&!MCsample)text1.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(2P) = %1.0f / %1.0f, #chi_{b}(2P)_{sig.} = %1.3f",rchib2P,rbgchib2P,TwoPsig);
    TLatex text2 = TLatex(xText,highestText*(0.95-2*deltaText),text);
    text2.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&PlotBothSteps&&!finalPlots&&!MCsample)text2.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P_1S,rbgchib3P_1S,ThreePsig_1S);
    TLatex text3 = TLatex(xText,highestText*(0.95-3*deltaText),text);
    text3.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&PlotBothSteps&&!finalPlots&&!MCsample)text3.Draw( "same" );
    sprintf(text,"Sig. #chi_{b}(2P)#rightarrow#Upsilon(2S)+#gamma#rightarrow#Upsilon(1S)+X+#gamma = %1.3f",Xsignificance);
    TLatex text3X = TLatex(xText,highestText*(0.95-4*deltaText),text);
    text3X.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&PlotBothSteps&&!finalPlots&&!MCsample)text3X.Draw( "same" );
    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",fmchib3P_1S,err_fmchib3P_1S);
    TLatex text0 = TLatex(xText,highestText*(0.95-5*deltaText),text)                                                                                                                            ;
    text0.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&PlotBothSteps&&!finalPlots&&!MCsample)text0.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"N_{tot} = %1.0f",data1S_masswindow->sumEntries());
    TLatex text4 = TLatex(xText,highestText*(0.95-6*deltaText),text);                                                                                                                                                                        ;
    text4.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//    if(DrawTextOnPlots)text4.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi21S);
    TLatex text6 = TLatex(xText+0.4,highestText*(0.95-0*deltaText),text)                                                                                                                                                                                          ;
    text6.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&OneSPlotModel)text6.Draw( "same" )                                                                                                                                                                                                                                                 ;

//    ChibCanvas1S->Modified();
    sprintf(saveName,"Figures/%s/InvMass_Y1S_%s.pdf",FitID,cutName_);
    ChibCanvas1S->SaveAs(saveName);




gamma2.setConstant(kFALSE);
alpha2.setConstant(kFALSE);
beta2.setConstant(kFALSE);
q02S.setConstant(kFALSE);
background_nevt2S.setConstant(kFALSE);
PhotonMassScale2P_2S.setConstant(kFALSE);
PhotonSigmaScale2P_2S.setConstant(kFALSE);

alpha_3J.setConstant(kTRUE);
PhotonMassScale2P.setConstant(kTRUE);
PhotonSigmaScale.setConstant(kTRUE);
PhotonMassScale.setConstant(kTRUE);
PhotonSigmaScale2P.setConstant(kTRUE);
alpha1.setConstant(kTRUE);
beta1.setConstant(kTRUE);
q01S.setConstant(kTRUE);
background_nevt1S.setConstant(kTRUE);
n2PnJ.setConstant(kTRUE);
n1PnJ.setConstant(kTRUE);
gamma1.setConstant(kTRUE);


    TCanvas* ChibCanvas2S = new TCanvas("#chi_{b} 2S invariant mass","#chi_{b} 2S invariant mass",1900, 800);
    ChibCanvas2S->Divide(1);
    ChibCanvas2S->SetFillColor(kWhite);
    ChibCanvas2S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    RooPlot* frame2S= invm2S.frame(invm_min1S,invm_max2S,FrameBins2S);
    frame2S->SetTitle("#chi_{b} invariant mass, Y2S decay");
    frame2S->GetXaxis()->SetLimits(invm_min1S,invm_max1S);
    frame2S->SetYTitle(plotYtitle);
    data2S->plotOn(frame2S,MarkerStyle(MarkerStyle_nS[2]),MarkerColor(MarkerColor_nS[2]),MarkerSize(MarkerSize_nS[2]));
    frame2S->SetTitleOffset(1.15,"X");
    frame2S->SetTitleOffset(0.825,"Y");

    if(OneSPlotModel){
    modelPdf2S.plotOn(frame2S,Range(invm_min2S,invm_max2S),LineWidth(linewidth),Normalization(totalEventsInFit2S,2));;
    chi22S = frame2S->chiSquare();

    modelPdf2S.plotOn(frame2S, Components("chib3P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("chib3P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("chib2P1_2S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("chib2P2_2S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("background2S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S,Range(invm_min2S,invm_max2S),LineColor(LineColor_nS[2]),LineWidth(linewidth),Normalization(totalEventsInFit2S,2));

    modelPdf2S.paramOn(frame2S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));
    }

    frame2S->SetTitle(0);
    frame2S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
    if(restrictMaximum) frame2S->SetMaximum(restrictMaximum_nPlots[1]);

    frame2S->Draw();

    xText=invm_min1S+0.1;
    highestText=frame2S->GetMaximum();
    deltaText=0.06;

    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",fmchib3P,err_fmchib3P);
    TLatex text5 = TLatex(xText,highestText*(0.95-1*deltaText),text)                                                                                                                            ;
    text5.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&!finalPlots)text5.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P,rbgchib3P,ThreePsig);
    TLatex text8 = TLatex(xText,highestText*(0.95-2*deltaText),text)                                                                                                                            ;
    text8.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&!finalPlots)text8.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"N_{tot} = %1.0f",data2S_masswindow->sumEntries());
    TLatex text9 = TLatex(xText,highestText*(0.95-3*deltaText),text);                                                                                                                                                                        ;
    text9.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//    if(DrawTextOnPlots)text9.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi22S);
    TLatex text10 = TLatex(xText,highestText*(0.95-0*deltaText),text)                                                                                                                                                                                          ;
    text10.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text10.Draw( "same" )                                                                                                                                                                                                                                                 ;


    ChibCanvas2S->Modified();
    sprintf(saveName,"Figures/%s/InvMass_Y2S_%s.pdf",FitID,cutName_);
    if(nState>1) ChibCanvas2S->SaveAs(saveName);


    n3PnJ_3S.setConstant(kFALSE);
    m_chib3P.setConstant(kFALSE);
    background_nevt3S.setConstant(kFALSE);

gamma2.setConstant(kTRUE);
alpha2.setConstant(kTRUE);
beta2.setConstant(kTRUE);
q02S.setConstant(kTRUE);
background_nevt2S.setConstant(kTRUE);
PhotonMassScale2P_2S.setConstant(kTRUE);
PhotonSigmaScale2P_2S.setConstant(kTRUE);

    TCanvas* ChibCanvas3S = new TCanvas("#chi_{b} 3S invariant mass","#chi_{b} 3S invariant mass",1900, 800);
    ChibCanvas3S->Divide(1);
    ChibCanvas3S->SetFillColor(kWhite);
    ChibCanvas3S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    RooPlot* frame3S= invm3S.frame(invm_min1S,invm_max3S,FrameBins3S);
    frame3S->SetTitle("#chi_{b} invariant mass, Y3S decay");
    frame3S->GetXaxis()->SetLimits(invm_min1S,invm_max1S);
    frame3S->SetYTitle(plotYtitle);
    data3S->plotOn(frame3S,MarkerStyle(MarkerStyle_nS[3]),MarkerColor(MarkerColor_nS[3]),MarkerSize(MarkerSize_nS[3]));
    frame3S->SetTitleOffset(1.15,"X");
    frame3S->SetTitleOffset(0.825,"Y");

    double chi23S;
    if(OneSPlotModel){
    modelPdf3S.plotOn(frame3S,Range(invm_min3S,invm_max3S),LineWidth(linewidth),Normalization(totalEventsInFit3S,2));;
    chi23S = frame3S->chiSquare();
    modelPdf3S.plotOn(frame3S, Components("chib3P1_3S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S,2));
    modelPdf3S.plotOn(frame3S, Components("chib3P2_3S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S,2));
    modelPdf3S.plotOn(frame3S, Components("background3S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S,2));
    modelPdf3S.plotOn(frame3S,Range(invm_min3S,invm_max3S),LineColor(LineColor_nS[3]),LineWidth(linewidth),Normalization(totalEventsInFit3S,2));

    modelPdf3S.paramOn(frame3S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));
    }

    frame3S->SetTitle(0);
    frame3S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
    if(restrictMaximum) frame3S->SetMaximum(restrictMaximum_nPlots[2]);

    frame3S->Draw();

    xText=invm_min1S+0.1;
    highestText=frame3S->GetMaximum();
    deltaText=0.06;

    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",fmchib3P,err_fmchib3P);
    TLatex text05 = TLatex(xText,highestText*(0.95-1*deltaText),text)                                                                                                                            ;
    text05.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&!finalPlots)text05.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P_3S,rbgchib3P_3S,ThreePsig_3S);
    TLatex text08 = TLatex(xText,highestText*(0.95-2*deltaText),text)                                                                                                                            ;
    text08.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots&&!finalPlots)text08.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"N_{tot} = %1.0f",data3S_masswindow->sumEntries());
    TLatex text09 = TLatex(xText,highestText*(0.95-3*deltaText),text);                                                                                                                                                                        ;
    text09.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//    if(DrawTextOnPlots)text09.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi23S);
    TLatex text010 = TLatex(xText,highestText*(0.95-0*deltaText),text)                                                                                                                                                                                          ;
    text010.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text010.Draw( "same" )                                                                                                                                                                                                                                                 ;


    ChibCanvas3S->Modified();
    sprintf(saveName,"Figures/%s/InvMass_Y3S_%s.pdf",FitID,cutName_);
    if(nState>2) ChibCanvas3S->SaveAs(saveName);






    TCanvas* ChibCanvas1S2S = new TCanvas("#chi_{b} 1S2S invariant mass","#chi_{b} 1S2S invariant mass",1900,800);
    ChibCanvas1S2S->Divide(1);
    ChibCanvas1S2S->SetFillColor(kWhite);
    ChibCanvas1S2S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    RooPlot* frame1S2S= invm1S.frame(FrameBins1S);
    frame1S2S->SetYTitle(plotYtitle);
    frame1S2S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    frame1S2S->GetXaxis()->SetLimits(invm_min1S,invm_max1S);
    frame1S2S->SetTitleOffset(1.15,"X");
    frame1S2S->SetTitleOffset(0.825,"Y");

    data1S->plotOn(frame1S2S,MarkerStyle(MarkerStyle_nS[1]),MarkerColor(MarkerColor_nS[1]),MarkerSize(MarkerSize_nS[1]));

    if(OneSPlotModel){
    modelPdf1S.plotOn(frame1S2S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));;
    modelPdf1S.plotOn(frame1S2S, Components("chib2P_sigInd_2Sin1S"), LineStyle(2),LineColor(LineColor_nP[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
    modelPdf1S.plotOn(frame1S2S, Components("chib1P_sig"), LineStyle(2),LineColor(LineColor_nP[1]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S, Components("chib2P_sig"), LineStyle(2),LineColor(LineColor_nP[2]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S, Components("chib3P_sig_1S"), LineStyle(2),LineColor(LineColor_nP[3]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S, Components("background1S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S,Range(invm_min1S,invm_max1S),LineColor(LineColor_nS[1]),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));

//    modelPdf1S.paramOn(frame1S2S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));
    }

    alpha_3J.setConstant();
    n_3J.setConstant();
    n1PnJ.setConstant();
    n2PnJ.setConstant();
    PhotonMassScale2P.setConstant();
    ratio_J2overJ1.setConstant();
    PhotonSigmaScale.setConstant();
    PhotonMassScale.setConstant();

    RooPlot* frame1S2S_= invm2S.frame(FrameBins2Sin1S);
    frame1S2S_->SetTitle("#chi_{b} invariant mass, Y1S decay");
    frame1S2S_->SetYTitle(plotYtitle);
    frame1S2S_->GetXaxis()->SetLimits(invm_min1S,invm_max1S);
    data2S->plotOn(frame1S2S_,MarkerStyle(MarkerStyle_nS[2]),MarkerColor(MarkerColor_nS[2]),MarkerSize(MarkerSize_nS[2]));

    if(OneSPlotModel){
    modelPdf2S.plotOn(frame1S2S_, Components("chib3P_sig"), LineStyle(2),LineColor(LineColor_nP[3]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame1S2S_, Components("chib2P_sig_2S"), LineStyle(2),LineColor(LineColor_nP[2]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame1S2S_, Components("background2S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame1S2S_,Range(invm_min2S,invm_max2S),LineWidth(linewidth),LineColor(LineColor_nS[2]),Normalization(totalEventsInFit2S,2));

//    modelPdf2S.paramOn(frame1S2S_, Layout(0.725,0.9875,0.275), Format("NE",AutoPrecision(2)));
    }

    RooPlot* frame1S2S__= invm3S.frame(FrameBins3Sin1S);
    frame1S2S__->SetTitle("#chi_{b} invariant mass, Y3S decay");
    frame1S2S__->SetYTitle(plotYtitle);
    frame1S2S__->SetTitleOffset(1.15,"X");
    frame1S2S__->SetTitleOffset(0.825,"Y");
    frame1S2S__->GetXaxis()->SetLimits(invm_min1S,invm_max1S);

    data3S->plotOn(frame1S2S__,MarkerStyle(MarkerStyle_nS[3]),MarkerColor(MarkerColor_nS[3]),MarkerSize(MarkerSize_nS[3]));
    if(OneSPlotModel){
    modelPdf3S.plotOn(frame1S2S__, Components("chib3P_sig_3S"), LineStyle(2),LineColor(LineColor_nP[3]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S,2));
    modelPdf3S.plotOn(frame1S2S__, Components("background3S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S,2));
    modelPdf3S.plotOn(frame1S2S__,Range(invm_min3S,invm_max3S),LineWidth(linewidth),LineColor(LineColor_nS[3]),Normalization(totalEventsInFit3S,2));

//    modelPdf3S.paramOn(frame1S2S_, Layout(0.725,0.9875,0.275), Format("NE",AutoPrecision(2)));
    }

    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(nS)^{PDG}} [GeV]");

    frame1S2S->GetXaxis()->SetTitle(invmName);
    frame1S2S->SetTitle(0);
    frame1S2S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
    if(restrictMaximum) frame1S2S->SetMaximum(restrictMaximum_nPlots[3]);

    frame1S2S->Draw();
    frame1S2S_->SetTitle(0);
    frame1S2S_->SetMinimum(PlotMinimum);
    frame1S2S_->Draw("same");
    frame1S2S__->SetTitle(0);
    frame1S2S__->SetMinimum(PlotMinimum);
    if(nState>2) frame1S2S__->Draw("same");


    xText=invm_max1S-0.8;
    highestText=frame1S->GetMaximum();
    deltaText=0.06;

    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",fmchib3P,err_fmchib3P);
    TLatex text14 = TLatex(xText,highestText*(0.95-0*deltaText),text)                                                                                                                            ;
    text14.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text14.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Y(3S) + #gamma: S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P_3S,rbgchib3P_3S,ThreePsig_3S);
    TLatex text111 = TLatex(xText,highestText*(0.95-3*deltaText),text)                                                                                                                            ;
    text111.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text111.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Y(2S) + #gamma: S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P,rbgchib3P,ThreePsig);
    TLatex text11 = TLatex(xText,highestText*(0.95-2*deltaText),text)                                                                                                                            ;
    text11.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text11.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Y(1S) + #gamma: S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P_1S,rbgchib3P_1S,ThreePsig_1S);
    TLatex text12 = TLatex(xText,highestText*(0.95-1*deltaText),text)                                                                                                                            ;
    text12.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text12.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Combined Significance: #chi_{b}(3P)_{sig.} = %1.3f",ThreePsigALL);
    TLatex text13 = TLatex(xText,highestText*(0.95-4*deltaText),text)                                                                                                                            ;
    text13.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)text13.Draw( "same" )                                                                                                                                                                                                                                                 ;

    ChibCanvas1S2S->Modified();
    sprintf(saveName,"Figures/%s/InvMass_Y1S2S3S_%s.pdf",FitID,cutName_);
    if(nState>1) ChibCanvas1S2S->SaveAs(saveName);





    alpha_3J.setConstant(kFALSE);
    n_3J.setConstant(kFALSE);
    n1PnJ.setConstant(kFALSE);
    n2PnJ.setConstant();
    PhotonMassScale2P.setConstant();
    PhotonMassScale3P.setConstant();
    ratio_J2overJ1.setConstant(kFALSE);
    PhotonSigmaScale.setConstant(kFALSE);
    PhotonMassScale.setConstant(kFALSE);
    alpha_3J_2P.setConstant();
    alpha_3J_3P.setConstant();
    PhotonSigmaScale2P.setConstant();
    PhotonSigmaScale3P.setConstant();
    n3PnJ_1S.setConstant();
    m_chib3P_1S.setConstant();
    m_chib3P.setConstant();
    n_3J_2P.setConstant();


    bool DrawZoom=false;
    if(DrawZoom){

    	double ZoomPlotMin=9.75;
    	double ZoomPlotMax=10.05;
    	double ScaleZoomBinsTo=5.;
        int FrameBins1SZoom=(ZoomPlotMax-ZoomPlotMin)*1000./ScaleZoomBinsTo;//400/binWidth


    TH1F  *data1SZommHist = new TH1F("data1SZommHist","",FrameBins1SZoom,ZoomPlotMin,ZoomPlotMax);//9.7,10.1);

    tree1S->Draw("invm1S>>data1SZommHist");

    int nBinsReplace=10;//20;
    int nRebin=2;//4;

    double invm1SCentre[FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin];
    double nEventsZoom[FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin];
    double err_invm1SCentre[FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin];
    double err_up_nEventsZoom[FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin];
    double err_low_nEventsZoom[FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin];

	  for(int invm1SBin=0;invm1SBin<FrameBins1SZoom-nBinsReplace;invm1SBin++){
		  invm1SCentre[invm1SBin]=data1SZommHist->GetXaxis()->GetBinCenter(invm1SBin+1);
		  nEventsZoom[invm1SBin] = data1SZommHist->GetBinContent(invm1SBin+1);
		  err_invm1SCentre[invm1SBin] = data1SZommHist->GetBinWidth(invm1SBin+1)/2;
		  err_up_nEventsZoom[invm1SBin] = data1SZommHist->GetBinError(invm1SBin+1);
		  err_low_nEventsZoom[invm1SBin] = data1SZommHist->GetBinError(invm1SBin+1);
//		  cout<<"invm1SBin "<<invm1SBin<<endl;
//		  cout<<"invm1SCentre[invm1SBin] "<<invm1SCentre[invm1SBin]<<endl;
//		  cout<<"nEventsZoom[invm1SBin] "<<nEventsZoom[invm1SBin]<<endl;
	  }

	  int nAdd=0;
	  for(int invm1SBin=FrameBins1SZoom-nBinsReplace;invm1SBin<FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin;invm1SBin++){
		  double invm1SCentreBuffer=0;
		  double nEventsZoomBuffer=0;
		  cout<<"invm1SBin "<<invm1SBin<<endl;
		  for(int replace=0;replace<nRebin;replace++){
			  invm1SCentreBuffer=invm1SCentreBuffer+data1SZommHist->GetXaxis()->GetBinCenter(invm1SBin+replace+nAdd+1);
			  nEventsZoomBuffer=nEventsZoomBuffer+data1SZommHist->GetBinContent(invm1SBin+replace+nAdd+1);
			  cout<<"invm1SCentreBuffer "<<invm1SCentreBuffer<<endl;
			  cout<<"nEventsZoomBuffer "<<nEventsZoomBuffer<<endl;
		  }
		  invm1SCentre[invm1SBin]=invm1SCentreBuffer/nRebin;
		  nEventsZoom[invm1SBin] = nEventsZoomBuffer/nRebin;
		  err_invm1SCentre[invm1SBin] = data1SZommHist->GetBinWidth(invm1SBin+1)*nRebin/2;
		  err_up_nEventsZoom[invm1SBin] = TMath::Sqrt(nEventsZoom[invm1SBin]);
		  err_low_nEventsZoom[invm1SBin] = TMath::Sqrt(nEventsZoom[invm1SBin]);
		  cout<<"invm1SCentre[invm1SBin] "<<invm1SCentre[invm1SBin]<<endl;
		  cout<<"nEventsZoom[invm1SBin] "<<nEventsZoom[invm1SBin]<<endl;
		  nAdd=nAdd+nRebin-1;
	  }

	  FrameBins1SZoom=(ZoomPlotMax-ZoomPlotMin)*1000./ScaleZoomBinsTo;

	TGraphAsymmErrors *zoomGraphinvm1S = new TGraphAsymmErrors(FrameBins1SZoom-nBinsReplace+nBinsReplace/nRebin,invm1SCentre,nEventsZoom,err_invm1SCentre,err_invm1SCentre,err_up_nEventsZoom,err_low_nEventsZoom);
	zoomGraphinvm1S->SetMarkerSize(MarkerSize_nS[1]);
	  zoomGraphinvm1S->SetMarkerStyle(MarkerStyle_nS[1]);
	  zoomGraphinvm1S->SetMarkerColor(MarkerColor_nS[1]);


    TCanvas* ZoomChibCanvas1S = new TCanvas("#chi_{b} 1S2S invariant mass zoom","#chi_{b} 1S2S invariant mass zoom",2100,800);
    ZoomChibCanvas1S->Divide(1);
    ZoomChibCanvas1S->SetFillColor(kWhite);
    ZoomChibCanvas1S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    RooPlot* Zoomframe1S= invm1S.frame(ZoomPlotMin,ZoomPlotMax,FrameBins1SZoom);
    Zoomframe1S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    data1S->plotOn(Zoomframe1S,MarkerColor(MarkerColor_nS[1]),MarkerStyle(MarkerStyle_nS[1]),MarkerSize(MarkerSize_nS[1]));
    sprintf(plotYtitle,"Events per %1.1f  MeV",ScaleZoomBinsTo);
    Zoomframe1S->SetYTitle(plotYtitle);
        modelPdf1S.plotOn(Zoomframe1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));;
       double chi21Szoom = Zoomframe1S->chiSquare();
    /////////// Two CB adventure //////////////////
       modelPdf1S.plotOn(Zoomframe1S, Components("chib3P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib3P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib1P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib1P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib2P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib2P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib2P1_2Sin1S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
       modelPdf1S.plotOn(Zoomframe1S, Components("chib2P2_2Sin1S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(Zoomframe1S, Components("background1S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(Zoomframe1S,Range(invm_min1S,invm_max1S),LineColor(LineColor_nS[1]),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));

//    modelPdf1S.paramOn(Zoomframe1S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

    Zoomframe1S->SetTitle(0);
    Zoomframe1S->SetMinimum(PlotMinimum);
//    Zoomframe1S->SetMaximum(75.*ScaleZoomBinsTo/5.);
    Zoomframe1S->Draw();

//    zoomGraphinvm1S->Draw("P");

    xText=9.925;
    highestText=Zoomframe1S->GetMaximum();
    deltaText=0.06;
    FontSize=0.0425;


    	cout<<"DRAW LATEX"<<endl;

        sprintf(text,"CMS preliminary");
        TLatex zoomtext4 = TLatex(xText,highestText*(0.95-1*deltaText),text);
        zoomtext4.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)zoomtext4.Draw( "same" )                                                                                                                                                                                                                                                 ;
        sprintf(text,"L_{int} = 4.7 fb^{-1}");
        TLatex zoomtext5 = TLatex(xText,highestText*(0.95-2*deltaText),text);
        zoomtext5.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)zoomtext5.Draw( "same" )                                                                                                                                                                                                                                                 ;
        FontSize=0.0325;

    sprintf(text,"#chi_{bj}(1P) width = %1.1f #pm %1.1f MeV",onePwidth*1000,onePwidtherr*1000);
    TLatex zoomtext = TLatex(xText,highestText*(0.95-3*deltaText),text);
    zoomtext.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)zoomtext.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"N_{#chi_{b2}(1P)} / N_{#chi_{b1}(1P)} = %1.3f #pm %1.3f",ratio_J2overJ1_1P.getVal(),ratio_J2overJ1_1P.getError());
    TLatex zoomtext2 = TLatex(xText,highestText*(0.95-4*deltaText),text);
    zoomtext2.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
    if(DrawTextOnPlots)zoomtext2.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi21Szoom);
     TLatex zoomtext3 = TLatex(xText,highestText*(0.95-2*deltaText),text)                                                                                                                                                                                          ;
     zoomtext3.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//     if(DrawTextOnPlots)zoomtext3.Draw( "same" )                                                                                                                                                                                                                                                 ;

    sprintf(saveName,"Figures/%s/InvMass_Zoom_Y1S_%s.pdf",FitID,cutName_);
    ZoomChibCanvas1S->SaveAs(saveName);
    delete ZoomChibCanvas1S;
    delete zoomGraphinvm1S;
    delete data1SZommHist;
	delete Zoomframe1S;
    }



    bool DrawNonEqu=true;
    if(DrawNonEqu){



        int nBins_buff;
        int nBins_buff2S;
        int nBins_buff3S;
        double MinDistancetoMax=0.1;

	double resolutionMin=0.01;
    double nResolution=1.;
    double BinBorders[1000];
    BinBorders[0]=invm_min1S;
    for(int iBin=1;iBin<1000;iBin++){
    	BinBorders[iBin]=BinBorders[iBin-1]+nResolution*(BinBorders[iBin-1]-Ymass1S)*PhotonSigmaScale3P.getVal();
    	if(BinBorders[iBin]-BinBorders[iBin-1]<resolutionMin){
    		BinBorders[iBin]=BinBorders[iBin-1]+resolutionMin;
    		cout<<"using min resolution = "<<resolutionMin<<" at "<<BinBorders[iBin]<<" GeV "<<endl;
    	}
    	if(BinBorders[iBin-1]>9.9&&BinBorders[iBin-1]<10.15){
    		BinBorders[iBin]=BinBorders[iBin-1]+2.5*nResolution*(BinBorders[iBin-1]-Ymass1S)*PhotonSigmaScale3P.getVal();
    	}
    	if(BinBorders[iBin-1]>bb_thresh){
    		BinBorders[iBin]=BinBorders[iBin-1]+2.5*nResolution*(BinBorders[iBin-1]-Ymass1S)*PhotonSigmaScale3P.getVal();
    	}
//Carlos additions:
		BinBorders[iBin]=BinBorders[iBin-1]+0.01;
    	if(BinBorders[iBin-1]>10.){
    		BinBorders[iBin]=BinBorders[iBin-1]+0.02;
    	}
    	if(BinBorders[iBin-1]>10.6){
    		BinBorders[iBin]=BinBorders[iBin-1]+0.04;
    	}
/////////
    	if(BinBorders[iBin]>invm_max1S-MinDistancetoMax) {
    		nBins_buff=iBin; BinBorders[iBin]=invm_max1S; break;
    	}
    }

    resolutionMin=0.02;
    double BinBorders2S[1000];
    BinBorders2S[0]=invm_min2S;
    BinBorders2S[1]=invm_min2S+0.01;
    for(int iBin=2;iBin<1000;iBin++){
    	BinBorders2S[iBin]=BinBorders2S[iBin-1]+nResolution*(BinBorders2S[iBin-1]-Ymass2S)*PhotonSigmaScale3P.getVal();
    	if(BinBorders2S[iBin]-BinBorders2S[iBin-1]<resolutionMin){
    		BinBorders2S[iBin]=BinBorders2S[iBin-1]+resolutionMin;
    		cout<<"using min resolution = "<<resolutionMin<<" at "<<BinBorders2S[iBin]<<" GeV "<<endl;
    	}
    	if(BinBorders2S[iBin-1]>10.5999){
    		BinBorders2S[iBin]=BinBorders2S[iBin-1]+0.1;
    	}
    	if(BinBorders2S[iBin]>invm_max2S) {
    		nBins_buff2S=iBin; BinBorders2S[iBin]=invm_max2S; break;
    	}
    }

    resolutionMin=0.02;
    double BinBorders3S[1000];
    BinBorders3S[0]=invm_min3S;
    BinBorders3S[1]=invm_min3S+0.01;
    for(int iBin=2;iBin<1000;iBin++){
    	BinBorders3S[iBin]=BinBorders3S[iBin-1]+nResolution*(BinBorders3S[iBin-1]-Ymass3S)*PhotonSigmaScale3P.getVal();
    	if(BinBorders3S[iBin]-BinBorders3S[iBin-1]<resolutionMin){
    		BinBorders3S[iBin]=BinBorders3S[iBin-1]+resolutionMin;
    		cout<<"using min resolution = "<<resolutionMin<<" at "<<BinBorders3S[iBin]<<" GeV "<<endl;
    	}
    	if(BinBorders3S[iBin-1]>10.5999){
    		BinBorders3S[iBin]=BinBorders3S[iBin-1]+0.1;
    	}
    	if(BinBorders3S[iBin]>invm_max3S) {
    		nBins_buff3S=iBin; BinBorders3S[iBin]=invm_max3S; break;
    	}
    }

    const int nBins=nBins_buff;
    const int nBins2S=nBins_buff2S;
    const int nBins3S=nBins_buff3S;
    TH1F  *data1SNonEquHist = new TH1F("data1SNonEquHist","",nBins,invm_min1S,invm_max1S);//9.7,10.1);
    TH1F  *data2SNonEquHist = new TH1F("data2SNonEquHist","",nBins2S,invm_min2S,invm_max2S);//9.7,10.1);
    TH1F  *data3SNonEquHist = new TH1F("data3SNonEquHist","",nBins3S,invm_min3S,invm_max3S);//9.7,10.1);

    double xBins[nBins+1];
    for(int iBin=0;iBin<nBins+1;iBin++){
    	xBins[iBin]=BinBorders[iBin];
    }

    double xBins2S[nBins+1];
    for(int iBin=0;iBin<nBins2S+1;iBin++){
    	xBins2S[iBin]=BinBorders2S[iBin];
    }

    double xBins3S[nBins+1];
    for(int iBin=0;iBin<nBins3S+1;iBin++){
    	xBins3S[iBin]=BinBorders3S[iBin];
    }

    data1SNonEquHist->GetXaxis()->Set(nBins, xBins);
    data2SNonEquHist->GetXaxis()->Set(nBins2S, xBins2S);
    data3SNonEquHist->GetXaxis()->Set(nBins3S, xBins3S);


    tree1S->Draw("invm1S>>data1SNonEquHist");
    tree2S->Draw("invm2S>>data2SNonEquHist");
    tree3S->Draw("invm3S>>data3SNonEquHist");

        double ScaleBinsTo=ScalePlot1S_NE/1000.;
        double ScaledContent;
        for(int iBin=1;iBin<nBins+1;iBin++){
        	ScaledContent=data1SNonEquHist->GetBinContent(iBin)*ScaleBinsTo/data1SNonEquHist->GetBinWidth(iBin);
        	data1SNonEquHist->SetBinContent(iBin,ScaledContent);
        }
        double ScaleBinsTo2S=ScalePlot2S_NE/1000.;
        for(int iBin=1;iBin<nBins2S+1;iBin++){
        	ScaledContent=data2SNonEquHist->GetBinContent(iBin)*ScaleBinsTo2S/data2SNonEquHist->GetBinWidth(iBin);
        	data2SNonEquHist->SetBinContent(iBin,ScaledContent);
        }
        double ScaleBinsTo3S=ScalePlot3S_NE/1000.;
        for(int iBin=1;iBin<nBins3S+1;iBin++){
        	ScaledContent=data3SNonEquHist->GetBinContent(iBin)*ScaleBinsTo3S/data3SNonEquHist->GetBinWidth(iBin);
        	data3SNonEquHist->SetBinContent(iBin,ScaledContent);
        }


        data1SNonEquHist->SetMarkerSize(MarkerSize_nS[1]);
        data1SNonEquHist->SetMarkerStyle(MarkerStyle_nS[1]);
        data1SNonEquHist->SetMarkerColor(LineColor_nS[1]);

        data2SNonEquHist->SetMarkerSize(MarkerSize_nS[2]);
        data2SNonEquHist->SetMarkerStyle(MarkerStyle_nS[2]);
        data2SNonEquHist->SetMarkerColor(LineColor_nS[2]);

        data3SNonEquHist->SetMarkerSize(MarkerSize_nS[3]);
        data3SNonEquHist->SetMarkerStyle(MarkerStyle_nS[3]);
        data3SNonEquHist->SetMarkerColor(LineColor_nS[3]);

        RooDataHist NonEquHist("NonEquHist", "NonEquHist",RooArgList(invm1S), data1SNonEquHist);
    	RooDataHist NonEquHist2S("NonEquHist2S", "NonEquHist2S",RooArgList(invm2S), data2SNonEquHist);
    	RooDataHist NonEquHist3S("NonEquHist3S", "NonEquHist3S",RooArgList(invm3S), data3SNonEquHist);

    TCanvas* NonEquChibCanvas1S = new TCanvas("#chi_{b} 1S2S invariant mass NonEqu","#chi_{b} 1S2S invariant mass NonEqu",1900,800);
    NonEquChibCanvas1S->Divide(1);
    NonEquChibCanvas1S->SetFillColor(kWhite);
    NonEquChibCanvas1S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

	int FrameBins1SNonEqu=(invm_max1S-invm_min1S)/ScaleBinsTo;

    RooPlot* NonEquframe1S= invm1S.frame(invm_min1S,invm_max1S,FrameBins1SNonEqu);
    NonEquframe1S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    sprintf(plotYtitle,"Events per %1.1f MeV",ScaleBinsTo*1000.);
    NonEquframe1S->SetYTitle(plotYtitle);
    NonEquframe1S->SetTitleOffset(1.15,"X");
    NonEquframe1S->SetTitleOffset(0.825,"Y");
//    data1S->plotOn(NonEquframe1S);
//    NonEquHist.plotOn(NonEquframe1S);
        modelPdf1S.plotOn(NonEquframe1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));;
       double chi21SNonEqu = NonEquframe1S->chiSquare();
    /////////// Two CB adventure //////////////////
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("chib2P1_2Sin1S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("chib2P2_2Sin1S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("chib1P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("chib1P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("chib2P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("chib2P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep2) modelPdf1S.plotOn(NonEquframe1S, Components("chib3P1_1S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep2) modelPdf1S.plotOn(NonEquframe1S, Components("chib3P2_1S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S, Components("background1S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
//       data1SNonEquHist->Draw("same,E");
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S,Range(PlotRangeSteps),LineColor(LineColor_nS[1]),LineWidth(linewidth),Normalization(totalEventsInFit1S_plot,2));

//       if(OneSPlotModel) modelPdf1S.paramOn(NonEquframe1S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

//    modelPdf1S.paramOn(NonEquframe1S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

    NonEquframe1S->SetTitle(0);
    NonEquframe1S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
    if(restrictMaximum) NonEquframe1S->SetMaximum(restrictMaximum_nPlots[4]);
    NonEquframe1S->Draw();

    data1SNonEquHist->Draw("same,E");

    NonEquframe1S->Draw("same");

//    NonEquGraphinvm1S->Draw("P");
//    data1SNonEquHist->Draw("same,E");
    xText=invm_max1S-0.4;
    highestText=NonEquframe1S->GetMaximum();
    deltaText=0.06;
    FontSize=0.0425;


    	cout<<"DRAW LATEX"<<endl;

        sprintf(text,"CMS preliminary");
        TLatex NonEqutext4 = TLatex(xText,highestText*(0.95-1*deltaText),text);
        NonEqutext4.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext4.Draw( "same" )                                                                                                                                                                                                                                                 ;
        sprintf(text,"L_{int} = 4.7 fb^{-1}");
        TLatex NonEqutext5 = TLatex(xText,highestText*(0.95-2*deltaText),text);
        NonEqutext5.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext5.Draw( "same" )                                                                                                                                                                                                                                                 ;
        FontSize=0.0325;

    sprintf(text,"#chi_{bj}(1P) width = %1.1f MeV",onePwidth*1000);
    TLatex NonEqutext = TLatex(xText,highestText*(0.95-3*deltaText),text);
    NonEqutext.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//    if(DrawTextOnPlots)NonEqutext.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"N_{#chi_{b2}(1P)} / N_{#chi_{b1}(1P)} = %1.3f #pm %1.3f",ratio_J2overJ1.getVal(),ratio_J2overJ1.getError());
    TLatex NonEqutext2 = TLatex(xText,highestText*(0.95-4*deltaText),text);
    NonEqutext2.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//    if(DrawTextOnPlots)NonEqutext2.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi21SNonEqu);
     TLatex NonEqutext3 = TLatex(xText,highestText*(0.95-3*deltaText),text)                                                                                                                                                                                          ;
     NonEqutext3.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
//     if(DrawTextOnPlots)NonEqutext3.Draw( "same" )                                                                                                                                                                                                                                                 ;

    sprintf(saveName,"Figures/%s/InvMass_NonEqu_Y1S_%s.pdf",FitID,cutName_);
    NonEquChibCanvas1S->SaveAs(saveName);





    TCanvas* NonEquChibCanvas2S = new TCanvas("#chi_{b} 2S2S invariant mass NonEqu","#chi_{b} 2S2S invariant mass NonEqu",1900,800);
    NonEquChibCanvas2S->Divide(1);
    NonEquChibCanvas2S->SetFillColor(kWhite);
    NonEquChibCanvas2S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

	int FrameBins2SNonEqu=(invm_max2S-invm_min1S)/ScaleBinsTo2S;

    RooPlot* NonEquframe2S= invm2S.frame(invm_min1S,invm_max2S,FrameBins2SNonEqu);
    NonEquframe2S->SetTitle("#chi_{b} invariant mass, Y2S decay");
    sprintf(plotYtitle,"Events per %1.1f MeV",ScaleBinsTo2S*1000.);
    NonEquframe2S->SetYTitle(plotYtitle);
    NonEquframe2S->SetTitleOffset(1.15,"X");
    NonEquframe2S->SetTitleOffset(0.825,"Y");
//    data2S->plotOn(NonEquframe2S);
//    NonEquHist.plotOn(NonEquframe2S);
        modelPdf2S.plotOn(NonEquframe2S,Range(invm_min2S,invm_max2S),LineWidth(linewidth),Normalization(totalEventsInFit2S,2));;
       double chi22SNonEqu = NonEquframe2S->chiSquare();
    /////////// Two CB adventure //////////////////
       modelPdf2S.plotOn(NonEquframe2S, Components("chib2P1_2S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S_plot,2));
       modelPdf2S.plotOn(NonEquframe2S, Components("chib2P2_2S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S_plot,2));
       modelPdf2S.plotOn(NonEquframe2S, Components("chib3P1"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S_plot,2));
       modelPdf2S.plotOn(NonEquframe2S, Components("chib3P2"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S_plot,2));
       modelPdf2S.plotOn(NonEquframe2S, Components("background2S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S_plot,2));
       modelPdf2S.plotOn(NonEquframe2S,Range(invm_min2S,invm_max2S),LineColor(LineColor_nS[2]),LineWidth(linewidth),Normalization(totalEventsInFit2S_plot,2));

//       if(OneSPlotModel) modelPdf2S.paramOn(NonEquframe2S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

//    modelPdf2S.paramOn(NonEquframe2S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

    NonEquframe2S->SetTitle(0);
    NonEquframe2S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
//    if(restrictMaximum)
    NonEquframe2S->SetMaximum(restrictMaximum_nPlots[5]);
    NonEquframe2S->Draw();

    data2SNonEquHist->Draw("same,E");

    NonEquframe2S->Draw("same");

//    NonEquGraphinvm2S->Draw("P");
//    data2SNonEquHist->Draw("same,E");
    xText=invm_max1S-0.4;
    highestText=NonEquframe2S->GetMaximum();
    deltaText=0.06;
    FontSize=0.0425;


    	cout<<"DRAW LATEX"<<endl;

        sprintf(text,"CMS preliminary");
        TLatex NonEqutext42S = TLatex(xText,highestText*(0.95-1*deltaText),text);
        NonEqutext42S.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext42S.Draw( "same" )                                                                                                                                                                                                                                                 ;
        sprintf(text,"L_{int} = 4.7 fb^{-1}");
        TLatex NonEqutext52S = TLatex(xText,highestText*(0.95-2*deltaText),text);
        NonEqutext52S.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext52S.Draw( "same" )                                                                                                                                                                                                                                                 ;
        FontSize=0.0325;

    sprintf(saveName,"Figures/%s/InvMass_NonEqu_Y2S_%s.pdf",FitID,cutName_);
    if(nState>1) NonEquChibCanvas2S->SaveAs(saveName);








    TCanvas* NonEquChibCanvas3S = new TCanvas("#chi_{b} 3S3S invariant mass NonEqu","#chi_{b} 3S3S invariant mass NonEqu",1900,800);
    NonEquChibCanvas3S->Divide(1);
    NonEquChibCanvas3S->SetFillColor(kWhite);
    NonEquChibCanvas3S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

	int FrameBins3SNonEqu=(invm_max3S-invm_min1S)/ScaleBinsTo3S;

    RooPlot* NonEquframe3S= invm3S.frame(invm_min1S,invm_max3S,FrameBins3SNonEqu);
    NonEquframe3S->SetTitle("#chi_{b} invariant mass, Y3S decay");
    sprintf(plotYtitle,"Events per %1.1f MeV",ScaleBinsTo3S*1000.);
    NonEquframe3S->SetYTitle(plotYtitle);
    NonEquframe3S->SetTitleOffset(1.15,"X");
    NonEquframe3S->SetTitleOffset(0.825,"Y");
//    data3S->plotOn(NonEquframe3S);
//    NonEquHist.plotOn(NonEquframe3S);
        modelPdf3S.plotOn(NonEquframe3S,Range(invm_min3S,invm_max3S),LineWidth(linewidth),Normalization(totalEventsInFit3S,2));;
       double chi23SNonEqu = NonEquframe3S->chiSquare();
    /////////// Two CB adventure //////////////////
       modelPdf3S.plotOn(NonEquframe3S, Components("chib3P1_3S"), LineStyle(2),LineColor(LineColor_nJ[1]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S_plot,2));
       modelPdf3S.plotOn(NonEquframe3S, Components("chib3P2_3S"), LineStyle(2),LineColor(LineColor_nJ[2]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S_plot,2));
       modelPdf3S.plotOn(NonEquframe3S, Components("background3S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S_plot,2));
       modelPdf3S.plotOn(NonEquframe3S,Range(invm_min3S,invm_max3S),LineColor(LineColor_nS[3]),LineWidth(linewidth),Normalization(totalEventsInFit3S_plot,2));

//       if(OneSPlotModel) modelPdf3S.paramOn(NonEquframe3S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

//    modelPdf3S.paramOn(NonEquframe3S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

    NonEquframe3S->SetTitle(0);
    NonEquframe3S->SetMinimum(PlotMinimum);
    if(MCsample) restrictMaximum=false;
    if(Optimize) restrictMaximum=true;
//    if(restrictMaximum)
    NonEquframe3S->SetMaximum(restrictMaximum_nPlots[6]);
    NonEquframe3S->Draw();

    data3SNonEquHist->Draw("same,E");

    NonEquframe3S->Draw("same");

//    NonEquGraphinvm3S->Draw("P");
//    data3SNonEquHist->Draw("same,E");
    xText=invm_max1S-0.4;
    highestText=NonEquframe3S->GetMaximum();
    deltaText=0.06;
    FontSize=0.0425;


    	cout<<"DRAW LATEX"<<endl;

        sprintf(text,"CMS preliminary");
        TLatex NonEqutext43S = TLatex(xText,highestText*(0.95-1*deltaText),text);
        NonEqutext43S.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext43S.Draw( "same" )                                                                                                                                                                                                                                                 ;
        sprintf(text,"L_{int} = 4.7 fb^{-1}");
        TLatex NonEqutext53S = TLatex(xText,highestText*(0.95-2*deltaText),text);
        NonEqutext53S.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext53S.Draw( "same" )                                                                                                                                                                                                                                                 ;
        FontSize=0.0325;

    sprintf(saveName,"Figures/%s/InvMass_NonEqu_Y3S_%s.pdf",FitID,cutName_);
    if(nState>2) NonEquChibCanvas3S->SaveAs(saveName);





// NonEQU-All-In-One:::

    TCanvas* NonEquChibCanvas1S2S3S = new TCanvas("#chi_{b} 1S2S3S1S2S3S invariant mass NonEqu","#chi_{b} 1S2S3S1S2S3S invariant mass NonEqu",1900,800);
    NonEquChibCanvas1S2S3S->Divide(1);
    NonEquChibCanvas1S2S3S->SetFillColor(kWhite);
    NonEquChibCanvas1S2S3S->cd(1);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);





    RooPlot* NonEquframe1S2S3S_1S= invm1S.frame(invm_min1S,invm_max1S,FrameBins1SNonEqu);
    NonEquframe1S2S3S_1S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    sprintf(plotYtitle,"Events per %1.1f MeV",ScaleBinsTo*1000.);
    NonEquframe1S2S3S_1S->SetYTitle(plotYtitle);
    NonEquframe1S2S3S_1S->SetTitleOffset(1.15,"X");
    NonEquframe1S2S3S_1S->SetTitleOffset(0.825,"Y");
        modelPdf1S.plotOn(NonEquframe1S2S3S_1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));;
    /////////// Two CB adventure //////////////////
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S2S3S_1S, Components("background1S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
       if(OneSPlotModel&&PlotStep1) modelPdf1S.plotOn(NonEquframe1S2S3S_1S,Range(PlotRangeSteps),LineColor(LineColor_nS[1]),LineWidth(linewidth),Normalization(totalEventsInFit1S_plot,2));

       NonEquframe1S2S3S_1S->SetTitle(0);
       NonEquframe1S2S3S_1S->SetMinimum(PlotMinimum);
       NonEquframe1S2S3S_1S->SetMaximum(restrictMaximum_nPlots[7]);



    RooPlot* NonEquframe1S2S3S_2S= invm2S.frame(invm_min1S,invm_max2S,FrameBins2SNonEqu);
    NonEquframe1S2S3S_2S->SetTitle("#chi_{b} invariant mass, Y2S decay");
    sprintf(plotYtitle,"Events per %1.1f MeV",ScaleBinsTo2S*1000.);
    NonEquframe1S2S3S_2S->SetYTitle(plotYtitle);
    NonEquframe1S2S3S_2S->SetTitleOffset(1.15,"X");
    NonEquframe1S2S3S_2S->SetTitleOffset(0.825,"Y");
        modelPdf2S.plotOn(NonEquframe1S2S3S_2S,Range(invm_min2S,invm_max2S),LineWidth(linewidth),Normalization(totalEventsInFit2S,2));;
    /////////// Two CB adventure //////////////////
        modelPdf2S.plotOn(NonEquframe1S2S3S_2S, Components("background2S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S_plot,2));
       modelPdf2S.plotOn(NonEquframe1S2S3S_2S,Range(invm_min2S,invm_max2S),LineColor(LineColor_nS[2]),LineWidth(linewidth),Normalization(totalEventsInFit2S_plot,2));

       NonEquframe1S2S3S_2S->SetTitle(0);
       NonEquframe1S2S3S_2S->SetMinimum(PlotMinimum);



       RooPlot* NonEquframe1S2S3S_3S= invm3S.frame(invm_min1S,invm_max3S,FrameBins3SNonEqu);
       NonEquframe1S2S3S_3S->SetTitle("#chi_{b} invariant mass, Y2S decay");
       sprintf(plotYtitle,"Events per %1.1f MeV",ScaleBinsTo3S*1000.);
       NonEquframe1S2S3S_3S->SetYTitle(plotYtitle);
       NonEquframe1S2S3S_3S->SetTitleOffset(1.15,"X");
       NonEquframe1S2S3S_3S->SetTitleOffset(0.825,"Y");
           modelPdf3S.plotOn(NonEquframe1S2S3S_3S,Range(invm_min3S,invm_max3S),LineWidth(linewidth),Normalization(totalEventsInFit3S,2));;
       /////////// Two CB adventure //////////////////
           modelPdf3S.plotOn(NonEquframe1S2S3S_3S, Components("background3S"), LineStyle(2),LineColor(LineColor_nJ[3]),LineWidth(linewidth),Range(invm_min3S,invm_max3S),Normalization(totalEventsInFit3S_plot,2));
          modelPdf3S.plotOn(NonEquframe1S2S3S_3S,Range(invm_min3S,invm_max3S),LineColor(LineColor_nS[3]),LineWidth(linewidth),Normalization(totalEventsInFit3S_plot,2));

          NonEquframe1S2S3S_3S->SetTitle(0);
          NonEquframe1S2S3S_3S->SetMinimum(PlotMinimum);


          sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(nS)^{PDG}} [GeV]");

          NonEquframe1S2S3S_1S->GetXaxis()->SetTitle(invmName);

    NonEquframe1S2S3S_1S->Draw();
    data1SNonEquHist->Draw("same,E");
    NonEquframe1S2S3S_1S->Draw("same");

    if(nState>1.5){
    NonEquframe1S2S3S_2S->Draw("same");
    data2SNonEquHist->Draw("same,E");
    NonEquframe1S2S3S_2S->Draw("same");
    }
    if(nState>2.5){
    NonEquframe1S2S3S_3S->Draw("same");
    data3SNonEquHist->Draw("same,E");
    NonEquframe1S2S3S_3S->Draw("same");
    }

    xText=invm_max1S-0.4;
    highestText=NonEquframe1S2S3S_1S->GetMaximum();
    deltaText=0.06;
    FontSize=0.0425;


    	cout<<"DRAW LATEX"<<endl;

        sprintf(text,"CMS preliminary");
        TLatex NonEqutext41S2S3S = TLatex(xText,highestText*(0.95-1*deltaText),text);
        NonEqutext41S2S3S.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext41S2S3S.Draw( "same" )                                                                                                                                                                                                                                                 ;
        sprintf(text,"L_{int} = 4.7 fb^{-1}");
        TLatex NonEqutext51S2S3S = TLatex(xText,highestText*(0.95-2*deltaText),text);
        NonEqutext51S2S3S.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
        if(DrawTextOnPlots)NonEqutext51S2S3S.Draw( "same" )                                                                                                                                                                                                                                                 ;
        FontSize=0.0325;

    sprintf(saveName,"Figures/%s/InvMass_NonEqu_Y1S2S3S_%s.pdf",FitID,cutName_);
    if(nState>1) NonEquChibCanvas1S2S3S->SaveAs(saveName);








    if(SaveAll){

            TH1F  *AN_Q1SHisto = new TH1F("AN_Q1SHisto","",int((invm_max1S-invm_min1S)/binWidth*1000),0,invm_max1S-invm_min1S);//binWidth MeV
            TH1F  *AN_Q2SHisto = new TH1F("AN_Q2SHisto","",int((invm_max1S-invm_min1S)/binWidth*1000),0,invm_max1S-invm_min1S);//binWidth MeV
            TH1F  *AN_Q3SHisto = new TH1F("AN_Q3SHisto","",int((invm_max1S-invm_min1S)/binWidth*1000),0,invm_max1S-invm_min1S);//binWidth MeV

            char PlotCutChar[1000];

            sprintf(PlotCutChar,"%s && TMath::Abs(Y1Smass_nSigma)<%f", FinalCutChar, nYsig);
            treeBeforeAllCuts->Draw("Q>>AN_Q1SHisto",PlotCutChar);
            sprintf(PlotCutChar,"%s && TMath::Abs(Y2Smass_nSigma)<%f", FinalCutChar, nYsig*Ymass2S/Ymass1S);
            treeBeforeAllCuts->Draw("Q>>AN_Q2SHisto",PlotCutChar);
            sprintf(PlotCutChar,"%s && TMath::Abs(Y3Smass_nSigma)<%f", FinalCutChar, nYsig*Ymass3S/Ymass2S);
            treeBeforeAllCuts->Draw("Q>>AN_Q3SHisto",PlotCutChar);

            TCanvas* PlotCanvas = new TCanvas("PlotCanvas","PlotCanvas",1600, 1200);
            PlotCanvas->SetFillColor(kWhite);
            gPad->SetLeftMargin(0.125);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.05);
            gPad->SetFillColor(kWhite);

            double MaxFact=1.2;
            double yOff=1.5;
            double PlotMinDist=1.;
            char saveMass[200];

        	  AN_Q1SHisto->SetMinimum(PlotMinDist); AN_Q1SHisto->SetMaximum(AN_Q1SHisto->GetMaximum()*MaxFact);
        	  MaxFact=1.;
        AN_Q1SHisto->GetYaxis()->SetTitleOffset(yOff); AN_Q1SHisto->GetXaxis()->SetTitle("m_{#mu#mu#gamma}-m_{#mu#mu} [GeV]");	AN_Q1SHisto->GetYaxis()->SetTitle("Counts per 12.5 MeV");

        AN_Q1SHisto->SetStats(0);
        AN_Q1SHisto->SetMarkerColor(MarkerColor_nS[1]);
        AN_Q1SHisto->SetMarkerStyle(MarkerStyle_nS[1]);
        AN_Q1SHisto->SetMarkerSize(MarkerSize_nS[1]);
        AN_Q1SHisto->Draw("E");
        AN_Q2SHisto->SetStats(0);
        AN_Q2SHisto->SetMarkerColor(MarkerColor_nS[2]);
        AN_Q2SHisto->SetMarkerStyle(MarkerStyle_nS[2]);
        AN_Q2SHisto->SetMarkerSize(MarkerSize_nS[2]);
        if(nState>1.5) AN_Q2SHisto->Draw("same,E");
        AN_Q3SHisto->SetStats(0);
        AN_Q3SHisto->SetMarkerColor(MarkerColor_nS[3]);
        AN_Q3SHisto->SetMarkerStyle(MarkerStyle_nS[3]);
        AN_Q3SHisto->SetMarkerSize(MarkerSize_nS[3]);
        if(nState>2.5) AN_Q3SHisto->Draw("same,E");


        bool addStandardQ=true;
        if(addStandardQ){
        	  cout<<"addStandardQ"<<endl;

              int ColorQ1P=632;
              int ColorQ2P=600;
              int ColorQ3P=416;
      const int nQLinesStandard=14;
  	  TLine* QLine[nQLinesStandard];

  	  double QvalForQlines[nQLinesStandard]={
  			  Q_MEAS_chib0_1P_1S,
  			  Q_MEAS_chib1_1P_1S,
  		      Q_MEAS_chib2_1P_1S,
  		      Q_MEAS_chib0_2P_1S,
  		      Q_MEAS_chib1_2P_1S,
  		      Q_MEAS_chib2_2P_1S,
  		      Q_MEAS_chib1_2P_2S,
  		      Q_MEAS_chib2_2P_2S,
  		      Q_MEAS_chib1_3P_1S,
  		      Q_MEAS_chib2_3P_1S,
  		      Q_MEAS_chib1_3P_2S,
  		      Q_MEAS_chib2_3P_2S,
  		      Q_MEAS_chib1_3P_3S,
  		      Q_MEAS_chib2_3P_3S
  	  };

  	  int ColorForQlines[nQLinesStandard]={ColorQ1P,ColorQ1P,ColorQ1P, ColorQ2P,ColorQ2P,ColorQ2P,ColorQ2P,ColorQ2P, ColorQ3P,ColorQ3P,ColorQ3P,ColorQ3P,ColorQ3P,ColorQ3P};

  	  int StyleForQlines[nQLinesStandard]={2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

  	  for(int iQLine=0;iQLine<nQLinesStandard;iQLine++){
  	  QLine[iQLine] = new TLine( 0,0,0,0 );
      QLine[iQLine]->SetLineWidth( 1 );
  	  QLine[iQLine]->SetLineStyle(StyleForQlines[iQLine]);
  	  QLine[iQLine]->SetLineColor(ColorForQlines[iQLine]);
  	  QLine[iQLine]->SetY1(0);
  	  QLine[iQLine]->SetY2(AN_Q1SHisto->GetMaximum());
  	  QLine[iQLine]->SetX1(QvalForQlines[iQLine]);
  	  QLine[iQLine]->SetX2(QvalForQlines[iQLine]);
  	  QLine[iQLine]->Draw("same");
  	  cout<<QvalForQlines[iQLine]<<endl;
  	  }

        }

        bool addReflectedQ=true;
        if(addReflectedQ){
            int ColorQ1P=627;
            int ColorQ2P=595;
            int ColorQ3P=411;
      	  cout<<"addReflectedQ"<<endl;
      const int nQLinesReflected=8;
  	  TLine* QLine[nQLinesReflected];

  	  double QvalForQlines[nQLinesReflected]={
  		  	Q_MEAS_chib1_2P_2Sin1S,
  		  	Q_MEAS_chib2_2P_2Sin1S,
  		  	Q_MEAS_chib1_3P_2Sin1S,
  		  	Q_MEAS_chib2_3P_2Sin1S,
  		  	Q_MEAS_chib1_3P_3Sin1S,
  		  	Q_MEAS_chib2_3P_3Sin1S,
  		  	Q_MEAS_chib1_3P_3Sin2S,
  		  	Q_MEAS_chib2_3P_3Sin2S
  	  };

  	  int ColorForQlines[nQLinesReflected]={ColorQ2P,ColorQ2P, ColorQ3P,ColorQ3P,ColorQ3P,ColorQ3P,ColorQ3P,ColorQ3P};

  	  int StyleForQlines[nQLinesReflected]={3, 3, 3, 3, 3, 3, 3, 3};

  	  for(int iQLine=0;iQLine<nQLinesReflected;iQLine++){
  	  QLine[iQLine] = new TLine( 0,0,0,0 );
      QLine[iQLine]->SetLineWidth( 1 );
  	  QLine[iQLine]->SetLineStyle(StyleForQlines[iQLine]);
  	  QLine[iQLine]->SetLineColor(ColorForQlines[iQLine]);
  	  QLine[iQLine]->SetY1(0);
  	  QLine[iQLine]->SetY2(AN_Q1SHisto->GetMaximum());
  	  QLine[iQLine]->SetX1(QvalForQlines[iQLine]);
  	  QLine[iQLine]->SetX2(QvalForQlines[iQLine]);
  	  QLine[iQLine]->Draw("same");
  	  cout<<QvalForQlines[iQLine]<<endl;
  	  }

        }


        bool addMirroredQ=true;
        if(addMirroredQ){
            int ColorQ1P=635;
            int ColorQ2P=603;
            int ColorQ3P=419;
        	  cout<<"addMirroredQ"<<endl;
      const int nQLinesMirrored=8;
  	  TLine* QLine[nQLinesMirrored];

  	  double QvalForQlines[nQLinesMirrored]={
  		  	Q_MEAS_chib1_1P_1Sfrom3S,
  		  	Q_MEAS_chib2_1P_1Sfrom3S,
  		  	Q_MEAS_chib1_1P_1Sfrom2S,
  		  	Q_MEAS_chib2_1P_1Sfrom2S,
  		  	Q_MEAS_chib1_2P_1Sfrom3S,
  		  	Q_MEAS_chib2_2P_1Sfrom3S,
  		  	Q_MEAS_chib1_2P_2Sfrom3S,
  		  	Q_MEAS_chib2_2P_2Sfrom3S
  	  };

  	  int ColorForQlines[nQLinesMirrored]={ColorQ1P,ColorQ1P,ColorQ1P,ColorQ1P, ColorQ2P,ColorQ2P,ColorQ2P,ColorQ2P};

  	  int StyleForQlines[nQLinesMirrored]={10,10,10,10,10,10,10,10};

  	  for(int iQLine=0;iQLine<nQLinesMirrored;iQLine++){
  	  QLine[iQLine] = new TLine( 0,0,0,0 );
      QLine[iQLine]->SetLineWidth( 1 );
  	  QLine[iQLine]->SetLineStyle(StyleForQlines[iQLine]);
  	  QLine[iQLine]->SetLineColor(ColorForQlines[iQLine]);
  	  QLine[iQLine]->SetY1(0);
  	  QLine[iQLine]->SetY2(AN_Q1SHisto->GetMaximum());
  	  QLine[iQLine]->SetX1(QvalForQlines[iQLine]);
  	  QLine[iQLine]->SetX2(QvalForQlines[iQLine]);
  	  QLine[iQLine]->Draw("same");
  	  cout<<QvalForQlines[iQLine]<<endl;
  	  }

        }


        TLegend* plotcompLegend=new TLegend(0.8,0.725,0.875,0.925);
        plotcompLegend->SetFillColor(0);
        plotcompLegend->SetTextFont(72);
        plotcompLegend->SetTextSize(0.04);
        plotcompLegend->SetBorderSize(0);
        char complegendentry[200];
        sprintf(complegendentry,"#Upsilon(1S)");
        plotcompLegend->AddEntry(AN_Q1SHisto,complegendentry,"elp");
        sprintf(complegendentry,"#Upsilon(2S)");
        if(nState>1.5) plotcompLegend->AddEntry(AN_Q2SHisto,complegendentry,"elp");
        sprintf(complegendentry,"#Upsilon(3S)");
        if(nState>2.5) plotcompLegend->AddEntry(AN_Q3SHisto,complegendentry,"elp");
        if(nState>1.5) plotcompLegend->Draw();

        PlotCanvas->Modified();
        PlotCanvas->SetLogy(false);
        sprintf(saveMass,"Figures/%s/AN_QnSHisto.pdf",FitID);
        PlotCanvas->SaveAs(saveMass);
    }

    delete NonEquChibCanvas1S2S3S;
	delete NonEquframe1S2S3S_1S;
	delete NonEquframe1S2S3S_2S;
	delete NonEquframe1S2S3S_3S;

    delete NonEquChibCanvas1S;
    delete data1SNonEquHist;
	delete NonEquframe1S;

    delete NonEquChibCanvas2S;
    delete data2SNonEquHist;
	delete NonEquframe2S;

    delete NonEquChibCanvas3S;
    delete data3SNonEquHist;
	delete NonEquframe3S;





	    }











    bool PlotPull=false;
	if(PlotPull){


		  TH1 * ModelHist =  modelPdf1S.createHistogram("invm1S",FrameBins1S) ;
		  TH1F  *DataHist = new TH1F("DataHist","",FrameBins1S,invm_min1S,invm_max1S);
		  tree1S->Draw("invm1S>>DataHist");

		  cout<<"ModelHist entries "<<ModelHist->GetSumOfWeights()<<endl;
		  ModelHist->Scale(totalEventsInFit1S/ModelHist->GetSumOfWeights());
		  cout<<"ModelHist entries "<<ModelHist->GetSumOfWeights()<<endl;

		  double deltaPull[FrameBins1S];
		  double errdeltaPull[FrameBins1S];
		  double FitValue;
		  double invmCenter[FrameBins1S];
		  double errinvmCenter[FrameBins1S];

		  for(int BinRun = 1; BinRun < FrameBins1S+1; ++BinRun) {

			  FitValue=ModelHist->GetBinContent(BinRun);

			  if(DataHist->GetBinError(BinRun)==0) deltaPull[BinRun-1]=100000;
			  else deltaPull[BinRun-1]=((DataHist->GetBinContent(BinRun)-FitValue)/DataHist->GetBinError(BinRun));
//			  cout<<deltaPull[BinRun-1]<<endl;

			  errdeltaPull[BinRun-1]=0;
			  invmCenter[BinRun-1]=DataHist->GetBinCenter(BinRun);
			  errinvmCenter[BinRun-1]=0;

		  }


			 TCanvas *PullCanvas = new TCanvas("PullCanvas","",1330,800);
			 PullCanvas->SetFillColor(kWhite);
			 PullCanvas->cd(1);
//			 PullCanvas->SetLeftMargin(0.2);
			 gPad->SetFillColor(kWhite);
//			 gPad->SetRightMargin(0.3);/////


					TH1F *PullHisto = new TH1F;
					PullHisto = PullCanvas->DrawFrame(invm_min1S,-4,invm_max1S,4);
					PullHisto->SetXTitle(invm1S.GetTitle());
					PullHisto->SetYTitle("( Data - Model ) / #sigma_{Data}");
					PullHisto->GetYaxis()->SetTitleOffset(1.1);

					  TLine* refLine[5];
					  double refLinePosition[5]={-2.,-1.,0.,1.,2.};
					  double refLineWidth[5]={0.5,1.,2.,1.,0.5};
					  for(int iLine=0;iLine<5;iLine++){
					  refLine[iLine] = new TLine( invm_min1S, refLinePosition[iLine], invm_max1S, refLinePosition[iLine] );
					  refLine[iLine]->SetLineWidth( refLineWidth[iLine] );
					  refLine[iLine]->SetLineStyle( 2 );
					  refLine[iLine]->SetLineColor( kBlack );
					  refLine[iLine]->Draw( "same" );
					  }

			  TGraphErrors *PullGraph = new TGraphErrors(FrameBins1S,invmCenter,deltaPull,errinvmCenter,errdeltaPull);
			  PullGraph->SetMarkerColor(kGreen+2);
			  PullGraph->SetMarkerStyle(20);
			  PullGraph->SetTitle(0);
			  PullGraph->Draw("P");

		    sprintf(saveName,"Figures/%s/InvMass_Pull_Y1S_%s.pdf",FitID,cutName_);
		    PullCanvas->SaveAs(saveName);
   		    PullCanvas->Close();

}













    cout<<"deleting pointers"<<endl;

    delete dataAfterBasicCuts;
    delete data_;
    delete nll1S_1P2P;
    delete nll1S_1P2P_;
    delete nll1S_BB;
    delete nll2S_2P;
    delete nll2S_BB;
    delete nll1S;
    delete nll2S;
    if(BkgMixerBool&&!useExistingMixFile){
    delete data1SExceptGammaPtCut;
    delete data2SExceptGammaPtCut;
    delete data3SExceptGammaPtCut;
    delete dataAfterBasicCutsExceptGammaPtCut;
    }

/*    delete data_______;
    delete data__;
    delete data___;
    delete data____;
    delete data_____;
    delete data______;
    delete data______vtx;
*/
    delete data1S;
    delete data2S;
    delete data3S;
    delete dataSB1S;
    delete dataSB2S;

    delete minuit;
    delete PhotonEnergyCanvas;
    delete PhotonPtCanvas;
    delete fQ;

	delete frame1S;
	delete frame2S;
	delete frame1S2S_;

	delete ChibCanvas1S;
    delete ChibCanvas2S;
    delete ChibCanvas1S2S;

    delete RooBkgToySet;
    delete RooBkgToyHist1S;
    delete RooBkgToyHist2S;

    cout<<"deleted pointers"<<endl;
    }

    sprintf(filename_,"Figures/%s/CutHistos_Y%dS.root",FitID,nState);
	TFile *CutHistos_file = new TFile(filename_,"RECREATE");

	char XTitle[200];
	if(nYsigCut  )sprintf(XTitle,"Y(1S) mass window, n#sigma_{m_{#mu#mu}}(|y|^{Y(1S)}) ");
	if(vtxProbCut)sprintf(XTitle,"Dimuon vertex probability cut");
	if(ctSigCut  )sprintf(XTitle,"c#tau significance cut, n#sigma");
	if(gammaptCut)sprintf(XTitle,"#gamma - p_{T} cut, GeV");
	if(gammaptDCut)sprintf(XTitle,"Param. D of Q-dep. #gamma - p_{T} cut, GeV");
	if(ctCut     )sprintf(XTitle,"c#tau cut, cm");
	if(RConvUpperCut  )sprintf(XTitle,"RConvMax cut");
	if(RConvLowerCut  )sprintf(XTitle,"RConvMin cut");
	if(pTCut  )sprintf(XTitle,"lower p_{T}^{Y} cut, GeV");
	if(vtxChi2ProbGammaCut  )sprintf(XTitle,"vertexChi2ProbGamma cut");
	if(vtxChi2ProbGammaLogCut  )sprintf(XTitle,"log10(vertexChi2ProbGamma) cut");
	if(Pi0Cut  )sprintf(XTitle,"rejected #pi^{0} mass window: |0.134 - m_{#gamma#gamma}| < x GeV");
	if(dzSigCut  )sprintf(XTitle,"dz/#sigma_{dz} cut");
	if(dzCut  )sprintf(XTitle,"dz cut");

	if(n_Cut>1.5){

	double yOffset=1.3;

    TCanvas* SummaryCanvas = new TCanvas("SummaryCanvas","SummaryCanvas",1600, 4224);
    SummaryCanvas->Divide(2,11);
    SummaryCanvas->SetFillColor(kWhite);
    SummaryCanvas->cd(1); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_1Psig->SetStats(0);hCut_1Psig->GetYaxis()->SetTitle("#chi_{b}(1P)_{sign.}");	hCut_1Psig->GetYaxis()->SetTitleOffset(yOffset); hCut_1Psig->GetXaxis()->SetTitle(XTitle);	hCut_1Psig->Draw();
    SummaryCanvas->cd(2); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_2Psig->SetStats(0);hCut_2Psig->GetYaxis()->SetTitle("#chi_{b}(2P)_{sign.}");	hCut_2Psig->GetYaxis()->SetTitleOffset(yOffset); hCut_2Psig->GetXaxis()->SetTitle(XTitle);	hCut_2Psig->Draw();

    SummaryCanvas->cd(3); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_1PsigLL->SetStats(0);hCut_1PsigLL->GetYaxis()->SetTitle("#chi_{b}(1P)_{sign.LL}");		hCut_1PsigLL->GetYaxis()->SetTitleOffset(yOffset); hCut_1PsigLL->GetXaxis()->SetTitle(XTitle);	hCut_1PsigLL->Draw();
    SummaryCanvas->cd(4); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_2PsigLL->SetStats(0);hCut_2PsigLL->GetYaxis()->SetTitle("#chi_{b}(2P)_{sign.LL}");		hCut_2PsigLL->GetYaxis()->SetTitleOffset(yOffset); hCut_2PsigLL->GetXaxis()->SetTitle(XTitle);	hCut_2PsigLL->Draw();


    SummaryCanvas->cd(5); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_1PSB->SetStats(0);hCut_1PSB->GetYaxis()->SetTitle("#chi_{b}(1P)_{S/B}");		hCut_1PSB->GetYaxis()->SetTitleOffset(yOffset); hCut_1PSB->GetXaxis()->SetTitle(XTitle);	hCut_1PSB->Draw();
    SummaryCanvas->cd(6); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_2PSB->SetStats(0);hCut_2PSB->GetYaxis()->SetTitle("#chi_{b}(2P)_{S/B}");		hCut_2PSB->GetYaxis()->SetTitleOffset(yOffset); hCut_2PSB->GetXaxis()->SetTitle(XTitle);	hCut_2PSB->Draw();

    SummaryCanvas->cd(7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_1PS->SetStats(0);hCut_1PS->GetYaxis()->SetTitle("#chi_{b}(1P)_{S}");			hCut_1PS->GetYaxis()->SetTitleOffset(yOffset); hCut_1PS->GetXaxis()->SetTitle(XTitle);	hCut_1PS->Draw();
    SummaryCanvas->cd(8); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_2PS->SetStats(0);hCut_2PS->GetYaxis()->SetTitle("#chi_{b}(2P)_{S}");			hCut_2PS->GetYaxis()->SetTitleOffset(yOffset); hCut_2PS->GetXaxis()->SetTitle(XTitle);	hCut_2PS->Draw();
    SummaryCanvas->cd(9); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_1PB->SetStats(0);hCut_1PB->GetYaxis()->SetTitle("#chi_{b}(1P)_{B}");			hCut_1PB->GetYaxis()->SetTitleOffset(yOffset); hCut_1PB->GetXaxis()->SetTitle(XTitle);	hCut_1PB->Draw();
    SummaryCanvas->cd(10); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_2PB->SetStats(0);hCut_2PB->GetYaxis()->SetTitle("#chi_{b}(2P)_{B}");			hCut_2PB->GetYaxis()->SetTitleOffset(yOffset); hCut_2PB->GetXaxis()->SetTitle(XTitle);	hCut_2PB->Draw();

    SummaryCanvas->cd(11); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_sig->SetStats(0);hCut_sig->GetYaxis()->SetTitle("(#chi_{b}(1P)_{sign.}+#chi_{b}(2P)_{sign.})/2");			hCut_sig->GetYaxis()->SetTitleOffset(yOffset); hCut_sig->GetXaxis()->SetTitle(XTitle);	hCut_sig->Draw();
    SummaryCanvas->cd(12); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_2PSover1PS->SetStats(0);hCut_2PSover1PS->GetYaxis()->SetTitle("#chi_{b}(2P)_{S}/#chi_{b}(1P)_{S}");			hCut_2PSover1PS->GetYaxis()->SetTitleOffset(yOffset); hCut_2PSover1PS->GetXaxis()->SetTitle(XTitle);	hCut_2PSover1PS->Draw();

    SummaryCanvas->cd(13); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3Psig->SetStats(0);hCut_3Psig->GetYaxis()->SetTitle("#chi_{b}(3P)_{sign.}");	hCut_3Psig->GetYaxis()->SetTitleOffset(yOffset); hCut_3Psig->GetXaxis()->SetTitle(XTitle);	hCut_3Psig->Draw();
    SummaryCanvas->cd(14); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3Psig_1S->SetStats(0);hCut_3Psig_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{sign.} (1S)");	hCut_3Psig_1S->GetYaxis()->SetTitleOffset(yOffset); hCut_3Psig_1S->GetXaxis()->SetTitle(XTitle);	hCut_3Psig_1S->Draw();
    SummaryCanvas->cd(15); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PSB->SetStats(0);hCut_3PSB->GetYaxis()->SetTitle("#chi_{b}(3P)_{S/B}");		hCut_3PSB->GetYaxis()->SetTitleOffset(yOffset); hCut_3PSB->GetXaxis()->SetTitle(XTitle);	hCut_3PSB->Draw();
    SummaryCanvas->cd(16); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PSB_1S->SetStats(0);hCut_3PSB_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{S/B} (1S)");		hCut_3PSB_1S->GetYaxis()->SetTitleOffset(yOffset); hCut_3PSB_1S->GetXaxis()->SetTitle(XTitle);	hCut_3PSB_1S->Draw();

    SummaryCanvas->cd(17); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PS->SetStats(0);hCut_3PS->GetYaxis()->SetTitle("#chi_{b}(3P)_{S}");			hCut_3PS->GetYaxis()->SetTitleOffset(yOffset); hCut_3PS->GetXaxis()->SetTitle(XTitle);	hCut_3PS->Draw();
    SummaryCanvas->cd(18); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PS_1S->SetStats(0);hCut_3PS_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{S} (1S)");			hCut_3PS_1S->GetYaxis()->SetTitleOffset(yOffset); hCut_3PS_1S->GetXaxis()->SetTitle(XTitle);	hCut_3PS_1S->Draw();
    SummaryCanvas->cd(19); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PB->SetStats(0);hCut_3PB->GetYaxis()->SetTitle("#chi_{b}(3P)_{B}");			hCut_3PB->GetYaxis()->SetTitleOffset(yOffset); hCut_3PB->GetXaxis()->SetTitle(XTitle);	hCut_3PB->Draw();
    SummaryCanvas->cd(20); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PB_1S->SetStats(0);hCut_3PB_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{B} (1S)");			hCut_3PB_1S->GetYaxis()->SetTitleOffset(yOffset); hCut_3PB_1S->GetXaxis()->SetTitle(XTitle);	hCut_3PB_1S->Draw();

    SummaryCanvas->cd(21); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_3PsigALL->SetStats(0);hCut_3PsigALL->GetYaxis()->SetTitle("#chi_{b}(3P)_{sign.comb.}");			hCut_3PsigALL->GetYaxis()->SetTitleOffset(yOffset); hCut_3PsigALL->GetXaxis()->SetTitle(XTitle);	hCut_3PsigALL->Draw();
    SummaryCanvas->cd(22); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hCut_PhotonEnergyScale->SetStats(0);hCut_PhotonEnergyScale->GetYaxis()->SetTitle("Photon Energy Scale");			hCut_PhotonEnergyScale->GetYaxis()->SetTitleOffset(yOffset); hCut_PhotonEnergyScale->GetXaxis()->SetTitle(XTitle);	hCut_PhotonEnergyScale->Draw();

    SummaryCanvas->Modified();
    char saveName[200];
    sprintf(saveName,"Figures/%s/SummaryCanvas_Y%dS.pdf",FitID,nState);
    SummaryCanvas->SaveAs(saveName);

    TCanvas* IndividualCanvas = new TCanvas("IndividualCanvas","IndividualCanvas",1600, 800);
    IndividualCanvas->SetFillColor(kWhite);
    hCut_1Psig->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1Psig->SetStats(0);hCut_1Psig->GetYaxis()->SetTitle("#chi_{b}(1P)_{sign.}");	hCut_1Psig->GetXaxis()->SetTitle(XTitle);	hCut_1Psig->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1Psig.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2Psig->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Psig->SetStats(0);hCut_2Psig->GetYaxis()->SetTitle("#chi_{b}(2P)_{sign.}");	hCut_2Psig->GetXaxis()->SetTitle(XTitle);	hCut_2Psig->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2Psig.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1PS->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1PS->SetStats(0);hCut_1PS->GetYaxis()->SetTitle("#chi_{b}(1P)_{S}");	hCut_1PS->GetXaxis()->SetTitle(XTitle);	hCut_1PS->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1PS.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2PS->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2PS->SetStats(0);hCut_2PS->GetYaxis()->SetTitle("#chi_{b}(2P)_{S}");	hCut_2PS->GetXaxis()->SetTitle(XTitle);	hCut_2PS->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PS.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1PB->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1PB->SetStats(0);hCut_1PB->GetYaxis()->SetTitle("#chi_{b}(1P)_{B}");	hCut_1PB->GetXaxis()->SetTitle(XTitle);	hCut_1PB->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1PB.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2PB->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2PB->SetStats(0);hCut_2PB->GetYaxis()->SetTitle("#chi_{b}(2P)_{B}");	hCut_2PB->GetXaxis()->SetTitle(XTitle);	hCut_2PB->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PB.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1PSB->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1PSB->SetStats(0);hCut_1PSB->GetYaxis()->SetTitle("#chi_{b}(1P)_{S/B}");	hCut_1PSB->GetXaxis()->SetTitle(XTitle);	hCut_1PSB->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1PSB.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2PSB->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2PSB->SetStats(0);hCut_2PSB->GetYaxis()->SetTitle("#chi_{b}(2P)_{S/B}");	hCut_2PSB->GetXaxis()->SetTitle(XTitle);	hCut_2PSB->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PSB.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_sig->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_sig->SetStats(0);hCut_sig->GetYaxis()->SetTitle("(#chi_{b}(1P)_{sign.} + #chi_{b}(2P)_{sign.})/2");	hCut_sig->GetXaxis()->SetTitle(XTitle);	hCut_sig->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/sig.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_covQual1->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_covQual1->SetStats(0);hCut_covQual1->GetYaxis()->SetTitle("covQual fit 1");	hCut_covQual1->GetXaxis()->SetTitle(XTitle);	hCut_covQual1->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/covQual1.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1Pmean->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1Pmean->SetStats(0);hCut_1Pmean->GetYaxis()->SetTitle("#chi_{b}(1P)_{mean}");	hCut_1Pmean->GetXaxis()->SetTitle(XTitle);	hCut_1Pmean->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1Pmean.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1Pwidth->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1Pwidth->SetStats(0);hCut_1Pwidth->GetYaxis()->SetTitle("#chi_{b}(1P)_{width}");	hCut_1Pwidth->GetXaxis()->SetTitle(XTitle);	hCut_1Pwidth->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1Pwidth.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1Pnevt->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1Pnevt->SetStats(0);hCut_1Pnevt->GetYaxis()->SetTitle("#chi_{b}(1P)_{nevt}");	hCut_1Pnevt->GetXaxis()->SetTitle(XTitle);	hCut_1Pnevt->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1Pnevt.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2Pmean->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Pmean->SetStats(0);hCut_2Pmean->GetYaxis()->SetTitle("#chi_{b}(2P)_{mean}");	hCut_2Pmean->GetXaxis()->SetTitle(XTitle);	hCut_2Pmean->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2Pmean.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2Pwidth->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Pwidth->SetStats(0);hCut_2Pwidth->GetYaxis()->SetTitle("#chi_{b}(2P)_{width}");	hCut_2Pwidth->GetXaxis()->SetTitle(XTitle);	hCut_2Pwidth->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2Pwidth.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);


    hCut_3Pwidth->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3Pwidth->SetStats(0);hCut_3Pwidth->GetYaxis()->SetTitle("#chi_{b}(3P)_{width}, 2S sample");	hCut_3Pwidth->GetXaxis()->SetTitle(XTitle);	hCut_3Pwidth->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3Pwidth_2S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3Pwidth_1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3Pwidth_1S->SetStats(0);hCut_3Pwidth_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{width}, 1S sample");	hCut_3Pwidth_1S->GetXaxis()->SetTitle(XTitle);	hCut_3Pwidth_1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3Pwidth_1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2Pwidth_2S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Pwidth_2S->SetStats(0);hCut_2Pwidth_2S->GetYaxis()->SetTitle("#chi_{b}(2P)_{width}, 2S sample");	hCut_2Pwidth_2S->GetXaxis()->SetTitle(XTitle);	hCut_2Pwidth_2S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2Pwidth_2S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2Pnevt->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Pnevt->SetStats(0);hCut_2Pnevt->GetYaxis()->SetTitle("#chi_{b}(2P)_{nevt}");	hCut_2Pnevt->GetXaxis()->SetTitle(XTitle);	hCut_2Pnevt->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2Pnevt.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_CBa->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_CBa->SetStats(0);hCut_CBa->GetYaxis()->SetTitle("#alpha");	hCut_CBa->GetXaxis()->SetTitle(XTitle);	hCut_CBa->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/CBa.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_CBn->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_CBn->SetStats(0);hCut_CBn->GetYaxis()->SetTitle("n");	hCut_CBn->GetXaxis()->SetTitle(XTitle);	hCut_CBn->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/CBn.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_BGnevt1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_BGnevt1S->SetStats(0);hCut_BGnevt1S->GetYaxis()->SetTitle("BG_{nevt1S}");	hCut_BGnevt1S->GetXaxis()->SetTitle(XTitle);	hCut_BGnevt1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/BGnevt1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_BGnevt2S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_BGnevt2S->SetStats(0);hCut_BGnevt2S->GetYaxis()->SetTitle("BG_{nevt2S}");	hCut_BGnevt2S->GetXaxis()->SetTitle(XTitle);	hCut_BGnevt2S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/BGnevt2S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_Ntotal1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_Ntotal1S->SetStats(0);hCut_Ntotal1S->GetYaxis()->SetTitle("N_{total1S}");	hCut_Ntotal1S->GetXaxis()->SetTitle(XTitle);	hCut_Ntotal1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/Ntotal1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_Ntotal2S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_Ntotal2S->SetStats(0);hCut_Ntotal2S->GetYaxis()->SetTitle("N_{total2S}");	hCut_Ntotal2S->GetXaxis()->SetTitle(XTitle);	hCut_Ntotal2S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/Ntotal2S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2PSover1PS->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2PSover1PS->SetStats(0);hCut_2PSover1PS->GetYaxis()->SetTitle("#chi_{b}(2P)_{S} / #chi_{b}(1P)_{S}");	hCut_2PSover1PS->GetXaxis()->SetTitle(XTitle);	hCut_2PSover1PS->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PSover1PS.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);



    hCut_3Psig->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3Psig->SetStats(0);hCut_3Psig->GetYaxis()->SetTitle("#chi_{b}(3P)_{sign.}, 2S sample");	hCut_3Psig->GetXaxis()->SetTitle(XTitle);	hCut_3Psig->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3Psig.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3Psig_1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3Psig_1S->SetStats(0);hCut_3Psig_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{sign.}, 1S sample");	hCut_3Psig_1S->GetXaxis()->SetTitle(XTitle);	hCut_3Psig_1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3Psig_1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PsigALL->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PsigALL->SetStats(0);hCut_3PsigALL->GetYaxis()->SetTitle("#chi_{b}(3P)_{sign.}, 1S, 2S sample combined");	hCut_3PsigALL->GetXaxis()->SetTitle(XTitle);	hCut_3PsigALL->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PsigALL.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);


    hCut_3PSB->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PSB->SetStats(0);hCut_3PSB->GetYaxis()->SetTitle("#chi_{b}(3P)_{S/B}, 2S sample");	hCut_3PSB->GetXaxis()->SetTitle(XTitle);	hCut_3PSB->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PSB.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PSB_1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PSB_1S->SetStats(0);hCut_3PSB_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{S/B}, 1S sample");	hCut_3PSB_1S->GetXaxis()->SetTitle(XTitle);	hCut_3PSB_1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PSB_1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);


    hCut_3PS->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PS->SetStats(0);hCut_3PS->GetYaxis()->SetTitle("#chi_{b}(3P)_{S}, 2S sample");	hCut_3PS->GetXaxis()->SetTitle(XTitle);	hCut_3PS->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PS.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PS_1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PS_1S->SetStats(0);hCut_3PS_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{S}, 1S sample");	hCut_3PS_1S->GetXaxis()->SetTitle(XTitle);	hCut_3PS_1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PS_1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PB->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PB->SetStats(0);hCut_3PB->GetYaxis()->SetTitle("#chi_{b}(3P)_{B}, 2S sample");	hCut_3PB->GetXaxis()->SetTitle(XTitle);	hCut_3PB->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PB.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PB_1S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PB_1S->SetStats(0);hCut_3PB_1S->GetYaxis()->SetTitle("#chi_{b}(3P)_{B}, 1S sample");	hCut_3PB_1S->GetXaxis()->SetTitle(XTitle);	hCut_3PB_1S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PB_1S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3Pmass->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3Pmass->SetStats(0);hCut_3Pmass->GetYaxis()->SetTitle("#chi_{b}(3P) mass");	hCut_3Pmass->GetXaxis()->SetTitle(XTitle);	hCut_3Pmass->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3Pmass.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_PhotonEnergyScale->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonEnergyScale->SetStats(0);hCut_PhotonEnergyScale->GetYaxis()->SetTitle("Photon Energy Scale");	hCut_PhotonEnergyScale->GetXaxis()->SetTitle(XTitle);	hCut_PhotonEnergyScale->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/PhotonEnergyScale.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_PhotonEnergyScale2P->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonEnergyScale2P->SetStats(0);hCut_PhotonEnergyScale2P->GetYaxis()->SetTitle("Photon Energy Scale 2P");	hCut_PhotonEnergyScale2P->GetXaxis()->SetTitle(XTitle);	hCut_PhotonEnergyScale2P->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/PhotonEnergyScale2P.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2PsigLL->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2PsigLL->SetStats(0);hCut_2PsigLL->GetYaxis()->SetTitle("#chi_{b}(2P)_{sign.LL}");	hCut_2PsigLL->GetXaxis()->SetTitle(XTitle);	hCut_2PsigLL->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PsigLL.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1PsigLL->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1PsigLL->SetStats(0);hCut_1PsigLL->GetYaxis()->SetTitle("#chi_{b}(1P)_{sign.LL}");	hCut_1PsigLL->GetXaxis()->SetTitle(XTitle);	hCut_1PsigLL->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1PsigLL.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_sigLL->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_sigLL->SetStats(0);hCut_sigLL->GetYaxis()->SetTitle("comb. LLratio-significance of #chi_{b}(1P) and #chi_{b}(2P)");	hCut_sigLL->GetXaxis()->SetTitle(XTitle);	hCut_sigLL->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/CombsigLL.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_sigLLhalf->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_sigLLhalf->SetStats(0);hCut_sigLLhalf->GetYaxis()->SetTitle("(#chi_{b}(1P)_{sign.LL} + #chi_{b}(2P)_{sign.LL})/2");	hCut_sigLLhalf->GetXaxis()->SetTitle(XTitle);	hCut_sigLLhalf->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/sigLL.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PsigLLhalf->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PsigLLhalf->SetStats(0);hCut_3PsigLLhalf->GetYaxis()->SetTitle("( #chi_{b}(3P)_{sign.LL} (1S) + #chi_{b}(3P)_{sign.LL} (2S) )/2");	hCut_3PsigLLhalf->GetXaxis()->SetTitle(XTitle);	hCut_3PsigLLhalf->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PsigBothHalf.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);


    hCut_2Psig_2S->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Psig_2S->SetStats(0);hCut_2Psig_2S->GetYaxis()->SetTitle("#chi_{b}(2P)_{sign.LL}, 2S sample");	hCut_2Psig_2S->GetXaxis()->SetTitle(XTitle);	hCut_2Psig_2S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PsigLL_2S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_2Pnevt_2S->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2Pnevt_2S->SetStats(0);hCut_2Pnevt_2S->GetYaxis()->SetTitle("#chi_{b}(2P)_{nevt}, 2S sample");	hCut_2Pnevt_2S->GetXaxis()->SetTitle(XTitle);	hCut_2Pnevt_2S->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2Pnevt_2S.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_PhotonSigmaScale->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonSigmaScale->SetStats(0);hCut_PhotonSigmaScale->GetYaxis()->SetTitle("#sigma/Q");	hCut_PhotonSigmaScale->GetXaxis()->SetTitle(XTitle);	hCut_PhotonSigmaScale->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/PhotonSigmaScale.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_PhotonSigmaScale2P->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonSigmaScale2P->SetStats(0);hCut_PhotonSigmaScale2P->GetYaxis()->SetTitle("#sigma/Q, 2P");	hCut_PhotonSigmaScale2P->GetXaxis()->SetTitle(XTitle);	hCut_PhotonSigmaScale2P->Draw("E1");
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/PhotonSigmaScale2P.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_3PMerr->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3PMerr->SetStats(0);hCut_3PMerr->GetYaxis()->SetTitle("Stat uncertainty on 3P mass meas.");	hCut_3PMerr->GetXaxis()->SetTitle(XTitle);	hCut_3PMerr->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3PmassErr.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

	}

	hCut_1Psig->Write();
	hCut_1PSB->Write();
	hCut_1PS->Write();
	hCut_1PB->Write();
	hCut_2Psig->Write();
	hCut_2PSB->Write();
	hCut_2PS->Write();
	hCut_2PB->Write();
	hCut_sig->Write();
	hCut_1Pmean ->Write();
	hCut_1Pwidth->Write();
	hCut_1Pnevt ->Write();
	hCut_2Pmean ->Write();
	hCut_2Pwidth->Write();
	hCut_2Pnevt ->Write();
	hCut_CBa    ->Write();
	hCut_CBn    ->Write();
	hCut_2PSover1PS->Write();

	hCut_3Pwidth->Write();
	hCut_3Pwidth_1S->Write();
	hCut_2Pwidth_2S->Write();

	hCut_covQual1->Write();
	hCut_BGnevt1S ->Write();
	hCut_BGnevt2S ->Write();
	hCut_Ntotal1S->Write();
	hCut_Ntotal2S->Write();

	hCut_3Psig->Write();
	hCut_3Psig_1S->Write();
	hCut_3PsigALL->Write();

	hCut_3PSB->Write();
	hCut_3PS->Write();
	hCut_3PB->Write();
	hCut_3Pmass->Write();

	hCut_3PSB_1S->Write();
	hCut_3PS_1S->Write();
	hCut_3PB_1S->Write();
	hCut_PhotonEnergyScale->Write();
	hCut_PhotonEnergyScale2P->Write();
	hCut_PhotonSigmaScale->Write();
	hCut_PhotonSigmaScale2P->Write();
	hCut_1PsigLL->Write();
	hCut_1PsigLL->Write();

	hCut_2Pnevt_2S->Write();
	hCut_sigLL->Write();
	hCut_2Psig_2S->Write();
	hCut_sigLLhalf->Write();
	hCut_3PsigLLhalf->Write();

	hCut_3PMerr->Write();

	CutHistos_file->Write();
    CutHistos_file->Close();

	return 0;
}


