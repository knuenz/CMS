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


double NormalizedIntegral(RooAbsPdf * function, RooRealVar& integrationVar, double lowerLimit, double upperLimit)
{
using namespace RooFit;
integrationVar.setRange("integralRange", lowerLimit, upperLimit);
RooAbsReal* integral = (*function).createIntegral(integrationVar,NormSet(integrationVar), Range("integralRange"));
double normalizedIntegralValue = integral->getVal();
return normalizedIntegralValue;
}
int main(int argc, char** argv) {

	Char_t *FitID = "Default"; //Storage Directory
	Char_t *cutName = "Default"; //Code Directory
  	Char_t *fileName = "Default";
  	int nState=1000;
  	char cutName_[500];
    int nToy=5000000;

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

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("FitID") != std::string::npos) {char* FitIDchar = argv[i]; char* FitIDchar2 = strtok (FitIDchar, "="); FitID = FitIDchar2; cout<<"FitID = "<<FitID<<endl;}
		    if(std::string(argv[i]).find("cutName") != std::string::npos) {char* cutNamechar = argv[i]; char* cutNamechar2 = strtok (cutNamechar, "="); cutName = cutNamechar2; cout<<"cutName = "<<cutName<<endl;}
		    if(std::string(argv[i]).find("fileName") != std::string::npos) {char* fileNamechar = argv[i]; char* fileNamechar2 = strtok (fileNamechar, "="); fileName = fileNamechar2; cout<<"fileName = "<<fileName<<endl;}
		    if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
		    if(std::string(argv[i]).find("nToy") != std::string::npos) {char* nToychar = argv[i]; char* nToychar2 = strtok (nToychar, "p"); nToy = atof(nToychar2); cout<<"nToy = "<<nToy<<endl;}

		    if(std::string(argv[i]).find("nYsigCut") != std::string::npos) {nYsigCut=true; cout<<"Optimizing nYsigCut"<<endl;}
		    if(std::string(argv[i]).find("vtxProbCut") != std::string::npos) {vtxProbCut=true; cout<<"Optimizing vtxProbCut"<<endl;}
		    if(std::string(argv[i]).find("ctSigCut") != std::string::npos) {ctSigCut=true; cout<<"Optimizing ctSigCut"<<endl;}
		    if(std::string(argv[i]).find("ctCut") != std::string::npos) {ctCut=true; cout<<"Optimizing ctCut"<<endl;}
		    if(std::string(argv[i]).find("gammaptCut") != std::string::npos) {gammaptCut=true; cout<<"Optimizing gammaptCut"<<endl;}
		    if(std::string(argv[i]).find("RConvUpperCut") != std::string::npos) {RConvUpperCut=true; cout<<"Optimizing RConvUpperCut"<<endl;}
		    if(std::string(argv[i]).find("RConvLowerCut") != std::string::npos) {RConvLowerCut=true; cout<<"Optimizing RConvLowerCut"<<endl;}
		    if(std::string(argv[i]).find("pTCut") != std::string::npos) {pTCut=true; cout<<"Optimizing pTCut"<<endl;}
		    if(std::string(argv[i]).find("vtxChi2ProbGammaCut") != std::string::npos) {vtxChi2ProbGammaCut=true; cout<<"Optimizing vtxChi2ProbGammaCut"<<endl;}
		    if(std::string(argv[i]).find("vtxChi2ProbGammaLogCut") != std::string::npos) {vtxChi2ProbGammaLogCut=true; cout<<"Optimizing vtxChi2ProbGammaLogCut"<<endl;}

		    if(std::string(argv[i]).find("useSBforBKGmodel=true") != std::string::npos) {useSBforBKGmodel=true; cout<<"use SB for BKG model"<<endl;}
		    if(std::string(argv[i]).find("SaveAll=true") != std::string::npos) {SaveAll=true; cout<<"SaveAll"<<endl;}
		    if(std::string(argv[i]).find("SetSignalToZero=true") != std::string::npos) {SetSignalToZero=true; cout<<"SetSignalToZero"<<endl;}
		    if(std::string(argv[i]).find("DrawTextOnPlots=true") != std::string::npos) {DrawTextOnPlots=true; cout<<"DrawTextOnPlots"<<endl;}
		    if(std::string(argv[i]).find("useLeftSB=true") != std::string::npos) {useLeftSB=true; cout<<"useLeftSB"<<endl;}
		    if(std::string(argv[i]).find("useRightSB=true") != std::string::npos) {useRightSB=true; cout<<"useRightSB"<<endl;}
		    if(std::string(argv[i]).find("alteredToy=true") != std::string::npos) {alteredToy=true; cout<<"alteredToy"<<endl;}

	  }


	  char dirstruct[200];
	  sprintf(dirstruct,"Figures/%s",FitID);
	  gSystem->mkdir(dirstruct);

//load RooFit library
	using namespace RooFit;
    gSystem->Load("libRooFit");
    gROOT->SetBatch(1);

//open file
    char filename_[500];
    sprintf(filename_,"rooDS_%s.root",fileName);
    TFile data_file(filename_);

//grab roodataset
    RooDataSet* data_=(RooDataSet*)data_file.Get("d");

//declare RooDataSet data variables
    char invmName[200];

    double Ymass1S;
    double Ymass2S;
    Ymass1S=3.096916;
    Ymass2S=10.02326;

    double invm_min1S;
    double invm_max1S;
    double invm_min2S;
    double invm_max2S;
    double invm_min3S;
    double invm_max3S;

    invm_min1S=3.1; invm_max1S=4;
    invm_min2S=3.1; invm_max2S=4;
    invm_min3S=3.1; invm_max3S=4;

    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{J/#psi^{PDG}} [GeV]");
    RooRealVar invm1S("invm1S",invmName,invm_min1S,invm_max1S);
    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(2S)^{PDG}} [GeV]");
    RooRealVar invm2S("invm2S",invmName,invm_min2S,invm_max2S);
    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(3S)^{PDG}} [GeV]");
    RooRealVar invm3S("invm3S",invmName,invm_min3S,invm_max3S);

    invm1S.setRange("TotalRange1S",invm_min1S,invm_max1S);
    invm2S.setRange("TotalRange2S",invm_min2S,invm_max2S);
    invm3S.setRange("TotalRange3S",invm_min3S,invm_max3S);

    RooRealVar jpsipt   = RooRealVar("jpsipt", "jpsipt",0,100);
    RooRealVar jpsimass   = RooRealVar("jpsimass", "jpsimass",1,5);
    RooRealVar gammapt   = RooRealVar("gammapt", "gammapt",0,100);
    RooRealVar deltaRChiJpsi  = RooRealVar("deltaRChiJpsi","deltaRChiJpsi",0,1);
    RooRealVar deltaRJpsig    = RooRealVar("deltaRJpsig","deltaRJpsig",0,1);
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
    RooRealVar Y1Smass_nSigma    = RooRealVar("Y1Smass_nSigma","Y1Smass_nSigma",-1000,1000);
    RooRealVar Y2Smass_nSigma    = RooRealVar("Y2Smass_nSigma","Y2Smass_nSigma",-1000,1000);
    RooRealVar Y3Smass_nSigma    = RooRealVar("Y3Smass_nSigma","Y3Smass_nSigma",-1000,1000);
    RooRealVar vertexChi2ProbGamma    = RooRealVar("vertexChi2ProbGamma","vertexChi2ProbGamma",0,100);

    RooRealVar jpsipx   = RooRealVar("jpsipx", "jpsipx",-100,100);
    RooRealVar jpsipy   = RooRealVar("jpsipy", "jpsipy",-100,100);
    RooRealVar jpsipz   = RooRealVar("jpsipz", "jpsipz",-100,100);
    RooRealVar gammapx   = RooRealVar("gammapx", "gammapx",-100,100);
    RooRealVar gammapy   = RooRealVar("gammapy", "gammapy",-100,100);
    RooRealVar gammapz   = RooRealVar("gammapz", "gammapz",-100,100);
    RooRealVar Q   = RooRealVar("Q", "Q",0,100);

    int n_Cut=1;
    double n_Cut_0;
    double deltaCut;
    int n_Cut_Default=100;

// nYsig
    if(nYsigCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.1;
        deltaCut=0.05;
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
        n_Cut_0=0.1;
        deltaCut=0.025;
    }
// ct
    if(ctCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.0001;
        deltaCut=0.00025;
    }
// gammapt
    if(gammaptCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0.5;
        deltaCut=0.02;
    }
// RConvUpperCut
    if(RConvUpperCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=10;
        deltaCut=0.2;
    }
// RConvLowerCut
    if(RConvLowerCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=1.5;
        deltaCut=0.025;
    }
// pT
    if(pTCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0;
        deltaCut=0.2;
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

    TH1F  *hCut_3Psig_1S = new TH1F("hCut_3Psig_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PSB_1S = new TH1F("hCut_3PSB_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PS_1S = new TH1F("hCut_3PS_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PB_1S = new TH1F("hCut_3PB_1S","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_3PsigALL = new TH1F("hCut_3PsigALL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_PhotonEnergyScale = new TH1F("hCut_PhotonEnergyScale","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_PhotonEnergyScale2P = new TH1F("hCut_PhotonEnergyScale2P","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_1PsigLL = new TH1F("hCut_1PsigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PsigLL = new TH1F("hCut_2PsigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_sigLLhalf = new TH1F("hCut_sigLLhalf","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_3PsigLLhalf = new TH1F("hCut_3PsigLLhalf","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);


    for(int inCut=0; inCut<n_Cut; inCut++){

///////////////////////////////////////////////
//////////// CUTS /////////////////////////////
///////////////////////////////////////////////
    	cout<<"Cut number "<<inCut<<endl;
        data_->Print();

// Y cuts

    double cut_rap=1000;
    double nYsig=2.5;
    double cut_vtxProb=0.01;
    double cut_ctSig=1000;
    double cut_ct=0.01;
    double cut_Ypt=0.;

// Gamma cuts

    double cut_gammapt = 1.1;
    double cut_gammapt1S = cut_gammapt;

    double cut_RconvMin = 2.5;
    double cut_RconvMax = 200;
    double cut_vtxChi2ProbGamma = 0.;//1e-5;//000005;
    double cut_vtxChi2ProbGammaLog = -1000.;//1e-5;//000005;

// Y - Gamma cuts



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

//Cut RooDataSet


    char DScutChar[1000];

    sprintf(DScutChar,"jpsipt>%f",cut_Ypt);
	cout<<DScutChar<<endl;
    RooDataSet* data_______=(RooDataSet*)data_->reduce(Cut(DScutChar));
    sprintf(DScutChar,"gammapt>%f",cut_gammapt);
	cout<<DScutChar<<endl;
    RooDataSet* data__=(RooDataSet*)data_______->reduce(Cut(DScutChar));
    data__->Print();
    sprintf(DScutChar,"Rconv > %f && Rconv < %f",cut_RconvMin,cut_RconvMax);
    cout<<DScutChar<<endl;
    RooDataSet* data___=(RooDataSet*)data__->reduce(Cut(DScutChar));
    data___->Print();
    sprintf(DScutChar,"jpsieta < %f && jpsieta > %f",cut_rap,-cut_rap);
    cout<<DScutChar<<endl;
    RooDataSet* data____=(RooDataSet*)data___->reduce(Cut(DScutChar));
    data____->Print();
    sprintf(DScutChar,"ctpv < %f &&  ctpv > %f",cut_ct,-cut_ct);
    cout<<DScutChar<<endl;
    RooDataSet* data_____=(RooDataSet*)data____->reduce(Cut(DScutChar));
    data_____->Print();

//    sprintf(DScutChar,"log10(vertexChi2ProbGamma) > %f",cut_vtxChi2ProbGamma);
    sprintf(DScutChar,"vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f",cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog);
    cout<<DScutChar<<endl;
    RooDataSet* data______vtx=(RooDataSet*)data_____->reduce(Cut(DScutChar));
    data______vtx->Print();

    sprintf(DScutChar,"jpsiVprob > %f",cut_vtxProb);
    cout<<DScutChar<<endl;
    RooDataSet* data______=(RooDataSet*)data______vtx->reduce(Cut(DScutChar));
    data______->Print();
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
        	LowerSBmin=2.8;
        	nSigSBlow=-3.5;
        	nSigSBhigh=3.5;
        	UpperSBmax=4;

        	LowerSBmin2S=8.7;
        	nSigSBlow2S=-3.5;
        	nSigSBhigh2S=3.5;
        	UpperSBmax2S=10.85;
        }

        if(useLeftSB&&!useRightSB){
        	LowerSBmin=2.8;
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

        double MassSignalMin=3;
        double MassSignalMax=3.2;

    sprintf(DScutChar,"jpsimass > %f && Y1Smass_nSigma < %f || jpsimass < %f && Y2Smass_nSigma > %f",LowerSBmin,nSigSBlow,UpperSBmax,nSigSBhigh);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB1S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    dataSB1S->Print();

    sprintf(DScutChar,"jpsimass > %f && Y1Smass_nSigma < %f || jpsimass < %f && Y2Smass_nSigma > %f",LowerSBmin2S,nSigSBlow2S,UpperSBmax2S,nSigSBhigh2S);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB2S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    dataSB2S->Print();


/////// Producing Signal dataset
    sprintf(DScutChar,"jpsimass > %f && jpsimass < %f",MassSignalMin,MassSignalMax);
    cout<<DScutChar<<endl;
    RooDataSet* data1S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    data1S->Print();

    RooDataSet* data2S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    data2S->Print();

    sprintf(DScutChar,"Y3Smass_nSigma > %f && Y3Smass_nSigma < %f && Y2Smass_nSigma > %f",-nYsig,nYsig,nYsig);
    cout<<DScutChar<<endl;
    RooDataSet* data3S=(RooDataSet*)data______->reduce(Cut(DScutChar));
//    data3S->Print();


    if(SetSignalToZero&&useSBforBKGmodel){
    	data1S=dataSB1S;
    	data2S=dataSB2S;
    }
    TTree* tree1S=(TTree*)data1S->tree();
    TTree* tree2S=(TTree*)data2S->tree();
    TTree* treeAllCuts=(TTree*)data______->tree();

    if(useSBforBKGmodel){
    	tree1S=(TTree*)dataSB1S->tree();
    	tree2S=(TTree*)dataSB2S->tree();
    }

    if(SetSignalToZero){
    	data1S=dataSB1S;
    	data2S=dataSB2S;
    }

    char treeName[200];
    sprintf(treeName,"Figures/%s/tree_Y1S.root",FitID);
    if(SaveAll) tree1S->SaveAs(treeName);
    sprintf(treeName,"Figures/%s/tree_Y2S.root",FitID);
    if(SaveAll) tree2S->SaveAs(treeName);
    sprintf(treeName,"Figures/%s/treeAllCuts.root",FitID);
    if(SaveAll) treeAllCuts->SaveAs(treeName);

    if(SaveAll){
    char saveMass[200];
    TH1F  *MassHisto = new TH1F("MassHisto","",100,2,5);
    treeAllCuts->Draw("jpsimass>>MassHisto");
    TCanvas* MassCanvas = new TCanvas("MassCanvas","MassCanvas",1600, 800);
    MassCanvas->SetFillColor(kWhite);
    MassHisto->GetYaxis()->SetTitleOffset(1.7); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); MassHisto->SetStats(0);MassHisto->GetXaxis()->SetTitle("m_{#mu#mu}");	MassHisto->Draw();
    MassCanvas->Modified();
    sprintf(saveMass,"Figures/%s/DimuonMass.pdf",FitID);
    MassCanvas->SaveAs(saveMass);

    TH1F  *MassHisto1S = new TH1F("MassHisto1S","",100,-20,20);
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

    }

    // Ernest's ToyMC for background estimation:

    double bb_thresh=3.585;
    double chib1Pmin=3.325;
    double chib1Pmax=9.95;
    double chib2Pmin=10.15;
    double chib2Pmax=10.325;
    double chic1min=3.25;

        int nHistBins=1000;
        float m_gamma = 0.;
        char DrawChar[500];

        sprintf(DrawChar,"invm1S>0");//chic1min
        if(useSBforBKGmodel) sprintf(DrawChar,"invm1S<10000");

        TH1F  *hYmass1S = new TH1F("hYmass1S","",nHistBins,2,5);
        TH1F  *hYmass_check1S = new TH1F("hYmass_check1S","",nHistBins,2,5);
        TH1F  *hYmass2S = new TH1F("hYmass2S","",nHistBins,2,5);
        TH1F  *hYmass_check2S = new TH1F("hYmass_check2S","",nHistBins,2,5);

       tree1S->Draw("jpsimass>>hYmass1S",DrawChar);
       tree2S->Draw("jpsimass>>hYmass2S");

        TH1F  *hGammaP1S = new TH1F("hGammaP1S","",nHistBins,0,100);
        TH1F  *hGammaP_check1S = new TH1F("hGammaP_check1S","",nHistBins,0,100);
        TH1F  *hGammaP2S = new TH1F("hGammaP2S","",nHistBins,0,100);
        TH1F  *hGammaP_check2S = new TH1F("hGammaP_check2S","",nHistBins,0,100);

        tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP1S",DrawChar);
        tree2S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP2S");

        TH1F  *hUpsP1S = new TH1F("hUpsP1S","",nHistBins,0,100);
        TH1F  *hUpsP_check1S = new TH1F("hUpsP_check1S","",nHistBins,0,100);
        TH1F  *hUpsP2S = new TH1F("hUpsP2S","",nHistBins,0,100);
        TH1F  *hUpsP_check2S = new TH1F("hUpsP_check2S","",nHistBins,0,100);

        tree1S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP1S",DrawChar);
        tree2S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP2S");

        TH1F  *hCosAlphaP1S = new TH1F("hCosAlphaP1S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check1S = new TH1F("hCosAlphaP_check1S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP2S = new TH1F("hCosAlphaP2S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check2S = new TH1F("hCosAlphaP_check2S","",nHistBins,-1,1);

        tree1S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP1S",DrawChar);
        tree2S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP2S");


        bool saveBkgToy=false;

        char bkgToyname[200];
        sprintf(bkgToyname,"%s/bkgToy.root",dirstruct);
        bool existingBKGfile=false;
        TFile *fQ = new TFile(bkgToyname,"RECREATE");
        TTree* tQ = new TTree("tQ","toy MC");


//        bool existingBKGfile=true;
//        TFile *fQ = new TFile("/Users/valentinknuenz/usr/local/workspace/Chib2011_/Chi2OniaGamma/analysis/bkgToy.root","READ");
//        TTree* tQ = (TTree*)fQ->Get("tQ");



        if(existingBKGfile) {
        	cout<<"using existing bkt-toyMC file"<<endl;
        }

        if(!existingBKGfile){
        	cout<<"producing new bkt-toyMC file"<<endl;

        hYmass1S->Write();
        hGammaP1S->Write();
        hUpsP1S->Write();
        hCosAlphaP1S->Write();


        float Q1S,Q2S;
        float M1S,M2S;

//        tQ->Branch("Q1S",&Q1S,"Q1S/F");
        tQ->Branch("invm1S",&M1S,"invm1S/F");
//        tQ->Branch("Q2S",&Q2S,"Q2S/F");
        tQ->Branch("invm2S",&M1S,"invm2S/F");

        TRandom3 *randOM = new TRandom3();

        TF1 *f1 = new TF1("f1","exp(x)",-1,1);

        for(int i=0;i<nToy;i++){

          float Ymass_1S = hYmass1S->GetRandom();
          float p_gamma1S = hGammaP1S->GetRandom();
          float p_Ups1S = hUpsP1S->GetRandom();

          if(!alteredToy) Ymass_1S=Ymass1S;

//          float cosAlpha1S = 2*randOM->Uniform()-1;
          float cosAlpha1S = hCosAlphaP1S->GetRandom();
//          cosAlpha1S = cosAlpha1S*randOM->Uniform(0.988,1);
          float e_gamma1S = sqrt(m_gamma*m_gamma+p_gamma1S*p_gamma1S);
          float e_Ups1S = sqrt(Ymass_1S*Ymass_1S+p_Ups1S*p_Ups1S);


//          Q1S = sqrt(m_gamma*m_gamma + Ymass_1S*Ymass_1S + 2*e_gamma1S*e_Ups1S - 2*p_gamma1S*p_Ups1S*cosAlpha1S) - m_gamma - Ymass_1S;
          M1S = sqrt(m_gamma*m_gamma + Ymass_1S*Ymass_1S + 2*e_gamma1S*e_Ups1S - 2*p_gamma1S*p_Ups1S*cosAlpha1S)-Ymass_1S+Ymass1S;

       	  hUpsP_check1S->Fill(p_Ups1S);
          hGammaP_check1S->Fill(p_gamma1S);
          hYmass_check1S->Fill(Ymass_1S);
          hCosAlphaP_check1S->Fill(cosAlpha1S);


          if (M1S<4.2) {
        	  tQ->Fill();
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

        tQ->Write();
    	cout<<"finalising bkg-toyMC file"<<endl;
        delete randOM;


        }

        if(SaveAll){
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

        }





        TH1F  *BkgToyHist1S = new TH1F("BkgToyHist1S","",100,3,4.5);
        tQ->Draw("invm1S>>BkgToyHist1S");
        TH1F  *BkgToyHist2S = new TH1F("BkgToyHist2S","",100,3,4.5);
        tQ->Draw("invm2S>>BkgToyHist2S");


        invm1S.setBins(50);
        invm2S.setBins(50);

        RooDataSet* RooBkgToySet = new RooDataSet("RooBkgToySet","RooBkgToySet",tQ,RooArgList(invm1S,invm2S));
        RooBkgToySet->Print();
        RooDataHist* RooBkgToyHist1S = new RooDataHist("RooBkgToyHist1S","RooBkgToyHist1S",RooArgList(invm1S),*RooBkgToySet);
        RooBkgToyHist1S->Print();
        RooHistPdf background1S = RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S,1);
        BkgToyHist1S->Write();
        RooDataHist* RooBkgToyHist2S = new RooDataHist("RooBkgToyHist2S","RooBkgToyHist2S",RooArgList(invm2S),*RooBkgToySet);
        RooBkgToyHist2S->Print();
        RooHistPdf background2S = RooHistPdf("background2S","background2S",RooArgSet(invm2S),*RooBkgToyHist2S,1);
        BkgToyHist2S->Write();

        fQ->Write();
        fQ->Close();


//        gSystem->Unlink(bkgToyname);


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



//declare fit variables

    RooRealVar m_chib1P= RooRealVar("Mean #chi_{b}(1P)","Mean #chi_{b}(1P)",9.888,9.85,9.925);
    RooRealVar m_chib2P= RooRealVar("Mean #chi_{b}(2P)","Mean #chi_{b}(2P)",10.25,10.225,10.275);

    RooRealVar sigma  = RooRealVar("#sigma","#sigma",0.01,0.005,0.02);
    RooRealVar alpha  = RooRealVar("#alpha","#alpha",0.654,0.,5.);
    RooRealVar n      = RooRealVar("n","n",2.58,.5,15.);
    RooRealVar sigma1  = RooRealVar("#sigma1","#sigma1",0.01,0.005,0.02);
    RooRealVar sigma2  = RooRealVar("#sigma2","#sigma2",0.01,0.005,0.02);
    RooRealVar sigma3  = RooRealVar("#sigma3","#sigma3",1.2094e-02);



    double m_chib3P_START=10.45;
    double m_chib3P_1S_START=10.45;

    RooRealVar m_chib3P= RooRealVar("Mean #chi_{b}(3P)","Mean #chi_{b}(3P)",m_chib3P_START,10.325,10.56);
    RooRealVar m_chib3P_1S= RooRealVar("Mean #chi_{b}(3P) in 1S sample","Mean #chi_{b}(3P) in 1S sample",m_chib3P_1S_START,10.325,10.56);

    m_chib3P.setVal(m_chib3P_START);
    m_chib3P_1S.setVal(m_chib3P_1S_START);



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
    double alpha_3J_START=0.74;
    double n_3J_START=3.;
    double PhotonMassScale_START=1.;
    double PhotonMassScale2P_START=1.;
    double PhotonSigmaScale_START=0.017;

    RooRealVar sigmaIndJ  = RooRealVar("#sigma_{J}","#sigma_{J}",sigmaInd_START,0.003,0.015);
    RooRealVar sigmaIndJ1  = RooRealVar("#sigma_{J1}","#sigma_{J}",sigmaInd_START,0.003,0.015);
    RooRealVar sigmaIndJ2  = RooRealVar("#sigma_{J2}","#sigma_{J}",sigmaInd_START,0.003,0.015);

    RooRealVar alpha_3J_J  = RooRealVar("#alpha","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar alpha_3J_J1  = RooRealVar("#alpha_{J1}","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar alpha_3J_J2  = RooRealVar("#alpha_{J2}","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar n_3J_J      = RooRealVar("n","CB_{n}",n_3J_START,.5,15.);
    RooRealVar n_3J_J1      = RooRealVar("n_{J1}","CB_{n}",n_3J_START,.5,15.);
    RooRealVar n_3J_J2      = RooRealVar("n_{J1}","CB_{n}",n_3J_START,.5,15.);
    RooRealVar PhotonMassScale= RooRealVar("PES","PES",PhotonMassScale_START,0.95,1.05);
    RooRealVar PhotonMassScale2P= RooRealVar("PES_{J2}","PhotonMassScale2P",PhotonMassScale2P_START,0.95,1.05);
    RooRealVar PhotonSigmaScale= RooRealVar("#sigma_{Q}/Q","PhotonSigmaScale",PhotonSigmaScale_START,0.01,0.025);
    RooRealVar PhotonSigmaScaleJ2= RooRealVar("#sigma_{Q}/Q, J2","PhotonSigmaScale",PhotonSigmaScale_START,0.01,0.025);

    RooRealVar n_3J      = RooRealVar("CB_{n}","CB_{n}",n_3J_START,.5,15.);
    RooRealVar alpha_3J  = RooRealVar("CB_{#alpha}","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar sigmaInd  = RooRealVar("#sigma_{J}","#sigma_{J}",sigmaInd_START,0.003,0.015);
    RooRealVar sigmaInd2P  = RooRealVar("sigmaInd2P","sigmaInd2P",sigmaInd2P_START,0.003,0.025);
    RooRealVar sigmaInd2P_2S  = RooRealVar("sigmaInd2P_2S","sigmaInd2P_2S",sigmaInd2P_2S_START,0.003,0.025);
    RooRealVar sigmaInd3P  = RooRealVar("sigmaInd3P","sigmaInd3P",sigmaInd3P_START,0.003,0.025);
    RooRealVar sigmaInd3P_1S  = RooRealVar("sigmaInd3P_1S","sigmaInd3P_1S",sigmaInd3P_1S_START,0.003,0.04);

    sigmaInd.setVal(sigmaInd_START);
    sigmaInd2P.setVal(sigmaInd2P_START);
    sigmaInd2P_2S.setVal(sigmaInd2P_2S_START);
    sigmaInd3P.setVal(sigmaInd3P_START);
    sigmaInd3P_1S.setVal(sigmaInd3P_1S_START);
    alpha_3J.setVal(alpha_3J_START);
    n_3J.setVal(n_3J_START);
    PhotonMassScale.setVal(PhotonMassScale_START);
    PhotonMassScale2P.setVal(PhotonMassScale2P_START);
    PhotonSigmaScale.setVal(PhotonSigmaScale_START);

//    n_3J.setConstant();
//    alpha_3J.setConstant();

    double ratio_J2overJ1_START=.5;
    RooRealVar ratio_J2overJ1= RooRealVar("ratio_J2overJ1","ratio_J2overJ1",ratio_J2overJ1_START);//,0.,5.);
    ratio_J2overJ1.setVal(ratio_J2overJ1_START);
	RooFormulaVar fracJ1("fracJ1", "1./(1.+@0)", RooArgList(ratio_J2overJ1));







    RooRealVar fractionJ0= RooRealVar("N_{#chi_{c0}}/N_{#chi_{cJ}}","#chi_{c0}/#Sigma_{J}(#chi_{cJ})",0.05,0.,0.1);
    RooRealVar fractionJ1= RooRealVar("N_{#chi_{c1}}/N_{#chi_{cJ}}","#chi_{c1}/#Sigma_{J}(#chi_{cJ})",0.4,0.,.8);

    RooRealVar m_chib1P0fix= RooRealVar("m_chib1P0fix","m_chib1P0fix",3.41475);
    RooRealVar m_chib1P1fix= RooRealVar("m_chib1P1fix","m_chib1P1fix",3.51066);
    RooRealVar m_chib1P2fix= RooRealVar("m_chib1P2fix","m_chib1P2fix",3.55620);

    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass1S,Ymass1S);

    sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass1S);

    RooFormulaVar m_chib1P0float("m_chib1P0float", MassScaleFormula, RooArgList(m_chib1P0fix, PhotonMassScale));
    RooFormulaVar m_chib1P1float("m_chib1P1float", MassScaleFormula, RooArgList(m_chib1P1fix, PhotonMassScale));
	RooFormulaVar m_chib1P2float("m_chib1P2float", MassScaleFormula, RooArgList(m_chib1P2fix, PhotonMassScale));

    RooFormulaVar w_chib1P0float("w_chib1P0float", SigmaScaleFormula, RooArgList(m_chib1P0float, PhotonSigmaScale));
    RooFormulaVar w_chib1P1float("w_chib1P1float", SigmaScaleFormula, RooArgList(m_chib1P1float, PhotonSigmaScale));
    RooFormulaVar w_chib1P2float("w_chib1P2float", SigmaScaleFormula, RooArgList(m_chib1P2float, PhotonSigmaScale));

	RooCBShape chib1P0      = RooCBShape ("chib1P0","chib1P0",invm1S,m_chib1P0float,w_chib1P0float,alpha_3J_J1,n_3J_J1);
	RooCBShape chib1P1      = RooCBShape ("chib1P1","chib1P1",invm1S,m_chib1P1float,w_chib1P1float,alpha_3J_J1,n_3J_J1);
    RooCBShape chib1P2      = RooCBShape ("chib1P2","chib1P2",invm1S,m_chib1P2float,w_chib1P2float,alpha_3J_J2,n_3J_J1);

    double n1PnJ_START=12000;
    RooRealVar n1PnJ= RooRealVar("N_{#chi_{cJ}}","#Sigma_{J}(#chi_{cJ})",n1PnJ_START,1000,20000);
    n1PnJ.setVal(n1PnJ_START);

    RooAddPdf chib1P_sigInd= RooAddPdf("chib1P_sigInd","chib1P_sigInd",RooArgList(chib1P0,chib1P1,chib1P2),RooArgList(fractionJ0,fractionJ1));

    RooFormulaVar chib1P_nevt("chib1P_nevt", "@1*@0+@2*@0+(1-@1-@2)*@0", RooArgList(n1PnJ,fractionJ0,fractionJ1));
    RooAddPdf chib1P_sig= RooAddPdf("chib1P_sig","chib1P_sig",RooArgList(chib1P_sigInd),RooArgList(chib1P_nevt));














    RooRealVar m_chib2P1fix= RooRealVar("m_chib2P1fix","m_chib2P1fix",10.25546);
    RooRealVar m_chib2P2fix= RooRealVar("m_chib2P2fix","m_chib2P2fix",10.26865);

	RooFormulaVar m_chib2P1float("m_chib2P1float", MassScaleFormula, RooArgList(m_chib2P1fix, PhotonMassScale));
	RooFormulaVar m_chib2P2float("m_chib2P2float", MassScaleFormula, RooArgList(m_chib2P2fix, PhotonMassScale));

    RooCBShape chib2P1      = RooCBShape ("chib2P1","chib2P1",invm1S,m_chib2P1float,sigmaInd2P,alpha_3J,n_3J);
    RooCBShape chib2P2      = RooCBShape ("chib2P2","chib2P2",invm1S,m_chib2P2float,sigmaInd2P,alpha_3J,n_3J);

    double n2PnJ_START=160;
    RooRealVar n2PnJ= RooRealVar("n2PJ","n2PJ",n2PnJ_START,0,450);
    n2PnJ.setVal(n2PnJ_START);
    RooAddPdf chib2P_sigInd= RooAddPdf("chib2P_sigInd","chib2P_sigInd",RooArgList(chib2P1,chib2P2),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt("chib2P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n2PnJ,fracJ1));
    RooAddPdf chib2P_sig= RooAddPdf("chib2P_sig","chib2P_sig",RooArgList(chib2P_sigInd),RooArgList(chib2P_nevt));


    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass2S,Ymass2S);

	RooFormulaVar m_chib2P1float_2S("m_chib2P1float_2S", MassScaleFormula, RooArgList(m_chib2P1fix, PhotonMassScale));
	RooFormulaVar m_chib2P2float_2S("m_chib2P2float_2S", MassScaleFormula, RooArgList(m_chib2P2fix, PhotonMassScale));


	RooCBShape chib2P1_2S      = RooCBShape ("chib2P1_2S","chib2P1_2S",invm2S,m_chib2P1float_2S,sigmaInd2P_2S,alpha_3J,n_3J);
    RooCBShape chib2P2_2S      = RooCBShape ("chib2P2_2S","chib2P2_2S",invm2S,m_chib2P2float_2S,sigmaInd2P_2S,alpha_3J,n_3J);

    double n2PnJ_2S_START=8;
    RooRealVar n2PnJ_2S= RooRealVar("n2PnJ_2S","n2PnJ_2S",n2PnJ_2S_START,0,200);
    n2PnJ_2S.setVal(n2PnJ_2S_START);
    RooAddPdf chib2P_sigInd_2S= RooAddPdf("chib2P_sigInd_2S","chib2P_sigInd_2S",RooArgList(chib2P1_2S,chib2P2_2S),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt_2S("chib2P_nevt_2S", "@1*@0+(1-@1)*@0", RooArgList(n2PnJ_2S,fracJ1));
    RooAddPdf chib2P_sig_2S= RooAddPdf("chib2P_sig_2S","chib2P_sig_2S",RooArgList(chib2P_sigInd_2S),RooArgList(chib2P_nevt_2S));


    sprintf(MassScaleFormula,"(@0-%f-0.012*(1-@2))*@1+%f",Ymass2S,Ymass2S);
	RooFormulaVar m_chib3P1float("m_chib3P1float", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale, fracJ1));
    sprintf(MassScaleFormula,"(@0-%f+0.012*@2)*@1+%f",Ymass2S,Ymass2S);
	RooFormulaVar m_chib3P2float("m_chib3P2float", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale, fracJ1));

    RooCBShape chib3P1      = RooCBShape ("chib3P1","chib3P1",invm2S,m_chib3P1float,sigmaInd3P,alpha_3J,n_3J);
    RooCBShape chib3P2      = RooCBShape ("chib3P2","chib3P2",invm2S,m_chib3P2float,sigmaInd3P,alpha_3J,n_3J);

    double n3PnJ_START=8;
    RooRealVar n3PnJ= RooRealVar("n3PnJ","n3PnJ",n3PnJ_START,0,100);
    n3PnJ.setVal(n3PnJ_START);
    RooAddPdf chib3P_sigInd= RooAddPdf("chib3P_sigInd","chib3P_sigInd",RooArgList(chib3P1,chib3P2),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt("chib3P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ,fracJ1));
    RooAddPdf chib3P_sig= RooAddPdf("chib3P_sig","chib3P_sig",RooArgList(chib3P_sigInd),RooArgList(chib3P_nevt));





    sprintf(MassScaleFormula,"(@0-%f-0.012*(1-@2))*@1+%f",Ymass1S,Ymass1S);
	RooFormulaVar m_chib3P1float_1S("m_chib3P1float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale, fracJ1));
    sprintf(MassScaleFormula,"(@0-%f+0.012*@2)*@1+%f",Ymass1S,Ymass1S);
	RooFormulaVar m_chib3P2float_1S("m_chib3P2float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale, fracJ1));

    RooCBShape chib3P1_1S      = RooCBShape ("chib3P1_1S","chib3P1_1S",invm1S,m_chib3P1float_1S,sigmaInd3P_1S,alpha_3J,n_3J);
    RooCBShape chib3P2_1S     = RooCBShape ("chib3P2_1S","chib3P2_1S",invm1S,m_chib3P2float_1S,sigmaInd3P_1S,alpha_3J,n_3J);

    double n3PnJ_1S_START=8;
    RooRealVar n3PnJ_1S= RooRealVar("n3PnJ_1S","n3PnJ_1S",n3PnJ_1S_START,0,250);
    n3PnJ_1S.setVal(n3PnJ_1S_START);
    RooAddPdf chib3P_sigInd_1S= RooAddPdf("chib3P_sigInd_1S","chib3P_sigInd_1S",RooArgList(chib3P1_1S,chib3P2_1S),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt_1S("chib3P_nevt_1S", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ_1S,fracJ1));
    RooAddPdf chib3P_sig_1S= RooAddPdf("chib3P_sig_1S","chib3P_sig_1S",RooArgList(chib3P_sigInd_1S),RooArgList(chib3P_nevt_1S));

////////////////////////////////////////////////////////


    RooRealVar background_nevt1S = RooRealVar("N_{bkg}1S","N_{bkg}1S",1500,0,20000);
    RooRealVar background_nevtSB1S = RooRealVar("background_nevtSB1S","background_nevtSB1S",1500,0,10000);
    RooAddPdf backgroundSB1S= RooAddPdf("backgroundSB1S","backgroundSB1S",RooArgList(background1S),RooArgList(background_nevtSB1S));
    RooAddPdf backgroundE1S= RooAddPdf("backgroundE1S","backgroundE1S",RooArgList(background1S),RooArgList(background_nevt1S));
    RooRealVar background_nevt2S = RooRealVar("N_{bkg}2S","N_{bkg}2S",1500,0,10000);
    RooRealVar background_nevtSB2S = RooRealVar("background_nevtSB2S","background_nevtSB2S",1500,0,10000);
    RooAddPdf backgroundSB2S= RooAddPdf("backgroundSB2S","backgroundSB2S",RooArgList(background2S),RooArgList(background_nevtSB2S));
    RooAddPdf backgroundE2S= RooAddPdf("backgroundE2S","backgroundE2S",RooArgList(background2S),RooArgList(background_nevt2S));


    sprintf(DScutChar,"invm1S > %f && invm1S < %f",invm_min1S,invm_max1S);
    RooDataSet* data1S_masswindow=(RooDataSet*)data1S->reduce(Cut(DScutChar));
    sprintf(DScutChar,"invm2S > %f && invm2S < %f",invm_min2S,invm_max2S);
    RooDataSet* data2S_masswindow=(RooDataSet*)data2S->reduce(Cut(DScutChar));



    RooAddPdf modelPdf1S= RooAddPdf("modelPdf1S","modelPdf1S",RooArgList(chib1P_sig,background1S),RooArgList(chib1P_nevt,background_nevt1S));
    RooAddPdf modelPdf1S_= RooAddPdf("modelPdf1S_","modelPdf1S_",RooArgList(chib1P_sig,background1S),RooArgList(chib1P_nevt,background_nevt1S));

    RooAddPdf modelPdf2S= RooAddPdf("modelPdf2S","modelPdf2S",RooArgList(chib3P_sig,chib2P_sig_2S,background2S),RooArgList(chib3P_nevt,chib2P_nevt_2S,background_nevt2S));

    RooAddPdf modelPdf1S_1P2P= RooAddPdf("modelPdf1S_1P2P","modelPdf1S_1P2P",RooArgList(chib1P_sig,chib2P_sig,background1S),RooArgList(chib1P_nevt,chib2P_nevt,background_nevt1S));
    RooAddPdf modelPdf1S_3P= RooAddPdf("modelPdf1S_3P","modelPdf1S_3P",RooArgList(chib3P_sig_1S,background1S),RooArgList(chib3P_nevt_1S,background_nevt1S));
    RooAddPdf modelPdf2S_red= RooAddPdf("modelPdf2S_red","modelPdf2S_red",RooArgList(chib3P_sig,chib2P_sig_2S,background2S),RooArgList(chib3P_nevt,chib2P_nevt_2S,background_nevt2S));




    invm1S.setRange("SBregion2",chib1Pmax,chib2Pmin);
    invm1S.setRange("SBregion31S",bb_thresh,invm_max1S);

    invm1S.setRange("Chib1Pregion",invm_min1S,bb_thresh);
    invm1S.setRange("Chib2Pregion",chib2Pmin,chib2Pmax);
    invm1S.setRange("Chib3Pregion",chib2Pmax,bb_thresh);

    invm1S.setRange("TotalBlindedRegion1",invm_min1S,chib2Pmax);
    invm1S.setRange("TotalBlindedRegion2",bb_thresh,invm_max1S);

    invm1S.setRange("Chib1P2Pregion",invm_min1S,chib2Pmax);

    invm2S.setRange("Y2SAll",invm_min2S,invm_max2S);
    invm2S.setRange("SBregion32S",bb_thresh,invm_max2S);
    invm2S.setRange("Chib3Pregion2S",invm_min2S,bb_thresh);

    char SBregion2Char[200];
    sprintf(SBregion2Char,"invm1S > %f && invm1S < %f",bb_thresh,invm_max1S);
    char SBregion31SChar[200];
    sprintf(SBregion31SChar,"invm1S > %f && invm1S < %f",invm_min1S,chic1min);
    char SBregion32SChar[200];
    sprintf(SBregion32SChar,"invm2S > %f && invm2S < %f",bb_thresh,invm_max2S);

    double Nbkg_1S_SB1=data1S->sumEntries(SBregion2Char);
    double Nbkg_1S_SB2=data1S->sumEntries(SBregion31SChar);
    double Nbkg_2S_SB1=data2S->sumEntries(SBregion32SChar);
    cout<<"background_nevt1S sumentriesrange = " << Nbkg_1S_SB1<<endl;
    cout<<"background_nevt1S sumentriesrange2 = " << Nbkg_1S_SB2<<endl;
    cout<<"background_nevt2S sumentriesrange = " << Nbkg_2S_SB1<<endl;


    double BkgSBInt1S;
    BkgSBInt1S = NormalizedIntegral(&backgroundSB1S, invm1S, bb_thresh,invm_max1S);
    cout<<BkgSBInt1S<<endl;
    BkgSBInt1S = BkgSBInt1S+ NormalizedIntegral(&backgroundSB1S, invm1S, invm_min1S,chic1min);
    cout<<BkgSBInt1S<<endl;
    background_nevt1S.setVal((Nbkg_1S_SB1+Nbkg_1S_SB2)/BkgSBInt1S);
    background_nevt1S.setConstant();
    cout<<"background normalization 1S = " << background_nevt1S.getVal()<<endl;

    double BkgSBInt2S;
    BkgSBInt2S = NormalizedIntegral(&backgroundSB2S, invm2S, bb_thresh, invm_max2S);
    cout<<BkgSBInt2S<<endl;
    background_nevt2S.setVal(Nbkg_2S_SB1/BkgSBInt2S);
    background_nevt2S.setConstant();
    cout<<"background normalization 2S = " << background_nevt2S.getVal()<<endl;


    double covQual_BG2S=0;
    double covQual_1P2P2S=0;
    double covQual_BG1S=0;
    double covQual_1P2P1S=0;


    double Nevt_buff=n3PnJ_1S.getVal();
    n3PnJ_1S.setVal(0);
    n3PnJ_1S.setConstant();
    m_chib3P.setConstant();
    sigmaInd3P_1S.setConstant();

    char SignormRegion[200];
    sprintf(SignormRegion,"Chib1Pregion");
    char OnePTwoPDatacut[200];
    sprintf(OnePTwoPDatacut,"invm1S<%f",chib2Pmax);



    RooDataSet* data1S_=(RooDataSet*)data1S->reduce(Cut(OnePTwoPDatacut));

        RooFitResult* Fit1P2P = modelPdf1S_.fitTo(*data1S,Save(1),Range(SignormRegion));
        Fit1P2P->Print();
        covQual_1P2P1S = Fit1P2P->covQual();

        n3PnJ_1S.setConstant(kFALSE);
        n3PnJ_1S.setVal(Nevt_buff);
        m_chib3P.setConstant(kFALSE);
        sigmaInd3P_1S.setConstant(kFALSE);

        sigma.setConstant();
        alpha.setConstant();
        n.setConstant();

//        sigmaInd.setConstant();
        sigmaInd2P.setConstant();
//        alpha_3J.setConstant();
//        n_3J.setConstant();
//        PhotonMassScale2P.setConstant();
        ratio_J2overJ1.setConstant();

//        PhotonMassScale.setVal(PhotonMassScale.getVal()-PhotonMassScale.getError());
//        PhotonMassScale.setConstant();


        background_nevt1S.setConstant();
        background_nevt2S.setConstant();

        delete Fit1P2P;
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

double deltaMJ0=(m_chib1P0float.getVal()-m_chib1P0fix.getVal())*1000;
double deltaMJ1=(m_chib1P1float.getVal()-m_chib1P1fix.getVal())*1000;
double deltaMJ2=(m_chib1P2float.getVal()-m_chib1P2fix.getVal())*1000;

m_chib1P.setVal(m_chib1P1float.getVal()*fracJ1.getVal()+m_chib1P2float.getVal()*(1-fracJ1.getVal()));
m_chib2P.setVal(m_chib2P1float.getVal()*fracJ1.getVal()+m_chib2P2float.getVal()*(1-fracJ1.getVal()));
cout<<"m_chib1P: "<<m_chib1P.getVal()<<endl;
cout<<"m_chib2P: "<<m_chib2P.getVal()<<endl;




        cout<<"chib1P_nevt: "<<chib1P_nevt.getVal()<<endl;
        cout<<"chib2P_nevt: "<<chib2P_nevt.getVal()<<endl;
        cout<<"chib3P_nevt_1S: "<<chib3P_nevt_1S.getVal()<<endl;
        cout<<"chib3P_nevt: "<<chib3P_nevt.getVal()<<endl;
        cout<<"chib2P_nevt_2S: "<<chib2P_nevt_2S.getVal()<<endl;


     double TotalRealEvents1S = data1S_masswindow->sumEntries();
     double totalEventsInFit1S;
     double TotalRealEvents2S = data2S_masswindow->sumEntries();
     double totalEventsInFit2S;

     totalEventsInFit1S = chib1P_nevt.getVal()+background_nevt1S.getVal();
     totalEventsInFit2S = chib3P_nevt.getVal()+chib2P_nevt_2S.getVal()+background_nevt2S.getVal();

     cout<< "N_tot_Fit1S =                            "<<totalEventsInFit1S <<endl;
     cout<< "N_tot_Samlpe1S =                         "<<TotalRealEvents1S <<endl;
     cout<< "Delta_N_tot1S = N_tot_Sample1S-N_tot_Fit1S = "<<TotalRealEvents1S-totalEventsInFit1S <<endl;

     cout<< "N_tot_Fit2S =                            "<<totalEventsInFit2S <<endl;
     cout<< "N_tot_Samlpe2S =                         "<<TotalRealEvents2S <<endl;
     cout<< "Delta_N_tot2S = N_tot_Sample2S-N_tot_Fit2S = "<<TotalRealEvents2S-totalEventsInFit2S <<endl;


// likelihood gaussian check

     double nSig=1.5;

     double fmchib1P = m_chib1P.getVal();
     double fmchib2P = m_chib2P.getVal();
     double fmchib3P = m_chib3P.getVal();
     double err_fmchib3P = m_chib3P.getError();
     double fmchib3P_1S = m_chib3P.getVal();//m_chib3P_1S.getVal();
     double err_fmchib3P_1S = m_chib3P.getError();//m_chib3P_1S.getError();


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
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib3P_sig, invm2S, (fmchib3P-Ymass2S)*PhotonMassScale.getVal()+Ymass2S-sigmaCalc, (fmchib3P-Ymass2S)*PhotonMassScale.getVal()+Ymass2S+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma3=sigmaCalc;

     sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib3P_sig_1S, invm1S, (fmchib3P_1S-Ymass1S)*PhotonMassScale2P.getVal()+Ymass1S-sigmaCalc, (fmchib3P_1S-Ymass1S)*PhotonMassScale2P.getVal()+Ymass1S+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma3_1S=sigmaCalc;



     cout<<"m_chib1P: "<<m_chib1P.getVal()<<endl;
     cout<<"m_chib2P: "<<m_chib2P.getVal()<<endl;
     cout<<"m_chib3P: "<<m_chib3P.getVal()<<endl;
     cout<<"m_chib3P_1S: "<<m_chib3P_1S.getVal()<<endl;
     cout<<"fsigma1: "<<fsigma1<<endl;
     cout<<"fsigma2: "<<fsigma2<<endl;
     cout<<"fsigma3: "<<fsigma3<<endl;
     cout<<"fsigma3_1S: "<<fsigma3_1S<<endl;

     double fnchib1P = chib1P_nevt.getVal();
     double fnchib2P = chib2P_nevt.getVal();
     double fnchib3P = chib3P_nevt.getVal();
     double fnchib3P_1S = chib3P_nevt_1S.getVal();
     double fnbg1S = background_nevt1S.getVal();
     double fnbg2S = background_nevt2S.getVal();
     double rchib1P,rchib2P,rchib3P,rbgchib1P,rbgchib2P,rbgchib3P,rsigchib1P,rsigchib2P,rsigchib3P,rchib3P_1S,rbgchib3P_1S,rsigchib3P_BityukovKrasnikov,rsigchib3P_BityukovKrasnikov_1S,rsigchib3P_ScL,rsigchib3P_ScL_1S;
     double rratchib1P,rratchib2P,rratchib3P,rratchib3P_1S,rsigchib3P_1S;
     rchib1P=fnchib1P*NormalizedIntegral(&chib1P_sig, invm1S, fmchib1P-nSig*fsigma1, fmchib1P+nSig*fsigma1);
     rchib2P=fnchib2P*NormalizedIntegral(&chib2P_sig, invm1S, fmchib2P-nSig*fsigma2, fmchib2P+nSig*fsigma2);
     rchib3P=fnchib3P*NormalizedIntegral(&chib3P_sig, invm2S, fmchib3P-nSig*fsigma3, fmchib3P+nSig*fsigma3);
     rchib3P_1S=fnchib3P_1S*NormalizedIntegral(&chib3P_sig_1S, invm1S, fmchib3P_1S-nSig*fsigma3_1S, fmchib3P_1S+nSig*fsigma3_1S);
     rbgchib1P=fnbg1S*NormalizedIntegral(&background1S, invm1S, fmchib1P-nSig*fsigma1, fmchib1P+nSig*fsigma1);
     rbgchib2P=fnbg1S*NormalizedIntegral(&background1S, invm1S, fmchib2P-nSig*fsigma2, fmchib2P+nSig*fsigma2);
     rbgchib3P=fnbg2S*NormalizedIntegral(&background2S, invm2S, fmchib3P-nSig*fsigma3, fmchib3P+nSig*fsigma3);
     rbgchib3P_1S=fnbg1S*NormalizedIntegral(&background1S, invm1S, fmchib3P_1S-nSig*fsigma3_1S, fmchib3P_1S+nSig*fsigma3_1S);
     rratchib1P=rchib1P/rbgchib1P;
     rratchib2P=rchib2P/rbgchib2P;
     rratchib3P=rchib3P/rbgchib3P;
     rratchib3P_1S=rchib3P_1S/rbgchib3P_1S;
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

     char MasswindowChar[200];
     TH1F  *hGammaE_1P = new TH1F("hGammaE_1P","",nHistBins,0,10);
     TH1F  *hGammaE_2P = new TH1F("hGammaE_2P","",nHistBins,0,10);
     TH1F  *hGammaE_3P = new TH1F("hGammaE_3P","",nHistBins,0,10);
     TH1F  *hGammaE_3P_1S = new TH1F("hGammaE_3P_1S","",nHistBins,0,10);

     sprintf(MasswindowChar,"invm1S < %f && invm1S > %f",fmchib1P+nSig*fsigma1,fmchib1P-nSig*fsigma1);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_1P",MasswindowChar);
     sprintf(MasswindowChar,"invm1S < %f && invm1S > %f",fmchib2P+nSig*fsigma2,fmchib2P-nSig*fsigma2);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_2P",MasswindowChar);
     sprintf(MasswindowChar,"invm1S < %f && invm1S > %f",fmchib3P+nSig*fsigma3_1S,fmchib3P-nSig*fsigma3_1S);
     tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_3P_1S",MasswindowChar);
     sprintf(MasswindowChar,"invm2S < %f && invm2S > %f",fmchib3P+nSig*fsigma3,fmchib3P-nSig*fsigma3);
     tree2S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaE_3P",MasswindowChar);

     hGammaE_1P->Print();
     hGammaE_2P->Print();
     hGammaE_3P->Print();
     hGammaE_3P_1S->Print();

     char saveName[200];

     TCanvas* PhotonEnergyCanvas = new TCanvas("PhotonEnergyCanvas","PhotonEnergyCanvas",1600, 1200);
     PhotonEnergyCanvas->SetFillColor(kWhite);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_1P->GetXaxis()->SetTitle("E_{#gamma} in 1P peak (1S sample)");	hGammaE_1P->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_1P_1S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_2P->GetXaxis()->SetTitle("E_{#gamma} in 2P peak (1S sample)");	hGammaE_2P->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_2P_1S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_3P->GetXaxis()->SetTitle("E_{#gamma} in 3P peak (2S sample)");	hGammaE_3P->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_3P_2S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_3P_1S->GetXaxis()->SetTitle("E_{#gamma} in 3P peak (1S sample)");	hGammaE_3P_1S->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy3P_1S.pdf",FitID);
     if(SaveAll) PhotonEnergyCanvas->SaveAs(saveName);

     double plotFrom=invm_min1S;
     double plotTo=invm_max1S;

     double linewidth = 1;
     double chi21S;
     double chi22S;
     double binWidth=5;
     int FrameBins1S=(plotTo-plotFrom)*1000/binWidth;


     int FrameBins2Sin1S=100./4.*44./binWidth;
     int FrameBins2S=100./4.*64./binWidth;
     int FrameBins3Sin1S=100./4.*32./binWidth;

     char plotYtitle[200];

     if(SetSignalToZero){
     n1PnJ.setVal(1e-8);
     n2PnJ.setVal(1e-8);
     n3PnJ.setVal(1e-8);
     n2PnJ_2S.setVal(1e-8);
     n3PnJ_1S.setVal(1e-8);
     }

//     TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",1600, 800);
     TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",2100, 800);
    ChibCanvas1S->Divide(1);
    ChibCanvas1S->SetFillColor(kWhite);
    ChibCanvas1S->cd(1);
//    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    double FontSize=0.0325;
    data1S->Print();
    RooPlot* frame1S= invm1S.frame(FrameBins1S);
    frame1S->SetTitle("#chi_{c} invariant mass, J/#psi decay");
    data1S->plotOn(frame1S);
    sprintf(plotYtitle,"Events per %1.1f MeV",binWidth);
    frame1S->SetYTitle(plotYtitle);
        modelPdf1S.plotOn(frame1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));;
    chi21S = frame1S->chiSquare();
    /////////// Two CB adventure //////////////////
    modelPdf1S.plotOn(frame1S, Components("chib1P0"), LineStyle(2),LineColor(kMagenta),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib1P1"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib1P2"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("background1S"), LineStyle(2),LineColor(kBlack),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));

//    modelPdf1S.paramOn(frame1S, Layout(0.1875,0.42,0.875), Format("NE",AutoPrecision(1)));
    modelPdf1S.paramOn(frame1S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));

    frame1S->SetTitle(0);
    frame1S->SetMinimum(0);
//    frame1S->SetMaximum(2000);
    frame1S->Draw();

    double xText=3.65;
    double highestText=frame1S->GetMaximum();
    double deltaText=0.06;

    char text[200];

    	cout<<"DRAW LATEX"<<endl;
    sprintf(text,"PES = %1.4f +- %1.4f",PhotonMassScale.getVal(),PhotonMassScale.getError());
    TLatex text1 = TLatex(xText,highestText*(0.95-0*deltaText),text);
    text1.SetTextSize(FontSize);
    text1.Draw( "same" );
    sprintf(text,"#Delta M_{#chi_{c0}} = %1.3f MeV",deltaMJ0);
    TLatex text4 = TLatex(xText,highestText*(0.95-1*deltaText),text);
    text4.SetTextSize(FontSize);
    text4.Draw( "same" );
    sprintf(text,"#Delta M_{#chi_{c1}} = %1.3f MeV",deltaMJ1);
    TLatex text40 = TLatex(xText,highestText*(0.95-2*deltaText),text);
    text40.SetTextSize(FontSize);
    text40.Draw( "same" );
    sprintf(text,"#Delta M_{#chi_{c2}} = %1.3f MeV",deltaMJ2);
    TLatex text400 = TLatex(xText,highestText*(0.95-3*deltaText),text);
    text400.SetTextSize(FontSize);
    text400.Draw( "same" );
    sprintf(text,"#chi^{2} / ndf = %1.3f",chi21S);
    TLatex text4000 = TLatex(xText,highestText*(0.95-4.5*deltaText),text);
    text4000.SetTextSize(FontSize);
    text4000.Draw( "same" );


//    ChibCanvas1S->Modified();
        sprintf(cutName_,"vtxProb%d_ctCut%d_nYsig%d_cut_gammapt%d_cut_RconvMin%d_cut_RconvMax%d_cut_Ypt%d_vtxChi2ProbGamma%d_vtxChi2ProbGammaLog%d",int(1000000*cut_vtxProb),int(1000000*cut_ct),int(1000000*nYsig),int(1000000*cut_gammapt),int(1000000*cut_RconvMin),int(1000000*cut_RconvMax),int(1000000*cut_Ypt),int(100000000*cut_vtxChi2ProbGamma),-int(10000*cut_vtxChi2ProbGammaLog));
    sprintf(saveName,"Figures/%s/InvMass_Y1S_%s.pdf",FitID,cutName_);
    ChibCanvas1S->SaveAs(saveName);




    cout<<"deleting pointers"<<endl;

    delete data_______;
    delete data__;
    delete data___;
    delete data____;
    delete data_____;
    delete data______;
    delete data______vtx;

    delete data1S;
    delete data2S;
    delete data3S;
    delete dataSB1S;
    delete dataSB2S;

    delete PhotonEnergyCanvas;
    delete fQ;

	delete frame1S;

	delete ChibCanvas1S;


    cout<<"deleted pointers"<<endl;
    }

    sprintf(filename_,"Figures/%s/CutHistos_Y%dS.root",FitID,nState);
	TFile *CutHistos_file = new TFile(filename_,"RECREATE");

	char XTitle[200];
	if(nYsigCut  )sprintf(XTitle,"Y(1S) mass window, n#sigma_{m_{#mu#mu}}(|y|^{Y(1S)}) ");
	if(vtxProbCut)sprintf(XTitle,"Dimuon vertex probability cut");
	if(ctSigCut  )sprintf(XTitle,"c#tau significance cut, n#sigma");
	if(gammaptCut)sprintf(XTitle,"#gamma - p_{T} cut, GeV");
	if(ctCut     )sprintf(XTitle,"c#tau cut, cm");
	if(RConvUpperCut  )sprintf(XTitle,"RConvMax cut");
	if(RConvLowerCut  )sprintf(XTitle,"RConvMin cut");
	if(pTCut  )sprintf(XTitle,"lower p_{T}^{Y} cut, GeV");
	if(vtxChi2ProbGammaCut  )sprintf(XTitle,"vertexChi2ProbGamma cut");
	if(vtxChi2ProbGammaLogCut  )sprintf(XTitle,"log10(vertexChi2ProbGamma) cut");

	if(n_Cut>1.5){

	double yOffset=1.7;

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

    hCut_CBa->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_CBa->SetStats(0);hCut_CBa->GetYaxis()->SetTitle("#alpha_{CB}");	hCut_CBa->GetXaxis()->SetTitle(XTitle);	hCut_CBa->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/CBa.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_CBn->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_CBn->SetStats(0);hCut_CBn->GetYaxis()->SetTitle("n_{CB}");	hCut_CBn->GetXaxis()->SetTitle(XTitle);	hCut_CBn->Draw();
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

    hCut_PhotonEnergyScale->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonEnergyScale->SetStats(0);hCut_PhotonEnergyScale->GetYaxis()->SetTitle("Photon Energy Scale");	hCut_PhotonEnergyScale->GetXaxis()->SetTitle(XTitle);	hCut_PhotonEnergyScale->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/PhotonEnergyScale.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_PhotonEnergyScale2P->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonEnergyScale2P->SetStats(0);hCut_PhotonEnergyScale2P->GetYaxis()->SetTitle("Photon Energy Scale 2P");	hCut_PhotonEnergyScale2P->GetXaxis()->SetTitle(XTitle);	hCut_PhotonEnergyScale2P->Draw();
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

    hCut_sigLL->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_sigLL->SetStats(0);hCut_sigLL->GetYaxis()->SetTitle("combined LLratio-significance of #chi_{b}(1P) and #chi_{b}(2P)");	hCut_sigLL->GetXaxis()->SetTitle(XTitle);	hCut_sigLL->Draw();
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
	hCut_1PsigLL->Write();
	hCut_1PsigLL->Write();

	hCut_2Pnevt_2S->Write();
	hCut_sigLL->Write();
	hCut_2Psig_2S->Write();
	hCut_sigLLhalf->Write();
	hCut_3PsigLLhalf->Write();

	CutHistos_file->Write();
    CutHistos_file->Close();

	return 0;
}


