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
    bool useAnalyticalBKG=false;
    bool Search2P=false;
    bool BkgMixerBool=false;
    bool BkgToy=false;
	bool useExistingMixFile=false;

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
		    if(std::string(argv[i]).find("useAnalyticalBKG=true") != std::string::npos) {useAnalyticalBKG=true; cout<<"useAnalyticalBKG"<<endl;}
		    if(std::string(argv[i]).find("Search2P=true") != std::string::npos) {Search2P=true; cout<<"Search2P"<<endl;}
		    if(std::string(argv[i]).find("BkgMixer=true") != std::string::npos) {BkgMixerBool=true; cout<<"BkgMixerBool"<<endl;}
		    if(std::string(argv[i]).find("BkgToy=true") != std::string::npos) {BkgToy=true; cout<<"BkgToy"<<endl;}
		    if(std::string(argv[i]).find("useExistingMixFile=true") != std::string::npos) {useExistingMixFile=true; cout<<"useExistingMixFile"<<endl;}

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

    invm_min1S=3.2; invm_max1S=3.8;
    invm_min2S=3.1; invm_max2S=4;
    invm_min3S=3.1; invm_max3S=4;
    if(Search2P)invm_max1S=4.1;

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

    RooRealVar vtxNsigmadz   = RooRealVar("vtxNsigmadz", "vtxNsigmadz",-15,15);
    RooRealVar vtxdz   = RooRealVar("vtxdz", "vtxdz",-15,15);
    RooRealVar vtxerrdz   = RooRealVar("vtxerrdz", "vtxerrdz",0,15);

    RooRealVar UpsIndex   	= RooRealVar("UpsIndex", "UpsIndex",0,10000000);
	RooRealVar GammaIndex   	= RooRealVar("GammaIndex", "GammaIndex",0,500);

	RooRealVar RunNb     = RooRealVar("RunNb", "RunNb",0,1e9);
	RooRealVar EventNb     = RooRealVar("EventNb", "EventNb",0,1e15);

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



    for(int inCut=0; inCut<n_Cut; inCut++){

///////////////////////////////////////////////
//////////// CUTS /////////////////////////////
///////////////////////////////////////////////
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

/*	 RooArgSet CutVars("CutVars");
	 CutVars.add(invm1S);
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
*/
    char DScutChar[1000];
//    sprintf(DScutChar,"jpsipt>%f && gammapt>%f && gammapt<%f && Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && gammapt>%f*Q+%f && TMath::Abs(0.5*TMath::Log((sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)+gammapz)/(sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)-gammapz)))<%f",cut_Ypt,cut_gammapt,cut_gammapt_max,cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,gammaptK,gammaptD,cut_gammaeta);
    sprintf(DScutChar,"jpsipt>%f && gammapt>%f && gammapt<%f && Rconv > %f && Rconv < %f && jpsieta < %f && jpsieta > %f && ctpv < %f &&  ctpv > %f && vertexChi2ProbGamma > %f && log10(vertexChi2ProbGamma) > %f && jpsiVprob > %f && gammapt>%f*Q+%f && (gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>%f && TMath::Abs(vtxNsigmadz) < %f &&TMath::Abs(vtxdz)<%f &&TMath::Abs(ctpvsig)<%f",cut_Ypt,cut_gammapt,cut_gammapt_max,cut_RconvMin,cut_RconvMax,cut_rap,-cut_rap,cut_ct,-cut_ct,cut_vtxChi2ProbGamma,cut_vtxChi2ProbGammaLog,cut_vtxProb,gammaptK,gammaptD,MINcosalphaCut,cut_dzSig, cut_dz, cut_ctSig);
   RooDataSet* dataAfterBasicCuts=(RooDataSet*)data_->reduce(SelectVars(RooArgSet(invm1S,invm2S,invm3S,Y1Smass_nSigma,Y2Smass_nSigma,Y3Smass_nSigma,jpsimass)),Cut(DScutChar));
    cout<<DScutChar<<endl;cout<<endl;
    dataAfterBasicCuts->Print();

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
        	nSigSBlow=3.0;
        	nSigSBhigh=3.2;
        	UpperSBmax=4;

        	LowerSBmin2S=8.7;
        	nSigSBlow2S=-3.5;
        	nSigSBhigh2S=3.5;
        	UpperSBmax2S=10.85;
        }

        if(useLeftSB&&!useRightSB){
        	LowerSBmin=2.8;
        	nSigSBlow=3.0;
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
        	nSigSBhigh=3.2;
        	UpperSBmax=4.;

        	LowerSBmin2S=1000;
        	nSigSBlow2S=-1000;
        	nSigSBhigh2S=3.5;
        	UpperSBmax2S=10.85;
        }

        double MassSignalMin=3;
        double MassSignalMax=3.2;

    sprintf(DScutChar,"jpsimass > %f && jpsimass < %f || jpsimass < %f && jpsimass > %f",LowerSBmin,nSigSBlow,UpperSBmax,nSigSBhigh);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB1S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    dataSB1S->Print();



/////// Producing Signal dataset
    sprintf(DScutChar,"jpsimass > %f && jpsimass < %f",MassSignalMin,MassSignalMax);
    cout<<DScutChar<<endl;
    RooDataSet* data1S=(RooDataSet*)dataAfterBasicCuts->reduce(Cut(DScutChar));
    data1S->Print();



    if(SetSignalToZero&&useSBforBKGmodel){
    	data1S=dataSB1S;
    }
    TTree* tree1S=(TTree*)data1S->tree();
    TTree* treeAllCuts=(TTree*)dataAfterBasicCuts->tree();

    if(useSBforBKGmodel){
    	tree1S=(TTree*)dataSB1S->tree();
    }

    if(SetSignalToZero){
    	data1S=dataSB1S;
    }

    char treeName[200];
    sprintf(treeName,"Figures/%s/tree_Y1S.root",FitID);
    if(SaveAll) tree1S->SaveAs(treeName);
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


    }

    /////////////////////////////////////////////
    ///////////// THE Background Saga ///////////
    /////////////////////////////////////////////

    RooAbsPdf* background1S;

    double bb_thresh=3.585;
    double chib1Pmin=3.325;
    double chic1min=3.25;

        int nHistBins=10000;
        float m_gamma = 0.;
        char DrawChar[500];

//        sprintf(DrawChar,"invm1S>0 && Q<1000");//chic1min
        sprintf(DrawChar,"invm1S<%f || invm1S>%f",chib1Pmin,bb_thresh);//chic1min
        if(useSBforBKGmodel) sprintf(DrawChar,"invm1S<10000");

        TH1F  *hYmass1S = new TH1F("hYmass1S","",nHistBins,2,5);
        TH1F  *hYmass_check1S = new TH1F("hYmass_check1S","",nHistBins,2,5);

       tree1S->Draw("jpsimass>>hYmass1S",DrawChar);

        TH1F  *hGammaP1S = new TH1F("hGammaP1S","",nHistBins,0,100);
        TH1F  *hGammaP_check1S = new TH1F("hGammaP_check1S","",nHistBins,0,100);

        tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP1S",DrawChar);

        TH1F  *hUpsP1S = new TH1F("hUpsP1S","",nHistBins,0,100);
        TH1F  *hUpsP_check1S = new TH1F("hUpsP_check1S","",nHistBins,0,100);

        tree1S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP1S",DrawChar);

        TH1F  *hCosAlphaP1S = new TH1F("hCosAlphaP1S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check1S = new TH1F("hCosAlphaP_check1S","",nHistBins,-1,1);

        tree1S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP1S",DrawChar);


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


    	double p_gammaInt=0;
        double p_gammaArray[nToy];
        double p_UpsInt=0;
        double p_UpsArray[nToy];
        double cosAlphaInt=0;
        double cosAlphaArray[nToy];


        float M1S;

        tQ->Branch("invm1S",&M1S,"invm1S/F");

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
        	  p_gammaInt+=p_gamma1S; p_gammaArray[i]=p_gamma1S;
        	  p_UpsInt+=p_Ups1S; p_UpsArray[i]=p_Ups1S;
        	  cosAlphaInt+=cosAlpha1S; cosAlphaArray[i]=cosAlpha1S;
          }

          else i--;

        }



        hCosAlphaP_check1S->Write();
        hUpsP_check1S->Write();
        hGammaP_check1S->Write();
        hYmass_check1S->Write();

        double avGamma;
        double avUps;
        double avAlpha;

        double avGammUps=0;
        double avGammAlpha=0;
        double avAlphaUps=0;

        for(int i=0;i<nToy;i++){
        avGamma=(p_gammaArray[i]-p_gammaInt/nToy);
        avUps=(p_UpsArray[i]-p_UpsInt/nToy);
        avAlpha=(cosAlphaArray[i]-cosAlphaInt/nToy);

        avGammUps+=avGamma*avUps;
        avGammAlpha+=avGamma*avAlpha;
        avAlphaUps+=avAlpha*avUps;
        }

        avGammUps=avGammUps/nToy;
        avGammAlpha=avGammAlpha/nToy;
        avAlphaUps=avAlphaUps/nToy;

        cout<<"avGammUps = "<<avGammUps/hGammaP1S->GetRMS()/hUpsP1S->GetRMS()<<endl;
        cout<<"avGammAlpha = "<<avGammAlpha/hGammaP1S->GetRMS()/hCosAlphaP1S->GetRMS()<<endl;
        cout<<"avAlphaUps = "<<avAlphaUps/hCosAlphaP1S->GetRMS()/hUpsP1S->GetRMS()<<endl;




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


        }





        TH1F  *BkgToyHist1S = new TH1F("BkgToyHist1S","",100,3,4.5);
        tQ->Draw("invm1S>>BkgToyHist1S");


        invm1S.setBins(50);
        invm2S.setBins(50);

        RooDataSet* RooBkgToySet = new RooDataSet("RooBkgToySet","RooBkgToySet",tQ,RooArgList(invm1S,invm2S));
        RooBkgToySet->Print();
        RooDataHist* RooBkgToyHist1S = new RooDataHist("RooBkgToyHist1S","RooBkgToyHist1S",RooArgList(invm1S),*RooBkgToySet);
        RooBkgToyHist1S->Print();
        BkgToyHist1S->Write();

        fQ->Write();
        fQ->Close();

        double q01S_Start=3.2;
        RooRealVar alpha1("alpha1","alpha1",0.6,0.,5.);
        RooRealVar beta1("beta1","beta1",-2.5,-10.,0.);
        RooRealVar q01S("q01S","q01S",q01S_Start,3.15,3.25);
        RooFormulaVar a1("a1","TMath::Abs(@0-@1)",RooArgList(invm1S,q01S));
        RooFormulaVar b1("b1","@0*(@1-@2)",RooArgList(beta1,invm1S,q01S));
        RooFormulaVar signum1("signum1","(TMath::Sign(-1.,@0-@1)+1)/2.",RooArgList(invm1S,q01S));

        if(BkgToy) background1S = new RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S,1);
        if(useAnalyticalBKG) background1S = new RooGenericPdf("background1S","signum1*pow(a1,alpha1)*exp(b1)",RooArgSet(signum1,a1,alpha1,b1));



        if(BkgMixerBool){
        	cout<<"Using BkgMixer tool to build background shape"<<endl;
     	    gROOT->ProcessLine(".L BkgMixer.C++");
     	    gROOT->ProcessLine(".L BkgMixer.C++");
     	    gROOT->ProcessLine(".L BkgMixer.C++");
     	    char MixerInputTree[200];
     	    int nStateToMix;
     	    char ExcludeSignal[200];
     	    int nBinsCosAlphaTry;

     	    //change accordingly (feed the macro with the tree of this iteration...)
     	    nStateToMix=1;
     	    nBinsCosAlphaTry=2500;
     	    sprintf(ExcludeSignal,"invm1S<%f || invm1S>%f",chib1Pmin,bb_thresh);
     	    sprintf(MixerInputTree,"BkgMixTree/tree_Jpsi_Pt0____________________________.root");
     	    BkgMixer(MixerInputTree,tree1S,nStateToMix,invm_min1S,invm_max1S,useExistingMixFile,nMix,dirstruct,ExcludeSignal,cut_gammapt,gammaptK,gammaptD,Ymass1S,nBinsCosAlphaTry);
     	    TH1F *hInvM1S=(TH1F*)hInvMnS->Clone("hInvM1S");
     	    hInvMnS->Print();
     	    hInvM1S->Print();

          RooDataHist* RooBkgToyHist1S_Mixer = new RooDataHist("RooBkgToyHist1S_Mixer","RooBkgToyHist1S_Mixer",RooArgList(invm1S),hInvM1S);
          RooBkgToyHist1S_Mixer->Print();
          background1S = new  RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S_Mixer,1);


          fMix->Close();

        }


//        gSystem->Unlink(bkgToyname);


        delete hYmass1S;
        delete hGammaP1S;
        delete hUpsP1S;
        delete hCosAlphaP1S;
        delete hCosAlphaP_check1S;
        delete hUpsP_check1S;
        delete hGammaP_check1S;
        delete hYmass_check1S;



//declare fit variables

    RooRealVar m_chib1P= RooRealVar("Mean #chi_{b}(1P)","Mean #chi_{b}(1P)",9.888,9.85,9.925);



/////////// TwoCB adventure //////////////////

    char MassScaleFormula[200];
    char SigmaScaleFormula[200];

    double alpha_3J_START=0.583;//0.74;
    double n_3J_START=3.;
    double PhotonMassScale_START=1.;
    double PhotonMassScale2P_START=1.;
    double PhotonSigmaScale_START=0.017;

    RooRealVar alpha_3J_J  = RooRealVar("#alpha","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar alpha_3J_J1  = RooRealVar("#alpha_{J1}","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar alpha_3J_J2  = RooRealVar("#alpha_{J2}","CB_{#alpha}",alpha_3J_START,0.,2.);
    RooRealVar n_3J_J      = RooRealVar("n","CB_{n}",n_3J_START,.5,15.);
    RooRealVar n_3J_J1      = RooRealVar("n_{J1}","CB_{n}",n_3J_START,.5,15.);
    RooRealVar n_3J_J2      = RooRealVar("n_{J2}","CB_{n}",n_3J_START,.5,15.);
    RooRealVar PhotonMassScale= RooRealVar("PES","PES",PhotonMassScale_START,0.95,1.05);
    RooRealVar PhotonMassScale2P= RooRealVar("PES_{J2}","PhotonMassScale2P",PhotonMassScale2P_START,0.95,1.05);
    RooRealVar PhotonSigmaScale= RooRealVar("#sigma_{Q}/Q","PhotonSigmaScale",PhotonSigmaScale_START,0.01,0.025);
    RooRealVar PhotonSigmaScaleJ2= RooRealVar("#sigma_{Q}/Q, J2","PhotonSigmaScale",PhotonSigmaScale_START,0.01,0.025);

    RooRealVar n_3J      = RooRealVar("CB_{n}","CB_{n}",n_3J_START,.5,15.);
    RooRealVar alpha_3J  = RooRealVar("CB_{#alpha}","CB_{#alpha}",alpha_3J_START,0.,2.);

    alpha_3J.setVal(alpha_3J_START);
    n_3J.setVal(n_3J_START);
    PhotonMassScale.setVal(PhotonMassScale_START);
    PhotonMassScale2P.setVal(PhotonMassScale2P_START);
    PhotonSigmaScale.setVal(PhotonSigmaScale_START);

//    n_3J.setConstant();
//    n_3J_J1.setConstant();
 //   n_3J_J2.setConstant();
    //    alpha_3J.setConstant();
//    alpha_3J_J1.setConstant();
//    alpha_3J_J2.setConstant();

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
	RooFormulaVar m_chib1P2float("m_chib1P2float", MassScaleFormula, RooArgList(m_chib1P2fix, PhotonMassScale2P));

    RooFormulaVar w_chib1P0float("w_chib1P0float", SigmaScaleFormula, RooArgList(m_chib1P0float, PhotonSigmaScale));
    RooFormulaVar w_chib1P1float("w_chib1P1float", SigmaScaleFormula, RooArgList(m_chib1P1float, PhotonSigmaScale));
    RooFormulaVar w_chib1P2float("w_chib1P2float", SigmaScaleFormula, RooArgList(m_chib1P2float, PhotonSigmaScaleJ2));//J2

	RooCBShape chib1P0      = RooCBShape ("chib1P0","chib1P0",invm1S,m_chib1P0float,w_chib1P0float,alpha_3J_J1,n_3J_J1);
	RooCBShape chib1P1      = RooCBShape ("chib1P1","chib1P1",invm1S,m_chib1P1float,w_chib1P1float,alpha_3J_J1,n_3J_J1);
    RooCBShape chib1P2      = RooCBShape ("chib1P2","chib1P2",invm1S,m_chib1P2float,w_chib1P2float,alpha_3J_J2,n_3J_J1);

    double n1PnJ_START=12000;
    RooRealVar n1PnJ= RooRealVar("N_{#chi_{cJ}}","#Sigma_{J}(#chi_{cJ})",n1PnJ_START,1000,20000);
    n1PnJ.setVal(n1PnJ_START);

    RooAddPdf chib1P_sigInd= RooAddPdf("chib1P_sigInd","chib1P_sigInd",RooArgList(chib1P0,chib1P1,chib1P2),RooArgList(fractionJ0,fractionJ1));

    RooFormulaVar chib1P_nevt("chib1P_nevt", "@1*@0+@2*@0+(1-@1-@2)*@0", RooArgList(n1PnJ,fractionJ0,fractionJ1));
    RooAddPdf chib1P_sig= RooAddPdf("chib1P_sig","chib1P_sig",RooArgList(chib1P_sigInd),RooArgList(chib1P_nevt));





////////////////////////////////////////////////////////


    RooRealVar background_nevt1S = RooRealVar("N_{bkg}1S","N_{bkg}1S",1500,0,20000);
    RooRealVar background_nevtSB1S = RooRealVar("background_nevtSB1S","background_nevtSB1S",1500,0,10000);
    RooAddPdf backgroundSB1S= RooAddPdf("backgroundSB1S","backgroundSB1S",RooArgList(*background1S),RooArgList(background_nevtSB1S));
    RooAddPdf backgroundE1S= RooAddPdf("backgroundE1S","backgroundE1S",RooArgList(*background1S),RooArgList(background_nevt1S));


    sprintf(DScutChar,"invm1S > %f && invm1S < %f",invm_min1S,invm_max1S);
    RooDataSet* data1S_masswindow=(RooDataSet*)data1S->reduce(Cut(DScutChar));



    RooAddPdf modelPdf1S= RooAddPdf("modelPdf1S","modelPdf1S",RooArgList(chib1P_sig,*background1S),RooArgList(chib1P_nevt,background_nevt1S));
    RooAddPdf modelPdf1S_= RooAddPdf("modelPdf1S_","modelPdf1S_",RooArgList(chib1P_sig,*background1S),RooArgList(chib1P_nevt,background_nevt1S));






    invm1S.setRange("SBregion31S",bb_thresh,invm_max1S);

    invm1S.setRange("Chib1Pregion",invm_min1S,bb_thresh);

    invm1S.setRange("TotalBlindedRegion2",bb_thresh,invm_max1S);

    invm1S.setRange("FullRegion",invm_min1S,invm_max1S);


    char SBregion2Char[200];
    sprintf(SBregion2Char,"invm1S > %f && invm1S < %f",bb_thresh,invm_max1S);
    char SBregion31SChar[200];
    sprintf(SBregion31SChar,"invm1S > %f && invm1S < %f",invm_min1S,invm_min1S);//chic1min);

    double Nbkg_1S_SB1=data1S->sumEntries(SBregion2Char);
    double Nbkg_1S_SB2=data1S->sumEntries(SBregion31SChar);
    cout<<"background_nevt1S sumentriesrange = " << Nbkg_1S_SB1<<endl;
    cout<<"background_nevt1S sumentriesrange2 = " << Nbkg_1S_SB2<<endl;


    double BkgSBInt1S;
    BkgSBInt1S = NormalizedIntegral(&backgroundSB1S, invm1S, bb_thresh,invm_max1S);
    cout<<BkgSBInt1S<<endl;
    BkgSBInt1S = BkgSBInt1S+ NormalizedIntegral(&backgroundSB1S, invm1S, invm_min1S,invm_min1S);//chic1min);
    cout<<BkgSBInt1S<<endl;
    background_nevt1S.setVal((Nbkg_1S_SB1+Nbkg_1S_SB2)/BkgSBInt1S);
    background_nevt1S.setConstant();
    cout<<"background normalization 1S = " << background_nevt1S.getVal()<<endl;


    double covQual_BG1S=0;
    double covQual_1P2P1S=0;



    char SignormRegion[200];

    if(useAnalyticalBKG) sprintf(SignormRegion,"FullRegion");
    if(!useAnalyticalBKG) sprintf(SignormRegion,"Chib1Pregion");

    if(useAnalyticalBKG) background_nevt1S.setConstant(kFALSE);
    cout<<"fitting in region: "<<SignormRegion<<endl;


        RooFitResult* Fit1P2P = modelPdf1S_.fitTo(*data1S,Save(1),Range(SignormRegion),PrintEvalErrors(-1),PrintLevel(1),Warnings(kFALSE));
        Fit1P2P->Print();
        covQual_1P2P1S = Fit1P2P->covQual();


        ratio_J2overJ1.setConstant();


        background_nevt1S.setConstant();

        delete Fit1P2P;


cout<<"deltaM 1P, J2: "<<m_chib1P2float.getVal()-m_chib1P2fix.getVal()<<endl;
cout<<"deltaM 1P, J1: "<<m_chib1P1float.getVal()-m_chib1P1fix.getVal()<<endl;


double deltaMJ0=(m_chib1P0float.getVal()-m_chib1P0fix.getVal())*1000;
double deltaMJ1=(m_chib1P1float.getVal()-m_chib1P1fix.getVal())*1000;
double deltaMJ2=(m_chib1P2float.getVal()-m_chib1P2fix.getVal())*1000;

m_chib1P.setVal(m_chib1P1float.getVal()*fracJ1.getVal()+m_chib1P2float.getVal()*(1-fracJ1.getVal()));
cout<<"m_chib1P: "<<m_chib1P.getVal()<<endl;




        cout<<"chib1P_nevt: "<<chib1P_nevt.getVal()<<endl;


     double TotalRealEvents1S = data1S_masswindow->sumEntries();
     double totalEventsInFit1S;

     totalEventsInFit1S = chib1P_nevt.getVal()+background_nevt1S.getVal();

     cout<< "N_tot_Fit1S =                            "<<totalEventsInFit1S <<endl;
     cout<< "N_tot_Samlpe1S =                         "<<TotalRealEvents1S <<endl;
     cout<< "Delta_N_tot1S = N_tot_Sample1S-N_tot_Fit1S = "<<TotalRealEvents1S-totalEventsInFit1S <<endl;


// likelihood gaussian check

     double nSig=1.5;

     double fmchib1P = m_chib1P.getVal();


     double OneSigmaCL = 0.682689492137;

     double sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib1P_sig, invm1S, fmchib1P-sigmaCalc, fmchib1P+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma1=sigmaCalc;




     cout<<"m_chib1P: "<<m_chib1P.getVal()<<endl;
     cout<<"fsigma1: "<<fsigma1<<endl;

     double fnchib1P = chib1P_nevt.getVal();
     double fnbg1S = background_nevt1S.getVal();
     double rchib1P,rchib2P,rchib3P,rbgchib1P,rbgchib2P,rbgchib3P,rsigchib1P,rsigchib2P,rsigchib3P,rchib3P_1S,rbgchib3P_1S,rsigchib3P_BityukovKrasnikov,rsigchib3P_BityukovKrasnikov_1S,rsigchib3P_ScL,rsigchib3P_ScL_1S;
     double rratchib1P,rratchib2P,rratchib3P,rratchib3P_1S,rsigchib3P_1S;
     rchib1P=fnchib1P*NormalizedIntegral(&chib1P_sig, invm1S, fmchib1P-nSig*fsigma1, fmchib1P+nSig*fsigma1);
     rbgchib1P=fnbg1S*NormalizedIntegral(&*background1S, invm1S, fmchib1P-nSig*fsigma1, fmchib1P+nSig*fsigma1);
     rratchib1P=rchib1P/rbgchib1P;
     rsigchib1P=rchib1P/sqrt(rchib1P+rbgchib1P);




     double plotFrom=invm_min1S;
     double plotTo=invm_max1S;
     if(Search2P) plotFrom=3.58;

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
      }

//     TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",1600, 800);
     TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",2100, 800);
    ChibCanvas1S->Divide(1);
    ChibCanvas1S->SetFillColor(kWhite);
    ChibCanvas1S->cd(1);
    if(Search2P) ChibCanvas1S->cd(1)->SetLogy(1);

//    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.3);/////
    gPad->SetFillColor(kWhite);

    double FontSize=0.0325;
    data1S->Print();
    RooPlot* frame1S= invm1S.frame(FrameBins1S);
    frame1S->SetTitle("#chi_{c} invariant mass, J/#psi decay");
    frame1S->GetXaxis()->SetLimits(plotFrom,plotTo);
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
    bool restrictMaximum=false;
    if(Search2P) restrictMaximum=true;
    if(restrictMaximum) frame1S->SetMaximum(150);
    frame1S->Draw();

    double xText=3.68;
    double highestText=frame1S->GetMaximum();
    double deltaText=0.06;

    char text[200];

    	cout<<"DRAW LATEX"<<endl;
    sprintf(text,"PES = %1.4f +- %1.4f",PhotonMassScale.getVal(),PhotonMassScale.getError());
    TLatex text1 = TLatex(xText,highestText*(0.95-0*deltaText),text);
    text1.SetTextSize(FontSize);
//    text1.Draw( "same" );
    sprintf(text,"#Delta M_{#chi_{c0}} = %1.3f MeV",deltaMJ0);
    TLatex text4 = TLatex(xText,highestText*(0.95-1*deltaText),text);
    text4.SetTextSize(FontSize);
//    text4.Draw( "same" );
    sprintf(text,"#Delta M_{#chi_{c1}} = %1.3f MeV",deltaMJ1);
    TLatex text40 = TLatex(xText,highestText*(0.95-2*deltaText),text);
    text40.SetTextSize(FontSize);
//    text40.Draw( "same" );
    sprintf(text,"#Delta M_{#chi_{c2}} = %1.3f MeV",deltaMJ2);
    TLatex text400 = TLatex(xText,highestText*(0.95-3*deltaText),text);
    text400.SetTextSize(FontSize);
//    text400.Draw( "same" );
    sprintf(text,"#chi^{2} / ndf = %1.3f",chi21S);
    TLatex text4000 = TLatex(xText,highestText*(0.95-0*deltaText),text);
    text4000.SetTextSize(FontSize);
    text4000.Draw( "same" );


//    ChibCanvas1S->Modified();
    char saveName[2000];
    sprintf(cutName_,"vtxProb%d_ctCut%d_nYsig%d_cut_gammapt%d_cut_RconvMin%d_cut_RconvMax%d_cut_Ypt%d_vtxChi2ProbGamma%d_vtxChi2ProbGammaLog%d",int(1000000*cut_vtxProb),int(1000000*cut_ct),int(1000000*nYsig),int(1000000*cut_gammapt),int(1000000*cut_RconvMin),int(1000000*cut_RconvMax),int(1000000*cut_Ypt),int(100000000*cut_vtxChi2ProbGamma),-int(10000*cut_vtxChi2ProbGammaLog));
    sprintf(saveName,"Figures/%s/InvMass_Y1S_%s.pdf",FitID,cutName_);
    ChibCanvas1S->SaveAs(saveName);
//    sprintf(saveName,"Figures/%s/InvMass_Y1S_Logy_%s.pdf",FitID,cutName_);
//    ChibCanvas1S->SaveAs(saveName);




    cout<<"deleting pointers"<<endl;

    delete dataAfterBasicCuts;

    delete data1S;
    delete dataSB1S;

    delete fQ;

	delete frame1S;

	delete ChibCanvas1S;


    cout<<"deleted pointers"<<endl;
    }



	return 0;
}


