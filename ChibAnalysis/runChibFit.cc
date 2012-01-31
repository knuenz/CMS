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
  	char cutName_[200];
    int nToy=5000000;

  	bool nYsigCut=false;
  	bool vtxProbCut=false;
  	bool ctSigCut=false;
  	bool gammaptCut=false;
  	bool ctCut=false;
  	bool RConvCut=false;
  	bool pTCut=false;
  	bool vtxChi2ProbGamma=false;

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
		    if(std::string(argv[i]).find("RConvCut") != std::string::npos) {RConvCut=true; cout<<"Optimizing RConvCut"<<endl;}
		    if(std::string(argv[i]).find("pTCut") != std::string::npos) {pTCut=true; cout<<"Optimizing pTCut"<<endl;}
		    if(std::string(argv[i]).find("vtxChi2ProbGamma") != std::string::npos) {vtxChi2ProbGamma=true; cout<<"Optimizing vtxChi2ProbGamma"<<endl;}

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
    Ymass1S=9.460;
    Ymass2S=10.023;

    double invm_min1S;
    double invm_max1S;
    double invm_min2S;
    double invm_max2S;
    double invm_min3S;
    double invm_max3S;

    invm_min1S=9.5; invm_max1S=11.1;
    invm_min2S=10; invm_max2S=11.1;
    invm_min3S=10.3; invm_max3S=11.1;

    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(1S)^{PDG}}");
    RooRealVar invm1S("invm1S",invmName,invm_min1S,invm_max1S);
    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(2S)^{PDG}}");
    RooRealVar invm2S("invm2S",invmName,invm_min2S,invm_max2S);
    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(3S)^{PDG}}");
    RooRealVar invm3S("invm3S",invmName,invm_min3S,invm_max3S);

    invm1S.setRange("TotalRange1S",invm_min1S,invm_max1S);
    invm2S.setRange("TotalRange2S",invm_min2S,invm_max2S);
    invm3S.setRange("TotalRange3S",invm_min3S,invm_max3S);

    RooRealVar jpsipt   = RooRealVar("jpsipt", "jpsipt",0,100);
    RooRealVar jpsimass   = RooRealVar("jpsimass", "jpsimass",8,12);
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
// RConvCut
    if(RConvCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0;
        deltaCut=0.07;
    }
// pT
    if(pTCut){
        n_Cut=n_Cut_Default;
        n_Cut_0=0;
        deltaCut=0.2;
    }
// pT
    if(vtxChi2ProbGamma){
        n_Cut=n_Cut_Default;
        n_Cut_0= 0.0000001;
        deltaCut=0.0000005;
    }

    TH1F  *hCut_1Psig = new TH1F("hCut_1Psig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1PSB = new TH1F("hCut_1PSB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1PS = new TH1F("hCut_1PS","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1PB = new TH1F("hCut_1PB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Psig = new TH1F("hCut_2Psig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PSB = new TH1F("hCut_2PSB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PS = new TH1F("hCut_2PS","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PB = new TH1F("hCut_2PB","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_sig = new TH1F("hCut_sig","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_covQual1 = new TH1F("hCut_covQual1","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

    TH1F  *hCut_1Pmean = new TH1F("hCut_1Pmean","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1Pwidth = new TH1F("hCut_1Pwidth","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_1Pnevt = new TH1F("hCut_1Pnevt","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pmean = new TH1F("hCut_2Pmean","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pwidth = new TH1F("hCut_2Pwidth","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2Pnevt = new TH1F("hCut_2Pnevt","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
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

    TH1F  *hCut_1PsigLL = new TH1F("hCut_1PsigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);
    TH1F  *hCut_2PsigLL = new TH1F("hCut_2PsigLL","",n_Cut,n_Cut_0,n_Cut_0+deltaCut*n_Cut);

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
    double cut_Ypt=0;

// Gamma cuts

    double cut_gammapt = 1.1;
    double cut_RconvMin = 2;
    double cut_RconvMax = 12;
    double cut_vtxChi2ProbGamma = 0.;//000005;

// Y - Gamma cuts



///////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////

    if(nYsigCut)    nYsig=n_Cut_0+inCut*deltaCut;
    if(vtxProbCut)    cut_vtxProb=n_Cut_0+inCut*deltaCut;
    if(ctSigCut)    cut_ctSig=n_Cut_0+inCut*deltaCut;
    if(gammaptCut)    cut_gammapt=n_Cut_0+inCut*deltaCut;
    if(ctCut)    	cut_ct=n_Cut_0+inCut*deltaCut;
    if(RConvCut)	{cut_RconvMin=n_Cut_0+inCut*deltaCut; cut_RconvMax=n_Cut_0+14-inCut*deltaCut;}
    if(pTCut)    	cut_Ypt=n_Cut_0+inCut*deltaCut;

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

    sprintf(DScutChar,"vertexChi2ProbGamma > %f",cut_vtxChi2ProbGamma);
    cout<<DScutChar<<endl;
    RooDataSet* data______vtx=(RooDataSet*)data_____->reduce(Cut(DScutChar));
    data______vtx->Print();

    sprintf(DScutChar,"jpsiVprob > %f",cut_vtxProb);
    cout<<DScutChar<<endl;
    RooDataSet* data______=(RooDataSet*)data______vtx->reduce(Cut(DScutChar));
    data______->Print();
/////// Producing SB dataset

        double LowerSBmin=8.7;
        double nSigSBlow=-3.5;
        double nSigSBhigh=300.5;
        double UpperSBmax=9.8;




    sprintf(DScutChar,"jpsimass > %f && Y1Smass_nSigma < %f || jpsimass < %f && Y1Smass_nSigma > %f",LowerSBmin,nSigSBlow,UpperSBmax,nSigSBhigh);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB1S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    dataSB1S->Print();

    sprintf(DScutChar,"jpsimass > %f && Y2Smass_nSigma < %f || jpsimass < %f && Y2Smass_nSigma > %f",LowerSBmin,nSigSBlow,UpperSBmax,nSigSBhigh);
    cout<<DScutChar<<endl;
    RooDataSet* dataSB2S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    dataSB2S->Print();


/////// Producing Signal dataset
    sprintf(DScutChar,"Y1Smass_nSigma > %f && Y1Smass_nSigma < %f ",-nYsig,nYsig);
    cout<<DScutChar<<endl;
    RooDataSet* data1S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    data1S->Print();

    sprintf(DScutChar,"Y2Smass_nSigma > %f && Y2Smass_nSigma < %f ",-nYsig,nYsig);
    cout<<DScutChar<<endl;
    RooDataSet* data2S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    data2S->Print();

    sprintf(DScutChar,"Y3Smass_nSigma > %f && Y3Smass_nSigma < %f ",-nYsig,nYsig);
    cout<<DScutChar<<endl;
    RooDataSet* data3S=(RooDataSet*)data______->reduce(Cut(DScutChar));
    data3S->Print();

    TTree* tree1S=(TTree*)data1S->tree();
    TTree* tree2S=(TTree*)data2S->tree();

    char treeName[200];
    sprintf(treeName,"tree_Y1S.root");
//    tree1S->SaveAs(treeName);
    sprintf(treeName,"tree_Y2S.root");
//    tree2S->SaveAs(treeName);


    // Ernest's ToyMC for background estimation:


        int nHistBins=100;
        float m_gamma = 0.;

        TH1F  *hYmass1S = new TH1F("hYmass1S","",nHistBins,8.7,11.3);
        TH1F  *hYmass_check1S = new TH1F("hYmass_check1S","",nHistBins,8.7,11.3);
        TH1F  *hYmass2S = new TH1F("hYmass2S","",nHistBins,8.7,11.3);
        TH1F  *hYmass_check2S = new TH1F("hYmass_check2S","",nHistBins,8.7,11.3);

       tree1S->Draw("jpsimass>>hYmass1S","invm1S<9.85 | invm1S>9.95&&invm1S<10.15|invm1S>10.3");
       tree2S->Draw("jpsimass>>hYmass2S");

        TH1F  *hGammaP1S = new TH1F("hGammaP1S","",nHistBins,0,10);
        TH1F  *hGammaP_check1S = new TH1F("hGammaP_check1S","",nHistBins,0,10);
        TH1F  *hGammaP2S = new TH1F("hGammaP2S","",nHistBins,0,10);
        TH1F  *hGammaP_check2S = new TH1F("hGammaP_check2S","",nHistBins,0,10);

        tree1S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP1S","invm1S<9.85 | invm1S>9.95&&invm1S<10.15|invm1S>10.3");
        tree2S->Draw("sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)>>hGammaP2S");

        TH1F  *hUpsP1S = new TH1F("hUpsP1S","",nHistBins,0,100);
        TH1F  *hUpsP_check1S = new TH1F("hUpsP_check1S","",nHistBins,0,100);
        TH1F  *hUpsP2S = new TH1F("hUpsP2S","",nHistBins,0,100);
        TH1F  *hUpsP_check2S = new TH1F("hUpsP_check2S","",nHistBins,0,100);

        tree1S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP1S","invm1S<9.85 | invm1S>9.95&&invm1S<10.15|invm1S>10.3");
        tree2S->Draw("sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hUpsP2S");

        TH1F  *hCosAlphaP1S = new TH1F("hCosAlphaP1S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check1S = new TH1F("hCosAlphaP_check1S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP2S = new TH1F("hCosAlphaP2S","",nHistBins,-1,1);
        TH1F  *hCosAlphaP_check2S = new TH1F("hCosAlphaP_check2S","",nHistBins,-1,1);

        tree1S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP1S","invm1S<9.85 | invm1S>9.95&&invm1S<10.15|invm1S>10.3");
        tree2S->Draw("(gammapx*jpsipx+gammapy*jpsipy+gammapz*jpsipz)/sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)/sqrt(jpsipx*jpsipx+jpsipy*jpsipy+jpsipz*jpsipz)>>hCosAlphaP2S");


        bool saveBkgToy=false;

        char bkgToyname[200];
        sprintf(bkgToyname,"%s/bkgToy.root",dirstruct);
        TFile *fQ = new TFile(bkgToyname,"RECREATE");

        bool existingBKGfile=true;
        if(fQ->Get("tQ")==NULL) existingBKGfile=false;

//        TTree tQ("tQ","toy MC");
        TTree* tQ = new TTree("tQ","toy MC");

        existingBKGfile=false; //// temporary, if you want to use existing TTree, also change RECREATE->UPDATE

        if(existingBKGfile) {
        	cout<<"using existing bkt-toyMC file"<<endl;
        	tQ=(TTree*)fQ->Get("tQ");
        }

        if(!existingBKGfile){
        	cout<<"producing new bkt-toyMC file"<<endl;

        hYmass1S->Write();
        hGammaP1S->Write();
        hUpsP1S->Write();
        hCosAlphaP1S->Write();
        hYmass2S->Write();
        hGammaP2S->Write();
        hUpsP2S->Write();
        hCosAlphaP2S->Write();


        float Q1S,Q2S;
        float M1S,M2S;

        tQ->Branch("Q1S",&Q1S,"Q1S/F");
        tQ->Branch("invm1S",&M1S,"invm1S/F");
        tQ->Branch("Q2S",&Q2S,"Q2S/F");
        tQ->Branch("invm2S",&M2S,"invm2S/F");

        TRandom3 *randOM = new TRandom3();

        for(int i=0;i<nToy;i++){

          float Ymass_1S = hYmass1S->GetRandom();
          float p_gamma1S = hGammaP1S->GetRandom();
          float p_Ups1S = hUpsP1S->GetRandom();

          float cosAlpha1S = 2*randOM->Uniform()-1;
          cosAlpha1S = hCosAlphaP1S->GetRandom();
          float e_gamma1S = sqrt(m_gamma*m_gamma+p_gamma1S*p_gamma1S);
          float e_Ups1S = sqrt(Ymass_1S*Ymass_1S+p_Ups1S*p_Ups1S);

          Q1S = sqrt(m_gamma*m_gamma + Ymass_1S*Ymass_1S + 2*e_gamma1S*e_Ups1S - 2*p_gamma1S*p_Ups1S*cosAlpha1S) - m_gamma - Ymass_1S;
          M1S = sqrt(m_gamma*m_gamma + Ymass_1S*Ymass_1S + 2*e_gamma1S*e_Ups1S - 2*p_gamma1S*p_Ups1S*cosAlpha1S)-Ymass_1S+Ymass1S;

       	  hUpsP_check1S->Fill(p_Ups1S);
          hGammaP_check1S->Fill(p_gamma1S);
          hYmass_check1S->Fill(Ymass_1S);
          hCosAlphaP_check1S->Fill(cosAlpha1S);

          float Ymass_2S = hYmass2S->GetRandom();
          float p_gamma2S = hGammaP2S->GetRandom();
          float p_Ups2S = hUpsP2S->GetRandom();

          float cosAlpha2S = 2*randOM->Uniform()-1;
          cosAlpha2S = hCosAlphaP2S->GetRandom();
          float e_gamma2S = sqrt(m_gamma*m_gamma+p_gamma2S*p_gamma2S);
          float e_Ups2S = sqrt(Ymass_2S*Ymass_2S+p_Ups2S*p_Ups2S);

          Q2S = sqrt(m_gamma*m_gamma + Ymass_2S*Ymass_2S + 2*e_gamma2S*e_Ups2S - 2*p_gamma2S*p_Ups2S*cosAlpha2S) - m_gamma - Ymass_2S;
          M2S = sqrt(m_gamma*m_gamma + Ymass_2S*Ymass_2S + 2*e_gamma2S*e_Ups2S - 2*p_gamma2S*p_Ups2S*cosAlpha2S)-Ymass_2S+Ymass2S;

       	  hUpsP_check2S->Fill(p_Ups2S);
          hGammaP_check2S->Fill(p_gamma2S);
          hYmass_check2S->Fill(Ymass_2S);
          hCosAlphaP_check2S->Fill(cosAlpha2S);

          if (Q1S<3 | Q2S<3) {
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


        TH1F  *BkgToyHist1S = new TH1F("BkgToyHist1S","",100,9.5,12);
        tQ->Draw("invm1S>>BkgToyHist1S");
        TH1F  *BkgToyHist2S = new TH1F("BkgToyHist2S","",100,9.5,12);
        tQ->Draw("invm2S>>BkgToyHist2S");


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
    RooRealVar m_chib3P= RooRealVar("Mean #chi_{b}(3P)","Mean #chi_{b}(3P)",10.45,10.3,10.56);

    RooRealVar m_chib3P_1S= RooRealVar("Mean #chi_{b}(3P) in 1S sample","Mean #chi_{b}(3P) in 1S sample",10.45,10.3,10.56);


    RooRealVar sigma  = RooRealVar("#sigma","#sigma",0.01,0.005,0.02);
    RooRealVar alpha  = RooRealVar("#alpha","#alpha",0.654,0.,5.);
    RooRealVar n      = RooRealVar("n","n",2.58,.5,15.);
    RooRealVar sigma1  = RooRealVar("#sigma1","#sigma1",0.01,0.005,0.02);
    RooRealVar sigma2  = RooRealVar("#sigma2","#sigma2",0.01,0.005,0.02);
    RooRealVar sigma3  = RooRealVar("#sigma3","#sigma3",1.2094e-02);




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

    RooRealVar sigmaInd  = RooRealVar("sigmaInd","sigmaInd",0.010,0.003,0.015);
    RooRealVar alpha_3J  = RooRealVar("#alpha_3J","#alpha_3J",0.654);//,0.,2.);
    RooRealVar n_3J      = RooRealVar("n_3J","n_3J",4.,.5,15.);
    RooRealVar PhotonMassScale= RooRealVar("PhotonMassScale","PhotonMassScale",1.,0.95,1.05);
    RooRealVar PhotonMassScale2P= RooRealVar("PhotonMassScale2P","PhotonMassScale2P",1.,0.95,1.05);
    RooRealVar fracJ1= RooRealVar("fracJ1","fracJ1",0.666666);


    RooRealVar m_chib1P1fix= RooRealVar("m_chib1P1fix","m_chib1P1fix",9.89278);
    RooRealVar m_chib1P2fix= RooRealVar("m_chib1P2fix","m_chib1P2fix",9.91221);

    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass1S,Ymass1S);

    RooFormulaVar m_chib1P1float("m_chib1P1float", MassScaleFormula, RooArgList(m_chib1P1fix, PhotonMassScale));
	RooFormulaVar m_chib1P2float("m_chib1P2float", MassScaleFormula, RooArgList(m_chib1P2fix, PhotonMassScale));

	RooCBShape chib1P1      = RooCBShape ("chib1P1","chib1P1",invm1S,m_chib1P1float,sigmaInd,alpha_3J,n_3J);
    RooCBShape chib1P2      = RooCBShape ("chib1P2","chib1P2",invm1S,m_chib1P2float,sigmaInd,alpha_3J,n_3J);

    RooRealVar n1PnJ= RooRealVar("n1PJ","n1PJ",200,0,350);
    RooAddPdf chib1P_sigInd= RooAddPdf("chib1P_sigInd","chib1P_sigInd",RooArgList(chib1P1,chib1P2),RooArgList(fracJ1));

    RooFormulaVar chib1P_nevt("chib1P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n1PnJ,fracJ1));
    RooAddPdf chib1P_sig= RooAddPdf("chib1P_sig","chib1P_sig",RooArgList(chib1P_sigInd),RooArgList(chib1P_nevt));

    RooRealVar m_chib2P1fix= RooRealVar("m_chib2P1fix","m_chib2P1fix",10.25546);
    RooRealVar m_chib2P2fix= RooRealVar("m_chib2P2fix","m_chib2P2fix",10.26865);

	RooFormulaVar m_chib2P1float("m_chib2P1float", MassScaleFormula, RooArgList(m_chib2P1fix, PhotonMassScale2P));
	RooFormulaVar m_chib2P2float("m_chib2P2float", MassScaleFormula, RooArgList(m_chib2P2fix, PhotonMassScale2P));

    RooCBShape chib2P1      = RooCBShape ("chib2P1","chib2P1",invm1S,m_chib2P1float,sigmaInd,alpha_3J,n_3J);
    RooCBShape chib2P2      = RooCBShape ("chib2P2","chib2P2",invm1S,m_chib2P2float,sigmaInd,alpha_3J,n_3J);

    RooRealVar n2PnJ= RooRealVar("n2PJ","n2PJ",160,0,350);
    RooAddPdf chib2P_sigInd= RooAddPdf("chib2P_sigInd","chib2P_sigInd",RooArgList(chib2P1,chib2P2),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt("chib2P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n2PnJ,fracJ1));
    RooAddPdf chib2P_sig= RooAddPdf("chib2P_sig","chib2P_sig",RooArgList(chib2P_sigInd),RooArgList(chib2P_nevt));


    sprintf(MassScaleFormula,"(@0-%f)*@1+%f",Ymass2S,Ymass2S);

	RooFormulaVar m_chib2P1float_2S("m_chib2P1float_2S", MassScaleFormula, RooArgList(m_chib2P1fix, PhotonMassScale));
	RooFormulaVar m_chib2P2float_2S("m_chib2P2float_2S", MassScaleFormula, RooArgList(m_chib2P2fix, PhotonMassScale));

	RooCBShape chib2P1_2S      = RooCBShape ("chib2P1_2S","chib2P1_2S",invm2S,m_chib2P1float_2S,sigmaInd,alpha_3J,n_3J);
    RooCBShape chib2P2_2S      = RooCBShape ("chib2P2_2S","chib2P2_2S",invm2S,m_chib2P2float_2S,sigmaInd,alpha_3J,n_3J);

    RooRealVar n2PnJ_2S= RooRealVar("n2PnJ_2S","n2PnJ_2S",8,0,200);
    RooAddPdf chib2P_sigInd_2S= RooAddPdf("chib2P_sigInd_2S","chib2P_sigInd_2S",RooArgList(chib2P1_2S,chib2P2_2S),RooArgList(fracJ1));

    RooFormulaVar chib2P_nevt_2S("chib2P_nevt_2S", "@1*@0+(1-@1)*@0", RooArgList(n2PnJ_2S,fracJ1));
    RooAddPdf chib2P_sig_2S= RooAddPdf("chib2P_sig_2S","chib2P_sig_2S",RooArgList(chib2P_sigInd_2S),RooArgList(chib2P_nevt_2S));


    sprintf(MassScaleFormula,"(@0-%f-0.006)*@1+%f",Ymass2S,Ymass2S);
	RooFormulaVar m_chib3P1float("m_chib3P1float", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale));
    sprintf(MassScaleFormula,"(@0-%f+0.006)*@1+%f",Ymass2S,Ymass2S);
	RooFormulaVar m_chib3P2float("m_chib3P2float", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale));

    RooCBShape chib3P1      = RooCBShape ("chib3P1","chib3P1",invm2S,m_chib3P1float,sigmaInd,alpha_3J,n_3J);
    RooCBShape chib3P2      = RooCBShape ("chib3P2","chib3P2",invm2S,m_chib3P2float,sigmaInd,alpha_3J,n_3J);

    RooRealVar n3PnJ= RooRealVar("n3PnJ","n3PnJ",8,0,40);
    RooAddPdf chib3P_sigInd= RooAddPdf("chib3P_sigInd","chib3P_sigInd",RooArgList(chib3P1,chib3P2),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt("chib3P_nevt", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ,fracJ1));
    RooAddPdf chib3P_sig= RooAddPdf("chib3P_sig","chib3P_sig",RooArgList(chib3P_sigInd),RooArgList(chib3P_nevt));





    sprintf(MassScaleFormula,"(@0-%f-0.006)*@1+%f",Ymass1S,Ymass1S);
	RooFormulaVar m_chib3P1float_1S("m_chib3P1float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale2P));
    sprintf(MassScaleFormula,"(@0-%f+0.006)*@1+%f",Ymass1S,Ymass1S);
	RooFormulaVar m_chib3P2float_1S("m_chib3P2float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale2P));

    RooCBShape chib3P1_1S      = RooCBShape ("chib3P1_1S","chib3P1_1S",invm1S,m_chib3P1float_1S,sigmaInd,alpha_3J,n_3J);
    RooCBShape chib3P2_1S     = RooCBShape ("chib3P2_1S","chib3P2_1S",invm1S,m_chib3P2float_1S,sigmaInd,alpha_3J,n_3J);

    RooRealVar n3PnJ_1S= RooRealVar("n3PnJ_1S","n3PnJ_1S",8,0,40);
    RooAddPdf chib3P_sigInd_1S= RooAddPdf("chib3P_sigInd_1S","chib3P_sigInd_1S",RooArgList(chib3P1_1S,chib3P2_1S),RooArgList(fracJ1));

    RooFormulaVar chib3P_nevt_1S("chib3P_nevt_1S", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ_1S,fracJ1));
    RooAddPdf chib3P_sig_1S= RooAddPdf("chib3P_sig_1S","chib3P_sig_1S",RooArgList(chib3P_sigInd_1S),RooArgList(chib3P_nevt_1S));

////////////////////////////////////////////////////////


    RooRealVar background_nevt1S = RooRealVar("N_{bkg}1S","N_{bkg}1S",1500,0,10000);
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



    RooAddPdf modelPdf1S= RooAddPdf("modelPdf1S","modelPdf1S",RooArgList(chib1P_sig,chib2P_sig,chib3P_sig_1S,background1S),RooArgList(chib1P_nevt,chib2P_nevt,chib3P_nevt_1S,background_nevt1S));

    RooAddPdf modelPdf2S= RooAddPdf("modelPdf2S","modelPdf2S",RooArgList(chib3P_sig,chib2P_sig_2S,background2S),RooArgList(chib3P_nevt,chib2P_nevt_2S,background_nevt2S));




    double bb_thresh=10.56;

    invm1S.setRange("SBregion1",invm_min1S,9.85);
    invm1S.setRange("SBregion2",9.95,10.15);
    invm1S.setRange("SBregion31S",bb_thresh,invm_max1S);

    invm1S.setRange("Chib1Pregion",invm_min1S,9.95);
    invm1S.setRange("Chib2Pregion",10.15,10.3);
    invm1S.setRange("Chib3Pregion",10.3,bb_thresh);

    invm1S.setRange("TotalBlindedRegion1",invm_min1S,10.3);
    invm1S.setRange("TotalBlindedRegion2",bb_thresh,invm_max1S);


    invm2S.setRange("Y2SAll",invm_min2S,invm_max2S);
    invm2S.setRange("SBregion32S",bb_thresh,invm_max2S);




    double covQual_BG1S=0;
    double covQual_1P2P1S=0;
    char BGnormRegion[200];
    sprintf(BGnormRegion,"SBregion2");

    int nFits=10;
    int nBGFits=0;
    for(int n=0;n<nFits;n++){
    RooFitResult* FitBG1S = backgroundSB1S.fitTo(*data1S,Range(BGnormRegion),Save(1));
    FitBG1S->Print();
    covQual_BG1S = FitBG1S->covQual();
    if(covQual_BG1S==3) break;
    nBGFits++;
    delete FitBG1S;
    }

    cout<<"background_nevt1S = " << background_nevtSB1S.getVal()<<endl;

    double N_backgroundSB21S;
    double N_backgroundSB31S;
    double BkgSBInt1S;
    BkgSBInt1S = NormalizedIntegral(&backgroundSB1S, invm1S, 9.95,10.15);
    cout<<BkgSBInt1S<<endl;
    background_nevt1S.setVal(background_nevtSB1S.getVal()/BkgSBInt1S);
    background_nevt1S.setConstant();
    N_backgroundSB21S=background_nevt1S.getVal();

    cout<<"background_nevt1S all = " << background_nevt1S.getVal()<<endl;
    cout<<"background_nevt1S SB = " << background_nevtSB1S.getVal()<<endl;


    sprintf(BGnormRegion,"SBregion31S");

    for(int n=0;n<nFits;n++){
    RooFitResult* FitBG = backgroundSB1S.fitTo(*data1S,Range(BGnormRegion),Save(1));
    FitBG->Print();
    covQual_BG1S = FitBG->covQual();
    if(covQual_BG1S==3) break;
    nBGFits++;
    delete FitBG;
    }
    BkgSBInt1S = NormalizedIntegral(&backgroundSB1S, invm1S, bb_thresh, invm_max1S);
    cout<<BkgSBInt1S<<endl;
    background_nevt1S.setVal(background_nevtSB1S.getVal()/BkgSBInt1S);
    background_nevt1S.setConstant();
    N_backgroundSB31S=background_nevt1S.getVal();

    cout<<"background_nevt1S all = " << background_nevt1S.getVal()<<endl;
    cout<<"background_nevt1S SB = " << background_nevtSB1S.getVal()<<endl;

    background_nevt1S.setVal((N_backgroundSB21S+N_backgroundSB31S)/2.);
    background_nevt1S.setConstant();

    cout<<"background_nevt1S all = " << background_nevt1S.getVal()<<endl;








    double covQual_BG2S=0;
    double covQual_1P2P2S=0;
    sprintf(BGnormRegion,"SBregion32S");

    if(nState==2){
    for(int n=0;n<nFits;n++){
    RooFitResult* FitBG2S = backgroundSB2S.fitTo(*data2S,Range(BGnormRegion),Save(1));
    FitBG2S->Print();
    covQual_BG2S = FitBG2S->covQual();
    if(covQual_BG2S==3) break;
    nBGFits++;
    delete FitBG2S;
    }
    }
    cout<<"background_nevt2S = " << background_nevtSB2S.getVal()<<endl;

    double N_backgroundSB22S;
    double N_backgroundSB32S;
    double BkgSBInt2S;
    BkgSBInt2S = NormalizedIntegral(&backgroundSB2S, invm2S, bb_thresh, invm_max2S);
    cout<<BkgSBInt2S<<endl;
    background_nevt2S.setVal(background_nevtSB2S.getVal()/BkgSBInt2S);
    background_nevt2S.setConstant();
    N_backgroundSB22S=background_nevt2S.getVal();

    cout<<"background_nevt2S all = " << background_nevt2S.getVal()<<endl;
    cout<<"background_nevt2S SB = " << background_nevtSB2S.getVal()<<endl;


    double Nevt_buff=n3PnJ_1S.getVal();
    n3PnJ_1S.setVal(0);
    n3PnJ_1S.setConstant();
    m_chib3P.setConstant();

    char SignormRegion[200];
    sprintf(SignormRegion,"Chib1Pregion,Chib2Pregion");

        RooFitResult* Fit1P2P = modelPdf1S.fitTo(*data1S,Save(1));
        Fit1P2P->Print();
        covQual_1P2P1S = Fit1P2P->covQual();

        n3PnJ_1S.setConstant(kFALSE);
        n3PnJ_1S.setVal(Nevt_buff);
        m_chib3P.setConstant(kFALSE);

        sigma.setConstant();
        alpha.setConstant();
        n.setConstant();

        sigmaInd.setConstant();
        alpha_3J.setConstant();
        n_3J.setConstant();
        PhotonMassScale.setConstant();
        PhotonMassScale2P.setConstant();
        fracJ1.setConstant();
        n1PnJ.setConstant();
        n2PnJ.setConstant();
        delete Fit1P2P;



if(nState==2){

		Nevt_buff=n3PnJ.getVal();
		n3PnJ.setVal(0.00001);
		n3PnJ.setConstant();
        m_chib3P.setConstant();

	    RooFitResult* Fit1P2P = modelPdf2S.fitTo(*data2S,Save(1));
        Fit1P2P->Print();
        covQual_1P2P2S = Fit1P2P->covQual();

        n3PnJ.setConstant(kFALSE);
        n3PnJ.setVal(Nevt_buff);
        m_chib3P.setConstant(kFALSE);
		n2PnJ_2S.setConstant();
        delete Fit1P2P;

}

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



		 RooAbsReal* nll1S=modelPdf1S.createNLL(*data1S);
		 RooAbsReal* nll2S=modelPdf2S.createNLL(*data2S);

         RooAddition simNLL = RooAddition("add","add",RooArgSet(*nll1S,*nll2S));
         double Nllbefore3P=simNLL.getVal();

		 RooMinuit* minuit = new RooMinuit(simNLL);

		 minuit->setStrategy(2);
		 minuit->setPrintEvalErrors(-1);
		 minuit->setPrintLevel(3);

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

         double Nllafter3P=simNLL.getVal();

         cout<<"minNllNoSigAll "<<Nllbefore3P<<endl;
         cout<<"minNll3PAll    "<<Nllafter3P<<endl;
         double logLikelihoodRatio3PAll=TMath::Log(TMath::Exp(Nllafter3P-Nllbefore3P));
         cout<<"logLikelihoodRatio3PAll "<<logLikelihoodRatio3PAll<<endl;
         double ThreePsigALL=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3PAll));
         cout<<"ThreePsigALL "<<ThreePsigALL<<endl;






        Nevt_buff=n3PnJ.getVal();
		double minNll3P=modelPdf2S.createNLL(*data2S)->getVal();
		n3PnJ.setVal(0);
		double minNllNoSig2S=modelPdf2S.createNLL(*data2S)->getVal();
		n3PnJ.setVal(Nevt_buff);

        cout<<"minNllNoSig2S "<<minNllNoSig2S<<endl;
        cout<<"minNll3P    "<<minNll3P<<endl;
        double logLikelihoodRatio3P=TMath::Log(TMath::Exp(minNll3P-minNllNoSig2S));
        cout<<"logLikelihoodRatio3P "<<logLikelihoodRatio3P<<endl;
        double ThreePsig=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3P));
        cout<<"ThreePsig "<<ThreePsig<<endl;

        Nevt_buff=n3PnJ_1S.getVal();
		double minNll3P_1S=modelPdf1S.createNLL(*data1S)->getVal();
		n3PnJ_1S.setVal(0);
		double minNllNoSig1S=modelPdf1S.createNLL(*data1S)->getVal();
		n3PnJ_1S.setVal(Nevt_buff);

		cout<<"minNllNoSig1S "<<minNllNoSig1S<<endl;
        cout<<"minNll3P_1S    "<<minNll3P_1S<<endl;
        double logLikelihoodRatio3P_1S=TMath::Log(TMath::Exp(minNll3P_1S-minNllNoSig1S));
        cout<<"logLikelihoodRatio3P_1S "<<logLikelihoodRatio3P_1S<<endl;
        double ThreePsig_1S=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio3P_1S));
        cout<<"ThreePsig_1S "<<ThreePsig_1S<<endl;

        Nevt_buff=n1PnJ.getVal();
		double minNll1P=modelPdf1S.createNLL(*data1S)->getVal();
		n1PnJ.setVal(0);
		double minNllNo1P1S=modelPdf1S.createNLL(*data1S)->getVal();
		n1PnJ.setVal(Nevt_buff);

		cout<<"minNllNo1P1S "<<minNllNo1P1S<<endl;
        cout<<"minNll1P    "<<minNll1P<<endl;
        double logLikelihoodRatio1P=TMath::Log(TMath::Exp(minNll1P-minNllNo1P1S));
        cout<<"logLikelihoodRatio1P "<<logLikelihoodRatio1P<<endl;
        double OnePsig=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio1P));
        cout<<"OnePsig "<<OnePsig<<endl;

        Nevt_buff=n2PnJ.getVal();
		double minNll2P=modelPdf1S.createNLL(*data1S)->getVal();
		n2PnJ.setVal(0);
		double minNllNo2P1S=modelPdf1S.createNLL(*data1S)->getVal();
		n2PnJ.setVal(Nevt_buff);

		cout<<"minNllNo2P1S "<<minNllNo2P1S<<endl;
        cout<<"minNll2P    "<<minNll2P<<endl;
        double logLikelihoodRatio2P=TMath::Log(TMath::Exp(minNll2P-minNllNo2P1S));
        cout<<"logLikelihoodRatio2P "<<logLikelihoodRatio2P<<endl;
        double TwoPsig=TMath::Sqrt(2*TMath::Abs(logLikelihoodRatio2P));
        cout<<"TwoPsig "<<TwoPsig<<endl;



        cout<<"chib1P_nevt: "<<chib1P_nevt.getVal()<<endl;
        cout<<"chib2P_nevt: "<<chib2P_nevt.getVal()<<endl;
        cout<<"chib3P_nevt_1S: "<<chib3P_nevt_1S.getVal()<<endl;
        cout<<"chib3P_nevt: "<<chib3P_nevt.getVal()<<endl;
        cout<<"chib2P_nevt_2S: "<<chib2P_nevt_2S.getVal()<<endl;


     double TotalRealEvents1S = data1S_masswindow->sumEntries();
     double totalEventsInFit1S;
     double TotalRealEvents2S = data2S_masswindow->sumEntries();
     double totalEventsInFit2S;

     totalEventsInFit1S = chib1P_nevt.getVal()+chib2P_nevt.getVal()+chib3P_nevt_1S.getVal()+background_nevt1S.getVal();
     totalEventsInFit2S = chib3P_nevt.getVal()+chib2P_nevt_2S.getVal()+background_nevt2S.getVal();

     cout<< "N_tot_Fit1S =                            "<<totalEventsInFit1S <<endl;
     cout<< "N_tot_Samlpe1S =                         "<<TotalRealEvents1S <<endl;
     cout<< "Delta_N_tot1S = N_tot_Sample1S-N_tot_Fit1S = "<<TotalRealEvents1S-totalEventsInFit1S <<endl;

     cout<< "N_tot_Fit2S =                            "<<totalEventsInFit2S <<endl;
     cout<< "N_tot_Samlpe2S =                         "<<TotalRealEvents2S <<endl;
     cout<< "Delta_N_tot2S = N_tot_Sample2S-N_tot_Fit2S = "<<TotalRealEvents2S-totalEventsInFit2S <<endl;


     double nSig=1.5;



     double fmchib1P = m_chib1P.getVal();
     double fmchib2P = m_chib2P.getVal();
     double fmchib3P = (m_chib3P.getVal()-0.006)*fracJ1.getVal()+(m_chib3P.getVal()+0.006)*(1-fracJ1.getVal());
     double fmchib3P_1S = m_chib3P.getVal();

     m_chib3P.setVal(fmchib3P);

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
     cout<<"m_chib3P: "<<m_chib3P.getVal()<<endl;
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
     PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_2P->GetXaxis()->SetTitle("E_{#gamma} in 2P peak (1S sample)");	hGammaE_2P->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_2P_1S.pdf",FitID);
     PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_3P->GetXaxis()->SetTitle("E_{#gamma} in 3P peak (2S sample)");	hGammaE_3P->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy_3P_2S.pdf",FitID);
     PhotonEnergyCanvas->SaveAs(saveName);
     gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);hGammaE_3P_1S->GetXaxis()->SetTitle("E_{#gamma} in 3P peak (1S sample)");	hGammaE_3P_1S->Draw();
     PhotonEnergyCanvas->Modified();
     sprintf(saveName,"Figures/%s/PhotonEnergy3P_1S.pdf",FitID);
     PhotonEnergyCanvas->SaveAs(saveName);

     hCut_1Psig->SetBinContent(inCut+1,rsigchib1P);
     if(rratchib1P<100) hCut_1PSB->SetBinContent(inCut+1,rratchib1P);
     hCut_1PS->SetBinContent(inCut+1,rchib1P);
     hCut_1PB->SetBinContent(inCut+1,rbgchib1P);
     hCut_2Psig->SetBinContent(inCut+1,rsigchib2P);
     if(rratchib2P<100) hCut_2PSB->SetBinContent(inCut+1,rratchib2P);
     hCut_2PS->SetBinContent(inCut+1,rchib2P);
     hCut_2PB->SetBinContent(inCut+1,rbgchib2P);
     hCut_sig->SetBinContent(inCut+1,(rsigchib1P+rsigchib2P)/2.);
     hCut_covQual1->SetBinContent(inCut+1,covQual_BG1S);
     hCut_1Pmean ->SetBinContent(inCut+1,fmchib1P);
     hCut_1Pwidth->SetBinContent(inCut+1,sigmaInd.getVal());
     hCut_1Pnevt ->SetBinContent(inCut+1,fnchib1P);
     hCut_2Pmean ->SetBinContent(inCut+1,fmchib2P);
     hCut_2Pwidth->SetBinContent(inCut+1,sigmaInd.getVal());
     hCut_2Pnevt ->SetBinContent(inCut+1,fnchib2P);
     hCut_CBa    ->SetBinContent(inCut+1,alpha_3J.getVal());
     hCut_CBn    ->SetBinContent(inCut+1,n_3J.getVal());
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

     hCut_3Psig_1S->SetBinContent(inCut+1,ThreePsig_1S);
     if(rratchib3P_1S<100) hCut_3PSB_1S->SetBinContent(inCut+1,rratchib3P_1S);
     hCut_3PS_1S->SetBinContent(inCut+1,rchib3P_1S);
     hCut_3PB_1S->SetBinContent(inCut+1,rbgchib3P_1S);

     hCut_3PsigALL->SetBinContent(inCut+1,ThreePsigALL);

     hCut_PhotonEnergyScale->SetBinContent(inCut+1,PhotonMassScale.getVal());

     hCut_1PsigLL->SetBinContent(inCut+1,OnePsig);
     hCut_2PsigLL->SetBinContent(inCut+1,TwoPsig);


     double linewidth = 1;
     double chi21S;
     double chi22S;
     int FrameBins1S=10./4.*64.;
     int FrameBins2Sin1S=10./4.*44.;
     int FrameBins2S=10./4.*64.;
     int FrameBins3Sin1S=10./4.*32.;




    TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",1600, 800);
    ChibCanvas1S->Divide(1);
    ChibCanvas1S->SetFillColor(kWhite);
    ChibCanvas1S->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(kWhite);

    RooPlot* frame1S= invm1S.frame(FrameBins1S);
    frame1S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    data1S->plotOn(frame1S);

        modelPdf1S.plotOn(frame1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));;
    chi21S = frame1S->chiSquare();
    /////////// Two CB adventure //////////////////
    modelPdf1S.plotOn(frame1S, Components("chib1P1"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib1P2"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib2P1"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib2P2"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib3P1_1S"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("chib3P2_1S"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S, Components("background1S"), LineStyle(2),LineColor(kBlack),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));

    frame1S->SetTitle(0);
    frame1S->SetMinimum(0);
    frame1S->Draw();

    double xText=10.5;
    char text[200];
    sprintf(text,"S/B #chi_{b}(1P) = %1.0f / %1.0f, #chi_{b}(1P)_{sig.} = %1.4f",rchib1P,rbgchib1P,rsigchib1P);
    TLatex text1 = TLatex(xText,frame1S->GetMaximum()*0.9,text);
    text1.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text1.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(2P) = %1.0f / %1.0f, #chi_{b}(2P)_{sig.} = %1.4f",rchib2P,rbgchib2P,rsigchib2P);
    TLatex text2 = TLatex(xText,frame1S->GetMaximum()*0.85,text);
    text2.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text2.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.4f",rchib3P_1S,rbgchib3P_1S,ThreePsig_1S);
    TLatex text3 = TLatex(xText,frame1S->GetMaximum()*0.8,text);
    text3.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text3.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",m_chib3P.getVal(),m_chib3P.getError());
    TLatex text0 = TLatex(xText,frame1S->GetMaximum()*0.75,text)                                                                                                                            ;
    text0.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text0.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"N_{tot} = %1.0f",data1S_masswindow->sumEntries());
    TLatex text4 = TLatex(xText,frame1S->GetMaximum()*0.65,text);                                                                                                                                                                        ;
    text4.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text4.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi21S);
    TLatex text6 = TLatex(xText,frame1S->GetMaximum()*0.6,text)                                                                                                                                                                                          ;
    text6.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text6.Draw( "same" )                                                                                                                                                                                                                                                 ;

    ChibCanvas1S->Modified();
    sprintf(cutName_,"vtxProb%d_ctCut%d_nYsig%d_cut_gammapt%d_cut_RconvMin%d_cut_Ypt%d",int(1000000*cut_vtxProb),int(1000000*cut_ct),int(1000000*nYsig),int(1000000*cut_gammapt),int(1000000*cut_RconvMin),int(1000000*cut_Ypt));
    sprintf(saveName,"Figures/%s/InvMass_Y1S_%s.pdf",FitID,cutName_);
    ChibCanvas1S->SaveAs(saveName);






    TCanvas* ChibCanvas2S = new TCanvas("#chi_{b} 2S invariant mass","#chi_{b} 2S invariant mass",1600, 800);
    ChibCanvas2S->Divide(1);
    ChibCanvas2S->SetFillColor(kWhite);
    ChibCanvas2S->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(kWhite);

    RooPlot* frame2S= invm2S.frame(invm_min1S,invm_max2S,FrameBins2S);
    frame2S->SetTitle("#chi_{b} invariant mass, Y2S decay");
    data2S->plotOn(frame2S,MarkerStyle(24));


    modelPdf2S.plotOn(frame2S,Range(invm_min2S,invm_max2S),LineWidth(linewidth),Normalization(totalEventsInFit2S,2));;
    chi22S = frame2S->chiSquare();

    modelPdf2S.plotOn(frame2S, Components("chib3P1"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("chib3P2"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("chib2P1_2S"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("chib2P2_2S"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S, Components("background2S"), LineStyle(2),LineColor(kBlack),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame2S,Range(invm_min2S,invm_max2S),LineColor(kGreen+2),LineWidth(linewidth),Normalization(totalEventsInFit2S,2));


    frame2S->SetTitle(0);
    frame2S->SetMinimum(0);
    frame2S->Draw();


    sprintf(text,"N_{tot} = %1.0f",data2S_masswindow->sumEntries());
    TLatex text9 = TLatex(xText,frame2S->GetMaximum()*0.75,text);                                                                                                                                                                        ;
    text9.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text9.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"#chi^{2}/ndf = %1.4f",chi22S);
    TLatex text10 = TLatex(xText,frame2S->GetMaximum()*0.7,text)                                                                                                                                                                                          ;
    text10.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text10.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",m_chib3P.getVal(),m_chib3P.getError());
    TLatex text5 = TLatex(xText,frame2S->GetMaximum()*0.9,text)                                                                                                                            ;
    text5.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text5.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.4f",rchib3P,rbgchib3P,ThreePsig);
    TLatex text8 = TLatex(xText,frame2S->GetMaximum()*0.85,text)                                                                                                                            ;
    text8.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text8.Draw( "same" )                                                                                                                                                                                                                                                 ;

    ChibCanvas2S->Modified();
    sprintf(cutName_,"vtxProb%d_ctCut%d_nYsig%d_cut_gammapt%d_cut_RconvMin%d_cut_Ypt%d",int(1000000*cut_vtxProb),int(1000000*cut_ct),int(1000000*nYsig),int(1000000*cut_gammapt),int(1000000*cut_RconvMin),int(1000000*cut_Ypt));
    sprintf(saveName,"Figures/%s/InvMass_Y2S_%s.pdf",FitID,cutName_);
    if(nState==2) ChibCanvas2S->SaveAs(saveName);







    TCanvas* ChibCanvas1S2S = new TCanvas("#chi_{b} 1S2S invariant mass","#chi_{b} 1S2S invariant mass",1600, 800);
    ChibCanvas1S2S->Divide(1);
    ChibCanvas1S2S->SetFillColor(kWhite);
    ChibCanvas1S2S->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(kWhite);

    RooPlot* frame1S2S= invm1S.frame(FrameBins1S);
    frame1S2S->SetTitle("#chi_{b} invariant mass, Y1S decay");
    data1S->plotOn(frame1S2S);

    modelPdf1S.plotOn(frame1S2S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));;
    modelPdf1S.plotOn(frame1S2S, Components("chib1P_sig"), LineStyle(2),LineColor(kMagenta),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S, Components("chib2P_sig"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S, Components("chib3P_sig_1S"), LineStyle(2),LineColor(kOrange),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S, Components("background1S"), LineStyle(2),LineColor(kBlack),LineWidth(linewidth),Range(invm_min1S,invm_max1S),Normalization(totalEventsInFit1S,2));
    modelPdf1S.plotOn(frame1S2S,Range(invm_min1S,invm_max1S),LineWidth(linewidth),Normalization(totalEventsInFit1S,2));

    RooPlot* frame1S2S_= invm2S.frame(FrameBins2Sin1S);
    frame1S2S_->SetTitle("#chi_{b} invariant mass, Y1S decay");
    data2S->plotOn(frame1S2S_,MarkerStyle(24));

    modelPdf2S.plotOn(frame1S2S_, Components("chib3P_sig"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame1S2S_, Components("chib2P_sig_2S"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame1S2S_, Components("background2S"), LineStyle(2),LineColor(kBlack),LineWidth(linewidth),Range(invm_min2S,invm_max2S),Normalization(totalEventsInFit2S,2));
    modelPdf2S.plotOn(frame1S2S_,Range(invm_min2S,invm_max2S),LineWidth(linewidth),LineColor(kGreen+2),Normalization(totalEventsInFit2S,2));

    RooPlot* frame1S2S__= invm3S.frame(FrameBins3Sin1S);
    frame1S2S__->SetTitle("#chi_{b} invariant mass, Y3S decay");
    data3S->plotOn(frame1S2S__,MarkerStyle(22));

    sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(nS)^{PDG}}");

    frame1S2S->GetXaxis()->SetTitle(invmName);
    frame1S2S->SetTitle(0);
    frame1S2S->SetMinimum(0);
    frame1S2S->Draw();
    frame1S2S_->SetTitle(0);
    frame1S2S_->SetMinimum(0);
    frame1S2S_->Draw("same");
    frame1S2S__->SetTitle(0);
    frame1S2S__->SetMinimum(0);
//    frame1S2S__->Draw("same");


    sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",m_chib3P.getVal(),m_chib3P.getError());
    TLatex text14 = TLatex(xText,frame1S->GetMaximum()*0.9,text)                                                                                                                            ;
    text14.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text14.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Y(2S) + #gamma: S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.4f",rchib3P,rbgchib3P,ThreePsig);
    TLatex text11 = TLatex(xText,frame1S->GetMaximum()*0.85,text)                                                                                                                            ;
    text11.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text11.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Y(1S) + #gamma: S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.4f",rchib3P_1S,rbgchib3P_1S,ThreePsig_1S);
    TLatex text12 = TLatex(xText,frame1S->GetMaximum()*0.8,text)                                                                                                                            ;
    text12.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text12.Draw( "same" )                                                                                                                                                                                                                                                 ;
    sprintf(text,"Combined Significance: #chi_{b}(3P)_{sig.} = %1.4f",ThreePsigALL);
    TLatex text13 = TLatex(xText,frame1S->GetMaximum()*0.75,text)                                                                                                                            ;
    text13.SetTextSize(0.025)                                                                                                                                                                                                                                             ;
    text13.Draw( "same" )                                                                                                                                                                                                                                                 ;


    ChibCanvas1S2S->Modified();
    sprintf(cutName_,"vtxProb%d_ctCut%d_nYsig%d_cut_gammapt%d_cut_RconvMin%d_cut_Ypt%d",int(1000000*cut_vtxProb),int(1000000*cut_ct),int(1000000*nYsig),int(1000000*cut_gammapt),int(1000000*cut_RconvMin),int(1000000*cut_Ypt));
    sprintf(saveName,"Figures/%s/InvMass_Y1S2S_%s.pdf",FitID,cutName_);
    if(nState==2) ChibCanvas1S2S->SaveAs(saveName);



    cout<<"deleting pointers"<<endl;

    delete data_______;
    delete data__;
    delete data___;
    delete data____;
    delete data_____;
    delete data______;
    delete data______vtx;

    delete dataSB1S;
    delete dataSB2S;
    delete data1S;
    delete data2S;
    delete data3S;

    delete minuit;
    delete PhotonEnergyCanvas;
    delete fQ;

	delete frame1S;
	delete frame2S;
	delete frame1S2S_;

	delete ChibCanvas1S;
    delete ChibCanvas2S;
    delete ChibCanvas1S2S;


    cout<<"deleted pointers"<<endl;
    }

    sprintf(filename_,"Figures/%s/CutHistos_Y%dS.root",FitID,nState);
	TFile *CutHistos_file = new TFile(filename_,"RECREATE");

	char XTitle[200];
	if(nYsigCut  )sprintf(XTitle,"Y(1S) mass window, n#sigma_{m_{#mu#mu}}(|y|^{Y(1S)}) ");
	if(vtxProbCut)sprintf(XTitle,"Vertex probability cut");
	if(ctSigCut  )sprintf(XTitle,"c#tau significance cut, n#sigma");
	if(gammaptCut)sprintf(XTitle,"#gamma - p_{T} cut, GeV");
	if(ctCut     )sprintf(XTitle,"c#tau cut, cm");
	if(RConvCut  )sprintf(XTitle,"RConv cut");
	if(pTCut  )sprintf(XTitle,"lower p_{T}^{Y} cut, GeV");

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

    hCut_3Pmass->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_3Pmass->SetStats(0);hCut_3Pmass->GetYaxis()->SetTitle("#chi_{b}(3P) mass");	hCut_3Pmass->GetXaxis()->SetTitle(XTitle);	hCut_3Pmass->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/3Pmass.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_PhotonEnergyScale->GetYaxis()->SetTitleOffset(yOffset); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_PhotonEnergyScale->SetStats(0);hCut_PhotonEnergyScale->GetYaxis()->SetTitle("Photon Energy Scale");	hCut_PhotonEnergyScale->GetXaxis()->SetTitle(XTitle);	hCut_PhotonEnergyScale->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/PhotonEnergyScale.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);


    hCut_2PsigLL->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_2PsigLL->SetStats(0);hCut_2PsigLL->GetYaxis()->SetTitle("#chi_{b}(2P)_{sign.LL}");	hCut_2PsigLL->GetXaxis()->SetTitle(XTitle);	hCut_2PsigLL->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/2PsigLL.pdf",FitID);
    IndividualCanvas->SaveAs(saveName);

    hCut_1PsigLL->GetYaxis()->SetTitleOffset(yOffset);gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite); hCut_1PsigLL->SetStats(0);hCut_1PsigLL->GetYaxis()->SetTitle("#chi_{b}(1P)_{sign.LL}");	hCut_1PsigLL->GetXaxis()->SetTitle(XTitle);	hCut_1PsigLL->Draw();
    IndividualCanvas->Modified();
    sprintf(saveName,"Figures/%s/1PsigLL.pdf",FitID);
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
	hCut_1PsigLL->Write();
	hCut_1PsigLL->Write();


	CutHistos_file->Write();
    CutHistos_file->Close();

	return 0;
}


