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
	  TFile* resultsFile = new TFile(filename, "RECREATE", "results");


	  double S_in_tree;
	  double S_out_tree;
	  double B_in_tree;
	  double S_tree;
	  double B_tree;
	  double SB_tree;
	  double Sig_tree;
	  double DeltaM_tree;
	  double Merr_tree;
	  double Width_tree;

	  TTree* Results = new TTree("Results","Results");
	  Results->Branch("S_in",     	   &S_in_tree,         "S_in/D");
	  Results->Branch("S_out",     	   &S_out_tree,         "S_out/D");
	  Results->Branch("B_in",    	   &B_in_tree,         "B_in/D");
	  Results->Branch("S",      	   &S_tree,         "S/D");
	  Results->Branch("B",      	   &B_tree,         "B/D");
	  Results->Branch("SB",      	   &SB_tree,         "SB/D");
	  Results->Branch("Sig",     	   &Sig_tree,         "Sig/D");
	  Results->Branch("DeltaM",        &DeltaM_tree,         "DeltaM/D");
	  Results->Branch("Merr",     	   &Merr_tree,         "Merr/D");
	  Results->Branch("Width",     	   &Width_tree,         "Width/D");
	  Results->Branch("Width",     	   &Width_tree,         "Width/D");

//load RooFit library
	using namespace RooFit;
  gSystem->Load("libRooFit");
  gROOT->SetBatch(1);
  gStyle->SetFrameBorderMode(0);

//declare RooDataSet data variables
  char invmName[200];

  double Ymass1S;
  double Ymass2S;
  Ymass1S=9.4603;
  Ymass2S=10.02326;

  double invm_min1S;
  double invm_max1S;
  double invm_min2S;
  double invm_max2S;

  invm_min1S=10.325; invm_max1S=10.58;
  invm_min2S=10; invm_max2S=11.1;

  double InvSpace=invm_max1S-invm_min1S;

  sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(1S)^{PDG}} [GeV]");
  RooRealVar invm1S("invm1S",invmName,invm_min1S,invm_max1S);
  sprintf(invmName,"m_{#mu#mu#gamma}-m_{#mu#mu}+m_{#Upsilon(2S)^{PDG}} [GeV]");
  RooRealVar invm2S("invm2S",invmName,invm_min2S,invm_max2S);

  invm1S.setRange("Y1SAll",invm_min1S,invm_max1S);


  char filename_[500];
  sprintf(filename_,"EventMixer/eventMixer1S.root");
  TFile EventMixerfile1S(filename_,"READ");
  sprintf(filename_,"EventMixer/eventMixer2S.root");
  TFile EventMixerfile2S(filename_,"READ");

  TH1F* hInvM1S;
  TH1F* hInvM2S;
  hInvM1S=(TH1F*)EventMixerfile1S.Get("hInvMnS");
  hInvM2S=(TH1F*)EventMixerfile2S.Get("hInvMnS");



  for(int iSig=1;iSig<niSig+1;iSig++){

	  int nSig=nSig0+(iSig-1)*delta_nSig;

	  for(int iBkg=1;iBkg<niBkg+1;iBkg++){

		  int nBkg=nBkg0+(iBkg-1)*delta_nBkg;

		  double PhotonSigmaScale3P_STARTStart=0.015;
		  double sigRatioEstimate=0.75;
		  double bkgRatioEstimate=0.35*PhotonSigmaScale3P_STARTStart/0.015;
		  double SoverBCheck=(nSig*sigRatioEstimate)/(nBkg*bkgRatioEstimate);
		  double SigCheck=(nSig*sigRatioEstimate)/TMath::Sqrt(nSig*sigRatioEstimate+nBkg*bkgRatioEstimate);

		  double minSig=1.5;
		  double maxSig=13.5;

		  if(SigCheck>maxSig||SigCheck<minSig) continue;

  double S_in_treeBuffer=0;
  double S_out_treeBuffer=0;
  double B_in_treeBuffer=0;
  double S_treeBuffer=0;
  double B_treeBuffer=0;
  double SB_treeBuffer=0;
  double Sig_treeBuffer=0;
  double DeltaM_treeBuffer=0;
  double Merr_treeBuffer=0;
  double Width_treeBuffer=0;

  for(int iToy=1;iToy<nToy+1;iToy++){

  RooDataHist* RooBkgToyHist1S_Mixer = new RooDataHist("RooBkgToyHist1S_Mixer","RooBkgToyHist1S_Mixer",RooArgList(invm1S),hInvM2S);
  RooAbsPdf* background1S;
  background1S = new  RooHistPdf("background1S","background1S",RooArgSet(invm1S),*RooBkgToyHist1S_Mixer,1);


  double alpha_3J_START=0.5;
  double n_3J_START=3.;
  double PhotonMassScale3P_START=0.985;
  double PhotonSigmaScale3P_START=PhotonSigmaScale3P_STARTStart;
  double ratio_J2overJ1_START=.5;
  double background_nevt1S_START=nBkg;
  double m_chib3P_START=10.5;
  double n3PnJ_1S_START=nSig;


  RooRealVar m_chib3P= RooRealVar("M_{#chi_{b}(3P)}","Mean #chi_{b}(3P)",m_chib3P_START,10.325,10.56);
  RooRealVar alpha_3J  = RooRealVar("#alpha","#alpha_3J",alpha_3J_START,0.,3.);//0.2,1.);
  RooRealVar n_3J      = RooRealVar("n","n_3J",n_3J_START,.5,15.);
  RooRealVar PhotonSigmaScale3P= RooRealVar("#sigma_{Q}/Q, 3P","PhotonSigmaScale",PhotonSigmaScale3P_START,0.0001,0.05);
  RooRealVar PhotonMassScale3P= RooRealVar("PES_{3P}","PhotonMassScale3P",PhotonMassScale3P_START,0.95,1.05);
  RooRealVar ratio_J2overJ1= RooRealVar("N_{#chi_{b2}}/N_{#chi_{b1}}","ratio_J2overJ1",ratio_J2overJ1_START,0.,5.);

  RooFormulaVar fracJ1("fracJ1", "1./(1.+@0)", RooArgList(ratio_J2overJ1));

  m_chib3P.setVal(m_chib3P_START);
  alpha_3J.setVal(alpha_3J_START);
  n_3J.setVal(n_3J_START);
  PhotonMassScale3P.setVal(PhotonMassScale3P_START);
  PhotonSigmaScale3P.setVal(PhotonSigmaScale3P_START);
  ratio_J2overJ1.setVal(ratio_J2overJ1_START);

  ratio_J2overJ1.setConstant();
  alpha_3J.setConstant();
  n_3J.setConstant();
  PhotonSigmaScale3P.setConstant();
  PhotonMassScale3P.setConstant();

  char MassScaleFormula[200];
  sprintf(MassScaleFormula,"(@0-%f-0.012*(1-@2))*@1+%f",Ymass1S,Ymass1S);
  RooFormulaVar m_chib3P1float_1S("m_chib3P1float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale3P, fracJ1));
  sprintf(MassScaleFormula,"(@0-%f+0.012*@2)*@1+%f",Ymass1S,Ymass1S);
  RooFormulaVar m_chib3P2float_1S("m_chib3P2float_1S", MassScaleFormula, RooArgList(m_chib3P, PhotonMassScale3P, fracJ1));

  char SigmaScaleFormula[200];
  sprintf(SigmaScaleFormula,"(@0-%f)*@1",Ymass1S);
  RooFormulaVar w_chib3P1float_1S("w_chib3P1float_1S", SigmaScaleFormula, RooArgList(m_chib3P1float_1S, PhotonSigmaScale3P));
  RooFormulaVar w_chib3P2float_1S("w_chib3P2float_1S", SigmaScaleFormula, RooArgList(m_chib3P2float_1S, PhotonSigmaScale3P));

  RooCBShape chib3P1_1S      = RooCBShape ("chib3P1_1S","chib3P1_1S",invm1S,m_chib3P1float_1S,w_chib3P1float_1S,alpha_3J,n_3J);
  RooCBShape chib3P2_1S     = RooCBShape ("chib3P2_1S","chib3P2_1S",invm1S,m_chib3P2float_1S,w_chib3P2float_1S,alpha_3J,n_3J);

  RooRealVar n3PnJ_1S= RooRealVar("N_{#chi_{b}(3P),1S}","n3PnJ_1S",n3PnJ_1S_START,0,100000);
  n3PnJ_1S.setVal(n3PnJ_1S_START);
  RooAddPdf chib3P_sigInd_1S= RooAddPdf("chib3P_sigInd_1S","chib3P_sigInd_1S",RooArgList(chib3P1_1S,chib3P2_1S),RooArgList(fracJ1));

  RooFormulaVar chib3P_nevt_1S("chib3P_nevt_1S", "@1*@0+(1-@1)*@0", RooArgList(n3PnJ_1S,fracJ1));
  RooAddPdf chib3P_sig_1S= RooAddPdf("chib3P_sig_1S","chib3P_sig_1S",RooArgList(chib3P_sigInd_1S),RooArgList(chib3P_nevt_1S));

  RooRealVar background_nevt1S = RooRealVar("N_{bkg}1S","N_{bkg}1S",background_nevt1S_START,0,100000);
  background_nevt1S.setConstant();

  RooAddPdf modelPdf1S_3P= RooAddPdf("modelPdf1S_3P","modelPdf1S_3P",RooArgList(chib3P_sig_1S,*background1S),RooArgList(chib3P_nevt_1S,background_nevt1S));

  RooDataSet* toySet =(RooDataSet*)modelPdf1S_3P.generate(RooArgSet(invm1S));

  toySet->Print();




	 RooAbsReal* nll1S=modelPdf1S_3P.createNLL(*toySet);

	 RooArgSet simNLL_Set2("simNLL_Set2");
	 simNLL_Set2.add(*nll1S);

	 RooAddition simNLL = RooAddition("add","add",simNLL_Set2);
     double Nllbefore3P=simNLL.getVal();

	 RooMinuit* minuit = new RooMinuit(simNLL);

	 minuit->setStrategy(2);
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



     double fmchib3P_1S=m_chib3P.getVal();
     double err_fmchib3P_1S=m_chib3P.getError();

     double OneSigmaCL = 0.682689492137;

     double sigmaCalc=0.005;
     for(int i = 1; i < 10000; i++){
    	 sigmaCalc+=i*0.000001;double Int1 = NormalizedIntegral(&chib3P_sig_1S, invm1S, (fmchib3P_1S-Ymass1S)*PhotonMassScale3P.getVal()+Ymass1S-sigmaCalc, (fmchib3P_1S-Ymass1S)*PhotonMassScale3P.getVal()+Ymass1S+sigmaCalc);if(Int1 > OneSigmaCL) break;
     }
     double fsigma3_1S=sigmaCalc;
     cout<<"fsigma3_1S = "<<fsigma3_1S<<endl;

     double fnchib3P_1S = chib3P_nevt_1S.getVal();
     double fnbg1S = background_nevt1S.getVal();

     double ThreePsig_1S=fnchib3P_1S/TMath::Sqrt(fnchib3P_1S+fnbg1S);

     double nSigma=1.5;

     double rchib3P_1S=fnchib3P_1S*NormalizedIntegral(&chib3P_sig_1S, invm1S, fmchib3P_1S-nSigma*fsigma3_1S, fmchib3P_1S+nSigma*fsigma3_1S);
     double rbgchib3P_1S=fnbg1S*NormalizedIntegral(&*background1S, invm1S, fmchib3P_1S-nSigma*fsigma3_1S, fmchib3P_1S+nSigma*fsigma3_1S);

     double width=PhotonSigmaScale3P_START*(m_chib3P_START-Ymass1S)*PhotonMassScale3P_START;

  double totalEventsInFit1S_plot = chib3P_nevt_1S.getVal()+background_nevt1S.getVal();


  if(iToy==1){
  TCanvas* ChibCanvas1S = new TCanvas("#chi_{b} 1S invariant mass","#chi_{b} 1S invariant mass",1900, 800);
  ChibCanvas1S->Divide(1);
  ChibCanvas1S->SetFillColor(kWhite);
  ChibCanvas1S->cd(1);
  gPad->SetRightMargin(0.3);/////
  gPad->SetFillColor(kWhite);

  double FontSize=0.0325;
  char PlotRangeSteps[200];
  char saveName[200];
  char plotYtitle[200];
  sprintf(PlotRangeSteps,"Y1SAll");
  double linewidth = 1;

  double BinWidth=15.;//MeV
  int FrameBins=InvSpace*1000./BinWidth;

  RooPlot* frame1S= invm1S.frame(FrameBins);
  frame1S->SetTitle("#chi_{b} invariant mass, Y1S decay");
  frame1S->GetXaxis()->SetLimits(invm_min1S,invm_max1S);
  toySet->plotOn(frame1S);
  sprintf(plotYtitle,"Events per %1.1f MeV",BinWidth);
  frame1S->SetYTitle(plotYtitle);
  frame1S->SetTitleOffset(1.15,"X");
  frame1S->SetTitleOffset(0.825,"Y");
  modelPdf1S_3P.plotOn(frame1S,Range(PlotRangeSteps),LineWidth(linewidth),Normalization(totalEventsInFit1S_plot,2));;
  double chi21S = frame1S->chiSquare();
  /////////// Two CB adventure //////////////////
  modelPdf1S_3P.plotOn(frame1S, Components("chib3P1_1S"), LineStyle(2),LineColor(kGreen),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
  modelPdf1S_3P.plotOn(frame1S, Components("chib3P2_1S"), LineStyle(2),LineColor(kRed),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
  modelPdf1S_3P.plotOn(frame1S, Components("background1S"), LineStyle(2),LineColor(kBlack),LineWidth(linewidth),Range(PlotRangeSteps),Normalization(totalEventsInFit1S_plot,2));
  modelPdf1S_3P.plotOn(frame1S,Range(PlotRangeSteps),LineWidth(linewidth),Normalization(totalEventsInFit1S_plot,2));

  modelPdf1S_3P.paramOn(frame1S, Layout(0.725,0.9875,0.9), Format("NE",AutoPrecision(2)));
  frame1S->SetTitle(0);
  frame1S->SetMinimum(0);
  frame1S->Draw();


  double xText=10.6;
  double highestText=frame1S->GetMaximum();
  double deltaText=0.06;
  char text[200];

  cout<<"DRAW LATEX"<<endl;
  sprintf(text,"S/B #chi_{b}(3P) = %1.0f / %1.0f, #chi_{b}(3P)_{sig.} = %1.3f",rchib3P_1S,rbgchib3P_1S,ThreePsig_1S);
  TLatex text3 = TLatex(xText,highestText*(0.95-6*deltaText),text);
  text3.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
  text3.Draw( "same" )                                                                                                                                                                                                                                                 ;
  sprintf(text,"M_{#chi_{b}(3P)} =  %1.4f #pm %1.4f",fmchib3P_1S,err_fmchib3P_1S);
  TLatex text0 = TLatex(xText,highestText*(0.95-7*deltaText),text)                                                                                                                            ;
  text0.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
  text0.Draw( "same" )                                                                                                                                                                                                                                                 ;
  sprintf(text,"S_{tot} = %d, B_{tot} = %d",nSig,nBkg);
  TLatex text4 = TLatex(xText,highestText*(0.95-8*deltaText),text);                                                                                                                                                                        ;
  text4.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
  text4.Draw( "same" )                                                                                                                                                                                                                                                 ;
  sprintf(text,"#chi^{2}/ndf = %1.4f",chi21S);
  TLatex text6 = TLatex(xText,highestText*(0.95-5*deltaText),text)                                                                                                                                                                                          ;
  text6.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
  text6.Draw( "same" )                                                                                                                                                                                                                                                 ;

  sprintf(saveName,"ToyFigures/%s/SamplePlots/Toy_InvMass_chib_S%d_B%d.pdf",FitID,nSig,nBkg);
  ChibCanvas1S->SaveAs(saveName);

  delete ChibCanvas1S;
  delete frame1S;

  }



  S_in_treeBuffer+=nSig;
  S_out_treeBuffer+=fnchib3P_1S;
  B_in_treeBuffer+=nBkg;
  S_treeBuffer+=rchib3P_1S;
  B_treeBuffer+=rbgchib3P_1S;
  SB_treeBuffer+=rchib3P_1S/rbgchib3P_1S;
  Sig_treeBuffer+=ThreePsig_1S;
  DeltaM_treeBuffer+=fmchib3P_1S-m_chib3P_START;
  Merr_treeBuffer+=err_fmchib3P_1S;
  Width_treeBuffer+=width;


     delete nll1S;
     delete minuit;
  delete background1S;
  delete RooBkgToyHist1S_Mixer;
  delete toySet;

  }

  S_in_tree=S_in_treeBuffer/nToy;
  S_out_tree=S_out_treeBuffer/nToy;
  B_in_tree=B_in_treeBuffer/nToy;
  S_tree=S_treeBuffer/nToy;
  B_tree=B_treeBuffer/nToy;
  SB_tree=SB_treeBuffer/nToy;
  Sig_tree=Sig_treeBuffer/nToy;
  DeltaM_tree=DeltaM_treeBuffer/nToy;
  Merr_tree=Merr_treeBuffer/nToy;
  Width_tree=Width_treeBuffer/nToy;

  Results->Fill();
}
}
  resultsFile->Write();



	return 0;
}

