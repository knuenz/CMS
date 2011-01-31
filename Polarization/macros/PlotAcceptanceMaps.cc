#include <iostream>
#include <sstream>

//J/Psi common vars
#include "commonVar.h"

//Fitting routine
#include "CompositeModelBuilder.h"



// RooFit Includes
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TClass.h"
#include "TLegend.h"
#include "TSystem.h"


int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;


  gSystem->mkdir("Plots/AcceptanceMaps");

   gStyle->SetPalette(1);
   gStyle->SetPadLeftMargin(0.20);
//   gStyle->SetPadRightMargin(0.02);
   gStyle->SetTitleFillColor(kWhite);


  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5); 
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);
  RooRealVar MCweight("MCweight","MCweight",1);
  RooRealVar JpsiRap("JpsiRap","#nu",-2.3,2.3);
  RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);

  RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
  varlist.add(MCweight);

  
  TFile *fInput = new TFile("jPsiFit_PR_FSR_GA07RE17smsm0NEW.root","UPDATE");

  Char_t *fileNameInPR = "/scratch/knuenz/Polarization/RootInput/TTree_final_notrigger_MCprompt_Jpsi_Fall10_folded_.root";
//  Char_t *fileNameInPR = "/scratch/knuenz/Polarization/RootInput/TTree_red_PR_pseudo.root";

  TFile* fInPR = new TFile(fileNameInPR);
  TTree* dataTreesPR = (TTree*)fInPR->Get("data");
  RooDataSet *dataPR = new RooDataSet("dataPR","Supplied Data Prompt",varlist,Import(*dataTreesPR),WeightVar(MCweight));


  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-3; ++yBin) {
	  for(int ptBin = 1; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {

 		  if(ptBin==7 && yBin==0)continue;


      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];

  	std::cout << cutString.str() << std::endl;


	  char reducePR[200];
	  sprintf(reducePR,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);// && JpsictErr > 0.0004 && JpsictErr < 0.999999

	  RooDataSet* dataPRbin = (RooDataSet*)dataPR->reduce(reducePR);

	  TH2F * promptMap_CS_MC = (TH2F*)dataPRbin->createHistogram("promptMap_CS_MC",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol)));
	  TH2F * promptMap_HX_MC = (TH2F*)dataPRbin->createHistogram("promptMap_HX_MC",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol)));

      char AccCSPRTitle[200];
      sprintf(AccCSPRTitle,"Prompt CS Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char AccHXPRTitle[200];
      sprintf(AccHXPRTitle,"Prompt HX Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char AccCSNPTitle[200];
      sprintf(AccCSNPTitle,"Non Prompt CS Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char AccHXNPTitle[200];
      sprintf(AccHXNPTitle,"Non Prompt HX Acc*eff Map %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char MCCSPRTitle[200];
      sprintf(MCCSPRTitle,"MC Prompt CS %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);
      char MCHXPRTitle[200];
      sprintf(MCHXPRTitle,"MC Prompt HX %1.1f < pT < %1.1f, %1.1f < |y| < %1.1f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);


      char FilenameAcc[200];


      char DirectoryPath[200];
         sprintf(DirectoryPath,"/pt%d_rapidity%d",ptBin+1,yBin+1);

         TDirectory *InputDirectory = (TDirectory*)fInput->GetDirectory(DirectoryPath);

         CompositeModelBuilder* modelHX = new CompositeModelBuilder("HX");
         CompositeModelBuilder* modelCS = new CompositeModelBuilder("CS");


         modelHX->setUseLifetime(false);
         modelHX->setUsePol(false);
         modelHX->setUseMass(false);
         modelHX->setUseBkg(false);
         modelHX->loadParameters(*InputDirectory);
         modelHX->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);

         modelCS->setUseLifetime(false);
         modelCS->setUsePol(false);
         modelCS->setUseMass(false);
         modelCS->setUseBkg(false);
         modelCS->loadParameters(*InputDirectory);
         modelCS->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);


         RooHistPdf* promptCShistPdf_ = modelCS->getAcceptanceMaps()->prompt();
         RooHistPdf* promptHXhistPdf_ = modelHX->getAcceptanceMaps()->prompt();
//         RooHistPdf* nonpromptCShistPdf_ = modelCS->getAcceptanceMaps()->nonPrompt();
//         RooHistPdf* nonpromptHXhistPdf_ = modelHX->getAcceptanceMaps()->nonPrompt();


    	TH1 * hist_CS_PR_model = promptCShistPdf_->createHistogram("hist_CS_PR_model",costh_CS,Binning(250),YVar(phi_CS,Binning(250))) ;

			RooRealVar norm1("norm1","norm1",1);
     		RooDataHist accmodelCShist("accmodelCShist","accmodelCShist", RooArgList(costh_CS,phi_CS),hist_CS_PR_model);
     		RooHistPdf accmodelCShistpdf("accmodelCShistpdf","accmodelCShistpdf",RooArgList(costh_CS,phi_CS),accmodelCShist,1);
     		RooExtendPdf accmodelCS("accmodelCS","accmodelCS",accmodelCShistpdf,norm1);

        TH1 * hist_HX_PR_model = promptHXhistPdf_->createHistogram("hist_HX_PR_model",costh_HX,Binning(250),YVar(phi_HX,Binning(250))) ;

         	RooDataHist accmodelHXhist("accmodelHXhist","accmodelHXhist", RooArgList(costh_HX,phi_HX),hist_HX_PR_model);
         	RooHistPdf accmodelHXhistpdf("accmodelHXhistpdf","accmodelHXhistpdf",RooArgList(costh_HX,phi_HX),accmodelHXhist,1);
         	RooExtendPdf accmodelHX("accmodelHX","accmodelHX",accmodelHXhistpdf,norm1);



   	  TH1 * promptCShistPdf = promptCShistPdf_->createHistogram("promptCShistPdf",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;
   	  TH1 * promptHXhistPdf = promptHXhistPdf_->createHistogram("promptHXhistPdf",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol))) ;
//   	  TH1 * nonpromptCShistPdf = nonpromptCShistPdf_->createHistogram("nonpromptCShistPdf",costh_CS,Binning(jpsi::kNbBinsCosT),YVar(phi_CS,Binning(jpsi::kNbBinsPhiPol))) ;
//   	  TH1 * nonpromptHXhistPdf = nonpromptHXhistPdf_->createHistogram("nonpromptHXhistPdf",costh_HX,Binning(jpsi::kNbBinsCosT),YVar(phi_HX,Binning(jpsi::kNbBinsPhiPol))) ;


   	TCanvas* AcceptanceHistCanvas = new TCanvas("AcceptanceHistCanvas","AcceptanceHistCanvas",1600, 700);
   	AcceptanceHistCanvas->Divide(2);  AcceptanceHistCanvas->SetFillColor(kWhite);
   	AcceptanceHistCanvas->cd(1) ; promptCShistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptCShistPdf->SetStats(0); promptCShistPdf->GetXaxis()->SetTitle("cos #theta_{CS}"); promptCShistPdf->GetYaxis()->SetTitle("#phi_{CS} [deg]"); promptCShistPdf->SetTitle(AccCSPRTitle); promptCShistPdf->Draw("colz");
   	AcceptanceHistCanvas->cd(2) ; promptHXhistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptHXhistPdf->SetStats(0); promptHXhistPdf->GetXaxis()->SetTitle("cos #theta_{HX}"); promptHXhistPdf->GetYaxis()->SetTitle("#phi_{HX} [deg]"); promptHXhistPdf->SetTitle(AccHXPRTitle); promptHXhistPdf->Draw("colz");

   	      sprintf(FilenameAcc,"Plots/AcceptanceMaps/AccHists_rapidity%d_pt%d.png",yBin+1,ptBin+1);

   	   AcceptanceHistCanvas->SaveAs(FilenameAcc);
   	AcceptanceHistCanvas->Close();


    TCanvas* AcceptanceCanvas2 = new TCanvas("AcceptanceCanvas2","AcceptanceCanvas2",1600, 1400);
    AcceptanceCanvas2->Divide(2,2);  AcceptanceCanvas2->SetFillColor(kWhite);
    AcceptanceCanvas2->cd(1) ; promptCShistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptCShistPdf->SetStats(0); promptCShistPdf->GetXaxis()->SetTitle("cos #theta_{CS}"); promptCShistPdf->GetYaxis()->SetTitle("#phi_{CS} [deg]"); promptCShistPdf->SetTitle(AccCSPRTitle); promptCShistPdf->Draw("colz");
    AcceptanceCanvas2->cd(2) ; promptHXhistPdf->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptHXhistPdf->SetStats(0); promptHXhistPdf->GetXaxis()->SetTitle("cos #theta_{HX}"); promptHXhistPdf->GetYaxis()->SetTitle("#phi_{HX} [deg]"); promptHXhistPdf->SetTitle(AccHXPRTitle); promptHXhistPdf->Draw("colz");
    AcceptanceCanvas2->cd(3) ; promptMap_CS_MC->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptMap_CS_MC->SetStats(0); promptMap_CS_MC->GetXaxis()->SetTitle("cos #theta_{CS}"); promptMap_CS_MC->GetYaxis()->SetTitle("#phi_{CS} [deg]"); promptMap_CS_MC->SetTitle(MCCSPRTitle); promptMap_CS_MC->Draw("colz");
    AcceptanceCanvas2->cd(4) ; promptMap_HX_MC->GetYaxis()->SetTitleOffset(2); gPad->SetFillColor(kWhite); promptMap_HX_MC->SetStats(0); promptMap_HX_MC->GetXaxis()->SetTitle("cos #theta_{HX}"); promptMap_HX_MC->GetYaxis()->SetTitle("#phi_{HX} [deg]"); promptMap_HX_MC->SetTitle(MCHXPRTitle); promptMap_HX_MC->Draw("colz");

    sprintf(FilenameAcc,"Plots/AcceptanceMaps/AccMaps_rapidity%d_pt%d.png",yBin+1,ptBin+1);

    AcceptanceCanvas2->SaveAs(FilenameAcc);

    AcceptanceCanvas2->Close();



    RooPlot* AccFramePR_costh_CS = new RooPlot;
    AccFramePR_costh_CS = costh_CS.frame() ;
	dataPRbin->plotOn(AccFramePR_costh_CS,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	accmodelCS.plotOn(AccFramePR_costh_CS,LineWidth(2),Normalization(1.0));

	RooPlot* AccFramePR_phi_CS = new RooPlot;
    AccFramePR_phi_CS = phi_CS.frame() ;
	dataPRbin->plotOn(AccFramePR_phi_CS,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	accmodelCS.plotOn(AccFramePR_phi_CS,LineWidth(2),Normalization(1.0));

    RooPlot* AccFramePR_costh_HX = new RooPlot;
    AccFramePR_costh_HX = costh_HX.frame() ;
 	dataPRbin->plotOn(AccFramePR_costh_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 	accmodelHX.plotOn(AccFramePR_costh_HX,LineWidth(2),Normalization(1.0));

 	RooPlot* AccFramePR_phi_HX = new RooPlot;
    AccFramePR_phi_HX = phi_HX.frame() ;
 	dataPRbin->plotOn(AccFramePR_phi_HX,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 	accmodelHX.plotOn(AccFramePR_phi_HX,LineWidth(2),Normalization(1.0));

    TCanvas* Acc1DCanvas = new TCanvas("Acc1DCanvas","Acc1DCanvas",1600,1000);
    Acc1DCanvas->Divide(2,2);  Acc1DCanvas->SetFillColor(kWhite);
    Acc1DCanvas->cd(1) ; gPad->SetFillColor(kWhite); AccFramePR_costh_CS->Draw();
    Acc1DCanvas->cd(2) ; gPad->SetFillColor(kWhite); AccFramePR_phi_CS->Draw();
    Acc1DCanvas->cd(3) ; gPad->SetFillColor(kWhite); AccFramePR_costh_HX->Draw();
    Acc1DCanvas->cd(4) ; gPad->SetFillColor(kWhite); AccFramePR_phi_HX->Draw();

    sprintf(FilenameAcc,"Plots/AcceptanceMaps/AccHistsCS1D_rapidity%d_pt%d.png",yBin+1,ptBin+1);

    Acc1DCanvas->SaveAs(FilenameAcc);
    Acc1DCanvas->Close();

    delete AcceptanceCanvas2;
    delete AcceptanceHistCanvas;
    delete Acc1DCanvas;
    delete AccFramePR_costh_CS;
    delete AccFramePR_phi_CS;
    delete AccFramePR_costh_HX;
    delete AccFramePR_phi_HX;

    }
  }
  
  delete fInput;

  return 0;
}
