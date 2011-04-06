#include <iostream>
#include <sstream>
#include <cstring>

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
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPolarizationPdf.h"
#include "RooGaussian.h"
#include "RooPolarizationConstraint.h"
#include "RooProdPdf.h"
#include "RooNLLVar.h"
#include "RooMinimizer.h"
#include "RooNumIntConfig.h"
#include "RooCmdArg.h"
#include "RooLinkedList.h"
#include "RooAbsArg.h"
#include "RooHistFunc.h"
#include "RooFormulaVar.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TText.h"


int main(int argc, char** argv) {
  using namespace JPsiPolarization;
  using namespace RooFit;


  char geomAcccharPR[200];
  char recoEffcharPR[200];
  char trigEffcharPR[200];
  char geomAcccharNP[200];
  char recoEffcharNP[200];
  char trigEffcharNP[200];

  char geomAcccharPRin[200];
  char recoEffcharPRin[200];
  char trigEffcharPRin[200];
  char geomAcccharNPin[200];
  char recoEffcharNPin[200];
  char trigEffcharNPin[200];
  char geomAcccharPRout[200];
  char recoEffcharPRout[200];
  char trigEffcharPRout[200];
  char geomAcccharNPout[200];
  char recoEffcharNPout[200];
  char trigEffcharNPout[200];

  char dirstruct[200];
  char cutstring[200];

  sprintf(cutstring,"AccEffCutUnRe70_");
  double realcut = 0.7;
  double cut=realcut;

  sprintf(dirstruct,"/scratch/knuenz/Polarization/RootInput/");
  sprintf(geomAcccharPR,"geomAccHistos_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected");
  sprintf(recoEffcharPR,"recoEffHistos_ATLASPT_20March2011_phiFolded_zeroBinsCorrected");
  sprintf(trigEffcharPR,"trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected");
  sprintf(geomAcccharNP,"geomAccHistos_NP_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected");
  sprintf(recoEffcharNP,"recoEffHistos_NP_ATLASPT_19March2011_phiFolded_zeroBinsCorrected");
  sprintf(trigEffcharNP,"trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected");

  sprintf(geomAcccharPRin,"%s%s.root",dirstruct,geomAcccharPR);
  sprintf(geomAcccharPRout,"%s%s%s.root",dirstruct,cutstring,geomAcccharPR);
  sprintf(geomAcccharNPin,"%s%s.root",dirstruct,geomAcccharNP);
  sprintf(geomAcccharNPout,"%s%s%s.root",dirstruct,cutstring,geomAcccharNP);
  sprintf(recoEffcharPRin,"%s%s.root",dirstruct,recoEffcharPR);
  sprintf(recoEffcharPRout,"%s%s%s.root",dirstruct,cutstring,recoEffcharPR);
  sprintf(recoEffcharNPin,"%s%s.root",dirstruct,recoEffcharNP);
  sprintf(recoEffcharNPout,"%s%s%s.root",dirstruct,cutstring,recoEffcharNP);
  sprintf(trigEffcharPRin,"%s%s.root",dirstruct,trigEffcharPR);
  sprintf(trigEffcharPRout,"%s%s%s.root",dirstruct,cutstring,trigEffcharPR);
  sprintf(trigEffcharNPin,"%s%s.root",dirstruct,trigEffcharNP);
  sprintf(trigEffcharNPout,"%s%s%s.root",dirstruct,cutstring,trigEffcharNP);


  TFile *geomAccFilePR = new TFile(geomAcccharPRin,"UPDATE");
  TFile *geomAccFileNP= new TFile(geomAcccharNPin,"UPDATE");
  TFile *recoEffFilePR = new TFile(recoEffcharPRin,"UPDATE");
  TFile *recoEffFileNP= new TFile(recoEffcharNPin,"UPDATE");
  TFile *trigEffFilePR = new TFile(trigEffcharPRin,"UPDATE");
  TFile *trigEffFileNP= new TFile(trigEffcharNPin,"UPDATE");

//  TFile *geomAccFilePRout = new TFile(geomAcccharPRout,"RECREATE");
//  TFile *geomAccFileNPout= new TFile(geomAcccharNPout,"RECREATE");
//  TFile *recoEffFilePRout = new TFile(recoEffcharPRout,"RECREATE");
//  TFile *recoEffFileNPout= new TFile(recoEffcharNPout,"RECREATE");
  TFile *trigEffFilePRout = new TFile(trigEffcharPRout,"RECREATE");
  TFile *trigEffFileNPout= new TFile(trigEffcharNPout,"RECREATE");


  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-3; ++yBin) {
 	  for(int ptBin =1; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {


  std::stringstream nameHXp,nameCSp;
  nameHXp << "hAcc2D_HX_pT" << ptBin+1 << "_rap" << yBin+1;
  nameCSp << "hAcc2D_CS_pT" << ptBin+1 << "_rap" << yBin+1;

  TH2F *AccMap_CS_PR,*AccMap_HX_PR;
  TH2F *RecoEffMap_CS_PR,*RecoEffMap_HX_PR;
  TH2F *TrigEffMap_CS_PR,*TrigEffMap_HX_PR;
  AccMap_CS_PR = (TH2F*)geomAccFilePR->Get(nameCSp.str().c_str());
  AccMap_HX_PR = (TH2F*)geomAccFilePR->Get(nameHXp.str().c_str());
  RecoEffMap_CS_PR = (TH2F*)recoEffFilePR->Get(nameCSp.str().c_str());
  RecoEffMap_HX_PR = (TH2F*)recoEffFilePR->Get(nameHXp.str().c_str());
  TrigEffMap_CS_PR = (TH2F*)trigEffFilePR->Get(nameCSp.str().c_str());
  TrigEffMap_HX_PR = (TH2F*)trigEffFilePR->Get(nameHXp.str().c_str());

  TH2F *AccMap_CS_NP,*AccMap_HX_NP;
  TH2F *RecoEffMap_CS_NP,*RecoEffMap_HX_NP;
  TH2F *TrigEffMap_CS_NP,*TrigEffMap_HX_NP;
  AccMap_CS_NP = (TH2F*)geomAccFileNP->Get(nameCSp.str().c_str());
  AccMap_HX_NP = (TH2F*)geomAccFileNP->Get(nameHXp.str().c_str());
  RecoEffMap_CS_NP = (TH2F*)recoEffFileNP->Get(nameCSp.str().c_str());
  RecoEffMap_HX_NP = (TH2F*)recoEffFileNP->Get(nameHXp.str().c_str());
  TrigEffMap_CS_NP = (TH2F*)trigEffFileNP->Get(nameCSp.str().c_str());
  TrigEffMap_HX_NP = (TH2F*)trigEffFileNP->Get(nameHXp.str().c_str());









  double RecoEffMap_CS_PR_double[RecoEffMap_CS_PR->GetNbinsX()+1][RecoEffMap_CS_PR->GetNbinsY()+1];
  double RecoEffMap_HX_PR_double[RecoEffMap_HX_PR->GetNbinsX()+1][RecoEffMap_HX_PR->GetNbinsY()+1];
  double TrigEffMap_CS_PR_double[TrigEffMap_CS_PR->GetNbinsX()+1][TrigEffMap_CS_PR->GetNbinsY()+1];
  double TrigEffMap_HX_PR_double[TrigEffMap_HX_PR->GetNbinsX()+1][TrigEffMap_HX_PR->GetNbinsY()+1];
  double RecoEffMap_CS_NP_double[RecoEffMap_CS_NP->GetNbinsX()+1][RecoEffMap_CS_NP->GetNbinsY()+1];
  double RecoEffMap_HX_NP_double[RecoEffMap_HX_NP->GetNbinsX()+1][RecoEffMap_HX_NP->GetNbinsY()+1];
  double TrigEffMap_CS_NP_double[TrigEffMap_CS_NP->GetNbinsX()+1][TrigEffMap_CS_NP->GetNbinsY()+1];
  double TrigEffMap_HX_NP_double[TrigEffMap_HX_NP->GetNbinsX()+1][TrigEffMap_HX_NP->GetNbinsY()+1];

  for(int costhbin=1;costhbin<RecoEffMap_CS_PR->GetNbinsX()+1;costhbin++){
      for(int phibin=1;phibin<RecoEffMap_CS_PR->GetNbinsY()+1;phibin++){
    	  RecoEffMap_CS_PR_double[costhbin][phibin]=RecoEffMap_CS_PR->GetBinContent(costhbin,phibin);
    	  RecoEffMap_HX_PR_double[costhbin][phibin]=RecoEffMap_HX_PR->GetBinContent(costhbin,phibin);
    	  TrigEffMap_CS_PR_double[costhbin][phibin]=TrigEffMap_CS_PR->GetBinContent(costhbin,phibin);
    	  TrigEffMap_HX_PR_double[costhbin][phibin]=TrigEffMap_HX_PR->GetBinContent(costhbin,phibin);
    	  RecoEffMap_CS_NP_double[costhbin][phibin]=RecoEffMap_CS_NP->GetBinContent(costhbin,phibin);
    	  RecoEffMap_HX_NP_double[costhbin][phibin]=RecoEffMap_HX_NP->GetBinContent(costhbin,phibin);
    	  TrigEffMap_CS_NP_double[costhbin][phibin]=TrigEffMap_CS_NP->GetBinContent(costhbin,phibin);
    	  TrigEffMap_HX_NP_double[costhbin][phibin]=TrigEffMap_HX_NP->GetBinContent(costhbin,phibin);
  }}

  int binfactor_costh;
  int binfactor_phi;
  int binscosth;
  int binsphi;

  TH2F *RecoEffMap_CS_PR_out = new TH2F(nameCSp.str().c_str(),nameCSp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *RecoEffMap_HX_PR_out = new TH2F(nameHXp.str().c_str(),nameHXp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *TrigEffMap_CS_PR_out = new TH2F(nameCSp.str().c_str(),nameCSp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *TrigEffMap_HX_PR_out = new TH2F(nameHXp.str().c_str(),nameHXp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *RecoEffMap_CS_NP_out = new TH2F(nameCSp.str().c_str(),nameCSp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *RecoEffMap_HX_NP_out = new TH2F(nameHXp.str().c_str(),nameHXp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *TrigEffMap_CS_NP_out = new TH2F(nameCSp.str().c_str(),nameCSp.str().c_str(),binscosth,-1,1,binsphi,0,90);
  TH2F *TrigEffMap_HX_NP_out = new TH2F(nameHXp.str().c_str(),nameHXp.str().c_str(),binscosth,-1,1,binsphi,0,90);

//  bool rebin(true);

//  if (rebin){
  cout<<"REBINNING"<<endl;
  binfactor_costh=2;
  binfactor_phi=2;
  binscosth=40;
  binsphi=36;



  for(int costhbin=1;costhbin<AccMap_CS_NP->GetNbinsX()+1;costhbin++){
      for(int phibin=1;phibin<AccMap_CS_NP->GetNbinsY()+1;phibin++){
    	  for(int i=1;i<binfactor_costh+1;i++){
        	  for(int j=1;j<binfactor_phi+1;j++){

    		  RecoEffMap_CS_PR_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,RecoEffMap_CS_PR_double[costhbin][phibin]);
    		  RecoEffMap_HX_PR_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,RecoEffMap_HX_PR_double[costhbin][phibin]);
    		  TrigEffMap_CS_PR_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,TrigEffMap_CS_PR_double[costhbin][phibin]);
    		  TrigEffMap_HX_PR_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,TrigEffMap_HX_PR_double[costhbin][phibin]);
    		  RecoEffMap_CS_NP_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,RecoEffMap_CS_NP_double[costhbin][phibin]);
    		  RecoEffMap_HX_NP_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,RecoEffMap_HX_NP_double[costhbin][phibin]);
    		  TrigEffMap_CS_NP_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,TrigEffMap_CS_NP_double[costhbin][phibin]);
    		  TrigEffMap_HX_NP_out->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,TrigEffMap_HX_NP_double[costhbin][phibin]);

        	  }
    	  }
      }
	  }




//  AccMap_CS_PR->Rebin2D(2,2);
//  AccMap_HX_PR->Rebin2D(2,2);
//  AccMap_CS_NP->Rebin2D(2,2);
//  AccMap_HX_NP->Rebin2D(2,2);

  AccMap_CS_PR->Multiply(RecoEffMap_CS_PR_out);
  AccMap_CS_PR->Multiply(TrigEffMap_CS_PR_out);
  AccMap_CS_NP->Multiply(RecoEffMap_CS_NP_out);
  AccMap_CS_NP->Multiply(TrigEffMap_CS_NP_out);
  AccMap_HX_PR->Multiply(RecoEffMap_HX_PR_out);
  AccMap_HX_PR->Multiply(TrigEffMap_HX_PR_out);
  AccMap_HX_NP->Multiply(RecoEffMap_HX_NP_out);
  AccMap_HX_NP->Multiply(TrigEffMap_HX_NP_out);

//  AccMap_HX_PR->Print("all");

//  double quasiNumEntriesCS=0;
//  double quasiNumEntriesCSnoot=0;
//  double quasiNumEntriesHX=0;
//  double quasiNumEntriesHXnoot=0;

  for(int costhbin=1;costhbin<TrigEffMap_CS_NP_out->GetNbinsX()+1;costhbin++){
      for(int phibin=1;phibin<TrigEffMap_CS_NP_out->GetNbinsY()+1;phibin++){

//    	  quasiNumEntriesCS=+AccMap_CS_PR->GetBinContent(costhbin,phibin);
//    	  quasiNumEntriesHX=+AccMap_HX_PR->GetBinContent(costhbin,phibin);

//    	  if(RecoEffMap_CS_PR->GetBinContent(costhbin,phibin)*TrigEffMap_CS_PR->GetBinContent(costhbin,phibin)*AccMap_CS_PR->GetBinContent(costhbin,phibin)<cut || RecoEffMap_CS_NP->GetBinContent(costhbin,phibin)*TrigEffMap_CS_NP->GetBinContent(costhbin,phibin)*AccMap_CS_NP->GetBinContent(costhbin,phibin)<cut)
		  if(AccMap_CS_PR->GetBinContent(costhbin,phibin)<cut || AccMap_CS_NP->GetBinContent(costhbin,phibin)<cut)
    	  {
    		  TrigEffMap_CS_PR_out->SetBinContent(costhbin,phibin,0);
    		  TrigEffMap_CS_NP_out->SetBinContent(costhbin,phibin,0);
    	  }
//    	  if(RecoEffMap_HX_PR->GetBinContent(costhbin,phibin)*TrigEffMap_HX_PR->GetBinContent(costhbin,phibin)*AccMap_HX_PR->GetBinContent(costhbin,phibin)<cut || RecoEffMap_HX_NP->GetBinContent(costhbin,phibin)*TrigEffMap_HX_NP->GetBinContent(costhbin,phibin)*AccMap_HX_NP->GetBinContent(costhbin,phibin)<cut)
    	  if(AccMap_HX_PR->GetBinContent(costhbin,phibin)<cut || AccMap_HX_NP->GetBinContent(costhbin,phibin)<cut)
    	  {
    		  TrigEffMap_HX_PR_out->SetBinContent(costhbin,phibin,0);
    		  TrigEffMap_HX_NP_out->SetBinContent(costhbin,phibin,0);
    	  }
      }
      }


//  }


//  geomAccFilePRout->Add(AccMap_CS_PR);
//  geomAccFilePRout->Add(AccMap_HX_PR);
//  geomAccFileNPout->Add(AccMap_CS_NP);
//  geomAccFileNPout->Add(AccMap_HX_NP);
//  recoEffFilePRout->Add(RecoEffMap_CS_PR);
//  recoEffFilePRout->Add(RecoEffMap_HX_PR);
//  recoEffFileNPout->Add(RecoEffMap_CS_NP);
//  recoEffFileNPout->Add(RecoEffMap_HX_NP);
  trigEffFilePRout->Add(TrigEffMap_CS_PR_out);
  trigEffFilePRout->Add(TrigEffMap_HX_PR_out);
  trigEffFileNPout->Add(TrigEffMap_CS_NP_out);
  trigEffFileNPout->Add(TrigEffMap_HX_NP_out);

 	  }}

//  geomAccFilePRout->Write();
//  geomAccFileNPout->Write();
//  recoEffFilePRout->Write();
//  recoEffFileNPout->Write();
  trigEffFilePRout->Write();
  trigEffFileNPout->Write();

//  geomAccFilePRout->Close();
//  geomAccFileNPout->Close();
//  recoEffFilePRout->Close();
//  recoEffFileNPout->Close();
  trigEffFilePRout->Close();
  trigEffFileNPout->Close();

//  geomAccFilePR->Close();
//  geomAccFileNP->Close();
//  recoEffFilePR->Close();
//  recoEffFileNP->Close();
  trigEffFilePR->Close();
  trigEffFileNP->Close();

/*
delete RecoEffMap_CS_PR_out;
delete RecoEffMap_HX_PR_out;
delete TrigEffMap_CS_PR_out;
delete TrigEffMap_HX_PR_out;
delete RecoEffMap_CS_NP_out;
delete RecoEffMap_HX_NP_out;
delete TrigEffMap_CS_NP_out;
delete TrigEffMap_HX_NP_out;
*/
//  delete geomAccFilePRout;
//  delete geomAccFileNPout;
//  delete recoEffFilePRout;
//  delete recoEffFileNPout;
  delete trigEffFilePRout;
  delete trigEffFileNPout;
//  delete geomAccFilePR;
//  delete geomAccFileNP;
//  delete recoEffFilePR;
//  delete recoEffFileNP;
  delete trigEffFilePR;
  delete trigEffFileNP;

  return 0;

}
