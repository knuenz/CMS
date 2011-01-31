
///////
// This program is meant to extract the background shapes from the data side bands. 
// Only accepts one file, weight doesn't matter.
// \author Lindsey Gray (UW Madison)
//////

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

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TObject.h"

int main(int argc, char** argv) {
  using namespace JPsiPolarization;

  if(argc != 4) {
    std::cout << "Usage: ./injectAcceptanceMaps /path/to/promptAcceptance.root /path/to/nonPromptAcceptance.root /path/to/promptRecoEff.root" << std::endl;
    std::cout << "Two files accepted as input." << std::endl;
    return 1;
  }  

  bool folded(true);


  RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
  RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5); 
  RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
//  RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",-180,180);
  RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);
//  RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",-180,180);

  TFile *output = new TFile("jPsiFit_PR_FSR_GA07RE17smsm0NEW.root","UPDATE");
  TFile *promptMaps = NULL, *nonPromptMaps = NULL, *promptRecoEff = NULL;
 
  std::cout << "Loading acceptance maps from: " << argv[1] << " and " << argv[2] << std::endl;
  std::cout << "Loading reco efficienciy maps from: " << argv[3] << std::endl;
  promptMaps = TFile::Open(argv[1],"UPDATE");
  nonPromptMaps = TFile::Open(argv[2],"UPDATE");
  promptRecoEff = TFile::Open(argv[3],"UPDATE");


  for(int yBin = 0; yBin < jpsi::kNbRapForPTBins-3; ++yBin) {
 	  for(int ptBin = 1; ptBin < jpsi::kNbPTBins[yBin+1]; ++ptBin) {

 		  if(ptBin==7 && yBin==0)continue;

      std::stringstream binName,cutString;
      binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
      cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
		<< "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];

      CompositeModelBuilder* modelHX = new CompositeModelBuilder("HX");
      CompositeModelBuilder* modelCS = new CompositeModelBuilder("CS");
           
      TDirectory *current; 
      if(!(current = output->GetDirectory(binName.str().c_str()))) {
	current = output->mkdir(binName.str().c_str());
      }
      
      modelHX->setUseLifetime(false);
      modelHX->setUsePol(false);
      modelHX->setUseMass(false);
      modelHX->setUseBkg(false);
      
      modelCS->setUseLifetime(false);
      modelCS->setUsePol(false);
      modelCS->setUseMass(false);
      modelCS->setUseBkg(false);
      
//      modelCS->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);

      TH2F *promptMap_CS,*nonPromptMap_CS;
      TH2F *promptMap_HX,*nonPromptMap_HX;
      TH2F *promptEff_CS,*promptEff_HX;


      std::stringstream nameHXp,nameCSp;
      std::stringstream nameHXnp,nameCSnp;      
      
      nameHXp << "hAcc2D_HX_pT" << ptBin+1 << "_rap" << yBin+1;
      nameCSp << "hAcc2D_CS_pT" << ptBin+1 << "_rap" << yBin+1;
      nameHXnp << "hAcc2D_HX_pT" << ptBin+1 << "_rap" << yBin+1 << "_NP";
      nameCSnp << "hAcc2D_CS_pT" << ptBin+1 << "_rap" << yBin+1 << "_NP";

      cout<<nameCSp.str().c_str()<<endl;

      promptMap_HX = (TH2F*)promptMaps->Get(nameHXp.str().c_str());
      nonPromptMap_HX = (TH2F*)nonPromptMaps->Get(nameHXp.str().c_str());
      promptEff_HX = (TH2F*)promptRecoEff->Get(nameHXp.str().c_str());

      promptMap_CS = (TH2F*)promptMaps->Get(nameCSp.str().c_str());
      nonPromptMap_CS= (TH2F*)nonPromptMaps->Get(nameCSp.str().c_str());
      promptEff_CS = (TH2F*)promptRecoEff->Get(nameCSp.str().c_str());

      //      TH2F *promptCorrectionCS, *promptCorrectionHX;

//      cout<<promptMap_CS->GetBinContent(10,10)<<"+-"<<promptMap_CS->GetBinError(10,10)<<endl;
//      cout<<promptEff_CS->GetBinContent(10,10)<<"+-"<<promptEff_CS->GetBinError(10,10)<<endl;

///////////////////////// REBIN HISTORGRAMS TO MULTIPLY ACC*EFF /////////////////////////////////////////////

      cout<<promptEff_CS->GetNbinsX()<<" "<<promptEff_CS->GetNbinsY()<<endl;

      double accPrompt_CS[promptMap_CS->GetNbinsX()+1][promptMap_CS->GetNbinsY()+1];
      double accPrompt_HX[promptMap_HX->GetNbinsX()+1][promptMap_HX->GetNbinsY()+1];
      double effPrompt_CS[promptEff_CS->GetNbinsX()+1][promptEff_CS->GetNbinsY()+1];
      double effPrompt_HX[promptEff_HX->GetNbinsX()+1][promptEff_HX->GetNbinsY()+1];

      for(int costhbin=1;costhbin<promptMap_CS->GetNbinsX()+1;costhbin++){
          for(int phibin=1;phibin<promptMap_CS->GetNbinsY()+1;phibin++){
        	  accPrompt_CS[costhbin][phibin]=promptMap_CS->GetBinContent(costhbin,phibin);
        	  accPrompt_HX[costhbin][phibin]=promptMap_HX->GetBinContent(costhbin,phibin);
      }

      }
      for(int costhbin=1;costhbin<promptEff_CS->GetNbinsX()+1;costhbin++){
           for(int phibin=1;phibin<promptEff_CS->GetNbinsY()+1;phibin++){
              effPrompt_CS[costhbin][phibin]=promptEff_CS->GetBinContent(costhbin,phibin);
              effPrompt_HX[costhbin][phibin]=promptEff_HX->GetBinContent(costhbin,phibin);
      }
      }

      int binfactor_costh;
      int binfactor_phi;
      int binscosth;
      int binsphi;

      if (promptMap_CS->GetNbinsX()!=promptEff_CS->GetNbinsX()){
      cout<<"REBINNING"<<endl;
      //////////// if there is a gemeinsamer teiler, change following 4 parameters!! //////////////////
      binfactor_costh=5;//promptMap_CS->GetNbinsX();
      binfactor_phi=5;//promptMap_CS->GetNbinsY();
      binscosth=promptMap_CS->GetNbinsX()*binfactor_costh;
      binsphi=promptMap_CS->GetNbinsY()*binfactor_phi;
      }
      else{
      cout<<"no REBINNING necessary"<<endl;
      binfactor_costh=1;
      binfactor_phi=1;
      binscosth=promptMap_CS->GetNbinsX()*binfactor_costh;
      binsphi=promptMap_CS->GetNbinsY()*binfactor_phi;
      }

      TH2F *promptMapRebin_CS = new TH2F("promptMapRebin_CS","promptMapRebin_CS",binscosth,-1,1,binsphi,0,90);
      TH2F *promptMapRebin_HX = new TH2F("promptMapRebin_HX","promptMapRebin_HX",binscosth,-1,1,binsphi,0,90);
      TH2F *promptEffRebin_CS = new TH2F("promptEffRebin_CS","promptEffRebin_CS",binscosth,-1,1,binsphi,0,90);
      TH2F *promptEffRebin_HX = new TH2F("promptEffRebin_HX","promptEffRebin_HX",binscosth,-1,1,binsphi,0,90);

      for(int costhbin=1;costhbin<promptMap_CS->GetNbinsX()+1;costhbin++){
          for(int phibin=1;phibin<promptMap_CS->GetNbinsY()+1;phibin++){
        	  for(int i=1;i<binfactor_costh+1;i++){
            	  for(int j=1;j<binfactor_phi+1;j++){

        		  promptMapRebin_CS->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,accPrompt_CS[costhbin][phibin]);
        		  promptMapRebin_HX->SetBinContent((costhbin-1)*binfactor_costh+i,(phibin-1)*binfactor_phi+j,accPrompt_HX[costhbin][phibin]);

            	  }
        	  }
          }
 	  }

/*      cout<<promptMap_CS->GetBinContent(1,1)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(1,1)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(1,2)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(1,3)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(1,4)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(1,5)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(2,1)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(2,2)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(2,3)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(2,4)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(2,5)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(2,6)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(5,5)<<endl;
      cout<<promptMapRebin_CS->GetBinContent(6,5)<<endl;

      promptMapRebin_CS->Print("all");
*/
//      promptMap_CS->Multiply(promptEff_CS);
//      promptMap_HX->Multiply(promptEff_HX);

//      cout<<"fffff"<<endl;
//      promptEff_CS->Print("all");
//      cout<<"fffff"<<endl;

 //     promptMapRebin_CS->Print("all");

//      cout<<promptEff_CS->GetNbinsX()<<" "<<promptEff_CS->GetNbinsY()<<endl;
//      cout<<promptMapRebin_CS->GetNbinsX()<<" "<<promptMapRebin_CS->GetNbinsY()<<endl;


      cout<<"Multiplying"<<endl;
      promptMapRebin_CS->Multiply(promptEff_CS);
      promptMapRebin_HX->Multiply(promptEff_HX);

      promptMap_CS->Multiply(promptEff_CS);
      promptMap_HX->Multiply(promptEff_HX);


//      cout<<promptMapRebin_CS->GetNbinsX()<<" "<<promptMapRebin_CS->GetNbinsY()<<endl;

//      promptMapRebin_CS->Print("all");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      
      modelCS->setPromptAccHist(promptMap_CS);
      modelCS->setNonPromptAccHist(nonPromptMap_CS);

      modelHX->setPromptAccHist(promptMap_HX);
      modelHX->setNonPromptAccHist(nonPromptMap_HX);

      modelCS->initModel(JpsiMass,Jpsict,costh_CS,phi_CS);
      modelHX->initModel(JpsiMass,Jpsict,costh_HX,phi_HX);

      modelCS->saveParameters(*current);
      modelHX->saveParameters(*current);
      
      delete modelHX;
      delete modelCS;      
    }
  }

  output->Write();
  output->Close();
  promptMaps->Close();
  nonPromptMaps->Close();
  
  delete promptMaps;
  delete nonPromptMaps;
  delete output;
  return 0;
}
