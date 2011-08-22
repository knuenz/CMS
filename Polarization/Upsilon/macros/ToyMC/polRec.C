#include "Riostream.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"


// WARNING!! match the pT-eta binning of the efficiency map with the borders
// of the regions of different pT-eta cuts



enum { loose, tight };

bool isMuonInAcceptance(int iCut, double pT, double eta){

  double etaBorderHLT[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts
  double pTBorderHLT[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};


  double minPT_HLT;
  bool decision = kFALSE;

  //loop over higher pT muon
  for(int iEta = 0; iEta < 3; iEta++){
    if(fabs(eta) > etaBorderHLT[iCut][iEta] && fabs(eta) < etaBorderHLT[iCut][iEta+1]){
      minPT_HLT = (pTBorderHLT[iCut][iEta+1]-pTBorderHLT[iCut][iEta]) / (etaBorderHLT[iCut][iEta+1]-etaBorderHLT[iCut][iEta]) * (fabs(eta) - etaBorderHLT[iCut][iEta]) + pTBorderHLT[iCut][iEta];
      break;
    }
    else if(fabs(eta) > etaBorderHLT[iCut][3])
      minPT_HLT = 1000.; //reject all events with |eta| > 2.4 (or 2.2, ...)
  }

  if(pT > minPT_HLT)
    decision = kTRUE;

  return decision;
}


const int nbin_pT     = 400;
const int nbin_abseta = 28;

/*
double singleLeptonEfficiency( double& pT, double& eta ) {
   double eff = 0.;
   if ( (   fabs(eta) < 0.8 && pT > 3.75
          || fabs(eta) > 0.8 && fabs(eta) < 1.6 && pT > 3.5

          || fabs(eta) > 1.6 && fabs(eta) < 2.4 && pT > 3.0  )
      ) eff = 1.;
   return eff;
}
*/

double singleLeptonEfficiency1( double& pT, double& eta, int FidCuts ) {

	double eff=1;

//	  if(FidCuts==1 && !isMuonInAcceptance( loose, pT, eta) ) eff = 0.;
//	  if(FidCuts==2 && !isMuonInAcceptance( tight, pT, eta) ) eff = 0.;

  return eff;
}


double singleLeptonEfficiency2( double& pT, double& eta ,int FidCuts) {

  const double mu_pT_min = 4.0;
  const double smoothcutpar = 3.0;
  double eff = 1. / ( 1. + exp(-smoothcutpar*(pT - mu_pT_min)) );
  if ( fabs(eta) > 2.4 ) eff = 0.;

//  if(FidCuts==1 && !isMuonInAcceptance( loose, pT, eta) ) eff = 0.;
//  if(FidCuts==2 && !isMuonInAcceptance( tight, pT, eta) ) eff = 0.;

  return eff;

}


/// Matt's efficiencies
double singleLeptonEfficiency3( double& pT, double& eta ,int FidCuts) {

  double c0 = 0.878;
  double c1 = 3.894;
  double c2 = 0.957;

  if ( fabs(eta) > 0.8 && fabs(eta) < 1.2 ) {
    c0 = 0.839;
    c1 = 3.860;
    c2 = 0.512;
  }

  if ( fabs(eta) > 1.2 && fabs(eta) < 1.6 ) {
    c0 = 0.882;
    c1 = 2.984;
    c2 = 0.405;
  }

  if ( fabs(eta) > 1.6 && fabs(eta) < 2.0 ) {
    c0 = 0.839;
    c1 = 2.280;
    c2 = 1.398;
  }

  if ( fabs(eta) > 2.0 && fabs(eta) < 2.4 ) {
    c0 = 0.713;
    c1 = 0.0;
    c2 = 0.0;
  }

  double eff = 0.5*c0 * (1. + TMath::Erf( (pT - c1) / ( sqrt(2.) * c2) ) );

//  if(FidCuts==1 && !isMuonInAcceptance( loose, pT, eta) ) eff = 0.;
//  if(FidCuts==2 && !isMuonInAcceptance( tight, pT, eta) ) eff = 0.;

  return eff;
}






/*
double singleLeptonEfficiency( double& pT, double& eta ) {
   double eff = 0.;
   if (  pT > 2.0 && fabs(eta) < 2.4 ) eff = 1.;
   return eff;
}
*/


void polRec(int nEff=1,
			int FidCuts=0,
			Char_t *dirstruct = "ToyDirectory_Default"){

  gROOT->Reset();

  char filename [500];

  // get ntuple of generated data from file

  sprintf(filename,"%s/genData.root",dirstruct);
  TFile* genFile = new TFile(filename,"READ");
  TTree* genData = (TTree*)genFile->Get("genData");


  // create output data file

  sprintf(filename,"%s/data.root",dirstruct);
  TFile* dataFile = new TFile(filename, "RECREATE", "dataFile");


  // output ntuple

  TTree* data = new TTree("data","data");


  // histograms to be output to the opened file

  // signal+background and background distributions in one frame.
  // chosen frame: PX (almost uniform distribution in phi -> less phi bins, less fluctuations)

  const int nbinth = 40;
  const int nbinph = 12;

  TH2D* total_PX      = new TH2D( "total_PX",      "", nbinth, -1., 1., nbinph, -180., 180. );
  TH2D* background_PX = new TH2D( "background_PX", "", nbinth, -1., 1., nbinph, -180., 180. );

  TH1D* backgroundFraction = new TH1D( "backgroundFraction", "", 1, 0., 1. );
  TH1D* isBG_distribution = new TH1D( "isBG_distribution", "", 2, 0., 2. );

  sprintf(filename,"%s/efficiency.root",dirstruct);
  TFile* efficiencyFile = new TFile(filename, "RECREATE", "efficiencyFile");

  TH2D* lepton_pT_vs_eta_rec = new TH2D( "lepton_pT_vs_eta_rec", "", nbin_abseta, 0.0, 5.6, nbin_pT, 0., 100. );
  TH2D* lepton_pT_vs_eta_gen = new TH2D( "lepton_pT_vs_eta_gen", "", nbin_abseta, 0.0, 5.6, nbin_pT, 0., 100. );
  TH2D* leptonEfficiency     = new TH2D( "leptonEfficiency",     "", nbin_abseta, 0.0, 5.6, nbin_pT, 0., 100. );

  lepton_pT_vs_eta_rec->Sumw2(); lepton_pT_vs_eta_gen->Sumw2(); leptonEfficiency->Sumw2();


  // structure of input ntuple

  TLorentzVector* lepP_gen = 0;        genData->SetBranchAddress("lepP",    &lepP_gen);
  TLorentzVector* lepN_gen = 0;        genData->SetBranchAddress("lepN",    &lepN_gen);

  double costh_CS;  genData->SetBranchAddress("costh_CS",     &costh_CS);
  double phi_CS;    genData->SetBranchAddress("phi_CS",       &phi_CS  );

  double costh_HX;  genData->SetBranchAddress("costh_HX",     &costh_HX);
  double phi_HX;    genData->SetBranchAddress("phi_HX",       &phi_HX  );

  double costh_PX;  genData->SetBranchAddress("costh_PX",     &costh_PX);
  double phi_PX;    genData->SetBranchAddress("phi_PX",       &phi_PX  );

  int isBG;         genData->SetBranchAddress("isBG",         &isBG        );


  // structure of output ntuple

  TLorentzVector* lepP = new TLorentzVector(0.,0.,0.,0.);
  data->Branch("lepP", "TLorentzVector", &lepP);
  TLorentzVector* lepN = new TLorentzVector(0.,0.,0.,0.);
  data->Branch("lepN", "TLorentzVector", &lepN);


  // for the calculation of the bin boundaries from the generated data themselves:

  double pTdilepton_min = 10000.;
  double pTdilepton_max = 0.;
  double rapdilepton_min = 10000;
  double rapdilepton_max = 0.;
  double massdilepton_min = 10000;
  double massdilepton_max = 0.;

  // loop over events in the input ntuple

  int numEvts = int( genData->GetEntries() );


  cout << endl;
  cout << "Reading " << numEvts << " dilepton events"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: " << endl;

  int n_step = numEvts/5;
  int n_step_=1;

  int rejected=0;

  for ( int evtNumber = 0; evtNumber < numEvts; evtNumber++ ) {

	    if ((evtNumber+1)%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    genData->GetEvent( evtNumber );


    double lepP_pT  = lepP_gen->Pt();
    double lepN_pT  = lepN_gen->Pt();

    double lepP_eta  = lepP_gen->PseudoRapidity();
    double lepN_eta  = lepN_gen->PseudoRapidity();

    lepton_pT_vs_eta_gen->Fill( fabs(lepP_eta), lepP_pT );
    lepton_pT_vs_eta_gen->Fill( fabs(lepN_eta), lepN_pT );

    bool isMuPinAcc = isMuonInAcceptance( FidCuts-1 , lepP_pT, lepP_eta) ;
    bool isMuNinAcc = isMuonInAcceptance( FidCuts-1 , lepN_pT, lepN_eta) ;

    bool FidCutRejection(false);
    if(FidCuts>0) if(!isMuPinAcc || !isMuNinAcc) FidCutRejection=true;

    double effP;
    double effN;
    if(nEff==1) {effP=singleLeptonEfficiency1( lepP_pT, lepP_eta, FidCuts ); effN=singleLeptonEfficiency1( lepN_pT, lepN_eta, FidCuts );}
    if(nEff==2) {effP=singleLeptonEfficiency2( lepP_pT, lepP_eta, FidCuts ); effN=singleLeptonEfficiency2( lepN_pT, lepN_eta, FidCuts );}
    if(nEff==3) {effP=singleLeptonEfficiency3( lepP_pT, lepP_eta, FidCuts ); effN=singleLeptonEfficiency3( lepN_pT, lepN_eta, FidCuts );}

    double rndmeffP = gRandom->Uniform(1.);
    double rndmeffN = gRandom->Uniform(1.);


    if ( rndmeffP < effP ) lepton_pT_vs_eta_rec->Fill( fabs(lepP_eta), lepP_pT );
    if ( rndmeffN < effN ) lepton_pT_vs_eta_rec->Fill( fabs(lepN_eta), lepN_pT );



    if ( rndmeffP > effP || rndmeffN > effN || FidCutRejection) {rejected++; continue; }

    total_PX->Fill( costh_PX, phi_PX );

    isBG_distribution->Fill( isBG );

    if ( isBG ) background_PX->Fill( costh_PX, phi_PX );

    *lepP = *lepP_gen;
    *lepN = *lepN_gen;

    data->Fill();

  } // end of RD ntuple loop

  cout << endl << endl;



  // background fraction

  double f_BG = isBG_distribution->GetMean();


  backgroundFraction->SetBinContent( 1, f_BG );

  double effFrac=(double(numEvts-rejected))/double(numEvts);


  cout<<"Effective Background Fraction:           "<<f_BG<<endl;
  cout<<"Fraction of Events Surviving Efficiency: "<<effFrac<<endl;
  cout<<"Surviving Signal Events:                 "<<isBG_distribution->GetEntries()*(1-f_BG)<<endl;


  // efficiency



  leptonEfficiency->Divide( lepton_pT_vs_eta_rec, lepton_pT_vs_eta_gen, 1., 1., "B" );

  lepton_pT_vs_eta_rec->Delete();
  lepton_pT_vs_eta_gen->Delete();

  // end

  genFile->Close();

  dataFile->Write();
  dataFile->Close();

  efficiencyFile->Write();
  efficiencyFile->Close();

}

