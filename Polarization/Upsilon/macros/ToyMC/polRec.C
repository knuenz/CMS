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

double singleLeptonEfficiency( double& pT, double& eta ) {

  const double mu_pT_min = 4.0;
  const double smoothcutpar = 3.0;
  double eff = 1. / ( 1. + exp(-smoothcutpar*(pT - mu_pT_min)) );
  if ( fabs(eta) > 2.4 ) eff = 0.;
  return eff;
}


/*
double singleLeptonEfficiency( double& pT, double& eta ) {

  return 1.;
}
*/

/*
double singleLeptonEfficiency( double& pT, double& eta ) {
   double eff = 0.;
   if (  pT > 2.0 && fabs(eta) < 2.4 ) eff = 1.;
   return eff;
}
*/


void polRec(int nEff=1,
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

    double effP = singleLeptonEfficiency( lepP_pT, lepP_eta );
    double effN = singleLeptonEfficiency( lepN_pT, lepN_eta );

    double rndmeffP = gRandom->Uniform(1.);
    double rndmeffN = gRandom->Uniform(1.);


    if ( rndmeffP < effP ) lepton_pT_vs_eta_rec->Fill( fabs(lepP_eta), lepP_pT );
    if ( rndmeffN < effN ) lepton_pT_vs_eta_rec->Fill( fabs(lepN_eta), lepN_pT );


    if ( rndmeffP > effP || rndmeffN > effN ) {rejected++; continue; }

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

  delete isBG_distribution;

  backgroundFraction->SetBinContent( 1, f_BG );

  cout<<"Effective Background Fraction:           "<<f_BG<<endl;
  cout<<"Fraction of Events Surviving Efficiency: "<<(double(numEvts-rejected))/double(numEvts)<<endl;

  cout<<"Surviving Signal Events:                 "<<(numEvts-rejected)*(1-f_BG)<<endl;
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

