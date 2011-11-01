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
#include "TH3.h"
#include "TFile.h"

//#include "genDefs.h"
//#include "effsAndCuts.h"

bool isMuonInAcceptance(int iCut, double pT, double eta);
double singleLeptonEfficiency( double& pT, double& eta, int nEff, TH1* hEff);
void EvaluateEffFileName(int nEff, char EffFileName [200], bool singleLeptonEff);
double DiLeptonEfficiency( double& Dilepton_pT, double& Dilepton_rap, int nDileptonEff, TH1* hDileptonEff);

void polRec(double rapdilepton_min = 1,
		double rapdilepton_max = 1,
		double pTdilepton_min = 1,
		double pTdilepton_max = 1,
		double mass_signal_peak  =  1,
		double mass_signal_sigma =  1,
		double n_sigmas_signal = 1,
		int nEff=1,
		int nDileptonEff=1,
		int FidCuts=0,
		Char_t *dirstruct = "ToyDirectory_Default",
		bool applyFidCuts=false,
		Char_t *effDir = "effDir_Default"){

  gROOT->Reset();


  //Get single Lepton Efficiency File name
    char EffFile[200];
    char EffFileName[200];
    EvaluateEffFileName(nEff,EffFileName,true);
    sprintf(EffFile,"%s/%s",effDir,EffFileName);

    TFile *fInEff = new TFile(EffFile);
    TH1* hEff=(TH1*) fInEff->Get("hEff_DATA_central");

  //Get DiLepton Efficiency File name

    EvaluateEffFileName(nDileptonEff,EffFileName,false);
    sprintf(EffFile,"%s/%s",effDir,EffFileName);

    TFile *fInDileptonEff = new TFile(EffFile);
    TH1* hDileptonEff=(TH1*) fInDileptonEff->Get("hEff_DATA_central");



	double mass_min = mass_signal_peak - n_sigmas_signal*mass_signal_sigma;
	double mass_max = mass_signal_peak + n_sigmas_signal*mass_signal_sigma;

  char filename [500];

  // get ntuple of generated data from file

  sprintf(filename,"%s/genData.root",dirstruct);
  TFile* genFile = new TFile(filename,"READ");
  TTree* genData = (TTree*)genFile->Get("genData");


  // create output data file

  sprintf(filename,"%s/data.root",dirstruct);
  TFile* dataFile = new TFile(filename, "RECREATE", "dataFile");


  // output ntuple

  TTree* data = new TTree("selectedData","selectedData");


  // histograms to be output to the opened file

  // background distributions

  const int nbinth   = 10;
  const int nbinph   = 10;
  const int nbinpT   =  7;
  const int nbinrap  =  2;
  const int nbinmass =  7;

  TH2D* background_costhphiPX = new TH2D( "background_costhphiPHX", "", nbinth,    -1.,    1.,
                                                                        nbinph,   180.,  180.  );
  TH3D* background_pTrapMass  = new TH3D( "background_pTrapMass",   "", nbinpT,   pTdilepton_min,  pTdilepton_max,
                                                                        nbinrap,  rapdilepton_min, rapdilepton_max,
                                                                        nbinmass, mass_min,        mass_max         );
  TH1D* background_fraction = new TH1D( "background_fraction", "", 1, 0., 1. );

  // this is a temporary histogram to calculate BG fraction after acceptence and efficiency cuts
  TH1D* isBGdistribution    = new TH1D( "isBGdistribution", "", 2, 0., 2. );


  // structure of input ntuple

  TLorentzVector* lepP_gen = 0;    genData->SetBranchAddress( "lepP",    &lepP_gen );
  TLorentzVector* lepN_gen = 0;    genData->SetBranchAddress( "lepN",    &lepN_gen );

  double costh_CS;  genData->SetBranchAddress( "costh_CS",     &costh_CS );
  double phi_CS;    genData->SetBranchAddress( "phi_CS",       &phi_CS   );
  double phith_CS;  genData->SetBranchAddress( "phith_CS",     &phith_CS );

  double costh_HX;  genData->SetBranchAddress( "costh_HX",     &costh_HX );
  double phi_HX;    genData->SetBranchAddress( "phi_HX",       &phi_HX   );
  double phith_HX;  genData->SetBranchAddress( "phith_HX",     &phith_HX );

  double costh_PX;  genData->SetBranchAddress( "costh_PX",     &costh_PX );
  double phi_PX;    genData->SetBranchAddress( "phi_PX",       &phi_PX   );
  double phith_PX;  genData->SetBranchAddress( "phith_PX",     &phith_PX );

  double mass;      genData->SetBranchAddress( "mass",         &mass     );
  double pT;        genData->SetBranchAddress( "pT",           &pT       );
  double rap;       genData->SetBranchAddress( "rap",          &rap      );

  int isBG;         genData->SetBranchAddress( "isBG",         &isBG     );


  // structure of output ntuple

  TLorentzVector* lepP = new TLorentzVector(0.,0.,0.,0.);  data->Branch( "lepP", "TLorentzVector", &lepP );
  TLorentzVector* lepN = new TLorentzVector(0.,0.,0.,0.);  data->Branch( "lepN", "TLorentzVector", &lepN );


  // loop over events in the input ntuple

  int numEvts = int( genData->GetEntries() );


  cout << endl;
  cout << "Reading " << numEvts << " dilepton events"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;

  int n_step = numEvts/5;
  int n_step_=1;
  int rejected=0;

  for ( int evtNumber = 0; evtNumber < numEvts; evtNumber++ ) {

	    if ((evtNumber+1)%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    genData->GetEvent( evtNumber );




    // select data in acceptance and apply efficiency

    double lepP_pT  = lepP_gen->Pt();
    double lepN_pT  = lepN_gen->Pt();
    double lepP_eta  = lepP_gen->PseudoRapidity();
    double lepN_eta  = lepN_gen->PseudoRapidity();

    bool isEventAccepted = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
                         * isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );

    if ( !isEventAccepted ) {rejected++; continue;}

    double effP = singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, hEff );
    double effN = singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, hEff );
    double DileptonEff = DiLeptonEfficiency( pT, rap, nDileptonEff, hDileptonEff );

    double rndmeffP = gRandom->Uniform(1.);
    double rndmeffN = gRandom->Uniform(1.);
    double rndmDileptoneff = gRandom->Uniform(1.);

    if ( rndmeffP > effP || rndmeffN > effN || rndmDileptoneff > DileptonEff ) {rejected++; continue;}

    // fill background histograms and output ntuple

    isBGdistribution->Fill( isBG );
    if ( isBG ) {
       background_costhphiPX->Fill( costh_PX, phi_PX );
       background_pTrapMass->Fill( pT, TMath::Abs(rap), mass );
    }

    *lepP = *lepP_gen;
    *lepN = *lepN_gen;

    data->Fill();

  } // end of RD ntuple loop

  cout << endl << endl;



  // background fraction

  double f_BG = isBGdistribution->GetMean();


  background_fraction->SetBinContent( 1, f_BG );

  double effFrac=(double(numEvts-rejected))/double(numEvts);


  cout<<"Effective Background Fraction:           "<<f_BG<<endl;
  cout<<"Fraction of Events Surviving Efficiency: "<<effFrac<<endl;
  cout<<"Surviving Signal Events:                 "<<isBGdistribution->GetEntries()*(1-f_BG)<<endl;

  // end

  genFile->Close();
  dataFile->Write();
 dataFile->Close();

}

