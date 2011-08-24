#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"

// number of random extractions for the sampling of the parameter space
// (the number of entries in the results ntuple will be smaller than this,
// due to efficiency < 1 of the sampling method and the rejected extractions of the burn-in period):
//const int n_sampledPoints = 12000;  // this INCLUDES the burn-in below.
                                    // 12000 is already perfectly good for obtaining central values and errors
                                    // more itetarions necessary for smooth 2D distributions

// number of discarded intial random extractions (burn-in period):
const int n_burnIn = 2000; // do not change this

// number of MC events for the preliminary integration of the
// parameter-independent parts of the likelihood
const double n_MCevents = 1000000;

// minimum value of the dilepton efficiency
const double min_dileptonEff = 0.01;

// Numerical inputs to calculate decay angles:
const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam_ >> Mprot_
const double Mprot_ = 0.9382720;
const double Mlepton_ = 0.10566;  // (muon)
const double gPI_ = TMath::Pi();
const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );


////////////////////////////////////////////////////////////////////////////////////////////////////
// Function randomly sampling the lambda values of the "next iteration"
// given the current values, according to a Gaussian "proposal pdf"
// including or not the positivity constraints of the angular distribution:

void extractFromProposalPDF( double& lth_candidate,     double& lph_candidate,     double& ltp_candidate,
                             double& lth,               double& lph,               double& ltp,
                             double& proposalWidth_lth, double& proposalWidth_lph, double& proposalWidth_ltp ) {

	double burnwidth=0.1;

  if ( proposalWidth_lth < 0. || proposalWidth_lph < 0. || proposalWidth_ltp < 0. ) {
       do {
           lth_candidate = gRandom->Gaus( lth, burnwidth ); // Gaussian proposal pdf with "large" sigma and positivity constraints
           lph_candidate = gRandom->Gaus( lph, burnwidth );
           ltp_candidate = gRandom->Gaus( ltp, burnwidth );
       }
       while ( TMath::Abs( lph_candidate ) > 0.5*( 1 + lth_candidate ) || lth_candidate*lth_candidate + 2.*ltp_candidate*ltp_candidate > 1
            || TMath::Abs( ltp_candidate ) > 0.5*( 1 - lph_candidate )
            || (  (1.+2.*lph_candidate)*(1.+2.*lph_candidate) + 2.*ltp_candidate*ltp_candidate > 1 && lph_candidate < -1./3. ) );
  }
  else {
      lth_candidate = gRandom->Gaus ( lth, proposalWidth_lth ); // Gaussian proposal pdf with possibly smaller sigmas
      lph_candidate = gRandom->Gaus ( lph, proposalWidth_lph ); // and no positivity constraints
      ltp_candidate = gRandom->Gaus ( ltp, proposalWidth_ltp );
  }
}

// starting values of the widths of the proposal functions to scan the parameter space.
// Negative = the burn-in period starts with "large" sigma and positivity constraints

const double proposalWidth_lth_CS_start = -1.;
const double proposalWidth_lph_CS_start = -1.;
const double proposalWidth_ltp_CS_start = -1.;

const double proposalWidth_lth_HX_start = -1.;
const double proposalWidth_lph_HX_start = -1.;
const double proposalWidth_ltp_HX_start = -1.;

const double proposalWidth_lth_PX_start = -1.;
const double proposalWidth_lph_PX_start = -1.;
const double proposalWidth_ltp_PX_start = -1.;


// calculation of the frame-invariants lambdatheta^star and lambdaphi^star:

void calcLambdastar( double& lthstar, double& lphstar,
                     double& lth,     double& lph,     double& ltp ) {

  double LamPlus = 0.25 * ( lth - lph + sqrt( pow( lth - lph, 2. ) + 4. * ltp*ltp ) );
  double LamMnus = 0.25 * ( lth - lph - sqrt( pow( lth - lph, 2. ) + 4. * ltp*ltp ) );

  lthstar = ( lth - 3.*LamPlus ) / ( 1. + LamPlus );
  lphstar = ( lph + LamPlus )    / ( 1. + LamPlus );

  double lphstarMnus = ( lph + LamMnus ) / ( 1. + LamMnus );

  if ( TMath::Abs(lphstarMnus) < TMath::Abs(lphstar) ) {
     lthstar = ( lth - 3.*LamMnus ) / ( 1. + LamMnus );
     lphstar = lphstarMnus;
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// main program /////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void polFit(int n_sampledPoints=1,
		Char_t *dirstruct = "ToyDirectory_Default"){

  gROOT->Reset();

  delete gRandom;
  gRandom = new TRandom3(0);  // better random generator


  // input files:

  // 1) file with the histogram of the single-lepton pT-eta efficiency
  char filename [500];

  sprintf(filename,"%s/efficiency.root",dirstruct);
  TFile* efficiencyFile = new TFile(filename);
  TH2D* leptonEfficiency = (TH2D*)efficiencyFile->Get("leptonEfficiency");

  // 2) file with the signal+background and background costheta-phi histograms(in the PX frame),
  //    a 2nd histogram with the background fraction,
  //    a 3rd histogram with the cell boundaries in dilepton variables
  //    and the ntuple of dilepton events

  sprintf(filename,"%s/data.root",dirstruct);
  TFile* dataFile = new TFile(filename);

  TH2D* total_PX = (TH2D*)dataFile->Get("total_PX");
  TH2D* background_PX = (TH2D*)dataFile->Get("background_PX");

  TH1D* backgroundFraction = (TH1D*)dataFile->Get("backgroundFraction");

  double  f_background = backgroundFraction->GetBinContent( 1 );


  // calculate background-subtracted ("signal") costh-phi distribution

  TH2D* signal_PX = (TH2D*)total_PX->Clone();
  signal_PX->SetName("signal_PX");

  total_PX->Sumw2(); background_PX->Sumw2(); signal_PX->Sumw2();

  total_PX->Scale( 1. / total_PX->Integral() );
  background_PX->Scale( 1. / background_PX->Integral() );

  signal_PX->Add( total_PX,                background_PX,
                  1. / (1.-f_background),  -f_background / (1.-f_background) );


  // input ntuple (dilepton events)

  TTree* data = (TTree*)dataFile->Get("data");

  TLorentzVector* lepP = 0;  data->SetBranchAddress( "lepP",  &lepP );  // lepton 4-vectors
  TLorentzVector* lepN = 0;  data->SetBranchAddress( "lepN",  &lepN );


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////    CALCULATIONS    /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Likelihood(lambdas) = Product_{i=1...n_events} PDF(lepton_momenta[i], lambdas)
  //  where PDF(lepton_momenta[i], lambdas) = probability distribution of the dilepton events
  //
  //  PDF[i] =     W_signal_theor(costheta[i], phi[i], lambdas) * epsilon(lepton_momenta[i])
  //                                / normalization
  //
  //  where epsilon[i] = product of two single-lepton efficiencies
  //
  //  W_signal_theor[i]  is proportional to   1 + lambdatheta    * costheta[i]^2
  //                                            + lambdaphi      * sintheta[i]^2 cos2phi[i]
  //                                            + lambdathetaphi * sin2theta[i] cosphi[i]
  //
  //  This is the PDF for the signal events. The background is subtracted preliminary.
  //  The PDF must be normalized so that
  //         Sum_{i=1...n_events} PDF[i] = constant independent of parameters.
  //  In general the normalization *depends* on the lambda parameters, so that it
  //  should be re-calculated (by re-scanning all events) many times during the "fit",
  //  each time that the values of the parameters change. Fortunately, in our case
  //  PDF[i] is a linear function of the lambdas, so that its normalization can
  //  be calculated as an analytical function of the lambdas, therefore *before*
  //  starting the scan of the parameter space.
  //  We preliminary calculate the parts of the PDF that do not depend on the lambdas
  //  event by event, filling arrays with one element per event:
  //
  //  w_1[i]        = epsilon[i] * 1,
  //  w_cth2[i]     = epsilon[i] * costheta[i]^2,
  //  w_sth2c2ph[i] = epsilon[i] * sintheta[i]^2 cos2phi[i],
  //  w_s2thcph[i]  = epsilon[i] * sin2theta[i] cosphi[i],       i = single *real* envent
  //
  //  To calculate the normalization of the PDF as a function of the parameters, we integrate
  //  numerically these PDF terms.
  //  To perform the integration we generate events according to the *efficiency-corrected*
  //  (not *acceptance* corrected) pT-eta-mass distribution (in the considered cell). Histograms
  //  for the three dilepton variables corrected by efficiency are prepared preliminarly while
  //  we scan over all events in the ntuple.
  //  The integration must be uniform over costh and phi in the acceptance region. This condition
  //  is realized by generating the decay distributions of the dileptons as uniform.
  //  We cannot use the data themselves (corrected by efficiency) to calculate the integrals,
  //  just because the dilepton in the data are not unpolarized.
  //
  //  sum_1        = Sum_{k} w_1[k]
  //  sum_cth2     = Sum_{k} w_cth2[k]
  //  sum_sth2c2ph = Sum_{k} w_sth2c2ph[k]
  //  sum_s2thcph  = Sum_{k} w_s2thcph[k]        k = single *generated* (unpolarized) envent
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////
  // Preliminary normalization calculations

  // calculation of the bin boundaries from the data themselves:

  double pTdilepton_min = 10000.;
  double pTdilepton_max = 0.;
  double rapdilepton_min = 10000;
  double rapdilepton_max = 0.;
  double massdilepton_min = 10000;
  double massdilepton_max = 0.;

  int n_events = int( data->GetEntries() );

  for ( int i_event = 1; i_event <= n_events; i_event++ ) {

    data->GetEvent( i_event-1 );

    TLorentzVector dilepton = *lepP + *lepN;

    // calculate minimum and maximum pT, y and mass of dilepton

    double pT   = dilepton.Pt();
    double rap  = TMath::Abs( dilepton.Rapidity() );
    double mass  = TMath::Abs( dilepton.M() );

    if ( pT > pTdilepton_max ) pTdilepton_max = pT;
    if ( pT < pTdilepton_min ) pTdilepton_min = pT;

    if ( rap > rapdilepton_max ) rapdilepton_max = rap;
    if ( rap < rapdilepton_min ) rapdilepton_min = rap;

    if ( mass > massdilepton_max ) massdilepton_max = mass;
    if ( mass < massdilepton_min ) massdilepton_min = mass;
  }


  // create histogram with |y|-pT efficiency-(not acceptance)-corrected distribution of the dilepton

  const int nbin_pT  = 10;
  const int nbin_rap = 10;
  const int nbin_mass = 10;

  TH1D* pT_effcorr = new TH1D( "pT_effcorr", "", nbin_pT, pTdilepton_min, pTdilepton_max );
  pT_effcorr->Sumw2();
  TH1D* rap_effcorr = new TH1D( "rap_effcorr", "", nbin_rap, rapdilepton_min, rapdilepton_max );
  rap_effcorr->Sumw2();
  TH1D* mass_effcorr = new TH1D( "mass_effcorr", "", nbin_mass, massdilepton_min, massdilepton_max );
  mass_effcorr->Sumw2();


  //////////////////////////////////////////////////////////////////////////////
  // create output file with ntuples of results
  // (distributions of angles and parameters)

  sprintf(filename,"%s/results.root",dirstruct);
  TFile* resultsFile = new TFile(filename, "RECREATE", "results");


  // output ntuple with angular distribution of background-subtracted events,
  // for the subsequent calculation of efficiency-corrected angular distributions


  TTree* angles = new TTree("angles","angles");

  double costh_CS;  angles->Branch("costh_CS",     &costh_CS,     "costh_CS/D");
  double phi_CS;    angles->Branch("phi_CS",       &phi_CS,       "phi_CS/D"  );
  double phith_CS;  angles->Branch("phithCS",      &phith_CS,     "phith_CS/D");

  double costh_HX;  angles->Branch("costh_HX",     &costh_HX,     "costh_HX/D");
  double phi_HX;    angles->Branch("phi_HX",       &phi_HX,       "phi_HX/D"  );
  double phith_HX;  angles->Branch("phith_HX",     &phith_HX,     "phith_HX/D");

  double costh_PX;  angles->Branch("costh_PX",     &costh_PX,     "costh_PX/D");
  double phi_PX;    angles->Branch("phi_PX",       &phi_PX,       "phi_PX/D"  );
  double phith_PX;  angles->Branch("phith_PX",     &phith_PX,     "phith_PX/D");

  double cosalpha;  angles->Branch("cosalpha",     &cosalpha,     "cosalpha/D");

  double epsilon;   angles->Branch("epsilon",      &epsilon,      "epsilon/D");


  // Read all events of the data ntuple and
  // -- subtract background
  // -- calculate the parameter-independent pieces of the event PDF and store them
  //    in large arrays, so that after this loop the the ntuple will not be accessed anymore
  // -- create output ntuple with decay angles in all frames for the background-subtracted
  //    events and the efficiency weight


  cout << endl << endl;
  cout << "Reading input ntuple (" << n_events << " dilepton events)," << endl;
  cout << "subtracting background"<< endl;
  cout << "and calculating parameter-independent PDF terms"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;


  int i_signalEvent = 0;

  int n_step = n_events/5;  // to visualize progress of the event scan (50 steps)
  int n_step_=1;

  // arrays to be filled, one element per event

  double* w_1           = new double[n_events];

  double* w_cth2_CS     = new double[n_events];
  double* w_sth2c2ph_CS = new double[n_events];
  double* w_s2thcph_CS  = new double[n_events];

  double* w_cth2_HX     = new double[n_events];
  double* w_sth2c2ph_HX = new double[n_events];
  double* w_s2thcph_HX  = new double[n_events];

  double* w_cth2_PX     = new double[n_events];
  double* w_sth2c2ph_PX = new double[n_events];
  double* w_s2thcph_PX  = new double[n_events];


  int n_background_tot = int( f_background * double ( n_events ) );
  int n_background_tmp = 0;


  // backgrond distribution TEST (goes into output file)
  // should be flat (except border fluctuations - it is a division)
  // indicating that the background distribution is well reproduced
  // (and subtracted) by the likelihood-ratio method

  TH2D* background_PX_test = (TH2D*)background_PX->Clone();
  background_PX_test->SetName("background_PX_test");
  background_PX_test->Reset();

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // loop over dilepton events in the ntuple

  for ( int i_event = 1; i_event <= n_events; i_event++ ) {


	    if (i_event%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    data->GetEvent( i_event-1 );


    // calculate efficiency ("epsilon") of the dilepton event as a product of the single-lepton
    // efficiencies, applying to them two Gaussian fluctuations (assumed to be independent)
    // using the errors in the efficiency histogram:

    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double max_pT = leptonEfficiency->GetYaxis()->GetXmax();
    if ( lepP_pT > max_pT || lepN_pT > max_pT ) continue;

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();


    double  epsilonP = leptonEfficiency->GetBinContent( leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepP_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepP_pT ) );
    double depsilonP = leptonEfficiency->GetBinError(  leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepP_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepP_pT ) );
    epsilonP = gRandom->Gaus (epsilonP, depsilonP);

    double  epsilonN = leptonEfficiency->GetBinContent( leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepN_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepN_pT ) );
    double depsilonN = leptonEfficiency->GetBinError(  leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepN_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepN_pT ) );
    epsilonN = gRandom->Gaus (epsilonN, depsilonN);

    epsilon = epsilonP * epsilonN;


    // calculation of decay angles in three polarization frames

    // dilepton 4-vector:

    TLorentzVector dilepton = *lepP + *lepN;

    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons:

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity:

    if ( dilepton.Rapidity() < 0 ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


    // positive lepton in the dilepton rest frame:

    TLorentzVector lepton_DILEP = *lepP;
    lepton_DILEP.Boost(lab_to_dilep);

    // CS frame angles:

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    TRotation rotation;
    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    costh_CS = lepton_DILEP_rotated.CosTheta();

    phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;

    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;

    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame angles:

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    costh_HX = lepton_DILEP_rotated.CosTheta();

    phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;

    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;

    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;


    // PERPENDICULAR HELICITY frame angles:

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    costh_PX = lepton_DILEP_rotated.CosTheta();

    phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;

    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;

    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


    // invariant polarization angle

    cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );


    // background subtraction: reject events having a costh_PX, phi_PH distribution
    // compatible with the corresponding background distribution.
    // Do this until a fraction f_background of the events (= n_background_tot) is rejected.

    if ( n_background_tmp < n_background_tot ) {

        int ibin_costh_PX = total_PX->GetXaxis()->FindBin( costh_PX );
        int ibin_phi_PX   = total_PX->GetYaxis()->FindBin( phi_PX );

        double likelihood_BG  = background_PX->GetBinContent( ibin_costh_PX, ibin_phi_PX );
        double likelihood_SIG =     signal_PX->GetBinContent( ibin_costh_PX, ibin_phi_PX );

        likelihood_BG  = gRandom->Uniform( likelihood_BG );
        likelihood_SIG = gRandom->Uniform( likelihood_SIG );

        if ( likelihood_BG > likelihood_SIG ) {
              ++n_background_tmp;
              background_PX_test->Fill( costh_PX, phi_PX );
              continue;
        }
    }


    if ( epsilon > min_dileptonEff ) {

      ++i_signalEvent;

      // store calculated values in the arrays

      w_1[i_signalEvent]           = epsilon;

      w_cth2_CS[i_signalEvent]     = epsilon * costh_CS*costh_CS;
      w_sth2c2ph_CS[i_signalEvent] = epsilon * (1. - costh_CS*costh_CS) * cos(gPI_/90.*phi_CS);
      w_s2thcph_CS[i_signalEvent]  = epsilon * 2.*costh_CS*sqrt(1.-costh_CS*costh_CS) * cos(gPI_/180.*phi_CS);

      w_cth2_HX[i_signalEvent]     = epsilon * costh_HX*costh_HX;
      w_sth2c2ph_HX[i_signalEvent] = epsilon * (1. - costh_HX*costh_HX) * cos(gPI_/90.*phi_HX);
      w_s2thcph_HX[i_signalEvent]  = epsilon * 2.*costh_HX*sqrt(1.-costh_HX*costh_HX) * cos(gPI_/180.*phi_HX);

      w_cth2_PX[i_signalEvent]     = epsilon * costh_PX*costh_PX;
      w_sth2c2ph_PX[i_signalEvent] = epsilon * (1. - costh_PX*costh_PX) * cos(gPI_/90.*phi_PX);
      w_s2thcph_PX[i_signalEvent]  = epsilon * 2.*costh_PX*sqrt(1.-costh_PX*costh_PX) * cos(gPI_/180.*phi_PX);


      // fill ntuple with decay angles

      angles->Fill();


      // fill pT, y and mass efficiency-corrected histograms

      pT_effcorr->Fill( dilepton.Pt(), 1./epsilon );
      rap_effcorr->Fill( TMath::Abs( dilepton.Rapidity() ), 1./epsilon );
      mass_effcorr->Fill( dilepton.M(), 1./epsilon );
    }


  } // end of loop over dilepton events in the ntuple

  cout << endl;

  double n_signalEvents = i_signalEvent;

  if ( n_background_tmp < n_background_tot ) cout << "Error: background subtraction incomplete" << endl;


  background_PX_test->Divide(background_PX);


  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////// Monte Carlo integration of the parameter-independent pieces of the PDF
  ////// *within the acceptance*, using the efficiency-corrected pT-y distribution

  n_step = n_MCevents/5;
  n_step_=1;

  cout << endl;
  cout << "MC integration of the parameter-inpependent pieces" << endl;
  cout << "(" << n_MCevents << " dilepton events)"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;


  double sum_1           = 0.;

  double sum_cth2_CS     = 0.;
  double sum_sth2c2ph_CS = 0.;
  double sum_s2thcph_CS  = 0.;

  double sum_cth2_HX     = 0.;
  double sum_sth2c2ph_HX = 0.;
  double sum_s2thcph_HX  = 0.;

  double sum_cth2_PX     = 0.;
  double sum_sth2c2ph_PX = 0.;
  double sum_s2thcph_PX  = 0.;

  int n_acceptedEvents = 0;


  // histograms with the angular dependencies of the parameter-independent PDF terms,
  // to calculate the event PDFs to be compared with the data distributions in the final plots

  int nbin_cth = 80.;
  int nbin_ph  = 72.;

  TH1D* PDF_1_vs_cth_CS       = new TH1D( "PDF_1_vs_cth_CS",       "", nbin_cth, -1., 1. );
  TH1D* PDF_cth2_vs_cth_CS    = new TH1D( "PDF_cth2_vs_cth_CS",    "", nbin_cth, -1., 1. );

  TH1D* PDF_1_vs_ph_CS        = new TH1D( "PDF_1_vs_ph_CS",        "", nbin_ph, -180., 180. );
  TH1D* PDF_c2ph_vs_ph_CS     = new TH1D( "PDF_c2ph_vs_ph_CS",     "", nbin_ph, -180., 180. );

  TH1D* PDF_1_vs_phth_CS      = new TH1D( "PDF_1_vs_phth_CS",      "", nbin_ph, -180., 180. );
  TH1D* PDF_cphth_vs_phth_CS  = new TH1D( "PDF_cphth_vs_phth_CS",  "", nbin_ph, -180., 180. );


  TH1D* PDF_1_vs_cth_HX       = new TH1D( "PDF_1_vs_cth_HX",       "", nbin_cth, -1., 1. );
  TH1D* PDF_cth2_vs_cth_HX    = new TH1D( "PDF_cth2_vs_cth_HX",    "", nbin_cth, -1., 1. );

  TH1D* PDF_1_vs_ph_HX        = new TH1D( "PDF_1_vs_ph_HX",        "", nbin_ph, -180., 180. );
  TH1D* PDF_c2ph_vs_ph_HX     = new TH1D( "PDF_c2ph_vs_ph_HX",     "", nbin_ph, -180., 180. );

  TH1D* PDF_1_vs_phth_HX      = new TH1D( "PDF_1_vs_phth_HX",      "", nbin_ph, -180., 180. );
  TH1D* PDF_cphth_vs_phth_HX  = new TH1D( "PDF_cphth_vs_phth_HX",  "", nbin_ph, -180., 180. );


  TH1D* PDF_1_vs_cth_PX       = new TH1D( "PDF_1_vs_cth_PX",       "", nbin_cth, -1., 1. );
  TH1D* PDF_cth2_vs_cth_PX    = new TH1D( "PDF_cth2_vs_cth_PX",    "", nbin_cth, -1., 1. );

  TH1D* PDF_1_vs_ph_PX        = new TH1D( "PDF_1_vs_ph_PX",        "", nbin_ph, -180., 180. );
  TH1D* PDF_c2ph_vs_ph_PX     = new TH1D( "PDF_c2ph_vs_ph_PX",     "", nbin_ph, -180., 180. );

  TH1D* PDF_1_vs_phth_PX      = new TH1D( "PDF_1_vs_phth_PX",      "", nbin_ph, -180., 180. );
  TH1D* PDF_cphth_vs_phth_PX  = new TH1D( "PDF_cphth_vs_phth_PX",  "", nbin_ph, -180., 180. );


  TH1D* PDF_1_vs_calpha       = new TH1D( "PDF_1_vs_calpha",       "", nbin_cth, -1., 1. );
  TH1D* PDF_calpha2_vs_calpha = new TH1D( "PDF_calpha2_vs_calpha", "", nbin_cth, -1., 1. );


  ///////////////// cycle of MC events ////////////////////////
  for(int i_MCevent = 1; i_MCevent <= n_MCevents; i_MCevent++){

	    if (i_MCevent%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    // generation of dilepton in the pp CM

    // pT, rapidity, mass:

    double pT = pT_effcorr->GetRandom();
    double rap = rap_effcorr->GetRandom();
    double mass = mass_effcorr->GetRandom();

    // pL:

    double rap_sign = gRandom->Uniform(-1., 1.); rap_sign /= TMath::Abs(rap_sign);
    rap *= rap_sign;
    double mT = sqrt( mass*mass + pT*pT );
    double pL1 = 0.5 *mT * exp(rap);
    double pL2 = - 0.5 *mT * exp(-rap);
    double pL = pL1 + pL2;

    // Phi:

    double Phi = 2. * gPI_ * gRandom->Uniform(1.);

    // 4-vector:

    TLorentzVector dilepton;
    dilepton.SetXYZM( pT * cos(Phi) , pT * sin(Phi), pL, mass );


    // random extraction of decay angles (generic reference frame)

    double costh_gen = -1. + 2. * gRandom->Uniform(1.);
    double sinth_gen = sqrt( 1. - costh_gen*costh_gen );
    double phi_gen = 2. * gPI_ * gRandom->Uniform(1.);

    // lepton momentum in the dilepton rest frame:

    double p_lepton_DILEP = sqrt( 0.25*mass*mass - Mlepton_*Mlepton_ );

    TLorentzVector lepton_DILEP;

    lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
                          p_lepton_DILEP * sinth_gen * sin(phi_gen),
                          p_lepton_DILEP * costh_gen,
                          Mlepton_ );


    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity

    if ( rap < 0 ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


    // step 1: transform (rotation) lepton momentum components from generation frame
    // to the frame with x,y,z axes as in the laboratory

    TVector3 oldZaxis = beam1_beam2_bisect;

    TVector3 oldYaxis = Yaxis;
    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

    TRotation rotation;
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
                     // transforms coordinates from the "old" frame to the "xyz" frame

    TLorentzVector lepton_DILEP_xyz = lepton_DILEP;

    lepton_DILEP_xyz.Transform(rotation);
                     // lepton_DILEP_xyz is the lepton in the dilepton rest frame
                     // wrt to the lab axes

    // CS frame

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame

    TVector3 lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    double costh_CS = lepton_DILEP_rotated.CosTheta();

    double phi_CS = lepton_DILEP_rotated.Phi() * 180. / gPI_;

    double phith_CS;

    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;

    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    double costh_HX = lepton_DILEP_rotated.CosTheta();

    double phi_HX = lepton_DILEP_rotated.Phi() * 180. / gPI_;

    double phith_HX;

    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;

    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;


    // PERPENDICULAR HELICITY frame

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    double costh_PX = lepton_DILEP_rotated.CosTheta();

    double phi_PX = lepton_DILEP_rotated.Phi() * 180. / gPI_;

    double phith_PX;

    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;

    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


    // invariant polarization angle

    double cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );



    // lepton 4-vectors in the LAB frame:

    TVector3 dilep_to_lab = dilepton.BoostVector();

    *lepP = lepton_DILEP_xyz;
    lepP->Boost(dilep_to_lab);
    lepN->SetPxPyPzE(-lepton_DILEP_xyz.Px(),-lepton_DILEP_xyz.Py(),-lepton_DILEP_xyz.Pz(),lepton_DILEP_xyz.E());
    lepN->Boost(dilep_to_lab);

    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double max_pT = leptonEfficiency->GetYaxis()->GetXmax();
    if ( lepP_pT > max_pT || lepN_pT > max_pT ) continue;

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();


    // efficiencies

    double  epsilonP = leptonEfficiency->GetBinContent( leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepP_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepP_pT ) );
    double depsilonP = leptonEfficiency->GetBinError(  leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepP_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepP_pT ) );
    epsilonP = gRandom->Gaus (epsilonP, depsilonP);

    double  epsilonN = leptonEfficiency->GetBinContent( leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepN_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepN_pT ) );
    double depsilonN = leptonEfficiency->GetBinError(  leptonEfficiency->GetXaxis()->FindBin( TMath::Abs( lepN_eta ) ),
                                                        leptonEfficiency->GetYaxis()->FindBin( lepN_pT ) );
    epsilonN = gRandom->Gaus (epsilonN, depsilonN);

    double epsilon = epsilonP * epsilonN;


    if ( epsilon > min_dileptonEff ) {

        //  update the sums

        ++n_acceptedEvents;

        sum_1           += epsilon;

        sum_cth2_CS     += epsilon * costh_CS*costh_CS;
        sum_sth2c2ph_CS += epsilon * (1.- costh_CS*costh_CS) * cos(gPI_/90.*phi_CS);
        sum_s2thcph_CS  += epsilon * 2.*costh_CS*sqrt(1.-costh_CS*costh_CS) * cos(gPI_/180.*phi_CS);

        sum_cth2_HX     += epsilon * costh_HX*costh_HX;
        sum_sth2c2ph_HX += epsilon * (1.- costh_HX*costh_HX) * cos(gPI_/90.*phi_HX);
        sum_s2thcph_HX  += epsilon * 2.*costh_HX*sqrt(1.-costh_HX*costh_HX) * cos(gPI_/180.*phi_HX);

        sum_cth2_PX     += epsilon * costh_PX*costh_PX;
        sum_sth2c2ph_PX += epsilon * (1.- costh_PX*costh_PX) * cos(gPI_/90.*phi_PX);
        sum_s2thcph_PX  += epsilon * 2.*costh_PX*sqrt(1.-costh_PX*costh_PX) * cos(gPI_/180.*phi_PX);

        // fill histograms for the PDF calculations

        PDF_1_vs_cth_CS->Fill( costh_CS, epsilon );
        PDF_cth2_vs_cth_CS->Fill( costh_CS, epsilon * costh_CS*costh_CS );
        PDF_1_vs_ph_CS->Fill( phi_CS, epsilon );
        PDF_c2ph_vs_ph_CS->Fill( phi_CS, epsilon * cos(gPI_/90.*phi_CS) );
        PDF_1_vs_phth_CS->Fill( phith_CS, epsilon );
        PDF_cphth_vs_phth_CS->Fill( phith_CS, epsilon * cos(gPI_/180.*phith_CS) );

        PDF_1_vs_cth_HX->Fill( costh_HX, epsilon );
        PDF_cth2_vs_cth_HX->Fill( costh_HX, epsilon * costh_HX*costh_HX );
        PDF_1_vs_ph_HX->Fill( phi_HX, epsilon );
        PDF_c2ph_vs_ph_HX->Fill( phi_HX, epsilon * cos(gPI_/90.*phi_HX) );
        PDF_1_vs_phth_HX->Fill( phith_HX, epsilon );
        PDF_cphth_vs_phth_HX->Fill( phith_HX, epsilon * cos(gPI_/180.*phith_HX) );

        PDF_1_vs_cth_PX->Fill( costh_PX, epsilon );
        PDF_cth2_vs_cth_PX->Fill( costh_PX, epsilon * costh_PX*costh_PX );
        PDF_1_vs_ph_PX->Fill( phi_PX, epsilon );
        PDF_c2ph_vs_ph_PX->Fill( phi_PX, epsilon * cos(gPI_/90.*phi_PX) );
        PDF_1_vs_phth_PX->Fill( phith_PX, epsilon );
        PDF_cphth_vs_phth_PX->Fill( phith_PX, epsilon * cos(gPI_/180.*phith_PX) );

        PDF_1_vs_calpha->Fill( cosalpha, epsilon );
        PDF_calpha2_vs_calpha->Fill( cosalpha, epsilon * cosalpha*cosalpha );

    }


  } // end loop of MC integration

  cout << endl;



  sum_1           *= n_signalEvents / n_acceptedEvents;

  sum_cth2_CS     *= n_signalEvents / n_acceptedEvents;
  sum_sth2c2ph_CS *= n_signalEvents / n_acceptedEvents;
  sum_s2thcph_CS  *= n_signalEvents / n_acceptedEvents;

  sum_cth2_HX     *= n_signalEvents / n_acceptedEvents;
  sum_sth2c2ph_HX *= n_signalEvents / n_acceptedEvents;
  sum_s2thcph_HX  *= n_signalEvents / n_acceptedEvents;

  sum_cth2_PX     *= n_signalEvents / n_acceptedEvents;
  sum_sth2c2ph_PX *= n_signalEvents / n_acceptedEvents;
  sum_s2thcph_PX  *= n_signalEvents / n_acceptedEvents;


  // input data file is no more used and is closed

  efficiencyFile->Close();
  dataFile->Close();


  // output ntuples for the CS, HX, PX (perpendicular helicity) frames
  // lth = lambdatheta, lph = lambdaphi, ltp = lambdathetaphi

  TTree* results_CS = new TTree("lambdaCS","lambdaCS");
  double lth_CS;        results_CS->Branch("lth",         &lth_CS,         "lth/D");
  double lph_CS;        results_CS->Branch("lph",         &lph_CS,         "lph/D");
  double ltp_CS;        results_CS->Branch("ltp",         &ltp_CS,         "ltp/D");
  double lthstar_CS;    results_CS->Branch("lthstar",     &lthstar_CS,     "lthstar/D");
  double lphstar_CS;    results_CS->Branch("lphstar",     &lphstar_CS,     "lphstar/D");
  double ltilde_CS;     results_CS->Branch("ltilde",      &ltilde_CS,      "ltilde/D");
  int    positivity_CS; results_CS->Branch("positivity",  &positivity_CS,  "positivity/I");

  TTree* results_HX = new TTree("lambdaHX","lambdaHX");
  double lth_HX;        results_HX->Branch("lth",         &lth_HX,         "lth/D");
  double lph_HX;        results_HX->Branch("lph",         &lph_HX,         "lph/D");
  double ltp_HX;        results_HX->Branch("ltp",         &ltp_HX,         "ltp/D");
  double lthstar_HX;    results_HX->Branch("lthstar",     &lthstar_HX,     "lthstar/D");
  double lphstar_HX;    results_HX->Branch("lphstar",     &lphstar_HX,     "lphstar/D");
  double ltilde_HX;     results_HX->Branch("ltilde",      &ltilde_HX,      "ltilde/D");
  int    positivity_HX; results_HX->Branch("positivity",  &positivity_HX,  "positivity/I");

  TTree* results_PX = new TTree("lambdaPX","lambdaPX");
  double lth_PX;        results_PX->Branch("lth",         &lth_PX,         "lth/D");
  double lph_PX;        results_PX->Branch("lph",         &lph_PX,         "lph/D");
  double ltp_PX;        results_PX->Branch("ltp",         &ltp_PX,         "ltp/D");
  double lthstar_PX;    results_PX->Branch("lthstar",     &lthstar_PX,     "lthstar/D");
  double lphstar_PX;    results_PX->Branch("lphstar",     &lphstar_PX,     "lphstar/D");
  double ltilde_PX;     results_PX->Branch("ltilde",      &ltilde_PX,      "ltilde/D");
  int    positivity_PX; results_PX->Branch("positivity",  &positivity_PX,  "positivity/I");


  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////// sampling of the parameter space ("fit") with the Markov chain

  // starting point for the "chain" of sampled values

  lth_CS = -0.5;//-0.5;
  lph_CS = 0.5;//0.5;
  ltp_CS = 0.0;

  lth_HX = 0.5;//0.5;
  lph_HX = 0.0;
  ltp_HX = 0.0;

  lth_PX = 0.5;//0.5;
  lph_PX = 0.0;
  ltp_PX = 0.0;

  double initiallikelihood = -1.e12;

  double loglikelihood_CS = initiallikelihood;  // intial (arbitrary) values
  double loglikelihood_HX = initiallikelihood;
  double loglikelihood_PX = initiallikelihood;

  // test histograms of the distributions after the burn-in period,
  // for a subsequent adjustment of the algorithm

  TH1D* test_lth_CS = new TH1D( "test_lth_CS", "", 3000, -1.5, 1.5 );
  TH1D* test_lph_CS = new TH1D( "test_lph_CS", "", 3000, -1.5, 1.5 );
  TH1D* test_ltp_CS = new TH1D( "test_ltp_CS", "", 3000, -1.0, 1.0 );

  TH1D* test_lth_HX = new TH1D( "test_lth_HX", "", 3000, -1.5, 1.5 );
  TH1D* test_lph_HX = new TH1D( "test_lph_HX", "", 3000, -1.5, 1.5 );
  TH1D* test_ltp_HX = new TH1D( "test_ltp_HX", "", 3000, -1.0, 1.0 );

  TH1D* test_lth_PX = new TH1D( "test_lth_PX", "", 3000, -1.5, 1.5 );
  TH1D* test_lph_PX = new TH1D( "test_lph_PX", "", 3000, -1.5, 1.5 );
  TH1D* test_ltp_PX = new TH1D( "test_ltp_PX", "", 3000, -1.0, 1.0 );


  cout << endl;
  cout << "Sampling of the parameter space (" << n_burnIn << "+"
                         << n_sampledPoints - n_burnIn << " iterations)" << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;

  // Initial values of the sigmas of the proposal functions.
  // Speed and efficiency of the algorithm are extremely sensistive to the
  // sigma of the proposal function. The intial sigma should be large enough, because
  // a large sigma helps finding qickly where the distribution is centred. However, if the
  // sigma is too large, the efficiency of the method is very small, i.e. many events are
  // rejected and it is necessary to increase a lot the number of extractions.
  // For these reasons, during the burn-in period a large sigma is used for each parameter
  // and after this period the sigma is adjusted looking at the distribution produced so far.
  // Optimal: sigma(proposal) ~0.3 * sigma(target)


  double proposalWidth_lth_CS = proposalWidth_lth_CS_start;
  double proposalWidth_lph_CS = proposalWidth_lph_CS_start;
  double proposalWidth_ltp_CS = proposalWidth_ltp_CS_start;

  double proposalWidth_lth_HX = proposalWidth_lth_HX_start;
  double proposalWidth_lph_HX = proposalWidth_lph_HX_start;
  double proposalWidth_ltp_HX = proposalWidth_ltp_HX_start;

  double proposalWidth_lth_PX = proposalWidth_lth_PX_start;
  double proposalWidth_lph_PX = proposalWidth_lph_PX_start;
  double proposalWidth_ltp_PX = proposalWidth_ltp_PX_start;

  int rejectPointCS=0;
  int rejectPointHX=0;
  int rejectPointPX=0;

  n_step = n_sampledPoints/5;  // visualize progress of the parameter sampling
  n_step_=1;

  for(int i_sampledPoint = 1; i_sampledPoint <= n_sampledPoints; i_sampledPoint++){

	    if (i_sampledPoint%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

     // initialize positivity condition to positive

     positivity_CS = 1;
     positivity_HX = 1;
     positivity_PX = 1;

     // random extraction of parameter values according to the Metropolis-Hastings
     // implementation of the Markov chain

     double lth_CS_candidate, lph_CS_candidate, ltp_CS_candidate;
     double lth_HX_candidate, lph_HX_candidate, ltp_HX_candidate;
     double lth_PX_candidate, lph_PX_candidate, ltp_PX_candidate;

     // sampling of the parameter space using Gaussian "proposal PDF"
     // centred in the "previous" parameter values and generating
     // "candidate" values for the next iteration

     extractFromProposalPDF( lth_CS_candidate,     lph_CS_candidate,     ltp_CS_candidate,
                             lth_CS,               lph_CS,               ltp_CS,
                             proposalWidth_lth_CS, proposalWidth_lph_CS, proposalWidth_ltp_CS );
     extractFromProposalPDF( lth_HX_candidate,     lph_HX_candidate,     ltp_HX_candidate,
                             lth_HX,               lph_HX,               ltp_HX,
                             proposalWidth_lth_HX, proposalWidth_lph_HX, proposalWidth_ltp_HX );
     extractFromProposalPDF( lth_PX_candidate,     lph_PX_candidate,     ltp_PX_candidate,
                             lth_PX,               lph_PX,               ltp_PX,
                             proposalWidth_lth_PX, proposalWidth_lph_PX, proposalWidth_ltp_PX );


     // calculation of the likelihoods to decide to keep or reject the candidates

     double loglikelihood_CS_candidate = 0; // calculate log(likelihood) to avoid that
     double loglikelihood_HX_candidate = 0; // the likelihood goes to zero in the product over all events
     double loglikelihood_PX_candidate = 0; // because of limited numerical precision

     const double min_likelihood = 1.e-100;


     // loop over arrays of selected dilepton events to calculate the likelihoods
     /////////////////////////////////////////////////////////////////////////////////////
     for ( int i = 0; i < n_signalEvents; i++ ) {

       double singleEventLikelihood_CS = (  w_1[i] + lth_CS_candidate * w_cth2_CS[i]
                                                   + lph_CS_candidate * w_sth2c2ph_CS[i]
                                                   + ltp_CS_candidate * w_s2thcph_CS[i]  )
                                      /  (  sum_1  + lth_CS_candidate * sum_cth2_CS
                                                   + lph_CS_candidate * sum_sth2c2ph_CS
                                                   + ltp_CS_candidate * sum_s2thcph_CS   ) ;
       if ( singleEventLikelihood_CS < min_likelihood ) singleEventLikelihood_CS = min_likelihood;
       loglikelihood_CS_candidate += log( singleEventLikelihood_CS );

       double singleEventLikelihood_HX = (  w_1[i] + lth_HX_candidate * w_cth2_HX[i]
                                                   + lph_HX_candidate * w_sth2c2ph_HX[i]
                                                   + ltp_HX_candidate * w_s2thcph_HX[i]  )
                                      /  (  sum_1  + lth_HX_candidate * sum_cth2_HX
                                                   + lph_HX_candidate * sum_sth2c2ph_HX
                                                   + ltp_HX_candidate * sum_s2thcph_HX   ) ;
       if ( singleEventLikelihood_HX < min_likelihood ) singleEventLikelihood_HX = min_likelihood;
       loglikelihood_HX_candidate += log( singleEventLikelihood_HX );

       double singleEventLikelihood_PX = (  w_1[i] + lth_PX_candidate * w_cth2_PX[i]
                                                   + lph_PX_candidate * w_sth2c2ph_PX[i]
                                                   + ltp_PX_candidate * w_s2thcph_PX[i]  )
                                      /  (  sum_1  + lth_PX_candidate * sum_cth2_PX
                                                   + lph_PX_candidate * sum_sth2c2ph_PX
                                                   + ltp_PX_candidate * sum_s2thcph_PX   ) ;
       if ( singleEventLikelihood_PX < min_likelihood ) singleEventLikelihood_PX = min_likelihood;
       loglikelihood_PX_candidate += log( singleEventLikelihood_PX );

     } // end of loop over arrays of selected dilepton events


     // apply Metropolis-Hastings algorithm: the candidate parameter values are
     // kept or rejected depending on the likelihood ratio wrt the "previous" values.
     // The likelihood ratio, if smaller than 1, represents the "probability" with which
     // the new values must be taken as next "good" values of the "chain" and as starting
     // point for the next iteration.
     // If the ratio is > 1 the new values are always taken as good.
     // The condition "likelihood ratio > or < 1" is translated into the condition
     // "log(likelihood) difference > or < 0"

     double loglikelihood_CS_difference = loglikelihood_CS_candidate - loglikelihood_CS;

//     cout<<i_sampledPoint<<"delta "<<loglikelihood_CS_difference<<" = "<<loglikelihood_CS_candidate<<" - "<<loglikelihood_CS<<endl;

     if(  loglikelihood_CS_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_CS_difference  ) {
         lth_CS = lth_CS_candidate; lph_CS = lph_CS_candidate; ltp_CS = ltp_CS_candidate;
         loglikelihood_CS = loglikelihood_CS_candidate;
         if ( i_sampledPoint > n_burnIn ) {  // reject first n_burnIn extractions in the chain

            if ( TMath::Abs( lph_CS ) > 0.5*( 1 + lth_CS ) || lth_CS*lth_CS + 2.*ltp_CS*ltp_CS > 1
            || TMath::Abs( ltp_CS ) > 0.5*( 1 - lph_CS )
            || (  (1.+2.*lph_CS)*(1.+2.*lph_CS) + 2.*ltp_CS*ltp_CS > 1 && lph_CS < -1./3. ) )
            positivity_CS = 0; // apply positivity constraint in a flag-variable of the ntuple

            calcLambdastar( lthstar_CS, lphstar_CS, lth_CS, lph_CS, ltp_CS );
            ltilde_CS = (lth_CS + 3.*lph_CS)/(1.-lph_CS);
            results_CS->Fill();
//            cout<<lth_CS_candidate<<endl;
         }
         else if ( i_sampledPoint > n_burnIn/2 ) { test_lth_CS->Fill( lth_CS );
                                                   test_lph_CS->Fill( lph_CS );
                                                   test_ltp_CS->Fill( ltp_CS ); }
     }
     else {rejectPointCS++; }
     if ( i_sampledPoint == n_burnIn +1 ) { proposalWidth_lth_CS = TMath::Max( 0.3 * test_lth_CS->GetRMS(), 0.001 );
                                            proposalWidth_lph_CS = TMath::Max( 0.3 * test_lph_CS->GetRMS(), 0.001 );
                                            proposalWidth_ltp_CS = TMath::Max( 0.3 * test_ltp_CS->GetRMS(), 0.001 ); }
     // fill test histograms in the second half of the burn-in period (when the algorithm
     // should already have found where the bulk of the distribution is) in order to estimate
     // the sigmas of the output distributions and to adjust the sigmas of the proposal functions
     // accordingly (optimal: between 1/4 and 1/3 of the target distributions)

     double loglikelihood_HX_difference = loglikelihood_HX_candidate - loglikelihood_HX;
     if(  loglikelihood_HX_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_HX_difference  ) {
         lth_HX = lth_HX_candidate; lph_HX = lph_HX_candidate; ltp_HX = ltp_HX_candidate;
         loglikelihood_HX = loglikelihood_HX_candidate;
         if ( i_sampledPoint > n_burnIn ) {

            if ( TMath::Abs( lph_HX ) > 0.5*( 1 + lth_HX ) || lth_HX*lth_HX + 2.*ltp_HX*ltp_HX > 1
            || TMath::Abs( ltp_HX ) > 0.5*( 1 - lph_HX )
            || (  (1.+2.*lph_HX)*(1.+2.*lph_HX) + 2.*ltp_HX*ltp_HX > 1 && lph_HX < -1./3. ) )
            positivity_HX = 0;

            calcLambdastar( lthstar_HX, lphstar_HX, lth_HX, lph_HX, ltp_HX );
            ltilde_HX = (lth_HX + 3.*lph_HX)/(1.-lph_HX);
            results_HX->Fill();
//            if(lph_HX>0) cout<<"lph_HX "<<lph_HX<<", lth_HX "<<lth_HX<<", ltp_HX "<<ltp_HX<<endl;
         }
         else if ( i_sampledPoint > n_burnIn/2 ) { test_lth_HX->Fill( lth_HX );
                                                   test_lph_HX->Fill( lph_HX );
                                                   test_ltp_HX->Fill( ltp_HX ); }
     }
     else {rejectPointHX++; }
     if ( i_sampledPoint == n_burnIn +1 ) { proposalWidth_lth_HX = TMath::Max( 0.3 * test_lth_HX->GetRMS(), 0.001 );
                                            proposalWidth_lph_HX = TMath::Max( 0.3 * test_lph_HX->GetRMS(), 0.001 );
                                            proposalWidth_ltp_HX = TMath::Max( 0.3 * test_ltp_HX->GetRMS(), 0.001 ); }


     double loglikelihood_PX_difference = loglikelihood_PX_candidate - loglikelihood_PX;
     if(  loglikelihood_PX_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_PX_difference  ) {
         lth_PX = lth_PX_candidate; lph_PX = lph_PX_candidate; ltp_PX = ltp_PX_candidate;
         loglikelihood_PX = loglikelihood_PX_candidate;
         if ( i_sampledPoint > n_burnIn ) {

            if ( TMath::Abs( lph_PX ) > 0.5*( 1 + lth_PX ) || lth_PX*lth_PX + 2.*ltp_PX*ltp_PX > 1
            || TMath::Abs( ltp_PX ) > 0.5*( 1 - lph_PX )
            || (  (1.+2.*lph_PX)*(1.+2.*lph_PX) + 2.*ltp_PX*ltp_PX > 1 && lph_PX < -1./3. ) )
            positivity_PX = 0;

            calcLambdastar( lthstar_PX, lphstar_PX, lth_PX, lph_PX, ltp_PX );
            ltilde_PX = (lth_PX + 3.*lph_PX)/(1.-lph_PX);
            results_PX->Fill();
         }
         else if ( i_sampledPoint > n_burnIn/2 ) { test_lth_PX->Fill( lth_PX );
                                                   test_lph_PX->Fill( lph_PX );
                                                   test_ltp_PX->Fill( ltp_PX ); }
     }
     else {rejectPointPX++; }
     if ( i_sampledPoint == n_burnIn +1 ) { proposalWidth_lth_PX = TMath::Max( 0.3 * test_lth_PX->GetRMS(), 0.001 );
                                            proposalWidth_lph_PX = TMath::Max( 0.3 * test_lph_PX->GetRMS(), 0.001 );
                                            proposalWidth_ltp_PX = TMath::Max( 0.3 * test_ltp_PX->GetRMS(), 0.001 ); }




  } // end of parameter sampling loop

  cout << endl << endl;

  cout<<"Fraction of rejected Points CS: "<<double(rejectPointCS)/double(n_sampledPoints)<<endl;
  cout<<"Fraction of rejected Points HX: "<<double(rejectPointHX)/double(n_sampledPoints)<<endl;
  cout<<"Fraction of rejected Points PX: "<<double(rejectPointPX)/double(n_sampledPoints)<<endl;

///// Extract numerical values of the lambda parameters from the result TTrees /////

  // extremes and binning of lambda plots
  const double l_min = -1.4;
  const double l_max =  1.4;
  const double l_step_1D = 0.02;


  TH1D* h_lth_CS = new TH1D( "h_lth_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_lph_CS = new TH1D( "h_lph_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_ltp_CS = new TH1D( "h_ltp_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_ltilde_CS = new TH1D( "h_ltilde_CS", "", int((l_max-l_min)/l_step_1D), l_min, 3*l_max );
  TH1D* h_lphstar_CS = new TH1D( "h_lphstar_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_lthstar_CS = new TH1D( "h_lthstar_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

  TH1D* h_lth_HX = new TH1D( "h_lth_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_lph_HX = new TH1D( "h_lph_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_ltp_HX = new TH1D( "h_ltp_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_ltilde_HX = new TH1D( "h_ltilde_HX", "", int((l_max-l_min)/l_step_1D), l_min, 3*l_max );
  TH1D* h_lphstar_HX = new TH1D( "h_lphstar_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_lthstar_HX = new TH1D( "h_lthstar_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

  TH1D* h_lth_PX = new TH1D( "h_lth_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_lph_PX = new TH1D( "h_lph_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_ltp_PX = new TH1D( "h_ltp_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_ltilde_PX = new TH1D( "h_ltilde_PX", "", int((l_max-l_min)/l_step_1D), l_min, 3*l_max );
  TH1D* h_lphstar_PX = new TH1D( "h_lphstar_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
  TH1D* h_lthstar_PX = new TH1D( "h_lthstar_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

  // loop over entries in the ntuples

  int n_entries = int( results_CS->GetEntries() );
  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
    results_CS->GetEvent( i_entry );
    h_lth_CS->Fill( lth_CS );
    h_lph_CS->Fill( lph_CS );
    h_ltp_CS->Fill( ltp_CS );
    h_ltilde_CS->Fill( ltilde_CS );
    h_lphstar_CS->Fill( lphstar_CS );
    h_lthstar_CS->Fill( lthstar_CS );
  }

  n_entries = int( results_HX->GetEntries() );
  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
    results_HX->GetEvent( i_entry );
    h_lth_HX->Fill( lth_HX );
    h_lph_HX->Fill( lph_HX );
    h_ltp_HX->Fill( ltp_HX );
    h_ltilde_HX->Fill( ltilde_HX );
    h_lphstar_HX->Fill( lphstar_HX );
    h_lthstar_HX->Fill( lthstar_HX );
  }

  n_entries = int( results_PX->GetEntries() );
  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
    results_PX->GetEvent( i_entry );
    h_lth_PX->Fill( lth_PX );
    h_lph_PX->Fill( lph_PX );
    h_ltp_PX->Fill( ltp_PX );
    h_ltilde_PX->Fill( ltilde_PX );
    h_lphstar_PX->Fill( lphstar_PX );
    h_lthstar_PX->Fill( lthstar_PX );
  }

  // end of loop over entries in the ntuples

//  cout << endl << endl;


  TTree* Results = new TTree("Results","Results");
  Results->Branch("lthCS",         &lth_CS,         "lthCS/D");
  Results->Branch("lphCS",         &lph_CS,         "lphCS/D");
  Results->Branch("ltpCS",         &ltp_CS,         "ltpCS/D");
  Results->Branch("lthstarCS",     &lthstar_CS,     "lthstarCS/D");
  Results->Branch("lphstarCS",     &lphstar_CS,     "lphstarCS/D");
  Results->Branch("ltildeCS",      &ltilde_CS,      "ltildeCS/D");
  Results->Branch("positivityCS",  &positivity_CS,  "positivityCS/I");
  Results->Branch("lthHX",         &lth_HX,         "lthHX/D");
  Results->Branch("lphHX",         &lph_HX,         "lphHX/D");
  Results->Branch("ltpHX",         &ltp_HX,         "ltpHX/D");
  Results->Branch("lthstarHX",     &lthstar_HX,     "lthstarHX/D");
  Results->Branch("lphstarHX",     &lphstar_HX,     "lphstarHX/D");
  Results->Branch("ltildeHX",      &ltilde_HX,      "ltildeHX/D");
  Results->Branch("positivityHX",  &positivity_HX,  "positivityHX/I");
  Results->Branch("lthPX",         &lth_PX,         "lthPX/D");
  Results->Branch("lphPX",         &lph_PX,         "lphPX/D");
  Results->Branch("ltpPX",         &ltp_PX,         "ltpPX/D");
  Results->Branch("lthstarPX",     &lthstar_PX,     "lthstarPX/D");
  Results->Branch("lphstarPX",     &lphstar_PX,     "lphstarPX/D");
  Results->Branch("ltildePX",      &ltilde_PX,      "ltildePX/D");
  Results->Branch("positivityPX",  &positivity_PX,  "positivityPX/I");

  double err_lth_CS;
  double err_lph_CS;
  double err_ltp_CS;
  double err_lthstar_CS;
  double err_lphstar_CS;
  double err_ltilde_CS;
  double err_lth_HX;
  double err_lph_HX;
  double err_ltp_HX;
  double err_lthstar_HX;
  double err_lphstar_HX;
  double err_ltilde_HX;
  double err_lth_PX;
  double err_lph_PX;
  double err_ltp_PX;
  double err_lthstar_PX;
  double err_lphstar_PX;
  double err_ltilde_PX;

  Results->Branch("err_lthCS",         &err_lth_CS,         "err_lthCS/D");
  Results->Branch("err_lphCS",         &err_lph_CS,         "err_lphCS/D");
  Results->Branch("err_ltpCS",         &err_ltp_CS,         "err_ltpCS/D");
  Results->Branch("err_lthstarCS",     &err_lthstar_CS,     "err_lthstarCS/D");
  Results->Branch("err_lphstarCS",     &err_lphstar_CS,     "err_lphstarCS/D");
  Results->Branch("err_ltildeCS",      &err_ltilde_CS,      "err_ltildeCS/D");
  Results->Branch("err_lthHX",         &err_lth_HX,         "err_lthHX/D");
  Results->Branch("err_lphHX",         &err_lph_HX,         "err_lphHX/D");
  Results->Branch("err_ltpHX",         &err_ltp_HX,         "err_ltpHX/D");
  Results->Branch("err_lthstarHX",     &err_lthstar_HX,     "err_lthstarHX/D");
  Results->Branch("err_lphstarHX",     &err_lphstar_HX,     "err_lphstarHX/D");
  Results->Branch("err_ltildeHX",      &err_ltilde_HX,      "err_ltildeHX/D");
  Results->Branch("err_lthPX",         &err_lth_PX,         "err_lthPX/D");
  Results->Branch("err_lphPX",         &err_lph_PX,         "err_lphPX/D");
  Results->Branch("err_ltpPX",         &err_ltp_PX,         "err_ltpPX/D");
  Results->Branch("err_lthstarPX",     &err_lthstar_PX,     "err_lthstarPX/D");
  Results->Branch("err_lphstarPX",     &err_lphstar_PX,     "err_lphstarPX/D");
  Results->Branch("err_ltildePX",      &err_ltilde_PX,      "err_ltildePX/D");

  lth_CS=h_lth_CS->GetMean();
  err_lth_CS=h_lth_CS->GetRMS();
  lph_CS=h_lph_CS->GetMean();
  err_lph_CS=h_lph_CS->GetRMS();
  ltp_CS=h_ltp_CS->GetMean();
  err_ltp_CS=h_ltp_CS->GetRMS();
  lthstar_CS=h_lthstar_CS->GetMean();
  err_lthstar_CS=h_lthstar_CS->GetRMS();
  lphstar_CS=h_lphstar_CS->GetMean();
  err_lphstar_CS=h_lphstar_CS->GetRMS();
  ltilde_CS=h_ltilde_CS->GetMean();
  err_ltilde_CS=h_ltilde_CS->GetRMS();

  lth_HX=h_lth_HX->GetMean();
  err_lth_HX=h_lth_HX->GetRMS();
  lph_HX=h_lph_HX->GetMean();
  err_lph_HX=h_lph_HX->GetRMS();
  ltp_HX=h_ltp_HX->GetMean();
  err_ltp_HX=h_ltp_HX->GetRMS();
  lthstar_HX=h_lthstar_HX->GetMean();
  err_lthstar_HX=h_lthstar_HX->GetRMS();
  lphstar_HX=h_lphstar_HX->GetMean();
  err_lphstar_HX=h_lphstar_HX->GetRMS();
  ltilde_HX=h_ltilde_HX->GetMean();
  err_ltilde_HX=h_ltilde_HX->GetRMS();

  lth_PX=h_lth_PX->GetMean();
  err_lth_PX=h_lth_PX->GetRMS();
  lph_PX=h_lph_PX->GetMean();
  err_lph_PX=h_lph_PX->GetRMS();
  ltp_PX=h_ltp_PX->GetMean();
  err_ltp_PX=h_ltp_PX->GetRMS();
  lthstar_PX=h_lthstar_PX->GetMean();
  err_lthstar_PX=h_lthstar_PX->GetRMS();
  lphstar_PX=h_lphstar_PX->GetMean();
  err_lphstar_PX=h_lphstar_PX->GetRMS();
  ltilde_PX=h_ltilde_PX->GetMean();
  err_ltilde_PX=h_ltilde_PX->GetRMS();


  Results->Fill();


  test_lth_CS->Delete();
  test_lph_CS->Delete();
  test_ltp_CS->Delete();

  test_lth_HX->Delete();
  test_lph_HX->Delete();
  test_ltp_HX->Delete();

  test_lth_PX->Delete();
  test_lph_PX->Delete();
  test_ltp_PX->Delete();

  resultsFile->Write();

} // end of main
