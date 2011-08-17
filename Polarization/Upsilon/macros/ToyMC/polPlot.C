#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TPostscript.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLatex.h"

// binning of costh - phi plots
const int nbin_cth = 80.;
const int nbin_ph  = 72.;

// extremes and binning of lambda plots
const double lth_min = -1.4;
const double lth_max =  1.4;
const double lth_step_1D = 0.02;
const double lth_step_2D = 0.02;

const double lph_min = -1.4;
const double lph_max =  1.4;
const double lph_step_1D = 0.02;
const double lph_step_2D = 0.02;

const double ltp_min = -1.4;
const double ltp_max =  1.4;
const double ltp_step_1D = 0.02;
const double ltp_step_2D = 0.02;


// function to calculate height of the contour of a 2D distribution for a certain confidence level
inline double contourHeight2D ( TH2D *h, double confidenceLevel ) {
  int Nx = h->GetXaxis()->GetNbins();
  int Ny = h->GetYaxis()->GetNbins();

  double totSum = h->GetSum();
  double targetSum = confidenceLevel * totSum;
  double maxHeight = h->GetMaximum();
  double step = 0.001*maxHeight;

  double tempHeight = 0.;
  double tempSum = totSum;

  while ( tempSum > targetSum && tempHeight < maxHeight ) {
    tempHeight += step;
    tempSum = 0.;
    for ( int ix = 0; ix < Nx; ix++ ) {
      for ( int iy = 0; iy < Ny; iy++ ) {
      	double binContent = h->GetBinContent(ix,iy);
        if ( binContent > tempHeight ) tempSum += binContent;
      }
    }
  }
  return tempHeight;
}

// function to set the 99% and 68% C.L. contours of a 2D histogram
inline void setContourHistogram ( TH2D *h ) {
  double cont0 = contourHeight2D( h, 0.99 );
  double cont1 = contourHeight2D( h, 0.68 );
  h->SetContour(2);
  h->SetContourLevel(0,cont0);
  h->SetContourLevel(1,cont1);
}



void polPlot(Char_t *dirstruct = "ToyDirectory_Default"){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  char filename [500];

  // input file:

  sprintf(filename,"%s/results.root",dirstruct);
  TFile* results = new TFile(filename);

  int n_entries, n_step;

  // input ntuples with PDFs of the anisotropy parameters

  TTree* lambdaCS = (TTree*)results->Get("lambdaCS");

  double lth_CS;        lambdaCS->SetBranchAddress("lth",         &lth_CS        );
  double lph_CS;        lambdaCS->SetBranchAddress("lph",         &lph_CS        );
  double ltp_CS;        lambdaCS->SetBranchAddress("ltp",         &ltp_CS        );
  double lthstar_CS;    lambdaCS->SetBranchAddress("lthstar",     &lthstar_CS    );
  double lphstar_CS;    lambdaCS->SetBranchAddress("lphstar",     &lphstar_CS    );
  double ltilde_CS;     lambdaCS->SetBranchAddress("ltilde",      &ltilde_CS     );
  int    positivity_CS; lambdaCS->SetBranchAddress("positivity",  &positivity_CS );

  TTree* lambdaHX = (TTree*)results->Get("lambdaHX");

  double lth_HX;        lambdaHX->SetBranchAddress("lth",         &lth_HX        );
  double lph_HX;        lambdaHX->SetBranchAddress("lph",         &lph_HX        );
  double ltp_HX;        lambdaHX->SetBranchAddress("ltp",         &ltp_HX        );
  double lthstar_HX;    lambdaHX->SetBranchAddress("lthstar",     &lthstar_HX    );
  double lphstar_HX;    lambdaHX->SetBranchAddress("lphstar",     &lphstar_HX    );
  double ltilde_HX;     lambdaHX->SetBranchAddress("ltilde",      &ltilde_HX     );
  int    positivity_HX; lambdaHX->SetBranchAddress("positivity",  &positivity_HX );

  TTree* lambdaPX = (TTree*)results->Get("lambdaPX");

  double lth_PX;        lambdaPX->SetBranchAddress("lth",         &lth_PX        );
  double lph_PX;        lambdaPX->SetBranchAddress("lph",         &lph_PX        );
  double ltp_PX;        lambdaPX->SetBranchAddress("ltp",         &ltp_PX        );
  double lthstar_PX;    lambdaPX->SetBranchAddress("lthstar",     &lthstar_PX    );
  double lphstar_PX;    lambdaPX->SetBranchAddress("lphstar",     &lphstar_PX    );
  double ltilde_PX;     lambdaPX->SetBranchAddress("ltilde",      &ltilde_PX     );
  int    positivity_PX; lambdaPX->SetBranchAddress("positivity",  &positivity_PX );


  // create CS histograms

  TH1D* h_lth_CS = new TH1D( "h_lth_CS", "", int((lth_max-lth_min)/lth_step_1D), lth_min, lth_max );
  TH1D* h_lph_CS = new TH1D( "h_lph_CS", "", int((lph_max-lph_min)/lph_step_1D), lph_min, lph_max );
  TH1D* h_ltp_CS = new TH1D( "h_ltp_CS", "", int((ltp_max-ltp_min)/ltp_step_1D), ltp_min, ltp_max );

  TH2D* h_lph_vs_lth_CS = new TH2D( "h_lph_vs_lth_CS", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
  TH2D* h_ltp_vs_lth_CS = new TH2D( "h_ltp_vs_lth_CS", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
  TH2D* h_ltp_vs_lph_CS = new TH2D( "h_ltp_vs_lph_CS", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );

  TH1D* h_ltilde_CS = new TH1D( "h_ltilde_CS", "", int((lth_max-lth_min)/lth_step_1D), lth_min, lth_max );
  TH2D* h_lphstar_vs_lthstar_CS = new TH2D( "h_lphstar_vs_lthstar_CS", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );


  // loop over entries in the CS ntuple

  n_entries = int( lambdaCS->GetEntries() );

  cout << endl;
  cout << "Reading distribution of CS parameters (" << n_entries << ") entries"<< endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << "Progress: ";

  n_step = n_entries/50;

  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

    if ( i_entry%n_step == 0 ) cout << "X";

    lambdaCS->GetEvent( i_entry );

    h_lth_CS->Fill( lth_CS );
    h_lph_CS->Fill( lph_CS );
    h_ltp_CS->Fill( ltp_CS );

    h_lph_vs_lth_CS->Fill( lth_CS, lph_CS );
    h_ltp_vs_lth_CS->Fill( lth_CS, ltp_CS );
    h_ltp_vs_lph_CS->Fill( lph_CS, ltp_CS );

    h_ltilde_CS->Fill( ltilde_CS );
    h_lphstar_vs_lthstar_CS->Fill( lthstar_CS, lphstar_CS );

  }
  // end of loop over entries in the CS ntuple

  cout << endl << endl;



  // create HX histograms

  TH1D* h_lth_HX = new TH1D( "h_lth_HX", "", int((lth_max-lth_min)/lth_step_1D), lth_min, lth_max );
  TH1D* h_lph_HX = new TH1D( "h_lph_HX", "", int((lph_max-lph_min)/lph_step_1D), lph_min, lph_max );
  TH1D* h_ltp_HX = new TH1D( "h_ltp_HX", "", int((ltp_max-ltp_min)/ltp_step_1D), ltp_min, ltp_max );

  TH2D* h_lph_vs_lth_HX = new TH2D( "h_lph_vs_lth_HX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
  TH2D* h_ltp_vs_lth_HX = new TH2D( "h_ltp_vs_lth_HX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
  TH2D* h_ltp_vs_lph_HX = new TH2D( "h_ltp_vs_lph_HX", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );

  TH1D* h_ltilde_HX = new TH1D( "h_ltilde_HX", "", int((lth_max-lth_min)/lth_step_1D), lth_min, lth_max );
  TH2D* h_lphstar_vs_lthstar_HX = new TH2D( "h_lphstar_vs_lthstar_HX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );

  // loop over entries in the HX ntuple

  n_entries = int( lambdaHX->GetEntries() );

  cout << endl;
  cout << "Reading distribution of HX parameters (" << n_entries << ") entries"<< endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << "Progress: ";

  n_step = n_entries/50;

  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

    if ( i_entry%n_step == 0 ) cout << "X";

    lambdaHX->GetEvent( i_entry );

    h_lth_HX->Fill( lth_HX );
    h_lph_HX->Fill( lph_HX );
    h_ltp_HX->Fill( ltp_HX );

    h_lph_vs_lth_HX->Fill( lth_HX, lph_HX );
    h_ltp_vs_lth_HX->Fill( lth_HX, ltp_HX );
    h_ltp_vs_lph_HX->Fill( lph_HX, ltp_HX );

    h_ltilde_HX->Fill( ltilde_HX );
    h_lphstar_vs_lthstar_HX->Fill( lthstar_HX, lphstar_HX );

  }
  // end of loop over entries in the HX ntuple

  cout << endl << endl;


  // create PX histograms

  TH1D* h_lth_PX = new TH1D( "h_lth_PX", "", int((lth_max-lth_min)/lth_step_1D), lth_min, lth_max );
  TH1D* h_lph_PX = new TH1D( "h_lph_PX", "", int((lph_max-lph_min)/lph_step_1D), lph_min, lph_max );
  TH1D* h_ltp_PX = new TH1D( "h_ltp_PX", "", int((ltp_max-ltp_min)/ltp_step_1D), ltp_min, ltp_max );

  TH2D* h_lph_vs_lth_PX = new TH2D( "h_lph_vs_lth_PX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
  TH2D* h_ltp_vs_lth_PX = new TH2D( "h_ltp_vs_lth_PX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
  TH2D* h_ltp_vs_lph_PX = new TH2D( "h_ltp_vs_lph_PX", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );

  TH1D* h_ltilde_PX = new TH1D( "h_ltilde_PX", "", int((lth_max-lth_min)/lth_step_1D), lth_min, lth_max );
  TH2D* h_lphstar_vs_lthstar_PX = new TH2D( "h_lphstar_vs_lthstar_PX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );

  // loop over entries in the PX ntuple

  n_entries = int( lambdaPX->GetEntries() );

  cout << endl;
  cout << "Reading distribution of PX parameters (" << n_entries << ") entries"<< endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << "Progress: ";

  n_step = n_entries/50;

  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

    if ( i_entry%n_step == 0 ) cout << "X";

    lambdaPX->GetEvent( i_entry );

    h_lth_PX->Fill( lth_PX );
    h_lph_PX->Fill( lph_PX );
    h_ltp_PX->Fill( ltp_PX );

    h_lph_vs_lth_PX->Fill( lth_PX, lph_PX );
    h_ltp_vs_lth_PX->Fill( lth_PX, ltp_PX );
    h_ltp_vs_lph_PX->Fill( lph_PX, ltp_PX );

    h_ltilde_PX->Fill( ltilde_PX );
    h_lphstar_vs_lthstar_PX->Fill( lthstar_PX, lphstar_PX );

  }
  // end of loop over entries in the PX ntuple

  cout << endl << endl;


  // canvas

  TCanvas *c1 = new TCanvas("c1", "c1", 10, 28, 580,571);
  c1->Range(-237.541,-66.47556,187.377,434.8609);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetLeftMargin(0.1354167);
  c1->SetRightMargin(0.01736111);
  c1->SetTopMargin(0.01841621);
  c1->SetBottomMargin(0.1325967);
  c1->SetFrameBorderMode(0);


  ////////////////////////////////////////////////////////////////
  // 2D plot lph vs lth

  TH2D* h_lph_vs_lth = (TH2D*)h_lph_vs_lth_CS->Clone();
  h_lph_vs_lth->SetName("h_lph_vs_lth");
  h_lph_vs_lth->Reset();

  h_lph_vs_lth->GetXaxis()->SetTitle("#lambda_{#vartheta}");
  h_lph_vs_lth->GetXaxis()->SetLabelOffset(0.028);
  h_lph_vs_lth->GetXaxis()->SetTitleSize(0.05);
  h_lph_vs_lth->GetXaxis()->SetTickLength(-0.03);
  h_lph_vs_lth->GetXaxis()->SetTitleOffset(1.15);
  h_lph_vs_lth->GetYaxis()->SetTitle("#lambda_{#varphi}");
  h_lph_vs_lth->GetYaxis()->SetLabelOffset(0.032);
  h_lph_vs_lth->GetYaxis()->SetTitleSize(0.05);
  h_lph_vs_lth->GetYaxis()->SetTickLength(-0.03);
  h_lph_vs_lth->GetYaxis()->SetTitleOffset(1.3);
  h_lph_vs_lth->Draw("");

  double x_background_ph_vs_th[4] = { lth_min+0.01, lth_min+0.01, lth_max-0.01, lth_max-0.01 };
  double y_background_ph_vs_th[4] = { lph_min+0.01, lph_max-0.01, lph_max-0.01, lph_min+0.01 };
  TPolyLine *background_ph_vs_th = new TPolyLine( 4, x_background_ph_vs_th, y_background_ph_vs_th );
  background_ph_vs_th->SetFillColor(kGray);
  background_ph_vs_th->SetLineStyle(0);
  background_ph_vs_th->Draw("f same");

  double x_triangle_ph_vs_th[3] = {-1.,  1., 1.};
  double y_triangle_ph_vs_th[3] = { 0., -1., 1.};
  TPolyLine *triangle_ph_vs_th = new TPolyLine( 3, x_triangle_ph_vs_th, y_triangle_ph_vs_th );
  triangle_ph_vs_th->SetFillColor(kWhite);
  triangle_ph_vs_th->SetLineStyle(0);
  triangle_ph_vs_th->Draw("f same");

  TLine* lh_p1 = new TLine( lth_min, 1., lth_max, 1. );
  lh_p1->SetLineWidth( 1 );
  lh_p1->SetLineStyle( 2 );
  lh_p1->SetLineColor( kGray+2 );
  lh_p1->Draw( "same" );

  TLine* lh_m1 = new TLine( lth_min, -1., lth_max, -1. );
  lh_m1->SetLineWidth( 1 );
  lh_m1->SetLineStyle( 2 );
  lh_m1->SetLineColor( kGray+2 );
  lh_m1->Draw( "same" );

  TLine* lh_0 = new TLine( lth_min, 0., lth_max, 0. );
  lh_0->SetLineWidth( 1 );
  lh_0->SetLineStyle( 2 );
  lh_0->SetLineColor( kGray+2 );
  lh_0->Draw( "same" );

  TLine* lv_p1 = new TLine( 1., lph_min, 1., lph_max );
  lv_p1->SetLineWidth( 1 );
  lv_p1->SetLineStyle( 2 );
  lv_p1->SetLineColor( kGray+2 );
  lv_p1->Draw( "same" );

  TLine* lv_m1 = new TLine( -1., lph_min, -1., lph_max );
  lv_m1->SetLineWidth( 1 );
  lv_m1->SetLineStyle( 2 );
  lv_m1->SetLineColor( kGray+2 );
  lv_m1->Draw( "same" );

  TLine* lv_0 = new TLine( 0., lph_min, 0., lph_max );
  lv_0->SetLineWidth( 1 );
  lv_0->SetLineStyle( 2 );
  lv_0->SetLineColor( kGray+2 );
  lv_0->Draw( "same" );


  // CS frame

  setContourHistogram ( h_lph_vs_lth_CS );
  h_lph_vs_lth_CS->SetLineColor( kBlue );
  h_lph_vs_lth_CS->SetLineWidth( 2 );
  h_lph_vs_lth_CS->SetLineStyle( 1 );
  h_lph_vs_lth_CS->Draw( "cont2, same" );

  // HX frame

  setContourHistogram ( h_lph_vs_lth_HX );
  h_lph_vs_lth_HX->SetLineColor( kRed );
  h_lph_vs_lth_HX->SetLineWidth( 2 );
  h_lph_vs_lth_HX->SetLineStyle( 1 );
  h_lph_vs_lth_HX->Draw( "cont2, same" );

  // PX frame

  setContourHistogram ( h_lph_vs_lth_PX );
  h_lph_vs_lth_PX->SetLineColor( kGreen+2 );
  h_lph_vs_lth_PX->SetLineWidth( 2 );
  h_lph_vs_lth_PX->SetLineStyle( 1 );
  h_lph_vs_lth_PX->Draw( "cont2, same" );


  sprintf(filename,"%s/lph_vs_lth.eps",dirstruct);
  c1->Print( filename );




  ////////////////////////////////////////////////////////////////
  // 2D plot lphstar vs lthstar

  TH2D* h_lphstar_vs_lthstar = (TH2D*)h_lphstar_vs_lthstar_CS->Clone();
  h_lphstar_vs_lthstar->SetName("h_lphstar_vs_lthstar");
  h_lphstar_vs_lthstar->Reset();

  h_lphstar_vs_lthstar->GetXaxis()->SetTitle("#lambda_{#vartheta}^{*}");
  h_lphstar_vs_lthstar->GetXaxis()->SetLabelOffset(0.028);
  h_lphstar_vs_lthstar->GetXaxis()->SetTitleSize(0.05);
  h_lphstar_vs_lthstar->GetXaxis()->SetTickLength(-0.03);
  h_lphstar_vs_lthstar->GetXaxis()->SetTitleOffset(1.15);
  h_lphstar_vs_lthstar->GetYaxis()->SetTitle("#lambda_{#varphi}^{*}");
  h_lphstar_vs_lthstar->GetYaxis()->SetLabelOffset(0.032);
  h_lphstar_vs_lthstar->GetYaxis()->SetTitleSize(0.05);
  h_lphstar_vs_lthstar->GetYaxis()->SetTickLength(-0.03);
  h_lphstar_vs_lthstar->GetYaxis()->SetTitleOffset(1.3);
  h_lphstar_vs_lthstar->Draw("");

  background_ph_vs_th->Draw("f same");
  triangle_ph_vs_th->Draw("f same");

  lh_p1->Draw( "same" );
  lh_m1->Draw( "same" );
  lh_0->Draw( "same" );
  lv_p1->Draw( "same" );
  lv_m1->Draw( "same" );
  lv_0->Draw( "same" );


  // CS frame

  setContourHistogram ( h_lphstar_vs_lthstar_CS );
  h_lphstar_vs_lthstar_CS->SetLineColor( kBlue );
  h_lphstar_vs_lthstar_CS->SetLineWidth( 2 );
  h_lphstar_vs_lthstar_CS->SetLineStyle( 1 );
  h_lphstar_vs_lthstar_CS->Draw( "cont2, same" );

  // HX frame

  setContourHistogram ( h_lphstar_vs_lthstar_HX );
  h_lphstar_vs_lthstar_HX->SetLineColor( kRed );
  h_lphstar_vs_lthstar_HX->SetLineWidth( 2 );
  h_lphstar_vs_lthstar_HX->SetLineStyle( 1 );
  h_lphstar_vs_lthstar_HX->Draw( "cont2, same" );

  // PX frame

  setContourHistogram ( h_lphstar_vs_lthstar_PX );
  h_lphstar_vs_lthstar_PX->SetLineColor( kGreen+2 );
  h_lphstar_vs_lthstar_PX->SetLineWidth( 2 );
  h_lphstar_vs_lthstar_PX->SetLineStyle( 1 );
  h_lphstar_vs_lthstar_PX->Draw( "cont2, same" );


  sprintf(filename,"%s/lphstar_vs_lthstar.eps",dirstruct);
  c1->Print( filename );




  ////////////////////////////////////////////////////////////////
  // 2D plot ltp vs lth

  TH2D* h_ltp_vs_lth = (TH2D*)h_ltp_vs_lth_CS->Clone();
  h_ltp_vs_lth->SetName("h_ltp_vs_lth");
  h_ltp_vs_lth->Reset();

  h_ltp_vs_lth->GetXaxis()->SetTitle("#lambda_{#vartheta}");
  h_ltp_vs_lth->GetXaxis()->SetLabelOffset(0.028);
  h_ltp_vs_lth->GetXaxis()->SetTitleSize(0.05);
  h_ltp_vs_lth->GetXaxis()->SetTickLength(-0.03);
  h_ltp_vs_lth->GetXaxis()->SetTitleOffset(1.15);
  h_ltp_vs_lth->GetYaxis()->SetTitle("#lambda_{#vartheta#varphi}");
  h_ltp_vs_lth->GetYaxis()->SetLabelOffset(0.032);
  h_ltp_vs_lth->GetYaxis()->SetTitleSize(0.05);
  h_ltp_vs_lth->GetYaxis()->SetTickLength(-0.03);
  h_ltp_vs_lth->GetYaxis()->SetTitleOffset(1.3);
  h_ltp_vs_lth->Draw("");

  double x_background_tp_vs_th[4] = { lth_min+0.01, lth_min+0.01, lth_max-0.01, lth_max-0.01 };
  double y_background_tp_vs_th[4] = { ltp_min+0.01, ltp_max-0.01, ltp_max-0.01, ltp_min+0.01 };
  TPolyLine *background_tp_vs_th = new TPolyLine( 4, x_background_tp_vs_th, y_background_tp_vs_th );
  background_tp_vs_th->SetFillColor(kGray);
  background_tp_vs_th->SetLineStyle(0);
  background_tp_vs_th->Draw("f same");

  TEllipse *ellipse_tp_vs_th = new TEllipse( 0., 0., 1., sqrt(2.)/2. );
  ellipse_tp_vs_th->SetFillColor(kWhite);
  ellipse_tp_vs_th->SetLineStyle(0);
  ellipse_tp_vs_th->Draw("f same");


  lh_p1 = new TLine( lth_min, sqrt(2.)/2., lth_max, sqrt(2.)/2. );
  lh_p1->SetLineWidth( 1 );
  lh_p1->SetLineStyle( 2 );
  lh_p1->SetLineColor( kGray+2 );
  lh_p1->Draw( "same" );

  lh_m1 = new TLine( lth_min, -sqrt(2.)/2., lth_max, -sqrt(2.)/2. );
  lh_m1->SetLineWidth( 1 );
  lh_m1->SetLineStyle( 2 );
  lh_m1->SetLineColor( kGray+2 );
  lh_m1->Draw( "same" );

  lh_0 = new TLine( lth_min, 0., lth_max, 0. );
  lh_0->SetLineWidth( 1 );
  lh_0->SetLineStyle( 2 );
  lh_0->SetLineColor( kGray+2 );
  lh_0->Draw( "same" );

  lv_p1 = new TLine( 1., ltp_min, 1., ltp_max );
  lv_p1->SetLineWidth( 1 );
  lv_p1->SetLineStyle( 2 );
  lv_p1->SetLineColor( kGray+2 );
  lv_p1->Draw( "same" );

  lv_m1 = new TLine( -1., ltp_min, -1., ltp_max );
  lv_m1->SetLineWidth( 1 );
  lv_m1->SetLineStyle( 2 );
  lv_m1->SetLineColor( kGray+2 );
  lv_m1->Draw( "same" );

  lv_0 = new TLine( 0., ltp_min, 0., ltp_max );
  lv_0->SetLineWidth( 1 );
  lv_0->SetLineStyle( 2 );
  lv_0->SetLineColor( kGray+2 );
  lv_0->Draw( "same" );


  // CS frame

  setContourHistogram ( h_ltp_vs_lth_CS );
  h_ltp_vs_lth_CS->SetLineColor( kBlue );
  h_ltp_vs_lth_CS->SetLineWidth( 2 );
  h_ltp_vs_lth_CS->SetLineStyle( 1 );
  h_ltp_vs_lth_CS->Draw( "cont2, same" );

  // HX frame

  setContourHistogram ( h_ltp_vs_lth_HX );
  h_ltp_vs_lth_HX->SetLineColor( kRed );
  h_ltp_vs_lth_HX->SetLineWidth( 2 );
  h_ltp_vs_lth_HX->SetLineStyle( 1 );
  h_ltp_vs_lth_HX->Draw( "cont2, same" );

  // PX frame

  setContourHistogram ( h_ltp_vs_lth_PX );
  h_ltp_vs_lth_PX->SetLineColor( kGreen+2 );
  h_ltp_vs_lth_PX->SetLineWidth( 2 );
  h_ltp_vs_lth_PX->SetLineStyle( 1 );
  h_ltp_vs_lth_PX->Draw( "cont2, same" );


  sprintf(filename,"%s/ltp_vs_lth.eps",dirstruct);
  c1->Print( filename );


  ////////////////////////////////////////////////////////////////
  // 2D plot ltp vs lth

  TH2D* h_ltp_vs_lph = (TH2D*)h_ltp_vs_lph_CS->Clone();
  h_ltp_vs_lph->SetName("h_ltp_vs_lph");
  h_ltp_vs_lph->Reset();

  h_ltp_vs_lph->GetXaxis()->SetTitle("#lambda_{#varphi}");
  h_ltp_vs_lph->GetXaxis()->SetLabelOffset(0.028);
  h_ltp_vs_lph->GetXaxis()->SetTitleSize(0.05);
  h_ltp_vs_lph->GetXaxis()->SetTickLength(-0.03);
  h_ltp_vs_lph->GetXaxis()->SetTitleOffset(1.15);
  h_ltp_vs_lph->GetYaxis()->SetTitle("#lambda_{#vartheta#varphi}");
  h_ltp_vs_lph->GetYaxis()->SetLabelOffset(0.032);
  h_ltp_vs_lph->GetYaxis()->SetTitleSize(0.05);
  h_ltp_vs_lph->GetYaxis()->SetTickLength(-0.03);
  h_ltp_vs_lph->GetYaxis()->SetTitleOffset(1.3);
  h_ltp_vs_lph->Draw("");

  double x_background_tp_vs_ph[4] = { lph_min+0.01, lph_min+0.01, lph_max-0.01, lph_max-0.01 };
  double y_background_tp_vs_ph[4] = { ltp_min+0.01, ltp_max-0.01, ltp_max-0.01, ltp_min+0.01 };
  TPolyLine *background_tp_vs_ph = new TPolyLine( 4, x_background_tp_vs_ph, y_background_tp_vs_ph );
  background_tp_vs_ph->SetFillColor(kGray);
  background_tp_vs_ph->SetLineStyle(0);
  background_tp_vs_ph->Draw("f same");

  TEllipse *ellipse_tp_vs_ph = new TEllipse( -0.5, 0., 0.5, sqrt(2.)/2. );
  ellipse_tp_vs_ph->SetFillColor(kWhite);
  ellipse_tp_vs_ph->SetLineStyle(0);
  ellipse_tp_vs_ph->Draw("f same");

  double x_triangle_tp_vs_th[3] = {-1./3., -1./3., 1.};
  double y_triangle_tp_vs_th[3] = { 2./3., -2./3., 0.};
  TPolyLine *triangle_tp_vs_th = new TPolyLine( 3, x_triangle_tp_vs_th, y_triangle_tp_vs_th );
  triangle_tp_vs_th->SetFillColor(kWhite);
  triangle_tp_vs_th->SetLineStyle(0);
  triangle_tp_vs_th->Draw("f same");

  lh_p1 = new TLine( lph_min, sqrt(2.)/2., lph_max, sqrt(2.)/2. );
  lh_p1->SetLineWidth( 1 );
  lh_p1->SetLineStyle( 2 );
  lh_p1->SetLineColor( kGray+2 );
  lh_p1->Draw( "same" );

  lh_m1 = new TLine( lph_min, -sqrt(2.)/2., lph_max, -sqrt(2.)/2. );
  lh_m1->SetLineWidth( 1 );
  lh_m1->SetLineStyle( 2 );
  lh_m1->SetLineColor( kGray+2 );
  lh_m1->Draw( "same" );

  lh_0 = new TLine( lph_min, 0., lph_max, 0. );
  lh_0->SetLineWidth( 1 );
  lh_0->SetLineStyle( 2 );
  lh_0->SetLineColor( kGray+2 );
  lh_0->Draw( "same" );

  lv_p1 = new TLine( 1., ltp_min, 1., ltp_max );
  lv_p1->SetLineWidth( 1 );
  lv_p1->SetLineStyle( 2 );
  lv_p1->SetLineColor( kGray+2 );
  lv_p1->Draw( "same" );

  lv_m1 = new TLine( -1., ltp_min, -1., ltp_max );
  lv_m1->SetLineWidth( 1 );
  lv_m1->SetLineStyle( 2 );
  lv_m1->SetLineColor( kGray+2 );
  lv_m1->Draw( "same" );

  lv_0 = new TLine( 0., ltp_min, 0., ltp_max );
  lv_0->SetLineWidth( 1 );
  lv_0->SetLineStyle( 2 );
  lv_0->SetLineColor( kGray+2 );
  lv_0->Draw( "same" );


  // CS frame

  setContourHistogram ( h_ltp_vs_lph_CS );
  h_ltp_vs_lph_CS->SetLineColor( kBlue );
  h_ltp_vs_lph_CS->SetLineWidth( 2 );
  h_ltp_vs_lph_CS->SetLineStyle( 1 );
  h_ltp_vs_lph_CS->Draw( "cont2, same" );

  // HX frame

  setContourHistogram ( h_ltp_vs_lph_HX );
  h_ltp_vs_lph_HX->SetLineColor( kRed );
  h_ltp_vs_lph_HX->SetLineWidth( 2 );
  h_ltp_vs_lph_HX->SetLineStyle( 1 );
  h_ltp_vs_lph_HX->Draw( "cont2, same" );

  // PX frame

  setContourHistogram ( h_ltp_vs_lph_PX );
  h_ltp_vs_lph_PX->SetLineColor( kGreen+2 );
  h_ltp_vs_lph_PX->SetLineWidth( 2 );
  h_ltp_vs_lph_PX->SetLineStyle( 1 );
  h_ltp_vs_lph_PX->Draw( "cont2, same" );


  sprintf(filename,"%s/ltp_vs_lph.eps",dirstruct);
  c1->Print( filename );




  ///////////////////////////////////////////////////////////////////////////////
  // plot PDFs for lambdatilde in three frames

  c1 = new TCanvas("c1", "c1", 10, 28, 588,563);
  c1->Range(-1.370833,-114.012,1.029167,745.8285);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetLeftMargin(0.1545139);
  c1->SetRightMargin(0.01215278);
  c1->SetTopMargin(0.01841621);
  c1->SetBottomMargin(0.1325967);
  c1->SetFrameBorderMode(0);

  TH1D* h_ltilde = (TH1D*)h_ltilde_CS->Clone();
  h_ltilde->SetName("h_ltilde");
  h_ltilde->Reset();

  h_ltilde_CS->Scale( 1./h_ltilde_CS->Integral() );
  h_ltilde_HX->Scale( 1./h_ltilde_HX->Integral() );
  h_ltilde_PX->Scale( 1./h_ltilde_PX->Integral() );

  double plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  double plotMax = plotborder * h_ltilde_CS->GetMaximum();
  if ( plotborder * h_ltilde_HX->GetMaximum() > plotMax ) plotMax = plotborder * h_ltilde_HX->GetMaximum();
  if ( plotborder * h_ltilde_PX->GetMaximum() > plotMax ) plotMax = plotborder * h_ltilde_PX->GetMaximum();

  h_ltilde->GetXaxis()->SetTitle("#tilde{#lambda}");
  h_ltilde->GetXaxis()->SetLabelOffset(0.028);
  h_ltilde->GetXaxis()->SetTitleSize(0.05);
  h_ltilde->GetXaxis()->SetTickLength(-0.03);
  h_ltilde->GetXaxis()->SetTitleOffset(1.20);
  h_ltilde->GetYaxis()->SetTitle("parameter PDF [a.u.]");
  h_ltilde->GetYaxis()->SetLabelOffset(0.032);
  h_ltilde->GetYaxis()->SetTitleSize(0.05);
  h_ltilde->GetYaxis()->SetTickLength(-0.03);
  h_ltilde->GetYaxis()->SetTitleOffset(1.55);
  h_ltilde->SetMinimum(0.);
  h_ltilde->SetMaximum(plotMax);
  h_ltilde->Draw("");

  // CS frame

  h_ltilde_CS->SetLineColor( kBlue );
  h_ltilde_CS->SetLineWidth( 2 );
  h_ltilde_CS->SetLineStyle( 1 );
  h_ltilde_CS->Draw( "L same" );

  // HX frame

  h_ltilde_HX->SetLineColor( kRed );
  h_ltilde_HX->SetLineWidth( 2 );
  h_ltilde_HX->SetLineStyle( 1 );
  h_ltilde_HX->Draw( "L same" );

  // PX frame

  h_ltilde_PX->SetLineColor( kGreen+2 );
  h_ltilde_PX->SetLineWidth( 2 );
  h_ltilde_PX->SetLineStyle( 1 );
  h_ltilde_PX->Draw( "L same" );


  sprintf(filename,"%s/ltilde.eps",dirstruct);
  c1->Print( filename );


  ///////////////////////////////////////////////////////////////////////////
  // extract best-fit values from the 1D histograms (taking mean)

  double lth_CS_best = h_lth_CS->GetMean();
  double lph_CS_best = h_lph_CS->GetMean();
  double ltp_CS_best = h_ltp_CS->GetMean();

  double lth_HX_best = h_lth_HX->GetMean();
  double lph_HX_best = h_lph_HX->GetMean();
  double ltp_HX_best = h_ltp_HX->GetMean();

  double lth_PX_best = h_lth_PX->GetMean();
  double lph_PX_best = h_lph_PX->GetMean();
  double ltp_PX_best = h_ltp_PX->GetMean();

  double ltilde_CS_best = h_ltilde_CS->GetMean();
  double ltilde_HX_best = h_ltilde_HX->GetMean();
  double ltilde_PX_best = h_ltilde_PX->GetMean();

  double ltilde_best = (ltilde_CS_best + ltilde_HX_best + ltilde_PX_best)/3.;


  // input ntuple with anguar distribution of background-subtracted dilepton events

  TTree* angles = (TTree*)results->Get("angles");

  double costh_CS;  angles->SetBranchAddress("costh_CS",     &costh_CS );
  double phi_CS;    angles->SetBranchAddress("phi_CS",       &phi_CS   );
  double phith_CS;  angles->SetBranchAddress("phithCS",      &phith_CS );

  double costh_HX;  angles->SetBranchAddress("costh_HX",     &costh_HX );
  double phi_HX;    angles->SetBranchAddress("phi_HX",       &phi_HX   );
  double phith_HX;  angles->SetBranchAddress("phith_HX",     &phith_HX );

  double costh_PX;  angles->SetBranchAddress("costh_PX",     &costh_PX );
  double phi_PX;    angles->SetBranchAddress("phi_PX",       &phi_PX   );
  double phith_PX;  angles->SetBranchAddress("phith_PX",     &phith_PX );

  double cosalpha;  angles->SetBranchAddress("cosalpha",     &cosalpha );
  double epsilon;   angles->SetBranchAddress("epsilon",      &epsilon  );


  // input histograms with angular dependencies of the parameter-independent PDF terms

  TH1D* PDF_1_vs_cth_CS = (TH1D*)results->Get("PDF_1_vs_cth_CS");
  TH1D* PDF_cth2_vs_cth_CS = (TH1D*)results->Get("PDF_cth2_vs_cth_CS");

  TH1D* PDF_1_vs_ph_CS = (TH1D*)results->Get("PDF_1_vs_ph_CS");
  TH1D* PDF_c2ph_vs_ph_CS = (TH1D*)results->Get("PDF_c2ph_vs_ph_CS");

  TH1D* PDF_1_vs_phth_CS = (TH1D*)results->Get("PDF_1_vs_phth_CS");
  TH1D* PDF_cphth_vs_phth_CS = (TH1D*)results->Get("PDF_cphth_vs_phth_CS");


  TH1D* PDF_1_vs_cth_HX = (TH1D*)results->Get("PDF_1_vs_cth_HX");
  TH1D* PDF_cth2_vs_cth_HX = (TH1D*)results->Get("PDF_cth2_vs_cth_HX");

  TH1D* PDF_1_vs_ph_HX = (TH1D*)results->Get("PDF_1_vs_ph_HX");
  TH1D* PDF_c2ph_vs_ph_HX = (TH1D*)results->Get("PDF_c2ph_vs_ph_HX");

  TH1D* PDF_1_vs_phth_HX = (TH1D*)results->Get("PDF_1_vs_phth_HX");
  TH1D* PDF_cphth_vs_phth_HX = (TH1D*)results->Get("PDF_cphth_vs_phth_HX");


  TH1D* PDF_1_vs_cth_PX = (TH1D*)results->Get("PDF_1_vs_cth_PX");
  TH1D* PDF_cth2_vs_cth_PX = (TH1D*)results->Get("PDF_cth2_vs_cth_PX");

  TH1D* PDF_1_vs_ph_PX = (TH1D*)results->Get("PDF_1_vs_ph_PX");
  TH1D* PDF_c2ph_vs_ph_PX = (TH1D*)results->Get("PDF_c2ph_vs_ph_PX");

  TH1D* PDF_1_vs_phth_PX = (TH1D*)results->Get("PDF_1_vs_phth_PX");
  TH1D* PDF_cphth_vs_phth_PX = (TH1D*)results->Get("PDF_cphth_vs_phth_PX");


  // create data histograms

  TH1D* DATA_cth_CS   = new TH1D( "DATA_cth_CS",  "", nbin_cth, -1., 1. );
  TH1D* DATA_ph_CS    = new TH1D( "DATA_ph_CS",   "", nbin_ph, -180., 180. );
  TH1D* DATA_phth_CS  = new TH1D( "DATA_phth_CS", "", nbin_ph, -180., 180. );

  TH1D* DATA_cth_HX   = new TH1D( "DATA_cth_HX",  "", nbin_cth, -1., 1. );
  TH1D* DATA_ph_HX    = new TH1D( "DATA_ph_HX",   "", nbin_ph, -180., 180. );
  TH1D* DATA_phth_HX  = new TH1D( "DATA_phth_HX", "", nbin_ph, -180., 180. );

  TH1D* DATA_cth_PX   = new TH1D( "DATA_cth_PX",  "", nbin_cth, -1., 1. );
  TH1D* DATA_ph_PX    = new TH1D( "DATA_ph_PX",   "", nbin_ph, -180., 180. );
  TH1D* DATA_phth_PX  = new TH1D( "DATA_phth_PX", "", nbin_ph, -180., 180. );


  // loop over entries in the ntuple of angles

  n_entries = int( angles->GetEntries() );

  cout << endl;
  cout << "Reading angular distributions (" << n_entries << ") entries"<< endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << "Progress: ";

  n_step = n_entries/50;

  for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

    if ( i_entry%n_step == 0 ) cout << "X";

    angles->GetEvent( i_entry );

    DATA_cth_CS->Fill( costh_CS );
    DATA_ph_CS->Fill( phi_CS );
    DATA_phth_CS->Fill( phith_CS );

    DATA_cth_HX->Fill( costh_HX );
    DATA_ph_HX->Fill( phi_HX );
    DATA_phth_HX->Fill( phith_HX );

    DATA_cth_PX->Fill( costh_PX );
    DATA_ph_PX->Fill( phi_PX );
    DATA_phth_PX->Fill( phith_PX );

  }
  // end of loop over entries in the ntuple of angles

  cout << endl << endl;


  // plot angular PDFs (calculated combining parameter-independent terms with best-fit parameters)
  // normalized to data distribution

  c1 = new TCanvas("c1", "c1", 10, 28, 588,563);
  c1->Range(-1.370833,-114.012,1.029167,745.8285);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetLeftMargin(0.1545139);
  c1->SetRightMargin(0.01215278);
  c1->SetTopMargin(0.01841621);
  c1->SetBottomMargin(0.1325967);
  c1->SetFrameBorderMode(0);


  // costh_CS

  TH1D* PDF_cth_CS = (TH1D*)PDF_1_vs_cth_CS->Clone();
  PDF_cth_CS->SetName("PDF_cth_CS");
  PDF_cth_CS->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., lth_CS_best );
  PDF_cth_CS->Scale( DATA_cth_CS->Integral() * PDF_cth_CS->GetNbinsX() / ( PDF_cth_CS->Integral()*DATA_cth_CS->GetNbinsX() ) );

  TH1D* PDF_cth_CS_p1 = (TH1D*)PDF_1_vs_cth_CS->Clone();
  PDF_cth_CS_p1->SetName("PDF_cth_CS_p1");
  PDF_cth_CS_p1->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., 1. );
  PDF_cth_CS_p1->Scale( DATA_cth_CS->Integral() * PDF_cth_CS_p1->GetNbinsX() / ( PDF_cth_CS_p1->Integral()*DATA_cth_CS->GetNbinsX() ) );

  TH1D* PDF_cth_CS_m1 = (TH1D*)PDF_1_vs_cth_CS->Clone();
  PDF_cth_CS_m1->SetName("PDF_cth_CS_m1");
  PDF_cth_CS_m1->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., -1. );
  PDF_cth_CS_m1->Scale( DATA_cth_CS->Integral() * PDF_cth_CS_m1->GetNbinsX() / ( PDF_cth_CS_m1->Integral()*DATA_cth_CS->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_cth_CS->GetMaximum();
  if ( plotborder * PDF_cth_CS_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_CS_p1->GetMaximum();
  if ( plotborder * PDF_cth_CS_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_CS_m1->GetMaximum();

  PDF_cth_CS->GetXaxis()->SetTitle("cos#vartheta_{CS}");
  PDF_cth_CS->GetXaxis()->SetLabelOffset(0.028);
  PDF_cth_CS->GetXaxis()->SetTitleSize(0.05);
  PDF_cth_CS->GetXaxis()->SetTickLength(-0.03);
  PDF_cth_CS->GetXaxis()->SetTitleOffset(1.20);
  PDF_cth_CS->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_cth_CS->GetYaxis()->SetLabelOffset(0.032);
  PDF_cth_CS->GetYaxis()->SetTitleSize(0.05);
  PDF_cth_CS->GetYaxis()->SetTickLength(-0.03);
  PDF_cth_CS->GetYaxis()->SetTitleOffset(1.55);
  PDF_cth_CS->SetMinimum(0.);
  PDF_cth_CS->SetLineColor(kRed);
  PDF_cth_CS->SetLineStyle(1);
  PDF_cth_CS->SetLineWidth(2);
  PDF_cth_CS->SetMaximum(plotMax);
  PDF_cth_CS->Draw("L");

  PDF_cth_CS_p1->SetLineColor(kRed);
  PDF_cth_CS_p1->SetLineStyle(2);
  PDF_cth_CS_p1->Draw("L same");

  PDF_cth_CS_m1->SetLineColor(kRed);
  PDF_cth_CS_m1->SetLineStyle(10);
  PDF_cth_CS_m1->Draw("L same");

  DATA_cth_CS->SetLineColor(kBlack);
  DATA_cth_CS->Draw("E same");

  sprintf(filename,"%s/fit_CS_costh.eps",dirstruct);
  c1->Print( filename );


  // phi_CS

  TH1D* PDF_ph_CS = (TH1D*)PDF_1_vs_ph_CS->Clone();
  PDF_ph_CS->SetName("PDF_ph_CS");
  PDF_ph_CS->Add( PDF_1_vs_ph_CS, PDF_c2ph_vs_ph_CS, 1., 2.*lph_CS_best/(3.+lth_CS_best) );
  PDF_ph_CS->Scale( DATA_ph_CS->Integral() * PDF_ph_CS->GetNbinsX() / ( PDF_ph_CS->Integral()*DATA_ph_CS->GetNbinsX() ) );

  TH1D* PDF_ph_CS_p1 = (TH1D*)PDF_1_vs_ph_CS->Clone();
  PDF_ph_CS_p1->SetName("PDF_ph_CS_p1");
  PDF_ph_CS_p1->Add( PDF_1_vs_ph_CS, PDF_c2ph_vs_ph_CS, 1., 0.5 );
  PDF_ph_CS_p1->Scale( DATA_ph_CS->Integral() * PDF_ph_CS_p1->GetNbinsX() / ( PDF_ph_CS_p1->Integral()*DATA_ph_CS->GetNbinsX() ) );

  TH1D* PDF_ph_CS_m1 = (TH1D*)PDF_1_vs_ph_CS->Clone();
  PDF_ph_CS_m1->SetName("PDF_ph_CS_m1");
  PDF_ph_CS_m1->Add( PDF_1_vs_ph_CS, PDF_c2ph_vs_ph_CS, 1., -0.5 );
  PDF_ph_CS_m1->Scale( DATA_ph_CS->Integral() * PDF_ph_CS_m1->GetNbinsX() / ( PDF_ph_CS_m1->Integral()*DATA_ph_CS->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_ph_CS->GetMaximum();
  if ( plotborder * PDF_ph_CS_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_CS_p1->GetMaximum();
  if ( plotborder * PDF_ph_CS_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_CS_m1->GetMaximum();

  PDF_ph_CS->GetXaxis()->SetTitle("#varphi_{CS}");
  PDF_ph_CS->GetXaxis()->SetLabelOffset(0.028);
  PDF_ph_CS->GetXaxis()->SetTitleSize(0.05);
  PDF_ph_CS->GetXaxis()->SetTickLength(-0.03);
  PDF_ph_CS->GetXaxis()->SetTitleOffset(1.20);
  PDF_ph_CS->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_ph_CS->GetYaxis()->SetLabelOffset(0.032);
  PDF_ph_CS->GetYaxis()->SetTitleSize(0.05);
  PDF_ph_CS->GetYaxis()->SetTickLength(-0.03);
  PDF_ph_CS->GetYaxis()->SetTitleOffset(1.55);
  PDF_ph_CS->SetMinimum(0.);
  PDF_ph_CS->SetLineColor(kRed);
  PDF_ph_CS->SetLineStyle(1);
  PDF_ph_CS->SetLineWidth(2);
  PDF_ph_CS->SetMaximum(plotMax);
  PDF_ph_CS->Draw("L");

  PDF_ph_CS_p1->SetLineColor(kRed);
  PDF_ph_CS_p1->SetLineStyle(2);
  PDF_ph_CS_p1->Draw("L same");

  PDF_ph_CS_m1->SetLineColor(kRed);
  PDF_ph_CS_m1->SetLineStyle(10);
  PDF_ph_CS_m1->Draw("L same");

  DATA_ph_CS->SetLineColor(kBlack);
  DATA_ph_CS->Draw("E same");

  sprintf(filename,"%s/fit_CS_phi.eps",dirstruct);
  c1->Print( filename );


  // phith_CS

  TH1D* PDF_phth_CS = (TH1D*)PDF_1_vs_phth_CS->Clone();
  PDF_phth_CS->SetName("PDF_phth_CS");
  PDF_phth_CS->Add( PDF_1_vs_phth_CS, PDF_cphth_vs_phth_CS, 1., sqrt(2.)*ltp_CS_best/(3.+lth_CS_best) );
  PDF_phth_CS->Scale( DATA_phth_CS->Integral() * PDF_phth_CS->GetNbinsX() / ( PDF_phth_CS->Integral()*DATA_phth_CS->GetNbinsX() ) );

  TH1D* PDF_phth_CS_p1 = (TH1D*)PDF_1_vs_phth_CS->Clone();
  PDF_phth_CS_p1->SetName("PDF_phth_CS_p1");
  PDF_phth_CS_p1->Add( PDF_1_vs_phth_CS, PDF_cphth_vs_phth_CS, 1., sqrt(2.)/4. );
  PDF_phth_CS_p1->Scale( DATA_phth_CS->Integral() * PDF_phth_CS_p1->GetNbinsX() / ( PDF_phth_CS_p1->Integral()*DATA_phth_CS->GetNbinsX() ) );

  TH1D* PDF_phth_CS_m1 = (TH1D*)PDF_1_vs_phth_CS->Clone();
  PDF_phth_CS_m1->SetName("PDF_phth_CS_m1");
  PDF_phth_CS_m1->Add( PDF_1_vs_phth_CS, PDF_cphth_vs_phth_CS, 1., -sqrt(2.)/4. );
  PDF_phth_CS_m1->Scale( DATA_phth_CS->Integral() * PDF_phth_CS_m1->GetNbinsX() / ( PDF_phth_CS_m1->Integral()*DATA_phth_CS->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_phth_CS->GetMaximum();
  if ( plotborder * PDF_phth_CS_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_CS_p1->GetMaximum();
  if ( plotborder * PDF_phth_CS_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_CS_m1->GetMaximum();

  PDF_phth_CS->GetXaxis()->SetTitle("#tilde{#varphi}_{CS}");
  PDF_phth_CS->GetXaxis()->SetLabelOffset(0.028);
  PDF_phth_CS->GetXaxis()->SetTitleSize(0.05);
  PDF_phth_CS->GetXaxis()->SetTickLength(-0.03);
  PDF_phth_CS->GetXaxis()->SetTitleOffset(1.20);
  PDF_phth_CS->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_phth_CS->GetYaxis()->SetLabelOffset(0.032);
  PDF_phth_CS->GetYaxis()->SetTitleSize(0.05);
  PDF_phth_CS->GetYaxis()->SetTickLength(-0.03);
  PDF_phth_CS->GetYaxis()->SetTitleOffset(1.55);
  PDF_phth_CS->SetMinimum(0.);
  PDF_phth_CS->SetLineColor(kRed);
  PDF_phth_CS->SetLineStyle(1);
  PDF_phth_CS->SetLineWidth(2);
  PDF_phth_CS->SetMaximum(plotMax);
  PDF_phth_CS->Draw("L");

  PDF_phth_CS_p1->SetLineColor(kRed);
  PDF_phth_CS_p1->SetLineStyle(2);
  PDF_phth_CS_p1->Draw("L same");

  PDF_phth_CS_m1->SetLineColor(kRed);
  PDF_phth_CS_m1->SetLineStyle(10);
  PDF_phth_CS_m1->Draw("L same");

  DATA_phth_CS->SetLineColor(kBlack);
  DATA_phth_CS->Draw("E same");

  sprintf(filename,"%s/fit_CS_phith.eps",dirstruct);
  c1->Print( filename );

  
  // costh_HX

  TH1D* PDF_cth_HX = (TH1D*)PDF_1_vs_cth_HX->Clone();
  PDF_cth_HX->SetName("PDF_cth_HX");
  PDF_cth_HX->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., lth_HX_best );
  PDF_cth_HX->Scale( DATA_cth_HX->Integral() * PDF_cth_HX->GetNbinsX() / ( PDF_cth_HX->Integral()*DATA_cth_HX->GetNbinsX() ) );

  TH1D* PDF_cth_HX_p1 = (TH1D*)PDF_1_vs_cth_HX->Clone();
  PDF_cth_HX_p1->SetName("PDF_cth_HX_p1");
  PDF_cth_HX_p1->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., 1. );
  PDF_cth_HX_p1->Scale( DATA_cth_HX->Integral() * PDF_cth_HX_p1->GetNbinsX() / ( PDF_cth_HX_p1->Integral()*DATA_cth_HX->GetNbinsX() ) );

  TH1D* PDF_cth_HX_m1 = (TH1D*)PDF_1_vs_cth_HX->Clone();
  PDF_cth_HX_m1->SetName("PDF_cth_HX_m1");
  PDF_cth_HX_m1->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., -1. );
  PDF_cth_HX_m1->Scale( DATA_cth_HX->Integral() * PDF_cth_HX_m1->GetNbinsX() / ( PDF_cth_HX_m1->Integral()*DATA_cth_HX->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_cth_HX->GetMaximum();
  if ( plotborder * PDF_cth_HX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_HX_p1->GetMaximum();
  if ( plotborder * PDF_cth_HX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_HX_m1->GetMaximum();

  PDF_cth_HX->GetXaxis()->SetTitle("cos#vartheta_{HX}");
  PDF_cth_HX->GetXaxis()->SetLabelOffset(0.028);
  PDF_cth_HX->GetXaxis()->SetTitleSize(0.05);
  PDF_cth_HX->GetXaxis()->SetTickLength(-0.03);
  PDF_cth_HX->GetXaxis()->SetTitleOffset(1.20);
  PDF_cth_HX->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_cth_HX->GetYaxis()->SetLabelOffset(0.032);
  PDF_cth_HX->GetYaxis()->SetTitleSize(0.05);
  PDF_cth_HX->GetYaxis()->SetTickLength(-0.03);
  PDF_cth_HX->GetYaxis()->SetTitleOffset(1.55);
  PDF_cth_HX->SetMinimum(0.);
  PDF_cth_HX->SetLineColor(kRed);
  PDF_cth_HX->SetLineStyle(1);
  PDF_cth_HX->SetLineWidth(2);
  PDF_cth_HX->SetMaximum(plotMax);
  PDF_cth_HX->Draw("L");

  PDF_cth_HX_p1->SetLineColor(kRed);
  PDF_cth_HX_p1->SetLineStyle(2);
  PDF_cth_HX_p1->Draw("L same");

  PDF_cth_HX_m1->SetLineColor(kRed);
  PDF_cth_HX_m1->SetLineStyle(10);
  PDF_cth_HX_m1->Draw("L same");

  DATA_cth_HX->SetLineColor(kBlack);
  DATA_cth_HX->Draw("E same");

  sprintf(filename,"%s/fit_HX_costh.eps",dirstruct);
  c1->Print( filename );


  // phi_HX

  TH1D* PDF_ph_HX = (TH1D*)PDF_1_vs_ph_HX->Clone();
  PDF_ph_HX->SetName("PDF_ph_HX");
  PDF_ph_HX->Add( PDF_1_vs_ph_HX, PDF_c2ph_vs_ph_HX, 1., 2.*lph_HX_best/(3.+lth_HX_best) );
  PDF_ph_HX->Scale( DATA_ph_HX->Integral() * PDF_ph_HX->GetNbinsX() / ( PDF_ph_HX->Integral()*DATA_ph_HX->GetNbinsX() ) );

  TH1D* PDF_ph_HX_p1 = (TH1D*)PDF_1_vs_ph_HX->Clone();
  PDF_ph_HX_p1->SetName("PDF_ph_HX_p1");
  PDF_ph_HX_p1->Add( PDF_1_vs_ph_HX, PDF_c2ph_vs_ph_HX, 1., 0.5 );
  PDF_ph_HX_p1->Scale( DATA_ph_HX->Integral() * PDF_ph_HX_p1->GetNbinsX() / ( PDF_ph_HX_p1->Integral()*DATA_ph_HX->GetNbinsX() ) );

  TH1D* PDF_ph_HX_m1 = (TH1D*)PDF_1_vs_ph_HX->Clone();
  PDF_ph_HX_m1->SetName("PDF_ph_HX_m1");
  PDF_ph_HX_m1->Add( PDF_1_vs_ph_HX, PDF_c2ph_vs_ph_HX, 1., -0.5 );
  PDF_ph_HX_m1->Scale( DATA_ph_HX->Integral() * PDF_ph_HX_m1->GetNbinsX() / ( PDF_ph_HX_m1->Integral()*DATA_ph_HX->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_ph_HX->GetMaximum();
  if ( plotborder * PDF_ph_HX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_HX_p1->GetMaximum();
  if ( plotborder * PDF_ph_HX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_HX_m1->GetMaximum();

  PDF_ph_HX->GetXaxis()->SetTitle("#varphi_{HX}");
  PDF_ph_HX->GetXaxis()->SetLabelOffset(0.028);
  PDF_ph_HX->GetXaxis()->SetTitleSize(0.05);
  PDF_ph_HX->GetXaxis()->SetTickLength(-0.03);
  PDF_ph_HX->GetXaxis()->SetTitleOffset(1.20);
  PDF_ph_HX->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_ph_HX->GetYaxis()->SetLabelOffset(0.032);
  PDF_ph_HX->GetYaxis()->SetTitleSize(0.05);
  PDF_ph_HX->GetYaxis()->SetTickLength(-0.03);
  PDF_ph_HX->GetYaxis()->SetTitleOffset(1.55);
  PDF_ph_HX->SetMinimum(0.);
  PDF_ph_HX->SetLineColor(kRed);
  PDF_ph_HX->SetLineStyle(1);
  PDF_ph_HX->SetLineWidth(2);
  PDF_ph_HX->SetMaximum(plotMax);
  PDF_ph_HX->Draw("L");

  PDF_ph_HX_p1->SetLineColor(kRed);
  PDF_ph_HX_p1->SetLineStyle(2);
  PDF_ph_HX_p1->Draw("L same");

  PDF_ph_HX_m1->SetLineColor(kRed);
  PDF_ph_HX_m1->SetLineStyle(10);
  PDF_ph_HX_m1->Draw("L same");

  DATA_ph_HX->SetLineColor(kBlack);
  DATA_ph_HX->Draw("E same");

  sprintf(filename,"%s/fit_HX_phi.eps",dirstruct);
  c1->Print( filename );


  // phith_HX

  TH1D* PDF_phth_HX = (TH1D*)PDF_1_vs_phth_HX->Clone();
  PDF_phth_HX->SetName("PDF_phth_HX");
  PDF_phth_HX->Add( PDF_1_vs_phth_HX, PDF_cphth_vs_phth_HX, 1., sqrt(2.)*ltp_HX_best/(3.+lth_HX_best) );
  PDF_phth_HX->Scale( DATA_phth_HX->Integral() * PDF_phth_HX->GetNbinsX() / ( PDF_phth_HX->Integral()*DATA_phth_HX->GetNbinsX() ) );

  TH1D* PDF_phth_HX_p1 = (TH1D*)PDF_1_vs_phth_HX->Clone();
  PDF_phth_HX_p1->SetName("PDF_phth_HX_p1");
  PDF_phth_HX_p1->Add( PDF_1_vs_phth_HX, PDF_cphth_vs_phth_HX, 1., sqrt(2.)/4. );
  PDF_phth_HX_p1->Scale( DATA_phth_HX->Integral() * PDF_phth_HX_p1->GetNbinsX() / ( PDF_phth_HX_p1->Integral()*DATA_phth_HX->GetNbinsX() ) );

  TH1D* PDF_phth_HX_m1 = (TH1D*)PDF_1_vs_phth_HX->Clone();
  PDF_phth_HX_m1->SetName("PDF_phth_HX_m1");
  PDF_phth_HX_m1->Add( PDF_1_vs_phth_HX, PDF_cphth_vs_phth_HX, 1., -sqrt(2.)/4. );
  PDF_phth_HX_m1->Scale( DATA_phth_HX->Integral() * PDF_phth_HX_m1->GetNbinsX() / ( PDF_phth_HX_m1->Integral()*DATA_phth_HX->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_phth_HX->GetMaximum();
  if ( plotborder * PDF_phth_HX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_HX_p1->GetMaximum();
  if ( plotborder * PDF_phth_HX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_HX_m1->GetMaximum();

  PDF_phth_HX->GetXaxis()->SetTitle("#tilde{#varphi}_{HX}");
  PDF_phth_HX->GetXaxis()->SetLabelOffset(0.028);
  PDF_phth_HX->GetXaxis()->SetTitleSize(0.05);
  PDF_phth_HX->GetXaxis()->SetTickLength(-0.03);
  PDF_phth_HX->GetXaxis()->SetTitleOffset(1.20);
  PDF_phth_HX->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_phth_HX->GetYaxis()->SetLabelOffset(0.032);
  PDF_phth_HX->GetYaxis()->SetTitleSize(0.05);
  PDF_phth_HX->GetYaxis()->SetTickLength(-0.03);
  PDF_phth_HX->GetYaxis()->SetTitleOffset(1.55);
  PDF_phth_HX->SetMinimum(0.);
  PDF_phth_HX->SetLineColor(kRed);
  PDF_phth_HX->SetLineStyle(1);
  PDF_phth_HX->SetLineWidth(2);
  PDF_phth_HX->SetMaximum(plotMax);
  PDF_phth_HX->Draw("L");

  PDF_phth_HX_p1->SetLineColor(kRed);
  PDF_phth_HX_p1->SetLineStyle(2);
  PDF_phth_HX_p1->Draw("L same");

  PDF_phth_HX_m1->SetLineColor(kRed);
  PDF_phth_HX_m1->SetLineStyle(10);
  PDF_phth_HX_m1->Draw("L same");

  DATA_phth_HX->SetLineColor(kBlack);
  DATA_phth_HX->Draw("E same");

  sprintf(filename,"%s/fit_HX_phith.eps",dirstruct);
  c1->Print( filename );



  // costh_PX

  TH1D* PDF_cth_PX = (TH1D*)PDF_1_vs_cth_PX->Clone();
  PDF_cth_PX->SetName("PDF_cth_PX");
  PDF_cth_PX->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., lth_PX_best );
  PDF_cth_PX->Scale( DATA_cth_PX->Integral() * PDF_cth_PX->GetNbinsX() / ( PDF_cth_PX->Integral()*DATA_cth_PX->GetNbinsX() ) );

  TH1D* PDF_cth_PX_p1 = (TH1D*)PDF_1_vs_cth_PX->Clone();
  PDF_cth_PX_p1->SetName("PDF_cth_PX_p1");
  PDF_cth_PX_p1->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., 1. );
  PDF_cth_PX_p1->Scale( DATA_cth_PX->Integral() * PDF_cth_PX_p1->GetNbinsX() / ( PDF_cth_PX_p1->Integral()*DATA_cth_PX->GetNbinsX() ) );

  TH1D* PDF_cth_PX_m1 = (TH1D*)PDF_1_vs_cth_PX->Clone();
  PDF_cth_PX_m1->SetName("PDF_cth_PX_m1");
  PDF_cth_PX_m1->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., -1. );
  PDF_cth_PX_m1->Scale( DATA_cth_PX->Integral() * PDF_cth_PX_m1->GetNbinsX() / ( PDF_cth_PX_m1->Integral()*DATA_cth_PX->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_cth_PX->GetMaximum();
  if ( plotborder * PDF_cth_PX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_PX_p1->GetMaximum();
  if ( plotborder * PDF_cth_PX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_PX_m1->GetMaximum();

  PDF_cth_PX->GetXaxis()->SetTitle("cos#vartheta_{PX}");
  PDF_cth_PX->GetXaxis()->SetLabelOffset(0.028);
  PDF_cth_PX->GetXaxis()->SetTitleSize(0.05);
  PDF_cth_PX->GetXaxis()->SetTickLength(-0.03);
  PDF_cth_PX->GetXaxis()->SetTitleOffset(1.20);
  PDF_cth_PX->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_cth_PX->GetYaxis()->SetLabelOffset(0.032);
  PDF_cth_PX->GetYaxis()->SetTitleSize(0.05);
  PDF_cth_PX->GetYaxis()->SetTickLength(-0.03);
  PDF_cth_PX->GetYaxis()->SetTitleOffset(1.55);
  PDF_cth_PX->SetMinimum(0.);
  PDF_cth_PX->SetLineColor(kRed);
  PDF_cth_PX->SetLineStyle(1);
  PDF_cth_PX->SetLineWidth(2);
  PDF_cth_PX->SetMaximum(plotMax);
  PDF_cth_PX->Draw("L");

  PDF_cth_PX_p1->SetLineColor(kRed);
  PDF_cth_PX_p1->SetLineStyle(2);
  PDF_cth_PX_p1->Draw("L same");

  PDF_cth_PX_m1->SetLineColor(kRed);
  PDF_cth_PX_m1->SetLineStyle(10);
  PDF_cth_PX_m1->Draw("L same");

  DATA_cth_PX->SetLineColor(kBlack);
  DATA_cth_PX->Draw("E same");

  sprintf(filename,"%s/fit_PX_costh.eps",dirstruct);
  c1->Print( filename );


  // phi_PX

  TH1D* PDF_ph_PX = (TH1D*)PDF_1_vs_ph_PX->Clone();
  PDF_ph_PX->SetName("PDF_ph_PX");
  PDF_ph_PX->Add( PDF_1_vs_ph_PX, PDF_c2ph_vs_ph_PX, 1., 2.*lph_PX_best/(3.+lth_PX_best) );
  PDF_ph_PX->Scale( DATA_ph_PX->Integral() * PDF_ph_PX->GetNbinsX() / ( PDF_ph_PX->Integral()*DATA_ph_PX->GetNbinsX() ) );

  TH1D* PDF_ph_PX_p1 = (TH1D*)PDF_1_vs_ph_PX->Clone();
  PDF_ph_PX_p1->SetName("PDF_ph_PX_p1");
  PDF_ph_PX_p1->Add( PDF_1_vs_ph_PX, PDF_c2ph_vs_ph_PX, 1., 0.5 );
  PDF_ph_PX_p1->Scale( DATA_ph_PX->Integral() * PDF_ph_PX_p1->GetNbinsX() / ( PDF_ph_PX_p1->Integral()*DATA_ph_PX->GetNbinsX() ) );

  TH1D* PDF_ph_PX_m1 = (TH1D*)PDF_1_vs_ph_PX->Clone();
  PDF_ph_PX_m1->SetName("PDF_ph_PX_m1");
  PDF_ph_PX_m1->Add( PDF_1_vs_ph_PX, PDF_c2ph_vs_ph_PX, 1., -0.5 );
  PDF_ph_PX_m1->Scale( DATA_ph_PX->Integral() * PDF_ph_PX_m1->GetNbinsX() / ( PDF_ph_PX_m1->Integral()*DATA_ph_PX->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_ph_PX->GetMaximum();
  if ( plotborder * PDF_ph_PX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_PX_p1->GetMaximum();
  if ( plotborder * PDF_ph_PX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_PX_m1->GetMaximum();

  PDF_ph_PX->GetXaxis()->SetTitle("#varphi_{PX}");
  PDF_ph_PX->GetXaxis()->SetLabelOffset(0.028);
  PDF_ph_PX->GetXaxis()->SetTitleSize(0.05);
  PDF_ph_PX->GetXaxis()->SetTickLength(-0.03);
  PDF_ph_PX->GetXaxis()->SetTitleOffset(1.20);
  PDF_ph_PX->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_ph_PX->GetYaxis()->SetLabelOffset(0.032);
  PDF_ph_PX->GetYaxis()->SetTitleSize(0.05);
  PDF_ph_PX->GetYaxis()->SetTickLength(-0.03);
  PDF_ph_PX->GetYaxis()->SetTitleOffset(1.55);
  PDF_ph_PX->SetMinimum(0.);
  PDF_ph_PX->SetLineColor(kRed);
  PDF_ph_PX->SetLineStyle(1);
  PDF_ph_PX->SetLineWidth(2);
  PDF_ph_PX->SetMaximum(plotMax);
  PDF_ph_PX->Draw("L");

  PDF_ph_PX_p1->SetLineColor(kRed);
  PDF_ph_PX_p1->SetLineStyle(2);
  PDF_ph_PX_p1->Draw("L same");

  PDF_ph_PX_m1->SetLineColor(kRed);
  PDF_ph_PX_m1->SetLineStyle(10);
  PDF_ph_PX_m1->Draw("L same");

  DATA_ph_PX->SetLineColor(kBlack);
  DATA_ph_PX->Draw("E same");

  sprintf(filename,"%s/fit_PX_phi.eps",dirstruct);
  c1->Print( filename );


  // phith_PX

  TH1D* PDF_phth_PX = (TH1D*)PDF_1_vs_phth_PX->Clone();
  PDF_phth_PX->SetName("PDF_phth_PX");
  PDF_phth_PX->Add( PDF_1_vs_phth_PX, PDF_cphth_vs_phth_PX, 1., sqrt(2.)*ltp_PX_best/(3.+lth_PX_best) );
  PDF_phth_PX->Scale( DATA_phth_PX->Integral() * PDF_phth_PX->GetNbinsX() / ( PDF_phth_PX->Integral()*DATA_phth_PX->GetNbinsX() ) );

  TH1D* PDF_phth_PX_p1 = (TH1D*)PDF_1_vs_phth_PX->Clone();
  PDF_phth_PX_p1->SetName("PDF_phth_PX_p1");
  PDF_phth_PX_p1->Add( PDF_1_vs_phth_PX, PDF_cphth_vs_phth_PX, 1., sqrt(2.)/4. );
  PDF_phth_PX_p1->Scale( DATA_phth_PX->Integral() * PDF_phth_PX_p1->GetNbinsX() / ( PDF_phth_PX_p1->Integral()*DATA_phth_PX->GetNbinsX() ) );

  TH1D* PDF_phth_PX_m1 = (TH1D*)PDF_1_vs_phth_PX->Clone();
  PDF_phth_PX_m1->SetName("PDF_phth_PX_m1");
  PDF_phth_PX_m1->Add( PDF_1_vs_phth_PX, PDF_cphth_vs_phth_PX, 1., -sqrt(2.)/4. );
  PDF_phth_PX_m1->Scale( DATA_phth_PX->Integral() * PDF_phth_PX_m1->GetNbinsX() / ( PDF_phth_PX_m1->Integral()*DATA_phth_PX->GetNbinsX() ) );

  plotborder = 1.1; // how much space to leave beyond the maximum of the plots

  plotMax = plotborder * PDF_phth_PX->GetMaximum();
  if ( plotborder * PDF_phth_PX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_PX_p1->GetMaximum();
  if ( plotborder * PDF_phth_PX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_PX_m1->GetMaximum();

  PDF_phth_PX->GetXaxis()->SetTitle("#tilde{#varphi}_{PX}");
  PDF_phth_PX->GetXaxis()->SetLabelOffset(0.028);
  PDF_phth_PX->GetXaxis()->SetTitleSize(0.05);
  PDF_phth_PX->GetXaxis()->SetTickLength(-0.03);
  PDF_phth_PX->GetXaxis()->SetTitleOffset(1.20);
  PDF_phth_PX->GetYaxis()->SetTitle("event PDF [a.u.]");
  PDF_phth_PX->GetYaxis()->SetLabelOffset(0.032);
  PDF_phth_PX->GetYaxis()->SetTitleSize(0.05);
  PDF_phth_PX->GetYaxis()->SetTickLength(-0.03);
  PDF_phth_PX->GetYaxis()->SetTitleOffset(1.55);
  PDF_phth_PX->SetMinimum(0.);
  PDF_phth_PX->SetLineColor(kRed);
  PDF_phth_PX->SetLineStyle(1);
  PDF_phth_PX->SetLineWidth(2);
  PDF_phth_PX->SetMaximum(plotMax);
  PDF_phth_PX->Draw("L");

  PDF_phth_PX_p1->SetLineColor(kRed);
  PDF_phth_PX_p1->SetLineStyle(2);
  PDF_phth_PX_p1->Draw("L same");

  PDF_phth_PX_m1->SetLineColor(kRed);
  PDF_phth_PX_m1->SetLineStyle(10);
  PDF_phth_PX_m1->Draw("L same");

  DATA_phth_PX->SetLineColor(kBlack);
  DATA_phth_PX->Draw("E same");

  sprintf(filename,"%s/fit_PX_phith.eps",dirstruct);
  c1->Print( filename );



}
