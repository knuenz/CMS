#include "TLorentzVector.h"
#include "TMath.h"

namespace onia{

  // beam energy in GeV
  const double pbeam = 3500.;
  // masses
  const double Mprot = 0.9382720;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );
  const double muMass = 0.105658;
  Int_t const kNbRapForPTBins = 2;
  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.6, 1.2}; //2 September 2011
//  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.3, 0.6}; //2 September 2011
  //study the negative and positive rapidity sides separately
  Int_t const kNbRapBins = kNbRapForPTBins;
  Double_t rapRange[2*kNbRapBins+1] = {-1.2, -0.6, -0., 0.6, 1.2};
//  Double_t rapRange[2*kNbRapBins+1] = {-0.6, -0.3, -0., 0.3, 0.6};

  Int_t const kNbPTMaxBins = 10;
/*  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 12, 12};//all y, y1
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {5., 6., 7., 8., 9., 10., 11., 12., 14., 16., 20., 25., 30., 50., 100.},//all rapidities
    {5., 6., 7., 8., 9., 10., 11., 12., 14., 16., 20., 30., 50.},//mid-rap
    {5., 6., 7., 8., 9., 10., 11., 12., 14., 16., 20., 30., 50.}};//most forward
*/

  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 10, 10};//all y, y1
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
	{5., 6., 7., 8., 9., 10., 12., 16., 20., 30., 50.},//all rapidities
    {5., 6., 7., 8., 9., 10., 12., 16., 20., 30., 50.},//mid-rap
    {5., 6., 7., 8., 9., 10., 12., 16., 20., 30., 50.}};//most forward
  // High ctau
/*    Int_t const kNbPTBins[kNbRapForPTBins+1] = {3, 3, 3};//all y, y1
  	  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
  			  {10., 13., 18., 30.},//mid-rap
  			  {10., 13., 18., 30.},//mid-rap
  			  {10., 13., 18., 30.}};//mid-rap
*/

// SingleMuPD binning
/*  Int_t const kNbPTBins[kNbRapForPTBins+1] = {2, 2, 2};//all y, y1
	  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
			  {20., 30., 50.},//all rapidities
			  {20., 30., 50.},//mid-rap
			  {20., 30., 50.}};//most forward
*/

  //number of reference frames
  Int_t const kNbFrames = 6;
  Char_t *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  //polarization variables
  Int_t const kNbPolVar = 2; //cosTheta, phi
  enum {cosThPol,phiPol};
  //cosTheta
  Int_t const kNbBinsCosT = 15;
  Double_t cosTMin = -1., cosTMax = 1.;
  //phi for pol. 
   Int_t const kNbBinsPhiPol = 15;
  Double_t phiPolMin = -180., phiPolMax = 180.;
  Double_t massMinL = 8.6; //left mass edge of the BG mass window
  Double_t massMaxR = 11.4; //right mass edge of the BG mass window
  Double_t nSigmaL = 4.; //nSigma to define the right edge of the L BG mass window
  Double_t nSigmaR = 3.; //nSigma to define the left edge of the R BG mass window

  //phase space limiting cuts:
  Double_t rapYPS = 1.2;
  Double_t JpsiCtauMax = 1000.; //effectively no cut on lifetime
  Double_t JpsiCtauMin = -1000; //effectively no cut on lifetime
  Double_t JpsiCtauSigMin=-1000; //-1000
  Double_t JpsiCtauSigMax=2.; //1000//2.
  /* Double_t nSigMass = 1000.; */
  Double_t nSigMass = 1.;
  /* Double_t polMassJpsi[kNbRapForPTBins+1] = {3.094, 3.094, 3.094, 3.092, 3.092, 3.092};//[all rap, rap bin 1-5] */
  /* Double_t sigmaMassJpsi[kNbRapForPTBins+1] = {0.037, 0.022, 0.029, 0.039, 0.045, 0.059};//[all rap, rap bin 1-5] */
  Double_t polMassOnia[kNbRapForPTBins+1] = {9.45, 9.45, 9.45};//[all rap, rap bin 1-n]
  Double_t sigmaMassOnia[kNbRapForPTBins+1] = {0.080, 0.065, 0.085};//[all rap, rap bin 1-n]

  //some make up to use the same colour and marker for each pT and rapidity bin
  //in every plotting macro:
//  Int_t colour_pT[kNbPTMaxBins+1] = {1, 2, 3, 4, 6, 7, 8, 49, 38, 46, 12, 40};
//  Int_t marker_pT[kNbPTMaxBins+1] = {20, 21, 25, 22, 23, 26, 27, 28, 29, 30, 20, 20};
  
  /* Int_t colour_rapForPTBins[kNbRapForPTBins+1] = {1, 30, 4, 2, 3, kMagenta+1}; */
  /* Int_t marker_rapForPTBins[kNbRapForPTBins+1] = {20, 21, 25, 20, 22, 29}; */
  Int_t colour_rapForPTBins[kNbRapForPTBins+1] = {1, 30, 4};
  Int_t marker_rapForPTBins[kNbRapForPTBins+1] = {20, 21, 25};

}
