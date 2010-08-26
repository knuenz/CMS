#include "TLorentzVector.h"
#include "TMath.h"

//enum {L1DoubleMuOpen, Mu0Track0Jpsi, Mu3Track0Jpsi, DoubleMu0, DoubleMu3};

// beam energy in GeV
const double pbeam = 3500.;
// masses
const double Mprot = 0.9382720;
const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );

//phase space limiting cuts:
Double_t etaPS = 2.4; //pseudo-rap cuts for muons
Double_t pTMuMin = 0.;
Double_t rapYPS = 2.3;

//pT bins
Int_t const kNbPTBins = 2;
Double_t pTRange[kNbPTBins+1] = {7.5,12.5,20};//{0., 2.5, 5.0, 7.5, 12.5, 20., 30.};
//need to be extracted from the real data:
Double_t pTWCentre[kNbPTBins] = {5,10};//{1.527, 3.561, 6.021, 9.540, 15.205, 23.566};
//rap bins
Int_t const kNbRapForPTBins = 1;
Double_t rapForPTRange[kNbRapForPTBins+1] = {0,0.9};//{0., 0.9, 1.5, 2.3};
Double_t pTWCentre_rap[kNbRapForPTBins][kNbPTBins] = {{5,10}};
//  {{1.529, 3.569, 5.991, 9.468, 15.206, 23.586},
//  {1.522, 3.561, 6.015, 9.536, 15.163, 23.618},
//   {1.533, 3.550, 6.085, 9.626, 15.267, 23.544}};
//number of reference frames
Int_t const kNbFrames = 2;
Char_t *frameLabel[kNbFrames] = {"CS", "HX"};
enum {CS,HX};
//polarization variables (3rd one was for testing purposes)
Int_t const kNbPolVar = 3; //cosTheta, phi, cos2Phi
enum {cosThPol,phiPol,cos2PhiPol};
//cosTheta 
Int_t const nBinsCosT = 40;
Double_t cosTMin = -1., cosTMax = 1.;
//phi for pol. 
Int_t const nBinsPhiPol = 36;
Double_t phiPolMin = 0., phiPolMax = 360.;
//cos2Phi
Int_t const nBinsCos2Phi = 40;
Double_t cos2PhiMin = -1., cos2PhiMax = 1.;

//study the negative and positive rapidity sides separately
Int_t const kNbRapBins = kNbRapForPTBins;
Double_t rapRange[2*kNbRapBins+1] = {-2.3,0,2.3};//{-2.3, -1.5, -0.9, 0., 0.9, 1.5, 2.3};

//some make up to use the same colour and marker for each pT and rapidity bin
//in every plotting macro:
Int_t colour_pT[kNbPTBins+1] = {1,2,3};//{1, 2, 3, 4, 6, 7, 8};
Int_t marker_pT[kNbPTBins+1] = {1,2,3};;//{20, 21, 22, 23, 20, 21, 29};
/* Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 6, 7, 7, 6, 4, 3, 2}; */
/* Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 29, 20, 24, 30, 23, 25, 24}; */
Int_t colour_rap[2*kNbRapBins+1] = {1,2,3};//{1, 2, 3, 4, 4, 3, 2};
Int_t marker_rap[2*kNbRapBins+1] = {1,2,3};//{20, 20, 21, 22, 23, 25, 24};
