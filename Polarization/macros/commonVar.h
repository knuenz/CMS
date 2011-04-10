#include "TLorentzVector.h"
#include "TMath.h"

namespace jpsi{

  //enum {L1DoubleMuOpen, Mu0Track0Jpsi, Mu3Track0Jpsi, DoubleMu0, DoubleMu3};

  // beam energy in GeV
  const double pbeam = 3500.;
  // masses
  const double Mprot = 0.9382720;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );
  const double muMass = 0.105658;
  //rap bins
  Int_t const kNbRapForPTBins = 5;
  /* Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.5, 1.9, 2.3}; */
  Double_t rapForPTRange[kNbRapForPTBins+1] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
  //pT bins (optimized to have at least 10.000 entries per bin)
  Int_t const kNbPTMaxBins = 12;
  Int_t const kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, 7,8,9,12,12};//all y, y1, y2, y3, y4, y5
  Double_t pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},//all rapidities
    {0., 6., 7., 8., 10., 15., 20., 30.},//mid-rap
    {0., 4., 6., 7., 8., 10., 15., 20., 30.},
    {0., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.},
    {0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 15., 20., 30.}};//most forward


  //the following values are dummy and need to be filled out at end of analysis
  Double_t pTWCentre_rap[kNbRapForPTBins+1][23] =
    {{0.5, 1.25, 1.75, 1.95, 2.5, 2.8, 3.1, 3.4, 3.9, 4.3, 4.6, 5.0, 5.5, 6.0, 6.5, 7.2, 7.8, 8.5, 9.6, 11.0, 14., 27.},
     //{7., 12.},
     //{4.5, 6.0, 7.3, 15.},
    		{3,6.5,7.5,9,12.5,17.5,25},
    					{2,5,6.5,7.5,9,12.5,17.5,25},
    {1.0, 2.5, 3.5, 4.5, 5.5, 6.7, 10.},
     {0.5, 1.5, 2.5, 3.5, 5.0, 13.},
     {0.5, 1.5, 2.5, 3.5, 5.0, 13.}};
  //need to be extracted from the real data:
  /* Double_t pTWCentre[kNbPTBins] = {1.527, 3.561, 6.021, 9.540, 15.205, 23.566}; */
  //number of reference frames
  Int_t const kNbFrames = 6;
  Char_t *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  //polarization variables (3rd one was for testing purposes)
  Int_t const kNbPolVar = 3; //cosTheta, phi, cos2Phi
  enum {cosThPol,phiPol,cos2PhiPol};
  //cosTheta
  Int_t const kNbBinsCosT = 40;
  Double_t cosTMin = -1., cosTMax = 1.;
  //phi for pol.
  Int_t const kNbBinsPhiPol = 36;
  //  Double_t phiPolMin = 0., phiPolMax = 360.;
  Double_t phiPolMin = -180., phiPolMax = 180.;
  //cos2Phi
  Int_t const kNbBinsCos2Phi = 40;
  Double_t cos2PhiMin = -1., cos2PhiMax = 1.;

  //study the negative and positive rapidity sides separately
  Int_t const kNbRapBins = kNbRapForPTBins;
  /* Double_t rapRange[2*kNbRapBins+1] = {-2.3, -1.9, -1.5, -0.9, 0., 0.9, 1.5, 1.9, 2.3}; */
  Double_t rapRange[2*kNbRapBins+1] = {-2.4, -2.1, -1.6, -1.2, -0.9, 0., 0.9, 1.2, 1.6, 2.1, 2.4};

  //phase space limiting cuts:
  Int_t const kNbEtaRegions = 3;
  Double_t etaPS[kNbEtaRegions] = {1.3, 2.2, 2.4}; //pseudo-rap cuts for muons
  Double_t pTMuMin[kNbEtaRegions] = {3.3, 0., 0.8};
  Double_t pMuMin[kNbEtaRegions] = {0., 2.9, 0.};
  Double_t rapYPS = 2.4;
  /* Double_t JpsiCtauMax = 0.100; //100 micron */
  Double_t JpsiCtauMax = 1000.; //effectively no cut on lifetime
  Double_t nSigMass = 2.;
  Double_t polMassJpsi[kNbRapForPTBins+1] = {3.092, 3.094, 3.094, 3.092, 3.092, 3.090};//[all rap, rap bin 1-3]
  Double_t sigmaMassJpsi[kNbRapForPTBins+1] = {0.042, 0.024, 0.035, 0.040, 0.045, 0.056};//[all rap, rap bin 1-3]
  int nSigBkgLow=2;
  int nSigBkgHigh=2;

  //some make up to use the same colour and marker for each pT and rapidity bin
  //in every plotting macro:
  Int_t colour_pT[kNbPTMaxBins+1] = {1, 2, 3, 4, 6, 7, 8, 49, 38, 46, 12, 40};
  Int_t marker_pT[kNbPTMaxBins+1] = {20, 21, 25, 22, 23, 26, 27, 28, 29, 30, 20, 20};
  /* Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 6, 7, 7, 6, 4, 3, 2}; */
  /* Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 29, 20, 24, 30, 23, 25, 24}; */
  /* Int_t colour_rap[2*kNbRapBins+1] = {1, 2, 3, 4, 4, 3, 2}; */
  /* Int_t marker_rap[2*kNbRapBins+1] = {20, 20, 21, 22, 23, 25, 24}; */

  Int_t colour_rapForPTBins[kNbRapForPTBins+1] = {0,kRed,kBlue};//{1, 30, 4, 2, 3, kMagenta+1};
  Int_t marker_rapForPTBins[kNbRapForPTBins+1] = {20, 20,21};//, 20, 25, 20, 22};

  //min number of entries to consider the bin
  //(in acceptance map, etc)
  Int_t minEntriesPerBin = 10;

  //////// number of events to generate for ToyMC study///////////////
  int numEventsBin [5][12]={
		  {8175, 10774, 25593, 58835, 65923, 16447, 6742},
		  {354,28142,25301,23313,31274,25281, 5467, 2179},
		  {50413, 52517, 64535, 58920, 43601, 50493, 37124,  7160,  2567},
		  {39791, 144638, 175916, 183508, 152802, 117774,  85012,  57282,  62223,  43722,   7880,   2815},
		  {29341,   78976,   65594,   51442,   37421,   26644,   18135,   12324,   13998,   10598,    1973,     681}};

  int numEventsBin_HLT_Mu0_TkMu0_OST_Jpsi [5][12]={
		  {4465,5884,14040,31820,35577,8833,3641},
		  {209,15803,14018,12942,17118,13759,3001,1158}
  };


  int numEventsPR [5][12]={
		  {0,0,0,0,0,5842,2178},
		  {0,0,0,0,14347,10733,2101,0}
  };

  int numEventsNP [5][12]={
		  {0,0,0,0,0,3572,1712},
		  {0,0,0,0,4181,4202,1182,0}
  };


  double inj_PR_CS_th [5][12]={
		  {0,0,0,0,0,-0.00922041,0.0189234},
		  {0,0,0,0,-0.0826305,-0.103903,-0.00231016,0}
  };
  double inj_PR_CS_ph [5][12]={
		  {0,0,0,0,0,0.0020492,0.10338},
		  {0,0,0,0,0.0938028,0.165187,0.156993,0}
  };
  double inj_PR_CS_thph [5][12]={
		  {0,0,0,0,0,-0.02629752,0.0402561},
		  {0,0,0,0,-0.0154505,-0.0748545,-0.0350425,0}
  };
  double inj_PR_HX_th [5][12]={
		  {0,0,0,0,0,-0.133236,-0.0504058},
		  {0,0,0,0,0.380172,0.388609,0.273877,0}
  };
  double inj_PR_HX_ph [5][12]={
		  {0,0,0,0,0,0.00963261,0.0395828},
		  {0,0,0,0,0.0260826,0.0210654,0.0741257,0}
  };
  double inj_PR_HX_thph [5][12]={
		  {0,0,0,0,0,0.01813,-0.0375903},
		  {0,0,0,0,0.0438046,0.0207189,0.0231536,0}
  };


  double inj_NP_CS_th [5][12]={
		  {0,0,0,0,0,0.141544,0.08739241},
		  {0,0,0,0,-0.256359,0.264646,-0.17409,0}
  };
  double inj_NP_CS_ph [5][12]={
		  {0,0,0,0,0,-0.0405901,0.0192388},
		  {0,0,0,0,0.132394,-0.146828,0.125612,0}
  };
  double inj_NP_CS_thph [5][12]={
		  {0,0,0,0,0,0.0367851,-0.0084211},
		  {0,0,0,0,-0.0672116,0.00945514,-0.0514772,0}
  };
  double inj_NP_HX_th [5][12]={
		  {0,0,0,0,0,-0.228897,-0.242696},
		  {0,0,0,0,0.548767,-0.290254,0.298969,0}
  };
  double inj_NP_HX_ph [5][12]={
		  {0,0,0,0,0,0.0617815,0.0248167},
		  {0,0,0,0,-0.0649005,0.0392934,-0.043394,0}
  };
  double inj_NP_HX_thph [5][12]={
		  {0,0,0,0,0,-0.0345403,0.0102396},
		  {0,0,0,0,0.0173524,0.0528877,0.0201468,0}
  };
/*
  numEventsBin_HLT_Mu0_TkMu0_OST_Jpsi Lindsey:
  rap =  1  pt =  6
  N_Events (l < .1 mm):  5842.0
  N_Events (l > .1 mm):  3572.0

  rap =  1  pt =  7
  N_Events (l < .1 mm):  2178.0
  N_Events (l > .1 mm):  1712.0

  rap =  2  pt =  5
  N_Events (l < .1 mm):  14347.0
  N_Events (l > .1 mm):  4181.0

  rap =  2  pt =  6
  N_Events (l < .1 mm):  10733.0
  N_Events (l > .1 mm):  4202.0

  rap =  2  pt =  7
  N_Events (l < .1 mm):  2101.0
  N_Events (l > .1 mm):  1182.0

  rap =  2  pt =  8
  N_Events (l < .1 mm):  728.0
  N_Events (l > .1 mm):  542.0
  */

// MC Fall10 Prompt numEvents in bins (noTrigger):
// from pt2 onwards: rap1: 21772,40993,83369,89817,20855,8551
// from pt2 onwards: rap1: 59939,38192,32495,42338,35506,7273,2778

  int numEventsBin_DoubleMu0 [2][8]={
		  {0,0,0,0,0,18837,8175},
		  {0,0,0,0,0,25713,6352,2596}
  };

  double Bfrac [2][8]={
		  {0.000000,0.000000,0.000000,0.000000,0.000000,0.486931,0.569504},
		  {0.000000,0.000000,0.000000,0.000000,0.000000,0.362768,0.471131,0.575429}
  };


  double scenarioLambdas [5][6]={//5 scenarios, parameters in order: lam_th_PR, lam_ph_PR,lam_thph_PR, lam_th_NP, lam_ph_NP,lam_thph_NP
  		{+1,0,0,-1,0,0},//{+0.5,0,0,-0.5,0,0},//REAL SCENARIO 1:
  		{0,0.5,0,0,-0.5,0},
  		{0,0,0,0,0,0},
  		{0.5,0,0,-0.5,0,0},// Original {0,-0.5,0,0,0.5,0},
  		{-0.5,0,0,0.5,0,0}//Original {-1,0,0,+1,0,0}
  };

  double fPinP_actual[2][8]={{0.000000,0.000000,0.000000,0.000000,0.000000,0.820462,0.772049},{0.000000,0.000000,0.000000,0.000000,0.000000,0.885443,0.839310,0.768876}};
  double fNPinP_actual[2][8]={{0.000000,0.000000,0.000000,0.000000,0.000000,0.179538,0.227951},{0.000000,0.000000,0.000000,0.000000,0.000000,0.114557,0.160690,0.231124}};
  double fPinNP_actual[2][8]={{0.000000,0.000000,0.000000,0.000000,0.000000,0.000144,0.000052},{0.000000,0.000000,0.000000,0.000000,0.000000,0.002540,0.000212,0.000039}};
  double fNPinNP_actual[2][8]={{0.000000,0.000000,0.000000,0.000000,0.000000,0.999856,0.999948},{0.000000,0.000000,0.000000,0.000000,0.000000,0.997460,0.999788,0.999961}};
  double fracPRregion[2][8]={{0.000000,0.000000,0.000000,0.000000,0.000000,0.625275,0.557571},{0.000000,0.000000,0.000000,0.000000,0.000000,0.718869,0.630030,0.552174}};

 /* 	ltheta_P   lphi_P     ltheta_NP   lphi_NP
  		+1         0           -1         0
  		 0         +0.5          0         -0.5
  		 0          0            0         0
  		 0         -0.5          0         +0.5
  		 -1         0           +1         0
*/

  char scratchLocation[200] = "/scratch/knuenz/Polarization/RootInput/ProjectClosure/";
  char WorkspaceDir[200] = "/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/";
}
