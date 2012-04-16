/*
 * upsilon_MCMassFit.C
 *
 *  Created on: Jan 30, 2012
 *      Author: valentinknuenz
 */

#include "rootIncludes.inc"
#include "commonVar.h"
#include "CBFunction.C"
#include "TMath.h"
#include "THStack.h"

Double_t peak_min = 8.9, peak_max = 10.6; //veto this range when fitting the continuum

const Double_t massPDG1S = 9.460;
const Double_t massPDG2S = 10.023;
const Double_t massPDG3S = 10.355;

Int_t const kNbSpecies = 3;
Char_t *specName[kNbSpecies] = {"1S", "2S", "3S"};
enum {UPS1S, UPS2S, UPS3S, BG};
Int_t colour[kNbSpecies+1] = {kRed-9,kGreen-3,kBlue-8,kGray};
TH1F *hMass;
Double_t binWidth;
TF1* fRECO;
Double_t fitParBG[3];
Double_t intCB[kNbSpecies]; //integral values of a CB with N=1 and fixed alpha, n, sigma, width

Double_t massMin[kNbSpecies], massMax[kNbSpecies];
TF1 *fUps1S, *fUps2S, *fUps3S, *fBG;
Double_t fracBG[kNbSpecies];
Double_t nY[kNbSpecies];

void GetHisto(Char_t *fileNameIn, Int_t iRapBin, Int_t iPTBin);
void FitSignalBG(Double_t nSigma, Int_t iRapBin, Int_t iPTBin);
Double_t fitPolyCrystal3(Double_t *x, Double_t *par);
Double_t fitContinuum(Double_t *x, Double_t *par);
Double_t DrawContinuum(Double_t *x, Double_t *par);
void DrawFit(Double_t nSigma, Int_t iRapBin, Int_t iPTBin);
void SaveCBParameters(Int_t iRapBin, Int_t iPTBin, Double_t alpha, Double_t n, Double_t alphaErr, Double_t nErr);
void SaveFitPars(Int_t iRapBin, Int_t iPTBin);
//==============================
void upsilon_MCMassFit(Int_t iRapBin = 0,
		      Int_t iPTBin = 0,
		      Double_t nSigma = 2.,
		      Char_t *fileNameIn = "RootFiles/selEvents_data_Ups_2Aug2011.root"){

  GetHisto(fileNameIn, iRapBin, iPTBin);
  if(hMass->GetEntries() < 200.){
    printf("\n\n\nskip processing this bin, because the number of entries is smaller than 200\n\n\n");
    return;
  }
  FitSignalBG(nSigma, iRapBin, iPTBin);
  DrawFit(nSigma, iRapBin, iPTBin);
}

//===============================
void DrawFit(Double_t nSigma, Int_t iRapBin, Int_t iPTBin){

  Char_t name[100];
  //prepare the drawing of the individual components:
  fBG->SetFillColor(colour[BG]);
  fBG->SetLineColor(colour[BG]);
  fBG->SetFillStyle(1001);
  fBG->SetNpx(1000);
  TH1 *hBG = fBG->GetHistogram();

  fUps1S->SetNpx(1000);
  fUps1S->SetFillColor(colour[UPS1S]);
  fUps1S->SetLineColor(colour[UPS1S]);
  fUps1S->SetFillStyle(1001);

  fUps2S->SetNpx(1000);
  fUps2S->SetFillColor(colour[UPS2S]);
  fUps2S->SetLineColor(colour[UPS2S]);
  fUps2S->SetFillStyle(1001);

  fUps3S->SetNpx(1000);
  fUps3S->SetFillColor(colour[UPS3S]);
  fUps3S->SetLineColor(colour[UPS3S]);
  fUps3S->SetFillStyle(1001);

  TH1 *hUps1S = fUps1S->GetHistogram();
  TH1 *hUps2S = fUps2S->GetHistogram();
  TH1 *hUps3S = fUps3S->GetHistogram();

  THStack *hStack = new THStack("hMass_Stack", "");
  hStack->Add(hBG);
  hStack->Add(hUps3S);
  hStack->Add(hUps2S);
  hStack->Add(hUps1S);
  hStack->Draw("same");

  hMass->Draw("same");
  fRECO->Draw("same");

  TLine *line[3];
  Double_t max[3] = {0.5, 0.5, 0.3};
  for(int iL = 0; iL < 3; iL++){
    line[iL]= new TLine(massMin[iL], 0.1, massMin[iL], max[iL]*hUps1S->GetMaximum());
    line[iL]->SetLineStyle(2); line[iL]->SetLineColor(colour[iL]);
    line[iL]->SetLineWidth(2); line[iL]->Draw();
    line[iL]->DrawLine(massMax[iL], 0.1, massMax[iL], max[iL]*hUps1S->GetMaximum());
  }

  if(iRapBin == 0) sprintf(name, "|y| < %1.1f", onia::rapYPS);
  else if(iRapBin == 1) sprintf(name, "|y| < %1.1f", onia::rapForPTRange[iRapBin]);
  else if(iRapBin > 1)  sprintf(name, "%1.1f < |y| < %1.1f", onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin]);
  TLatex *tex = new TLatex(8.5, hStack->GetMaximum(), name);
  tex->SetTextSize(0.04); tex->Draw();
  if(iPTBin == 0) sprintf(name, "all p_{T}");
  //  else if(iPTBin == 1) sprintf(name, "p_{T} < %1.1f GeV", onia::pTRange[iRapBin][iPTBin]);
  else sprintf(name, "%1.1f < p_{T} < %1.1f", onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin]);
  tex->DrawLatex(8.5, 0.94*hStack->GetMaximum(), name);

  sprintf(name, "frac(BG) in #pm %1.1f#sigma:", nSigma);
  tex->DrawLatex(8.5, 0.86*hStack->GetMaximum(), name);
  sprintf(name, "%1.2f, %1.2f, %1.2f", fracBG[0], fracBG[1], fracBG[2]);
  tex->DrawLatex(8.5, 0.80*hStack->GetMaximum(), name);

  sprintf(name, "Figures/massFit_rap%d_pT%d.pdf", iRapBin, iPTBin);  gPad->Print(name);
}

//===============================
void FitSignalBG(Double_t nSigma, Int_t iRapBin, Int_t iPTBin){

  gStyle->SetOptFit(kTRUE);
  gStyle->SetOptStat(kFALSE);

  Char_t name[100];
  //1.) perform the fit to the continuum, using the sidebands
  sprintf(name, "c1_rap%d_pT%d", iRapBin, iPTBin);
  TCanvas *c1 = new TCanvas(name);
  sprintf(name, "Events in %1.0f MeV", 1000.*binWidth);
  hMass->SetYTitle(name);
  hMass->SetTitleOffset(1.7, "y");
  hMass->Draw();

  //starting values for fit:
  Double_t a = 1, b = 1, c = 1;
  Double_t range_min = 8.6, range_max = 11.4;
  Double_t normY1S = hMass->GetMaximum() / binWidth;
  Double_t normY2S = 0;
  Double_t normY3S = 0;
  Double_t sigma1S = 0.1, mean1S = 9.45;
  Double_t alpha = 1.33, n = 6.6; //CB-tail parameters

  printf("will be fitting the continuum between %1.2f < M < %1.2f\n", range_min, range_max);
  sprintf(name, "fCont");
  TF1* fCONT = new TF1(name, fitContinuum, range_min, range_max, 3);
//  fCONT->SetParameters(a, b, c);
//  hMass->Fit(fCONT, "0", "", range_min, range_max);
//  fCONT = hMass->GetFunction(name);

  Double_t chisqrd = fCONT->GetChisquare();
  Int_t NDF = fCONT->GetNDF();
  Double_t redChi2 = chisqrd/NDF;
  printf("\nChisqrd = %1.3f, NDF = %d, Chisqrd/NDF = %1.3f, Prob = %1.3f\n\n", chisqrd, NDF, redChi2, TMath::Prob(chisqrd,NDF));

  fCONT->GetParameters(fitParBG);

  //2.) fit the peaks on the top of a fixed continuum
  sprintf(name,"fPeaks");
  Int_t const npar = 10;
  fRECO = new TF1(name, fitPolyCrystal3, peak_min, peak_max, npar);
  fRECO->SetParNames("normY1S", "mass_Ups1S", "sigma_Ups1S", "normY2S", "normY3S", "n", "alpha", "a", "b", "c");
  fRECO->SetParameters(normY1S, mean1S, sigma1S, normY2S, normY3S, n, alpha, a, b, c);
  fRECO->FixParameter(7,a);
  fRECO->FixParameter(8,b);
  fRECO->FixParameter(9,c);

  fRECO->FixParameter(3,0);
  fRECO->FixParameter(4,0);

  //fix alpha and n from the fit to all bins
  if(iPTBin > 0 && iRapBin > 0){

    Char_t fileName[100];
    if(iRapBin == 0)
      sprintf(fileName, "tmpFiles/CBParameters.root");
    else
      sprintf(fileName, "tmpFiles/CBParameters_rap%d.root", iRapBin);
    TFile *fIn = new TFile(fileName);
    TTree *treeIn = (TTree *) gDirectory->Get("CBPars");
    Double_t alphaAll, nAll;
    TBranch *b_alphaAll, *b_nAll;
    treeIn->SetBranchAddress("alphaAll", &alphaAll, &b_alphaAll);
    treeIn->SetBranchAddress("nAll", &nAll, &b_nAll);
    Long64_t iEntry = treeIn->LoadTree(0);
    treeIn->GetEntry(0);
    printf("alpha and n from File: %1.3f, %1.3f\n", alphaAll, nAll);
    fRECO->FixParameter(5, nAll);
    fRECO->FixParameter(6, alphaAll);
  }

  hMass->Fit(fRECO, "0", "", peak_min, peak_max);

  fRECO = hMass->GetFunction(name);
  fRECO->SetLineWidth(1);

  Double_t fitParTot[npar];
  fRECO->GetParameters(fitParTot);
  normY1S = fitParTot[0];   printf("normY1S = %1.3e\n", normY1S);
  mean1S = fitParTot[1];
  sigma1S = fitParTot[2];
  normY2S = 0;
  normY3S = 0;
  n = fitParTot[5];  printf("n = %1.3f\n", n);
  alpha = fitParTot[6];  printf("alpha = %1.3f\n", alpha);

  chisqrd = fRECO->GetChisquare();
  NDF = fRECO->GetNDF();
  redChi2 = chisqrd/NDF;
  printf("\nChisqrd = %1.3f, NDF = %d, Chisqrd/NDF = %1.3f, Prob = %1.3e\n", chisqrd, NDF, redChi2, TMath::Prob(chisqrd,NDF));

  //save alpha and n parameters if fit is
  //for integrated bins in y and pT:
  if(iPTBin == 0)
    SaveCBParameters(iRapBin, iPTBin, alpha, n, fRECO->GetParError(6), fRECO->GetParError(5));

  Double_t mean2S = mean1S*(massPDG2S/massPDG1S);
  Double_t mean3S = mean1S*(massPDG3S/massPDG1S);
  Double_t sigma2S = sigma1S*(massPDG2S/massPDG1S);
  Double_t sigma3S = sigma1S*(massPDG3S/massPDG1S);

  printf("=========================================\n");
  printf("Calculate the number of Y's in the sample\n");
  printf("=========================================\n");

  TF1 *CB[kNbSpecies];
  Double_t intCBFit[kNbSpecies]; //integral values of a CB with N=Nfit and fixed alpha, n, sigma, width

  for(int iUps = 0; iUps < kNbSpecies; iUps++){
    sprintf(name, "CB_%d", iUps);
    CB[iUps] = new TF1(name, CBFunction, range_min, range_max, 5);
    CB[iUps]->SetParameter(0, 1.);
    if(iUps == 0){
      CB[iUps]->FixParameter(1, mean1S);
      CB[iUps]->FixParameter(2, sigma1S);
    }
    else if(iUps == 1){
      CB[iUps]->FixParameter(1, mean2S);
      CB[iUps]->FixParameter(2, sigma2S);
    }
    else if(iUps == 2){
      CB[iUps]->FixParameter(1, mean3S);
      CB[iUps]->FixParameter(2, sigma3S);
    }
    CB[iUps]->FixParameter(3, alpha);
    CB[iUps]->FixParameter(4, n);
    intCB[iUps] = CB[iUps]->Integral(range_min, range_max);
  }

  nY[UPS1S] = normY1S * intCB[UPS1S];
  nY[UPS2S] = normY2S * intCB[UPS2S];
  nY[UPS3S] = normY3S * intCB[UPS3S];

  printf("rapidity bin %d, pT bin %d\n", iRapBin, iPTBin);
  printf("the integral of the fitted CB for the %s is: %1.3e --> #%s = %1.3e\n", specName[UPS1S], intCB[UPS1S], specName[UPS1S], nY[UPS1S]);
  printf("the integral of the fitted CB for the %s is: %1.3e --> #%s = %1.3e\n", specName[UPS2S], intCB[UPS2S], specName[UPS2S], nY[UPS2S]);
  printf("the integral of the fitted CB for the %s is: %1.3e --> #%s = %1.3e\n", specName[UPS3S], intCB[UPS3S], specName[UPS3S], nY[UPS3S]);

  //calculate the fraction of BG in a given mass interval
  massMin[UPS1S] = mean1S - nSigma*sigma1S;
  massMin[UPS2S] = mean2S - nSigma*sigma2S;
  massMin[UPS3S] = mean3S - nSigma*sigma3S;
  massMax[UPS1S] = mean1S + nSigma*sigma1S;
  massMax[UPS2S] = mean2S + nSigma*sigma2S;
  massMax[UPS3S] = mean3S + nSigma*sigma3S;
  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++){
    printf("integrating histos between %1.3f and %1.3f GeV (+- %1.1f sigma window)\n",
	   massMin[iSpecies], massMax[iSpecies], nSigma);
  }

  fUps1S =  new TF1("fUps1S", CBFunction, range_min, range_max, 5);
  fUps1S->FixParameter(0, normY1S * binWidth);
  fUps1S->FixParameter(1, mean1S);
  fUps1S->FixParameter(2, sigma1S);
  fUps1S->FixParameter(3, alpha);
  fUps1S->FixParameter(4, n);

  fUps2S = new TF1("fUps2S", CBFunction, range_min, range_max, 5);
  fUps2S->FixParameter(0, normY2S * binWidth);
  fUps2S->FixParameter(1, mean2S);
  fUps2S->FixParameter(2, sigma2S);
  fUps2S->FixParameter(3, alpha);
  fUps2S->FixParameter(4, n);

  fUps3S = new TF1("fUps3S", CBFunction, range_min, range_max, 5);
  fUps3S->FixParameter(0, normY3S * binWidth);
  fUps3S->FixParameter(1, mean3S);
  fUps3S->FixParameter(2, sigma3S);
  fUps3S->FixParameter(3, alpha);
  fUps3S->FixParameter(4, n);

  Double_t nUps[kNbSpecies];
  nUps[UPS1S] = fUps1S->Integral(massMin[UPS1S], massMax[UPS1S]);
  nUps[UPS2S] = fUps2S->Integral(massMin[UPS2S], massMax[UPS2S]);
  nUps[UPS3S] = fUps3S->Integral(massMin[UPS3S], massMax[UPS3S]);

  fBG = new TF1("fBG", DrawContinuum, range_min, range_max, 3);
  for(int iPar = 0; iPar < 3; iPar++)
    fBG->FixParameter(iPar, fitParBG[iPar]);

  Double_t nBG[kNbSpecies];
  nBG[UPS1S] = 0;
  nBG[UPS2S] = 0;
  nBG[UPS3S] = 0;

  if(iRapBin == 0 || iPTBin == 0)
    printf("rapidity bin %d, pTBin %d\n", iRapBin, iPTBin);
  else{
    printf("rapidity bin %d (%1.1f - %1.1f), pT bin %d (%1.0f - %1.0f GeV/c)\n",
	   iRapBin, onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin],
	   iPTBin, onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin]);
  }

  FILE *fInfo = fopen("tmpFiles/statistics.txt", "a");
  if(iRapBin == 0 || iPTBin == 0)
    fprintf(fInfo, "rapidity bin %d, pTBin %d\n", iRapBin, iPTBin);
  else{
    fprintf(fInfo, "rapidity bin %d (%1.1f - %1.1f), pT bin %d (%1.0f - %1.0f GeV/c)\n",
	   iRapBin, onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin],
	   iPTBin, onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin]);
  }
  for(int iSpecies = 0; iSpecies < kNbSpecies; iSpecies++){
    fracBG[iSpecies] = 0.01;
    // printf("nUps = %1.3f\n", nUps[iSpecies]);
    // printf("nBG = %1.3f\n", nBG[iSpecies]);
    printf("%s: fraction of BG in a +- %1.1f sigma window is %1.3f\n",
	   specName[iSpecies], nSigma, fracBG[iSpecies]);

    fprintf(fInfo, "integral of the fitted CB for the %s is (norm factor= %1.3e): %7.0f\n", specName[iSpecies], intCB[iSpecies], nY[iSpecies]);
  }
  fclose(fInfo);

  if(iPTBin > 0 && iRapBin > 0)
    SaveFitPars(iRapBin, iPTBin);

}

//==============================
void GetHisto(Char_t *fileNameIn, Int_t iRapBin, Int_t iPTBin){

  TFile *fin = new TFile(fileNameIn);
  Char_t name[100];
  sprintf(name, "Reco_Onia_mass_rap%d_pT%d", iRapBin, iPTBin);
  hMass = (TH1F*) gDirectory->Get(name);
  hMass->Rebin(2);
  binWidth = hMass->GetBinWidth(1); //valid only for an equal bin histogram!
  printf("binwidth = %1.2e\n", binWidth);
}

//==============================
void SaveCBParameters(Int_t iRapBin, Int_t iPTBin, Double_t alpha, Double_t n, Double_t alphaErr, Double_t nErr){

  Char_t name[100];
  if(iPTBin == 0 && iRapBin == 0)
    sprintf(name, "tmpFiles/CBParameters.root");
  else if(iPTBin == 0)
    sprintf(name, "tmpFiles/CBParameters_rap%d.root", iRapBin);
  else{
    printf("<SaveCBParameters> can only be called for pT = 0!\n");
    exit(0);
  }
  TFile *fOut = new TFile(name, "RECREATE");
  TTree *treeOut = new TTree("CBPars", "");
  Double_t alphaAll = alpha, nAll = n;
  Double_t alphaAllErr = alphaErr, nAllErr = nErr;
  treeOut->Branch("alphaAll", &alphaAll, "alphaAll/D");
  treeOut->Branch("nAll", &nAll, "nAll/D");
  treeOut->Branch("alphaAllErr", &alphaAllErr, "alphaAllErr/D");
  treeOut->Branch("nAllErr", &nAllErr, "nAllErr/D");
  treeOut->Fill();
  treeOut->Write();
  fOut->Close();
}

//=========================
void SaveFitPars(Int_t iRapBin, Int_t iPTBin){

  Char_t name[100];
  sprintf(name, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
  TFile *fOut = new TFile(name, "RECREATE");
  TTree *treeOut = new TTree("massFitParameters", "");
  Int_t bufsize = 32000; //default = 32000
  Int_t splitlevel = 0; //recommended by R. Brun
  treeOut->Branch("fUps1S", "TF1", &fUps1S, bufsize, splitlevel);
  treeOut->Branch("fUps2S", "TF1", &fUps2S, bufsize, splitlevel);
  treeOut->Branch("fUps3S", "TF1", &fUps3S, bufsize, splitlevel);
  treeOut->Branch("fBG", "TF1", &fBG, bufsize, splitlevel);
  treeOut->Fill();

  treeOut->Write();
  fOut->Close();
}

//=========================
Double_t fitPolyCrystal3(Double_t *x, Double_t *par){

  Double_t normY1S = par[0];
  Double_t mean1S = par[1];
  Double_t sigma1S = par[2];
  Double_t normY2S = par[3];
  Double_t normY3S = par[4];
  Double_t n = par[5];
  Double_t alpha = par[6];
  Double_t a = par[7];
  Double_t b = par[8];
  Double_t c = par[9];

  Double_t mean2S = mean1S*(massPDG2S/massPDG1S);
  Double_t mean3S = mean1S*(massPDG3S/massPDG1S);
  Double_t sigma2S = sigma1S*(massPDG2S/massPDG1S);
  Double_t sigma3S = sigma1S*(massPDG3S/massPDG1S);

  Double_t poly2 = a + b*x[0] + c*x[0]*x[0];

  Double_t CB1 = 0.;
  if(((x[0] - mean1S)/sigma1S) > -alpha)
    CB1 = TMath::Exp(-(pow(x[0]-mean1S,2)/(2.*sigma1S*sigma1S)));
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB1 = A*pow(B - (x[0]-mean1S)/sigma1S, -n);
  }

  Double_t CB2 = 0.;
  if(((x[0] - mean2S)/sigma2S) > -alpha)
    CB2 = TMath::Exp(-(pow(x[0]-mean2S,2)/(2.*sigma2S*sigma2S)));
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB2= A*pow(B - (x[0]-mean2S)/sigma2S, -n);
  }
  Double_t CB3 = 0.;
  if(((x[0] - mean3S)/sigma3S) > -alpha)
    CB3 = TMath::Exp(-(pow(x[0]-mean3S,2)/(2.*sigma3S*sigma3S)));
  else{
    Double_t A = pow(n / TMath::Abs(alpha),n) * TMath::Exp(-alpha*alpha/2.);
    Double_t B = n / TMath::Abs(alpha) - TMath::Abs(alpha);
    CB3 = A*pow(B - (x[0]-mean3S)/sigma3S, -n);
  }

  Double_t result = normY1S*CB1;

  result *= binWidth; //correct for the bin width
  return result;
}

//==================================
Double_t fitContinuum(Double_t *x, Double_t *par){

  if (x[0] > peak_min && x[0] < peak_max) {
    TF1::RejectPoint();
    return 0;
  }
  Double_t result = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  result *= binWidth; //correct for the bin width
  return result;
}

//==================================
Double_t DrawContinuum(Double_t *x, Double_t *par){

  Double_t result = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
  result *= binWidth; //correct for the bin width
  return result;
}


