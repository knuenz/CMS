#include "rootIncludes.inc"
#include "calcPol.C"
#include "TH3D.h"
#include "TRandom3.h"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};

double contamination2Sin1S;
double contamination1Sin2S;
double contamination3Sin2S;
double contamination2Sin3S;
double backgroundFrac;
double err_backgroundFrac;
double SignalEvents;

double TrimEventContent(Int_t iRapBin = 1,
		      Int_t iPTBin = 1,
		      Double_t fracL = 0.5, Double_t nSigma = 2., 
		      Int_t nUpsState=0,//[0]... 1S, [1]... 2S, [2]... 3S
		      bool UpsMC=false,
		      bool f_BG_zero=false,
		      bool ProjectLSBdata=false,
		      bool ProjectRSBdata=false,
		      bool CombineSignalPeaks=false,
		      bool Y1Sto2S_SB=false,
		      bool LeftSided=false,
		      bool RightSided=false,
		      bool MassScan=false,
		      bool adjustOverlapBorders=true
		      ){

  printf("\n\n\nfracL = %1.3f, nSigma = %1.1f, iState = %d, rap %d, pT %d\n", fracL, nSigma, nUpsState, iRapBin, iPTBin);

  Char_t name[100], title[100];
  Char_t fileNameIn[100];
  sprintf(fileNameIn, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
  //==============================
  //read inputs from input file:
  TFile *fIn = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");
  TH2D *hBG_cosThetaPhiLR[onia::kNbFrames][2];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
    hBG_cosThetaPhiLR[iFrame][L] = (TH2D *) gDirectory->Get(name);
    sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
    hBG_cosThetaPhiLR[iFrame][R] = (TH2D *) gDirectory->Get(name);
  }
  //==============================

  //definition of output variables 
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "AllStates_%1.2fSigma_FracLSB%dPercent/data_%dSUps_rap%d_pT%d.root", nSigma, int(fracL*100), nUpsState+1, iRapBin, iPTBin);
  TFile *fOut = new TFile(fileNameOut, "RECREATE");
  gStyle->SetPadRightMargin(0.2);
  TTree *treeOut = treeIn->CloneTree(0);
  // treeOut->SetName("data");
  TH2D *hBG_cosThetaPhi[onia::kNbFrames];
  // TH2D *hBG_cosThetaPhiSignal[onia::kNbFrames];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    // //book the histo for the signal
    // sprintf(name, "total_%s", onia::frameLabel[iFrame]);
    // sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    // hBG_cosThetaPhiSignal[iFrame] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
    // 					  onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    // hBG_cosThetaPhiSignal[iFrame]->Sumw2();
    //copy the L and R sideband histos into one output BG histogram
    hBG_cosThetaPhiLR[iFrame][L]->Scale(fracL/hBG_cosThetaPhiLR[iFrame][L]->Integral());
    hBG_cosThetaPhiLR[iFrame][R]->Scale((1.-fracL)/hBG_cosThetaPhiLR[iFrame][R]->Integral());
    sprintf(name, "background_costhphi%s", onia::frameLabel[iFrame]);
    hBG_cosThetaPhi[iFrame] = (TH2D *) hBG_cosThetaPhiLR[iFrame][L]->Clone(name);
    hBG_cosThetaPhi[iFrame]->Add(hBG_cosThetaPhiLR[iFrame][R]);
  }

  //==========================================================
  //reading fit parameters to establish signal mass window
  //as well as the L and R sideband window for the 3D BG histo
  //==========================================================
  fIn->cd();
  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
  TF1 *fUps[kNbSpecies], *fBG = 0;
  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
  treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
  treeFitPar->SetBranchAddress("fBG", &fBG);
  treeFitPar->LoadTree(0);
  treeFitPar->GetEntry(0);


  Double_t mass[kNbSpecies], sigma[kNbSpecies];
  for(int iState = 0; iState < kNbSpecies; iState++){
    mass[iState] = fUps[iState]->GetParameter(1);
    sigma[iState] = fUps[iState]->GetParameter(2);
  }
  printf("1S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS1S], sigma[UPS1S]);
  printf("2S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS2S], sigma[UPS2S]);
  printf("3S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS3S], sigma[UPS3S]);
  Double_t poleMass = mass[nUpsState], massMin, massMax;
  massMin = poleMass - nSigma*sigma[nUpsState];
  massMax = poleMass + nSigma*sigma[nUpsState];

  if(LeftSided){
	  massMin = poleMass - nSigma*sigma[nUpsState];
	  massMax = poleMass;
  }

  if(RightSided){
	  massMin = poleMass;
	  massMax = poleMass + nSigma*sigma[nUpsState];
  }


  if(adjustOverlapBorders){
  if( nUpsState==2 && mass[2]-nSigma*sigma[2]<mass[1]+nSigma*sigma[1] ){
	  cout<<"adjusting lower border of Y(3S) mass window due to overlap"<<endl;
	  massMin=(sigma[2]*mass[1]+sigma[1]*mass[2])/(sigma[1]+sigma[2]);
  }
  if( nUpsState==1 && mass[1]-nSigma*sigma[1]<mass[0]+nSigma*sigma[0] ){
	  cout<<"adjusting lower border of Y(2S) mass window due to overlap"<<endl;
	  massMin=(sigma[1]*mass[0]+sigma[0]*mass[1])/(sigma[1]+sigma[0]);
  }
  if( nUpsState==1 && mass[2]-nSigma*sigma[2]<mass[1]+nSigma*sigma[1] ){
	  cout<<"adjusting upper border of Y(2S) mass window due to overlap"<<endl;
	  massMax=(sigma[2]*mass[1]+sigma[1]*mass[2])/(sigma[1]+sigma[2]);
  }
  if( nUpsState==0 && mass[1]-nSigma*sigma[1]<mass[0]+nSigma*sigma[0] ){
	  cout<<"adjusting upper border of Y(1S) mass window due to overlap"<<endl;
	  massMax=(sigma[1]*mass[0]+sigma[0]*mass[1])/(sigma[1]+sigma[0]);
  }
  }


  if( nUpsState==0){
  contamination2Sin1S=fUps[1]->Integral(massMin, massMax)/(fUps[0]->Integral(massMin, massMax)+fUps[1]->Integral(massMin, massMax));
  cout<<"Contamination of the 1S sample by 2S events = "<<contamination2Sin1S*100<<" %"<<endl;
  }
  if( nUpsState==1){
  contamination1Sin2S=fUps[0]->Integral(massMin, massMax)/(fUps[0]->Integral(massMin, massMax)+fUps[1]->Integral(massMin, massMax));
  contamination3Sin2S=fUps[2]->Integral(massMin, massMax)/(fUps[2]->Integral(massMin, massMax)+fUps[1]->Integral(massMin, massMax));
  cout<<"Contamination of the 2S sample by 1S events = "<<contamination1Sin2S*100<<" %"<<endl;
  cout<<"Contamination of the 2S sample by 3S events = "<<contamination3Sin2S*100<<" %"<<endl;
  }
  if( nUpsState==2){
  contamination2Sin3S=fUps[1]->Integral(massMin, massMax)/(fUps[2]->Integral(massMin, massMax)+fUps[1]->Integral(massMin, massMax));
  cout<<"Contamination of the 3S sample by 2S events = "<<contamination2Sin3S*100<<" %"<<endl;
  }

  if(CombineSignalPeaks) {
	  massMin=9.15;
	  massMax=10.65;
  }

  if(Y1Sto2S_SB){
	  massMin = mass[0] + 3.*sigma[0];
	  massMax = mass[1] - 4.*sigma[1];
  }

////////////////// Background mass scan
const int nMassScan=12;
int nMassScanCurrent;
//double MassScanBorders[nMassScan+1] = {8.6, 8.95, 9.3, 9.45, 9.6, 9.85, 10.0125, 10.175, 10.3425, 10.51, 10.8, 11.1, 11.4};
double MassScanBorders[nMassScan+1] = {8.6, 8.95, 9.3, 9.45, 9.6, 9.85, 10.0125, 10.175, 10.3425, 10.51, 10.8, 11.1, 11.4};

for(int i=1;i<nMassScan+1;i++){
if(nSigma<i+0.5) {nMassScanCurrent=i; break;}
}

double BuffMinL=onia::massMinL;
double BuffMaxL=mass[UPS1S] - onia::nSigmaL*sigma[UPS1S];
double BuffMinR=mass[UPS3S] + onia::nSigmaR*sigma[UPS3S];
double BuffMaxR=onia::massMaxR;

double MassScanMin=MassScanBorders[nMassScanCurrent-1];
double MassScanMax=MassScanBorders[nMassScanCurrent];

if(nSigma==1){
	MassScanMin=BuffMinL;
	MassScanMax=BuffMinL+1./3.*(BuffMaxL-BuffMinL);
}
if(nSigma==2){
	MassScanMin=BuffMinL+1./3.*(BuffMaxL-BuffMinL);
	MassScanMax=BuffMinL+2./3.*(BuffMaxL-BuffMinL);
}
if(nSigma==3){
	MassScanMin=BuffMinL+2./3.*(BuffMaxL-BuffMinL);
	MassScanMax=BuffMaxL;
}
if(nSigma==4){
	MassScanMin=BuffMinL;
	MassScanMax=BuffMaxL;
}

if(nSigma==5){
	MassScanMin=BuffMinR;
	MassScanMax=BuffMinR+1./3.*(BuffMaxR-BuffMinR);
}
if(nSigma==6){
	MassScanMin=BuffMinR+1./3.*(BuffMaxR-BuffMinR);
	MassScanMax=BuffMinR+2./3.*(BuffMaxR-BuffMinR);
}
if(nSigma==7){
	MassScanMin=BuffMinR+2./3.*(BuffMaxR-BuffMinR);
	MassScanMax=BuffMaxR;
}
if(nSigma==8){
	MassScanMin=BuffMinR;
	MassScanMax=BuffMaxR;
}
if(nSigma>8){
	BuffMaxL=mass[UPS1S] - 7*sigma[UPS1S];
	BuffMinR=mass[UPS3S] + onia::nSigmaR*sigma[UPS3S];
}
if(nSigma==9){
	MassScanMin=BuffMinL;
	MassScanMax=BuffMinL+1./3.*(BuffMaxL-BuffMinL);
}
if(nSigma==10){
	MassScanMin=BuffMinL+1./3.*(BuffMaxL-BuffMinL);
	MassScanMax=BuffMinL+2./3.*(BuffMaxL-BuffMinL);
}
if(nSigma==11){
	MassScanMin=BuffMinL+2./3.*(BuffMaxL-BuffMinL);
	MassScanMax=BuffMaxL;
}
if(nSigma==12){
	MassScanMin=BuffMinL;
	MassScanMax=BuffMaxL;
}


  if(MassScan){
  massMin=MassScanMin;
  massMax=MassScanMax;
  }
/////////////////////////////////////

  printf("--> signal mass window: %1.3f < M < %1.3f GeV\n", massMin, massMax);

  //calculate the fraction of BG under the signal mass window
  Double_t nBG = fBG->Integral(massMin, massMax);
//  Double_t nSignal = fUps[nUpsState]->Integral(massMin, massMax);
  Double_t nSignal = fUps[0]->Integral(massMin, massMax)+fUps[1]->Integral(massMin, massMax)+fUps[2]->Integral(massMin, massMax);
  Double_t fracBG = nBG / (nBG + nSignal);
  sprintf(name, ";;fraction of BG in %1.1f sigma window", nSigma);
  TH1D *hFracBG = new TH1D("background_fraction", name, 1, 0., 1.);
  if(!UpsMC&&!f_BG_zero) hFracBG->SetBinContent(1, fracBG);
  if(UpsMC||f_BG_zero)  hFracBG->SetBinContent(1, 0.001);

 /* fBG->ReleaseParameter(1);
  cout<<"fBG "<<fBG->Integral(massMin, massMax) / (fBG->Integral(massMin, massMax) + fUps[nUpsState]->Integral(massMin, massMax))<<endl;
  double fBGPar0=fBG->GetParameter(1);
  double err_fBGPar0=fBG->GetParError(1);
  cout<<"fBGPar0 "<<fBGPar0<<endl;
  cout<<"err_fBGPar0 "<<err_fBGPar0<<endl;

  fBG->FixParameter(1,fBGPar0+100*err_fBGPar0);
  cout<<"fBGPar0 "<<fBG->GetParameter(1)<<endl;
  cout<<"fBG "<<fBG->Integral(massMin, massMax) / (fBG->Integral(massMin, massMax) + fUps[nUpsState]->Integral(massMin, massMax))<<endl;
  fBG->SetParameter(1,fBGPar0);
  cout<<"fBG "<<fBG->Integral(massMin, massMax) / (fBG->Integral(massMin, massMax) + fUps[nUpsState]->Integral(massMin, massMax))<<endl;
*/

  if(Y1Sto2S_SB){
cout<<"Y1Sto2S_SB fBG = "<<fBG->Integral(massMin, massMax) / (fBG->Integral(massMin, massMax) + fUps[0]->Integral(massMin, massMax)+fUps[1]->Integral(massMin, massMax))<<endl;
  }

  backgroundFrac=fracBG;

  //calculate the L and R mass windows:
  Double_t massMinBG[2], massMaxBG[2];
  massMinBG[L] = onia::massMinL;
  massMaxBG[L] = mass[UPS1S] - onia::nSigmaL*sigma[UPS1S];
  massMinBG[R] = mass[UPS3S] + onia::nSigmaR*sigma[UPS3S];
  massMaxBG[R] = onia::massMaxR;
  printf("--> L mass window: %1.3f < M < %1.3f GeV\n", massMinBG[L], massMaxBG[L]);
  printf("--> R mass window: %1.3f < M < %1.3f GeV\n", massMinBG[R], massMaxBG[R]);
//cout<<"here?"<<endl;
  bool PseudoBin=false;

  if(CombineSignalPeaks) {
	  massMin=9.15;
	  massMax=10.65;
  }

  if(massMaxBG[L] > massMin){
    printf("the right sideband window is LARGER than the left signal window!!!!\n\n\n\n");
    massMaxBG[L]=9.;
    massMinBG[L]=8.6;
    massMin=9.;
    PseudoBin=true;
//    exit(0);
  }
  if(massMinBG[R] < massMax){
    printf("the left sideband window is SMALLER than the right signal window!!!!\n\n\n\n");
    massMaxBG[R]=11.4;
    massMinBG[R]=11.;
    massMax=10.;
    PseudoBin=true;
//    exit(0);
  }

  if(ProjectLSBdata){
	  massMin=massMinBG[L];
	  massMax=massMaxBG[L];
  }
  if(ProjectRSBdata){
	  massMin=massMinBG[R];
	  massMax=massMaxBG[R];
  }

  if(CombineSignalPeaks) {
	  massMin=9.15;
	  massMax=10.65;
  }

  double meanMass=0;
  double MassDistCurrent=0.001;
  TH1D *hMassScanInfo = new TH1D("hMassScanInfo", "hMassScanInfo", 2, 0., 1.);

  if(MassScan){
  massMin=MassScanMin;
  massMax=MassScanMax;

  double IntCurrent = fBG->Integral(massMin, massMax);
  for(int i=1;i<1000000;i++){
	  if(fBG->Integral(massMin, massMin+i*MassDistCurrent)>IntCurrent/2.) {meanMass=massMin+i*MassDistCurrent; break;}
  }

  hMassScanInfo->SetBinContent(1,meanMass);
  hMassScanInfo->SetBinContent(2,1-fracBG);

  }

  if(Y1Sto2S_SB){
	  massMin = mass[0] + 3.*sigma[0];
	  massMax = mass[1] - 4.*sigma[1];

	  double IntCurrent = fBG->Integral(massMin, massMax);
	  for(int i=1;i<1000000;i++){
		  if(fBG->Integral(massMin, massMin+i*MassDistCurrent)>IntCurrent/2.) {meanMass=massMin+i*MassDistCurrent; break;}
	  }
	  hMassScanInfo->SetBinContent(1,meanMass);
	  hMassScanInfo->SetBinContent(2,1-fracBG);
  }

  //calculate central fracL

  double fracLCentral;

  if(!f_BG_zero){
  double mean_LSB;
  double mean_RSB;
  double mean_nS;

  double MassDist=0.001;
  double IntLSB = fBG->Integral(massMinBG[L], massMaxBG[L]);
  for(int i=1;i<1000000;i++){
	  if(fBG->Integral(massMinBG[L], massMinBG[L]+i*MassDist)>IntLSB/2.) {mean_LSB=massMinBG[L]+i*MassDist; break;}
  }
  double IntRSB = fBG->Integral(massMinBG[R], massMaxBG[R]);
  for(int i=1;i<1000000;i++){
	  if(fBG->Integral(massMinBG[R], massMinBG[R]+i*MassDist)>IntRSB/2.) {mean_RSB=massMinBG[R]+i*MassDist; break;}
  }
  double IntSig = fBG->Integral(massMin, massMax);
  for(int i=1;i<1000000;i++){
	  if(fBG->Integral(massMin, massMin+i*MassDist)>IntSig/2.) {mean_nS=massMin+i*MassDist; break;}
  }

  fracLCentral=1-(mean_nS-mean_LSB)/(mean_RSB-mean_LSB);
  cout<<"Median LSB: "<<mean_LSB<<endl;
  cout<<"Median RSB: "<<mean_RSB<<endl;
  cout<<"Median signal region: "<<mean_nS<<endl;
  cout<<"Central FracL: "<<fracLCentral<<endl;


  }



  //build the 3D (pT, |y|, M) histos for the L and R mass sideband 
  TH3D *hBG_pTRapMass[2];
  hBG_pTRapMass[L] = new TH3D("hBG_pTRapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]", 
			      7, onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin],
			      2, onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin],
			      7, massMin, massMax);//*signal* mass window!
  hBG_pTRapMass[L]->Sumw2();
  //
  hBG_pTRapMass[R] = new TH3D("hBG_pTRapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]", 
			      7, onia::pTRange[iRapBin][iPTBin-1], onia::pTRange[iRapBin][iPTBin],
			      2, onia::rapForPTRange[iRapBin-1], onia::rapForPTRange[iRapBin],
			      7, massMin, massMax);//*signal* mass window!
  hBG_pTRapMass[R]->Sumw2();

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();


  Double_t onia_mass;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 10000 == 0)
      cout << "entry " << iEntry << " out of " << treeIn->GetEntries() << endl;

    *onia = *(lepP) + *(lepN);
    onia_mass = onia->M();

	  if(UpsMC||f_BG_zero) hBG_pTRapMass[L]->Fill(onia->Pt(), TMath::Abs(onia->Rapidity()), gRandom->Uniform(massMin, massMax));
	  if(UpsMC||f_BG_zero) hBG_pTRapMass[R]->Fill(onia->Pt(), TMath::Abs(onia->Rapidity()), gRandom->Uniform(massMin, massMax));

	  if(!UpsMC&&!f_BG_zero){
    if(onia_mass  > massMinBG[L] && onia_mass < massMaxBG[L])
      hBG_pTRapMass[L]->Fill(onia->Pt(), TMath::Abs(onia->Rapidity()), fBG->GetRandom(massMin, massMax));
    else if(onia_mass > massMinBG[R] && onia_mass < massMaxBG[R])
      hBG_pTRapMass[R]->Fill(onia->Pt(), TMath::Abs(onia->Rapidity()), fBG->GetRandom(massMin, massMax));
	  }

    if(onia_mass > massMin && onia_mass < massMax){
      treeOut->Fill(); //stores TLorenzVectors of the two muons

      // //store now the cosTheta and phi distributions of the signal window:
      // calcPol(*lepP, *lepN);

      // for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
      // 	hCosThetaPhiSignal[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
    }


  }



/*  Char_t name2[100];
  sprintf(name2, "Reco_Onia_mass_rap%d_pT%d", iRapBin, iPTBin);
  TH1F* hMass = (TH1F*) gDirectory->Get(name2);
  hMass->Rebin(2);
  double binWidth = hMass->GetBinWidth(1); //valid only for an equal bin histogram!
  printf("binwidth = %1.2e\n", binWidth);
*/
  double nY[3];
  double nBG_;
  double binWidth = 0.02; //You have to manually change that parameter from the mass fit!!!

    nY[2] = fUps[2]->Integral(massMin, massMax)/binWidth;
    nY[1] = fUps[1]->Integral(massMin, massMax)/binWidth;
    nY[0] = fUps[0]->Integral(massMin, massMax)/binWidth;
    nBG_ = fBG->Integral(massMin, massMax)/binWidth;
    double nAll=nY[0]+nY[1]+nY[2]+nBG_;

    printf("1 sigma region num background = %1.3f\n",nBG_);
    printf("1 sigma region num 1S = %1.3f\n",nY[0]);
    printf("1 sigma region num 2S = %1.3f\n",nY[1]);
    printf("1 sigma region num 3S = %1.3f\n",nY[2]);
    printf("1 sigma region num all = %1.3f\n",nAll);
    printf("1 sigma region calcNumAll - numEventsInDataSet = %1.3f\n",nAll-treeOut->GetEntries());

    double tempBfrac=nBG_/(treeOut->GetEntries());
    printf("1 sigma region background fraction = %1.3f\n",tempBfrac);
    printf("Used           background fraction = %1.3f\n",backgroundFrac);

    SignalEvents=nY[nUpsState];

    double BackgroundEvents=SignalEvents*(1-backgroundFrac)/backgroundFrac;

    err_backgroundFrac=TMath::Power(SignalEvents*BackgroundEvents,0.5)/TMath::Power(SignalEvents+BackgroundEvents,1.5);

  //now, add the L and R 3D (pT, |y|, M) histos
  TH3D *hBG_pTRapMassSum;
  hBG_pTRapMass[L]->Scale(fracL/hBG_pTRapMass[L]->Integral());
  hBG_pTRapMass[R]->Scale((1.-fracL)/hBG_pTRapMass[R]->Integral());
  sprintf(name, "background_pTrapMass");
  hBG_pTRapMassSum = (TH3D*) hBG_pTRapMass[L]->Clone(name);
  hBG_pTRapMassSum->Add(hBG_pTRapMass[R]);

  //write the output
  fOut->cd();
  treeOut->Write();

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    hBG_cosThetaPhi[iFrame]->Write();
    hBG_cosThetaPhiLR[iFrame][L]->Write();
    hBG_cosThetaPhiLR[iFrame][R]->Write();
    // hCosThetaPhiSignal[iFrame]->Write();
  }
  hBG_pTRapMassSum->Write();
  hFracBG->Write();
  if(MassScan){
	  hMassScanInfo->Write();
  }
  fOut->Close();
  fIn->Close();

  return fracLCentral;
}
