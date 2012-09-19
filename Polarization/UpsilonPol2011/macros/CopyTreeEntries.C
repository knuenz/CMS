#include "rootIncludes.inc"
#include "TLorentzVector.h"
#include "calcPol.C"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};
//=====================================================
void CopyTreeEntries(Int_t iRapBin = 1, 
		     Int_t iPTBin = 1,
		     Char_t *fileNameIn = "RootFiles/selEvents_data_Ups_noCowboys_3Aug2011.root",
		     bool UpsMC=false,
		     bool DoCPUconsumingPlots=false){

  Char_t name[100], title[100];  
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "tmpFiles/data_Ups_rap%d_pT%d.root", iRapBin, iPTBin);
  printf("updating file %s\n", fileNameOut);

  //==============================
  TFile *fIn = new TFile(fileNameIn);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) gDirectory->Get("selectedData");

  //==============================

  //==============================
  //definition of output variables 
  TFile *fOut = new TFile(fileNameOut, "UPDATE");
  gStyle->SetPadRightMargin(0.2);
  TTree *treeOut = treeIn->CloneTree(0);
  TH2D *hBG_cosThetaPhi[onia::kNbFrames][2];
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    //book the 2D (cosTheta, phi) histos for the L and R mass sideband
    sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    hBG_cosThetaPhi[iFrame][L] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
			       onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhi[iFrame][L]->Sumw2();
    //
    sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
    sprintf(title, ";cos#theta_{%s};#phi_{%s} [deg]", onia::frameLabel[iFrame], onia::frameLabel[iFrame]);
    hBG_cosThetaPhi[iFrame][R] = new TH2D(name, title, onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, 
			       onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
    hBG_cosThetaPhi[iFrame][R]->Sumw2();
  }

  //==============================
  //reading info from input file:
  //===============================
  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
  if(gDirectory->Get("massFitParameters")==NULL){
    printf("\n\n\nskip processing this bin.\n\n\n");
    return;
  }

  TF1 *fUps[kNbSpecies], *fBG = 0;
  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
  treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
  treeFitPar->LoadTree(0);
  treeFitPar->GetEntry(0);
  fUps[0]->Print();

  Double_t mass1S = fUps[UPS1S]->GetParameter(1);
  Double_t sigma1S = fUps[UPS1S]->GetParameter(2);
  Double_t mass2S = fUps[UPS2S]->GetParameter(1);
  Double_t sigma2S = fUps[UPS2S]->GetParameter(2);
  Double_t mass3S = fUps[UPS3S]->GetParameter(1);
  Double_t sigma3S = fUps[UPS3S]->GetParameter(2);
  printf("1S: mass = %1.3f, sigma = %1.3f\n", mass1S, sigma1S);
  printf("3S: mass = %1.3f, sigma = %1.3f\n", mass3S, sigma3S);
  Double_t massMin[2], massMax[2];
  massMin[L] = onia::massMinL;
  massMax[L] = mass1S - onia::nSigmaL*sigma1S;
  massMin[R] = mass3S + onia::nSigmaR*sigma3S;
  massMax[R] = onia::massMaxR;
  printf("--> L mass window: %1.3f < M < %1.3f GeV\n", massMin[L], massMax[L]);
  printf("--> R mass window: %1.3f < M < %1.3f GeV\n", massMin[R], massMax[R]);

  lepP = 0; lepN = 0;
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();

  bool CMSprelim=true;
  bool PlotSpecials=false;
  if(DoCPUconsumingPlots) PlotSpecials=true;

  if(iRapBin==0&&iPTBin==0&&PlotSpecials){
	  char CutChar[200];
	  char savename[200];
	  TCanvas *c2;




	  double DMmin=8.6;
	  double DMmax=11.4;
	  double DMbinwidth1=0.04;
	  double DMbinwidth2=0.04;
	  int DMnBins1=(DMmax-DMmin)/DMbinwidth1;
	  int DMnBins2=(DMmax-DMmin)/DMbinwidth2;

	  TH1D* dimumassRap1   = new TH1D( "dimumassRap1", "dimumassRap1", DMnBins1, DMmin,DMmax);
	  TH1D* dimumassRap2   = new TH1D( "dimumassRap2", "dimumassRap2", DMnBins2, DMmin,DMmax);

	  char AbsRapChar[1000];
	  char PtChar[1000];
	  char DimuonMassChar[1000];

	  sprintf(AbsRapChar,"TMath::Abs(1/2*TMath::Log((TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())+(lepP->Pz()+lepN->Pz())*(lepP->Pz()+lepN->Pz())+2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))+(lepP->Pz()+lepN->Pz()))/(TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())+(lepP->Pz()+lepN->Pz())*(lepP->Pz()+lepN->Pz())+2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))-(lepP->Pz()+lepN->Pz()))))");
	  sprintf(PtChar,"TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py()))");
	  sprintf(DimuonMassChar,"TMath::Sqrt(2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))");

	  char DMDrawChar[2000];
	  char DMCutChar[2000];


//		TLegend* DMLegend=new TLegend(0.53,0.75,0.95,0.9);
		TLegend* DMLegend=new TLegend(0.68,0.585,0.95,0.725);
		DMLegend->SetFillColor(0);
//		DMLegend->SetTextFont(72);
		DMLegend->SetTextSize(0.0345);
		DMLegend->SetBorderSize(1);
//		DMLegend->SetMargin(0.135);
		char DMLegendEntry[200];


	  double ptMin=10;
	  double ptMax=50;
	  double rapMin=0;
	  double rapMax=0.6;
	  sprintf(DMDrawChar,"%s>>dimumassRap1",DimuonMassChar);
	  sprintf(DMCutChar,"%s>%f && %s<%f && %s>%f && %s<%f",PtChar,ptMin,PtChar,ptMax,AbsRapChar,rapMin,AbsRapChar,rapMax);

	  treeIn->Draw(DMDrawChar,DMCutChar);

	  ptMin=10;
	  ptMax=50;
	  rapMin=0.6;
	  rapMax=1.2;
	  sprintf(DMDrawChar,"%s>>dimumassRap2",DimuonMassChar);
	  sprintf(DMCutChar,"%s>%f && %s<%f && %s>%f && %s<%f",PtChar,ptMin,PtChar,ptMax,AbsRapChar,rapMin,AbsRapChar,rapMax);

	  treeIn->Draw(DMDrawChar,DMCutChar);

	  char DMyTitle[200];
	  sprintf(DMyTitle,"Counts per %d MeV",int(1000*DMbinwidth1));
	  c2 = new TCanvas("c2","c2",1200,1100);
	  gStyle->SetPalette(1);
 	  gPad->SetFillColor(kWhite);

 	  double MargLeft=0.125;//0.175;
 	  double MargRight=0.025;//0.1;
 	  double MargTop=0.05;//0.175;

 	  double MarkerSizeTom=1.75;

 	  gPad->SetLeftMargin(MargLeft);
      gPad->SetRightMargin(MargRight);
      gPad->SetTopMargin(MargTop);

 	  dimumassRap1->SetYTitle(DMyTitle);
	  dimumassRap1->SetXTitle("Dimuon mass [GeV]");
	  dimumassRap1->SetTitle(0);
	  dimumassRap1->SetStats(0);
	  dimumassRap1->GetYaxis()->SetTitleOffset(1.25);
	  dimumassRap1->SetMarkerStyle(25);
	  dimumassRap1->SetMarkerSize(MarkerSizeTom);
	  dimumassRap1->SetMarkerColor(kBlack);
	  dimumassRap1->SetLineColor(kBlack);
	  dimumassRap2->SetTitle(0);
	  dimumassRap2->SetStats(0);
	  dimumassRap2->SetMarkerColor(kRed);
	  dimumassRap2->SetMarkerStyle(20);
	  dimumassRap2->SetLineColor(kRed);
	  dimumassRap2->SetMarkerSize(MarkerSizeTom);


	  double SF_DM=1000;

	    for(int iX=0;iX<dimumassRap2->GetNbinsX()+1;iX++){
	        	int globalBin=dimumassRap2->GetBin(iX);
	        	dimumassRap2->SetBinError(iX,dimumassRap2->GetBinError(iX)/SF_DM);
	        	dimumassRap2->SetBinContent(iX,dimumassRap2->GetBinContent(iX)/SF_DM);
	        }
	    for(int iX=0;iX<dimumassRap1->GetNbinsX()+1;iX++){
	        	int globalBin=dimumassRap1->GetBin(iX);
	        	dimumassRap1->SetBinError(iX,dimumassRap1->GetBinError(iX)/SF_DM);
	        	dimumassRap1->SetBinContent(iX,dimumassRap1->GetBinContent(iX)/SF_DM);
	        }


	  dimumassRap1->Draw("E1");
	  dimumassRap2->Draw("same,E1");
	  dimumassRap1->Draw("same,chist");
	  dimumassRap2->Draw("same,chist");


//		sprintf(DMLegendEntry,"|#it{y}| < 0.6, per %d MeV",int(1000*DMbinwidth1));
		sprintf(DMLegendEntry,"|#it{y}| < 0.6");
		DMLegend->AddEntry(dimumassRap1,DMLegendEntry,"p");
//		sprintf(DMLegendEntry,"0.6 < |#it{y}| < 1.2, per %d MeV",int(1000*DMbinwidth2));
		sprintf(DMLegendEntry,"0.6 < |#it{y}| < 1.2");
		DMLegend->AddEntry(dimumassRap2,DMLegendEntry,"p");

		DMLegend->Draw();

		double xCentralsDM=10.475;
		double xCentralsDM_shift=0.;

		 cout<<"DRAW CMS preliminary Latex"<<endl;
	 double CentralsFontSize=0.035;
	 char text[200];
	 sprintf(text,"CMS preliminary");
	 if(!CMSprelim) sprintf(text,"CMS");
	 cout<<text<<endl;
	 TLatex *CentralsText1DM = new TLatex(xCentralsDM+xCentralsDM_shift,dimumassRap1->GetMaximum()*0.975,text);
	 CentralsText1DM->SetTextSize(CentralsFontSize);
	 CentralsText1DM->Draw( "same" );
	 sprintf(text,"L = 4.9 fb^{-1}");
	 cout<<text<<endl;
	 TLatex *CentralsText2DM = new TLatex(xCentralsDM+xCentralsDM_shift,dimumassRap1->GetMaximum()*0.825,text);
//	 if(!CMSprelim) CentralsText2DM = new TLatex(xCentralsDM+0.15,dimumassRap2->GetMaximum()*0.65,text);
	 CentralsText2DM->SetTextSize(CentralsFontSize);
	 CentralsText2DM->Draw( "same" );
	 sprintf(text,"pp    #sqrt{s} = 7 TeV");
	 cout<<text<<endl;
	 TLatex *CentralsText3DM = new TLatex(xCentralsDM+xCentralsDM_shift,dimumassRap1->GetMaximum()*0.9,text);
//	 if(!CMSprelim) CentralsText3DM = new TLatex(xCentralsDM+0.15,dimumassRap2->GetMaximum()*0.725,text);
	 CentralsText3DM->SetTextSize(CentralsFontSize);
	 CentralsText3DM->Draw( "same" );

	   char abcdef[200];
	   sprintf(abcdef,"a)");
	   TLatex *tex_abcdef = new TLatex(8.95,dimumassRap2->GetMaximum()*0.95,abcdef);
	   tex_abcdef->SetTextSize(CentralsFontSize*1.25);
//	   tex_abcdef->Draw( "same" );

	   sprintf(abcdef,"x 10^{3}");
	   tex_abcdef = new TLatex(8.5,dimumassRap1->GetMaximum()*1.07,abcdef);
	   	   tex_abcdef->SetTextSize(CentralsFontSize*1.05);
	   	   tex_abcdef->Draw( "same" );

	  c2->SetFrameBorderMode(0);
  	  sprintf(savename,"Figures/dimuonMass.pdf");
  	  c2->SaveAs(savename);
  	  sprintf(savename,"Figures/dimuonMass.C");
  	  c2->SaveAs(savename);
  	  sprintf(savename,"Figures/dimuonMass.jpg");
  	  c2->SaveAs(savename);



		TLegend* PTLegend;

		  for(int iState=1;iState<4;iState++){


	  double MassMin;
	  double MassMax;
	  if(iState==1){
	  MassMin=9.373;//mass1S-onia::nSigMass*sigma1S;
	  MassMax=9.526;//mass1S+onia::nSigMass*sigma1S;
	  // sigma 76.5 MeV
	  // mean 9.4495 GeV
	  }
	  if(iState==2){
	  MassMin=9.931;//mass2S-onia::nSigMass*sigma2S;
	  MassMax=10.093;//mass2S+onia::nSigMass*sigma2S;
	  // sigma 81 MeV
	  // mean 10.012 GeV
	  }
	  if(iState==3){
	  MassMin=10.260;//mass3S-onia::nSigMass*sigma3S;
	  MassMax=10.428;//mass3S+onia::nSigMass*sigma3S;
	  // sigma 84 MeV
	  // mean 10.344 GeV
	  }

		PTLegend=new TLegend(0.68,0.55,0.95,0.7);
		PTLegend->SetFillColor(0);
//		PTLegend->SetTextFont(72);
		PTLegend->SetTextSize(0.0345);
		PTLegend->SetBorderSize(1);
		char PTLegendEntry[200];


	  double PTplotmin=10;//10
	  double PTplotmax=79.99;//79.99
	  double PTbinwidth=0.5;//0.5
	  int PTnBins=(PTplotmax-PTplotmin)/PTbinwidth;

	  TH1D* ptRap1   = new TH1D( "ptRap1", "ptRap1", PTnBins, PTplotmin,PTplotmax);
	  TH1D* ptRap2   = new TH1D( "ptRap2", "ptRap2", PTnBins, PTplotmin,PTplotmax);

	  char PTDrawChar[2000];
	  char PTCutChar[2000];

	  rapMin=0.;
	  rapMax=0.6;
	  sprintf(PTDrawChar,"%s>>ptRap1",PtChar);
	  sprintf(PTCutChar,"%s>%f && %s<%f && %s>%f && %s<%f",AbsRapChar,rapMin,AbsRapChar,rapMax,DimuonMassChar,MassMin,DimuonMassChar,MassMax);

	  treeIn->Draw(PTDrawChar,PTCutChar);

	  rapMin=0.6;
	  rapMax=1.2;
	  sprintf(PTDrawChar,"%s>>ptRap2",PtChar);
	  sprintf(PTCutChar,"%s>%f && %s<%f && %s>%f && %s<%f",AbsRapChar,rapMin,AbsRapChar,rapMax,DimuonMassChar,MassMin,DimuonMassChar,MassMax);

	  treeIn->Draw(PTDrawChar,PTCutChar);

	  char PTyTitle[200];
	  sprintf(PTyTitle,"Counts per %d MeV",int(PTbinwidth*1000));
	  c2 = new TCanvas("c2","c2",1200,1100);
	  gStyle->SetPalette(1);
	  c2->SetFillColor(kWhite);
      gPad->SetLeftMargin(MargLeft);
      gPad->SetRightMargin(MargRight);
      gPad->SetTopMargin(MargTop);
	  c2->SetFrameBorderMode(0);

	  double histMin=0.5;
	  double histMax;
	  if(iState==1) histMax=2e4;
	  if(iState==2) histMax=2e4;
	  if(iState==3) histMax=2e4;

	   TH1F *PThist = new TH1F;
	   PThist = c2->DrawFrame(PTplotmin-4,histMin,PTplotmax,histMax);

	   PThist->SetXTitle("Dimuon #it{p}_{T} [GeV]");
	   PThist->SetYTitle(PTyTitle);
	   PThist->GetYaxis()->SetTitleOffset(1.25);

 	  gPad->SetFillColor(kWhite);
//      ptRap1->SetYTitle(PTyTitle);
//     ptRap1->SetXTitle("p_{T}(#mu#mu) [GeV]");
      ptRap1->SetTitle(0);
      ptRap1->SetStats(0);
//      ptRap1->GetYaxis()->SetTitleOffset(1.75);
      ptRap2->SetTitle(0);
      ptRap2->SetStats(0);
      ptRap1->SetMarkerStyle(25);
      ptRap1->SetMarkerColor(kBlack);
      ptRap1->SetLineColor(kBlack);
      ptRap2->SetMarkerStyle(20);
      ptRap2->SetMarkerColor(kRed);
      ptRap2->SetLineColor(kRed);
      ptRap2->SetMarkerSize(MarkerSizeTom);
      ptRap1->SetMarkerSize(MarkerSizeTom);



      ptRap1->Draw("same,E1");
	  ptRap2->Draw("same,E1");


		sprintf(PTLegendEntry,"|#it{y}| < 0.6");
		PTLegend->AddEntry(ptRap1,PTLegendEntry,"p");
		sprintf(PTLegendEntry,"0.6 < |#it{y}| < 1.2");
		PTLegend->AddEntry(ptRap2,PTLegendEntry,"p");


		TLine *PtLine;
		  for(int iPt=5;iPt<onia::kNbPTBins[1]+1;iPt++){
		  PtLine= new TLine( onia::pTRange[0][iPt], 0,  onia::pTRange[0][iPt] ,histMax);
		  PtLine->SetLineWidth( 2 );
		  PtLine->SetLineStyle( 2 );
		  PtLine->SetLineColor( kBlack );
		  PtLine->Draw();
		  }


		PTLegend->Draw();

		double xCentralsPt=55;
		double xCentralsPt_shift=0.;
		 cout<<"DRAW CMS preliminary Latex"<<endl;
	 double CentralsFontSize=0.035;
	 char text[200];
	 sprintf(text,"CMS preliminary");
	 if(!CMSprelim) sprintf(text,"CMS");
     cout<<text<<endl;
	 TLatex *CentralsText1Pt = new TLatex(xCentralsPt+xCentralsPt_shift,6e3,text);
	 CentralsText1Pt->SetTextSize(CentralsFontSize);
	 CentralsText1Pt->Draw( "same" );
	 sprintf(text,"L = 4.9 fb^{-1}");
	 cout<<text<<endl;
	 TLatex *CentralsText2Pt = new TLatex(xCentralsPt+xCentralsPt_shift,1.5e3,text);
//	 if(!CMSprelim) CentralsText2Pt = new TLatex(xCentralsPt+2.5,4e2,text);
	 CentralsText2Pt->SetTextSize(CentralsFontSize);
	 CentralsText2Pt->Draw( "same" );
	 sprintf(text,"pp    #sqrt{s} = 7 TeV");
	 cout<<text<<endl;
	 TLatex *CentralsText3Pt = new TLatex(xCentralsPt+xCentralsPt_shift,3e3,text);
//	 if(!CMSprelim) CentralsText3Pt = new TLatex(xCentralsPt+2.5,8e2,text);
	 CentralsText3Pt->SetTextSize(CentralsFontSize);
	 CentralsText3Pt->Draw( "same" );

	 sprintf(text,"#Upsilon(%dS)",iState);
	 cout<<text<<endl;
	 TLatex *CentralsText4Pt = new TLatex(37,6e3,text);
	 CentralsText4Pt->SetTextSize(CentralsFontSize);
	 CentralsText4Pt->Draw( "same" );


	   if(iState==1) sprintf(abcdef,"b)");
	   if(iState==2) sprintf(abcdef,"c)");
	   if(iState==3) sprintf(abcdef,"d)");
	   tex_abcdef = new TLatex(23.5,7e3,abcdef);
	   tex_abcdef->SetTextSize(CentralsFontSize*1.25);
//	   tex_abcdef->Draw( "same" );

  	  sprintf(savename,"Figures/pTdist_Ups%dS.pdf",iState);
  	  c2->SetLogy(true);
  	  c2->SaveAs(savename);
  	  sprintf(savename,"Figures/pTdist_Ups%dS.C",iState);
  	  c2->SaveAs(savename);
  	  sprintf(savename,"Figures/pTdist_Ups%dS.jpg",iState);
  	  c2->SaveAs(savename);

		  }



	  /// produce some more interesting plots

	  	  cout<<"Plot dimuon mass vs. pT / rap"<<endl;

	      int nHistBins=100;

	  	  TH2D* massrap   = new TH2D( "massrap", "massrap", nHistBins,8.6,11.4,104,0,1.3);
	  	  TH2D* masspt   = new TH2D( "masspt", "masspt", nHistBins,8.6,11.4,110,0,55);

	  	  treeIn->Draw("TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())):TMath::Sqrt(2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))>>masspt","","colz");
	  	  treeIn->Draw("TMath::Abs(1/2*TMath::Log((TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())+(lepP->Pz()+lepN->Pz())*(lepP->Pz()+lepN->Pz())+2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))+(lepP->Pz()+lepN->Pz()))/(TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())+(lepP->Pz()+lepN->Pz())*(lepP->Pz()+lepN->Pz())+2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))-(lepP->Pz()+lepN->Pz())))):TMath::Sqrt(2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))>>massrap","","colz");

	  	  c2 = new TCanvas("c2","c2",1200,1000);
	  	  massrap->SetYTitle("|y|(#mu#mu)");
	  	  massrap->SetXTitle("M(#mu#mu) [GeV]");
	  	  gStyle->SetPalette(1);
	  	  gPad->SetFillColor(kWhite);
	  	  massrap->SetTitle(0);
	  	  massrap->SetStats(0);
	      gPad->SetLeftMargin(0.15);
	      massrap->GetYaxis()->SetTitleOffset(1.5);
	  	  massrap->Draw("colz");

	  	  sprintf(savename,"Figures/massrap.pdf");
	  	  c2->SaveAs(savename);

	  	  c2 = new TCanvas("c2","c2",1200,1000);
	  	  masspt->SetYTitle("#it{p}_{T}(#mu#mu) [GeV]");
	  	  masspt->SetXTitle("M(#mu#mu) [GeV]");
	  	  gStyle->SetPalette(1);
	  	  gPad->SetFillColor(kWhite);
	  	  masspt->SetTitle(0);
	  	  masspt->SetStats(0);
	      gPad->SetLeftMargin(0.15);
	      masspt->GetYaxis()->SetTitleOffset(1.5);
	  	  masspt->Draw("colz");

		  c2->SetFrameBorderMode(0);
	  	  sprintf(savename,"Figures/masspt.pdf");
	  	  c2->SaveAs(savename);


	  	  char MasswindowChar[2000];
	      TH1F  *hMeanPt = new TH1F("hMeanPt","",nHistBins,8.6,11.4);
	      TH1F  *hMeanRap = new TH1F("hMeanRap","",nHistBins,8.6,11.4);
	      double massrange=11.4-8.6;

	      for(int i=1;i<nHistBins+1;i++){

	    	  TH1D  *hMeanPtBuffer=masspt->ProjectionY("hMeanPtBuffer", i, i,"e");
	    	  TH1D  *hMeanRapBuffer=massrap->ProjectionY("hMeanRapBuffer", i, i,"e");


	      hMeanPt->SetBinContent(i,hMeanPtBuffer->GetMean());
	      hMeanPt->SetBinError(i,hMeanPtBuffer->GetMeanError());
	      hMeanRap->SetBinContent(i,hMeanRapBuffer->GetMean());
	      hMeanRap->SetBinError(i,hMeanRapBuffer->GetMeanError());
	      }

	  	  c2 = new TCanvas("c2","c2",1200,1000);
	  	hMeanPt->SetYTitle("mean #it{p}_{T}(#mu#mu) [GeV]");
	  	hMeanPt->SetXTitle("M(#mu#mu) [GeV]");
	  	  gStyle->SetPalette(1);
	  	  gPad->SetFillColor(kWhite);
	  	hMeanPt->SetTitle(0);
	  	hMeanPt->SetStats(0);
	      gPad->SetLeftMargin(0.15);
	      gPad->SetRightMargin(0.15);
	      hMeanPt->GetYaxis()->SetTitleOffset(1.5);
	  	hMeanPt->Draw("E1");

		  c2->SetFrameBorderMode(0);
	  	  sprintf(savename,"Figures/meanpt_vs_mass.pdf");
	  	  c2->SaveAs(savename);

	  	  c2 = new TCanvas("c2","c2",1200,1000);
	  	hMeanRap->SetYTitle("mean |y|(#mu#mu)");
	  	hMeanRap->SetXTitle("M(#mu#mu) [GeV]");
	  	  gStyle->SetPalette(1);
	  	  gPad->SetFillColor(kWhite);
	  	hMeanRap->SetTitle(0);
	  	hMeanRap->SetStats(0);
	      gPad->SetLeftMargin(0.15);
	      gPad->SetRightMargin(0.15);
	      hMeanRap->GetYaxis()->SetTitleOffset(1.5);
	  	hMeanRap->Draw("E1");

		  c2->SetFrameBorderMode(0);
	  	  sprintf(savename,"Figures/meanrap_vs_mass.pdf");
	  	  c2->SaveAs(savename);

	  for(int iState=1;iState<4;iState++){

  TH2D* rapPt   = new TH2D( "rapPt", "rapPt", 52,-1.3,1.3,110,0,55);

  double MassMin;
  double MassMax;
  if(iState==1){
  MassMin=9.373;//mass1S-onia::nSigMass*sigma1S;
  MassMax=9.526;//mass1S+onia::nSigMass*sigma1S;
  // sigma 76.5 MeV
  // mean 9.4495 GeV
  }
  if(iState==2){
  MassMin=9.931;//mass2S-onia::nSigMass*sigma2S;
  MassMax=10.093;//mass2S+onia::nSigMass*sigma2S;
  // sigma 81 MeV
  // mean 10.012 GeV
  }
  if(iState==3){
  MassMin=10.260;//mass3S-onia::nSigMass*sigma3S;
  MassMax=10.428;//mass3S+onia::nSigMass*sigma3S;
  // sigma 84 MeV
  // mean 10.344 GeV
  }

  cout<<"Plotting rap-Pt for Ups"<<iState<<"S"<<endl;
  cout<<"MassMin for rap-Pt plot = "<<MassMin<<endl;
  cout<<"MassMax for rap-Pt plot = "<<MassMax<<endl;

  sprintf(CutChar,"TMath::Sqrt(2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))>%f&&TMath::Sqrt(2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))<%f",MassMin,MassMax);

  treeIn->Draw("TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())):1/2*TMath::Log((TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())+(lepP->Pz()+lepN->Pz())*(lepP->Pz()+lepN->Pz())+2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))+(lepP->Pz()+lepN->Pz()))/(TMath::Sqrt((lepP->Px()+lepN->Px())*(lepP->Px()+lepN->Px())+(lepP->Py()+lepN->Py())*(lepP->Py()+lepN->Py())+(lepP->Pz()+lepN->Pz())*(lepP->Pz()+lepN->Pz())+2*0.105658*0.105658+2*(lepP->Energy()*lepN->Energy()-(lepP->Px()*lepN->Px()+lepP->Py()*lepN->Py()+lepP->Pz()*lepN->Pz())))-(lepP->Pz()+lepN->Pz())))>>rapPt",CutChar,"colz");

  c2 = new TCanvas("c2","c2",1200,1500);
  rapPt->SetYTitle("#it{p}_{T}(#mu#mu) [GeV]");
  rapPt->SetXTitle("y(#mu#mu)");
  gStyle->SetPalette(1);
//  gPad->SetTopMargin(0.1);
  gPad->SetFillColor(kWhite);
  rapPt->SetTitle(0);
  rapPt->SetStats(0);
  gPad->SetLeftMargin(0.15);
  rapPt->GetYaxis()->SetTitleOffset(1.5);
  rapPt->Draw("colz");

  TLine* rapPtLine;

  for(int iRap=0;iRap<onia::kNbRapForPTBins+1;iRap++){
	  rapPtLine= new TLine( -onia::rapForPTRange[iRap], onia::pTRange[0][0], -onia::rapForPTRange[iRap], onia::pTRange[0][onia::kNbPTBins[iRap]] );
	  rapPtLine->SetLineWidth( 2 );
	  rapPtLine->SetLineStyle( 1 );
	  rapPtLine->SetLineColor( kWhite );
	  rapPtLine->Draw();
	  rapPtLine= new TLine( onia::rapForPTRange[iRap], onia::pTRange[0][0], onia::rapForPTRange[iRap], onia::pTRange[0][onia::kNbPTBins[iRap]] );
	  rapPtLine->SetLineWidth( 2 );
	  rapPtLine->SetLineStyle( 1 );
	  rapPtLine->SetLineColor( kWhite );
	  rapPtLine->Draw();

  for(int iPt=0;iPt<onia::kNbPTBins[iRap+1]+1;iPt++){
  rapPtLine= new TLine( -onia::rapForPTRange[onia::kNbRapForPTBins], onia::pTRange[0][iPt], onia::rapForPTRange[onia::kNbRapForPTBins], onia::pTRange[0][iPt] );
  rapPtLine->SetLineWidth( 2 );
  rapPtLine->SetLineStyle( 1 );
  rapPtLine->SetLineColor( kWhite );
  rapPtLine->Draw();
  }}
  sprintf(savename,"Figures/rapPt_Ups%dS.pdf",iState);
  c2->SaveAs(savename);

  delete rapPt;
  delete rapPtLine;
	  }


	  delete massrap;
	  delete masspt;
	  delete c2;


  }



  Double_t onia_mass;
  Int_t index;
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 100000 == 0)
      cout << "entry " << iEntry << " out of " << treeIn->GetEntries() << endl;



    double ProjectPtMin;
    double ProjectPtMax;
    double ProjectRapMin;
    double ProjectRapMax;

    if(iPTBin>0){
    	ProjectPtMin=onia::pTRange[iRapBin][iPTBin-1];
    	ProjectPtMax=onia::pTRange[iRapBin][iPTBin];
    }
    if(iPTBin==0){
    	ProjectPtMin=onia::pTRange[iRapBin][5];//[0]
    	ProjectPtMax=onia::pTRange[iRapBin][onia::kNbPTBins[iRapBin]];
    }
    if(iRapBin>0){
    	ProjectRapMin=onia::rapForPTRange[iRapBin-1];
    	ProjectRapMax=onia::rapForPTRange[iRapBin];
    }
    if(iRapBin==0){
    	ProjectRapMin=onia::rapForPTRange[0];
    	ProjectRapMax=onia::rapForPTRange[onia::kNbRapForPTBins];
    }

    if(UpsMC&&iPTBin<6){
    	ProjectPtMin=onia::pTRange[iRapBin][5];//[0]
    	ProjectPtMax=onia::pTRange[iRapBin][onia::kNbPTBins[iRapBin]];
    }

    *onia = *(lepP) + *(lepN);
    if(onia->Pt() > ProjectPtMin && onia->Pt() < ProjectPtMax &&
       TMath::Abs(onia->Rapidity()) > ProjectRapMin && TMath::Abs(onia->Rapidity()) < ProjectRapMax){

      treeOut->Fill(); //stores TLorenzVectors of the two muons in the given pT and rap cell

      if(UpsMC){
      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    		hBG_cosThetaPhi[iFrame][L]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
    		hBG_cosThetaPhi[iFrame][R]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
      }
	  }

      //store now the cosTheta and phi distributions of the BG:
      onia_mass = onia->M();
      if(onia_mass > massMin[L] && onia_mass < massMax[L])
	index = L;
      else if(onia_mass > massMin[R] && onia_mass < massMax[R])
	index = R;
      else 
	continue;

      calcPol(*lepP, *lepN);

      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
    	  if(!UpsMC) hBG_cosThetaPhi[iFrame][index]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
    }



  }
  fOut->cd();

  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
    hBG_cosThetaPhi[iFrame][L]->Write();
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
    hBG_cosThetaPhi[iFrame][R]->Write();
  treeOut->Write();  
  fOut->Close();
}
