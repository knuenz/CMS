#include "Riostream.h"
#include "TEfficiency.h"


enum { loose, tight };

bool isMuonInAcceptance(int iCut, double pT, double eta){

	//iCut=x correspinding to FidCuts=x+1 in scripts
	//iCut=-1 (FidCuts=0) - no cuts -> decision is always true
	//iCut=0  (FidCuts=1) - LOOSE cuts
	//iCut=1  (FidCuts=2) - TIGHT cuts
	//iCut=2  (FidCuts=3) - LOOSE cuts with wheel cut
	//iCut=3  (FidCuts=4) - TIGHT cuts with wheel cut
	//iCut=6  (FidCuts=7) - Matts Cuts (AN2011-130)
	//iCut=7  (FidCuts=8) - New LOOSE cuts Feb12
	//iCut=8  (FidCuts=9) - New TIGHT cuts Feb12
	//iCut=9  (FidCuts=10) - New TIGHT cuts Feb12, without an eta<1.6 cut
	//iCut=10  (FidCuts=11) - TPV cuts
	//iCut=11  (FidCuts=12) - TPV cuts + 500 MeV
	//iCut=12  (FidCuts=13) - TPV cuts + 0.2-0.3 eta cut
	//iCut=13  (FidCuts=14) - Chib cuts

	  Double_t etaBorderHLT[2][4] = {{0., 1.2, 1.6, 2.1}, {0., 1.2, 1.6, 2.1}}; //LOOSE, TIGHT cuts, Tracker muons
	  Double_t pTBorderHLT[2][4] = {{3.5, 3.5, 2.0, 2.0}, {3.8, 3.8, 2.0, 2.0}};
//  double etaBorderHLT[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts, Global Muons
//  double pTBorderHLT[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};
//  double etaBorderHLT[2][4] = {{0., 0.8, 1.2, 2.1}, {0., 0.8, 1.2, 2.1}}; //LOOSE, TIGHT cuts, Jozko Design
//  double pTBorderHLT[2][4] = {{4.6, 3.5, 2.75, 2.0}, {5.2, 4.0, 2.75, 2.0}};

  double minPT_HLT;
  bool decision = kTRUE;

  if(iCut==-1) return decision;
  decision = kFALSE;

  if(iCut==0 || iCut==1){
	  //loop over higher pT muon
	  for(int iEta = 0; iEta < 3; iEta++){
		  if(TMath::Abs(eta) > etaBorderHLT[iCut][iEta] && TMath::Abs(eta) < etaBorderHLT[iCut][iEta+1]){
			  minPT_HLT = (pTBorderHLT[iCut][iEta+1]-pTBorderHLT[iCut][iEta]) / (etaBorderHLT[iCut][iEta+1]-etaBorderHLT[iCut][iEta]) * (TMath::Abs(eta) - etaBorderHLT[iCut][iEta]) + pTBorderHLT[iCut][iEta];
			  break;
		 }
	  }
                  if(TMath::Abs(eta) > 1.6)
                          minPT_HLT = 1000.; //reject all events with |eta| > 1.6 (or 2.2, ...)

	  if(pT > minPT_HLT)
		  decision = kTRUE;
	  if(pT < 2.5)
		  decision = kFALSE;
  }

  if(iCut==2 || iCut==3){
          //loop over higher pT muon
          for(int iEta = 0; iEta < 3; iEta++){
                  if(TMath::Abs(eta) > etaBorderHLT[iCut-2][iEta] && TMath::Abs(eta) < etaBorderHLT[iCut-2][iEta+1]){
                          minPT_HLT = (pTBorderHLT[iCut-2][iEta+1]-pTBorderHLT[iCut-2][iEta]) / (etaBorderHLT[iCut-2][iEta+1]-etaBorderHLT[iCut-2][iEta]) * (TMath::Abs(eta) - etaBorderHLT[iCut-2][iEta]) + pTBorderHLT[iCut-2][iEta];
                          break;
                  }
           }
                           if(TMath::Abs(eta) > 1.2)
                          minPT_HLT = 1000.; //reject all events with |eta| > 1.6 (or 2.2, ...)
 
          if(pT > minPT_HLT)
                  decision = kTRUE;
	  if(TMath::Abs(eta) < 0.3 && TMath::Abs(eta) > 0.2) decision = kFALSE;
	  if(pT < 2.5)
		  decision = kFALSE;
                  
  }             

  if(iCut==6){
	  if(TMath::Abs(eta)<0.8 && pT>3.75) decision=kTRUE;
	  if(TMath::Abs(eta)>0.8 && TMath::Abs(eta)<1.6 && pT>3.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.6 && TMath::Abs(eta)<2.4 && pT>3.0) decision=kTRUE;
  }

  if(iCut==7){
	  if(TMath::Abs(eta)<1.2 && pT>4.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<1.4 && pT>4.) decision=kTRUE;
	  if(TMath::Abs(eta)>1.4 && TMath::Abs(eta)<1.6 && pT>3.) decision=kTRUE;
  }

  if(iCut==8){
	  if(TMath::Abs(eta)<1.4 && pT>4.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.4 && TMath::Abs(eta)<1.6 && pT>3.5) decision=kTRUE;
  }

  if(iCut==9){
	  if(TMath::Abs(eta)<1.4 && pT>4.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.4 && pT>3.5) decision=kTRUE;
  }

  if(iCut==10){
	  if(TMath::Abs(eta)<1.2 && pT>4.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<1.4 && pT>3.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.4 && TMath::Abs(eta)<1.6 && pT>3.) decision=kTRUE;
  }

  if(iCut==11){
	  if(TMath::Abs(eta)<1.2 && pT>5.0) decision=kTRUE;
	  if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<1.4 && pT>4.0) decision=kTRUE;
	  if(TMath::Abs(eta)>1.4 && TMath::Abs(eta)<1.6 && pT>3.5) decision=kTRUE;
  }

  if(iCut==12){
	  if(TMath::Abs(eta)<1.2 && pT>4.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.2 && TMath::Abs(eta)<1.4 && pT>3.5) decision=kTRUE;
	  if(TMath::Abs(eta)>1.4 && TMath::Abs(eta)<1.6 && pT>3.) decision=kTRUE;
	  if(TMath::Abs(eta)<0.3 && TMath::Abs(eta)>0.2) decision=kFALSE;
  }

  if(iCut==13){
	  if(TMath::Abs(eta)<1.3 && pT>3.3) decision=kTRUE;
	  if(TMath::Abs(eta)>1.3 && TMath::Abs(eta)<2.2 && pT>2.9) decision=kTRUE;
  }

  return decision;
}


void EvaluateEffFileName(int nEff, char EffFileName [200], bool singleLeptonEff) {

	if(singleLeptonEff) {

		if(nEff==101 || nEff==1001 || nEff==10001) sprintf(EffFileName,"EfficiencyProductDimuon0Jpsi_TrkCuts_2Nov2011.root");
		if(nEff==102 || nEff==1002 || nEff==10002) sprintf(EffFileName,"EfficiencyProductDimuon0Jpsi_TrkCuts_2Nov2011_handCorrected.root");
		if(nEff==103 || nEff==1003 || nEff==10003) sprintf(EffFileName,"EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Tracker80.root");
		if(nEff==104) sprintf(EffFileName,"ScaledMCEff2_EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Tracker80.root");
		if(nEff==105) sprintf(EffFileName,"singleMuTruthEff_18Jan2012_40GeVrap1_2pT100GeV_EtaCut_FineBins200MeV.root");
		if(nEff==106) sprintf(EffFileName,"singleMuTruthEff_1Feb2012_40GeVrap1_2pT100GeV_EtaCut_FineBins200MeV_L1DoubleMu0_HighQ.root");
		if(nEff==107) sprintf(EffFileName,"EfficiencyProductDimuon0Jpsi_MuonID-MC_MuonQualRunA_L1L2L3Run1_Trk80Cuts_29Jan2012.root");
		if(nEff==108 || nEff==1008 || nEff==10008) sprintf(EffFileName,"singleMuonEff_combinedMC_fineBins_Trk80Cuts_3Feb2011.root");
		if(nEff==109 || nEff==1009 || nEff==10009) sprintf(EffFileName,"EfficiencyProduct_combinedMC_fineBins_Trk80Cuts_3Feb2011.root");
		if(nEff==110 || nEff==1010 || nEff==10010) sprintf(EffFileName,"EfficiencyProductDimuon0Jpsi_combinedMC_fineBins_Trk80Cuts_6Mar2011.root");
		if(nEff==111 || nEff==1011 || nEff==10011) sprintf(EffFileName,"EfficiencyFactorized_Dimuon0Jpsi_combined_DATA_MC_Trk80Cuts_14Mar2012.root");

		if(nEff==1020) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_Central.root");
		if(nEff==1021) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_pTshift_plus.root");
		if(nEff==1022) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_pTshift_minus.root");
		if(nEff==1023) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_pTscale_plus.root");
		if(nEff==1024) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_pTscale_minus.root");
		if(nEff==1025) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_effshift_plus.root");
		if(nEff==1026) sprintf(EffFileName,"ParametrizedFactDataEff_16Mar_effshift_minus.root");

		if(nEff==1030) sprintf(EffFileName,"ParametrizedFactDataEff_May19_Central.root");
		if(nEff==1031) sprintf(EffFileName,"ParametrizedFactDataEff_May19_pTshift_plus.root");
		if(nEff==1032) sprintf(EffFileName,"ParametrizedFactDataEff_May19_pTshift_minus.root");
		if(nEff==1033) sprintf(EffFileName,"ParametrizedFactDataEff_May19_pTscale_plus.root");
		if(nEff==1034) sprintf(EffFileName,"ParametrizedFactDataEff_May19_pTscale_minus.root");
		if(nEff==1035) sprintf(EffFileName,"ParametrizedFactDataEff_May19_effshift_plus.root");
		if(nEff==1036) sprintf(EffFileName,"ParametrizedFactDataEff_May19_effshift_minus.root");

	}

	if(!singleLeptonEff) {

		if(nEff==201) sprintf(EffFileName,"DimuonVtxEff_Dimuon0Jpsi_cosTheta_Phi_TrkCuts80_CS_01Dec2011.root");
		if(nEff==211) sprintf(EffFileName,"DimuonVtxEff_Dimuon0Jpsi_cosTheta_Phi_TrkCuts80_HX_01Dec2011.root");
		if(nEff==222) sprintf(EffFileName,"DimuVtxModuleDimuon10JpsiBarrel_onePair_cosTheta_phi_PHX_seagulls_Trk80Cuts_22March2011_corrected.root");

		if(nEff==301 || nEff==311 || nEff==321) sprintf(EffFileName,"rhoFactor_Ups1S_MCTruthEff_22Jan2011.root");
		if(nEff==302 || nEff==312 || nEff==322) sprintf(EffFileName,"rhoFactor_Ups1S_MCTruthEff_23Jan2011.root");
		if(nEff==303 || nEff==313 || nEff==323) sprintf(EffFileName,"rhoFactor_Ups1S_SingleMuEff_7Feb2012.root");
		if(nEff==304 || nEff==314 || nEff==324) sprintf(EffFileName,"rhoFactor_Ups1S_ProdSingleMuEff_7Feb2012.root");
		if(nEff==305 || nEff==315 || nEff==325) sprintf(EffFileName,"rhoFactor_Ups1S_MCTruthEff_21March2012_finalPTBins_TPVCuts.root");

	}


}

double EvaluateRhoFactor( double& costh, double& phi, int nEff, TFile* fInRhoFactor, double rap, double pT) {

	double eff=1;
	if(nEff==1) return eff;

	int pTbin;
	int rapBin;
	const int nRhoPtBins=12;
	const int nRhorapBins=2;
	Double_t rapRangeRho[nRhorapBins+1] = {0.,0.6,1.2};
	Double_t pTRangeRho[nRhoPtBins+1] = {5., 6., 7., 8., 9., 10., 12., 16., 20., 30., 50.};

	if(pT>pTRangeRho[nRhoPtBins]){eff=0;return eff;}

	for(int i=1;i<nRhoPtBins+1;i++){
//		cout<<i;
		if(pT>pTRangeRho[i-1]&&pT<pTRangeRho[i]) {pTbin=i; break; }
	}
	for(int i=1;i<nRhorapBins+1;i++){
		if(TMath::Abs(rap)>rapRangeRho[i-1]&&TMath::Abs(rap)<rapRangeRho[i]) {rapBin=i; break; }
	}
//	cout<<"rap "<<rap<<" bin "<<rapBin<<" pT "<<pT<<" bin "<<pTbin<<endl;


	if(nEff>300){
		char EffType[200];
		if(nEff>300 && nEff<311) sprintf(EffType,"hRho_pol_totEff_CS_rap%d_pT%d",rapBin,pTbin);
		if(nEff>310 && nEff<321) sprintf(EffType,"hRho_pol_totEff_HX_rap%d_pT%d",rapBin,pTbin);
		if(nEff>320 && nEff<331) sprintf(EffType,"hRho_pol_totEff_PHX_rap%d_pT%d",rapBin,pTbin);

	  TH1* hEff=(TH1*) fInRhoFactor->Get(EffType);
	  Int_t binX = hEff->GetXaxis()->FindBin(TMath::Abs(costh));
	  Int_t binY = hEff->GetYaxis()->FindBin(phi);
	  eff = hEff->GetBinContent(binX, binY);
//	  cout<<eff<<endl;
	  return eff;
	}


}


double DiLeptonEfficiency( double& costh, double& phi, int nEff, TFile* fInDileptonEff, bool MCeff) {


	double eff=1;
	if(nEff==1) return eff;
	if(nEff==2) return eff/2.;


	if(nEff>200){
		char EffType[200];
		if(MCeff) sprintf(EffType,"hEff_MC_central");
		else sprintf(EffType,"hEff_DATA_central");

	  TH1* hEff=(TH1*) fInDileptonEff->Get(EffType);
	  Int_t binX = hEff->GetXaxis()->FindBin(TMath::Abs(costh));
	  Int_t binY = hEff->GetYaxis()->FindBin(phi);
	  eff = hEff->GetBinContent(binX, binY);
	  return eff;
	}


}

double singleLeptonEfficiency( double& pT, double& eta, int nEff, TFile* fInEff, TH2D* hEvalEff,bool MCeff, TEfficiency* TEff) {

	//nEff=1 all muons have eff=1
	//nEff=2 see algorithm below
	//nEff=3 Matt's efficiencies
	//nEff=4 I suggest this to be the final efficiencies from the TnP studies, not yet implemented

	double eff;
	char EffType[200];

	if(nEff==105 || nEff==106){
		sprintf(EffType,"totEff_MCTRUTH_pT_eta");


/*		  int binX = hEvalEff->GetXaxis()->FindBin(TMath::Abs(eta));
		  int binY = hEvalEff->GetYaxis()->FindBin(pT);

		  double effHist = hEvalEff->GetBinContent(binX, binY);
//		  cout<<"HistEff "<<eff<<endl;
*/
		  Int_t globalBin = TEff->FindFixBin(TMath::Abs(eta), pT);
		  eff = TEff->GetEfficiency(globalBin);


		return eff;
	}

	if(MCeff) sprintf(EffType,"hEff_MC_central");
	else sprintf(EffType,"hEff_DATA_central");

	if(nEff > 100 && nEff < 1000){//binned efficiencies: 1XX

	  TH1* hEff=(TH1*) fInEff->Get(EffType);
	  Int_t binX = hEff->GetXaxis()->FindBin(TMath::Abs(eta));
	  Int_t binY = hEff->GetYaxis()->FindBin(pT);

	  eff = hEff->GetBinContent(binX, binY);

	  return eff;
	}


	if(nEff > 1000){//linear interpolated efficiencies: 1XXnEff (1D) 1XXXnEff (2D)

	  int binX = hEvalEff->GetXaxis()->FindBin(TMath::Abs(eta));
	  int binY = hEvalEff->GetYaxis()->FindBin(pT);

	  eff = hEvalEff->GetBinContent(binX, binY);

  	  return eff;
	}


  if(nEff==1) return 1;

  const double mu_pT_min = 4.0;
  const double smoothcutpar = 3.0;
  eff = 1. / ( 1. + exp(-smoothcutpar*(pT - mu_pT_min)) );
  if ( TMath::Abs(eta) > 2.4 ) eff = 0.;

  if(nEff==2) return eff;

  double c0 = 0.878;
  double c1 = 3.894;
  double c2 = 0.957;

  if ( TMath::Abs(eta) > 0.8 && TMath::Abs(eta) < 1.2 ) {
    c0 = 0.839;
    c1 = 3.860;
    c2 = 0.512;
  }

  if ( TMath::Abs(eta) > 1.2 && TMath::Abs(eta) < 1.6 ) {
    c0 = 0.882;
    c1 = 2.984;
    c2 = 0.405;
  }

  if ( TMath::Abs(eta) > 1.6 && TMath::Abs(eta) < 2.0 ) {
    c0 = 0.839;
    c1 = 2.280;
    c2 = 1.398;
  }

  if ( TMath::Abs(eta) > 2.0 && TMath::Abs(eta) < 2.4 ) {
    c0 = 0.713;
    c1 = 0.0;
    c2 = 0.0;
  }

  eff = 0.5*c0 * (1. + TMath::Erf( (pT - c1) / ( sqrt(2.) * c2) ) );

  if(nEff==3) return eff;


//implement our TnP efficiencies here, and calculate it as eff:

  double b = -2.821;
  double k =  1.341;
  double a = -0.098;


  if ( TMath::Abs(eta) > 0.2 && TMath::Abs(eta) < 0.3 ) {
  b = -0.644;
  k =  3.460;
  a = -0.145;
  }

  if ( TMath::Abs(eta) > 0.3 && TMath::Abs(eta) < 0.8 ) {
  b = -0.551;
  k =  3.260;
  a = -0.071;
  }

  if ( TMath::Abs(eta) > 0.8 && TMath::Abs(eta) < 1.2 ) {
  b = -0.006;
  k =  4.884;
  a = -0.098;
  }

  if ( TMath::Abs(eta) > 1.2 && TMath::Abs(eta) < 1.4 ) {
  b = -0.942;
  k =  1.778;
  a = -0.117;
  }

  if ( TMath::Abs(eta) > 1.4 && TMath::Abs(eta) < 2.1 ) {
  b = -0.834;
  k =  1.908;
  a = -0.111;
  }

  eff = a + TMath::Erf( b + pT/k );

  if ( TMath::Abs(eta) > 2.1 ) eff = 0.;

  if(nEff==4) return eff;


//Eff Ilse 21.October 2011

  b = -1.9598;
  k =  1.6968;
  a = -0.0860;


  if ( TMath::Abs(eta) > 0.2 && TMath::Abs(eta) < 0.3 ) {
  b = -0.0153;
  k =  6.0333;
  a = -0.1160;
  }

  if ( TMath::Abs(eta) > 0.3 && TMath::Abs(eta) < 0.8 ) {
  b =  -0.7084;
  k =  2.8787;
  a = -0.0629;
  }

  if ( TMath::Abs(eta) > 0.8 && TMath::Abs(eta) < 1.2 ) {
  b = -0.2856;
  k =  3.7530;
  a = -0.0839;
  }

  if ( TMath::Abs(eta) > 1.2 && TMath::Abs(eta) < 1.4 ) {
  b =  -0.2846;
  k =  2.5955;
  a = -0.0773;
  }

  if ( TMath::Abs(eta) > 1.4 && TMath::Abs(eta) < 2.1 ) {
  b = -0.7710;
  k =  1.9250;
  a = -0.0848;
  }

  eff = a + TMath::Erf( b + pT/k );

  if ( TMath::Abs(eta) > 2.1 ) eff = 0.;

  if(nEff==5) return eff;


//Eff Ilse 20.October 2011, Trk Muon Cuts


  b = -1.9803;
  k =  1.6726;
  a = -0.0890;


  if ( TMath::Abs(eta) > 0.2 && TMath::Abs(eta) < 0.3 ) {
  b = -1.7473;
  k =  1.9034;
  a = -0.1741;
  }

  if ( TMath::Abs(eta) > 0.3 && TMath::Abs(eta) < 0.6 ) {
  b = -1.8522;
  k =  1.7584;
  a = -0.0662;
  }

  if ( TMath::Abs(eta) > 0.6 && TMath::Abs(eta) < 0.8 ) {
  b = -2.3682;
  k =  1.4189;
  a = -0.0694;
  }

  if ( TMath::Abs(eta) > 0.8 && TMath::Abs(eta) < 1.2 ) {
  b = -1.4904;
  k =  1.8760;
  a = -0.1203;
  }

  if ( TMath::Abs(eta) > 1.2 && TMath::Abs(eta) < 1.6 ) {
  b = -1.7625;
  k =  1.2212;
  a = -0.0832;
  }

  if ( TMath::Abs(eta) > 1.6 && TMath::Abs(eta) < 2.1 ) {
  b = -0.9992;
  k =  1.6830;
  a = -0.1081;
  }

  eff = a + TMath::Erf( b + pT/k );

  if ( TMath::Abs(eta) > 2.1 ) eff = 0.;

  if(nEff==6) return eff;


}