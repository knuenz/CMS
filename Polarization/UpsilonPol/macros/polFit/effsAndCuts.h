

enum { loose, tight };

bool isMuonInAcceptance(int iCut, double pT, double eta){

	//iCut=x correspinding to FidCuts=x+1 in scripts
	//iCut=-1 (FidCuts=0) - no cuts -> decision is always true
	//iCut=0  (FidCuts=1) - LOOSE cuts
	//iCut=1  (FidCuts=2) - TIGHT cuts
	//iCut=2  (FidCuts=3) - simple 6 GeV cut (including |eta|<2.4 cut)
	//iCut=3  (FidCuts=4) - simple 4.6 GeV cut (including |eta|<2.4 cut)
	//iCut=4  (FidCuts=5) - for |eta|<0.3: 6 GeV cut, for |eta|>0.3: 4.6 GeV cut (including |eta|<2.4 cut)

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
		  else if(TMath::Abs(eta) > etaBorderHLT[iCut][3])
			  minPT_HLT = 1000.; //reject all events with |eta| > 2.4 (or 2.2, ...)
	  }

	  if(pT > minPT_HLT)
		  decision = kTRUE;

  }

  if(iCut==2 && pT>6 && TMath::Abs(eta) < 2.4) decision = kTRUE;

  if(iCut==3 && pT>4.6 && TMath::Abs(eta) < 2.4) decision = kTRUE;

  if(iCut==4) {
	  if(TMath::Abs(eta) < 0.3 && pT>6) decision = kTRUE;
	  if(TMath::Abs(eta) > 0.3 && pT>4.6 && TMath::Abs(eta) < 2.4) decision = kTRUE;
  }

  return decision;
}


void EvaluateEffFileName(int nEff, char EffFileName [200], bool singleLeptonEff) {

	if(singleLeptonEff) {

		if(nEff==101) sprintf(EffFileName,"EfficiencyProductDimuon0Jpsi_TrkCuts_20Oct2011.root");

	}

	if(!singleLeptonEff) {

		if(nEff==201) sprintf(EffFileName,"DimuVtxModule_onePair_pair_pt_pair_absrapidity_prescaleModule_Tr.root");

	}


}


double DiLeptonEfficiency( double& pT, double& rap, int nEff, TH1* hEff) {


	double eff=1;
	if(nEff==1) return eff;
	if(nEff==2) return eff/2.;

	if(nEff>100){
	  Int_t binX = hEff->GetXaxis()->FindBin(TMath::Abs(rap));
	  Int_t binY = hEff->GetYaxis()->FindBin(pT);
	  eff = hEff->GetBinContent(binX, binY);
	  return eff;
	}


}

double singleLeptonEfficiency( double& pT, double& eta, int nEff, TH1* hEff) {

	//nEff=1 all muons have eff=1
	//nEff=2 see algorithm below
	//nEff=3 Matt's efficiencies
	//nEff=4 I suggest this to be the final efficiencies from the TnP studies, not yet implemented

	double eff;

	if(nEff>100){
	  Int_t binX = hEff->GetXaxis()->FindBin(TMath::Abs(eta));
	  Int_t binY = hEff->GetYaxis()->FindBin(pT);
	  eff = hEff->GetBinContent(binX, binY);
//	  printf("%f",eff);
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
