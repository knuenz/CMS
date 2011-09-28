

enum { loose, tight };

bool isMuonInAcceptance(int iCut, double pT, double eta){

	//iCut=x correspinding to FidCuts=x+1 in scripts
	//iCut=-1 (FidCuts=0) - no cuts -> decision is always true
	//iCut=0  (FidCuts=1) - LOOSE cuts
	//iCut=1  (FidCuts=2) - TIGHT cuts
	//iCut=2  (FidCuts=3) - simple 6 GeV cut (including |eta|<2.4 cut)
	//iCut=3  (FidCuts=4) - simple 4.6 GeV cut (including |eta|<2.4 cut)
	//iCut=4  (FidCuts=5) - for |eta|<0.3: 6 GeV cut, for |eta|>0.3: 4.6 GeV cut (including |eta|<2.4 cut)

  double etaBorderHLT[2][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2}}; //LOOSE, TIGHT cuts
  double pTBorderHLT[2][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0}};

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




double singleLeptonEfficiency( double& pT, double& eta, int nEff) {

	//nEff=1 all muons have eff=1
	//nEff=2 see algorithm below
	//nEff=3 Matt's efficiencies
	//nEff=4 I suggest this to be the final efficiencies from the TnP studies, not yet implemented

	double eff;

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

//
//
//

  if(nEff==4) return eff;

}
