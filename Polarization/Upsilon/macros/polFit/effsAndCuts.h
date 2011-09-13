

enum { loose, tight };

//int fidCut = loose;

// Hermine's fiducial cuts

bool isMuonInAcceptance(int iCut, double pT, double eta){

  double etaBorderHLT[4][4] = {{0., 1.1, 1.4, 2.4}, {0., 1.2, 1.3, 2.2},{0., 1.1, 1.4, 2.4},{0., 0.75, 1.5, 2.25}}; //LOOSE, TIGHT cuts
  double pTBorderHLT[4][4] = {{4.6, 4.0, 2.75, 2.0}, {5.2, 4.7, 3.3, 3.0},{6.,6.,6.,6.},{9.,5.,3.,1.}};

  double minPT_HLT;
  bool decision = kTRUE;

	if(iCut==-1) return decision;
	decision = kFALSE;

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

  return decision;
}

/*
// Matt's fiducial cuts
bool isMuonInAcceptance(int iCut, double pT, double eta){

   bool decision = kFALSE;
   if ( (   TMath::Abs(eta) < 0.8 && pT > 3.75
          || TMath::Abs(eta) > 0.8 && TMath::Abs(eta) < 1.6 && pT > 3.5

          || TMath::Abs(eta) > 1.6 && TMath::Abs(eta) < 2.4 && pT > 3.0  )
      ) decision = kTRUE;

  return decision;
}
*/

/// Matt's single-muon efficiencies
double singleLeptonEfficiency( double& pT, double& eta, int nEff) {

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

}

// simple efficiency
/*
double singleLeptonEfficiency( double& pT, double& eta ) {

  const double mu_pT_min = 5.0; // point of 50% efficiency
  const double smoothcutpar = 3.0; // smaller -> smoother
  return 1. / ( 1. + exp(-smoothcutpar*(pT - mu_pT_min)) );
}
*/
