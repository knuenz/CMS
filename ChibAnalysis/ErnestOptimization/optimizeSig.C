//using namespace RooFit;

int optimizeSig(TString varset="", TString what="chib"){

  TString cuts;

  TFile *f;

  if (what=="chib") f = new TFile("/afs/cern.ch/work/k/knuenz/CMSSW_4_2_3/src/ErnestOptimization/rooDS_hltbest.root");
  else if (what=="chic") f = new TFile("/afs/cern.ch/work/e/eaguiloc/CMSSW_4_2_3/src/TTree_ChicMay03.root");
  else return 0;

  TTree *t = (TTree*)f->Get("atree");

  TString vars[30];
  vector<float> varcuts[30];
  float inicuts[30];
  int nvars;

  float margin;
  int maxcount,nbins;

  ifstream in("variables"+varset+".txt");
  in >> margin >> maxcount >> nbins;
  in >> cuts;
  in >> nvars;

  for(int i=0 ; i<nvars && i<30 ; i++) {
    float min,max;
    in >> vars[i] >> inicuts[i] >> min >> max;
    if(min==0 && max==0){
      for(int ibin=0;ibin<nbins;ibin++) varcuts[i].push_back(inicuts[i]);
      continue;
    }
    TH1F *hVar = new TH1F("hVar","vars",1000000,min,max);
    t->Draw(vars[i]+">>hVar",cuts+"&&(invm1S>9.7&&invm1S<9.95)||(invm1S>10.1&&invm1S<10.35)");
    float sum = 0;
    int ibin=0;
    float interval = (hVar->GetEntries()-hVar->GetBinContent(0)-hVar->GetBinContent(1000001))/(1.*nbins); 
    for (int j=hVar->GetNbinsX() ; j>0 && ibin<nbins ; j--){
      sum += hVar->GetBinContent(j);
      if (sum>=interval){
	sum -= interval;
	varcuts[i].push_back(hVar->GetBinLowEdge(j));
	ibin++;
      }
    }
    for (int j=1 ; j<=hVar->GetNbinsX() ; j++){
      if (hVar->GetBinContent(j)>0){
	cout << "& " << vars[i] << " >= " << hVar->GetBinLowEdge(j) << " &" << endl;
	if (varcuts[i].size()<nbins) varcuts[i].push_back(hVar->GetBinLowEdge(j));
	else varcuts[i][nbins-1] = hVar->GetBinLowEdge(j);
	j=hVar->GetNbinsX();
      }
    }
    if (varcuts[i].size()<nbins) {
      float last = varcuts[i][varcuts[i].size()-1];
      int missing = nbins-varcuts[i].size();
      for (int j=0 ; j<missing ; j++) varcuts[i].push_back(last);
    }
    delete hVar;
  }



  float ii[30];
  for (int i=0 ; i<nvars ; i++){
    ii[i]=nbins-1;
  }
  //cout << buildstr(nvars,nbins,vars,varcuts,ii) << endl;
  //pair<float,float> maxpair = SIG(cuts+buildstr(nvars,nbins,vars,varcuts,ii),t);

  for (int i=0 ; i<nvars ; i++){
    for(int j=0 ; j<nbins ; j++) if(varcuts[i][j]<inicuts[i]){
      ii[i] = j-1+(varcuts[i][j-1] - inicuts[i])/(varcuts[i][j-1] - varcuts[i][j]);
      j=nbins;
    }
  }  

  cout << buildstr(nvars,nbins,vars,varcuts,ii) << endl;

  pair<float,float> startpair = SIG(cuts+buildstr(nvars,nbins,vars,varcuts,ii),t);

  pair<float,float> bestpair;
  bestpair.first = startpair.first;
  bestpair.second = startpair.second;

  float best[30],jj[30];
  for (int i=0 ; i<nvars ; i++) {
    best[i] = ii[i];
    jj[i] = ii[i];
  }

  TString outname = "optimizeSig"+varset+".txt";

  ofstream out(outname);
  out << cuts+buildstr(nvars,nbins,vars,varcuts,ii) << endl;
  out << " Sig = " << bestpair.first << endl;

  TRandom3 *r = new TRandom3(0);

  int counter=maxcount;

  for (int j=0 ; j<100000 ; j++){

    for (int i=0 ; i<nvars ; i++) jj[i] = ii[i];

    if (counter>=maxcount){
      bool findich=true;
      
      for (int i=0 ; i<nvars ; i++) {
	ii[i] = best[i];
	jj[i] = best[i];
      }

      int ich = nvars*r->Uniform();
      while(findich) {
	if (varcuts[ich][0]!=varcuts[ich][nbins-1]) findich=false;
	else ich = nvars*r->Uniform();
	if(ich==nvars) ich=0; 
      }
      //for (int ich=0 ; ich<nvars ; ich++){
      jj[ich] = TMath::Max(0.,TMath::Min(1.*(nbins-1),ii[ich] + 2*r->Uniform()));
      
      findich=true;
      
      int ich2 = nvars*r->Uniform();
      while(findich) {
	if (varcuts[ich2][0]!=varcuts[ich2][nbins-1]&&ich2!=ich) findich=false;
	else ich2 = nvars*r->Uniform();
	if(ich2==nvars) ich2=0; 
      }
      //for (int ich2=0 ; ich2<nvars ; ich2++){
      jj[ich2] = TMath::Max(0.,TMath::Min(1.*(nbins-1),ii[ich2] - 2*r->Uniform()));
    }
    counter++;

    pair<float,float> nextpair=SIG(cuts+buildstr(nvars,nbins,vars,varcuts,jj),t);

    if(nextpair.first>bestpair.first&&nextpair.second<=bestpair.second&&nextpair.first>0) {
      counter=0;
      for (int i=0 ; i<nvars ; i++) {
	best[i] = jj[i];
	ii[i] = jj[i];
      }
      bestpair.first = nextpair.first;
      bestpair.second = nextpair.second;
      cout << "\t\t\t\tfound " << bestpair.first << " " << bestpair.second << endl;
      out << cuts + buildstr(nvars,nbins,vars,varcuts,jj) << endl;
      out << " Sig = " << bestpair.first << endl;
    } else if ((nextpair.second-bestpair.second)/bestpair.second > margin*r->Uniform()&&nextpair.first>=bestpair.first&&nextpair.first>0) {
      for (int i=0 ; i<nvars ; i++) ii[i] = jj[i];
      cout << "\t\t\t\tmargin " << nextpair.first << " " << nextpair.second << endl;
    }
  }
  
  return 0;
  
}

TString buildstr(int nvars, int nbins, TString *vars, vector<float> *varcuts, float *ii){
  
  TString result = "";

  for (int i=0 ; i<nvars ; i++){
    int iii = ii[i];
    float x = ii[i]-iii;
    stringstream kk;
    float kkflt;
    vector<float> vcuts = varcuts[i];
    if (iii>=nbins-1) kkflt = vcuts[nbins-1];
    else if (iii<0) kkflt = vcuts[0];
    else {
      kkflt = vcuts[iii]*(1-x) + vcuts[iii+1]*x;
    }
    kk << kkflt;
    TString kkstr = kk.str();
    result += " && " + vars[i] + " >= " + kkstr;
  }

  return result;

}

pair<float,float> SIG(TString cuts="",TTree *t,TString what="chib"){
  if (what=="chib"){
    TH1F *hMass = new TH1F("hMass","",30,9.5,11);
    t->Draw("invm1S>>hMass",cuts+"&&invm1S>9.5&&invm1S<11");
    
    pair<float,float> result;
    
    float nbkg = 6./7.*(hMass->GetBinContent(5)+hMass->GetBinContent(6)+hMass->GetBinContent(10)+hMass->GetBinContent(11))+(hMass->GetBinContent(12)+hMass->GetBinContent(13)+hMass->GetBinContent(17));
    
    float npk = hMass->GetBinContent(7)+hMass->GetBinContent(8)+hMass->GetBinContent(9)+hMass->GetBinContent(14)+hMass->GetBinContent(15)+hMass->GetBinContent(16);
    
    result.first = (npk - nbkg)/sqrt(npk);
    
    result.second = hMass->GetEntries() - npk - (hMass->GetBinContent(17)+hMass->GetBinContent(18)+hMass->GetBinContent(19)+hMass->GetBinContent(20)+hMass->GetBinContent(21)+hMass->GetBinContent(22));
    
    delete hMass;

    return result;

  }else{
    TH1F *hMass = new TH1F("hMass","",30,3.1,4);
    t->Draw("invm1S>>hMass",cuts+"&&invm1S>3.1&&invm1S<4");
    
    pair<float,float> result;
    
    float nbkg = 7/4.*(hMass->GetBinContent(9)+hMass->GetBinContent(17)+hMass->GetBinContent(18)+hMass->GetBinContent(19));
    
    float npk = hMass->GetBinContent(10)+hMass->GetBinContent(11)+hMass->GetBinContent(12)+hMass->GetBinContent(13)+hMass->GetBinContent(14)+hMass->GetBinContent(15)+hMass->GetBinContent(16);
    
    result.first = (npk - nbkg)/sqrt(npk);
    
    result.second = hMass->GetEntries() - npk;
    
    delete hMass;

    return result;

  }

};
