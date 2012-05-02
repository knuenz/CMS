#include <iostream>
#include <sstream>
#include <iomanip>


// RooFit Includes
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooArgusBG.h"
#include "RooAbsReal.h"
#include "RooChebychev.h"
#include "RooBinning.h"
#include "RooMinuit.h"
#include "RooAddition.h"

//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TFormula.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
//#include "ProfileLikelihoodCalculator.h"
//#include "HypoTestResult.h"

TH1F *hInvMnS;
TFile *fMix;

void BkgMixer(
		Char_t *MixerInputTree = "MixerInputTree",
		TTree* t=NULL,
		int nStateToMix=999,
		double rangeMin=0,
		double rangeMax=0,
		bool useExistingMixFile=true,
		int nMix=0,
		Char_t *dirstruct = "dirstruct",
		Char_t *ExcludeSignal = "ExcludeSignal",
		double cut_gammapt=1.,
		double gammaptK=1.,
		double gammaptD=1.,
		double DimuonMassPDG=1.,
		int nBinsCosAlphaTry=1000
		){

//	cout<<"root BkgMixer.C"<<endl;
//	return;

	hInvMnS = new TH1F("hInvMnS","",100,rangeMin-0.1,rangeMax+0.1);
        char bkgMixname[200];

        	if(!useExistingMixFile){
        		cout<<"Producing new BkgMixer file..."<<endl;

        	  int n_evt=nMix;

        	  char MixerInputTreeFileName[200];
        	  sprintf(MixerInputTreeFileName,"%s",MixerInputTree);
//        	  TFile *f = new TFile(MixerInputTreeFileName);
        	  sprintf(bkgMixname,"%s/eventMixer%dS.root",dirstruct,nStateToMix);
        	  fMix = new TFile(bkgMixname,"recreate");

        	  char gammaPtCutChar[200];
        	  sprintf(gammaPtCutChar,"gammapt>%f && gammapt>%f *Q+%f",cut_gammapt,gammaptK,gammaptD);
        	  TTree *t1 = (TTree*)t->CopyTree(ExcludeSignal);
        	  TTree *t2 = (TTree*)t1->CopyTree(gammaPtCutChar);




              const int nBinsCosAlphaDataHist=25;

              TH1F *hWNonEqu = new TH1F("hWNonEqu","; weights for cos#alpha",nBinsCosAlphaDataHist,-1,1);
              TH1F *hCosAlphaNonEqu = new TH1F("hCosAlphaNonEqu",";cos#alpha",nBinsCosAlphaDataHist,-1,1);
              TH1F *hCosAlphaDataNonEqu = new TH1F("hCosAlphaDataNonEqu",";cos#alpha",nBinsCosAlphaDataHist,-1,1);

              cout<<"Begin weighting of cosAlpha histo"<<endl;

              TH1F *hCosAlphaDataBuffer;

              double BinBordersCosAlpha[nBinsCosAlphaDataHist+1];
              int n_t2 = t2->GetEntries(ExcludeSignal);
              int targetNumEventsPerBin=n_t2/nBinsCosAlphaDataHist;

              BinBordersCosAlpha[0]=-1.;
              int iRealBin=1;
              for(int iBins=1;iBins<nBinsCosAlphaTry+1;iBins++){
              hCosAlphaDataBuffer = new TH1F("hCosAlphaDataBuffer",";cos#alpha",50,-1,1);
              double cosAlphaBufferCut=-1+2.*double(iBins)/double(nBinsCosAlphaTry);
              char AlphaCut[200];
              sprintf(AlphaCut,"(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz) < %f",cosAlphaBufferCut);
              t2->Draw("(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz)>>hCosAlphaDataBuffer",AlphaCut);
//              cout<<"hCosAlphaDataBuffer->GetEntries() "<<hCosAlphaDataBuffer->GetEntries()<<endl;

              if(hCosAlphaDataBuffer->GetEntries()>iRealBin*targetNumEventsPerBin){
//              	cout<<"hCosAlphaDataBuffer->GetEntries() = "<<hCosAlphaDataBuffer->GetEntries()<<endl;
              	cout<<"iCosAlphaBin = "<<iRealBin<<" / "<<nBinsCosAlphaDataHist<<endl;
              	BinBordersCosAlpha[iRealBin]=cosAlphaBufferCut-1./double(nBinsCosAlphaTry);
//              	cout<<"BinBordersCosAlpha[iRealBin] = "<<BinBordersCosAlpha[iRealBin]<<endl;
              	iRealBin++;
              }
              if(iRealBin>nBinsCosAlphaDataHist-1) break;
              delete hCosAlphaDataBuffer;
              }
              BinBordersCosAlpha[nBinsCosAlphaDataHist]=1.;

              double xBinsCosAlpha[nBinsCosAlphaDataHist+1];

              for(int iBinCosAlpha=0;iBinCosAlpha<nBinsCosAlphaDataHist+1;iBinCosAlpha++){
              	xBinsCosAlpha[iBinCosAlpha]=BinBordersCosAlpha[iBinCosAlpha];
              	cout<<iBinCosAlpha<<" xBinsCosAlpha[iBinCosAlpha]"<<xBinsCosAlpha[iBinCosAlpha]<<endl;
              }


              hCosAlphaDataNonEqu->GetXaxis()->Set(nBinsCosAlphaDataHist, xBinsCosAlpha);
              hCosAlphaNonEqu->GetXaxis()->Set(nBinsCosAlphaDataHist, xBinsCosAlpha);
              hWNonEqu->GetXaxis()->Set(nBinsCosAlphaDataHist, xBinsCosAlpha);

              t2->Draw("(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz)>>hCosAlphaDataNonEqu");




        	  double jpsimass,jpsipx,jpsipy,jpsipz,gammapt_,gammapx,gammapy,gammapz;

        	  t->SetBranchAddress("jpsimass",&jpsimass);
        	  t->SetBranchAddress("jpsipx",&jpsipx);
        	  t->SetBranchAddress("jpsipy",&jpsipy);
        	  t->SetBranchAddress("jpsipz",&jpsipz);
        	  t->SetBranchAddress("gammapt",&gammapt_);
        	  t->SetBranchAddress("gammapx",&gammapx);
        	  t->SetBranchAddress("gammapy",&gammapy);
        	  t->SetBranchAddress("gammapz",&gammapz);

//              RooRealVar invm1S_write = RooRealVar("invm1S_write", "invm1S_write",0,100);
//              RooRealVar gammapt_write = RooRealVar("gammapt_write", "gammapt_write",0,100);
//              RooRealVar cosalpha_write = RooRealVar("cosalpha_write", "cosalpha_write",-1,1);
//              RooRealVar Q_write = RooRealVar("Q_write", "Q_write",0,100);
//              RooArgSet argSet = RooArgSet();
//              argSet.add(invm1S_write);
//              argSet.add(gammapt_write);
//              argSet.add(cosalpha_write);
//              argSet.add(Q_write);
//
//              RooDataSet rds = RooDataSet("d","d",argSet);


        	  int i_evt=0;
        	  int bin;

        	  TRandom3 *r = new TRandom3();

        	  int n = t->GetEntries();

        	  for(int iIter=1;iIter<3;iIter++){

        		  if(iIter==1) cout<<"Start mixing "<<nStateToMix<<"S + gamma..."<<endl;

        		  if(iIter==2){
            		for(int iDivide=1;iDivide<nBinsCosAlphaDataHist+1;iDivide++){
//            			cout<<"dataHist "<<hCosAlphaDataNonEqu->GetBinContent(iDivide)<<endl;
//            			cout<<" mixHist "<<hCosAlphaNonEqu->GetBinContent(iDivide)<<endl;
            		hWNonEqu->SetBinContent(iDivide,hCosAlphaDataNonEqu->GetBinContent(iDivide)/hCosAlphaNonEqu->GetBinContent(iDivide));
            		}
            		cout<<"Finish weighting of cosAlpha histo"<<endl;
        		  }

        		  i_evt=0;

        	  while(i_evt<n_evt){
        	    double m_jpsimass,m_jpsipx,m_jpsipy,m_jpsipz;
        	    int i=n*r->Uniform();
        	    if (i==n) i=0;
        	    t->GetEntry(i);
        	    m_jpsimass = jpsimass;
        	    m_jpsipx = jpsipx;
        	    m_jpsipy = jpsipy;
        	    m_jpsipz = jpsipz;
        	    double m_jpsip = sqrt(m_jpsipx*m_jpsipx + m_jpsipy*m_jpsipy + m_jpsipz*m_jpsipz);
        	    double m_jpsie = sqrt(m_jpsimass*m_jpsimass + m_jpsip*m_jpsip);
        	    int j=n*r->Uniform();
        	    if (j==n) j=0;
        	    if (j!=i){
        	      t->GetEntry(j);
        	      double gammae = sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz);
        	      double chibe = gammae + m_jpsie;
        	      double chibpx = gammapx + m_jpsipx;
        	      double chibpy = gammapy + m_jpsipy;
        	      double chibpz = gammapz + m_jpsipz;
        	      double chibmass = sqrt(chibe*chibe - chibpx*chibpx - chibpy*chibpy - chibpz*chibpz);
        	      double Q = chibmass - m_jpsimass;
        	      if (gammapt_>cut_gammapt&&gammapt_>gammaptK*Q+gammaptD){
        		double cosalpha = (m_jpsipx*gammapx + m_jpsipy*gammapy + m_jpsipz*gammapz)/gammae/m_jpsip;
        		double chibmass_corr = Q + DimuonMassPDG;
        		if(iIter==1) hCosAlphaNonEqu->Fill(cosalpha);
        		if(iIter==2) {
         	    	bin = hWNonEqu->FindBin(cosalpha);
              	    hInvMnS->Fill(chibmass_corr,hWNonEqu->GetBinContent(bin));
//              	    cout<<"chibmass_corr "<<chibmass_corr<<endl;
//            	    cout<<"bin "<<bin<<endl;
//            	    cout<<"hWNonEqu->GetBinContent(bin) "<<hWNonEqu->GetBinContent(bin)<<endl;
        		}

       		if (i_evt%10000==9999) cout << i_evt+1 << endl;
        		i_evt ++;
        	      }
        	    }
        	  }
        	  }

          	  hInvMnS->Write();
          	  hCosAlphaNonEqu->Write();
          	  hCosAlphaDataNonEqu->Write();
          	  hWNonEqu->Write();


        	}

        	else {
        		cout<<"Using existing BkgMixer file..."<<endl;
          	    sprintf(bkgMixname,"%s/eventMixer%dS.root",dirstruct,nStateToMix);
        		fMix = new TFile(bkgMixname,"update");
              	hInvMnS=(TH1F*)fMix->Get("hInvMnS");
        	}

        	hInvMnS->Print();


//        	hInvMnS->Rebin(2);



}
