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

TH1F *hInvMnS;
TFile *fMix;

void BkgMixer(
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
		int nBinsCosAlphaTry=1000,
		double GammaEtaBorder=1000.,
		double cut_gammapt_Midrap=1000.,
		double cut_gammapt_Forward=1000.,
		double cut_Ypt=1000.,
		double cut_gammaeta=1000.
		){


	hInvMnS = new TH1F("hInvMnS","",100,rangeMin-0.1,rangeMax+0.1);
        char bkgMixname[200];

        	if(!useExistingMixFile){
        		cout<<"Producing new BkgMixer file..."<<endl;

        	  int n_evt=nMix;

        	  sprintf(bkgMixname,"%s/eventMixer%dS.root",dirstruct,nStateToMix);
        	  fMix = new TFile(bkgMixname,"recreate");


        	  char gammaPtCutChar[200];
//        	  sprintf(gammaPtCutChar,"gammapt>%f && gammapt>%f *Q+%f",cut_gammapt,gammaptK,gammaptD);
//        	  sprintf(gammaPtCutChar,"jpsipt>%f&&TMath::Abs(gammaeta)<%f&&gammapt>%f && gammapt>%f *Q+%f",cut_Ypt,cut_gammaeta,cut_gammapt,gammaptK,gammaptD);
        	  sprintf(gammaPtCutChar,"");
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

              TH1F *hCosAlphaDataBuffer2 = new TH1F("hCosAlphaDataBuffer2",";cos#alpha",nBinsCosAlphaTry,-1,1);
              t2->Draw("(jpsipx*gammapx + jpsipy*gammapy + jpsipz*gammapz)/sqrt(gammapx*gammapx + gammapy*gammapy + gammapz*gammapz)/sqrt(jpsipx*jpsipx + jpsipy*jpsipy + jpsipz*jpsipz)>>hCosAlphaDataBuffer2");


              Double_t xq[nBinsCosAlphaDataHist-1];  // position where to compute the quantiles in [0,1]
              Double_t yq[nBinsCosAlphaDataHist-1];  // array to contain the quantiles
              for (Int_t i=0;i<nBinsCosAlphaDataHist-1;i++) xq[i] = Float_t(i+1)/nBinsCosAlphaDataHist;
              hCosAlphaDataBuffer2->GetQuantiles(nBinsCosAlphaDataHist-1,yq,xq);

              BinBordersCosAlpha[0]=-1.;
              for (Int_t i=0;i<nBinsCosAlphaDataHist-1;i++){
            	  BinBordersCosAlpha[i+1]=yq[i];
                	cout<<"iCosAlphaBin = "<<i+1<<" / "<<nBinsCosAlphaDataHist<<" --> "<<BinBordersCosAlpha[i+1]<<endl;
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




        	  double jpsimass,jpsipx,jpsipy,jpsipz,gammapt_,gammapx,gammapy,gammapz,jpsipt_,gammaeta_;

        	  t->SetBranchAddress("jpsimass",&jpsimass);
        	  t->SetBranchAddress("jpsipx",&jpsipx);
        	  t->SetBranchAddress("jpsipy",&jpsipy);
        	  t->SetBranchAddress("jpsipz",&jpsipz);
        	  t->SetBranchAddress("jpsipt",&jpsipt_);
        	  t->SetBranchAddress("gammapt",&gammapt_);
        	  t->SetBranchAddress("gammapx",&gammapx);
        	  t->SetBranchAddress("gammapy",&gammapy);
        	  t->SetBranchAddress("gammapz",&gammapz);
        	  t->SetBranchAddress("gammaeta",&gammaeta_);


    	      cout<<"Not Cutting in BkgMixer.C"<<endl;

        	  int i_evt=0;
        	  int bin;

        	  TRandom3 *r = new TRandom3();

        	  int n = t->GetEntries();

        	  for(int iIter=1;iIter<3;iIter++){

        		  if(iIter==1) cout<<"Start mixing "<<nStateToMix<<"S + gamma..."<<endl;

        		  if(iIter==2){
            		for(int iDivide=1;iDivide<nBinsCosAlphaDataHist+1;iDivide++){
                		hWNonEqu->SetBinContent(iDivide,hCosAlphaDataNonEqu->GetBinContent(iDivide)/hCosAlphaNonEqu->GetBinContent(iDivide));
                		if(nStateToMix==2) hWNonEqu->SetBinContent(iDivide, 1);
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
        	      double Absgammaeta_=TMath::Abs(gammaeta_);//TMath::Abs(0.5*TMath::Log((TMath::Sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)+gammapz)/(TMath::Sqrt(gammapx*gammapx+gammapy*gammapy+gammapz*gammapz)-gammapz)));

//        	      if (gammapt_>cut_gammapt&&gammapt_>gammaptK*Q+gammaptD){

//        	    	  if (jpsipt_>cut_Ypt&&Absgammaeta_<cut_gammaeta&&gammapt_>cut_gammapt&&gammapt_>gammaptK*Q+gammaptD){


//        	    		  if(jpsipt_>-100){


        	     double cosalpha = (m_jpsipx*gammapx + m_jpsipy*gammapy + m_jpsipz*gammapz)/gammae/m_jpsip;
        		double chibmass_corr = Q + DimuonMassPDG;
        		if(iIter==1) hCosAlphaNonEqu->Fill(cosalpha);
        		if(iIter==2) {
         	    	bin = hWNonEqu->FindBin(cosalpha);
              	    hInvMnS->Fill(chibmass_corr,hWNonEqu->GetBinContent(bin));
        		}

       		if (i_evt%10000==9999) cout << i_evt+1 << endl;
        		i_evt ++;
//        	      }
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





}
