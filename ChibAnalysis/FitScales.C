/*
 * SmearEff.C
 *
 *  Created on: Dec 4, 2011
 *      Author: valentinknuenz
 */

#include "Riostream.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TAttLine.h"
#include "TColor.h"


void FitScales(){

	   gStyle->SetPalette(1,0);
	   gStyle->SetPadBottomMargin(0.12);
	   gStyle->SetPadLeftMargin(0.13);
	   gStyle->SetPadRightMargin(0.15);

	   gStyle->SetTickLength(-0.02, "xyz");
	   gStyle->SetLabelOffset(0.02, "x");
	   gStyle->SetLabelOffset(0.02, "y");
	   gStyle->SetTitleOffset(1.3, "x");
	   gStyle->SetTitleOffset(1.4, "y");
	   gStyle->SetTitleFillColor(kWhite);


	   double M_chic1=3.51066;
	   double M_chic2=3.55620;
	   double M_chib1P1=9.89278;
	   double M_chib1P2=9.91221;
	   double M_chib2P1=10.25546;
	   double M_chib2P2=10.26865;

	   double ratio2over1=0.5;
	   double fracJ1=1/(1+ratio2over1);

	   double M_Jpsi=3.096916;
	   double M_Y1S=9.4603;
	   double M_Y2S=10.02326;


	   double Q[4];


	  char graphName[200];

	  sprintf(graphName,"PESgraph");

	  double PES_Q[4]={0.98494,0.9841,0.9856,0.9835};
	  double err_PES_Q[4]={0.00058,0.0033,0.0014,0.003};

	  double PWS_Q[4]={1.629,1.70,1.74,1.71};
	  double err_PWS_Q[4]={0.058,0.27,0.11,0.32};

	   Q[0]=(M_chic1-M_Jpsi)*PES_Q[0];
	   Q[1]=(M_chib1P1*fracJ1+M_chib1P2*(1-fracJ1)-M_Y1S)*PES_Q[1];
	   Q[2]=(M_chic2-M_Jpsi)*PES_Q[2];
	   Q[3]=(M_chib2P1*fracJ1+M_chib2P2*(1-fracJ1)-M_Y1S)*PES_Q[3];

	   cout<<Q[0]<<endl;
	   cout<<Q[1]<<endl;
	   cout<<Q[2]<<endl;
	   cout<<Q[3]<<endl;
	   cout<<PES_Q[0]<<endl;
	   cout<<PES_Q[1]<<endl;
	   cout<<PES_Q[2]<<endl;
	   cout<<PES_Q[3]<<endl;

	   cout<<"$Slashchi_{c1}$&"<<Q[0]<<"&"<<PES_Q[0]<<"$Slashpm$ "<<err_PES_Q[0]<<"SlashSlashSlashhline"<<endl;
	   cout<<"$Slashchi_{b}(1P)$&"<<Q[1]<<"&"<<PES_Q[1]<<"$Slashpm$ "<<err_PES_Q[1]<<"SlashSlashSlashhline"<<endl;
	   cout<<"$Slashchi_{c2}$&"<<Q[2]<<"&"<<PES_Q[2]<<"$Slashpm$ "<<err_PES_Q[2]<<"SlashSlashSlashhline"<<endl;
	   cout<<"$Slashchi_{b}(2P)$&"<<Q[3]<<"&"<<PES_Q[3]<<"$Slashpm$ "<<err_PES_Q[3]<<"SlashSlashSlashhline"<<endl;

	   cout<<"$Slashchi_{c1}$&"<<Q[0]<<"&"<<PWS_Q[0]<<"$Slashpm$ "<<err_PWS_Q[0]<<"SlashSlashSlashhline"<<endl;
	   cout<<"$Slashchi_{b}(1P)$&"<<Q[1]<<"&"<<PWS_Q[1]<<"$Slashpm$ "<<err_PWS_Q[1]<<"SlashSlashSlashhline"<<endl;
	   cout<<"$Slashchi_{c2}$&"<<Q[2]<<"&"<<PWS_Q[2]<<"$Slashpm$ "<<err_PWS_Q[2]<<"SlashSlashSlashhline"<<endl;
	   cout<<"$Slashchi_{b}(2P)$&"<<Q[3]<<"&"<<PWS_Q[3]<<"$Slashpm$ "<<err_PWS_Q[3]<<"SlashSlashSlashhline"<<endl;


	   char text[200];


	      double plotQmin=0;
	      double plotQmax=1.1;
	      double plotPESmin=0.975;
	      double plotPESmax=0.995;
	      double plotPWSmin=1;
	      double plotPWSmax=2.5;

	      double FontSize=0.0215;
	      double xText=plotQmin+(plotQmax-plotQmin)*0.025;

	      char FitOptions[200];
	      sprintf(FitOptions,"EFNR");

	    TCanvas* ScaleCanvas = new TCanvas("ScaleCanvas","ScaleCanvas",1000, 800);
	    ScaleCanvas->SetFillColor(kWhite);
	    ScaleCanvas->cd(1);
	    ScaleCanvas->SetLeftMargin(0.2);
	    gPad->SetFillColor(kWhite);


		  TGraphErrors *PESgraph = new TGraphErrors(4,Q,PES_Q,0,err_PES_Q);
		  TGraphErrors *PWSgraph = new TGraphErrors(4,Q,PWS_Q,0,err_PWS_Q);


		  /////////////////////
		  bool DrawExt=true;


			TH1F *plotHisto = new TH1F;
			plotHisto = ScaleCanvas->DrawFrame(plotQmin,plotPESmin,plotQmax,plotPESmax);
			plotHisto->SetXTitle("Q [GeV]");
			plotHisto->SetYTitle("Photon energy scale");
			plotHisto->GetYaxis()->SetTitleOffset(2);


		  PESgraph->SetMarkerColor(kGreen-2);
		  PESgraph->SetMarkerStyle(21);
		  PESgraph->SetTitle(0);
		  PESgraph->Draw("P");
		  PESgraph->SaveAs("ScalePlots/PESgraph.root");


	      TF1* f1PES = new TF1("f1PES","pol1",plotQmin,plotQmax);
	      PESgraph->Fit("f1PES",FitOptions);
	      f1PES->SetLineWidth(0.4);
	      f1PES->SetLineColor(kRed);
	      f1PES->SetLineStyle(1);
//	      f1PES->Draw("same");

	      double PES_p0 = f1PES->GetParameter(0);
	      double err_PES_p0 = f1PES->GetParError(0);
	      double PES_p1 = f1PES->GetParameter(1);
	      double err_PES_p1 = f1PES->GetParError(1);
	      double PES_chi2=f1PES->GetChisquare();
	      double PES_NDF=f1PES->GetNDF();
	      double PES_BIC=PES_chi2+2*TMath::Log(4);

	      double highest=plotPESmax-(plotPESmax-plotPESmin)*0.05;

	      sprintf(text,"#color[2]{Fitting pol1:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", PES_p0, err_PES_p0, PES_p1, err_PES_p1,PES_chi2,PES_NDF,PES_BIC);
	      TLatex textPES1 = TLatex(xText,highest,text);
	      textPES1.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
	      textPES1.Draw( "same" )                                                                                                                                                                                                                                                 ;

	      TF1* f1PESconst = new TF1("f1PESconst","pol0",plotQmin,plotQmax);
	      PESgraph->Fit("f1PESconst",FitOptions);
	      f1PESconst->SetLineWidth(0.4);
	      f1PESconst->SetLineStyle(1);
	      f1PESconst->SetLineColor(kBlue);
	      f1PESconst->Draw("same");

	      double PES_const_p0 = f1PESconst->GetParameter(0);
	      double err_PES_const_p0 = f1PESconst->GetParError(0);
	      double PES_const_chi2=f1PESconst->GetChisquare();
	      double PES_const_NDF=f1PESconst->GetNDF();
	      double PES_const_BIC=PES_const_chi2+1*TMath::Log(4);

	      TF1* f1PES_m = new TF1("f1PES_m","pol0",plotQmin,plotQmax);
	      f1PES_m->SetParameter(0,PES_const_p0+err_PES_const_p0);
	      f1PES_m->SetLineWidth(0.2);
	      f1PES_m->SetLineColor(kBlue);
	      f1PES_m->SetLineStyle(2);
	      if(DrawExt) f1PES_m->Draw("same");
	      TF1* f1PES_p = new TF1("f1PES_p","pol0",plotQmin,plotQmax);
	      f1PES_p->SetParameter(0,PES_const_p0-err_PES_const_p0);
	      f1PES_p->SetLineWidth(0.2);
	      f1PES_p->SetLineColor(kBlue);
	      f1PES_p->SetLineStyle(2);
	      if(DrawExt) f1PES_p->Draw("same");


	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.1;

	      sprintf(text,"#color[4]{Fitting constant:} PES = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", PES_const_p0, err_PES_const_p0, PES_const_chi2, PES_const_NDF,PES_const_BIC);
	      TLatex textPES2 = TLatex(xText,highest,text);
	      textPES2.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
	      textPES2.Draw( "same" )                                                                                                                                                                                                                                                 ;



	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.95; TLatex textPESchi1 = TLatex(Q[0]-0.05,highest,"#chi_{c1}"); textPESchi1.SetTextSize(FontSize); textPESchi1.Draw( "same" );
	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.9; TLatex textPESchi2 = TLatex(Q[1]-0.025,highest,"#chi_{b}(1P)"); textPESchi2.SetTextSize(FontSize); textPESchi2.Draw( "same" );
	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.95; TLatex textPESchi3 = TLatex(Q[2]+0.025,highest,"#chi_{c2}"); textPESchi3.SetTextSize(FontSize); textPESchi3.Draw( "same" );
	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.9; TLatex textPESchi4 = TLatex(Q[3]-0.025,highest,"#chi_{b}(2P)"); textPESchi4.SetTextSize(FontSize); textPESchi4.Draw( "same" );

		  PESgraph->Draw("P");
	      ScaleCanvas->SaveAs("ScalePlots/PESfit.pdf");





			plotHisto = ScaleCanvas->DrawFrame(plotQmin,plotPWSmin,plotQmax,plotPWSmax);
			plotHisto->SetXTitle("Q [GeV]");
			plotHisto->SetYTitle("Q resolution #sigma_{Q}/Q [%]");
			plotHisto->GetYaxis()->SetTitleOffset(2);


		  PWSgraph->SetMarkerColor(kGreen-2);
		  PWSgraph->SetMarkerStyle(21);
		  PWSgraph->SetTitle(0);
		  PWSgraph->Draw("P");
		  PWSgraph->SaveAs("ScalePlots/PWSgraph.root");


	      TF1* f1PWS = new TF1("f1PWS","pol1",plotQmin,plotQmax);
	      PWSgraph->Fit("f1PWS",FitOptions);
	      f1PWS->SetLineWidth(0.4);
	      f1PWS->SetLineStyle(1);
	      f1PWS->SetLineColor(kRed);
//	      f1PWS->Draw("same");

	      double PWS_p0 = f1PWS->GetParameter(0);
	      double err_PWS_p0 = f1PWS->GetParError(0);
	      double PWS_p1 = f1PWS->GetParameter(1);
	      double err_PWS_p1 = f1PWS->GetParError(1);
	      double PWS_chi2=f1PWS->GetChisquare();
	      double PWS_NDF=f1PWS->GetNDF();
	      double PWS_BIC=PWS_chi2+2*TMath::Log(4);

	      highest=plotPWSmax-(plotPWSmax-plotPWSmin)*0.05;

	      sprintf(text,"#color[2]{Fitting pol1:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", PWS_p0, err_PWS_p0, PWS_p1, err_PWS_p1,PWS_chi2,PWS_NDF,PWS_BIC);
	      TLatex textPWS1 = TLatex(xText,highest,text);
	      textPWS1.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
	      textPWS1.Draw( "same" )                                                                                                                                                                                                                                                 ;

	      TF1* f1PWSconst = new TF1("f1PWSconst","pol0",plotQmin,plotQmax);
	      PWSgraph->Fit("f1PWSconst",FitOptions);
	      f1PWSconst->SetLineWidth(0.4);
	      f1PWSconst->SetLineStyle(1);
	      f1PWSconst->SetLineColor(kBlue);
	      f1PWSconst->Draw("same");

	      double PWS_const_p0 = f1PWSconst->GetParameter(0);
	      double err_PWS_const_p0 = f1PWSconst->GetParError(0);
	      double PWS_const_chi2=f1PWSconst->GetChisquare();
	      double PWS_const_NDF=f1PWSconst->GetNDF();
	      double PWS_const_BIC=PWS_const_chi2+1*TMath::Log(4);

	      TF1* f1PWS_m = new TF1("f1PWS_m","pol0",plotQmin,plotQmax);
	      f1PWS_m->SetParameter(0,PWS_const_p0+err_PWS_const_p0);
	      f1PWS_m->SetLineWidth(0.2);
	      f1PWS_m->SetLineColor(kBlue);
	      f1PWS_m->SetLineStyle(2);
	      if(DrawExt) f1PWS_m->Draw("same");
	      TF1* f1PWS_p = new TF1("f1PWS_p","pol0",plotQmin,plotQmax);
	      f1PWS_p->SetParameter(0,PWS_const_p0-err_PWS_const_p0);
	      f1PWS_p->SetLineWidth(0.2);
	      f1PWS_p->SetLineColor(kBlue);
	      f1PWS_p->SetLineStyle(2);
	      if(DrawExt) f1PWS_p->Draw("same");

	      highest=plotPWSmax-(plotPWSmax-plotPWSmin)*0.1;

	      sprintf(text,"#color[4]{Fitting constant:} #sigma_{Q}/Q = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", PWS_const_p0, err_PWS_const_p0, PWS_const_chi2, PWS_const_NDF,PWS_const_BIC);
	      TLatex textPWS2 = TLatex(xText,highest,text);
	      textPWS2.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
	      textPWS2.Draw( "same" )                                                                                                                                                                                                                                                 ;

	      highest=plotPWSmax-(plotPWSmax-plotPWSmin)*0.95; TLatex textPWSchi1 = TLatex(Q[0]-0.05,highest,"#chi_{c1}"); textPWSchi1.SetTextSize(FontSize); textPWSchi1.Draw( "same" );
	      highest=plotPWSmax-(plotPWSmax-plotPWSmin)*0.9; TLatex textPWSchi2 = TLatex(Q[1]-0.025,highest,"#chi_{b}(1P)"); textPWSchi2.SetTextSize(FontSize); textPWSchi2.Draw( "same" );
	      highest=plotPWSmax-(plotPWSmax-plotPWSmin)*0.95; TLatex textPWSchi3 = TLatex(Q[2]+0.025,highest,"#chi_{c2}"); textPWSchi3.SetTextSize(FontSize); textPWSchi3.Draw( "same" );
	      highest=plotPWSmax-(plotPWSmax-plotPWSmin)*0.9; TLatex textPWSchi4 = TLatex(Q[3]-0.025,highest,"#chi_{b}(2P)"); textPWSchi4.SetTextSize(FontSize); textPWSchi4.Draw( "same" );

		  PWSgraph->Draw("P");
	      ScaleCanvas->SaveAs("ScalePlots/PWSfit.pdf");


}

Double_t funcPES(Double_t *x, Double_t *par){

  Double_t result = par[0]*x+par[1];

  return result;
}

