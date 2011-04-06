/*
 * GenScript.cc
 *
 *  Created on: 20.03.2011
 *      Author: valentinknuenz
 */

#include <iostream>
#include <sstream>
#include <cstring>

#include "commonVar.h"
#include "CompositeModelBuilder.h"

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolarizationPdf.h"
#include "RooGenericPdf.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TLegend.h"
#include "TList.h"


void Plot(RooAbsData &plotset, RooRealVar &plotvar, RooAbsPdf &plotpdf, bool prompt){
	using namespace RooFit;
	using namespace std;


 	RooPlot* plot = new RooPlot;
 	plot = plotvar.frame(50) ;
 	plotset.plotOn(plot,DataError(RooAbsData::SumW2),MarkerSize(0.4));
 	plotpdf.plotOn(plot,LineWidth(2),Normalization(1.0));

    TCanvas* plotcanvas = new TCanvas("plotcanvas","plotcanvas",1600,1000);
    plotcanvas->SetFillColor(kWhite);
    plotcanvas->cd(1) ; gPad->SetFillColor(kWhite); plot->Draw();

    const char* PlotFilenamePart=plotvar.GetName();
    const char* PlotFilenameTitle=plotset.GetTitle();
    const char* PlotFilenameName=plotset.GetName();

	char contributionname[200];
    if(prompt) sprintf(contributionname,"PRregion");
    if(!prompt) sprintf(contributionname,"NPregion");

    char dirstruct[200];
    sprintf(dirstruct,"Plots/ClosureTest/%s",PlotFilenameName);
    gSystem->mkdir(dirstruct);

	char PlotFilename[200];
    sprintf(PlotFilename,"%s/%s_%s_%s.png",dirstruct,PlotFilenamePart,contributionname,PlotFilenameTitle);

    plotcanvas->SaveAs(PlotFilename);
    plotcanvas->Close();

    delete plotcanvas;
    delete plot;
}



int main(int argc, char** argv) {
	using namespace std;
	using namespace RooFit;



cout<<"Hello Project Closure"<<endl;

//////////////////// ROOREALVARS ///////////////////////

RooRealVar JpsiMass("JpsiMass","M [GeV]",2.7,3.5);
RooRealVar JpsiRap("JpsiRap","#nu",-2.4,2.4);
RooRealVar Jpsict("Jpsict","l_{J/#psi} [mm]",-1,2.5);
RooRealVar JpsictErr("JpsictErr","Error on l_{J/#psi} [mm]",1e-3,(1-1e-6));
RooRealVar JpsiPt("JpsiPt","pT [GeV]",0,40);
RooRealVar costh_CS("costh_CS","cos #theta_{CS}",-1,1);
RooRealVar phi_CS("phi_CS","#phi_{CS} [deg]",0,90);
RooRealVar costh_HX("costh_HX","cos#theta_{HX}",-1,1);
RooRealVar phi_HX("phi_HX","#phi_{HX} [deg]",0,90);

RooRealVar HLT_DoubleMu0("HLT_DoubleMu0","HLT_DoubleMu0",0.5,1.5);
RooRealVar JpsiVprob("JpsiVprob","JpsiVprob",0.01,100000);

RooRealVar MCType("MCType","MCType",0,2);//0=PR,1=NP,2=BK
RooRealVar MCweight("MCweight","MCweight",1);
RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
varlist.add(JpsictErr);
varlist.add(JpsiVprob);
varlist.add(HLT_DoubleMu0);
varlist.add(MCType);
varlist.add(MCweight);

//////////////////// INPUT FILES ///////////////////////

//gSystem->Rename("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempPRNP.root","/scratch/knuenz/Polarization/RootInput/ProjectClosure/temptempPRNP.root");
//gSystem->Rename("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempGEN.root","/scratch/knuenz/Polarization/RootInput/ProjectClosure/temptempGEN.root");

//TFile *PRfile = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_prep_PR_.root","UPDATE");
//TFile *NPfile = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_prep_NP_.root","UPDATE");


TFile *PRfile = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_HXth_plus1.root","UPDATE");
TFile *NPfile = new TFile("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_HXth_plus1.root","UPDATE");

int cutVal=70;

TFile *PRgeomAccfile = new TFile("/scratch/knuenz/Polarization/RootInput/geomAccHistos_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root","UPDATE");
TFile *PRrecoEfffile = new TFile("/scratch/knuenz/Polarization/RootInput/recoEffHistos_ATLASPT_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
TFile *PRtrigEfffile = new 
TFile("/scratch/knuenz/Polarization/RootInput/AccEffCutUnRe70_trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root","UPDATE");

TFile *NPgeomAccfile = new TFile("/scratch/knuenz/Polarization/RootInput/geomAccHistos_NP_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root","UPDATE");
TFile *NPrecoEfffile = new TFile("/scratch/knuenz/Polarization/RootInput/recoEffHistos_NP_ATLASPT_19March2011_phiFolded_zeroBinsCorrected.root","UPDATE");
TFile *NPtrigEfffile = new 
TFile("/scratch/knuenz/Polarization/RootInput/AccEffCutUnRe70_trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root","UPDATE");

/*
TFile *PRgeomAccfile = new TFile("2_7_PRgeomAccATLASpT24Mar.root","UPDATE");
TFile *PRrecoEfffile = new TFile("2_7_PRrecoEffATLASpT24Mar.root","UPDATE");
TFile *PRtrigEfffile = new TFile("2_7_PRtrigEffATLASpT24Mar.root","UPDATE");

TFile *NPgeomAccfile = new TFile("2_7_NPgeomAccATLASpT24Mar.root","UPDATE");
TFile *NPrecoEfffile = new TFile("2_7_NPrecoEffATLASpT24Mar.root","UPDATE");
TFile *NPtrigEfffile = new TFile("2_7_NPtrigEffATLASpT24Mar.root","UPDATE");
*/


/*
TFile *PRgeomAccfile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
TFile *PRrecoEfffile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
TFile *PRtrigEfffile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");

TFile *NPgeomAccfile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
TFile *NPrecoEfffile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
TFile *NPtrigEfffile = new TFile("/scratch/knuenz/Polarization/RootInput/flatHisto.root","UPDATE");
*/

TFile *fileforwhateverrootorcplusplusneedsme = new TFile("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/hollodrio.root","RECREATE");

////////////////////////////////////////////////////////
//////////////////// ADJUSTMENTS ///////////////////////
////////////////////////////////////////////////////////

bool scenIn(false);
double pTin;
double yin;
double scen;
double iter;
bool plot(false);

for( int i=0;i < argc; ++i ) {
  if(std::string(argv[i]).find("pT") != std::string::npos) {char* pTchar = argv[i]; char* pTchar2 = strtok (pTchar, "p"); pTin = atof(pTchar2); cout<<"pTin = "<<pTin<<endl;}
  if(std::string(argv[i]).find("rap") != std::string::npos) {char* ychar = argv[i]; char* ychar2 = strtok (ychar, "r"); yin = atof(ychar2); cout<<"yin = "<<yin<<endl;}
  if(std::string(argv[i]).find("scen") != std::string::npos) {char* scenchar = argv[i]; char* scenchar2 = strtok (scenchar, "s"); scen = atof(scenchar2); cout<<"scen = "<<scen<<endl;}
  if(std::string(argv[i]).find("iter") != std::string::npos) {char* iterchar = argv[i]; char* iterchar2 = strtok (iterchar, "i"); iter = atof(iterchar2); cout<<"iter = "<<iter<<endl;}
  if(std::string(argv[i]).find("--plot") != std::string::npos) plot=true;

}

char JOBNAME[200];
sprintf(JOBNAME,"AccEffCutUnRe%d",cutVal);
char GenID[200];
sprintf(GenID,"GEN%1.0f_%s_scen%1.0f",iter,JOBNAME,scen);
cout<<JOBNAME<<endl;
int scenario = scen;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


//////////////////// LOOP IT ///////////////////////

TTree* dataTreePR = (TTree*)PRfile->Get("data");
TTree* dataTreeNP = (TTree*)NPfile->Get("data");

RooDataSet *dataPR = new RooDataSet("dataPR","Prompt MC",varlist,Import(*dataTreePR),WeightVar(MCweight));
RooDataSet *dataNP = new RooDataSet("dataNP","Non Prompt MC",varlist,Import(*dataTreeNP),WeightVar(MCweight));

dataPR->Print();
dataNP->Print();

//dataPR->reduce("JpsiVprob > 0.01");
//dataNP->reduce("JpsiVprob > 0.01");

//dataPR->Print();
//dataNP->Print();

int generations=iter;

char filenamePRNP[200];
char filenameGEN[200];

//	  gSystem->Rename("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root","/scratch/knuenz/Polarization/RootInput/ProjectClosure/temptempYPT.root");

for(int yBin =yin-1; yBin < jpsi::kNbRapForPTBins+yin-5; ++yBin) {
	  for(int ptBin =pTin-1; ptBin < jpsi::kNbPTBins[yBin+1]+pTin-6-yin; ++ptBin) {

		  for(int generation = 1; generation < generations+1; ++generation) {



//////////////////// REDUCE DATASETS AND READ IN HISTOGRAMS ///////////////////////

	 std::stringstream binName,cutString;
	 binName << "pt" << ptBin+1 << "_rapidity" << yBin+1;
	 cutString << "JpsiPt > " << jpsi::pTRange[yBin+1][ptBin] << " && JpsiPt < " << jpsi::pTRange[yBin+1][ptBin+1] << " && "
	 << "abs(JpsiRap) > " << jpsi::rapForPTRange[yBin] << " && abs(JpsiRap) < " << jpsi::rapForPTRange[yBin+1];
	 std::cout << cutString.str() << std::endl;

	 char reduce[200];
	 sprintf(reduce,"JpsiPt > %f && JpsiPt < %f && abs(JpsiRap) > %f && abs(JpsiRap) < %f",jpsi::pTRange[yBin+1][ptBin],jpsi::pTRange[yBin+1][ptBin+1],jpsi::rapForPTRange[yBin],jpsi::rapForPTRange[yBin+1]);// && JpsictErr > 0.0004 && JpsictErr < 0.999999

	 RooDataSet* dataPRbin = (RooDataSet*)dataPR->reduce(reduce);
	 RooDataSet* dataNPbin = (RooDataSet*)dataNP->reduce(reduce);

	 char thisplotname[200];
	 sprintf(thisplotname,"massplot");

	 cout<<"Extract maps from data"<<endl;

     TH2F *AccMap_CS_PR,*AccMap_HX_PR;
     TH2F *RecoEffMap_CS_PR,*RecoEffMap_HX_PR;
     TH2F *TrigEffMap_CS_PR,*TrigEffMap_HX_PR;
     TH2F *AccMap_CS_NP,*AccMap_HX_NP;
     TH2F *RecoEffMap_CS_NP,*RecoEffMap_HX_NP;
     TH2F *TrigEffMap_CS_NP,*TrigEffMap_HX_NP;

     std::stringstream nameHXp,nameCSp;
     nameHXp << "hAcc2D_HX_pT" << ptBin+1 << "_rap" << yBin+1;
     nameCSp << "hAcc2D_CS_pT" << ptBin+1 << "_rap" << yBin+1;
     AccMap_CS_PR = (TH2F*)PRgeomAccfile->Get(nameCSp.str().c_str());
     AccMap_HX_PR = (TH2F*)PRgeomAccfile->Get(nameHXp.str().c_str());
     RecoEffMap_CS_PR = (TH2F*)PRrecoEfffile->Get(nameCSp.str().c_str());
     RecoEffMap_HX_PR = (TH2F*)PRrecoEfffile->Get(nameHXp.str().c_str());
     TrigEffMap_CS_PR = (TH2F*)PRtrigEfffile->Get(nameCSp.str().c_str());
     TrigEffMap_HX_PR = (TH2F*)PRtrigEfffile->Get(nameHXp.str().c_str());
     AccMap_CS_NP = (TH2F*)NPgeomAccfile->Get(nameCSp.str().c_str());
     AccMap_HX_NP = (TH2F*)NPgeomAccfile->Get(nameHXp.str().c_str());
     RecoEffMap_CS_NP = (TH2F*)NPrecoEfffile->Get(nameCSp.str().c_str());
     RecoEffMap_HX_NP = (TH2F*)NPrecoEfffile->Get(nameHXp.str().c_str());
     TrigEffMap_CS_NP = (TH2F*)NPtrigEfffile->Get(nameCSp.str().c_str());
     TrigEffMap_HX_NP = (TH2F*)NPtrigEfffile->Get(nameHXp.str().c_str());

	 cout<<"Extracted maps from data"<<endl;

/*     TH2F *Corr_CS_PR = (TH2F*) AccMap_CS_PR->Clone("Corr_CS_PR");
     TH2F *Corr_HX_PR = (TH2F*) AccMap_HX_PR->Clone("Corr_HX_PR");
     TH2F *Corr_CS_NP = (TH2F*) AccMap_CS_NP->Clone("Corr_CS_NP");
     TH2F *Corr_HX_NP = (TH2F*) AccMap_HX_NP->Clone("Corr_HX_NP");

    RecoEffMap_CS_PR->Multiply(TrigEffMap_CS_PR);
    Corr_CS_PR->Multiply(RecoEffMap_CS_PR);
    RecoEffMap_HX_PR->Multiply(TrigEffMap_HX_PR);
    Corr_HX_PR->Multiply(RecoEffMap_HX_PR);
    RecoEffMap_CS_NP->Multiply(TrigEffMap_CS_NP);
    Corr_CS_NP->Multiply(RecoEffMap_CS_NP);
    RecoEffMap_HX_NP->Multiply(TrigEffMap_HX_NP);
    Corr_HX_NP->Multiply(RecoEffMap_HX_NP);
*/


//////////////////// CUT TURN-ONS //////////////////////////


/*	  double cut=0.7;


	  for(int costhbin=1;costhbin<TrigEffMap_CS_PR->GetNbinsX()+1;costhbin++){
	      for(int phibin=1;phibin<TrigEffMap_CS_PR->GetNbinsY()+1;phibin++){
	    	  if(TrigEffMap_CS_PR->GetBinContent(costhbin,phibin)<cut) TrigEffMap_CS_PR->SetBinContent(costhbin,phibin,0);
	    	  if(TrigEffMap_CS_NP->GetBinContent(costhbin,phibin)<cut) TrigEffMap_CS_NP->SetBinContent(costhbin,phibin,0);
	    	  if(TrigEffMap_HX_PR->GetBinContent(costhbin,phibin)<cut) TrigEffMap_HX_PR->SetBinContent(costhbin,phibin,0);
	    	  if(TrigEffMap_HX_NP->GetBinContent(costhbin,phibin)<cut) TrigEffMap_HX_NP->SetBinContent(costhbin,phibin,0);
	      }}
*/




//////////////////// BUILD TEMPLATES ///////////////////////

    JpsiMass.setBins(100);
    Jpsict.setBins(100);
    JpsictErr.setBins(100);
    costh_CS.setBins(AccMap_CS_PR->GetNbinsX());
    costh_HX.setBins(AccMap_CS_PR->GetNbinsX());
    phi_CS.setBins(AccMap_CS_PR->GetNbinsY());
    phi_HX.setBins(AccMap_CS_PR->GetNbinsY());

    int linsmooth = 0;

    RooDataHist* PRmassHist = new RooDataHist("PRmassHist","PRmassHist",RooArgSet(JpsiMass),*dataPRbin);
    RooDataHist* PRlifetimeHist = new RooDataHist("PRlifetimeHist","PRlifetimeHist",RooArgSet(Jpsict),*dataPRbin);
    RooDataHist* PRlifetimeErrHist = new RooDataHist("PRlifetimeErrHist","PRlifetimeErrHist",RooArgSet(JpsictErr),*dataPRbin);

    RooHistPdf* PRmassHistPdf = new RooHistPdf("PRmassHistPdf","PRmassHistPdf",RooArgSet(JpsiMass),*PRmassHist,1);
    RooHistPdf* PRlifetimeHistPdf = new RooHistPdf("PRlifetimeHistPdf","PRlifetimeHistPdf",RooArgSet(Jpsict),*PRlifetimeHist,1);
    RooHistPdf* PRlifetimeErrHistPdf = new RooHistPdf("PRlifetimeErrHistPdf","PRlifetimeErrHistPdf",RooArgSet(JpsictErr),*PRlifetimeErrHist,1);

    RooDataHist* NPmassHist = new RooDataHist("NPmassHist","NPmassHist",RooArgSet(JpsiMass),*dataNPbin);
    RooDataHist* NPlifetimeHist = new RooDataHist("NPlifetimeHist","NPlifetimeHist",RooArgSet(Jpsict),*dataNPbin);
    RooDataHist* NPlifetimeErrHist = new RooDataHist("NPlifetimeErrHist","NPlifetimeErrHist",RooArgSet(JpsictErr),*dataNPbin);

    RooHistPdf* NPmassHistPdf = new RooHistPdf("NPmassHistPdf","NPmassHistPdf",RooArgSet(JpsiMass),*NPmassHist,1);
    RooHistPdf* NPlifetimeHistPdf = new RooHistPdf("NPlifetimeHistPdf","NPlifetimeHistPdf",RooArgSet(Jpsict),*NPlifetimeHist,1);
    RooHistPdf* NPlifetimeErrHistPdf = new RooHistPdf("NPlifetimeErrHistPdf","NPlifetimeErrHistPdf",RooArgSet(JpsictErr),*NPlifetimeErrHist,1);
/*
    RooDataHist* PR_CSHist = new RooDataHist("PR_CSHist","PR_CSHist",RooArgList(costh_CS,phi_CS),Corr_CS_PR);
    RooDataHist* PR_HXHist = new RooDataHist("PR_HXHist","PR_HXHist",RooArgList(costh_HX,phi_HX),Corr_HX_PR);
    RooDataHist* NP_CSHist = new RooDataHist("NP_CSHist","NP_CSHist",RooArgList(costh_CS,phi_CS),Corr_CS_NP);
    RooDataHist* NP_HXHist = new RooDataHist("NP_HXHist","NP_HXHist",RooArgList(costh_HX,phi_HX),Corr_HX_NP);

    RooHistPdf* PR_CSHistPdf = new RooHistPdf("PR_CSHistPdf","PR_CSHistPdf",RooArgSet(costh_CS,phi_CS),*PR_CSHist,1);
    RooHistPdf* PR_HXHistPdf = new RooHistPdf("PR_HXHistPdf","PR_HXHistPdf",RooArgSet(costh_HX,phi_HX),*PR_HXHist,1);
    RooHistPdf* NP_CSHistPdf = new RooHistPdf("NP_CSHistPdf","NP_CSHistPdf",RooArgSet(costh_CS,phi_CS),*NP_CSHist,1);
    RooHistPdf* NP_HXHistPdf = new RooHistPdf("NP_HXHistPdf","NP_HXHistPdf",RooArgSet(costh_HX,phi_HX),*NP_HXHist,1);
*/
	 cout<<"Mass and Lifetime built"<<endl;


    RooDataHist* PR_CS_geomAccHist = new RooDataHist("PR_CS_geomAccHist","PR_CS_geomAccHist",RooArgList(costh_CS,phi_CS),AccMap_CS_PR);
    RooDataHist* PR_HX_geomAccHist = new RooDataHist("PR_HX_geomAccHist","PR_HX_geomAccHist",RooArgList(costh_HX,phi_HX),AccMap_HX_PR);
    RooDataHist* NP_CS_geomAccHist = new RooDataHist("NP_CS_geomAccHist","NP_CS_geomAccHist",RooArgList(costh_CS,phi_CS),AccMap_CS_NP);
    RooDataHist* NP_HX_geomAccHist = new RooDataHist("NP_HX_geomAccHist","NP_HX_geomAccHist",RooArgList(costh_HX,phi_HX),AccMap_HX_NP);

    RooHistPdf* PR_CS_geomAccHistPdf = new RooHistPdf("PR_CS_geomAccHistPdf","PR_CS_geomAccHistPdf",RooArgSet(costh_CS,phi_CS),*PR_CS_geomAccHist,linsmooth);
    RooHistPdf* PR_HX_geomAccHistPdf = new RooHistPdf("PR_HX_geomAccHistPdf","PR_HX_geomAccHistPdf",RooArgSet(costh_HX,phi_HX),*PR_HX_geomAccHist,linsmooth);
    RooHistPdf* NP_CS_geomAccHistPdf = new RooHistPdf("NP_CS_geomAccHistPdf","NP_CS_geomAccHistPdf",RooArgSet(costh_CS,phi_CS),*NP_CS_geomAccHist,linsmooth);
    RooHistPdf* NP_HX_geomAccHistPdf = new RooHistPdf("NP_HX_geomAccHistPdf","NP_HX_geomAccHistPdf",RooArgSet(costh_HX,phi_HX),*NP_HX_geomAccHist,linsmooth);

    RooDataHist* PR_CS_recoEffHist = new RooDataHist("PR_CS_recoEffHist","PR_CS_recoEffHist",RooArgList(costh_CS,phi_CS),RecoEffMap_CS_PR);
    RooDataHist* PR_HX_recoEffHist = new RooDataHist("PR_HX_recoEffHist","PR_HX_recoEffHist",RooArgList(costh_HX,phi_HX),RecoEffMap_HX_PR);
    RooDataHist* NP_CS_recoEffHist = new RooDataHist("NP_CS_recoEffHist","NP_CS_recoEffHist",RooArgList(costh_CS,phi_CS),RecoEffMap_CS_NP);
    RooDataHist* NP_HX_recoEffHist = new RooDataHist("NP_HX_recoEffHist","NP_HX_recoEffHist",RooArgList(costh_HX,phi_HX),RecoEffMap_HX_NP);

    RooHistPdf* PR_CS_recoEffHistPdf = new RooHistPdf("PR_CS_recoEffHistPdf","PR_CS_recoEffHistPdf",RooArgSet(costh_CS,phi_CS),*PR_CS_recoEffHist,linsmooth);
    RooHistPdf* PR_HX_recoEffHistPdf = new RooHistPdf("PR_HX_recoEffHistPdf","PR_HX_recoEffHistPdf",RooArgSet(costh_HX,phi_HX),*PR_HX_recoEffHist,linsmooth);
    RooHistPdf* NP_CS_recoEffHistPdf = new RooHistPdf("NP_CS_recoEffHistPdf","NP_CS_recoEffHistPdf",RooArgSet(costh_CS,phi_CS),*NP_CS_recoEffHist,linsmooth);
    RooHistPdf* NP_HX_recoEffHistPdf = new RooHistPdf("NP_HX_recoEffHistPdf","NP_HX_recoEffHistPdf",RooArgSet(costh_HX,phi_HX),*NP_HX_recoEffHist,linsmooth);

    RooDataHist* PR_CS_trigEffHist = new RooDataHist("PR_CS_trigEffHist","PR_CS_trigEffHist",RooArgList(costh_CS,phi_CS),TrigEffMap_CS_PR);
    RooDataHist* PR_HX_trigEffHist = new RooDataHist("PR_HX_trigEffHist","PR_HX_trigEffHist",RooArgList(costh_HX,phi_HX),TrigEffMap_HX_PR);
    RooDataHist* NP_CS_trigEffHist = new RooDataHist("NP_CS_trigEffHist","NP_CS_trigEffHist",RooArgList(costh_CS,phi_CS),TrigEffMap_CS_NP);
    RooDataHist* NP_HX_trigEffHist = new RooDataHist("NP_HX_trigEffHist","NP_HX_trigEffHist",RooArgList(costh_HX,phi_HX),TrigEffMap_HX_NP);

    RooHistPdf* PR_CS_trigEffHistPdf = new RooHistPdf("PR_CS_trigEffHistPdf","PR_CS_trigEffHistPdf",RooArgSet(costh_CS,phi_CS),*PR_CS_trigEffHist,linsmooth);
    RooHistPdf* PR_HX_trigEffHistPdf = new RooHistPdf("PR_HX_trigEffHistPdf","PR_HX_trigEffHistPdf",RooArgSet(costh_HX,phi_HX),*PR_HX_trigEffHist,linsmooth);
    RooHistPdf* NP_CS_trigEffHistPdf = new RooHistPdf("NP_CS_trigEffHistPdf","NP_CS_trigEffHistPdf",RooArgSet(costh_CS,phi_CS),*NP_CS_trigEffHist,linsmooth);
    RooHistPdf* NP_HX_trigEffHistPdf = new RooHistPdf("NP_HX_trigEffHistPdf","NP_HX_trigEffHistPdf",RooArgSet(costh_HX,phi_HX),*NP_HX_trigEffHist,linsmooth);

	 cout<<"Costh/phi built"<<endl;

//////////////////// BUILD MODEL ///////////////////////

    RooRealVar PR_lam_th_CS("PR_lam_th_CS","PR_lam_th_CS",0);
    RooRealVar PR_lam_ph_CS("PR_lam_ph_CS","PR_lam_ph_CS",0);
    RooRealVar PR_lam_thph_CS("PR_lam_thph_CS","PR_lam_thph_CS",0);
    RooRealVar PR_lam_th_HX("PR_lam_th_HX","PR_lam_th_HX",0);
    RooRealVar PR_lam_ph_HX("PR_lam_ph_HX","PR_lam_ph_HX",0);
    RooRealVar PR_lam_thph_HX("PR_lam_thph_HX","PR_lam_thph_HX",0);

    RooRealVar NP_lam_th_CS("NP_lam_th_CS","NP_lam_th_CS",0);
    RooRealVar NP_lam_ph_CS("NP_lam_ph_CS","NP_lam_ph_CS",0);
    RooRealVar NP_lam_thph_CS("NP_lam_thph_CS","NP_lam_thph_CS",0);
    RooRealVar NP_lam_th_HX("NP_lam_th_HX","NP_lam_th_HX",0);
    RooRealVar NP_lam_ph_HX("NP_lam_ph_HX","NP_lam_ph_HX",0);
    RooRealVar NP_lam_thph_HX("NP_lam_thph_HX","NP_lam_thph_HX",0);

    RooPolarizationPdf* PR_CS_PolPdf = new RooPolarizationPdf("PR_CS_PolPdf","PR_CS_PolPdf",costh_CS,phi_CS,PR_lam_th_CS,PR_lam_ph_CS,PR_lam_thph_CS);
    RooPolarizationPdf* PR_HX_PolPdf = new RooPolarizationPdf("PR_HX_PolPdf","PR_HX_PolPdf",costh_HX,phi_HX,PR_lam_th_HX,PR_lam_ph_HX,PR_lam_thph_HX);
    RooPolarizationPdf* NP_CS_PolPdf = new RooPolarizationPdf("NP_CS_PolPdf","NP_CS_PolPdf",costh_CS,phi_CS,NP_lam_th_CS,NP_lam_ph_CS,NP_lam_thph_CS);
    RooPolarizationPdf* NP_HX_PolPdf = new RooPolarizationPdf("NP_HX_PolPdf","NP_HX_PolPdf",costh_HX,phi_HX,NP_lam_th_HX,NP_lam_ph_HX,NP_lam_thph_HX);


    RooArgList prodListPR;
//    prodListPR.add(*PRmassHistPdf);
//    prodListPR.add(*PRlifetimeHistPdf);
//    prodListPR.add(*PRlifetimeErrHistPdf);

    prodListPR.add(*PR_CS_geomAccHistPdf);
    prodListPR.add(*PR_CS_recoEffHistPdf);
    prodListPR.add(*PR_CS_trigEffHistPdf);
    prodListPR.add(*PR_CS_PolPdf);
    prodListPR.add(*PR_HX_geomAccHistPdf);
    prodListPR.add(*PR_HX_recoEffHistPdf);
    prodListPR.add(*PR_HX_trigEffHistPdf);
    prodListPR.add(*PR_HX_PolPdf);

    RooProdPdf* PRpdf = new RooProdPdf("PRpdf","PRpdf",prodListPR);

    RooArgList prodListNP;
//    prodListNP.add(*NPmassHistPdf);
//    prodListNP.add(*NPlifetimeHistPdf);
//    prodListNP.add(*NPlifetimeErrHistPdf);

    prodListNP.add(*NP_CS_geomAccHistPdf);
    prodListNP.add(*NP_CS_recoEffHistPdf);
    prodListNP.add(*NP_CS_trigEffHistPdf);
    prodListNP.add(*NP_CS_PolPdf);
    prodListNP.add(*NP_HX_geomAccHistPdf);
    prodListNP.add(*NP_HX_recoEffHistPdf);
    prodListNP.add(*NP_HX_trigEffHistPdf);
    prodListNP.add(*NP_HX_PolPdf);

    RooProdPdf* NPpdf = new RooProdPdf("NPpdf","NPpdf",prodListNP);


//////////////////// DEFINE POLARIZATION AND B-FRACTION ///////////////////////
	 cout<<"Define Polarization"<<endl;


    PR_lam_th_CS.setVal(jpsi::scenarioLambdas[scenario-1][0]);
    PR_lam_ph_CS.setVal(jpsi::scenarioLambdas[scenario-1][1]);
    PR_lam_thph_CS.setVal(jpsi::scenarioLambdas[scenario-1][2]);
    NP_lam_th_CS.setVal(jpsi::scenarioLambdas[scenario-1][3]);
    NP_lam_ph_CS.setVal(jpsi::scenarioLambdas[scenario-1][4]);
    NP_lam_thph_CS.setVal(jpsi::scenarioLambdas[scenario-1][5]);

    PR_lam_th_HX.setVal(jpsi::scenarioLambdas[scenario-1][0]);
    PR_lam_ph_HX.setVal(jpsi::scenarioLambdas[scenario-1][1]);
    PR_lam_thph_HX.setVal(jpsi::scenarioLambdas[scenario-1][2]);
    NP_lam_th_HX.setVal(jpsi::scenarioLambdas[scenario-1][3]);
    NP_lam_ph_HX.setVal(jpsi::scenarioLambdas[scenario-1][4]);
    NP_lam_thph_HX.setVal(jpsi::scenarioLambdas[scenario-1][5]);

    cout<<jpsi::scenarioLambdas[scenario-1][0]<<endl;
    cout<<jpsi::scenarioLambdas[scenario-1][1]<<endl;
    cout<<jpsi::scenarioLambdas[scenario-1][2]<<endl;
    cout<<jpsi::scenarioLambdas[scenario-1][3]<<endl;
    cout<<jpsi::scenarioLambdas[scenario-1][4]<<endl;
    cout<<jpsi::scenarioLambdas[scenario-1][5]<<endl;


    double numEvents = jpsi::numEventsBin_DoubleMu0[yBin][ptBin];
    double numEventsPRregion = numEvents*jpsi::fracPRregion[yBin][ptBin];
    double numEventsNPregion = numEvents*(1-jpsi::fracPRregion[yBin][ptBin]);

    double numEventsPRinPR = numEventsPRregion*jpsi::fPinP_actual[yBin][ptBin];
    double numEventsPRinNP = numEventsNPregion*jpsi::fPinNP_actual[yBin][ptBin];
    double numEventsNPinPR = numEventsPRregion*jpsi::fNPinP_actual[yBin][ptBin];
    double numEventsNPinNP = numEventsNPregion*jpsi::fNPinNP_actual[yBin][ptBin];

//////////////////// GENERATE DATASET ///////////////////////

	RooAbsData *thisBinGenPRinPR = PRpdf->generate(RooArgSet(/*JpsiMass,Jpsict,JpsictErr,*/costh_CS,phi_CS,costh_HX,phi_HX),numEventsPRinPR);
	RooAbsData *thisBinGenNPinPR = NPpdf->generate(RooArgSet(/*JpsiMass,Jpsict,JpsictErr,*/costh_CS,phi_CS,costh_HX,phi_HX),numEventsNPinPR);

	RooAbsData *thisBinGenPRinNP = PRpdf->generate(RooArgSet(/*JpsiMass,Jpsict,JpsictErr,*/costh_CS,phi_CS,costh_HX,phi_HX),numEventsPRinNP);
	RooAbsData *thisBinGenNPinNP = NPpdf->generate(RooArgSet(/*JpsiMass,Jpsict,JpsictErr,*/costh_CS,phi_CS,costh_HX,phi_HX),numEventsNPinNP);

	cout<<"Generation finished"<<endl;
	char title[200];
	sprintf(title,"PRinPR_scen%d",scenario);
	thisBinGenPRinPR->SetTitle(title);
	sprintf(title,"NPinPR_scen%d",scenario);
	thisBinGenNPinPR->SetTitle(title);
	sprintf(title,"PRinNP_scen%d",scenario);
	thisBinGenPRinNP->SetTitle(title);
	sprintf(title,"NPinNP_scen%d",scenario);
	thisBinGenNPinNP->SetTitle(title);

	thisBinGenPRinPR->SetName(JOBNAME);
	thisBinGenPRinNP->SetName(JOBNAME);
	thisBinGenNPinPR->SetName(JOBNAME);
	thisBinGenNPinNP->SetName(JOBNAME);


	cout<<"Start Plotting"<<endl;


//	Plot(*thisBinGenPR,JpsiMass,*PRpdf,true);
//  Plot(*thisBinGenPR,Jpsict,*PRpdf,true);
//	Plot(*thisBinGenPR,JpsictErr,*PRpdf,true);
//	Plot(*thisBinGenNP,JpsiMass,*NPpdf,false);
//	Plot(*thisBinGenNP,Jpsict,*NPpdf,false);
//	Plot(*thisBinGenNP,JpsictErr,*NPpdf,false);

	if(plot){
	Plot(*thisBinGenPRinPR,costh_CS,*PRpdf,true);
    Plot(*thisBinGenPRinPR,costh_HX,*PRpdf,true);
    Plot(*thisBinGenPRinPR,phi_CS,*PRpdf,true);
    Plot(*thisBinGenPRinPR,phi_HX,*PRpdf,true);

    Plot(*thisBinGenNPinPR,costh_CS,*NPpdf,false);
    Plot(*thisBinGenNPinPR,costh_HX,*NPpdf,false);
    Plot(*thisBinGenNPinPR,phi_CS,*NPpdf,false);
    Plot(*thisBinGenNPinPR,phi_HX,*NPpdf,false);

	Plot(*thisBinGenPRinNP,costh_CS,*PRpdf,true);
    Plot(*thisBinGenPRinNP,costh_HX,*PRpdf,true);
    Plot(*thisBinGenPRinNP,phi_CS,*PRpdf,true);
    Plot(*thisBinGenPRinNP,phi_HX,*PRpdf,true);

    Plot(*thisBinGenNPinNP,costh_CS,*NPpdf,false);
    Plot(*thisBinGenNPinNP,costh_HX,*NPpdf,false);
    Plot(*thisBinGenNPinNP,phi_CS,*NPpdf,false);
    Plot(*thisBinGenNPinNP,phi_HX,*NPpdf,false);
	}
//////////////////// DELETE POINTERS (MODEL) ///////////////////////


	delete PRmassHist;
	delete PRlifetimeHist;
	delete PRlifetimeErrHist;
	delete PRmassHistPdf;
	delete PRlifetimeHistPdf;
	delete PRlifetimeErrHistPdf;
	delete NPmassHist;
	delete NPlifetimeHist;
	delete NPlifetimeErrHist;
	delete NPmassHistPdf;
	delete NPlifetimeHistPdf;
	delete NPlifetimeErrHistPdf;
/*	delete PR_CSHist;
	delete PR_HXHist;
	delete NP_CSHist;
	delete NP_HXHist;
	delete PR_CSHistPdf;
	delete PR_HXHistPdf;
	delete NP_CSHistPdf;
	delete NP_HXHistPdf;
*/
	delete PR_CS_geomAccHist;
	delete PR_HX_geomAccHist;
	delete NP_CS_geomAccHist;
	delete NP_HX_geomAccHist;
	delete PR_CS_geomAccHistPdf;
	delete PR_HX_geomAccHistPdf;
	delete NP_CS_geomAccHistPdf;
	delete NP_HX_geomAccHistPdf;
	delete PR_CS_recoEffHist;
	delete PR_HX_recoEffHist;
	delete NP_CS_recoEffHist;
	delete NP_HX_recoEffHist;
	delete PR_CS_recoEffHistPdf;
	delete PR_HX_recoEffHistPdf;
	delete NP_CS_recoEffHistPdf;
	delete NP_HX_recoEffHistPdf;
	delete PR_CS_trigEffHist;
	delete PR_HX_trigEffHist;
	delete NP_CS_trigEffHist;
	delete NP_HX_trigEffHist;
	delete PR_CS_trigEffHistPdf;
	delete PR_HX_trigEffHistPdf;
	delete NP_CS_trigEffHistPdf;
	delete NP_HX_trigEffHistPdf;

	delete PR_CS_PolPdf;
	delete PR_HX_PolPdf;
	delete NP_CS_PolPdf;
	delete NP_HX_PolPdf;
	delete PRpdf;
	delete NPpdf;


//////////////////// ADD CONSTANT VARIABLES ///////////////////////

    TTree* thisBinGenPRinPRtree = (TTree*)thisBinGenPRinPR->tree();
    TTree* thisBinGenNPinPRtree = (TTree*)thisBinGenNPinPR->tree();

    TTree* thisBinGenPRinNPtree = (TTree*)thisBinGenPRinNP->tree();
    TTree* thisBinGenNPinNPtree = (TTree*)thisBinGenNPinNP->tree();

    thisBinGenPRinPRtree->SetName("data");
    thisBinGenNPinPRtree->SetName("data");
    thisBinGenPRinNPtree->SetName("data");
    thisBinGenNPinNPtree->SetName("data");


    Float_t MCType_;
    Float_t pT_pseudo;
    Float_t y_pseudo;
    Float_t generation_;
    Float_t Jpsict_;

     TBranch *MCTypeBranchPRinPR = thisBinGenPRinPRtree->Branch("MCType",&MCType_,"MCType/F");
     TBranch *y_pseudoBranchPRinPR = thisBinGenPRinPRtree->Branch("JpsiRap", &y_pseudo, "JpsiRap/F");
     TBranch *pT_pseudoBranchPRinPR = thisBinGenPRinPRtree->Branch("JpsiPt",&pT_pseudo,"JpsiPt/F");
     TBranch *generationPRinPR = thisBinGenPRinPRtree->Branch("generation",&generation_,"generation/F");
     TBranch *lifetimeBranchPRinPR = thisBinGenPRinPRtree->Branch("Jpsict",&Jpsict_,"Jpsict/F");

     TBranch *MCTypeBranchNPinPR = thisBinGenNPinPRtree->Branch("MCType",&MCType_,"MCType/F");
     TBranch *y_pseudoBranchNPinPR = thisBinGenNPinPRtree->Branch("JpsiRap", &y_pseudo, "JpsiRap/F");
     TBranch *pT_pseudoBranchNPinPR = thisBinGenNPinPRtree->Branch("JpsiPt",&pT_pseudo,"JpsiPt/F");
     TBranch *generationNPinPR = thisBinGenNPinPRtree->Branch("generation",&generation_,"generation/F");
     TBranch *lifetimeBranchNPinPR = thisBinGenNPinPRtree->Branch("Jpsict",&Jpsict_,"Jpsict/F");

     TBranch *MCTypeBranchPRinNP = thisBinGenPRinNPtree->Branch("MCType",&MCType_,"MCType/F");
     TBranch *y_pseudoBranchPRinNP = thisBinGenPRinNPtree->Branch("JpsiRap", &y_pseudo, "JpsiRap/F");
     TBranch *pT_pseudoBranchPRinNP = thisBinGenPRinNPtree->Branch("JpsiPt",&pT_pseudo,"JpsiPt/F");
     TBranch *generationPRinNP = thisBinGenPRinNPtree->Branch("generation",&generation_,"generation/F");
     TBranch *lifetimeBranchPRinNP = thisBinGenPRinNPtree->Branch("Jpsict",&Jpsict_,"Jpsict/F");

     TBranch *MCTypeBranchNPinNP = thisBinGenNPinNPtree->Branch("MCType",&MCType_,"MCType/F");
     TBranch *y_pseudoBranchNPinNP = thisBinGenNPinNPtree->Branch("JpsiRap", &y_pseudo, "JpsiRap/F");
     TBranch *pT_pseudoBranchNPinNP = thisBinGenNPinNPtree->Branch("JpsiPt",&pT_pseudo,"JpsiPt/F");
     TBranch *generationNPinNP = thisBinGenNPinNPtree->Branch("generation",&generation_,"generation/F");
     TBranch *lifetimeBranchNPinNP = thisBinGenNPinNPtree->Branch("Jpsict",&Jpsict_,"Jpsict/F");

    Int_t numEntriesPRinPR = (Int_t)thisBinGenPRinPRtree->GetEntries();
    Int_t numEntriesNPinPR = (Int_t)thisBinGenNPinPRtree->GetEntries();

    Int_t numEntriesPRinNP = (Int_t)thisBinGenPRinNPtree->GetEntries();
    Int_t numEntriesNPinNP = (Int_t)thisBinGenNPinNPtree->GetEntries();

    for(Int_t j=0; j<numEntriesPRinPR; ++j)
	{
      thisBinGenPRinPRtree->GetEntry(j);
      MCType_=0;
      MCTypeBranchPRinPR->Fill();
      pT_pseudo = (jpsi::pTRange[yBin+1][ptBin]+jpsi::pTRange[yBin+1][ptBin+1])/2;
	  pT_pseudoBranchPRinPR->Fill();
      y_pseudo=(jpsi::rapForPTRange[yBin]+jpsi::rapForPTRange[yBin+1])/2;
      y_pseudoBranchPRinPR->Fill();
      generation_=generation;
      generationPRinPR->Fill();
      Jpsict_=0;
      lifetimeBranchPRinPR->Fill();
	}

    for(Int_t j=0; j<numEntriesNPinPR; ++j)
 	{
      thisBinGenNPinPRtree->GetEntry(j);
      MCType_=1;
      MCTypeBranchNPinPR->Fill();
      pT_pseudo = (jpsi::pTRange[yBin+1][ptBin]+jpsi::pTRange[yBin+1][ptBin+1])/2;
  	  pT_pseudoBranchNPinPR->Fill();
      y_pseudo=(jpsi::rapForPTRange[yBin]+jpsi::rapForPTRange[yBin+1])/2;
      y_pseudoBranchNPinPR->Fill();
      generation_=generation;
      generationNPinPR->Fill();
      Jpsict_=0;
      lifetimeBranchNPinPR->Fill();
 	}

    for(Int_t j=0; j<numEntriesPRinNP; ++j)
	{
      thisBinGenPRinNPtree->GetEntry(j);
      MCType_=0;
      MCTypeBranchPRinNP->Fill();
      pT_pseudo = (jpsi::pTRange[yBin+1][ptBin]+jpsi::pTRange[yBin+1][ptBin+1])/2;
	  pT_pseudoBranchPRinNP->Fill();
      y_pseudo=(jpsi::rapForPTRange[yBin]+jpsi::rapForPTRange[yBin+1])/2;
      y_pseudoBranchPRinNP->Fill();
      generation_=generation;
      generationPRinNP->Fill();
      Jpsict_=1;
      lifetimeBranchPRinNP->Fill();
	}

    for(Int_t j=0; j<numEntriesNPinNP; ++j)
 	{
      thisBinGenNPinNPtree->GetEntry(j);
      MCType_=1;
      MCTypeBranchNPinNP->Fill();
      pT_pseudo = (jpsi::pTRange[yBin+1][ptBin]+jpsi::pTRange[yBin+1][ptBin+1])/2;
  	  pT_pseudoBranchNPinNP->Fill();
      y_pseudo=(jpsi::rapForPTRange[yBin]+jpsi::rapForPTRange[yBin+1])/2;
      y_pseudoBranchNPinNP->Fill();
      generation_=generation;
      generationNPinNP->Fill();
      Jpsict_=1;
      lifetimeBranchNPinNP->Fill();
 	}
//////////////////// MERGE PR AND NP ///////////////////////

    TList list_prnp;

    list_prnp.Add(thisBinGenPRinPRtree);
    list_prnp.Add(thisBinGenNPinPRtree);
    list_prnp.Add(thisBinGenPRinNPtree);
    list_prnp.Add(thisBinGenNPinNPtree);

    sprintf(filenamePRNP,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/%s_tempPRNP.root",JOBNAME);
    TFile *f = new TFile(filenamePRNP,"RECREATE");

    sprintf(filenameGEN,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/%s_tempGEN.root",JOBNAME);

    TTree::MergeTrees(&list_prnp);

    f->Write();
    f->Close();

	delete f;


//////////////////// MERGE GENERATIONS ///////////////////////

	TFile *gen0 = new TFile(filenameGEN,"UPDATE");
	if(gen0->Get("data")==NULL){
	cout<<"NULL"<<endl;
	gen0->Close();
	delete gen0;
	gSystem->CopyFile(filenamePRNP,filenameGEN,true);
	}
	else{
	cout<<"merge Gen "<<generation<<endl;
	gen0->Close();
	delete gen0;
    TFile *prnp = new TFile(filenamePRNP,"UPDATE");
	TTree* tree_prnp_existing = (TTree*)prnp->Get("data");

	TFile *gen1 = new TFile(filenameGEN,"UPDATE");
	TTree* tree_gen_existing = (TTree*)gen1->Get("data");

	TFile *gen = new TFile(filenameGEN,"RECREATE");
	TList list_gen;
	list_gen.Add(tree_gen_existing);
	list_gen.Add(tree_prnp_existing);
	TTree::MergeTrees(&list_gen);
	gen->Write();
	gen->Close();
	delete gen;

	prnp->Close();
	delete prnp;
	gen1->Close();
	delete gen1;
	}

//////////////////// MERGE PT AND Y ///////////////////////

/*    TFile *pty0 = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root","UPDATE");
	if(pty0->Get("data")==NULL){
	cout<<"NULL"<<endl;
	pty0->Close();
	delete pty0;
    gSystem->CopyFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempPRNP.root","/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root",true);
	}
	else{
	cout<<"merge y/pt"<<endl;
	pty0->Close();
	delete pty0;
    TFile *prnp = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempPRNP.root","UPDATE");
	TTree* tree_prnp_existing = (TTree*)prnp->Get("data");

    TFile *pty1 = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root","UPDATE");
    TTree* tree_pty_existing = (TTree*)pty1->Get("data");

    TFile *pty = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root","RECREATE");
    TList list_pty;
	list_pty.Add(tree_pty_existing);
	list_pty.Add(tree_prnp_existing);
	TTree::MergeTrees(&list_pty);
	pty->Write();
	pty->Close();
	delete pty;

	prnp->Close();
	delete prnp;
	pty1->Close();
	delete pty1;
	}
*/

//////////////////// END ITERATION LOOP ///////////////////////

	  }


		  char dirstruct2[200];
		  sprintf(dirstruct2,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/scenario%d",scenario);
		  gSystem->mkdir(dirstruct2);
		  char outputname[200];
		  sprintf(outputname,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/saveTrees/TTree_%s_rap%d_pT%d.root",GenID,yBin+1,ptBin+1);


		  gSystem->Rename(filenameGEN,outputname);
		  gSystem->Rename(filenamePRNP,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/temptempPRNP.root");

//////////////////// END PT LOOP ///////////////////////


}


//////////////////// END Y LOOP ///////////////////////

//////////////////// MERGE GENERATIONS ///////////////////////
/*
	TFile *gen0 = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempGEN.root","UPDATE");
	if(gen0->Get("data")==NULL){
	cout<<"NULL"<<endl;
	gen0->Close();
	delete gen0;
	gSystem->CopyFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root","/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempGEN.root",true);
	}
	else{
	cout<<"merge Gen "<<generation<<endl;
	gen0->Close();
	delete gen0;
	TFile *pty2 = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempYPT.root","UPDATE");
	TTree* tree_pty2_existing = (TTree*)pty2->Get("data");

	TFile *gen1 = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempGEN.root","UPDATE");
	TTree* tree_gen_existing = (TTree*)gen1->Get("data");

	TFile *gen = new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/tempGEN.root","RECREATE");
	TList list_gen;
	list_gen.Add(tree_gen_existing);
	list_gen.Add(tree_pty2_existing);
	TTree::MergeTrees(&list_gen);
	gen->Write();
	gen->Close();
	delete gen;

	pty2->Close();
	delete pty2;
	gen1->Close();
	delete gen1;
}
*/
//////////////////// CHANGE FORMAT ///////////////////////


//////////////////// SAVE ///////////////////////


//////////////////// MOVE ON TO AN OTHER GENERATION :) ///////////////////////

}





delete dataPR;
delete dataNP;

PRfile->Close();
NPfile->Close();

PRgeomAccfile->Close();
PRrecoEfffile->Close();
PRtrigEfffile->Close();

NPgeomAccfile->Close();
NPrecoEfffile->Close();
NPtrigEfffile->Close();

fileforwhateverrootorcplusplusneedsme->Close();

	return 0;
}


