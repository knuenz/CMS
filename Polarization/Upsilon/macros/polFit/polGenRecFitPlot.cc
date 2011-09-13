#include "../../interface/rootIncludes.inc"
#include "../../interface/commonVar.h"
#include "ToyMC.h"
#include "effsAndCuts.h"

#include "TSystem.h"
#include "TROOT.h"
#include "polGen.C"
#include "polRec.C"
#include "polFit.C"
#include "polPlot.C"

#include <time.h>

//====================================

int main(int argc, char** argv) {


	int nGenerations=999;
	int polScenSig=999;
	int frameSig=999;
	int polScenBkg=999;
	int frameBkg=999;
	int rapBinMin=999;
	int rapBinMax=999;
	int ptBinMin=999;
	int ptBinMax=999;
	int nEff=999;
	int FidCuts=999;
	int nSample=999;
	int ConstEvents=999;
	int nSkipGen=999;
	int ThisGen=999;

	bool ConstEvents_(false);
	bool gen(false);
	bool rec(false);
	bool fit(false);
	bool plot(false);
	bool RealData(false);

	Char_t *storagedir = "Default"; //Storage Directory
	Char_t *basedir = "Default"; //Code Directory
  	Char_t *JobID = "Default";
	Char_t *realdatadir = "Default"; //Storage Directory
	Char_t *TreeID = "ToyMC"; //Storage Directory

	  for( int i=0;i < argc; ++i ) {
	    if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
	    if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
	    if(std::string(argv[i]).find("rapBinMin") != std::string::npos) {char* rapBinMinchar = argv[i]; char* rapBinMinchar2 = strtok (rapBinMinchar, "p"); rapBinMin = atof(rapBinMinchar2); cout<<"rapBinMin = "<<rapBinMin<<endl;}
	    if(std::string(argv[i]).find("rapBinMax") != std::string::npos) {char* rapBinMaxchar = argv[i]; char* rapBinMaxchar2 = strtok (rapBinMaxchar, "p"); rapBinMax = atof(rapBinMaxchar2); cout<<"rapBinMax = "<<rapBinMax<<endl;}
	    if(std::string(argv[i]).find("nGenerations") != std::string::npos) {char* nGenerationschar = argv[i]; char* nGenerationschar2 = strtok (nGenerationschar, "p"); nGenerations = atof(nGenerationschar2); cout<<"nGenerations = "<<nGenerations<<endl;}
	    if(std::string(argv[i]).find("frameSig") != std::string::npos) {char* framecharSig = argv[i]; char* framecharSig2 = strtok (framecharSig, "p"); frameSig = atof(framecharSig2); cout<<"frameSig = "<<frameSig<<endl;}
	    if(std::string(argv[i]).find("polScenSig") != std::string::npos) {char* polScencharSig = argv[i]; char* polScencharSig2 = strtok (polScencharSig, "p"); polScenSig = atof(polScencharSig2); cout<<"polScenSig = "<<polScenSig<<endl;}
	    if(std::string(argv[i]).find("frameBkg") != std::string::npos) {char* framecharBkg = argv[i]; char* framecharBkg2 = strtok (framecharBkg, "p"); frameBkg = atof(framecharBkg2); cout<<"frameBkg = "<<frameBkg<<endl;}
	    if(std::string(argv[i]).find("polScenBkg") != std::string::npos) {char* polScencharBkg = argv[i]; char* polScencharBkg2 = strtok (polScencharBkg, "p"); polScenBkg = atof(polScencharBkg2); cout<<"polScenBkg = "<<polScenBkg<<endl;}
	    if(std::string(argv[i]).find("nEff") != std::string::npos) {char* nEffchar = argv[i]; char* nEffchar2 = strtok (nEffchar, "p"); nEff = atof(nEffchar2); cout<<"nEff = "<<nEff<<endl;}
	    if(std::string(argv[i]).find("FidCuts") != std::string::npos) {char* FidCutschar = argv[i]; char* FidCutschar2 = strtok (FidCutschar, "p"); FidCuts = atof(FidCutschar2); cout<<"FidCuts = "<<FidCuts<<endl;}
	    if(std::string(argv[i]).find("nSample") != std::string::npos) {char* nSamplechar = argv[i]; char* nSamplechar2 = strtok (nSamplechar, "p"); nSample = atof(nSamplechar2); cout<<"nSample = "<<nSample<<endl;}
	    if(std::string(argv[i]).find("ConstEvents") != std::string::npos) {char* ConstEventschar = argv[i]; char* ConstEventschar2 = strtok (ConstEventschar, "p"); ConstEvents = atof(ConstEventschar2); cout<<"ConstEvents = "<<ConstEvents<<endl;}
	    if(std::string(argv[i]).find("nSkipGen") != std::string::npos) {char* nSkipGenchar = argv[i]; char* nSkipGenchar2 = strtok (nSkipGenchar, "p"); nSkipGen = atof(nSkipGenchar2); cout<<"nSkipGen = "<<nSkipGen<<endl;}
	    if(std::string(argv[i]).find("ThisGen") != std::string::npos) {char* ThisGenchar = argv[i]; char* ThisGenchar2 = strtok (ThisGenchar, "p"); ThisGen = atof(ThisGenchar2); cout<<"ThisGen = "<<ThisGen<<endl;}

	    if(std::string(argv[i]).find("UseConstEv=true") != std::string::npos) {ConstEvents_=true; cout<<"use constant number of reconstructed events"<<endl;}
	    if(std::string(argv[i]).find("gen=true") != std::string::npos) {gen=true; cout<<"run polGen.C"<<endl;}
	    if(std::string(argv[i]).find("rec=true") != std::string::npos) {rec=true; cout<<"run polRec.C"<<endl;}
	    if(std::string(argv[i]).find("fit=true") != std::string::npos) {fit=true; cout<<"run polFit.C"<<endl;}
	    if(std::string(argv[i]).find("plot=true") != std::string::npos) {plot=true; cout<<"run polPlot.C"<<endl;}

	    if(std::string(argv[i]).find("JobID") != std::string::npos) {char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<<"JobID = "<<JobID<<endl;}
	    if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}
	    if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}

	    if(std::string(argv[i]).find("realdatadir") != std::string::npos) {RealData=true; cout<<"Fit to real data"<<endl; char* realdatadirchar = argv[i]; char* realdatadirchar2 = strtok (realdatadirchar, "="); realdatadir = realdatadirchar2; cout<<"realdatadir = "<<realdatadir<<endl;}
	    if(std::string(argv[i]).find("TreeID") != std::string::npos) {char* TreeIDchar = argv[i]; char* TreeIDchar2 = strtok (TreeIDchar, "="); TreeID = TreeIDchar2; cout<<"TreeID = "<<TreeID<<endl;}

	    }

  		double lambda_theta_sig_ = ToyMC::ScenarioSig[0][polScenSig-1];
  		double lambda_phi_sig_ = ToyMC::ScenarioSig[1][polScenSig-1];
  		double lambda_thetaphi_sig_ = ToyMC::ScenarioSig[2][polScenSig-1];

  		double lambda_theta_bkg_ = ToyMC::ScenarioBkg[0][polScenBkg-1];
  		double lambda_phi_bkg_ = ToyMC::ScenarioBkg[1][polScenBkg-1];
  		double lambda_thetaphi_bkg_ = ToyMC::ScenarioBkg[2][polScenBkg-1];

		double mass_signal_peak  =  9.5;
		double mass_signal_sigma =  0.1;
		double n_sigmas_signal = 3.;


  	  	Char_t *OutputDirectory;
  	  	Char_t *TreeBinID;
  	  	double f_BG;
  	 	int n_events;
  		char basestruct[200],substruct[200], dirstruct[200], rapptstruct[200], filenameFrom[200], filenameTo[200] , tmpfilename[200], TreeBinID_[200];

  		sprintf(basestruct,"%s/%s",storagedir,JobID);gSystem->mkdir(basestruct);
  		sprintf(substruct,"%s/Sig_frame%dscen%d_Bkg_frame%dscen%d",basestruct,frameSig,polScenSig,frameBkg,polScenBkg); if(!RealData) gSystem->mkdir(substruct);

  		time_t seconds; seconds = time (NULL); double time_0=seconds; double time_1;


		int iRap = rapBinMin;
  	    int iPt = ptBinMin;
  		sprintf(TreeBinID_,"%s_rap%d_pT%d",TreeID,iRap,iPt);TreeBinID=TreeBinID_;

  	    double ptlow=onia::pTRange[iRap][iPt-1];
  	    double pthigh=onia::pTRange[iRap][iPt];
  	    double raplow=onia::rapForPTRange[iRap-1];
  	    double raphigh=onia::rapForPTRange[iRap];

 	    sprintf(rapptstruct,"%s/rap%d_pT%d",substruct,iRap,iPt); if(!RealData) gSystem->mkdir(rapptstruct);

 	    if(gen) {sprintf(filenameFrom,"%s/polGen.C",substruct);					sprintf(filenameTo,"%s/polGen.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	    if(rec) {sprintf(filenameFrom,"%s/polRec.C",substruct);					sprintf(filenameTo,"%s/polRec.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	    if(fit) {sprintf(filenameFrom,"%s/polFit.C",substruct);					sprintf(filenameTo,"%s/polFit.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	    if(plot) {sprintf(filenameFrom,"%s/polPlot.C",substruct);				sprintf(filenameTo,"%s/polPlot.C",rapptstruct); 						gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}

 	    sprintf(filenameFrom,"%s/effsAndCuts.h",substruct);					sprintf(filenameTo,"%s/effsAndCuts.h",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);

 	    gSystem->cd(rapptstruct);

 	    if(gen) gROOT->ProcessLine(".L polGen.C+");
 	    if(rec) gROOT->ProcessLine(".L polRec.C+");
 	    if(fit) gROOT->ProcessLine(".L polFit.C+");
 	    if(plot) gROOT->ProcessLine(".L polPlot.C+");


/// Extract number of signal and background events to be generated, as well as f_BG to be generated to result in desired effective f_BG:

 	    int nTargetEvents;
 	    nTargetEvents = ToyMC::numEvents[iRap-1][iPt-1];

 	    if(gen){

			int numEvCheck = 500000;
			f_BG = ToyMC::fracBackground[iRap-1][iPt-1];
			sprintf(tmpfilename,"%s/data.root",rapptstruct);
			TFile* dataFile = new TFile(tmpfilename, "READ");

			if(dataFile->Get("isBGdistribution")==NULL){
				OutputDirectory=rapptstruct;
				polGen(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,numEvCheck,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,-999,OutputDirectory);
				if(rec)polRec(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,nEff,FidCuts,OutputDirectory, true);
				sprintf(tmpfilename,"%s/genData.root",rapptstruct);			gSystem->Unlink(tmpfilename);
				sprintf(tmpfilename,"%s/GenResults.root",rapptstruct);		gSystem->Unlink(tmpfilename);
			}

			sprintf(tmpfilename,"%s/data.root",rapptstruct);
			dataFile = new TFile(tmpfilename, "READ");
			TH1D* isBG_distribution = (TH1D*)dataFile->Get("isBGdistribution");

			double sigFact = isBG_distribution->GetBinContent(1)/(numEvCheck*(1-f_BG));
			double bkgFact = isBG_distribution->GetBinContent(2)/(numEvCheck*f_BG);

			dataFile->Close();


			if(ConstEvents_) nTargetEvents = ConstEvents;

			n_events = nTargetEvents/sigFact+(nTargetEvents/(1-f_BG)-nTargetEvents)/bkgFact;

			f_BG = (n_events-nTargetEvents/sigFact)/n_events;

 	    }
/// Start actual Generation and Fits:

		int iGen = ThisGen;

  		seconds = time (NULL); time_1=seconds;
  		sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen); if(!RealData) gSystem->mkdir(dirstruct);
  		OutputDirectory=dirstruct;
  		if(RealData) OutputDirectory=basestruct;

		if(gen)polGen(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,n_events,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,iGen,OutputDirectory);
		if(rec)polRec(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,nEff,FidCuts,OutputDirectory, false);
  		if(fit)polFit(nSample,FidCuts, nEff, OutputDirectory, realdatadir, TreeBinID, RealData);
  		if(plot)polPlot(OutputDirectory, TreeBinID, RealData);

	 	sprintf(tmpfilename,"%s/genData.root",dirstruct);			gSystem->Unlink(tmpfilename);
  		sprintf(tmpfilename,"%s/data.root",dirstruct);				gSystem->Unlink(tmpfilename);
  		sprintf(tmpfilename,"%s/efficiency.root",dirstruct);		gSystem->Unlink(tmpfilename);

  		seconds = time (NULL);

  		if(fit) cout<<"Proccessing time for this generation: "<<seconds-time_1<<" s"<<" (Corresponding to "<<(seconds-time_1)/60<<" min)"<<endl;
  		if(fit) cout<<"Per signal event in final sample:     "<<(seconds-time_1)/nTargetEvents*1000<<" ms"<<endl;

 	    if(gen)  {sprintf(tmpfilename,"%s/polGen.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	    if(rec)  {sprintf(tmpfilename,"%s/polRec.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	    if(fit)  {sprintf(tmpfilename,"%s/polFit.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	    if(plot)  {sprintf(tmpfilename,"%s/polPlot.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	    if(fit && RealData)  {sprintf(tmpfilename,"%s/polFit_C.d",basestruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",basestruct);			gSystem->Unlink(tmpfilename);}
 	    if(plot && RealData)  {sprintf(tmpfilename,"%s/polPlot_C.d",basestruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",basestruct);			gSystem->Unlink(tmpfilename);}

 	    sprintf(tmpfilename,"%s/effsAndCuts.h",rapptstruct);			gSystem->Unlink(tmpfilename);

return 0;
  	}


