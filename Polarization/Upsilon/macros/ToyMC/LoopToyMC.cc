#include "../../interface/rootIncludes.inc"
#include "../../interface/commonVar.h"
#include "ToyMC.h"

#include "TSystem.h"
#include "TROOT.h"
#include "polGen.C"
#include "polRec.C"
#include "polFit.C"

#include <time.h>

//====================================

int main(int argc, char** argv) {


	int nGenerations=1;
	int polScenSig=1;
	int frameSig=1;
	int polScenBkg=1;
	int frameBkg=1;
	int rapBinMin=1;
	int rapBinMax=1;
	int ptBinMin=1;
	int ptBinMax=1;
	int nEff=1;
	int FidCuts=1;
	int nSample=1;
	int ConstEvents=1;
	int nSkipGen=1;
	int ThisGen=1;

	bool ConstEvents_(false);
	bool gen(false);
	bool rec(false);
	bool fit(false);
	bool plot(false);

	Char_t *storagedir = "/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC"; //Storage Directory
	Char_t *basedir = "/Users/valentinknuenz/usr/local/workspace/Upsilon"; //Code Directory
  	Char_t *JobID = "Default";

	  for( int i=0;i < argc; ++i ) {
	    if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
	    if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
	    if(std::string(argv[i]).find("rapBinMin") != std::string::npos) {char* rapBinMinchar = argv[i]; char* rapBinMinchar2 = strtok (rapBinMinchar, "p"); rapBinMin = atof(rapBinMinchar2); cout<<"rapBinMin = "<<rapBinMin<<endl;}
	    if(std::string(argv[i]).find("rapBinMax") != std::string::npos) {char* rapBinMaxchar = argv[i]; char* rapBinMaxchar2 = strtok (rapBinMaxchar, "p"); rapBinMax = atof(rapBinMaxchar2); cout<<"rapBinMax = "<<rapBinMax<<endl;}
	    if(std::string(argv[i]).find("nGenerations") != std::string::npos) {char* nGenerationschar = argv[i]; char* nGenerationschar2 = strtok (nGenerationschar, "p"); nGenerations = atof(nGenerationschar2); cout<<"nGenerations = "<<nGenerations<<endl;}
	    if(std::string(argv[i]).find("Sigframe") != std::string::npos) {char* framecharSig = argv[i]; char* framecharSig2 = strtok (framecharSig, "p"); frameSig = atof(framecharSig2); cout<<"frameSig = "<<frameSig<<endl;}
	    if(std::string(argv[i]).find("polScenSig") != std::string::npos) {char* polScencharSig = argv[i]; char* polScencharSig2 = strtok (polScencharSig, "p"); polScenSig = atof(polScencharSig2); cout<<"polScenSig = "<<polScenSig<<endl;}
	    if(std::string(argv[i]).find("Bkgframe") != std::string::npos) {char* framecharBkg = argv[i]; char* framecharBkg2 = strtok (framecharBkg, "p"); frameBkg = atof(framecharBkg2); cout<<"frameBkg = "<<frameBkg<<endl;}
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
	  }

  		double lambda_theta_sig_ = ToyMC::ScenarioSig[0][polScenSig-1];
  		double lambda_phi_sig_ = ToyMC::ScenarioSig[1][polScenSig-1];
  		double lambda_thetaphi_sig_ = ToyMC::ScenarioSig[2][polScenSig-1];

  		double lambda_theta_bkg_ = ToyMC::ScenarioBkg[0][polScenBkg-1];
  		double lambda_phi_bkg_ = ToyMC::ScenarioBkg[1][polScenBkg-1];
  		double lambda_thetaphi_bkg_ = ToyMC::ScenarioBkg[2][polScenBkg-1];


  	  	Char_t *ToyDirectory;
  	  	double f_BG;
  	 	int n_events;
  		char basestruct[200],substruct[200], dirstruct[200], rapptstruct[200], filenameFrom[200], filenameTo[200] , tmpfilename[200];

  		sprintf(basestruct,"%s/%s",storagedir,JobID);gSystem->mkdir(basestruct);

  		sprintf(substruct,"%s/Sig_frame%dscen%d_Bkg_frame%dscen%d",basestruct,frameSig,polScenSig,frameBkg,polScenBkg);gSystem->mkdir(substruct);

  		time_t seconds; seconds = time (NULL); double time_0=seconds; double time_1;


		int iRap = rapBinMin;
  	    int iPt = ptBinMin;

  	    double ptlow=onia::pTRange[iRap][iPt-1];
  	    double pthigh=onia::pTRange[iRap][iPt];
  	    double raplow=onia::rapForPTRange[iRap-1];
  	    double raphigh=onia::rapForPTRange[iRap];

 	    sprintf(rapptstruct,"%s/rap%d_pT%d",substruct,iRap,iPt);gSystem->mkdir(rapptstruct);

 	    if(gen) {sprintf(filenameFrom,"%s/polGen.C",substruct);					sprintf(filenameTo,"%s/polGen.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	    if(rec) {sprintf(filenameFrom,"%s/polRec.C",substruct);					sprintf(filenameTo,"%s/polRec.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	    if(fit) {sprintf(filenameFrom,"%s/polFit.C",substruct);					sprintf(filenameTo,"%s/polFit.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}

 	    gSystem->cd(rapptstruct);

 	    if(gen) gROOT->ProcessLine(".L polGen.C+");
 	    if(rec) gROOT->ProcessLine(".L polRec.C+");
 	    if(fit) gROOT->ProcessLine(".L polFit.C+");


/// Extract number of signal and background events to be generated, as well as f_BG to be generated to result in desired effective f_BG:

 	    int numEvCheck = 500000;
		f_BG = ToyMC::fracBackground[iRap-1][iPt-1];
 	    sprintf(tmpfilename,"%s/data.root",rapptstruct);
		TFile* dataFile = new TFile(tmpfilename, "READ");

		if(dataFile->Get("isBG_distribution")==NULL){
			ToyDirectory=rapptstruct;
			if(gen)polGen(raplow,raphigh,ptlow,pthigh,numEvCheck,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,-999,ToyDirectory);
			if(rec)polRec(nEff,FidCuts,ToyDirectory);
  			sprintf(tmpfilename,"%s/genData.root",rapptstruct);			gSystem->Unlink(tmpfilename);
  			sprintf(tmpfilename,"%s/efficiency.root",rapptstruct);		gSystem->Unlink(tmpfilename);
  			sprintf(tmpfilename,"%s/GenResults.root",rapptstruct);		gSystem->Unlink(tmpfilename);
		}

 	    sprintf(tmpfilename,"%s/data.root",rapptstruct);
		dataFile = new TFile(tmpfilename, "READ");
		TH1D* isBG_distribution = (TH1D*)dataFile->Get("isBG_distribution");

		double sigFact = isBG_distribution->GetBinContent(1)/(numEvCheck*(1-f_BG));
		double bkgFact = isBG_distribution->GetBinContent(2)/(numEvCheck*f_BG);

		dataFile->Close();

		int nTargetEvents = ToyMC::numEvents[iRap-1][iPt-1];
		if(ConstEvents_) nTargetEvents = ConstEvents;

		n_events = nTargetEvents/sigFact+(nTargetEvents/(1-f_BG)-nTargetEvents)/bkgFact;

		f_BG = (n_events-nTargetEvents/sigFact)/n_events;


/// Start actual Generation and Fits:

		int iGen = ThisGen;

  		seconds = time (NULL); time_1=seconds;
  		sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen);gSystem->mkdir(dirstruct);
  		ToyDirectory=dirstruct;

		if(gen)polGen(raplow,raphigh,ptlow,pthigh,n_events,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,iGen,ToyDirectory);
		if(rec)polRec(nEff,FidCuts,ToyDirectory);
  		if(fit)polFit(nSample,ToyDirectory);

	 	sprintf(tmpfilename,"%s/genData.root",dirstruct);			gSystem->Unlink(tmpfilename);
  		sprintf(tmpfilename,"%s/data.root",dirstruct);				gSystem->Unlink(tmpfilename);
  		sprintf(tmpfilename,"%s/efficiency.root",dirstruct);		gSystem->Unlink(tmpfilename);

  		seconds = time (NULL);

  		if(fit) cout<<"Proccessing time for this generation: "<<seconds-time_1<<" s"<<" (Corresponding to "<<(seconds-time_1)/60<<" min)"<<endl;
  		if(fit) cout<<"Per signal event in final sample:     "<<(seconds-time_1)/nTargetEvents*1000<<" ms"<<endl;

 	    if(gen)  {sprintf(tmpfilename,"%s/polGen.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	    if(rec)  {sprintf(tmpfilename,"%s/polRec.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	    if(fit)  {sprintf(tmpfilename,"%s/polFit.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}



return 0;
  	}


