#include "../../interface/rootIncludes.inc"
#include "../../interface/commonVar.h"
#include "ToyMC.h"


#include <time.h>

//====================================

void runToyMC(){

	bool gen(true);
	bool rec(true);
	bool fit(true);
	bool plot(false);

// Define Output Directory Structure

	char storagedir [200]= "/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC"; //Storage Directory
	char basedir [200]= "/Users/valentinknuenz/usr/local/workspace/Upsilon"; //Code Directory
  	Char_t *JobID = "FidTry";

  	// Define Kinematic Bins to Test

  	  	int rapBinMin = 1;
  	    int rapBinMax = 1;
  	 	int ptBinMin = 1;
  	 	int ptBinMax = 1;

  	 // Define Polarization Scenario

  	 	int polScenSig=3;
  	  	int polScenBkg=3;
  	  	int frameSig=1;
  	  	int frameBkg=1;

  	// Define Basic Variables Needed For Toys

  	  	int nGenerations=1;
  	  	int nEff = 1;//1...Eff=const=1, 2...oldEff, 3...newEff
  	  	int FidCuts = 1;//0...no cuts, 1... loose cuts, 2... tight cuts
  	  	int nSample = 6000;//2000 burn-in iterations included

  	  	bool ConstEvents_(true);
  	  	int ConstEvents=20000;

  	//////////////////////////////////////////////////////////////////////////////////////
  	///CVS

  	  	int nSkipGen=0;

  	  	Char_t *ToyDirectory;
  	  	double f_BG;
  	 	int n_events;
  		char basestruct[200],substruct[200], dirstruct[200], rapptstruct[200], filenameFrom[200], filenameTo[200] , tmpfilename[200];

  		sprintf(basestruct,"%s/%s",storagedir,JobID);gSystem->mkdir(basestruct);

  		sprintf(substruct,"%s/Sig_frame%dscen%d_Bkg_frame%dscen%d",basestruct,frameSig,polScenSig,frameBkg,polScenBkg);gSystem->mkdir(substruct);

  		time_t seconds; seconds = time (NULL); double time_0=seconds; double time_1;


  	    for(int iRap = rapBinMin; iRap < rapBinMax+1; iRap++){
  	        for(int iPt = ptBinMin; iPt < ptBinMax+1; iPt++){

 	        	sprintf(rapptstruct,"%s/rap%d_pT%d",substruct,iRap,iPt);gSystem->mkdir(rapptstruct);

 	     		if(gen) {sprintf(filenameFrom,"%s/macros/ToyMC/polGen.C",basedir);					sprintf(filenameTo,"%s/polGen.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	     		if(rec) {sprintf(filenameFrom,"%s/macros/ToyMC/polRec.C",basedir);					sprintf(filenameTo,"%s/polRec.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	     		if(fit) {sprintf(filenameFrom,"%s/macros/ToyMC/polFit.C",basedir);					sprintf(filenameTo,"%s/polFit.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	     		if(plot){sprintf(filenameFrom,"%s/macros/ToyMC/polPlot.C",basedir);					sprintf(filenameTo,"%s/polPlot.C",rapptstruct); 						gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
 	     		sprintf(filenameFrom,"%s/interface/commonVar.h",basedir);							sprintf(filenameTo,"%s/commonVar.h",rapptstruct); 						gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);
 	     		sprintf(filenameFrom,"%s/macros/ToyMC/ToyMC.h",basedir);							sprintf(filenameTo,"%s/ToyMC.h",rapptstruct); 						gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);

 	     		gSystem->cd(rapptstruct);

 	     		if(gen) gROOT->ProcessLine(".L polGen.C+");
 	     		if(rec) gROOT->ProcessLine(".L polRec.C+");
 	     		if(fit) gROOT->ProcessLine(".L polFit.C+");
 	     		if(plot)gROOT->ProcessLine(".L polPlot.C+");

  	        	sprintf(tmpfilename,"%s/EffFraction.txt",rapptstruct); //gSystem->Unlink("AccEffFraction.txt");

  	        	/// Extract number of signal and background events to be generated, as well as f_BG to be generated to result in desired effective f_BG:
  	        	int numEvCheck = 500000;
				f_BG = ToyMC::fracBackground[iRap-1][iPt-1];
 	        	sprintf(tmpfilename,"%s/data.root",rapptstruct);
				TFile* dataFile = new TFile(tmpfilename, "READ");
				if(dataFile->Get("isBG_distribution")==NULL){
					ToyDirectory=rapptstruct;
					if(gen)polGen(iRap,iPt,numEvCheck,f_BG,polScenSig,polScenBkg,frameSig,frameBkg,-999,ToyDirectory);
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


				/// Start actual Toys:
				for(int iGen = 1; iGen < nGenerations+1; iGen++){

  					seconds = time (NULL); time_1=seconds;
  					sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen);gSystem->mkdir(dirstruct);
  					ToyDirectory=dirstruct;

  					if(gen)polGen(iRap,iPt,n_events,f_BG,polScenSig,polScenBkg,frameSig,frameBkg,iGen+nSkipGen,ToyDirectory);
  					if(rec)polRec(nEff,FidCuts,ToyDirectory);
  					if(fit)polFit(nSample,ToyDirectory);

  					if(plot)polPlot(ToyDirectory);

//	 					sprintf(tmpfilename,"%s/genData.root",dirstruct);			gSystem->Unlink(tmpfilename);
//  					sprintf(tmpfilename,"%s/data.root",dirstruct);				gSystem->Unlink(tmpfilename);
//  					sprintf(tmpfilename,"%s/efficiency.root",dirstruct);		gSystem->Unlink(tmpfilename);

  					seconds = time (NULL);

  					if(fit) cout<<"Proccessing time for this generation: "<<seconds-time_1<<" s"<<" (Corresponding to "<<(seconds-time_1)/60<<" min)"<<endl;
  					if(fit) cout<<"Per signal event in final sample:     "<<(seconds-time_1)/ToyMC::numEvents[iRap-1][iPt-1]*1000<<" ms"<<endl;
  				}// end iGen

 	     		if(gen)  {sprintf(tmpfilename,"%s/polGen.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	     		if(rec)  {sprintf(tmpfilename,"%s/polRec.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	     		if(fit)  {sprintf(tmpfilename,"%s/polFit.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	     		if(plot) {sprintf(tmpfilename,"%s/polPlot.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
 	     		sprintf(tmpfilename,"%s/commonVar.h",rapptstruct);					gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/ToyMC.h",rapptstruct);			gSystem->Unlink(tmpfilename);

  	        }// end iPt
  	    }// end iRap

  		  seconds = time (NULL);
  		  cout<<"Total time needed for ToyMC test: "<<seconds-time_0<<" s"<<" (Corresponding to "<<(seconds-time_0)/3600<<" h)"<<endl;


  	}


