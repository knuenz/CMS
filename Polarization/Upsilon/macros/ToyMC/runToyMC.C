//#include "/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/Upsilon/interface/rootIncludes.inc"
//#include "/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/Upsilon/interface/commonVar.h"
//#include "/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/Upsilon/macros/ToyMC/ToyMC.h"
#include "/Users/valentinknuenz/usr/local/workspace/Upsilon/interface/rootIncludes.inc"
#include "/Users/valentinknuenz/usr/local/workspace/Upsilon/interface/commonVar.h"
#include "/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC/ToyMC.h"


#include <time.h>

//====================================

void runToyMC(){

	bool gen(true);
	bool recfit(true);
	bool plot(false);

// Define Output Directory Structure

//	char storagedir [200]= "/scratch/knuenz/Polarization/Upsilon/ToyMC"; //Storage Directory
//	char basedir [200]= "/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/Upsilon"; //Code Directory
	char storagedir [200]= "/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC"; //Storage Directory
	char basedir [200]= "/Users/valentinknuenz/usr/local/workspace/Upsilon"; //Code Directory
  	Char_t *JobID = "Test";

// Define Kinematic Bins to Test

  	int rapBinMin = 2;
    int rapBinMax = 2;
 	int ptBinMin = 8;
 	int ptBinMax = 8;

 // Define Polarization Scenario

 	int polScenSig=3;
  	int polScenBkg=3;
  	int frameSig=1;
  	int frameBkg=1;

// Define Basic Variables Needed For Toys

  	int nGenerations=1;
  	int nEff = 1;
  	int nSample = 12000;

//////////////////////////////////////////////////////////////////////////////////////

  	int nSkipGen=4;

  	Char_t *ToyDirectory;
  	double f_BG;
 	int n_events;
	char basestruct[200],substruct[200], dirstruct[200], rapptstruct[200], filenameFrom[200], filenameTo[200] , tmpfilename[200];

	sprintf(basestruct,"%s/%s",storagedir,JobID);gSystem->mkdir(basestruct);

	sprintf(filenameFrom,"%s/macros/ToyMC/polGen.C",basedir);					sprintf(filenameTo,"%s/polGen.C",basestruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);
	sprintf(filenameFrom,"%s/macros/ToyMC/polRec.C",basedir);					sprintf(filenameTo,"%s/polRec.C",basestruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);
	sprintf(filenameFrom,"%s/macros/ToyMC/polFit.C",basedir);					sprintf(filenameTo,"%s/polFit.C",basestruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);
	sprintf(filenameFrom,"%s/macros/ToyMC/polPlot.C",basedir);					sprintf(filenameTo,"%s/polPlot.C",basestruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);

	gSystem->cd(basestruct);

	if(gen){gROOT->ProcessLine(".L polGen.C+");
			if(recfit){	gROOT->ProcessLine(".L polRec.C+");
						gROOT->ProcessLine(".L polFit.C+");}
			}
	if(plot)gROOT->ProcessLine(".L polPlot.C+");

	sprintf(substruct,"%s/Sig_frame%dscen%d_Bkg_frame%dscen%d",basestruct,frameSig,polScenSig,frameBkg,polScenBkg);gSystem->mkdir(substruct);

	time_t seconds; seconds = time (NULL); double time_0=seconds; double time_1;


    for(int iRap = rapBinMin; iRap < rapBinMax+1; iRap++){
        for(int iPt = ptBinMin; iPt < ptBinMax+1; iPt++){

        	sprintf(rapptstruct,"%s/rap%d_pT%d",substruct,iRap,iPt);gSystem->mkdir(rapptstruct);

        	f_BG = ToyMC::fracBackground[iRap-1][iPt-1];
        	n_events = ToyMC::numEvents[iRap-1][iPt-1]/((1-f_BG)*ToyMC::EffCorrFrac[iRap-1][iPt-1]);


			for(int iGen = 1; iGen < nGenerations+1; iGen++){

				seconds = time (NULL); time_1=seconds;
				sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen);gSystem->mkdir(dirstruct);
				ToyDirectory=dirstruct;

				if(gen){polGen(iRap,iPt,n_events,f_BG,polScenSig,polScenBkg,frameSig,frameBkg,iGen+nSkipGen,ToyDirectory);
						if(recfit){	polRec(nEff,ToyDirectory);
									polFit(nSample,ToyDirectory);}
						}
				if(plot)polPlot(ToyDirectory);

				sprintf(tmpfilename,"%s/genData.root",dirstruct);		gSystem->Unlink(tmpfilename);
				sprintf(tmpfilename,"%s/data.root",dirstruct);		gSystem->Unlink(tmpfilename);
				sprintf(tmpfilename,"%s/efficiency.root",dirstruct);		gSystem->Unlink(tmpfilename);

				seconds = time (NULL);

				if(gen) cout<<"Proccessing time for this generation: "<<seconds-time_1<<" s"<<" (Corresponding to "<<(seconds-time_1)/60<<" min)"<<endl;
				if(gen) cout<<"Per signal event in final sample:     "<<(seconds-time_1)/ToyMC::numEvents[iRap-1][iPt-1]*1000<<" ms"<<endl;
			}
        }
    }

	  seconds = time (NULL);
	  cout<<"Total time needed for ToyMC test: "<<seconds-time_0<<" s"<<" (Corresponding to "<<(seconds-time_0)/3600<<" h)"<<endl;



}


