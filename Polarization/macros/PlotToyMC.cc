#include <iostream>
#include <sstream>
#include <cstring>

#include "commonVar.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphErrors.h"

#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TLegend.h"

#include "RooWorkspace.h"

int main(int argc, char** argv) {
//	using namespace JPsiPolarization;
	using namespace std;
	using namespace RooFit;

	TFile* f= new TFile("/scratch/knuenz/Polarization/RootInput/ProjectClosure/scenario1/CLOSURE_scen1_rap2_pT8_gen1_CS-CS.root");
	RooWorkspace* space=(RooWorkspace*)f->Get("CLOSURE_scen1_rap2_pT8_gen1_CS");
	space->Print();

	char JobID[200];
	sprintf(JobID,"multi_try");
	char miniCrit[200];
	char miniCrit2[200];

	bool TwoCrit(false);
	bool logicOr(false);

	bool all(false);

	bool normFitted(false);
  	bool FileSearch(false);
  	bool noCut(false);

	int critCovVal = 3;

	double edmcondition = 10e100;

	double promptlambda_phi_injected_CS[3][9];
  	double promptlambda_theta_injected_CS[3][9];
  	double promptlambda_thetaphi_injected_CS[3][9];

  	double promptlambda_phi_injected_HX[3][9];
  	double promptlambda_theta_injected_HX[3][9];
  	double promptlambda_thetaphi_injected_HX[3][9];

  	double promptlambda_theta_injected = 0;

  	bool RealDataPRNPstats(false);
  	bool prompt(true);

  	bool saveedm(true);
  	bool serial(false);


//////////////// CLOSURE //////////////////////

  	bool Closure(false);
  	double scen=100;
  	bool tilde(false);
  	bool workspace(false);
  	char JOBNAME[200];
  	bool onlyOneComponent(false);
	bool JOBNAMEinj(false);

	  int generations = 50;

for( int i=0;i < argc; ++i ) {
	    if(std::string(argv[i]).find("--Criterium1=migrad") != std::string::npos) sprintf(miniCrit,"migrad");
	    if(std::string(argv[i]).find("--Criterium1=hesse") != std::string::npos) sprintf(miniCrit,"hesse");
	    if(std::string(argv[i]).find("--Criterium1=minos") != std::string::npos) sprintf(miniCrit,"minos");
	    if(std::string(argv[i]).find("--Criterium1=improve") != std::string::npos) sprintf(miniCrit,"improve");
	    if(std::string(argv[i]).find("--Criterium2=migrad") != std::string::npos) sprintf(miniCrit2,"migrad");
	 	if(std::string(argv[i]).find("--Criterium2=hesse") != std::string::npos) sprintf(miniCrit2,"hesse");
	 	if(std::string(argv[i]).find("--Criterium2=minos") != std::string::npos) sprintf(miniCrit2,"minos");
	    if(std::string(argv[i]).find("--Criterium2=improve") != std::string::npos) sprintf(miniCrit2,"improve");
	 	if(std::string(argv[i]).find("--TwoCrit") != std::string::npos) TwoCrit=true;
	 	if(std::string(argv[i]).find("--OneCrit") != std::string::npos) TwoCrit=false;
	 	if(std::string(argv[i]).find("--logicOr") != std::string::npos) logicOr=true;
	 	if(std::string(argv[i]).find("--logicAnd") != std::string::npos) logicOr=false;
	 	if(std::string(argv[i]).find("--all") != std::string::npos) all=true;
	 	if(std::string(argv[i]).find("--normFitted") != std::string::npos) normFitted=true;
	 	if(std::string(argv[i]).find("multi") != std::string::npos) sprintf(JobID,argv[i]);
	 	if(std::string(argv[i]).find("less") != std::string::npos) sprintf(JobID,argv[i]);
	 	if(std::string(argv[i]).find("--lambda_theta_plus05") != std::string::npos) {promptlambda_theta_injected = +0.5;}
	 	if(std::string(argv[i]).find("--lambda_theta_minus05") != std::string::npos) {promptlambda_theta_injected = -0.5;}
	 	if(std::string(argv[i]).find("--FileSearch") != std::string::npos) FileSearch=true;
	 	if(std::string(argv[i]).find("--noCut") != std::string::npos) noCut=true;
	 	if(std::string(argv[i]).find("--RealDataVal-PR") != std::string::npos) {RealDataPRNPstats=true;prompt=true;}
	 	if(std::string(argv[i]).find("--RealDataVal-NP") != std::string::npos) {RealDataPRNPstats=true;prompt=false;}
	 	if(std::string(argv[i]).find("--gen200") != std::string::npos) {generations = 200;}
	 	if(std::string(argv[i]).find("--gen1000") != std::string::npos) {generations = 1000;}
	 	if(std::string(argv[i]).find("--gen500") != std::string::npos) {generations = 500;}
	 	if(std::string(argv[i]).find("--RealDataVal-NP") != std::string::npos) {RealDataPRNPstats=true;prompt=false;}
	 	if(std::string(argv[i]).find("serial") != std::string::npos) {sprintf(JobID,argv[i]);serial=true;}
///////////////////// CLOSURE /////////////////////////
	 	if(std::string(argv[i]).find("--Closure") != std::string::npos) {Closure=true;}
	 	if(std::string(argv[i]).find("generations") != std::string::npos) {char* generationschar = argv[i]; char* generationschar2 = strtok (generationschar, "g"); generations = atof(generationschar2); cout<<"generations = "<<generations<<endl;}
	 	if(std::string(argv[i]).find("scen") != std::string::npos) {char* scenchar = argv[i]; char* scenchar2 = strtok (scenchar, "s"); scen = atof(scenchar2); cout<<"scen = "<<scen<<endl;}
	 	if(std::string(argv[i]).find("--PR") != std::string::npos) {prompt=true;}
	 	if(std::string(argv[i]).find("--NP") != std::string::npos) {prompt=false;}
	 	if(std::string(argv[i]).find("--tilde") != std::string::npos) {tilde=true;}
	 	if(std::string(argv[i]).find("--workspace") != std::string::npos) {workspace=true;}
	 	if(std::string(argv[i]).find("--onlyOne") != std::string::npos) {onlyOneComponent=true;}
		if(std::string(argv[i]).find("AccEffCut") != std::string::npos) {sprintf(JOBNAME,argv[i]);JOBNAMEinj=true;}

	  }


	int scenario=scen;

	if(Closure){
		if(!JOBNAMEinj) sprintf(JOBNAME,"UNSMOOTHnlb200_new");
		if(prompt)sprintf(JobID,"CLOSURE_%s_GEN%d_scen%d_PR",JOBNAME,generations,scenario);
		if(!prompt)sprintf(JobID,"CLOSURE_%s_GEN%d_scen%d_NP",JOBNAME,generations,scenario);

	}

	if(!TwoCrit) sprintf(miniCrit2,miniCrit);

	cout<<"Criterium 1 = "<<miniCrit<<endl;
	cout<<"Criterium 2 = "<<miniCrit2<<endl;
	cout<<"TwoCrit? "<<TwoCrit<<endl;
	cout<<"logicOr? "<<logicOr<<endl;
	cout<<"all? "<<all<<endl;
	cout<<"normFitted? "<<normFitted<<endl;
	cout<<"JobID = "<<JobID<<endl;
//	cout<<"promptlambda_theta_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_theta_injected_CS[yBinstart][ptBinstart]<<endl;
//	cout<<"promptlambda_theta_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_theta_injected_HX[yBinstart][ptBinstart]<<endl;
	cout<<"FileSearch? "<<FileSearch<<endl;

////////////	if TwoCrit = true:
//////////// 	if logicOr = true, one of both Results must fullfill covQual() == critCovQual
//////////// 	if logicOr = false, both Results must fullfill covQual() == critCovQual
////////////	Parameter results: if TwoCrit = true, second Result is valid

	if(all) TwoCrit=false;

	char miniCritLogic[200];
	if(logicOr) sprintf(miniCritLogic,"|");
	if(!logicOr) sprintf(miniCritLogic,"&");

	char dirStruct[200];
	sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC/%s",JobID);

	gSystem->mkdir("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC");
	gSystem->mkdir(dirStruct);

	if(TwoCrit) sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC/%s/%s%s%s",JobID,miniCrit,miniCritLogic,miniCrit2);
	else sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC/%s/%s",JobID,miniCrit);
	if(all) sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC/%s/%s_all",JobID,miniCrit);
	if(noCut) sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC/%s/%s_noCut",JobID,miniCrit);
	if(TwoCrit && noCut) sprintf(dirStruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Plots/ToyMC/%s/%s%s%s_noCut",JobID,miniCrit,miniCritLogic,miniCrit2);

	gSystem->mkdir(dirStruct);

	cout<<dirStruct<<endl;

	  gSystem->mkdir("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Results/ToyMC");

	  gStyle->SetTitleFillColor(kWhite);


///////// VARIABLE PLOTTING ADJUSTMENT AND OUTLIER CONFIGURATION ///////////////

	    double borders=10;//cut
	    double borders_val=2;//cut

	  	double cut;
	  	if(noCut) cut = 5000000;
	  	else cut = 10;
	  	double cutval=10;

/////////////////////////////////////////////////////////////////////////////////

	  char line[200];

	  char outputfilename2[200];
	  sprintf(outputfilename2,"%s/parametersToyMCTotal.txt",dirStruct);
	  FILE *outputFile2 = fopen(outputfilename2,"w");

	  int CONV[2][6][13][100]={NULL};

	  double mean_phi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errmean_phi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double mean_theta[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errmean_theta[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double mean_thetaphi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errmean_thetaphi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double sigma_phi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errsigma_phi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double sigma_theta[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errsigma_theta[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double sigma_thetaphi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errsigma_thetaphi[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};


	  double mean_phi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errmean_phi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double mean_theta_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errmean_theta_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double mean_thetaphi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errmean_thetaphi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double sigma_phi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errsigma_phi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double sigma_theta_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errsigma_theta_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	  double sigma_thetaphi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  double errsigma_thetaphi_val[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};


	  	for(int frameDex = 0; frameDex < 2; frameDex++) {
  	for(int yBinstart = 1; yBinstart < 3; yBinstart++) {
  		for(int ptBinstart = 1; ptBinstart < 9; ptBinstart++) {

	  mean_phi[frameDex][yBinstart][ptBinstart]=-9999;
	  sigma_phi[frameDex][yBinstart][ptBinstart]=-9999;
	  mean_theta[frameDex][yBinstart][ptBinstart]=-9999;
	  sigma_theta[frameDex][yBinstart][ptBinstart]=-9999;
	  mean_thetaphi[frameDex][yBinstart][ptBinstart]=-9999;
	  sigma_thetaphi[frameDex][yBinstart][ptBinstart]=-9999;
	  mean_phi_val[frameDex][yBinstart][ptBinstart]=-9999;
	  sigma_phi_val[frameDex][yBinstart][ptBinstart]=-9999;
	  mean_theta_val[frameDex][yBinstart][ptBinstart]=-9999;
	  sigma_theta_val[frameDex][yBinstart][ptBinstart]=-9999;
	  mean_thetaphi_val[frameDex][yBinstart][ptBinstart]=-9999;
	  sigma_thetaphi_val[frameDex][yBinstart][ptBinstart]=-9999;

  		}}}

	  double overflow[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};
	  int MissingGes[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	    int convCountCheck[2][6][13]={NULL};

	    int convCount[2][6][13]={NULL};

	    bool nofile(false);

		  int missingGenerations=0;

	    	for(int yBinstart = 1; yBinstart < 2; yBinstart++) {
	    		for(int ptBinstart = 6; ptBinstart < 9; ptBinstart++) {
//3,9 for mid, intermediate bins

	    	if (yBinstart==1 && ptBinstart==8) continue;
	    //	if (yBinstart==2 && ptBinstart==7) continue;
	    //	if (yBinstart==2 && ptBinstart==8) continue;

	    			cout<<"rap"<<yBinstart<<"_pt"<<ptBinstart<<endl;

	    		//	yBinstart=3;
	    		//	ptBinstart=9;

	    			nofile=false;


		char inputfilename[200];
		sprintf(inputfilename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/CRAB/%s/RooFitResult_rap%d_pt%d.root",JobID,yBinstart,ptBinstart);
		if(serial) sprintf(inputfilename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/SERIAL/%s_res/RooFitResult_rap%d_pt%d.root",JobID,yBinstart,ptBinstart);
		//if(Closure && !workspace) sprintf(inputfilename,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/scenario%d/RooFitResult_scen%d_rap%d_pT%d.root",scenario,scenario,yBinstart,ptBinstart);
		if(Closure && !workspace) sprintf(inputfilename,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/PlotResults/%s/RooFitResult_%s_scen%d_rap%d_pt%d.root",JOBNAME,JOBNAME,scenario,yBinstart,ptBinstart);

		TFile *resfile = new TFile(inputfilename,"UPDATE");




		char outputfilename[200];
		sprintf(outputfilename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/Results/ToyMC/parametersToyMC_rap%d_pt%d.txt",yBinstart,ptBinstart);
		FILE *outputFile = fopen(outputfilename,"w");

		  int ptBin=0;
		  int yBin;

		  bool RealMC(false);


		  int CONV[2][6][13][101]={NULL};
		  int generation=0;
		  char zeile[200];
		  int i=0;
		  int fitresIndex=1000;
		  int MinosresIndex=1000;
		  double edmCondition=0.01;

		  double edm;
		  double nPrompt[2][6][13][101];
		  double errnPrompt[2][6][13][101];
		  double promptlambda_phi[2][6][13][101];
		  double errpromptlambda_phi[2][6][13][101];
		  double promptlambda_theta[2][6][13][101];
		  double errpromptlambda_theta[2][6][13][101];
		  double promptlambda_thetaphi[2][6][13][101];
		  double errpromptlambda_thetaphi[2][6][13][101];
		  double edm1[2][6][13][101];
		  double edm2[2][6][13][101];
		  double status1[2][6][13][101];
		  double status2[2][6][13][101];
		  double covQualy1[2][6][13][101];
		  double covQualy2[2][6][13][101];


		 	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = 0;
			promptlambda_theta_injected_CS[yBinstart][ptBinstart] = promptlambda_theta_injected;
			promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = 0;
			promptlambda_phi_injected_HX[yBinstart][ptBinstart] = 0;
			promptlambda_theta_injected_HX[yBinstart][ptBinstart] = promptlambda_theta_injected;
			promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = 0;


		  if(RealDataPRNPstats){
		  if(prompt){
		  	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = jpsi::inj_PR_CS_ph[yBinstart-1][ptBinstart-1];
		    promptlambda_theta_injected_CS[yBinstart][ptBinstart] = jpsi::inj_PR_CS_th[yBinstart-1][ptBinstart-1];
		    promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = jpsi::inj_PR_CS_thph[yBinstart-1][ptBinstart-1];
		    promptlambda_phi_injected_HX[yBinstart][ptBinstart] = jpsi::inj_PR_HX_ph[yBinstart-1][ptBinstart-1];
		    promptlambda_theta_injected_HX[yBinstart][ptBinstart] = jpsi::inj_PR_HX_th[yBinstart-1][ptBinstart-1];
		    promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = jpsi::inj_PR_HX_thph[yBinstart-1][ptBinstart-1];
		  }
		  if(!prompt){
		  	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = jpsi::inj_NP_CS_ph[yBinstart-1][ptBinstart-1];
		    promptlambda_theta_injected_CS[yBinstart][ptBinstart] = jpsi::inj_NP_CS_th[yBinstart-1][ptBinstart-1];
		    promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = jpsi::inj_NP_CS_thph[yBinstart-1][ptBinstart-1];
		    promptlambda_phi_injected_HX[yBinstart][ptBinstart] = jpsi::inj_NP_HX_ph[yBinstart-1][ptBinstart-1];
		    promptlambda_theta_injected_HX[yBinstart][ptBinstart] = jpsi::inj_NP_HX_th[yBinstart-1][ptBinstart-1];
		    promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = jpsi::inj_NP_HX_thph[yBinstart-1][ptBinstart-1];
		  }
		  }

//////// Start rolling over the file, row by row //////////////////




///////////////////////// END OF FILE READING ///////////////////////////////////////////////////////////

if(Closure){

	  if(prompt){
	  	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][1];
	    promptlambda_theta_injected_CS[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][0];
	    promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][2];
	    promptlambda_phi_injected_HX[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][1];
	    promptlambda_theta_injected_HX[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][0];
	    promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][2];
	    if(tilde){
		  	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = (jpsi::scenarioLambdas[scenario-1][0]+3*jpsi::scenarioLambdas[scenario-1][1])/(1-jpsi::scenarioLambdas[scenario-1][1]);
		  	promptlambda_phi_injected_HX[yBinstart][ptBinstart] = (jpsi::scenarioLambdas[scenario-1][0]+3*jpsi::scenarioLambdas[scenario-1][1])/(1-jpsi::scenarioLambdas[scenario-1][1]);
	    	  }
	  }
	  if(!prompt){
	  	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][4];
	    promptlambda_theta_injected_CS[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][3];
	    promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][5];
	    promptlambda_phi_injected_HX[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][4];
	    promptlambda_theta_injected_HX[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][3];
	    promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = jpsi::scenarioLambdas[scenario-1][5];
	    if(tilde){
	    	promptlambda_phi_injected_CS[yBinstart][ptBinstart] = (jpsi::scenarioLambdas[scenario-1][3]+3*jpsi::scenarioLambdas[scenario-1][4])/(1-jpsi::scenarioLambdas[scenario-1][4]);
	    	promptlambda_phi_injected_HX[yBinstart][ptBinstart] = (jpsi::scenarioLambdas[scenario-1][3]+3*jpsi::scenarioLambdas[scenario-1][4])/(1-jpsi::scenarioLambdas[scenario-1][4]);
	    	  }
	  }


	  cout<<"promptlambda_theta_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_theta_injected_CS[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_phi_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_phi_injected_CS[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_theta_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_theta_injected_HX[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_phi_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_phi_injected_HX[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart]<<endl;




	    for(int frameDex = 0; frameDex<2; ++frameDex){ missingGenerations=0;
				for(int generation = 1; generation < generations+1; generation++) {
					cout<<"GEN"<<generation<<" (yBin= "<<yBinstart<<" , ptBin= "<<ptBinstart<<" frame= "<<frameDex<<")"<<endl;

					  RooFitResult* roofitres_Closure_migrad;
					  RooFitResult* roofitres_Closure_hesse;

					  bool onlyCrit1converged(false);
					  char resname[200];
					  char res2name[200];

					  char framechar[200];
					  if(frameDex==0) sprintf(framechar,"CS");
					  if(frameDex==1) sprintf(framechar,"HX");

					  sprintf(resname,"migrad_pol_fitresult_rap%d_pt%d_%s_gen%d",yBinstart,ptBinstart,framechar,generation);
					  sprintf(res2name,"hesse_pol_fitresult_rap%d_pt%d_%s_gen%d",yBinstart,ptBinstart,framechar,generation);

					  char Closure_inputfilename[200];
					  if(frameDex==0)sprintf(Closure_inputfilename,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/scenario%d/CLOSURE_scen%d_rap%d_pT%d_gen%d_CS-CS.root",scenario,scenario,yBinstart,ptBinstart,generation);
					  if(frameDex==1)sprintf(Closure_inputfilename,"/scratch/knuenz/Polarization/RootInput/ProjectClosure/scenario%d/CLOSURE_scen%d_rap%d_pT%d_gen%d_HX-HX.root",scenario,scenario,yBinstart,ptBinstart,generation);

					  TFile *Closure_resfile;
					  RooWorkspace* wss;
					  if(workspace){
						  Closure_resfile = new TFile(Closure_inputfilename,"UPDATE");

						  char workspacename[200];
						  if(frameDex==0)sprintf(workspacename,"CLOSURE_scen%d_rap%d_pT%d_gen%d_CS",scenario,yBinstart,ptBinstart,generation);
						  if(frameDex==1)sprintf(workspacename,"CLOSURE_scen%d_rap%d_pT%d_gen%d_HX",scenario,yBinstart,ptBinstart,generation);

						  wss=(RooWorkspace*)Closure_resfile->Get(workspacename);
					  }

					  if(resfile->Get(resname) == NULL) {
					  missingGenerations++;
					  nofile=true;
					  nPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  errnPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
			          promptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  promptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  edm1[frameDex][yBinstart][ptBinstart][generation]=-9999;
					  edm2[frameDex][yBinstart][ptBinstart][generation]=-9999;

/*					  mean_phi[frameDex][yBinstart][ptBinstart]=-9999;
					  sigma_phi[frameDex][yBinstart][ptBinstart]=-9999;
					  mean_theta[frameDex][yBinstart][ptBinstart]=-9999;
					  sigma_theta[frameDex][yBinstart][ptBinstart]=-9999;
					  mean_thetaphi[frameDex][yBinstart][ptBinstart]=-9999;
					  sigma_thetaphi[frameDex][yBinstart][ptBinstart]=-9999;
					  mean_phi_val[frameDex][yBinstart][ptBinstart]=-9999;
					  sigma_phi_val[frameDex][yBinstart][ptBinstart]=-9999;
					  mean_theta_val[frameDex][yBinstart][ptBinstart]=-9999;
					  sigma_theta_val[frameDex][yBinstart][ptBinstart]=-9999;
					  mean_thetaphi_val[frameDex][yBinstart][ptBinstart]=-9999;
					  sigma_thetaphi_val[frameDex][yBinstart][ptBinstart]=-9999;
*/
					  }
					  else {

						 if(workspace){
					  roofitres_Closure_migrad = (RooFitResult*)wss->genobj(resname);
					  roofitres_Closure_hesse = (RooFitResult*)wss->genobj(res2name);
						 }
						 else{
					  roofitres_Closure_migrad = (RooFitResult*)resfile->Get(resname);
					  roofitres_Closure_hesse = (RooFitResult*)resfile->Get(res2name);

					  roofitres_Closure_hesse->Print();
						 }

					  cout<<"status  migrad: "<<roofitres_Closure_migrad->status()<<endl;
					  cout<<"covQual migrad: "<<roofitres_Closure_migrad->covQual()<<endl;
					  cout<<"status  minos: "<<roofitres_Closure_hesse->status()<<endl;
					  cout<<"covQual minos: "<<roofitres_Closure_hesse->covQual()<<endl;


					  RooArgList params = roofitres_Closure_hesse->floatParsFinal();
					//  if(!TwoCrit || all) {/*cout<<"Values from 1st"<<endl;*/ params = roofitres->floatParsFinal();}

					  edm1[frameDex][yBinstart][ptBinstart][generation]=roofitres_Closure_migrad->edm();
					  edm2[frameDex][yBinstart][ptBinstart][generation]=roofitres_Closure_hesse->edm();
					  status1[frameDex][yBinstart][ptBinstart][generation]=roofitres_Closure_migrad->status();
					  status2[frameDex][yBinstart][ptBinstart][generation]=roofitres_Closure_hesse->status();
					  covQualy1[frameDex][yBinstart][ptBinstart][generation]=roofitres_Closure_migrad->covQual();
					  covQualy2[frameDex][yBinstart][ptBinstart][generation]=roofitres_Closure_hesse->covQual();

					  bool criterium(false);
					//  if(TwoCrit && logicOr && roofitres_Closure_migrad->covQual() == critCovVal && roofitres_Closure_migrad->edm() < edmcondition || TwoCrit && logicOr && roofitres_Closure_minos->covQual() == critCovVal && roofitres_Closure_minos->edm() < edmcondition) {criterium=true; if(roofitres_Closure_minos->covQual() != critCovVal) onlyCrit1converged=true;;};
				//	  if(TwoCrit && !logicOr && roofitres->covQual() == critCovVal && roofitres->edm() < edmcondition && roofitres2->covQual() == critCovVal && roofitres2->edm() < edmcondition ) {criterium=true;}
				//	  if (!TwoCrit && roofitres->edm() < edmcondition && roofitres->covQual() == critCovVal/* || roofitres->status() == 0*/) {criterium=true;}
					//  if(all) criterium=true;

					//  if(onlyCrit1converged) params = roofitres_Closure_migrad->floatParsFinal();

					  if(TwoCrit && logicOr && roofitres_Closure_migrad->covQual() == critCovVal || TwoCrit && logicOr && roofitres_Closure_hesse->covQual() == critCovVal) {criterium=true;cout<<"Criterium satisfied"<<endl;}

					  RooRealVar* lambdaphi_;
					  RooRealVar* lambdatheta_;
					  RooRealVar* lambdathetaphi_;

if(!onlyOneComponent){
					  if(prompt){
					  lambdaphi_ =  (RooRealVar*)params.at(5);
					  lambdatheta_ =  (RooRealVar*)params.at(1);
					  lambdathetaphi_ =  (RooRealVar*)params.at(3);
					  }
					  if(!prompt){
					  lambdaphi_ =  (RooRealVar*)params.at(4);
					  lambdatheta_ =  (RooRealVar*)params.at(0);
					  lambdathetaphi_ =  (RooRealVar*)params.at(2);
					  }
}

if(onlyOneComponent){

					  lambdaphi_ =  (RooRealVar*)params.at(2);
					  lambdatheta_ =  (RooRealVar*)params.at(0);
					  lambdathetaphi_ =  (RooRealVar*)params.at(1);
}

					  if(criterium){
						//  cout<<"status  "<<roofitres->status()<<endl;
						//  cout<<"covqual "<<roofitres->covQual()<<endl;

						  convCount[frameDex][yBinstart][ptBinstart]++;
						  CONV[frameDex][yBinstart][ptBinstart][generation]=1;

			/*
						  nPrompt[frameDex][yBinstart][ptBinstart][generation]=nPrompt_->getVal();
						  if(nPrompt_->getAsymErrorHi()>fabs(nPrompt_->getAsymErrorLo())) errnPrompt[frameDex][yBinstart][ptBinstart][generation]=nPrompt_->getAsymErrorHi();
						  else errnPrompt[frameDex][yBinstart][ptBinstart][generation]=fabs(nPrompt_->getAsymErrorLo());
						  if(errnPrompt[frameDex][yBinstart][ptBinstart][generation]==0) errnPrompt[frameDex][yBinstart][ptBinstart][generation]=nPrompt_->getError();
			*/
						  promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=lambdaphi_->getVal();
						  if(lambdaphi_->getAsymErrorHi()>fabs(lambdaphi_->getAsymErrorLo())) errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=lambdaphi_->getAsymErrorHi();
						  else errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=fabs(lambdaphi_->getAsymErrorLo());
						  if(errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]==0) errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=lambdaphi_->getError();

						  promptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=lambdatheta_->getVal();
						  if(lambdatheta_->getAsymErrorHi()>fabs(lambdatheta_->getAsymErrorLo())) errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=lambdatheta_->getAsymErrorHi();
						  else errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=fabs(lambdatheta_->getAsymErrorLo());
						  if(errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]==0) errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=lambdatheta_->getError();

						  promptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=lambdathetaphi_->getVal();
						  if(lambdathetaphi_->getAsymErrorHi()>fabs(lambdathetaphi_->getAsymErrorLo())) errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=lambdathetaphi_->getAsymErrorHi();
						  else errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=fabs(lambdathetaphi_->getAsymErrorLo());
						  if(errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]==0) errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=lambdathetaphi_->getError();

				//	  cout<<"+"<<lambdathetaphi_->getAsymErrorHi()<<"-"<<lambdathetaphi_->getAsymErrorLo()<<"+-"<<lambdathetaphi_->getError()<<endl;
				//	  cout<<errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]<<endl;

						  cout<<"lambda_phi"<<promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]<<endl;



					  }

					  else{
						  nPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
						  errnPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
						  promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
						  errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
			              promptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
						  errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
						  promptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;
						  errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;

					  }
					  }


					  if(workspace){
					  Closure_resfile->Close();
					  delete Closure_resfile;
					  }

					 // delete roofitres_Closure_migrad;
					 // delete roofitres_Closure_hesse;


				}

				cout<<"MISSING GENERATIONS!!! = "<<missingGenerations<<" (yBin= "<<yBinstart<<" , ptBin= "<<ptBinstart<<" frame= "<<frameDex<<")"<<endl;
				MissingGes[frameDex][yBinstart][ptBinstart]=missingGenerations;
	    }

	    //if(nofile) continue;
		  //cout<<"ddddddddd"<<endl;

}



if(!Closure){


	  cout<<"promptlambda_theta_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_theta_injected_CS[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_phi_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_phi_injected_CS[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart] = "<<promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_theta_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_theta_injected_HX[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_phi_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_phi_injected_HX[yBinstart][ptBinstart]<<endl;
	  cout<<"promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart] = "<<promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart]<<endl;




		  RooFitResult* roofitres;
		  RooFitResult* roofitres2;


		    for(int frameDex = 0; frameDex<2; ++frameDex){
					for(int generation = 1; generation < generations+0.5; generation++) {


		  bool onlyCrit1converged(false);
		  char resname[200];
		  char res2name[200];

		  if(frameDex==0){sprintf(resname,"cs_%s_rap%d_pt%d_gen%d",miniCrit,yBinstart,ptBinstart,generation); sprintf(res2name,"cs_%s_rap%d_pt%d_gen%d",miniCrit2,yBinstart,ptBinstart,generation);}
		  if(frameDex==1){sprintf(resname,"hx_%s_rap%d_pt%d_gen%d",miniCrit,yBinstart,ptBinstart,generation); sprintf(res2name,"hx_%s_rap%d_pt%d_gen%d",miniCrit2,yBinstart,ptBinstart,generation);}

		  roofitres = (RooFitResult*) resfile->Get(resname);
		  roofitres2 = (RooFitResult*) resfile->Get(res2name);


		  if(resfile->Get(resname) == NULL) {
			//  cout<<"nofile"<<endl;
		  nofile=true;
		  nPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  errnPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
          promptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  promptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  edm1[frameDex][yBinstart][ptBinstart][generation]=-9999;
		  edm2[frameDex][yBinstart][ptBinstart][generation]=-9999;

		  mean_phi[frameDex][yBinstart][ptBinstart]=-9999;
		  sigma_phi[frameDex][yBinstart][ptBinstart]=-9999;
		  mean_theta[frameDex][yBinstart][ptBinstart]=-9999;
		  sigma_theta[frameDex][yBinstart][ptBinstart]=-9999;
		  mean_thetaphi[frameDex][yBinstart][ptBinstart]=-9999;
		  sigma_thetaphi[frameDex][yBinstart][ptBinstart]=-9999;
		  mean_phi_val[frameDex][yBinstart][ptBinstart]=-9999;
		  sigma_phi_val[frameDex][yBinstart][ptBinstart]=-9999;
		  mean_theta_val[frameDex][yBinstart][ptBinstart]=-9999;
		  sigma_theta_val[frameDex][yBinstart][ptBinstart]=-9999;
		  mean_thetaphi_val[frameDex][yBinstart][ptBinstart]=-9999;
		  sigma_thetaphi_val[frameDex][yBinstart][ptBinstart]=-9999;

		  }
		  else {


		  roofitres = (RooFitResult*) resfile->Get(resname);
		  roofitres2 = (RooFitResult*) resfile->Get(res2name);


		  RooArgList params = roofitres2->floatParsFinal();
		  if(!TwoCrit || all) {/*cout<<"Values from 1st"<<endl;*/ params = roofitres->floatParsFinal();}

		  edm1[frameDex][yBinstart][ptBinstart][generation]=roofitres->edm();
		  edm2[frameDex][yBinstart][ptBinstart][generation]=roofitres2->edm();
		  status1[frameDex][yBinstart][ptBinstart][generation]=roofitres->status();
		  status2[frameDex][yBinstart][ptBinstart][generation]=roofitres2->status();
		  covQualy1[frameDex][yBinstart][ptBinstart][generation]=roofitres->covQual();
		  covQualy2[frameDex][yBinstart][ptBinstart][generation]=roofitres2->covQual();



		  bool criterium(false);
		  if(TwoCrit && logicOr && roofitres->covQual() == critCovVal && roofitres->edm() < edmcondition || TwoCrit && logicOr && roofitres2->covQual() == critCovVal && roofitres2->edm() < edmcondition) {criterium=true; if(roofitres2->covQual() != critCovVal) onlyCrit1converged=true;;};
		  if(TwoCrit && !logicOr && roofitres->covQual() == critCovVal && roofitres->edm() < edmcondition && roofitres2->covQual() == critCovVal && roofitres2->edm() < edmcondition ) {criterium=true;}
		  if (!TwoCrit && roofitres->edm() < edmcondition && roofitres->covQual() == critCovVal/* || roofitres->status() == 0*/) {criterium=true;}
		  if(all) criterium=true;

		  if(onlyCrit1converged) params = roofitres->floatParsFinal();

		  RooRealVar* nPrompt_;
		  RooRealVar* lambdaphi_;
		  RooRealVar* lambdatheta_;
		  RooRealVar* lambdathetaphi_;

		  if(normFitted){
		  nPrompt_ =  (RooRealVar*)params.at(0);
		  lambdaphi_ =  (RooRealVar*)params.at(1);
		  lambdatheta_ =  (RooRealVar*)params.at(2);
		  lambdathetaphi_ =  (RooRealVar*)params.at(3);
		  }

		  if(!normFitted){
		  lambdaphi_ =  (RooRealVar*)params.at(0);
		  lambdatheta_ =  (RooRealVar*)params.at(1);
		  lambdathetaphi_ =  (RooRealVar*)params.at(2);
		  }

		  if(criterium){
			//  cout<<"status  "<<roofitres->status()<<endl;
			//  cout<<"covqual "<<roofitres->covQual()<<endl;

			  convCount[frameDex][yBinstart][ptBinstart]++;
			  CONV[frameDex][yBinstart][ptBinstart][generation]=1;

/*
			  nPrompt[frameDex][yBinstart][ptBinstart][generation]=nPrompt_->getVal();
			  if(nPrompt_->getAsymErrorHi()>fabs(nPrompt_->getAsymErrorLo())) errnPrompt[frameDex][yBinstart][ptBinstart][generation]=nPrompt_->getAsymErrorHi();
			  else errnPrompt[frameDex][yBinstart][ptBinstart][generation]=fabs(nPrompt_->getAsymErrorLo());
			  if(errnPrompt[frameDex][yBinstart][ptBinstart][generation]==0) errnPrompt[frameDex][yBinstart][ptBinstart][generation]=nPrompt_->getError();
*/
			  promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=lambdaphi_->getVal();
			  if(lambdaphi_->getAsymErrorHi()>fabs(lambdaphi_->getAsymErrorLo())) errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=lambdaphi_->getAsymErrorHi();
			  else errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=fabs(lambdaphi_->getAsymErrorLo());
			  if(errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]==0) errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=lambdaphi_->getError();

			  promptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=lambdatheta_->getVal();
			  if(lambdatheta_->getAsymErrorHi()>fabs(lambdatheta_->getAsymErrorLo())) errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=lambdatheta_->getAsymErrorHi();
			  else errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=fabs(lambdatheta_->getAsymErrorLo());
			  if(errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]==0) errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=lambdatheta_->getError();

			  promptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=lambdathetaphi_->getVal();
			  if(lambdathetaphi_->getAsymErrorHi()>fabs(lambdathetaphi_->getAsymErrorLo())) errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=lambdathetaphi_->getAsymErrorHi();
			  else errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=fabs(lambdathetaphi_->getAsymErrorLo());
			  if(errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]==0) errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=lambdathetaphi_->getError();

	//	  cout<<"+"<<lambdathetaphi_->getAsymErrorHi()<<"-"<<lambdathetaphi_->getAsymErrorLo()<<"+-"<<lambdathetaphi_->getError()<<endl;
	//	  cout<<errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]<<endl;




		  }

		  else{
			  nPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
			  errnPrompt[frameDex][yBinstart][ptBinstart][generation]=-9999;
			  promptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
			  errpromptlambda_phi[frameDex][yBinstart][ptBinstart][generation]=-9999;
              promptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
			  errpromptlambda_theta[frameDex][yBinstart][ptBinstart][generation]=-9999;
			  promptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;
			  errpromptlambda_thetaphi[frameDex][yBinstart][ptBinstart][generation]=-9999;

		  }
		  }


					}}


			if(nofile) fprintf(outputFile2, "Bin rap%d_pt%d not yet delivered\n",yBinstart,ptBinstart);

		    if(nofile) continue;


}



///////////////////////// END OF FILE READING ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// FIT THE RESULTS WITH GAUSSIAN ////////////////////////////////////////////////





	    RooRealVar *promptlambda_phi_ = new RooRealVar("promptlambda_phi_","promptlambda_phi_",-cut,cut);
	    RooRealVar *promptlambda_theta_ = new RooRealVar("promptlambda_theta_","promptlambda_theta_",-cut,cut);
	    RooRealVar *promptlambda_thetaphi_ = new RooRealVar("promptlambda_thetaphi_","promptlambda_thetaphi_",-cut,cut);

	    RooRealVar *promptlambda_phi_val = new RooRealVar("promptlambda_phi_val","promptlambda_phi_val",-cutval,cutval);
		RooRealVar *promptlambda_theta_val = new RooRealVar("promptlambda_theta_val","promptlambda_theta_val",-cutval,cutval);
		RooRealVar *promptlambda_thetaphi_val = new RooRealVar("promptlambda_thetaphi_val","promptlambda_thetaphi_val",-cutval,cutval);

		double cutedm = 0.0001;
		RooRealVar *edm1_ = new RooRealVar("edm1_","edm1_",0,cutedm);
		RooRealVar *edm2_ = new RooRealVar("edm2_","edm2_",0,cutedm);

	    double meanstart=0;
	    double sigmastart=1;

	    RooRealVar mean("mean","mean",meanstart,-100,100);
	    RooRealVar sigma("sigma","sigma",sigmastart,0,1000);
	    RooGaussian gauss_phi("gauss_phi","gauss_phi",*promptlambda_phi_,mean,sigma);
	    RooGaussian gauss_theta("gauss_theta","gauss_theta",*promptlambda_theta_,mean,sigma);
	    RooGaussian gauss_thetaphi("gauss_thetaphi","gauss_thetaphi",*promptlambda_thetaphi_,mean,sigma);

	    RooGaussian gauss_phi_val("gauss_phi_val","gauss_phi_val",*promptlambda_phi_val,mean,sigma);
	    RooGaussian gauss_theta_val("gauss_theta_val","gauss_theta_val",*promptlambda_theta_val,mean,sigma);
	    RooGaussian gauss_thetaphi_val("gauss_thetaphi_val","gauss_thetaphi_val",*promptlambda_thetaphi_val,mean,sigma);

	    TCanvas* promptlambdaCanvas;
	    promptlambdaCanvas = new TCanvas("promptlambdaCanvas","promptlambdaCanvas",3800,3200);
	  	promptlambdaCanvas->Divide(2,3);  promptlambdaCanvas->SetFillColor(kWhite);

	  	TCanvas* promptlambdaCanvas_val;
	  	promptlambdaCanvas_val = new TCanvas("promptlambdaCanvas_val","promptlambdaCanvas_val",3800,3200);
	  	promptlambdaCanvas_val->Divide(2,3);  promptlambdaCanvas_val->SetFillColor(kWhite);

	    TCanvas* edmCanvas;
//	    edmCanvas = new TCanvas("edmCanvas","edmCanvas",1900,3200);
	    if(TwoCrit) {edmCanvas = new TCanvas("edmCanvas","edmCanvas",3200,1900); edmCanvas->Divide(2,2);}
	    else {edmCanvas = new TCanvas("edmCanvas","edmCanvas",1600,1900); edmCanvas->Divide(2);}
	    edmCanvas->SetFillColor(kWhite);

	    int PlotBins=20;//generations/5;
	    if(generations==1000) PlotBins=generations/20;


	    char Filename[200];
		sprintf(Filename,"%s/ToyMCPlot_rapidity%d_pt%d.png",dirStruct,yBinstart,ptBinstart);

//	    	for(int yBin = 1; yBin < 6; yBin++) {
//	    		for(int ptBin = 1; ptBin < 13; ptBin++) {

	    yBin=yBinstart;
	  	ptBin=ptBinstart;


	    		    for(int frameDex = 0; frameDex<2; ++frameDex){





	    		    	RooArgSet varlist(*promptlambda_phi_,*promptlambda_theta_,*promptlambda_thetaphi_,*promptlambda_phi_val,*promptlambda_theta_val,*promptlambda_thetaphi_val,*edm1_,*edm2_);
	    		    	RooDataSet* data = new RooDataSet("data","A sample",varlist);
	    		    	RooDataSet* dataRef = new RooDataSet("dataRef","A Ref sample",varlist);

	    if(CONV[frameDex][yBin][ptBin][generation]) fprintf(outputFile, "\n");
	    for(int generation = 1; generation < generations+1; generation++) {

	    	bool plotpull(true);

	    promptlambda_phi_val->setVal(promptlambda_phi[frameDex][yBin][ptBin][generation]);
	    promptlambda_theta_val->setVal(promptlambda_theta[frameDex][yBin][ptBin][generation]);
	    promptlambda_thetaphi_val->setVal(promptlambda_thetaphi[frameDex][yBin][ptBin][generation]);

	    if(frameDex==0){
	    promptlambda_phi_->setVal(-(promptlambda_phi_injected_CS[yBinstart][ptBinstart]-promptlambda_phi[frameDex][yBin][ptBin][generation])/errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    promptlambda_theta_->setVal(-(promptlambda_theta_injected_CS[yBinstart][ptBinstart]-promptlambda_theta[frameDex][yBin][ptBin][generation])/errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    promptlambda_thetaphi_->setVal(-(promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart]-promptlambda_thetaphi[frameDex][yBin][ptBin][generation])/errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
	    }
	    if(frameDex==1){
	    promptlambda_phi_->setVal(-(promptlambda_phi_injected_HX[yBinstart][ptBinstart]-promptlambda_phi[frameDex][yBin][ptBin][generation])/errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    promptlambda_theta_->setVal(-(promptlambda_theta_injected_HX[yBinstart][ptBinstart]-promptlambda_theta[frameDex][yBin][ptBin][generation])/errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    promptlambda_thetaphi_->setVal(-(promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart]-promptlambda_thetaphi[frameDex][yBin][ptBin][generation])/errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
	    }

	    edm1_->setVal(edm1[frameDex][yBin][ptBin][generation]);
	    edm2_->setVal(edm2[frameDex][yBin][ptBin][generation]);

	    if (edm1[frameDex][yBin][ptBin][generation]>cutedm) edm1_->setVal(cutedm-0.01);
	    if (edm2[frameDex][yBin][ptBin][generation]>cutedm) edm2_->setVal(cutedm-0.01);

	    double th_valerr=promptlambda_theta_->getVal();
	    double ph_valerr=promptlambda_phi_->getVal();
	    double thph_valerr=promptlambda_thetaphi_->getVal();

/*	    if(CONV[frameDex][yBin][ptBin][generation]<0.5){
		cout<<"NotConverged:"<<endl;
		cout<<" (generation) "<<generation<<" frame "<<frameDex<<" TH ERROR: "<<(promptlambda_theta[frameDex][yBin][ptBin][generation]-promptlambda_theta_injected)/errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<"; "<<promptlambda_theta[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<endl;
		cout<<" (generation) "<<generation<<" frame "<<frameDex<<" PH ERROR: "<<(promptlambda_phi[frameDex][yBin][ptBin][generation]-promptlambda_phi_injected)/errpromptlambda_phi[frameDex][yBin][ptBin][generation]<<"; "<<promptlambda_phi[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_phi[frameDex][yBin][ptBin][generation]<<endl;
		cout<<" (generation) "<<generation<<" frame "<<frameDex<<" TP ERROR: "<<(promptlambda_thetaphi[frameDex][yBin][ptBin][generation]-promptlambda_thetaphi_injected)/errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]<<"; "<<promptlambda_thetaphi[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]<<endl;
		cout<<miniCrit<<" status "<<status1[frameDex][yBin][ptBin][generation]<<" covQual "<<covQualy1[frameDex][yBin][ptBin][generation]<<", edm "<<edm1[frameDex][yBin][ptBin][generation]<<endl;
		cout<<miniCrit2<<" status "<<status2[frameDex][yBin][ptBin][generation]<<" covQual "<<covQualy2[frameDex][yBin][ptBin][generation]<<", edm "<<edm2[frameDex][yBin][ptBin][generation]<<endl;
	    }
*/	    if(CONV[frameDex][yBin][ptBin][generation]){
		dataRef->add(varlist);
		if(fabs(th_valerr)>=cut || fabs(ph_valerr)>=cut || fabs(thph_valerr)>=cut){
		cout<<"OUTLIER:"<<endl;
//		cout<<" (generation) "<<generation<<" frame "<<frameDex<<" TH ERROR: "<<(promptlambda_theta[frameDex][yBin][ptBin][generation]-promptlambda_theta_injected_CS[yBinstart][ptBinstart])/errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<"; "<<promptlambda_theta[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<endl;
//		cout<<" (generation) "<<generation<<" frame "<<frameDex<<" PH ERROR: "<<(promptlambda_phi[frameDex][yBin][ptBin][generation]-promptlambda_phi_injected_CS[yBinstart][ptBinstart])/errpromptlambda_phi[frameDex][yBin][ptBin][generation]<<"; "<<promptlambda_phi[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_phi[frameDex][yBin][ptBin][generation]<<endl;
//		cout<<" (generation) "<<generation<<" frame "<<frameDex<<" TP ERROR: "<<(promptlambda_thetaphi[frameDex][yBin][ptBin][generation]-promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart])/errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]<<"; "<<promptlambda_thetaphi[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]<<endl;
		cout<<miniCrit<<" status "<<status1[frameDex][yBin][ptBin][generation]<<" covQual "<<covQualy1[frameDex][yBin][ptBin][generation]<<", edm "<<edm1[frameDex][yBin][ptBin][generation]<<endl;
		cout<<miniCrit2<<" status "<<status2[frameDex][yBin][ptBin][generation]<<" covQual "<<covQualy2[frameDex][yBin][ptBin][generation]<<", edm "<<edm2[frameDex][yBin][ptBin][generation]<<endl;

		cout<<" "<<endl;
		}
	    }
	    if(CONV[frameDex][yBin][ptBin][generation] && fabs(th_valerr) < cut && fabs(ph_valerr) < cut && fabs(thph_valerr) < cut  ){
//	    cout<<"lambda_theta "<<generation<<" "<<promptlambda_theta[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<endl;
//	    cout<<promptlambda_theta[frameDex][yBin][ptBin][generation]/errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<endl;
	    data->add(varlist);
		if(frameDex==0) fprintf(outputFile, "CS rapidity%d_pt%d generation%d promptlambda_phi %1.4f +/- %1.4f\n",yBin,ptBin,generation,promptlambda_phi[frameDex][yBin][ptBin][generation],errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
		if(frameDex==1) fprintf(outputFile, "HX rapidity%d_pt%d generation%d promptlambda_phi %1.4f +/- %1.4f\n",yBin,ptBin,generation,promptlambda_phi[frameDex][yBin][ptBin][generation],errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    }
	    }
		if(CONV[frameDex][yBin][ptBin][generation]) fprintf(outputFile, "\n");
	    for(int generation = 1; generation < 100; generation++) {
	    	if(CONV[frameDex][yBin][ptBin][generation]){
	    if(frameDex==0) fprintf(outputFile, "CS rapidity%d_pt%d generation%d promptlambda_theta %1.4f +/- %1.4f\n",yBin,ptBin,generation,promptlambda_theta[frameDex][yBin][ptBin][generation],errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    if(frameDex==1) fprintf(outputFile, "HX rapidity%d_pt%d generation%d promptlambda_theta %1.4f +/- %1.4f\n",yBin,ptBin,generation,promptlambda_theta[frameDex][yBin][ptBin][generation],errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    	}}
	    if(CONV[frameDex][yBin][ptBin][generation]) fprintf(outputFile, "\n");
	    for(int generation = 1; generation < 100; generation++) {
	    	if(CONV[frameDex][yBin][ptBin][generation]){
	    if(frameDex==0) fprintf(outputFile, "CS rapidity%d_pt%d generation%d promptlambda_thetaphi %1.4f +/- %1.4f\n",yBin,ptBin,generation,promptlambda_thetaphi[frameDex][yBin][ptBin][generation],errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
	    if(frameDex==1) fprintf(outputFile, "HX rapidity%d_pt%d generation%d promptlambda_thetaphi %1.4f +/- %1.4f\n",yBin,ptBin,generation,promptlambda_thetaphi[frameDex][yBin][ptBin][generation],errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
	    	}}

	    data->Print();
	    dataRef->Print();

	    cout<<"mean_phi"<<data->mean(*promptlambda_phi_)<<endl;
	    cout<<"sigma_phi"<<data->sigma(*promptlambda_phi_)<<endl;
	    cout<<"Overflow = "<<dataRef->sumEntries()-data->sumEntries()<<endl;

	    overflow[frameDex][yBin][ptBin]=dataRef->sumEntries()-data->sumEntries();
	//    cout<<data->mean()<<endl;

/*	    	    			mean.setVal(meanstart);
	    	    			sigma.setVal(sigmastart);


	    		    		for(int i=0;i<5;i++){
	    		    RooFitResult* phi_result = gauss_phi.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Minos(0),RooFit::Strategy(2));
	    		    phi_result->Print(); cout<<"Fit number "<<i+1<<endl;
	    		    if(sigma.getVal()<1.5)continue;
	    		    		}
*/

	    	    	mean.setVal(data->meanVar(*promptlambda_phi_)->getVal());
	    	    	mean.setError(data->meanVar(*promptlambda_phi_)->getError());
	    	    	sigma.setVal(data->rmsVar(*promptlambda_phi_)->getVal());
	    	    	sigma.setError(data->rmsVar(*promptlambda_phi_)->getError());

	    		    mean_phi[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    		    errmean_phi[frameDex][yBinstart][ptBinstart]=mean.getError();
	    		    sigma_phi[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    		    errsigma_phi[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    		    RooPlot* promptlambda_phi_frame = new RooPlot;
	    		    promptlambda_phi_frame = promptlambda_phi_->frame(-borders,borders,PlotBins) ;
	    		    promptlambda_phi_frame->SetTitle(0);
	    		    if(frameDex==0) {promptlambda_phi_frame->SetTitleSize(0.06,"X"); if(!tilde) promptlambda_phi_frame->SetXTitle("z(#lambda_{#phi CS})"); else promptlambda_phi_frame->SetXTitle("z(#tilde{#lambda}_{CS})");}
	    		    if(frameDex==1) {promptlambda_phi_frame->SetTitleSize(0.06,"X"); if(!tilde) promptlambda_phi_frame->SetXTitle("z(#lambda_{#phi HX})"); else promptlambda_phi_frame->SetXTitle("z(#tilde{#lambda}_{HX})");}
	    		    data->plotOn(promptlambda_phi_frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    		    gauss_phi.plotOn(promptlambda_phi_frame,LineWidth(2),Normalization(1.0));
	    		    gauss_phi.paramOn(promptlambda_phi_frame,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));

/*	    			mean.setVal(meanstart);
	    			sigma.setVal(sigmastart);

	    			for(int i=0;i<5;i++){
	    			RooFitResult* theta_result = gauss_theta.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Strategy(2));
	    			theta_result->Print(); cout<<"Fit number "<<i+1<<endl;
	    		    if(sigma.getVal()<1.5)continue;
	    		    		}
*/
	    			mean.setVal(data->meanVar(*promptlambda_theta_)->getVal());
	    			mean.setError(data->meanVar(*promptlambda_theta_)->getError());
	    			sigma.setVal(data->rmsVar(*promptlambda_theta_)->getVal());
	    			sigma.setError(data->rmsVar(*promptlambda_theta_)->getError());

	    			mean_theta[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    			errmean_theta[frameDex][yBinstart][ptBinstart]=mean.getError();
	    			sigma_theta[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    			errsigma_theta[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    			RooPlot* promptlambda_theta_frame = new RooPlot;
	    			promptlambda_theta_frame = promptlambda_theta_->frame(-borders,borders,PlotBins) ;
	    			promptlambda_theta_frame->SetTitle(0);
	    			if(frameDex==0) {promptlambda_theta_frame->SetTitleSize(0.06,"X"); promptlambda_theta_frame->SetXTitle("z(#lambda_{#theta CS})");}
	    			if(frameDex==1) {promptlambda_theta_frame->SetTitleSize(0.06,"X"); promptlambda_theta_frame->SetXTitle("z(#lambda_{#theta HX})");}
	    			data->plotOn(promptlambda_theta_frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    			gauss_theta.plotOn(promptlambda_theta_frame,LineWidth(2),Normalization(1.0));
	    			gauss_theta.paramOn(promptlambda_theta_frame,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));

/*	    			mean.setVal(meanstart);
	    			sigma.setVal(sigmastart);


	    			for(int i=0;i<5;i++){
	    			RooFitResult* thetaphi_result = gauss_thetaphi.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Strategy(2));
	    			thetaphi_result->Print(); cout<<"Fit number "<<i+1<<endl;
	    		    if(sigma.getVal()<1.5)continue;
	    		    		}
*/

	    			mean.setVal(data->meanVar(*promptlambda_thetaphi_)->getVal());
	    			mean.setError(data->meanVar(*promptlambda_thetaphi_)->getError());
	    			sigma.setVal(data->rmsVar(*promptlambda_thetaphi_)->getVal());
	    			sigma.setError(data->rmsVar(*promptlambda_thetaphi_)->getError());

	    			mean_thetaphi[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    			errmean_thetaphi[frameDex][yBinstart][ptBinstart]=mean.getError();
	    			sigma_thetaphi[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    			errsigma_thetaphi[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    			cout<<"meancheck"<<mean_thetaphi[frameDex][yBinstart][ptBinstart]<<endl;

	    			RooPlot* promptlambda_thetaphi_frame = new RooPlot;
	    			promptlambda_thetaphi_frame = promptlambda_thetaphi_->frame(-borders,borders,PlotBins) ;
	    			promptlambda_thetaphi_frame->SetTitle(0); //promptlambda_thetaphi_frame->SetTitleOffset(0.005);
	    			if(frameDex==0) {promptlambda_thetaphi_frame->SetTitleSize(0.06,"X"); promptlambda_thetaphi_frame->SetXTitle("z(#lambda_{#theta#phi CS})");}
	    			if(frameDex==1) {promptlambda_thetaphi_frame->SetTitleSize(0.06,"X"); promptlambda_thetaphi_frame->SetXTitle("z(#lambda_{#theta#phi HX})");}
	    			data->plotOn(promptlambda_thetaphi_frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    			gauss_thetaphi.plotOn(promptlambda_thetaphi_frame,LineWidth(2),Normalization(1.0));
	    			gauss_thetaphi.paramOn(promptlambda_thetaphi_frame,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));



	    		    if(frameDex==0){
	    		    promptlambdaCanvas->cd(3) ; gPad->SetFillColor(kWhite); promptlambda_phi_frame->Draw();
	    		    promptlambdaCanvas->cd(1) ; gPad->SetFillColor(kWhite); promptlambda_theta_frame->Draw();
	    		    promptlambdaCanvas->cd(5) ; gPad->SetFillColor(kWhite); promptlambda_thetaphi_frame->Draw();
	    		    }
	    		    if(frameDex==1){
	    		   	promptlambdaCanvas->cd(4) ; gPad->SetFillColor(kWhite); promptlambda_phi_frame->Draw();
	    		   	promptlambdaCanvas->cd(2) ; gPad->SetFillColor(kWhite); promptlambda_theta_frame->Draw();
	    		   	promptlambdaCanvas->cd(6) ; gPad->SetFillColor(kWhite); promptlambda_thetaphi_frame->Draw();
	    		   	}





		    	  mean.setVal(dataRef->meanVar(*promptlambda_thetaphi_val)->getVal());
		    	  mean.setError(dataRef->meanVar(*promptlambda_thetaphi_val)->getError());
		    	  sigma.setVal(dataRef->rmsVar(*promptlambda_thetaphi_val)->getVal());
		    	  sigma.setError(dataRef->rmsVar(*promptlambda_thetaphi_val)->getError());

	    		  mean_thetaphi_val[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    		  errmean_thetaphi_val[frameDex][yBinstart][ptBinstart]=mean.getError();
	    		  sigma_thetaphi_val[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    		  errsigma_thetaphi_val[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    		  RooPlot* promptlambda_thetaphi_frame_val = new RooPlot;
	    		  if(frameDex==0) promptlambda_thetaphi_frame_val = promptlambda_thetaphi_val->frame(promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart]-borders_val,promptlambda_thetaphi_injected_CS[yBinstart][ptBinstart]+borders_val,PlotBins) ;
	    		  if(frameDex==1) promptlambda_thetaphi_frame_val = promptlambda_thetaphi_val->frame(promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart]-borders_val,promptlambda_thetaphi_injected_HX[yBinstart][ptBinstart]+borders_val,PlotBins) ;
	    		  promptlambda_thetaphi_frame_val->SetTitle(0); //promptlambda_thetaphi_frame->SetTitleOffset(0.005);
	    		  if(frameDex==0) {promptlambda_thetaphi_frame_val->SetTitleSize(0.06,"X"); promptlambda_thetaphi_frame_val->SetXTitle("#lambda_{#theta#phi CS}");}
	    		  if(frameDex==1) {promptlambda_thetaphi_frame_val->SetTitleSize(0.06,"X"); promptlambda_thetaphi_frame_val->SetXTitle("#lambda_{#theta#phi HX}");}
	    		  dataRef->plotOn(promptlambda_thetaphi_frame_val,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    		  gauss_thetaphi_val.plotOn(promptlambda_thetaphi_frame_val,LineWidth(2),Normalization(1.0));
	    		  gauss_thetaphi_val.paramOn(promptlambda_thetaphi_frame_val,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));

	    		  mean.setVal(dataRef->meanVar(*promptlambda_theta_val)->getVal());
	    		  mean.setError(dataRef->meanVar(*promptlambda_theta_val)->getError());
	    		  sigma.setVal(dataRef->rmsVar(*promptlambda_theta_val)->getVal());
	    		  sigma.setError(dataRef->rmsVar(*promptlambda_theta_val)->getError());

	    		  mean_theta_val[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    		  errmean_theta_val[frameDex][yBinstart][ptBinstart]=mean.getError();
	    		  sigma_theta_val[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    		  errsigma_theta_val[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    		  RooPlot* promptlambda_theta_frame_val = new RooPlot;
	    		  if(frameDex==0) promptlambda_theta_frame_val = promptlambda_theta_val->frame(promptlambda_theta_injected_CS[yBinstart][ptBinstart]-borders_val,promptlambda_theta_injected_CS[yBinstart][ptBinstart]+borders_val,PlotBins) ;
	    		  if(frameDex==1) promptlambda_theta_frame_val = promptlambda_theta_val->frame(promptlambda_theta_injected_HX[yBinstart][ptBinstart]-borders_val,promptlambda_theta_injected_HX[yBinstart][ptBinstart]+borders_val,PlotBins) ;
	    		  promptlambda_theta_frame_val->SetTitle(0); //promptlambda_theta_frame->SetTitleOffset(0.005);
	    		  if(frameDex==0) {promptlambda_theta_frame_val->SetTitleSize(0.06,"X"); promptlambda_theta_frame_val->SetXTitle("#lambda_{#theta CS}");}
	    		  if(frameDex==1) {promptlambda_theta_frame_val->SetTitleSize(0.06,"X"); promptlambda_theta_frame_val->SetXTitle("#lambda_{#theta HX}");}
	    		  dataRef->plotOn(promptlambda_theta_frame_val,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    		  gauss_theta_val.plotOn(promptlambda_theta_frame_val,LineWidth(2),Normalization(1.0));
	    		  gauss_theta_val.paramOn(promptlambda_theta_frame_val,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));

	    		  mean.setVal(dataRef->meanVar(*promptlambda_phi_val)->getVal());
	    		  mean.setError(dataRef->meanVar(*promptlambda_phi_val)->getError());
	    		  sigma.setVal(dataRef->rmsVar(*promptlambda_phi_val)->getVal());
	    		  sigma.setError(dataRef->rmsVar(*promptlambda_phi_val)->getError());

	    		  mean_phi_val[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    		  errmean_phi_val[frameDex][yBinstart][ptBinstart]=mean.getError();
	    		  sigma_phi_val[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    		  errsigma_phi_val[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    		  RooPlot* promptlambda_phi_frame_val = new RooPlot;
	    		  if(frameDex==0) promptlambda_phi_frame_val = promptlambda_phi_val->frame(promptlambda_phi_injected_CS[yBinstart][ptBinstart]-borders_val,promptlambda_phi_injected_CS[yBinstart][ptBinstart]+borders_val,PlotBins) ;
	    		  if(frameDex==1) promptlambda_phi_frame_val = promptlambda_phi_val->frame(promptlambda_phi_injected_HX[yBinstart][ptBinstart]-borders_val,promptlambda_phi_injected_HX[yBinstart][ptBinstart]+borders_val,PlotBins) ;
	    		  promptlambda_phi_frame_val->SetTitle(0); //promptlambda_phi_frame->SetTitleOffset(0.005);
	    		  if(frameDex==0) {promptlambda_phi_frame_val->SetTitleSize(0.06,"X");if(!tilde) promptlambda_phi_frame_val->SetXTitle("#lambda_{#phi CS}"); else promptlambda_phi_frame_val->SetXTitle("#tilde{#lambda}_{CS}");}
	    		  if(frameDex==1) {promptlambda_phi_frame_val->SetTitleSize(0.06,"X");if(!tilde) promptlambda_phi_frame_val->SetXTitle("#lambda_{#phi HX}"); else promptlambda_phi_frame_val->SetXTitle("#tilde{#lambda}_{HX}");}
	    		  dataRef->plotOn(promptlambda_phi_frame_val,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    		  gauss_phi_val.plotOn(promptlambda_phi_frame_val,LineWidth(2),Normalization(1.0));
	    		  gauss_phi_val.paramOn(promptlambda_phi_frame_val,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));


	    		if(frameDex==0){
	    			promptlambdaCanvas_val->cd(3) ; gPad->SetFillColor(kWhite); promptlambda_phi_frame_val->Draw();
	    			promptlambdaCanvas_val->cd(1) ; gPad->SetFillColor(kWhite); promptlambda_theta_frame_val->Draw();
	    			promptlambdaCanvas_val->cd(5) ; gPad->SetFillColor(kWhite); promptlambda_thetaphi_frame_val->Draw();
	    		}
	    		if(frameDex==1){
	    			promptlambdaCanvas_val->cd(4) ; gPad->SetFillColor(kWhite); promptlambda_phi_frame_val->Draw();
	    			promptlambdaCanvas_val->cd(2) ; gPad->SetFillColor(kWhite); promptlambda_theta_frame_val->Draw();
	    			promptlambdaCanvas_val->cd(6) ; gPad->SetFillColor(kWhite); promptlambda_thetaphi_frame_val->Draw();
	    		}


	    		char edmtitlename[200];


	    		  RooPlot* edm1frame = new RooPlot;
	    		  edm1frame = edm1_->frame(0,cutedm,50) ;
	    		  edm1frame->SetTitle(0);
	    		  if(frameDex==0) {edm1frame->SetTitleSize(0.06,"X"); sprintf(edmtitlename,"edm CS %s",miniCrit); edm1frame->SetXTitle(edmtitlename);}
	    		  if(frameDex==1) {edm1frame->SetTitleSize(0.06,"X"); sprintf(edmtitlename,"edm HX %s",miniCrit); edm1frame->SetXTitle(edmtitlename);}
	    		  dataRef->plotOn(edm1frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));

	    		  RooPlot* edm2frame = new RooPlot;
	    		  edm2frame = edm2_->frame(0,cutedm,50) ;
	    		  edm2frame->SetTitle(0);
	    		  if(frameDex==0) {edm2frame->SetTitleSize(0.06,"X"); sprintf(edmtitlename,"edm CS %s",miniCrit2); edm2frame->SetXTitle(edmtitlename);}
	    		  if(frameDex==1) {edm2frame->SetTitleSize(0.06,"X"); sprintf(edmtitlename,"edm HX %s",miniCrit2); edm2frame->SetXTitle(edmtitlename);}
	    		  dataRef->plotOn(edm2frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));

	    		  edmCanvas->SetLogx(1);

	    		  if(frameDex==0) {edmCanvas->cd(1) ; gPad->SetFillColor(kWhite); edm1frame->Draw(); if(TwoCrit) edmCanvas->cd(3) ; gPad->SetFillColor(kWhite); edm2frame->Draw(); }
	    		  if(frameDex==1) {edmCanvas->cd(2) ; gPad->SetFillColor(kWhite); edm1frame->Draw(); if(TwoCrit) edmCanvas->cd(4) ; gPad->SetFillColor(kWhite); edm2frame->Draw(); }

	    		    	delete data;
	    		    	delete dataRef;

	    		    }

//	    		}}


	    		    for(int frameDex = 0; frameDex<2; frameDex++){
	    		    for(int generation = 1; generation < generations+1; generation++) {
	    		//    	if(CONV[frameDex][yBin][ptBin][generation]) convCountCheck[frameDex][yBin][ptBin]++;
	    		//    	if(frameDex==0) cout<<"CS generation "<<generation<<": "<<CONV[frameDex][yBin][ptBin][generation]<<endl;
	    		//    	if(frameDex==1) cout<<"HX generation "<<generation<<": "<<CONV[frameDex][yBin][ptBin][generation]<<endl;

	    		    }}
	    		    cout<<convCount[0][yBinstart][ptBinstart]<<endl;
	    		    cout<<convCount[1][yBinstart][ptBinstart]<<endl;

//	    		    cout<<convCount[0][1][5]<<endl;//<<convCount[0][1][6]<<convCount[0][2][5]<<convCount[0][2][6]<<endl;
//	    		    cout<<convCount[1][1][5]<<endl;//<<convCount[1][1][6]<<convCount[1][2][5]<<convCount[1][2][6]<<endl;

	    		    cout<<"saving canvasses"<<endl;

	    		    promptlambdaCanvas->SaveAs(Filename);
	    		    promptlambdaCanvas->Close();

	    			sprintf(Filename,"%s/ToyMCPlot_val_rapidity%d_pt%d.png",dirStruct,yBinstart,ptBinstart);

	    		    promptlambdaCanvas_val->SaveAs(Filename);
	    		    promptlambdaCanvas_val->Close();

	    			sprintf(Filename,"%s/ToyMCPlot_edm_rapidity%d_pt%d.png",dirStruct,yBinstart,ptBinstart);

	    			if(saveedm) edmCanvas->SaveAs(Filename);
	    			edmCanvas->Close();



	    		    fprintf(outputFile, "\n");
	    			fprintf(outputFile, "CS convergence %d\n",convCountCheck[0][yBin][ptBin]);
	    			fprintf(outputFile, "HX convergence %d\n",convCountCheck[1][yBin][ptBin]);

	    			fprintf(outputFile2, "\n");
	    			fprintf(outputFile2, "Rapidity %d pT %d\n",yBin,ptBin);
	    			fprintf(outputFile2, "CS convergence %d\n",convCountCheck[0][yBin][ptBin]);
	    			fprintf(outputFile2, "HX convergence %d\n",convCountCheck[1][yBin][ptBin]);



	    			if(RealMC){

	    				    		for(int frameDex = 0; frameDex<2; frameDex++){
	    				    			 for(int yBin = 1; yBin < 3; yBin++) {
	    					    			 for(int ptBin = 1; ptBin < 9; ptBin++) {
	    				    			    		    	if(CONV[frameDex][yBin][ptBin][1]){
	    				    			    		    	if(frameDex==0) {
	    				    			    		    		fprintf(outputFile2, "CS rapidity%d_pt%d promptlambda_theta %1.4f +/- %1.4f\n",yBin,ptBin,promptlambda_theta[frameDex][yBin][ptBin][generation],errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    				    			    		    		fprintf(outputFile2, "CS rapidity%d_pt%d promptlambda_phi %1.4f +/- %1.4f\n",yBin,ptBin,promptlambda_phi[frameDex][yBin][ptBin][generation],errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    				    			    		    		fprintf(outputFile2, "CS rapidity%d_pt%d promptlambda_thetaphi %1.4f +/- %1.4f\n",yBin,ptBin,promptlambda_thetaphi[frameDex][yBin][ptBin][generation],errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);


	    				    			    		    	}
	    			  	    			    		        if(frameDex==1) {
	    			  	    			    		        	fprintf(outputFile2, "HX rapidity%d_pt%d promptlambda_theta %1.4f +/- %1.4f\n",yBin,ptBin,promptlambda_theta[frameDex][yBin][ptBin][generation],errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    			  	    			    		        	fprintf(outputFile2, "HX rapidity%d_pt%d promptlambda_phi %1.4f +/- %1.4f\n",yBin,ptBin,promptlambda_phi[frameDex][yBin][ptBin][generation],errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    			  	    			    		            fprintf(outputFile2, "HX rapidity%d_pt%d promptlambda_thetaphi %1.4f +/- %1.4f\n",yBin,ptBin,promptlambda_thetaphi[frameDex][yBin][ptBin][generation],errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);

	    			  	    			    		        }
	    				    			    		    	}
	    				    			    		    }}}



			    		    		fprintf(outputFile2, "\nlambda_theta_CS\nlambda_theta_HX\nlambda_phi_CS\nlambda_phi_HX\nlambda_thetaphi_CS\nlambda_thetaphi_HX\n(Errors in same order)\n\n");


	    				    		for(int var = 1; var<7; var++){


	    				    		for(int frameDex = 0; frameDex<2; frameDex++){
    			    		    		fprintf(outputFile2, "{");
	    				    			for(int yBin = 1; yBin < 3; yBin++) {
	    				    				 if(yBin!=1) fprintf(outputFile2, ",");
	    				    				 fprintf(outputFile2, "{");
	    					    			 for(int ptBin = 2; ptBin < 9; ptBin++) {
		 	     			    		    	 if(yBin==1 && ptBin>2 && ptBin<8) fprintf(outputFile2, ",");
	    					    				 if(yBin==2 && ptBin!=2) fprintf(outputFile2, ",");
	    				    			 		 if(yBin==1 && ptBin==8) continue;
	    					    				 if(CONV[frameDex][yBin][ptBin][1]){
	    				    			    		    	if (var==1)	fprintf(outputFile2, "%1.4f",promptlambda_theta[frameDex][yBin][ptBin][generation]);
	    				    			    		    	if (var==2)	fprintf(outputFile2, "%1.4f",promptlambda_phi[frameDex][yBin][ptBin][generation]);
	    				    			    		    	if (var==3)	fprintf(outputFile2, "%1.4f",promptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
	    				    			    		    	if (var==4)	fprintf(outputFile2, "%1.4f",errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    				    			    		    	if (var==5)	fprintf(outputFile2, "%1.4f",errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    				    			    		    	if (var==6)	fprintf(outputFile2, "%1.4f",errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);

	    				    			    	 }
	    					    			else {if(var<3.5) fprintf(outputFile2, "9999"); if(var>3.5) fprintf(outputFile2, "0");}}
	     			    		    		fprintf(outputFile2, "}");}
	    			    		    		fprintf(outputFile2, "}");
	    			    		    		if (var!=3 && var!= 6)fprintf(outputFile2, ",");
	    			    		    		if (var==3 && frameDex==0)fprintf(outputFile2, ",");
	    			    		    		if (var==6 && frameDex==0)fprintf(outputFile2, ",");
	    			    		    		fprintf(outputFile2, "\n");}
											if(var==3) fprintf(outputFile2, "\n");}


	    				    		fprintf(outputFile2, "\n");

	    				    		}


	    				    	//}

	    		    fclose(outputFile);

	    		    resfile->Close();
	    		    delete resfile;

	    		}}

			cout<<"meancheck1"<<mean_thetaphi[0][1][6]<<endl;
			cout<<"meancheck1"<<mean_thetaphi[1][1][6]<<endl;


			fprintf(outputFile2, "convergence counter CS = {%d,%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d,%d,%d}\n",convCountCheck[0][1][2],convCountCheck[0][1][3],convCountCheck[0][1][4],convCountCheck[0][1][5],convCountCheck[0][1][6],convCountCheck[0][1][7],convCountCheck[0][2][2],convCountCheck[0][2][3],convCountCheck[0][2][4],convCountCheck[0][2][5],convCountCheck[0][2][6],convCountCheck[0][2][7],convCountCheck[0][2][8]);
			fprintf(outputFile2, "convergence counter HX = {%d,%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d,%d,%d}\n",convCountCheck[1][1][2],convCountCheck[1][1][3],convCountCheck[1][1][4],convCountCheck[1][1][5],convCountCheck[1][1][6],convCountCheck[1][1][7],convCountCheck[1][2][2],convCountCheck[1][2][3],convCountCheck[1][2][4],convCountCheck[1][2][5],convCountCheck[1][2][6],convCountCheck[1][2][7],convCountCheck[1][2][8]);
			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "mean_phi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_phi[0][1][2],mean_phi[0][1][3],mean_phi[0][1][4],mean_phi[0][1][5],mean_phi[0][1][6],mean_phi[0][1][7],mean_phi[0][2][2],mean_phi[0][2][3],mean_phi[0][2][4],mean_phi[0][2][5],mean_phi[0][2][6],mean_phi[0][2][7],mean_phi[0][2][8]);
			fprintf(outputFile2, "errmean_phi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_phi[0][1][2],errmean_phi[0][1][3],errmean_phi[0][1][4],errmean_phi[0][1][5],errmean_phi[0][1][6],errmean_phi[0][1][7],errmean_phi[0][2][2],errmean_phi[0][2][3],errmean_phi[0][2][4],errmean_phi[0][2][5],errmean_phi[0][2][6],errmean_phi[0][2][7],errmean_phi[0][2][8]);
			fprintf(outputFile2, "mean_phi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_phi[1][1][2],mean_phi[1][1][3],mean_phi[1][1][4],mean_phi[1][1][5],mean_phi[1][1][6],mean_phi[1][1][7],mean_phi[1][2][2],mean_phi[1][2][3],mean_phi[1][2][4],mean_phi[1][2][5],mean_phi[1][2][6],mean_phi[1][2][7],mean_phi[1][2][8]);
			fprintf(outputFile2, "errmean_phi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_phi[1][1][2],errmean_phi[1][1][3],errmean_phi[1][1][4],errmean_phi[1][1][5],errmean_phi[1][1][6],errmean_phi[1][1][7],errmean_phi[1][2][2],errmean_phi[1][2][3],errmean_phi[1][2][4],errmean_phi[1][2][5],errmean_phi[1][2][6],errmean_phi[1][2][7],errmean_phi[1][2][8]);
			fprintf(outputFile2, "sigma_phi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_phi[0][1][2],sigma_phi[0][1][3],sigma_phi[0][1][4],sigma_phi[0][1][5],sigma_phi[0][1][6],sigma_phi[0][1][7],sigma_phi[0][2][2],sigma_phi[0][2][3],sigma_phi[0][2][4],sigma_phi[0][2][5],sigma_phi[0][2][6],sigma_phi[0][2][7],sigma_phi[0][2][8]);
			fprintf(outputFile2, "errsigma_phi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_phi[0][1][2],errsigma_phi[0][1][3],errsigma_phi[0][1][4],errsigma_phi[0][1][5],errsigma_phi[0][1][6],errsigma_phi[0][1][7],errsigma_phi[0][2][2],errsigma_phi[0][2][3],errsigma_phi[0][2][4],errsigma_phi[0][2][5],errsigma_phi[0][2][6],errsigma_phi[0][2][7],errsigma_phi[0][2][8]);
			fprintf(outputFile2, "sigma_phi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_phi[1][1][2],sigma_phi[1][1][3],sigma_phi[1][1][4],sigma_phi[1][1][5],sigma_phi[1][1][6],sigma_phi[1][1][7],sigma_phi[1][2][2],sigma_phi[1][2][3],sigma_phi[1][2][4],sigma_phi[1][2][5],sigma_phi[1][2][6],sigma_phi[1][2][7],sigma_phi[1][2][8]);
			fprintf(outputFile2, "errsigma_phi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_phi[1][1][2],errsigma_phi[1][1][3],errsigma_phi[1][1][4],errsigma_phi[1][1][5],errsigma_phi[1][1][6],errsigma_phi[1][1][7],errsigma_phi[1][2][2],errsigma_phi[1][2][3],errsigma_phi[1][2][4],errsigma_phi[1][2][5],errsigma_phi[1][2][6],errsigma_phi[1][2][7],errsigma_phi[1][2][8]);
			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "mean_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_theta[0][1][2],mean_theta[0][1][3],mean_theta[0][1][4],mean_theta[0][1][5],mean_theta[0][1][6],mean_theta[0][1][7],mean_theta[0][2][2],mean_theta[0][2][3],mean_theta[0][2][4],mean_theta[0][2][5],mean_theta[0][2][6],mean_theta[0][2][7],mean_theta[0][2][8]);
			fprintf(outputFile2, "errmean_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_theta[0][1][2],errmean_theta[0][1][3],errmean_theta[0][1][4],errmean_theta[0][1][5],errmean_theta[0][1][6],errmean_theta[0][1][7],errmean_theta[0][2][2],errmean_theta[0][2][3],errmean_theta[0][2][4],errmean_theta[0][2][5],errmean_theta[0][2][6],errmean_theta[0][2][7],errmean_theta[0][2][8]);
			fprintf(outputFile2, "mean_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_theta[1][1][2],mean_theta[1][1][3],mean_theta[1][1][4],mean_theta[1][1][5],mean_theta[1][1][6],mean_theta[1][1][7],mean_theta[1][2][2],mean_theta[1][2][3],mean_theta[1][2][4],mean_theta[1][2][5],mean_theta[1][2][6],mean_theta[1][2][7],mean_theta[1][2][8]);
			fprintf(outputFile2, "errmean_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_theta[1][1][2],errmean_theta[1][1][3],errmean_theta[1][1][4],errmean_theta[1][1][5],errmean_theta[1][1][6],errmean_theta[1][1][7],errmean_theta[1][2][2],errmean_theta[1][2][3],errmean_theta[1][2][4],errmean_theta[1][2][5],errmean_theta[1][2][6],errmean_theta[1][2][7],errmean_theta[1][2][8]);
			fprintf(outputFile2, "sigma_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_theta[0][1][2],sigma_theta[0][1][3],sigma_theta[0][1][4],sigma_theta[0][1][5],sigma_theta[0][1][6],sigma_theta[0][1][7],sigma_theta[0][2][2],sigma_theta[0][2][3],sigma_theta[0][2][4],sigma_theta[0][2][5],sigma_theta[0][2][6],sigma_theta[0][2][7],sigma_theta[0][2][8]);
			fprintf(outputFile2, "errsigma_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_theta[0][1][2],errsigma_theta[0][1][3],errsigma_theta[0][1][4],errsigma_theta[0][1][5],errsigma_theta[0][1][6],errsigma_theta[0][1][7],errsigma_theta[0][2][2],errsigma_theta[0][2][3],errsigma_theta[0][2][4],errsigma_theta[0][2][5],errsigma_theta[0][2][6],errsigma_theta[0][2][7],errsigma_theta[0][2][8]);
			fprintf(outputFile2, "sigma_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_theta[1][1][2],sigma_theta[1][1][3],sigma_theta[1][1][4],sigma_theta[1][1][5],sigma_theta[1][1][6],sigma_theta[1][1][7],sigma_theta[1][2][2],sigma_theta[1][2][3],sigma_theta[1][2][4],sigma_theta[1][2][5],sigma_theta[1][2][6],sigma_theta[1][2][7],sigma_theta[1][2][8]);
			fprintf(outputFile2, "errsigma_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_theta[1][1][2],errsigma_theta[1][1][3],errsigma_theta[1][1][4],errsigma_theta[1][1][5],errsigma_theta[1][1][6],errsigma_theta[1][1][7],errsigma_theta[1][2][2],errsigma_theta[1][2][3],errsigma_theta[1][2][4],errsigma_theta[1][2][5],errsigma_theta[1][2][6],errsigma_theta[1][2][7],errsigma_theta[1][2][8]);
			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "mean_thetaphi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_thetaphi[0][1][2],mean_thetaphi[0][1][3],mean_thetaphi[0][1][4],mean_thetaphi[0][1][5],mean_thetaphi[0][1][6],mean_thetaphi[0][1][7],mean_thetaphi[0][2][2],mean_thetaphi[0][2][3],mean_thetaphi[0][2][4],mean_thetaphi[0][2][5],mean_thetaphi[0][2][6],mean_thetaphi[0][2][7],mean_thetaphi[0][2][8]);
			fprintf(outputFile2, "errmean_thetaphi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_thetaphi[0][1][2],errmean_thetaphi[0][1][3],errmean_thetaphi[0][1][4],errmean_thetaphi[0][1][5],errmean_thetaphi[0][1][6],errmean_thetaphi[0][1][7],errmean_thetaphi[0][2][2],errmean_thetaphi[0][2][3],errmean_thetaphi[0][2][4],errmean_thetaphi[0][2][5],errmean_thetaphi[0][2][6],errmean_thetaphi[0][2][7],errmean_thetaphi[0][2][8]);
			fprintf(outputFile2, "mean_thetaphi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_thetaphi[1][1][2],mean_thetaphi[1][1][3],mean_thetaphi[1][1][4],mean_thetaphi[1][1][5],mean_thetaphi[1][1][6],mean_thetaphi[1][1][7],mean_thetaphi[1][2][2],mean_thetaphi[1][2][3],mean_thetaphi[1][2][4],mean_thetaphi[1][2][5],mean_thetaphi[1][2][6],mean_thetaphi[1][2][7],mean_thetaphi[1][2][8]);
			fprintf(outputFile2, "errmean_thetaphi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_thetaphi[1][1][2],errmean_thetaphi[1][1][3],errmean_thetaphi[1][1][4],errmean_thetaphi[1][1][5],errmean_thetaphi[1][1][6],errmean_thetaphi[1][1][7],errmean_thetaphi[1][2][2],errmean_thetaphi[1][2][3],errmean_thetaphi[1][2][4],errmean_thetaphi[1][2][5],errmean_thetaphi[1][2][6],errmean_thetaphi[1][2][7],errmean_thetaphi[1][2][8]);
			fprintf(outputFile2, "sigma_thetaphi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_thetaphi[0][1][2],sigma_thetaphi[0][1][3],sigma_thetaphi[0][1][4],sigma_thetaphi[0][1][5],sigma_thetaphi[0][1][6],sigma_thetaphi[0][1][7],sigma_thetaphi[0][2][2],sigma_thetaphi[0][2][3],sigma_thetaphi[0][2][4],sigma_thetaphi[0][2][5],sigma_thetaphi[0][2][6],sigma_thetaphi[0][2][7],sigma_thetaphi[0][2][8]);
			fprintf(outputFile2, "errsigma_thetaphi_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_thetaphi[0][1][2],errsigma_thetaphi[0][1][3],errsigma_thetaphi[0][1][4],errsigma_thetaphi[0][1][5],errsigma_thetaphi[0][1][6],errsigma_thetaphi[0][1][7],errsigma_thetaphi[0][2][2],errsigma_thetaphi[0][2][3],errsigma_thetaphi[0][2][4],errsigma_thetaphi[0][2][5],errsigma_thetaphi[0][2][6],errsigma_thetaphi[0][2][7],errsigma_thetaphi[0][2][8]);
			fprintf(outputFile2, "sigma_thetaphi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_thetaphi[1][1][2],sigma_thetaphi[1][1][3],sigma_thetaphi[1][1][4],sigma_thetaphi[1][1][5],sigma_thetaphi[1][1][6],sigma_thetaphi[1][1][7],sigma_thetaphi[1][2][2],sigma_thetaphi[1][2][3],sigma_thetaphi[1][2][4],sigma_thetaphi[1][2][5],sigma_thetaphi[1][2][6],sigma_thetaphi[1][2][7],sigma_thetaphi[1][2][8]);
			fprintf(outputFile2, "errsigma_thetaphi_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_thetaphi[1][1][2],errsigma_thetaphi[1][1][3],errsigma_thetaphi[1][1][4],errsigma_thetaphi[1][1][5],errsigma_thetaphi[1][1][6],errsigma_thetaphi[1][1][7],errsigma_thetaphi[1][2][2],errsigma_thetaphi[1][2][3],errsigma_thetaphi[1][2][4],errsigma_thetaphi[1][2][5],errsigma_thetaphi[1][2][6],errsigma_thetaphi[1][2][7],errsigma_thetaphi[1][2][8]);


			fprintf(outputFile2, "pull       CS\n");

			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_theta[0][1][2],errmean_theta[0][1][2],mean_theta[0][1][3],errmean_theta[0][1][3],mean_theta[0][1][4],errmean_theta[0][1][4],mean_theta[0][1][5],errmean_theta[0][1][5],mean_theta[0][1][6],errmean_theta[0][1][6],mean_theta[0][1][7],errmean_theta[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_theta[0][1][2],errsigma_theta[0][1][2],sigma_theta[0][1][3],errsigma_theta[0][1][3],sigma_theta[0][1][4],errsigma_theta[0][1][4],sigma_theta[0][1][5],errsigma_theta[0][1][5],sigma_theta[0][1][6],errsigma_theta[0][1][6],sigma_theta[0][1][7],errsigma_theta[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_phi[0][1][2],errmean_phi[0][1][2],mean_phi[0][1][3],errmean_phi[0][1][3],mean_phi[0][1][4],errmean_phi[0][1][4],mean_phi[0][1][5],errmean_phi[0][1][5],mean_phi[0][1][6],errmean_phi[0][1][6],mean_phi[0][1][7],errmean_phi[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_phi[0][1][2],errsigma_phi[0][1][2],sigma_phi[0][1][3],errsigma_phi[0][1][3],sigma_phi[0][1][4],errsigma_phi[0][1][4],sigma_phi[0][1][5],errsigma_phi[0][1][5],sigma_phi[0][1][6],errsigma_phi[0][1][6],sigma_phi[0][1][7],errsigma_phi[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_thetaphi[0][1][2],errmean_thetaphi[0][1][2],mean_thetaphi[0][1][3],errmean_thetaphi[0][1][3],mean_thetaphi[0][1][4],errmean_thetaphi[0][1][4],mean_thetaphi[0][1][5],errmean_thetaphi[0][1][5],mean_thetaphi[0][1][6],errmean_thetaphi[0][1][6],mean_thetaphi[0][1][7],errmean_thetaphi[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_thetaphi[0][1][2],errsigma_thetaphi[0][1][2],sigma_thetaphi[0][1][3],errsigma_thetaphi[0][1][3],sigma_thetaphi[0][1][4],errsigma_thetaphi[0][1][4],sigma_thetaphi[0][1][5],errsigma_thetaphi[0][1][5],sigma_thetaphi[0][1][6],errsigma_thetaphi[0][1][6],sigma_thetaphi[0][1][7],errsigma_thetaphi[0][1][7]);

			fprintf(outputFile2, "pull       HX\n");

			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_theta[1][1][2],errmean_theta[1][1][2],mean_theta[1][1][3],errmean_theta[1][1][3],mean_theta[1][1][4],errmean_theta[1][1][4],mean_theta[1][1][5],errmean_theta[1][1][5],mean_theta[1][1][6],errmean_theta[1][1][6],mean_theta[1][1][7],errmean_theta[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_theta[1][1][2],errsigma_theta[1][1][2],sigma_theta[1][1][3],errsigma_theta[1][1][3],sigma_theta[1][1][4],errsigma_theta[1][1][4],sigma_theta[1][1][5],errsigma_theta[1][1][5],sigma_theta[1][1][6],errsigma_theta[1][1][6],sigma_theta[1][1][7],errsigma_theta[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_phi[1][1][2],errmean_phi[1][1][2],mean_phi[1][1][3],errmean_phi[1][1][3],mean_phi[1][1][4],errmean_phi[1][1][4],mean_phi[1][1][5],errmean_phi[1][1][5],mean_phi[1][1][6],errmean_phi[1][1][6],mean_phi[1][1][7],errmean_phi[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_phi[1][1][2],errsigma_phi[1][1][2],sigma_phi[1][1][3],errsigma_phi[1][1][3],sigma_phi[1][1][4],errsigma_phi[1][1][4],sigma_phi[1][1][5],errsigma_phi[1][1][5],sigma_phi[1][1][6],errsigma_phi[1][1][6],sigma_phi[1][1][7],errsigma_phi[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_thetaphi[1][1][2],errmean_thetaphi[1][1][2],mean_thetaphi[1][1][3],errmean_thetaphi[1][1][3],mean_thetaphi[1][1][4],errmean_thetaphi[1][1][4],mean_thetaphi[1][1][5],errmean_thetaphi[1][1][5],mean_thetaphi[1][1][6],errmean_thetaphi[1][1][6],mean_thetaphi[1][1][7],errmean_thetaphi[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_thetaphi[1][1][2],errsigma_thetaphi[1][1][2],sigma_thetaphi[1][1][3],errsigma_thetaphi[1][1][3],sigma_thetaphi[1][1][4],errsigma_thetaphi[1][1][4],sigma_thetaphi[1][1][5],errsigma_thetaphi[1][1][5],sigma_thetaphi[1][1][6],errsigma_thetaphi[1][1][6],sigma_thetaphi[1][1][7],errsigma_thetaphi[1][1][7]);

			fprintf(outputFile2, "param       CS\n");

			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_theta_val[0][1][2],errmean_theta_val[0][1][2],mean_theta_val[0][1][3],errmean_theta_val[0][1][3],mean_theta_val[0][1][4],errmean_theta_val[0][1][4],mean_theta_val[0][1][5],errmean_theta_val[0][1][5],mean_theta_val[0][1][6],errmean_theta_val[0][1][6],mean_theta_val[0][1][7],errmean_theta_val[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_theta_val[0][1][2],errsigma_theta_val[0][1][2],sigma_theta_val[0][1][3],errsigma_theta_val[0][1][3],sigma_theta_val[0][1][4],errsigma_theta_val[0][1][4],sigma_theta_val[0][1][5],errsigma_theta_val[0][1][5],sigma_theta_val[0][1][6],errsigma_theta_val[0][1][6],sigma_theta_val[0][1][7],errsigma_theta_val[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_phi_val[0][1][2],errmean_phi_val[0][1][2],mean_phi_val[0][1][3],errmean_phi_val[0][1][3],mean_phi_val[0][1][4],errmean_phi_val[0][1][4],mean_phi_val[0][1][5],errmean_phi_val[0][1][5],mean_phi_val[0][1][6],errmean_phi_val[0][1][6],mean_phi_val[0][1][7],errmean_phi_val[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_phi_val[0][1][2],errsigma_phi_val[0][1][2],sigma_phi_val[0][1][3],errsigma_phi_val[0][1][3],sigma_phi_val[0][1][4],errsigma_phi_val[0][1][4],sigma_phi_val[0][1][5],errsigma_phi_val[0][1][5],sigma_phi_val[0][1][6],errsigma_phi_val[0][1][6],sigma_phi_val[0][1][7],errsigma_phi_val[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_thetaphi_val[0][1][2],errmean_thetaphi_val[0][1][2],mean_thetaphi_val[0][1][3],errmean_thetaphi_val[0][1][3],mean_thetaphi_val[0][1][4],errmean_thetaphi_val[0][1][4],mean_thetaphi_val[0][1][5],errmean_thetaphi_val[0][1][5],mean_thetaphi_val[0][1][6],errmean_thetaphi_val[0][1][6],mean_thetaphi_val[0][1][7],errmean_thetaphi_val[0][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_thetaphi_val[0][1][2],errsigma_thetaphi_val[0][1][2],sigma_thetaphi_val[0][1][3],errsigma_thetaphi_val[0][1][3],sigma_thetaphi_val[0][1][4],errsigma_thetaphi_val[0][1][4],sigma_thetaphi_val[0][1][5],errsigma_thetaphi_val[0][1][5],sigma_thetaphi_val[0][1][6],errsigma_thetaphi_val[0][1][6],sigma_thetaphi_val[0][1][7],errsigma_thetaphi_val[0][1][7]);

			fprintf(outputFile2, "param       HX\n");

			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_theta_val[1][1][2],errmean_theta_val[1][1][2],mean_theta_val[1][1][3],errmean_theta_val[1][1][3],mean_theta_val[1][1][4],errmean_theta_val[1][1][4],mean_theta_val[1][1][5],errmean_theta_val[1][1][5],mean_theta_val[1][1][6],errmean_theta_val[1][1][6],mean_theta_val[1][1][7],errmean_theta_val[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_theta_val[1][1][2],errsigma_theta_val[1][1][2],sigma_theta_val[1][1][3],errsigma_theta_val[1][1][3],sigma_theta_val[1][1][4],errsigma_theta_val[1][1][4],sigma_theta_val[1][1][5],errsigma_theta_val[1][1][5],sigma_theta_val[1][1][6],errsigma_theta_val[1][1][6],sigma_theta_val[1][1][7],errsigma_theta_val[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_phi_val[1][1][2],errmean_phi_val[1][1][2],mean_phi_val[1][1][3],errmean_phi_val[1][1][3],mean_phi_val[1][1][4],errmean_phi_val[1][1][4],mean_phi_val[1][1][5],errmean_phi_val[1][1][5],mean_phi_val[1][1][6],errmean_phi_val[1][1][6],mean_phi_val[1][1][7],errmean_phi_val[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_phi_val[1][1][2],errsigma_phi_val[1][1][2],sigma_phi_val[1][1][3],errsigma_phi_val[1][1][3],sigma_phi_val[1][1][4],errsigma_phi_val[1][1][4],sigma_phi_val[1][1][5],errsigma_phi_val[1][1][5],sigma_phi_val[1][1][6],errsigma_phi_val[1][1][6],sigma_phi_val[1][1][7],errsigma_phi_val[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",mean_thetaphi_val[1][1][2],errmean_thetaphi_val[1][1][2],mean_thetaphi_val[1][1][3],errmean_thetaphi_val[1][1][3],mean_thetaphi_val[1][1][4],errmean_thetaphi_val[1][1][4],mean_thetaphi_val[1][1][5],errmean_thetaphi_val[1][1][5],mean_thetaphi_val[1][1][6],errmean_thetaphi_val[1][1][6],mean_thetaphi_val[1][1][7],errmean_thetaphi_val[1][1][7]);
			fprintf(outputFile2, "%1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f & %1.3f       %1.3f\n",sigma_thetaphi_val[1][1][2],errsigma_thetaphi_val[1][1][2],sigma_thetaphi_val[1][1][3],errsigma_thetaphi_val[1][1][3],sigma_thetaphi_val[1][1][4],errsigma_thetaphi_val[1][1][4],sigma_thetaphi_val[1][1][5],errsigma_thetaphi_val[1][1][5],sigma_thetaphi_val[1][1][6],errsigma_thetaphi_val[1][1][6],sigma_thetaphi_val[1][1][7],errsigma_thetaphi_val[1][1][7]);

			fprintf(outputFile2, "sigmaCS_rap1\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_theta[0][1][2],sigma_theta[0][1][3],sigma_theta[0][1][4],sigma_theta[0][1][5],sigma_theta[0][1][6],sigma_theta[0][1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_phi[0][1][2],sigma_phi[0][1][3],sigma_phi[0][1][4],sigma_phi[0][1][5],sigma_phi[0][1][6],sigma_phi[0][1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_thetaphi[0][1][2],sigma_thetaphi[0][1][3],sigma_thetaphi[0][1][4],sigma_thetaphi[0][1][5],sigma_thetaphi[0][1][6],sigma_thetaphi[0][1][7]);

			fprintf(outputFile2, "mean_valCS_rap1\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_theta_val[0][1][2]-promptlambda_theta_injected_CS[1][2],mean_theta_val[0][1][3]-promptlambda_theta_injected_CS[1][3],mean_theta_val[0][1][4]-promptlambda_theta_injected_CS[1][4],mean_theta_val[0][1][5]-promptlambda_theta_injected_CS[1][5],mean_theta_val[0][1][6]-promptlambda_theta_injected_CS[1][6],mean_theta_val[0][1][7]-promptlambda_theta_injected_CS[1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_phi_val[0][1][2]-promptlambda_phi_injected_CS[1][2],mean_phi_val[0][1][3]-promptlambda_phi_injected_CS[1][3],mean_phi_val[0][1][4]-promptlambda_phi_injected_CS[1][4],mean_phi_val[0][1][5]-promptlambda_phi_injected_CS[1][5],mean_phi_val[0][1][6]-promptlambda_phi_injected_CS[1][6],mean_phi_val[0][1][7]-promptlambda_phi_injected_CS[1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_thetaphi_val[0][1][2]-promptlambda_thetaphi_injected_CS[1][2],mean_thetaphi_val[0][1][3]-promptlambda_thetaphi_injected_CS[1][3],mean_thetaphi_val[0][1][4]-promptlambda_thetaphi_injected_CS[1][4],mean_thetaphi_val[0][1][5]-promptlambda_thetaphi_injected_CS[1][5],mean_thetaphi_val[0][1][6]-promptlambda_thetaphi_injected_CS[1][6],mean_thetaphi_val[0][1][7]-promptlambda_thetaphi_injected_CS[1][7]);

			fprintf(outputFile2, "sigmaCS_rap2\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_theta[0][2][2],sigma_theta[0][2][3],sigma_theta[0][2][4],sigma_theta[0][2][5],sigma_theta[0][2][6],sigma_theta[0][2][7],sigma_theta[0][2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_phi[0][2][2],sigma_phi[0][2][3],sigma_phi[0][2][4],sigma_phi[0][2][5],sigma_phi[0][2][6],sigma_phi[0][2][7],sigma_phi[0][2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_thetaphi[0][2][2],sigma_thetaphi[0][2][3],sigma_thetaphi[0][2][4],sigma_thetaphi[0][2][5],sigma_thetaphi[0][2][6],sigma_thetaphi[0][2][7],sigma_thetaphi[0][2][8]);

			fprintf(outputFile2, "mean_valCS_rap2\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_theta_val[0][2][2]-promptlambda_theta_injected_CS[2][2],mean_theta_val[0][2][3]-promptlambda_theta_injected_CS[2][3],mean_theta_val[0][2][4]-promptlambda_theta_injected_CS[2][4],mean_theta_val[0][2][5]-promptlambda_theta_injected_CS[2][5],mean_theta_val[0][2][6]-promptlambda_theta_injected_CS[2][6],mean_theta_val[0][2][7]-promptlambda_theta_injected_CS[2][7],mean_theta_val[0][2][8]-promptlambda_theta_injected_CS[2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_phi_val[0][2][2]-promptlambda_phi_injected_CS[2][3],mean_phi_val[0][2][3]-promptlambda_phi_injected_CS[2][3],mean_phi_val[0][2][4]-promptlambda_phi_injected_CS[2][4],mean_phi_val[0][2][5]-promptlambda_phi_injected_CS[2][5],mean_phi_val[0][2][6]-promptlambda_phi_injected_CS[2][6],mean_phi_val[0][2][7]-promptlambda_phi_injected_CS[2][7],mean_phi_val[0][2][8]-promptlambda_phi_injected_CS[2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_thetaphi_val[0][2][2]-promptlambda_thetaphi_injected_CS[2][3],mean_thetaphi_val[0][2][3]-promptlambda_thetaphi_injected_CS[2][3],mean_thetaphi_val[0][2][4]-promptlambda_thetaphi_injected_CS[2][4],mean_thetaphi_val[0][2][5]-promptlambda_thetaphi_injected_CS[2][5],mean_thetaphi_val[0][2][6]-promptlambda_thetaphi_injected_CS[2][6],mean_thetaphi_val[0][2][7]-promptlambda_thetaphi_injected_CS[2][7],mean_thetaphi_val[0][2][8]-promptlambda_thetaphi_injected_CS[2][8]);

			fprintf(outputFile2, "sigmaHX_rap1\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_theta[1][1][2],sigma_theta[1][1][3],sigma_theta[1][1][4],sigma_theta[1][1][5],sigma_theta[1][1][6],sigma_theta[1][1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_phi[1][1][2],sigma_phi[1][1][3],sigma_phi[1][1][4],sigma_phi[1][1][5],sigma_phi[1][1][6],sigma_phi[1][1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_thetaphi[1][1][2],sigma_thetaphi[1][1][3],sigma_thetaphi[1][1][4],sigma_thetaphi[1][1][5],sigma_thetaphi[1][1][6],sigma_thetaphi[1][1][7]);

			fprintf(outputFile2, "mean_valHX_rap1\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_theta_val[1][1][2]-promptlambda_theta_injected_HX[1][2],mean_theta_val[1][1][3]-promptlambda_theta_injected_HX[1][3],mean_theta_val[1][1][4]-promptlambda_theta_injected_HX[1][4],mean_theta_val[1][1][5]-promptlambda_theta_injected_HX[1][5],mean_theta_val[1][1][6]-promptlambda_theta_injected_HX[1][6],mean_theta_val[1][1][7]-promptlambda_theta_injected_HX[1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_phi_val[1][1][2]-promptlambda_phi_injected_HX[1][2],mean_phi_val[1][1][3]-promptlambda_phi_injected_HX[1][3],mean_phi_val[1][1][4]-promptlambda_phi_injected_HX[1][4],mean_phi_val[1][1][5]-promptlambda_phi_injected_HX[1][5],mean_phi_val[1][1][6]-promptlambda_phi_injected_HX[1][6],mean_phi_val[1][1][7]-promptlambda_phi_injected_HX[1][7]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_thetaphi_val[1][1][2]-promptlambda_thetaphi_injected_HX[1][2],mean_thetaphi_val[1][1][3]-promptlambda_thetaphi_injected_HX[1][3],mean_thetaphi_val[1][1][4]-promptlambda_thetaphi_injected_HX[1][4],mean_thetaphi_val[1][1][5]-promptlambda_thetaphi_injected_HX[1][5],mean_thetaphi_val[1][1][6]-promptlambda_thetaphi_injected_HX[1][6],mean_thetaphi_val[1][1][7]-promptlambda_thetaphi_injected_HX[1][7]);

			fprintf(outputFile2, "sigmaHX_rap2\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_theta[1][2][2],sigma_theta[1][2][3],sigma_theta[1][2][4],sigma_theta[1][2][5],sigma_theta[1][2][6],sigma_theta[1][2][7],sigma_theta[1][2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_phi[1][2][2],sigma_phi[1][2][3],sigma_phi[1][2][4],sigma_phi[1][2][5],sigma_phi[1][2][6],sigma_phi[1][2][7],sigma_phi[1][2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",sigma_thetaphi[1][2][2],sigma_thetaphi[1][2][3],sigma_thetaphi[1][2][4],sigma_thetaphi[1][2][5],sigma_thetaphi[1][2][6],sigma_thetaphi[1][2][7],sigma_thetaphi[1][2][8]);

			fprintf(outputFile2, "mean_valHX_rap2\n");

			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_theta_val[1][2][2]-promptlambda_theta_injected_HX[2][2],mean_theta_val[1][2][3]-promptlambda_theta_injected_HX[2][3],mean_theta_val[1][2][4]-promptlambda_theta_injected_HX[2][4],mean_theta_val[1][2][5]-promptlambda_theta_injected_HX[2][5],mean_theta_val[1][2][6]-promptlambda_theta_injected_HX[2][6],mean_theta_val[1][2][7]-promptlambda_theta_injected_HX[2][7],mean_theta_val[1][2][8]-promptlambda_theta_injected_HX[2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_phi_val[1][2][2]-promptlambda_phi_injected_HX[2][2],mean_phi_val[1][2][3]-promptlambda_phi_injected_HX[2][3],mean_phi_val[1][2][4]-promptlambda_phi_injected_HX[2][4],mean_phi_val[1][2][5]-promptlambda_phi_injected_HX[2][5],mean_phi_val[1][2][6]-promptlambda_phi_injected_HX[2][6],mean_phi_val[1][2][7]-promptlambda_phi_injected_HX[2][7],mean_phi_val[1][2][8]-promptlambda_phi_injected_HX[2][8]);
			fprintf(outputFile2, "%1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f\n",mean_thetaphi_val[1][2][2]-promptlambda_thetaphi_injected_HX[2][2],mean_thetaphi_val[1][2][3]-promptlambda_thetaphi_injected_HX[2][3],mean_thetaphi_val[1][2][4]-promptlambda_thetaphi_injected_HX[2][4],mean_thetaphi_val[1][2][5]-promptlambda_thetaphi_injected_HX[2][5],mean_thetaphi_val[1][2][6]-promptlambda_thetaphi_injected_HX[2][6],mean_thetaphi_val[1][2][7]-promptlambda_thetaphi_injected_HX[2][7],mean_thetaphi_val[1][2][8]-promptlambda_thetaphi_injected_HX[2][8]);



//			fprintf(outputFile2, " & $|y|<0.9$, $15<pT<20$ & $|y|<0.9$, $20<pT<30$ & $0.9<|y|<1.2$ , $10<pT<15$ & $0.9<|y|<1.2$ , $10<pT<15$\n");
//			fprintf(outputFile2, "\\mu(\\lambda_{\\theta CS} - \\lambda^{Truth}_{\\theta CS} & %1.4f & %1.4f & %1.4f & %1.4f\n",mean_theta_val[0][1][2]-promptlambda_theta_injected_CS[1][2]);

			bool tableallbins(false);
			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "SLASHbegin{table}[!h]\nSLASHcentering SLASHcaption{Mean of the distribution of the standard score $z(SLASHlambda_{i})=(SLASHlambda_{i}-SLASHlambda^{Truth}_{i})/SLASHsigma(SLASHlambda_{i})$}\nSLASHbegin{tabular}{|cc|ccc|}\nSLASHhline\n");
			if(!tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHmu[z(SLASHlambda_{SLASHvartheta})]$ & $SLASHmu[z(SLASHlambda_{SLASHvarphi})]$ &  $SLASHmu[z(SLASHlambda_{SLASHvartheta SLASHvarphi})]$SLASHSLASH\n");
			if(tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHmu[z(SLASHlambda_{SLASHvartheta})]$ & $SLASHmu[z(SLASHtilde{SLASHlambda})]$ &  $SLASHmu[z(SLASHlambda_{SLASHvartheta SLASHvarphi})]$SLASHSLASH\n");
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{HX frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][1][2],errmean_theta[1][1][2],mean_phi[1][1][2],errmean_phi[1][1][2],mean_thetaphi[1][1][2],errmean_thetaphi[1][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][1][3],errmean_theta[1][1][3],mean_phi[1][1][3],errmean_phi[1][1][3],mean_thetaphi[1][1][3],errmean_thetaphi[1][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][1][4],errmean_theta[1][1][4],mean_phi[1][1][4],errmean_phi[1][1][4],mean_thetaphi[1][1][4],errmean_thetaphi[1][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][1][5],errmean_theta[1][1][5],mean_phi[1][1][5],errmean_phi[1][1][5],mean_thetaphi[1][1][5],errmean_thetaphi[1][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][1][6],errmean_theta[1][1][6],mean_phi[1][1][6],errmean_phi[1][1][6],mean_thetaphi[1][1][6],errmean_thetaphi[1][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][1][7],errmean_theta[1][1][7],mean_phi[1][1][7],errmean_phi[1][1][7],mean_thetaphi[1][1][7],errmean_thetaphi[1][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][2],errmean_theta[1][2][2],mean_phi[1][2][2],errmean_phi[1][2][2],mean_thetaphi[1][2][2],errmean_thetaphi[1][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][3],errmean_theta[1][2][3],mean_phi[1][2][3],errmean_phi[1][2][3],mean_thetaphi[1][2][3],errmean_thetaphi[1][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][4],errmean_theta[1][2][4],mean_phi[1][2][4],errmean_phi[1][2][4],mean_thetaphi[1][2][4],errmean_thetaphi[1][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][5],errmean_theta[1][2][5],mean_phi[1][2][5],errmean_phi[1][2][5],mean_thetaphi[1][2][5],errmean_thetaphi[1][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][6],errmean_theta[1][2][6],mean_phi[1][2][6],errmean_phi[1][2][6],mean_thetaphi[1][2][6],errmean_thetaphi[1][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][7],errmean_theta[1][2][7],mean_phi[1][2][7],errmean_phi[1][2][7],mean_thetaphi[1][2][7],errmean_thetaphi[1][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[1][2][8],errmean_theta[1][2][8],mean_phi[1][2][8],errmean_phi[1][2][8],mean_thetaphi[1][2][8],errmean_thetaphi[1][2][8]);
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{CS frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][1][2],errmean_theta[0][1][2],mean_phi[0][1][2],errmean_phi[0][1][2],mean_thetaphi[0][1][2],errmean_thetaphi[0][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][1][3],errmean_theta[0][1][3],mean_phi[0][1][3],errmean_phi[0][1][3],mean_thetaphi[0][1][3],errmean_thetaphi[0][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][1][4],errmean_theta[0][1][4],mean_phi[0][1][4],errmean_phi[0][1][4],mean_thetaphi[0][1][4],errmean_thetaphi[0][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][1][5],errmean_theta[0][1][5],mean_phi[0][1][5],errmean_phi[0][1][5],mean_thetaphi[0][1][5],errmean_thetaphi[0][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][1][6],errmean_theta[0][1][6],mean_phi[0][1][6],errmean_phi[0][1][6],mean_thetaphi[0][1][6],errmean_thetaphi[0][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][1][7],errmean_theta[0][1][7],mean_phi[0][1][7],errmean_phi[0][1][7],mean_thetaphi[0][1][7],errmean_thetaphi[0][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][2],errmean_theta[0][2][2],mean_phi[0][2][2],errmean_phi[0][2][2],mean_thetaphi[0][2][2],errmean_thetaphi[0][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][3],errmean_theta[0][2][3],mean_phi[0][2][3],errmean_phi[0][2][3],mean_thetaphi[0][2][3],errmean_thetaphi[0][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][4],errmean_theta[0][2][4],mean_phi[0][2][4],errmean_phi[0][2][4],mean_thetaphi[0][2][4],errmean_thetaphi[0][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][5],errmean_theta[0][2][5],mean_phi[0][2][5],errmean_phi[0][2][5],mean_thetaphi[0][2][5],errmean_thetaphi[0][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][6],errmean_theta[0][2][6],mean_phi[0][2][6],errmean_phi[0][2][6],mean_thetaphi[0][2][6],errmean_thetaphi[0][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][7],errmean_theta[0][2][7],mean_phi[0][2][7],errmean_phi[0][2][7],mean_thetaphi[0][2][7],errmean_thetaphi[0][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta[0][2][8],errmean_theta[0][2][8],mean_phi[0][2][8],errmean_phi[0][2][8],mean_thetaphi[0][2][8],errmean_thetaphi[0][2][8]);
			fprintf(outputFile2, "SLASHhline\n");
			fprintf(outputFile2, "SLASHend{tabular}\n");
			fprintf(outputFile2, "SLASHlabel{tab:syst_acceptance}\n");
			fprintf(outputFile2, "SLASHend{table}\n");

			fprintf(outputFile2, "\n");


			fprintf(outputFile2, "SLASHbegin{table}[!h]\nSLASHcentering SLASHcaption{R.m.s. of the distribution of the standard score $z(SLASHlambda_{i})=(SLASHlambda_{i}-SLASHlambda^{Truth}_{i})/SLASHsigma(SLASHlambda_{i})$}\nSLASHbegin{tabular}{|cc|ccc|}\nSLASHhline\n");
			if(!tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHsigma[z(SLASHlambda_{SLASHvartheta})]$ & $SLASHsigma[z(SLASHlambda_{SLASHvarphi})]$ &  $SLASHsigma[z(SLASHlambda_{SLASHvartheta SLASHvarphi})]$SLASHSLASH\n");
			if(tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHsigma[z(SLASHlambda_{SLASHvartheta})]$ & $SLASHsigma[z(SLASHtilde{SLASHlambda})]$ &  $SLASHsigma[z(SLASHlambda_{SLASHvartheta SLASHvarphi})]$SLASHSLASH\n");
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{HX frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][1][2],errsigma_theta[1][1][2],sigma_phi[1][1][2],errsigma_phi[1][1][2],sigma_thetaphi[1][1][2],errsigma_thetaphi[1][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][1][3],errsigma_theta[1][1][3],sigma_phi[1][1][3],errsigma_phi[1][1][3],sigma_thetaphi[1][1][3],errsigma_thetaphi[1][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][1][4],errsigma_theta[1][1][4],sigma_phi[1][1][4],errsigma_phi[1][1][4],sigma_thetaphi[1][1][4],errsigma_thetaphi[1][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][1][5],errsigma_theta[1][1][5],sigma_phi[1][1][5],errsigma_phi[1][1][5],sigma_thetaphi[1][1][5],errsigma_thetaphi[1][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][1][6],errsigma_theta[1][1][6],sigma_phi[1][1][6],errsigma_phi[1][1][6],sigma_thetaphi[1][1][6],errsigma_thetaphi[1][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][1][7],errsigma_theta[1][1][7],sigma_phi[1][1][7],errsigma_phi[1][1][7],sigma_thetaphi[1][1][7],errsigma_thetaphi[1][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][2],errsigma_theta[1][2][2],sigma_phi[1][2][2],errsigma_phi[1][2][2],sigma_thetaphi[1][2][2],errsigma_thetaphi[1][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][3],errsigma_theta[1][2][3],sigma_phi[1][2][3],errsigma_phi[1][2][3],sigma_thetaphi[1][2][3],errsigma_thetaphi[1][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][4],errsigma_theta[1][2][4],sigma_phi[1][2][4],errsigma_phi[1][2][4],sigma_thetaphi[1][2][4],errsigma_thetaphi[1][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][5],errsigma_theta[1][2][5],sigma_phi[1][2][5],errsigma_phi[1][2][5],sigma_thetaphi[1][2][5],errsigma_thetaphi[1][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][6],errsigma_theta[1][2][6],sigma_phi[1][2][6],errsigma_phi[1][2][6],sigma_thetaphi[1][2][6],errsigma_thetaphi[1][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][7],errsigma_theta[1][2][7],sigma_phi[1][2][7],errsigma_phi[1][2][7],sigma_thetaphi[1][2][7],errsigma_thetaphi[1][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[1][2][8],errsigma_theta[1][2][8],sigma_phi[1][2][8],errsigma_phi[1][2][8],sigma_thetaphi[1][2][8],errsigma_thetaphi[1][2][8]);
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{CS frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][1][2],errsigma_theta[0][1][2],sigma_phi[0][1][2],errsigma_phi[0][1][2],sigma_thetaphi[0][1][2],errsigma_thetaphi[0][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][1][3],errsigma_theta[0][1][3],sigma_phi[0][1][3],errsigma_phi[0][1][3],sigma_thetaphi[0][1][3],errsigma_thetaphi[0][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][1][4],errsigma_theta[0][1][4],sigma_phi[0][1][4],errsigma_phi[0][1][4],sigma_thetaphi[0][1][4],errsigma_thetaphi[0][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][1][5],errsigma_theta[0][1][5],sigma_phi[0][1][5],errsigma_phi[0][1][5],sigma_thetaphi[0][1][5],errsigma_thetaphi[0][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][1][6],errsigma_theta[0][1][6],sigma_phi[0][1][6],errsigma_phi[0][1][6],sigma_thetaphi[0][1][6],errsigma_thetaphi[0][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][1][7],errsigma_theta[0][1][7],sigma_phi[0][1][7],errsigma_phi[0][1][7],sigma_thetaphi[0][1][7],errsigma_thetaphi[0][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][2],errsigma_theta[0][2][2],sigma_phi[0][2][2],errsigma_phi[0][2][2],sigma_thetaphi[0][2][2],errsigma_thetaphi[0][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][3],errsigma_theta[0][2][3],sigma_phi[0][2][3],errsigma_phi[0][2][3],sigma_thetaphi[0][2][3],errsigma_thetaphi[0][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][4],errsigma_theta[0][2][4],sigma_phi[0][2][4],errsigma_phi[0][2][4],sigma_thetaphi[0][2][4],errsigma_thetaphi[0][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][5],errsigma_theta[0][2][5],sigma_phi[0][2][5],errsigma_phi[0][2][5],sigma_thetaphi[0][2][5],errsigma_thetaphi[0][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][6],errsigma_theta[0][2][6],sigma_phi[0][2][6],errsigma_phi[0][2][6],sigma_thetaphi[0][2][6],errsigma_thetaphi[0][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][7],errsigma_theta[0][2][7],sigma_phi[0][2][7],errsigma_phi[0][2][7],sigma_thetaphi[0][2][7],errsigma_thetaphi[0][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta[0][2][8],errsigma_theta[0][2][8],sigma_phi[0][2][8],errsigma_phi[0][2][8],sigma_thetaphi[0][2][8],errsigma_thetaphi[0][2][8]);
			fprintf(outputFile2, "SLASHhline\n");
			fprintf(outputFile2, "SLASHend{tabular}\n");
			fprintf(outputFile2, "SLASHlabel{tab:syst_acceptance}\n");
			fprintf(outputFile2, "SLASHend{table}\n");

			fprintf(outputFile2, "\n");


			fprintf(outputFile2, "SLASHbegin{table}[!h]\nSLASHcentering SLASHcaption{Mean deviation of fitted parameters from injected values , $SLASHDeltaSLASHlambda_{i} = SLASHmu(SLASHlambda_{i} - SLASHlambda^{Truth}_{i}$)}\nSLASHbegin{tabular}{|cc|ccc|}\nSLASHhline\n");
			if(!tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHDeltaSLASHlambda_{SLASHvartheta}$ & $SLASHDeltaSLASHlambda_{SLASHvarphi}$ &  $SLASHDeltaSLASHlambda_{SLASHvartheta SLASHvarphi}$SLASHSLASH\n");
			if(tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHDeltaSLASHlambda_{SLASHvartheta}$ & $SLASHDeltaSLASHtilde{SLASHlambda}$ &  $SLASHDeltaSLASHlambda_{SLASHvartheta SLASHvarphi}$SLASHSLASH\n");
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{HX frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][1][2]-promptlambda_theta_injected_HX[1][2],errmean_theta_val[1][1][2],mean_phi_val[1][1][2]-promptlambda_phi_injected_HX[1][2],errmean_phi_val[1][1][2],mean_thetaphi_val[1][1][2]-promptlambda_thetaphi_injected_HX[1][2],errmean_thetaphi_val[1][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][1][3]-promptlambda_theta_injected_HX[1][3],errmean_theta_val[1][1][3],mean_phi_val[1][1][3]-promptlambda_phi_injected_HX[1][3],errmean_phi_val[1][1][3],mean_thetaphi_val[1][1][3]-promptlambda_thetaphi_injected_HX[1][3],errmean_thetaphi_val[1][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][1][4]-promptlambda_theta_injected_HX[1][4],errmean_theta_val[1][1][4],mean_phi_val[1][1][4]-promptlambda_phi_injected_HX[1][4],errmean_phi_val[1][1][4],mean_thetaphi_val[1][1][4]-promptlambda_thetaphi_injected_HX[1][4],errmean_thetaphi_val[1][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][1][5]-promptlambda_theta_injected_HX[1][5],errmean_theta_val[1][1][5],mean_phi_val[1][1][5]-promptlambda_phi_injected_HX[1][5],errmean_phi_val[1][1][5],mean_thetaphi_val[1][1][5]-promptlambda_thetaphi_injected_HX[1][5],errmean_thetaphi_val[1][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][1][6]-promptlambda_theta_injected_HX[1][6],errmean_theta_val[1][1][6],mean_phi_val[1][1][6]-promptlambda_phi_injected_HX[1][6],errmean_phi_val[1][1][6],mean_thetaphi_val[1][1][6]-promptlambda_thetaphi_injected_HX[1][6],errmean_thetaphi_val[1][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][1][7]-promptlambda_theta_injected_HX[1][7],errmean_theta_val[1][1][7],mean_phi_val[1][1][7]-promptlambda_phi_injected_HX[1][7],errmean_phi_val[1][1][7],mean_thetaphi_val[1][1][7]-promptlambda_thetaphi_injected_HX[1][7],errmean_thetaphi_val[1][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][2]-promptlambda_theta_injected_HX[2][2],errmean_theta_val[1][2][2],mean_phi_val[1][2][2]-promptlambda_phi_injected_HX[2][2],errmean_phi_val[1][2][2],mean_thetaphi_val[1][2][2]-promptlambda_thetaphi_injected_HX[2][2],errmean_thetaphi_val[1][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][3]-promptlambda_theta_injected_HX[2][3],errmean_theta_val[1][2][3],mean_phi_val[1][2][3]-promptlambda_phi_injected_HX[2][3],errmean_phi_val[1][2][3],mean_thetaphi_val[1][2][3]-promptlambda_thetaphi_injected_HX[2][3],errmean_thetaphi_val[1][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][4]-promptlambda_theta_injected_HX[2][4],errmean_theta_val[1][2][4],mean_phi_val[1][2][4]-promptlambda_phi_injected_HX[2][4],errmean_phi_val[1][2][4],mean_thetaphi_val[1][2][4]-promptlambda_thetaphi_injected_HX[2][4],errmean_thetaphi_val[1][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][5]-promptlambda_theta_injected_HX[2][5],errmean_theta_val[1][2][5],mean_phi_val[1][2][5]-promptlambda_phi_injected_HX[2][5],errmean_phi_val[1][2][5],mean_thetaphi_val[1][2][5]-promptlambda_thetaphi_injected_HX[2][5],errmean_thetaphi_val[1][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][6]-promptlambda_theta_injected_HX[2][6],errmean_theta_val[1][2][6],mean_phi_val[1][2][6]-promptlambda_phi_injected_HX[2][6],errmean_phi_val[1][2][6],mean_thetaphi_val[1][2][6]-promptlambda_thetaphi_injected_HX[2][6],errmean_thetaphi_val[1][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][7]-promptlambda_theta_injected_HX[2][7],errmean_theta_val[1][2][7],mean_phi_val[1][2][7]-promptlambda_phi_injected_HX[2][7],errmean_phi_val[1][2][7],mean_thetaphi_val[1][2][7]-promptlambda_thetaphi_injected_HX[2][7],errmean_thetaphi_val[1][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[1][2][8]-promptlambda_theta_injected_HX[2][8],errmean_theta_val[1][2][8],mean_phi_val[1][2][8]-promptlambda_phi_injected_HX[2][8],errmean_phi_val[1][2][8],mean_thetaphi_val[1][2][8]-promptlambda_thetaphi_injected_HX[2][8],errmean_thetaphi_val[1][2][8]);
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{CS frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][1][2]-promptlambda_theta_injected_CS[1][2],errmean_theta_val[0][1][2],mean_phi_val[0][1][2]-promptlambda_phi_injected_CS[1][2],errmean_phi_val[0][1][2],mean_thetaphi_val[0][1][2]-promptlambda_thetaphi_injected_CS[1][2],errmean_thetaphi_val[0][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][1][3]-promptlambda_theta_injected_CS[1][3],errmean_theta_val[0][1][3],mean_phi_val[0][1][3]-promptlambda_phi_injected_CS[1][3],errmean_phi_val[0][1][3],mean_thetaphi_val[0][1][3]-promptlambda_thetaphi_injected_CS[1][3],errmean_thetaphi_val[0][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][1][4]-promptlambda_theta_injected_CS[1][4],errmean_theta_val[0][1][4],mean_phi_val[0][1][4]-promptlambda_phi_injected_CS[1][4],errmean_phi_val[0][1][4],mean_thetaphi_val[0][1][4]-promptlambda_thetaphi_injected_CS[1][4],errmean_thetaphi_val[0][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][1][5]-promptlambda_theta_injected_CS[1][5],errmean_theta_val[0][1][5],mean_phi_val[0][1][5]-promptlambda_phi_injected_CS[1][5],errmean_phi_val[0][1][5],mean_thetaphi_val[0][1][5]-promptlambda_thetaphi_injected_CS[1][5],errmean_thetaphi_val[0][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][1][6]-promptlambda_theta_injected_CS[1][6],errmean_theta_val[0][1][6],mean_phi_val[0][1][6]-promptlambda_phi_injected_CS[1][6],errmean_phi_val[0][1][6],mean_thetaphi_val[0][1][6]-promptlambda_thetaphi_injected_CS[1][6],errmean_thetaphi_val[0][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][1][7]-promptlambda_theta_injected_CS[1][7],errmean_theta_val[0][1][7],mean_phi_val[0][1][7]-promptlambda_phi_injected_CS[1][7],errmean_phi_val[0][1][7],mean_thetaphi_val[0][1][7]-promptlambda_thetaphi_injected_CS[1][7],errmean_thetaphi_val[0][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][2]-promptlambda_theta_injected_CS[2][2],errmean_theta_val[0][2][2],mean_phi_val[0][2][2]-promptlambda_phi_injected_CS[2][2],errmean_phi_val[0][2][2],mean_thetaphi_val[0][2][2]-promptlambda_thetaphi_injected_CS[2][2],errmean_thetaphi_val[0][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][3]-promptlambda_theta_injected_CS[2][3],errmean_theta_val[0][2][3],mean_phi_val[0][2][3]-promptlambda_phi_injected_CS[2][3],errmean_phi_val[0][2][3],mean_thetaphi_val[0][2][3]-promptlambda_thetaphi_injected_CS[2][3],errmean_thetaphi_val[0][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][4]-promptlambda_theta_injected_CS[2][4],errmean_theta_val[0][2][4],mean_phi_val[0][2][4]-promptlambda_phi_injected_CS[2][4],errmean_phi_val[0][2][4],mean_thetaphi_val[0][2][4]-promptlambda_thetaphi_injected_CS[2][4],errmean_thetaphi_val[0][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][5]-promptlambda_theta_injected_CS[2][5],errmean_theta_val[0][2][5],mean_phi_val[0][2][5]-promptlambda_phi_injected_CS[2][5],errmean_phi_val[0][2][5],mean_thetaphi_val[0][2][5]-promptlambda_thetaphi_injected_CS[2][5],errmean_thetaphi_val[0][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][6]-promptlambda_theta_injected_CS[2][6],errmean_theta_val[0][2][6],mean_phi_val[0][2][6]-promptlambda_phi_injected_CS[2][6],errmean_phi_val[0][2][6],mean_thetaphi_val[0][2][6]-promptlambda_thetaphi_injected_CS[2][6],errmean_thetaphi_val[0][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][7]-promptlambda_theta_injected_CS[2][7],errmean_theta_val[0][2][7],mean_phi_val[0][2][7]-promptlambda_phi_injected_CS[2][7],errmean_phi_val[0][2][7],mean_thetaphi_val[0][2][7]-promptlambda_thetaphi_injected_CS[2][7],errmean_thetaphi_val[0][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",mean_theta_val[0][2][8]-promptlambda_theta_injected_CS[2][8],errmean_theta_val[0][2][8],mean_phi_val[0][2][8]-promptlambda_phi_injected_CS[2][8],errmean_phi_val[0][2][8],mean_thetaphi_val[0][2][8]-promptlambda_thetaphi_injected_CS[2][8],errmean_thetaphi_val[0][2][8]);			fprintf(outputFile2, "SLASHhline\n");
			fprintf(outputFile2, "SLASHend{tabular}\n");
			fprintf(outputFile2, "SLASHlabel{tab:syst_acceptance}\n");
			fprintf(outputFile2, "SLASHend{table}\n");


			fprintf(outputFile2, "\n");


			fprintf(outputFile2, "SLASHbegin{table}[!h]\nSLASHcentering SLASHcaption{r.m.s. of the extracted $SLASHlambda_{i}$ distributions, $SLASHsigma(SLASHlambda_{i})$}\nSLASHbegin{tabular}{|cc|ccc|}\nSLASHhline\n");
			if(!tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHsigma(SLASHlambda_{SLASHvartheta})$ & $SLASHsigma(SLASHlambda_{SLASHvarphi})$ &  $SLASHsigma(SLASHlambda_{SLASHvartheta SLASHvarphi})$SLASHSLASH\n");
			if(tilde)fprintf(outputFile2, "$|y|$ & $p_{T}$ [GeV] & $SLASHsigma(SLASHlambda_{SLASHvartheta})$ & $SLASHsigma(SLASHtilde{SLASHlambda})$ &  $SLASHsigma(SLASHlambda_{SLASHvartheta SLASHvarphi})$SLASHSLASH\n");
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{HX frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][1][2],errsigma_theta_val[1][1][2],sigma_phi_val[1][1][2],errsigma_phi_val[1][1][2],sigma_thetaphi_val[1][1][2],errsigma_thetaphi_val[1][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][1][3],errsigma_theta_val[1][1][3],sigma_phi_val[1][1][3],errsigma_phi_val[1][1][3],sigma_thetaphi_val[1][1][3],errsigma_thetaphi_val[1][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][1][4],errsigma_theta_val[1][1][4],sigma_phi_val[1][1][4],errsigma_phi_val[1][1][4],sigma_thetaphi_val[1][1][4],errsigma_thetaphi_val[1][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][1][5],errsigma_theta_val[1][1][5],sigma_phi_val[1][1][5],errsigma_phi_val[1][1][5],sigma_thetaphi_val[1][1][5],errsigma_thetaphi_val[1][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][1][6],errsigma_theta_val[1][1][6],sigma_phi_val[1][1][6],errsigma_phi_val[1][1][6],sigma_thetaphi_val[1][1][6],errsigma_thetaphi_val[1][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][1][7],errsigma_theta_val[1][1][7],sigma_phi_val[1][1][7],errsigma_phi_val[1][1][7],sigma_thetaphi_val[1][1][7],errsigma_thetaphi_val[1][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][2],errsigma_theta_val[1][2][2],sigma_phi_val[1][2][2],errsigma_phi_val[1][2][2],sigma_thetaphi_val[1][2][2],errsigma_thetaphi_val[1][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][3],errsigma_theta_val[1][2][3],sigma_phi_val[1][2][3],errsigma_phi_val[1][2][3],sigma_thetaphi_val[1][2][3],errsigma_thetaphi_val[1][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][4],errsigma_theta_val[1][2][4],sigma_phi_val[1][2][4],errsigma_phi_val[1][2][4],sigma_thetaphi_val[1][2][4],errsigma_thetaphi_val[1][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][5],errsigma_theta_val[1][2][5],sigma_phi_val[1][2][5],errsigma_phi_val[1][2][5],sigma_thetaphi_val[1][2][5],errsigma_thetaphi_val[1][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][6],errsigma_theta_val[1][2][6],sigma_phi_val[1][2][6],errsigma_phi_val[1][2][6],sigma_thetaphi_val[1][2][6],errsigma_thetaphi_val[1][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][7],errsigma_theta_val[1][2][7],sigma_phi_val[1][2][7],errsigma_phi_val[1][2][7],sigma_thetaphi_val[1][2][7],errsigma_thetaphi_val[1][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[1][2][8],errsigma_theta_val[1][2][8],sigma_phi_val[1][2][8],errsigma_phi_val[1][2][8],sigma_thetaphi_val[1][2][8],errsigma_thetaphi_val[1][2][8]);
			fprintf(outputFile2, "SLASHhline SLASHmulticolumn{5}{|c|}{CS frame}SLASHSLASH SLASHhline SLASHrule{0pt}{4mm}\n");
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][1][2],errsigma_theta_val[0][1][2],sigma_phi_val[0][1][2],errsigma_phi_val[0][1][2],sigma_thetaphi_val[0][1][2],errsigma_thetaphi_val[0][1][2]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][1][3],errsigma_theta_val[0][1][3],sigma_phi_val[0][1][3],errsigma_phi_val[0][1][3],sigma_thetaphi_val[0][1][3],errsigma_thetaphi_val[0][1][3]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][1][4],errsigma_theta_val[0][1][4],sigma_phi_val[0][1][4],errsigma_phi_val[0][1][4],sigma_thetaphi_val[0][1][4],errsigma_thetaphi_val[0][1][4]);
			if(tableallbins)fprintf(outputFile2, "0.0--0.9 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][1][5],errsigma_theta_val[0][1][5],sigma_phi_val[0][1][5],errsigma_phi_val[0][1][5],sigma_thetaphi_val[0][1][5],errsigma_thetaphi_val[0][1][5]);
			fprintf(outputFile2, "0.0--0.9 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][1][6],errsigma_theta_val[0][1][6],sigma_phi_val[0][1][6],errsigma_phi_val[0][1][6],sigma_thetaphi_val[0][1][6],errsigma_thetaphi_val[0][1][6]);
			fprintf(outputFile2, "0.0--0.9 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][1][7],errsigma_theta_val[0][1][7],sigma_phi_val[0][1][7],errsigma_phi_val[0][1][7],sigma_thetaphi_val[0][1][7],errsigma_thetaphi_val[0][1][7]);
			fprintf(outputFile2, "SLASHhline\n");
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 4--6   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][2],errsigma_theta_val[0][2][2],sigma_phi_val[0][2][2],errsigma_phi_val[0][2][2],sigma_thetaphi_val[0][2][2],errsigma_thetaphi_val[0][2][2]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 6--7   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][3],errsigma_theta_val[0][2][3],sigma_phi_val[0][2][3],errsigma_phi_val[0][2][3],sigma_thetaphi_val[0][2][3],errsigma_thetaphi_val[0][2][3]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 7--8   &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][4],errsigma_theta_val[0][2][4],sigma_phi_val[0][2][4],errsigma_phi_val[0][2][4],sigma_thetaphi_val[0][2][4],errsigma_thetaphi_val[0][2][4]);
			if(tableallbins)fprintf(outputFile2, "0.9--1.2 & 8--10  &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][5],errsigma_theta_val[0][2][5],sigma_phi_val[0][2][5],errsigma_phi_val[0][2][5],sigma_thetaphi_val[0][2][5],errsigma_thetaphi_val[0][2][5]);
			fprintf(outputFile2, "0.9--1.2 & 10--15 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][6],errsigma_theta_val[0][2][6],sigma_phi_val[0][2][6],errsigma_phi_val[0][2][6],sigma_thetaphi_val[0][2][6],errsigma_thetaphi_val[0][2][6]);
			fprintf(outputFile2, "0.9--1.2 & 15--20 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][7],errsigma_theta_val[0][2][7],sigma_phi_val[0][2][7],errsigma_phi_val[0][2][7],sigma_thetaphi_val[0][2][7],errsigma_thetaphi_val[0][2][7]);
			fprintf(outputFile2, "0.9--1.2 & 20--30 &  $%1.4f SLASHpm %1.4f$  & $%1.4f SLASHpm %1.4f$  &  $%1.4f SLASHpm %1.4f$ SLASHSLASH\n",sigma_theta_val[0][2][8],errsigma_theta_val[0][2][8],sigma_phi_val[0][2][8],errsigma_phi_val[0][2][8],sigma_thetaphi_val[0][2][8],errsigma_thetaphi_val[0][2][8]);
			fprintf(outputFile2, "SLASHhline\n");
			fprintf(outputFile2, "SLASHend{tabular}\n");
			fprintf(outputFile2, "SLASHlabel{tab:syst_acceptance}\n");
			fprintf(outputFile2, "SLASHend{table}\n");








			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "mean_phi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_phi_val[0][1][2],mean_phi_val[0][1][3],mean_phi_val[0][1][4],mean_phi_val[0][1][5],mean_phi_val[0][1][6],mean_phi_val[0][1][7],mean_phi_val[0][2][2],mean_phi_val[0][2][3],mean_phi_val[0][2][4],mean_phi_val[0][2][5],mean_phi_val[0][2][6],mean_phi_val[0][2][7],mean_phi_val[0][2][8]);
			fprintf(outputFile2, "errmean_phi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_phi_val[0][1][2],errmean_phi_val[0][1][3],errmean_phi_val[0][1][4],errmean_phi_val[0][1][5],errmean_phi_val[0][1][6],errmean_phi_val[0][1][7],errmean_phi_val[0][2][2],errmean_phi_val[0][2][3],errmean_phi_val[0][2][4],errmean_phi_val[0][2][5],errmean_phi_val[0][2][6],errmean_phi_val[0][2][7],errmean_phi_val[0][2][8]);
			fprintf(outputFile2, "mean_phi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_phi_val[1][1][2],mean_phi_val[1][1][3],mean_phi_val[1][1][4],mean_phi_val[1][1][5],mean_phi_val[1][1][6],mean_phi_val[1][1][7],mean_phi_val[1][2][2],mean_phi_val[1][2][3],mean_phi_val[1][2][4],mean_phi_val[1][2][5],mean_phi_val[1][2][6],mean_phi_val[1][2][7],mean_phi_val[1][2][8]);
			fprintf(outputFile2, "errmean_phi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_phi_val[1][1][2],errmean_phi_val[1][1][3],errmean_phi_val[1][1][4],errmean_phi_val[1][1][5],errmean_phi_val[1][1][6],errmean_phi_val[1][1][7],errmean_phi_val[1][2][2],errmean_phi_val[1][2][3],errmean_phi_val[1][2][4],errmean_phi_val[1][2][5],errmean_phi_val[1][2][6],errmean_phi_val[1][2][7],errmean_phi_val[1][2][8]);
			fprintf(outputFile2, "sigma_phi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_phi_val[0][1][2],sigma_phi_val[0][1][3],sigma_phi_val[0][1][4],sigma_phi_val[0][1][5],sigma_phi_val[0][1][6],sigma_phi_val[0][1][7],sigma_phi_val[0][2][2],sigma_phi_val[0][2][3],sigma_phi_val[0][2][4],sigma_phi_val[0][2][5],sigma_phi_val[0][2][6],sigma_phi_val[0][2][7],sigma_phi_val[0][2][8]);
			fprintf(outputFile2, "errsigma_phi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_phi_val[0][1][2],errsigma_phi_val[0][1][3],errsigma_phi_val[0][1][4],errsigma_phi_val[0][1][5],errsigma_phi_val[0][1][6],errsigma_phi_val[0][1][7],errsigma_phi_val[0][2][2],errsigma_phi_val[0][2][3],errsigma_phi_val[0][2][4],errsigma_phi_val[0][2][5],errsigma_phi_val[0][2][6],errsigma_phi_val[0][2][7],errsigma_phi_val[0][2][8]);
			fprintf(outputFile2, "sigma_phi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_phi_val[1][1][2],sigma_phi_val[1][1][3],sigma_phi_val[1][1][4],sigma_phi_val[1][1][5],sigma_phi_val[1][1][6],sigma_phi_val[1][1][7],sigma_phi_val[1][2][2],sigma_phi_val[1][2][3],sigma_phi_val[1][2][4],sigma_phi_val[1][2][5],sigma_phi_val[1][2][6],sigma_phi_val[1][2][7],sigma_phi_val[1][2][8]);
			fprintf(outputFile2, "errsigma_phi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_phi_val[1][1][2],errsigma_phi_val[1][1][3],errsigma_phi_val[1][1][4],errsigma_phi_val[1][1][5],errsigma_phi_val[1][1][6],errsigma_phi_val[1][1][7],errsigma_phi_val[1][2][2],errsigma_phi_val[1][2][3],errsigma_phi_val[1][2][4],errsigma_phi_val[1][2][5],errsigma_phi_val[1][2][6],errsigma_phi_val[1][2][7],errsigma_phi_val[1][2][8]);
			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "mean_theta_valCS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_theta_val[0][1][2],mean_theta_val[0][1][3],mean_theta_val[0][1][4],mean_theta_val[0][1][5],mean_theta_val[0][1][6],mean_theta_val[0][1][7],mean_theta_val[0][2][2],mean_theta_val[0][2][3],mean_theta_val[0][2][4],mean_theta_val[0][2][5],mean_theta_val[0][2][6],mean_theta_val[0][2][7],mean_theta_val[0][2][8]);
			fprintf(outputFile2, "errmean_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_theta_val[0][1][2],errmean_theta_val[0][1][3],errmean_theta_val[0][1][4],errmean_theta_val[0][1][5],errmean_theta_val[0][1][6],errmean_theta_val[0][1][7],errmean_theta_val[0][2][2],errmean_theta_val[0][2][3],errmean_theta_val[0][2][4],errmean_theta_val[0][2][5],errmean_theta_val[0][2][6],errmean_theta_val[0][2][7],errmean_theta_val[0][2][8]);
			fprintf(outputFile2, "mean_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_theta_val[1][1][2],mean_theta_val[1][1][3],mean_theta_val[1][1][4],mean_theta_val[1][1][5],mean_theta_val[1][1][6],mean_theta_val[1][1][7],mean_theta_val[1][2][2],mean_theta_val[1][2][3],mean_theta_val[1][2][4],mean_theta_val[1][2][5],mean_theta_val[1][2][6],mean_theta_val[1][2][7],mean_theta_val[1][2][8]);
			fprintf(outputFile2, "errmean_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_theta_val[1][1][2],errmean_theta_val[1][1][3],errmean_theta_val[1][1][4],errmean_theta_val[1][1][5],errmean_theta_val[1][1][6],errmean_theta_val[1][1][7],errmean_theta_val[1][2][2],errmean_theta_val[1][2][3],errmean_theta_val[1][2][4],errmean_theta_val[1][2][5],errmean_theta_val[1][2][6],errmean_theta_val[1][2][7],errmean_theta_val[1][2][8]);
			fprintf(outputFile2, "sigma_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_theta_val[0][1][2],sigma_theta_val[0][1][3],sigma_theta_val[0][1][4],sigma_theta_val[0][1][5],sigma_theta_val[0][1][6],sigma_theta_val[0][1][7],sigma_theta_val[0][2][2],sigma_theta_val[0][2][3],sigma_theta_val[0][2][4],sigma_theta_val[0][2][5],sigma_theta_val[0][2][6],sigma_theta_val[0][2][7],sigma_theta_val[0][2][8]);
			fprintf(outputFile2, "errsigma_theta_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_theta_val[0][1][2],errsigma_theta_val[0][1][3],errsigma_theta_val[0][1][4],errsigma_theta_val[0][1][5],errsigma_theta_val[0][1][6],errsigma_theta_val[0][1][7],errsigma_theta_val[0][2][2],errsigma_theta_val[0][2][3],errsigma_theta_val[0][2][4],errsigma_theta_val[0][2][5],errsigma_theta_val[0][2][6],errsigma_theta_val[0][2][7],errsigma_theta_val[0][2][8]);
			fprintf(outputFile2, "sigma_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_theta_val[1][1][2],sigma_theta_val[1][1][3],sigma_theta_val[1][1][4],sigma_theta_val[1][1][5],sigma_theta_val[1][1][6],sigma_theta_val[1][1][7],sigma_theta_val[1][2][2],sigma_theta_val[1][2][3],sigma_theta_val[1][2][4],sigma_theta_val[1][2][5],sigma_theta_val[1][2][6],sigma_theta_val[1][2][7],sigma_theta_val[1][2][8]);
			fprintf(outputFile2, "errsigma_theta_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_theta_val[1][1][2],errsigma_theta_val[1][1][3],errsigma_theta_val[1][1][4],errsigma_theta_val[1][1][5],errsigma_theta_val[1][1][6],errsigma_theta_val[1][1][7],errsigma_theta_val[1][2][2],errsigma_theta_val[1][2][3],errsigma_theta_val[1][2][4],errsigma_theta_val[1][2][5],errsigma_theta_val[1][2][6],errsigma_theta_val[1][2][7],errsigma_theta_val[1][2][8]);
			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "mean_thetaphi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_thetaphi_val[0][1][2],mean_thetaphi_val[0][1][3],mean_thetaphi_val[0][1][4],mean_thetaphi_val[0][1][5],mean_thetaphi_val[0][1][6],mean_thetaphi_val[0][1][7],mean_thetaphi_val[0][2][2],mean_thetaphi_val[0][2][3],mean_thetaphi_val[0][2][4],mean_thetaphi_val[0][2][5],mean_thetaphi_val[0][2][6],mean_thetaphi_val[0][2][7],mean_thetaphi_val[0][2][8]);
			fprintf(outputFile2, "errmean_thetaphi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_thetaphi_val[0][1][2],errmean_thetaphi_val[0][1][3],errmean_thetaphi_val[0][1][4],errmean_thetaphi_val[0][1][5],errmean_thetaphi_val[0][1][6],errmean_thetaphi_val[0][1][7],errmean_thetaphi_val[0][2][2],errmean_thetaphi_val[0][2][3],errmean_thetaphi_val[0][2][4],errmean_thetaphi_val[0][2][5],errmean_thetaphi_val[0][2][6],errmean_thetaphi_val[0][2][7],errmean_thetaphi_val[0][2][8]);
			fprintf(outputFile2, "mean_thetaphi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",mean_thetaphi_val[1][1][2],mean_thetaphi_val[1][1][3],mean_thetaphi_val[1][1][4],mean_thetaphi_val[1][1][5],mean_thetaphi_val[1][1][6],mean_thetaphi_val[1][1][7],mean_thetaphi_val[1][2][2],mean_thetaphi_val[1][2][3],mean_thetaphi_val[1][2][4],mean_thetaphi_val[1][2][5],mean_thetaphi_val[1][2][6],mean_thetaphi_val[1][2][7],mean_thetaphi_val[1][2][8]);
			fprintf(outputFile2, "errmean_thetaphi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errmean_thetaphi_val[1][1][2],errmean_thetaphi_val[1][1][3],errmean_thetaphi_val[1][1][4],errmean_thetaphi_val[1][1][5],errmean_thetaphi_val[1][1][6],errmean_thetaphi_val[1][1][7],errmean_thetaphi_val[1][2][2],errmean_thetaphi_val[1][2][3],errmean_thetaphi_val[1][2][4],errmean_thetaphi_val[1][2][5],errmean_thetaphi_val[1][2][6],errmean_thetaphi_val[1][2][7],errmean_thetaphi_val[1][2][8]);
			fprintf(outputFile2, "sigma_thetaphi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_thetaphi_val[0][1][2],sigma_thetaphi_val[0][1][3],sigma_thetaphi_val[0][1][4],sigma_thetaphi_val[0][1][5],sigma_thetaphi_val[0][1][6],sigma_thetaphi_val[0][1][7],sigma_thetaphi_val[0][2][2],sigma_thetaphi_val[0][2][3],sigma_thetaphi_val[0][2][4],sigma_thetaphi_val[0][2][5],sigma_thetaphi_val[0][2][6],sigma_thetaphi_val[0][2][7],sigma_thetaphi_val[0][2][8]);
			fprintf(outputFile2, "errsigma_thetaphi_val_CS = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_thetaphi_val[0][1][2],errsigma_thetaphi_val[0][1][3],errsigma_thetaphi_val[0][1][4],errsigma_thetaphi_val[0][1][5],errsigma_thetaphi_val[0][1][6],errsigma_thetaphi_val[0][1][7],errsigma_thetaphi_val[0][2][2],errsigma_thetaphi_val[0][2][3],errsigma_thetaphi_val[0][2][4],errsigma_thetaphi_val[0][2][5],errsigma_thetaphi_val[0][2][6],errsigma_thetaphi_val[0][2][7],errsigma_thetaphi_val[0][2][8]);
			fprintf(outputFile2, "sigma_thetaphi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",sigma_thetaphi_val[1][1][2],sigma_thetaphi_val[1][1][3],sigma_thetaphi_val[1][1][4],sigma_thetaphi_val[1][1][5],sigma_thetaphi_val[1][1][6],sigma_thetaphi_val[1][1][7],sigma_thetaphi_val[1][2][2],sigma_thetaphi_val[1][2][3],sigma_thetaphi_val[1][2][4],sigma_thetaphi_val[1][2][5],sigma_thetaphi_val[1][2][6],sigma_thetaphi_val[1][2][7],sigma_thetaphi_val[1][2][8]);
			fprintf(outputFile2, "errsigma_thetaphi_val_HX = {%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}\n",errsigma_thetaphi_val[1][1][2],errsigma_thetaphi_val[1][1][3],errsigma_thetaphi_val[1][1][4],errsigma_thetaphi_val[1][1][5],errsigma_thetaphi_val[1][1][6],errsigma_thetaphi_val[1][1][7],errsigma_thetaphi_val[1][2][2],errsigma_thetaphi_val[1][2][3],errsigma_thetaphi_val[1][2][4],errsigma_thetaphi_val[1][2][5],errsigma_thetaphi_val[1][2][6],errsigma_thetaphi_val[1][2][7],errsigma_thetaphi_val[1][2][8]);

			fprintf(outputFile2, "\n");



			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",mean_phi[0][1][2],mean_phi[0][1][3],mean_phi[0][1][4],mean_phi[0][1][5],mean_phi[0][1][6],mean_phi[0][1][7],mean_phi[0][2][2],mean_phi[0][2][3],mean_phi[0][2][4],mean_phi[0][2][5],mean_phi[0][2][6],mean_phi[0][2][7],mean_phi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",mean_theta[0][1][2],mean_theta[0][1][3],mean_theta[0][1][4],mean_theta[0][1][5],mean_theta[0][1][6],mean_theta[0][1][7],mean_theta[0][2][2],mean_theta[0][2][3],mean_theta[0][2][4],mean_theta[0][2][5],mean_theta[0][2][6],mean_theta[0][2][7],mean_theta[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",mean_thetaphi[0][1][2],mean_thetaphi[0][1][3],mean_thetaphi[0][1][4],mean_thetaphi[0][1][5],mean_thetaphi[0][1][6],mean_thetaphi[0][1][7],mean_thetaphi[0][2][2],mean_thetaphi[0][2][3],mean_thetaphi[0][2][4],mean_thetaphi[0][2][5],mean_thetaphi[0][2][6],mean_thetaphi[0][2][7],mean_thetaphi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",mean_phi[1][1][2],mean_phi[1][1][3],mean_phi[1][1][4],mean_phi[1][1][5],mean_phi[1][1][6],mean_phi[1][1][7],mean_phi[1][2][2],mean_phi[1][2][3],mean_phi[1][2][4],mean_phi[1][2][5],mean_phi[1][2][6],mean_phi[1][2][7],mean_phi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",mean_theta[1][1][2],mean_theta[1][1][3],mean_theta[1][1][4],mean_theta[1][1][5],mean_theta[1][1][6],mean_theta[1][1][7],mean_theta[1][2][2],mean_theta[1][2][3],mean_theta[1][2][4],mean_theta[1][2][5],mean_theta[1][2][6],mean_theta[1][2][7],mean_theta[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",mean_thetaphi[1][1][2],mean_thetaphi[1][1][3],mean_thetaphi[1][1][4],mean_thetaphi[1][1][5],mean_thetaphi[1][1][6],mean_thetaphi[1][1][7],mean_thetaphi[1][2][2],mean_thetaphi[1][2][3],mean_thetaphi[1][2][4],mean_thetaphi[1][2][5],mean_thetaphi[1][2][6],mean_thetaphi[1][2][7],mean_thetaphi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",sigma_phi[0][1][2],sigma_phi[0][1][3],sigma_phi[0][1][4],sigma_phi[0][1][5],sigma_phi[0][1][6],sigma_phi[0][1][7],sigma_phi[0][2][2],sigma_phi[0][2][3],sigma_phi[0][2][4],sigma_phi[0][2][5],sigma_phi[0][2][6],sigma_phi[0][2][7],sigma_phi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",sigma_theta[0][1][2],sigma_theta[0][1][3],sigma_theta[0][1][4],sigma_theta[0][1][5],sigma_theta[0][1][6],sigma_theta[0][1][7],sigma_theta[0][2][2],sigma_theta[0][2][3],sigma_theta[0][2][4],sigma_theta[0][2][5],sigma_theta[0][2][6],sigma_theta[0][2][7],sigma_theta[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",sigma_thetaphi[0][1][2],sigma_thetaphi[0][1][3],sigma_thetaphi[0][1][4],sigma_thetaphi[0][1][5],sigma_thetaphi[0][1][6],sigma_thetaphi[0][1][7],sigma_thetaphi[0][2][2],sigma_thetaphi[0][2][3],sigma_thetaphi[0][2][4],sigma_thetaphi[0][2][5],sigma_thetaphi[0][2][6],sigma_thetaphi[0][2][7],sigma_thetaphi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",sigma_phi[1][1][2],sigma_phi[1][1][3],sigma_phi[1][1][4],sigma_phi[1][1][5],sigma_phi[1][1][6],sigma_phi[1][1][7],sigma_phi[1][2][2],sigma_phi[1][2][3],sigma_phi[1][2][4],sigma_phi[1][2][5],sigma_phi[1][2][6],sigma_phi[1][2][7],sigma_phi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",sigma_theta[1][1][2],sigma_theta[1][1][3],sigma_theta[1][1][4],sigma_theta[1][1][5],sigma_theta[1][1][6],sigma_theta[1][1][7],sigma_theta[1][2][2],sigma_theta[1][2][3],sigma_theta[1][2][4],sigma_theta[1][2][5],sigma_theta[1][2][6],sigma_theta[1][2][7],sigma_theta[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",sigma_thetaphi[1][1][2],sigma_thetaphi[1][1][3],sigma_thetaphi[1][1][4],sigma_thetaphi[1][1][5],sigma_thetaphi[1][1][6],sigma_thetaphi[1][1][7],sigma_thetaphi[1][2][2],sigma_thetaphi[1][2][3],sigma_thetaphi[1][2][4],sigma_thetaphi[1][2][5],sigma_thetaphi[1][2][6],sigma_thetaphi[1][2][7],sigma_thetaphi[1][2][8]);
			fprintf(outputFile2, "{{%d,%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d,%d,%d}},\n",convCount[0][1][2],convCount[0][1][3],convCount[0][1][4],convCount[0][1][5],convCount[0][1][6],convCount[0][1][7],convCount[0][2][2],convCount[0][2][3],convCount[0][2][4],convCount[0][2][5],convCount[0][2][6],convCount[0][2][7],convCount[0][2][8]);
			fprintf(outputFile2, "{{%d,%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d,%d,%d}},\n",convCount[1][1][2],convCount[1][1][3],convCount[1][1][4],convCount[1][1][5],convCount[1][1][6],convCount[1][1][7],convCount[1][2][2],convCount[1][2][3],convCount[1][2][4],convCount[1][2][5],convCount[1][2][6],convCount[1][2][7],convCount[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",overflow[0][1][2],overflow[0][1][3],overflow[0][1][4],overflow[0][1][5],overflow[0][1][6],overflow[0][1][7],overflow[0][2][2],overflow[0][2][3],overflow[0][2][4],overflow[0][2][5],overflow[0][2][6],overflow[0][2][7],overflow[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}}\n",overflow[1][1][2],overflow[1][1][3],overflow[1][1][4],overflow[1][1][5],overflow[1][1][6],overflow[1][1][7],overflow[1][2][2],overflow[1][2][3],overflow[1][2][4],overflow[1][2][5],overflow[1][2][6],overflow[1][2][7],overflow[1][2][8]);

			fprintf(outputFile2, " ");

			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errmean_phi[0][1][2],errmean_phi[0][1][3],errmean_phi[0][1][4],errmean_phi[0][1][5],errmean_phi[0][1][6],errmean_phi[0][1][7],errmean_phi[0][2][2],errmean_phi[0][2][3],errmean_phi[0][2][4],errmean_phi[0][2][5],errmean_phi[0][2][6],errmean_phi[0][2][7],errmean_phi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errmean_theta[0][1][2],errmean_theta[0][1][3],errmean_theta[0][1][4],errmean_theta[0][1][5],errmean_theta[0][1][6],errmean_theta[0][1][7],errmean_theta[0][2][2],errmean_theta[0][2][3],errmean_theta[0][2][4],errmean_theta[0][2][5],errmean_theta[0][2][6],errmean_theta[0][2][7],errmean_theta[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errmean_thetaphi[0][1][2],errmean_thetaphi[0][1][3],errmean_thetaphi[0][1][4],errmean_thetaphi[0][1][5],errmean_thetaphi[0][1][6],errmean_thetaphi[0][1][7],errmean_thetaphi[0][2][2],errmean_thetaphi[0][2][3],errmean_thetaphi[0][2][4],errmean_thetaphi[0][2][5],errmean_thetaphi[0][2][6],errmean_thetaphi[0][2][7],errmean_thetaphi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errmean_phi[1][1][2],errmean_phi[1][1][3],errmean_phi[1][1][4],errmean_phi[1][1][5],errmean_phi[1][1][6],errmean_phi[1][1][7],errmean_phi[1][2][2],errmean_phi[1][2][3],errmean_phi[1][2][4],errmean_phi[1][2][5],errmean_phi[1][2][6],errmean_phi[1][2][7],errmean_phi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errmean_theta[1][1][2],errmean_theta[1][1][3],errmean_theta[1][1][4],errmean_theta[1][1][5],errmean_theta[1][1][6],errmean_theta[1][1][7],errmean_theta[1][2][2],errmean_theta[1][2][3],errmean_theta[1][2][4],errmean_theta[1][2][5],errmean_theta[1][2][6],errmean_theta[1][2][7],errmean_theta[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errmean_thetaphi[1][1][2],errmean_thetaphi[1][1][3],errmean_thetaphi[1][1][4],errmean_thetaphi[1][1][5],errmean_thetaphi[1][1][6],errmean_thetaphi[1][1][7],errmean_thetaphi[1][2][2],errmean_thetaphi[1][2][3],errmean_thetaphi[1][2][4],errmean_thetaphi[1][2][5],errmean_thetaphi[1][2][6],errmean_thetaphi[1][2][7],errmean_thetaphi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errsigma_phi[0][1][2],errsigma_phi[0][1][3],errsigma_phi[0][1][4],errsigma_phi[0][1][5],errsigma_phi[0][1][6],errsigma_phi[0][1][7],errsigma_phi[0][2][2],errsigma_phi[0][2][3],errsigma_phi[0][2][4],errsigma_phi[0][2][5],errsigma_phi[0][2][6],errsigma_phi[0][2][7],errsigma_phi[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errsigma_theta[0][1][2],errsigma_theta[0][1][3],errsigma_theta[0][1][4],errsigma_theta[0][1][5],errsigma_theta[0][1][6],errsigma_theta[0][1][7],errsigma_theta[0][2][2],errsigma_theta[0][2][3],errsigma_theta[0][2][4],errsigma_theta[0][2][5],errsigma_theta[0][2][6],errsigma_theta[0][2][7],errsigma_theta[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errsigma_thetaphi[1][1][2],errsigma_thetaphi[1][1][3],errsigma_thetaphi[1][1][4],errsigma_thetaphi[1][1][5],errsigma_thetaphi[1][1][6],errsigma_thetaphi[1][1][7],errsigma_thetaphi[1][2][2],errsigma_thetaphi[1][2][3],errsigma_thetaphi[1][2][4],errsigma_thetaphi[1][2][5],errsigma_thetaphi[1][2][6],errsigma_thetaphi[1][2][7],errsigma_thetaphi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errsigma_phi[1][1][2],errsigma_phi[1][1][3],errsigma_phi[1][1][4],errsigma_phi[1][1][5],errsigma_phi[1][1][6],errsigma_phi[1][1][7],errsigma_phi[1][2][2],errsigma_phi[1][2][3],errsigma_phi[1][2][4],errsigma_phi[1][2][5],errsigma_phi[1][2][6],errsigma_phi[1][2][7],errsigma_phi[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errsigma_theta[1][1][2],errsigma_theta[1][1][3],errsigma_theta[1][1][4],errsigma_theta[1][1][5],errsigma_theta[1][1][6],errsigma_theta[1][1][7],errsigma_theta[1][2][2],errsigma_theta[1][2][3],errsigma_theta[1][2][4],errsigma_theta[1][2][5],errsigma_theta[1][2][6],errsigma_theta[1][2][7],errsigma_theta[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",errsigma_thetaphi[0][1][2],errsigma_thetaphi[0][1][3],errsigma_thetaphi[0][1][4],errsigma_thetaphi[0][1][5],errsigma_thetaphi[0][1][6],errsigma_thetaphi[0][1][7],errsigma_thetaphi[0][2][2],errsigma_thetaphi[0][2][3],errsigma_thetaphi[0][2][4],errsigma_thetaphi[0][2][5],errsigma_thetaphi[0][2][6],errsigma_thetaphi[0][2][7],errsigma_thetaphi[0][2][8]);
			fprintf(outputFile2, "{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},\n");
			fprintf(outputFile2, "{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},\n");
			fprintf(outputFile2, "{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},\n");
			fprintf(outputFile2, "{{0,0,0,0,0,0},{0,0,0,0,0,0,0}}\n");

			fprintf(outputFile2, "\n");

			fprintf(outputFile2, "MissingGenerations(rap1_pt6 CS) = %d",MissingGes[0][1][6]);
			fprintf(outputFile2, "MissingGenerations(rap1_pt6 HX) = %d",MissingGes[1][1][6]);
			fprintf(outputFile2, "MissingGenerations(rap1_pt7 CS) = %d",MissingGes[0][1][7]);
			fprintf(outputFile2, "MissingGenerations(rap1_pt7 HX) = %d",MissingGes[1][1][7]);
			fprintf(outputFile2, "MissingGenerations(rap2_pt6 CS) = %d",MissingGes[0][2][6]);
			fprintf(outputFile2, "MissingGenerations(rap2_pt6 HX) = %d",MissingGes[1][2][6]);
			fprintf(outputFile2, "MissingGenerations(rap2_pt7 CS) = %d",MissingGes[0][2][7]);
			fprintf(outputFile2, "MissingGenerations(rap2_pt7 HX) = %d",MissingGes[1][2][7]);
			fprintf(outputFile2, "MissingGenerations(rap2_pt8 CS) = %d",MissingGes[0][2][8]);
			fprintf(outputFile2, "MissingGenerations(rap2_pt8 HX) = %d",MissingGes[1][2][8]);

		    fclose(outputFile2);

			cout<<" "<<endl;
			//if(missingGenerations>0) cout<<"WARNING:::MISSING_GENERATIONS_TOTAL:::"<<missingGenerations<<endl;

			cout<<"MissingGenerations(rap1_pt6 CS) = "<<MissingGes[0][1][6]<<endl;
			cout<<"MissingGenerations(rap1_pt6 HX) = "<<MissingGes[1][1][6]<<endl;
			cout<<"MissingGenerations(rap1_pt7 CS) = "<<MissingGes[0][1][7]<<endl;
			cout<<"MissingGenerations(rap1_pt7 HX) = "<<MissingGes[1][1][7]<<endl;
			cout<<"MissingGenerations(rap2_pt6 CS) = "<<MissingGes[0][2][6]<<endl;
			cout<<"MissingGenerations(rap2_pt6 HX) = "<<MissingGes[1][2][6]<<endl;
			cout<<"MissingGenerations(rap2_pt7 CS) = "<<MissingGes[0][2][7]<<endl;
			cout<<"MissingGenerations(rap2_pt7 HX) = "<<MissingGes[1][2][7]<<endl;
			cout<<"MissingGenerations(rap2_pt8 CS) = "<<MissingGes[0][2][8]<<endl;
			cout<<"MissingGenerations(rap2_pt8 HX) = "<<MissingGes[1][2][8]<<endl;

		    cout<<"here..."<<endl;



////////////////// shouldn't do that... got no time now!!! /////////////
/*		    mean_phi[0][1][5]=-9999;
		    mean_theta[0][1][5]=-9999;
		    mean_thetaphi[0][1][5]=-9999;
		    mean_phi_val[0][1][5]=-9999;
		    mean_theta_val[0][1][5]=-9999;
		    mean_thetaphi_val[0][1][5]=-9999;
		    mean_phi[1][1][5]=-9999;
		    mean_theta[1][1][5]=-9999;
		    mean_thetaphi[1][1][5]=-9999;
		    mean_phi_val[1][1][5]=-9999;
		    mean_theta_val[1][1][5]=-9999;
		    mean_thetaphi_val[1][1][5]=-9999;
*//////////////////////////////////////////////////////////////////////////

			cout<<"meancheck1"<<mean_thetaphi[0][1][6]<<endl;
			cout<<"meancheck1"<<mean_thetaphi[1][1][6]<<endl;

		    double plotvalues[100][2][12];
		    double errplotvalues[100][2][12];

			    for(int yBin=0; yBin<2; yBin++){
				    for(int ptBin=0; ptBin<10; ptBin++){
				    	plotvalues[0][yBin][ptBin]=mean_phi[0][yBin+1][ptBin+2];
				    	plotvalues[1][yBin][ptBin]=mean_theta[0][yBin+1][ptBin+2];
				    	plotvalues[2][yBin][ptBin]=mean_thetaphi[0][yBin+1][ptBin+2];
				    	plotvalues[3][yBin][ptBin]=mean_phi[1][yBin+1][ptBin+2];
				    	plotvalues[4][yBin][ptBin]=mean_theta[1][yBin+1][ptBin+2];
				    	plotvalues[5][yBin][ptBin]=mean_thetaphi[1][yBin+1][ptBin+2];
				    	plotvalues[6][yBin][ptBin]=sigma_phi[0][yBin+1][ptBin+2];
				    	plotvalues[7][yBin][ptBin]=sigma_theta[0][yBin+1][ptBin+2];
				    	plotvalues[8][yBin][ptBin]=sigma_thetaphi[0][yBin+1][ptBin+2];
				    	plotvalues[9][yBin][ptBin]=sigma_phi[1][yBin+1][ptBin+2];
				    	plotvalues[10][yBin][ptBin]=sigma_theta[1][yBin+1][ptBin+2];
				    	plotvalues[11][yBin][ptBin]=sigma_thetaphi[1][yBin+1][ptBin+2];
				    	plotvalues[12][yBin][ptBin]=convCount[0][yBin+1][ptBin+2];
				    	plotvalues[13][yBin][ptBin]=convCount[1][yBin+1][ptBin+2];
				    	plotvalues[14][yBin][ptBin]=overflow[0][yBin+1][ptBin+2];
				    	plotvalues[15][yBin][ptBin]=overflow[1][yBin+1][ptBin+2];
				    	plotvalues[16][yBin][ptBin]=mean_phi_val[0][yBin+1][ptBin+2];
				    	plotvalues[17][yBin][ptBin]=mean_theta_val[0][yBin+1][ptBin+2];
				    	plotvalues[18][yBin][ptBin]=mean_thetaphi_val[0][yBin+1][ptBin+2];
				    	plotvalues[19][yBin][ptBin]=mean_phi_val[1][yBin+1][ptBin+2];
				    	plotvalues[20][yBin][ptBin]=mean_theta_val[1][yBin+1][ptBin+2];
				    	plotvalues[21][yBin][ptBin]=mean_thetaphi_val[1][yBin+1][ptBin+2];
				    	plotvalues[22][yBin][ptBin]=sigma_phi_val[0][yBin+1][ptBin+2];
				    	plotvalues[23][yBin][ptBin]=sigma_theta_val[0][yBin+1][ptBin+2];
				    	plotvalues[24][yBin][ptBin]=sigma_thetaphi_val[0][yBin+1][ptBin+2];
				    	plotvalues[25][yBin][ptBin]=sigma_phi_val[1][yBin+1][ptBin+2];
				    	plotvalues[26][yBin][ptBin]=sigma_theta_val[1][yBin+1][ptBin+2];
				    	plotvalues[27][yBin][ptBin]=sigma_thetaphi_val[1][yBin+1][ptBin+2];
				    	plotvalues[28][yBin][ptBin]=mean_phi_val[0][yBin+1][ptBin+2]-promptlambda_phi_injected_CS[yBin+1][ptBin+2];
				    	plotvalues[29][yBin][ptBin]=mean_theta_val[0][yBin+1][ptBin+2]-promptlambda_theta_injected_CS[yBin+1][ptBin+2];
				    	plotvalues[30][yBin][ptBin]=mean_thetaphi_val[0][yBin+1][ptBin+2]-promptlambda_thetaphi_injected_CS[yBin+1][ptBin+2];
				    	plotvalues[31][yBin][ptBin]=mean_phi_val[1][yBin+1][ptBin+2]-promptlambda_phi_injected_HX[yBin+1][ptBin+2];
				    	plotvalues[32][yBin][ptBin]=mean_theta_val[1][yBin+1][ptBin+2]-promptlambda_theta_injected_HX[yBin+1][ptBin+2];
				    	plotvalues[33][yBin][ptBin]=mean_thetaphi_val[1][yBin+1][ptBin+2]-promptlambda_thetaphi_injected_HX[yBin+1][ptBin+2];
				    	plotvalues[34][yBin][ptBin]=MissingGes[0][yBin+1][ptBin+2];
				    	plotvalues[35][yBin][ptBin]=MissingGes[1][yBin+1][ptBin+2];
				    	plotvalues[36][yBin][ptBin]=MissingGes[0][yBin+1][ptBin+2]+convCount[0][yBin+1][ptBin+2];
				    	plotvalues[37][yBin][ptBin]=MissingGes[1][yBin+1][ptBin+2]+convCount[1][yBin+1][ptBin+2];

				    	errplotvalues[0][yBin][ptBin]=errmean_phi[0][yBin+1][ptBin+2];
				    	errplotvalues[1][yBin][ptBin]=errmean_theta[0][yBin+1][ptBin+2];
				    	errplotvalues[2][yBin][ptBin]=errmean_thetaphi[0][yBin+1][ptBin+2];
				    	errplotvalues[3][yBin][ptBin]=errmean_phi[1][yBin+1][ptBin+2];
				    	errplotvalues[4][yBin][ptBin]=errmean_theta[1][yBin+1][ptBin+2];
				    	errplotvalues[5][yBin][ptBin]=errmean_thetaphi[1][yBin+1][ptBin+2];
				    	errplotvalues[6][yBin][ptBin]=errsigma_phi[0][yBin+1][ptBin+2];
				    	errplotvalues[7][yBin][ptBin]=errsigma_theta[0][yBin+1][ptBin+2];
				    	errplotvalues[8][yBin][ptBin]=errsigma_thetaphi[0][yBin+1][ptBin+2];
				    	errplotvalues[9][yBin][ptBin]=errsigma_phi[1][yBin+1][ptBin+2];
				    	errplotvalues[10][yBin][ptBin]=errsigma_theta[1][yBin+1][ptBin+2];
				    	errplotvalues[11][yBin][ptBin]=errsigma_thetaphi[1][yBin+1][ptBin+2];
				    	errplotvalues[12][yBin][ptBin]=0;
				    	errplotvalues[13][yBin][ptBin]=0;
				    	errplotvalues[14][yBin][ptBin]=0;
				    	errplotvalues[15][yBin][ptBin]=0;
				    	errplotvalues[16][yBin][ptBin]=sigma_phi_val[0][yBin+1][ptBin+2];
				    	errplotvalues[17][yBin][ptBin]=sigma_theta_val[0][yBin+1][ptBin+2];
				    	errplotvalues[18][yBin][ptBin]=sigma_thetaphi_val[0][yBin+1][ptBin+2];
				    	errplotvalues[19][yBin][ptBin]=sigma_phi_val[1][yBin+1][ptBin+2];
				    	errplotvalues[20][yBin][ptBin]=sigma_theta_val[1][yBin+1][ptBin+2];
				    	errplotvalues[21][yBin][ptBin]=sigma_thetaphi_val[1][yBin+1][ptBin+2];
				    	errplotvalues[22][yBin][ptBin]=errsigma_phi_val[0][yBin+1][ptBin+2];
				    	errplotvalues[23][yBin][ptBin]=errsigma_theta_val[0][yBin+1][ptBin+2];
				    	errplotvalues[24][yBin][ptBin]=errsigma_thetaphi_val[0][yBin+1][ptBin+2];
				    	errplotvalues[25][yBin][ptBin]=errsigma_phi_val[1][yBin+1][ptBin+2];
				    	errplotvalues[26][yBin][ptBin]=errsigma_theta_val[1][yBin+1][ptBin+2];
				    	errplotvalues[27][yBin][ptBin]=errsigma_thetaphi_val[1][yBin+1][ptBin+2];
				    	errplotvalues[28][yBin][ptBin]=errmean_phi_val[0][yBin+1][ptBin+2];//sigma_phi_val[0][yBin+1][ptBin+2];
				    	errplotvalues[29][yBin][ptBin]=errmean_theta_val[0][yBin+1][ptBin+2];//sigma_theta_val[0][yBin+1][ptBin+2];
				    	errplotvalues[30][yBin][ptBin]=errmean_thetaphi_val[0][yBin+1][ptBin+2];//sigma_thetaphi_val[0][yBin+1][ptBin+2];
				    	errplotvalues[31][yBin][ptBin]=errmean_phi_val[1][yBin+1][ptBin+2];//sigma_phi_val[1][yBin+1][ptBin+2];
				    	errplotvalues[32][yBin][ptBin]=errmean_theta_val[1][yBin+1][ptBin+2];//sigma_theta_val[1][yBin+1][ptBin+2];
				    	errplotvalues[33][yBin][ptBin]=errmean_thetaphi_val[1][yBin+1][ptBin+2];//sigma_thetaphi_val[1][yBin+1][ptBin+2];
				    	errplotvalues[34][yBin][ptBin]=0;
				    	errplotvalues[35][yBin][ptBin]=0;
				    	errplotvalues[36][yBin][ptBin]=0;
				    	errplotvalues[37][yBin][ptBin]=0;
				    }
		    }


	//		char dirStruct[200];
	//		sprintf(dirStruct,"Plots/ParameterPlots/multiGen000");

			gSystem->mkdir(dirStruct);

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

			double plotvalues_[jpsi::kNbPTMaxBins];
			double errplotvalues_[jpsi::kNbPTMaxBins];
			double ptCentre[jpsi::kNbRapForPTBins][jpsi::kNbPTMaxBins]={
					{3,6.5,7.5,9,12.5,17.5,25},
					{2,5,6.5,7.5,9,12.5,17.5,25}
			};
			double ptMean[jpsi::kNbPTMaxBins]={-9999,-9999,-9999,-9999,-9999,-9999};
			double errptMean[jpsi::kNbPTMaxBins]={NULL};

			char legendrap[200];
			char Filename[200];
			char yAxis[200];
			double yMin;
			double yMax;

		//	double plotvalues[2][12];
		//	double errplotvalues[2][12];

			for(int i=1;i<39;i++){


		////////////////////////////// FILL INS /////////////////////////////
/*				double plotvalues[22][2][12]={
						{{-0.07240,0.11677,-0.14533,-0.29752,-0.04619,-9999.00000},{0.09032,0.11540,-0.02027,-0.07337,-0.13592,-0.10115,0.09064}},
						{{0.11646,-0.22619,-0.10316,0.14406,0.04314,-9999.00000},{0.00101,0.04654,0.34408,0.20649,0.24520,0.15611,0.18207}},
						{{-0.01986,0.01301,-0.01714,0.04978,-0.13467,-9999.00000},{0.03394,-0.12039,0.08061,-0.10259,0.11117,0.01947,0.13219}},
						{{0.16535,0.03471,0.09044,-0.19271,-0.31755,-9999.00000},{0.06379,0.14039,0.21751,0.17525,0.13289,-0.01600,0.03353}},
						{{-0.13676,0.11479,-0.16688,0.29657,-0.20683,-9999.00000},{0.05664,-0.18328,-0.05448,-0.00938,0.00961,-0.08935,0.00800}},
						{{-0.06855,0.21371,0.17602,-0.14533,-0.14850,-9999.00000},{0.12266,-0.02246,-0.04664,-0.00511,-0.04197,0.04235,0.08434}},
						{{0.84285,1.12189,0.87579,1.06521,1.06482,-9999.00000},{0.87995,0.91128,1.03967,0.89443,0.91656,1.01205,0.92472}},
						{{0.95640,0.92226,1.03505,0.89825,0.96664,-9999.00000},{1.03002,1.06554,0.92142,0.92426,0.99351,1.08098,1.00305}},
						{{1.05581,0.98838,1.09055,0.69054,0.90197,-9999.00000},{1.04265,1.08439,1.09067,1.09411,1.07211,1.01879,0.85151}},
						{{0.94977,0.99739,0.87443,1.04632,0.99787,-9999.00000},{1.02570,1.11737,1.16007,1.22825,1.06840,1.11417,0.88095}},
						{{0.94789,1.01579,1.00604,0.88485,1.11341,-9999.00000},{0.93295,1.03066,1.04311,0.86823,1.05470,1.04909,1.01993}},
						{{1.16568,0.89397,0.99386,1.20574,0.95077,-9999.00000},{1.01179,0.91177,0.94019,0.88844,0.83012,1.04642,1.09326}},
						{{50,50,50,50,50,0},{50,50,50,50,50,50,50}},
						{{50,50,50,50,50,0},{50,50,50,50,50,50,50}},
						{{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000}},
						{{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},{0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000}}

			};

			double errplotvalues[22][2][12]={


					{{0.11993,0.14486,0.13455,0.12094,0.09807,0.11071},{0.16557,0.17942,0.11901,0.09310,0.09926,0.11261,0.10169}},
					 {{0.11919,0.15864,0.12385,0.15063,0.15057,0.00000},{0.12443,0.12886,0.14702,0.12648,0.12961,0.14311,0.13076}},
									{{0.13524,0.13042,0.14636,0.12702,0.13669,0.00000},{0.14565,0.15068,0.13029,0.13070,0.14049,0.15285,0.14184}},
									{{0.14930,0.13976,0.15420,0.09765,0.12755,0.00000},{0.14744,0.15333,0.15422,0.15471,0.15160,0.14406,0.12041}},
									{{0.13431,0.14104,0.12365,0.14795,0.14111,0.00000},{0.14504,0.15800,0.16404,0.17367,0.15108,0.15755,0.12457}},
									{{0.13404,0.14364,0.14226,0.12513,0.15744,0.00000},{0.13193,0.14575,0.14750,0.12278,0.14914,0.14835,0.14423}},
									{{0.16483,0.12642,0.14054,0.17049,0.13445,0.00000},{0.14307,0.12893,0.13295,0.12563,0.11739,0.14797,0.15459}},
									{{0.08424,0.11212,0.08754,0.10643,0.10633,0.00000},{0.08796,0.09104,0.10389,0.08943,0.09162,0.10115,0.09244}},
									{{0.09562,0.09221,0.10341,0.08981,0.09663,0.00000},{0.10296,0.10656,0.09212,0.09243,0.09932,0.10807,0.10028}},
									{{0.11652,0.08937,0.09936,0.12053,0.09501,0.00000},{0.10108,0.09116,0.09399,0.08882,0.08301,0.10462,0.10920}},
									{{0.09493,0.09963,0.08743,0.10459,0.09971,0.00000},{0.10255,0.11171,0.11598,0.12274,0.10682,0.11133,0.08804}},
									{{0.09477,0.10151,0.10056,0.08847,0.11128,0.00000},{0.09329,0.10301,0.10424,0.08676,0.10545,0.10483,0.10196}},
									{{0.10556,0.09882,0.10892,0.06903,0.09020,0.00000},{0.10418,0.10842,0.10904,0.10934,0.10712,0.10186,0.08513}},
									{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
									{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
									{{0,0,0,0,0,0},{0,0,0,0,0,0,0}},
									{{0,0,0,0,0,0},{0,0,0,0,0,0,0}}




			};

*/


				borders_val=1;
		/////// mean /////
			if(i<6.5) {yMin=-1.5;yMax=1.5;}
		////// sigma /////
			if(i>6.5 && i<12.5) {yMin=0;yMax=3;}
		////// Conv counter ////
			if(i>12.5 && i<14.5) {yMin=1;yMax=generations*1.3;}
		////// Overflow counter ////
			if(i>14.5 && i<16.5) {yMin=0;yMax=10;}
		////// mean val ////
			if(i==17) {yMin=promptlambda_phi_injected_CS[1][6]-borders_val;yMax=promptlambda_phi_injected_CS[1][6]+borders_val;}            //{yMin=-1.5;yMax=1.5;}
			if(i==18) {yMin=promptlambda_theta_injected_CS[1][6]-borders_val;yMax=promptlambda_theta_injected_CS[1][6]+borders_val;}        //{yMin=-1.5;yMax=1.5;}
			if(i==19) {yMin=promptlambda_thetaphi_injected_CS[1][6]-borders_val;yMax=promptlambda_thetaphi_injected_CS[1][6]+borders_val;}  //{yMin=-1.5;yMax=1.5;}
			if(i==20) {yMin=promptlambda_phi_injected_HX[1][6]-borders_val;yMax=promptlambda_phi_injected_HX[1][6]+borders_val;}            //{yMin=-1.5;yMax=1.5;}
			if(i==21) {yMin=promptlambda_theta_injected_HX[1][6]-borders_val;yMax=promptlambda_theta_injected_HX[1][6]+borders_val;}        //{yMin=-1.5;yMax=1.5;}
			if(i==22) {yMin=promptlambda_thetaphi_injected_HX[1][6]-borders_val;yMax=promptlambda_thetaphi_injected_HX[1][6]+borders_val;}  //{yMin=-1.5;yMax=1.5;}
		////// sigma val ////
			if(i>22 && i<29) {yMin=0;yMax=0.5;}
		////// delta mean val ////
			if(i>28.5 && i<34.5) {yMin=-0.5;yMax=0.5;}
		////// MissingGen counter ////
			if(i>34.5 && i<36.5) {yMin=1;yMax=generations*1.3;}
		////// MissingGen + conv counter ////
			if(i>36.5 && i<38.5) {yMin=1;yMax=generations*1.3;}

			if(i==1) if(!tilde) sprintf(yAxis,"#mu[z(#lambda_{#phi CS})]"); else sprintf(yAxis,"#mu[z(#tilde{#lambda}_{CS})]");
			if(i==2) sprintf(yAxis,"#mu[z(#lambda_{#theta CS})]");
			if(i==3) sprintf(yAxis,"#mu[z(#lambda_{#theta #phi CS})]");
			if(i==4) if(!tilde) sprintf(yAxis,"#mu[z(#lambda_{#phi HX})]"); else sprintf(yAxis,"#mu[z(#tilde{#lambda}_{HX})]");
			if(i==5) sprintf(yAxis,"#mu[z(#lambda_{#theta HX})]");
			if(i==6) sprintf(yAxis,"#mu[z(#lambda_{#theta #phi HX})]");
			if(i==7) sprintf(yAxis,"#sigma[z(#lambda_{#phi CS})]");
			if(i==8) sprintf(yAxis,"#sigma[z(#lambda_{#theta CS})]");
			if(i==9) sprintf(yAxis,"#sigma[z(#lambda_{#theta #phi CS})]");
			if(i==10) sprintf(yAxis,"#sigma[z(#lambda_{#phi HX})]");
			if(i==11) sprintf(yAxis,"#sigma[z(#lambda_{#theta HX})]");
			if(i==12) sprintf(yAxis,"#sigma[z(#lambda_{#theta #phi HX})]");
			if(i==13) sprintf(yAxis,"Converged fits CS");
			if(i==14) sprintf(yAxis,"Converged fits HX");
			if(i==15) sprintf(yAxis,"Overflow CS");
			if(i==16) sprintf(yAxis,"Overflow HX");
			if(i==17) if(!tilde) sprintf(yAxis,"#mu[#lambda_{#phi CS}]"); else sprintf(yAxis,"#mu[#tilde{#lambda}_{CS}]");
			if(i==18) sprintf(yAxis,"#mu[#lambda_{#theta CS}]");
			if(i==19) sprintf(yAxis,"#mu[#lambda_{#theta #phi CS}]");
			if(i==20) if(!tilde) sprintf(yAxis,"#mu[#lambda_{#phi HX}]"); else sprintf(yAxis,"#mu[#tilde{#lambda}_{HX}]");
			if(i==21) sprintf(yAxis,"#mu[#lambda_{#theta HX}]");
			if(i==22) sprintf(yAxis,"#mu[#lambda_{#theta #phi HX}]");
			if(i==23) if(!tilde) sprintf(yAxis,"#sigma[#lambda_{#phi CS}]"); else sprintf(yAxis,"#sigma[#tilde{#lambda}_{CS}]");
			if(i==24) sprintf(yAxis,"#sigma[#lambda_{#theta CS}]");
			if(i==25) sprintf(yAxis,"#sigma[#lambda_{#theta #phi CS}]");
			if(i==26) if(!tilde) sprintf(yAxis,"#sigma[#lambda_{#phi HX}]"); else sprintf(yAxis,"#sigma[#tilde{#lambda}_{HX}]");
			if(i==27) sprintf(yAxis,"#sigma[#lambda_{#theta HX}]");
			if(i==28) sprintf(yAxis,"#sigma[#lambda_{#theta #phi HX}]");

			if(i==29) if(!tilde) sprintf(yAxis,"#Delta#lambda_{#phi CS}"); else sprintf(yAxis,"#Delta#tilde{#lambda}_{CS}");
			if(i==30) sprintf(yAxis,"#Delta#lambda_{#theta CS}");
			if(i==31) sprintf(yAxis,"#Delta#lambda_{#theta #phi CS}");
			if(i==32) if(!tilde) sprintf(yAxis,"#Delta#lambda_{#phi HX}"); else sprintf(yAxis,"#Delta#tilde{#lambda}_{HX}");
			if(i==33) sprintf(yAxis,"#Delta#lambda_{#theta HX}");
			if(i==34) sprintf(yAxis,"#Delta#lambda_{#theta #phi HX}");
			if(i==35) sprintf(yAxis,"Missing Generations CS");
			if(i==36) sprintf(yAxis,"Missing Generations HX");
			if(i==37) sprintf(yAxis,"Missing Generations + converged CS");
			if(i==38) sprintf(yAxis,"Missing Generations + converged HX");


			if(i==1) if(!tilde) sprintf(Filename,"%s/pT_ToyMCmean_phi_CS.png",dirStruct);  else sprintf(Filename,"%s/pT_ToyMCmean_tilde_CS.png",dirStruct);
			if(i==2) sprintf(Filename,"%s/pT_ToyMCmean_theta_CS.png",dirStruct);
			if(i==3) sprintf(Filename,"%s/pT_ToyMCmean_thetaphi_CS.png",dirStruct);
			if(i==4) if(!tilde) sprintf(Filename,"%s/pT_ToyMCmean_phi_HX.png",dirStruct);  else sprintf(Filename,"%s/pT_ToyMCmean_tilde_HX.png",dirStruct);
			if(i==5) sprintf(Filename,"%s/pT_ToyMCmean_theta_HX.png",dirStruct);
			if(i==6) sprintf(Filename,"%s/pT_ToyMCmean_thetaphi_HX.png",dirStruct);
			if(i==7) if(!tilde) sprintf(Filename,"%s/pT_ToyMCsigma_phi_CS.png",dirStruct); else sprintf(Filename,"%s/pT_ToyMCsigma_tilde_CS.png",dirStruct);
			if(i==8) sprintf(Filename,"%s/pT_ToyMCsigma_theta_CS.png",dirStruct);
			if(i==9) sprintf(Filename,"%s/pT_ToyMCsigma_thetaphi_CS.png",dirStruct);
			if(i==10) if(!tilde) sprintf(Filename,"%s/pT_ToyMCsigma_phi_HX.png",dirStruct);  else sprintf(Filename,"%s/pT_ToyMCsigma_tilde_HX.png",dirStruct);
			if(i==11) sprintf(Filename,"%s/pT_ToyMCsigma_theta_HX.png",dirStruct);
			if(i==12) sprintf(Filename,"%s/pT_ToyMCsigma_thetaphi_HX.png",dirStruct);
			if(i==13) sprintf(Filename,"%s/pT_ToyMCconvCS.png",dirStruct);
			if(i==14) sprintf(Filename,"%s/pT_ToyMCconvHX.png",dirStruct);
			if(i==15) sprintf(Filename,"%s/pT_ToyMCoverflowCS.png",dirStruct);
			if(i==16) sprintf(Filename,"%s/pT_ToyMCoverflowHX.png",dirStruct);
			if(i==17) if(!tilde) sprintf(Filename,"%s/pT_ToyMCval_mean_phi_CS.png",dirStruct);  else sprintf(Filename,"%s/pT_ToyMCval_mean_tilde_CS.png",dirStruct);
			if(i==18) sprintf(Filename,"%s/pT_ToyMCval_mean_theta_CS.png",dirStruct);
			if(i==19) sprintf(Filename,"%s/pT_ToyMCval_mean_thetaphi_CS.png",dirStruct);
			if(i==20) if(!tilde) sprintf(Filename,"%s/pT_ToyMCval_mean_phi_HX.png",dirStruct); else sprintf(Filename,"%s/pT_ToyMCval_mean_tilde_HX.png",dirStruct);
			if(i==21) sprintf(Filename,"%s/pT_ToyMCval_mean_theta_HX.png",dirStruct);
			if(i==22) sprintf(Filename,"%s/pT_ToyMCval_mean_thetaphi_HX.png",dirStruct);
			if(i==23) if(!tilde) sprintf(Filename,"%s/pT_ToyMCval_sigma_phi_CS.png",dirStruct); else sprintf(Filename,"%s/pT_ToyMCval_sigma_tilde_CS.png",dirStruct);
			if(i==24) sprintf(Filename,"%s/pT_ToyMCval_sigma_theta_CS.png",dirStruct);
			if(i==25) sprintf(Filename,"%s/pT_ToyMCval_sigma_thetaphi_CS.png",dirStruct);
			if(i==26) if(!tilde) sprintf(Filename,"%s/pT_ToyMCval_sigma_phi_HX.png",dirStruct); else sprintf(Filename,"%s/pT_ToyMCval_sigma_tilde_HX.png",dirStruct);
			if(i==27) sprintf(Filename,"%s/pT_ToyMCval_sigma_theta_HX.png",dirStruct);
			if(i==28) sprintf(Filename,"%s/pT_ToyMCval_sigma_thetaphi_HX.png",dirStruct);
			if(i==29) if(!tilde) sprintf(Filename,"%s/pT_ToyMCval_deltamean_phi_CS.png",dirStruct); else sprintf(Filename,"%s/pT_ToyMCval_deltamean_tilde_CS.png",dirStruct);
			if(i==30) sprintf(Filename,"%s/pT_ToyMCval_deltamean_theta_CS.png",dirStruct);
			if(i==31) sprintf(Filename,"%s/pT_ToyMCval_deltamean_thetaphi_CS.png",dirStruct);
			if(i==32) if(!tilde) sprintf(Filename,"%s/pT_ToyMCval_deltamean_phi_HX.png",dirStruct); else sprintf(Filename,"%s/pT_ToyMCval_deltamean_tilde_HX.png",dirStruct);
			if(i==33) sprintf(Filename,"%s/pT_ToyMCval_deltamean_theta_HX.png",dirStruct);
			if(i==34) sprintf(Filename,"%s/pT_ToyMCval_deltamean_thetaphi_HX.png",dirStruct);
			if(i==35) sprintf(Filename,"%s/pT_MissingGenerations_CS.png",dirStruct);
			if(i==36) sprintf(Filename,"%s/pT_MissingGenerations_HX.png",dirStruct);
			if(i==37) sprintf(Filename,"%s/pT_MissingGenerations_andconverved_CS.png",dirStruct);
			if(i==38) sprintf(Filename,"%s/pT_MissingGenerations_andconverved_HX.png",dirStruct);

		/////////////////////////////////////////////////////////////////////////////////////////////////////

			TLegend* plotLegend=new TLegend(0.7125,0.75,0.95,0.9);
			plotLegend->SetFillColor(kWhite);
			plotLegend->SetTextFont(72);
			plotLegend->SetTextSize(0.035);
			plotLegend->SetBorderSize(1);


			TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

			plotCanvas->SetFillColor(kWhite);
			plotCanvas->SetGrid();
			plotCanvas->GetFrame()->SetFillColor(kWhite);
			plotCanvas->GetFrame()->SetBorderSize(0);
			plotCanvas->SetRightMargin(0.05) ;

			TH1F *plotHisto = new TH1F;
		//	TH1F *
			plotHisto = plotCanvas->DrawFrame(10,yMin,30,yMax);
			plotHisto->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
			plotHisto->SetYTitle(yAxis);
			plotHisto->GetYaxis()->SetTitleOffset(1.5);

			for(int yBin = 0; yBin < 2; ++yBin) {
			    	for(int ptBin = 0; ptBin < 10; ++ptBin) {

				    				   ptMean[ptBin]=ptCentre[yBin][ptBin]+0.25*yBin;
			  					   errptMean[ptBin]=0;
			  					 plotvalues_[ptBin+1]=plotvalues[i-1][yBin][ptBin];
			  					errplotvalues_[ptBin+1]=errplotvalues[i-1][yBin][ptBin];
			  				//	cout<<"i="<<i<<" "<<plotvalues[i-1][yBin][ptBin]<<endl;
			  					  }

			TGraphErrors *plotGraph = new TGraphErrors(jpsi::kNbPTBins[yBin]+1,ptMean,plotvalues_,errptMean,errplotvalues_);
			plotGraph->SetMarkerColor(jpsi::colour_rapForPTBins[yBin+1]);
		//	plotGraph->SetLineColor(jpsi::colour_rapForPTBins[yBin+1]);
			plotGraph->SetMarkerStyle(jpsi::marker_rapForPTBins[yBin+1]);
			if(yBin==0) sprintf(legendrap,"|y| < 0.9");
			if(yBin==1) sprintf(legendrap,"0.9 < |y| < 1.2");
			plotLegend->AddEntry(plotGraph,legendrap,"ple");
			plotGraph->Draw("P");
			//delete plotGraph;

			  		}




			  		plotLegend->Draw();

			  		plotCanvas->SaveAs(Filename);

			  		plotCanvas->Close();

			  		delete plotLegend;
			  		delete plotCanvas;
		//	  		delete plotHisto;
			}






  return 0;
}


