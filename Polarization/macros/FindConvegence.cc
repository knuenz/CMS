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

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"

int main(int argc, char** argv) {
//	using namespace JPsiPolarization;
	using namespace std;
	using namespace RooFit;

	bool doFit(false);

	  gSystem->mkdir("Plots/ToyMC");
	  gSystem->mkdir("Results/ToyMC");

	  gStyle->SetTitleFillColor(kWhite);

	  char line[200];

	  char outputfilename2[200];
	  sprintf(outputfilename2,"Results/ToyMC/parametersToyMCTotal.txt");
	  FILE *outputFile2 = fopen(outputfilename2,"w");

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

	  double overflow[2][jpsi::kNbRapForPTBins+1][jpsi::kNbPTMaxBins+1]={NULL};

	    int convCountCheck[2][6][13]={NULL};


	    	for(int yBinstart = 1; yBinstart < 2; yBinstart++) {
	    		for(int ptBinstart = 2; ptBinstart < 3; ptBinstart++) {
//3,9 for mid, intermediate bins
	    			 		  if(yBinstart==1 && ptBinstart==8) continue;

//	  int yBinstart=1;
//	  int ptBinstart=2;

		char inputfilename[200];
//		sprintf(inputfilename,"terminalToyMC_rap%d_pt%d.txt",yBinstart,ptBinstart);
		sprintf(inputfilename,"terminalToyMC_rap1_pt2.txt");

		FILE *inputFile = fopen(inputfilename,"r");

		char outputfilename[200];
		sprintf(outputfilename,"Results/ToyMC/parametersToyMC_rap%d_pt%d.txt",yBinstart,ptBinstart);
//		sprintf(outputfilename,"Results/ToyMC/parametersRealMC.txt");
		FILE *outputFile = fopen(outputfilename,"w");



		  if (inputFile == NULL)
		  {printf("\nQuelle \"%s\" falsch - kann Datei nicht šffnen\n", inputfilename);return -1;}

		  int ptBin=0;
		  int yBin;

		  bool RealMC(true);

		  bool isCS;
		  bool isHX;
		  bool converged(false);
		  bool errormatrixuncertainty(false);
		  bool minos(true);
		  int convCount[2][6][13]={NULL};
		  int CONV[2][6][13][100]={NULL};
		  int generation=0;
		  char zeile[200];
		  int i=0;
		  int fitresIndex=1000;
		  int MinosresIndex=1000;

		  double edmCondition=0.001;

		  double edm;
		  double nPrompt[2][6][13][100];
		  double errnPrompt[2][6][13][100];
		  double promptlambda_phi[2][6][13][100];
		  double errpromptlambda_phi[2][6][13][100];
		  double promptlambda_theta[2][6][13][100];
		  double errpromptlambda_theta[2][6][13][100];
		  double promptlambda_thetaphi[2][6][13][100];
		  double errpromptlambda_thetaphi[2][6][13][100];

//////// Start rolling over the file, row by row //////////////////

	 do{
	  fgets(zeile, 200, inputFile);


	  char* str1 = zeile;
	  char* str2 = "JpsiPt >";
	  char* str3 = "EDM=";
	  char* str4 = "Fitting CS";
	  char* str5 = "Fitting HX";
	  char* str6 = "abs(JpsiRap) > 0";
	  char* str7 = "abs(JpsiRap) > 0.9";
	  char* str8 = "abs(JpsiRap) > 1.2";
	  char* str9 = "abs(JpsiRap) > 1.6";
	  char* str10 = "abs(JpsiRap) > 2.1";
	  char* str11 = "JpsiPt > 0";
	  char* str12 = "JpsiPt > 1";
	  char* str13 = "JpsiPt > 2";
	  char* str14 = "JpsiPt > 3";
	  char* str15 = "JpsiPt > 4";
	  char* str16 = "JpsiPt > 5";
	  char* str17 = "JpsiPt > 6";
	  char* str18 = "JpsiPt > 7";
	  char* str19 = "JpsiPt > 8";
	  char* str20 = "JpsiPt > 10";
	  char* str21 = "JpsiPt > 15";
	  char* str22 = "JpsiPt > 20";
	  char* str23 = "JpsiPt > 30";
	  char* str24 = "MIGRAD FAILS TO FIND IMPROVEMENT";
	  char* str25 = "GENERATION";
	  char* str26 = "STATUS=";
	  char* str27 = "FROM MIGRAD    STATUS=CONVERGED";
	  char* str28 = "ERROR MATRIX ACCURATE";
	  char* str29 = "Floating Parameter    FinalValue +/-  Error";
	  char* str30 = "ERROR MATRIX UNCERTAINTY 100.0 per cent";
	  char* str31 = "FROM MINOS     STATUS=SUCCESSFUL";
	  char* str32 = "FROM MIGRAD    STATUS=FAILED";


	  char* result2 = strstr( str1, str2 );
	  char* result3 = strstr( str1, str3 );
	  char* result4 = strstr( str1, str4 );
	  char* result5 = strstr( str1, str5 );
	  char* result6 = strstr( str1, str6 );
	  char* result7 = strstr( str1, str7 );
	  char* result8 = strstr( str1, str8 );
	  char* result9 = strstr( str1, str9 );
	  char* result10 = strstr( str1, str10 );
	  char* result11 = strstr( str1, str11 );
	  char* result12 = strstr( str1, str12 );
	  char* result13 = strstr( str1, str13 );
	  char* result14 = strstr( str1, str14 );
	  char* result15 = strstr( str1, str15 );
	  char* result16 = strstr( str1, str16 );
	  char* result17 = strstr( str1, str17 );
	  char* result18 = strstr( str1, str18 );
	  char* result19 = strstr( str1, str19 );
	  char* result20 = strstr( str1, str20 );
	  char* result21 = strstr( str1, str21 );
	  char* result22 = strstr( str1, str22 );
	  char* result23 = strstr( str1, str23 );
	  char* result24 = strstr( str1, str24 );
	  char* result25 = strstr( str1, str25 );
	  char* result26 = strstr( str1, str26 );
	  char* result27 = strstr( str1, str27 );
	  char* result28 = strstr( str1, str28 );
	  char* result29 = strstr( str1, str29 );
	  char* result30 = strstr( str1, str30 );
	  char* result31 = strstr( str1, str31 );
	  char* result32 = strstr( str1, str32 );

	  if(result6!=0) yBin=1;
	  if(result7!=0) yBin=2;
	  if(result8!=0) yBin=3;
	  if(result9!=0) yBin=4;
	  if(result10!=0) yBin=5;

	  if(result11!=0) ptBin=1;

	  if(yBin==1 && result17!=0) ptBin=2; if(yBin==2 | yBin==3 && result15!=0) ptBin=2; if(yBin==4 | yBin==5 && result12!=0) ptBin=2;

	  if(yBin==1 && result18!=0) ptBin=3; if(yBin==2 && result17!=0) ptBin=3; if(yBin==3 && result16!=0) ptBin=3; if(yBin==4 | yBin==5 && result13!=0) ptBin=3;

	  if(yBin==1 && result19!=0) ptBin=4; if(yBin==2 && result18!=0) ptBin=4; if(yBin==3 && result17!=0) ptBin=4; if(yBin==4 | yBin==5 && result14!=0) ptBin=4;

	  if(yBin==1 && result20!=0) ptBin=5; if(yBin==2 && result19!=0) ptBin=5; if(yBin==3 && result18!=0) ptBin=5; if(yBin==4 | yBin==5 && result15!=0) ptBin=5;

	  if(yBin==1 && result21!=0) ptBin=6; if(yBin==2 && result20!=0) ptBin=6; if(yBin==3 && result19!=0) ptBin=6; if(yBin==4 | yBin==5 && result16!=0) ptBin=6;

	  if(yBin==1 && result22!=0) ptBin=7; if(yBin==2 && result21!=0) ptBin=7; if(yBin==3 && result20!=0) ptBin=7; if(yBin==4 | yBin==5 && result17!=0) ptBin=7;

	  if(yBin==2 && result22!=0) ptBin=8; if(yBin==3 && result21!=0) ptBin=8; if(yBin==4 | yBin==5 && result18!=0) ptBin=8;

	  if(yBin==3 && result22!=0) ptBin=9; if(yBin==4 | yBin==5 && result19!=0) ptBin=9;

	  if(yBin==4 | yBin==5 && result20!=0) ptBin=10;

	  if(yBin==4 | yBin==5 && result21!=0) ptBin=11;

	  if(yBin==4 | yBin==5 && result22!=0) ptBin=12;

//	  if(result4!=0 | result5!=0 | result2!=0) {cout<<zeile<<endl;}
	  if(result25!=0) { generation++; cout<<"Generation "<<generation<<endl; }//cout<<result25<<endl; }
	  if(result2!=0) { cout<<"rap"<<yBin<<"_pT"<<ptBin<<endl;}// cout<<NonConvCount<<endl; }

	  if(result4!=0) {isCS=true;isHX=false; cout<<"CS"<<endl; }
	  if(result5!=0) {isCS=false;isHX=true; cout<<"HX"<<endl; }

//	  if(result4!=0 | result5!=0) {CONV[0][yBin][ptBin][generation]=0; CONV[1][yBin][ptBin][generation]=0; }

/////////////// Now you know the exact bin, and if you fit CS or HX and generation /////////////////////////
///////////////////// by ptBin , yBin and isCS, isHX ////////////////////////////////////////

///////////////// convergence criteria ////////////////////////////////

	  edm=1000;

	  if(result3!=0){
		 // cout<<zeile<<endl;
		  char* part = strtok (zeile,"=");
		  part = strtok (NULL, "=");
		  edm = atof(part);
		  }

	  if(result30!=0){errormatrixuncertainty=true;}
	  if(converged/* && edm<=edmCondition && !errormatrixuncertainty*/) { cout<<"converged"<<endl;if(isCS) { convCount[0][yBin][ptBin]++; CONV[0][yBin][ptBin][generation]=1;} if(isHX) {convCount[1][yBin][ptBin]++; CONV[1][yBin][ptBin][generation]=1;}}
	  converged=false;
//	  if(result27!=0) {converged=true; }
	  if(result32!=0) {converged=true; } /// FAIL...

	  errormatrixuncertainty=false;



//////////////// read out results /////////////////////////////////////////

	  if(result29!=0) fitresIndex=0;
	  if(fitresIndex==2) {
			  char* part = strtok (zeile," ");
			  part = strtok (NULL, " ");
			  if(isCS) nPrompt[0][yBin][ptBin][generation] = atof(part);
			  if(isHX) nPrompt[1][yBin][ptBin][generation] = atof(part);
			  part = strtok (NULL, " ");
			  part = strtok (NULL, " ");
			  if(isCS) errnPrompt[0][yBin][ptBin][generation] = atof(part);
			  if(isHX) errnPrompt[1][yBin][ptBin][generation] = atof(part);
			  }
	  if(fitresIndex==3) {
	  			  char* part = strtok (zeile," ");
	  			  part = strtok (NULL, " ");
	  			  if(isCS) promptlambda_phi[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) promptlambda_phi[1][yBin][ptBin][generation] = atof(part);
				  part = strtok (NULL, " ");
				  part = strtok (NULL, " ");
				  if(isCS) errpromptlambda_phi[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) errpromptlambda_phi[1][yBin][ptBin][generation] = atof(part);
	  			  }
	  if(fitresIndex==4) {
	  			  char* part = strtok (zeile," ");
	  			  part = strtok (NULL, " ");
	  			  if(isCS) promptlambda_theta[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) promptlambda_theta[1][yBin][ptBin][generation] = atof(part);
				  part = strtok (NULL, " ");
				  part = strtok (NULL, " ");
				  if(isCS) errpromptlambda_theta[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) errpromptlambda_theta[1][yBin][ptBin][generation] = atof(part);
	  			  }
	  if(fitresIndex==5) {
	  			  char* part = strtok (zeile," ");
	  			  part = strtok (NULL, " ");
	  			  if(isCS) promptlambda_thetaphi[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) promptlambda_thetaphi[1][yBin][ptBin][generation] = atof(part);
				  part = strtok (NULL, " ");
				  part = strtok (NULL, " ");
				  if(isCS) errpromptlambda_thetaphi[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) errpromptlambda_thetaphi[1][yBin][ptBin][generation] = atof(part);
	  			  }

//////////////////// GER MINOS ERRORS //////////////////////////////////////////////

	  if(minos){
		  if(result31!=0) MinosresIndex=0;
		  if(MinosresIndex==4) {
			      char* part = strtok (zeile," ");
			      part = strtok (NULL, " ");
			      part = strtok (NULL, " ");
			      if(isCS) nPrompt[0][yBin][ptBin][generation] = atof(part);
			      if(isHX) nPrompt[1][yBin][ptBin][generation] = atof(part);
			      part = strtok (NULL, " ");
			      part = strtok (NULL, " ");
			      double negErr = atof(part);
			      part = strtok (NULL, " ");
			      double posErr = atof(part);
			      double MinosError = posErr;
			      if (fabs(negErr)>posErr) MinosError = negErr;
			      if(isCS) errnPrompt[0][yBin][ptBin][generation] = MinosError;
			      if(isHX) errnPrompt[1][yBin][ptBin][generation] = MinosError;
		  }
		  if(MinosresIndex==5) {
	  		      char* part = strtok (zeile," ");
	  			  part = strtok (NULL, " ");
	  			  part = strtok (NULL, " ");
	  			  if(isCS) promptlambda_phi[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) promptlambda_phi[1][yBin][ptBin][generation] = atof(part);
				  part = strtok (NULL, " ");
				  part = strtok (NULL, " ");
				  double negErr = atof(part);
				  part = strtok (NULL, " ");
				  double posErr = atof(part);
				  double MinosError = posErr;
				  if (fabs(negErr)>posErr) MinosError = negErr;
				  if(isCS) errpromptlambda_phi[0][yBin][ptBin][generation] = MinosError;
	  			  if(isHX) errpromptlambda_phi[1][yBin][ptBin][generation] = MinosError;
	  			  }
		  if(MinosresIndex==6) {
	  			  char* part = strtok (zeile," ");
	  			  part = strtok (NULL, " ");
	  			  part = strtok (NULL, " ");
	  			  if(isCS) promptlambda_theta[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) promptlambda_theta[1][yBin][ptBin][generation] = atof(part);
				  part = strtok (NULL, " ");
				  part = strtok (NULL, " ");
				  double negErr = atof(part);
				  part = strtok (NULL, " ");
				  double posErr = atof(part);
				  double MinosError = posErr;
				  if (fabs(negErr)>posErr) MinosError = negErr;
				  if(isCS) errpromptlambda_theta[0][yBin][ptBin][generation] = MinosError;
	  			  if(isHX) errpromptlambda_theta[1][yBin][ptBin][generation] = MinosError;
	  			  }
		  if(MinosresIndex==7) {
	  			  char* part = strtok (zeile," ");
	  			  part = strtok (NULL, " ");
	  			  part = strtok (NULL, " ");
	  			  if(isCS) promptlambda_thetaphi[0][yBin][ptBin][generation] = atof(part);
	  			  if(isHX) promptlambda_thetaphi[1][yBin][ptBin][generation] = atof(part);
				  part = strtok (NULL, " ");
				  part = strtok (NULL, " ");
				  double negErr = atof(part);
				  part = strtok (NULL, " ");
				  double posErr = atof(part);
				  double MinosError = posErr;
				  if (fabs(negErr)>posErr) MinosError = negErr;
				  if(isCS) errpromptlambda_thetaphi[0][yBin][ptBin][generation] = MinosError;
	  			  if(isHX) errpromptlambda_thetaphi[1][yBin][ptBin][generation] = MinosError;
	  			  }

	  }

	  fitresIndex++;
	  MinosresIndex++;
	  i++;
	} while(i!=3000000);
	  fclose(inputFile);


///////////////////////// END OF FILE READING ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// FIT THE RESULTS WITH GAUSSIAN ////////////////////////////////////////////////

	  	double cut=10;

	    RooRealVar *promptlambda_phi_ = new RooRealVar("promptlambda_phi_","promptlambda_phi_",-cut,cut);
	    RooRealVar *promptlambda_theta_ = new RooRealVar("promptlambda_theta_","promptlambda_theta_",-cut,cut);
	    RooRealVar *promptlambda_thetaphi_ = new RooRealVar("promptlambda_thetaphi_","promptlambda_thetaphi_",-cut,cut);

	    RooRealVar mean("mean","mean",0,-2,2);
	    RooRealVar sigma("sigma","sigma",1,0.01,10);
	    RooGaussian gauss_phi("gauss_phi","gauss_phi",*promptlambda_phi_,mean,sigma);
	    RooGaussian gauss_theta("gauss_theta","gauss_theta",*promptlambda_theta_,mean,sigma);
	    RooGaussian gauss_thetaphi("gauss_thetaphi","gauss_thetaphi",*promptlambda_thetaphi_,mean,sigma);


	    TCanvas* promptlambdaCanvas;
	    promptlambdaCanvas = new TCanvas("promptlambdaCanvas","promptlambdaCanvas",3800,3200);
	  	promptlambdaCanvas->Divide(2,3);  promptlambdaCanvas->SetFillColor(kWhite);


	    int PlotBins=25;
	    double borders=cut;

	    char Filename[200];
		sprintf(Filename,"Plots/ToyMC/ToyMCPlot_rapidity%d_pt%d.png",yBin,ptBin);

//	    	for(int yBin = 1; yBin < 6; yBin++) {
//	    		for(int ptBin = 1; ptBin < 13; ptBin++) {

	    yBin=yBinstart;
	  	ptBin=ptBinstart;


	    		    for(int frameDex = 0; frameDex<2; ++frameDex){





	    		    	RooArgSet varlist(*promptlambda_phi_,*promptlambda_theta_,*promptlambda_thetaphi_);
	    		    	RooDataSet* data = new RooDataSet("data","A sample",varlist);
	    		    	RooDataSet* dataRef = new RooDataSet("dataRef","A Ref sample",varlist);

	    if(CONV[frameDex][yBin][ptBin][generation]) fprintf(outputFile, "\n");
	    for(int generation = 1; generation < 100; generation++) {

/*	    promptlambda_phi_->setVal(promptlambda_phi[frameDex][yBin][ptBin][generation]);
	    promptlambda_phi_->setError(errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    promptlambda_theta_->setVal(promptlambda_theta[frameDex][yBin][ptBin][generation]);
	    promptlambda_theta_->setError(errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    promptlambda_thetaphi_->setVal(promptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
	    promptlambda_thetaphi_->setError(errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);
*/
	    promptlambda_phi_->setVal(promptlambda_phi[frameDex][yBin][ptBin][generation]/errpromptlambda_phi[frameDex][yBin][ptBin][generation]);
	    promptlambda_theta_->setVal(promptlambda_theta[frameDex][yBin][ptBin][generation]/errpromptlambda_theta[frameDex][yBin][ptBin][generation]);
	    promptlambda_thetaphi_->setVal(promptlambda_thetaphi[frameDex][yBin][ptBin][generation]/errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation]);

	    double th_valerr=promptlambda_theta[frameDex][yBin][ptBin][generation]/errpromptlambda_theta[frameDex][yBin][ptBin][generation];
	    double ph_valerr=promptlambda_phi[frameDex][yBin][ptBin][generation]/errpromptlambda_phi[frameDex][yBin][ptBin][generation];
	    double thph_valerr=promptlambda_thetaphi[frameDex][yBin][ptBin][generation]/errpromptlambda_thetaphi[frameDex][yBin][ptBin][generation];

	    if(CONV[frameDex][yBin][ptBin][generation]){
		dataRef->add(varlist);
		if(fabs(th_valerr)>cut)cout<<th_valerr<<" (generation) "<<generation<<promptlambda_theta[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<endl;
	    }
	    if(CONV[frameDex][yBin][ptBin][generation] && fabs(th_valerr) < cut && fabs(ph_valerr) < cut && fabs(thph_valerr) < cut  ){
//	    cout<<"lambda_theta "<<generation<<endl;//promptlambda_theta[frameDex][yBin][ptBin][generation]<<" +- "<<errpromptlambda_theta[frameDex][yBin][ptBin][generation]<<endl;
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
	    		    	if(doFit){

	    		    		for(int i=0;i<5;i++){
	    		    RooFitResult* phi_result = gauss_phi.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Minos(1),RooFit::Strategy(2));
	    		    phi_result->Print(); cout<<"Fit number "<<i+1<<endl;
	    		    if(sigma.getVal()<1.5)continue;
	    		    		}

//	    		    phi_result = gauss_phi.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Strategy(2));
//	    		    phi_result->Print(); cout<<"Fit number "<<i+1<<endl;

	    		    mean_phi[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    		    errmean_phi[frameDex][yBinstart][ptBinstart]=mean.getError();
	    		    sigma_phi[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    		    errsigma_phi[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    		    RooPlot* promptlambda_phi_frame = new RooPlot;
	    		    promptlambda_phi_frame = promptlambda_phi_->frame(-borders,borders,PlotBins) ;
	    		    promptlambda_phi_frame->SetTitle(0);
	    		    if(frameDex==0) {promptlambda_phi_frame->SetTitleSize(0.06,"X"); promptlambda_phi_frame->SetXTitle("#lambda_{#phi CS}/err(#lambda_{#phi CS})");}
	    		    if(frameDex==1) {promptlambda_phi_frame->SetTitleSize(0.06,"X"); promptlambda_phi_frame->SetXTitle("#lambda_{#phi HX}/err(#lambda_{#phi HX})");}
	    		    data->plotOn(promptlambda_phi_frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    		    gauss_phi.plotOn(promptlambda_phi_frame,LineWidth(2),Normalization(1.0));
	    			cout<<"heere"<<endl;
	    		    gauss_phi.paramOn(promptlambda_phi_frame,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));

	    			for(int i=0;i<5;i++){
	    			RooFitResult* theta_result = gauss_theta.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Strategy(2));
	    			theta_result->Print(); cout<<"Fit number "<<i+1<<endl;
	    		    if(sigma.getVal()<1.5)continue;
	    		    		}

	    			mean_theta[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    			errmean_theta[frameDex][yBinstart][ptBinstart]=mean.getError();
	    			sigma_theta[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    			errsigma_theta[frameDex][yBinstart][ptBinstart]=sigma.getError();

	    			RooPlot* promptlambda_theta_frame = new RooPlot;
	    			promptlambda_theta_frame = promptlambda_theta_->frame(-borders,borders,PlotBins) ;
	    			promptlambda_theta_frame->SetTitle(0);
	    			if(frameDex==0) {promptlambda_theta_frame->SetTitleSize(0.06,"X"); promptlambda_theta_frame->SetXTitle("#lambda_{#theta CS}/err(#lambda_{#theta CS})");}
	    			if(frameDex==1) {promptlambda_theta_frame->SetTitleSize(0.06,"X"); promptlambda_theta_frame->SetXTitle("#lambda_{#theta HX}/err(#lambda_{#theta HX})");}
	    			data->plotOn(promptlambda_theta_frame,DataError(RooAbsData::SumW2),MarkerSize(0.4));
	    			gauss_theta.plotOn(promptlambda_theta_frame,LineWidth(2),Normalization(1.0));
	    			gauss_theta.paramOn(promptlambda_theta_frame,Format("NEU",AutoPrecision(2)),Parameters(RooArgList(mean,sigma)),Layout(0.6,0.9,0.9));

	    			for(int i=0;i<5;i++){
	    			RooFitResult* thetaphi_result = gauss_thetaphi.fitTo(*data,Save(true),RooFit::Minimizer("Minuit","migrad"),RooFit::Strategy(2));
	    			thetaphi_result->Print(); cout<<"Fit number "<<i+1<<endl;
	    		    if(sigma.getVal()<1.5)continue;
	    		    		}

	    			mean_thetaphi[frameDex][yBinstart][ptBinstart]=mean.getVal();
	    			errmean_thetaphi[frameDex][yBinstart][ptBinstart]=mean.getError();
	    			sigma_thetaphi[frameDex][yBinstart][ptBinstart]=sigma.getVal();
	    			errsigma_thetaphi[frameDex][yBinstart][ptBinstart]=sigma.getError();


	    			RooPlot* promptlambda_thetaphi_frame = new RooPlot;
	    			promptlambda_thetaphi_frame = promptlambda_thetaphi_->frame(-borders,borders,PlotBins) ;
	    			promptlambda_thetaphi_frame->SetTitle(0); //promptlambda_thetaphi_frame->SetTitleOffset(0.005);
	    			if(frameDex==0) {promptlambda_thetaphi_frame->SetTitleSize(0.06,"X"); promptlambda_thetaphi_frame->SetXTitle("#lambda_{#theta#phi CS}/err(#lambda_{#theta#phi CS})");}
	    			if(frameDex==1) {promptlambda_thetaphi_frame->SetTitleSize(0.06,"X"); promptlambda_thetaphi_frame->SetXTitle("#lambda_{#theta#phi HX}/err(#lambda_{#theta#phi HX})");}
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


	    		    	}


	    		    	delete data;
	    		    }

//	    		}}


	    		    for(int frameDex = 0; frameDex<2; frameDex++){
	    		    for(int generation = 1; generation < 101; generation++) {
	    		    	if(CONV[frameDex][yBin][ptBin][generation]) convCountCheck[frameDex][yBin][ptBin]++;
	    		    	if(frameDex==0) cout<<"CS generation "<<generation<<": "<<CONV[frameDex][yBin][ptBin][generation]<<endl;
	    		    	if(frameDex==1) cout<<"HX generation "<<generation<<": "<<CONV[frameDex][yBin][ptBin][generation]<<endl;

	    		    }}
	    		    cout<<convCountCheck[0][yBinstart][ptBinstart]<<endl;
	    		    cout<<convCountCheck[1][yBinstart][ptBinstart]<<endl;

//	    		    cout<<convCount[0][1][5]<<endl;//<<convCount[0][1][6]<<convCount[0][2][5]<<convCount[0][2][6]<<endl;
//	    		    cout<<convCount[1][1][5]<<endl;//<<convCount[1][1][6]<<convCount[1][2][5]<<convCount[1][2][6]<<endl;

    		    	if(doFit){

	    		    promptlambdaCanvas->SaveAs(Filename);
	    		    promptlambdaCanvas->Close();
    		    	}

	    		    fprintf(outputFile, "\n");
	    			fprintf(outputFile, "EDM condition %1.5f\n",edmCondition);
	    			fprintf(outputFile, "CS convergence %d\n",convCountCheck[0][yBin][ptBin]);
	    			fprintf(outputFile, "HX convergence %d\n",convCountCheck[1][yBin][ptBin]);

	    			fprintf(outputFile2, "\n");
	    			fprintf(outputFile2, "Rapidity %d pT %d\n",yBin,ptBin);
	    			fprintf(outputFile2, "EDM condition %1.5f\n",edmCondition);
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

	    		}}



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
			fprintf(outputFile2, "{{%d,%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d,%d,%d}},\n",convCountCheck[0][1][2],convCountCheck[0][1][3],convCountCheck[0][1][4],convCountCheck[0][1][5],convCountCheck[0][1][6],convCountCheck[0][1][7],convCountCheck[0][2][2],convCountCheck[0][2][3],convCountCheck[0][2][4],convCountCheck[0][2][5],convCountCheck[0][2][6],convCountCheck[0][2][7],convCountCheck[0][2][8]);
			fprintf(outputFile2, "{{%d,%d,%d,%d,%d,%d},{%d,%d,%d,%d,%d,%d,%d}},\n",convCountCheck[1][1][2],convCountCheck[1][1][3],convCountCheck[1][1][4],convCountCheck[1][1][5],convCountCheck[1][1][6],convCountCheck[1][1][7],convCountCheck[1][2][2],convCountCheck[1][2][3],convCountCheck[1][2][4],convCountCheck[1][2][5],convCountCheck[1][2][6],convCountCheck[1][2][7],convCountCheck[1][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}},\n",overflow[0][1][2],overflow[0][1][3],overflow[0][1][4],overflow[0][1][5],overflow[0][1][6],overflow[0][1][7],overflow[0][2][2],overflow[0][2][3],overflow[0][2][4],overflow[0][2][5],overflow[0][2][6],overflow[0][2][7],overflow[0][2][8]);
			fprintf(outputFile2, "{{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f},{%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%1.5f}}\n",overflow[1][1][2],overflow[1][1][3],overflow[1][1][4],overflow[1][1][5],overflow[1][1][6],overflow[1][1][7],overflow[1][2][2],overflow[1][2][3],overflow[1][2][4],overflow[1][2][5],overflow[1][2][6],overflow[1][2][7],overflow[1][2][8]);


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


		    fclose(outputFile2);

  return 0;
}
