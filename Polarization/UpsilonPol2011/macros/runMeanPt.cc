#include <iostream>
#include <string>
#include <sstream>
using namespace std;


#include "rootIncludes.inc"
#include "commonVar.h"
#include "calcMeanPt.C"

#include "TROOT.h"

//void calcMeanPt(Int_t iRapBin, Int_t iPTBin, Double_t nSigma, Int_t nUpsState, Int_t output);

//====================================
int main(int argc, char** argv){

	Double_t nSigma = 2.;
	Double_t fracL = 0.5;

	  for( int i=0;i < argc; ++i ) {
		    if(std::string(argv[i]).find("FracLSB") != std::string::npos) {char* fracLchar = argv[i]; char* fracLchar2 = strtok (fracLchar, "p"); fracL = atof(fracLchar2); fracL=fracL/100; cout<<"fracLSB = "<<fracL<<endl;}
		    if(std::string(argv[i]).find("nSigma") != std::string::npos) {char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "p"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl;}
	    }

	gROOT->ProcessLine(".L calcMeanPt.C+");

	char filename[200];
	sprintf(filename,"AllStates_%1.2fSigma_FracLSB%dPercent/meanPtlll.txt", nSigma, int(fracL*100));
	gSystem->Unlink(filename);


	  int numEvents=0;

		  for(int iState = 0; iState < 3; iState++){

			  cout<<"Upsilon("<<iState+1<<"):"<<endl;

  cout<<endl;
  cout<<endl;

  cout<<"double ptCentre[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,0, fracL);
    	  if(iPT<max-1) {cout<<", ";  }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{";}
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;

  cout<<"int numEvents[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,1, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};"; }
      }
    }
//  }
  cout<<endl;
  cout<<endl;

  cout<<"double fracBackground[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,2, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;
/*
  cout<<"double massMean[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,3, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;

  cout<<"double fracSignal[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,4, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;
*/

  cout<<"double meanRap[nRapBins][nPtBins]={{";

//  for(int iState = iState_-1; iState < iState_; iState++){
    for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
      Int_t max = onia::kNbPTBins[iRap]+1;
      for(int iPT = 1; iPT < max; iPT++){
    	  calcMeanPt(iRap, iPT, nSigma, iState,5, fracL);
    	  if(iPT<max-1) {cout<<", "; }
    	  else if(iPT==max-1 && iRap < onia::kNbRapForPTBins) {cout<<"},{"; }
    	  else {cout<<"}};";}
      }
    }
//  }
  cout<<endl;
  cout<<endl;

  int numEventsTotal=0;
  for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){
    Int_t max = onia::kNbPTBins[iRap]+1;
    for(int iPT = 1; iPT < max; iPT++){
                 numEventsTotal+=numEvents_[iState][iRap][iPT];
    }}
  cout<<"Total Number of Signal Events: "<<numEventsTotal<<endl<<endl;



		  }





			FILE *NumFile;

			bool SaveTables=true;
		    if(SaveTables){


		    char framerap[200];

			char NumFileName[200];
			sprintf(NumFileName,"AllStates_%1.2fSigma_FracLSB%dPercent/meanPt.tex", nSigma, int(fracL*100));
			NumFile = fopen(NumFileName,"w");

			fprintf(NumFile, "\n");
			fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

			fprintf(NumFile, "\n\n\n\n");

		    	int nTables=1;
			    for(int iTab=1; iTab<nTables+1;iTab++){

					fprintf(NumFile, "\n\n\n\n");

					if(iTab==1){
						fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Estimated mean $p_{T}$ and $|y|$, for each kinematic bin of the $\\Upsilon(nS)$ analyses}\n \\begin{tabular}{|c|cccccc|}\n\\hline\n");
						fprintf(NumFile, "$p_{T}$ [GeV] & $\\hat{p_{T}}^{\\Upsilon(1S)} [GeV] $ & $\\hat{|y|}^{\\Upsilon(1S)}$ & $\\hat{p_{T}}^{\\Upsilon(2S)} [GeV] $ & $\\hat{|y|}^{\\Upsilon(2S)}$ & $\\hat{p_{T}}^{\\Upsilon(3S)} [GeV] $ & $\\hat{|y|}^{\\Upsilon(3S)}$\\\\\n");
					}
					if(iTab==2){
						fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Estimated number of signal events and background fraction, as a function of $p_{T}$}\n \\begin{tabular}{|c|cccccc|}\n\\hline\n");
						fprintf(NumFile, "$p_{T}$ [GeV] & $N_{\\Upsilon(1S)}$ & $f_{BG}^{\\Upsilon(1S)} [\\%] $ & $N_{\\Upsilon(2S)}$ & $f_{BG}^{\\Upsilon(2S)} [\\%]$& $N_{\\Upsilon(3S)}$ & $f_{BG}^{\\Upsilon(3S)} [\\%]$ \\\\\n");
					}



				int rap=0;
			    for(int rapBin = 1; rapBin < 3; rapBin++) {

			    	sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{$%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);

						int pt=0;
							for(int ptBin = 6; ptBin < 11; ptBin++) {

								if(iTab==1){
									fprintf(NumFile, "%1.0f--%1.0f   &  $%1.3f  $  & $%1.3f $  &  $%1.3f $ &  $%1.3f $ &  $%1.3f $ &  $%1.3f $\\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],meanPt[0][rapBin][ptBin] ,meanRap[0][rapBin][ptBin],meanPt[1][rapBin][ptBin] ,meanRap[1][rapBin][ptBin],meanPt[2][rapBin][ptBin] ,meanRap[2][rapBin][ptBin]);
								}
								if(iTab==2){
									fprintf(NumFile, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%d $ &  $%1.1f $ &  $%d$ &  $%1.1f $\\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],numEvents_[0][rapBin][ptBin] ,BackgroundFrac[0][rapBin][ptBin],numEvents_[1][rapBin][ptBin] ,BackgroundFrac[1][rapBin][ptBin],numEvents_[2][rapBin][ptBin] ,BackgroundFrac[2][rapBin][ptBin]);
								}


								pt++;
							}



							rap++;

			    }//end rapBin




				fprintf(NumFile, "\\hline\n");
				fprintf(NumFile, "\\end{tabular}\n");
				if(iTab==1) fprintf(NumFile, "\\label{tab:MeanKinematics}\n");
				if(iTab==2) fprintf(NumFile, "\\label{tab:NumEventsFracBG}\n");
				fprintf(NumFile, "\\end{table}\n");
				fprintf(NumFile, "\n");

			    }//end iTab



		    }

				fprintf(NumFile, "\\end{document}");

				fclose(NumFile);



  return 0;

}

