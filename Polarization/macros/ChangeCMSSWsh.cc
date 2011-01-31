#include <iostream>
#include <sstream>
#include <cstring>

//#include "commonVar.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

int main(int argc, char** argv) {
	using namespace std;
	using namespace RooFit;

	char* executable;

	  for( unsigned int arg = 1; arg < argc; ++arg ) {

	      executable = argv[arg];
	      cout<<"Changing executable cmsRun to "<<executable<<endl;

	  }

	  char line[200];

	  cout<<"Finding CrabJob:"<<endl;


	    char* part;

	    bool endisnear(false);
	    bool percentprob(false);

		char inputfilename[300];
		sprintf(inputfilename,"create.txt");
		FILE *inputFile = fopen(inputfilename,"rw");

		  if (inputFile == NULL)
		  {printf("\nQuelle \"%s\" falsch - kann Datei nicht šffnen\n", inputfilename);return -1;}

		  char zeile[500];

		  int i=0;
		  char* inputfilename2;
		  char* inputfilename3;
		  char* outputfilename;
		  char* outputfilename2;

	 do{
	  fgets(zeile, 200, inputFile);

	  char* str35 = "working directory";

	  char* result35 = strstr( zeile, str35 );

	  if(result35!=0) {
	  part = strtok (zeile," ");
	  for(int j=0;j<11;j++){
	  part = strtok (NULL, "////");
	  }

	  cout<<part<<endl;
	  inputfilename2=Form("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/CRAB/%s/job/CMSSW.sh",part);
//	  outputfilename=Form("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/CRAB/%s/job/CMSSW_2.sh",part);
	  outputfilename2=Form("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/CRAB/CMSSW_2.sh");
//	  inputfilename3=Form("/afs/hephy.at/scratch/k/knuenz/CMSSW_3_8_1/src/JPsiPolarizationSave3/CRAB/%s/job/CMSSW_3.sh",part);


	  }





	  i++;
	} while(i!=1000);
	  fclose(inputFile);

	  cout<<"Changing executable in "<<inputfilename2<<endl;

		inputFile = fopen(inputfilename2,"rw");


		  ofstream myfile;
		  myfile.open (outputfilename2);

		FILE *outputFile = fopen(outputfilename,"w");

		  if (inputFile == NULL)
		  {printf("\nQuelle \"%s\" falsch - kann Datei nicht šffnen\n", inputfilename2);return -1;}

		  i=0;

	  do{
	  	  fgets(zeile, 200, inputFile);

	  	  char* str32 = "executable=cmsRun";
	  	  char* str33 = "func_exit";
	  	  char* str34 = "# END OF PROCESS THE PRODUCED RESULTS";

	  	  char* result32 = strstr( zeile, str32 );
	  	  char* result33 = strstr( zeile, str33 );
	  	  char* result34 = strstr( zeile, str34 );


	  	  if(result32!=0) {cout<<"found line"<<endl; sprintf(zeile,"executable=$CMSSW_BASE/src/JPsiPolarizationSave3/data/%s\n",executable);}

		  myfile << zeile;

	  	 fprintf(outputFile, zeile);

	  	if(result34!=0) endisnear=true;

	  		  if(result33!=0 && endisnear) i=999;



	  	  i++;
	  	} while(i!=1000);
	  	  fclose(inputFile);
	  	  fclose(outputFile);
		  myfile.close();

			cout<<"Saving changed file"<<endl;

	  	endisnear=false;

	  	 ofstream outfile;
	     outfile.open (inputfilename2);
		 inputFile = fopen(outputfilename2,"rw");

			  if (inputFile == NULL)
			  {printf("\nQuelle \"%s\" falsch - kann Datei nicht šffnen\n", inputfilename2);return -1;}

			  i=0;

		  do{
		  	  fgets(zeile, 200, inputFile);
				  outfile << zeile;

		  		char* str33 = "func_exit";
		  	  	char* str34 = "# END OF PROCESS THE PRODUCED RESULTS";
		  	  	char* result33 = strstr( zeile, str33 );
		  	  	char* result34 = strstr( zeile, str34 );

		  	  	if(result34!=0) endisnear=true;

		  		  if(result33!=0 && endisnear) i=999;

		  	  i++;
		  	} while(i!=1000);
		  	  fclose(inputFile);
		      outfile.close();

	  return 0;
	}
