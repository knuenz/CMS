
//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"

int main(int argc, char** argv) {

	bool CS=false;
	bool HX=false;

	int rap;
	int pt;
	int iter;

	for( int i=0;i < argc; ++i ) {
    if(std::string(argv[i]).find("--rap1") != std::string::npos) rap=1;
    if(std::string(argv[i]).find("--rap2") != std::string::npos) rap=2;
    if(std::string(argv[i]).find("--rap3") != std::string::npos) rap=3;
    if(std::string(argv[i]).find("--rap4") != std::string::npos) rap=4;
    if(std::string(argv[i]).find("--rap5") != std::string::npos) rap=5;
    if(std::string(argv[i]).find("--pt1") != std::string::npos) pt=1;
    if(std::string(argv[i]).find("--pt2") != std::string::npos) pt=2;
    if(std::string(argv[i]).find("--pt3") != std::string::npos) pt=3;
    if(std::string(argv[i]).find("--pt4") != std::string::npos) pt=4;
    if(std::string(argv[i]).find("--pt5") != std::string::npos) pt=5;
    if(std::string(argv[i]).find("--pt6") != std::string::npos) pt=6;
    if(std::string(argv[i]).find("--pt7") != std::string::npos) pt=7;
    if(std::string(argv[i]).find("--pt8") != std::string::npos) pt=8;
    if(std::string(argv[i]).find("--pt9") != std::string::npos) pt=9;
    if(std::string(argv[i]).find("--pt10") != std::string::npos) pt=10;
    if(std::string(argv[i]).find("--pt11") != std::string::npos) pt=11;
    if(std::string(argv[i]).find("--pt12") != std::string::npos) pt=12;

    if(std::string(argv[i]).find("--iter16") != std::string::npos) iter=16;
    if(std::string(argv[i]).find("--iter22") != std::string::npos) iter=22;
    if(std::string(argv[i]).find("--iter23") != std::string::npos) iter=23;
    if(std::string(argv[i]).find("--iter24") != std::string::npos) iter=24;


    if(std::string(argv[i]).find("--CS") != std::string::npos) CS=true;
    if(std::string(argv[i]).find("--HX") != std::string::npos) HX=true;

  }

		TFile* fInData = new TFile("/scratch/knuenz/Polarization/RootInput/trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root");
//		TFile* fInData = new TFile("DataGen_eff_plus05.root");

		char histname[200];

//if(CS){
		 sprintf(histname,"hAcc2D_CS_pT%d_rap%d",pt,rap);
	//	sprintf(histname,"GenData_HX_rap1_pt6_generation6__costh_HX_phi_HX");


//}
//if(HX){
//		 sprintf(histname,"hGen_HX_pT%d_rap%d",pt,rap);
//}
		TObject* accmap= fInData->Get(histname);
		accmap->Print("all");
				return 0;
}
