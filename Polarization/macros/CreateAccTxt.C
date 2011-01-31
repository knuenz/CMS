#include "TFile.h"
#include "TH2.h"

void CreateAccTxt(){



		TFile* fInData = new TFile("/scratch/knuenz/Polarization/RootInput/geomAcc_WithFSR_uniform_8Dec2010_merged.root");

		TH2F* accmap= fInData->Get("hGen_HX_pT10_rap5");
		accmap->SaveAs("accMap.xml");
				return;
}
