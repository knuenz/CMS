#include "TTree.h"

Char_t *fileNameData = "/scratch/knuenz/Polarization/RootInput/TTree_red_PR.root";

void ChangePolarization(){

	TFile* fInData = new TFile(fileNameData);

	TTree *treeData = (TTree*)fInData->Get("data");

	treeData->Print();

	double lambdaTh = 0.5;


	TTree *fChain;
	fChain=treeData;

	Long64_t nentries = fChain->GetEntries();
	cout<<nentries<<endl;
	Long64_t nb = 0;

	Double_t        JpsiMass;
	Double_t        JpsiPt;
	Double_t        JpsiRap;
	Double_t        Jpsict;
	Double_t        JpsictErr;
	Double_t        JpsiType_idx;
	Double_t        MCType_idx;
	Double_t        costh_CS;
	Double_t        phi_CS;
	Double_t        costh_HX;
	Double_t        phi_HX;
	Double_t        MCweight;

	TBranch        *b_JpsiMass;   //!
	TBranch        *b_JpsiPt;   //!
	TBranch        *b_JpsiRap;   //!
	TBranch        *b_Jpsict;   //!
	TBranch        *b_JpsictErr;   //!
	TBranch        *b_JpsiType_idx;   //!
	TBranch        *b_MCType_idx;   //!
	TBranch        *b_costh_CS;   //!
	TBranch        *b_phi_CS;   //!
	TBranch        *b_costh_HX;   //!
	TBranch        *b_phi_HX;   //!
	TBranch        *b_MCweight;   //!

	fChain->SetBranchAddress("JpsiMass", &JpsiMass, &b_JpsiMass);
	fChain->SetBranchAddress("JpsiPt", &JpsiPt, &b_JpsiPt);
	fChain->SetBranchAddress("JpsiRap", &JpsiRap, &b_JpsiRap);
	fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
	fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
	fChain->SetBranchAddress("JpsiType_idx", &JpsiType_idx, &b_JpsiType_idx);
	fChain->SetBranchAddress("MCType_idx", &MCType_idx, &b_MCType_idx);
	fChain->SetBranchAddress("costh_CS", &costh_CS, &b_costh_CS);
	fChain->SetBranchAddress("phi_CS", &phi_CS, &b_phi_CS);
	fChain->SetBranchAddress("costh_HX", &costh_HX, &b_costh_HX);
	fChain->SetBranchAddress("phi_HX", &phi_HX, &b_phi_HX);
	fChain->SetBranchAddress("MCweight", &MCweight, &b_MCweight);

	double newMCweight;

	double preweight=0.851;
	double pseudopreweight=7.414;

    TNtuple* data = new TNtuple("data","data","JpsiMass:JpsiPt:JpsiRap:Jpsict:JpsictErr:costh_CS:phi_CS:costh_HX:phi_HX:JpsiType_idx:MCType_idx:MCweight");
//    TNtuple* data = new TNtuple("data","data","JpsiMass:JpsiPt:JpsiRap:Jpsict:JpsictErr:costh_CS_prime:phi_CS_prime:costh_HX_prime:phi_HX_prime:JpsiType_idx:MCType_idx:MCweight");

    Double_t phiFolded[2];
    Double_t thetaAdjusted[2];
    Double_t thisPhi[2];
    Double_t thisCosTh[2];

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
   		    	 if(jentry % 10000 == 0) printf("event %d\n", (Int_t) jentry);

   	   nb = fChain->GetEntry(jentry);

   	thisPhi[0]=phi_CS-180;
   	thisPhi[1]=phi_HX-180;
   	thisCosTh[0]=costh_CS;
   	thisCosTh[1]=costh_HX;

	    for (int iFrame=0; iFrame<2;iFrame++) {

    	phiFolded [iFrame] = thisPhi[iFrame];
     	thetaAdjusted [iFrame] = thisCosTh[iFrame];
     	  if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.)
     	    phiFolded[iFrame]= phiFolded[iFrame]*(-1);
     	  else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
     	    phiFolded[iFrame] = 180. - thisPhi[iFrame];
     	    thetaAdjusted [iFrame]= thetaAdjusted [iFrame]*(-1);
     	  }
     	  else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
     	    phiFolded [iFrame]= 180. + thisPhi[iFrame];
     	    thetaAdjusted [iFrame]= thetaAdjusted [iFrame]*(-1);
     	  }

	    }


   	newMCweight=preweight;//*(1+lambdaTh*TMath::Power(costh_CS,2));

 //  	data->Fill(JpsiMass,JpsiPt,JpsiRap,Jpsict,JpsictErr,costh_CS,phi_CS,costh_HX,phi_HX,JpsiType_idx,MCType_idx,newMCweight);
   	data->Fill(JpsiMass,JpsiPt,JpsiRap,Jpsict,JpsictErr,thetaAdjusted[0],phiFolded[0],thetaAdjusted[1],phiFolded[1],JpsiType_idx,MCType_idx,newMCweight);


  	 if(jentry % 10000 == 0) printf("New weight %1.3f\n", newMCweight);

    }


    data->SaveAs("/scratch/knuenz/Polarization/RootInput/TTree_red_PR_folded.root");

	data->Print();


	return;
}
