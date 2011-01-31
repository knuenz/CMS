//'.include /afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include';
#include "RooDataSet.h"
#include <map>
//#include <vector>

void MergeRooDataSet(){
using namespace std;
    Jpsi_Mass = new RooRealVar("JpsiMass","J/psi mass",0,50,"GeV/c^{2}");
    Jpsi_Pt = new RooRealVar("JpsiPt","J/psi pt",0,100,"GeV/c");
    Jpsi_Rap = new RooRealVar("JpsiRap","J/psi eta",-10,10);
    Jpsi_ct = new RooRealVar("Jpsict","J/psi ctau",-10,10,"mm");
    Jpsi_ctErr = new RooRealVar("JpsictErr","J/psi ctau error",-10,10,"mm");
    costh_CS = new RooRealVar("costh_CS","costh_CS",-1.,1.);
    phi_CS = new RooRealVar("phi_CS","phi_CS",-180.,180.,"deg");
    costh_HX = new RooRealVar("costh_HX","costh_HX",-1.,1.);
    phi_HX = new RooRealVar("phi_HX","phi_HX",-180.,180.,"deg");
    costh_PHX = new RooRealVar("costh_PHX","costh_PHX",-1.,1.);
    phi_PHX = new RooRealVar("phi_PHX","phi_PHX",-180.,180.,"deg");
    costh_sGJ = new RooRealVar("costh_sGJ","costh_sGJ",-1.,1.);
    phi_sGJ = new RooRealVar("phi_sGJ","phi_sGJ",-180.,180.,"deg");
    costh_GJ1 = new RooRealVar("costh_GJ1","costh_GJ1",-1.,1.);
    phi_GJ1 = new RooRealVar("phi_GJ1","phi_GJ1",-180.,180.,"deg");
    costh_GJ2 = new RooRealVar("costh_GJ2","costh_GJ2",-1.,1.);
    phi_GJ2 = new RooRealVar("phi_GJ2","phi_GJ2",-180.,180.,"deg");
    muPos_Pt = new RooRealVar("muPos_Pt","muPos_Pt",-1000.,1000.);
    muNeg_Pt = new RooRealVar("muNeg_Pt","muNeg_Pt",-1000.,1000.);
    muPos_Eta = new RooRealVar("muPos_Eta","muPos_Eta",-10.,10.);
    muNeg_Eta = new RooRealVar("muNeg_Eta","muNeg_Eta",-10.,10.);
    muPos_Phi = new RooRealVar("muPos_Phi","muPos_Phi",-4.,4.);
    muNeg_Phi = new RooRealVar("muNeg_Phi","muNeg_Phi",-4.,4.);

    RooArgList varlist(*Jpsi_Mass,*Jpsi_ct,*Jpsi_Pt,*Jpsi_Rap);
    varlist.add(*Jpsi_ctErr);
    varlist.add(*costh_CS); varlist.add(*phi_CS); //varlist.add(*phi_CS_prime);//varlist.add(*costh_CS_prime);
    varlist.add(*costh_HX); varlist.add(*phi_HX);// varlist.add(*costh_HX_prime);  varlist.add(*phi_HX_prime);
    varlist.add(*costh_PHX); varlist.add(*phi_PHX);// varlist.add(*costh_PHX_prime); varlist.add(*phi_PHX_prime);
    varlist.add(*costh_sGJ); varlist.add(*phi_sGJ);// varlist.add(*costh_sGJ_prime); varlist.add(*phi_sGJ_prime);
    varlist.add(*costh_GJ1); varlist.add(*phi_GJ1);// varlist.add(*costh_GJ1_prime); varlist.add(*phi_GJ1_prime);
    varlist.add(*costh_GJ2); varlist.add(*phi_GJ2);// varlist.add(*costh_GJ2_prime); varlist.add(*phi_GJ2_prime);
    varlist.add(*muPos_Pt); varlist.add(*muPos_Eta); varlist.add(*muPos_Phi); varlist.add(*muNeg_Pt); varlist.add(*muNeg_Eta); varlist.add(*muNeg_Phi);

	TFile* fIn1 = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/RooDataSet_Nov04_RunB_folded_1.root");
	TFile* fIn2 = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/RooDataSet_Nov04_RunB_folded_2.root");
	TFile* fIn3 = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/RooDataSet_Nov04_RunB_1.root");
	TFile* fIn4 = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/RooDataSet_Nov04_RunB_2.root");

//	map<string,RooDataSet*> datasets1; datasets2["ds4"]=
//	map<string,RooDataSet*> datasets2;

	RooDataSet* data1=(RooDataSet*)fIn1->Get("data");
	RooDataSet* data2=(RooDataSet*)fIn2->Get("data");
	RooDataSet* data3=(RooDataSet*)fIn3->Get("data");
	RooDataSet* data4=(RooDataSet*)fIn4->Get("data");

	cout<<"data1 entries "<<data1->numEntries()<<endl;
	cout<<"data2 entries "<<data2->numEntries()<<endl;
	cout<<"data3 entries "<<data3->numEntries()<<endl;
	cout<<"data4 entries "<<data4->numEntries()<<endl;

	cout<<"data1 +data2 entries "<<data1->numEntries()+data2->numEntries()<<endl;
	cout<<"data3 +data4 entries "<<data3->numEntries()+data4->numEntries()<<endl;


	data1->append(*data2);
	data3->append(*data4);

//	RooDataSet* merge1 = new RooDataSet("merge1","merge1",varlist,Import(*data1,*data2));


//	data1.append(data2);
//	data3.append(data4);

	cout<<"data12 merged entries "<<data1->numEntries()<<endl;

	cout<<"data34 merged entries "<<data3->numEntries()<<endl;


return;
}
