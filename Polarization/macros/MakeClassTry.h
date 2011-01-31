//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 13 11:41:58 2010 by ROOT version 5.22/00d
// from TTree data/CMSSW Quarkonia J/psi Polarization+Trigger Tree
// found on file: /tmp_mnt/scratch/knuenz/Polarization/RootInput/TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root
//////////////////////////////////////////////////////////

#ifndef MakeClassTry_h
#define MakeClassTry_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MakeClassTry {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           eventNb;
   Int_t           runNb;
   Int_t           lumiBlock;
   Int_t           nPriVtx;
   Int_t           JpsiType;
 //TLorentzVector  *JpsiP;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          fP_fUniqueID;
   UInt_t          fP_fBits;
   Double_t        fP_fX;
   Double_t        fP_fY;
   Double_t        fP_fZ;
   Double_t        fE;
   Int_t           JpsiCharge;
   Double_t        Jpsict;
   Double_t        JpsictErr;
   Double_t        JpsiVprob;
 //TLorentzVector  *muPosP;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          fP_fUniqueID;
   UInt_t          fP_fBits;
   Double_t        fP_fX;
   Double_t        fP_fY;
   Double_t        fP_fZ;
   Double_t        fE;
 //TLorentzVector  *muNegP;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          fP_fUniqueID;
   UInt_t          fP_fBits;
   Double_t        fP_fX;
   Double_t        fP_fY;
   Double_t        fP_fZ;
   Double_t        fE;
   Int_t           HLT_Mu0_Track0_Jpsi;
   Int_t           HLT_Mu3_Track0_Jpsi;
   Int_t           HLT_Mu5_Track0_Jpsi;
   Int_t           HLT_Mu0_TkMu0_Jpsi;
   Int_t           HLT_Mu3_TkMu0_Jpsi;
   Int_t           HLT_Mu5_TkMu0_Jpsi;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;
   Int_t           HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;
   Int_t           HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;
   Int_t           HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;
   Int_t           HLT_Mu3_Track3_Jpsi;
   Int_t           HLT_Mu3_Track3_Jpsi_v2;
   Int_t           HLT_Mu3_Track5_Jpsi_v1;
   Int_t           HLT_Mu3_Track5_Jpsi_v2;
   Int_t           HLT_DoubleMu0;
   Int_t           HLT_DoubleMu0_Quarkonium_v1;
   Int_t           HLT_DoubleMu0_Quarkonium_LS_v1;
   Int_t           HLT_L1DoubleMuOpen;
   Int_t           HLT_L1DoubleMuOpen_Tight;
   Int_t           HLT_DoubleMu3;
   Int_t           HLT_Mu3;
   Int_t           HLT_Mu5;
   Int_t           HLT_Mu7;
   Int_t           HLT_Mu9;
   Int_t           HLT_Mu11;

   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_nPriVtx;   //!
   TBranch        *b_JpsiType;   //!
   TBranch        *b_JpsiP_fUniqueID;   //!
   TBranch        *b_JpsiP_fBits;   //!
   TBranch        *b_JpsiP_fP_fUniqueID;   //!
   TBranch        *b_JpsiP_fP_fBits;   //!
   TBranch        *b_JpsiP_fP_fX;   //!
   TBranch        *b_JpsiP_fP_fY;   //!
   TBranch        *b_JpsiP_fP_fZ;   //!
   TBranch        *b_JpsiP_fE;   //!
   TBranch        *b_JpsiCharge;   //!
   TBranch        *b_Jpsict;   //!
   TBranch        *b_JpsictErr;   //!
   TBranch        *b_JpsiVprob;   //!
   TBranch        *b_muPosP_fUniqueID;   //!
   TBranch        *b_muPosP_fBits;   //!
   TBranch        *b_muPosP_fP_fUniqueID;   //!
   TBranch        *b_muPosP_fP_fBits;   //!
   TBranch        *b_muPosP_fP_fX;   //!
   TBranch        *b_muPosP_fP_fY;   //!
   TBranch        *b_muPosP_fP_fZ;   //!
   TBranch        *b_muPosP_fE;   //!
   TBranch        *b_muNegP_fUniqueID;   //!
   TBranch        *b_muNegP_fBits;   //!
   TBranch        *b_muNegP_fP_fUniqueID;   //!
   TBranch        *b_muNegP_fP_fBits;   //!
   TBranch        *b_muNegP_fP_fX;   //!
   TBranch        *b_muNegP_fP_fY;   //!
   TBranch        *b_muNegP_fP_fZ;   //!
   TBranch        *b_muNegP_fE;   //!
   TBranch        *b_HLT_Mu0_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu3_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu5_Track0_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu3_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu5_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;   //!
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;   //!
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu3_Track5_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu3_Track5_Jpsi_v2;   //!
   TBranch        *b_HLT_DoubleMu0;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_v1;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_LS_v1;   //!
   TBranch        *b_HLT_L1DoubleMuOpen;   //!
   TBranch        *b_HLT_L1DoubleMuOpen_Tight;   //!
   TBranch        *b_HLT_DoubleMu3;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu7;   //!
   TBranch        *b_HLT_Mu9;   //!
   TBranch        *b_HLT_Mu11;   //!

   MakeClassTry(TTree *tree=0);
   virtual ~MakeClassTry();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MakeClassTry_cxx
MakeClassTry::MakeClassTry(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp_mnt/scratch/knuenz/Polarization/RootInput/TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root");
      if (!f) {
         f = new TFile("/tmp_mnt/scratch/knuenz/Polarization/RootInput/TTree_pol_noTriggerFilter_Run2010A-Nov4ReReco_v1_06Dec2010.root");
      }
      tree = (TTree*)gDirectory->Get("data");

   }
   Init(tree);
}

MakeClassTry::~MakeClassTry()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MakeClassTry::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MakeClassTry::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MakeClassTry::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("nPriVtx", &nPriVtx, &b_nPriVtx);
   fChain->SetBranchAddress("JpsiType", &JpsiType, &b_JpsiType);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_JpsiP_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_JpsiP_fBits);
   fChain->SetBranchAddress("fP.fUniqueID", &fP_fUniqueID, &b_JpsiP_fP_fUniqueID);
   fChain->SetBranchAddress("fP.fBits", &fP_fBits, &b_JpsiP_fP_fBits);
   fChain->SetBranchAddress("fP.fX", &fP_fX, &b_JpsiP_fP_fX);
   fChain->SetBranchAddress("fP.fY", &fP_fY, &b_JpsiP_fP_fY);
   fChain->SetBranchAddress("fP.fZ", &fP_fZ, &b_JpsiP_fP_fZ);
   fChain->SetBranchAddress("fE", &fE, &b_JpsiP_fE);
   fChain->SetBranchAddress("JpsiCharge", &JpsiCharge, &b_JpsiCharge);
   fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
   fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
   fChain->SetBranchAddress("JpsiVprob", &JpsiVprob, &b_JpsiVprob);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_muPosP_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_muPosP_fBits);
   fChain->SetBranchAddress("fP.fUniqueID", &fP_fUniqueID, &b_muPosP_fP_fUniqueID);
   fChain->SetBranchAddress("fP.fBits", &fP_fBits, &b_muPosP_fP_fBits);
   fChain->SetBranchAddress("fP.fX", &fP_fX, &b_muPosP_fP_fX);
   fChain->SetBranchAddress("fP.fY", &fP_fY, &b_muPosP_fP_fY);
   fChain->SetBranchAddress("fP.fZ", &fP_fZ, &b_muPosP_fP_fZ);
   fChain->SetBranchAddress("fE", &fE, &b_muPosP_fE);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_muNegP_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_muNegP_fBits);
   fChain->SetBranchAddress("fP.fUniqueID", &fP_fUniqueID, &b_muNegP_fP_fUniqueID);
   fChain->SetBranchAddress("fP.fBits", &fP_fBits, &b_muNegP_fP_fBits);
   fChain->SetBranchAddress("fP.fX", &fP_fX, &b_muNegP_fP_fX);
   fChain->SetBranchAddress("fP.fY", &fP_fY, &b_muNegP_fP_fY);
   fChain->SetBranchAddress("fP.fZ", &fP_fZ, &b_muNegP_fP_fZ);
   fChain->SetBranchAddress("fE", &fE, &b_muNegP_fE);
   fChain->SetBranchAddress("HLT_Mu0_Track0_Jpsi", &HLT_Mu0_Track0_Jpsi, &b_HLT_Mu0_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_Track0_Jpsi", &HLT_Mu3_Track0_Jpsi, &b_HLT_Mu3_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_Track0_Jpsi", &HLT_Mu5_Track0_Jpsi, &b_HLT_Mu5_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_Jpsi", &HLT_Mu0_TkMu0_Jpsi, &b_HLT_Mu0_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_Jpsi", &HLT_Mu3_TkMu0_Jpsi, &b_HLT_Mu3_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_Jpsi", &HLT_Mu5_TkMu0_Jpsi, &b_HLT_Mu5_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi", &HLT_Mu0_TkMu0_OST_Jpsi, &b_HLT_Mu0_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi", &HLT_Mu3_TkMu0_OST_Jpsi, &b_HLT_Mu3_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi", &HLT_Mu5_TkMu0_OST_Jpsi, &b_HLT_Mu5_TkMu0_OST_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3);
   fChain->SetBranchAddress("HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3", &HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3, &b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi", &HLT_Mu3_Track3_Jpsi, &b_HLT_Mu3_Track3_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v2", &HLT_Mu3_Track3_Jpsi_v2, &b_HLT_Mu3_Track3_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v1", &HLT_Mu3_Track5_Jpsi_v1, &b_HLT_Mu3_Track5_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu3_Track5_Jpsi_v2", &HLT_Mu3_Track5_Jpsi_v2, &b_HLT_Mu3_Track5_Jpsi_v2);
   fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &HLT_DoubleMu0_Quarkonium_v1, &b_HLT_DoubleMu0_Quarkonium_v1);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_LS_v1", &HLT_DoubleMu0_Quarkonium_LS_v1, &b_HLT_DoubleMu0_Quarkonium_LS_v1);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, &b_HLT_L1DoubleMuOpen);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen_Tight", &HLT_L1DoubleMuOpen_Tight, &b_HLT_L1DoubleMuOpen_Tight);
   fChain->SetBranchAddress("HLT_DoubleMu3", &HLT_DoubleMu3, &b_HLT_DoubleMu3);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   Notify();
}

Bool_t MakeClassTry::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MakeClassTry::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MakeClassTry::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MakeClassTry_cxx
