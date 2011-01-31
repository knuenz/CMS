//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 28 14:54:21 2010 by ROOT version 5.26/00
// from TTree data/new data
// found on file: tree_promptJPsi_v4.root
//////////////////////////////////////////////////////////

#ifndef PolData_h
#define PolData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class PolData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        JpsiMass;
   Double_t        JpsiPt;
   Double_t        JpsiRap;
   Double_t        JpsiPx;
   Double_t        JpsiPy;
   Double_t        JpsiPz;
   Double_t        Jpsict;
   Double_t        JpsictErr;
   Int_t           JpsiType_idx;
   Char_t          JpsiType_lbl[3];
   Double_t        muPosPx;
   Double_t        muPosPy;
   Double_t        muPosPz;
   Double_t        muNegPx;
   Double_t        muNegPy;
   Double_t        muNegPz;
   /* Double_t        JpsiMass_Gen; */
   /* Double_t        JpsiPt_Gen; */
   /* Double_t        JpsiRap_Gen; */
   /* Double_t        JpsiPx_Gen; */
   /* Double_t        JpsiPy_Gen; */
   /* Double_t        JpsiPz_Gen; */
   /* Double_t        muPosPx_Gen; */
   /* Double_t        muPosPy_Gen; */
   /* Double_t        muPosPz_Gen; */
   /* Double_t        muNegPx_Gen; */
   /* Double_t        muNegPy_Gen; */
   /* Double_t        muNegPz_Gen; */
   Int_t        HLT_Mu0_Track0_Jpsi;
   Int_t        HLT_Mu3_Track0_Jpsi;
   Int_t        HLT_Mu5_Track0_Jpsi;
   Int_t        HLT_Mu0_TkMu0_Jpsi;
   Int_t        HLT_Mu3_TkMu0_Jpsi;
   Int_t        HLT_Mu5_TkMu0_Jpsi;
   Int_t        HLT_Mu0_TkMu0_OST_Jpsi;
   Int_t        HLT_Mu3_TkMu0_OST_Jpsi;
   Int_t        HLT_Mu5_TkMu0_OST_Jpsi;
   Int_t        HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;
   Int_t        HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;
   Int_t        HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;
   Int_t        HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
   Int_t        HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;
   Int_t        HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;
   Int_t        HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;
   Int_t        HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;
   Int_t        HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;
   Int_t        HLT_DoubleMu0;
   Int_t        HLT_DoubleMu0_Quarkonium_v1;
   Int_t        HLT_DoubleMu3;
   Int_t        HLT_L1DoubleMuOpen;
   Int_t        HLT_L1DoubleMuOpen_Tight;
   Int_t        HLT_Mu3;
   Int_t        HLT_Mu5;
   Int_t        HLT_Mu7;
   Int_t        HLT_Mu9;
   Int_t        HLT_Mu11;
   Double_t        eventNb;
   Double_t        runNb;
   Double_t        lumiBlock;
   /* Int_t           MCType_idx; */
   /* Char_t          MCType_lbl[3]; */
   Double_t        costh_HX;
   Double_t        costh_CS;
   Double_t        phi_HX;
   Double_t        phi_CS;

   // List of branches
   TBranch        *b_JpsiMass;   //!
   TBranch        *b_JpsiPt;   //!
   TBranch        *b_JpsiRap;   //!
   TBranch        *b_JpsiPx;   //!
   TBranch        *b_JpsiPy;   //!
   TBranch        *b_JpsiPz;   //!
   TBranch        *b_Jpsict;   //!
   TBranch        *b_JpsictErr;   //!
   TBranch        *b_JpsiType_idx;   //!
   TBranch        *b_JpsiType_lbl;   //!
   TBranch        *b_muPosPx;   //!
   TBranch        *b_muPosPy;   //!
   TBranch        *b_muPosPz;   //!
   TBranch        *b_muNegPx;   //!
   TBranch        *b_muNegPy;   //!
   TBranch        *b_muNegPz;   //!
   /* TBranch        *b_JpsiMass_Gen;   //! */
   /* TBranch        *b_JpsiPt_Gen;   //! */
   /* TBranch        *b_JpsiRap_Gen;   //! */
   /* TBranch        *b_JpsiPx_Gen;   //! */
   /* TBranch        *b_JpsiPy_Gen;   //! */
   /* TBranch        *b_JpsiPz_Gen;   //! */
   /* TBranch        *b_muPosPx_Gen;   //! */
   /* TBranch        *b_muPosPy_Gen;   //! */
   /* TBranch        *b_muPosPz_Gen;   //! */
   /* TBranch        *b_muNegPx_Gen;   //! */
   /* TBranch        *b_muNegPy_Gen;   //! */
   /* TBranch        *b_muNegPz_Gen;   //! */
   TBranch        *b_HLT_Mu0_Track0_Jpsi;
   TBranch        *b_HLT_Mu3_Track0_Jpsi;
   TBranch        *b_HLT_Mu5_Track0_Jpsi;
   TBranch        *b_HLT_Mu0_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu3_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu5_TkMu0_Jpsi;   //!
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi;
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi;
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi;
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1;
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v1;
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v1;
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2;
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v2;
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v2;
   TBranch        *b_HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3;
   TBranch        *b_HLT_Mu3_TkMu0_OST_Jpsi_Tight_v3;
   TBranch        *b_HLT_Mu5_TkMu0_OST_Jpsi_Tight_v3;
   TBranch        *b_HLT_DoubleMu0;   //!
   TBranch        *b_HLT_DoubleMu0_Quarkonium_v1;   //!
   TBranch        *b_HLT_DoubleMu3;   //!
   TBranch        *b_HLT_L1DoubleMuOpen;   //!
   TBranch        *b_HLT_L1DoubleMuOpen_Tight;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu7;   //!
   TBranch        *b_HLT_Mu9;   //!
   TBranch        *b_HLT_Mu11;   //!
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   /* TBranch        *b_MCType_idx;   //! */
   /* TBranch        *b_MCType_lbl;   //! */
   TBranch        *b_costh_HX;   //!
   TBranch        *b_costh_CS;   //!
   TBranch        *b_phi_HX;   //!
   TBranch        *b_phi_CS;   //!

   PolData(TTree *tree=0);
   virtual ~PolData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t selDimuType, Bool_t writeOutEvents);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PolData_cxx
PolData::PolData(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/hwoehri/TTree_pol_Mu0Track0Jpsi_MCprompt.root");
      if (!f) {
         f = new TFile("/tmp/hwoehri/TTree_pol_Mu0Track0Jpsi_MCprompt.root");
      }
      tree = (TTree*)gDirectory->Get("data");

   }
   Init(tree);
}

PolData::~PolData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PolData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PolData::LoadTree(Long64_t entry)
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

void PolData::Init(TTree *tree)
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

   fChain->SetBranchAddress("JpsiMass", &JpsiMass, &b_JpsiMass);
   fChain->SetBranchAddress("JpsiPt", &JpsiPt, &b_JpsiPt);
   fChain->SetBranchAddress("JpsiRap", &JpsiRap, &b_JpsiRap);
   fChain->SetBranchAddress("JpsiPx", &JpsiPx, &b_JpsiPx);
   fChain->SetBranchAddress("JpsiPy", &JpsiPy, &b_JpsiPy);
   fChain->SetBranchAddress("JpsiPz", &JpsiPz, &b_JpsiPz);
   fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
   fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
   fChain->SetBranchAddress("JpsiType_idx", &JpsiType_idx, &b_JpsiType_idx);
   fChain->SetBranchAddress("JpsiType_lbl", JpsiType_lbl, &b_JpsiType_lbl);
   fChain->SetBranchAddress("muPosPx", &muPosPx, &b_muPosPx);
   fChain->SetBranchAddress("muPosPy", &muPosPy, &b_muPosPy);
   fChain->SetBranchAddress("muPosPz", &muPosPz, &b_muPosPz);
   fChain->SetBranchAddress("muNegPx", &muNegPx, &b_muNegPx);
   fChain->SetBranchAddress("muNegPy", &muNegPy, &b_muNegPy);
   fChain->SetBranchAddress("muNegPz", &muNegPz, &b_muNegPz);
   /* fChain->SetBranchAddress("JpsiMass_Gen", &JpsiMass_Gen, &b_JpsiMass_Gen); */
   /* fChain->SetBranchAddress("JpsiPt_Gen", &JpsiPt_Gen, &b_JpsiPt_Gen); */
   /* fChain->SetBranchAddress("JpsiRap_Gen", &JpsiRap_Gen, &b_JpsiRap_Gen); */
   /* fChain->SetBranchAddress("JpsiPx_Gen", &JpsiPx_Gen, &b_JpsiPx_Gen); */
   /* fChain->SetBranchAddress("JpsiPy_Gen", &JpsiPy_Gen, &b_JpsiPy_Gen); */
   /* fChain->SetBranchAddress("JpsiPz_Gen", &JpsiPz_Gen, &b_JpsiPz_Gen); */
   /* fChain->SetBranchAddress("muPosPx_Gen", &muPosPx_Gen, &b_muPosPx_Gen); */
   /* fChain->SetBranchAddress("muPosPy_Gen", &muPosPy_Gen, &b_muPosPy_Gen); */
   /* fChain->SetBranchAddress("muPosPz_Gen", &muPosPz_Gen, &b_muPosPz_Gen); */
   /* fChain->SetBranchAddress("muNegPx_Gen", &muNegPx_Gen, &b_muNegPx_Gen); */
   /* fChain->SetBranchAddress("muNegPy_Gen", &muNegPy_Gen, &b_muNegPy_Gen); */
   /* fChain->SetBranchAddress("muNegPz_Gen", &muNegPz_Gen, &b_muNegPz_Gen); */
   fChain->SetBranchAddress("HLT_Mu0_Track0_Jpsi", &HLT_Mu0_Track0_Jpsi, &b_HLT_Mu0_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_Track0_Jpsi", &HLT_Mu3_Track0_Jpsi, &b_HLT_Mu3_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu5_Track0_Jpsi", &HLT_Mu3_Track0_Jpsi, &b_HLT_Mu3_Track0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu0_TkMu0_Jpsi", &HLT_Mu0_TkMu0_Jpsi, &b_HLT_Mu0_TkMu0_Jpsi);
   fChain->SetBranchAddress("HLT_Mu3_TkMu0_Jpsi", &HLT_Mu3_TkMu0_Jpsi, &b_HLT_Mu3_TkMu0_Jpsi);
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
   fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
   fChain->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &HLT_DoubleMu0_Quarkonium_v1, &b_HLT_DoubleMu0_Quarkonium_v1);
   fChain->SetBranchAddress("HLT_DoubleMu3", &HLT_DoubleMu3, &b_HLT_DoubleMu3);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, &b_HLT_L1DoubleMuOpen);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen_Tight", &HLT_L1DoubleMuOpen_Tight, &b_HLT_L1DoubleMuOpen_Tight);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu7", &HLT_Mu7, &b_HLT_Mu7);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Mu11", &HLT_Mu11, &b_HLT_Mu11);
   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   /* fChain->SetBranchAddress("MCType_idx", &MCType_idx, &b_MCType_idx); */
   /* fChain->SetBranchAddress("MCType_lbl", MCType_lbl, &b_MCType_lbl); */
   fChain->SetBranchAddress("costh_HX", &costh_HX, &b_costh_HX);
   fChain->SetBranchAddress("costh_CS", &costh_CS, &b_costh_CS);
   fChain->SetBranchAddress("phi_HX", &phi_HX, &b_phi_HX);
   fChain->SetBranchAddress("phi_CS", &phi_CS, &b_phi_CS);
   Notify();
}

Bool_t PolData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PolData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PolData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PolData_cxx
