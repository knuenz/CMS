#ifndef PolData_h
#define PolData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"

class PolData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *onia, *muPos, *muNeg;
   Int_t           eventNb;
   Int_t           runNb;
   Int_t           lumiBlock;
   Int_t           nPriVtx;
   Int_t           JpsiType;
 //TLorentzVector  *JpsiP;
   /* UInt_t          fUniqueID; */
   /* UInt_t          fBits; */
   /* UInt_t          fP_fUniqueID; */
   /* UInt_t          fP_fBits; */
   /* Double_t        fP_fX; */
   /* Double_t        fP_fY; */
   /* Double_t        fP_fZ; */
   /* Double_t        fE; */
   Int_t           JpsiCharge;
   Double_t        Jpsict;
   Double_t        JpsictErr;
   Double_t        JpsiVprob;
   Double_t        JpsiDistM1;
   Double_t        JpsiDphiM1;
   Double_t        JpsiDrM1;
   Double_t        JpsiDistM2;
   Double_t        JpsiDphiM2;
   Double_t        JpsiDrM2;
 //TLorentzVector  *muPosP;
   /* UInt_t          fUniqueID; */
   /* UInt_t          fBits; */
   /* UInt_t          fP_fUniqueID; */
   /* UInt_t          fP_fBits; */
   /* Double_t        fP_fX; */
   /* Double_t        fP_fY; */
   /* Double_t        fP_fZ; */
   /* Double_t        fE; */
 //TLorentzVector  *muNegP;
   /* UInt_t          fUniqueID; */
   /* UInt_t          fBits; */
   /* UInt_t          fP_fUniqueID; */
   /* UInt_t          fP_fBits; */
   /* Double_t        fP_fX; */
   /* Double_t        fP_fY; */
   /* Double_t        fP_fZ; */
   /* Double_t        fE; */
   Int_t           HLT_DoubleMu3_Jpsi_v1;
   Int_t           HLT_DoubleMu3_Jpsi_v1_PreScale;
   Int_t           HLT_DoubleMu3_Jpsi_v2;
   Int_t           HLT_DoubleMu3_Jpsi_v2_PreScale;
   Int_t           HLT_Dimuon6p5_Jpsi_v1;
   Int_t           HLT_Dimuon6p5_Jpsi_v1_PreScale;
   Int_t           HLT_Dimuon6p5_Jpsi_Displaced_v1;
   Int_t           HLT_Dimuon6p5_Jpsi_Displaced_v1_PreScale;
   Int_t           HLT_Dimuon6p5_Barrel_Jpsi_v1;
   Int_t           HLT_Dimuon6p5_Barrel_Jpsi_v1_PreScale;
   Int_t           HLT_DoubleMu3_Quarkonium_v1;
   Int_t           HLT_DoubleMu3_Quarkonium_v1_PreScale;
   Int_t           HLT_DoubleMu3_Quarkonium_v2;
   Int_t           HLT_DoubleMu3_Quarkonium_v2_PreScale;
   Int_t           HLT_Dimuon6p5_Barrel_PsiPrime_v1;
   Int_t           HLT_Dimuon6p5_Barrel_PsiPrime_v1_PreScale;
   Int_t           HLT_DoubleMu3_Upsilon_v1;
   Int_t           HLT_DoubleMu3_Upsilon_v1_PreScale;
   Int_t           HLT_Dimuon0_Barrel_Upsilon_v1;
   Int_t           HLT_Dimuon0_Barrel_Upsilon_v1_PreScale;
   Int_t           HLT_DoubleMu2_Bs_v3;
   Int_t           HLT_DoubleMu2_Bs_v3_PreScale;
   Int_t           HLT_DoubleMu2_Bs_v4;
   Int_t           HLT_DoubleMu2_Bs_v4_PreScale;
   Int_t           HLT_DoubleMu2_Bs_v5;
   Int_t           HLT_DoubleMu2_Bs_v5_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_v1;
   Int_t           HLT_Dimuon0_Jpsi_v1_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_v2;
   Int_t           HLT_Dimuon0_Jpsi_v2_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_v3;
   Int_t           HLT_Dimuon0_Jpsi_v3_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_v1;
   Int_t           HLT_Dimuon0_Upsilon_v1_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_v2;
   Int_t           HLT_Dimuon0_Upsilon_v2_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_v3;
   Int_t           HLT_Dimuon0_Upsilon_v3_PreScale;
   Int_t           HLT_Dimuon4_Bs_Barrel_v2;
   Int_t           HLT_Dimuon4_Bs_Barrel_v2_PreScale;
   Int_t           HLT_Dimuon4_Bs_Barrel_v3;
   Int_t           HLT_Dimuon4_Bs_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon4_Bs_Barrel_v4;
   Int_t           HLT_Dimuon4_Bs_Barrel_v4_PreScale;
   Int_t           HLT_Dimuon4_Bs_Barrel_v5;
   Int_t           HLT_Dimuon4_Bs_Barrel_v5_PreScale;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v1;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v1_PreScale;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v2;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v2_PreScale;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v3;
   Int_t           HLT_Dimuon5_Upsilon_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon6_Bs_v1;
   Int_t           HLT_Dimuon6_Bs_v1_PreScale;
   Int_t           HLT_Dimuon6_Bs_v2;
   Int_t           HLT_Dimuon6_Bs_v2_PreScale;
   Int_t           HLT_Dimuon6_Bs_v3;
   Int_t           HLT_Dimuon6_Bs_v3_PreScale;
   Int_t           HLT_Dimuon6_Bs_v4;
   Int_t           HLT_Dimuon6_Bs_v4_PreScale;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v1;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v1_PreScale;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v2;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v2_PreScale;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v3;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v3_PreScale;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v4;
   Int_t           HLT_Dimuon7_LowMass_Displaced_v4_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_Displaced_v1;
   Int_t           HLT_Dimuon7_Jpsi_Displaced_v1_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_Displaced_v2;
   Int_t           HLT_Dimuon7_Jpsi_Displaced_v2_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_Displaced_v3;
   Int_t           HLT_Dimuon7_Jpsi_Displaced_v3_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v1;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v1_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v2;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v2_PreScale;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v3;
   Int_t           HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon7_PsiPrime_v1;
   Int_t           HLT_Dimuon7_PsiPrime_v1_PreScale;
   Int_t           HLT_Dimuon7_PsiPrime_v2;
   Int_t           HLT_Dimuon7_PsiPrime_v2_PreScale;
   Int_t           HLT_Dimuon7_PsiPrime_v3;
   Int_t           HLT_Dimuon7_PsiPrime_v3_PreScale;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v1;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v1_PreScale;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v2;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v2_PreScale;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v3;
   Int_t           HLT_Dimuon10_Jpsi_Barrel_v3_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v1;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v1_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v2;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v2_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v3;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v3_PreScale;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v4;
   Int_t           HLT_Dimuon0_Jpsi_Muon_v4_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v1;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v1_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v2;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v2_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v3;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v3_PreScale;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v4;
   Int_t           HLT_Dimuon0_Upsilon_Muon_v4_PreScale;
   Int_t           HLT_Mu5_L2Mu2_v1;
   Int_t           HLT_Mu5_L2Mu2_v1_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v1;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v1_PreScale;
   Int_t           HLT_Mu5_L2Mu2_v2;
   Int_t           HLT_Mu5_L2Mu2_v2_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v2;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v2_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v3;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v3_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v4;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v4_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v5;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v5_PreScale;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v6;
   Int_t           HLT_Mu5_L2Mu2_Jpsi_v6_PreScale;
   Int_t           HLT_Mu3_Track3_Jpsi_v4;
   Int_t           HLT_Mu3_Track3_Jpsi_v4_PreScale;
   Int_t           HLT_Mu7_Track5_Jpsi_v1;
   Int_t           HLT_Mu7_Track5_Jpsi_v1_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v1;
   Int_t           HLT_Mu7_Track7_Jpsi_v1_PreScale;
   Int_t           HLT_Mu3_Track3_Jpsi_v5;
   Int_t           HLT_Mu3_Track3_Jpsi_v5_PreScale;
   Int_t           HLT_Mu5_Track2_Jpsi_v1;
   Int_t           HLT_Mu5_Track2_Jpsi_v1_PreScale;
   Int_t           HLT_Mu7_Track5_Jpsi_v2;
   Int_t           HLT_Mu7_Track5_Jpsi_v2_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v2;
   Int_t           HLT_Mu7_Track7_Jpsi_v2_PreScale;
   Int_t           HLT_Mu5_Track2_Jpsi_v2;
   Int_t           HLT_Mu5_Track2_Jpsi_v2_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v3;
   Int_t           HLT_Mu7_Track7_Jpsi_v3_PreScale;
   Int_t           HLT_Mu5_Track2_Jpsi_v4;
   Int_t           HLT_Mu5_Track2_Jpsi_v4_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v5;
   Int_t           HLT_Mu7_Track7_Jpsi_v5_PreScale;
   Int_t           HLT_Mu5_Track2_Jpsi_v5;
   Int_t           HLT_Mu5_Track2_Jpsi_v5_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v6;
   Int_t           HLT_Mu7_Track7_Jpsi_v6_PreScale;
   Int_t           HLT_Mu5_Track2_Jpsi_v6;
   Int_t           HLT_Mu5_Track2_Jpsi_v6_PreScale;
   Int_t           HLT_Mu7_Track7_Jpsi_v7;
   Int_t           HLT_Mu7_Track7_Jpsi_v7_PreScale;

   Int_t HLT_Dimuon5_Upsilon_Barrel_v5;
   Int_t HLT_Dimuon7_Upsilon_Barrel_v1;
   Int_t HLT_Dimuon9_Upsilon_Barrel_v1;
   Int_t HLT_Dimuon7_Upsilon_Barrel_v4;
   Int_t HLT_Dimuon9_Upsilon_Barrel_v4;


   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_nPriVtx;   //!
   TBranch        *b_JpsiType;   //!
   /* TBranch        *b_JpsiP_fUniqueID;   //! */
   /* TBranch        *b_JpsiP_fBits;   //! */
   /* TBranch        *b_JpsiP_fP_fUniqueID;   //! */
   /* TBranch        *b_JpsiP_fP_fBits;   //! */
   /* TBranch        *b_JpsiP_fP_fX;   //! */
   /* TBranch        *b_JpsiP_fP_fY;   //! */
   /* TBranch        *b_JpsiP_fP_fZ;   //! */
   /* TBranch        *b_JpsiP_fE;   //! */
   TBranch        *b_JpsiCharge;   //!
   TBranch        *b_Jpsict;   //!
   TBranch        *b_JpsictErr;   //!
   TBranch        *b_JpsiVprob;   //!
   TBranch        *b_JpsiDistM1;   //!
   TBranch        *b_JpsiDphiM1;   //!
   TBranch        *b_JpsiDrM1;   //!
   TBranch        *b_JpsiDistM2;   //!
   TBranch        *b_JpsiDphiM2;   //!
   TBranch        *b_JpsiDrM2;   //!
   /* TBranch        *b_muPosP_fUniqueID;   //! */
   /* TBranch        *b_muPosP_fBits;   //! */
   /* TBranch        *b_muPosP_fP_fUniqueID;   //! */
   /* TBranch        *b_muPosP_fP_fBits;   //! */
   /* TBranch        *b_muPosP_fP_fX;   //! */
   /* TBranch        *b_muPosP_fP_fY;   //! */
   /* TBranch        *b_muPosP_fP_fZ;   //! */
   /* TBranch        *b_muPosP_fE;   //! */
   /* TBranch        *b_muNegP_fUniqueID;   //! */
   /* TBranch        *b_muNegP_fBits;   //! */
   /* TBranch        *b_muNegP_fP_fUniqueID;   //! */
   /* TBranch        *b_muNegP_fP_fBits;   //! */
   /* TBranch        *b_muNegP_fP_fX;   //! */
   /* TBranch        *b_muNegP_fP_fY;   //! */
   /* TBranch        *b_muNegP_fP_fZ;   //! */
   /* TBranch        *b_muNegP_fE;   //! */
   TBranch        *b_HLT_DoubleMu3_Jpsi_v1;   //!
   TBranch        *b_HLT_DoubleMu3_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_DoubleMu3_Jpsi_v2;   //!
   TBranch        *b_HLT_DoubleMu3_Jpsi_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon6p5_Jpsi_v1;   //!
   TBranch        *b_HLT_Dimuon6p5_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon6p5_Jpsi_Displaced_v1;   //!
   TBranch        *b_HLT_Dimuon6p5_Jpsi_Displaced_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon6p5_Barrel_Jpsi_v1;   //!
   TBranch        *b_HLT_Dimuon6p5_Barrel_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_DoubleMu3_Quarkonium_v1;   //!
   TBranch        *b_HLT_DoubleMu3_Quarkonium_v1_PreScale;   //!
   TBranch        *b_HLT_DoubleMu3_Quarkonium_v2;   //!
   TBranch        *b_HLT_DoubleMu3_Quarkonium_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon6p5_Barrel_PsiPrime_v1;   //!
   TBranch        *b_HLT_Dimuon6p5_Barrel_PsiPrime_v1_PreScale;   //!
   TBranch        *b_HLT_DoubleMu3_Upsilon_v1;   //!
   TBranch        *b_HLT_DoubleMu3_Upsilon_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Barrel_Upsilon_v1;   //!
   TBranch        *b_HLT_Dimuon0_Barrel_Upsilon_v1_PreScale;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v3;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v3_PreScale;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v4;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v4_PreScale;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v5;   //!
   TBranch        *b_HLT_DoubleMu2_Bs_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v1;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v2;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v3;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v1;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v2;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v3;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v2;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v4;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v4_PreScale;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v5;   //!
   TBranch        *b_HLT_Dimuon4_Bs_Barrel_v5_PreScale;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v1;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v2;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon5_Upsilon_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v1;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v2;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v3;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v4;   //!
   TBranch        *b_HLT_Dimuon6_Bs_v4_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v1;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v2;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v3;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v4;   //!
   TBranch        *b_HLT_Dimuon7_LowMass_Displaced_v4_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_Displaced_v1;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_Displaced_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_Displaced_v2;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_Displaced_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_Displaced_v3;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_Displaced_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v1;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v2;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v1;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v2;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v3;   //!
   TBranch        *b_HLT_Dimuon7_PsiPrime_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v1;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v2;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v3;   //!
   TBranch        *b_HLT_Dimuon10_Jpsi_Barrel_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v1;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v2;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v3;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v4;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_Muon_v4_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v1;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v1_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v2;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v2_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v3;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v3_PreScale;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v4;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_v4_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_v1;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_v1_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_v2;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_v2_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v2_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v3;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v3_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v4;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v4_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v5;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v5_PreScale;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v6;   //!
   TBranch        *b_HLT_Mu5_L2Mu2_Jpsi_v6_PreScale;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v4;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v4_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track5_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu7_Track5_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v5;   //!
   TBranch        *b_HLT_Mu3_Track3_Jpsi_v5_PreScale;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v1;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v1_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track5_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu7_Track5_Jpsi_v2_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v2_PreScale;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v2;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v2_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v3;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v3_PreScale;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v4;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v4_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v5;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v5_PreScale;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v5;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v5_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v6;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v6_PreScale;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v6;   //!
   TBranch        *b_HLT_Mu5_Track2_Jpsi_v6_PreScale;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v7;   //!
   TBranch        *b_HLT_Mu7_Track7_Jpsi_v7_PreScale;   //!


   TBranch *b_HLT_Dimuon5_Upsilon_Barrel_v5;
   TBranch *b_HLT_Dimuon7_Upsilon_Barrel_v1;
   TBranch *b_HLT_Dimuon9_Upsilon_Barrel_v1;
   TBranch *b_HLT_Dimuon7_Upsilon_Barrel_v4;
   TBranch *b_HLT_Dimuon9_Upsilon_Barrel_v4;


   PolData(TTree *tree=0);
   virtual ~PolData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t selDimuType, Bool_t rejectCowboys, Int_t FidCuts, bool UpsMC, bool RequestTrigger);
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTree_Onia2MuMu_V8_PromptReco_v4.root");
      if (!f) {
         f = new TFile("TTree_Onia2MuMu_V8_PromptReco_v4.root");
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
   /* fChain->SetMakeClass(1); */
   onia = 0;
   muNeg = 0;
   muPos = 0;

   fChain->SetBranchAddress("JpsiP", &onia);
   fChain->SetBranchAddress("muNegP", &muNeg);
   fChain->SetBranchAddress("muPosP", &muPos);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("nPriVtx", &nPriVtx, &b_nPriVtx);
   fChain->SetBranchAddress("JpsiType", &JpsiType, &b_JpsiType);
   /* fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_JpsiP_fUniqueID); */
   /* fChain->SetBranchAddress("fBits", &fBits, &b_JpsiP_fBits); */
   /* fChain->SetBranchAddress("fP.fUniqueID", &fP_fUniqueID, &b_JpsiP_fP_fUniqueID); */
   /* fChain->SetBranchAddress("fP.fBits", &fP_fBits, &b_JpsiP_fP_fBits); */
   /* fChain->SetBranchAddress("fP.fX", &fP_fX, &b_JpsiP_fP_fX); */
   /* fChain->SetBranchAddress("fP.fY", &fP_fY, &b_JpsiP_fP_fY); */
   /* fChain->SetBranchAddress("fP.fZ", &fP_fZ, &b_JpsiP_fP_fZ); */
   /* fChain->SetBranchAddress("fE", &fE, &b_JpsiP_fE); */
   fChain->SetBranchAddress("JpsiCharge", &JpsiCharge, &b_JpsiCharge);
   fChain->SetBranchAddress("Jpsict", &Jpsict, &b_Jpsict);
   fChain->SetBranchAddress("JpsictErr", &JpsictErr, &b_JpsictErr);
   fChain->SetBranchAddress("JpsiVprob", &JpsiVprob, &b_JpsiVprob);
   fChain->SetBranchAddress("JpsiDistM1", &JpsiDistM1, &b_JpsiDistM1);
   fChain->SetBranchAddress("JpsiDphiM1", &JpsiDphiM1, &b_JpsiDphiM1);
   fChain->SetBranchAddress("JpsiDrM1", &JpsiDrM1, &b_JpsiDrM1);
   fChain->SetBranchAddress("JpsiDistM2", &JpsiDistM2, &b_JpsiDistM2);
   fChain->SetBranchAddress("JpsiDphiM2", &JpsiDphiM2, &b_JpsiDphiM2);
   fChain->SetBranchAddress("JpsiDrM2", &JpsiDrM2, &b_JpsiDrM2);
   /* fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_muPosP_fUniqueID); */
   /* fChain->SetBranchAddress("fBits", &fBits, &b_muPosP_fBits); */
   /* fChain->SetBranchAddress("fP.fUniqueID", &fP_fUniqueID, &b_muPosP_fP_fUniqueID); */
   /* fChain->SetBranchAddress("fP.fBits", &fP_fBits, &b_muPosP_fP_fBits); */
   /* fChain->SetBranchAddress("fP.fX", &fP_fX, &b_muPosP_fP_fX); */
   /* fChain->SetBranchAddress("fP.fY", &fP_fY, &b_muPosP_fP_fY); */
   /* fChain->SetBranchAddress("fP.fZ", &fP_fZ, &b_muPosP_fP_fZ); */
   /* fChain->SetBranchAddress("fE", &fE, &b_muPosP_fE); */
   /* fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_muNegP_fUniqueID); */
   /* fChain->SetBranchAddress("fBits", &fBits, &b_muNegP_fBits); */
   /* fChain->SetBranchAddress("fP.fUniqueID", &fP_fUniqueID, &b_muNegP_fP_fUniqueID); */
   /* fChain->SetBranchAddress("fP.fBits", &fP_fBits, &b_muNegP_fP_fBits); */
   /* fChain->SetBranchAddress("fP.fX", &fP_fX, &b_muNegP_fP_fX); */
   /* fChain->SetBranchAddress("fP.fY", &fP_fY, &b_muNegP_fP_fY); */
   /* fChain->SetBranchAddress("fP.fZ", &fP_fZ, &b_muNegP_fP_fZ); */
   /* fChain->SetBranchAddress("fE", &fE, &b_muNegP_fE); */
   fChain->SetBranchAddress("HLT_DoubleMu3_Jpsi_v1", &HLT_DoubleMu3_Jpsi_v1, &b_HLT_DoubleMu3_Jpsi_v1);
   fChain->SetBranchAddress("HLT_DoubleMu3_Jpsi_v1_PreScale", &HLT_DoubleMu3_Jpsi_v1_PreScale, &b_HLT_DoubleMu3_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu3_Jpsi_v2", &HLT_DoubleMu3_Jpsi_v2, &b_HLT_DoubleMu3_Jpsi_v2);
   fChain->SetBranchAddress("HLT_DoubleMu3_Jpsi_v2_PreScale", &HLT_DoubleMu3_Jpsi_v2_PreScale, &b_HLT_DoubleMu3_Jpsi_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Jpsi_v1", &HLT_Dimuon6p5_Jpsi_v1, &b_HLT_Dimuon6p5_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Jpsi_v1_PreScale", &HLT_Dimuon6p5_Jpsi_v1_PreScale, &b_HLT_Dimuon6p5_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Jpsi_Displaced_v1", &HLT_Dimuon6p5_Jpsi_Displaced_v1, &b_HLT_Dimuon6p5_Jpsi_Displaced_v1);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Jpsi_Displaced_v1_PreScale", &HLT_Dimuon6p5_Jpsi_Displaced_v1_PreScale, &b_HLT_Dimuon6p5_Jpsi_Displaced_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Barrel_Jpsi_v1", &HLT_Dimuon6p5_Barrel_Jpsi_v1, &b_HLT_Dimuon6p5_Barrel_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Barrel_Jpsi_v1_PreScale", &HLT_Dimuon6p5_Barrel_Jpsi_v1_PreScale, &b_HLT_Dimuon6p5_Barrel_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu3_Quarkonium_v1", &HLT_DoubleMu3_Quarkonium_v1, &b_HLT_DoubleMu3_Quarkonium_v1);
   fChain->SetBranchAddress("HLT_DoubleMu3_Quarkonium_v1_PreScale", &HLT_DoubleMu3_Quarkonium_v1_PreScale, &b_HLT_DoubleMu3_Quarkonium_v1_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu3_Quarkonium_v2", &HLT_DoubleMu3_Quarkonium_v2, &b_HLT_DoubleMu3_Quarkonium_v2);
   fChain->SetBranchAddress("HLT_DoubleMu3_Quarkonium_v2_PreScale", &HLT_DoubleMu3_Quarkonium_v2_PreScale, &b_HLT_DoubleMu3_Quarkonium_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Barrel_PsiPrime_v1", &HLT_Dimuon6p5_Barrel_PsiPrime_v1, &b_HLT_Dimuon6p5_Barrel_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon6p5_Barrel_PsiPrime_v1_PreScale", &HLT_Dimuon6p5_Barrel_PsiPrime_v1_PreScale, &b_HLT_Dimuon6p5_Barrel_PsiPrime_v1_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu3_Upsilon_v1", &HLT_DoubleMu3_Upsilon_v1, &b_HLT_DoubleMu3_Upsilon_v1);
   fChain->SetBranchAddress("HLT_DoubleMu3_Upsilon_v1_PreScale", &HLT_DoubleMu3_Upsilon_v1_PreScale, &b_HLT_DoubleMu3_Upsilon_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Barrel_Upsilon_v1", &HLT_Dimuon0_Barrel_Upsilon_v1, &b_HLT_Dimuon0_Barrel_Upsilon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Barrel_Upsilon_v1_PreScale", &HLT_Dimuon0_Barrel_Upsilon_v1_PreScale, &b_HLT_Dimuon0_Barrel_Upsilon_v1_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v3", &HLT_DoubleMu2_Bs_v3, &b_HLT_DoubleMu2_Bs_v3);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v3_PreScale", &HLT_DoubleMu2_Bs_v3_PreScale, &b_HLT_DoubleMu2_Bs_v3_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v4", &HLT_DoubleMu2_Bs_v4, &b_HLT_DoubleMu2_Bs_v4);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v4_PreScale", &HLT_DoubleMu2_Bs_v4_PreScale, &b_HLT_DoubleMu2_Bs_v4_PreScale);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v5", &HLT_DoubleMu2_Bs_v5, &b_HLT_DoubleMu2_Bs_v5);
   fChain->SetBranchAddress("HLT_DoubleMu2_Bs_v5_PreScale", &HLT_DoubleMu2_Bs_v5_PreScale, &b_HLT_DoubleMu2_Bs_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v1", &HLT_Dimuon0_Jpsi_v1, &b_HLT_Dimuon0_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v1_PreScale", &HLT_Dimuon0_Jpsi_v1_PreScale, &b_HLT_Dimuon0_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v2", &HLT_Dimuon0_Jpsi_v2, &b_HLT_Dimuon0_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v2_PreScale", &HLT_Dimuon0_Jpsi_v2_PreScale, &b_HLT_Dimuon0_Jpsi_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v3", &HLT_Dimuon0_Jpsi_v3, &b_HLT_Dimuon0_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v3_PreScale", &HLT_Dimuon0_Jpsi_v3_PreScale, &b_HLT_Dimuon0_Jpsi_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v1", &HLT_Dimuon0_Upsilon_v1, &b_HLT_Dimuon0_Upsilon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v1_PreScale", &HLT_Dimuon0_Upsilon_v1_PreScale, &b_HLT_Dimuon0_Upsilon_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v2", &HLT_Dimuon0_Upsilon_v2, &b_HLT_Dimuon0_Upsilon_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v2_PreScale", &HLT_Dimuon0_Upsilon_v2_PreScale, &b_HLT_Dimuon0_Upsilon_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v3", &HLT_Dimuon0_Upsilon_v3, &b_HLT_Dimuon0_Upsilon_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v3_PreScale", &HLT_Dimuon0_Upsilon_v3_PreScale, &b_HLT_Dimuon0_Upsilon_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v2", &HLT_Dimuon4_Bs_Barrel_v2, &b_HLT_Dimuon4_Bs_Barrel_v2);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v2_PreScale", &HLT_Dimuon4_Bs_Barrel_v2_PreScale, &b_HLT_Dimuon4_Bs_Barrel_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v3", &HLT_Dimuon4_Bs_Barrel_v3, &b_HLT_Dimuon4_Bs_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v3_PreScale", &HLT_Dimuon4_Bs_Barrel_v3_PreScale, &b_HLT_Dimuon4_Bs_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v4", &HLT_Dimuon4_Bs_Barrel_v4, &b_HLT_Dimuon4_Bs_Barrel_v4);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v4_PreScale", &HLT_Dimuon4_Bs_Barrel_v4_PreScale, &b_HLT_Dimuon4_Bs_Barrel_v4_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v5", &HLT_Dimuon4_Bs_Barrel_v5, &b_HLT_Dimuon4_Bs_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon4_Bs_Barrel_v5_PreScale", &HLT_Dimuon4_Bs_Barrel_v5_PreScale, &b_HLT_Dimuon4_Bs_Barrel_v5_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v1", &HLT_Dimuon5_Upsilon_Barrel_v1, &b_HLT_Dimuon5_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v1_PreScale", &HLT_Dimuon5_Upsilon_Barrel_v1_PreScale, &b_HLT_Dimuon5_Upsilon_Barrel_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v2", &HLT_Dimuon5_Upsilon_Barrel_v2, &b_HLT_Dimuon5_Upsilon_Barrel_v2);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v2_PreScale", &HLT_Dimuon5_Upsilon_Barrel_v2_PreScale, &b_HLT_Dimuon5_Upsilon_Barrel_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3", &HLT_Dimuon5_Upsilon_Barrel_v3, &b_HLT_Dimuon5_Upsilon_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3_PreScale", &HLT_Dimuon5_Upsilon_Barrel_v3_PreScale, &b_HLT_Dimuon5_Upsilon_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v1", &HLT_Dimuon6_Bs_v1, &b_HLT_Dimuon6_Bs_v1);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v1_PreScale", &HLT_Dimuon6_Bs_v1_PreScale, &b_HLT_Dimuon6_Bs_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v2", &HLT_Dimuon6_Bs_v2, &b_HLT_Dimuon6_Bs_v2);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v2_PreScale", &HLT_Dimuon6_Bs_v2_PreScale, &b_HLT_Dimuon6_Bs_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v3", &HLT_Dimuon6_Bs_v3, &b_HLT_Dimuon6_Bs_v3);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v3_PreScale", &HLT_Dimuon6_Bs_v3_PreScale, &b_HLT_Dimuon6_Bs_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v4", &HLT_Dimuon6_Bs_v4, &b_HLT_Dimuon6_Bs_v4);
   fChain->SetBranchAddress("HLT_Dimuon6_Bs_v4_PreScale", &HLT_Dimuon6_Bs_v4_PreScale, &b_HLT_Dimuon6_Bs_v4_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v1", &HLT_Dimuon7_LowMass_Displaced_v1, &b_HLT_Dimuon7_LowMass_Displaced_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v1_PreScale", &HLT_Dimuon7_LowMass_Displaced_v1_PreScale, &b_HLT_Dimuon7_LowMass_Displaced_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v2", &HLT_Dimuon7_LowMass_Displaced_v2, &b_HLT_Dimuon7_LowMass_Displaced_v2);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v2_PreScale", &HLT_Dimuon7_LowMass_Displaced_v2_PreScale, &b_HLT_Dimuon7_LowMass_Displaced_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v3", &HLT_Dimuon7_LowMass_Displaced_v3, &b_HLT_Dimuon7_LowMass_Displaced_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v3_PreScale", &HLT_Dimuon7_LowMass_Displaced_v3_PreScale, &b_HLT_Dimuon7_LowMass_Displaced_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v4", &HLT_Dimuon7_LowMass_Displaced_v4, &b_HLT_Dimuon7_LowMass_Displaced_v4);
   fChain->SetBranchAddress("HLT_Dimuon7_LowMass_Displaced_v4_PreScale", &HLT_Dimuon7_LowMass_Displaced_v4_PreScale, &b_HLT_Dimuon7_LowMass_Displaced_v4_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced_v1", &HLT_Dimuon7_Jpsi_Displaced_v1, &b_HLT_Dimuon7_Jpsi_Displaced_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced_v1_PreScale", &HLT_Dimuon7_Jpsi_Displaced_v1_PreScale, &b_HLT_Dimuon7_Jpsi_Displaced_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced_v2", &HLT_Dimuon7_Jpsi_Displaced_v2, &b_HLT_Dimuon7_Jpsi_Displaced_v2);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced_v2_PreScale", &HLT_Dimuon7_Jpsi_Displaced_v2_PreScale, &b_HLT_Dimuon7_Jpsi_Displaced_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced_v3", &HLT_Dimuon7_Jpsi_Displaced_v3, &b_HLT_Dimuon7_Jpsi_Displaced_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_Displaced_v3_PreScale", &HLT_Dimuon7_Jpsi_Displaced_v3_PreScale, &b_HLT_Dimuon7_Jpsi_Displaced_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v1", &HLT_Dimuon7_Jpsi_X_Barrel_v1, &b_HLT_Dimuon7_Jpsi_X_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v1_PreScale", &HLT_Dimuon7_Jpsi_X_Barrel_v1_PreScale, &b_HLT_Dimuon7_Jpsi_X_Barrel_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v2", &HLT_Dimuon7_Jpsi_X_Barrel_v2, &b_HLT_Dimuon7_Jpsi_X_Barrel_v2);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v2_PreScale", &HLT_Dimuon7_Jpsi_X_Barrel_v2_PreScale, &b_HLT_Dimuon7_Jpsi_X_Barrel_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v3", &HLT_Dimuon7_Jpsi_X_Barrel_v3, &b_HLT_Dimuon7_Jpsi_X_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale", &HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale, &b_HLT_Dimuon7_Jpsi_X_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v1", &HLT_Dimuon7_PsiPrime_v1, &b_HLT_Dimuon7_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v1_PreScale", &HLT_Dimuon7_PsiPrime_v1_PreScale, &b_HLT_Dimuon7_PsiPrime_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v2", &HLT_Dimuon7_PsiPrime_v2, &b_HLT_Dimuon7_PsiPrime_v2);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v2_PreScale", &HLT_Dimuon7_PsiPrime_v2_PreScale, &b_HLT_Dimuon7_PsiPrime_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v3", &HLT_Dimuon7_PsiPrime_v3, &b_HLT_Dimuon7_PsiPrime_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v3_PreScale", &HLT_Dimuon7_PsiPrime_v3_PreScale, &b_HLT_Dimuon7_PsiPrime_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v1", &HLT_Dimuon10_Jpsi_Barrel_v1, &b_HLT_Dimuon10_Jpsi_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v1_PreScale", &HLT_Dimuon10_Jpsi_Barrel_v1_PreScale, &b_HLT_Dimuon10_Jpsi_Barrel_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v2", &HLT_Dimuon10_Jpsi_Barrel_v2, &b_HLT_Dimuon10_Jpsi_Barrel_v2);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v2_PreScale", &HLT_Dimuon10_Jpsi_Barrel_v2_PreScale, &b_HLT_Dimuon10_Jpsi_Barrel_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3, &b_HLT_Dimuon10_Jpsi_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3_PreScale", &HLT_Dimuon10_Jpsi_Barrel_v3_PreScale, &b_HLT_Dimuon10_Jpsi_Barrel_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v1", &HLT_Dimuon0_Jpsi_Muon_v1, &b_HLT_Dimuon0_Jpsi_Muon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v1_PreScale", &HLT_Dimuon0_Jpsi_Muon_v1_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v2", &HLT_Dimuon0_Jpsi_Muon_v2, &b_HLT_Dimuon0_Jpsi_Muon_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v2_PreScale", &HLT_Dimuon0_Jpsi_Muon_v2_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v3", &HLT_Dimuon0_Jpsi_Muon_v3, &b_HLT_Dimuon0_Jpsi_Muon_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v3_PreScale", &HLT_Dimuon0_Jpsi_Muon_v3_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v4", &HLT_Dimuon0_Jpsi_Muon_v4, &b_HLT_Dimuon0_Jpsi_Muon_v4);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon_v4_PreScale", &HLT_Dimuon0_Jpsi_Muon_v4_PreScale, &b_HLT_Dimuon0_Jpsi_Muon_v4_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v1", &HLT_Dimuon0_Upsilon_Muon_v1, &b_HLT_Dimuon0_Upsilon_Muon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v1_PreScale", &HLT_Dimuon0_Upsilon_Muon_v1_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v1_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v2", &HLT_Dimuon0_Upsilon_Muon_v2, &b_HLT_Dimuon0_Upsilon_Muon_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v2_PreScale", &HLT_Dimuon0_Upsilon_Muon_v2_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v2_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v3", &HLT_Dimuon0_Upsilon_Muon_v3, &b_HLT_Dimuon0_Upsilon_Muon_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v3_PreScale", &HLT_Dimuon0_Upsilon_Muon_v3_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v3_PreScale);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v4", &HLT_Dimuon0_Upsilon_Muon_v4, &b_HLT_Dimuon0_Upsilon_Muon_v4);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_v4_PreScale", &HLT_Dimuon0_Upsilon_Muon_v4_PreScale, &b_HLT_Dimuon0_Upsilon_Muon_v4_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v1", &HLT_Mu5_L2Mu2_v1, &b_HLT_Mu5_L2Mu2_v1);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v1_PreScale", &HLT_Mu5_L2Mu2_v1_PreScale, &b_HLT_Mu5_L2Mu2_v1_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v1", &HLT_Mu5_L2Mu2_Jpsi_v1, &b_HLT_Mu5_L2Mu2_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v1_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v1_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v2", &HLT_Mu5_L2Mu2_v2, &b_HLT_Mu5_L2Mu2_v2);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_v2_PreScale", &HLT_Mu5_L2Mu2_v2_PreScale, &b_HLT_Mu5_L2Mu2_v2_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v2", &HLT_Mu5_L2Mu2_Jpsi_v2, &b_HLT_Mu5_L2Mu2_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v2_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v2_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v2_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v3", &HLT_Mu5_L2Mu2_Jpsi_v3, &b_HLT_Mu5_L2Mu2_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v3_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v3_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v3_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v4", &HLT_Mu5_L2Mu2_Jpsi_v4, &b_HLT_Mu5_L2Mu2_Jpsi_v4);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v4_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v4_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v4_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v5", &HLT_Mu5_L2Mu2_Jpsi_v5, &b_HLT_Mu5_L2Mu2_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v5_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v5_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v5_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v6", &HLT_Mu5_L2Mu2_Jpsi_v6, &b_HLT_Mu5_L2Mu2_Jpsi_v6);
   fChain->SetBranchAddress("HLT_Mu5_L2Mu2_Jpsi_v6_PreScale", &HLT_Mu5_L2Mu2_Jpsi_v6_PreScale, &b_HLT_Mu5_L2Mu2_Jpsi_v6_PreScale);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v4", &HLT_Mu3_Track3_Jpsi_v4, &b_HLT_Mu3_Track3_Jpsi_v4);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v4_PreScale", &HLT_Mu3_Track3_Jpsi_v4_PreScale, &b_HLT_Mu3_Track3_Jpsi_v4_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v1", &HLT_Mu7_Track5_Jpsi_v1, &b_HLT_Mu7_Track5_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v1_PreScale", &HLT_Mu7_Track5_Jpsi_v1_PreScale, &b_HLT_Mu7_Track5_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v1", &HLT_Mu7_Track7_Jpsi_v1, &b_HLT_Mu7_Track7_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v1_PreScale", &HLT_Mu7_Track7_Jpsi_v1_PreScale, &b_HLT_Mu7_Track7_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v5", &HLT_Mu3_Track3_Jpsi_v5, &b_HLT_Mu3_Track3_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Mu3_Track3_Jpsi_v5_PreScale", &HLT_Mu3_Track3_Jpsi_v5_PreScale, &b_HLT_Mu3_Track3_Jpsi_v5_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v1", &HLT_Mu5_Track2_Jpsi_v1, &b_HLT_Mu5_Track2_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v1_PreScale", &HLT_Mu5_Track2_Jpsi_v1_PreScale, &b_HLT_Mu5_Track2_Jpsi_v1_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v2", &HLT_Mu7_Track5_Jpsi_v2, &b_HLT_Mu7_Track5_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu7_Track5_Jpsi_v2_PreScale", &HLT_Mu7_Track5_Jpsi_v2_PreScale, &b_HLT_Mu7_Track5_Jpsi_v2_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v2", &HLT_Mu7_Track7_Jpsi_v2, &b_HLT_Mu7_Track7_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v2_PreScale", &HLT_Mu7_Track7_Jpsi_v2_PreScale, &b_HLT_Mu7_Track7_Jpsi_v2_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v2", &HLT_Mu5_Track2_Jpsi_v2, &b_HLT_Mu5_Track2_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v2_PreScale", &HLT_Mu5_Track2_Jpsi_v2_PreScale, &b_HLT_Mu5_Track2_Jpsi_v2_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v3", &HLT_Mu7_Track7_Jpsi_v3, &b_HLT_Mu7_Track7_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v3_PreScale", &HLT_Mu7_Track7_Jpsi_v3_PreScale, &b_HLT_Mu7_Track7_Jpsi_v3_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v4", &HLT_Mu5_Track2_Jpsi_v4, &b_HLT_Mu5_Track2_Jpsi_v4);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v4_PreScale", &HLT_Mu5_Track2_Jpsi_v4_PreScale, &b_HLT_Mu5_Track2_Jpsi_v4_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v5", &HLT_Mu7_Track7_Jpsi_v5, &b_HLT_Mu7_Track7_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v5_PreScale", &HLT_Mu7_Track7_Jpsi_v5_PreScale, &b_HLT_Mu7_Track7_Jpsi_v5_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v5", &HLT_Mu5_Track2_Jpsi_v5, &b_HLT_Mu5_Track2_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v5_PreScale", &HLT_Mu5_Track2_Jpsi_v5_PreScale, &b_HLT_Mu5_Track2_Jpsi_v5_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v6", &HLT_Mu7_Track7_Jpsi_v6, &b_HLT_Mu7_Track7_Jpsi_v6);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v6_PreScale", &HLT_Mu7_Track7_Jpsi_v6_PreScale, &b_HLT_Mu7_Track7_Jpsi_v6_PreScale);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v6", &HLT_Mu5_Track2_Jpsi_v6, &b_HLT_Mu5_Track2_Jpsi_v6);
   fChain->SetBranchAddress("HLT_Mu5_Track2_Jpsi_v6_PreScale", &HLT_Mu5_Track2_Jpsi_v6_PreScale, &b_HLT_Mu5_Track2_Jpsi_v6_PreScale);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v7", &HLT_Mu7_Track7_Jpsi_v7, &b_HLT_Mu7_Track7_Jpsi_v7);
   fChain->SetBranchAddress("HLT_Mu7_Track7_Jpsi_v7_PreScale", &HLT_Mu7_Track7_Jpsi_v7_PreScale, &b_HLT_Mu7_Track7_Jpsi_v7_PreScale);

   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v5", &HLT_Dimuon5_Upsilon_Barrel_v5,&b_HLT_Dimuon5_Upsilon_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v1", &HLT_Dimuon7_Upsilon_Barrel_v1,&b_HLT_Dimuon7_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v1", &HLT_Dimuon9_Upsilon_Barrel_v1,&b_HLT_Dimuon9_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v4", &HLT_Dimuon7_Upsilon_Barrel_v4,&b_HLT_Dimuon7_Upsilon_Barrel_v4);
   fChain->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v4", &HLT_Dimuon9_Upsilon_Barrel_v4,&b_HLT_Dimuon9_Upsilon_Barrel_v4);


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


