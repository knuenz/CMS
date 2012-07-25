#define PolData_cxx
#include "PolData.h"
#include "commonVar.h"
//#include "isMuonInAcceptance.C"
#include "effsAndCuts.h"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

//some statistics
TH1F *Reco_StatEv;
TH1F *Reco_Onia_mass[onia::kNbPTMaxBins+1][onia::kNbRapForPTBins+1];
TH2F *Reco_Onia_rap_pT;
TH1F *Reco_Onia_pt[onia::kNbRapForPTBins+1];
TH1F *Reco_Onia_rap[onia::kNbPTMaxBins+1];
TTree *treeOut;
TLorentzVector *lepP, *lepN;
Int_t RunID;

//polarization histos:
// TH2F *Reco2D_Onia_pol_pT_rap[onia::kNbFrames][onia::kNbPTMaxBins+1][onia::kNbRapForPTBins+1];

enum {LOOSE,TIGHT};//set of muon fiducial cuts
//==============================================
void PolData::Loop(Int_t selDimuType, Bool_t rejectCowboys, Int_t FidCuts, bool UpsMC, bool RequestTrigger, bool selectSOFT, bool selectTIGHT, bool selectMIXED, bool selectNOTMIXED)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  Long64_t cutAtRecEvent = nentries;
  Long64_t countRecEvent = 0;
  Long64_t nb = 0;
  printf("number of entries = %d\n", (Int_t) nentries);

//  nentries=2000000;
  //loop over the events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<1000000;jentry++) {
//	   cout<<onia::kNbRapForPTBins<<endl;

    if(jentry % 100000 == 0) printf("event %d\n", (Int_t) jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    //if we process MC, we must ensure that we only consider
    //reconstructed events
    if(onia->Pt() > 990.)
      continue;

    if(JpsiVprob < 0.01)
      continue;

//		double cut_nPriVtx=7.5;
//		if(nPriVtx>cut_nPriVtx) continue;

//    Run1: 165088 - 172868 (1a + 2a) --> 53%
//    Run2: 175971 - 178379 (1b + 2a) --> 23%
//    Run3: 178380 - 180252 (1b + 2b) --> 15%

    if(runNb>=165088 && runNb<=172868) RunID=1;
    if(runNb>=175971 && runNb<=178379) RunID=2;
    if(runNb>=178380 && runNb<=180252) RunID=3;

    //if(RunID!=3 && onia->Pt() > 10.) continue;


/*    cout<<"muPosPglobalchi2"<<			muPosPglobalchi2           <<endl;
    cout<<"muNegPglobalchi2"<<			muNegPglobalchi2           <<endl;
    cout<<"muPosPglobalMuonHits"<<		muPosPglobalMuonHits       <<endl;
    cout<<"muNegPglobalMuonHits"<<		muNegPglobalMuonHits       <<endl;
    cout<<"muPosPMuonMatchedStations"<<	muPosPMuonMatchedStations   <<endl;
    cout<<"muNegPMuonMatchedStations"<<	muNegPMuonMatchedStations   <<endl;
    cout<<"ismuPosTMOneStationTight"<<	ismuPosTMOneStationTight   <<endl;
    cout<<"ismuNegTMOneStationTight"<<	ismuNegTMOneStationTight   <<endl;
*/
    if(selectTIGHT){
    if(
      !(
          (muPosPglobalchi2>-1 && muNegPglobalchi2>-1) &&
          (muPosPglobalchi2<10 && muNegPglobalchi2<10)&&
          (muPosPglobalMuonHits>0 && muNegPglobalMuonHits>0)
//          && (muPosPMuonMatchedStations>1, muNegPMuonMatchedStations>1)
      )
    ) continue;
    if(
      !(
          ismuPosTMOneStationTight==1 && ismuNegTMOneStationTight==1
      )
    ) continue;
//    if(JpsiDrM1 < 0.5) continue;
    }

    if(selectSOFT){
    if(
      !(
          ismuPosTMOneStationTight==1 && ismuNegTMOneStationTight==1
      )
    ) continue;
    }

    if(selectMIXED){
    	if(!(ismuPosTMOneStationTight==1&&ismuNegTMOneStationTight==1))
    	  continue;
    	if(muPosPglobalMuonHits>0)
    	  if( !(muPosPglobalchi2<10) )
    	    continue;
    	if(muNegPglobalMuonHits>0)
    	  if( !(muNegPglobalchi2<10) )
    	    continue;
    	if(muPosPglobalMuonHits<=0) continue;
    	if(muNegPglobalMuonHits<=0) continue;

    }

    if(selectNOTMIXED){
    	if(!(ismuPosTMOneStationTight==1&&ismuNegTMOneStationTight==1))
    	    continue;
      if((muPosPglobalMuonHits>0 && muPosPglobalchi2<10) && (muNegPglobalMuonHits>0 && muNegPglobalchi2<10))
        continue;
    }

//    if(muPosPglobalOK!=1||muNegPglobalOK!=1) continue;

    Reco_StatEv->Fill(0.5);//count all events

    //check the trigger flag: 0... no trigger, 1 ... triggered+matched, 3 ... triggered (HLT_DoubleMu0)
    //for a full list of accessible triggers, check out "PolData.h"
    //and https://espace.cern.ch/cms-quarkonia/onia-polarization/L1%20%20HLT/unprescaledTriggersVsRun.aspx
    //for the run ranges per HLT path during unprescaled running periods

    Int_t trigDecision = -99;
    Int_t trigPtDecision = -99;
    
/*    if(//HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
       //HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
       //HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
       HLT_Dimuon5_Upsilon_Barrel_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
       HLT_Dimuon5_Upsilon_Barrel_v2 == 1 || //  1e33: 166346
       HLT_Dimuon5_Upsilon_Barrel_v3 == 1)   //1.4e33: 167078 - 167913
      trigDecision = 1;
*/

  //Upsilon trigger paths:
   if(//HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
      //HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
      //HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
      HLT_Dimuon5_Upsilon_Barrel_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
      HLT_Dimuon5_Upsilon_Barrel_v2 == 1 || //  1e33: 166346
      HLT_Dimuon5_Upsilon_Barrel_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 2)
      HLT_Dimuon5_Upsilon_Barrel_v5 == 1 || //2E33 (no cowboys)
      HLT_Dimuon7_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ)
      HLT_Dimuon7_Upsilon_Barrel_v4 == 1 || //5E33 (becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v4 == 1) //5E33
      trigDecision = 1;

   if(trigDecision != 1 && RequestTrigger)
     continue;


   if(//HLT_DoubleMu3_Quarkonium_v1 == 1 ||   //  5e32: 160404 - 161176
      //HLT_DoubleMu3_Upsilon_v1 == 1 ||      //  5e32: 161216 - 163261
      //HLT_Dimuon0_Barrel_Upsilon_v1 == 1 || //  5e32: 163269 - 163869
      HLT_Dimuon5_Upsilon_Barrel_v1 == 1 && onia->Pt() < 5.5 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
      HLT_Dimuon5_Upsilon_Barrel_v2 == 1 && onia->Pt() < 5.5  || //  1e33: 166346
      HLT_Dimuon5_Upsilon_Barrel_v3 == 1 && onia->Pt() < 5.5  ||  //1.4e33: 167078 - 167913 (prescale of 2)
      HLT_Dimuon5_Upsilon_Barrel_v5 == 1 && onia->Pt() < 5.5  || //2E33 (no cowboys)
      HLT_Dimuon7_Upsilon_Barrel_v1 == 1 && onia->Pt() < 7.5  || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v1 == 1 && onia->Pt() < 9.5  || //3E33 (L1_DoubleMu0_HighQ)
      HLT_Dimuon7_Upsilon_Barrel_v4 == 1 && onia->Pt() < 7.5  || //5E33 (becomes inactive for Linst >= 5E33)
      HLT_Dimuon9_Upsilon_Barrel_v4 == 1 && onia->Pt() < 9.5 ) //5E33
      trigPtDecision = 1;

   if(trigPtDecision == 1)
     continue;




    //reject processing of events where the dimuon type (GG, GT or TT)
    //does not correspond to the chosen one
    // if(selDimuType < 3 && JpsiType != selDimuType)
    //   continue;
    // else if(selDimuType == 3 && JpsiType > 1) //only GG or GT
    //   continue;

    Reco_StatEv->Fill(1.5); //count all events

    Double_t onia_mass = onia->M();
    Double_t onia_pt = onia->Pt();
    Double_t onia_P = onia->P();
    Double_t onia_eta = onia->PseudoRapidity();
    Double_t onia_rap = onia->Rapidity();
    Double_t onia_phi = onia->Phi();
    Double_t onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

    if(TMath::Abs(onia_rap) > onia::rapYPS)
      continue;
    Reco_StatEv->Fill(2.5);

    Double_t etaMuPos = muPos->PseudoRapidity();
    Double_t etaMuNeg = muNeg->PseudoRapidity();
    Double_t pTMuPos = muPos->Pt();
    Double_t pTMuNeg = muNeg->Pt();
    Double_t pMuPos = muPos->P();
    Double_t pMuNeg = muNeg->P();


	Double_t deltaPhi = muNeg->Phi() - muPos->Phi();
	if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
	else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

	if(rejectCowboys)
		  if(deltaPhi < 0.) continue;

	//if(rejectCowboys)
	//	if((muNeg->Phi() - muPos->Phi()) < 0.)
	//		continue;

    Reco_StatEv->Fill(3.5);

    //select events with a cut on the lifetime to reject NP J/psis:
    if(Jpsict > onia::JpsiCtauMax)
      continue;
    if(Jpsict < onia::JpsiCtauMin)
      continue;

    if(TMath::Abs(Jpsict/JpsictErr)<onia::JpsiCtauSigMin) continue;
    if(TMath::Abs(Jpsict/JpsictErr)>onia::JpsiCtauSigMax) continue;

    Reco_StatEv->Fill(4.5);

    if(onia_mass < 8.4 || onia_mass > 11.6) //all Upsilons triggered
      continue;
    // if(onia_mass < (9.45 - 2.*0.080) || onia_mass > (9.45 + 2.*0.080)) //Ups(1S) only
    //   continue;

    Reco_StatEv->Fill(5.5);

    //apply different fiducial cuts for the muon matched to the HLT leg or the TM leg
    //for simplicity we decided to apply the stronger cuts to the higher pT muon and not
    //strictly to the HLT muon
    Bool_t muonsInAcc = kFALSE;
    if(isMuonInAcceptance(FidCuts-1, pTMuPos, etaMuPos) && isMuonInAcceptance(FidCuts-1, pTMuNeg, etaMuNeg))
      muonsInAcc = kTRUE;
    if(!muonsInAcc)
      continue;

   Reco_StatEv->Fill(5.5);

   //reject furthermore all events in which one of the muons has a pT smaller than 3 GeV/c
   // if(pTMuPos < 3.0 || pTMuNeg < 3.0)
   //   continue;

   Reco_StatEv->Fill(6.5);

   //set up the pT and y indices
   Int_t rapForPTIndex = -1;
   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
     if(TMath::Abs(onia_rap) > onia::rapForPTRange[iRap] && 
	TMath::Abs(onia_rap) < onia::rapForPTRange[iRap+1]){
       rapForPTIndex = iRap+1;
       break;
     }
   }
    Int_t pTIndex = -1;
    for(int iPT = 0; iPT < onia::kNbPTBins[rapForPTIndex]; iPT++){
      if(onia_pt > onia::pTRange[rapForPTIndex][iPT] && onia_pt < onia::pTRange[rapForPTIndex][iPT+1]){
	pTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIntegratedPTIndex = -1;
    for(int iPT = 0; iPT < onia::kNbPTBins[0]; iPT++){
      if(onia_pt > onia::pTRange[0][iPT] && onia_pt < onia::pTRange[0][iPT+1]){
	rapIntegratedPTIndex = iPT+1;
	break;
      }
    }
    Int_t rapIndex = -1;
    for(int iRap = 0; iRap < 2*onia::kNbRapBins; iRap++){
      if(onia_rap > onia::rapRange[iRap] && onia_rap < onia::rapRange[iRap+1]){
    	rapIndex = iRap;
    	break;
      }
    }

    Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
    Reco_Onia_mass[0][0]->Fill(onia_mass);
    Reco_Onia_rap[0]->Fill(onia_rap);
    Reco_Onia_pt[0]->Fill(onia_pt);

    if(rapIntegratedPTIndex >= 0){
      Reco_Onia_mass[rapIntegratedPTIndex][0]->Fill(onia_mass);
      Reco_Onia_rap[rapIntegratedPTIndex]->Fill(onia_rap);
    }
    if(rapForPTIndex > 0){
      Reco_Onia_mass[0][rapForPTIndex]->Fill(onia_mass);
      Reco_Onia_pt[rapForPTIndex]->Fill(onia_pt);
    }
    //continue only if we have events within the bins we are interested in
    if(rapIndex < 0){
      // printf("rapIndex %d, rap(onia) = %f\n", rapIndex, onia_rap);
      continue;
    }
    if(rapForPTIndex < 1){
      // printf("rapForPTIndex %d, rap(onia) = %f\n", rapForPTIndex, onia_rap);
      continue;
    }
    if(pTIndex < 1){
      // printf("pTIndex %d, pT(onia) = %f\n", pTIndex, onia_pt);
      continue;
    }

    Reco_Onia_mass[pTIndex][rapForPTIndex]->Fill(onia_mass);

/*    if(Jpsict > onia::JpsiCtauMax)
      continue;
    if(Jpsict < onia::JpsiCtauMin)
      continue;
*/
//    if(TMath::Abs(Jpsict/JpsictErr)<3.) continue;
//    if(TMath::Abs(Jpsict/JpsictErr)>onia::JpsiCtauSigMax) continue;

    //remaining of the events will be used for the analysis
    countRecEvent++;

    lepP = muPos;
    lepN = muNeg;
    treeOut->Fill();
  }
  printf("nb. of rec. events is %d of a total of %d events\n", (Int_t) countRecEvent, (Int_t) nentries);
}
