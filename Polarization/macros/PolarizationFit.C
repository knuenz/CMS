gROOT->Reset();
#include "RootIncludes.h"
Int_t SavePlotIndex = 1; //0...don't save plots, 1... save plots (.png)
Int_t RealDataIndex = 1; //0...fit to MC, 1...fit to Data
Int_t PromptIndex = 1; //1...fit Prompt Jpsi, 0...fit Non Prompt Jpsi
Int_t AngularIndex = 1; //0...no polarization fits
Char_t *fileName2 = "RootInput/RooDataSet_pol_Mu0Track0Jpsi_dataR_1_600.root"; //analyser_files/DataSet_OniaCS.root"; // Data-Sample
Char_t *fileName3 = "RootInput/PolTree.root";//TNtuple with pol vars (MC Model Sample)
Char_t *fileName4 = "RootInput/PolTreeTest.root";//TNtuple with pol vars (MC Test Sample)
Char_t *fileNamePrompt = "RootInput/TTree_pol_Mu0Track0Jpsi_MCpromptR.root";//RooDataSet Prompt MC
Char_t *fileNameMinBias = "RootInput/TTree_pol_Mu0Track0Jpsi_MCMinBias.root";//RooDataSet Background and NonPrompt MC
Char_t *fileNameAcceptance = "RootInput/accHistos_HLT_Mu0Track0Jpsi_30Aug2010.root";//RootFile Acceptance maps


#include "GlobalVars.h"

void ReadDataSets();
void Legend();
void DrawAndSave();
void DrawFinal();
void DefineMassModel();
void DefinePromptDecayLengthModel();
void ReduceDataSetToBin();
void FitMass();
void PreFitPrompt();
void ProduceSideBandTemplates();
void AnalyticalBKandNPModels();
void CalcSigOvBkg();
void FitGlobalAndAnalytical();
void Plot();
void FillArrays();
void FitPol();
void CalcAndPrintResults();
void Chi2();



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void PolarizationFit(){

	ReadDataSets();

	Legend();

	DefineMassModel();

	DefinePromptDecayLengthModel();


	char outputfilename2[200];
	sprintf(outputfilename2,"Results/results_PolarizationFit.txt");
	printf("output filename is: %s\n", outputfilename2);
	outputFile2 = fopen(outputfilename2,"w");

	fprintf(outputFile2, "Number of bins in pT = %1.0f\n",numberofpTbins);
	fprintf(outputFile2, "Number of bins in Rapidity = %1.0f\n",numberofRapbins);
	fprintf(outputFile2, "Number of bins in costh and phi = %1.0f and %1.0f\n\n",Angular_Binning_costh,Angular_Binning_phi);

		numberofbins=numberofpTbins*numberofRapbins;
		nbin = 1;
		npT=NpTBinStart;

		while(npT<=numberofpTbins+NpTBinStart-1) {

			nRap=NRapBinStart;

			pTMin=pTRange[npT-1];
			pTMax=pTRange[npT];

			while(nRap<=numberofRapbins+NRapBinStart-1) {

				RapMin=rapForPTRange[nRap-1];
				RapMax=rapForPTRange[nRap];


	ReduceDataSetToBin();

	FitMass();

	CalcSigOvBkg();

	PreFitPrompt();

	ProduceSideBandTemplates();

	AnalyticalBKandNPModels();

	FitGlobalAndAnalytical();

	CalcAndPrintResults();

	FitPol();

	Plot();

	Chi2();

	DrawAndSave();

	FillArrays();


				nbin++;
				nRap++;
			}

			npT++;

		}

	fclose(outputFile2);

	DrawFinal();


	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void ReadDataSets(){

	fAcc = new TFile(fileNameAcceptance);
	   Char_t name[100];
	   for(int iFrame = 0; iFrame < 2; iFrame++){
	   for(int iPTBin = 1; iPTBin < kNbPTBins+1; iPTBin++){
	     for(int iRapBin = 1; iRapBin < kNbRapForPTBins+1; iRapBin++){
	       sprintf(name, "hAcc2D_Onia_%s_pT%d_rap%d", frameLabel[iFrame], iPTBin, iRapBin);
	       hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin] = (TH2F *) fAcc->Get(name);//gDirectory->Get(name);
//	       hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Rebin2D(rebinCosTh, rebinPhi);
	       //normalise to 1:
//	       Double_t integral = hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->GetIntegral();
//	      hAcc2D_pol_pT_rap[iFrame][iPTBin][iRapBin]->Scale(1./integral);
	  }
	}
	   }


	char reduceMC[200];

	fIn2 = new TFile(fileName2);
	realdata000 = (RooDataSet*)fIn2->Get("data");
	realdata000->Print();


	RooBinning rbcosth(-1,1);
	rbcosth.addUniform(Angular_Binning_costh,-1,1);
	costh_CS.setBinning(rbcosth);

	RooBinning rbphi(0,360);
	rbphi.addUniform(Angular_Binning_phi,0,360);
	phi_CS.setBinning(rbphi);

//	costh_HX.setBinning(rbcosth);
//	phi_HX.setBinning(rbphi);



	fIn3 = new TFile(fileName3);
	fIn4 = new TFile(fileName4);


	treeData = (TTree*)fIn3->Get("PolTree");
	treeDataTest = (TTree*)fIn4->Get("PolTreeTest");


	RooArgSet varlist(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX,MCweight);
	varlist.add(MCType_idx);
	varlist.add(JpsiType_idx);

	RooArgSet varlist2(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX);
	varlist2.add(MCType_idx);
	varlist2.add(JpsiType_idx);


//	MCPROMPTALLVar = (RooDataSet*)fInPrompt->Get("data");
//	MCMINBIASALLVar = (RooDataSet*)fInMinBias->Get("data");

//	TTree* MCPROMPTALLVartree = MCPROMPTALLVar->tree();
//	TTree* MCMINBIASALLVartree = MCMINBIASALLVar->tree();

	if (ValIndex == 0){

	fInPrompt = new TFile(fileNamePrompt);
	fInMinBias = new TFile(fileNameMinBias);

	MCPROMPTALLVartree = (TTree*)fInPrompt->Get("data");
	MCMINBIASALLVartree = (TTree*)fInMinBias->Get("data");

	MCPROMPT = new RooDataSet("MCPROMPT","MCPROMPT",MCPROMPTALLVartree,varlist2,0);
	MCMINBIAS = new RooDataSet("MCMINBIAS","MCMINBIAS",MCMINBIASALLVartree,varlist2,0);

	MCNONPROMPT = (RooDataSet*)MCMINBIAS->reduce("MCType_idx > 1.5");
	MCBACKGROUND = (RooDataSet*)MCMINBIAS->reduce("MCType_idx > 0.5 && MCType_idx < 1.5");
	}
	else if (ValIndex == 1){

	fInPromptVal = new TFile(fileNamePromptVal);
	fInNonPromptVal = new TFile(fileNameNonPromptVal);
	fInBackgroundVal = new TFile(fileNameBackgroundVal);

	MCPROMPT = (RooDataSet*)fInPromptVal->Get("MCPROMPT");
	MCNONPROMPT = (RooDataSet*)fInNonPromptVal->Get("MCMINBIAS");
	MCBACKGROUND = (RooDataSet*)fInNonPromptVal->Get("MCMINBIAS");

	}

	MCMINBIAS->Print();
	MCPROMPT->Print();
	MCNONPROMPT->Print();
	MCBACKGROUND->Print();

	MCModelDataSet = new RooDataSet("MCModelDataSet","MCModelDataSet",treeData,varlist,0,"MCweight");//RooArgSet(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX,MCweight,MCType_idx,JpsiType_idx));
	MCModelDataSet->Print();
	MCTestDataSet = new RooDataSet("MCTestDataSet","MCTestDataSet",treeDataTest,varlist,0,"MCweight");//RooArgSet(JpsiMass,Jpsict,JpsiPt,JpsiRap,costh_CS,phi_CS,costh_HX,phi_HX,MCweight,MCType_idx,JpsiType_idx));
	MCTestDataSet->Print();


	sprintf(reduceMassData,"JpsiMass > %f && JpsiMass < %f && Jpsict < %f && Jpsict > %f && JpsiType < 1.5",JpsiMassMin,JpsiMassMax,JpsictMax,JpsictMin);
	realdata00 = (RooDataSet*)realdata000->reduce(reduceMassData);
	realdata00->Print();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DefineMassModel(){

	RooRealVar lambda("lambda","lambda",-1,-3,0.1);
	RooRealVar NbkgMass("NbkgMass","NbkgMass",1,0,150000);

	RooExponential expo("expo","expo",JpsiMass,lambda);
	RooExtendPdf eexpo("eexpo","eexpo",expo,NbkgMass);

	RooRealVar CBm("CBm","CBm",3.1,3.05,3.15);
	RooRealVar CBs("CBs","CBs",0.02,0.008,0.1);
	RooRealVar CBa("CBa","CBa",0.5,0.,3);
	RooRealVar CBn("CBn","CBn",10,1,30);

	RooCBShape CBMass("CBMass","CBMass",JpsiMass,CBm,CBs,CBa,CBn);
	RooGaussian gauss("gauss","gauss",JpsiMass,CBm,CBs);

	RooRealVar fracgauss("fracgauss","fracgauss",0.5,0,1);

	RooAddPdf gaussCB("gaussCB","gaussCB",RooArgList(gauss,CBMass),fracgauss);


	RooRealVar NsigMassCB("NsigMassCB","NsigMassCB",1,0,200000);
	RooExtendPdf eCBMass("eCBMass","eCBMass",gaussCB,NsigMassCB);

	RooAddPdf sumJpsiMass2("sumJpsiMass2","sumJpsiMass2",RooArgList(eCBMass,eexpo));


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DefinePromptDecayLengthModel(){

// 		Gauss 1
		RooRealVar ctmeanSig1("ctmeanSig1","ctmeanSig1",0);
		RooRealVar ctsigmaSig1("ctsigmaSig1","ctsigmaSig1",0.2,0.02,1);
		RooRealVar ctcoefSig1("ctcoefSig1","ctcoefSig1",1,0,100000);
		RooGaussModel ctgaussSig1("ctgaussSig1","ctgaussSig1",Jpsict,ctmeanSig1,ctsigmaSig1);
		RooExtendPdf ctegaussSig1("ctegaussSig1","ctegaussSig1",ctgaussSig1,ctcoefSig1);
//		Gauss 2 (same mean as Gauss 1)
		RooRealVar ctsigmaSig2("ctsigmaSig2","ctsigmaSig2",0.2,0.02,0.75);
		RooRealVar ctcoefSig2("ctcoefSig2","ctcoefSig2",1,0,100000);
		RooGaussModel ctgaussSig2("ctgaussSig2","ctgaussSig2",Jpsict,ctmeanSig1,ctsigmaSig2);
		RooExtendPdf ctegaussSig2("ctegaussSig2","ctegaussSig2",ctgaussSig2,ctcoefSig2);
//		Gauss 3
		RooRealVar ctmeanSig3("ctmeanSig3","ctmeanSig3",0,-0.025,0.025);
		RooRealVar ctsigmaSig3("ctsigmaSig3","ctsigmaSig3",0.5,0.05,1.2);
		RooRealVar ctcoefSig3("ctcoefSig3","ctcoefSig3",1,0,100000);
		RooGaussModel ctgaussSig3("ctgaussSig3","ctgaussSig3",Jpsict,ctmeanSig1,ctsigmaSig3);
		RooExtendPdf ctegaussSig3("ctegaussSig3","ctegaussSig3",ctgaussSig3,ctcoefSig3);

		RooAddModel sumJpsictPrompt("sumJpsictPrompt","sumJpsictPrompt",RooArgList(ctgaussSig1,ctgaussSig2,ctgaussSig3),RooArgList(ctcoefSig1,ctcoefSig2,ctcoefSig3));
		RooRealVar ctcoefPrompt("ctcoefPrompt","ctcoefPrompt",1,0,100000);
		RooExtendPdf esumJpsictPrompt("esumJpsictPrompt","esumJpsictPrompt",sumJpsictPrompt,ctcoefPrompt);


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ReduceDataSetToBin(){

		char reduceMCTestBin[200];
		sprintf(reduceMCTestBin,"abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiPt < %f && JpsiPt > %f",RapMin,RapMax,pTMax,pTMin);
		char reduceMCBin[200];
		sprintf(reduceMCBin,"abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiPt < %f && JpsiPt > %f",RapMin,RapMax,pTMax,pTMin);
		char reduceDataBin[200];
		sprintf(reduceDataBin,"abs(JpsiRap) > %f && abs(JpsiRap) < %f && JpsiPt < %f && JpsiPt > %f",RapMin,RapMax,pTMax,pTMin);

		reddata = (RooDataSet*)MCModelDataSet->reduce(reduceMCBin);

		reddataMCtest = (RooDataSet*)MCTestDataSet->reduce(reduceMCTestBin);

		realdata = (RooDataSet*)realdata00->reduce(reduceDataBin);

		realhData = new RooDataHist("realhData","realhData",RooArgSet(JpsiMass,Jpsict),*realdata);

		dataMCtesthist = new RooDataHist("dataMCtesthist","dataMCtesthist",RooArgSet(JpsiMass,Jpsict),*reddataMCtest);

		redMCPROMPT= (RooDataSet*)MCPROMPT->reduce(reduceMCTestBin);
		redMCNONPROMPT= (RooDataSet*)MCNONPROMPT->reduce(reduceMCTestBin);
		redMCBACKGROUND= (RooDataSet*)MCBACKGROUND->reduce(reduceMCTestBin);


		redMCPROMPThist = new RooDataHist("redMCPROMPThist","redMCPROMPThist",RooArgSet(JpsiMass,Jpsict),*redMCPROMPT);
		redMCNONPROMPThist = new RooDataHist("redMCNONPROMPThist","redMCNONPROMPThist",RooArgSet(JpsiMass,Jpsict),*redMCNONPROMPT);
		redMCBACKGROUNDhist = new RooDataHist("redMCBACKGROUNDhist","redMCBACKGROUNDhist",RooArgSet(JpsiMass,Jpsict),*redMCBACKGROUND);


		hltIndex0PR = (RooDataSet*)reddata->reduce("MCType_idx < 0.5");
		hltIndex0NP = (RooDataSet*)reddata->reduce("MCType_idx > 0.5 && MCType_idx < 1.5");
		SignalDataSet = (RooDataSet*)reddata->reduce("MCType_idx < 1.5");
		hltIndex0BK = (RooDataSet*)reddata->reduce("MCType_idx > 1.5 && MCType_idx < 2.5");

		hltIndex0PRMCtest = (RooDataSet*)reddataMCtest->reduce("MCType_idx < 0.5");
		hltIndex0NPMCtest = (RooDataSet*)reddataMCtest->reduce("MCType_idx > 0.5 && MCType_idx < 1.5");
		hltIndex0BKMCtest = (RooDataSet*)reddataMCtest->reduce("MCType_idx > 1.5 && MCType_idx < 2.5");
		hltIndex0PRhistMCtest = new RooDataHist("hltIndex0PRhistMCtest","",RooArgSet(JpsiMass,Jpsict),*hltIndex0PRMCtest);
		hltIndex0NPhistMCtest = new RooDataHist("hltIndex0NPhistMCtest","",RooArgSet(JpsiMass,Jpsict),*hltIndex0NPMCtest);
		hltIndex0BKhistMCtest = new RooDataHist("hltIndex0BKhistMCtest","",RooArgSet(JpsiMass,Jpsict),*hltIndex0BKMCtest);

		hltIndex0NPhist = new RooDataHist("hltIndex0NPhist","hltIndex0NPhist",RooArgSet(Jpsict),*hltIndex0NP,1);
		hltIndex0PRhist = new RooDataHist("hltIndex0PRhist","hltIndex0PRhist",RooArgSet(Jpsict),*hltIndex0PR,1);
		hltIndex0BKhist = new RooDataHist("hltIndex0BKhist","hltIndex0BKhist",RooArgSet(Jpsict),*hltIndex0BK,1);

		hData = new RooDataHist("hData","hData",RooArgSet(Jpsict,JpsiMass),*reddata);
		hDataMCtest = new RooDataHist("hDataMCtest","hDataMCtest",RooArgSet(Jpsict,JpsiMass),*reddataMCtest);


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitMass(){


	if(RealDataIndex == 0){//fit to MC
		fitresMass = sumJpsiMass2.fitTo(*hDataMCtest,Save(1),Minos(0),SumW2Error(kFALSE));

	}

	else if(RealDataIndex == 1){
		fitresMass = sumJpsiMass2.fitTo(*realhData,Save(1),Minos(0),SumW2Error(kFALSE));
	}

		NsigTot   = NsigMassCB.getVal();
		errNsigTot = NsigMassCB.getError();
		NbkgTot   = NbkgMass.getVal();
		errNbkgTot = NbkgMass.getError();

		JpsiMassframeMCTruthtestCB = JpsiMass.frame() ;
		SignalDataSet->plotOn(JpsiMassframeMCTruthtestCB,DataError(RooAbsData::SumW2));
		eCBMass.plotOn(JpsiMassframeMCTruthtestCB,Normalization(1.0));
		eCBMass.plotOn(JpsiMassframeMCTruthtestCB, RooFit::Components(RooArgList(gauss)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		eCBMass.plotOn(JpsiMassframeMCTruthtestCB, RooFit::Components(RooArgList(CBMass)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));


		errfracgauss = fracgauss.getError();
		fracgaussdouble = fracgauss.getVal();

		resol = CBs.getVal();
		errresol = CBs.getError();
//		resol = sqrt(fracgaussdouble*sigmaSig1double*sigmaSig1double + (1-fracgaussdouble)*sigmaSig2double*sigmaSig2double);
//		errresol = (0.5/resol)*sqrt(pow(sigmaSig1double*fracgaussdouble*errsigmaSig1,2) + pow(sigmaSig2double*(1-fracgaussdouble)*errsigmaSig2,2) + pow(0.5*(sigmaSig1double*sigmaSig1double - sigmaSig2double*sigmaSig2double)*errfracgauss,2));
//		meanSig1double = meanSig1.getVal();
		CBmean = CBm.getVal();
		errCBmean = CBm.getError();
		JpsiMass3SigmaMin = CBmean-3.*resol;
		JpsiMass3SigmaMax = CBmean+3.*resol;
		JpsiMass6SigmaMin = CBmean-6.*resol;
		JpsiMass6SigmaMax = CBmean+6.*resol;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CalcSigOvBkg(){

		JpsiMass.setRange("3sigmaRange", JpsiMass3SigmaMin, JpsiMass3SigmaMax);
		RooAbsReal* integralSignal =  eCBMass.createIntegral(RooArgSet(JpsiMass),NormSet(JpsiMass),Range("3sigmaRange"));
		integralNsigRange = integralSignal->getVal();
		RooAbsReal* integralbkg =  eexpo.createIntegral(RooArgSet(JpsiMass),NormSet(RooArgSet(JpsiMass)),Range("3sigmaRange"));
		integralNbkgRange = integralbkg->getVal();

		NsigRange = NsigTot*integralNsigRange;
		errNsigTotscaled = errNsigTot*NsigRange/NsigTot;

		NbkgRange = NbkgTot*integralNbkgRange;
		errNbkgTotscaled = errNbkgTot*NbkgRange/NbkgTot;
		sigOvbkg;
		errsigOvbkg;

		if(NbkgRange > 0){
			sigOvbkg = NsigRange / NbkgRange;
			errsigOvbkg = sigOvbkg * sqrt(pow(errNsigTot/NsigRange,2) + pow(errNbkgTot/NbkgRange,2));
			printf("sig/bkg = %f", sigOvbkg);
			printf(" +- %f\n\n", errsigOvbkg);
		  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PreFitPrompt(){

	fitresctPrompt = esumJpsictPrompt.fitTo(*redMCPROMPThist,Save(1),Minos(0),SumW2Error(kFALSE));

	JpsictframePromptAnal = Jpsict.frame(-1,2.5) ;
	redMCPROMPThist->plotOn(JpsictframePromptAnal,DataError(RooAbsData::SumW2));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal,Normalization(1.0));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal, RooFit::Components(RooArgList(ctgaussSig1)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal, RooFit::Components(RooArgList(ctgaussSig2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal, RooFit::Components(RooArgList(ctgaussSig3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

	JpsictframePromptAnalChi2 = Jpsict.frame(-1,2.5,1000) ;
	redMCPROMPThist->plotOn(JpsictframePromptAnalChi2,DataError(RooAbsData::SumW2));
	esumJpsictPrompt.plotOn(JpsictframePromptAnalChi2,Normalization(1.0));


	ctmeanSig1double = ctmeanSig1.getVal();
	ctsigmaSig1double = ctsigmaSig1.getVal();
	ctsigmaSig2double = ctsigmaSig2.getVal();
	ctsigmaSig3double = ctsigmaSig3.getVal();
	ctcoefSig1double = ctcoefSig1.getVal();
	ctcoefSig2double = ctcoefSig2.getVal();
	ctcoefSig3double = ctcoefSig3.getVal();

	ctcoefSigdoubleall=ctcoefSig1double+ctcoefSig2double+ctcoefSig3double;
	ctcoefSig1doublenorm=ctcoefSig1double/ctcoefSigdoubleall;
	ctcoefSig2doublenorm=ctcoefSig2double/ctcoefSigdoubleall;
	ctcoefSig3doublenorm=ctcoefSig3double/ctcoefSigdoubleall;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DrawAndSave(){



		char titlebackgroundlog[200];
		sprintf(titlebackgroundlog,"Jpsi c-tau background log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlebackgroundanallog[200];
		sprintf(titlebackgroundanallog,"Jpsi c-tau background anal log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titlesideband[200];
		sprintf(titlesideband,"Jpsi c-tau sideband, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlesidebandlog[200];
		sprintf(titlesidebandlog,"Jpsi c-tau sideband log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titleMass[200];
		sprintf(titleMass,"Jpsi mass, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlectlog[200];
		sprintf(titlectlog,"Jpsi c-tau global log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlect[200];
		sprintf(titlect,"Jpsi c-tau global, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);


		char titlecosth_CS[200];
		sprintf(titlecosth_CS,"costh_CS, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlecosth_HX[200];
		sprintf(titlecosth_HX,"costh_HX, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titlephi_CS[200];
		sprintf(titlephi_CS,"phi_CS, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlephi_HX[200];
		sprintf(titlephi_HX,"phi_HX, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titlecosth_phi_CS[200];
		sprintf(titlecosth_phi_CS,"costh_phi_CS, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlecosth_phi_HX[200];
		sprintf(titlecosth_phi_HX,"costh_phi_HX, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titlecosth_phi_CSPrompt[200];
		sprintf(titlecosth_phi_CSPrompt,"costh_phi_CSPrompt, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlecosth_phi_HXPrompt[200];
		sprintf(titlecosth_phi_HXPrompt,"costh_phi_HXPrompt, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titlecosth_phi_CSNonPrompt[200];
		sprintf(titlecosth_phi_CSNonPrompt,"costh_phi_CSNonPrompt, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlecosth_phi_HXNonPrompt[200];
		sprintf(titlecosth_phi_HXNonPrompt,"costh_phi_HXNonPrompt, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titlecosth_phi_CSBackground[200];
		sprintf(titlecosth_phi_CSBackground,"costh_phi_CSBackground, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
		char titlecosth_phi_HXBackground[200];
		sprintf(titlecosth_phi_HXBackground,"costh_phi_HXBackground, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char titleprompt[200];
		char titlepromptanal[200];
		char titlepromptlog[200];
		char titlepromptanallog[200];
		char titlenonprompt[200];
		char titlenonpromptanal[200];
		char titlenonpromptlog[200];
		char titlenonpromptanallog[200];
		char titlebackground[200];
		char titlebackgroundanal[200];
		char titlemasstestCB[200];

		TCanvas *c4;
		TCanvas *c5;
		TCanvas *c6;
		TCanvas *c7;
		TCanvas *c8;
		TCanvas *c9;
		TCanvas *c10;
		TCanvas *c11;
		TCanvas *c12;
		TCanvas *c13;
		TCanvas *c14;
		TCanvas *c15;
		TCanvas *c101;
		TCanvas *c102;

		TCanvas *SBRcomparisonCanvas;
		TCanvas *SBRcomparisonlogCanvas;
		TCanvas *SBLcomparisonCanvas;
		TCanvas *SBLcomparisonlogCanvas;

		if(RealDataIndex == 0){

			sprintf(titleprompt,"Jpsi c-tau prompt, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlepromptanal,"Jpsi c-tau prompt anal, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlepromptlog,"Jpsi c-tau prompt log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlepromptanallog,"Jpsi c-tau prompt anal log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlenonprompt,"Jpsi c-tau non prompt, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlenonpromptanal,"Jpsi c-tau non prompt anal, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlenonpromptlog,"Jpsi c-tau non prompt log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlenonpromptanallog,"Jpsi c-tau non prompt anal log, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlebackground,"Jpsi c-tau background, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlebackgroundanal,"Jpsi c-tau background anal, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);
			sprintf(titlemasstestCB,"Mass Signal Comparison, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		if (AnalyticIndex == 1){

	    c4 = new TCanvas(titleprompt,titleprompt,1200, 600);
	    c4->cd(1) ; JpsictframePrompt->Draw();

	    c5 = new TCanvas(titlepromptlog,titlepromptlog,1200, 600);
	    c5->SetLogy(1);
	    JpsictframePrompt->SetMinimum(1);
	    c5->cd(1) ; JpsictframePrompt->Draw();

	 	c6 = new TCanvas(titlepromptanal,titlepromptanal,1200, 600);
	 	c6->cd(1) ; JpsictframePromptAnal->Draw();
	 	c6->cd(1) ; PromptLegendAnal->Draw();

	 	c7 = new TCanvas(titlepromptanallog,titlepromptanallog,1200, 600);
	 	c7->SetLogy(1);
	 	JpsictframePromptAnal->SetMinimum(1);
	 	c7->cd(1) ; JpsictframePromptAnal->Draw();
	 	c7->cd(1) ; PromptLegendAnal->Draw();

	    c8 = new TCanvas(titlenonprompt,titlenonprompt,1200, 600);
	    c8->cd(1) ; JpsictframeNonPrompt->Draw();

	    c9 = new TCanvas(titlenonpromptlog,titlenonpromptlog,1200, 600);
	    c9->SetLogy(1);
	    JpsictframeNonPrompt->SetMinimum(1);
	    c9->cd(1) ; JpsictframeNonPrompt->Draw();

	    c10 = new TCanvas(titlenonpromptanal,titlenonpromptanal,1200, 600);
	    c10->cd(1) ; JpsictframeNonPromptAnal->Draw();
	    c10->cd(1) ; NonPromptLegendAnal->Draw();

	    c11 = new TCanvas(titlenonpromptanallog,titlenonpromptanallog,1200, 600);
	    c11->SetLogy(1);
	    JpsictframeNonPromptAnal->SetMinimum(1);
	    c11->cd(1) ; JpsictframeNonPromptAnal->Draw();
	    c11->cd(1) ; NonPromptLegendAnal->Draw();

	    c12 = new TCanvas(titlebackground,titlebackground,1200, 600);
	    c12->cd(1) ; JpsictframeBackground->Draw();

	    c13 = new TCanvas(titlebackgroundlog,titlebackgroundlog,1200, 600);
	    c13->SetLogy(1);
	    JpsictframeBackground->SetMinimum(1);
	    c13->cd(1) ; JpsictframeBackground->Draw();

	    c14 = new TCanvas(titlebackgroundanal,titlebackgroundanal,1200, 600);
	    c14->cd(1) ; JpsictframeBackgroundAnal->Draw();
	    c14->cd(1) ; BackgroundLegendAnal->Draw();

	    c15 = new TCanvas(titlebackgroundanallog,titlebackgroundanallog,1200, 600);
		c15->SetLogy(1);
		JpsictframeBackgroundAnal->SetMinimum(1);
		c15->cd(1) ; JpsictframeBackgroundAnal->Draw();
		c15->cd(1) ; BackgroundLegendAnal->Draw();
		}
		}

		else if(RealDataIndex == 1){

		if (SBComparisonIndex == 1){

		c101 = new TCanvas(titlesideband,titlesideband,1200, 600);
		c101->cd(1) ; JpsictframeSB->Draw();
		c101->cd(1) ; SBLegend->Draw();

		c102 = new TCanvas(titlesidebandlog,titlesidebandlog,1200, 600);
		c102->SetLogy(1);
		JpsictframeSB->SetMinimum(1);
		c102->cd(1) ; JpsictframeSB->Draw();
		c102->cd(1) ; SBLegend->Draw();

		SBRcomparisonCanvas = new TCanvas("right","right",1200, 600);
		SBRcomparisonCanvas->cd(1) ; JpsictframeBackgroundAnalRealDataSBR->Draw();
		SBRcomparisonCanvas->cd(1) ; BackgroundLegendAnalDataSBR->Draw();

		SBRcomparisonlogCanvas = new TCanvas("right log","right log",1200, 600);
		SBRcomparisonlogCanvas->SetLogy(1);
		JpsictframeBackgroundAnalRealDataSBR->SetMinimum(1);
		SBRcomparisonlogCanvas->cd(1) ; JpsictframeBackgroundAnalRealDataSBR->Draw();
		SBRcomparisonlogCanvas->cd(1) ; BackgroundLegendAnalDataSBR->Draw();

		SBLcomparisonCanvas = new TCanvas("left","left",1200, 600);
		SBLcomparisonCanvas->cd(1) ; JpsictframeBackgroundAnalRealDataSBL->Draw();
		SBLcomparisonCanvas->cd(1) ; BackgroundLegendAnalDataSBL->Draw();

		SBLcomparisonlogCanvas = new TCanvas("left log","left log",1200, 600);
		SBLcomparisonlogCanvas->SetLogy(1);
		JpsictframeBackgroundAnalRealDataSBL->SetMinimum(1);
		SBLcomparisonlogCanvas->cd(1) ; JpsictframeBackgroundAnalRealDataSBL->Draw();
		SBLcomparisonlogCanvas->cd(1) ; BackgroundLegendAnalDataSBL->Draw();

		}

		}
		TCanvas *c1;
		c1 = new TCanvas(titleMass,titleMass,1200, 600);
		c1->cd(1) ; JpsiMassframe->Draw() ;
		c1->cd(1) ; MassLegend->Draw();

		TCanvas *c10000;

		if (RealDataIndex == 0){

		c10000 = new TCanvas(titlemasstestCB,titlemasstestCB,1200, 600);
		c10000->cd(1)->SetLogy(1);
		JpsiMassframeMCTruthtestCB->SetMinimum(1);
		c10000->cd(1) ; JpsiMassframeMCTruthtestCB->Draw() ;
		}



	    TCanvas *c2 = new TCanvas(titlectlog,titlectlog,1200, 600);
	    c2->SetLogy(1);
	    Jpsictframe->SetMinimum(1);
	    c2->cd(1) ; Jpsictframe->Draw();
	    c2->cd(1) ; GlobalLegend->Draw();

	    TCanvas *c3 = new TCanvas(titlect,titlect,1200, 600);
	    c3->cd(1) ; Jpsictframe->Draw();
	    c3->cd(1) ; GlobalLegend->Draw();

		char title_CS_2D[200];
		sprintf(title_CS_2D,"Collins soper 2D Histo, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		char title_CS_1D[200];
		sprintf(title_CS_1D,"Collins soper 1D Histo, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

		TCanvas *thetaphi2_CS = new TCanvas(title_CS_1D,title_CS_1D,1200, 600);
		thetaphi2_CS->Divide(2);
		thetaphi2_CS->cd(1) ; phi_CSframe->Draw();
		thetaphi2_CS->cd(2) ; costh_CSframe->Draw();


			TCanvas *thetaphi_CS = new TCanvas(title_CS_2D,title_CS_2D,1200, 600);
		thetaphi_CS->Divide(2);
		thetaphi_CS->cd(1) ; datahistoangular_CS->Draw("colz");
		thetaphi_CS->cd(2) ; polfunchisto_CS->Draw("colz");

		TCanvas *thetaphi_CS_acc = new TCanvas("thetaphi_CS_acc","thetaphi_CS_acc",1200, 600);
		thetaphi_CS_acc->Divide(2);
		thetaphi_CS_acc->cd(1) ; acceptancehist_CS->Draw("lego");
		thetaphi_CS_acc->cd(2) ; polfunchist_CS->Draw("lego");

		TCanvas *thetaphi_HX_acc = new TCanvas("thetaphi_HX_acc","thetaphi_HX_acc",1200, 600);
		thetaphi_HX_acc->Divide(2);
		thetaphi_HX_acc->cd(1) ; acceptancehist_HX->Draw("lego");
		thetaphi_HX_acc->cd(2) ; polfunchist_HX->Draw("lego");

/*		TCanvas *thetaphi_Data = new TCanvas("thetaphi_Data","thetaphi_Data",1200, 600);
		thetaphi_Data->Divide(2);
		thetaphi_Data->cd(1) ; AngularRealData_hist_CS_TH1->Draw("lego");
		thetaphi_Data->cd(2) ; AngularRealData_hist_HX_TH1->Draw("lego");
*/

	    char title_HX_2D[200];
	    sprintf(title_HX_2D,"Helicity 2D Histo, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);

	    char title_HX_1D[200];
	    sprintf(title_HX_1D,"Helicity 1D Histo, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);


	    TCanvas *thetaphi2_HX = new TCanvas(title_HX_2D,title_HX_2D,1200, 600);
	    thetaphi2_HX->Divide(2);
	    thetaphi2_HX->cd(1) ; phi_HXframe->Draw();
	    thetaphi2_HX->cd(2) ; costh_HXframe->Draw();

	    TCanvas *thetaphi_HX = new TCanvas(title_HX_1D,title_HX_1D,1200, 600);
	    thetaphi_HX->Divide(2);
	    thetaphi_HX->cd(1) ; datahistoangular_HX->Draw("colz");
	    thetaphi_HX->cd(2) ; polfunchisto_HX->Draw("colz");




//////////////////////////////////////////////////// SAVE //////////////////////////////////////////////////////////////////////////////////


		cout << " " << endl;

		char oFile1[200];
		sprintf(oFile1,"MyFigures/Binning/MassPlot_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile2[200];
		sprintf(oFile2,"MyFigures/Binning/ctPlot_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile3[200];
		sprintf(oFile3,"MyFigures/Binning/ctPlot_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile4[200];
		sprintf(oFile4,"MyFigures/Binning/ctPlot_Prompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile5[200];
		sprintf(oFile5,"MyFigures/Binning/ctPlot_Prompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile6[200];
		sprintf(oFile6,"MyFigures/Binning/ctPlot_Prompt_Anal_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile7[200];
		sprintf(oFile7,"MyFigures/Binning/ctPlot_Prompt_Anal_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile8[200];
		sprintf(oFile8,"MyFigures/Binning/ctPlot_NonPrompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile9[200];
		sprintf(oFile9,"MyFigures/Binning/ctPlot_NonPrompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile10[200];
		sprintf(oFile10,"MyFigures/Binning/ctPlot_NonPrompt_Anal_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile11[200];
		sprintf(oFile11,"MyFigures/Binning/ctPlot_NonPrompt_Anal_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile12[200];
		sprintf(oFile12,"MyFigures/Binning/ctPlot_Background_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile13[200];
		sprintf(oFile13,"MyFigures/Binning/ctPlot_Background_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile14[200];
		sprintf(oFile14,"MyFigures/Binning/ctPlot_Background_Anal_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile15[200];
		sprintf(oFile15,"MyFigures/Binning/ctPlot_Background_Anal_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1DATA[200];
		sprintf(oFile1DATA,"MyFigures/Binning/DATA_MassPlot_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile2DATA[200];
		sprintf(oFile2DATA,"MyFigures/Binning/DATA_ctPlot_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile3DATA[200];
		sprintf(oFile3DATA,"MyFigures/Binning/DATA_ctPlot_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile101[200];
		sprintf(oFile101,"MyFigures/Binning/DATA_ctPlot_Sideband_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile102[200];
		sprintf(oFile102,"MyFigures/Binning/DATA_ctPlot_Sideband_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f_log.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1001[200];
		sprintf(oFile1001,"MyFigures/Binning/Angular_costh_CS_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1002[200];
		sprintf(oFile1002,"MyFigures/Binning/Angular_phi_CS_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1003a[200];
		sprintf(oFile1003a,"MyFigures/Binning/Angular_costh_phi_CS_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1003b[200];
		sprintf(oFile1003b,"MyFigures/Binning/Angular_costh_phi_CSPrompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1003c[200];
		sprintf(oFile1003c,"MyFigures/Binning/Angular_costh_phi_CSNonPrompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1003d[200];
		sprintf(oFile1003d,"MyFigures/Binning/Angular_costh_phi_CSBackground_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1004[200];
		sprintf(oFile1004,"MyFigures/Binning/Angular_costh_HX_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1005[200];
		sprintf(oFile1005,"MyFigures/Binning/Angular_phi_HX_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1006a[200];
		sprintf(oFile1006a,"MyFigures/Binning/Angular_costh_phi_HX_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1006b[200];
		sprintf(oFile1006b,"MyFigures/Binning/Angular_costh_phi_HXPrompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1006c[200];
		sprintf(oFile1006c,"MyFigures/Binning/Angular_costh_phi_HXNonPrompt_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile1006d[200];
		sprintf(oFile1006d,"MyFigures/Binning/Angular_costh_phi_HXBackground_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile10000[200];
		sprintf(oFile10000,"MyFigures/Binning/MassSignalTest_pT_%1.2f-%1.2f,Rap_%1.2f-%1.2f.png", pTMin, pTMax, RapMin, RapMax);


		if (SavePlotIndex == 1){

	    if(RealDataIndex == 0){
	    	c1->SaveAs(oFile1);
	    	c1->Close();
	    	c2->SaveAs(oFile2);
	    	c2->Close();
	    	c3->SaveAs(oFile3);
	    	c3->Close();
	    	c10000->SaveAs(oFile10000);
	    	c10000->Close();

	        if (AnalyticIndex == 1){

	    	c4->SaveAs(oFile4);
	    	c4->Close();
	    	c5->SaveAs(oFile5);
	    	c5->Close();
	    	c6->SaveAs(oFile6);
	    	c6->Close();
	    	c7->SaveAs(oFile7);
	    	c7->Close();
	    	c8->SaveAs(oFile8);
	    	c8->Close();
	    	c9->SaveAs(oFile9);
	    	c9->Close();
	    	c10->SaveAs(oFile10);
	    	c10->Close();
	    	c11->SaveAs(oFile11);
	    	c11->Close();
	    	c12->SaveAs(oFile12);
	    	c12->Close();
	    	c13->SaveAs(oFile13);
	    	c13->Close();
	    	c14->SaveAs(oFile14);
	    	c14->Close();
	    	c15->SaveAs(oFile15);
	    	c15->Close();
	        }
	    }
	    else if(RealDataIndex == 1){
	    	c1->SaveAs(oFile1DATA);
	    	c1->Close();
	    	c2->SaveAs(oFile2DATA);
	    	c2->Close();
	    	c3->SaveAs(oFile3DATA);
	    	c3->Close();
	   if (SBComparisonIndex == 1){
	    	c101->SaveAs(oFile101);
	    	c101->Close();
	    	c102->SaveAs(oFile102);
	    	c102->Close();
	    	SBRcomparisonCanvas->SaveAs("MyFigures/Binning/ProjectNewNew/SBrightanal.png");
			SBRcomparisonCanvas->Close();
			SBRcomparisonlogCanvas->SaveAs("MyFigures/Binning/ProjectNewNew/SBrightanallog.png");
			SBRcomparisonlogCanvas->Close();
			SBLcomparisonCanvas->SaveAs("MyFigures/Binning/ProjectNewNew/SBleftanal.png");
			SBLcomparisonCanvas->Close();
			SBLcomparisonlogCanvas->SaveAs("MyFigures/Binning/ProjectNewNew/SBleftanallog.png");
			SBLcomparisonlogCanvas->Close();
	   }
	    }

		char oFile_CS_1D[200];
		sprintf(oFile_CS_1D,"MyFigures/Binning/Pol_CS_1D_pT_%1.0f-%1.0f,Rap_%1.0f-%1.0f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile_CS_2D[200];
		sprintf(oFile_CS_2D,"MyFigures/Binning/Pol_CS_2D_pT_%1.0f-%1.0f,Rap_%1.0f-%1.0f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile_HX_1D[200];
		sprintf(oFile_HX_1D,"MyFigures/Binning/Pol_HX_1D_pT_%1.0f-%1.0f,Rap_%1.0f-%1.0f.png", pTMin, pTMax, RapMin, RapMax);
		char oFile_HX_2D[200];
		sprintf(oFile_HX_2D,"MyFigures/Binning/Pol_HX_2D_pT_%1.0f-%1.0f,Rap_%1.0f-%1.0f.png", pTMin, pTMax, RapMin, RapMax);

			thetaphi_CS->SaveAs(oFile_CS_2D);
		thetaphi2_CS->SaveAs(oFile_CS_1D);
			thetaphi_HX->SaveAs(oFile_HX_2D);
		thetaphi2_HX->SaveAs(oFile_HX_1D);

		}

		cout<<"DrawAndSave Complete"<<endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Legend(){


		legendMassDashed = realdata000->createHistogram("legendMassDashed",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassDashed->SetLineColor(kBlack) ;
		legendMassDashed->SetLineStyle(kDashed) ;
		legendMassDashed->SetLineWidth(2) ;

		legendMassDashedred = realdata000->createHistogram("legendMassDashedred",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassDashedred->SetLineColor(kRed) ;
		legendMassDashedred->SetLineStyle(kDashed) ;
		legendMassDashedred->SetLineWidth(2) ;

		legendMassDashedgreen = realdata000->createHistogram("legendMassDashedgreen",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassDashedgreen->SetLineColor(kGreen) ;
		legendMassDashedgreen->SetLineStyle(kDashed) ;
		legendMassDashedgreen->SetLineWidth(2) ;

		legendMass = realdata000->createHistogram("legendMass",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMass->SetLineColor(kBlue) ;
		legendMass->SetLineWidth(2) ;

		legendMassdata = realdata000->createHistogram("legendMassdata",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdata->SetLineStyle(kSolid) ;
		legendMassdata->SetLineWidth(1) ;
		legendMassdata->SetMarkerStyle(8) ;

		legendMassdatared = realdata000->createHistogram("legendMassdatared",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdatared->SetMarkerColor(kRed) ;
		legendMassdatared->SetLineStyle(kSolid) ;
		legendMassdatared->SetLineWidth(1) ;
		legendMassdatared->SetMarkerStyle(8) ;

		legendMassdatablue = realdata000->createHistogram("legendMassdatablue",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdatablue->SetMarkerColor(kBlue) ;
		legendMassdatablue->SetLineStyle(kSolid) ;
		legendMassdatablue->SetLineWidth(1) ;
		legendMassdatablue->SetMarkerStyle(8) ;

		legendMassdataStyle = realdata000->createHistogram("legendMassdataStyle",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdataStyle->SetLineStyle(kSolid) ;
		legendMassdataStyle->SetLineWidth(1) ;
		legendMassdataStyle->SetMarkerStyle(21) ;

		legendMassdataredStyle = realdata000->createHistogram("legendMassdataredStyle",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdataredStyle->SetMarkerColor(kRed) ;
		legendMassdataredStyle->SetLineStyle(kSolid) ;
		legendMassdataredStyle->SetLineWidth(1) ;
		legendMassdataredStyle->SetMarkerStyle(21) ;

		legendMassdatablueStyle = realdata000->createHistogram("legendMassdatablueStyle",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdatablueStyle->SetMarkerColor(kBlue) ;
		legendMassdatablueStyle->SetLineStyle(kSolid) ;
		legendMassdatablueStyle->SetLineWidth(1) ;
		legendMassdatablueStyle->SetMarkerStyle(21) ;

		legendMassdatagreenStyle = realdata000->createHistogram("legendMassdatagreenStyle",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdatagreenStyle->SetMarkerColor(kGreen) ;
		legendMassdatagreenStyle->SetLineStyle(kSolid) ;
		legendMassdatagreenStyle->SetLineWidth(1) ;
		legendMassdatagreenStyle->SetMarkerStyle(21) ;

		legendMassdatayellowStyle = realdata000->createHistogram("legendMassdatayellowStyle",JpsiMass,Binning(100),YVar(Jpsict,Binning(100))) ;
		legendMassdatablueStyle->SetMarkerColor(kYellow) ;
		legendMassdatablueStyle->SetLineStyle(kSolid) ;
		legendMassdatablueStyle->SetLineWidth(1) ;
		legendMassdatablueStyle->SetMarkerStyle(21) ;


	  MassLegend=new TLegend(0.65,0.65,0.92,0.8);
	  MassLegend->SetTextFont(72);
      MassLegend->SetTextSize(0.04);
  	if(RealDataIndex == 0){
      MassLegend->AddEntry(legendMassdata,"MC","ple");
  	}
	else if(RealDataIndex == 1){
	  MassLegend->AddEntry(legendMassdata,"Real data","ple");
	}
	  MassLegend->AddEntry(legendMass,"Total mass fit","l");
	  MassLegend->AddEntry(legendMassDashedred,"Background","l");
	  MassLegend->AddEntry(legendMassDashed,"Signal components","l");



	  PromptLegend=new TLegend(0.65,0.65,0.92,0.8);
	  PromptLegend->SetTextFont(72);
	  PromptLegend->SetTextSize(0.04);
	  PromptLegend->AddEntry(legendMassdata,"MC prompt","ple");
	  PromptLegend->AddEntry(legendMass,"Fit","l");

	  PromptLegendAnal=new TLegend(0.65,0.65,0.92,0.8);
	  PromptLegendAnal->SetTextFont(72);
	  PromptLegendAnal->SetTextSize(0.04);
	  PromptLegendAnal->AddEntry(legendMassdata,"MC prompt","ple");
	  PromptLegendAnal->AddEntry(legendMass,"Analytical fit","l");
	  PromptLegendAnal->AddEntry(legendMassDashed,"Individual components","l");

	  NonPromptLegend=new TLegend(0.65,0.65,0.92,0.8);
	  NonPromptLegend->SetTextFont(72);
	  NonPromptLegend->SetTextSize(0.04);
	  NonPromptLegend->AddEntry(legendMassdata,"MC non prompt","ple");
	  NonPromptLegend->AddEntry(legendMass,"Fit","l");

	  NonPromptLegendAnal=new TLegend(0.65,0.65,0.92,0.8);
	  NonPromptLegendAnal->SetTextFont(72);
	  NonPromptLegendAnal->SetTextSize(0.04);
	  NonPromptLegendAnal->AddEntry(legendMassdata,"MC non prompt","ple");
	  NonPromptLegendAnal->AddEntry(legendMass,"Analytical fit","l");
	  NonPromptLegendAnal->AddEntry(legendMassDashed,"Individual components","l");


	  BackgroundLegend=new TLegend(0.65,0.65,0.92,0.8);
	  BackgroundLegend->SetTextFont(72);
	  BackgroundLegend->SetTextSize(0.04);
	  BackgroundLegend->AddEntry(legendMassdata,"MC background","ple");
	  BackgroundLegend->AddEntry(legendMass,"Fit","l");

	  BackgroundLegendAnal=new TLegend(0.65,0.65,0.92,0.8);
	  BackgroundLegendAnal->SetTextFont(72);
	  BackgroundLegendAnal->SetTextSize(0.04);
	  BackgroundLegendAnal->AddEntry(legendMassdata,"MC background","ple");
	  BackgroundLegendAnal->AddEntry(legendMass,"Analytical fit","l");
	  BackgroundLegendAnal->AddEntry(legendMassDashed,"Individual components","l");

	  BackgroundLegendAnalDataSBR=new TLegend(0.65,0.65,0.92,0.8);
	  BackgroundLegendAnalDataSBR->SetTextFont(72);
	  BackgroundLegendAnalDataSBR->SetTextSize(0.04);
	  BackgroundLegendAnalDataSBR->AddEntry(legendMassdata,"Real data (SB right)","ple");
	  BackgroundLegendAnalDataSBR->AddEntry(legendMass,"Analytical fit","l");
	  BackgroundLegendAnalDataSBR->AddEntry(legendMassDashed,"Individual components","l");

	  BackgroundLegendAnalDataSBL=new TLegend(0.65,0.65,0.92,0.8);
	  BackgroundLegendAnalDataSBL->SetTextFont(72);
	  BackgroundLegendAnalDataSBL->SetTextSize(0.04);
	  BackgroundLegendAnalDataSBL->AddEntry(legendMassdata,"Real data (SB left)","ple");
	  BackgroundLegendAnalDataSBL->AddEntry(legendMass,"Analytical fit","l");
	  BackgroundLegendAnalDataSBL->AddEntry(legendMassDashed,"Individual components","l");

	  GlobalLegend=new TLegend(0.65,0.58,0.92,0.8);
	  GlobalLegend->SetTextFont(72);
	  GlobalLegend->SetTextSize(0.04);


	  if(RealDataIndex == 0){
	  GlobalLegend->AddEntry(legendMassdata,"MC components","ple");
	  }
	  else if(RealDataIndex == 1){
	  GlobalLegend->AddEntry(legendMassdata,"Real data","ple");
	  }
	  GlobalLegend->AddEntry(legendMass,"Global fit","l");
	  GlobalLegend->AddEntry(legendMassDashedred,"Prompt fit","l");
	  GlobalLegend->AddEntry(legendMassDashed,"Background fit","l");
	  GlobalLegend->AddEntry(legendMassDashedgreen,"Non prompt fit","l");

	  SBLegend=new TLegend(0.65,0.65,0.92,0.8);
	  SBLegend->SetTextFont(72);
	  SBLegend->SetTextSize(0.04);
	  SBLegend->AddEntry(legendMassdata,"SB left (scaled)","ple");
	  SBLegend->AddEntry(legendMassdatared,"SB right","ple");

	  char titlelegendrap1[200];
	  if (rapForPTRange[0]<0.001){
	  sprintf(titlelegendrap1,"|y| < %1.2f", rapForPTRange[1]);
	  }
	  else if (rapForPTRange[0]>0.001){
	  sprintf(titlelegendrap1,"%1.2f < |y| < %1.2f", rapForPTRange[0], rapForPTRange[1]);
	  }
	  if (numberofRapbins == 2 |numberofRapbins == 3){
	  char titlelegendrap2[200];
	  sprintf(titlelegendrap2,"%1.2f < |y| < %1.2f", rapForPTRange[1],rapForPTRange[2]);
	  }
	  if (numberofRapbins == 3){
	  char titlelegendrap3[200];
	  sprintf(titlelegendrap3,"%1.2f < |y| < %1.2f", rapForPTRange[2],rapForPTRange[3]);
	  }
	  char titlelegendrap1MCTRUTH[200];
	  if (rapForPTRange[1]<0.001){
	  sprintf(titlelegendrap1MCTRUTH,"|y| < %1.2f - MC Truth", rapForPTRange[1]);
	  }
	  else if (rapForPTRange[1]>0.001){
	  sprintf(titlelegendrap1MCTRUTH,"%1.2f < |y| < %1.2f - MC Truth", rapForPTRange[0], rapForPTRange[1]);
	  }
	  if (numberofRapbins == 2){
	  char titlelegendrap2MCTRUTH[200];
	  sprintf(titlelegendrap2MCTRUTH,"%1.2f < |y| < %1.2f - MC Truth", rapForPTRange[1],rapForPTRange[2]);
	  }
	  if (numberofRapbins == 3){
	  char titlelegendrap3MCTRUTH[200];
	  sprintf(titlelegendrap3MCTRUTH,"%1.2f < |y| < %1.2f - MC Truth", rapForPTRange[2],rapForPTRange[3]);
	  }
	  BFractionLegend=new TLegend(0.15,0.6,0.45,0.8);
	  BFractionLegend->SetTextFont(72);
	  BFractionLegend->SetTextSize(0.03);
	  BFractionLegend->SetFillColor(kWhite);
	  BFractionLegend->SetLineColor(kWhite);
	  if (RealDataIndex == 1){
	BFractionLegend->AddEntry(legendMassdataStyle,"CDF: |y|<0.6, #sqrt{s}=1.96 TeV","p");
	  }
	if (numberofRapbins==1){
	  BFractionLegend->AddEntry(legendMassdatablue,titlelegendrap1,"p");
	if (RealDataIndex == 0){
		BFractionLegend->AddEntry(legendMassdatagreenStyle,titlelegendrap1MCTRUTH,"p");
	}
	}
	 else if (numberofRapbins==2){
	BFractionLegend->AddEntry(legendMassdatablue,titlelegendrap1,"p");
	BFractionLegend->AddEntry(legendMassdatared,titlelegendrap2,"p");
	if (RealDataIndex == 0){
			BFractionLegend->AddEntry(legendMassdatagreenStyle,titlelegendrap1MCTRUTH,"p");
			BFractionLegend->AddEntry(legendMassdatayellowStyle,titlelegendrap2MCTRUTH,"p");
		}
	 }
	 else if (numberofRapbins==3){
	BFractionLegend->AddEntry(legendMassdatablue,titlelegendrap1,"p");
	BFractionLegend->AddEntry(legendMassdatared,titlelegendrap2,"p");
	BFractionLegend->AddEntry(legendMassdata,titlelegendrap3,"p");
	if (RealDataIndex == 0){
			BFractionLegend->AddEntry(legendMassdatagreenStyle,titlelegendrap1MCTRUTH,"p");
			BFractionLegend->AddEntry(legendMassdatayellowStyle,titlelegendrap2MCTRUTH,"p");
			BFractionLegend->AddEntry(legendMassdataStyle,titlelegendrap3MCTRUTH,"p");
		}
	 }

	cout<<"Legend Complete"<<endl;


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DrawFinal(){


			const int kNbPoints = 26;
		  Double_t bFraction[kNbPoints] = {0.094, 0.092, 0.085, 0.100,
						   0.091, 0.101, 0.099, 0.109,
						   0.112, 0.113, 0.133, 0.116,
						   0.126, 0.131, 0.147, 0.141,
						   0.156, 0.169, 0.182, 0.208,
						   0.227, 0.250, 0.279, 0.337,
						   0.397, 0.464};
		  Double_t bFrac_Stat[kNbPoints] = {0.010, 0.006, 0.006, 0.005,
						   0.005, 0.005, 0.005, 0.005,
						   0.005, 0.005, 0.005, 0.005,
						   0.006, 0.006, 0.007, 0.005,
						   0.006, 0.007, 0.007, 0.006,
						   0.009, 0.011, 0.012, 0.019,
						   0.025, 0.045};
		  Double_t bFrac_Syst[kNbPoints] = {0.012, 0.010, 0.009, 0.011,
						   0.010, 0.009, 0.008, 0.007,
						   0.008, 0.007, 0.007, 0.007,
						   0.007, 0.007, 0.008, 0.006,
						   0.007, 0.007, 0.008, 0.009,
						   0.007, 0.008, 0.008, 0.009,
						   0.009, 0.014};
		  Double_t bFracErr[kNbPoints];
		  for (int newn=0;newn<kNbPoints;newn++){

			bFracErr[newn] = bFrac_Stat[newn]+bFrac_Syst[newn];
		   }
		  Double_t bFracErrpT[kNbPoints] = {0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07};
		  Double_t pTMid[kNbPoints] = {1.38, 1.63, 1.87, 2.13, 2.38, 2.62,
						  2.87, 3.12, 3.38, 3.62, 3.87, 4.12,
						  4.38, 4.62, 4.88, 5.24, 5.75, 6.24,
						  6.74, 7.45, 8.46, 9.46, 10.8, 12.8,
						  15.2, 18.3};

		double BFractionplotmax=0.6;


		if (RealDataIndex==0){
		BFractionplotmax=BFRACTION+ERRBFRACTION+0.05;
		}
		if (BFractionplotmax>1){
			BFractionplotmax=1;
}
		   TCanvas *BfractionCanvas = new TCanvas("fraction of J/#psi from B hadrons","fraction of J/#psi from B hadrons",1000,700);

			  TH1F *hr = BfractionCanvas->DrawFrame(0,0,pTMax,BFractionplotmax);
			  hr->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
			  hr->SetYTitle("fraction of J/#psi from B hadrons");
			  BfractionCanvas->GetFrame()->SetBorderSize(6);
			  BFractionLegend->Draw();

		   TGraphErrors *gr = new TGraphErrors(numberofpTbins,pTmean,yrap1,errpTmean,eyrap1);
		   gr->SetMarkerColor(kBlue);
		   gr->SetMarkerStyle(8);
		   gr->Draw("P");

		   if (RealDataIndex == 0){

			  TGraphErrors *grMCTRUTH = new TGraphErrors(numberofpTbins,pTmean,yrap1MCTRUTH,errpTmean,eyrap1MCTRUTH);
			  grMCTRUTH->SetMarkerColor(kGreen);
			  grMCTRUTH->SetMarkerStyle(21);
			  grMCTRUTH->Draw("P");
		   }

		   if (RealDataIndex == 1){
			  TGraphErrors *grCDF = new TGraphErrors(kNbPoints,pTMid,bFraction,bFracErrpT,bFracErr);
			  grCDF->SetMarkerColor(kBlack);
			  grCDF->SetMarkerStyle(21);
			  grCDF->Draw("P");
		   }

		   if (numberofRapbins==2){

		   TGraphErrors *gr2 = new TGraphErrors(numberofpTbins,pTmean,yrap2,errpTmean,eyrap2);
		   gr2->SetMarkerColor(kRed);
		   gr2->SetMarkerStyle(8);
		   gr2->Draw("P");

		   if (RealDataIndex == 0){

			  TGraphErrors *gr2MCTRUTH = new TGraphErrors(numberofpTbins,pTmean,yrap2MCTRUTH,errpTmean,eyrap2MCTRUTH);
			  gr2MCTRUTH->SetMarkerColor(kYellow);
			  gr2MCTRUTH->SetMarkerStyle(21);
			  gr2MCTRUTH->Draw("P");
		   }

		   }

		   else if (numberofRapbins==3){

		TGraphErrors *gr2 = new TGraphErrors(numberofpTbins,pTmean,yrap2,errpTmean,eyrap2);
		gr2->SetMarkerColor(kRed);
		gr2->SetMarkerStyle(8);
		gr2->Draw("P");

		TGraphErrors *gr3 = new TGraphErrors(numberofpTbins,pTmean,yrap3,errpTmean,eyrap3);
		gr3->SetMarkerColor(kBlack);
		gr3->SetMarkerStyle(8);
		gr3->Draw("P");

		if (RealDataIndex == 0){

		TGraphErrors *gr2MCTRUTH = new TGraphErrors(numberofpTbins,pTmean,yrap2MCTRUTH,errpTmean,eyrap2MCTRUTH);
		gr2MCTRUTH->SetMarkerColor(kRed);
		gr2MCTRUTH->SetMarkerStyle(21);
		gr2MCTRUTH->Draw("P");

		TGraphErrors *gr3MCTRUTH = new TGraphErrors(numberofpTbins,pTmean,yrap3MCTRUTH,errpTmean,eyrap3MCTRUTH);
		gr3MCTRUTH->SetMarkerColor(kBlack);
		gr3MCTRUTH->SetMarkerStyle(8);
		gr3MCTRUTH->Draw("P");

		}

		   }
		 if (RealDataIndex==0){
		   char oFileBFraction[200];
		   sprintf(oFileBFraction,"MyFigures/Binning/BFraction_MC_NpT%1.0f_NRap%1.0f_RapMin%1.2f_RapMax%1.2f.png", numberofpTbins, numberofRapbins, RapMin, RapMax);
		 }
		 if (RealDataIndex==1){
				   char oFileBFraction[200];
				   sprintf(oFileBFraction,"MyFigures/Binning/BFraction_DATA_NpT%1.0f_NRap%1.0f_RapMin%1.2f_RapMax%1.2f.png", numberofpTbins, numberofRapbins, RapMin, RapMax);
				 }
				   if (SavePlotIndex == 1){
			   BfractionCanvas->SaveAs(oFileBFraction);
					   BfractionCanvas->Close();
				   }



		TCanvas *PolThetaCanvas_CS = new TCanvas("PolThetaCanvas_CS","PolThetaCanvas_CS",1000,700);

		TH1F *PolThetaCanvas_CS_TH1 = PolThetaCanvas_CS->DrawFrame(0,-1.3,pTMax,1.3);
		PolThetaCanvas_CS_TH1->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
		PolThetaCanvas_CS_TH1->SetYTitle("#lambda_{#theta}_{CS}");
		PolThetaCanvas_CS->GetFrame()->SetBorderSize(6);

		TGraphErrors *PolThetaCanvas_CS_gr = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_1,errpTmean,eypol_CS_1);
		PolThetaCanvas_CS_gr->SetMarkerColor(kBlue);
		PolThetaCanvas_CS_gr->SetMarkerStyle(8);
		PolThetaCanvas_CS_gr->Draw("P");

		if (numberofRapbins==2){

		TGraphErrors *PolThetaCanvas_CS_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_1_rap2,errpTmean,eypol_CS_1_rap2);
		PolThetaCanvas_CS_gr_rap2->SetMarkerColor(kRed);
		PolThetaCanvas_CS_gr_rap2->SetMarkerStyle(8);
		PolThetaCanvas_CS_gr_rap2->Draw("P");

		}

		if (numberofRapbins==3){

		TGraphErrors *PolThetaCanvas_CS_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_1_rap2,errpTmean,eypol_CS_1_rap2);
		PolThetaCanvas_CS_gr_rap2->SetMarkerColor(kRed);
		PolThetaCanvas_CS_gr_rap2->SetMarkerStyle(8);
		PolThetaCanvas_CS_gr_rap2->Draw("P");

		TGraphErrors *PolThetaCanvas_CS_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_1_rap3,errpTmean,eypol_CS_1_rap3);
		PolThetaCanvas_CS_gr_rap3->SetMarkerColor(kBlack);
		PolThetaCanvas_CS_gr_rap3->SetMarkerStyle(8);
		PolThetaCanvas_CS_gr_rap3->Draw("P");

		}

		TCanvas *PolPhiCanvas_CS = new TCanvas("PolPhiCanvas_CS","PolPhiCanvas_CS",1000,700);

		TH1F *PolPhiCanvas_CS_TH1 = PolPhiCanvas_CS->DrawFrame(0,-1.3,pTMax,1.3);
		PolPhiCanvas_CS_TH1->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
		PolPhiCanvas_CS_TH1->SetYTitle("#lambda_{#Phi}_{CS}");
		PolPhiCanvas_CS->GetFrame()->SetBorderSize(6);

		TGraphErrors *PolPhiCanvas_CS_gr = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_2,errpTmean,eypol_CS_2);
		PolPhiCanvas_CS_gr->SetMarkerColor(kBlue);
		PolPhiCanvas_CS_gr->SetMarkerStyle(8);
		PolPhiCanvas_CS_gr->Draw("P");

		if (numberofRapbins==2){

		TGraphErrors *PolPhiCanvas_CS_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_2_rap2,errpTmean,eypol_CS_2_rap2);
		PolPhiCanvas_CS_gr_rap2->SetMarkerColor(kRed);
		PolPhiCanvas_CS_gr_rap2->SetMarkerStyle(8);
		PolPhiCanvas_CS_gr_rap2->Draw("P");

		}

		if (numberofRapbins==3){

		TGraphErrors *PolPhiCanvas_CS_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_2_rap2,errpTmean,eypol_CS_2_rap2);
		PolPhiCanvas_CS_gr_rap2->SetMarkerColor(kRed);
		PolPhiCanvas_CS_gr_rap2->SetMarkerStyle(8);
		PolPhiCanvas_CS_gr_rap2->Draw("P");

		TGraphErrors *PolPhiCanvas_CS_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_2_rap3,errpTmean,eypol_CS_2_rap3);
		PolPhiCanvas_CS_gr_rap3->SetMarkerColor(kBlack);
		PolPhiCanvas_CS_gr_rap3->SetMarkerStyle(8);
		PolPhiCanvas_CS_gr_rap3->Draw("P");

		}



		/*
		TCanvas *PolThetaPhiCanvas_CS = new TCanvas("PolThetaPhiCanvas_CS","PolThetaPhiCanvas_CS",1000,700);

		TH1F *PolThetaPhiCanvas_CS_TH1 = PolThetaPhiCanvas_CS->DrawFrame(0,-1.3,pTMax,1.3);
		PolThetaPhiCanvas_CS_TH1->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
		PolThetaPhiCanvas_CS_TH1->SetYTitle("#lambda_{#Theta#Phi}_{CS}");
		PolThetaPhiCanvas_CS->GetFrame()->SetBorderSize(6);

		TGraphErrors *PolThetaPhiCanvas_CS_gr = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_3,errpTmean,eypol_CS_3);
		PolThetaPhiCanvas_CS_gr->SetMarkerColor(kBlue);
		PolThetaPhiCanvas_CS_gr->SetMarkerStyle(8);
		PolThetaPhiCanvas_CS_gr->Draw("P");


		if (numberofRapbins==2){

		TGraphErrors *PolThetaPhiCanvas_CS_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_3_rap2,errpTmean,eypol_CS_3_rap2);
		PolThetaPhiCanvas_CS_gr_rap2->SetMarkerColor(kRed);
		PolThetaPhiCanvas_CS_gr_rap2->SetMarkerStyle(8);
		PolThetaPhiCanvas_CS_gr_rap2->Draw("P");

		}

		if (numberofRapbins==3){

		TGraphErrors *PolThetaPhiCanvas_CS_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_3_rap2,errpTmean,eypol_CS_3_rap2);
		PolThetaPhiCanvas_CS_gr_rap2->SetMarkerColor(kRed);
		PolThetaPhiCanvas_CS_gr_rap2->SetMarkerStyle(8);
		PolThetaPhiCanvas_CS_gr_rap2->Draw("P");

		TGraphErrors *PolThetaPhiCanvas_CS_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_CS_3_rap3,errpTmean,eypol_CS_3_rap3);
		PolThetaPhiCanvas_CS_gr_rap3->SetMarkerColor(kBlack);
		PolThetaPhiCanvas_CS_gr_rap3->SetMarkerStyle(8);
		PolThetaPhiCanvas_CS_gr_rap3->Draw("P");

		}
		*/

		TCanvas *PolThetaCanvas_HX = new TCanvas("PolThetaCanvas_HX","PolThetaCanvas_HX",1000,700);

		TH1F *PolThetaCanvas_HX_TH1 = PolThetaCanvas_HX->DrawFrame(0,-1.3,pTMax,1.3);
		PolThetaCanvas_HX_TH1->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
		PolThetaCanvas_HX_TH1->SetYTitle("#lambda_{#theta}_{HX}");
		PolThetaCanvas_HX->GetFrame()->SetBorderSize(6);

		TGraphErrors *PolThetaCanvas_HX_gr = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1,errpTmean,eypol_HX_1);
		PolThetaCanvas_HX_gr->SetMarkerColor(kBlue);
		PolThetaCanvas_HX_gr->SetMarkerStyle(8);
		PolThetaCanvas_HX_gr->Draw("P");


		if (numberofRapbins==2){

		TGraphErrors *PolThetaCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1_rap2,errpTmean,eypol_HX_1_rap2);
		PolThetaCanvas_HX_gr_rap2->SetMarkerColor(kRed);
		PolThetaCanvas_HX_gr_rap2->SetMarkerStyle(8);
		PolThetaCanvas_HX_gr_rap2->Draw("P");

		}

		if (numberofRapbins==3){

		TGraphErrors *PolThetaCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1_rap2,errpTmean,eypol_HX_1_rap2);
		PolThetaCanvas_HX_gr_rap2->SetMarkerColor(kRed);
		PolThetaCanvas_HX_gr_rap2->SetMarkerStyle(8);
		PolThetaCanvas_HX_gr_rap2->Draw("P");

		TGraphErrors *PolThetaCanvas_HX_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1_rap3,errpTmean,eypol_HX_1_rap3);
		PolThetaCanvas_HX_gr_rap3->SetMarkerColor(kBlack);
		PolThetaCanvas_HX_gr_rap3->SetMarkerStyle(8);
		PolThetaCanvas_HX_gr_rap3->Draw("P");

		}




		TCanvas *PolPhiCanvas_HX = new TCanvas("PolPhiCanvas_HX","PolPhiCanvas_HX",1000,700);

		TH1F *PolPhiCanvas_HX_TH1 = PolPhiCanvas_HX->DrawFrame(0,-1.3,pTMax,1.3);
		PolPhiCanvas_HX_TH1->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
		PolPhiCanvas_HX_TH1->SetYTitle("#lambda_{#Phi}_{HX}");
		PolPhiCanvas_HX->GetFrame()->SetBorderSize(6);

		TGraphErrors *PolPhiCanvas_HX_gr = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_2,errpTmean,eypol_HX_2);
		PolPhiCanvas_HX_gr->SetMarkerColor(kBlue);
		PolPhiCanvas_HX_gr->SetMarkerStyle(8);
		PolPhiCanvas_HX_gr->Draw("P");

		if (numberofRapbins==2){

		TGraphErrors *PolPhiCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_2_rap2,errpTmean,eypol_HX_2_rap2);
		PolPhiCanvas_HX_gr_rap2->SetMarkerColor(kRed);
		PolPhiCanvas_HX_gr_rap2->SetMarkerStyle(8);
		PolPhiCanvas_HX_gr_rap2->Draw("P");

		}

		if (numberofRapbins==3){

		TGraphErrors *PolPhiCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_2_rap2,errpTmean,eypol_HX_2_rap2);
		PolPhiCanvas_HX_gr_rap2->SetMarkerColor(kRed);
		PolPhiCanvas_HX_gr_rap2->SetMarkerStyle(8);
		PolPhiCanvas_HX_gr_rap2->Draw("P");

		TGraphErrors *PolPhiCanvas_HX_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_2_rap3,errpTmean,eypol_HX_2_rap3);
		PolPhiCanvas_HX_gr_rap3->SetMarkerColor(kBlack);
		PolPhiCanvas_HX_gr_rap3->SetMarkerStyle(8);
		PolPhiCanvas_HX_gr_rap3->Draw("P");

		}

		/*
		TCanvas *PolThetaPhiCanvas_HX = new TCanvas("PolThetaPhiCanvas_HX","PolThetaPhiCanvas_HX",1000,700);

		TH1F *PolThetaPhiCanvas_HX_TH1 = PolThetaPhiCanvas_HX->DrawFrame(0,-1.3,pTMax,1.3);
		PolThetaPhiCanvas_HX_TH1->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
		PolThetaPhiCanvas_HX_TH1->SetYTitle("#lambda_{#Theta#Phi}_{HX}");
		PolThetaPhiCanvas_HX->GetFrame()->SetBorderSize(6);

		TGraphErrors *PolThetaPhiCanvas_HX_gr = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_3,errpTmean,eypol_HX_3);
		PolThetaPhiCanvas_HX_gr->SetMarkerColor(kBlue);
		PolThetaPhiCanvas_HX_gr->SetMarkerStyle(8);
		PolThetaPhiCanvas_HX_gr->Draw("P");


		if (numberofRapbins==2){

		TGraphErrors *PolThetaPhiCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_3_rap2,errpTmean,eypol_HX_3_rap2);
		PolThetaPhiCanvas_HX_gr_rap2->SetMarkerColor(kRed);
		PolThetaPhiCanvas_HX_gr_rap2->SetMarkerStyle(8);
		PolThetaPhiCanvas_HX_gr_rap2->Draw("P");

		}

		if (numberofRapbins==3){

		TGraphErrors *PolThetaPhiCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_3_rap2,errpTmean,eypol_HX_3_rap2);
		PolThetaPhiCanvas_HX_gr_rap2->SetMarkerColor(kRed);
		PolThetaPhiCanvas_HX_gr_rap2->SetMarkerStyle(8);
		PolThetaPhiCanvas_HX_gr_rap2->Draw("P");

		TGraphErrors *PolThetaPhiCanvas_HX_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_3_rap3,errpTmean,eypol_HX_3_rap3);
		PolThetaPhiCanvas_HX_gr_rap3->SetMarkerColor(kBlack);
		PolThetaPhiCanvas_HX_gr_rap3->SetMarkerStyle(8);
		PolThetaPhiCanvas_HX_gr_rap3->Draw("P");

		}
		*/


/*		TCanvas *InvariantlambdaCanvas_CS = new TCanvas("InvariantlambdaCanvas_CS","InvariantlambdaCanvas_CS",1000,700);

				TH1F *InvariantlambdaCanvas_CS_TH1 = InvariantlambdaCanvas_CS->DrawFrame(0,-1.3,pTMax,1.3);
				InvariantlambdaCanvas_CS->SetXTitle("p^{J/#psi}_{T} [GeV/c]");
				InvariantlambdaCanvas_CS->SetYTitle("#lambda_{CS}_{Invariant}");
				InvariantlambdaCanvas_CS->GetFrame()->SetBorderSize(6);

				TGraphErrors *InvariantlambdaCanvas_CS_gr = new TGraphErrors(numberofpTbins,pTmean,invariantlambda_CS,errpTmean,errinvariantlambda_CS);
				InvariantlambdaCanvas_CS_gr->SetMarkerColor(kBlue);
				InvariantlambdaCanvas_CS_gr->SetMarkerStyle(8);
				InvariantlambdaCanvas_CS_gr->Draw("P");


				if (numberofRapbins==2){

				TGraphErrors *PolThetaCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1_rap2,errpTmean,eypol_HX_1_rap2);
				PolThetaCanvas_HX_gr_rap2->SetMarkerColor(kRed);
				PolThetaCanvas_HX_gr_rap2->SetMarkerStyle(8);
				PolThetaCanvas_HX_gr_rap2->Draw("P");

				}

				if (numberofRapbins==3){

				TGraphErrors *PolThetaCanvas_HX_gr_rap2 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1_rap2,errpTmean,eypol_HX_1_rap2);
				PolThetaCanvas_HX_gr_rap2->SetMarkerColor(kRed);
				PolThetaCanvas_HX_gr_rap2->SetMarkerStyle(8);
				PolThetaCanvas_HX_gr_rap2->Draw("P");

				TGraphErrors *PolThetaCanvas_HX_gr_rap3 = new TGraphErrors(numberofpTbins,pTmean,ypol_HX_1_rap3,errpTmean,eypol_HX_1_rap3);
				PolThetaCanvas_HX_gr_rap3->SetMarkerColor(kBlack);
				PolThetaCanvas_HX_gr_rap3->SetMarkerStyle(8);
				PolThetaCanvas_HX_gr_rap3->Draw("P");

				}
*/






		if (SavePlotIndex == 1){

		char oFilePolTheta_CS[200];
		sprintf(oFilePolTheta_CS,"MyFigures/Binning/CS_LambdaTheta_NpT%1.0f_NRap%1.0f_pTMin%1.0f_pTMax%1.0f_RapMin%1.0f_RapMax%1.0f.png", numberofpTbins, numberofRapbins,pTMin, pTMax, RapMin, RapMax);

		PolThetaCanvas_CS->SaveAs(oFilePolTheta_CS);


		char oFilePolPhi_CS[200];
		sprintf(oFilePolPhi_CS,"MyFigures/Binning/CS_LambdaPhi_NpT%1.0f_NRap%1.0f_pTMin%1.0f_pTMax%1.0f_RapMin%1.0f_RapMax%1.0f.png", numberofpTbins, numberofRapbins,pTMin, pTMax, RapMin, RapMax);

		PolPhiCanvas_CS->SaveAs(oFilePolPhi_CS);

		char oFilePolThetaPhi_CS[200];
		sprintf(oFilePolThetaPhi_CS,"MyFigures/Binning/CS_LambdaThetaPhi_NpT%1.0f_NRap%1.0f_pTMin%1.2f_pTMax%1.2f_RapMin%1.2f_RapMax%1.2f.png", numberofpTbins, numberofRapbins,pTMin, pTMax, RapMin, RapMax);

		//PolThetaPhiCanvas_CS->SaveAs(oFilePolThetaPhi_CS);


		char oFilePolTheta_HX[200];
		sprintf(oFilePolTheta_HX,"MyFigures/Binning/HX_LambdaTheta_NpT%1.0f_NRap%1.0f_pTMin%1.0f_pTMax%1.0f_RapMin%1.0f_RapMax%1.0f.png", numberofpTbins, numberofRapbins,pTMin, pTMax, RapMin, RapMax);

		PolThetaCanvas_HX->SaveAs(oFilePolTheta_HX);


		char oFilePolPhi_HX[200];
		sprintf(oFilePolPhi_HX,"MyFigures/Binning/HX_LambdaPhi_NpT%1.0f_NRap%1.0f_pTMin%1.0f_pTMax%1.0f_RapMin%1.0f_RapMax%1.0f.png", numberofpTbins, numberofRapbins,pTMin, pTMax, RapMin, RapMax);

		PolPhiCanvas_HX->SaveAs(oFilePolPhi_HX);

		char oFilePolThetaPhi_HX[200];
		sprintf(oFilePolThetaPhi_HX,"MyFigures/Binning/HX_LambdaThetaPhi_NpT%1.0f_NRap%1.0f_pTMin%1.2f_pTMax%1.2f_RapMin%1.2f_RapMax%1.2f.png", numberofpTbins, numberofRapbins,pTMin, pTMax, RapMin, RapMax);

		//PolThetaPhiCanvas_HX->SaveAs(oFilePolThetaPhi_HX);
		}
		cout<<"DrawFinal Complete"<<endl;


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ProduceSideBandTemplates(){

			sprintf(reducestr,"JpsiMass > %f & %f > JpsiMass",JpsiMass3SigmaMin,JpsiMass3SigmaMax); //  & JpsiMass > 3.05 & 3.15 > JpsiMass  & MCType == MCType::PR | MCType == MCType::NP
			redreddata = (RooDataSet*)reddata->reduce(reducestr);
			hDatared = new RooDataHist("hDatared","",RooArgSet(JpsiMass,Jpsict),*redreddata);
			hltIndex0PRred = (RooDataSet*)redreddata->reduce("MCType_idx < 0.5");
			hltIndex0NPred = (RooDataSet*)redreddata->reduce("MCType_idx > 0.5 && MCType_idx < 1.5");
			hltIndex0BKred = (RooDataSet*)redreddata->reduce("MCType_idx > 1.5 && MCType_idx < 2.5");
			hltIndex0PRredcut = (RooDataSet*)hltIndex0PRred->reduce("Jpsict < 0.25");
			hltIndex0NPredcut = (RooDataSet*)hltIndex0NPred->reduce("Jpsict < 0.25");
			hltIndex0PRredhist = new RooDataHist("hltIndex0PRredhist","",RooArgSet(JpsiMass,Jpsict),*hltIndex0PRred);
			hltIndex0BKredhist = new RooDataHist("hltIndex0BKredhist","",RooArgSet(JpsiMass,Jpsict),*hltIndex0BKred);

			redreddataMCtest = (RooDataSet*)reddataMCtest->reduce(reducestr);
			hDataredMCtest= new RooDataHist("hDataredMCtest","",RooArgSet(JpsiMass,Jpsict),*redreddataMCtest);
			hltIndex0PRredMCtest = (RooDataSet*)redreddataMCtest->reduce("MCType_idx < 0.5");
			hltIndex0NPredMCtest = (RooDataSet*)redreddataMCtest->reduce("MCType_idx > 0.5 && MCType_idx < 1.5");
			hltIndex0BKredMCtest = (RooDataSet*)redreddataMCtest->reduce("MCType_idx > 1.5 && MCType_idx < 2.5");
			hltIndex0PRredcutMCtest = (RooDataSet*)hltIndex0PRredMCtest->reduce("Jpsict < 0.25");
			hltIndex0NPredcutMCtest = (RooDataSet*)hltIndex0NPredMCtest->reduce("Jpsict < 0.25");
			hltIndex0PRredhistMCtest = new RooDataHist("hltIndex0PRredhistMCtest","",RooArgSet(JpsiMass,Jpsict),*hltIndex0PRredMCtest);
			hltIndex0NPredhistMCtest = new RooDataHist("hltIndex0NPredhistMCtest","",RooArgSet(JpsiMass,Jpsict),*hltIndex0NPredMCtest);
			hltIndex0BKredhistMCtest = new RooDataHist("hltIndex0BKredhistMCtest","",RooArgSet(JpsiMass,Jpsict),*hltIndex0BKredMCtest);

			redredMCNONPROMPT = (RooDataSet*)redMCNONPROMPT->reduce(reducestr);
			redredMCBACKGROUND = (RooDataSet*)redMCBACKGROUND->reduce(reducestr);


			sprintf(reducestr6,"%f < JpsiMass  &  JpsiMass < %f | %f < JpsiMass  &  JpsiMass < %f", JpsiMass6SigmaMin,JpsiMass3SigmaMin,JpsiMass3SigmaMax,JpsiMass6SigmaMax);
			realdataBkg = (RooDataSet*)realdata->reduce(reducestr6);
			realdataBkghist = new RooDataHist("realdataBkghist","realdataBkghist",RooArgSet(Jpsict,JpsiMass),*realdataBkg,1);

			sprintf(reducestr5,"%f < JpsiMass  &  JpsiMass < %f", JpsiMass3SigmaMax,JpsiMass6SigmaMax);
			realdataSBR = (RooDataSet*)realdata->reduce(reducestr5);
			realdataSBRhistfit = new RooDataHist("realdataSBRhistfit","realdataSBRhistfit",RooArgSet(Jpsict),*realdataSBR,1);
			NSBR = realdataSBR->sumEntries();
			sprintf(reducestr,"%f < JpsiMass  &  JpsiMass < %f", JpsiMass6SigmaMin,JpsiMass3SigmaMin);
			realdataSBL = (RooDataSet*)realdata->reduce(reducestr);
			realdataSBLhistfit = new RooDataHist("realdataSBLhistfit","realdataSBLhistfit",RooArgSet(Jpsict),*realdataSBL,1);
			NSBL = realdataSBL->sumEntries();
			if (NSBL>0){
			SBscale=NSBR/NSBL;
			}
			else if (NSBL==0){
			SBscale = 1.0;
			}
			else if (NSBL<0){
			SBscale = 1.0;
			}

			realdataSBLscale0 = realdataSBL->createHistogram("realdataSBLscale0",Jpsict,Binning(100)) ;
			realdataSBLscale0->Scale(SBscale);
			realdataSBLscale = new RooDataHist("realdataSBLscale","",Jpsict,realdataSBLscale0);


			RooBinning rb(-1,2.5);
			rb.addUniform(40,-1,2.5);
			Jpsict.setBinning(rb);

			realdataSBLhist = new RooDataHist("realdataSBLhist","realdataSBLhist",RooArgSet(Jpsict),*realdataSBL,1);
			RooHistPdf histrealSBL("histrealSBL","histrealSBL",Jpsict,*realdataSBLhist,5);
			RooRealVar coefctrealSBL("coefctrealSBL","coefctrealSBL",0.5);

			realdataSBRhist = new RooDataHist("realdataSBRhist","realdataSBRhist",RooArgSet(Jpsict),*realdataSBR,1);
			RooHistPdf histrealSBR("histrealSBR","histrealSBR",Jpsict,*realdataSBRhist,5);
			RooRealVar coefctrealSBR("coefctrealSBR","coefctrealSBR",0.5);
			RooAddPdf histrealSB("histrealSB","histrealSB",RooArgList(histrealSBL,histrealSBR),RooArgList(coefctrealSBL,coefctrealSBR));

			RooBinning rb(-1,2.5);
			rb.addUniform(100,-1,2.5);
			Jpsict.setBinning(rb);

			MCdataBkg = (RooDataSet*)reddataMCtest->reduce(reducestr6);
			MCdataBkghist = new RooDataHist("MCdataBkghist","MCdataBkghist",RooArgSet(Jpsict,JpsiMass),*MCdataBkg,1);
			MCdataSBL = (RooDataSet*)reddata->reduce(reducestr);
			MCdataSBR = (RooDataSet*)reddata->reduce(reducestr5);

			RooBinning rb(-1,2.5);
			rb.addUniform(50,-1,2.5);
			Jpsict.setBinning(rb);

			MCdataBkghistctau = new RooDataHist("MCdataBkghist","MCdataBkghist",RooArgSet(Jpsict),*MCdataBkg,1);
			RooHistPdf histMCSB("histMCSB","histMCSB",Jpsict,*MCdataBkghistctau,2);

			RooBinning rb(-1,2.5);
			rb.addUniform(100,-1,2.5);
			Jpsict.setBinning(rb);

			sprintf(reduceSignalwindowData,"%f < JpsiMass  &  JpsiMass < %f", JpsiMass3SigmaMin,JpsiMass3SigmaMax);

			redrealdata = (RooDataSet*)realdata->reduce(reduceSignalwindowData);
			redrealhData = new RooDataHist("redrealhData","redrealhData",RooArgSet(JpsiMass,Jpsict),*redrealdata);


			RooBinning rb(-1,2.5);
			rb.addUniform(100,-1,2.5);
			Jpsict.setBinning(rb);

			hltIndex0PRhistbin = new RooDataHist("hltIndex0PRhistbin","hltIndex0PRhistbin",RooArgSet(Jpsict),*hltIndex0PR,1);
			RooHistPdf histPrompt("histPrompt","histPrompt",Jpsict,*hltIndex0PRhistbin,2);
			RooRealVar coefctPrompt("coefctPrompt","coefctPrompt",1,0,200000);
			RooExtendPdf ehistPrompt("ehistPrompt","ehistPrompt",histPrompt,coefctPrompt);

			RooBinning rb(-1,2.5);
			rb.addUniform(35,-1,2.5);
			Jpsict.setBinning(rb);

			hltIndex0BKhistbin = new RooDataHist("hltIndex0BKhistbin","hltIndex0BKhistbin",RooArgSet(Jpsict),*redredMCBACKGROUND,1);
			RooHistPdf histBackground("histBackground","histBackground",Jpsict,*hltIndex0BKhistbin,2);
			RooRealVar coefctBackground("coefctBackground","coefctBackground",NbkgRange);
			RooExtendPdf ehistBackground("ehistBackground","ehistBackground",histBackground,coefctBackground);

			hltIndex0NPhistbin = new RooDataHist("hltIndex0NPhistbin","hltIndex0NPhistbin",RooArgSet(Jpsict),*redredMCNONPROMPT,1);
			RooHistPdf histNonPrompt("histNonPrompt","histNonPrompt",Jpsict,*hltIndex0NPhistbin,2);
			RooRealVar coefctNonPrompt("coefctNonPrompt","coefctNonPrompt",100,1,100000);
			RooExtendPdf ehistNonPrompt("ehistNonPrompt","ehistNonPrompt",histNonPrompt,coefctNonPrompt);

			RooBinning rb(-1,2.5);
			rb.addUniform(100,-1,2.5);
			Jpsict.setBinning(rb);

			RooExtendPdf ehistMCSB("ehistMCSB","ehistMCSB",histMCSB,coefctBackground);
			RooExtendPdf ehistrealSB("ehistrealSB","ehistrealSB",histrealSB,coefctBackground);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AnalyticalBKandNPModels(){

///////////////////////////////// FIX THE ANALYTICAL PROMPT MODEL PARAMETERS FOR BK AND NP FITS //////////////////////////////////////

// Gauss 1
	RooRealVar ctmeanSig1fixnotfix("ctmeanSig1fixnotfix","ctmeanSig1fixnotfix",0,-1,1);//0.0292);
	RooRealVar ctmeanSig1fix("ctmeanSig1fix","ctmeanSig1fix",0);
	RooRealVar ctsigmaSig1fix("ctsigmaSig1fix","ctsigmaSig1fix",ctsigmaSig1double);
	RooRealVar ctcoefSig1fix("ctcoefSig1fix","ctcoefSig1fix",ctcoefSig1doublenorm);
	RooGaussModel ctgaussSig1fixnotfix("ctgaussSig1fixnotfix","ctgaussSig1fixnotfix",Jpsict,ctmeanSig1fixnotfix,ctsigmaSig1fix);
	RooGaussModel ctgaussSig1fix("ctgaussSig1fix","ctgaussSig1fix",Jpsict,ctmeanSig1fix,ctsigmaSig1fix);
// Gauss 2 (same mean as Gauss 1)
	RooRealVar ctsigmaSig2fix("ctsigmaSig2fix","ctsigmaSig2fix",ctsigmaSig2double);
	RooRealVar ctcoefSig2fix("ctcoefSig2fix","ctcoefSig2fix",ctcoefSig2doublenorm);
	RooGaussModel ctgaussSig2fix("ctgaussSig2fix","ctgaussSig2fix",Jpsict,ctmeanSig1fix,ctsigmaSig2fix);
	RooGaussModel ctgaussSig2fixnotfix("ctgaussSig2fixnotfix","ctgaussSig2fixnotfix",Jpsict,ctmeanSig1fixnotfix,ctsigmaSig2fix);
//Gauss 3
//	RooRealVar ctmeanSig3fix("ctmeanSig3fix","ctmeanSig3fix",0,-0.1,0.3);
	RooRealVar ctsigmaSig3fix("ctsigmaSig3fix","ctsigmaSig3fix",ctsigmaSig3double);
	RooRealVar ctcoefSig3fix("ctcoefSig3fix","ctcoefSig3fix",ctcoefSig3doublenorm);
	RooGaussModel ctgaussSig3fix("ctgaussSig3fix","ctgaussSig3fix",Jpsict,ctmeanSig1fix,ctsigmaSig3fix);
	RooGaussModel ctgaussSig3fixnotfix("ctgaussSig3fixnotfix","ctgaussSig3fixnotfix",Jpsict,ctmeanSig1fixnotfix,ctsigmaSig3fix);



	RooAddModel sumJpsictPromptfix("sumJpsictPromptfix","sumJpsictPromptfix",RooArgList(ctgaussSig1fix,ctgaussSig2fix,ctgaussSig3fix),RooArgList(ctcoefSig1fix,ctcoefSig2fix,ctcoefSig3fix));
	RooAddModel sumJpsictPromptfixnotfix("sumJpsictPromptfixnotfix","sumJpsictPromptfixnotfix",RooArgList(ctgaussSig1fixnotfix,ctgaussSig2fixnotfix,ctgaussSig3fixnotfix),RooArgList(ctcoefSig1fix,ctcoefSig2fix,ctcoefSig3fix));

	RooRealVar ctcoefPromptfix("ctcoefPromptfix","ctcoefPromptfix",1,0,100000);
	RooExtendPdf esumJpsictPromptfix("esumJpsictPromptfix","esumJpsictPromptfix",sumJpsictPromptfix,ctcoefPromptfix);
	RooExtendPdf esumJpsictPromptfixnotfix("esumJpsictPromptfixnotfix","esumJpsictPromptfixnotfix",sumJpsictPromptfixnotfix,ctcoefPromptfix);

//////////////////////////////// NON PROMPT ANALYTICAL MODEL /////////////////////////////////////////////////////////////////////
	// Decay Model 2
		RooRealVar tau2NonPrompt("tau2NonPrompt","tau2NonPrompt",1,0,5);//0.39);
		RooDecay NonPrompt2("NonPrompt2","NonPrompt2",Jpsict,tau2NonPrompt,sumJpsictPromptfixnotfix,RooDecay::SingleSided);
		RooRealVar coefctNonPrompt2("coefctNonPrompt2","coefctNonPrompt2",0.5,0,1);//0.78);
		RooExtendPdf eNonPrompt2("eNonPrompt2","eNonPrompt2",NonPrompt2,coefctNonPrompt2);
	// Decay Model 3
		RooRealVar tau3NonPrompt("tau3NonPrompt","tau3NonPrompt",1,0,5);//0.0273);
		RooDecay NonPrompt3("NonPrompt3","NonPrompt3",Jpsict,tau3NonPrompt,sumJpsictPromptfixnotfix,RooDecay::Flipped);
		RooRealVar coefctNonPrompt3("coefctNonPrompt3","coefctNonPrompt3",1,0,100000);
		RooExtendPdf eNonPrompt3("eNonPrompt3","eNonPrompt3",NonPrompt3,coefctNonPrompt3);
	// Decay Model 4
		RooRealVar tau4NonPrompt("tau4NonPrompt","tau4NonPrompt",1,0,5);
		RooDecay NonPrompt4("NonPrompt4","NonPrompt4",Jpsict,tau4NonPrompt,sumJpsictPromptfixnotfix,RooDecay::DoubleSided);
		RooRealVar coefctNonPrompt4("coefctNonPrompt4","coefctNonPrompt4",1,0,100000);
		RooExtendPdf eNonPrompt4("eNonPrompt4","eNonPrompt4",NonPrompt4,coefctNonPrompt4);

		RooAddPdf NonPromptctTOT("NonPromptctTOT","NonPromptctTOT",RooArgList(NonPrompt2,NonPrompt3),RooArgList(coefctNonPrompt2));

////////////////////////////////// BACKGROUND ANALYTICAL MODEL ////////////////////////////////////////////////////////////////////////

	// Decay Model 2
		RooRealVar tau2("tau2","tau2",0.168,0.03,0.5);
		RooDecay bkg2("bkg2","bkg2",Jpsict,tau2,sumJpsictPromptfix,RooDecay::SingleSided);
		RooRealVar coefctbkg2("coefctbkg2","coefctbkg2",0.43,0.2,0.6);
		RooExtendPdf ebkg2("ebkg2","ebkg2",bkg2,coefctbkg2);
	// Decay Model 3
		RooRealVar tau3("tau3","tau3",0.128,0.03,0.2);
		RooDecay bkg3("bkg3","bkg3",Jpsict,tau3,sumJpsictPromptfix,RooDecay::Flipped);
		RooRealVar coefctbkg3("coefctbkg3","coefctbkg3",0.34,0.2,0.6);
		RooExtendPdf ebkg3("ebkg3","ebkg3",bkg3,coefctbkg3);
	// Decay Model 4
		RooRealVar tau4("tau4","tau4",0.787,0.05,0.99);
		RooDecay bkg4("bkg4","bkg4",Jpsict,tau4,sumJpsictPromptfix,RooDecay::DoubleSided);
		RooRealVar coefctbkg4("coefctbkg4","coefctbkg4",0.5,0,1);
		RooExtendPdf ebkg4("ebkg4","ebkg4",bkg4,coefctbkg4);

		RooAddPdf bkgctTOT("bkgctTOT","bkgctTOT",RooArgList(bkg2,bkg3,bkg4),RooArgList(coefctbkg2,coefctbkg3));

		RooExtendPdf ebkgctTOT("ebkgctTOT","ebkgctTOT",bkgctTOT,coefctBackground);
		RooRealVar coefctNonPromptanal("coefctNonPromptanal","coefctNonPromptanal",100,50,100000);
		RooExtendPdf eNonPromptctTOT("eNonPromptctTOT","eNonPromptctTOT",NonPromptctTOT,coefctNonPromptanal);


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitGlobalAndAnalytical(){


	if(RealDataIndex == 0){//fit to MC

		RooAddPdf sumJpsict("sumJpsict","sumJpsict",RooArgList(esumJpsictPromptfix,ehistNonPrompt,ehistMCSB));

		if (AnalyticIndex == 1){

		fitresctPromptTemplate = ehistPrompt.fitTo(*hltIndex0PRredhistMCtest,Save(1),Minos(0),SumW2Error(kFALSE));
		fitresctNonPromptTemplate = ehistNonPrompt.fitTo(*hltIndex0NPredhistMCtest,Save(1),Minos(0),SumW2Error(kFALSE));
		fitresctBackgroundTemplate = histMCSB.fitTo(*hltIndex0BKredhistMCtest,Save(1),Minos(0),SumW2Error(kFALSE));

		fitresctNonPrompt = NonPromptctTOT.fitTo(*hltIndex0NPredhistMCtest,Save(1),Minos(0),SumW2Error(kFALSE));
		fitresctNonPrompt->Print();

		JpsictframeNonPromptAnal = Jpsict.frame(-1,2.5) ;
		hltIndex0NPredhistMCtest->plotOn(JpsictframeNonPromptAnal,DataError(RooAbsData::SumW2));
		NonPromptctTOT.plotOn(JpsictframeNonPromptAnal,Normalization(1.0));
		NonPromptctTOT.plotOn(JpsictframeNonPromptAnal, RooFit::Components(RooArgList(NonPrompt2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		NonPromptctTOT.plotOn(JpsictframeNonPromptAnal, RooFit::Components(RooArgList(NonPrompt3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

		JpsictframeNonPromptAnalChi2 = Jpsict.frame(-1,2.5,1000) ;
		hltIndex0NPredhistMCtest->plotOn(JpsictframeNonPromptAnalChi2,DataError(RooAbsData::SumW2));
		NonPromptctTOT.plotOn(JpsictframeNonPromptAnalChi2,Normalization(1.0));


		fitresctBackground = bkgctTOT.fitTo(*hltIndex0BKredhistMCtest,Save(1),Minos(0),SumW2Error(kFALSE));
		}

		fitresGlobal = sumJpsict.fitTo(*hDataredMCtest,Save(1),Minos(0),SumW2Error(kFALSE));

	}



	else if (RealDataIndex == 1){

		RooAddPdf sumJpsict("sumJpsict","sumJpsict",RooArgList(esumJpsictPromptfix,ehistNonPrompt,ehistrealSB));
		fitresGlobal = sumJpsict.fitTo(*redrealhData,Save(1),Minos(0),SumW2Error(kFALSE));

		if (SBComparisonIndex ==1){

		fitresctBackgroundRealDataSBR = bkgctTOT.fitTo(*realdataSBRhistfit,Save(1),Minos(0),SumW2Error(kFALSE));

		JpsictframeBackgroundAnalRealDataSBR = Jpsict.frame(-1,2.5) ;
		realdataSBRhistfit->plotOn(JpsictframeBackgroundAnalRealDataSBR,DataError(RooAbsData::SumW2));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBR,Normalization(1.0));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBR, RooFit::Components(RooArgList(bkg2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBR, RooFit::Components(RooArgList(bkg3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBR, RooFit::Components(RooArgList(bkg4)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

		JpsictframeBackgroundAnalRealDataSBRchi2 = Jpsict.frame(-1,2.5) ;
		realdataSBRhistfit->plotOn(JpsictframeBackgroundAnalRealDataSBRchi2,DataError(RooAbsData::SumW2));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBRchi2,Normalization(1.0));

		cout << "chi2/ndof (BackgroundctAnalRealDataSLR) = " << JpsictframeBackgroundAnalRealDataSBRchi2->chiSquare() << endl;

		fitresctBackgroundRealDataSBR->Print();

		fitresctBackgroundRealDataSBL = bkgctTOT.fitTo(*realdataSBLhistfit,Save(1),Minos(0),SumW2Error(kFALSE));

		JpsictframeBackgroundAnalRealDataSBL = Jpsict.frame(-1,2.5) ;
		realdataSBLhistfit->plotOn(JpsictframeBackgroundAnalRealDataSBL,DataError(RooAbsData::SumW2));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBL,Normalization(1.0));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBL, RooFit::Components(RooArgList(bkg2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBL, RooFit::Components(RooArgList(bkg3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBL, RooFit::Components(RooArgList(bkg4)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

		JpsictframeBackgroundAnalRealDataSBLchi2 = Jpsict.frame(-1,2.5) ;
		realdataSBLhistfit->plotOn(JpsictframeBackgroundAnalRealDataSBLchi2,DataError(RooAbsData::SumW2));
		bkgctTOT.plotOn(JpsictframeBackgroundAnalRealDataSBLchi2,Normalization(1.0));

		cout << "chi2/ndof (BackgroundctAnalRealDataSLL) = " << JpsictframeBackgroundAnalRealDataSBLchi2->chiSquare() << endl;

		fitresctBackgroundRealDataSBL->Print();
		}
	}


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Plot(){


	JpsiMassframe = JpsiMass.frame() ;

	if (RealDataIndex == 0){


		hDataMCtest->plotOn(JpsiMassframe,DataError(RooAbsData::SumW2));
		sumJpsiMass2.plotOn(JpsiMassframe,Normalization(1.0)) ;
		sumJpsiMass2.plotOn(JpsiMassframe, RooFit::Components(RooArgList(CBMass)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		sumJpsiMass2.plotOn(JpsiMassframe, RooFit::Components(RooArgList(gauss)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		sumJpsiMass2.plotOn(JpsiMassframe, RooFit::Components(RooArgList(eexpo)), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));

		JpsiMassframeChi2 = JpsiMass.frame(2.7,3.5) ;
		hDataMCtest->plotOn(JpsiMassframeChi2,DataError(RooAbsData::SumW2));
		sumJpsiMass2.plotOn(JpsiMassframeChi2,Normalization(1.0)) ;
	}

///////////////// DATA PLOTS //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	else if(RealDataIndex == 1){
		realhData->plotOn(JpsiMassframe,DataError(RooAbsData::SumW2));
		sumJpsiMass2.plotOn(JpsiMassframe,Normalization(1.0)) ;
		sumJpsiMass2.plotOn(JpsiMassframe, RooFit::Components(RooArgList(CBMass)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		sumJpsiMass2.plotOn(JpsiMassframe, RooFit::Components(RooArgList(gauss)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
		sumJpsiMass2.plotOn(JpsiMassframe, RooFit::Components(RooArgList(eexpo)), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));

		JpsiMassframeChi2 = JpsiMass.frame(2.7,3.5,1000) ;
		realhData->plotOn(JpsiMassframeChi2,DataError(RooAbsData::SumW2));
		sumJpsiMass2.plotOn(JpsiMassframeChi2,Normalization(1.0)) ;
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////------------------------------------------------////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////--------c-Tau Draw--------///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////-----------------------------------------------/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	if(RealDataIndex == 0){

	if (AnalyticIndex == 1){


	JpsictframePrompt = Jpsict.frame(-1,2.5) ;
	hltIndex0PRredhistMCtest->plotOn(JpsictframePrompt,DataError(RooAbsData::SumW2));
	ehistPrompt.plotOn(JpsictframePrompt,Normalization(1.0));

	JpsictframePromptAnal = Jpsict.frame(-1,2.5) ;
	hltIndex0PRhist->plotOn(JpsictframePromptAnal,DataError(RooAbsData::SumW2));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal,Normalization(1.0));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal, RooFit::Components(RooArgList(ctgaussSig1)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal, RooFit::Components(RooArgList(ctgaussSig2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	esumJpsictPrompt.plotOn(JpsictframePromptAnal, RooFit::Components(RooArgList(ctgaussSig3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

	JpsictframePromptAnalChi2 = Jpsict.frame(-1,2.5,1000) ;
	hltIndex0PRhist->plotOn(JpsictframePromptAnalChi2,DataError(RooAbsData::SumW2));
	esumJpsictPrompt.plotOn(JpsictframePromptAnalChi2,Normalization(1.0));

	JpsictframeBackground = Jpsict.frame(-1,2.5) ;
	hltIndex0BKredhistMCtest->plotOn(JpsictframeBackground,DataError(RooAbsData::SumW2));
	histMCSB.plotOn(JpsictframeBackground,Normalization(1.0));

	JpsictframeBackgroundAnal = Jpsict.frame(-1,2.5) ;
	hltIndex0BKredhistMCtest->plotOn(JpsictframeBackgroundAnal,DataError(RooAbsData::SumW2));
	bkgctTOT.plotOn(JpsictframeBackgroundAnal,Normalization(1.0));
	bkgctTOT.plotOn(JpsictframeBackgroundAnal, RooFit::Components(RooArgList(bkg2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	bkgctTOT.plotOn(JpsictframeBackgroundAnal, RooFit::Components(RooArgList(bkg3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	bkgctTOT.plotOn(JpsictframeBackgroundAnal, RooFit::Components(RooArgList(bkg4)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

	JpsictframeBackgroundAnalChi2 = Jpsict.frame(-1,2.5,1000) ;
	hltIndex0BKredhistMCtest->plotOn(JpsictframeBackgroundAnalChi2,DataError(RooAbsData::SumW2));
	bkgctTOT.plotOn(JpsictframeBackgroundAnalChi2,Normalization(1.0));

	JpsictframeNonPrompt = Jpsict.frame(-1,2.5) ;
	hltIndex0NPredhistMCtest->plotOn(JpsictframeNonPrompt,DataError(RooAbsData::SumW2));
	ehistNonPrompt.plotOn(JpsictframeNonPrompt,Normalization(1.0));

	JpsictframeNonPromptAnal = Jpsict.frame(-1,2.5) ;
	hltIndex0NPredhistMCtest->plotOn(JpsictframeNonPromptAnal,DataError(RooAbsData::SumW2));
	NonPromptctTOT.plotOn(JpsictframeNonPromptAnal,Normalization(1.0));
	NonPromptctTOT.plotOn(JpsictframeNonPromptAnal, RooFit::Components(RooArgList(NonPrompt2)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	NonPromptctTOT.plotOn(JpsictframeNonPromptAnal, RooFit::Components(RooArgList(NonPrompt3)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));
	NonPromptctTOT.plotOn(JpsictframeNonPromptAnal, RooFit::Components(RooArgList(NonPrompt4)), RooFit::LineStyle(kDashed),RooFit::LineColor(kBlack));

	JpsictframeNonPromptAnalChi2 = Jpsict.frame(-1,2.5,1000) ;
	hltIndex0NPredhistMCtest->plotOn(JpsictframeNonPromptAnalChi2,DataError(RooAbsData::SumW2));
	NonPromptctTOT.plotOn(JpsictframeNonPromptAnalChi2,Normalization(1.0));


	}
	}


	else if(RealDataIndex == 1){
	JpsictframeDataChi2 = Jpsict.frame(-1,2.5) ;
	redrealhData->plotOn(JpsictframeDataChi2,DataError(RooAbsData::SumW2));
	sumJpsict.plotOn(JpsictframeDataChi2,Normalization(1.0));

	JpsictframeSB = Jpsict.frame(-1,2.5) ;
	realdataSBLscale->plotOn(JpsictframeSB,DataError(RooAbsData::SumW2));
	realdataSBR->plotOn(JpsictframeSB,DataError(RooAbsData::SumW2),MarkerColor(2));

	JpsictframeChi2 = Jpsict.frame(-1,2.5,1000) ;
	realhData->plotOn(JpsictframeChi2,DataError(RooAbsData::SumW2));
	sumJpsict.plotOn(JpsictframeChi2,Normalization(1.0));

	}


	JpsictframeChi2 = Jpsict.frame(-1,2.5) ;
	hDataredMCtest->plotOn(JpsictframeChi2,DataError(RooAbsData::SumW2));
	sumJpsict.plotOn(JpsictframeChi2,Normalization(1.0));

	chi2ctGlobal;
	if (RealDataIndex == 0){
	cout << "chi2/ndof (ctGlobalMC) = " << JpsictframeChi2->chiSquare() << endl;
	chi2ctGlobal=JpsictframeChi2->chiSquare();
	}

///////////////// MC PLOTS ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Jpsictframe = Jpsict.frame(-1,2.5) ;

	if(RealDataIndex == 0){

			hltIndex0PRredMCtest->plotOn(Jpsictframe,DataError(RooAbsData::SumW2),RooFit::LineColor(kRed));
			hltIndex0NPredMCtest->plotOn(Jpsictframe,DataError(RooAbsData::SumW2),RooFit::LineColor(kGreen));
			hltIndex0BKredMCtest->plotOn(Jpsictframe,DataError(RooAbsData::SumW2),RooFit::LineColor(kBlack));
			hDataredMCtest->plotOn(Jpsictframe,DataError(RooAbsData::SumW2));

			sumJpsict.plotOn(Jpsictframe,Normalization(1.0));


	sumJpsict.plotOn(Jpsictframe, RooFit::Components(RooArgList(esumJpsictPromptfix)), RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
	sumJpsict.plotOn(Jpsictframe, RooFit::Components(RooArgList(ehistNonPrompt)), RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
	sumJpsict.plotOn(Jpsictframe, RooFit::Components(RooArgList(ehistMCSB)), RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));


	}

///////////////// DATA PLOTS //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	else if(RealDataIndex == 1){
		redrealhData->plotOn(Jpsictframe,DataError(RooAbsData::SumW2));

		sumJpsict.plotOn(Jpsictframe,Normalization(1.0),RooFit::LineColor(kBlue));

		sumJpsict.plotOn(Jpsictframe, RooFit::Components(RooArgList(esumJpsictPromptfix)), RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
		sumJpsict.plotOn(Jpsictframe, RooFit::Components(RooArgList(ehistNonPrompt)), RooFit::LineColor(kGreen),RooFit::LineStyle(kDashed));
		sumJpsict.plotOn(Jpsictframe, RooFit::Components(RooArgList(ehistrealSB)), RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));
	}

}

void FillArrays(){

	xploterror=(pTMax-pTMin)/2;



	nArray=npT-NpTBinStart;

	pTmean[nArray]=pTWCentre_rap[nRap-1][npT-1];
	errpTmean[nArray]=0.1;





	if (AngularIndex == 1){

	lambda_theta_CS_double = lambda_theta_CS.getVal();
	lambda_phi_CS_double = lambda_phi_CS.getVal();
	lambda_thetaphi_CS_double = lambda_thetaphi_CS.getVal();
	errlambda_theta_CS_double = lambda_theta_CS.getError();
	errlambda_phi_CS_double = lambda_phi_CS.getError();
	errlambda_thetaphi_CS_double = lambda_thetaphi_CS.getError();

	lambda_theta_HX_double = lambda_theta_HX.getVal();
	lambda_phi_HX_double = lambda_phi_HX.getVal();
	lambda_thetaphi_HX_double = lambda_thetaphi_HX.getVal();
	errlambda_theta_HX_double = lambda_theta_HX.getError();
	errlambda_phi_HX_double = lambda_phi_HX.getError();
	errlambda_thetaphi_HX_double = lambda_thetaphi_HX.getError();

	invariantlambda_CS=(lambda_theta_CS_double+3*lambda_phi_CS_double)/(1-lambda_phi_CS_double);
	invariantlambda_HX=(lambda_theta_HX_double+3*lambda_phi_HX_double)/(1-lambda_phi_HX_double);
	invariantF_CS=(1+lambda_theta_CS_double+2*lambda_phi_CS_double)/(3-lambda_theta_CS_double);
	invariantF_HX=(1+lambda_theta_HX_double+2*lambda_phi_HX_double)/(3-lambda_theta_HX_double);



	printf("\n%1.1f < pT < %1.1f GeV\n", pTMin, pTMax);
	printf("%1.1f < Rap < %1.1f\n", RapMin, RapMax);

	cout<<"POLARIZATION PARAMETERS:"<<endl;

	printf("\lambda_theta_CS = %1.3f", lambda_theta_CS_double);
	printf(" +- %1.3f \n", errlambda_theta_CS_double);
	printf("\lambda_phi_CS = %1.3f", lambda_phi_CS_double);
	printf(" +- %1.3f \n", errlambda_phi_CS_double);
	printf("\lambda_theta_HX = %1.3f", lambda_theta_HX_double);
	printf(" +- %1.3f \n", errlambda_theta_HX_double);
	printf("\lambda_phi_HX = %1.3f", lambda_phi_HX_double);
	printf(" +- %1.3f \n\n", errlambda_phi_HX_double);


	cout<<"INVARIANT POLARIZATION PARAMETERS:"<<endl;
	printf("\CS_lambda = %1.2f\n", invariantlambda_CS);
	printf("\HX_lambda = %1.2f\n", invariantlambda_HX);
	printf("\CS_F = %1.2f\n", invariantF_CS);
	printf("\HX_F = %1.2f\n", invariantF_HX);

	fprintf(outputFile2, "POLARIZATION PARAMETERS:\n");
	fprintf(outputFile2, "\lambda_theta_CS = %1.3f", lambda_theta_CS_double);
	fprintf(outputFile2, " +- %1.3f \n", errlambda_theta_CS_double);
	fprintf(outputFile2, "\lambda_phi_CS = %1.3f", lambda_phi_CS_double);
	fprintf(outputFile2, " +- %1.3f \n", errlambda_phi_CS_double);
	fprintf(outputFile2, "\lambda_theta_HX = %1.3f", lambda_theta_HX_double);
	fprintf(outputFile2, " +- %1.3f \n", errlambda_theta_HX_double);
	fprintf(outputFile2, "\lambda_phi_HX = %1.3f", lambda_phi_HX_double);
	fprintf(outputFile2, " +- %1.3f \n\n", errlambda_phi_HX_double);

	fprintf(outputFile2, "INVARIANT POLARIZATION PARAMETERS:\n");

	fprintf(outputFile2, "\CS_lambda = %1.2f\n", invariantlambda_CS);
	fprintf(outputFile2, "\HX_lambda = %1.2f\n", invariantlambda_HX);
	fprintf(outputFile2, "\CS_F = %1.2f\n", invariantF_CS);
	fprintf(outputFile2, "\HX_F = %1.2f\n", invariantF_HX);

	errinvariantlambda_CS=0.1;
	errinvariantlambda_HX=0.1;
	errinvariantF_CS=0.1;
	errinvariantF_HX=0.1;


	if (nRap==1){

	ypol_CS_1[nArray]=lambda_theta_CS_double;
	eypol_CS_1[nArray]=errlambda_theta_CS_double;

	ypol_CS_2[nArray]=lambda_phi_CS_double;
	eypol_CS_2[nArray]=errlambda_phi_CS_double;

	ypol_CS_3[nArray]=lambda_thetaphi_CS_double;
	eypol_CS_3[nArray]=errlambda_thetaphi_CS_double;


	ypol_HX_1[nArray]=lambda_theta_HX_double;
	eypol_HX_1[nArray]=errlambda_theta_HX_double;

	ypol_HX_2[nArray]=lambda_phi_HX_double;
	eypol_HX_2[nArray]=errlambda_phi_HX_double;

	ypol_HX_3[nArray]=lambda_thetaphi_HX_double;
	eypol_HX_3[nArray]=errlambda_thetaphi_HX_double;

	invariantlambda_CS_rap1[nArray]=invariantlambda_CS;
	invariantlambda_HX_rap1[nArray]=invariantlambda_HX;
	invariantF_CS_rap1[nArray]=invariantF_CS;
	invariantF_HX_rap1[nArray]=invariantF_HX;

	errinvariantlambda_CS_rap1[nArray]=errinvariantlambda_CS;
	errinvariantlambda_HX_rap1[nArray]=errinvariantlambda_HX;
	errinvariantF_CS_rap1[nArray]=errinvariantF_CS;
	errinvariantF_HX_rap1[nArray]=errinvariantF_HX;
	}

	if (nRap==2){

	ypol_CS_1_rap2[nArray]=lambda_theta_CS_double;
	eypol_CS_1_rap2[nArray]=errlambda_theta_CS_double;

	ypol_CS_2_rap2[nArray]=lambda_phi_CS_double;
	eypol_CS_2_rap2[nArray]=errlambda_phi_CS_double;

	ypol_CS_3_rap2[nArray]=lambda_thetaphi_CS_double;
	eypol_CS_3_rap2[nArray]=errlambda_thetaphi_CS_double;


	ypol_HX_1_rap2[nArray]=lambda_theta_HX_double;
	eypol_HX_1_rap2[nArray]=errlambda_theta_HX_double;

	ypol_HX_2_rap2[nArray]=lambda_phi_HX_double;
	eypol_HX_2_rap2[nArray]=errlambda_phi_HX_double;

	ypol_HX_3_rap2[nArray]=lambda_thetaphi_HX_double;
	eypol_HX_3_rap2[nArray]=errlambda_thetaphi_HX_double;
	}

	if (nRap==3){

	ypol_CS_1_rap3[nArray]=lambda_theta_CS_double;
	eypol_CS_1_rap3[nArray]=errlambda_theta_CS_double;

	ypol_CS_2_rap3[nArray]=lambda_phi_CS_double;
	eypol_CS_2_rap3[nArray]=errlambda_phi_CS_double;

	ypol_CS_3_rap3[nArray]=lambda_thetaphi_CS_double;
	eypol_CS_3_rap3[nArray]=errlambda_thetaphi_CS_double;


	ypol_HX_1_rap3[nArray]=lambda_theta_HX_double;
	eypol_HX_1_rap3[nArray]=errlambda_theta_HX_double;

	ypol_HX_2_rap3[nArray]=lambda_phi_HX_double;
	eypol_HX_2_rap3[nArray]=errlambda_phi_HX_double;

	ypol_HX_3_rap3[nArray]=lambda_thetaphi_HX_double;
	eypol_HX_3_rap3[nArray]=errlambda_thetaphi_HX_double;
	}


	}


	if (nRap==1){

	yrap1[nArray]=BFRACTION;
	eyrap1[nArray]=ERRBFRACTION;

	if (RealDataIndex == 0){

	yrap1MCTRUTH[nArray]=bfractionMCTRUTH;
	eyrap1MCTRUTH[nArray]=0.0001;
	}
	}

	else if (nRap==2){


	yrap2[nArray]=BFRACTION;
	eyrap2[nArray]=ERRBFRACTION;

	if (RealDataIndex == 0){

	yrap2MCTRUTH[nArray]=bfractionMCTRUTH;
	eyrap2MCTRUTH[nArray]=0.0001;
				}
	}

	else if (nRap==3){


	yrap3[nArray]=BFRACTION;
	eyrap3[nArray]=ERRBFRACTION;

	if (RealDataIndex == 0){

	yrap3MCTRUTH[nArray]=bfractionMCTRUTH;
	eyrap3MCTRUTH[nArray]=0.0001;
				}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CalcAndPrintResults(){

	fprintf(outputFile2, "npT = %1.0f, ",npT);
			fprintf(outputFile2, "pTMin = %1.2f, ",pTMin);
			fprintf(outputFile2, "pTMax = %1.2f\n",pTMax);

			fprintf(outputFile2, "nRap = %1.0f, ",nRap);

			fprintf(outputFile2, "RapMin = %1.2f, ",RapMin);
			fprintf(outputFile2, "RapMax = %1.2f\n",RapMax);

			printf("\nNEW BIN\n");
			numberofprocessedbins = (npT-1)*numberofRapbins+nRap;
			printf("Processing pT bin number %1.0f and Rap bin number %1.0f (%1.0f, %1.0f)\n in total: processing bin number %1.0f out of %1.0f\n", npT, nRap, numberofpTbins, numberofRapbins, numberofprocessedbins, numberofbins);
			printf("nbin = %1.0f, npT = %1.0f, nRap = %1.0f\n", nbin, npT, nRap);
			printf("pTMin und Max = %f %f\n", pTMin, pTMax);
			printf("RapMin und Max = %f %f\n\n", RapMin, RapMax);


			sprintf(outputfilename,"Results/results_pT%1.2f-%1.2f_Rap%1.2f-%1.2f.txt", pTMin, pTMax, RapMin, RapMax);
			printf("output filename is: %s\n", outputfilename);
			outputFile = fopen(outputfilename,"w");

		printf("\nJpsiMass3SigmaMin = %f\n", JpsiMass3SigmaMin);
		printf("\JpsiMass3SigmaMax = %f\n", JpsiMass3SigmaMax);


		printf("\nNsigTot = %f", NsigTot);
		printf(" +- %f \n", errNsigTot);
		printf("NbkgTot = %f", NbkgTot);
		printf(" +- %f \n\n", errNbkgTot);

		printf("fitted mean value is %f", CBmean);
		printf("+- %f\n\n", errCBmean);
		printf("fitted sigma value is %f", resol);
		printf("+- %f\n\n", errresol);

			printf("Nsig in 3 Sigma Range = %f", NsigRange);
			printf(" +- %f \n\n", errNsigTotscaled);

			printf("Nbkg in 3 Sigma Range = %f", NbkgRange);
			printf(" +- %f \n\n", errNbkgTotscaled);
			printf("sig/bkg = %f", sigOvbkg);
			printf(" +- %f\n\n", errsigOvbkg);





			if( outputFile != NULL ){
						fprintf(outputFile, "\nJpsiMass3SigmaMin = %f/n", JpsiMass3SigmaMin);
						fprintf(outputFile, "\JpsiMass3SigmaMax = %f/n", JpsiMass3SigmaMax);
						fprintf(outputFile, "\nNsigTot = %f", NsigTot);
						fprintf(outputFile, " +- %f \n", errNsigTot);
						fprintf(outputFile, "NbkgTot = %f", NbkgTot);
						fprintf(outputFile, " +- %f \n\n", errNbkgTot);
						fprintf(outputFile, "fitted mean value is %f", CBmean);
						fprintf(outputFile, "+- %f\n\n", errCBmean);
						fprintf(outputFile, "fitted sigma value is %f", resol);
						fprintf(outputFile, "+- %f\n\n", errresol);

						fprintf(outputFile, "Nsig in 3 Sigma Range = %f", NsigRange);
						fprintf(outputFile, " +- %f \n\n", errNsigTotscaled);
						fprintf(outputFile, "Nbkg in 3 Sigma Range = %f", NbkgRange);
						fprintf(outputFile, " +- %f \n\n", errNbkgTotscaled);
						fprintf(outputFile, "sig/bkg = %f", sigOvbkg);
						fprintf(outputFile, " +- %f\n\n", errsigOvbkg);
			}



		    fitresMass->Print();
		    fitresctPrompt->Print();
		    if (RealDataIndex == 1){
		        fitresGlobal->Print();
		    }
		    else if(RealDataIndex == 0){

		    if (AnalyticIndex == 1){

		    fitresGlobal->Print();

		    fitresctPromptTemplate->Print();
		    fitresctNonPromptTemplate->Print();
		    fitresctBackgroundTemplate->Print();
		    fitresctNonPrompt->Print();
		    fitresctBackground->Print();
		    }

		    }






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////------------------------------------------------////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////--------c-Tau Print--------//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////-----------------------------------------------/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



		  	NbkgTot= coefctBackground.getVal();
			errNbkgTot = coefctBackground.getError();
			NPromptSig1   = ctcoefSig1.getVal();
			NPromptSig2   = ctcoefSig2.getVal();
			NPromptSig3   = ctcoefSig3.getVal();

			errNPromptSig1   = ctcoefSig1.getError();
			errNPromptSig2   = ctcoefSig2.getError();
			errNPromptSig3   = ctcoefSig3.getError();

			NPromptTot   = NPromptSig1+NPromptSig2+NPromptSig3;
			errNPromptTot = sqrt(pow(errNPromptSig1,2)+pow(errNPromptSig2,2)+pow(errNPromptSig3,2));

			NPromptTottemp   = coefctPrompt.getVal();
			errNPromptTottemp = coefctPrompt.getError();



			NPromptTotfix   = ctcoefPromptfix.getVal();
			errNPromptTotfix   = ctcoefPromptfix.getError();
			NNonPromptTot   = coefctNonPrompt.getVal();
			errNNonPromptTot = coefctNonPrompt.getError();



			NPromptSig1fix   = ctcoefSig1fix.getVal();
			NPromptSig2fix   = ctcoefSig2fix.getVal();
			NPromptSig3fix   = ctcoefSig3fix.getVal();

			errNPromptSig1fix   = ctcoefSig1fix.getError();
			errNPromptSig2fix   = ctcoefSig2fix.getError();
			errNPromptSig3fix   = ctcoefSig3fix.getError();

			NPromptTotfixadd   = NPromptSig1fix+NPromptSig2fix+NPromptSig3fix;
			errNPromptTotfixadd = sqrt(pow(errNPromptSig1,2)+pow(errNPromptSig2,2)+pow(errNPromptSig3,2));

			NNonPrompt2   = coefctNonPrompt2.getVal();
			NNonPrompt3   = coefctNonPrompt3.getVal();
			NNonPromptges   = NNonPrompt2+NNonPrompt3;

			NNonPromptTotanal  = coefctNonPromptanal.getVal();
			errNNonPromptTotanal = coefctNonPromptanal.getError();

			ctresol = sqrt(ctcoefSig1doublenorm*ctsigmaSig1double*ctsigmaSig1double + ctcoefSig2doublenorm*ctsigmaSig2double*ctsigmaSig2double+ctcoefSig3doublenorm*ctsigmaSig3double*ctsigmaSig3double);
			cout<<"ctresol = "<<ctresol<<endl;
			DecayLengthCut = 3*ctresol;

			fprintf(outputFile2, "Prompt ct sigma = %f", ctresol);
			fprintf(outputFile2, "DecayLengthCut = %f", DecayLengthCut);


			LifeTimeCutMin = JpsictMin;
			LifeTimeCutMax = DecayLengthCut;
			Jpsict.setRange("LifeTimeCut",LifeTimeCutMin,LifeTimeCutMax);
			RooAbsReal* integralPromptFit =  ehistPrompt.createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("LifeTimeCut"));
			RooAbsReal* integralNonPromptFit =  ehistNonPrompt.createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("LifeTimeCut"));
			RooAbsReal* integralBkgFit =  eexpo.createIntegral(RooArgSet(Jpsict),NormSet(Jpsict),Range("LifeTimeCut"));
			NPromptRange = integralPromptFit->getVal();
			NNonPromptRange = integralNonPromptFit->getVal();
			NBkgRange = integralBkgFit->getVal();

			// B FRACTION


			NPROMPT=NPromptTotfix;
			ERRNPROMPT=errNPromptTotfix;
			NNONPROMPT=NNonPromptTot;
			ERRNNONPROMPT=errNNonPromptTot;


			printf("\nNbkgTot = %f", NbkgTot);
			printf(" +- %f \n", errNbkgTotscaled);
			printf("NPromptTot = %f", NPROMPT);
			printf(" +- %f \n", ERRNPROMPT);
			printf("NNonPromptTot = %f", NNONPROMPT);
			printf(" +- %f \n", ERRNNONPROMPT);


			printf("\nFraction of Prompt Jpsi's in Range = %f \n", NPromptRange);
			NPromptRange *= NPROMPT;
			errNPromptTotscaled = ERRNPROMPT*NPromptRange/NPROMPT;
			printf("NPromptRange %f", NPromptRange);
			printf(" +- %f\n", errNPromptTotscaled);
			printf("\nFraction B Feed Down in Range = %f \n", NNonPromptRange);
			NNonPromptRange *= NNONPROMPT;
			errNNonPromptTotscaled = ERRNNONPROMPT*NNonPromptRange/NNONPROMPT;
			printf("NNonPromptRange (=B-Feed Down Contamination) %f", NNonPromptRange);
			printf(" +- %f\n\n", errNNonPromptTotscaled);
			NPromptoverNNonPrompt = NPromptRange/NNonPromptRange;
			printf("Prompt over NonPrompt in Range = %f\n\n", NPromptoverNNonPrompt);

			NBkgRange = NBkgRange*NbkgTot;
			SigOvBkgRange = (NPromptRange+NNonPromptRange)/NBkgRange;
			printf("NBkgRange = %f\n\n", NBkgRange);
			printf("SigOvBkgRange = %f\n\n", SigOvBkgRange);

			BFRACTION=NNONPROMPT/(NPROMPT+NNONPROMPT);
			BFRACTIONRANGE=NNonPromptRange/(NPromptRange+NNonPromptRange);
			ERRBFRACTION=(1/(NPROMPT+NNONPROMPT)+NNONPROMPT/pow(NPROMPT+NNONPROMPT,2))*ERRNNONPROMPT+NNONPROMPT/pow(NPROMPT+NNONPROMPT,2)*ERRNPROMPT;
			ERRBFRACTIONRANGE=(1/(NPromptRange+NNonPromptRange)+NNonPromptRange/pow(NPromptRange+NNonPromptRange,2))*errNNonPromptTotscaled+NNonPromptRange/pow(NPromptRange+NNonPromptRange,2)*errNPromptTotscaled;
			printf("B-fraction = %f", BFRACTION);
			printf(" +- %f\n", ERRBFRACTION);
			printf("B-fraction in Range = %f", BFRACTIONRANGE);
			printf(" +- %f\n", ERRBFRACTIONRANGE);


			fprintf(outputFile2, "NbkgTot = %f", NbkgTot);
			fprintf(outputFile2, " +- %f \n", errNbkgTotscaled);
			fprintf(outputFile2, "NPromptTot = %f", NPROMPT);
			fprintf(outputFile2, " +- %f \n", ERRNPROMPT);
			fprintf(outputFile2, "NNonPromptTot = %f", NNONPROMPT);
			fprintf(outputFile2, " +- %f \n", ERRNNONPROMPT);

			fprintf(outputFile2, "\nS/B = %f", sigOvbkg);
			fprintf(outputFile2, " +- %f\n\n", errsigOvbkg);
			fprintf(outputFile2, "B-fraction = %f", BFRACTION);
			fprintf(outputFile2, " +- %f\n", ERRBFRACTION);
			fprintf(outputFile2, "B-fraction in Range = %f", BFRACTIONRANGE);
			fprintf(outputFile2, " +- %f\n", ERRBFRACTIONRANGE);

			fprintf(outputFile2, "\nDecayLengthCut = %f\n", DecayLengthCut);


			if( outputFile != NULL ){

			reddatasum = reddata->sumEntries();
			fprintf(outputFile, "Number of MC Events in Bin = %f\n", reddatasum);

			fprintf(outputFile, "\nNbkgTot = %f", NbkgTot);
			fprintf(outputFile, " +- %f \n", errNbkgTotscaled);
			fprintf(outputFile, "NPromptTot = %f", NPROMPT);
			fprintf(outputFile, " +- %f \n", ERRNPROMPT);
			fprintf(outputFile, "NNonPromptTot = %f", NNONPROMPT);
			fprintf(outputFile, " +- %f \n", ERRNNONPROMPT);
			fprintf(outputFile, "B-fraction = %f", BFRACTION);
			fprintf(outputFile, " +- %f\n", ERRBFRACTION);
			fprintf(outputFile, "B-fraction in Range = %f", BFRACTIONRANGE);
			fprintf(outputFile, " +- %f\n", ERRBFRACTIONRANGE);






			realdata000sum = realdata000->sumEntries();
				realdata00sum = realdata00->sumEntries();
				realdatasum = realdata->sumEntries();
				reddatasum = reddata->sumEntries();
					reddataMCtestsum = reddataMCtest->sumEntries();
					redreddatasum = redreddata->sumEntries();


					hltIndex0PRredsum = hltIndex0PRred->sumEntries();
						hltIndex0NPredsum = hltIndex0NPred->sumEntries();
						hltIndex0BKredsum = hltIndex0BKred->sumEntries();
						hltIndex0PRredcutsum = hltIndex0PRredcut->sumEntries();
						hltIndex0NPredcutsum = hltIndex0NPredcut->sumEntries();
						hltIndex0PRredMCtestsum = hltIndex0PRredMCtest->sumEntries();
						hltIndex0NPredMCtestsum = hltIndex0NPredMCtest->sumEntries();
						hltIndex0BKredMCtestsum = hltIndex0BKredMCtest->sumEntries();
						hltIndex0PRredcutMCtestsum = hltIndex0PRredcutMCtest->sumEntries();
						hltIndex0NPredcutMCtestsum = hltIndex0NPredcutMCtest->sumEntries();
						realdataBkgsum = realdataBkg->sumEntries();
						realdataSBLsum = realdataSBL->sumEntries();
						realdataSBRsum = realdataSBR->sumEntries();

						MCdataBkgsum = MCdataBkg->sumEntries();
							MCdataSBLsum = MCdataSBL->sumEntries();
							MCdataSBRsum = MCdataSBR->sumEntries();
							redrealdatasum = redrealdata->sumEntries();
							hltIndex0PRsum = hltIndex0PR->sumEntries();
							hltIndex0NPsum = hltIndex0NP->sumEntries();
							hltIndex0BKsum = hltIndex0BK->sumEntries();
							hltIndex0PRMCtestsum = hltIndex0PRMCtest->sumEntries();
							hltIndex0NPMCtestsum = hltIndex0NPMCtest->sumEntries();
							hltIndex0BKMCtestsum = hltIndex0BKMCtest->sumEntries();
							hltIndex0PRredMCtestsum = hltIndex0PRredMCtest->sumEntries();
							hltIndex0NPredMCtestsum = hltIndex0NPredMCtest->sumEntries();
							hltIndex0BKredMCtestsum = hltIndex0BKredMCtest->sumEntries();
							hltIndex0PRredcutMCtestsum = hltIndex0PRredcutMCtest->sumEntries();
							hltIndex0NPredcutMCtestsum = hltIndex0NPredcutMCtest->sumEntries();



							sigovbkgMCTRUTH = (hltIndex0PRredMCtestsum+hltIndex0NPredMCtestsum)/hltIndex0BKredMCtestsum;
							bfractionMCTRUTH = hltIndex0NPredMCtestsum/(hltIndex0PRredMCtestsum+hltIndex0NPredMCtestsum);
							bfractionMCTRUTHrange = hltIndex0NPredcutMCtestsum/(hltIndex0PRredcutMCtestsum+hltIndex0NPredcutMCtestsum);


									fprintf(outputFile, "realdata000sum = %f\n", realdata000sum);
									fprintf(outputFile, "realdata00sum = %f\n", realdata00sum);
									fprintf(outputFile, "realdatasum = %f\n", realdatasum);
									fprintf(outputFile, "reddatasum = %f\n", reddatasum);


							fprintf(outputFile2, "S/B MCTRUTH = %f\n", sigovbkgMCTRUTH);
							fprintf(outputFile2, "B-fraction MCTRUTH = %f\n", bfractionMCTRUTH);
							fprintf(outputFile2, "B-fraction in Range MCTRUTH = %f\n", bfractionMCTRUTHrange);


							fprintf(outputFile2, "\nNumber of MC events in bin = %1.2f resp. (test) %1.2f\n",reddatasum,reddataMCtestsum);
							fprintf(outputFile2, "Number of real data events in bin = %1.2f\n\n",realdatasum);

								fprintf(outputFile, "reddataMCtestsum = %f\n", reddataMCtestsum);
								fprintf(outputFile, "hltIndex0PRsum = %f\n", hltIndex0PRsum);
								fprintf(outputFile, "hltIndex0NPsum = %f\n", hltIndex0NPsum);
								fprintf(outputFile, "hltIndex0BKsum = %f\n", hltIndex0BKsum);
								fprintf(outputFile, "hltIndex0PRMCtestsum = %f\n", hltIndex0PRMCtestsum);
								fprintf(outputFile, "hltIndex0NPMCtestsum = %f\n", hltIndex0NPMCtestsum);
								fprintf(outputFile, "hltIndex0BKMCtestsum = %f\n", hltIndex0BKMCtestsum);
								fprintf(outputFile, "redreddatasum = %f\n", redreddatasum);
								fprintf(outputFile, "hltIndex0PRredsum = %f\n", hltIndex0PRredsum);
								fprintf(outputFile, "hltIndex0NPredsum = %f\n", hltIndex0NPredsum);
								fprintf(outputFile, "hltIndex0BKredsum = %f\n", hltIndex0BKredsum);
								fprintf(outputFile, "hltIndex0PRredcutsum = %f\n", hltIndex0PRredcutsum);
								fprintf(outputFile, "hltIndex0NPredcutsum = %f\n", hltIndex0NPredcutsum);
								fprintf(outputFile, "hltIndex0PRredMCtestsum = %f\n", hltIndex0PRredMCtestsum);
								fprintf(outputFile, "hltIndex0NPredMCtestsum = %f\n", hltIndex0NPredMCtestsum);
								fprintf(outputFile, "hltIndex0BKredMCtestsum = %f\n", hltIndex0BKredMCtestsum);
								fprintf(outputFile, "hltIndex0PRredcutMCtestsum = %f\n", hltIndex0PRredcutMCtestsum);
								fprintf(outputFile, "hltIndex0NPredcutMCtestsum = %f\n", hltIndex0NPredcutMCtestsum);
								fprintf(outputFile, "realdataBkgsum = %f\n", realdataBkgsum);
								fprintf(outputFile, "realdataSBLsum = %f\n", realdataSBLsum);
								fprintf(outputFile, "realdataSBRsum = %f\n", realdataSBRsum);
								fprintf(outputFile, "MCdataBkgsum = %f\n", MCdataBkgsum);
								fprintf(outputFile, "MCdataSBLsum = %f\n", MCdataSBLsum);
								fprintf(outputFile, "MCdataSBRsum = %f\n", MCdataSBRsum);
								fprintf(outputFile, "redrealdatasum = %f\n", redrealdatasum);

			}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitPol(){



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////------------------------------------------------////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////--------POLARIZATION CUTS--------////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////-----------------------------------------------/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			AngularPRhistMCtest_CS = new RooDataHist("AngularPRhistMCtest_CS","",RooArgSet(costh_CS,phi_CS),*hltIndex0PRMCtest);
			AngularNPhistMCtest_CS = new RooDataHist("AngularNPhistMCtest_CS","",RooArgSet(costh_CS,phi_CS),*hltIndex0NPMCtest);
			AngularBKhistMCtest_CS = new RooDataHist("AngularBKhistMCtest_CS","",RooArgSet(costh_CS,phi_CS),*hltIndex0BKMCtest);
			AngularPRhistMCtest_HX = new RooDataHist("AngularPRhistMCtest_HX","",RooArgSet(costh_HX,phi_HX),*hltIndex0PRMCtest);
			AngularNPhistMCtest_HX = new RooDataHist("AngularNPhistMCtest_HX","",RooArgSet(costh_HX,phi_HX),*hltIndex0NPMCtest);
			AngularBKhistMCtest_HX = new RooDataHist("AngularBKhistMCtest_HX","",RooArgSet(costh_HX,phi_HX),*hltIndex0BKMCtest);

//			RooHistPdf AngularHistPdf_CS_Bkg("AngularHistPdf_CS_Bkg","AngularHistPdf_CS_Bkg",RooArgSet(costh_CS,phi_CS),*AngularBKhistMCtest_CS,1);

//			RooHistPdf AngularHistPdf_HX_Bkg("AngularHistPdf_HX_Bkg","AngularHistPdf_HX_Bkg",RooArgSet(costh_HX,phi_HX),*AngularBKhistMCtest_HX,1);



			char reduceDecayLengthCut[200];

			if (PromptIndex == 1){
			sprintf(reduceDecayLengthCut,"Jpsict < %f",DecayLengthCut);
			}
			else if (PromptIndex == 0){
			sprintf(reduceDecayLengthCut,"Jpsict > %f",DecayLengthCut);
			}

			AngularMCtest = (RooDataSet*)reddataMCtest->reduce(reduceDecayLengthCut);
			AngularMCtest_hist_CS = new RooDataHist("AngularMCtest_hist_CS","",RooArgSet(costh_CS,phi_CS),*AngularMCtest);
			AngularMCtest_hist_HX = new RooDataHist("AngularMCtest_hist_HX","",RooArgSet(costh_HX,phi_HX),*AngularMCtest);


			AngularRealData = (RooDataSet*)redrealdata->reduce(reduceDecayLengthCut);


						RooBinning rbcosth(-1,1);
						rbcosth.addUniform(Angular_Binning_costh,-1,1);
						costh_CS.setBinning(rbcosth);

						RooBinning rbphi(0,360);
						rbphi.addUniform(Angular_Binning_phi,0,360);
						phi_CS.setBinning(rbphi);

						costh_HX.setBinning(rbcosth);
						phi_HX.setBinning(rbphi);



			AngularRealData_hist_CS = new RooDataHist("AngularRealData_hist_CS","AngularRealData_hist_CS",RooArgSet(costh_CS,phi_CS),*AngularRealData);
			AngularRealData_hist_HX = new RooDataHist("AngularRealData_hist_HX","AngularRealData_hist_HX",RooArgSet(costh_HX,phi_HX),*AngularRealData);

//			AngularRealData_hist_CS_TH1 = AngularRealData_hist_CS->createHistogram("AngularRealData_hist_CS",costh_CS,Binning(Angular_Binning_costh),YVar(phi_CS,Binning(Angular_Binning_phi))) ;
//			AngularRealData_hist_HX_TH1 = AngularRealData_hist_HX->createHistogram("AngularRealData_hist_HX",costh_HX,Binning(Angular_Binning_costh),YVar(phi_HX,Binning(Angular_Binning_phi))) ;


			AngularRealData_hist_CS->Print();
			AngularRealData_hist_HX->Print();

			RooBinning rbcosth(-1,1);
			rbcosth.addUniform(100,-1,1);
			costh_CS.setBinning(rbcosth);

			RooBinning rbphi(0,360);
			rbphi.addUniform(100,0,360);
			phi_CS.setBinning(rbphi);

			costh_HX.setBinning(rbcosth);
			phi_HX.setBinning(rbphi);


			RooDataHist *AcceptanceHist_CS = new RooDataHist("AcceptanceHist_CS","AcceptanceHist_CS",RooArgSet(costh_CS,phi_CS),hAcc2D_pol_pT_rap[CS][npT][nRap]);

			RooHistPdf AcceptanceHistPdf_CS("AcceptanceHistPdf_CS","AcceptanceHistPdf_CS",RooArgSet(costh_CS,phi_CS),*AcceptanceHist_CS,1);

			RooDataHist *AcceptanceHist_HX = new RooDataHist("AcceptanceHist_HX","AcceptanceHist_HX",RooArgSet(costh_HX,phi_HX),hAcc2D_pol_pT_rap[HX][npT][nRap]);

			RooHistPdf AcceptanceHistPdf_HX("AcceptanceHistPdf_HX","AcceptanceHistPdf_HX",RooArgSet(costh_HX,phi_HX),*AcceptanceHist_HX,1);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////------------------------------------------------////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////--------POLARIZATION FITS--------////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////-----------------------------------------------/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////--------Collins soper--------////////////////////////////////////////////////////////////////////////////////////////////////////////





			RooRealVar lambda_theta_CS("lambda_theta_CS","lambda_theta_CS",0,-2,2);
			RooRealVar lambda_phi_CS("lambda_phi_CS","lambda_phi_CS",0,-2,2);
			RooRealVar lambda_thetaphi_CS("lambda_thetaphi_CS","lambda_thetaphi_CS",0);


			RooRealVar SigOvBkgRangeVar("SigOvBkgRangeVar","SigOvBkgRangeVar",SigOvBkgRange);

			RooGenericPdf polfunc_CS("polfunc_CS","polfunc_CS","3/(4*pi*(3+lambda_theta_CS))*(1+lambda_theta_CS*pow(costh_CS,2)+lambda_phi_CS*(1-pow(costh_CS,2))*cos(2*(phi_CS/180*pi))+lambda_thetaphi_CS*sin(2*acos(costh_CS))*cos(phi_CS/180*pi))",RooArgSet(lambda_theta_CS,lambda_phi_CS,lambda_thetaphi_CS,costh_CS,phi_CS)) ;

//			RooGenericPdf polfunc_CS("polfunc_CS","polfunc_CS","1/(3+lambda_theta_CS))*(1+lambda_theta_CS*pow(costh_CS,2))",RooArgSet(lambda_theta_CS,costh_CS)) ;



//			RooDataHist *AcceptanceHist_CS = new RooDataHist("AcceptanceHist_CS","AcceptanceHist_CS",RooArgSet(costh_CS,phi_CS),hAcc2D_pol_pT_rap[CS][npT][nRap]);

//			RooHistPdf AcceptanceHistPdf_CS("AcceptanceHistPdf_CS","AcceptanceHistPdf_CS",RooArgSet(costh_CS,phi_CS),*AcceptanceHist_CS,3);

			RooProdPdf CSPDF("CSPDF","CSPDF",RooArgList(polfunc_CS,AcceptanceHistPdf_CS));
		//	RooAddPdf CSPDF("CSPDF","CSPDF",RooArgList(AngularHistPdf_CS_Bkg,CSPDF_Prod),SigOvBkgRangeVar);



			phi_CS.setRange("phi_CSFitRange",60,115);

			if (AngularIndex == 1){

			if (RealDataIndex == 0){

			fitresthetaphi_CS = CSPDF.fitTo(*AngularMCtest_hist_CS,Save(1),Minos(0),SumW2Error(kTRUE));

			}


			else if (RealDataIndex == 1){

			fitresthetaphi_CS = CSPDF.fitTo(*AngularRealData_hist_CS,Save(1),Minos(0),SumW2Error(kTRUE));

			}
			}

			char titlepolfunchisto_CS[200];
			sprintf(titlepolfunchisto_CS,"titlepolfunchisto_CS, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);



			sprintf(titlehFitAccCorr2D_CS, "hFitAccCorr2D_CS_pT%d_rap%d", npT, nRap);
			sprintf(titlehFit2D_CS, "hFit2D_CS_pT%d_rap%d",  npT, nRap);
			sprintf(titlehAcc2D_CS, "hAcc2D_CS_pT%d_rap%d",  npT, nRap);
			sprintf(titlehData2D_CS, "hData2D_CS_pT%d_rap%d",  npT, nRap);

			acceptancehist_CS= AcceptanceHistPdf_CS.createHistogram(titlehAcc2D_CS,costh_CS,Binning(Angular_Binning_costh),YVar(phi_CS,Binning(Angular_Binning_phi))) ;
			polfunchist_CS= polfunc_CS.createHistogram(titlehFit2D_CS,costh_CS,Binning(Angular_Binning_costh),YVar(phi_CS,Binning(Angular_Binning_phi))) ;

			polfunchisto_CS = CSPDF.createHistogram(titlehFitAccCorr2D_CS,costh_CS,Binning(Angular_Binning_costh),YVar(phi_CS,Binning(Angular_Binning_phi))) ;
			if (RealDataIndex == 0){
			datahistoangular_CS = redreddataMCtest->createHistogram("datahistoangular_CS",costh_CS,Binning(Angular_Binning_costh),YVar(phi_CS,Binning(Angular_Binning_phi))) ;
			}
			else if (RealDataIndex == 1){
			datahistoangular_CS = AngularRealData->createHistogram(titlehData2D_CS,costh_CS,Binning(Angular_Binning_costh),YVar(phi_CS,Binning(Angular_Binning_phi))) ;
			}


			phi_CSframe = phi_CS.frame(Angular_Binning_phi) ;
			if (RealDataIndex == 0){
			AngularMCtest_hist_CS->plotOn(phi_CSframe,DataError(RooAbsData::SumW2));
			}
			else if (RealDataIndex == 1){
			AngularRealData->plotOn(phi_CSframe,DataError(RooAbsData::SumW2));
			}
			CSPDF.plotOn(phi_CSframe,Normalization(1.0));
			CSPDF.paramOn(phi_CSframe,Layout(0.25,0.95,0.85),Format("NEU",FixedPrecision(3)),Parameters(RooArgSet(lambda_phi_CS,lambda_theta_CS))) ;
//			CSPDF.plotOn(phi_CSframe, RooFit::Components(RooArgList(polfunc_CS)), RooFit::LineColor(kRed),RooFit::LineStyle(kDashed));
//			CSPDF.plotOn(phi_CSframe, RooFit::Components(RooArgList(AcceptanceHistPdf_CS)), RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));
//			polfunc_CS.plotOn(phi_CSframe,Normalization(1.0));


			costh_CSframe = costh_CS.frame(Angular_Binning_costh) ;
			if (RealDataIndex == 0){
			AngularMCtest_hist_CS->plotOn(costh_CSframe,DataError(RooAbsData::SumW2));
			}
			else if (RealDataIndex == 1){
			AngularRealData->plotOn(costh_CSframe,DataError(RooAbsData::SumW2));
			}
			CSPDF.plotOn(costh_CSframe,Normalization(1.0));
			CSPDF.paramOn(costh_CSframe,Layout(0.25,0.95,0.85),Format("NEU",FixedPrecision(3)),Parameters(RooArgSet(lambda_phi_CS,lambda_theta_CS))) ;

//			acceptancehist_CS->SaveAs("RootOutput/TH2F_AccMap_CS.root");
//			polfunchist_CS->SaveAs("RootOutput/TH2F_Fit_CS.root");
//			polfunchisto_CS->SaveAs("RootOutput/TH2F_FitAccCorr_CS.root");
//			datahistoangular_CS->SaveAs("RootOutput/TH2F_Data_CS.root");


			fitresthetaphi_CS->Print();



//////////////////////--------Helicity--------////////////////////////////////////////////////////////////////////////////////////////////////////////


			RooRealVar lambda_theta_HX("lambda_theta_HX","lambda_theta_HX",0,-2,2);
			RooRealVar lambda_phi_HX("lambda_phi_HX","lambda_phi_HX",0,-2,2);
			RooRealVar lambda_thetaphi_HX("lambda_thetaphi_HX","lambda_thetaphi_HX",0);


			RooGenericPdf polfunc_HX("polfunc_HX","polfunc_HX","3/(4*pi*(3+lambda_theta_HX))*(1+lambda_theta_HX*pow(costh_HX,2)+lambda_phi_HX*(1-pow(costh_HX,2))*cos(2*(phi_HX/180*pi))+lambda_thetaphi_HX*sin(2*acos(costh_HX))*cos(phi_HX/180*pi))",RooArgSet(lambda_theta_HX,lambda_phi_HX,lambda_thetaphi_HX,costh_HX,phi_HX)) ;

//			RooDataHist *AcceptanceHist_HX = new RooDataHist("AcceptanceHist_HX","AcceptanceHist_HX",RooArgSet(costh_HX,phi_HX),hAcc2D_pol_pT_rap[HX][npT][nRap]);

//			RooHistPdf AcceptanceHistPdf_HX("AcceptanceHistPdf_HX","AcceptanceHistPdf_HX",RooArgSet(costh_HX,phi_HX),*AcceptanceHist_HX,3);
//			RooGenericPdf polfunc_HX_sign("polfunc_HX_sign","polfunc_HX_sign","sign(3/(4*pi*(3+lambda_theta_HX))*(1+lambda_theta_HX*pow(costh_HX,2)+lambda_phi_HX*(1-pow(costh_HX,2))*cos(2*(phi_HX/180*pi))+lambda_thetaphi_HX*sin(2*acos(costh_HX))*cos(phi_HX/180*pi)))+1",RooArgSet(lambda_theta_HX,lambda_phi_HX,lambda_thetaphi_HX,costh_HX,phi_HX)) ;

			RooProdPdf HXPDF("HXPDF","HXPDF",RooArgList(AcceptanceHistPdf_HX,polfunc_HX));
//			RooAddPdf HXPDF("HXPDF","HXPDF",RooArgList(AngularHistPdf_HX_Bkg,HXPDF_Prod),SigOvBkgRangeVar);



			costh_HX.setRange("costh_HXFitRange",-0.5,0.5);

			if (AngularIndex == 1){

			if (RealDataIndex == 0){

			fitresthetaphi_HX = HXPDF.fitTo(*AngularMCtest_hist_HX,Save(1),Minos(0),SumW2Error(kTRUE));

			}


			else if (RealDataIndex == 1){

			fitresthetaphi_HX = HXPDF.fitTo(*AngularRealData_hist_HX,Save(1),Minos(0),SumW2Error(kTRUE));//(*AngularRealData_hist_HX,Save(1),Minos(0),SumW2Error(kFALSE));

			}
			}

			char titlepolfunchisto_HX[200];
			sprintf(titlepolfunchisto_HX,"titlepolfunchisto_HX, pT = %1.2fGeV-%1.2fGeV, Rap = %1.2f-%1.2f", pTMin, pTMax, RapMin, RapMax);


			sprintf(titlehFitAccCorr2D_HX, "hFitAccCorr2D_HX_pT%d_rap%d", npT, nRap);
			sprintf(titlehFit2D_HX, "hFit2D_HX_pT%d_rap%d", npT, nRap);
			sprintf(titlehAcc2D_HX, "hAcc2D_HX_pT%d_rap%d", npT, nRap);
			sprintf(titlehData2D_HX, "hData2D_HX_pT%d_rap%d", npT, nRap);

			acceptancehist_HX= AcceptanceHistPdf_HX.createHistogram(titlehAcc2D_HX,costh_HX,Binning(Angular_Binning_costh),YVar(phi_HX,Binning(Angular_Binning_phi))) ;
			polfunchist_HX= polfunc_HX.createHistogram(titlehFit2D_HX,costh_HX,Binning(Angular_Binning_costh),YVar(phi_HX,Binning(Angular_Binning_phi))) ;

			polfunchisto_HX = HXPDF.createHistogram(titlehFitAccCorr2D_HX,costh_HX,Binning(Angular_Binning_costh),YVar(phi_HX,Binning(Angular_Binning_phi))) ;
			if (RealDataIndex == 0){
			datahistoangular_HX = redreddataMCtest->createHistogram("datahistoangular_HX",costh_HX,Binning(Angular_Binning_costh),YVar(phi_HX,Binning(Angular_Binning_phi))) ;
			}
			else if (RealDataIndex == 1){
			datahistoangular_HX = AngularRealData->createHistogram(titlehData2D_HX,costh_HX,Binning(Angular_Binning_costh),YVar(phi_HX,Binning(Angular_Binning_phi))) ;
			}

			phi_HXframe = phi_HX.frame(Angular_Binning_phi) ;
			if (RealDataIndex == 0){
			AngularMCtest_hist_HX->plotOn(phi_HXframe,DataError(RooAbsData::SumW2));
			}
			else if (RealDataIndex == 1){
			AngularRealData->plotOn(phi_HXframe,DataError(RooAbsData::SumW2));
			}
			HXPDF.plotOn(phi_HXframe,Normalization(1.0));
			HXPDF.paramOn(phi_HXframe,Layout(0.25,0.95,0.85),Format("NEU",FixedPrecision(3)),Parameters(RooArgSet(lambda_phi_HX,lambda_theta_HX))) ;

			costh_HXframe = costh_HX.frame(Angular_Binning_costh) ;
			if (RealDataIndex == 0){
			AngularMCtest_hist_HX->plotOn(costh_HXframe,DataError(RooAbsData::SumW2));
			}
			else if (RealDataIndex == 1){
			AngularRealData->plotOn(costh_HXframe,DataError(RooAbsData::SumW2));
			}
			HXPDF.plotOn(costh_HXframe,Normalization(1.0));
			HXPDF.paramOn(costh_HXframe,Layout(0.25,0.95,0.85),Format("NEU",FixedPrecision(3)),Parameters(RooArgSet(lambda_phi_HX,lambda_theta_HX))) ;

//			acceptancehist_HX->SaveAs("RootOutput/TH2F_AccMap_HX.root");
//			polfunchist_HX->SaveAs("RootOutput/TH2F_Fit_HX.root");
//			polfunchisto_HX->SaveAs("RootOutput/TH2F_FitAccCorr_HX.root");
//			datahistoangular_HX->SaveAs("RootOutput/TH2F_Data_HX.root");

			char HermineFileName[200];
			sprintf(HermineFileName,"RootOutput/PolarizationFitOutput_pT%d_rap%d.root", npT, nRap);

			TFile fOut(HermineFileName, "RECREATE");
			fOut.cd();

			acceptancehist_CS->Write();
			polfunchist_CS->Write();
			polfunchisto_CS->Write();
			datahistoangular_CS->Write();
			acceptancehist_HX->Write();
			polfunchist_HX->Write();
			polfunchisto_HX->Write();
			datahistoangular_HX->Write();

			fOut.Close();


			fitresthetaphi_HX->Print();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Chi2(){


	cout << "chi2/ndof (Mass) = " << JpsiMassframeChi2->chiSquare() << endl;
	chi2Mass=JpsiMassframeChi2->chiSquare();
  	if(RealDataIndex == 0){

	cout << "chi2/ndof (ctGlobalMC) = " << JpsictframeChi2->chiSquare() << endl;
	chi2ctGlobal=JpsictframeChi2->chiSquare();

	printf("nbin = %1.0f, npT = %1.0f, nRap = %1.0f\n", nbin, npT, nRap);







    if (AnalyticIndex == 1){

	cout << "chi2/ndof (Promptct) = " << JpsictframePrompt->chiSquare() << endl;
	chi2ctPrompt=JpsictframePrompt->chiSquare();
	cout << "chi2/ndof (NonPromptct) = " << JpsictframeNonPrompt->chiSquare() << endl;
	chi2ctNonPrompt=JpsictframeNonPrompt->chiSquare();
	cout << "chi2/ndof (Backgroundct) = " << JpsictframeBackground->chiSquare() << endl;
	chi2ctBackground=JpsictframeBackground->chiSquare();
	cout << "chi2/ndof (PromptctAnal) = " << JpsictframePromptAnalChi2->chiSquare() << endl;
	chi2ctPromptAnal=JpsictframePromptAnalChi2->chiSquare();

    }
  	}
 	else if(RealDataIndex == 1){
    cout << "chi2/ndof (ctGlobalData) = " << JpsictframeDataChi2->chiSquare() << endl;
	chi2ctGlobalData=JpsictframeDataChi2->chiSquare();
 	}

	fprintf(outputFile, "chi2-Mass = %f\n", chi2Mass);

	chi2_phi_HX=phi_HXframe->chiSquare();
	chi2_costh_HX=costh_HXframe->chiSquare();
	chi2_phi_CS=phi_CSframe->chiSquare();
	chi2_costh_CS=costh_CSframe->chiSquare();

	printf("\nchi2 costh_CS = %1.2f\n", chi2_costh_CS);
	printf("\nchi2 phi_CS = %1.2f\n", chi2_phi_CS);
	printf("\nchi2 costh_HX = %1.2f\n", chi2_costh_HX);
	printf("\nchi2 phi_HX = %1.2f\n", chi2_phi_HX);

	fprintf(outputFile2, "\nchi2 costh_CS = %1.2f\n", chi2_costh_CS);
	fprintf(outputFile2, "\nchi2 phi_CS = %1.2f\n", chi2_phi_CS);
	fprintf(outputFile2, "\nchi2 costh_HX = %1.2f\n", chi2_costh_HX);
	fprintf(outputFile2, "\nchi2 phi_HX = %1.2f\n", chi2_phi_HX);




	if (RealDataIndex == 0){


	fprintf(outputFile, "chi2-Global = %f\n", chi2ctGlobal);

	if (AnalyticIndex == 1){

		fprintf(outputFile, "chi2-Prompt = %f\n", chi2ctPrompt);
		fprintf(outputFile, "chi2-PromptAnal = %f\n", chi2ctPromptAnal);
		fprintf(outputFile, "chi2-NonPrompt = %f\n", chi2ctNonPrompt);
		fprintf(outputFile, "chi2-Background = %f\n", chi2ctBackground);

		}
		}
	else if (RealDataIndex == 1){

	fprintf(outputFile, "chi2-Global Data = %f\n", chi2ctGlobalData);
	fprintf(outputFile, "Numbers of events in sideband left = %f\n", NSBL);
	fprintf(outputFile, "Numbers of events in sideband right = %f\n", NSBR);

	}
	fclose(outputFile);


}
