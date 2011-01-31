/*
 * GenTest.C
 *
 *  Created on: 02.11.2010
 *      Author: valentinknuenz
 */
'.include /afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.25.02-cms6/include'

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooConstVar.h"
#include "RooCmdArg.h"
#include "RooGExpModel.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TRandom.h"
#include "RooFactoryWSTool.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooAbsData.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "TH1.h"
#include "RooBinning.h"
#include "RooUniformBinning.h"
#include "RooGenProdProj.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooNumConvPdf.h"
#include "RooHistPdf.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "TLegend.h"
#include "TMath.h"


gSystem->Load("libRooFit") ;
gSystem->Load("libRooFitCore") ;


using namespace RooFit ;



Char_t *fileName2 = "/scratch/knuenz/Polarization/RootInput/TTree_pol_Mu0Track0Jpsi_MCprompt.root";

RooRealVar costh_CS_Gen("costh_CS_Gen","cos #theta_{CS}",-1,1);
RooRealVar phi_CS_Gen("phi_CS_Gen","#phi_{CS} [deg]",0,360);
RooRealVar costh_HX_Gen("costh_HX_Gen","cos#theta_{HX}",-1,1);
RooRealVar phi_HX_Gen("phi_HX_Gen","#phi_{HX} [deg]",0,360);


void GenTest(){

	fIn2 = new TFile(fileName2);
	RooDataSet* data = (RooDataSet*)fIn2->Get("data");
	data.Print();



}
