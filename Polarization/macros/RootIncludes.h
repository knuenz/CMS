#include "commonVar.h"
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

gStyle->SetPalette(1);
gSystem->Load("libRooFit") ;
gSystem->Load("libRooFitCore") ;


using namespace RooFit ;
