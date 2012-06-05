#!/bin/sh

todir=/afs/hephy.at/scratch/k/knuenz/UserCode/Knuenz/Polarization/UpsilonPol
comment=Developments


cd ${todir}/macros/polFit

cvs add runDataFits.sh
cvs add runToyMC.sh
cvs add PlotDataFits.sh
cvs add PlotToyMC.sh
cvs add polGenRecFitPlot.cc
cvs add polRapPtPlot.cc
cvs add Makefile
cvs add polGen.C
cvs add polRec.C
cvs add polFit.C
cvs add polPlot.C
cvs add ToyMC.h
cvs add effsAndCuts.h
cvs add EvalSyst.sh
cvs add PlotResults.sh
cvs add AverageSyst.sh
cvs add PlotFinalResults.cc
cvs add EvaluateSyst.cc
cvs add AverageSystematics.cc
cvs add ScaleMCtruth.C
cvs add ScaleChi2Func.h
cvs add PlotCombinedToyMC.sh
cvs add ChangeTGraph.cc
cvs add ChangeTGraph.sh
cvs add RePlotDataFits.sh
cvs add PlotCentrals.sh
cvs add MattsResults.cc


cvs ci -m "${comment}" 

cd ${todir}/latex

cvs add ToyResults.tex
cvs add PullSummaryResults.tex
cvs add ParameterSummaryResults.tex
cvs add IndividualFitResults.tex
cvs add DataResults_vs_RapPt.tex
cvs add cosThetaPhi_BG.tex
cvs add cosThetaPhi.tex
cvs add massFits.tex
cvs add README.tex
cvs add README.pdf
cp ${fromdir}/latex/Systematics.tex ${todir}/latex/
cp ${fromdir}/latex/FinalDataResults.tex ${todir}/latex/

cvs ci -m "${comment}" 


cd ${todir}/interface

cvs add commonVar.h
cvs add rootIncludes.inc

cvs ci -m "${comment}" 


cd ${todir}/macros/polFit/Smearing

cvs add ScaledMCEff.tex

cvs ci -m "${comment}" 


cd ${todir}/macros

cvs add calcMeanPt.C
cvs add calcPol.C
cvs add CBFunction.C
cvs add CopyTreeEntries.C
cvs add PlotCosThetaPhiBG.cc
cvs add PlotCosThetaPhiDistribution.cc
cvs add PolData.C
cvs add PolData.h
cvs add runCopyTreeEntries.cc
cvs add runData.cc
cvs add runMassFit.cc
cvs add runMeanPt.cc
cvs add runTrimEventContent.cc
cvs add TrimEventContent.C
cvs add upsilon_2StepFit.C
cvs add Makefile
cvs add PrepareTTree.sh
cvs add upsilon_MCMassFit.C
cvs add PrepareTTree.sh

cvs ci -m "${comment}" 

cd ${todir}

cvs add addToysCVS.sh
cvs add CPafs.sh

cvs ci -m "${comment}" 

pwd