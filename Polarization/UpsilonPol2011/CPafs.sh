#!/bin/sh

#fromdir=/Users/valentinknuenz/usr/local/workspace/UpsilonPol
fromdir=/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol
todir=/afs/hephy.at/scratch/k/knuenz/UserCode/Knuenz/Polarization/UpsilonPol2011
#todir=/afs/hephy.at/user/k/knuenz/public/forLuca

mkdir ${todir}
mkdir ${todir}/macros
mkdir ${todir}/macros/polFit
mkdir ${todir}/macros/polFit/EffFiles
mkdir ${todir}/macros/polFit/MattRes
mkdir ${todir}/macros/polFit/Smearing
mkdir ${todir}/macros/polFit/Systematics
mkdir ${todir}/latex
mkdir ${todir}/interface

#cp ${fromdir}/macros/polFit/EffFiles/* ${todir}/macros/polFit/EffFiles/
#cp ${fromdir}/macros/polFit/MattRes/* ${todir}/macros/polFit/MattRes/
#cp -r ${fromdir}/macros/polFit/Systematics ${todir}/macros/polFit/

cp ${fromdir}/macros/polFit/runToyMC.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/runDataFits.sh ${todir}/macros/polFit/

cp ${fromdir}/macros/polFit/PlotDataFits.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotToyMC.sh ${todir}/macros/polFit/

cp ${fromdir}/macros/polFit/EvalSyst.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotResults.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotCentrals.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/AverageSyst.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/ChangeTGraph.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/RePlotDataFits.sh ${todir}/macros/polFit/

cp ${fromdir}/macros/polFit/AlterPPD.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/AlterPPD.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/FitBGratio.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/FitRho.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/FitTGraph.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/Chi2FitBGratio.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/Chi2FitBGSBratio.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/Make_fLSB_varPlot.C ${todir}/macros/polFit/

cp ${fromdir}/macros/polFit/cpDataFiguresToAFS.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/cpIndividuals.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/cpLessFiguresToAFS.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/cpToyMCFiguresToAFS.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/MergeDataFiles.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotCombinedToyMC.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotSingleDataFits.sh ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/RePlotDataFits.sh ${todir}/macros/polFit/


cp ${fromdir}/macros/polFit/polGenRecFitPlot.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/polRapPtPlot.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotFinalResults.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/EvaluateSyst.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/AverageSystematics.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/MattsResults.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/Makefile ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/ScaleMCtruth.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/ScaleChi2Func.h ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/ChangeTGraph.cc ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/PlotCombinedToyMC.sh ${todir}/macros/polFit/

cp ${fromdir}/macros/polFit/polFit.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/polGen.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/polPlot.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/polRec.C ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/effsAndCuts.h ${todir}/macros/polFit/
cp ${fromdir}/macros/polFit/ToyMC.h ${todir}/macros/polFit/

cp ${fromdir}/latex/ToyResults.tex ${todir}/latex/
cp ${fromdir}/latex/PullSummaryResults.tex ${todir}/latex/
cp ${fromdir}/latex/ParameterSummaryResults.tex ${todir}/latex/
cp ${fromdir}/latex/IndividualFitResults.tex ${todir}/latex/
cp ${fromdir}/latex/DataResults_vs_RapPt.tex ${todir}/latex/
cp ${fromdir}/latex/FinalDataResults.tex ${todir}/latex/
cp ${fromdir}/latex/Systematics.tex ${todir}/latex/

cp ${fromdir}/latex/README.tex ${todir}/latex
cp ${fromdir}/latex/README.pdf ${todir}/latex

cp ${fromdir}/latex/cosThetaPhi_BG.tex ${todir}/latex/
cp ${fromdir}/latex/cosThetaPhi.tex ${todir}/latex/
cp ${fromdir}/latex/massFits.tex ${todir}/latex/
cp ${fromdir}/latex/ScaledMCEff.tex ${todir}/latex/

cp ${fromdir}/interface/commonVar.h ${todir}/interface/
cp ${fromdir}/interface/rootIncludes.inc ${todir}/interface/

cp ${fromdir}/macros/calcMeanPt.C ${todir}/macros/
cp ${fromdir}/macros/calcPol.C ${todir}/macros/
cp ${fromdir}/macros/CBFunction.C ${todir}/macros/
cp ${fromdir}/macros/CopyTreeEntries.C ${todir}/macros/
cp ${fromdir}/macros/PlotCosThetaPhiBG.cc ${todir}/macros/
cp ${fromdir}/macros/PlotCosThetaPhiDistribution.cc ${todir}/macros/
cp ${fromdir}/macros/PolData.C ${todir}/macros/
cp ${fromdir}/macros/PolData.h ${todir}/macros/
cp ${fromdir}/macros/runCopyTreeEntries.cc ${todir}/macros/
cp ${fromdir}/macros/runData.cc ${todir}/macros/
cp ${fromdir}/macros/runMassFit.cc ${todir}/macros/
cp ${fromdir}/macros/runMeanPt.cc ${todir}/macros/
cp ${fromdir}/macros/runTrimEventContent.cc ${todir}/macros/
cp ${fromdir}/macros/TrimEventContent.C ${todir}/macros/
cp ${fromdir}/macros/upsilon_2StepFit.C ${todir}/macros/
cp ${fromdir}/macros/upsilon_MCMassFit.C ${todir}/macros/
cp ${fromdir}/macros/Makefile ${todir}/macros/
cp ${fromdir}/macros/PrepareTTree.sh ${todir}/macros/

cp ${fromdir}/addToysCVS.sh ${todir}/
cp ${fromdir}/CPafs.sh ${todir}/

