#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir

########## INPUTS ##########

for nState in 1 3;do

JobID=July4
additionalName=_ICHEPnew_Ups${nState}S #SystBkgNew for centrals

PlotMatt=1
PlotMattForICHEP=1
PlotCompare=1

PlotAsymm=0
PlotFinalData=1
PlotSystematics=0
PlotLegend=0
PlotBrazilian=1
FitGraph=0
DrawLatexStuff=1
DrawPreliminary=1;
MultiPanelPlots=1
MPCentralsWithTotalSystID=June20Centrals_CentralsFromAlteredPPDJune20_1SigmaStatError #May11Centrals_CentralsFromAlteredPPDMay17_1SigmaStatErrorRelative #May11CentralWithTotalSyst #May11CentralWithTotalSyst May11NewTTreeCentralWithTotalSyst Apr27CentralWithTotalSyst
PlotAlteredPPDResults=1

#UncommentIFplotSystOverview
#PlotFinalData=0
#PlotSystematics=1
#PlotLegend=1
#PlotBrazilian=0
#FitGraph=0
#DrawLatexStuff=1
#MultiPanelPlots=0

#UncommentIFplotCDFComparison
#PlotMatt=1
#PlotCompare=1
#PlotAsymm=0
#PlotFinalData=1
#PlotSystematics=0
#PlotLegend=0
#PlotBrazilian=0
#FitGraph=0
#DrawLatexStuff=1
#MultiPanelPlots=1

DefaultID=Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT_AlteredPPD_June20_BKGlinPLUSRestSquaredGauss_5nRand
CompareID1=Data_TheGreatRun_10B_May11_NewestCentrals
CompareID2=Data_TheGreatRun_10B_May11_NewestCentrals_AlteredPPD_May17_BKGlinPLUSRestSquaredGauss_10nRand
CompareID3=BG0_Mar19_HighCtauSigCheck3p25
CompareID4=BG0_Mar19_HighCtauSigCheck3p5
nComp=0

nSystematics=0
#nSystematics=7

SystID1Base=TheGreatRun_BKGmodel
SystID1Specify=BestSyst_28p_SQRT12
SystID1Title=BKGmodel

SystID2Base=TheGreatRun_Rho
SystID2Specify=ConstSyst
SystID2Title=RhoFactor

SystID3Base=TheGreatRun_Param
SystID3Specify=BestSyst
SystID3Title=Parametrization

SystID4Base=TheGreatRun_TnP
SystID4Specify=BestSyst
SystID4Title=TnP_model

SystID5Base=TheGreatRun_FrameworkII
SystID5Specify=BestSyst_Bkg
SystID5Title=FrameworkIII

SystID6Base=TheGreatRun_FrameworkII
SystID6Specify=BestSyst_Sig_NoUnpol
SystID6Title=FrameworkII

SystID7Base=TheGreatRun_FrameworkI
SystID7Specify=BestSyst
SystID7Title=FrameworkI


ptBinMin=6
ptBinMax=10



########################################

cd ${homedir}

touch PlotFinalResults.cc
make

mkdir FinalResults
mkdir FinalResults/${JobID}

JobIDDir=FinalResults/${JobID}

mkdir ${JobIDDir}


./PlotFinalResults PlotAlteredPPDResults=${PlotAlteredPPDResults} MultiPanelPlots=${MultiPanelPlots} DrawLatexStuff=${DrawLatexStuff} ${MPCentralsWithTotalSystID}=MPCentralsWithTotalSystID ${DefaultID}=DefaultID ${CompareID1}=CompareID1 ${CompareID2}=CompareID2 ${CompareID3}=CompareID3 ${CompareID4}=CompareID4 ${JobID}=JobID ${SystID1Base}=SystID1Base ${SystID1Specify}=SystID1Specify ${SystID1Title}=SystID1Title ${SystID2Base}=SystID2Base ${SystID2Specify}=SystID2Specify ${SystID2Title}=SystID2Title ${basedir}=basedir ${storagedir}=storagedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nSystematics}nSystematics ${nComp}nComp ${nState}nState ${SystID3Base}=SystID3Base ${SystID3Specify}=SystID3Specify ${SystID3Title}=SystID3Title ${SystID4Base}=SystID4Base ${SystID4Specify}=SystID4Specify ${SystID4Title}=SystID4Title ${SystID5Base}=SystID5Base ${SystID5Specify}=SystID5Specify ${SystID5Title}=SystID5Title ${SystID6Base}=SystID6Base ${SystID6Specify}=SystID6Specify ${SystID6Title}=SystID6Title ${SystID7Base}=SystID7Base ${SystID7Specify}=SystID7Specify ${SystID7Title}=SystID7Title ${SystID8Base}=SystID8Base ${SystID8Specify}=SystID8Specify ${SystID8Title}=SystID8Title PlotMatt=${PlotMatt} PlotAsymm=${PlotAsymm} PlotCompare=${PlotCompare} PlotFinalData=${PlotFinalData} PlotSystematics=${PlotSystematics} PlotLegend=${PlotLegend} PlotBrazilian=${PlotBrazilian} FitGraph=${FitGraph} DrawPreliminary=${DrawPreliminary} PlotMattForICHEP=${PlotMattForICHEP}

cd ${homedir}/FinalResults/${JobID}/${nState}SUps
if [ ${PlotFinalData} -eq 1 ]
then
cp ${basedir}/latex/FinalDataResults.tex .
pdflatex FinalDataResults.tex
mv FinalDataResults.pdf FinalDataResults${additionalName}.pdf
fi
if [ ${PlotSystematics} -eq 1 ]
then
cp ${basedir}/latex/Systematics.tex .
pdflatex Systematics.tex
mv Systematics.pdf Systematics${additionalName}.pdf
fi
rm *.aux
rm *.log
rm *.tex

#rm -r Figures${additionalName}
mkdir Figures${additionalName}
mv Figures/* Figures${additionalName}/

done

