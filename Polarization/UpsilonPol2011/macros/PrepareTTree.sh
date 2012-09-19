#!/bin/sh

########## INPUTS ##########

JobID=MCclosure_Ups1S_July25 #MCclosure_Ups3S_Sept8 #July27_CowboyFix #MCclosure_Ups3S_July25 #MCclosure_Ups1S_July25 #May20Centrals #HighCtau3_BGratioTest #MCclosure #FinalData_CtauSig2 #_NoCtauCut #

for FidCuts in 11;do 				#defines the set of cuts to be used, see macros/polFit/effsAndCuts.h

for FracLSB in 25;do				#in percent; the left mass sideband will be weighted accordingly with respect to the right mass sideband
for nSigma in 10.00;do				#needed in 2 decimal accuracy (x.yz); this value decides in which dimuon mass region the data is projected

UpsMC=1
UpsMCstate=1 #enters only runMassFit in MC closure tests
RequestTrigger=1 ###set1
f_BG_zero=0
ProjectLSBdata=0
ProjectRSBdata=0
CombineSignalPeaks=0
Y1Sto2S_SB=0
LeftSided=0
RightSided=0
MassScan=0
DoCPUconsumingPlots=1 ###set1 for centrals
adjustOverlapBorders=1 ###set1 for centrals
#SELECTION
selectSOFT=0
selectTIGHT=0
selectMIXED=0
selectNOTMIXED=0

#following flags decide if the step is executed (1) or not (0):
execute_runData=0					#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_runMassFit=0				#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_runCopyTreeEntries=0		#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_PlotCosThetaPhiBG=0  		#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_runTrimEventContent=1		#For each set of cuts, you can choose different values for FracLSB and nSigma
execute_PlotCosThetaPhiDistribution=0	#For each set of cuts, you can choose different values for FracLSB and nSigma
execute_runMeanPt=0					#Optional; This macro calculates the mean pT for the plotting (one has to manually copy the contents of meanPt.txt to ToyMC.h, different mean-pT of the three Y states are not yet implemented)

#OLDtree:
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Upsi_Onia2MuMu_v20_PromptReco_v4.root
#inputTree2=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Upsi_Onia2MuMu_v20_PromptReco_v5.root
#inputTree3=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Upsi_Onia2MuMu_v20_PromptReco_v6.root
#inputTree4=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Upsi_Onia2MuMu_v20_RunB_PromptReco_v1_full.root
#NEWtree (Linlin 3 days before pre-approval):
inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Upsi.root
#MIXEDtree:
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_13June2012_PromptRecoBV1.root
#inputTree2=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_13June2012_PromptRecoV4.root
#inputTree3=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_13June2012_PromptRecoV5.root
#inputTree4=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_13June2012_PromptRecoV6.root
#NEW MIXED TTREE
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_16June2012_PromptRecoBV1_v1.root
#inputTree2=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_16June2012_PromptRecoV4_v1.root
#inputTree3=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_16June2012_PromptRecoV5_v1.root
#inputTree4=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v20_Upsi_16June2012_PromptRecoV6_v1.root

#MCUps1Stree:
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_Onia2MuMu_v10_UpsiLowPt_PGun_HLT_1E33-3E33.root
#SingleMuonTriggertree:
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/TTree_SingleMu_Onia2MuMu_v20_PromptRecoAB.root

#ILSE MC
#inputTree1=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/TTree_Ups1S_pt0-30_20June2012.root 
#inputTree1=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/TTree_Onia2MuMu_Upsi2s_measuredPt30_50_17July2012.root
#inputTree2=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/TTree_Onia2MuMu_Upsi2s_measuredPt_2April2012.root
inputTree1=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/TTree_Onia2MuMu_Upsi3s_measuredPt_2April2012.root

#ILSE MC Ups3S High pT
#inputTree1=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/Generation/Tree_MC_Y3S_measuredPt28-52_1.root
inputTree2=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/Generation/Tree_MC_Y3S_measuredPt28-52_2.root
inputTree3=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/Generation/Tree_MC_Y3S_measuredPt28-52_3.root
inputTree4=/scratch/ikratsch/Polarization/Upsilon/InputFiles/MC/Generation/Tree_MC_Y3S_measuredPt28-52_4.root 

#HERMINE NEW 1S (realistic pT shape):
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/onia2MuMu_tree_Ups1S_mergeAllBins_1.root
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/onia2MuMu_tree_Ups1S_pT6.root
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/onia2MuMu_tree_Ups1S_pT7.root
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/onia2MuMu_tree_Ups1S_pT8.root
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/onia2MuMu_tree_Ups1S_pT9.root
#inputTree1=/scratch/knuenz/Polarization/RootInput/Upsilon/onia2MuMu_tree_Ups1S_pT10.root



#In case of more input Files: define inputTreeX and adapt the line starting with ./runData, implemented up to 4 Files



############################

rejectCowboys=true

CutDir=DataFiles/SetOfCuts${FidCuts}_${JobID}
mkdir DataFiles
mkdir ${CutDir}
mkdir ${CutDir}/Figures
mkdir ${CutDir}/tmpFiles
mkdir ${CutDir}/PDF
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/Figures
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/Figures/Ups1S
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/Figures/Ups2S
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/Figures/Ups3S
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/PDF

cp Makefile ${CutDir}/Makefile

cp runData.cc ${CutDir}/runData.cc
cp PolData.C ${CutDir}/PolData.C
cp PolData.h ${CutDir}/PolData.h
cp polFit/effsAndCuts.h ${CutDir}/effsAndCuts.h
cp ../interface/commonVar.h ${CutDir}/commonVar.h
cp ../interface/rootIncludes.inc ${CutDir}/rootIncludes.inc

cp runMassFit.cc ${CutDir}/runMassFit.cc
cp runMassFitMC.cc ${CutDir}/runMassFitMC.cc
cp upsilon_2StepFit.C ${CutDir}/upsilon_2StepFit.C
cp upsilon_MCMassFit.C ${CutDir}/upsilon_MCMassFit.C
cp CBFunction.C ${CutDir}/CBFunction.C

cp runCopyTreeEntries.cc ${CutDir}/runCopyTreeEntries.cc
cp CopyTreeEntries.C ${CutDir}/CopyTreeEntries.C
cp calcPol.C ${CutDir}/calcPol.C

cp PlotCosThetaPhiBG.cc ${CutDir}/PlotCosThetaPhiBG.cc
cp PlotCosThetaPhiDistribution.cc ${CutDir}/PlotCosThetaPhiDistribution.cc

cp runTrimEventContent.cc ${CutDir}/runTrimEventContent.cc
cp TrimEventContent.C ${CutDir}/TrimEventContent.C

cp runMeanPt.cc ${CutDir}/runMeanPt.cc
cp calcMeanPt.C ${CutDir}/calcMeanPt.C

cp ../latex/massFits.tex ${CutDir}/massFits.tex
cp ../latex/cosThetaPhi_BG.tex ${CutDir}/cosThetaPhi_BG.tex
cp ../latex/cosThetaPhi.tex ${CutDir}/cosThetaPhi.tex

cd ${CutDir}

make


if [ ${execute_runData} -eq 1 ]
then
./runData rejectCowboys=${rejectCowboys} ${inputTree1}=inputTree1 ${inputTree2}=inputTree2 ${inputTree3}=inputTree3 ${inputTree4}=inputTree4 ${FidCuts}FidCuts UpsMC=${UpsMC} RequestTrigger=${RequestTrigger} selectSOFT=${selectSOFT} selectTIGHT=${selectTIGHT} selectMIXED=${selectMIXED} selectNOTMIXED=${selectNOTMIXED}
fi

if [ ${execute_runMassFit} -eq 1 ]
then
if [ ${UpsMC} -eq 0 ]
then
./runMassFit
fi
if [ ${UpsMC} -eq 1 ]
then
./runMassFitMC ${UpsMCstate}=nState
fi
pdflatex massFits.tex
mv massFits.pdf PDF/massFits.pdf
fi


if [ ${execute_runCopyTreeEntries} -eq 1 ]
then
./runCopyTreeEntries UpsMC=${UpsMC} DoCPUconsumingPlots=${DoCPUconsumingPlots}
fi

if [ ${execute_PlotCosThetaPhiBG} -eq 1 ]
then
./PlotCosThetaPhiBG
pdflatex cosThetaPhi_BG.tex
mv cosThetaPhi_BG.pdf PDF/cosThetaPhi_BG.pdf
fi

if [ ${execute_runTrimEventContent} -eq 1 ]
then
./runTrimEventContent ${FracLSB}=FracLSB ${nSigma}=nSigma UpsMC=${UpsMC} f_BG_zero=${f_BG_zero} ProjectLSBdata=${ProjectLSBdata} ProjectRSBdata=${ProjectRSBdata} CombineSignalPeaks=${CombineSignalPeaks} Y1Sto2S_SB=${Y1Sto2S_SB} LeftSided=${LeftSided} RightSided=${RightSided} MassScan=${MassScan} adjustOverlapBorders=${adjustOverlapBorders}
fi

if [ ${execute_PlotCosThetaPhiDistribution} -eq 1 ]
then
./PlotCosThetaPhiDistribution 1nState AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent=DataPath
./PlotCosThetaPhiDistribution 2nState AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent=DataPath
./PlotCosThetaPhiDistribution 3nState AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent=DataPath
cp cosThetaPhi.tex AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/cosThetaPhi.tex
cd AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent
pdflatex "\newcommand\TreeBinID{Ups1S}\input{cosThetaPhi.tex}"
mv cosThetaPhi.pdf PDF/cosThetaPhi_1SUps.pdf
pdflatex "\newcommand\TreeBinID{Ups2S}\input{cosThetaPhi.tex}"
mv cosThetaPhi.pdf PDF/cosThetaPhi_2SUps.pdf
pdflatex "\newcommand\TreeBinID{Ups3S}\input{cosThetaPhi.tex}"
mv cosThetaPhi.pdf PDF/cosThetaPhi_3SUps.pdf
cd ..
fi

if [ ${execute_runMeanPt} -eq 1 ]
then
./runMeanPt ${FracLSB}=FracLSB ${nSigma}=nSigma | tee AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/meanPt.txt
fi



rm runData
rm runMassFit
rm runCopyTreeEntries
rm PlotCosThetaPhiBG
rm runTrimEventContent
rm runMeanPt
rm *.aux
rm *.log
rm *.tex
rm *.so
rm *.d

cd ..
cd ..

mkdir tmp
mv ${CutDir}/*.cc tmp/
mv ${CutDir}/*.C tmp/
mv ${CutDir}/*.h tmp/
mv ${CutDir}/*.inc tmp/
mv ${CutDir}/Makefile tmp/


done
done
done

