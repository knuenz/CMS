#!/bin/sh

# Reference vtxChi2ProbGammaLog nYsig ct gammapt vtxProb pTCut RConvUpper RConvLower Pi0Cut dzSigCut dzCut gammaptD ctSigCut

for cutVar in Reference;do

FitID=${cutVar}

Additionals=Apr23_chib_NewExtremeCuts1Scuts_WithX_
fileName=April9 #April5BugFix #Mar29refitted #April5 # #MC1P2P # # #Mar29refitted #vtxFixYMCQ # MC1P #
#fileName=Def
SaveAll=false

useSBforBKGmodel=false
useLeftSB=false
useRightSB=true

#BACKGROUND MODEL CHOICE:
useAnalyticalBKG=false

BkgMixer=true
useExistingMixFile=true
nMix=5000000

BkgToy=false
nToy=400000

PlotStep1=true
PlotStep2=true
finalPlots=false

FixStuffForStep2=true

MCsample=false

alteredToy=false
SetSignalToZero=false
DrawTextOnPlots=true
OneSPlotModel=true
restrictMaximum=true
NarrowOptimization=true
determine3PparametersDuringStep1=false

nState=1

####### RUN #########

make
mkdir tmp
cp runChibFit runChibFit_${Additionals}${cutVar}
rm BkgMixer_C*


./runChibFit_${Additionals}${cutVar} determine3PparametersDuringStep1=${determine3PparametersDuringStep1} NarrowOptimization=${NarrowOptimization} restrictMaximum=${restrictMaximum} alteredToy=${alteredToy} SetSignalToZero=${SetSignalToZero} DrawTextOnPlots=${DrawTextOnPlots} useLeftSB=${useLeftSB} useRightSB=${useRightSB} ${nState}nState ${Additionals}${FitID}=FitID ${cutVar}Cut ${fileName}=fileName ${nToy}nToy useSBforBKGmodel=${useSBforBKGmodel} SaveAll=${SaveAll} PlotStep1=${PlotStep1} PlotStep2=${PlotStep2} finalPlots=${finalPlots} FixStuffForStep2=${FixStuffForStep2} MCsample=${MCsample} useAnalyticalBKG=${useAnalyticalBKG} BkgMixer=${BkgMixer} BkgToy=${BkgToy} useExistingMixFile=${useExistingMixFile} ${nMix}nMix OneSPlotModel=${OneSPlotModel} | tee Terminal_${Additionals}${cutVar}.txt

#rm Figures/${Additionals}${cutVar}/bkgToy.root
mv Terminal_${Additionals}${cutVar}.txt tmp/Terminal_${Additionals}${cutVar}.txt
rm runChibFit_${Additionals}${cutVar}
rm BkgMixer_C*

done