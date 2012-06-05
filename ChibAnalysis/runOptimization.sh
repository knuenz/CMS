#!/bin/sh

# Reference vtxChi2ProbGammaLog nYsig ct gammapt vtxProb pTCut RConvUpper RConvLower Pi0Cut dzSigCut dzCut gammaptD ctSigCut

for cutVar in Reference;do

FitID=${cutVar}

Additionals=June5_chib_ErnestcutsJune4_NewBkgMixer_nYsigCorr_Max11_WithX_
fileName=April9
SaveAll=false

#BACKGROUND MODEL CHOICE:
useAnalyticalBKG=false

BkgMixer=true
useExistingMixFile=true #OldMixer
useNewExistingMixFile=true #NewMixer, if set true, the old mixer is inactive
nMix=5000000

BkgToy=false
nToy=1000000


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
determine3PparametersDuringStep1=true
BkgToy1SBkgMixer2S=false #if true: BkgToy=true, BkgMixer=true, useExistingMixFile=true
useSBforBKGmodel=false
useLeftSB=false
useRightSB=true


nState=3

####### RUN #########

touch runChibFit.cc
make
mkdir tmp
cp runChibFit runChibFit_${Additionals}${cutVar}
rm BkgMixer_C*


./runChibFit_${Additionals}${cutVar} determine3PparametersDuringStep1=${determine3PparametersDuringStep1} NarrowOptimization=${NarrowOptimization} restrictMaximum=${restrictMaximum} alteredToy=${alteredToy} SetSignalToZero=${SetSignalToZero} DrawTextOnPlots=${DrawTextOnPlots} useLeftSB=${useLeftSB} useRightSB=${useRightSB} ${nState}nState ${Additionals}${FitID}=FitID ${cutVar}Cut ${fileName}=fileName ${nToy}nToy useSBforBKGmodel=${useSBforBKGmodel} SaveAll=${SaveAll} PlotStep1=${PlotStep1} PlotStep2=${PlotStep2} finalPlots=${finalPlots} FixStuffForStep2=${FixStuffForStep2} MCsample=${MCsample} useAnalyticalBKG=${useAnalyticalBKG} BkgMixer=${BkgMixer} BkgToy=${BkgToy} useExistingMixFile=${useExistingMixFile} ${nMix}nMix OneSPlotModel=${OneSPlotModel} BkgToy1SBkgMixer2S=${BkgToy1SBkgMixer2S} useNewExistingMixFile=${useNewExistingMixFile} | tee Terminal_${Additionals}${cutVar}.txt

#rm Figures/${Additionals}${cutVar}/bkgToy.root
mv Terminal_${Additionals}${cutVar}.txt tmp/Terminal_${Additionals}${cutVar}.txt
rm runChibFit_${Additionals}${cutVar}
rm BkgMixer_C*

done
