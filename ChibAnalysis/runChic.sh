#!/bin/sh

# Reference vtxChi2ProbGammaLog nYsig ct gammapt vtxProb pTCut RConvUpper RConvLower

for cutVar in Reference;do

FitID=${cutVar}

Additionals=Apr21_chic_Step1_
fileName=ChicApr19
SaveAll=true

#BACKGROUND MODEL CHOICE:
useAnalyticalBKG=true

BkgMixer=false
useExistingMixFile=true
nMix=1000000

BkgToy=false
nToy=200000

useSBforBKGmodel=false
useLeftSB=true
useRightSB=true

alteredToy=false
SetSignalToZero=false
DrawTextOnPlots=true

Search2P=false

make
mkdir tmp
cp runChicFit runChicFit_${Additionals}${cutVar}
rm BkgMixer_C*

nState=1

./runChicFit_${Additionals}${cutVar} alteredToy=${alteredToy} SetSignalToZero=${SetSignalToZero} DrawTextOnPlots=${DrawTextOnPlots} useLeftSB=${useLeftSB} useRightSB=${useRightSB} ${nState}nState ${Additionals}${FitID}=FitID ${cutVar}Cut ${fileName}=fileName ${nToy}nToy useSBforBKGmodel=${useSBforBKGmodel} SaveAll=${SaveAll} useAnalyticalBKG=${useAnalyticalBKG} Search2P=${Search2P} BkgMixer=${BkgMixer} BkgToy=${BkgToy} useExistingMixFile=${useExistingMixFile} ${nMix}nMix | tee Terminal_${Additionals}${cutVar}.txt

#rm Figures/${Additionals}${cutVar}/bkgToy.root
mv Terminal_${Additionals}${cutVar}.txt tmp/Terminal_${Additionals}${cutVar}.txt
rm runChicFit_${Additionals}${cutVar}
rm BkgMixer_C*

done