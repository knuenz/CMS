#!/bin/sh

# Reference vtxChi2ProbGammaLog nYsig ct gammapt vtxProb pTCut RConvUpper RConvLower

for cutVar in Reference;do

FitID=${cutVar}

Additionals=tmp_Chib_Foldername_
fileName=vtxFixYMCQ
nToy=200000
SaveAll=false

useSBforBKGmodel=false
useLeftSB=false
useRightSB=true

alteredToy=false
SetSignalToZero=false
DrawTextOnPlots=true


make
mkdir tmp
cp runChibFit runChibFit_${Additionals}${cutVar}

nState=2

./runChibFit_${Additionals}${cutVar} alteredToy=${alteredToy} SetSignalToZero=${SetSignalToZero} DrawTextOnPlots=${DrawTextOnPlots} useLeftSB=${useLeftSB} useRightSB=${useRightSB} ${nState}nState ${Additionals}${FitID}=FitID ${cutVar}Cut ${fileName}=fileName ${nToy}nToy useSBforBKGmodel=${useSBforBKGmodel} SaveAll=${SaveAll} | tee Terminal_${Additionals}${cutVar}.txt

rm Figures/${Additionals}${cutVar}/bkgToy.root
mv Terminal_${Additionals}${cutVar}.txt tmp/Terminal_${Additionals}${cutVar}.txt
rm runChibFit_${Additionals}${cutVar}

done