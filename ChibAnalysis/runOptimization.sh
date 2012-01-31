#!/bin/sh

# Reference vtxChi2ProbGamma nYsig ct gammapt vtxProb pTCut RConv

for cutVar in Reference;do

FitID=${cutVar}

Additionals=January31_NoPi0RejA_
fileName=NoPi0RejA
nToy=20000000

make
mkdir tmp
cp runChibFit runChibFit_${Additionals}${cutVar}

nState=2

./runChibFit_${Additionals}${cutVar} ${nState}nState ${Additionals}${FitID}=FitID ${cutVar}Cut ${fileName}=fileName ${nToy}nToy | tee Terminal_${Additionals}${cutVar}.txt

rm Figures/${Additionals}${cutVar}/bkgToy.root
mv Terminal_${Additionals}${cutVar}.txt tmp/Terminal_${Additionals}${cutVar}.txt
rm runChibFit_${Additionals}${cutVar}

done