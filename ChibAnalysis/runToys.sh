#!/bin/sh


FitID=Apr20_10Toys50Bins

nToy=10
niSig=50
niBkg=50

nSig0=10
delta_nSig=15
nBkg0=10
delta_nBkg=50
 
run=1
plot=1

####### RUN #########

make
mkdir ToyFigures
mkdir ToyFigures/${FitID}
mkdir ToyFigures/${FitID}/SamplePlots

if [ ${run} -eq 1 ]
then
./minErrorToys ${FitID}=FitID ${nToy}nToy ${niBkg}niBkg ${niSig}niSig ${nSig0}nSig0 ${nBkg0}nBkg0 ${delta_nSig}delta_nSig ${delta_nBkg}delta_nBkg
fi
if [ ${plot} -eq 1 ]
then
./plotToys ${FitID}=FitID ${nToy}nToy ${niBkg}niBkg ${niSig}niSig ${nSig0}nSig0 ${nBkg0}nBkg0 ${delta_nSig}delta_nSig ${delta_nBkg}delta_nBkg
fi

