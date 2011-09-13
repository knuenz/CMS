#!/bin/sh

########## INPUTS ##########

datadir=/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/RootFiles
storagedir=/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/polFit
basedir=/Users/valentinknuenz/usr/local/workspace/Upsilon

TreeID=Fake

JobID=ThirdDataFitTest

rapBinMin=1
rapBinMax=1
ptBinMin=5
ptBinMax=8

nEff=3
FidCuts=1

nSample=6000

########################################

polScenSig=3
polScenBkg=3
frameSig=1
frameBkg=1
nGenerations=1
ConstEvents=15000
nSkipGen=0
UseConstEv=true

gen=false
rec=false
fit=true
plot=true

touch polGenRecFitPlot.cc
make

ScenDir=Default_ScenDir
mkdir ${storagedir}
mkdir ${storagedir}/${JobID}

cp ${basedir}/macros/polFit/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/polPlot.C

cd ${storagedir}/${JobID}

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

nGen_=${nSkipGen}
nGen_=$((nGen_+1))
while [ $nGen_ -le ${nGenerations} ]
do


cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_rap${rap_}_pt${pT_}


./polGenRecFitPlot_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir

rm polGenRecFitPlot_rap${rap_}_pt${pT_}

nGen_=$((nGen_+1))
done
pT_=$((pT_+1))
done
rap_=$((rap_+1))
done

#rm polGen.C 
#rm polRec.C 
#rm polFit.C 
#rm polPlot.C 

cp ${basedir}/macros/polFit/runDataFits.sh .
