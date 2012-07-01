#!/bin/sh

########## INPUTS ##########

nState=1

JobID=Toy_TheGreatRun_ARC_June23_1S_PhiRejection8th

nGenerations=50

rapBinMin=1
rapBinMax=1
ptBinMin=8
ptBinMax=8

polScenSig=3
polScenBkg=3
frameSig=2
frameBkg=1

nEff=1030
UseMCeff=true
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=1

FidCuts=11

nSample=10000
nSkipGen=0

#GENERATION SETTINGS
ConstEvents=10000
UseConstEv=true

UseDifferingEff=false
nEffRec=1030
UseMCReceff=true
nDileptonEffRec=1
UseMCDileptonReceff=true
nRecRhoFactor=1

gen=true
rec=true
fit=true
plot=false
deletePseudoData=true

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
nSigma=2
NewAccCalc=false


########################################

cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/ToyMC #please define the directory storagedir in the file macros/polFit/storagedir



touch polGenRecFitPlot.cc
make

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir ${storagedir}
mkdir ${storagedir}/${JobID}
mkdir ${storagedir}/${JobID}/${ScenDir}
cp ${basedir}/macros/polFit/polGenRecFitPlot ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/${ScenDir}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/${ScenDir}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/${ScenDir}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/${ScenDir}/polPlot.C

cd ${storagedir}/${JobID}/${ScenDir}
cp ${basedir}/macros/polFit/runToyMC.sh .

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

plot=${plot}

cp ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_}

if [ ${nGen_} -eq 1 ]
then
plot=true
fi

./polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${nEffRec}nRecEff ${nDileptonEffRec}nRecDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} UseDifferingEff=${UseDifferingEff} UseMCeff=${UseMCeff} UseMCReceff=${UseMCReceff} UseMCDileptoneff=${UseMCDileptoneff} UseMCDileptonReceff=${UseMCDileptonReceff}  ${nRhoFactor}nRhoFactor ${nRecRhoFactor}nRecRhoFactor ${MPValgo}MPValgo ${nSigma}nSigma ${nState}nState NewAccCalc=${NewAccCalc} deletePseudoData=${deletePseudoData}

rm polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_}

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

