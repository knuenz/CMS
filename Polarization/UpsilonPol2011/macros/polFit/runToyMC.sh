#!/bin/sh

########## INPUTS ##########

nState=2

JobID=Toy_TowardsPRL_Aug03_AmapsTest_1030t32104o105

nGenerations=50

rapBinMin=1
rapBinMax=1
ptBinMin=10
ptBinMax=10

polScenSig=3
polScenBkg=3
frameSig=1
frameBkg=1

nEff=1030
UseMCeff=false
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=1

useAmapApproach=true
nAmap=32104                     #frame/state/sigma/ID ( ID= 2 digits )
nDenominatorAmap=105 		    #the number here corresponds to the same notation as nEff
 
FidCuts=11

nSample=10000
nSkipGen=0

#GENERATION SETTINGS
ConstEvents=15000
UseConstEv=true

UseDifferingEff=false
nEffRec=1101
UseMCReceff=false
nDileptonEffRec=1
UseMCDileptonReceff=true
nRecRhoFactor=1

gen=true
rec=true
fit=true
plot=false
deletePseudoData=true

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
nSigma=1
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

./polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${nEffRec}nRecEff ${nDileptonEffRec}nRecDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} UseDifferingEff=${UseDifferingEff} UseMCeff=${UseMCeff} UseMCReceff=${UseMCReceff} UseMCDileptoneff=${UseMCDileptoneff} UseMCDileptonReceff=${UseMCDileptonReceff}  ${nRhoFactor}nRhoFactor ${nRecRhoFactor}nRecRhoFactor ${MPValgo}MPValgo ${nSigma}nSigma ${nState}nState NewAccCalc=${NewAccCalc} deletePseudoData=${deletePseudoData} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap

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

