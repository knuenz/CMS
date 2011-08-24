#!/bin/sh

########## INPUTS ##########

storagedir=/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC
basedir=/Users/valentinknuenz/usr/local/workspace/Upsilon

JobID=ToyTest

nGenerations=50;

rapBinMin=1;
rapBinMax=2;
ptBinMin=1;
ptBinMax=8;

polScenSig=3;
polScenBkg=3;
frameSig=1;
frameBkg=1;

nEff=3;
FidCuts=0;

nSample=6000;
ConstEvents=15000;
nSkipGen=0;

UseConstEv=true;

gen=true;
rec=true;
fit=true;

########################################

touch LoopToyMC.cc
make

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}
mkdir ${storagedir}/${JobID}
mkdir ${storagedir}/${JobID}/${ScenDir}
cp ${basedir}/macros/ToyMC/LoopToyMC ${storagedir}/${JobID}/${ScenDir}/LoopToyMC
cp ${basedir}/macros/ToyMC/polGen.C ${storagedir}/${JobID}/${ScenDir}/polGen.C
cp ${basedir}/macros/ToyMC/polRec.C ${storagedir}/${JobID}/${ScenDir}/polRec.C
cp ${basedir}/macros/ToyMC/polFit.C ${storagedir}/${JobID}/${ScenDir}/polFit.C

cd ${storagedir}/${JobID}/${ScenDir}

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


cp ${storagedir}/${JobID}/${ScenDir}/LoopToyMC ${storagedir}/${JobID}/${ScenDir}/LoopToyMC_rap${rap_}_pt${pT_}_Gen${nGen_}


./LoopToyMC_rap${rap_}_pt${pT_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit}

rm LoopToyMC_rap${rap_}_pt${pT_}_Gen${nGen_}

nGen_=$((nGen_+1))
done
pT_=$((pT_+1))
done
rap_=$((rap_+1))
done

rm polGen.C 
rm polRec.C 
rm polFit.C 

cp ${basedir}/macros/ToyMC/runToyMC.sh .