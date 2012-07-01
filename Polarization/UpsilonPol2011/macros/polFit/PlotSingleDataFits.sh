#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
datadir_Start=${basedir}/macros/DataFiles

########## INPUTS ##########

#fracL=0 #in percent
nSigma=1.00 #needed in 2 decimal accuracy (x.yz)

for nState in 1;do

JobID=Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT_AlteredPPD_June20_BKGlinPLUSRestSquaredGauss_5nRand #Data_TheGreatRun_10B_May11_NewestCentrals

rapBinMin=1
rapBinMax=1
ptBinMin=10
ptBinMax=10

#NOTE: take care of nFits!
nFits=50

########################################

FidCuts=11
nEff=1030
UseMCeff=true
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=1
nSample=60000
nSkipGen=0
DataID=_FinalData_CtauSig2_RequestTrigger
MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2


useCentralFracL=1

if [ $useCentralFracL -eq 1 ]
then
if [ $nState -eq 1 ]
then
fracL=72 #72
fi
if [ $nState -eq 2 ]
then
fracL=46 #46
fi
if [ $nState -eq 3 ]
then
fracL=30 #30
fi
fi

datadir=${datadir_Start}/SetOfCuts${FidCuts}${DataID}/AllStates_${nSigma}Sigma_FracLSB${fracL}Percent

TreeID=${nState}SUps

cd ${homedir}


polScenSig=3
polScenBkg=3
frameSig=1
frameBkg=1
ConstEvents=15000
UseConstEv=true
nGenerations=${nFits}
gen=false
rec=false
fit=false
plot=true
NewAccCalc=false

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
cp ${basedir}/macros/polFit/runDataFits.sh .


rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

nGen_=${nSkipGen}
nGen_=$((nGen_+1))
nMaxGen=$((nGenerations+nSkipGen))
while [ $nGen_ -le $nMaxGen ]
do

resultfilename=resultsMerged_${nState}SUps_rap${rap_}_pT${pT_}.root
nActualGen=$((nGen_-nSkipGen))

nGen_=$((nGen_+1))
done


cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_rap${rap_}_pt${pT_}
./polGenRecFitPlot_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${nState}nState ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo scalePlots=true NewAccCalc=${NewAccCalc}
rm polGenRecFitPlot_rap${rap_}_pt${pT_}


pT_=$((pT_+1))
done
rap_=$((rap_+1))
done


done

rm *.so
rm *.d
#rm polGenRecFitPlot

mkdir ../tmp
#mv *.C ../tmp

