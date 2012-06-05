#!/bin/sh

########## INPUTS ##########

nState=3

JobID=Data_TheGreatRun_10B_May20_NewestCentrals_BGmodel_PLUS_28_46_30

nGenerations=10
MergeFiles=1

rapBinMin=2
rapBinMax=2
ptBinMin=2
ptBinMax=2


########################################

TreeID=${nState}SUps

##############################################

cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir

mkdir ${storagedir}
mkdir ${storagedir}/${JobID}

cd ${storagedir}/${JobID}
pwd

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


resultfilename=resultsMerged_${nState}SUps_rap${rap_}_pT${pT_}.root


nFitResultName=results_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.root
echo ${nFitResultName}

if [ $MergeFiles -eq 1 ]
then

if [ $nGen_ -eq 1 ]
then
cp ${nFitResultName} ${resultfilename}
fi

if [ $nGen_ -ge 2 ]
then
mv ${resultfilename} BUFFER_${resultfilename}
hadd -f ${resultfilename} BUFFER_${resultfilename} ${nFitResultName}
rm BUFFER_${resultfilename}
fi

fi



nGen_=$((nGen_+1))
done


FinalDestinationName=results_${nState}SUps_rap${rap_}_pT${pT_}.root

if [ $MergeFiles -eq 1 ]
then
mv ${FinalDestinationName} tmp_${FinalDestinationName}
mv ${resultfilename} ${FinalDestinationName}
fi

pT_=$((pT_+1))
done
rap_=$((rap_+1))
done


