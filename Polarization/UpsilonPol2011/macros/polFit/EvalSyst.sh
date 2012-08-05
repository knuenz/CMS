#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir

########## INPUTS ##########

for nState in 1 2 3;do
for JobID1 in Data_CowboyFix_July16;do #this is the systematic error!!!

SystID=MassSigmaDep

JobID2=MCclosure_CombinedStates_MCtruthFineEta_1Sig   #this is the default!!!

ptBinMin=1
ptBinMax=10
statErrConsideration=false

########################################

cd ${homedir}

touch EvaluateSyst.cc
make

mkdir Systematics
mkdir Systematics/${SystID}

SystDir=Systematics/${SystID}/${JobID1}_TO_${JobID2}

mkdir ${SystDir}


./EvaluateSyst ${JobID1}=JobID1 ${JobID2}=JobID2 ${SystDir}=SystDir ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nState}nState statErrConsideration=${statErrConsideration}


done
done

