#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir

########## INPUTS ##########

for nState in 1 2 3;do

SystID=TheGreatRun_DeltaTilde

JobID=/scratch/knuenz/Polarization/Upsilon/Data/Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT_AlteredPPD_June20_BKGlinPLUSRestSquaredGauss_5nRand

ptBinMin=1
ptBinMax=10

########################################

cd ${homedir}

touch ChangeTGraph.cc
make

mkdir Systematics
mkdir Systematics/${SystID}

SystDir=Systematics/${SystID}/ChangedTGraph

mkdir ${SystDir}


./ChangeTGraph ${JobID}=JobID1 ${SystID}=SystID ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nState}nState


done


