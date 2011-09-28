#!/bin/sh

########## INPUTS ##########



for FidCuts in 1;do 				#defines the set of cuts to be used, see macros/polFit/effsAndCuts.h

for FracLSB in 50;do				#in percent; the left mass sideband will be weighted accordingly with respect to the right mass sideband
for nSigma in 2.00;do				#needed in 2 decimal accuracy (x.yz); this value decides in which dimuon mass region the data is projected

#following flags decide if the step is executed (1) or not (0):
execute_runData=1 					#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_runMassFit=1 				#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_runCopyTreeEntries=1		#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_PlotCosThetaPhiBG=1			#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_runTrimEventContent=1		#For each set of cuts, you can choose different values for FracLSB and nSigma
execute_runMeanPt=1					#Optional; This macro calculates the mean pT for the plotting (one has to manually copy the contents of meanPt.txt to ToyMC.h, different mean-pT of the three Y states are not yet implemented)

inputTree1=/afs/hephy.at/scratch/k/knuenz/TreeBuffer/TTree_Onia2MuMu_V9_PromptReco_v4.root
inputTree2=/afs/hephy.at/scratch/k/knuenz/TreeBuffer/TTree_Onia2MuMu_V9_May10ReReco_v1.root
#In case of more input Files: define inputTreeX and adapt the line starting with ./runData, implemented up to 4 Files



############################

rejectCowboys=true

CutDir=DataFiles/SetOfCuts${FidCuts}
mkdir DataFiles
mkdir ${CutDir}
mkdir ${CutDir}/Figures
mkdir ${CutDir}/tmpFiles
mkdir ${CutDir}/PDF
mkdir ${CutDir}/AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent

cp Makefile ${CutDir}/Makefile

cp runData.cc ${CutDir}/runData.cc
cp PolData.C ${CutDir}/PolData.C
cp PolData.h ${CutDir}/PolData.h
cp polFit/effsAndCuts.h ${CutDir}/effsAndCuts.h
cp ../interface/commonVar.h ${CutDir}/commonVar.h
cp ../interface/rootIncludes.inc ${CutDir}/rootIncludes.inc

cp runMassFit.cc ${CutDir}/runMassFit.cc
cp upsilon_2StepFit.C ${CutDir}/upsilon_2StepFit.C
cp CBFunction.C ${CutDir}/CBFunction.C

cp runCopyTreeEntries.cc ${CutDir}/runCopyTreeEntries.cc
cp CopyTreeEntries.C ${CutDir}/CopyTreeEntries.C
cp calcPol.C ${CutDir}/calcPol.C

cp PlotCosThetaPhiBG.cc ${CutDir}/PlotCosThetaPhiBG.cc

cp runTrimEventContent.cc ${CutDir}/runTrimEventContent.cc
cp TrimEventContent.C ${CutDir}/TrimEventContent.C

cp runMeanPt.cc ${CutDir}/runMeanPt.cc
cp calcMeanPt.C ${CutDir}/calcMeanPt.C

cp ../latex/massFits.tex ${CutDir}/massFits.tex
cp ../latex/cosThetaPhi_BG.tex ${CutDir}/cosThetaPhi_BG.tex

cd ${CutDir}

make


if [ ${execute_runData} -eq 1 ]
then
./runData rejectCowboys=${rejectCowboys} ${inputTree1}=inputTree1 ${inputTree2}=inputTree2 ${FidCuts}FidCuts
fi

if [ ${execute_runMassFit} -eq 1 ]
then
./runMassFit
pdflatex massFits.tex
mv massFits.pdf PDF/massFits.pdf
fi


if [ ${execute_runCopyTreeEntries} -eq 1 ]
then
./runCopyTreeEntries
fi

if [ ${execute_PlotCosThetaPhiBG} -eq 1 ]
then
./PlotCosThetaPhiBG
pdflatex cosThetaPhi_BG.tex
mv cosThetaPhi_BG.pdf PDF/cosThetaPhi_BG.pdf
fi

if [ ${execute_runTrimEventContent} -eq 1 ]
then
./runTrimEventContent ${FracLSB}=FracLSB ${nSigma}=nSigma
fi

if [ ${execute_runMeanPt} -eq 1 ]
then
./runMeanPt ${FracLSB}=FracLSB ${nSigma}=nSigma | tee AllStates_${nSigma}Sigma_FracLSB${FracLSB}Percent/meanPt.txt
fi



rm runData
rm runMassFit
rm runCopyTreeEntries
rm PlotCosThetaPhiBG
rm runTrimEventContent
rm runMeanPt
rm *.aux
rm *.log
rm *.tex
rm *.so
rm *.d

cd ..
cd ..

mkdir tmp
mv ${CutDir}/*.cc tmp/
mv ${CutDir}/*.C tmp/
mv ${CutDir}/*.h tmp/
mv ${CutDir}/*.inc tmp/
mv ${CutDir}/Makefile tmp/


done
done
done

