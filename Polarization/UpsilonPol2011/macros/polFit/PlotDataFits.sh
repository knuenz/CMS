#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir

########## INPUTS ##########

for JobID in MCclosure_July25_Ups1S_MCTnPparam;do

#Take Care of Mean pT in ToyMC.h
for nState in 1;do

ptBinMin=1
ptBinMax=10

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2

additionalName=MPV${MPValgo}

############################


TreeID=${nState}SUps

frameSig=1
polScenSig=3

frameBkg=1
polScenBkg=3

nGenerations=1


rapBinMin=1 #don't change
rapBinMax=2 #don't change

Jobdir=${storagedir}/${JobID}

mkdir ${basedir}/macros/polFit/FiguresData
mkdir ${Jobdir}
mkdir ${Jobdir}/Figures
mkdir ${Jobdir}/Figures/${TreeID}
mkdir ${Jobdir}/Figures/${TreeID}/Figures
mkdir ${basedir}/macros/polFit/FiguresData/${JobID}
mkdir ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}

cd ${basedir}/macros/polFit

touch polRapPtPlot.cc

make

cd ${Jobdir}

cp ${basedir}/macros/polFit/polRapPtPlot polRapPtPlot_${TreeID}
echo 'copy finished'

for nSigma in 1;do #3 2 1

cd ${Jobdir}

./polRapPtPlot_${TreeID} ${nSigma}nSigma ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}Sigframe ${polScenSig}polScen ${MPValgo}MPValgo ${nGenerations}nGenerations ${TreeID}=TreeID realdata ${Jobdir}=dirstruct ${nState}nState

mv ${Jobdir}/TGraphResults_${TreeID}_temp.root ${Jobdir}/TGraphResults_${TreeID}.root 
cp ${Jobdir}/TGraphResults_${TreeID}.root ${Jobdir}/TGraphResults_${TreeID}_${nSigma}sigma.root 

cd ${Jobdir}/Figures/${TreeID}

cp ${basedir}/latex/DataResults_vs_RapPt.tex .
cp ${basedir}/latex/IndividualFitResults.tex ../../.
mv ${Jobdir}/ToyNumericalResults.tex .

pdflatex ToyNumericalResults.tex
mv ToyNumericalResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/DataNumericalResults_${additionalName}.pdf
rm *.aux
rm *.log


pdflatex DataResults_vs_RapPt.tex
mv DataResults_vs_RapPt.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/DataResults_vs_RapPt_${additionalName}.pdf

rm *.aux
rm *.log
rm DataResults_vs_RapPt.tex

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

cd ${Jobdir}/Figures/${TreeID}

filename=../lph_vs_lth_${TreeID}_rap${rap_}_pT${pT_}.pdf
if test -s "$filename"
then
cd ../..
	pdflatex "\newcommand\TreeBinID{${TreeID}_rap${rap_}_pT${pT_}}\input{IndividualFitResults.tex}"
	mv IndividualFitResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/IndividualFitResults_rap${rap_}pt${pT_}_${additionalName}.pdf
fi


pT_=$((pT_+1))
done
rap_=$((rap_+1))
done


done


cd ${Jobdir}


done

rm polRapPtPlot_${TreeID}
rm IndividualFitResults.tex
rm *.aux
rm *.log

done