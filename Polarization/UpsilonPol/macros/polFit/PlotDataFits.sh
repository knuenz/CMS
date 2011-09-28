#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`pwd` #please define differently, if you have storage issues

########## INPUTS ##########

for JobID in DefaultID;do

for nState in 1 2 3;do


ptBinMin=1
ptBinMax=12

############################


TreeID=${nState}SUps

frameSig=1
polScenSig=3

frameBkg=1
polScenBkg=3

nGenerations=1

additionalName=

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

make

cd ${Jobdir}

cp ${basedir}/macros/polFit/polRapPtPlot .


./polRapPtPlot ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}Sigframe ${polScenSig}polScen ${nGenerations}nGenerations ${TreeID}=TreeID realdata ${Jobdir}=dirstruct

cd ${Jobdir}/Figures/${TreeID}

cp ${basedir}/latex/DataResults_vs_RapPt.tex .
cp ${basedir}/latex/IndividualFitResults.tex ../../.
cp ${Jobdir}/ToyNumericalResults.tex .

pdflatex ToyNumericalResults.tex
mv ToyNumericalResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/ToyNumericalResults_${additionalName}.pdf
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
	mv IndividualFitResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/IndividualFitResults_rap${rap_}pt${pT_}.pdf
fi


pT_=$((pT_+1))
done
rap_=$((rap_+1))
done

rm *.aux
rm *.log
rm *.tex

done


cd ${Jobdir}

rm polRapPtPlot

done
