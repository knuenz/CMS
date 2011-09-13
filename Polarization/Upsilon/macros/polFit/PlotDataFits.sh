#!/bin/sh

########## INPUTS ##########

storagedir=/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/polFit
basedir=/Users/valentinknuenz/usr/local/workspace/Upsilon

for JobID in ThirdDataFitTest;do

TreeID=Fake

ptBinMin=1
ptBinMax=8

############################

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
mkdir ${basedir}/macros/polFit/FiguresData/${JobID}

make

cd ${storagedir}/${JobID}

cp ${basedir}/macros/polFit/polRapPtPlot .


./polRapPtPlot ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}Sigframe ${polScenSig}polScen ${nGenerations}nGenerations ${TreeID}=TreeID realdata ${Jobdir}=dirstruct

cp ${basedir}/latex/DataResults_vs_RapPt.tex .
cp ${basedir}/latex/IndividualFitResults.tex .

pdflatex ToyNumericalResults.tex
mv ToyNumericalResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/ToyNumericalResults_${additionalName}.pdf
rm *.aux
rm *.log


pdflatex DataResults_vs_RapPt.tex
mv DataResults_vs_RapPt.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/DataResults_vs_RapPt_${additionalName}.pdf

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

filename=Figures/lph_vs_lth_${TreeID}_rap${rap_}_pT${pT_}.pdf
if test -s "$filename"
then
	pdflatex "\newcommand\TreeBinID{${TreeID}_rap${rap_}_pT${pT_}}\input{IndividualFitResults.tex}"
	mv IndividualFitResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/IndividualFitResults_rap${rap_}pt${pT_}.pdf
fi


pT_=$((pT_+1))
done
rap_=$((rap_+1))
done

rm *.aux
rm *.log
rm *.tex

cd ..
rm polRapPtPlot

done
