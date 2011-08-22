#!/bin/sh

########## INPUTS ##########

storagedir=/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC
basedir=/Users/valentinknuenz/usr/local/workspace/Upsilon

JobID=LowStatTry

rapBinMin=1#don't change
rapBinMax=2#don't change
ptBinMin=8
ptBinMax=8

frameSig=1
for polScenSig in 3;do

frameBkg=1
polScenBkg=3

nGenerations=50

additionalName=


############################

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir ${basedir}/macros/ToyMC/Figures
mkdir ${basedir}/macros/ToyMC/Figures/${JobID}
mkdir ${basedir}/macros/ToyMC/Figures/${JobID}/${ScenDir}
 

cd ${storagedir}/${JobID}
mkdir ${ScenDir}
cp ${basedir}/macros/ToyMC/polToyMCPlot.cc .
cp ${basedir}/macros/ToyMC/Makefile .
cp ${basedir}/interface/commonVar.h .
cp ${basedir}/interface/rootIncludes.inc .
cp ${basedir}/macros/ToyMC/ToyMC.h .
make


./polToyMCPlot ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}Sigframe ${polScenSig}polScen ${nGenerations}nGenerations ${ScenDir}

cp ${basedir}/latex/PullSummaryResults.tex ${ScenDir}/PullSummaryResults_${ScenDir}.tex
cp ${basedir}/latex/ParameterSummaryResults.tex ${ScenDir}/ParameterSummaryResults_${ScenDir}.tex
cp ${basedir}/latex/ToyResults.tex ${ScenDir}/ToyResults_${ScenDir}.tex

pdflatex ToyNumericalResults_${ScenDir}.tex
mv ToyNumericalResults_${ScenDir}.pdf ${basedir}/macros/ToyMC/Figures/${JobID}/${ScenDir}/ToyNumericalResults_${ScenDir}_${additionalName}.pdf
rm *.aux
rm *.log

cd ${ScenDir}
pdflatex PullSummaryResults_${ScenDir}.tex
pdflatex ParameterSummaryResults_${ScenDir}.tex

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

pdflatex "\newcommand\rappt{rap${rap_}pt${pT_}}\input{ToyResults_${ScenDir}.tex}"
mv ToyResults_${ScenDir}.pdf ${basedir}/macros/ToyMC/Figures/${JobID}/${ScenDir}/ToyResults_${ScenDir}_rap${rap_}pt${pT_}.pdf

pT_=$((pT_+1))
done
rap_=$((rap_+1))
done

mv PullSummaryResults_${ScenDir}.pdf ${basedir}/macros/ToyMC/Figures/${JobID}/${ScenDir}/PullSummaryResults_${ScenDir}_${additionalName}.pdf
mv ParameterSummaryResults_${ScenDir}.pdf ${basedir}/macros/ToyMC/Figures/${JobID}/${ScenDir}/ParameterSummaryResults_${ScenDir}_${additionalName}.pdf

rm *.aux
rm *.log
rm *.tex

cd ..
rm polToyMCPlot.cc
rm polToyMCPlot
rm Makefile
rm rootIncludes.inc
rm commonVar.h
rm ToyMC.h

done