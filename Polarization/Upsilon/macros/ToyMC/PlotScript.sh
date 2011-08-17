#!/bin/sh

########## INPUTS ##########

#storagedir=/scratch/knuenz/Polarization/Upsilon/ToyMC
#basedir=/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/Upsilon
storagedir=/Users/valentinknuenz/usr/local/workspace/Upsilon/macros/ToyMC
basedir=/Users/valentinknuenz/usr/local/workspace/Upsilon

JobID=Test

rapBinMin=2
rapBinMax=2
ptBinMin=8
ptBinMax=8

frameSig=1
polScenSig=3

frameBkg=1
polScenBkg=3

nGenerations=5

############################

cd ${storagedir}/${JobID}

cp ${basedir}/macros/ToyMC/polToyMCPlot.cc .
cp ${basedir}/macros/ToyMC/Makefile .
make

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

./polToyMCPlot ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}Sigframe ${polScenSig}polScen ${nGenerations}nGenerations ${ScenDir}

cp ${basedir}/latex/PullSummaryResults.tex ${ScenDir}/PullSummaryResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.tex
cp ${basedir}/latex/ParameterSummaryResults.tex ${ScenDir}/ParameterSummaryResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.tex
cp ${basedir}/latex/ToyResults.tex ${ScenDir}/ToyResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.tex


cd ${ScenDir}
pdflatex PullSummaryResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.tex
pdflatex ParameterSummaryResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.tex

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

pdflatex "\newcommand\rappt{rap${rap_}pt${pT_}}\input{ToyResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.tex}"
mv ToyResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.pdf ../ToyResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations_rap${rap_}pt${pT_}.pdf

pT_=$((pT_+1))
done
rap_=$((rap_+1))
done

mv PullSummaryResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.pdf ../.
mv ParameterSummaryResults_Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}_${nGenerations}Generations.pdf ../.

cd ..
rm polToyMCPlot.cc
rm polToyMCPlot
rm Makefile
