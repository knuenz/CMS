#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
storagedir=`pwd` #please define differently, if you have storage issues
datadir_Start=${basedir}/macros/DataFiles

########## INPUTS ##########

fracL=50 #in percent
nSigma=2.00 #needed in 2 decimal accuracy (x.yz)

for nState in 1;do

JobID=DefaultID #Please define nSigma and fracL yourself in the JobID, if needed

rapBinMin=1
rapBinMax=2
ptBinMin=1
ptBinMax=12

nEff=6	 				#nEff > 100 corresponds to binned efficiencies (specify file name in effsAndCuts.h and copy files to basedir/macros/polFit/EffFiles)
nDileptonEff=1			#nDileptonEff > 200 corresponds to binned efficiencies (specify file name in effsAndCuts.h and copy files to basedir/macros/polFit/EffFiles)
FidCuts=1

nSample=6000

nFits=1
nSkipGen=0

########################################


datadir=${datadir_Start}/SetOfCuts${FidCuts}/AllStates_${nSigma}Sigma_FracLSB${fracL}Percent

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
fit=true
plot=true

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
if [ $nSkipGen -ge 0 ]
then
if [ $nActualGen -eq 1 ]
then
cp results_${nState}SUps_rap${rap_}_pT${pT_}.root ${resultfilename}
fi
fi


cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_rap${rap_}_pt${pT_}
./polGenRecFitPlot_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir
rm polGenRecFitPlot_rap${rap_}_pt${pT_}

mv Figures/fit_CS_costh_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_CS_costh_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_CS_phi_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_CS_phi_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_CS_phith_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_CS_phith_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_HX_costh_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_HX_costh_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_HX_phi_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_HX_phi_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_HX_phith_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_HX_phith_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_PX_costh_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_PX_costh_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_PX_phi_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_PX_phi_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/fit_PX_phith_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/fit_PX_phith_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/lph_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/lph_vs_lth_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/lphstar_vs_lthstar_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/lphstar_vs_lthstar_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/ltp_vs_lph_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/ltp_vs_lph_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/ltp_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/ltp_vs_lth_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf
mv Figures/ltilde_${nState}SUps_rap${rap_}_pT${pT_}.pdf Figures/ltilde_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.pdf



mv results_${nState}SUps_rap${rap_}_pT${pT_}.root results_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.root

if [ $nGen_ -eq 1 ]
then
cp results_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.root ${resultfilename}
fi

if [ $nGen_ -ge 2 ]
then
mv ${resultfilename} BUFFER_${resultfilename}
hadd -f ${resultfilename} BUFFER_${resultfilename} results_Fit${nGen_}_${nState}SUps_rap${rap_}_pT${pT_}.root
rm BUFFER_${resultfilename}
fi


nGen_=$((nGen_+1))
done

mv ${resultfilename} results_${nState}SUps_rap${rap_}_pT${pT_}.root

cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_rap${rap_}_pt${pT_}
./polGenRecFitPlot_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir
rm polGenRecFitPlot_rap${rap_}_pt${pT_}


pT_=$((pT_+1))
done
rap_=$((rap_+1))
done


done

rm *.so
rm *.d
rm polGenRecFitPlot

mkdir ../tmp
mv *.C ../tmp

