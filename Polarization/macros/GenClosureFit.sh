#!/bin/sh

#INITIAL SETTINGS
cutt=60
TREENAME=AccEffCutUnRe${cutt}
JOBNAME=AccEffCutUnRe${cutt}
dirstruct=/scratch/knuenz/Polarization/RootInput/
mkdir ${dirstruct}ProjectClosure/${JOBNAME}
cp polarizationFit.py polarizationFit${JOBNAME}.py
cp polarizationFitSimple.py polarizationFitSimple${JOBNAME}.py

#SETTINGS:::MAPS
geomAccPR=${dirstruct}geomAccHistos_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root
recoEffPR=${dirstruct}recoEffHistos_ATLASPT_20March2011_phiFolded_zeroBinsCorrected.root
#trigEffPR=${dirstruct}trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root
geomAccNP=${dirstruct}geomAccHistos_NP_WithFSR_PTandEtaSmeared_uniform_ATLASPT_PythiaRap_18March2011_merged_phiFolded_zeroBinsCorrected.root
recoEffNP=${dirstruct}recoEffHistos_NP_ATLASPT_19March2011_phiFolded_zeroBinsCorrected.root
#trigEffNP=${dirstruct}trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root

trigEffPR=${dirstruct}AccEffCutUnRe${cutt}_trigEffHistos_ATLASPT_DoubleMu0_20March2011_phiFolded_zeroBinsCorrected.root
trigEffNP=${dirstruct}AccEffCutUnRe${cutt}_trigEffHistos_NP_ATLASPT_DoubleMu0_19March2011_phiFolded_zeroBinsCorrected.root

#geomAccPR=${dirstruct}2_7_PRgeomAccATLASpT24Mar.root
#recoEffPR=${dirstruct}2_7_PRrecoEffATLASpT24Mar.root
#trigEffPR=${dirstruct}2_7_PRtrigEffATLASpT24Mar.root
#geomAccNP=${dirstruct}2_7_NPgeomAccATLASpT24Mar.root
#recoEffNP=${dirstruct}2_7_NPrecoEffATLASpT24Mar.root
#trigEffNP=${dirstruct}2_7_NPtrigEffATLASpT24Mar.root

#geomAccPR=${dirstruct}flatHisto.root
#recoEffPR=${dirstruct}flatHisto.root
#trigEffPR=${dirstruct}flatHisto.root
#geomAccNP=${dirstruct}flatHisto.root
#recoEffNP=${dirstruct}flatHisto.root
#trigEffNP=${dirstruct}flatHisto.root

#SETTINGS::BIN
generations=200
for rap_ in 1;do
for pT_ in 6 7 8;do

if [ $pT_ -eq 8 ] && [ $rap_ -eq 1 ]; then
continue
fi

for scenario in 4;do
mkdir ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}
cp ${dirstruct}ProjectClosure/saveTrees/TTree_GEN${generations}_${TREENAME}_scen${scenario}_rap${rap_}_pT${pT_}.root ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/TTree_GEN${generations}_${TREENAME}_scen${scenario}_rap${rap_}_pT${pT_}.root
cp polarizationFit${JOBNAME}.py polarizationFit${JOBNAME}_rap${rap_}_pt${pT_}_scen${scenario}.py
cp polarizationFitSimple${JOBNAME}.py polarizationFitSimple${JOBNAME}_${rap_}_pt${pT_}_scen${scenario}.py
rm RooFitResult_${JOBNAME}_scen${scenario}_rap${rap_}_pt${pT_}.root

for generation in $(seq ${generations});do #1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50;do

python polarizationFit${JOBNAME}_rap${rap_}_pt${pT_}_scen${scenario}.py --workspaceName=CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/TTree_GEN${generations}_${TREENAME}_scen${scenario}_rap${rap_}_pT${pT_}.root --fitFrame=CS --testBin=${rap_},${pT_} --noBackground --gen=${generation}
python polarizationFit${JOBNAME}_rap${rap_}_pt${pT_}_scen${scenario}.py --workspaceName=CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/TTree_GEN${generations}_${TREENAME}_scen${scenario}_rap${rap_}_pT${pT_}.root --fitFrame=HX --testBin=${rap_},${pT_} --noBackground --gen=${generation}

python polarizationFitSimple${JOBNAME}_${rap_}_pt${pT_}_scen${scenario}.py --workspaceName=CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS --fitFrame=CS --testBin=${rap_},${pT_} --acceptanceMap=${geomAccPR},${geomAccNP} --recoEfficiencyMap=${recoEffPR},${recoEffNP} --trigEfficiencyMap=${trigEffPR},${trigEffNP} --noBackground --BfracTruth --lambdaPhiSub=lambda_tilde --gen=${generation} --scenario=${scenario}  --JOBNAME=${JOBNAME} | tee terminal_CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS.txt
python polarizationFitSimple${JOBNAME}_${rap_}_pt${pT_}_scen${scenario}.py --workspaceName=CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX --fitFrame=HX --testBin=${rap_},${pT_} --acceptanceMap=${geomAccPR},${geomAccNP} --recoEfficiencyMap=${recoEffPR},${recoEffNP} --trigEfficiencyMap=${trigEffPR},${trigEffNP} --noBackground --BfracTruth --lambdaPhiSub=lambda_tilde --gen=${generation} --scenario=${scenario}  --JOBNAME=${JOBNAME} | tee terminal_CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX.txt

mv CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS.root ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/CLOSURE_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS.root
mv CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX.root ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/CLOSURE_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX.root
mv CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS-CS.root ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/CLOSURE_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS-CS.root
mv CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX-HX.root ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/CLOSURE_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX-HX.root
mv terminal_CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS.txt ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/terminal_CLOSURE_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_CS.txt
mv terminal_CLOSURE_${JOBNAME}_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX.txt ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/terminal_CLOSURE_scen${scenario}_rap${rap_}_pT${pT_}_gen${generation}_HX.txt

done

rm polarizationFit${JOBNAME}_rap${rap_}_pt${pT_}_scen${scenario}.py
rm polarizationFitSimple${JOBNAME}_${rap_}_pt${pT_}_scen${scenario}.py
cp RooFitResult_${JOBNAME}_scen${scenario}_rap${rap_}_pt${pT_}.root ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/RooFitResult_${JOBNAME}_scen${scenario}_rap${rap_}_pt${pT_}.root
mkdir ${dirstruct}ProjectClosure/PlotResults/${JOBNAME}
cp ${dirstruct}ProjectClosure/${JOBNAME}/scenario${scenario}/RooFitResult_${JOBNAME}_scen${scenario}_rap${rap_}_pt${pT_}.root ${dirstruct}ProjectClosure/PlotResults/${JOBNAME}/RooFitResult_${JOBNAME}_scen${scenario}_rap${rap_}_pt${pT_}.root

done
done
done
