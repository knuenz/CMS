This README file intends to explain the ToyMC framework
-------------------------------------------------------

Summary:
--------
-runToyMC.sh:
This script steers all generation, reconstruction fitting and the storage of the result
files. This script can be run in parallel for any scenario, bins, efficiencies, fiducial
cuts... Before strating toy-tests in parallel, wait at least 5 seconds for the previous
test to complete the compilation.

-PlotScript.sh:
When all fits, bins, scenarios are finished, or if one wants to get a snapshot of the
results in between, run this script. It plots the results that runToyMC.C produces (3
summary pdfs and one pdf for each bin showing all parameter and pull distributions).


The base directory 'basedir' contains following substructures:
latex: Latex files for visualization of plots and numerical results
interface: commonVar.h and rootIncludes.inc, containing the standard configuration of the
analysis (e.g. ranges of kinematic cells)
macros/ToyMC: All relevant files for the ToyMC tests that need to be changed to customize
the test are in this directory and are explained below.


ToyMC.h
-------
This file contains all relevant information that is needed for the tests. Polarization scenario
definitions, Number of events estimated from data, background fractions estimated from data,
Color and Marker settings.


polGen.C
--------
Macro needed to generate pseudo sample 
Output:
-genData.root, containing the pseudo sample on generator level
-GenResults.root, contaning information about the actual injected polarization in all frames,
needed for plotting the results


polRec.C
--------
Macro needed to alter pseudo sample according to defined efficiencies and fiducial cuts
Input:
-genData.root
Output:
-data.root, containing the eventually used sample
-efficiency.root, contaning 2D efficiency histograms (muon pT vs. |eta|)


polFit.C
--------
This macro conducts the fit 
Input:
-data.root
-efficiency.root
Output:
-results.root, containing the posterior distributions of all parameters, and a TTree containing the
numerical results (mean and rms of the posterior distributions)


polPlot.C
--------
This macro can plot the individually generated data sets vs. the fit result. It is not implemented in
the standard ToyMC generation and only needed for debugging purposes
Input:
-results.root
Output:
-several plots


polToyMCPlot.cc
---------------
This is the source file for the executable that plots the ToyMC results and produces latex tables,
steered by PlotScript.sh.

LoopToyMC.cc
---------------
This is the source file for the executable that generates, reconstructs and fits the data sets,
steered by runToyMC.sh.

Makefile
--------
This file steers the compilation of polToyMCPlot.cc and LoopToyMC.cc.


runToyMC.sh
----------
This is the steering script for generation, reconstruction and fit of the pseudo samples. It can be
started by 'sh runToyMC.sh'. In this file the customizable settings for the ToyMC tests can be
changed:
gen................ boolean, if true, polGen.C is executed
rec................ boolean, if true, polRec.C is executed
fit................ boolean, if true, polFit.C is executed
storagedir......... has to be set to the directory where all datasets results and figures will be
stored
basedir............ base directory of the code, containing the substructure 'latex', 'interface',
'macros/ToyMC'
JobID.............. Name to be specified for a certain test
rapBinMin.......... test will be conducted from this rapidity bin...
rapBinMax.......... ...to this rapidity bin
ptBinMin........... test will be conducted from this pT bin...
ptBinMax........... ...to this pT bin
polScenSig......... polarization scenario Signal (see ToyMC.h)
polScenBkg......... polarization scenario Background (see ToyMC.h)
frameSig........... natural polarization frame signal (1...CS, 2...HX, 3...PX)
frameBkg........... natural polarization frame background (1...CS, 2...HX, 3...PX)
nGenerations....... number of pseudo samples to be generated and fit
nEff............... defines the efficiency to be used (see polRec.C, 1...Eff=const=1, 2...oldEff, 3...newEff)
FidCuts............ defines, which set of cuts will be used (0...no cuts, 1... loose cuts, 2... tight cuts)
nSample............ numer of iterations in the algorithm (see polFit.C, 2000 burn-in iterations are included in
this number)
ConstEvents........ see below
UseConstEv......... boolean, if true it generates ConstEvents events, if false, it uses the number of events
stored in ToyMC.h
nSkipGen........... if the fit for one bin is split on several machines, this variable has to be used, in order
not to overwrite existing results 


PlotScript.h
------------
This is the script that steers the plotting of the ToyMC results (providing the executable polToyMCPlot
the correct inputs, executing latex to produce summary pdf files and pdf files containing tables of numerical
results)
The output is saved in the directory 'basedir/Figures/JobID/....pdf'. This script can be executed by
'sh PlotScript.sh'. Following settings need to be costomized:
storagedir......... same as above
basedir............ same as above
JobID.............. same as above
ptBinMin........... can be same as above, this just defines the minpT of the plots
ptBinMax........... can be same as above, this just defines the maxpT of the plots
frameSig........... same as above
polScenSig......... same as above
frameBkg........... same as above
polScenBkg......... same as above
nGenerations....... same as above
additionalName..... can be set optionally, if one wants to add a certain extention to the pdf summary files
