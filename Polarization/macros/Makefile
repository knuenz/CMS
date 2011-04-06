CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

CXX+=-I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include
LIBS+=-L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/lib

FIT=AcceptanceMaps.o LifetimeModel.o MassModel.o PolarizationModel.o \
	CompositeModelBuilder.o RooPolarizationConstraint.o RooPolarizationPdf.o 
#DATA=

%.o : %.cc
	$(CXX) -c $<


all: GenpolarizationFit injectAcceptanceMaps PlotAcceptanceMaps CreateAccTxt PlotSingleParameter ChangeBin PlotToyMC GenScript SysCheck ChangeAccEff

RooPolarizationDict.cc: RooPolarizationPdf.h RooPolarizationConstraint.h LinkDef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

GenpolarizationFit: GenpolarizationFit.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
	
injectAcceptanceMaps: injectAcceptanceMaps.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

PlotAcceptanceMaps: PlotAcceptanceMaps.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

CreateAccTxt: CreateAccTxt.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
					
PlotSingleParameter: PlotSingleParameter.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

ChangeBin: ChangeBin.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

PlotToyMC: PlotToyMC.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit	
	
GenScript: GenScript.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
	
SysCheck: SysCheck.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

ChangeAccEff: ChangeAccEff.cc $(FIT) RooPolarizationDict.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
									
clean: 
	rm GenpolarizationFit injectAcceptanceMaps PlotAcceptanceMaps CreateAccTxt PlotSingleParameter ChangeBin PlotToyMC GenScript SysCheck ChangeAccEff*.o
