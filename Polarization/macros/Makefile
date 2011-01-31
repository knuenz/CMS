CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

CXX+=-I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include
LIBS+=-L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/lib

FIT=AcceptanceMaps.o LifetimeModel.o MassModel.o PolarizationModel.o \
	CompositeModelBuilder.o

#DATA=

%.o : %.cc
	$(CXX) -c $<

all: polarizationFit GenpolarizationFit injectAcceptanceMaps extractFSRShape \
	extractPromptLifetimeShape extractBackgroundShapes \
	PlotMass PlotLifetime PlotPolarization \
	PlotFinalMass PlotFinalLifetime PlotParameters PlotAcceptanceMaps\
	extractNonPromptLifetimeShape CreateAccTxt PrepareFitTree FindConvergence\
	PlotNLL PlotSingleParameter


polarizationFit: polarizationFit.cc $(FIT) 
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

GenpolarizationFit: GenpolarizationFit.cc $(FIT) 
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractFSRShape: extractFSRShape.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

injectAcceptanceMaps: injectAcceptanceMaps.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractPromptLifetimeShape: extractPromptLifetimeShape.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractNonPromptLifetimeShape: extractNonPromptLifetimeShape.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

extractBackgroundShapes: extractBackgroundShapes.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

PlotLifetime: PlotLifetime.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
	
PlotMass: PlotMass.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
	
PlotPolarization: PlotPolarization.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit	
	
PlotFinalLifetime: PlotFinalLifetime.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

PlotFinalMass: PlotFinalMass.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
		
PlotParameters: PlotParameters.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
			
PlotAcceptanceMaps: PlotAcceptanceMaps.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

CreateAccTxt: CreateAccTxt.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
	
PrepareFitTree: PrepareFitTree.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
		
FindConvergence: FindConvegence.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
					
PlotNLL: PlotNLL.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
		
PlotSingleParameter: PlotSingleParameter.cc $(FIT)
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

							
clean: 
	rm polarizationFit GenpolarizationFit injectAcceptanceMaps extractFSRShape \
	extractPromptLifetimeShape extractBackgroundShapes \
	extractNonPromptLifetimeShape PlotMass PlotLifetime \
	PlotFinalMass PlotFinalLifetime PlotParameters\
	PlotPolarization PlotAcceptanceMaps CreateAccTxt PrepareFitTree FindConvergence\
	PlotNLL PlotSingleParameter*.o
