CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

#CXX+=-I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include
#LIBS+=-L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/lib

#DATA=

%.o : %.cc
	$(CXX) -c $<


all: runChibFit runChicFit minErrorToys plotToys

runChibFit: runChibFit.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

runChicFit: runChicFit.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

minErrorToys: minErrorToys.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

plotToys: plotToys.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
									
clean: 
	rm runChibFit runChicFit minErrorToys plotToys*.o
