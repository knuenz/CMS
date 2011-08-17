CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

CXX+=-I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include
LIBS+=-L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/lib

#DATA=

%.o : %.cc
	$(CXX) -c $<

all: polToyMCPlot

polToyMCPlot: polToyMCPlot.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
									
clean: 
	rm polToyMCPlot*.o
