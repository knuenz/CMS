CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

%.o : %.cc
	$(CXX) -c $<

all: polRapPtPlot polGenRecFitPlot

polRapPtPlot: polRapPtPlot.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

polGenRecFitPlot: polGenRecFitPlot.cc
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
									
clean: 
	rm polRapPtPlot polGenRecFitPlot*.o
