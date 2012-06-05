CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

%.o : %.cc
	$(CXX) -c $<

all: runMeanPt runData runMassFit runCopyTreeEntries PlotCosThetaPhiBG runTrimEventContent PlotCosThetaPhiDistribution runMassFitMC

runMeanPt: runMeanPt.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runData: runData.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runMassFit: runMassFit.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runMassFitMC: runMassFitMC.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runCopyTreeEntries: runCopyTreeEntries.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

PlotCosThetaPhiBG: PlotCosThetaPhiBG.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

runTrimEventContent: runTrimEventContent.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit

PlotCosThetaPhiDistribution: PlotCosThetaPhiDistribution.cc 
	$(CXX) $^ -o $@ $(LIBS) -lFoam -lMinuit
									
clean: 
	rm runMeanPt runData runMassFit runCopyTreeEntries PlotCosThetaPhiBG runTrimEventContent PlotCosThetaPhiDistribution runMassFitMC*.o
