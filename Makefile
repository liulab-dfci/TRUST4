CXX = g++
CXXFLAGS= -g #-Wall #-O3
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = main.o #BaseReads.o Alignment.o 

all: bcr

bcr: $(OBJECTS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) $(LINKFLAGS)

main.o: main.cpp AlignAlgo.hpp ReadFiles.hpp kseq.h SeqSet.hpp KmerIndex.hpp SimpleVector.hpp defs.h StatsTests.h poa.hpp
#Alignment.o: Alignment.cpp Alignment.h SimpleVector.h defs.h StatsTests.h KmerTree.h ReadSet.h KmerIndex.h poa.h

clean:
	rm -f *.o *.gch lec
