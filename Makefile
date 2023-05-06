CXX = g++
CXXFLAGS= -O3 -Wall -g #-pg -g #-Wall #-O3
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = main.o #BaseReads.o Alignment.o 

#asan=1
ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

all: trust4 bam-extractor fastq-extractor annotator

trust4: main.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

bam-extractor: BamExtractor.o
	if [ ! -f ./samtools-0.1.19/libbam.a ] ; \
	        then \
		                cd samtools-0.1.19 ; make ;\
	fi ;
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS) -lbam

fastq-extractor: FastqExtractor.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

annotator: Annotator.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

main.o: main.cpp AlignAlgo.hpp ReadFiles.hpp kseq.h SeqSet.hpp KmerIndex.hpp SimpleVector.hpp defs.h KmerCode.hpp KmerCount.hpp
BamExtractor.o: BamExtractor.cpp alignments.hpp defs.h SeqSet.hpp
FastqExtractor.o: FastqExtractor.cpp ReadFiles.hpp defs.h SeqSet.hpp BarcodeCorrector.hpp SimpleVector.hpp ReadFormatter.hpp BarcodeTranslator.hpp
Annotator.o: Annotator.cpp AlignAlgo.hpp ReadFiles.hpp kseq.h SeqSet.hpp KmerIndex.hpp SimpleVector.hpp defs.h KmerCode.hpp KmerCount.hpp FileLineReader.hpp
#Alignment.o: Alignment.cpp Alignment.h SimpleVector.h defs.h StatsTests.h KmerTree.h ReadSet.h KmerIndex.h poa.h

clean:
	rm -f *.o *.gch trust4 bam-extractor annotator fastq-extractor
