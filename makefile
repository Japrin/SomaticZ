SAMTOOLS_DIR=	/home/zhoudonger/01bin/repos/samtools/samtools-0.1.6
CXXFLAGS=		-I$(SAMTOOLS_DIR) -O3 -Wall -c
CXX=			g++
OBJS=			somatic_detector.o sniper_pileup.o main.o
LIBPATH=		$(SAMTOOLS_DIR)

all: MySomaticDetector

MySomaticDetector: $(OBJS)
	$(CXX) $(OBJS) -o MySomaticDetector -L$(LIBPATH) -lz -lbam
sniper_pileup.o: sniper_pileup.h
	$(CXX) $(CXXFLAGS) sniper_pileup.cpp -o sniper_pileup.o -L$(LIBPATH) -lz -lbam
somatic_detector.o: somatic_detector.h
	$(CXX) $(CXXFLAGS) somatic_detector.cpp -o somatic_detector.o -L$(LIBPATH) -lz -lbam
main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o main.o -L$(LIBPATH) -lz -lbam
clean:
	rm *.o MySomaticDetector
