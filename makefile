CC=g++ #//PAR
CXX=g++ #//PAR
CFLAGS = -Wall -Wconversion -O3 -fPIC -fopenmp #UNBOOST
#CFLAGS = -Wall -Wconversion -O3 -fPIC -fopenmp -I boost #BOOST
SHVER = 1
OS = $(shell uname)

#LIBS = blas/blas.a
#LIBS = -lblas

all: trainst train predict input

train: utils.o StripedArray.o FitBit.o FitBitSlicedA.o FitBitSlicedAB.o PULogistic.o PUCorbettianAdditive.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSlicedAB.o HMMProblemSlicedA.o HMMProblemComp.o HMMProblemSt.o InputUtil.o trainhmm.cpp
	$(CXX) $(CFLAGS) -o trainhmm trainhmm.cpp utils.o FitBit.o FitBitSlicedA.o FitBitSlicedAB.o InputUtil.o PULogistic.o PUCorbettianAdditive.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSlicedAB.o HMMProblemSlicedA.o StripedArray.o HMMProblemComp.o

trainst: utilsSt.o StripedArray.o FitBitSt.o HMMProblemSt.o InputUtilSt.o trainhmmst.cpp
	$(CXX) $(CFLAGS) -o trainhmmst trainhmmst.cpp utilsSt.o FitBitSt.o InputUtilSt.o HMMProblemSt.o 

predict: utils.o StripedArray.o FitBit.o FitBitSlicedA.o FitBitSlicedAB.o  PULogistic.o PUCorbettianAdditive.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSlicedAB.o HMMProblemSlicedA.o HMMProblemComp.o HMMProblemSt.o InputUtil.o predicthmm.cpp
	$(CXX) $(CFLAGS) -o predicthmm predicthmm.cpp utils.o FitBit.o FitBitSlicedA.o FitBitSlicedAB.o InputUtil.o PULogistic.o PUCorbettianAdditive.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSlicedAB.o HMMProblemSlicedA.o StripedArray.o HMMProblemComp.o

input: utils.o StripedArray.o InputUtil.o inputconvert.cpp
	$(CXX) $(CFLAGS) -o inputconvert inputconvert.cpp utils.o StripedArray.o InputUtil.o

utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) -c -o utils.o utils.cpp

utilsSt.o: utilsSt.cpp utilsSt.h
	$(CXX) $(CFLAGS) -c -o utilsSt.o utilsSt.cpp

StripedArray.o: StripedArray.cpp StripedArray.h
	$(CXX) $(CFLAGS) -c -o StripedArray.o StripedArray.cpp

InputUtil.o: InputUtil.cpp InputUtil.h
	$(CXX) $(CFLAGS) -c -o InputUtil.o InputUtil.cpp

InputUtilSt.o: InputUtilSt.cpp InputUtilSt.h
	$(CXX) $(CFLAGS) -c -o InputUtilSt.o InputUtilSt.cpp

FitBit.o: FitBit.cpp FitBit.h
	$(CXX) $(CFLAGS) -c -o FitBit.o FitBit.cpp

FitBitSt.o: FitBitSt.cpp FitBitSt.h
	$(CXX) $(CFLAGS) -c -o FitBitSt.o FitBitSt.cpp

FitBitSlicedA.o: FitBitSlicedA.cpp FitBitSlicedA.h
	$(CXX) $(CFLAGS) -c -o FitBitSlicedA.o FitBitSlicedA.cpp

FitBitSlicedAB.o: FitBitSlicedAB.cpp FitBitSlicedAB.h
	$(CXX) $(CFLAGS) -c -o FitBitSlicedAB.o FitBitSlicedAB.cpp

PULogistic.o: PULogistic.cpp PULogistic.h ParameterUnion.h
	$(CXX) $(CFLAGS) -c -o PULogistic.o PULogistic.cpp

PUCorbettianAdditive.o: PUCorbettianAdditive.cpp PUCorbettianAdditive.h ParameterUnion.h
	$(CXX) $(CFLAGS) -c -o PUCorbettianAdditive.o PUCorbettianAdditive.cpp

HMMProblem.o: HMMProblem.cpp HMMProblem.h
	$(CXX) $(CFLAGS) -c -o HMMProblem.o HMMProblem.cpp
HMMProblemPiGK.o: HMMProblemPiGK.cpp HMMProblemPiGK.h
	$(CXX) $(CFLAGS) -c -o HMMProblemPiGK.o HMMProblemPiGK.cpp
HMMProblemPiGKww.o: HMMProblemPiGKww.cpp HMMProblemPiGKww.h
	$(CXX) $(CFLAGS) -c -o HMMProblemPiGKww.o HMMProblemPiGKww.cpp
HMMProblemAGK.o: HMMProblemAGK.cpp HMMProblemAGK.h
	$(CXX) $(CFLAGS) -c -o HMMProblemAGK.o HMMProblemAGK.cpp
HMMProblemPiAGK.o: HMMProblemPiAGK.cpp HMMProblemPiAGK.h
	$(CXX) $(CFLAGS) -c -o HMMProblemPiAGK.o HMMProblemPiAGK.cpp
HMMProblemPiABGK.o: HMMProblemPiABGK.cpp HMMProblemPiABGK.h
	$(CXX) $(CFLAGS) -c -o HMMProblemPiABGK.o HMMProblemPiABGK.cpp
HMMProblemSlicedAB.o: HMMProblemSlicedAB.cpp HMMProblemSlicedAB.h
	$(CXX) $(CFLAGS) -c -o HMMProblemSlicedAB.o HMMProblemSlicedAB.cpp
HMMProblemSlicedA.o: HMMProblemSlicedA.cpp HMMProblemSlicedA.h
	$(CXX) $(CFLAGS) -c -o HMMProblemSlicedA.o HMMProblemSlicedA.cpp
HMMProblemComp.o: HMMProblemComp.cpp HMMProblemComp.h
	$(CXX) $(CFLAGS) -c -o HMMProblemComp.o HMMProblemComp.cpp
#HMMProblemEloK.o: HMMProblemEloK.cpp HMMProblemEloK.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemEloK.o HMMProblemEloK.cpp
#HMMProblemSt.o: HMMProblemSt.cpp HMMProblemSt.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemSt.o HMMProblemSt.cpp


clean:
#make -C blas clean
	rm -f *.o trainhmm trainhmmst predicthmm inputconvert

tidy:
	rm -f *.o
