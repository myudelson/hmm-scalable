CC=g++ #//PAR
CXX=g++ #//PAR
#CXX ?= $(CXX) $(CFLAGS)
CFLAGS = -Wall -Wconversion -O3 -fPIC -fopenmp #UNBOOST
#CFLAGS = -Wall -Wconversion -O3 -fPIC -fopenmp -I boost #BOOST
SHVER = 1
OS = $(shell uname)

#LIBS = blas/blas.a
#LIBS = -lblas

all: train predict input tidy

# lib: linear.o tron.o blas/blas.a
#	if [ "$(OS)" = "Darwin" ]; then \
#		SHARED_LIB_FLAG="-dynamiclib -Wl,-install_name,liblinear.so.$(SHVER)"; \
#	else \
#		SHARED_LIB_FLAG="-shared -Wl,-soname,liblinear.so.$(SHVER)"; \
#	fi; \
#	$(CXX) $${SHARED_LIB_FLAG} linear.o tron.o blas/blas.a -o liblinear.so.$(SHVER)


#train: tron.o linear.o utils.o StripedArray.o FitBit.o HMMProblem.o HMMProblemKT.o HMMProblemPiG.o HMMProblemPiGK.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o trainhmm.cpp blas/blas.a
#	$(CXX) $(CFLAGS) -o trainhmm trainhmm.cpp tron.o linear.o utils.o HMMProblem.o HMMProblemKT.o HMMProblemPiG.o HMMProblemPiGK.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o $(LIBS)

train: utils.o StripedArray.o FitBit.o FitBitSliced.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSliced.o InputUtil.o trainhmm.cpp
	$(CXX) $(CFLAGS) -o trainhmm trainhmm.cpp utils.o FitBit.o FitBitSliced.o InputUtil.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSliced.o StripedArray.o

predict: utils.o StripedArray.o FitBit.o FitBitSliced.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSliced.o InputUtil.o predicthmm.cpp
	$(CXX) $(CFLAGS) -o predicthmm predicthmm.cpp utils.o FitBit.o FitBitSliced.o InputUtil.o HMMProblem.o HMMProblemPiGK.o HMMProblemPiGKww.o HMMProblemAGK.o HMMProblemPiAGK.o HMMProblemPiABGK.o HMMProblemSliced.o StripedArray.o

input: utils.o StripedArray.o InputUtil.o inputconvert.cpp
	$(CXX) $(CFLAGS) -o inputconvert inputconvert.cpp utils.o StripedArray.o InputUtil.o

#tron.o: tron.cpp tron.h
#	$(CXX) $(CFLAGS) -c -o tron.o tron.cpp

#linear.o: linear.cpp linear.h
#	$(CXX) $(CFLAGS) -c -o linear.o linear.cpp

#blas/blas.a: blas/*.c blas/*.h
#	make -C blas OPTFLAGS='$(CFLAGS)' CC='$(CC)';

utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) -c -o utils.o utils.cpp

StripedArray.o: StripedArray.cpp StripedArray.h
	$(CXX) $(CFLAGS) -c -o StripedArray.o StripedArray.cpp

InputUtil.o: InputUtil.cpp InputUtil.h
	$(CXX) $(CFLAGS) -c -o InputUtil.o InputUtil.cpp

FitBit.o: FitBit.cpp FitBit.h
	$(CXX) $(CFLAGS) -c -o FitBit.o FitBit.cpp

FitBitSliced.o: FitBitSliced.cpp FitBitSliced.h
	$(CXX) $(CFLAGS) -c -o FitBitSliced.o FitBitSliced.cpp

HMMProblem.o: HMMProblem.cpp HMMProblem.h
	$(CXX) $(CFLAGS) -c -o HMMProblem.o HMMProblem.cpp 
#HMMProblemKT.o: HMMProblemKT.cpp HMMProblemKT.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemKT.o HMMProblemKT.cpp
#HMMProblemPiG.o: HMMProblemPiG.cpp HMMProblemPiG.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemPiG.o HMMProblemPiG.cpp
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
HMMProblemSliced.o: HMMProblemSliced.cpp HMMProblemSliced.h
	$(CXX) $(CFLAGS) -c -o HMMProblemSliced.o HMMProblemSliced.cpp

clean:
	#make -C blas clean
	rm -f *.o trainhmm predicthmm inputconvert

tidy:
	rm -f *.o

