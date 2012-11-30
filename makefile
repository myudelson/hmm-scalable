CXX ?= $(CXX) $(CFLAGS)
CFLAGS = -Wall -Wconversion -O3 -fPIC
#LIBS = blas/blas.a
SHVER = 1
OS = $(shell uname)
#LIBS = -lblas

all: train predict

# lib: linear.o tron.o blas/blas.a
#	if [ "$(OS)" = "Darwin" ]; then \
#		SHARED_LIB_FLAG="-dynamiclib -Wl,-install_name,liblinear.so.$(SHVER)"; \
#	else \
#		SHARED_LIB_FLAG="-shared -Wl,-soname,liblinear.so.$(SHVER)"; \
#	fi; \
#	$(CXX) $${SHARED_LIB_FLAG} linear.o tron.o blas/blas.a -o liblinear.so.$(SHVER)


#train: tron.o linear.o utils.o StripedArray.o FitBit.o HMMProblem.o HMMProblemKT.o HMMProblemPiG.o HMMProblemPiGK.o HMMProblemAGK.o HMMProblemPiAGK.o trainhmm.cpp blas/blas.a
#	$(CXX) $(CFLAGS) -o trainhmm trainhmm.cpp tron.o linear.o utils.o HMMProblem.o HMMProblemKT.o HMMProblemPiG.o HMMProblemPiGK.o HMMProblemAGK.o HMMProblemPiAGK.o $(LIBS)

train: utils.o StripedArray.o FitBit.o HMMProblem.o HMMProblemPiGK.o HMMProblemAGK.o trainhmm.cpp
	$(CXX) $(CFLAGS) -o trainhmm trainhmm.cpp utils.o FitBit.o HMMProblem.o HMMProblemPiGK.o HMMProblemAGK.o


predict: predicthmm.o utils.o 
	$(CXX) $(CFLAGS) -o predicthmm predicthmm.o utils.o

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

FitBit.o: FitBit.cpp FitBit.h
	$(CXX) $(CFLAGS) -c -o FitBit.o FitBit.cpp

HMMProblem.o: HMMProblem.cpp HMMProblem.h
	$(CXX) $(CFLAGS) -c -o 	HMMProblem.o HMMProblem.cpp 
#HMMProblemKT.o: HMMProblemKT.cpp HMMProblemKT.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemKT.o HMMProblemKT.cpp
#HMMProblemPiG.o: HMMProblemPiG.cpp HMMProblemPiG.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemPiG.o HMMProblemPiG.cpp
HMMProblemPiGK.o: HMMProblemPiGK.cpp HMMProblemPiGK.h
	$(CXX) $(CFLAGS) -c -o HMMProblemPiGK.o HMMProblemPiGK.cpp
HMMProblemAGK.o: HMMProblemAGK.cpp HMMProblemAGK.h
	$(CXX) $(CFLAGS) -c -o HMMProblemAGK.o HMMProblemAGK.cpp
#HMMProblemPiAGK.o: HMMProblemPiAGK.cpp HMMProblemPiAGK.h
#	$(CXX) $(CFLAGS) -c -o HMMProblemPiAGK.o HMMProblemPiAGK.cpp

clean:
	make -C blas clean
	rm -f *.o trainhmm predicthmm