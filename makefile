
CXX = g++
CXXFLAGS = -std=c++11 #-Wall -O3 #-D__PURE_INTEL_C99_HEADERS__ 
LFLAGS = -lopennn -ltinyxml2 #-fopenmp #-ltensorflow

OBJS = main.o acsf.o atoms.o symmfunc.o cutofffunc.o nnp.o neuralnetwork.o symmfuncscaler.o logger.o

all: $(OBJS)
	$(CXX) $(OBJS) -o nnp.x $(LFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

acsf.o: acsf.cpp acsf.h
	$(CXX) $(CXXFLAGS) -c acsf.cpp

atoms.o: atoms.cpp atoms.h
	$(CXX) $(CXXFLAGS) -c atoms.cpp

symmfunc.o: symmfunc.cpp symmfunc.h
	$(CXX) $(CXXFLAGS) -c symmfunc.cpp

cutofffunc.o: cutofffunc.cpp cutofffunc.h
	$(CXX) $(CXXFLAGS) -c cutofffunc.cpp

nnp.o: nnp.cpp nnp.h
	$(CXX) $(CXXFLAGS) -c nnp.cpp

neuralnetwork.o: neuralnetwork.cpp neuralnetwork.h
	$(CXX) $(CXXFLAGS) -c neuralnetwork.cpp

symmfuncscaler.o: symmfuncscaler.cpp symmfuncscaler.h
	$(CXX) $(CXXFLAGS) -c symmfuncscaler.cpp

logger.o: logger.cpp logger.h
	$(CXX) $(CXXFLAGS) -c logger.cpp

clean:
	rm -f *.o nnp.x
