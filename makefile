
CXX = g++
CXXFLAGS = -Wall -O2 #-D__PURE_INTEL_C99_HEADERS__  # -std=c++11
LFLAGS = -lopennn -ltinyxml2 #-pg #-fopenmp #-ltensorflow

EXE = nnp-pred.x
OBJS = nnp-pred.o atom.o distance.o descriptor.o structure.o symmetry_function.o cutoff_function.o neural_network_potential.o neural_network.o symmetry_function_scaler.o logger.o

all: $(OBJS)
	$(CXX) $(OBJS) -o $(EXE) $(LFLAGS)

nnp-pred.o: nnp-pred.cpp
	$(CXX) $(CXXFLAGS) -c nnp-pred.cpp

atom.o: atom.cpp atom.h
	$(CXX) $(CXXFLAGS) -c atom.cpp

descriptor.o: descriptor.cpp descriptor.h
	$(CXX) $(CXXFLAGS) -c descriptor.cpp

structure.o: structure.cpp structure.h
	$(CXX) $(CXXFLAGS) -c structure.cpp

symmetry_function.o: symmetry_function.cpp symmetry_function.h
	$(CXX) $(CXXFLAGS) -c symmetry_function.cpp

cutoff_function.o: cutoff_function.cpp cutoff_function.h
	$(CXX) $(CXXFLAGS) -c cutoff_function.cpp

neural_network_potential.o: neural_network_potential.cpp neural_network_potential.h
	$(CXX) $(CXXFLAGS) -c neural_network_potential.cpp

neural_network.o: neural_network.cpp neural_network.h
	$(CXX) $(CXXFLAGS) -c neural_network.cpp

symmetry_function_scaler.o: symmetry_function_scaler.cpp symmetry_function_scaler.h
	$(CXX) $(CXXFLAGS) -c symmetry_function_scaler.cpp

logger.o: logger.cpp logger.h
	$(CXX) $(CXXFLAGS) -c logger.cpp

distance.o: distance.cpp distance.h
	$(CXX) $(CXXFLAGS) -c distance.cpp

clean:
	rm -f *.o $(EXE)
