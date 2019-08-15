
CXX = g++
CXXFLAGS = #-Wall #-pg  #-O3 #-D__PURE_INTEL_C99_HEADERS__  # -std=c++11
LFLAGS = -lopennn -ltinyxml2 #-pg #-fopenmp #-ltensorflow

OBJS = main.o descriptor.o structure.o symmetry_function.o cutoff_function.o neural_network_potential.o neural_network.o symmetry_function_scaler.o log.o

all: $(OBJS)
	$(CXX) $(OBJS) -o nnp.x $(LFLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

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

log.o: log.cpp log.h
	$(CXX) $(CXXFLAGS) -c log.cpp

clean:
	rm -f *.o nnp.x
