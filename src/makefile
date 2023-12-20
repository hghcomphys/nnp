CXX = g++
CXXFLAGS = -Wall -O2 -I ~/opennn/opennn #-D__PURE_INTEL_C99_HEADERS__ # -std=c++11
LFLAGS = -lopennn -L ~/opennn/opennn #-ltinyxml2 #-pg #-fopenmp #-ltensorflow

EXE = nnp-pred.x

OBJS = nnp-pred.o atom.o distance.o descriptor.o structure.o symmetry_function.o \
	cutoff_function.o neural_network_potential.o neural_network.o symmetry_function_scaler.o logger.o

DEPS = atom.h descriptor.h structure.h symmetry_function.h cutoff_function.h \
	neural_network_potential.h neural_network.h symmetry_function_scaler.h logger.h distance.h

all: $(OBJS)
	$(CXX) $(OBJS) -o $(EXE) $(LFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f $(OBJS) $(EXE)

