## C++ Implementation of Neural Network Potential

This repository is a c++ implementation of high-dimensional neural network interatomic potential (NNP) used for force field molecular dynamics simulation. At this moment, only predict mode is implemented. 
This means that if one has trained potential data files including symmetry functions and weights 
for the neural network, current repo can predict energy and force components 
for atom(s) for given atomic structure.

Dependencies:
- Open-source neural network library [OpenNN](http://www.opennn.net/documentation/)

How to compile:
```
$ make all
```

How to run:
```
$ ./nnp-pred.x path_to_nnp_potential_directory
```
*Example: ./nnp-pred.x potentials/water*