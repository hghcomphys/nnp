## NNP - C++ Implementation of Neural Network Potential

This repository is an experimental c++ implementation of high-dimensional neural network interatomic potential (HDNNP) used for force field molecular dynamics simulation. At this moment, only the _predict_ mode is implemented. 
This means that if one has a trained potential data files including symmetry functions and weights 
for the neural network, the excutable file `nnp-pred.x` will predict energy and force components of atom(s) for a given input structure file.

Dependencies:
- Open-source neural network library [OpenNN](http://www.opennn.net/documentation/)

How to compile:
```
make all
```

How to run:
```
./nnp-pred.x [path_to_nnp_potential] [input_file]
```
*Example: ./nnp-pred.x ../potentials/water ../example/input.data*
