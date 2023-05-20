## NNP - C++ Implementation of Neural Network Potential

This repository presents an experimental implementation in C++ of a high-dimensional neural network interatomic potential (HDNNP) used in force field molecular dynamics simulations. 
Currently, only the _predict_ mode has been implemented. This means that if we possess trained potential data files containing symmetry functions and weights for the neural network, the executable file `nnp-pred.x` will be able to predict the energy and force components of atom(s) for a provided input structure file.

This repository is intended solely for educational purposes. 

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
