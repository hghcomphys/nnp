## NNP - C++ Implementation of Neural Network Potential

This repository presents an implementation in C++ of a high-dimensional neural network interatomic potential (HDNNP) used in force field molecular dynamics simulations. 
Currently, only the _predict_ mode has been implemented. This means that if we possess trained potential data files containing symmetry functions and weights for the neural network, the executable file `nnp-pred.x` will be able to predict the energy and force components of atom(s) for a provided input structure file.

This repository is intended solely for educational purposes. 

### How to Build
#### Dependencies
- Neural network library [OpenNN](https://github.com/hghcomphys/OpenNN) (forked repo)
```bash
$ git clone git@github.com:hghcomphys/OpenNN.git ~/opennn && cd ~/opennn
$ cmake .
$ make
```

- C++ XML parser [TinyXML-2](https://github.com/leethomason/tinyxml2)
```bash
$ sudo apt install libtinyxml2-dev
```

#### NNP code
```bash
$ git git@github.com:hghcomphys/nnp.git ~/nnp && cd ~/nnp/src
$ make 
```

### Example 
Here we predict the total energy and force components for water molecules in the simulation box. 
```bash
$ nnp-pred.x ../potentials/water ../example/input.data
```
