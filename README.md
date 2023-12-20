## NNP - A C++ Implementation of Neural Network Potential

This repository presents an experimental implementation in C++ of a high-dimensional neural network interatomic potential (HDNNP) used in force field molecular dynamics simulations. 
Currently, only the _predict_ mode has been implemented. This means that if we possess trained potential data files containing symmetry functions and weights for the neural network, the executable file `nnp-pred.x` will be able to predict the energy and force components of atom(s) for a provided input structure file.

This code base is intended solely for educational purposes. 
Please see [RuNNer](https://www.uni-goettingen.de/de/560580.html) and 
[N2P2](https://github.com/CompPhysVienna/n2p2) packages for developing HDNNPs.


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

### Examples 
Here we predict the total energy and force components for 222 water molecules in the simulation box. 
```bash
$ nnp-pred.x ../potentials/water ../examples/input.data

[INFO ] NNP directory: ../potentials/water//
[INFO ] Input: ../examples/input.data
[INFO ] Output: ../examples/input.data.pred
[INFO ] Read ../examples/input.data (222 atoms)
[INFO ] Element H: 148
[INFO ] Element O: 74
[INFO ] Cell (PBC)
[INFO ] Read ../potentials/water//input.nn
[INFO ] number_of_elements: 2
[INFO ] Element: H
[INFO ] Element: O
[INFO ] ACSF(H)
[INFO ] ACSF(O)
[INFO ] global_hidden_layers_short: 2
[INFO ] global_nodes_short: 25
[INFO ] global_nodes_short: 25
[INFO ] global_activation_short: t
[INFO ] global_activation_short: t
[INFO ] global_activation_short: l
[INFO ] Read ../potentials/water//nnp-scaling.log.0000 (H)
[INFO ] Read ../potentials/water//nnp-scaling.log.0000 (O)
[INFO ] Read scaling.data (H)
[INFO ] Read scaling.data (O)
[INFO ] Neural Network (H)
[INFO ] Neural Network (O)
[INFO ] Read ../potentials/water//weights.001.data
[INFO ] Weight matrix in layer(1): (25, 27)
[INFO ] Biases vector in layer(1): (25)
[INFO ] Weight matrix in layer(2): (25, 25)
[INFO ] Biases vector in layer(2): (25)
[INFO ] Weight matrix in layer(3): (1, 25)
[INFO ] Biases vector in layer(3): (1)
[INFO ] Initialize weights for Neural network (H)
[INFO ] Read ../potentials/water//weights.008.data
[INFO ] Weight matrix in layer(1): (25, 30)
[INFO ] Biases vector in layer(1): (25)
[INFO ] Weight matrix in layer(2): (25, 25)
[INFO ] Biases vector in layer(2): (25)
[INFO ] Weight matrix in layer(3): (1, 25)
[INFO ] Biases vector in layer(3): (1)
[INFO ] Initialize weights for Neural network (O)
[INFO ] Calculate NNP energy
[INFO ] Calculating NNP forces for entire structure ... (it may take a while)
[INFO ] NNP force calculation was successfully done
[INFO ] Write structure into ../examples/input.data.pred
```
