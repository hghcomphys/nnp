# C++ Implementation of Neural Network Potential

This repository is a c++ implementation of high-dimensional neural network interatomic potential used for classical molecular dynamics simulation. At this moment, only predict mode is implemented. 
This means that if a trained potential data including symmetry functions and weights 
for the neural network are provided, this repo can predict total energy and force components 
on each atom for given atomic structure.

Dependencies:
- Open-source neural network library [OpenNN](http://www.opennn.net/documentation/)

How to compile:
```
make all
```