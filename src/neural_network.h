/*
  neural_network.h: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Hossein Ghorbanfekr

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//
// Neural Network
//

#ifndef NNP_NEURALNETWORK_H
#define NNP_NEURALNETWORK_H

#include <vector>
#include "opennn.h"
// using namespace OpenNN;

class NeuralNetwork {
public:
    NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize, int outputsSize);
    NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize);
    ~NeuralNetwork();
    void setLayersActivationFunction(const std::vector<std::string>& activationFucntionsType);
    void readParameters(const std::string& fullPathFileName);
    int getNumberOfInputs();
    int getNumberOfOutputs();
    int getNumberOfLayers();
    int getNumberOfHiddenLayers();
    double calculateEnergy(double *descriptorValues, int descriptorSize);
    OpenNN::Vector<double> calculateJacobian(double *descriptorValues, int descriptorSize);

private:
    OpenNN::MultilayerPerceptron* multilayerPerceptron;
};


#endif //NNP_NEURALNETWORK_H
