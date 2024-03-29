/*
  neural_network.cpp: This file is part of Free Molecular Dynamics

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

#include "neural_network.h"
#include "logger.h"
#include <fstream>
#include <sstream>

/* ----------------------------------------------------------------------
   setup for class NeuralNetwork
------------------------------------------------------------------------- */
NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int> &hiddenLayersSize, int outputsSize)
{
    const int numberOfHiddenLayers = hiddenLayersSize.size();
    OpenNN::Vector<int> layersSize(numberOfHiddenLayers + 2);

    // set inputs, hidden, and output layers sizes
    layersSize[0] = inputsSize; // input layer
    for (int i = 0; i < numberOfHiddenLayers; i++)
        layersSize[i + 1] = hiddenLayersSize[i];        // number of perceptrons at each hidden layer
    layersSize[numberOfHiddenLayers + 1] = outputsSize; // Output layer

    // create multilayer perceptron
    multilayerPerceptron = new OpenNN::MultilayerPerceptron(layersSize);
}

NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int> &hiddenLayersSize) : NeuralNetwork(inputsSize, hiddenLayersSize, 1) {}

NeuralNetwork::~NeuralNetwork()
{
    delete multilayerPerceptron;
}

int NeuralNetwork::getNumberOfInputs()
{
    return multilayerPerceptron->get_inputs_number();
}

int NeuralNetwork::getNumberOfOutputs()
{
    return multilayerPerceptron->get_outputs_number();
}

int NeuralNetwork::getNumberOfLayers()
{
    return multilayerPerceptron->get_layers_number();
}

int NeuralNetwork::getNumberOfHiddenLayers()
{
    return getNumberOfLayers() - 1; /*exclude the output layer*/
}

void NeuralNetwork::setLayersActivationFunction(const std::vector<std::string> &activationFunctionsType)
{
    // set activation function from file into a list
    OpenNN::Vector<OpenNN::Perceptron::ActivationFunction> activationFunctions;
    for (auto each : activationFunctionsType)
    {
        if (each == "t")
            activationFunctions.push_back(OpenNN::Perceptron::ActivationFunction::HyperbolicTangent);
        else if (each == "l")
            activationFunctions.push_back(OpenNN::Perceptron::ActivationFunction::Linear);
        else
            throw runtime_error((Log(ERROR) << "Unknown activation function type " + each).toString());
    }
    if (activationFunctions.size() != getNumberOfLayers())
        throw runtime_error((Log(ERROR) << "Inconsistent number of given activation functions").toString());

    // set activation functions for layers
    multilayerPerceptron->set_layers_activation_function(activationFunctions);
}

void NeuralNetwork::readParameters(const std::string &filename)
{
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error((Log(ERROR) << "Unable to open script file " + filename).toString());

    OpenNN::Vector<OpenNN::Matrix<double>> layers_synaptic_weights;
    OpenNN::Vector<OpenNN::Vector<double>> layers_biases;
    for (int i = 0; i < getNumberOfLayers(); i++)
    {
        const int nrow = multilayerPerceptron->get_layer(i).get_perceptrons_number();
        const int ncol = multilayerPerceptron->get_layers_inputs_number()[i];

        // synaptic weights
        layers_synaptic_weights.push_back(OpenNN::Matrix<double>(nrow, ncol, 0.0));
        Log(INFO) << "Weight matrix in layer(" << i + 1 << "): (" << nrow << ", " << ncol << ")";

        // biases
        layers_biases.push_back(OpenNN::Vector<double>(nrow, 0.0));
        Log(INFO) << "Biases vector in layer(" << i + 1 << "): (" << nrow << ")";
    }

    char cSpace = ' ';
    std::string line, dummy;

    while (std::getline(inFile, line))
    {
        std::stringstream ss(line);
        std::string sIndvStr;
        while (std::getline(ss, sIndvStr, cSpace))
        {

            if (sIndvStr == "#")
                continue; // skip comments

            double weight;
            std::string weightType;
            ss >> weight >> weightType;
            // Log(DEBUG) << weight << " " << weightType;

            if (weightType == "a")
            {
                int index, dummy, inputIndex, layerIndex, neuronIndex;
                ss >> index >> dummy >> inputIndex >> layerIndex >> neuronIndex;
                // Log(DEBUG) << weight << " a " << inputIndex << " " << layerIndex << " " << neuronIndex;

                layers_synaptic_weights[layerIndex - 1](neuronIndex - 1, inputIndex - 1) = weight;
                // Log(DEBUG) << layers_synaptic_weights[layerIndex-1][neuronIndex-1, inputIndex-1] << " " << weight;
            }
            else if (weightType == "b")
            {
                int index, layerIndex, neuronIndex;
                ss >> index >> layerIndex >> neuronIndex;
                // Log(DEBUG) << weight << " b " << layerIndex << " " << neuronIndex;

                layers_biases[layerIndex - 1][neuronIndex - 1] = weight;
                // Log(DEBUG) << layers_biases[layerIndex-1][neuronIndex-1] << " " << weight;
            }
        }
    }
    inFile.close();

    // set all weights and biases
    multilayerPerceptron->set_layers_synaptic_weights(layers_synaptic_weights);
    multilayerPerceptron->set_layers_biases(layers_biases);
}

double NeuralNetwork::calculateEnergy(double *descriptorValues, int descriptorSize)
{
    const int numberOfInputs = getNumberOfInputs();
    if (descriptorSize != numberOfInputs)
        throw std::runtime_error(
            (Log(ERROR) << "Unexpected size of inputs (" << numberOfInputs << "), expected "
                        << descriptorSize << "(energy calculation)")
                .toString());

    // TODO: optimize the conversion
    OpenNN::Vector<double> inputs(numberOfInputs);
    for (int i = 0; i < numberOfInputs; i++)
        inputs[i] = descriptorValues[i];

    return multilayerPerceptron->calculate_outputs(inputs)[0]; // last layer has one output perceptron
}

OpenNN::Vector<double> NeuralNetwork::calculateJacobian(double *descriptorValues, int descriptorSize)
{
    const int numberOfInputs = getNumberOfInputs();
    if (descriptorSize != numberOfInputs)
        throw std::runtime_error(
            (Log(ERROR) << "Unexpected size of inputs (" << numberOfInputs << "), expected "
                        << descriptorSize << "(force calculation)")
                .toString());

    // TODO: optimize the conversion
    OpenNN::Vector<double> inputs(getNumberOfInputs(), 1);
    for (int i = 0; i < getNumberOfInputs(); i++)
        inputs[i] = descriptorValues[i];

    return multilayerPerceptron->calculate_Jacobian(inputs).to_vector();
}
