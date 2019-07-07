//
// Neural Network
//

#include "neuralnetwork.h"
#include <fstream>
#include <sstream>

#include <iostream>
using namespace std;

/* ----------------------------------------------------------------------
   setup for class NeuralNetwork 
------------------------------------------------------------------------- */

NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize, int outputsSize)
{
    const int numberOfHiddenLayers = hiddenLayersSize.size();
    OpenNN::Vector<int> layersSize(numberOfHiddenLayers+2);

    layersSize[0] = inputsSize;  //input layer
    for (int i=0; i<numberOfHiddenLayers; i++)
            layersSize[i+1] = hiddenLayersSize[i]; // number of perceptrons at each hidden layer
    layersSize[numberOfHiddenLayers+1] = outputsSize;  //Output layer

    // for(int i=0; i< layersSize.size(); i++)
    //     cout << layersSize[i] << " ";
    // cout << endl;

    neuralNetwork.set(layersSize);
}

NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize): NeuralNetwork(inputsSize, hiddenLayersSize, 1) {}

const OpenNN::MultilayerPerceptron& NeuralNetwork::getNeuralNetwork() const { return neuralNetwork; }

int NeuralNetwork::getNumberOfInputs() const { return neuralNetwork.get_inputs_number(); }

int NeuralNetwork::getNumberOfOutputs() const { return neuralNetwork.get_outputs_number(); }

int NeuralNetwork::getNumberOfHiddenLayers() const { return neuralNetwork.get_layers_number()-1; } 