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

NeuralNetwork::NeuralNetwork(const int& numberOfInputs, const std::vector<int>& hiddenLayersSize)
{
    const int numberOfHiddenLayers = hiddenLayersSize.size();
    OpenNN::Vector<int> layersSize(numberOfHiddenLayers+2);
    layersSize[0] = numberOfInputs;  //input layer
    for (int i=0; i<numberOfHiddenLayers; i++)
            layersSize[i+1] = hiddenLayersSize[i];
    layersSize[numberOfHiddenLayers+1] = 1;  //Output layer

    for(int i=0; i< layersSize.size(); i++)
        cout << layersSize[i] << " ";
    cout << endl;

    OpenNN::MultilayerPerceptron mlp(layersSize);
    cout << mlp.get_inputs_number() << endl;
    cout << mlp.get_layers_number() << endl;  //hidden layers + output layer 
    cout << mlp.get_outputs_number() << endl;

    OpenNN::Vector<double> input(1);
    input[0] = 0.5;
    OpenNN::Vector<double> output = mlp.calculate_outputs( input );
    cout << "output = " << output << endl;

    OpenNN::Matrix<double> jacobian = mlp.calculate_Jacobian ( input );
    cout << "jacobian = " << jacobian << endl;
}

NeuralNetwork::~NeuralNetwork() {}

// int NeuralNetwork::getNumberOfInputs() const { return numberOfInputs; }

// int NeuralNetwork::getNumberOfOutputs() const { return numberOfOutputs; }

// const std::vector<int>& NeuralNetwork::getNumberOfHiddenLayers() const { return numberOfHiddenLayers; }
