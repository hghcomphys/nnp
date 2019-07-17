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

NeuralNetwork::~NeuralNetwork() {}

const OpenNN::MultilayerPerceptron& NeuralNetwork::getNeuralNetwork() const { return neuralNetwork; }

int NeuralNetwork::getNumberOfInputs() const { return neuralNetwork.get_inputs_number(); }

int NeuralNetwork::getNumberOfOutputs() const { return neuralNetwork.get_outputs_number(); }

int NeuralNetwork::getNumberOfLayers() const { return neuralNetwork.get_layers_number(); }

int NeuralNetwork::getNumberOfHiddenLayers() const { return getNumberOfLayers()-1; /*exclude the output layer*/} 

void NeuralNetwork::setLayersActivationFunction(const std::vector<std::string>& activationFucntionsType) 
{
    OpenNN::Vector<OpenNN::Perceptron::ActivationFunction> activationFunctions;
    for (auto each: activationFucntionsType) {
        if( each == "t" )
             activationFunctions.push_back( OpenNN::Perceptron::ActivationFunction::HyperbolicTangent );
        else if ( each == "l" )
            activationFunctions.push_back( OpenNN::Perceptron::ActivationFunction::Linear );
        else
            throw runtime_error("Unknown activation function type " + each);
    }
    if( activationFunctions.size() != getNumberOfLayers() )
        throw runtime_error("Inconsistent number of given activation functions");
}

void NeuralNetwork::readParameters(const std::string& filename) 
{
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error("Unable to open script file " + filename);

    OpenNN::Vector<OpenNN::Matrix<double>> layers_synaptic_weights;
    OpenNN::Vector<OpenNN::Vector< double >>  layers_biases;
    for (int i=0; i<getNumberOfLayers(); i++) {

        const int nrow = neuralNetwork.get_layer(i).get_perceptrons_number();
        const int ncol = neuralNetwork.get_layers_inputs_number()[i];
        
        // synamptic weights
        layers_synaptic_weights.push_back( OpenNN::Matrix<double>(nrow, ncol) );
        cout << "weight matrix in layer(" << i+1 << "): (" << nrow << ", " << ncol << ")" << endl;

        // biases
        layers_biases.push_back( OpenNN::Vector<double>(nrow) );
        cout << "biases vector in layer(" << i+1 << "): (" << nrow << ")" << endl;
    }

    char cSpace = ' ';
    std::string line, dummy;
    double ddummy;
    int idummy;

    while ( std::getline(inFile, line) ) {
        std::stringstream ss(line);
        std::string sIndvStr;
        while ( std::getline(ss, sIndvStr, cSpace) ) {

            if (sIndvStr == "#")
                continue; // skip comments

            double weight;
            std::string weightType;
            ss >> weight >> weightType;

            if (weightType == "a") {

                int index, l_s, n_s, l_e, n_e;
                ss >> index >> l_s >> n_s >> l_e >> n_e;
                // std::cout << weight << " a " << index << " " << l_s << " " << n_s << " " << l_e << " " << n_e << std::endl;

                layers_synaptic_weights[l_s][n_s-1, n_e-1] = weight;
                // cout << layers_synaptic_weights[l_s][n_s, n_e] << " " << weight << endl;
            }

            if (weightType == "b") {
                
                int index, l_s, n_s;
                ss >> index >> l_s >> n_s;
                // std::cout << weight << " b " << index << " " << l_s << " " << n_s <<  std::endl;
                
                layers_biases[l_s-1][n_s-1] = weight;
                // cout << layers_biases[l_s-1][n_s-1] << " " << bias << endl;
            }
        }           
    }

    neuralNetwork.set_layers_synaptic_weights( layers_synaptic_weights );
    neuralNetwork.set_layers_biases( layers_biases );
}

double NeuralNetwork::calculateEnergy(const std::vector<double>& descriptorValues) {

    if (descriptorValues.size() != getNumberOfInputs())
        throw std::runtime_error("Unexpected size of inputs");

    // TODO: optimize the conversion
    OpenNN::Vector<double> inputs( getNumberOfInputs() );
    for(int i=0; i<getNumberOfInputs(); i++ )
        inputs[i] = descriptorValues[i];

    return neuralNetwork.calculate_outputs( inputs )[0]; // last layer has one output perceptron
}