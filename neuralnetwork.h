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
    NeuralNetwork();
    // ~NeuralNetwork();
    NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize, int outputsSize);
    NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize);
    const OpenNN::MultilayerPerceptron& getNeuralNetwork() const;
    int getNumberOfInputs() const;
    int getNumberOfOutputs() const;
    int getNumberOfHiddenLayers() const;

private:
    OpenNN::MultilayerPerceptron neuralNetwork;
    
};


#endif //NNP_NEURALNETWORK_H
