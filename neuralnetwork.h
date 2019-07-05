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
    NeuralNetwork::NeuralNetwork(const int& numberOfInputs, const std::vector<int>& hiddenLayersSize);
    ~NeuralNetwork();
    // int getNumberOfInputs() const;
    // int getNumberOfOutputs() const;
    // const std::vector<int>& getNumberOfHiddenLayers() const;

private:
    OpenNN::MultilayerPerceptron neuralNetwork;
};


#endif //NNP_NEURALNETWORK_H
