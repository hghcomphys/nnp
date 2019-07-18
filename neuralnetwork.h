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
    NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize, int outputsSize);
    NeuralNetwork::NeuralNetwork(int inputsSize, const std::vector<int>& hiddenLayersSize);
    ~NeuralNetwork();
    void setLayersActivationFunction(const std::vector<std::string>& activationFucntionsType);
    void readParameters(const std::string& fullPathFileName);
    const OpenNN::MultilayerPerceptron& getNeuralNetwork() const;
    int getNumberOfInputs() const;
    int getNumberOfOutputs() const;
    int getNumberOfLayers() const;
    int getNumberOfHiddenLayers() const;
    double calculateEnergy(const std::vector<double>& descriptorValues);

private:
    // OpenNN::Vector<int> layersSize; 
    OpenNN::MultilayerPerceptron neuralNetwork;
};


#endif //NNP_NEURALNETWORK_H
