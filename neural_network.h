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
