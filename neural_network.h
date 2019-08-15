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
    // const OpenNN::MultilayerPerceptron* getPerceptron() const;
    int getNumberOfInputs() const;
    int getNumberOfOutputs() const;
    int getNumberOfLayers() const;
    int getNumberOfHiddenLayers() const;
    double calculateEnergy(const std::vector<double>& descriptorValues);
    OpenNN::Vector<double> calculateJacobian(const std::vector<double>& descriptorValues);

private:
    // OpenNN::Vector<int> layersSize; 
    OpenNN::MultilayerPerceptron* multilayerPerceptron;
};


#endif //NNP_NEURALNETWORK_H