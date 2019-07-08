//
// Created by hossein on 6/13/19.
//

#ifndef NNP_NEURALNETWORKPOTENTIAL_H
#define NNP_NEURALNETWORKPOTENTIAL_H

#include <vector>
#include "acsf.h"
#include "neuralnetwork.h"


class NeuralNetworkPotential {
public:
    NeuralNetworkPotential(const std::string& directory);
    ~NeuralNetworkPotential();
    void initilize();
    ACSF& getDescriptorForElement(const std::string& element);
    NeuralNetwork& NeuralNetworkPotential::getNeuralNetworkForElement(const std::string& element);
    int getNumberOfElements() const;
    const std::vector<std::string>& getElements() const;
    const std::vector<std::vector<double>>& getDescriptorValuesForElement(const std::string& element);
    void calculateDescriptor(Atoms &configuration);

private:
    std::string directory;
    std::vector<std::string> elements;
    std::vector<ACSF> descriptors;
    std::vector<int> hiddenLayersSize;
    std::vector<std::string> activationFunctionTypes;
    std::vector<NeuralNetwork> neuralNetworks;
    int NeuralNetworkPotential::getIndexForElement(const std::string& element) const;
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
