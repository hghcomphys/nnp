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
    NeuralNetworkPotential();
    ~NeuralNetworkPotential();
    void addDescriptor(const ACSF& descriptor);
    void initilize();
    void initilize(const std::string& fileName);
    ACSF& getDescriptorForElement(const std::string& element);
    int getTotalNumberOfDescriptors();
    int getNumberOfElements();
    std::vector<std::string> &getElements();
    void calculateDescriptor(Atoms &configuration);
    std::vector<std::vector<double>>& getDescriptorValuesForElement(const std::string& element);

private:
    std::string scriptFileName;
    std::vector<ACSF> descriptors;
    std::vector<std::string> elements;
    std::vector<NeuralNetwork> neuralnetworks;
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
