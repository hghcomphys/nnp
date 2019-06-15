//
// Created by hossein on 6/13/19.
//

#ifndef NNP_NEURALNETWORKPOTENTIAL_H
#define NNP_NEURALNETWORKPOTENTIAL_H

#include <vector>
#include "acsf.h"


class NeuralNetworkPotential {
public:
    NeuralNetworkPotential();
    ~NeuralNetworkPotential();
    void addDescriptor(const ACSF& descriptor);
    void readScript();
    void readScript(const std::string& fileName);
    ACSF& getDescriptorForElement(const std::string element);
    int getTotalNumberOfDescriptors();

private:
    std::string scriptFileName;
    std::vector<ACSF> descriptors;
    int number_of_elements;
    std::vector<std::string> elements;
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
