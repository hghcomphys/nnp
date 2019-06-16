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
    ACSF& getDescriptorForElement(const std::string& element);
    int getTotalNumberOfDescriptors();
    int getNumberOfElements();
    std::vector<std::string> &getElements();
    void calculate(Atoms &configuration);
    std::vector<double>& getDescriptorValuesForElement(const std::string& element);

private:
    std::string scriptFileName;
    std::vector<ACSF> descriptors;
    std::vector<std::string> elements;
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
