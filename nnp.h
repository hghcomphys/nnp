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

private:
    std::string scriptFileName;
    std::vector<ACSF> descriptors;
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
