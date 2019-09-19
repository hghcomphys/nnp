//
// Neural Network Potential
//

#ifndef NNP_NEURALNETWORKPOTENTIAL_H
#define NNP_NEURALNETWORKPOTENTIAL_H

#include <vector>
#include "descriptor.h"
#include "symmetry_function_scaler.h"
#include "neural_network.h"

class NeuralNetworkPotential {
public:
    NeuralNetworkPotential();
    ~NeuralNetworkPotential();
    void readSetupFiles();
    void readSetupFiles(const std::string& directory);
    ACSF& getDescriptorForElement(const std::string& element);
    SymmeryFunctionsScaler& getScalerForElement(const std::string& element);
    NeuralNetwork& getNeuralNetworkForElement(const std::string& element);
    int getNumberOfElements();
    const std::vector<std::vector<double>>& getDescriptorValuesForElement(const std::string& element);
    void calculateEnergy(AtomicStructure& structure, Atom *atom_i);
    void caculateTotalEnergy(AtomicStructure& structure);
    void calculateForce(AtomicStructure& structure, Atom *atom);
    void calculateForce(AtomicStructure& structure);

// private:
    std::vector<std::string> elements;
    std::vector<ACSF> descriptors;
    std::vector<SymmeryFunctionsScaler> scalers;
    std::vector<int> hiddenLayersSize;
    std::vector<std::string> activationFunctionTypes;
    std::vector<NeuralNetwork*> neuralNetworks;
    int getIndexForElement(const std::string& element);
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
