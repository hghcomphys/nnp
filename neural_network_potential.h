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
    NeuralNetwork* getNeuralNetworkForElement(const std::string& element);
    int getNumberOfElements() const;
    const std::vector<std::string>& getElements() const;
    const std::vector<std::vector<double>>& getDescriptorValuesForElement(const std::string& element);
    double calculateEnergy(AtomicStructure& structure, int AtomIndex);
    double caculateTotalEnergy(AtomicStructure& structure);
    std::vector<double> calculateForce(AtomicStructure& structure, int atomIndex);

private:
    std::vector<std::string> elements;
    std::vector<ACSF> descriptors;
    std::vector<SymmeryFunctionsScaler> scalers;
    std::vector<int> hiddenLayersSize;
    std::vector<std::string> activationFunctionTypes;
    std::vector<NeuralNetwork*> neuralNetworks;
    int getIndexForElement(const std::string& element) const;
};


#endif //NNP_NEURALNETWORKPOTENTIAL_H
