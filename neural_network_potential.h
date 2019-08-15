//
// Neural Network Potential
//

#ifndef NNP_NEURALNETWORKPOTENTIAL_H
#define NNP_NEURALNETWORKPOTENTIAL_H

#include <vector>
#include "acsf.h"
#include "symmfuncscaler.h"
#include "neuralnetwork.h"


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
    double calculateEnergy(Atoms& configuration, int AtomIndex);
    double caculateTotalEnergy(Atoms& configuration);
    std::vector<double> calculateForce(Atoms& configuration, int atomIndex);

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