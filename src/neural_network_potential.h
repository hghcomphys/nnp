/*
  base.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Hossein Ghorbanfekr

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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
