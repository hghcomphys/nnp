//
// Created by hossein on 6/13/19.
//

#include "nnp.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>


/* ----------------------------------------------------------------------
   setup for Neural Network Potential
------------------------------------------------------------------------- */

NeuralNetworkPotential::NeuralNetworkPotential(const std::string& directory): directory(directory) { readSetupFiles(); }

int NeuralNetworkPotential::getNumberOfElements() const { return elements.size(); }

const std::vector<std::string>& NeuralNetworkPotential::getElements() const { return elements; }

int NeuralNetworkPotential::getIndexForElement(const std::string& element) const {
    // TODO: optimize the finding algorithm
    int index;
    for(index = 0; index<element.size(); index++)
        if (elements[index] == element)
            break;
    return index;
}

ACSF& NeuralNetworkPotential::getDescriptorForElement(const std::string& element) { 
    return descriptors[getIndexForElement(element)];
}

NeuralNetwork& NeuralNetworkPotential::getNeuralNetworkForElement(const std::string& element) { 
    return neuralNetworks[getIndexForElement(element)];
}

int elementToAtomicNumber(const std::string& element) 
{
    // TODO: extend to all elements
    if(element == "H")
        return 1;
    else if (element == "C")
        return 6; 
    else if (element == "O")
        return 8;
}

void NeuralNetworkPotential::readSetupFiles()
{
    const std::string filename = directory + "input.nn";
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error("Unable to open script file");

    char cSpace = ' ';
    std::string line, dummy;
    double ddummy;
    int idummy;
    int number_of_elements = 0;
    int number_of_hidden_layers = 0;
    
    while ( std::getline(inFile, line) ) {
        std::stringstream ss(line);
        std::string sIndvStr;
        while ( std::getline(ss, sIndvStr, cSpace) ) {

            if (sIndvStr == "number_of_elements") {
                ss >> number_of_elements;
                // std::cout << sIndvStr << ' ' << number_of_elements << std::endl;
            }
            else if (sIndvStr == "elements") {

                if (number_of_elements == 0)
                        throw std::runtime_error("Number of elemets is zero");

                for (int i=0; i<number_of_elements; i++) {
                    ss >> dummy;
                    elements.push_back(dummy);
                    // std::cout << elements[i] << std::endl;
                }
                for (auto &element: elements) {
                    descriptors.push_back( ACSF(element) );
                    std::cout << "ACSF(" << element << ")" << std::endl;
                }
            }
            else if (sIndvStr == "symfunction_short")
            {
                int sfType;
                std::string centralElement;
                std::string neighborElement1, neighborElement2;
                std::vector<double> params;
                ss >> centralElement >> sfType;
                switch (sfType)
                {
                case 2:
                    ss >> neighborElement1;
                    for (int i=0; i<3; i++) {
                        ss >> ddummy; //eta >> rshift >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addTwoBodySF(new G2(params), neighborElement1 );
                    // std::cout << "add G2(<" << centralElement << ">, " << neighborElement1 << "): ";
                    // for (auto& p: params) 
                    //     std::cout << p << ' ';
                    // std::cout << std::endl;
                    break;

                case 3:
                    ss >> neighborElement1 >> neighborElement2;
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> rshift >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G4(params), neighborElement1, neighborElement2 );
                    // std::cout << "add G4(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    // for (auto& p: params) 
                    //     std::cout << p << ' ';
                    // std::cout << std::endl;
                    break;

                case 9:
                    ss >> neighborElement1 >> neighborElement2;
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> rshift >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G5(params), neighborElement1, neighborElement2 );
                    // std::cout << "add G5(" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    // for (auto& p: params) 
                    //     std::cout << p << ' ';
                    // std::cout << std::endl;
                    break;
                
                default:
                    break;
                }
            }
            else if (sIndvStr == "global_hidden_layers_short") {
                ss >> number_of_hidden_layers;
                // std::cout << sIndvStr << ' ' << number_of_elements << std::endl;
            } 
            else if (sIndvStr == "global_nodes_short") {

                if (number_of_hidden_layers == 0)
                        throw std::runtime_error("Number of hidden layers is zero");

                for (int i=0; i<number_of_hidden_layers; i++) {
                    ss >> idummy;
                    hiddenLayersSize.push_back(idummy);
                    // std::cout << hiddenLayersSize[i] << std::endl;
                }
            }
            else if (sIndvStr == "global_activation_short") {

                const int number_of_layers = number_of_hidden_layers + 1; // plus output layer
                for (int i=0; i<number_of_layers; i++) {
                    ss >> dummy;
                    activationFunctionTypes.push_back(dummy);
                    // std::cout << activationFunctionTypes[i] << std::endl;
                }
            }
        }           
    }

    // create neural network for each element
    for (auto &element: elements) {
        neuralNetworks.push_back( NeuralNetwork(getDescriptorForElement(element).getTotalNumberOfSF(), hiddenLayersSize) );
        std::cout << "Neural Network (" << element << "):" << std::endl;
    }

    // initilize neural network for each element
    for (auto &element: elements) {
            char filename[16];
            sprintf(filename, "weights.%3.3d.data", elementToAtomicNumber(element));
            const std::string fullPathFileName = directory + std::string(filename);
            std::cout << fullPathFileName << std::endl;
            getNeuralNetworkForElement(element).readParameters(fullPathFileName);
    }
}

// void NeuralNetworkPotential::initilize() { initilize(); }

// void NeuralNetworkPotential::calculateDescriptor(Atoms &configuration) {
//     for (auto& descriptor: descriptors)
//         descriptor.calculate(configuration);
// }

// const std::vector<std::vector<double>>& NeuralNetworkPotential::getDescriptorValuesForElement(const std::string& element) {
//     return getDescriptorForElement(element).getValues();
// }

double NeuralNetworkPotential::calculateEnergy(Atoms& configuration, int atomIndex) {
    Atom atom = configuration.getListOfAtoms()[atomIndex];
    std::vector<double> descriptorValues = getDescriptorForElement(atom.getElement()).calculateSF(configuration, atomIndex);
    return getNeuralNetworkForElement(atom.getElement()).calculateEnergy(descriptorValues);
}

double NeuralNetworkPotential::caculateTotalEnergy(Atoms &configuration) {
    double totalEnergy = 0.0;
    for(auto atom: configuration.getListOfAtoms())
        totalEnergy += calculateEnergy(configuration, atom.getIndex());
    return totalEnergy;
}