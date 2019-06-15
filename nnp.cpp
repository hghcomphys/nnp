//
// Created by hossein on 6/13/19.
//

#include "nnp.h"
#include <fstream>
#include <sstream>
#include <iostream>

/* ----------------------------------------------------------------------
   setup for Neural Network Potential
------------------------------------------------------------------------- */

NeuralNetworkPotential::NeuralNetworkPotential(): number_of_elements(0) {}

NeuralNetworkPotential::~NeuralNetworkPotential() {}

void NeuralNetworkPotential::addDescriptor(const ACSF& descriptor) {
    descriptors.push_back( descriptor );
}

ACSF& NeuralNetworkPotential::getDescriptorForElement(const std::string element) {
    for (auto& descriptor: descriptors)
        if (descriptor.getCentralElement() == element)
            return descriptor;
}

int NeuralNetworkPotential::getTotalNumberOfDescriptors() {
    int totDes = 0;
    for (auto& descriptor: descriptors)
        totDes += descriptor.getTotalNumberOfSF();
    return totDes;
}

void NeuralNetworkPotential::readScript(const std::string& fileName)
{
    scriptFileName = fileName;
    std::ifstream inFile(fileName);
    if (!inFile)
        throw std::runtime_error("Unable to open script file");

    char cSpace = ' ';
    std::string line, dummy;
    double ddummy;
    
    while ( std::getline(inFile, line) ) {
        std::stringstream ss(line);
        std::string sIndvStr;
        while ( std::getline(ss, sIndvStr, cSpace) ) {

            if (sIndvStr == "number_of_elements") {
                ss >> number_of_elements;
                std::cout << sIndvStr << ' ' << number_of_elements << std::endl;
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
                    std::cout << "new ACSF( " << element << " )" << std::endl;
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
                    getDescriptorForElement(centralElement).addTwoBodySF(new G1(params), neighborElement1 );
                    std::cout << "add G1(" << centralElement << ", " << neighborElement1 << "): ";
                    for (auto& p: params) 
                        std::cout << p << ' ';
                    std::cout << std::endl;
                    break;

                case 3:
                    ss >> neighborElement1 >> neighborElement2;
                    params.push_back(0.0); // TODO:: what is this param for?
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> rshift >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G4(params), neighborElement1, neighborElement2 );
                    std::cout << "add G4(" << centralElement << ", " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        std::cout << p << ' ';
                    std::cout << std::endl;
                    break;
                
                default:
                    break;
                }
            }
            
        }
                
    }
}

void NeuralNetworkPotential::readScript() { readScript("input.nn"); }

