//
// Created by hossein on 6/13/19.
//

#include "nnp.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#define REORDER_SYMMETRY_FUNCTIONS // this flag is temporary

/* ----------------------------------------------------------------------
   setup for Neural Network Potential
------------------------------------------------------------------------- */

NeuralNetworkPotential::NeuralNetworkPotential() {}

NeuralNetworkPotential::~NeuralNetworkPotential() {}

void NeuralNetworkPotential::readSetupFiles(const std::string& directory)
{
    const std::string filename = directory + "input.nn";
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error("Unable to open script file " + filename);

    char cSpace = ' ';
    std::string line;

    int number_of_elements = 0;
    int number_of_hidden_layers = 0;
    
    while ( std::getline(inFile, line) ) {

        std::stringstream ss(line);
        std::string sIndvStr;

        while ( std::getline(ss, sIndvStr, cSpace) ) {

            std::string dummy;
            double ddummy;
            int idummy;

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
#ifndef REORDER_SYMMETRY_FUNCTIONS
            // add symmetryc functions from input.nn
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
                    getDescriptorForElement(centralElement).addTwoBodySF(new G2(params), neighborElement1);
                    // std::cout << "add G2(<" << centralElement << ">, " << neighborElement1 << "): ";
                    // for (auto& p: params) 
                    //     std::cout << p << ' ';
                    // std::cout << std::endl;
                    break;

                case 3:
                    ss >> neighborElement1 >> neighborElement2;
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> lambda >> zeta >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G4(params), neighborElement1, neighborElement2);
                    // std::cout << "add G4(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    // for (auto& p: params) 
                    //     std::cout << p << ' ';
                    // std::cout << std::endl;
                    break;

                case 9:
                    ss >> neighborElement1 >> neighborElement2;
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> lambda >> zeta >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G5(params), neighborElement1, neighborElement2);
                    // std::cout << "add G5(" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    // for (auto& p: params) 
                    //     std::cout << p << ' ';
                    // std::cout << std::endl;
                    break;
                
                default:
                    throw std::runtime_error("Unexpected symmetry function in " + filename);
                    break;
                }
            }
#endif
        }           
    }

#ifdef REORDER_SYMMETRY_FUNCTIONS
    // add symmetryc functions from scaling log
    // TODO: improve design
    for (auto &element: elements) 
    {
        const std::string filenameSF = directory + "nnp-scaling.log.0000";
        std::ifstream inFileSF(filenameSF);
        if (!inFileSF)
            throw std::runtime_error("Unable to open file " + filenameSF);

        while ( std::getline(inFileSF, line) ) {

            int sfIndex, sfType;
            std::stringstream ss(line);
            std::string centralElement;
            std::string neighborElement1, neighborElement2;
            std::vector<double> params;
            double eta, zeta, lambda, rshift, rcutoff;

            ss >> sfIndex >> centralElement >> sfType;        
            if ( element == centralElement) 
            {
                switch (sfType)
                {
                case 2:
                    ss >> neighborElement1 >> eta >> rshift >> rcutoff;
                    params.push_back(eta); 
                    params.push_back(rshift); 
                    params.push_back(rcutoff);
                    getDescriptorForElement(centralElement).addTwoBodySF(new G2(params), neighborElement1);
                    std::cout << "add G2(<" << centralElement << ">, " << neighborElement1 << "): ";
                    for (auto& p: params) 
                        std::cout << p << ' ';
                    std::cout << std::endl;
                    break;

                case 3:
                    ss >> neighborElement1 >> neighborElement2 >> eta >> rshift >> lambda >> zeta >> rcutoff;
                    params.push_back(eta);
                    params.push_back(lambda);
                    params.push_back(zeta);
                    params.push_back(rcutoff);
                    // params.push_back(rshift);
                    getDescriptorForElement(centralElement).addThreeBodySF(new G4(params), neighborElement1, neighborElement2);
                    std::cout << "add G4(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        std::cout << p << ' ';
                    std::cout << std::endl;
                break;

                case 9:
                     ss >> neighborElement1 >> neighborElement2 >> eta >> rshift >> lambda >> zeta >> rcutoff;
                    params.push_back(eta);
                    params.push_back(lambda);
                    params.push_back(zeta);
                    params.push_back(rcutoff);
                    // params.push_back(rshift);
                    getDescriptorForElement(centralElement).addThreeBodySF(new G5(params), neighborElement1, neighborElement2);
                    std::cout << "add G5(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        std::cout << p << ' ';
                    std::cout << std::endl;
                    break;
                
                default:
                    throw std::runtime_error("Unexpected symmetry function in " + filenameSF);
                    break;
                }
            }       
        } 
        inFileSF.close();  
    }
#endif

    // reading scaling data
    // TODO: improve design
    const std::string filenameScale = directory + "scaling.data";
    std::ifstream inFileScale(filenameScale);
    if (!inFileScale)
        throw std::runtime_error("Unable to open file " + filenameScale);

    while ( std::getline(inFileScale, line) ) {

        if( line[0] == '#' ) continue;  // ignore comments

        std::stringstream ss(line);
        double sfMin, sfMax, sfMean, sfSigma;
        int elementIndex, sfIndex;

        ss >> elementIndex >> sfIndex >> sfMin >> sfMax >> sfMean >> sfSigma;
        // std::cout << elements[elementIndex-1] << " " << sfIndex << " " << sfMin << " " << sfMax << " " << sfMean << " " << sfSigma << "\n";
        getDescriptorForElement(elements[elementIndex-1]).addScaler( Scaler(sfMin, sfMax, sfMean, sfSigma) );
    } 
    inFileScale.close(); 
    // enable scaling
    for (auto &element: elements)
        getDescriptorForElement(element).scaleSymmetryFunctions();

    // create neural network for each element
    for (auto &element: elements) {
        const int numberOfInputs = getDescriptorForElement(element).getTotalNumberOfSF();
        neuralNetworks.push_back( NeuralNetwork(numberOfInputs, hiddenLayersSize) );
        std::cout << "Neural Network (" << element << "):" << std::endl;
    }

    // initilize neural network for each element
    for (auto &element: elements) {
            
            // read weights
            char filename[16];
            sprintf(filename, "weights.%3.3d.data", Atom::getAtomicNumber(element));
            const std::string fullPathFileName = directory + std::string(filename);
            getNeuralNetworkForElement(element).readParameters(fullPathFileName);
            std::cout << fullPathFileName << std::endl;

            // set activation functions
            getNeuralNetworkForElement(element).setLayersActivationFunction(activationFunctionTypes);
    }
}

void NeuralNetworkPotential::readSetupFiles() { readSetupFiles(""); }

int NeuralNetworkPotential::getNumberOfElements() const { return elements.size(); }

const std::vector<std::string>& NeuralNetworkPotential::getElements() const { return elements; }

int NeuralNetworkPotential::getIndexForElement(const std::string& element) const {
    // TODO: optimize the finding algorithm
    int index;
    for(index = 0; index<elements.size(); index++)
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

double NeuralNetworkPotential::calculateEnergy(Atoms& configuration, int atomIndex) {
    Atom& atom = configuration.getListOfAtoms()[atomIndex];
    std::vector<double> descriptorValues = getDescriptorForElement(atom.getElement()).calculateSF(configuration, atomIndex);
    return getNeuralNetworkForElement(atom.getElement()).calculateEnergy(descriptorValues);
}

double NeuralNetworkPotential::caculateTotalEnergy(Atoms &configuration) {
    double totalEnergy = 0.0;
    for(auto atom: configuration.getListOfAtoms())
        totalEnergy += calculateEnergy(configuration, atom.getIndex());
    return totalEnergy;
}