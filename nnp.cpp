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

NeuralNetworkPotential::NeuralNetworkPotential() {}

NeuralNetworkPotential::~NeuralNetworkPotential() {}

void NeuralNetworkPotential::addDescriptor(const ACSF& descriptor) {
    descriptors.push_back( descriptor );
}

void NeuralNetworkPotential::readScript(const std::string& fileName)
{
    scriptFileName = fileName;
    std::ifstream inFile(fileName);
    if (!inFile) {
        throw std::runtime_error("Unable to open script file");
        //exit(1);   // call system to stop
    }

    char cSpace = ' ';
    std::string line;
    while ( std::getline(inFile, line) ) {
        std::stringstream ss(line);
        std::string sIndvStr;
        while ( std::getline(ss, sIndvStr, cSpace) ) {
            std::cout << sIndvStr << std::endl;
            break;
        }
    }
}

void NeuralNetworkPotential::readScript() { readScript("input.nn"); }

