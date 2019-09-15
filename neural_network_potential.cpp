//
// Neural Network Potential
//

#include "neural_network_potential.h"
#include "logger.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

#define REORDER_SYMMETRY_FUNCTIONS  // this flag is temporary

/* ----------------------------------------------------------------------
   setup for Neural Network Potential
------------------------------------------------------------------------- */
NeuralNetworkPotential::NeuralNetworkPotential() {}

NeuralNetworkPotential::~NeuralNetworkPotential() {
    for (auto each: neuralNetworks)
        delete each;    
    neuralNetworks.clear();
}

void NeuralNetworkPotential::readSetupFiles(const std::string& directory)
{
    const std::string filename = directory + "input.nn";
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error( (Log(ERROR) << "Unable to open script file " + filename).toString() );
    Log(INFO) << "Read " << filename;

    char cSpace = ' ';
    std::string line;

    int number_of_elements = 0;
    int number_of_hidden_layers = 0;
    
    while ( std::getline(inFile, line) ) 
    {
        std::stringstream ss(line);
        std::string sIndvStr;

        while ( std::getline(ss, sIndvStr, cSpace) ) 
        {
            std::string dummy;
            double ddummy;
            int idummy;

            if (sIndvStr == "number_of_elements") {
                ss >> number_of_elements;
                Log(INFO) << sIndvStr + ": " << number_of_elements;
            }
            else if (sIndvStr == "elements") 
            {
                if (number_of_elements == 0)
                        throw std::runtime_error( (Log(ERROR) << "Number of elemets is zero").toString() );

                for (int i=0; i<number_of_elements; i++) {
                    ss >> dummy;
                    elements.push_back(dummy);
                    Log(INFO) << "Element: " << elements[i];
                }
                for (auto &element: elements) {
                    descriptors.push_back( ACSF(element) );  // add descriptor
                    scalers.push_back( SymmeryFunctionsScaler() );  // add scaler
                    Log(INFO) << "ACSF(" << element << ")";// << std::endl;
                }
            }
            else if (sIndvStr == "global_hidden_layers_short") 
            {
                ss >> number_of_hidden_layers;
                Log(INFO) << sIndvStr + ": " << number_of_elements;
            } 
            else if (sIndvStr == "global_nodes_short") 
            {
                if (number_of_hidden_layers == 0)
                        throw std::runtime_error( (Log(ERROR) << "Number of hidden layers is zero").toString() );

                for (int i=0; i<number_of_hidden_layers; i++) {
                    ss >> idummy;
                    hiddenLayersSize.push_back(idummy);
                    Log(INFO) << sIndvStr + ": " << hiddenLayersSize[i];
                }
            }
            else if (sIndvStr == "global_activation_short") 
            {
                const int number_of_layers = number_of_hidden_layers + 1; // plus output layer
                for (int i=0; i<number_of_layers; i++) {
                    ss >> dummy;
                    activationFunctionTypes.push_back(dummy);
                    Log(INFO) << sIndvStr + ": " << activationFunctionTypes[i];
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
                Log log(DEBUG);

                ss >> centralElement >> sfType;
                switch (sfType)
                {
                case 2:
                    // two-body symmetry functions
                    ss >> neighborElement1;
                    for (int i=0; i<3; i++) {
                        ss >> ddummy; //eta >> rshift >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addTwoBodySF(new G2(params), neighborElement1);
                    // log
                    log << "G2(<" << centralElement << ">, " << neighborElement1 << "): ";
                    for (auto& p: params) 
                        log << p << " ";
                    break;

                case 3:
                    // three-body symmetry functions (narrow)
                    ss >> neighborElement1 >> neighborElement2;
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> lambda >> zeta >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G4(params), neighborElement1, neighborElement2);
                    // log
                    log << "G4(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        log << p << " ";
                    break;

                case 9:
                    // three-body symmetry functions (wide)
                    ss >> neighborElement1 >> neighborElement2;
                    for (int i=0; i<4; i++) {
                        ss >> ddummy; //eta >> lambda >> zeta >> rcutoff;
                        params.push_back(ddummy);
                    }
                    getDescriptorForElement(centralElement).addThreeBodySF(new G5(params), neighborElement1, neighborElement2);
                    // log
                    log << "G5(" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        log << p << " ";
                    break;
                
                default:
                    throw std::runtime_error(  (Log(ERROR) << "Unexpected symmetry function in " + filename).toString() );
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
            throw std::runtime_error( (Log(ERROR) << "Unable to open file " + filenameSF).toString() );
        Log(INFO) << "Read " << filenameSF << " (" + element + ")";

        while ( std::getline(inFileSF, line) ) {

            int sfIndex, sfType;
            std::stringstream ss(line);
            std::string centralElement;
            std::string neighborElement1, neighborElement2;
            std::vector<double> params;
            double eta, zeta, lambda, rshift, rcutoff;

            // read symmetry function parameters
            ss >> sfIndex >> centralElement >> sfType;        
            if ( element == centralElement) 
            {
                Log log(DEBUG);
                switch (sfType)
                {
                case 2:
                    // two-body symmetry functions
                    ss >> neighborElement1 >> eta >> rshift >> rcutoff;
                    params.push_back(eta); 
                    params.push_back(rshift); 
                    params.push_back(rcutoff);
                    getDescriptorForElement(centralElement).addTwoBodySF(new G2(params), neighborElement1);
                    // log 
                    log << "G2(<" << centralElement << ">, " << neighborElement1 << "): ";
                    for (auto& p: params) 
                        log << p << " ";
                    break;

                case 3:
                    // three-body symmetry functions (narrow)
                    ss >> neighborElement1 >> neighborElement2 >> eta >> rshift >> lambda >> zeta >> rcutoff;
                    params.push_back(eta);
                    params.push_back(lambda);
                    params.push_back(zeta);
                    params.push_back(rcutoff);
                    // params.push_back(rshift);        
                    // log
                    getDescriptorForElement(centralElement).addThreeBodySF(new G4(params), neighborElement1, neighborElement2);
                    log << "G4(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        log << p << " ";
                    break;

                case 9:
                    // three-body symmetry functions (wide)
                    ss >> neighborElement1 >> neighborElement2 >> eta >> rshift >> lambda >> zeta >> rcutoff;
                    params.push_back(eta);
                    params.push_back(lambda);
                    params.push_back(zeta);
                    params.push_back(rcutoff);
                    // params.push_back(rshift);
                    // log
                    getDescriptorForElement(centralElement).addThreeBodySF(new G5(params), neighborElement1, neighborElement2);
                    log << "G5(<" << centralElement << ">, " << neighborElement1 << ", " << neighborElement2 << "): ";
                    for (auto& p: params) 
                        log << p << " ";
                    break;
                
                default:
                    throw std::runtime_error(
                        (Log(ERROR) << "Unexpected symmetry function type (" << sfType << ") in "+ filenameSF).toString()
                    );
                    break;
                }
            }       
        } 
        inFileSF.close();  
    }
#endif

    // reading scaling data
    for (int index=0; index<getNumberOfElements(); index++) {
        const std::string filename = "scaling.data";
        const std::string element = elements[index];
        getScalerForElement(element).readScaling(directory + filename , index+1);
        Log(INFO) << "Read " + filename + " (" + element + ")";
    }
   
    // create neural network for each element
    for (auto &element: elements) {
        const int numberOfInputs = getDescriptorForElement(element).getTotalNumberOfSF();
        neuralNetworks.push_back( new NeuralNetwork(numberOfInputs, hiddenLayersSize) );
        Log(INFO) << "Neural Network (" + element + ")";
    }

    // initilize neural network for each element
    for (auto &element: elements) {

            // read weights
            char filename[32];
            sprintf(filename, "weights.%3.3d.data", getAtomicNumber(element));
            const std::string fullPathFileName = directory + std::string(filename);
            Log(INFO) << "Read " << fullPathFileName;

            // read parameters file into neural network
            getNeuralNetworkForElement(element).readParameters(fullPathFileName);

            // set activation functions
            getNeuralNetworkForElement(element).setLayersActivationFunction(activationFunctionTypes);
    
            // log info
            Log(INFO) << "Initialize weights for Neural network (" + element + ")";
    }

}

void NeuralNetworkPotential::readSetupFiles() { readSetupFiles(""); }

int NeuralNetworkPotential::getNumberOfElements() const { return elements.size(); }

const std::vector<std::string>& NeuralNetworkPotential::getElements() const { return elements; }

int NeuralNetworkPotential::getIndexForElement(const std::string& element) const 
{
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

SymmeryFunctionsScaler& NeuralNetworkPotential::getScalerForElement(const std::string& element) {
    return scalers[getIndexForElement(element)];
}

NeuralNetwork&  NeuralNetworkPotential::getNeuralNetworkForElement(const std::string& element) {
    return *neuralNetworks[getIndexForElement(element)];
}

double NeuralNetworkPotential::calculateEnergy(AtomicStructure& structure, int atomIndex) 
{
    Atom& atom = structure.getAtom(atomIndex);
    std::vector<double> descriptorValues = getDescriptorForElement(atom.element).calculate(structure, atomIndex);
    std::vector<double> scaledDescriptorValues = getScalerForElement(atom.element).scale(descriptorValues);
    return getNeuralNetworkForElement(atom.element).calculateEnergy(scaledDescriptorValues);
}

double NeuralNetworkPotential::caculateTotalEnergy(AtomicStructure& structure) 
{
    double totalEnergy = 0.0;
    for(int i=0; i<structure.numberOfAtoms; i++)
         totalEnergy += calculateEnergy(structure, structure.listOfAtoms[i]->index);
   
    return totalEnergy;
}

std::vector<double> NeuralNetworkPotential::calculateForce(AtomicStructure& structure, int atomIndex) 
{
    Atom *atom_i = &structure.getAtom(atomIndex);
    std::vector<double> force({0, 0, 0});

    // list of index of atoms
    // const auto& listofAtomIndex = structure.getListOfAtomIndex();

    // sum over all atoms
    for (int j=0; j<structure.numberOfAtoms; j++) 
    {       
        // refer to atom
        Atom *atom_j = structure.listOfAtoms[j];

        // TODO: improve 
        const double rij = structure.distance(atom_i, atom_j);
        if ( rij > getDescriptorForElement(atom_j->element).getGlobalCutOffRadius() ) continue; // TODO: fix it!

        // gradient of neural network respect to symmetry functions
        const std::vector<double>& descriptorValues = getDescriptorForElement(atom_j->element).calculate(structure, atom_j->index);
        const std::vector<double>& scaledDescriptorValues = getScalerForElement(atom_j->element).scale(descriptorValues);
        const OpenNN::Vector<double>&  networkGradient = getNeuralNetworkForElement(atom_j->element).calculateJacobian(scaledDescriptorValues);
        const std::vector<double>& scalingFactors = getScalerForElement(atom_j->element).getScalingFactors();          

        // gradient of symmetry functions respect to atomic positions
        const std::vector<std::vector<double>>& descriptorGradient = getDescriptorForElement(atom_j->element).gradient(structure, atom_j->index, atom_i->index);
       
        // sum over symmetry functions
        for (int n=0; n<descriptorValues.size(); n++ ) 
        {
            for (int d=0; d<3; d++)
                force[d] -=  scalingFactors[n] * networkGradient[n] * descriptorGradient[n][d];          
        }
    }

    // return force vector applied on atomIndex
    return force;
}