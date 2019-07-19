#include <iostream>
#include <sstream>
#include "acsf.h"
#include "atoms.h"
#include "nnp.h"

#include "neuralnetwork.h"

// #include <tensorflow/c/c_api.h>
// #include "opennn.h"

using namespace std;
// using namespace OpenNN;

int main()
{
    try {

        // -----------------------------------------------
        // make atomic structure and read data from a file
        // -----------------------------------------------
        Atoms configuration;
        configuration.readFileFormatRuNNer();

        cout << "Number of atoms: " << configuration.getNumberOfAtoms() << endl;
        cout << "Is PBC: " << configuration.isPBC() << endl;

        // Note: xyz and cell are in Angstrom but internal unit is in Bohr!
        // configuration.readFileFormatXYZ("water12.xyz");
        // double cell[9] = {4, 0, 0, 0, 4, 0, 0, 0, 4}; // Angstrom
        // configuration.setCell(cell);

        //    Atom atom(1, 2, 3, "X", 0);
        //    cout << atom.element << " " << atom.x << " " << atom.z << " " << atom.z <<  " " << atom.index << endl;

        // auto atoms = configuration.getListOfAtoms();
        // for (auto index: configuration.getListOfIndexForElement("H"))
        //     cout << index << " " << atoms[index].getZ() << " " << atoms[index].getElement() << endl;
        // cout << configuration.getListOfIndexForElement("H").size() << endl;

        //    for(auto &atom: ac.atoms) {
        //        cout << atom.element << " " << atom.x << " " << atom.y << " " << atom.z <<  " " << atom.index << endl;
        //         break;
        //    }


        // ----------------------------------------------
        // create a descriptor and add symmetry functions
        // ----------------------------------------------
        // ACSF descriptor("O");
        // cout << descriptor.getCentralElement() << endl;

        // // symfunction_short <element-central> 2 <element-neighbor> <eta> <rshift> <rcutoff>
        // descriptor.addTwoBodySF(new G1({0.001, 0.0, 12.00}), "H");

        // // symfunction_short <element-central> 3 <element-neighbor1> <element-neighbor2> <eta> <lambda> <zeta> <rcutoff> <<rshift>
        // descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "O", "O" );
        // descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "H", "H" );
        // descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "O", "H" );
        // descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "H", "O" );

        // cout << "Two-body SF: " << descriptor.getNumberOfTwoBodySF() << endl;
        // cout << "Three-body SF: " << descriptor.getNumberOfThreeBodySF() << endl;
        // cout << "Total number of SF: " << descriptor.getTotalNumberOfSF() << endl;


        // -------------------------------------------------
        // calculate symmetry fucntions for a given strucure
        // -------------------------------------------------
        // descriptor.calculate(configuration);
        // for (auto &sf: descriptor.getValues())
        //     cout << sf << " ";
        // cout << endl;
        // cout << ac.atoms.size() << endl;
        // cout << acsf.descriptors[0]->calculate(ac) << endl;


        // -------------------------------------------------
        // make neural network potential
        // ------------------------------------------------

        NeuralNetworkPotential nnp;
        nnp.readSetupFiles();

        cout << "Number of elements: " << nnp.getNumberOfElements() << endl;
        for (auto& element: nnp.getElements()) 
        {
            cout << "Element: " << element << endl;
            cout << "Two-body SF: " << nnp.getDescriptorForElement(element).getNumberOfTwoBodySF() << endl;
            cout << "Three-body SF: " << nnp.getDescriptorForElement(element).getNumberOfThreeBodySF() << endl;
            cout << "Total SF: " << nnp.getDescriptorForElement(element).getTotalNumberOfSF() << endl;
            cout << "NN number of inputs: "  << nnp.getNeuralNetworkForElement(element).getNumberOfInputs() << endl;
            cout << "NN number of hidden layers: "  << nnp.getNeuralNetworkForElement(element).getNumberOfHiddenLayers() << endl;
        }

        for (int index=0; index<3; index++) 
        {
            Atom atom = configuration.getListOfAtoms()[index];
            // symmetry functions
            // cout << "SF--> ";
            // for(auto sf: nnp.getDescriptorForElement(atom.getElement()).calculateSF(configuration, atom.getIndex()))
            //     cout << sf << " ";
            // cout << "\n";
            // energy
            cout << "Atom[" << atom.getElement() << "]:(" << atom.getX() << ", " << atom.getY() << ", " << atom.getZ() << ") "
            << "energy: " << nnp.calculateEnergy(configuration, atom.getIndex()) << endl;
        }    
        cout << "Total energy: " << nnp.caculateTotalEnergy(configuration) << endl;

        // for (auto &atoms: nnp.getDescriptorForElement(element).getValues() ) {
        //     for (auto &sf: atoms)
        //         cout << sf << " ";
        //     cout << endl;
        //     break;
        // }
        // // cout << endl; // << endl;

        // for (auto &index: nnp.getDescriptorForElement(element).getIndex2() )
        //     cout << index << " ";
        // cout << endl;

        // for (auto &index: nnp.getDescriptorForElement(element).getIndex3() )
        //     cout << index << " ";
        // cout << endl;


        // -----------------------------------------------
        // Nerural Network
        // -----------------------------------------------
        // NeuralNetwork nn(10, {25, 25});
        // cout << "size of inputs: " << nn.getNumberOfInputs() << endl;
        // cout << "size of outputs: " << nn.getNumberOfOutputs() << endl;
        // cout << "number of hidden layers: " << nn.getNumberOfHiddenLayers() << endl;  //hidden layers + output layer 

        // OpenNN::Vector<double> input(1);
        // input[0] = 0.5;
        // OpenNN::Vector<double> output = nn.getNeuralNetwork().calculate_outputs( input );
        // cout << "output = " << output << endl;

        // OpenNN::Matrix<double> jacobian = nn.getNeuralNetwork().calculate_Jacobian ( input );
        // cout << "jacobian = " << jacobian << endl;

        // std::cout << "TensorFlow Version: " << TF_Version() << std::endl;
    }

    catch (runtime_error e)
    {
        cout << "Runtime Error: " << e.what() << "!" << endl;
    }

}