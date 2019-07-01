#include <iostream>
#include <sstream>
#include "acsf.h"
#include "atoms.h"
#include "nnp.h"

using namespace std;
//using namespace NNP_SF;

int main()
{
    try {

        // -----------------------------------------------
        // make atomic structure and read data from a file
        // -----------------------------------------------
        Atoms configuration;

        // Note: xyz and cell are in Angstrom but internal unit is in Bohr!
        // configuration.readFileFormatXYZ("water12.xyz");
        // double cell[9] = {4, 0, 0, 0, 4, 0, 0, 0, 4}; // Angstrom
        // configuration.setCell(cell);

        configuration.readFileFormatRuNNer();

        cout << "Number of atoms: " << configuration.getNumberOfAtoms() << endl;
        cout << "Is PBC: " << configuration.isPBC() << endl;

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
        nnp.readScript();
        nnp.calculate(configuration);

        cout << "Number of elements: " << nnp.getNumberOfElements() << endl;
        cout << "Total number of descriptors: " << nnp.getTotalNumberOfDescriptors() << endl;
        for (auto& element: nnp.getElements()) {
            cout << "Element: " << element << endl;
            cout << "Two-body SF: " << nnp.getDescriptorForElement(element).getNumberOfTwoBodySF() << endl;
            cout << "Three-body SF: " << nnp.getDescriptorForElement(element).getNumberOfThreeBodySF() << endl;
            cout << "Total SF: " << nnp.getDescriptorForElement(element).getTotalNumberOfSF() << endl;

            for (auto &sf: nnp.getDescriptorForElement(element).getValues())
                cout << sf << " ";
            cout << endl;
        }
    }

    catch (runtime_error e)
    {
        cout << "Runtime Error: " << e.what() << '!' << endl;
    }

}