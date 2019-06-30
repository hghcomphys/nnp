#include <iostream>
#include <sstream>
#include "acsf.h"
#include "atoms.h"
#include "nnp.h"

using namespace std;
//using namespace NNP_SF;

#define descriptor descriptor_H

int main()
{
    try {

        ACSF descriptor_H("H"), descriptor_O("O");
        cout << descriptor.getCentralElement() << endl;

        // symfunction_short <element-central> 2 <element-neighbor> <eta> <rshift> <rcutoff>

        // # radial H H
        // symfunction_short H 2 H 0.001 0.0 12.00
        // symfunction_short H 2 H 0.01  0.0 12.00
        // symfunction_short H 2 H 0.03  0.0 12.00
        // symfunction_short H 2 H 0.06  0.0 12.00
        // symfunction_short H 2 H 0.15  1.9 12.00
        // symfunction_short H 2 H 0.30  1.9 12.00
        // symfunction_short H 2 H 0.60  1.9 12.00
        // symfunction_short H 2 H 1.50  1.9 12.00
        descriptor_H.addTwoBodySF(new G1({0.001, 0.0, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({0.01, 0.0, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({0.03, 0.0, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({0.06, 0.0, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({0.15, 1.9, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({0.30, 1.9, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({0.60, 1.9, 12.00}), "H");
        // descriptor_H.addTwoBodySF(new G1({1.50, 1.9, 12.00}), "H");

        // # radial H O / O H
        // symfunction_short H 2 O 0.001 0.0 12.00
        // symfunction_short H 2 O 0.01  0.0 12.00
        // symfunction_short H 2 O 0.03  0.0 12.00
        // symfunction_short H 2 O 0.06  0.0 12.00
        // symfunction_short H 2 O 0.15  0.9 12.00
        // symfunction_short H 2 O 0.30  0.9 12.00
        // symfunction_short H 2 O 0.60  0.9 12.00
        // symfunction_short H 2 O 1.50  0.9 12.00
        // descriptor_H.addTwoBodySF(new G1({0.001, 0.0, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({0.01, 0.0, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({0.03, 0.0, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({0.06, 0.0, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({0.15, 0.9, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({0.30, 0.9, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({0.60, 0.9, 12.00}), "O");
        // descriptor_H.addTwoBodySF(new G1({1.50, 0.9, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.001, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.01, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.03, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.06, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.15, 0.9, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.30, 0.9, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.60, 0.9, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({1.50, 0.9, 12.00}), "O");

        // # radial O O
        // symfunction_short O 2 O 0.001 0.0 12.00
        // symfunction_short O 2 O 0.01  0.0 12.00
        // symfunction_short O 2 O 0.03  0.0 12.00
        // symfunction_short O 2 O 0.06  0.0 12.00
        // symfunction_short O 2 O 0.15  4.0 12.00
        // symfunction_short O 2 O 0.30  4.0 12.00
        // symfunction_short O 2 O 0.60  4.0 12.00
        // symfunction_short O 2 O 1.50  4.0 12.00
        // descriptor_O.addTwoBodySF(new G1({0.001, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.01, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.03, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.06, 0.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.15, 4.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.30, 4.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({0.60, 4.0, 12.00}), "O");
        // descriptor_O.addTwoBodySF(new G1({1.50, 4.0, 12.00}), "O");

        // symfunction_short <element-central> 3 <element-neighbor1> <element-neighbor2> <eta> <lambda> <zeta> <rcutoff> <<rshift>
        // symfunction_short O 3 H H 0.07  1.0 1.0  12.00000
        // symfunction_short H 3 O H 0.07  1.0 1.0  12.00000
        // symfunction_short O 3 H H 0.07 -1.0 1.0  12.00000
        // symfunction_short H 3 O H 0.07 -1.0 1.0  12.00000
        descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "O", "O" );
        descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "H", "H" );
        descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "O", "H" );
        descriptor.addThreeBodySF( new G4({0, 0.07,  1.0, 1.0,  12.00000}), "H", "O" );

        cout << "Two-body SF: " << descriptor.getNumberOfTwoBodySF() << endl;
        cout << "Three-body SF: " << descriptor.getNumberOfThreeBodySF() << endl;
        cout << "Total number of SF: " << descriptor.getTotalNumberOfSF() << endl;


        Atoms configuration;
        configuration.readFileFormatXYZ("water12.xyz");
        double cell[9] = {4, 0, 0, 0, 4, 0, 0, 0, 4};
        // configuration.setCell(cell);
        cout << "Number of atoms: " << configuration.getNumberOfAtoms() << endl;
        cout << "Is PBC: " << configuration.isPBC() << endl;

        //    Atom atom(1, 2, 3, "X", 0);
        //    cout << atom.element << " " << atom.x << " " << atom.z << " " << atom.z <<  " " << atom.index << endl;

        // auto atoms = configuration.getListOfAtoms();
        // for (auto index: configuration.getListOfIndexForElement("O"))
        //     cout << index << " " << atoms[index].getIndex() << endl;
        // cout << configuration.getListOfIndexForElement("H").size() << endl;

        //    for(auto &atom: ac.atoms) {
        //        cout << atom.element << " " << atom.x << " " << atom.y << " " << atom.z <<  " " << atom.index << endl;
        //         break;
        //    }

        descriptor.calculate(configuration);
        for (auto &sf: descriptor.getValues())
            cout << sf << " ";
        cout << endl;
        // cout << ac.atoms.size() << endl;
        // cout << acsf.descriptors[0]->calculate(ac) << endl;

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