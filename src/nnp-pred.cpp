/*
  nnp-pred.cpp: This file is part of Free Molecular Dynamics

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

#include <iostream>
#include <sstream>
#include "structure.h"
#include "neural_network_potential.h"

using namespace std;

int main(int argc, char *argv[])
{
    try
    {

        // string variables for nnp-pred.x input arguments
        string potentialPath, inStructureFile, outStructureFile;

        // read input argument
        if (argc == 1)
        {
            cerr << "Expected (1) path to NNP potential and (2) input structure filename\n"
                 << "        Example run: ./nnp-pred.x potentials/water example/input.data";
            exit(1);
        }
        else
        {
            potentialPath = string(argv[1]) + "/";
            cout << "NNP directory: " << potentialPath << endl;

            inStructureFile = string(argv[2]);
            cout << "Input: " << inStructureFile << endl;

            outStructureFile = inStructureFile + ".pred";
            cout << "Output: " << outStructureFile << endl;
        }

        // -----------------------------------------------
        // make atomic structure and read data from a file
        // -----------------------------------------------
        AtomicStructure structure;
        structure.readFileFormatRuNNer(inStructureFile);
        structure.calculateTableOfDistances();

        // -------------------------------------------------
        // make neural network potential
        // ------------------------------------------------
        NeuralNetworkPotential nnp;
        nnp.readSetupFiles(potentialPath);

        // // calculate NNP energy and forces for specific atoms
        // for (int i=0; i<10; i++)
        // {
        //     // get specific atom in structure file
        //     Atom atom = structure.getAtom(i);

        //     // calculate NNP energy and force for specific atom
        //     nnp.calculateEnergy(structure, &atom);
        //     nnp.calculateForce(structure, &atom);

        //     // print out atom (with overwritten NNP energy and forces)
        //     cout << atom.toString() << endl;
        // }

        // Or, calculate NNP energy and forces for all atoms in given structure
        nnp.calculateTotalEnergy(structure);
        // nnp.calculateForce(structure);

        cout << "Total energy: " << structure.totalEnergy << endl;

        // write (calculated NNP energy and forces) into RuNNer data format
        // structure.writeFileFormatRunner(outStructureFile);
    }

    catch (runtime_error e)
    {
        cerr << "Runtime Error: " << e.what() << "!" << endl;
    }
}