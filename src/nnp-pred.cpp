#include <iostream>
#include <sstream>
#include "structure.h"
#include "neural_network_potential.h"
#include "logger.h"

using namespace std;

int main(int argc,char* argv[])
{
    try {

        // set defaults
        string potentialPath = "potentials/water/";
        string inStructureFile = "example/input.data";
        
        // read input argument
        switch (argc) {
            case 1:
                Log(WARN) << "Expected (1) path to NNP potential and (2) input filename!";
                Log(WARN) << "'./nnp-pred potentials/water example/input.data' will be used instead"; 
                break;
            default:
                potentialPath = string(argv[1]) + "/";
                Log(INFO) << "NNP directory: " << potentialPath;

                 string inStructureFile = "example/input.data";
                Log(INFO) << "Input: " << inStructureFile;
        }

        string outStructureFile = inStructureFile + ".pred";
        Log(INFO) << "Output: " << outStructureFile;

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

        // calculate NNP energy and forces for specific atoms 
        for (int i=0; i<10; i++) 
        {
            // get specific atom in structure file
            Atom atom = structure.getAtom(i);

            // calculate NNP energy and force for specific atom
            nnp.calculateEnergy(structure, &atom); 
            nnp.calculateForce(structure, &atom);

            // print out atom (with overwritten NNP energy and forces)
            cout << atom.toString() << endl;
        }       

        // Or, calculate NNP energy and forces for all atoms in given structure 
        nnp.caculateTotalEnergy(structure); 
        nnp.calculateForce(structure);

        // write (calculated NNP energy and forces) into RuNNer data format
        structure.writeFileFormatRunner(outStructureFile);
    }

    catch (runtime_error e)
    {
        cout << "Runtime Error: " << e.what() << "!" << endl;
    }

}