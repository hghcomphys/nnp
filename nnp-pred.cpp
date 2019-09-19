#include <iostream>
#include <sstream>
#include "structure.h"
#include "neural_network_potential.h"
#include "logger.h"

using namespace std;

int main(int argc,char* argv[])
{
    try {

        // set path to NNP potential files
        string potentialPath = "potentials/water/";

        // read input argument
        if(argc==1) {
            Log(WARN) << "Path to NNP potential directory has to be given!";
            Log(WARN) << " [./nnp-pred potentials/water] will be used instead"; 
        }    
        else
        {
            potentialPath = string(argv[1]) + "/";
            Log(INFO) << "NNP directory: " << potentialPath;
        }

        // -----------------------------------------------
        // make atomic structure and read data from a file
        // -----------------------------------------------
        AtomicStructure structure;
        structure.readFileFormatRuNNer("example/input.data");
        structure.calculateTableOfDistances();

        // -------------------------------------------------
        // make neural network potential
        // ------------------------------------------------
        NeuralNetworkPotential nnp;
        nnp.readSetupFiles(potentialPath);

        // calculate NNP energy and forces specific atoms 
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

        // calculate NNP energy and forces for all atoms in given structure 
        nnp.caculateTotalEnergy(structure); 
        nnp.calculateForce(structure);

        // write (calculated NNP energy and forces) into RuNNer data format
        structure.writeFileFormatRunner("example/input.data.pred");
    }

    catch (runtime_error e)
    {
        cout << "Runtime Error: " << e.what() << "!" << endl;
    }

}