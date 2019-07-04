//
// Neural Network
//

#include "neuralnetwork.h"
#include <fstream>
#include <sstream>

#include <iostream>
using namespace std;

#include "opennn.h"
// using namespace OpenNN;

/* ----------------------------------------------------------------------
   setup for class NeuralNetwork 
------------------------------------------------------------------------- */

NeuralNetwork::NeuralNetwork() {

        OpenNN::MultilayerPerceptron mlp(1, 3, 1);
        cout << "number of inputs: " <<mlp.get_inputs_number() << endl;
        cout << "number of outputs: " <<mlp.get_outputs_number() << endl;
        cout << mlp.get_outputs_number() << endl;

        OpenNN::Vector<double> input(1);
        input[0] = 0.5;
        OpenNN::Vector<double> output = mlp.calculate_outputs( input );
        cout << "output = " << output << endl;

        OpenNN::Matrix<double> jacobian = mlp.calculate_Jacobian ( input );
        cout << "jacobian = " << jacobian << endl;

}

NeuralNetwork::~NeuralNetwork() {}
