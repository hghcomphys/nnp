#include <iostream>
#include <sstream>

#include "atoms.h"
#include "acsf.h"

using namespace std;
//using namespace NNP_SF;

int main()
{
     ACSF descriptor;
     // Descriptor des;
     descriptor.addTwoBodySymmetricFunction(new G0({6.0}));
     descriptor.addTwoBodySymmetricFunction( new G1({6.0, 2.3, 4.2}) );
     descriptor.addThreeBodySymmetricFunction( new G4({6.0, 0.01, 0.01, 2, 1}) );
     descriptor.addThreeBodySymmetricFunction( new G5({6.0, 0.01, 0.01, 2, 1}) );

     // cout << dsc.two_body_symmetric_functions.size() << endl;
//     cout << dsc.two_body_symmetric_functions[0]->function(2.6) << endl;
    // cout << dsc.two_body_symmetric_functions[1]->function(2.6) << endl;
    // cout << dsc.three_body_symmetric_functions[0]->function(2.6, 2.3, 5.6) << endl;
    // cout << dsc.three_body_symmetric_functions[1]->function(2.6, 2.3, 5.6) << endl;

//    G1 sf1({6.0, 2.3, 4.2});
//    cout << "G1 = " << sf1.function(2.6) << endl;
//
//    G4 sf4({6.0, 0.01, 0.01, 2, 1});
//    cout << "G4 = " << sf4.function(2.6, 2.3, 5.6) << endl;
//
//    G5 sf5({6.0, 0.01, 0.01, 2, 1});
//    cout << "G5 = " << sf5.function(2.6, 2.3, 5.6) << endl;

//    Atom atom(1, 2, 3, "X", 0);
//    cout << atom.element << " " << atom.x << " " << atom.z << " " << atom.z <<  " " << atom.index << endl;

    AtomicConfiguration configuration;
    configuration.read_xyz("ions.xyz");

//    for(auto &atom: ac.atoms) {
//        cout << atom.element << " " << atom.x << " " << atom.y << " " << atom.z <<  " " << atom.index << endl;
//         break;
//    }
//    cout << ac.atoms.size() << endl;

//    double cell[9] = {10, 0, 0, 0, 10, 0, 0, 0, 10};
//    ac.set_cell(cell);
//    for (int d=0; d<9; d++)
//        cout << ac.cell[d] << endl;

     for (auto &sf: descriptor.calculate(configuration))
          cout << sf << endl;
    // cout << ac.atoms.size() << endl;
    // cout << acsf.descriptors[0]->calculate(ac) << endl;
}