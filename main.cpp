#include <iostream>
#include <sstream>

#include "symfn.h"
#include "atoms.h"

using namespace std;
//using namespace NNP_SF;

int main()
{
    ACSF acsf;

    //    G0 sf0({6.0});
    acsf.add(TG0, {6.0});
    acsf.add(TG1, {6.0, 2.3, 4.2});
    acsf.add(TG4, {6.0, 0.01, 0.01, 2, 1});
    acsf.add(TG5, {6.0, 0.01, 0.01, 2, 1});

    cout << acsf.list_of_symmetric_functions.size() << endl;
    cout << acsf.list_of_symmetric_functions[0]->descriptor(2.6) << endl;
    cout << acsf.list_of_symmetric_functions[1]->descriptor(2.6) << endl;
    cout << acsf.list_of_symmetric_functions[2]->descriptor(2.6, 2.3, 5.6) << endl;
    cout << acsf.list_of_symmetric_functions[3]->descriptor(2.6, 2.3, 5.6) << endl;

//    G1 sf1({6.0, 2.3, 4.2});
//    cout << "G1 = " << sf1.function(2.6) << endl;
//
//    G4 sf4({6.0, 0.01, 0.01, 2, 1});
//    cout << "G4 = " << sf4.function(2.6, 2.3, 5.6) << endl;
//
//    G5 sf5({6.0, 0.01, 0.01, 2, 1});
//    cout << "G5 = " << sf5.function(2.6, 2.3, 5.6) << endl;


//    AtomicConfiguration ac;
//    ac.read_xyz("water.xyz");
//    for(auto &atom: ac.atoms) {
//        cout << atom.element << " " << atom.x << " " << atom.y << " " << atom.z << endl;
//    }
//    cout << ac.atoms.size() << endl;
//
//    double cell[9] = {10, 0, 0, 0, 10, 0, 0, 0, 10};
//    ac.set_cell(cell);
//    for (int d=0; d<9; d++)
//        cout << ac.cell[d] << endl;

}