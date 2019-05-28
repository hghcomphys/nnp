#include <iostream>
#include <sstream>

#include "symfn.h"
using namespace std;


int main()
{
    G0 sf0(6.0);
    cout << "G0 = " << sf0.function(2.6) << endl;

    G1 sf1(6.0, 2.3, 4.2);
    cout << "G1 = " << sf1.function(2.6) << endl;

    G4 sf4(6.0, 0.01, 0.01, 2, 1);
    cout << "G4 = " << sf4.function(2.6, 2.3, 5.6) << endl;

    G5 sf5(6.0, 0.01, 0.01, 2, 1);
    cout << "G5 = " << sf5.function(2.6, 2.3, 5.6) << endl;
}