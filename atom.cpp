//
// Atom
//

#include "atom.h"

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(double x, double y, double z, std::string element, int index): 
    x(x), y(y), z(z), element(element), index(index), fx(0.0), fy(0.0), fz(0.0) {}

Atom::Atom(double x, double y, double z, std::string element, int index, double fx, double fy, double fz):
    x(x), y(y), z(z), element(element), index(index), fx(fx), fy(fy), fz(fz) {}

int Atom::getAtomicNumber(const std::string& element) 
{
    // TODO: extend to all elements
    if(element == "H")
        return 1;
    else if (element == "C")
        return 6; 
    else if (element == "O")
        return 8;
}
