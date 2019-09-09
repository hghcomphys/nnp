//
// Atom
//

#include "atom.h"

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(int index, const std::string& element, double position[3]): index(index), element(element)
{
    setPosition(position);
}

Atom::Atom(int index, const std::string& element, double position[3], double force[3]):
   Atom(index, element, position) 
{
    setForce(force);
}

void Atom::setPosition(double position[3])
{
    x = position[0];
    y = position[1];
    z = position[2];
}

void Atom::setForce(double force[3])
{
    fx = force[0];
    fy = force[1];
    fz = force[2];
}

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
