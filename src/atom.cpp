//
// Atom
//

#include "atom.h"
#define ATOM_BUFSIZE 200

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(int index, const std::string& element, const double position[3]): 
    index(index), element(element), fx(0), fy(0), fz(0), charge(0), energy(0)
{
    setPosition(position);
}

Atom::Atom(int index, const std::string& element, const double position[3], const double force[3], double charge, double energy): 
    index(index), element(element), charge(charge), energy(energy)
{
    setPosition(position);
    setForce(force);
}

void Atom::setPosition(const double position[3])
{
    x = position[0];
    y = position[1];
    z = position[2];
}

void Atom::setForce(const double force[3])
{
    fx = force[0];
    fy = force[1];
    fz = force[2];
}

bool Atom::isElement(const std::string& element)
{
    if ( this->element == element )
        return true;
    else 
        return false;
}

std::string Atom::toString()
{
    char buff[ATOM_BUFSIZE];
    sprintf(buff, "Atom(index=%d, element=%s, position=(%f, %f, %f), force=(%f, %f, %f), energy=%f)",
        index, element.c_str(), x, y, z, fx, fy, fz, energy);
    return std::string(buff);
}
   
int getAtomicNumber(const std::string& element) 
{
    // TODO: extend to all elements
    if ( element == "H" )
        return 1;
    else if ( element == "C" )
        return 6; 
    else if ( element == "O" )
        return 8;
}
