//
// Atom
//

#include "atom.h"
#include <cstring>

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(int index, const char *element, const double position[3]): 
    index(index), fx(0), fy(0), fz(0), charge(0), energy(0)
{
    setElement(element);
    setPosition(position);
}

Atom::Atom(int index, const char *element, const double position[3], const double force[3], double charge, double energy):
    index(index), charge(charge), energy(energy)
{
    setElement(element);
    setPosition(position);
    setForce(force);
}

void Atom::setElement(const char *element)
{
    strcpy(this->element, element);
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

bool Atom::isElement(const char *element)
{
    if ( std::strcmp(this->element, element) )
        return true;
    else 
        return false;
}

std::string Atom::toString()
{
    char buff[ATOM_BUFSIZE];
    sprintf(buff, "Atom(index=%d, element=%s, position=(%f, %f, %f), force=(%f, %f, %f))",
        index, element, x, y, z, fx, fy, fz);
    return std::string(buff);
}
   
int Atom::getAtomicNumber(const std::string& elem) 
{
    // TODO: extend to all elements
    if ( elem == "H" )
        return 1;
    else if ( elem == "C" )
        return 6; 
    else if ( elem == "O" )
        return 8;
}
