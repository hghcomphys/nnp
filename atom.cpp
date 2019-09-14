//
// Atom
//

#include "atom.h"
#include <cstring>

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(int index, const char *element, const double position[3]): index(index)
{
    setElement(element);
    setPosition(position);
    fx = fy = fz = 0.0;
    charge = energy = 0.0;
}

Atom::Atom(int index, const char *element, const double position[3], const double force[3], double charge, double energy):
    Atom(index, element, position) 
{
    setForce(force);
    this->charge = charge;
    this->energy = energy;
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
    if (std::strcmp(this->element, element))
        return true;
    else 
        return false;
}

std::string Atom::toString()
{
    char buff[200];
    sprintf(buff, "Atom(index=%d, element=%s, position=(%f, %f, %f), force=(%f, %f, %f))",
        index, element, x, y, z, fx, fy, fz);
    return std::string(buff);
}

int Atom::getAtomicNumber(const char *element) 
{
    // TODO: extend to all elements
    if (strcmp(element, "H"))
        return 1;
    else if (strcmp(element, "C"))
        return 6; 
    else if (strcmp(element, "O"))
        return 8;
}
