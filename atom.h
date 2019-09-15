
//
// Atom
//

#ifndef NNP_ATOM_H
#define NNP_ATOM_H

#include <string>

class Atom {
public:
    Atom(int index, const std::string& element, const double position[3]);
    Atom(int index, const std::string& element, const double position[3], const double force[3], double charge=0.0, double energy=0.0);
    void setPosition(const double position[3]);
    void setForce(const double force[3]);
    bool isElement(const std::string& element);
    std::string toString();

// private:
    int index;
    double x, y, z;
    double fx, fy, fz;
    double charge, energy;
    std::string element;
};

int getAtomicNumber(const std::string& element);

#endif //NNP_ATOM_H