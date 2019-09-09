
//
// Atom
//

#ifndef NNP_ATOM_H
#define NNP_ATOM_H

#include <string>

class Atom {
public:
    Atom(int index, const std::string& element, double position[3]);
    Atom(int index, const std::string& element, double position[3], double force[3]);
    void setPosition(double position[3]);
    void setForce(double force[3]);
    static int getAtomicNumber(const std::string& element);
// private:
    int index;
    double x, y, z;
    double fx, fy, fz;
    std::string element;
};

#endif //NNP_ATOM_H