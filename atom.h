
//
// Atom
//

#pragma once

#include <string>

class Atom {
public:
    Atom(double x, double y, double z, std::string element, int index);
    Atom(double x, double y, double z, std::string element, int index, double fx, double fy, double fz);
    static int getAtomicNumber(const std::string& element);
public: /*used to be private*/
    int index;
    double x, y, z;
    double fx, fy, fz;
    std::string element;
};
