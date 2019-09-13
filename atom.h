
//
// Atom
//

#ifndef NNP_ATOM_H
#define NNP_ATOM_H

#define ELEMENT_BUFSIZE 16
#include <string>

class Atom {
public:
    Atom(int index, const char *element, const double position[3]);
    Atom(int index, const char *element, const double position[3], const double force[3]);
    void setPosition(const double position[3]);
    void setForce(const double force[3]);
    void setElement(const char *element);
    bool isElement(const char *element);
    static int getAtomicNumber(const char *element);
    std::string toString();
    
// private:
    int index;
    double x, y, z;
    double fx, fy, fz;
    char element[ELEMENT_BUFSIZE];
};

#endif //NNP_ATOM_H