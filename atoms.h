//
// Created by hossein on 5/30/19.
//

#ifndef NNP_ATOMS_H
#define NNP_ATOMS_H

#include <string>
#include <vector>


class Atom {
public:
    double x, y, z;
    std::string element;
    Atom(double x, double y, double z, std::string element);
};


class Atoms {
public:
    std::vector<Atom> atoms;
    void read_xyz(std::string filename);
};


#endif //NNP_ATOMS_H
