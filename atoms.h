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


class AtomicConfiguration {
private:
    bool is_atom;
    bool is_cell;

public:
    std::vector<Atom> atoms;
    double cell[9];
    AtomicConfiguration();
    void read_xyz(std::string filename);
    void set_cell(double cell[]);
};


#endif //NNP_ATOMS_H
