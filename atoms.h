//
// Atomic Configuration
//

#ifndef NNP_ATOMS_H
#define NNP_ATOMS_H

#include <string>
#include <vector>


class Atom {
public:
    int index;
    double x, y, z;
    std::string element;
    Atom(double x, double y, double z, std::string element, int index);
};


class Atoms {
    bool is_atom;
    bool is_cell;
    void apply_pbc(double &dx, double &dy, double &dz);
public:
    std::vector<Atom> atoms;
    double cell[9];
    Atoms();
    ~Atoms();
    void read_xyz(std::string filename);
    void set_cell(double cell[]);
    double distance(Atom &atom_i, Atom &atom_j);
};


#endif //NNP_ATOMS_H
