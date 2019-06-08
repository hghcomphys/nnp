//
// Atomic Configuration
//

#ifndef NNP_ATOMS_H
#define NNP_ATOMS_H

#include <string>
#include <vector>


class Atom {
public:
    Atom(double x, double y, double z, std::string element, int index);
    double getX();
    double getY();
    double getZ();
    double getIndex();
    std::string getElement();
private:
    int index;
    double x, y, z;
    std::string element;
};


class Atoms {
public:
    Atoms();
    ~Atoms();
    std::vector<Atom>& getAtoms();
    int getNumberOfAtoms();
    void readXYZ(std::string filename);
    void setCell(double cell[]);
    double distance(Atom &atom_i, Atom &atom_j);
private:
    bool isAtom;
    bool isCell;
    void applyPBC(double &dx, double &dy, double &dz);
    double cell[9];
    std::vector<Atom> atoms;
};


#endif //NNP_ATOMS_H
