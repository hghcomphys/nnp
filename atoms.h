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
    void addAtom(const Atom& atom);
    std::vector<Atom>& getListOfAtoms();
    int getNumberOfAtoms();
    void readFileFormatXYZ(std::string filename);
    void setCell(double cell[]);
    double distance(Atom &atom_i, Atom &atom_j);
    std::vector<int> getListOfIndexForElement(const std::string &element);
    // const Atom& operator[] (unsigned int i) const;
    bool isPBC();
private:
    bool isAtom;
    bool isCell;
    void applyPBC(double &dx, double &dy, double &dz);
    double cell[9];
    std::vector<Atom> listOfAtoms;
    int atomIndex; //TODO: improve the method for indexing atoms
};


#endif //NNP_ATOMS_H
