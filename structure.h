//
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include <string>
#include <vector>

class Atom {
public:
    Atom(double x, double y, double z, std::string element, int index);
    Atom(double x, double y, double z, std::string element, int index, double fx, double fy, double fz);
    static int getAtomicNumber(const std::string& element);
    // double getX();
    // double getY();
    // double getZ();
    // double getIndex();
    // std::string getElement(); 
    // double getFx() const;
    // double getFy() const;
    // double getFz() const; 

public:
    int index;
    double x, y, z;
    double fx, fy, fz;
    std::string element;
};


class AtomicStructure {
public:
    AtomicStructure();
    ~AtomicStructure();
    void addAtom(const Atom& atom);
    std::vector<Atom>& getListOfAtoms();
    int getNumberOfAtoms();
    int getNumberOfAtomsForElement(const std::string& element);
    void readFileFormatXYZ(const std::string& filename);
    void readFileFormatRuNNer(const std::string& filename);
    void readFileFormatRuNNer();
    void setCell(double cell[]);
    double distance(Atom &atom_i, Atom &atom_j, double drij[3]);
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


#endif //NNP_STRUCTURE_H
