//
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include "atom.h"
#include <vector>

class AtomicStructure {
public:
    AtomicStructure();
    ~AtomicStructure();
    int getNumberOfAtoms();
    int getNumberOfAtomsForElement(const std::string& element);
    void readFileFormatXYZ(const std::string& filename);
    void readFileFormatRuNNer(const std::string& filename);
    void readFileFormatRuNNer();
    void setCell(double cell[9]);
    bool isPBC() const;
    double distance(Atom &atom_i, Atom &atom_j, double drij[3]);
    double distance(Atom &atom_i, Atom &atom_j);
    // const Atom& operator[] (unsigned int i) const;
    // std::vector<Atom*> getListOfAtoms();
    std::vector<int> getListOfAtomIndexForElement(const std::string &element);
    std::vector<int> getListOfAtomIndex();
    inline Atom& getAtom(int index) { return *listOfAtoms[index]; }

private:
    bool isAtom;
    bool isCell;
    double cell[9];
    int atomIndex; //TODO: improve the method for indexing atoms
    std::vector<Atom*> listOfAtoms;
    void addAtom(Atom* atom);
    void applyPBC(double &dx, double &dy, double &dz);
};


#endif //NNP_STRUCTURE_H
