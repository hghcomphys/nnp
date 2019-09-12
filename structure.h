//
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include "atom.h"
#include "distance.h"
#include <vector>
// #include <map>

class AtomicStructure {
public:
    AtomicStructure();
    ~AtomicStructure();
    int getNumberOfAtoms();
    int getNumberOfAtomsForElement(const char * element);
    void readFileFormatXYZ(const char * filename);
    void readFileFormatRuNNer(const char * filename);
    void readFileFormatRuNNer();
    void setCell(const double cell[9]);
    bool isPBC();
    double distance(Atom &atom_i, Atom &atom_j, double drij[3]);
    double distance(Atom &atom_i, Atom &atom_j);
    std::vector<int> getListOfAtomIndexForElement(const char* element); //TODO: improve design
    std::vector<int> getListOfAtomIndex(); // TODO: improve design
    inline Atom& getAtom(int index) { return *listOfAtoms[index]; }
    void calculateTableOfDistances();
    Distance** tableOfDistances;

// private:
    bool isAtom;
    bool isCell;
    double cell[9];
    int numberOfAtoms;
    Atom **listOfAtoms;
    void applyPBC(double &dx, double &dy, double &dz);
};

#endif //NNP_STRUCTURE_H
