//
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include "atom.h"
#include "distance.h"
#include <vector>
#include <map>
#include <string>

class AtomicStructure {
public:
    AtomicStructure();
    ~AtomicStructure();
    void readFileFormatRuNNer(const char * filename);
    void readFileFormatRuNNer();
    void setCell(const double cell[9]);
    double distance(Atom &atom_i, Atom &atom_j, double drij[3]);
    double distance(Atom &atom_i, Atom &atom_j);
    void calculateTableOfDistances();
    void applyPBC(double & dx, double & dy, double & dz);
    inline Atom& getAtom(int index) { return *listOfAtoms[index]; }
    Atom ** getListOfAtoms(); 
    Atom ** getListOfAtomForElement(const char * element); 
    int getNumberOfAtomsForElement(const char * element);

private:
    bool isAtom;
    bool isCell;
    bool isTableOfDistances;
    int numberOfAtoms;
    double cell[9];
    Atom ** listOfAtoms;
    std::map<std::string, int> numberOfAtomsForElement;
    std::map<std::string, Atom**> listOfAtomsForElement; 
    Distance ** tableOfDistances;
    void addAtom(Atom * atom);
};

#endif //NNP_STRUCTURE_H
