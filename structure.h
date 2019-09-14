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
    void readFileFormatRuNNer(const char *filename);
    void readFileFormatRuNNer();
    int getNumberOfAtomsForElement(const char *element);
    inline Atom & getAtom(int atomIndex) { return *listOfAtoms[atomIndex]; }
    inline int getNumberOfAtoms() { return numberOfAtoms; }
    Atom **getListOfAtoms(); 
    Atom **getListOfAtomsForElement(const char *element);
    void setCell(const double cell[9]);
    double distance(Atom *atom_i, Atom *atom_j, double drij[3]);
    double distance(Atom *atom_i, Atom *atom_j);
    void calculateTableOfDistances(double globalCutOffRadius = 12.0);
    Distance **getTableOfDistances();
   
private:
    bool isAtom, isCell;
    bool isTableOfDistances;
    int numberOfAtoms;
    double cell[9];
    // TODO: merge together listOfAtoms and ListOfAtomsForElement
    Atom **listOfAtoms;
    std::map<std::string, int> numberOfAtomsForElement;
    std::map<std::string, Atom**> listOfAtomsForElement; 
    Distance **tableOfDistances;
    void applyPBC(double &dx, double &dy, double &dz);
};

#endif //NNP_STRUCTURE_H
