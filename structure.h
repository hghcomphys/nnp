//
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include "atom.h"
#include "distance.h"
#include <map>
#include <string>

class AtomicStructure {
public:
    AtomicStructure();
    ~AtomicStructure();
    void writeFileFormatRunner(const std::string& filename);
    void readFileFormatRuNNer(const std::string& filename);
    void readFileFormatRuNNer();
    Atom& getAtom(int atomIndex);
    void setCell(const double cell[9]);
    double distance(Atom *atom_i, Atom *atom_j, double drij[3]);
    double distance(Atom *atom_i, Atom *atom_j);
    void calculateTableOfDistances(double globalCutOffRadius = 12.0);

// private:
    bool isAtom, isCell;
    bool isTableOfDistances;
    int numberOfAtoms;
    double cell[9];
    double totalCharge, totalEnergy;
    // TODO: merge together listOfAtoms and ListOfAtomsForElement
    Atom **atoms;
    std::map<std::string, int> numberOfAtomsForElement;
    std::map<std::string, Atom**> atomsForElement; 
    Distance **tableOfDistances;
    void applyPBC(double& dx, double& dy, double& dz);
};

inline Atom& AtomicStructure::getAtom(int atomIndex) 
{ 
    return *atoms[atomIndex]; 
}

#endif //NNP_STRUCTURE_H