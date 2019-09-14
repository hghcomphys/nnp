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
    void writeFileFormatRunner(const char *filename);
    void readFileFormatRuNNer(const char *filename);
    void readFileFormatRuNNer();
    int getNumberOfAtomsForElement(const char *element);
    Atom &getAtom(int atomIndex);
    int getNumberOfAtoms();
    Atom **getListOfAtoms(); 
    Atom **getListOfAtomsForElement(const char *element);
    void setCell(const double cell[9]);
    double distance(Atom *atom_i, Atom *atom_j, double drij[3]);
    double distance(Atom *atom_i, Atom *atom_j);
    void calculateTableOfDistances(double globalCutOffRadius = 12.0);
    Distance **getTableOfDistances();

// private:
    bool isAtom, isCell;
    bool isTableOfDistances;
    int numberOfAtoms;
    double cell[9];
    double totalCharge, totalEnergy;
    // TODO: merge together listOfAtoms and ListOfAtomsForElement
    Atom **listOfAtoms;
    std::map<std::string, int> numberOfAtomsForElement;
    std::map<std::string, Atom**> listOfAtomsForElement; 
    Distance **tableOfDistances;
    void applyPBC(double &dx, double &dy, double &dz);
};

inline Atom &AtomicStructure::getAtom(int atomIndex) 
{ 
    return *listOfAtoms[atomIndex]; 
}

inline int AtomicStructure::getNumberOfAtoms() 
{ 
    return numberOfAtoms; 
}

inline Distance **AtomicStructure::getTableOfDistances()
{
    return tableOfDistances;
}

inline Atom **AtomicStructure::getListOfAtomsForElement(const char *element)
{
    // TODO:: throw error when element is not found
    // return vector's pointer
    return listOfAtomsForElement[element];
}

inline Atom **AtomicStructure::getListOfAtoms()
{
    // return vector's pointer
    return listOfAtoms;
}

inline int  AtomicStructure::getNumberOfAtomsForElement(const char *element)
{
    return numberOfAtomsForElement[element];
}

#endif //NNP_STRUCTURE_H