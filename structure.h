//
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include "atom.h"
#include <vector>

class Distance {
public:
    double dr;
    double drVec[3];
    double inv_dr;
    Distance();
    Distance(double r, double vec[3]); 
    void set_drVec(double vec[3], double factor=1.0);
    void set(double r, double vec[3], double factor=1.0);
};


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
    std::vector<int> getListOfAtomIndexForElement(const std::string &element); //TODO: improve design
    std::vector<int> getListOfAtomIndex(); // TODO: improve design
    inline Atom& getAtom(int index) { return *listOfAtoms[index]; }
    void calculateTableOfDistances();
    Distance** tableOfDistances;

private:
    bool isAtom;
    bool isCell;
    double cell[9];
    int atomIndex; //TODO: improve the method for indexing atoms
    std::vector<Atom*> listOfAtoms;
    void addAtom(Atom* atom);
    void applyPBC(double &dx, double &dy, double &dz);
};

inline void Distance::set_drVec(double vec[3], double factor)
{
    for (int i=0; i<3; i++)
        drVec[i] = vec[i] * factor;
}


#endif //NNP_STRUCTURE_H
