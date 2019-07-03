//
// Atomic Centered Symmetry Functions (ACSF)
//

#include <vector>
#include "acsf.h"
#include "atoms.h"
#include <iostream>

/* ----------------------------------------------------------------------
   setup for ACSF
------------------------------------------------------------------------- */
ACSF::ACSF(std::string element): centralElement(element) {}

ACSF::~ACSF()
{
    /* free the allocated memory for two-body symmetry fucntion*/
    for (TwoBodySymmetryFunction *each: listOfTwoBodySF)
        delete each;
    listOfTwoBodySF.clear();

    /* free the allocated memory for three-body symmetry fucntion*/
    for (ThreeBodySymmetryFunction *each: listOfThreeBodySF)
        delete each;
    listOfThreeBodySF.clear();
}

void ACSF::addTwoBodySF(TwoBodySymmetryFunction *symmetryFunction, const std::string& neighborElement1) {
    listOfTwoBodySFindex.push_back(getTotalNumberOfSF());
    listOfTwoBodySF.push_back(symmetryFunction);
    listOfTwoBodyNeighborElement.push_back(neighborElement1);
}

void ACSF::addThreeBodySF(ThreeBodySymmetryFunction *symmetryFunction,
                            const std::string& neighborElement1, const std::string& neighborElement2) {
    listOfThreeBodySFindex.push_back(getTotalNumberOfSF());
    listOfThreeBodySF.push_back(symmetryFunction);
    listOfThreeBodyNeighborElement1.push_back(neighborElement1);
    listOfThreeBodyNeighborElement2.push_back(neighborElement2);
    
}

int ACSF::getNumberOfTwoBodySF() const { return listOfTwoBodySF.size(); }

int ACSF::getNumberOfThreeBodySF() const { return listOfThreeBodySF.size(); }

int ACSF::getTotalNumberOfSF() const { return getNumberOfTwoBodySF() + getNumberOfThreeBodySF(); }

TwoBodySymmetryFunction& ACSF::getTwoBodySF(int index) const {
    // TODO: index error
    return *(listOfTwoBodySF[index]);
}

ThreeBodySymmetryFunction& ACSF::getThreeBodySF(int index) const {
    // TODO: index error
    return *(listOfThreeBodySF[index]);
}

void ACSF::calculate(Atoms &configuration)
{
    // TODO: optimization
    const int n_2b = listOfTwoBodySF.size();
    const int n_3b = listOfThreeBodySF.size();
    const int nrow =  configuration.getNumberOfAtomsForElement(centralElement);
    const int ncol =  n_2b + n_3b;

    values.resize(nrow, std::vector<double>(ncol));
    for (auto& col: values) // initialize to zero
        std::fill(col.begin(), col.end(), 0.0);

    int elementIndex = 0;
    auto atoms = configuration.getListOfAtoms();
    auto listOfIndexForCentralElement = configuration.getListOfIndexForElement(centralElement);
    for(int i: listOfIndexForCentralElement) 
    {
        // central element
        Atom& atom_i = atoms[i];

        // symmetry function values for each atom with central element type
        auto& elementValues = values[elementIndex++];

        // Loop over all two-body symmetry functions
        for (int n=0; n<n_2b; n++) 
        {
            for(int j: configuration.getListOfIndexForElement(listOfTwoBodyNeighborElement[n])) {
                    Atom& atom_j = atoms[j];
                    if (atom_j.getIndex() == atom_i.getIndex()) continue;
                    const double rij = configuration.distance(atom_i, atom_j);
                    elementValues[listOfTwoBodySFindex[n]] += listOfTwoBodySF[n]->function(rij);
                } 
        }

        // Loop over all tree-body symmetry functions
        for (int n=0; n<n_3b; n++) 
        {
            // first neighbors
            for(int j: configuration.getListOfIndexForElement(listOfThreeBodyNeighborElement1[n])) {
                    
                    Atom& atom_j = atoms[j];
                    if (atom_j.getIndex() == atom_i.getIndex()) continue;
                    
                    double drij[3];
                    const double rij = configuration.distance(atom_i, atom_j, drij);

                    // second neighbors
                    for(int k: configuration.getListOfIndexForElement(listOfThreeBodyNeighborElement2[n])) {
                        
                        Atom& atom_k = atoms[k];
                        if (atom_k.getIndex() == atom_i.getIndex()) continue;

                        double drik[3];
                        const double rik = configuration.distance(atom_i, atom_k, drik);
                        const double rjk = configuration.distance(atom_j, atom_k);

                        // cosine of angle between k--<i>--j atoms
                        double cost = 0;
                        for (int d=0; d<3; d++)
                            cost += drij[d] * drik[d];
                        const double inv_r = 1.0 / rij / rik;
                        cost *= inv_r;
                        // std::cout << cost << std::endl;

                        elementValues[listOfThreeBodySFindex[n]] += listOfThreeBodySF[n]->function(rij, rik, rjk, cost);
                    }
            }
        }
    }
}

std::vector<std::vector<double>>& ACSF::getValues() { return values; }

std::string ACSF::getCentralElement() { return centralElement; }
