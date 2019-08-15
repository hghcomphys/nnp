//
// Atomic Centered Symmetry Functions (ACSF)
//

#include <vector>
#include "descriptor.h"
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
    listOfTwoBodySF.push_back(symmetryFunction);
    listOfTwoBodyNeighborElement.push_back(neighborElement1);
}

void ACSF::addThreeBodySF(ThreeBodySymmetryFunction *symmetryFunction,
        const std::string& neighborElement1, const std::string& neighborElement2) {
    listOfThreeBodySF.push_back(symmetryFunction);
    listOfThreeBodyNeighborElement1.push_back(neighborElement1);
    listOfThreeBodyNeighborElement2.push_back(neighborElement2);
}

int ACSF::getNumberOfTwoBodySF() const { return listOfTwoBodySF.size(); }

int ACSF::getNumberOfThreeBodySF() const { return listOfThreeBodySF.size(); }

int ACSF::getTotalNumberOfSF() const { return getNumberOfTwoBodySF() + getNumberOfThreeBodySF(); }

const std::string& ACSF::getCentralElement() const { return centralElement; }


TwoBodySymmetryFunction& ACSF::getTwoBodySF(int index) const {
    // TODO: index error
    return *(listOfTwoBodySF[index]);
}

ThreeBodySymmetryFunction& ACSF::getThreeBodySF(int index) const {
    // TODO: index error
    return *(listOfThreeBodySF[index]);
}

std::vector<std::vector<double>> ACSF::calculate(AtomicStructure &structure)
{
    std::vector<std::vector<double>> values;
    const auto listOfIndexForCentralElement = structure.getListOfIndexForElement(centralElement);
    
    // loop over all atoms
    for(auto atomIndex: listOfIndexForCentralElement)
        values.push_back( calculate(structure, atomIndex) );

    return values;
}

std::vector<double> ACSF::calculate(AtomicStructure &structure, int atomIndex)
{
    // TODO: optimization
    const int n_2b = getNumberOfTwoBodySF();
    const int n_3b = getNumberOfThreeBodySF();

    std::vector<double> values(n_2b + n_3b);
    std::fill(values.begin(), values.end(), 0.0); // initialize to zero

    auto atoms = structure.getListOfAtoms();

    // central element
    Atom& atom_i = atoms[atomIndex];

    // Loop over all two-body symmetry functions
    for (int n=0; n<n_2b; n++) 
    {
        for(int j: structure.getListOfIndexForElement(listOfTwoBodyNeighborElement[n])) {
                Atom& atom_j = atoms[j];
                if (atom_j.getIndex() == atom_i.getIndex()) continue;
                const double rij = structure.distance(atom_i, atom_j);
                values[n] += listOfTwoBodySF[n]->function(rij);
            } 
    }

    // Loop over all tree-body symmetry functions
    for (int n=0; n<n_3b; n++) 
    {
        // first neighbors
        for(int j: structure.getListOfIndexForElement(listOfThreeBodyNeighborElement1[n])) {
                
            Atom& atom_j = atoms[j];
            if (atom_j.getIndex() == atom_i.getIndex()) continue;  

            double drij[3];
            const double rij = structure.distance(atom_i, atom_j, drij);

            // second neighbors
            for(int k: structure.getListOfIndexForElement(listOfThreeBodyNeighborElement2[n])) {
                
                Atom& atom_k = atoms[k];
                if (atom_k.getIndex() == atom_i.getIndex()) continue;
                if (atom_k.getIndex() <= atom_j.getIndex()) continue;

                double drik[3];
                const double rik = structure.distance(atom_i, atom_k, drik);
                const double rjk = structure.distance(atom_j, atom_k);

                // cosine of angle between k--<i>--j atoms
                double cost = 0;
                for (int d=0; d<3; d++)
                    cost += drij[d] * drik[d];
                const double inv_r = 1.0 / rij / rik;
                cost *= inv_r;

                values[n+n_2b] += listOfThreeBodySF[n]->function(rij, rik, rjk, cost);
            }
        }
    }
 
    // return array of symmetry function values
    return values;
}

bool isInList(const std::vector<int>& list, int item) 
{
    bool find = false;
    for(auto each: list)
        if ( item == each ) {
            find = true;
            break;
        }
    return find;
}

double calculateCosine(double rij, double rik, double drij[3], double drik[3]) 
{
    // cosine of angle between k--<i>--j atoms
    double cost = 0;
    for (int d=0; d<3; d++)
        cost += drij[d] * drik[d];
    return cost / rij / rik;
}

std::vector<std::vector<double>> ACSF::gradient(AtomicStructure &structure, int atomIndex_i, int atomIndex_ip) 
{
    // TODO: optimization
    const int n_2b = getNumberOfTwoBodySF();
    const int n_3b = getNumberOfThreeBodySF();

    std::vector<std::vector<double>> values(n_2b + n_3b, std::vector<double>(3));
    for (auto& value: values)
        std::fill(value.begin(), value.end(), 0.0); // initialize to zero

    // list of atoms in structural configuration
    auto atoms = structure.getListOfAtoms();

    // atoms with index i & j
    Atom& atom_i = atoms[atomIndex_i];

    // Loop over all two-body symmetry functions
    for (int n=0; n<n_2b; n++) 
    { 
        // list of atom index for a specific element
        const std::vector<int>& listOfAtomIndexForElement = structure.getListOfIndexForElement(listOfTwoBodyNeighborElement[n]);
        
        // check whether gradient is respect to atom itself or neighbor atom
        if ( atomIndex_ip == atomIndex_i )
        {
            for (int j: listOfAtomIndexForElement) {

                if ( j == atomIndex_i ) continue;  
                Atom& atom_j = atoms[j];
                
                double drij[3];
                const double rij = structure.distance(atom_i, atom_j, drij);

                const std::vector<double>& gradient = listOfTwoBodySF[n]->gradient_ii(rij, drij);
                for (int d=0; d<3; d++)
                    values[n][d] += gradient[d]; 
            }
        } 
        else
        {
            if ( isInList(listOfAtomIndexForElement, atomIndex_ip) ) {
                
                Atom& atom_j = atoms[atomIndex_ip];
                
                double drij[3];
                const double rij = structure.distance(atom_i, atom_j, drij);

                const std::vector<double>& gradient = listOfTwoBodySF[n]->gradient_ij(rij, drij);
                for (int d=0; d<3; d++)
                    values[n][d] = gradient[d]; 
            }
        }
    }

    // Loop over all tree-body symmetry functions
    for (int n=0; n<n_3b; n++) 
    {
        // list of atom index for specific elements
        const std::vector<int>& listOfAtomIndexForElement1 = structure.getListOfIndexForElement(listOfThreeBodyNeighborElement1[n]);
        const std::vector<int>& listOfAtomIndexForElement2 = structure.getListOfIndexForElement(listOfThreeBodyNeighborElement2[n]);

        // check whether the gradient is respect to atom itself or other atoms
        if ( atomIndex_ip == atomIndex_i )
        {
            // first neighbors
            for(int j: listOfAtomIndexForElement1) {
                    
                Atom& atom_j = atoms[j];
                if (atom_j.getIndex() == atom_i.getIndex()) continue;  

                double drij[3];
                const double rij = structure.distance(atom_i, atom_j, drij);

                // second neighbors
                for(int k: listOfAtomIndexForElement2) {
                    
                    Atom& atom_k = atoms[k];
                    if (atom_k.getIndex() == atom_i.getIndex()) continue;
                    if (atom_k.getIndex() <= atom_j.getIndex()) continue;

                    double drik[3], drjk[3];
                    const double rik = structure.distance(atom_i, atom_k, drik);
                    const double rjk = structure.distance(atom_j, atom_k, drjk);

                    // cosine of angle between k--<i>--j atoms
                    const double cost = calculateCosine(rij, rik, drij, drik);

                    const std::vector<double>& gradient = listOfThreeBodySF[n]->gradient_ii(rij, rik, rjk, cost, drij, drik, drjk);
                    for (int d=0; d<3; d++)
                        values[n+n_2b][d] += gradient[d];
                }
            }
        } 
        else
        {
            if ( isInList(listOfAtomIndexForElement1, atomIndex_ip) ) {

                Atom& atom_j = atoms[atomIndex_ip];
                if (atom_j.getIndex() == atom_i.getIndex()) continue;  

                double drij[3];
                const double rij = structure.distance(atom_i, atom_j, drij);

                // second neighbors
                for(int k: listOfAtomIndexForElement2) {
                    
                    Atom& atom_k = atoms[k];
                    if (atom_k.getIndex() == atom_i.getIndex()) continue;
                    if (atom_k.getIndex() <= atom_j.getIndex()) continue;
                    
                    double drik[3], drjk[3];
                    const double rik = structure.distance(atom_i, atom_k, drik);
                    const double rjk = structure.distance(atom_j, atom_k, drjk);

                    // cosine of angle between k--<i>--j atoms
                    const double cost = calculateCosine(rij, rik, drij, drik);

                    const std::vector<double>& gradient = listOfThreeBodySF[n]->gradient_ij(rij, rik, rjk, cost, drij, drik, drjk);
                    for (int d=0; d<3; d++)
                        values[n+n_2b][d] += gradient[d];
                }
            }

            if ( isInList(listOfAtomIndexForElement2, atomIndex_ip) ) {

                for(int j: listOfAtomIndexForElement1) {
                
                    Atom& atom_j = atoms[j];
                    if (atom_j.getIndex() == atom_i.getIndex()) continue;  

                    double drij[3];
                    const double rij = structure.distance(atom_i, atom_j, drij);
  
                    Atom& atom_k = atoms[atomIndex_ip];
                    if (atom_k.getIndex() == atom_i.getIndex()) continue;
                    if (atom_k.getIndex() <= atom_j.getIndex()) continue;

                    double drik[3], drjk[3];
                    const double rik = structure.distance(atom_i, atom_k, drik);
                    const double rjk = structure.distance(atom_j, atom_k, drjk);

                    // cosine of angle between k--<i>--j atoms
                    const double cost = calculateCosine(rij, rik, drij, drik);

                    const std::vector<double>& gradient = listOfThreeBodySF[n]->gradient_ik(rij, rik, rjk, cost, drij, drik, drjk);
                    for (int d=0; d<3; d++)
                        values[n+n_2b][d] += gradient[d];
                }
            }
        }   
    }
 
    // return array of gradient (vector) of symmetry functions
    return values;
}

