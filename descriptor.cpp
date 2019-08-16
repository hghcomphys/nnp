//
// Atomic Centered Symmetry Functions (ACSF)
//

#include <vector>
#include "descriptor.h"
#include <iostream>

#include "log.h"

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
    const auto& listOfIndexForCentralElement = structure.getListOfAtomIndexForElement(centralElement);
    
    // loop over all atoms
    for(int atomIndex: listOfIndexForCentralElement)
        values.push_back( calculate(structure, atomIndex) );

    return values;
}

double calculateCosine(double rij, double rik, double drij[3], double drik[3]) 
{
    // cosine of angle between k--<i>--j atoms
    double cost = 0;
    for (int d=0; d<3; d++)
        cost += drij[d] * drik[d];
    return cost / rij / rik;
}

std::vector<double> ACSF::calculate(AtomicStructure &structure, int atomIndex)
{
    // TODO: optimization
    const int n_2b = getNumberOfTwoBodySF();
    const int n_3b = getNumberOfThreeBodySF();

    std::vector<double> values(n_2b + n_3b);
    std::fill(values.begin(), values.end(), 0.0); // initialize to zero

    // central element
    Atom& atom_i = structure.getAtom(atomIndex);

    // list of index for all atoms
    const auto& listOfAtomIndex = structure.getListOfAtomIndex();

    // global cutoff radius
    const double rcGlobal = getGlobalCutOffRadius();
    // Log(INFO) << "Global cutoff radius: " << rcGlobal;

    // ===================================================================

    // loop over fist neighbors
    for(auto j: listOfAtomIndex) 
    {
        Atom& atom_j = structure.getAtom(j);
        if ( atom_j.index == atom_i.index ) continue;

        double drij[3];
        const double rij = structure.distance(atom_i, atom_j, drij);
        if ( rij > rcGlobal ) continue;

        // Loop over all two-body symmetry functions
        for (int n=0; n<n_2b; n++) 
        {
            if ( atom_j.element == listOfTwoBodyNeighborElement[n] )
            {
                values[n] += listOfTwoBodySF[n]->function(rij);
            }
        }

        // loop over second neighbors
        for(auto k: listOfAtomIndex) 
        {
            Atom& atom_k = structure.getAtom(k);
            if (atom_k.index == atom_i.index) continue;
            if (atom_k.index <= atom_j.index) continue;

            double drik[3];
            const double rik = structure.distance(atom_i, atom_k, drik);
            const double rjk = structure.distance(atom_j, atom_k);
            if ( rik > rcGlobal || rjk > rcGlobal ) continue;

            // Loop over all tree-body symmetry functions
            for (int n=0; n<n_3b; n++) 
            {
                if ( atom_j.element == listOfThreeBodyNeighborElement1[n] && 
                    atom_k.element == listOfThreeBodyNeighborElement2[n] )
                {
                    // cosine of angle between k--<i>--j atoms
                    const double cost = calculateCosine(rij, rik, drij, drik);
                
                    // add value of symmetry function
                    values[n+n_2b] += listOfThreeBodySF[n]->function(rij, rik, rjk, cost);
                } 
            }
        }
    } 

    // ===================================================================

    // // Loop over all two-body symmetry functions
    // for (int n=0; n<n_2b; n++) 
    // {
    //     const auto& listOfAtomIndexForElement = structure.getListOfAtomIndexForElement(listOfTwoBodyNeighborElement[n]);
    //     for(auto j: listOfAtomIndexForElement) 
    //     {
    //         Atom& atom_j = structure.getAtom(j);
    //         if (atom_j.index == atom_i.index) continue;

    //         const double rij = structure.distance(atom_i, atom_j);
    //         if ( rij > rcGlobal ) continue;

    //         values[n] += listOfTwoBodySF[n]->function(rij);
    //     } 
    // }
    
    // // loop over all tree-body symmetry functions
    // for (int n=0; n<n_3b; n++) 
    // {
    //     const auto& listOfIndexForElement1 = structure.getListOfAtomIndexForElement(listOfThreeBodyNeighborElement1[n]);
    //     const auto& listOfIndexForElement2 = structure.getListOfAtomIndexForElement(listOfThreeBodyNeighborElement2[n]);

    //     // first neighbors
    //     for(int j: listOfIndexForElement1) {
                
    //         Atom& atom_j = structure.getAtom(j);
    //         if (atom_j.index == atom_i.index) continue;  

    //         double drij[3];
    //         const double rij = structure.distance(atom_i, atom_j, drij);
    //         if ( rij > rcGlobal ) continue;

    //         // second neighbors
    //         for(int k: listOfIndexForElement2) 
    //         {
    //             Atom& atom_k = structure.getAtom(k);
    //             if (atom_k.index == atom_i.index) continue;
    //             if (atom_k.index <= atom_j.index) continue;

    //             double drik[3];
    //             const double rik = structure.distance(atom_i, atom_k, drik);
    //             const double rjk = structure.distance(atom_j, atom_k);
    //             if ( rik > rcGlobal || rjk > rcGlobal ) continue;

    //             // cosine of angle between k--<i>--j atoms
    //             const double cost = calculateCosine(rij, rik, drij, drik);
    //             values[n+n_2b] += listOfThreeBodySF[n]->function(rij, rik, rjk, cost);
    //         }
    //     }
    // }

    // ===================================================================
 
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

std::vector<std::vector<double>> ACSF::gradient(AtomicStructure &structure, int atomIndex_i, int atomIndex_ip) 
{
    // TODO: optimization
    const int n_2b = getNumberOfTwoBodySF();
    const int n_3b = getNumberOfThreeBodySF();

    std::vector<std::vector<double>> values(n_2b + n_3b, std::vector<double>(3));
    for (auto& value: values)
        std::fill(value.begin(), value.end(), 0.0); // initialize to zero

    // list of index of atom
    const auto& listOfAtomIndex = structure.getListOfAtomIndex();

    // atoms with index i & j
    Atom& atom_i = structure.getAtom(atomIndex_i);

    // Loop over all two-body symmetry functions
    for (int n=0; n<n_2b; n++) 
    { 
        // list of atom index for a specific element
        const auto& listOfAtomIndexForElement = structure.getListOfAtomIndexForElement(listOfTwoBodyNeighborElement[n]);
        
        // check whether gradient is respect to atom itself or neighbor atom
        if ( atomIndex_ip == atomIndex_i )
        {
            for (int j: listOfAtomIndexForElement) {

                if ( j == atomIndex_i ) continue;  
                Atom& atom_j = structure.getAtom(j);
                
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
                
                Atom& atom_j = structure.getAtom(atomIndex_ip);
                
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
        const auto& listOfAtomIndexForElement1 = structure.getListOfAtomIndexForElement(listOfThreeBodyNeighborElement1[n]);
        const auto& listOfAtomIndexForElement2 = structure.getListOfAtomIndexForElement(listOfThreeBodyNeighborElement2[n]);

        // check whether the gradient is respect to atom itself or other atoms
        if ( atomIndex_ip == atomIndex_i )
        {
            // first neighbors
            for(int j: listOfAtomIndexForElement1) {
                    
                Atom& atom_j = structure.getAtom(j);
                if (atom_j.index == atom_i.index) continue;  

                double drij[3];
                const double rij = structure.distance(atom_i, atom_j, drij);

                // second neighbors
                for(int k: listOfAtomIndexForElement2) {
                    
                    Atom& atom_k = structure.getAtom(k);
                    if (atom_k.index == atom_i.index) continue;
                    if (atom_k.index <= atom_j.index) continue;

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

                Atom& atom_j = structure.getAtom(atomIndex_ip);
                if (atom_j.index == atom_i.index) continue;  

                double drij[3];
                const double rij = structure.distance(atom_i, atom_j, drij);

                // second neighbors
                for(int k: listOfAtomIndexForElement2) {
                    
                    Atom& atom_k = structure.getAtom(k);
                    if (atom_k.index == atom_i.index) continue;
                    if (atom_k.index <= atom_j.index) continue;
                    
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
                
                    Atom& atom_j = structure.getAtom(j);
                    if (atom_j.index == atom_i.index) continue;  

                    double drij[3];
                    const double rij = structure.distance(atom_i, atom_j, drij);
  
                    Atom& atom_k = structure.getAtom(atomIndex_ip);
                    if (atom_k.index == atom_i.index) continue;
                    if (atom_k.index <= atom_j.index) continue;

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

double ACSF::getGlobalCutOffRadius() const 
{
    double maxValue = 0.0;

    // loop over two-body symmetry functions
    for (auto iter: listOfTwoBodySF)
        if( iter->getCutoffRadius() > maxValue )
            maxValue = iter->getCutoffRadius();

    // loop over three-body symmetry functions
    for (auto iter: listOfTwoBodySF)
        if( iter->getCutoffRadius() > maxValue )
            maxValue = iter->getCutoffRadius();

    return maxValue;
}