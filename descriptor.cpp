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

void ACSF::addTwoBodySF(TwoBodySymmetryFunction *symmetryFunction, const std::string& neighborElement1) 
{
    listOfTwoBodySF.push_back(symmetryFunction);
    listOfTwoBodyNeighborElement.push_back(neighborElement1);
}

void ACSF::addThreeBodySF(ThreeBodySymmetryFunction *symmetryFunction,
        const std::string& neighborElement1, const std::string& neighborElement2) 
{
    listOfThreeBodySF.push_back(symmetryFunction);
    listOfThreeBodyNeighborElement1.push_back(neighborElement1);
    listOfThreeBodyNeighborElement2.push_back(neighborElement2);
}

int ACSF::getNumberOfTwoBodySF() const 
{ 
    return listOfTwoBodySF.size(); 
}

int ACSF::getNumberOfThreeBodySF() const 
{ 
    return listOfThreeBodySF.size();
}

int ACSF::getTotalNumberOfSF() const 
{ 
    return getNumberOfTwoBodySF() + getNumberOfThreeBodySF(); 
}

const std::string& ACSF::getCentralElement() const 
{ 
    return centralElement; 
}

TwoBodySymmetryFunction& ACSF::getTwoBodySF(int index) const 
{
    // TODO: index error
    return *(listOfTwoBodySF[index]);
}

ThreeBodySymmetryFunction& ACSF::getThreeBodySF(int index) const 
{
    // TODO: index error
    return *(listOfThreeBodySF[index]);
}

std::vector<std::vector<double>> ACSF::calculate(AtomicStructure &structure)
{
    std::vector<std::vector<double>> values;
    Atom **listOfAtomsForCentralElement = structure.listOfAtomsForElement[centralElement.c_str()];
    const int numberOfAtomsForCentralElement = structure.numberOfAtomsForElement[centralElement.c_str()];
    
    // loop over all atoms
    for(int i=0; i<numberOfAtomsForCentralElement; i++)
        values.push_back( calculate(structure, listOfAtomsForCentralElement[i]->index) );

    return values;
}

std::vector<double> ACSF::calculate(AtomicStructure &structure, int atomIndex)
{
    // TODO: optimization
    const int n_2b = getNumberOfTwoBodySF();
    const int n_3b = getNumberOfThreeBodySF();

    // initialize values for array of symmetry functions
    std::vector<double> values(n_2b + n_3b);
    std::fill(values.begin(), values.end(), 0.0); // initialize to zero

    // get central atom i
    Atom *atom_i = &(structure.getAtom(atomIndex));

    // global cutoff radius
    const double rcGlobal = getGlobalCutOffRadius();

    // ===================================================================

    // // list of index for all atoms
    // // const auto& listOfAtomIndex = structure.getListOfAtomIndex();
    // const int nAtoms = structure.getNumberOfAtoms(); 

    // // loop over fist neighbors
    // for(int j=0; j<nAtoms; j++) 
    // {
    //     // get neighbor atom j
    //     const Atom& atom_j = structure.getAtom(j);
    //     if ( atom_j.index == atom_i.index ) continue;

    //     // get distance between atom i & j from table of distances
    //     const Distance& distance_ij = structure.tableOfDistances[atom_i.index][atom_j.index];
    //     if ( distance_ij.dr > rcGlobal ) continue;

    //     // Loop over all two-body symmetry functions
    //     for (int n=0; n<n_2b; n++) 
    //     {
    //         if ( atom_j.element == listOfTwoBodyNeighborElement[n] )
    //         {
    //             values[n] += listOfTwoBodySF[n]->function(distance_ij.dr);
    //         }
    //     }

    //     // loop over second neighbors
    //     for(int k=0; k<nAtoms; k++) 
    //     {
    //         const Atom& atom_k = structure.getAtom(k);
    //         if ( atom_k.index == atom_i.index ) continue;
    //         if ( atom_k.index <= atom_j.index ) continue;

    //         // get distance between atom i & k from table of distances
    //         const Distance& distance_ik = structure.tableOfDistances[atom_i.index][atom_k.index];
    //         if ( distance_ik.dr > rcGlobal ) continue;

    //         // get distance between atom j & k from table of distances
    //         const Distance& distance_jk = structure.tableOfDistances[atom_j.index][atom_k.index];
    //         if ( distance_jk.dr > rcGlobal ) continue;

    //         // Loop over all tree-body symmetry functions
    //         for (int n=0; n<n_3b; n++) 
    //         {
    //             if ( atom_j.element == listOfThreeBodyNeighborElement1[n] && 
    //                 atom_k.element == listOfThreeBodyNeighborElement2[n] )
    //             {
    //                 // cosine of angle between k--<i>--j atoms
    //                 const double cost = calculateCosine(distance_ij.dr, distance_ik.dr, distance_ij.drVec, distance_ik.drVec);
                    
    //                 // add value of symmetry function
    //                 values[n+n_2b] += listOfThreeBodySF[n]->function(distance_ij.dr, distance_ik.dr, distance_jk.dr, cost);
    //             } 
    //         }
    //     }
    // } 

    // ===================================================================

    // Loop over all two-body symmetry functions
    for (int n=0; n<n_2b; n++) 
    {
        Atom **listOfAtomsForNeighborElement = structure.listOfAtomsForElement[listOfTwoBodyNeighborElement[n].c_str()];
        const int numberOfAtomsForNeighborElement = structure.numberOfAtomsForElement[listOfTwoBodyNeighborElement[n].c_str()];
       
        for(int j=0; j<numberOfAtomsForNeighborElement; j++) 
        {
            Atom *atom_j = listOfAtomsForNeighborElement[j];
            if (atom_j->index == atom_i->index) continue;

            // const double rij = structure.distance(atom_i, atom_j);
            Distance& distance_ij = structure.tableOfDistances[atom_i->index][atom_j->index];
            if ( distance_ij.dr > rcGlobal ) continue;

            values[n] += listOfTwoBodySF[n]->function(distance_ij.dr);
        } 
    }
    
    // loop over all tree-body symmetry functions
    for (int n=0; n<n_3b; n++) 
    {
        Atom **listOfAtomsForNeighborElement1 = structure.listOfAtomsForElement[listOfThreeBodyNeighborElement1[n].c_str()];
        const int numberOfAtomsForNeighborElement1 = structure.numberOfAtomsForElement[listOfThreeBodyNeighborElement1[n].c_str()];
     
        Atom **listOfAtomsForNeighborElement2 = structure.listOfAtomsForElement[listOfThreeBodyNeighborElement2[n].c_str()];
        const int numberOfAtomsForNeighborElement2 = structure.numberOfAtomsForElement[listOfThreeBodyNeighborElement2[n].c_str()];
     
        // first neighbors
        for(int j=0; j<numberOfAtomsForNeighborElement1; j++) 
        {        
            Atom *atom_j = listOfAtomsForNeighborElement1[j];
            if (atom_j->index == atom_i->index) continue;  

            // double drij[3];
            // const double rij = structure.distance(atom_i, atom_j, drij);
            Distance& distance_ij = structure.tableOfDistances[atom_i->index][atom_j->index];
            if ( distance_ij.dr > rcGlobal ) continue;

            // second neighbors
            for(int k=0; k<numberOfAtomsForNeighborElement2; k++) 
            {
                Atom *atom_k = listOfAtomsForNeighborElement2[k];
                if (atom_k->index == atom_i->index) continue;
                if (atom_k->index <= atom_j->index) continue;

                // double drik[3];
                // const double rik = structure.distance(atom_i, atom_k, drik);
                Distance& distance_ik = structure.tableOfDistances[atom_i->index][atom_k->index];
                if ( distance_ik.dr > rcGlobal ) continue;

                // const double rjk = structure.distance(atom_j, atom_k);
                Distance& distance_jk = structure.tableOfDistances[atom_j->index][atom_k->index];
                if ( distance_jk.dr > rcGlobal ) continue;

                // cosine of angle between k--<i>--j atoms
                const double cost = calculateCosine(distance_ij.dr, distance_ik.dr, distance_ij.drVec, distance_ik.drVec);
                
                // add value of symmetry function
                values[n+n_2b] += listOfThreeBodySF[n]->function(distance_ij.dr, distance_ik.dr, distance_jk.dr, cost);
            }
        }
    }

    // ===================================================================
 
    // return array of symmetry function values
    return values;
}

inline bool isInList(const std::vector<int>& list, int item) 
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

    // // global cutoff radius
    // const double rcGlobal = getGlobalCutOffRadius();
    // // Log(INFO) << "Global cutoff radius: " << rcGlobal;

    // // list of index of atom
    // const auto& listOfAtomIndex = structure.getListOfAtomIndex();

    // // get central atom i
    // Atom& atom_i = structure.getAtom(atomIndex_i);

    // // Loop over all two-body symmetry functions
    // for (int n=0; n<n_2b; n++) 
    // { 
    //     // list of atom index for a specific element
    //     const auto& listOfAtomIndexForElement = structure.getListOfAtomIndexForElement(listOfTwoBodyNeighborElement[n]);
        
    //     // check whether gradient is respect to atom itself or neighbor atom
    //     if ( atomIndex_ip == atomIndex_i )
    //     {
    //         for (int j: listOfAtomIndexForElement) {

    //             if ( j == atomIndex_i ) continue;  
    //             Atom& atom_j = structure.getAtom(j);

    //             // get distance between atom i & j from table of distances
    //             Distance& distance_ij = structure.tableOfDistances[atom_i.index][j];
    //             // if ( distance_ij.dr > rcGlobal ) continue;

    //             const std::vector<double>& gradient = listOfTwoBodySF[n]->gradient_ii(distance_ij.dr, distance_ij.drVec);
    //             for (int d=0; d<3; d++)
    //                 values[n][d] += gradient[d]; 
    //         }
    //     } 
    //     else
    //     {
    //         if ( isInList(listOfAtomIndexForElement, atomIndex_ip) ) {
                
    //             Atom& atom_j = structure.getAtom(atomIndex_ip);
                
    //             // get distance between atom i & j from table of distances
    //             Distance& distance_ij = structure.tableOfDistances[atom_i.index][atom_j.index];
    //             if ( distance_ij.dr > rcGlobal ) continue;

    //             const std::vector<double>& gradient = listOfTwoBodySF[n]->gradient_ij(distance_ij.dr, distance_ij.drVec);
    //             for (int d=0; d<3; d++)
    //                 values[n][d] = gradient[d]; 
    //         }
    //     }
    // }

    // // Loop over all tree-body symmetry functions
    // for (int n=0; n<n_3b; n++) 
    // {
    //     // list of atom index for specific elements
    //     const auto& listOfAtomIndexForElement1 = structure.getListOfAtomIndexForElement(listOfThreeBodyNeighborElement1[n]);
    //     const auto& listOfAtomIndexForElement2 = structure.getListOfAtomIndexForElement(listOfThreeBodyNeighborElement2[n]);

    //     // check whether the gradient is respect to atom itself or other atoms
    //     if ( atomIndex_ip == atomIndex_i )
    //     {
    //         // first neighbors
    //         for(int j: listOfAtomIndexForElement1) {
                    
    //             Atom& atom_j = structure.getAtom(j);
    //             if (atom_j.index == atom_i.index) continue;  

    //             // get distance between atom i & j from table of distances
    //             Distance& distance_ij = structure.tableOfDistances[atom_i.index][atom_j.index];
    //             if ( distance_ij.dr > rcGlobal ) continue;

    //             // second neighbors
    //             for(int k: listOfAtomIndexForElement2) {
                    
    //                 Atom& atom_k = structure.getAtom(k);
    //                 if (atom_k.index == atom_i.index) continue;
    //                 if (atom_k.index <= atom_j.index) continue;

    //                 // get distance between atom i & k from table of distances
    //                 Distance& distance_ik = structure.tableOfDistances[atom_i.index][k];
    //                 if ( distance_ik.dr > rcGlobal ) continue;

    //                 // get distance between atom j & k from table of distances
    //                 Distance& distance_jk = structure.tableOfDistances[j][k];
    //                 if ( distance_jk.dr > rcGlobal ) continue;

    //                 // cosine of angle between k--<i>--j atoms
    //                 const double cost = calculateCosine(distance_ij.dr, distance_ik.dr, distance_ij.drVec, distance_ik.drVec);
            
    //                 // add gradient vector 
    //                 const std::vector<double>& gradient = listOfThreeBodySF[n]->gradient_ii(distance_ij.dr, distance_ik.dr, distance_jk.dr, 
    //                         cost, distance_ij.drVec, distance_ik.drVec, distance_jk.drVec);
    //                 for (int d=0; d<3; d++)
    //                     values[n+n_2b][d] += gradient[d];
    //             }
    //         }
    //     } 
    //     else
    //     {
    //         if ( isInList(listOfAtomIndexForElement1, atomIndex_ip) ) {

    //             Atom& atom_j = structure.getAtom(atomIndex_ip);
    //             if (atom_j.index == atom_i.index) continue;  

    //             // get distance between atom i & j from table of distances
    //             Distance& distance_ij = structure.tableOfDistances[atom_i.index][atom_j.index];
    //             if ( distance_ij.dr > rcGlobal ) continue;

    //             // second neighbors
    //             for(int k: listOfAtomIndexForElement2) {
                    
    //                 Atom& atom_k = structure.getAtom(k);
    //                 if (atom_k.index == atom_i.index) continue;
    //                 if (atom_k.index <= atom_j.index) continue;
                    
    //                 // get distance between atom i & k from table of distances
    //                 Distance& distance_ik = structure.tableOfDistances[atom_i.index][k];
    //                 if ( distance_ik.dr > rcGlobal ) continue;

    //                 // get distance between atom j & k from table of distances
    //                 Distance& distance_jk = structure.tableOfDistances[atom_j.index][k];
    //                 if ( distance_jk.dr > rcGlobal ) continue;

    //                 // cosine of angle between k--<i>--j atoms
    //                 const double cost = calculateCosine(distance_ij.dr, distance_ik.dr, distance_ij.drVec, distance_ik.drVec);
            
    //                 // add gradient vector 
    //                 const std::vector<double>& gradient = listOfThreeBodySF[n]->gradient_ij(distance_ij.dr, distance_ik.dr, distance_jk.dr, 
    //                         cost, distance_ij.drVec, distance_ik.drVec, distance_jk.drVec);
    //                 for (int d=0; d<3; d++)
    //                     values[n+n_2b][d] += gradient[d];
    //             }
    //         }

    //         if ( isInList(listOfAtomIndexForElement2, atomIndex_ip) ) {

    //             for(int j: listOfAtomIndexForElement1) {
                
    //                 Atom& atom_j = structure.getAtom(j);
    //                 if (atom_j.index == atom_i.index) continue;  

    //                 // get distance between atom i & j from table of distances
    //                 Distance& distance_ij = structure.tableOfDistances[atom_i.index][atom_j.index];
    //                 if ( distance_ij.dr > rcGlobal ) continue;
  
    //                 Atom& atom_k = structure.getAtom(atomIndex_ip);
    //                 if (atom_k.index == atom_i.index) continue;
    //                 if (atom_k.index <= atom_j.index) continue;

    //                 // get distance between atom i & k from table of distances
    //                 Distance& distance_ik = structure.tableOfDistances[atom_i.index][atom_k.index];
    //                 if ( distance_ik.dr > rcGlobal ) continue;

    //                 // get distance between atom j & k from table of distances
    //                 Distance& distance_jk = structure.tableOfDistances[j][atom_k.index];
    //                 if ( distance_jk.dr > rcGlobal ) continue;

    //                 // cosine of angle between k--<i>--j atoms
    //                 const double cost = calculateCosine(distance_ij.dr, distance_ik.dr, distance_ij.drVec, distance_ik.drVec);
            
    //                 // add gradient vector 
    //                 const std::vector<double>& gradient = listOfThreeBodySF[n]->gradient_ik(distance_ij.dr, distance_ik.dr, distance_jk.dr, 
    //                         cost, distance_ij.drVec, distance_ik.drVec, distance_jk.drVec);
    //                 for (int d=0; d<3; d++)
    //                     values[n+n_2b][d] += gradient[d];
    //             }
    //         }
    //     }   
    // }
 
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