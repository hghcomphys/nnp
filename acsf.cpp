//
// Atomic Centered Symmetry Functions (ACSF)
//

#include <vector>
#include "acsf.h"
#include "atoms.h"

/* ----------------------------------------------------------------------
   setup for ACSF
------------------------------------------------------------------------- */
ACSF::ACSF() {}

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

void ACSF::addTwoBodySymmetryFunction(TwoBodySymmetryFunction *symmetry_function) {
    listOfTwoBodySF.push_back(symmetry_function);
}

void ACSF::addThreeBodySymmetryFunction(ThreeBodySymmetryFunction *symmetry_function) {
    listOfThreeBodySF.push_back(symmetry_function);
}

int ACSF::getNumberOfTwoBodySF() const { return listOfTwoBodySF.size(); }

int ACSF::getNumberOfThreeBodySF() const { return listOfThreeBodySF.size(); }

int ACSF::getNumberOfSF() const { return getNumberOfTwoBodySF() + getNumberOfThreeBodySF(); }

TwoBodySymmetryFunction& ACSF::getTwoBodySF(const int index) const {
    // TODO: index error
    return *(listOfTwoBodySF[index]);
}

ThreeBodySymmetryFunction& ACSF::getThreeBodySF(const int index) const {
    // TODO: index error
    return *(listOfThreeBodySF[index]);
}

void ACSF::calculate(Atoms &configuration)
{
    // TODO: optimization
    int n_2b = listOfTwoBodySF.size();
    int n_3b = listOfThreeBodySF.size();

    std::vector<double> results(n_2b+n_3b);
    std::fill(results.begin(), results.end(), 0.0);
    // std::cout << results.size() << std::endl;

    for(Atom &atom_i: configuration.getAtoms()) {
        for(Atom &atom_j: configuration.getAtoms()) {

            /*Two-body symmetry functions*/
            if (atom_i.getIndex() == atom_j.getIndex()) continue;
            double rij = configuration.distance(atom_i, atom_j);

            for (int n=0; n<n_2b; n++) {
                results[n] += listOfTwoBodySF[n]->function(rij);
            }

            /*Three-body symmetry functions*/
            for(Atom &atom_k: configuration.getAtoms()) {

                if (atom_i.getIndex() == atom_k.getIndex()) continue;
                if (atom_j.getIndex() == atom_k.getIndex()) continue;

                double rik = configuration.distance(atom_i, atom_k);
                double rjk = configuration.distance(atom_j, atom_k);

                for (int n=0; n<n_3b; n++) {
                    results[n+n_2b] += listOfThreeBodySF[n]->function(rij, rik, rjk);
                }
            }
        }
    }

    // set results to values
    // TODO: improve assignment
    values.clear();
    for (auto it = results.begin(); it != results.end(); it++)
        values.push_back(*it);
}

std::vector<double>& ACSF::getValues() { return values; }