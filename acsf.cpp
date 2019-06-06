//
// Atomic Centered Symmetry Functions (ACSF)
//

#include <vector>
#include "acsf.h"
#include "symmetry_function.h"
#include "atoms.h"

/* ----------------------------------------------------------------------
   setup for ACSF
------------------------------------------------------------------------- */
ACSF::ACSF() {}

ACSF::~ACSF() {

    /* free the allocated memory for two-body symmetry fucntion*/
    for (auto *each: two_body_symmetry_functions)
        delete each;
    two_body_symmetry_functions.clear();

    /* free the allocated memory for three-body symmetry fucntion*/
    for (auto *each: three_body_symmetry_functions)
        delete each;
    three_body_symmetry_functions.clear();
}

void ACSF::addTwoBodySymmetryFunction(TwoBodySymmetryFunction *symmetry_function) {
    two_body_symmetry_functions.push_back(symmetry_function);
}

void ACSF::addThreeBodySymmetryFunction(ThreeBodySymmetryFunction *symmetry_function) {
    three_body_symmetry_functions.push_back(symmetry_function);
}

std::vector<double> ACSF::calculate(Atoms &configuration) 
{
    // TODO: optimization
    int n_2b = two_body_symmetry_functions.size();
    int n_3b = three_body_symmetry_functions.size();

    std::vector<double> results(n_2b+n_3b);
    std::fill(results.begin(), results.end(), 0.0);
    // std::cout << results.size() << std::endl;

    for(Atom &atom_i: configuration.atoms) {
        for(Atom &atom_j: configuration.atoms) {

            /*Two-body symmetry functions*/
            if (atom_i.index == atom_j.index) continue;
            double rij = configuration.distance(atom_i, atom_j);

            for (int n=0; n<n_2b; n++) {
                results[n] += two_body_symmetry_functions[n]->function(rij);
            }

            /*Three-body symmetry functions*/
            for(Atom &atom_k: configuration.atoms) { 

                if (atom_i.index == atom_k.index) continue;
                if (atom_j.index == atom_k.index) continue;

                double rik = configuration.distance(atom_i, atom_k);
                double rjk = configuration.distance(atom_j, atom_k);

                for (int n=0; n<n_3b; n++) {
                    results[n+n_2b] += three_body_symmetry_functions[n]->function(rij, rik, rjk);
                }
            }
        }
    }
    return results;
}
