//
// Atomic Centered Symmetric Functions
//

#include "acsf.h"
#include "atoms.h"
#include <cmath>
#include <vector>
#include <iostream>

//using namespace NNP_SF;

/* ----------------------------------------------------------------------
   setup for base symmetric function base
------------------------------------------------------------------------- */

SymmetricFunction::SymmetricFunction(double cutoff_radius): cutoff_radius(cutoff_radius) {}

double SymmetricFunction::cutoff_function(double r) {
    return ( cos(M_PI*r/cutoff_radius) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetric function
------------------------------------------------------------------------- */

G0::G0(std::vector<double> p): TwoBodySymmetricFunction(p[0]) {}

double G0::function(double rij) {
    return cutoff_function(rij);
}

/* ----------------------------------------------------------------------
   setup for G1 symmetric function
------------------------------------------------------------------------- */

G1::G1(std::vector<double> p): TwoBodySymmetricFunction(p[0]), eta(p[1]), rs(p[2]) {}

double G1::function(double rij) {
    return exp( -eta * (rij-rs) * (rij-rs) ) * cutoff_function(rij);
}


/* ----------------------------------------------------------------------
   setup for G4 symmetric function
------------------------------------------------------------------------ */

G4::G4(std::vector<double> p): ThreeBodySymmetricFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G4::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoff_function(rij) * cutoff_function(rjk) * cutoff_function(rik);
}

/* ----------------------------------------------------------------------
   setup for G5 symmetric function
------------------------------------------------------------------------- */

G5::G5(std::vector<double> p): ThreeBodySymmetricFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G5::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoff_function(rij) * cutoff_function(rjk);
}

/* ----------------------------------------------------------------------
   setup for ACSF
------------------------------------------------------------------------- */
ACSF::ACSF() {}

ACSF::~ACSF() {

    /* free the allocated memory for two-body symmetric fucntion*/
    for (auto *each: two_body_symmetric_functions)
        delete each;
    two_body_symmetric_functions.clear();

    /* free the allocated memory for three-body symmetric fucntion*/
    for (auto *each: three_body_symmetric_functions)
        delete each;
    three_body_symmetric_functions.clear();
}

void ACSF::addTwoBodySymmetricFunction(TwoBodySymmetricFunction *symmetric_function) {
    two_body_symmetric_functions.push_back(symmetric_function);
}

void ACSF::addThreeBodySymmetricFunction(ThreeBodySymmetricFunction *symmetric_function) {
    three_body_symmetric_functions.push_back(symmetric_function);
}

std::vector<double> ACSF::calculate(AtomicConfiguration &configuration) 
{
    // TODO: optimization
    int n_2b = two_body_symmetric_functions.size();
    int n_3b = three_body_symmetric_functions.size();
    std::cout << n_2b << " " << n_3b << "" << std::endl;

    std::vector<double> results(n_2b+n_3b);
    std::fill(results.begin(), results.end(), 0.0);
    std::cout << results.size() << std::endl;

    for(Atom &atom_i: configuration.atoms) {
        for(Atom &atom_j: configuration.atoms) {

            /*Two-body symmetric functions*/
            if (atom_i.index == atom_j.index) continue;
            double rij = configuration.distance(atom_i, atom_j);

            for (int n=0; n<n_2b; n++) {
                results[n] += two_body_symmetric_functions[n]->function(rij);
            }

            /*Three-body symmetric functions*/
            for(Atom &atom_k: configuration.atoms) { 

                if (atom_i.index == atom_k.index) continue;
                if (atom_j.index == atom_k.index) continue;

                double rik = configuration.distance(atom_i, atom_k);
                double rjk = configuration.distance(atom_j, atom_k);

                for (int n=0; n<n_3b; n++) {
                    results[n+n_2b] += three_body_symmetric_functions[n]->function(rij, rik, rjk);
                }
            }
        }
    }
    return results;
}
