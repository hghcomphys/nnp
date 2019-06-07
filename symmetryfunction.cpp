//
// Symmetry Functions
//

#include <cmath>
#include <vector>
#include "symmetry_function.h"

/* ----------------------------------------------------------------------
   setup for base symmetry function base
------------------------------------------------------------------------- */
SymmetryFunction::SymmetryFunction(double cutoff_radius): cutoff_radius(cutoff_radius) {}

double SymmetryFunction::cutoff_function(double r) {
    // if ( r > cutoff_radius ) return 0;
    return ( cos(M_PI*r/cutoff_radius) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetry function
------------------------------------------------------------------------- */
G0::G0(std::vector<double> p): TwoBodySymmetryFunction(p[0]) {}

double G0::function(double rij) {
    if ( rij > cutoff_radius ) return 0;
    return cutoff_function(rij);
}

/* ----------------------------------------------------------------------
   setup for G1 symmetry function
------------------------------------------------------------------------- */
G1::G1(std::vector<double> p): TwoBodySymmetryFunction(p[0]), eta(p[1]), rs(p[2]) {}

double G1::function(double rij) {
    if ( rij > cutoff_radius ) return 0;
    return exp( -eta * (rij-rs) * (rij-rs) ) * cutoff_function(rij);
}

/* ----------------------------------------------------------------------
   setup for G4 symmetry function
------------------------------------------------------------------------ */
G4::G4(std::vector<double> p): ThreeBodySymmetryFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G4::function(double rij, double rik, double rjk)
{
    if ( rij > cutoff_radius || rik > cutoff_radius || rjk > cutoff_radius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoff_function(rij) * cutoff_function(rjk) * cutoff_function(rik);
}

/* ----------------------------------------------------------------------
   setup for G5 symmetry function
------------------------------------------------------------------------- */
G5::G5(std::vector<double> p): ThreeBodySymmetryFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G5::function(double rij, double rik, double rjk)
{
    if ( rij > cutoff_radius || rik > cutoff_radius || rjk > cutoff_radius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoff_function(rij) * cutoff_function(rjk);
}