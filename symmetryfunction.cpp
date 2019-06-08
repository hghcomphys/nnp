//
// Symmetry Functions
//

#include <cmath>
#include <vector>
#include "symmetryfunction.h"

/* ----------------------------------------------------------------------
   setup for base symmetry function base
------------------------------------------------------------------------- */
SymmetryFunction::SymmetryFunction(double cutoffRadius): cutoffRadius(cutoffRadius) {}

double SymmetryFunction::getCutoffRadius() { return cutoffRadius; }

double SymmetryFunction::cutoffFunction(double r) {
    // if ( r > cutoffRadius ) return 0;
    return ( cos(M_PI*r/cutoffRadius) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetry function
------------------------------------------------------------------------- */
G0::G0(std::vector<double> p): TwoBodySymmetryFunction(p[0]) {}

double G0::function(double rij) {
    if ( rij > cutoffRadius ) return 0;
    return cutoffFunction(rij);
}

/* ----------------------------------------------------------------------
   setup for G1 symmetry function
------------------------------------------------------------------------- */
G1::G1(std::vector<double> p): TwoBodySymmetryFunction(p[0]), eta(p[1]), rs(p[2]) {}

double G1::function(double rij) {
    if ( rij > cutoffRadius ) return 0;
    return exp( -eta * (rij-rs) * (rij-rs) ) * cutoffFunction(rij);
}

/* ----------------------------------------------------------------------
   setup for G4 symmetry function
------------------------------------------------------------------------ */
G4::G4(std::vector<double> p): ThreeBodySymmetryFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G4::function(double rij, double rik, double rjk)
{
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoffFunction(rij) * cutoffFunction(rjk) * cutoffFunction(rik);
}

/* ----------------------------------------------------------------------
   setup for G5 symmetry function
------------------------------------------------------------------------- */
G5::G5(std::vector<double> p): ThreeBodySymmetryFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G5::function(double rij, double rik, double rjk)
{
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoffFunction(rij) * cutoffFunction(rjk);
}