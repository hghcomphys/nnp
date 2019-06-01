//
// Atomic Centered Symmetric Functions
//

#include "acsf.h"
#include <cmath>

//using namespace NNP_SF;

/* ----------------------------------------------------------------------
   setup for base symmetric function base
------------------------------------------------------------------------- */

ACSF::ACSF(double cutoff_radius): cutoff_radius(cutoff_radius) {}

double ACSF::cutoff_function(double r) {
    return ( cos(M_PI*r/cutoff_radius) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetric function
------------------------------------------------------------------------- */

G0::G0(std::vector<double> p): ACSF(p[0]) {}

double G0::function(double rij) {
    return cutoff_function(rij);
}

double G0::function(double rij, double rik, double jk) { return 0; }

double G0::calculate() {}

/* ----------------------------------------------------------------------
   setup for G1 symmetric function
------------------------------------------------------------------------- */

G1::G1(std::vector<double> p): ACSF(p[0]), eta(p[1]), rs(p[2]) {}

double G1::function(double rij) {
    return exp( -eta * (rij-rs) * (rij-rs) ) * cutoff_function(rij);
}

double G1::function(double rij, double rik, double jk) { return 0; }

double G1::calculate() {}

/* ----------------------------------------------------------------------
   setup for G4 symmetric function
------------------------------------------------------------------------ */

G4::G4(std::vector<double> p): ACSF(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G4::function(double rij) { return 0; }

double G4::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoff_function(rij) * cutoff_function(rjk) * cutoff_function(rik);
}

double G4::calculate() {}


/* ----------------------------------------------------------------------
   setup for G5 symmetric function
------------------------------------------------------------------------- */

G5::G5(std::vector<double> p): ACSF(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G5::function(double rij) { return 0; }

double G5::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoff_function(rij) * cutoff_function(rjk);
}

double G5::calculate() {}
