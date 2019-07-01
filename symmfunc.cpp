//
// Symmetry Functions
//

#include <cmath>
#include <vector>
#include <stdexcept>
#include "symmfunc.h"
#include <iostream>

/* ----------------------------------------------------------------------
   setup for base symmetry function base
------------------------------------------------------------------------- */
SymmetryFunction::SymmetryFunction(double cutoffRadius): cutoffRadius(cutoffRadius) {
    cutoffFunction.setCutoffRadius(cutoffRadius);
}

double SymmetryFunction::getCutoffRadius() { return cutoffRadius; }


/* ----------------------------------------------------------------------
   setup for G0 symmetry function
------------------------------------------------------------------------- */
G0::G0(std::vector<double> p): TwoBodySymmetryFunction(p[0]) {
    if ( p.size()!=1 ) 
        throw std::runtime_error("Expected rcutoff argument");
}

double G0::function(double rij) {
    if ( rij > cutoffRadius ) return 0;
    return cutoffFunction.fc(rij);
}

/* ----------------------------------------------------------------------
   setup for G2 symmetry function
------------------------------------------------------------------------- */
G2::G2(std::vector<double> p): eta(p[0]), rshift(p[1]), TwoBodySymmetryFunction(p[2]) {
    if ( p.size()<3 ) 
        throw std::runtime_error("Expected eta, rshift, and rcutoff arguments");
}

double G2::function(double rij) {
    if ( rij > cutoffRadius ) return 0;
    return exp( -eta * (rij-rshift) * (rij-rshift) ) * cutoffFunction.fc(rij);
}

/* ----------------------------------------------------------------------
   setup for G4 symmetry function
------------------------------------------------------------------------ */
G4::G4(std::vector<double> p): eta(p[0]), lambda(p[1]), zeta(p[2]), ThreeBodySymmetryFunction(p[3]), rshift(0) {
    if ( p.size()<4 ) 
        throw std::runtime_error("Expected eta, lambda, zeta, and rcutoff arguments");
}

double G4::function(double rij, double rik, double rjk, double cost)
{
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lambda*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
}

/* ----------------------------------------------------------------------
   setup for G5 symmetry function
------------------------------------------------------------------------- */
G5::G5(std::vector<double> p): eta(p[0]), lambda(p[1]), zeta(p[2]), ThreeBodySymmetryFunction(p[3]), rshift(0) {
     if ( p.size()<4 ) 
        throw std::runtime_error("Expected eta, lambda, zeta, and rcutoff arguments");
}

double G5::function(double rij, double rik, double rjk, double cost)
{
    if ( rij > cutoffRadius || rik > cutoffRadius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lambda*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoffFunction.fc(rij) * cutoffFunction.fc(rik);
}