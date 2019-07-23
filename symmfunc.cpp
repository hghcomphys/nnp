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

std::vector<double>  G0::gradient_ij(double rij, double drij[3]) {
    if ( rij > cutoffRadius ) 
        return std::vector<double>({0.0, 0.0, 0.0});

    std::vector<double> result(3);
    const double temp = -cutoffFunction.dfc(rij) / rij;
    for (int d=0; d<3; d++)
        result[d] = drij[d] * temp;
    return result;
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

std::vector<double>  G2::gradient_ij(double rij, double drij[3]) {
    if ( rij > cutoffRadius ) 
        return std::vector<double>({0.0, 0.0, 0.0});

    std::vector<double> result(3);
    const double rp = rij - rshift;
    const double exptemp = exp( -eta * rp * rp );
    const double temp = ( 2.0 * eta * rp * cutoffFunction.fc(rij) - cutoffFunction.dfc(rij) ) * exptemp / rij ;
    for (int d=0; d<3; d++)
        result[d] = drij[d] * temp;
    return result;
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
    const double res =  pow(2.0, 1.0-zeta) * pow(1.0+lambda*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
}

std::vector<double> G4::gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) 
{
    // TODO: optimize performance
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) 
        return std::vector<double>({0.0, 0.0, 0.0});

    const double coef = pow(2.0, 1.0-zeta);
    const double inv_rij = 1.0 / rij;
    const double inv_rik = 1.0 / rik;
    const double inv_rjk = 1.0 / rjk;

    const double term1 = pow(1.0+lambda*cost, zeta);
    const double coef1 = -lambda * zeta * pow(1.0+lambda*cost, zeta-1) * inv_rij;
    double dterm1[3];
    for (int d=0; d<3; d++)
        dterm1[d] = coef1 * ( drik[d] * inv_rik + cost * drij[d] * inv_rij );

    const double term2 = exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    const double coef2 = 2.0 * eta * term2;
    double dterm2[3];
    for (int d=0; d<3; d++)
        dterm2[d] =  coef2 * ( drij[d] - drjk[d] );  

    const double term3 = cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
    const double coef3 = cutoffFunction.fc(rik);
    double dterm3[3];
    for (int d=0; d<3; d++)
        dterm3[d] = coef3 * ( cutoffFunction.fc(rij) * cutoffFunction.dfc(rjk) * drjk[d] * inv_rjk -
            cutoffFunction.dfc(rij) * cutoffFunction.fc(rjk) * drij[d] * inv_rij );

    std::vector<double> result(3);
    for (int d=0; d<3; d++) 
        result[d] = coef * ( dterm1[d] * term2 * term3 + term1 * dterm2[d] * term3 + term1 * term2 * dterm3[d]);
    return result;
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

std::vector<double> G5::gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) 
{
    // TODO: Has to be implemented
}
