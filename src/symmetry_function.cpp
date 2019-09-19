/*
  base.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Hossein Ghorbanfekr

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//
// Symmetry Functions
//

#include <cmath>
#include <vector>
#include <stdexcept>
#include "symmetry_function.h"
#include "logger.h"

/* ----------------------------------------------------------------------
   setup for base symmetry function base
------------------------------------------------------------------------- */
SymmetryFunction::SymmetryFunction(double cutoffRadius): cutoffRadius(cutoffRadius) 
{
    cutoffFunction.setCutoffRadius(cutoffRadius);
}


/* ----------------------------------------------------------------------
   setup for G0 symmetry function
------------------------------------------------------------------------- */
G0::G0(std::vector<double> p): TwoBodySymmetryFunction(p[0]) 
{
    if ( p.size()!=1 ) 
        throw std::runtime_error( (Log(ERROR) << "Expected rcutoff argument").toString() );
}

double G0::function(double rij) 
{
    if ( rij > cutoffRadius ) return 0;
    return cutoffFunction.fc(rij);
}

void  G0::gradient_ii(double rij, double drij[3]) 
{
    if ( rij > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }
    // gradient of symmytry function of atom i respect to itself
    const double temp = cutoffFunction.dfc(rij) / rij;
    for (int d=0; d<3; d++)
        gradientValue[d] = drij[d] * temp;
}

void  G0::gradient_ij(double rij, double drij[3]) 
{
    if ( rij > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }
    // gradient of symmytry function of atom i respect to other atom j 
    gradient_ii(rij, drij);
    for (int d=0; d<3; d++)
        gradientValue[d] *= -1.0;
}

/* ----------------------------------------------------------------------
   setup for G2 symmetry function
------------------------------------------------------------------------- */
G2::G2(std::vector<double> p): eta(p[0]), rshift(p[1]), TwoBodySymmetryFunction(p[2]) 
{
    if ( p.size()<3 ) 
        throw std::runtime_error( (Log(ERROR) << "Expected eta, rshift, and rcutoff arguments").toString() );
}

double G2::function(double rij) 
{
    if ( rij > cutoffRadius ) return 0;
        return exp( -eta * (rij-rshift) * (rij-rshift) ) * cutoffFunction.fc(rij);
}

void G2::gradient_ii(double rij, double drij[3]) 
{
    if ( rij > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }
    // gradient of symmytry function of atom i respect to itself
    const double rp = rij - rshift;
    const double temp = ( cutoffFunction.dfc(rij) - 2.0 * eta * rp * cutoffFunction.fc(rij) ) * exp( -eta * rp * rp ) / rij ;
    for (int d=0; d<3; d++)
        gradientValue[d] = drij[d] * temp;
}

void G2::gradient_ij(double rij, double drij[3]) 
{
    if ( rij > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }
    // gradient of symmytry function of atom i respect to other atom j
    gradient_ii(rij, drij);
    for (int d=0; d<3; d++)
        gradientValue[d] *= -1.0;
}


/* ----------------------------------------------------------------------
   setup for G4 symmetry function
------------------------------------------------------------------------ */
G4::G4(std::vector<double> p): eta(p[0]), lambda(p[1]), zeta(p[2]), ThreeBodySymmetryFunction(p[3]), rshift(0) 
{
    if ( p.size()<4 ) 
        throw std::runtime_error( (Log(ERROR) << "Expected eta, lambda, zeta, and rcutoff arguments").toString() );
}

double G4::function(double rij, double rik, double rjk, double cost)
{
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) return 0;
    const double res =  pow(2.0, 1.0-zeta) * pow(1.0+lambda*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
}

void G4::gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) 
{
    // TODO: optimize performance
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }

    // gradient of atom i respect to j 
    const double inv_rij = 1.0 / rij;
    const double inv_rik = 1.0 / rik;
    const double inv_rjk = 1.0 / rjk;

    const double term1 = pow(1.0+lambda*cost, zeta);
    const double coef1 = lambda * zeta * pow(1.0+lambda*cost, zeta-1);
    double dterm1[3];
    for (int d=0; d<3; d++)
        dterm1[d] = coef1 * ( (drij[d] + drik[d]) * inv_rij * inv_rik - cost * (drij[d] * inv_rij * inv_rij + drik[d] * inv_rik * inv_rik) );

    const double term2 = exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    const double coef2 = -2.0 * eta * term2;
    double dterm2[3];
    for (int d=0; d<3; d++)
        dterm2[d] =  coef2 * ( drij[d] + drik[d] );  

    const double term3 = cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
    const double coef3 = cutoffFunction.fc(rjk);
    double dterm3[3];
    for (int d=0; d<3; d++)
        dterm3[d] = coef3 * ( cutoffFunction.dfc(rij) * cutoffFunction.fc(rik) * drij[d] * inv_rij + cutoffFunction.fc(rij) * cutoffFunction.dfc(rik) * drik[d] * inv_rik );

    const double coef = pow(2.0, 1.0-zeta);
    for (int d=0; d<3; d++) 
        gradientValue[d] = coef * ( dterm1[d] * term2 * term3 + term1 * dterm2[d] * term3 + term1 * term2 * dterm3[d] );
}


void G4::gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) 
{
    // TODO: optimize performance
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }

    // gradient of atom i respect to j
    const double inv_rij = 1.0 / rij;
    const double inv_rik = 1.0 / rik;
    const double inv_rjk = 1.0 / rjk;

    const double term1 = pow(1.0+lambda*cost, zeta);
    const double coef1 = lambda * zeta * pow(1.0+lambda*cost, zeta-1) * inv_rij;
    double dterm1[3];
    for (int d=0; d<3; d++)
        dterm1[d] = coef1 * ( -drik[d] * inv_rik  + cost * drij[d] * inv_rij );

    const double term2 = exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    const double coef2 = -2.0 * eta * term2;
    double dterm2[3];
    for (int d=0; d<3; d++)
        dterm2[d] =  coef2 * ( -drij[d] + drjk[d] );  

    const double term3 = cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
    const double coef3 = cutoffFunction.fc(rik);
    double dterm3[3];
    for (int d=0; d<3; d++)
        dterm3[d] = coef3 * ( -cutoffFunction.dfc(rij) * cutoffFunction.fc(rjk) * drij[d] * inv_rij + cutoffFunction.fc(rij) * cutoffFunction.dfc(rjk) * drjk[d] * inv_rjk );

    const double coef = pow(2.0, 1.0-zeta);
    for (int d=0; d<3; d++) 
        gradientValue[d] = coef * ( dterm1[d] * term2 * term3 + term1 * dterm2[d] * term3 + term1 * term2 * dterm3[d] );
}


void G4::gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) 
{
    // TODO: optimize performance
    if ( rij > cutoffRadius || rik > cutoffRadius || rjk > cutoffRadius ) 
    {
        for (int d=0; d<3; d++)
            gradientValue[d] = 0.0;
        return;
    }

    // gradient of atom i respect to j
    const double inv_rij = 1.0 / rij;
    const double inv_rik = 1.0 / rik;
    const double inv_rjk = 1.0 / rjk;

    const double term1 = pow(1.0+lambda*cost, zeta);
    const double coef1 = lambda * zeta * pow(1.0+lambda*cost, zeta-1) * inv_rik;
    double dterm1[3];
    for (int d=0; d<3; d++)
        dterm1[d] = coef1 * ( -drij[d] * inv_rij  + cost * drik[d] * inv_rik );

    const double term2 = exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    const double coef2 = -2.0 * eta * term2;
    double dterm2[3];
    for (int d=0; d<3; d++)
        dterm2[d] =  coef2 * ( -drik[d] - drjk[d] );  

    const double term3 = cutoffFunction.fc(rij) * cutoffFunction.fc(rik) * cutoffFunction.fc(rjk);
    const double coef3 = cutoffFunction.fc(rij);
    double dterm3[3];
    for (int d=0; d<3; d++)
        dterm3[d] = coef3 * ( -cutoffFunction.dfc(rik) * cutoffFunction.fc(rjk) * drik[d] * inv_rik - cutoffFunction.fc(rik) * cutoffFunction.dfc(rjk) * drjk[d] * inv_rjk );

    const double coef = pow(2.0, 1.0-zeta);
    for (int d=0; d<3; d++) 
        gradientValue[d] = coef * ( dterm1[d] * term2 * term3 + term1 * dterm2[d] * term3 + term1 * term2 * dterm3[d] );
}


/* ----------------------------------------------------------------------
   setup for G5 symmetry function
------------------------------------------------------------------------- */
G5::G5(std::vector<double> p): eta(p[0]), lambda(p[1]), zeta(p[2]), ThreeBodySymmetryFunction(p[3]), rshift(0) 
{
     if ( p.size()<4 ) 
        throw std::runtime_error( (Log(ERROR) << "Expected eta, lambda, zeta, and rcutoff arguments").toString() );
}

double G5::function(double rij, double rik, double rjk, double cost) 
{
    if ( rij > cutoffRadius || rik > cutoffRadius ) return 0;
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lambda*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoffFunction.fc(rij) * cutoffFunction.fc(rik);
}

