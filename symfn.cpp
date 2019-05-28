//
// Created by hossein on 5/27/19.
//

#include <cmath>
#include "symfn.h"

# define M_PI   3.14159265358979323846  /* pi */


/* ----------------------------------------------------------------------
   setup for base symmetric function
------------------------------------------------------------------------- */

double SymmetricFunction::fn_cutoff(double r)
{
    return ( cos(M_PI*r/r_cutoff) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetric function
------------------------------------------------------------------------- */

double G0::calculate()
{

}

/* ----------------------------------------------------------------------
   setup for G1 symmetric function
------------------------------------------------------------------------- */

double G1::function(double rij)
{
    return exp( -eta * (rij-rs) * (rij-rs) ) * fn_cutoff(rij);
}

double G1::calculate()
{

}

/* ----------------------------------------------------------------------
   setup for G4 symmetric function
------------------------------------------------------------------------- */

double G4::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * fn_cutoff(rij) * fn_cutoff(rjk) * fn_cutoff(rik);
}

double G4::calculate()
{

}


/* ----------------------------------------------------------------------
   setup for G5 symmetric function
------------------------------------------------------------------------- */

double G5::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * fn_cutoff(rij) * fn_cutoff(rjk);
}

double G5::calculate()
{

}