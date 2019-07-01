//
// Created by hossein on 6/8/19.
//

#include <cmath>
#include "cutofffunc.h"

/* ----------------------------------------------------------------------
   setup for base cutoff function
------------------------------------------------------------------------- */
CutoffFunction::CutoffFunction() {}

void CutoffFunction::setCutoffRadius(double cutoffRadius) {
    rc = cutoffRadius;
    inv_rc = 1.0 / cutoffRadius;
}

// TODO: other types of cutoff function
double CutoffFunction::fc(double r) 
{
    if ( r > rc ) return 0;
    
    // COS TYPE
    // return ( cos(M_PI*r*inv_rc) + 1.0 ) * 0.5; 
    
    // TANH TYPE
    double const tanh1 = tanh(1.0 - r * inv_rc);
    return tanh1 * tanh1 * tanh1;
}