//
// Cutoff Function
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
    // return ( cos(M_PI * r * inv_rc) + 1.0 ) * 0.5; 
    
    // TANH TYPE
    double const tmp = tanh(1.0 - r * inv_rc);
    return tmp * tmp * tmp;
}

// TODO: other types of cutoff function
double CutoffFunction::dfc(double r) 
{
    if ( r > rc ) return 0;
    
    // COS TYPE
    // return -M_PI_2 * inv_rc * sin(M_PI * r * inv_rc);
    
    // TANH TYPE
    double temp = tanh(1.0 - r * inv_rc);
    temp *= temp;
    return 3.0 * temp * (temp - 1.0) * inv_rc;
}