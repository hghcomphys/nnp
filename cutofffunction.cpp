//
// Created by hossein on 6/8/19.
//

#include <cmath>
#include "cutofffunction.h"

/* ----------------------------------------------------------------------
   setup for base cutoff function
------------------------------------------------------------------------- */
CutoffFunction::CutoffFunction() {}

void CutoffFunction::setCutoffRadius(double cutoffRadius) {
    rc = cutoffRadius;
    inv_rc = 1.0/cutoffRadius;
}

// TODO: other types of cutoff function
double CutoffFunction::fc(double r) {
    if ( r > rc ) return 0;
    return ( cos(M_PI*r*inv_rc) + 1.0 ) * 0.5;
}