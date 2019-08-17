//
// Cutoff Function
//

#include "cutoff_function.h"

/* ----------------------------------------------------------------------
   setup for base cutoff function
------------------------------------------------------------------------- */
CutoffFunction::CutoffFunction() {}

void CutoffFunction::setCutoffRadius(double cutoffRadius) {
    rc = cutoffRadius;
    inv_rc = 1.0 / cutoffRadius;
}
