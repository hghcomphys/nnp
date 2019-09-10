//
// Distance
//

#include "distance.h"

/* ----------------------------------------------------------------------
   setup for Distance
------------------------------------------------------------------------- */
Distance::Distance(): dr(0), drVec{0, 0 ,0}, inv_dr(0) {};

void Distance::set(double r, double vec[3], double factor) 
{
    dr = r;
    // inv_dr = 1.0 / r;
    set_drVec(vec, factor);
}

Distance::Distance(double r, double vec[3]) {
    set(r, vec); 
}