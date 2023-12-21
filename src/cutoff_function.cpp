/*
  cutoff_function.cpp: This file is part of Free Molecular Dynamics

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
// Cutoff Function
//

#include "cutoff_function.h"

/* ----------------------------------------------------------------------
   setup for base cutoff function
------------------------------------------------------------------------- */
CutoffFunction::CutoffFunction() {}

void CutoffFunction::setCutoffRadius(double cutoffRadius)
{
    rc = cutoffRadius;
    inv_rc = 1.0 / cutoffRadius;
}

// TODO: other types of cutoff function
// TODO: using inline can reduce computational time
double CutoffFunction::fc(double r)
{
    if (r > rc)
        return 0;

    // COS TYPE
    // return ( cos(M_PI * r * inv_rc) + 1.0 ) * 0.5;

    // TANH TYPE
    double const tmp = tanh(1.0 - r * inv_rc);
    return tmp * tmp * tmp;
}

// TODO: other types of cutoff function
// TODO: using inline can reduce computational time
double CutoffFunction::dfc(double r)
{
    if (r > rc)
        return 0;

    // COS TYPE
    // return -M_PI_2 * inv_rc * sin(M_PI * r * inv_rc);

    // TANH TYPE
    double tmp = tanh(1.0 - r * inv_rc);
    tmp *= tmp;
    return 3.0 * tmp * (tmp - 1.0) * inv_rc;
}
