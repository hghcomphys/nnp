/*
  distance.cpp: This file is part of Free Molecular Dynamics

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
// Distance
//

#include "distance.h"

/* ----------------------------------------------------------------------
   setup for Distance
------------------------------------------------------------------------- */
Distance::Distance() : dr(0), drVec{0, 0, 0}, inv_dr(0){};

void Distance::set(double r, const double vec[3], double factor)
{
  dr = r;
  // inv_dr = 1.0 / r;
  set_drVec(vec, factor);
}

Distance::Distance(double r, const double vec[3])
{
  set(r, vec);
}