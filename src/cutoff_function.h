/*
  cutoff_function.h: This file is part of Free Molecular Dynamics

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

#ifndef NNP_CUTOFFFUNCTION_H
#define NNP_CUTOFFFUNCTION_H

#include <cmath>

class CutoffFunction
{
public:
  CutoffFunction();
  void setCutoffRadius(double cutoffRadius);
  double fc(double r);
  double dfc(double r);

private:
  double rc;
  double inv_rc;
};

#endif // NNP_CUTOFFFUNCTION_H
