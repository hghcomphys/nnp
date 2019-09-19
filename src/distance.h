/*
  distance.h: This file is part of Free Molecular Dynamics

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

#ifndef NNP_DISTANCE_H
#define NNP_DISTANCE_H

class Distance {
public:
    Distance();
    Distance(double r, const double vec[3]); 
    void set_drVec(const double vec[3], double factor=1.0);
    void set(double r, const double vec[3], double factor=1.0);
//private:
    double dr;
    double inv_dr;
    double drVec[3];
};

inline void Distance::set_drVec(const double vec[3], double factor)
{
    for (int d=0; d<3; d++)
        drVec[d] = vec[d] * factor;
}

#endif //NNP_DISTANCE_H
