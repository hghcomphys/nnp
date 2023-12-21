/*
  atom.h: This file is part of Free Molecular Dynamics

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
// Atom
//

#ifndef NNP_ATOM_H
#define NNP_ATOM_H

#include <string>

class Atom
{
public:
  Atom(int index, const std::string &element, const double position[3]);
  Atom(int index, const std::string &element, const double position[3], const double force[3], double charge = 0.0, double energy = 0.0);
  void setPosition(const double position[3]);
  void setForce(const double force[3]);
  bool isElement(const std::string &element);
  std::string toString();

  // private:
  int index;
  double x, y, z;
  double fx, fy, fz;
  double charge, energy;
  std::string element;
};

int getAtomicNumber(const std::string &element);

#endif // NNP_ATOM_H