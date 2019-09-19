/*
  descriptor.h: This file is part of Free Molecular Dynamics

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
// Atomic Centered Symmetry Functions (ACSF)
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "structure.h"
#include "symmetry_function.h"
#include "symmetry_function_scaler.h"

class ACSF {
public:
    ACSF(std::string element);
    ~ACSF();
    void addTwoBodySF(TwoBodySymmetryFunction *symmetryFunction, const std::string& neighborElement1); /*add two-body symmetry function*/
    void addThreeBodySF(ThreeBodySymmetryFunction *symmetryFunction, const std::string& neighborElement1, const std::string& neibghorElement2); /*add three-body symmetry function*/
    TwoBodySymmetryFunction& getTwoBodySF(int index);
    ThreeBodySymmetryFunction& getThreeBodySF(int index);
    int getNumberOfTwoBodySF();
    int getNumberOfThreeBodySF();
    int getTotalNumberOfSF();
    double getGlobalCutOffRadius();
    void calculate(AtomicStructure& structure, Atom *atom, double *descriptorValues);
    void gradient(AtomicStructure &structure, Atom *atom_i, Atom *atom_ip, double **gradientValues, int gradientSize);

// private:
    std::string centralElement;
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<std::string> listOfTwoBodyNeighborElement;
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    std::vector<std::string> listOfThreeBodyNeighborElement1, listOfThreeBodyNeighborElement2;
};

inline double calculateCosine(double rij, double rik, const double drij[3], const double drik[3]) 
{
    // cosine of angle between k--<i>--j atoms
    double cost = 0;
    for (int d=0; d<3; d++)
        cost += drij[d] * drik[d];
    return cost / rij / rik;
}

#endif //NNP_ACSF_H
