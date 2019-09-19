/*
  structure.h: This file is part of Free Molecular Dynamics

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
// Atomic Structure
//

#ifndef NNP_STRUCTURE_H
#define NNP_STRUCTURE_H

#include "atom.h"
#include "distance.h"
#include <map>
#include <string>

class AtomicStructure {
public:
    AtomicStructure();
    ~AtomicStructure();
    void writeFileFormatRunner(const std::string& filename);
    void readFileFormatRuNNer(const std::string& filename);
    void readFileFormatRuNNer();
    Atom& getAtom(int atomIndex);
    void setCell(const double cell[9]);
    double distance(Atom *atom_i, Atom *atom_j, double drij[3]);
    double distance(Atom *atom_i, Atom *atom_j);
    void calculateTableOfDistances(double globalCutOffRadius = 12.0);

// private:
    bool isAtom, isCell;
    bool isTableOfDistances;
    int numberOfAtoms;
    double cell[9];
    double totalCharge, totalEnergy;
    // TODO: merge together listOfAtoms and ListOfAtomsForElement
    Atom **atoms;
    std::map<std::string, int> numberOfAtomsForElement;
    std::map<std::string, Atom**> atomsForElement; 
    Distance **tableOfDistances;
    void applyPBC(double& dx, double& dy, double& dz);
};

inline Atom& AtomicStructure::getAtom(int atomIndex) 
{ 
    return *atoms[atomIndex]; 
}

#endif //NNP_STRUCTURE_H