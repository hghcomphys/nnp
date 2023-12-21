/*
  symmetry_function.h: This file is part of Free Molecular Dynamics

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
// Symmetry Functions
//

#ifndef NNP_SYMMETRY_FUNCTION_H
#define NNP_SYMMETRY_FUNCTION_H

#include <vector>
#include "cutoff_function.h"

class SymmetryFunction
{
public:
    SymmetryFunction(double cutoffRadius);
    // protected:
    CutoffFunction cutoffFunction;
    double cutoffRadius;
    double gradientValue[3];
};

// Two-body symmetry function base class
class TwoBodySymmetryFunction : public SymmetryFunction
{
public:
    TwoBodySymmetryFunction(double cutoffRadius) : SymmetryFunction(cutoffRadius){};
    virtual double function(double rij){};
    virtual void gradient_ii(double rij, double drij[3]){};
    virtual void gradient_ij(double rij, double drij[3]){};
};

// Three-body symmetry function base class
class ThreeBodySymmetryFunction : public SymmetryFunction
{
public:
    ThreeBodySymmetryFunction(double cutoffRadius) : SymmetryFunction(cutoffRadius){};
    virtual double function(double rij, double rik, double rjk, double cost){};
    virtual void gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]){};
    virtual void gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]){};
    virtual void gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]){};
};

// ------------------------- two-body symmetry function derived classes

class G0 : public TwoBodySymmetryFunction
{
public:
    G0(std::vector<double> p);
    double function(double rij);
    void gradient_ii(double rij, double drij[3]);
    void gradient_ij(double rij, double drij[3]);
};

class G2 : public TwoBodySymmetryFunction
{
public:
    G2(std::vector<double> p);
    double function(double rij);
    void gradient_ii(double rij, double drij[3]);
    void gradient_ij(double rij, double drij[3]);

private:
    double eta, rshift;
};

// ------------------------- three-body symmetry function derived classes

class G4 : public ThreeBodySymmetryFunction
{
public:
    G4(std::vector<double> p);
    double function(double rij, double rik, double rjk, double cost);
    void gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]);
    void gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]);
    void gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]);

private:
    double eta, zeta, lambda, rshift;
};

class G5 : public ThreeBodySymmetryFunction
{
private:
    double eta, zeta, lambda, rshift;

public:
    G5(std::vector<double> p);
    double function(double rij, double rik, double rjk, double cost);
    void gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]){};
    void gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]){};
    void gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]){};
};

#endif // NNP_SYMMETRY_FUNCTION_H