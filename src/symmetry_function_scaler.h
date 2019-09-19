/*
  base.c: This file is part of Free Molecular Dynamics

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
// Symmetric Funciton Scaler
//

#ifndef NNP_SYMMETRYFUNCTIONSCALER_H
#define NNP_SYMMETRYFUNCTIONSCALER_H

#include <vector>
#include <string>


class Scaler {
public:
    Scaler(double min, double max, double mean, double sigma);
    double scale(double value);
    void setScalingFactor();

// private:
    double min, max, mean, sigma;
    double smin, smax;
    double scalingFactor;
};


class SymmeryFunctionsScaler {
public:
    SymmeryFunctionsScaler();
    void addScaler(const Scaler& newScaler);
    int getNumberOfScalers();
    void readScaling(const std::string& filename, int elementIndex);
    void scale(double *descriptorValues, int descriptorSize);
    void setMaxNumberOfWarnings(double number);
    int getMaxNumberOfWarnings();
    std::vector<double> getScalingFactors();

private:
    std::vector<Scaler> listOfScalers;
    int numberOfWarnings, maxNumberOfWarnings;
};


#endif //NNP_SYMMETRYFUNCTIONSCALER_H
