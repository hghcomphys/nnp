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

#include "symmetry_function_scaler.h"
#include "logger.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

/* ----------------------------------------------------------------------
   setup for scaler 
------------------------------------------------------------------------- */
Scaler::Scaler(double min, double max, double mean, double sigma):
    min(min), max(max), mean(mean), sigma(sigma), smin(0.0), smax(1.0) 
{ 
    setScalingFactor(); 
}

void Scaler::setScalingFactor() 
{ 
    scalingFactor = (smax - smin) / (max - min); 
}

double Scaler::scale(double value) 
{
    return smin +  scalingFactor * (value - mean) ; 
}

/* ----------------------------------------------------------------------
   setup for symmetry function scaler
------------------------------------------------------------------------- */
SymmeryFunctionsScaler::SymmeryFunctionsScaler(): numberOfWarnings(0), maxNumberOfWarnings(100) {}

void SymmeryFunctionsScaler::addScaler(const Scaler& newScaler) 
{ 
    listOfScalers.push_back(newScaler); 
}

int SymmeryFunctionsScaler::getNumberOfScalers() 
{ 
    return listOfScalers.size(); 
}

void SymmeryFunctionsScaler::readScaling(const std::string& filename, int elementIndex) 
{
    std::string line;
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " + filename).toString() );

    while ( std::getline(inFile, line) ) 
    {
        if( line[0] == '#' ) continue;  // skip comments

        std::stringstream ss(line);
        double min, max, mean, sigma=0.0;
        int elemIndex, index;

        ss >> elemIndex >> index >> min >> max >> mean >> sigma;
        Log(DEBUG) << "SF-scaler (" << index << "): " << min << " " << max << " " << mean << " " << sigma;

        if (elementIndex == elemIndex)
            addScaler( Scaler(min, max, mean, sigma) );
    } 
    inFile.close(); 
}

void SymmeryFunctionsScaler::scale(double *descriptorValues, int descriptorSize) 
{
    const int numberOfScalers = getNumberOfScalers();
    if (descriptorSize != numberOfScalers) 
        throw std::runtime_error(
            (Log(ERROR) << "Inconsistent number of symmetry functions (" << descriptorSize
            << ") and scalers (" << numberOfScalers << ")").toString()
            );

    for(int i=0; i<numberOfScalers; i++) {

        Scaler& sc = listOfScalers[i];
        const double value = descriptorValues[i];

        if ( value > sc.max || value < sc.min ) 
        {
            // increase number of warrning
            numberOfWarnings++;
            Log(WARN) << "Exceed symmetry function scaler (index=" << i+1 << ")";
              
            // raise error if it exceeds maximum number of warnings
            if (numberOfWarnings > maxNumberOfWarnings)
                throw std::runtime_error(
                    (Log(ERROR) << "Exceed maximum number of symmetry function scaler warnings ("
                    << maxNumberOfWarnings << ")").toString()
                    );
        }
        descriptorValues[i] = sc.scale(value);
    }  
}

std::vector<double> SymmeryFunctionsScaler::getScalingFactors() 
{
    std::vector<double> scalingFactors;
    for (const Scaler& each: listOfScalers)
        scalingFactors.push_back( each.scalingFactor );
    return scalingFactors;
}

void SymmeryFunctionsScaler::setMaxNumberOfWarnings(double number) 
{ 
    maxNumberOfWarnings=number; 
}

int SymmeryFunctionsScaler::getMaxNumberOfWarnings() 
{ 
    return maxNumberOfWarnings; 
}