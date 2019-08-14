//
// Symmetric Funciton Scaler
//

#include <cmath>
#include "symmfuncscaler.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <iostream>


/* ----------------------------------------------------------------------
   setup for scaler
------------------------------------------------------------------------- */
Scaler::Scaler(double min, double max, double mean, double sigma):
            min(min), max(max), mean(mean), sigma(sigma),
            smin(0.0), smax(1.0) { setScalingFactor(); }

double Scaler::getMin() const { return min; }
double Scaler::getMax() const { return max; }
double Scaler::getMean() const { return mean; }
double Scaler::getSigma() const { return sigma; }
double Scaler::getScalingMinValue() const { return smin; }
double Scaler::getScalingMaxValue() const { return smax; }
void Scaler::setScalingMinValue(double value) { smin = value; }
void Scaler::setScalingMaxValue(double value) { smax = value; }
double Scaler::getScalingFactor() const { return scalingFactor; }

void Scaler::setScalingFactor() { 
    scalingFactor = (smax - smin) / (max - min); 
}

double Scaler::scale(double value) const {
    return smin +  scalingFactor * (value - mean) ; 
}


/* ----------------------------------------------------------------------
   setup for fymmetry function scaler
------------------------------------------------------------------------- */
SymmeryFunctionsScaler::SymmeryFunctionsScaler(): numberOfWarnings(0), maxNumberOfWarnings(100) {}

void SymmeryFunctionsScaler::addScaler(const Scaler& newScaler) { 
    listOfScalers.push_back(newScaler); 
}

int SymmeryFunctionsScaler::getNumberOfScalers() const { 
    return listOfScalers.size(); 
}

const std::vector<Scaler>& SymmeryFunctionsScaler::getListOfScalers() const { 
    return listOfScalers; 
}

void SymmeryFunctionsScaler::readScaling(const std::string& filename, int elementIndex) 
{
    std::string line;
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error("Unable to open file " + filename);

    while ( std::getline(inFile, line) ) {

        if( line[0] == '#' ) continue;  // skip comments

        std::stringstream ss(line);
        double min, max, mean, sigma=0.0;
        int elemIndex, index;

        ss >> elemIndex >> index >> min >> max >> mean >> sigma;
        // std::cout << << index << " " << min << " " << max << " " << mean << " " << sigma << "\n";

        if (elementIndex == elemIndex)
            addScaler( Scaler(min, max, mean, sigma) );
    } 
    inFile.close(); 
}


std::vector<double> SymmeryFunctionsScaler::scale(const std::vector<double>& values) 
{
    if (values.size() != getNumberOfScalers()) {
        std::cout << values.size() << " " << getNumberOfScalers() << "\n";
        throw std::runtime_error("Inconsistent number of symmetry functions and scalers");
    }

    std::vector<double> scaledValues( getNumberOfScalers() );
    for(int i=0; i<getNumberOfScalers(); i++) {

        const Scaler& sc = listOfScalers[i];
        const double value = values[i];

        if ( value > sc.getMax() || value < sc.getMin() ) {
            numberOfWarnings++;
            // std::cout << "Atom:" <<  atom_i.getIndex() << ":" << i+1 << ": " << (values[i] - sc.sfMin) / (sc.sfMax - sc.sfMin) << "\n";
            
            if (numberOfWarnings > maxNumberOfWarnings)
                throw std::runtime_error("Exceeds maximum number of symmetry function scaler warnings");
        }
        scaledValues[i] = sc.scale(value);
    }  
    return scaledValues;
}

std::vector<double> SymmeryFunctionsScaler::getScalingFactors() const {
    std::vector<double> scalingFactors;
    for (const Scaler& each: listOfScalers)
        scalingFactors.push_back( each.getScalingFactor() );
    return scalingFactors;
}

void SymmeryFunctionsScaler::setMaxNumberOfWarnings(double number) { 
    maxNumberOfWarnings=number; 
}

int SymmeryFunctionsScaler::getMaxNumberOfWarnings() const { 
    return maxNumberOfWarnings; 
}