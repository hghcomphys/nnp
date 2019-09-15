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
    std::vector<double> scale(const std::vector<double>& values);
    void setMaxNumberOfWarnings(double number);
    int getMaxNumberOfWarnings();
    std::vector<double> getScalingFactors();

private:
    std::vector<Scaler> listOfScalers;
    int numberOfWarnings, maxNumberOfWarnings;
};


#endif //NNP_SYMMETRYFUNCTIONSCALER_H
