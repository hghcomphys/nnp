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
    double getMin() const;
    double getMax() const;
    double getMean() const;
    double getSigma() const;
    double getScalingMinValue() const;
    double getScalingMaxValue() const;
    void setScalingMinValue(double value);
    void setScalingMaxValue(double value);
    double scale(double value) const;

private:
    double min, max, mean, sigma;
    double smin, smax;
};


class SymmeryFunctionsScaler {
public:
    SymmeryFunctionsScaler();
    void addScaler(const Scaler& newScaler);
    int getNumberOfScalers() const;
    const std::vector<Scaler>& getListOfScalers() const;
    void readScaling(const std::string& filename, int elementIndex);
    std::vector<double> scale(const std::vector<double>& values);
    void setMaxNumberOfWarnings(double number);
    int getMaxNumberOfWarnings() const;


private:
    std::vector<Scaler> listOfScalers;
    int numberOfWarnings, maxNumberOfWarnings;
};


#endif //NNP_SYMMETRYFUNCTIONSCALER_H
