
//
// Symmetry Functions
//

#ifndef NNP_SYMMETRY_FUNCTION_H
#define NNP_SYMMETRY_FUNCTION_H

#include <vector>
#include "cutofffunc.h"

class SymmetryFunction {
public:
    SymmetryFunction(double cutoffRadius);
    double getCutoffRadius();
    CutoffFunction cutoffFunction;
protected:
    double cutoffRadius;
};

// Two-body
class TwoBodySymmetryFunction: public SymmetryFunction { 
public:
    TwoBodySymmetryFunction(double cutoffRadius): SymmetryFunction(cutoffRadius) {};
    virtual double function(double rij) {}; 
};

// Three-body
class ThreeBodySymmetryFunction: public SymmetryFunction {
public:
    ThreeBodySymmetryFunction(double cutoffRadius): SymmetryFunction(cutoffRadius) {};
    virtual double function(double rij, double rik, double jk) {}; 
};


class G0 : public TwoBodySymmetryFunction {
public:
    G0(std::vector<double> p);
    double function(double rij);
};


class G1 : public TwoBodySymmetryFunction {
public:
    G1(std::vector<double> p);
    double function(double rij);
private:
    double eta, rshift;
};


class G4 : public ThreeBodySymmetryFunction {
public:
    G4(std::vector<double> p);
    double function(double rij, double rik, double rjk);
private:
    double cost, eta, zeta, lambda;
};


class G5 : public ThreeBodySymmetryFunction {
private:
    double cost, eta, zeta, lambda;
public:
    G5(std::vector<double> p);
    double function(double rij, double rik, double rjk);
};


#endif //NNP_SYMMETRY_FUNCTION_H