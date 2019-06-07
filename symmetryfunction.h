
//
// Symmetry Functions
//

#ifndef NNP_SYMMETRY_FUNCTION_H
#define NNP_SYMMETRY_FUNCTION_H

#include <vector>

//namespace NNP_SF {

//typedef enum  {
//    TG0, TG1, TG4, TG5
//}  SymmetryFunctionType;

class SymmetryFunction {
public:
    double cutoff_radius;
    double cutoff_function(double r);
    SymmetryFunction(double cutoff_radius);
};

// Two-body
class TwoBodySymmetryFunction: public SymmetryFunction { 
public:
    TwoBodySymmetryFunction(double cutoff_radius): SymmetryFunction(cutoff_radius) {};
    virtual double function(double rij) {}; 
};

// Three-body
class ThreeBodySymmetryFunction: public SymmetryFunction {
public:
    ThreeBodySymmetryFunction(double cutoff_radius): SymmetryFunction(cutoff_radius) {};
    virtual double function(double rij, double rik, double jk) {}; 
};


class G0 : public TwoBodySymmetryFunction {
public:
    G0(std::vector<double> p);
    double function(double rij);
};


class G1 : public TwoBodySymmetryFunction {
private:
    double eta, rs;
public:
    G1(std::vector<double> p);
    double function(double rij);
};


class G4 : public ThreeBodySymmetryFunction {
private:
    double cost, eta, zeta, lamb;
public:
    G4(std::vector<double> p);
    double function(double rij, double rik, double rjk);
};


class G5 : public ThreeBodySymmetryFunction {
private:
    double cost, eta, zeta, lamb;
public:
    G5(std::vector<double> p);
    double function(double rij, double rik, double rjk);
};

//}

#endif //NNP_SYMMETRY_FUNCTION_H