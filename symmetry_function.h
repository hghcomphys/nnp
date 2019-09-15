
//
// Symmetry Functions
//

#ifndef NNP_SYMMETRY_FUNCTION_H
#define NNP_SYMMETRY_FUNCTION_H

#include <vector>
#include "cutoff_function.h"

class SymmetryFunction {
public:
    SymmetryFunction(double cutoffRadius);
    double getCutoffRadius();
    CutoffFunction cutoffFunction;
protected:
    double cutoffRadius;
};

// Two-body symmetry function base class
class TwoBodySymmetryFunction: public SymmetryFunction { 
public:
    double gradientValue[3];
    TwoBodySymmetryFunction(double cutoffRadius): SymmetryFunction(cutoffRadius) {};
    virtual double function(double rij) {}; 
    virtual void gradient_ii(double rij, double drij[3]) {}; 
    virtual void gradient_ij(double rij, double drij[3]) {};
};

// Three-body symmetry function base class
class ThreeBodySymmetryFunction: public SymmetryFunction {
public:
    double gradientValue[3];
    ThreeBodySymmetryFunction(double cutoffRadius): SymmetryFunction(cutoffRadius) {};
    virtual double function(double rij, double rik, double rjk, double cost) {}; 
    virtual void gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) {}; 
    virtual void gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) {}; 
    virtual void gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) {}; 
};

// ------------------------- two-body symmetry function derived classes

class G0 : public TwoBodySymmetryFunction {
public:
    G0(std::vector<double> p);
    double function(double rij);
    void gradient_ii(double rij, double drij[3]);
    void gradient_ij(double rij, double drij[3]);
};


class G2 : public TwoBodySymmetryFunction {
public:
    G2(std::vector<double> p);
    double function(double rij);
    void gradient_ii(double rij, double drij[3]);
    void gradient_ij(double rij, double drij[3]);
private:
    double eta, rshift;
};

// ------------------------- three-body symmetry function derived classes

class G4 : public ThreeBodySymmetryFunction {
public:
    G4(std::vector<double> p);
    double function(double rij, double rik, double rjk, double cost);
    void gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]); 
    void gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]); 
    void gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]); 
private:
    double eta, zeta, lambda, rshift;
};


class G5 : public ThreeBodySymmetryFunction {
private:
    double eta, zeta, lambda, rshift;
public:
    G5(std::vector<double> p);
    double function(double rij, double rik, double rjk, double cost);
    void gradient_ii(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) {}; 
    void gradient_ij(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) {}; 
    void gradient_ik(double rij, double rik, double rjk, double cost, double drij[3], double drik[3], double drjk[3]) {}; 
};


#endif //NNP_SYMMETRY_FUNCTION_H