//
// Atomic Centered Symmetry Functions
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "atoms.h"

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

/*Two-body*/
class TwoBodySymmetryFunction: public SymmetryFunction { 
public:
    TwoBodySymmetryFunction(double cutoff_radius): SymmetryFunction(cutoff_radius) {};
    virtual double function(double rij) {}; 
};

 /*Three-body*/
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


class ACSF {
public:
    std::vector<TwoBodySymmetryFunction *> two_body_symmetry_functions; /* factory method */
    std::vector<ThreeBodySymmetryFunction *> three_body_symmetry_functions; /* factory method */
    ACSF();
    ~ACSF();
    void addTwoBodySymmetryFunction(TwoBodySymmetryFunction *symmetry_function); /*add two-body symmetry function*/
    void addThreeBodySymmetryFunction(ThreeBodySymmetryFunction *symmetry_function); /*add three-body symmetry function*/
    std::vector<double> calculate(Atoms &configuration);
};

//}

#endif //NNP_ACSF_H
