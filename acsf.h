//
// Atomic Centered Symmetric Functions
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "atoms.h"

//namespace NNP_SF {

//typedef enum  {
//    TG0, TG1, TG4, TG5
//}  SymmetricFunctionType;


class SymmetricFunction {
public:
    double cutoff_radius;
    double cutoff_function(double r);
    SymmetricFunction(double cutoff_radius);
};

/*Two-body*/
class TwoBodySymmetricFunction: public SymmetricFunction { 
public:
    TwoBodySymmetricFunction(double cutoff_radius): SymmetricFunction(cutoff_radius) {};
    virtual double function(double rij) {}; 
};

 /*Three-body*/
class ThreeBodySymmetricFunction: public SymmetricFunction {
public:
    ThreeBodySymmetricFunction(double cutoff_radius): SymmetricFunction(cutoff_radius) {};
    virtual double function(double rij, double rik, double jk) {}; 
};


class G0 : public TwoBodySymmetricFunction {
public:
    G0(std::vector<double> p);
    double function(double rij);
};


class G1 : public TwoBodySymmetricFunction {
private:
    double eta, rs;
public:
    G1(std::vector<double> p);
    double function(double rij);
};


class G4 : public ThreeBodySymmetricFunction {
private:
    double cost, eta, zeta, lamb;
public:
    G4(std::vector<double> p);
    double function(double rij, double rik, double rjk);
};


class G5 : public ThreeBodySymmetricFunction {
private:
    double cost, eta, zeta, lamb;
public:
    G5(std::vector<double> p);
    double function(double rij, double rik, double rjk);
};


class ACSF {
public:
    std::vector<TwoBodySymmetricFunction *> two_body_symmetric_functions; /* factory method */
    std::vector<ThreeBodySymmetricFunction *> three_body_symmetric_functions; /* factory method */
    ACSF();
    ~ACSF();
    void addTwoBodySymmetricFunction(TwoBodySymmetricFunction *symmetric_function); /*add two-body symmettic function*/
    void addThreeBodySymmetricFunction(ThreeBodySymmetricFunction *symmetric_function); /*add three-body symmettic function*/
    std::vector<double> calculate(AtomicConfiguration &configuration);
};

//}

#endif //NNP_ACSF_H
