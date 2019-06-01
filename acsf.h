//
// Created by hossein on 5/27/19.
//

#ifndef NNP_SYMMETRIC_FUNCTIONS_H
#define NNP_SYMMETRIC_FUNCTIONS_H

#include <vector>

//namespace NNP_SF {

typedef enum  {
    TG0, TG1, TG4, TG5
}  SymmetricFunctionType;


class SymmetricFunction {
public:
    double cutoff_radius;
    double cutoff_function(double r);
    SymmetricFunction(double cutoff_radius);
    virtual double descriptor(double rij) {}; /*two-body*/
    virtual double descriptor(double rij, double rik, double jk) {}; /*three-body*/
    virtual double calculate() {};
};


class G0 : public SymmetricFunction {
public:
    G0(std::vector<double> p);
    double descriptor(double rij);
    double descriptor(double rij, double rik, double jk);
    double calculate();
};


class G1 : public SymmetricFunction {
private:
    double eta, rs;
public:
    G1(std::vector<double> p);
    double descriptor(double rij);
    double descriptor(double rij, double rik, double jk);
    double calculate();
};


class G4 : public SymmetricFunction {
private:
    double cost, eta, zeta, lamb;
public:
    G4(std::vector<double> p);
    double descriptor(double rij);
    double descriptor(double rij, double rik, double rjk);
    double calculate();
};


class G5 : public SymmetricFunction {
private:
    double cost, eta, zeta, lamb;
public:
    G5(std::vector<double> p);
    double descriptor(double rij);
    double descriptor(double rij, double rik, double rjk);
    double calculate();
};


class ACSF {
public:
    /* factory method */
    std::vector<SymmetricFunction *> symmetric_functions;
    void add(SymmetricFunction *symmetric_function);
    ACSF();
    ~ACSF();
};

//}

#endif //NNP_SYMMETRIC_FUNCTIONS_H
