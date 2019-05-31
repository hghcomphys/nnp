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

protected:
    virtual double function(double r);
    virtual double calculate();

public:
    double r_cutoff;
    double fn_cutoff(double r);
    SymmetricFunction(double r_cutoff);
};


class G0 : SymmetricFunction {
public:
    G0(double r_cutoff);
    double function(double r);
    double calculate();
};


class G1 : SymmetricFunction {
private:
    double eta, rs;

public:
    G1(double r_cutoff, double eta, double rs);
    double function(double r);
    double calculate();
};


class G4 : SymmetricFunction {
private:
    double cost, eta, zeta, lamb;

public:
    G4(double r_cutoff, double cost, double eta, double zeta, double lamb);
    double function(double rij, double rik, double rjk);
    double calculate();
};


class G5 : SymmetricFunction {
private:
    double cost, eta, zeta, lamb;

public:
    G5(double r_cutoff, double cost, double eta, double zeta, double lamb);
    double function(double rij, double rik, double rjk);
    double calculate();
};


class ACSF {
public:
    /* factory method */
    std::vector<SymmetricFunction *> list_of_symmetric_functions;
    void add_symmetric_function(SymmetricFunctionType select);
};

//}

#endif //NNP_SYMMETRIC_FUNCTIONS_H
