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
    virtual double function();
    virtual double calculate();

public:
    double cutoff_radius;
    double cutoff_function(double r);
    SymmetricFunction(double cutoff_radius);
};


class G0 : SymmetricFunction {
public:
    G0(std::vector<double> p);
    double function(double r);
    double calculate();
};


class G1 : SymmetricFunction {
private:
    double eta, rs;

public:
    G1(std::vector<double> p);
    double function(double r);
    double calculate();
};


class G4 : SymmetricFunction {
private:
    double cost, eta, zeta, lamb;

public:
    G4(std::vector<double> p);
    double function(double rij, double rik, double rjk);
    double calculate();
};


class G5 : SymmetricFunction {
private:
    double cost, eta, zeta, lamb;

public:
    G5(std::vector<double> p);
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
