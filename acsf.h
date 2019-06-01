//
// Atomic Centered Symmetric Functions
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>

//namespace NNP_SF {

//typedef enum  {
//    TG0, TG1, TG4, TG5
//}  SymmetricFunctionType;


class ACSF {
public:
    double cutoff_radius;
    double cutoff_function(double r);
    ACSF(double cutoff_radius);
    virtual double function(double rij) {}; /*two-body*/
    virtual double function(double rij, double rik, double jk) {}; /*three-body*/
    virtual double calculate() {};
};


class G0 : public ACSF {
public:
    G0(std::vector<double> p);
    double function(double rij);
    double function(double rij, double rik, double jk);
    double calculate();
};


class G1 : public ACSF {
private:
    double eta, rs;
public:
    G1(std::vector<double> p);
    double function(double rij);
    double function(double rij, double rik, double jk);
    double calculate();
};


class G4 : public ACSF {
private:
    double cost, eta, zeta, lamb;
public:
    G4(std::vector<double> p);
    double function(double rij);
    double function(double rij, double rik, double rjk);
    double calculate();
};


class G5 : public ACSF {
private:
    double cost, eta, zeta, lamb;
public:
    G5(std::vector<double> p);
    double function(double rij);
    double function(double rij, double rik, double rjk);
    double calculate();
};


//}

#endif //NNP_ACSF_H
