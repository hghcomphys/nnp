//
// Created by hossein on 5/27/19.
//

#ifndef NNP_SYMMETRIC_FUNCTIONS_H
#define NNP_SYMMETRIC_FUNCTIONS_H



class SymmetricFunction {

public:
    double r_cutoff;
    double fn_cutoff(double r);
    SymmetricFunction(double r_cutoff): r_cutoff(r_cutoff) {};
    virtual double function(double r) {};
    virtual double calculate() {};

};


class G0 : public SymmetricFunction {

public:
    G0(double r_cutoff): SymmetricFunction(r_cutoff) {};
    double function(double r) { return fn_cutoff(r); };
    double calculate();

};


class G1 : public SymmetricFunction {

private:
    double eta, rs;

public:
    G1(double r_cutoff, double eta, double rs): SymmetricFunction(r_cutoff), eta(eta), rs(rs) {};
    double function(double r);
    double calculate();
};


class G4 : public SymmetricFunction {

private:
    double cost, eta, zeta, lamb;

public:
    G4(double r_cutoff, double cost, double eta, double zeta, double lamb): SymmetricFunction(r_cutoff),
                                                                            cost(cost), eta(eta),
                                                                            zeta(zeta), lamb(lamb) {};
    double function(double rij, double rik, double rjk);
    double calculate();
};



class G5 : public SymmetricFunction {

private:
    double cost, eta, zeta, lamb;

public:
    G5(double r_cutoff, double cost, double eta, double zeta, double lamb): SymmetricFunction(r_cutoff),
                                                                            cost(cost), eta(eta),
                                                                            zeta(zeta), lamb(lamb) {};
    double function(double rij, double rik, double rjk);
    double calculate();
};

#endif //NNP_SYMMETRIC_FUNCTIONS_H
