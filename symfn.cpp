//
// Created by hossein on 5/27/19.
//

#include <cmath>
#include "symfn.h"

//using namespace NNP_SF;

/* ----------------------------------------------------------------------
   setup for base symmetric function base
------------------------------------------------------------------------- */

SymmetricFunction::SymmetricFunction(double r_cutoff): r_cutoff(r_cutoff) {}

double SymmetricFunction::function(double r) { return 0; }

double SymmetricFunction::calculate() { return 0; }

double SymmetricFunction::fn_cutoff(double r)
{
    return ( cos(M_PI*r/r_cutoff) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetric function
------------------------------------------------------------------------- */

G0::G0(double r_cutoff): SymmetricFunction(r_cutoff) {}

double G0::function(double r) { return fn_cutoff(r); }

double G0::calculate() { return 0; }

/* ----------------------------------------------------------------------
   setup for G1 symmetric function
------------------------------------------------------------------------- */

G1::G1(double r_cutoff, double eta, double rs): SymmetricFunction(r_cutoff), eta(eta), rs(rs) {}

double G1::function(double rij)
{
    return exp( -eta * (rij-rs) * (rij-rs) ) * fn_cutoff(rij);
}

double G1::calculate() { return 0; }

/* ----------------------------------------------------------------------
   setup for G4 symmetric function
------------------------------------------------------------------------ */

G4::G4(double r_cutoff, double cost, double eta, double zeta, double lamb):
        SymmetricFunction(r_cutoff), cost(cost), eta(eta), zeta(zeta), lamb(lamb) {}

double G4::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * fn_cutoff(rij) * fn_cutoff(rjk) * fn_cutoff(rik);
}

double G4::calculate() { return 0; }


/* ----------------------------------------------------------------------
   setup for G5 symmetric function
------------------------------------------------------------------------- */

G5::G5(double r_cutoff, double cost, double eta, double zeta, double lamb): SymmetricFunction(r_cutoff),
                                                                        cost(cost), eta(eta),
                                                                        zeta(zeta), lamb(lamb) {};
double G5::function(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * fn_cutoff(rij) * fn_cutoff(rjk);
}

double G5::calculate() { return 0; }

/* ----------------------------------------------------------------------
   setup for symmetric function
------------------------------------------------------------------------- */

SymmetricFunction *ACSF::make_SymmetricFunction(SymmetricFunctionType select)
{
    switch(select)
    {
        case SymmetricFunctionType::TG0 :
            //return new G0();
            break;
        case SymmetricFunctionType::TG1 :
            //return new G1();
            break;
        case SymmetricFunctionType::TG4 :
            //return new G4();
            break;
        case SymmetricFunctionType::TG5 :
            //return new G5();
            break;
    }
}