//
// Created by hossein on 5/27/19.
//

#include <cmath>
#include "symfn.h"

//using namespace NNP_SF;

/* ----------------------------------------------------------------------
   setup for base symmetric function base
------------------------------------------------------------------------- */

SymmetricFunction::SymmetricFunction(double cutoff_radius): cutoff_radius(cutoff_radius) {}

double SymmetricFunction::cutoff_function(double r)
{
    return ( cos(M_PI*r/cutoff_radius) + 1.0 ) * 0.5;
}

/* ----------------------------------------------------------------------
   setup for G0 symmetric function
------------------------------------------------------------------------- */

G0::G0(std::vector<double> p): SymmetricFunction(p[0]) {}

double G0::descriptor(double rij)
{
    return cutoff_function(rij);
}

double G0::descriptor(double rij, double rik, double jk) { return 0; }

double G0::calculate() {}

/* ----------------------------------------------------------------------
   setup for G1 symmetric function
------------------------------------------------------------------------- */

G1::G1(std::vector<double> p): SymmetricFunction(p[0]), eta(p[1]), rs(p[2]) {}

double G1::descriptor(double rij)
{
    return exp( -eta * (rij-rs) * (rij-rs) ) * cutoff_function(rij);
}

double G1::descriptor(double rij, double rik, double jk) { return 0; }

double G1::calculate() {}

/* ----------------------------------------------------------------------
   setup for G4 symmetric function
------------------------------------------------------------------------ */

G4::G4(std::vector<double> p): SymmetricFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G4::descriptor(double rij) { return 0; }

double G4::descriptor(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik + rjk*rjk) );
    return res * cutoff_function(rij) * cutoff_function(rjk) * cutoff_function(rik);
}

double G4::calculate() {}


/* ----------------------------------------------------------------------
   setup for G5 symmetric function
------------------------------------------------------------------------- */

G5::G5(std::vector<double> p): SymmetricFunction(p[0]), cost(p[1]), eta(p[2]), zeta(p[3]), lamb(p[4]) {}

double G5::descriptor(double rij) { return 0; }

double G5::descriptor(double rij, double rik, double rjk)
{
    double res =  pow(2.0, 1.0-zeta) * pow(1.0+lamb*cost, zeta) * exp( -eta * (rij*rij + rik*rik) );
    return res * cutoff_function(rij) * cutoff_function(rjk);
}

double G5::calculate() {}

/* ----------------------------------------------------------------------
   setup for symmetric function
------------------------------------------------------------------------- */

//ACSF::ACSF() {}

void ACSF::add(SymmetricFunctionType select, std::vector<double> p)
{
    switch(select)
    {
        case SymmetricFunctionType::TG0 :
            list_of_symmetric_functions.push_back( new G0(p) );
            break;
        case SymmetricFunctionType::TG1 :
            list_of_symmetric_functions.push_back( new G1(p) );
            break;
        case SymmetricFunctionType::TG4 :
            list_of_symmetric_functions.push_back( new G4(p) );
            break;
        case SymmetricFunctionType::TG5 :
            list_of_symmetric_functions.push_back( new G5(p) );
            break;
    }
}

//ACSF::~ACSF()
//{
////    for (auto &sm: list_of_symmetric_functions)
////        delete sm;
//}