//
// Cutoff Function
//

#ifndef NNP_CUTOFFFUNCTION_H
#define NNP_CUTOFFFUNCTION_H

#include <cmath>

class CutoffFunction {
public:
    CutoffFunction();
    void setCutoffRadius(double cutoffRadius);
    double fc(double r);
    double dfc(double r);

private:
    double rc;
    double inv_rc;
};

// TODO: other types of cutoff function
// TODO: using inline can reduce computational time
inline double CutoffFunction::fc(double r) 
{
    if ( r > rc ) return 0;
    
    // COS TYPE
    // return ( cos(M_PI * r * inv_rc) + 1.0 ) * 0.5; 
    
    // TANH TYPE
    double const tmp = tanh(1.0 - r * inv_rc);
    return tmp * tmp * tmp;
}

// TODO: other types of cutoff function
// TODO: using inline can reduce computational time
inline double CutoffFunction::dfc(double r) 
{
    if ( r > rc ) return 0;
    
    // COS TYPE
    // return -M_PI_2 * inv_rc * sin(M_PI * r * inv_rc);
    
    // TANH TYPE
    double tmp = tanh(1.0 - r * inv_rc);
    tmp *= tmp;
    return 3.0 * tmp * (tmp - 1.0) * inv_rc;
}

#endif //NNP_CUTOFFFUNCTION_H
