//
// Atomic Centered Symmetry Functions (ACSF)
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "symmetry_function.h"
#include "atoms.h"


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


#endif //NNP_ACSF_H
