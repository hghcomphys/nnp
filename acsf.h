//
// Atomic Centered Symmetry Functions (ACSF)
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "symmetryfunction.h"
#include "atoms.h"


class ACSF {
public:
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    ACSF();
    ~ACSF();
    void addTwoBodySymmetryFunction(TwoBodySymmetryFunction *symmetry_function); /*add two-body symmetry function*/
    void addThreeBodySymmetryFunction(ThreeBodySymmetryFunction *symmetry_function); /*add three-body symmetry function*/
    std::vector<double> calculate(Atoms &configuration);
};


#endif //NNP_ACSF_H
