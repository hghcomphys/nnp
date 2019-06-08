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
    ACSF();
    ~ACSF();
    void addTwoBodySymmetryFunction(TwoBodySymmetryFunction *symmetry_function); /*add two-body symmetry function*/
    void addThreeBodySymmetryFunction(ThreeBodySymmetryFunction *symmetry_function); /*add three-body symmetry function*/
    TwoBodySymmetryFunction& getTwoBodySF(const int index) const;
    ThreeBodySymmetryFunction& getThreeBodySF(const int index) const;
    int getNumberOfTwoBodySF() const;
    int getNumberOfThreeBodySF() const;
    int getNumberOfSF() const;
    void calculate(Atoms &configuration);
    std::vector<double>& getValues();
private:
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    std::vector<double> values;
};


#endif //NNP_ACSF_H
