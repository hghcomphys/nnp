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
    void addTwoBodySF(TwoBodySymmetryFunction *symmetry_function); /*add two-body symmetry function*/
    void addThreeBodySF(ThreeBodySymmetryFunction *symmetry_function); /*add three-body symmetry function*/
    TwoBodySymmetryFunction& getTwoBodySF(int index) const;
    ThreeBodySymmetryFunction& getThreeBodySF(int index) const;
    int getNumberOfTwoBodySF() const;
    int getNumberOfThreeBodySF() const;
    int getTotalNumberOfSF() const;
    void calculate(Atoms &configuration);
    std::vector<double>& getValues();
private:
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    std::vector<double> values;
};


#endif //NNP_ACSF_H
