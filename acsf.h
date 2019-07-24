//
// Atomic Centered Symmetry Functions (ACSF)
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "symmfunc.h"
#include "atoms.h"
#include "symmfuncscaler.h"


class ACSF {
public:
    ACSF(std::string element);
    ~ACSF();
    void addTwoBodySF(TwoBodySymmetryFunction *symmetryFunction, const std::string& neighborElement1); /*add two-body symmetry function*/
    void addThreeBodySF(ThreeBodySymmetryFunction *symmetryFunction,
                            const std::string& neighborElement1, const std::string& neibghorElement2); /*add three-body symmetry function*/
    TwoBodySymmetryFunction& getTwoBodySF(int index) const;
    ThreeBodySymmetryFunction& getThreeBodySF(int index) const;
    int getNumberOfTwoBodySF() const;
    int getNumberOfThreeBodySF() const;
    int getTotalNumberOfSF() const;
    const std::string& getCentralElement() const;
    std::vector<std::vector<double>> calculate(Atoms &configuration);
    std::vector<double> calculate(Atoms &configuration, int atomIndex);
    std::vector<std::vector<double>> gradient(Atoms &configuration, int atomIndex_i, int atomIndex_j);
    
private:
    std::string centralElement;
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<std::string> listOfTwoBodyNeighborElement;
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    std::vector<std::string> listOfThreeBodyNeighborElement1, listOfThreeBodyNeighborElement2;
};


#endif //NNP_ACSF_H
