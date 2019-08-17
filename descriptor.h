//
// Atomic Centered Symmetry Functions (ACSF)
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "structure.h"
#include "symmetry_function.h"
#include "symmetry_function_scaler.h"

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
    std::vector<std::vector<double>> calculate(AtomicStructure &structure);
    std::vector<double> calculate(AtomicStructure &structure, int atomIndex);
    std::vector<std::vector<double>> gradient(AtomicStructure &structure, int atomIndex_i, int atomIndex_j);
    double getGlobalCutOffRadius() const;

private:
    std::string centralElement;
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<std::string> listOfTwoBodyNeighborElement;
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    std::vector<std::string> listOfThreeBodyNeighborElement1, listOfThreeBodyNeighborElement2;
};

inline double calculateCosine(double rij, double rik, double drij[3], double drik[3]) 
{
    // cosine of angle between k--<i>--j atoms
    double cost = 0;
    for (int d=0; d<3; d++)
        cost += drij[d] * drik[d];
    return cost / rij / rik;
}

#endif //NNP_ACSF_H
