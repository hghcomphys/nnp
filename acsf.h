//
// Atomic Centered Symmetry Functions (ACSF)
//

#ifndef NNP_ACSF_H
#define NNP_ACSF_H

#include <vector>
#include "symmfunc.h"
#include "atoms.h"


class Scaler {
public:
    Scaler(double sfMin, double sfMax, double sfMean, double sfSigma): sfMin(sfMin), sfMax(sfMax), sfMean(sfMean), sfSigma(sfSigma) {};
    ~Scaler() {};
    double sfMin, sfMax, sfMean, sfSigma;
};


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
    std::vector<std::vector<double>> calculateSF(Atoms &configuration);
    std::vector<double> calculateSF(Atoms &configuration, int atomIndex);

    void addScaler(const Scaler& newScaler) {  listOfScalers.push_back(newScaler); }
    int getNumberOfScalers() { return listOfScalers.size(); }
    void scaleSymmetryFunctions() { isScale = true; }

private:
    std::string centralElement;
    std::vector<TwoBodySymmetryFunction *> listOfTwoBodySF; /* factory method */
    std::vector<std::string> listOfTwoBodyNeighborElement;
    std::vector<ThreeBodySymmetryFunction *> listOfThreeBodySF; /* factory method */
    std::vector<std::string> listOfThreeBodyNeighborElement1, listOfThreeBodyNeighborElement2;
    std::vector<int> listOfTwoBodySFindex, listOfThreeBodySFindex;
    std::vector<std::vector<double>> values;

    std::vector<Scaler> listOfScalers;
    bool isScale;
};


#endif //NNP_ACSF_H
