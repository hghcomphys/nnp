//
// Cutoff Function
//

#ifndef NNP_CUTOFFFUNCTION_H
#define NNP_CUTOFFFUNCTION_H


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


#endif //NNP_CUTOFFFUNCTION_H
