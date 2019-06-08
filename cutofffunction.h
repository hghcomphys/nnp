//
// Created by hossein on 6/8/19.
//

#ifndef NNP_CUTOFFFUNCTION_H
#define NNP_CUTOFFFUNCTION_H


class CutoffFunction {
public:
    CutoffFunction();
    void setCutoffRadius(double cutoffRadius);
    double fc(double r);

private:
    double rc;
    double inv_rc;
};


#endif //NNP_CUTOFFFUNCTION_H