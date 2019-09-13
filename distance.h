//
// Distance
//

#ifndef NNP_DISTANCE_H
#define NNP_DISTANCE_H

class Distance {
public:
    Distance();
    Distance(double r, const double vec[3]); 
    void set_drVec(const double vec[3], double factor=1.0);
    void set(double r, const double vec[3], double factor=1.0);
//private:
    double dr;
    double inv_dr;
    double drVec[3];
};

inline void Distance::set_drVec(const double vec[3], double factor)
{
    for (int d=0; d<3; d++)
        drVec[d] = vec[d] * factor;
}

#endif //NNP_DISTANCE_H
