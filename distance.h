//
// Distance
//

#ifndef NNP_DISTANCE_H
#define NNP_DISTANCE_H

class Distance {
public:
    Distance();
    Distance(double r, double vec[3]); 
    void set_drVec(double vec[3], double factor=1.0);
    void set(double r, double vec[3], double factor=1.0);
//private:
    double dr;
    double inv_dr;
    double drVec[3];
};

inline void Distance::set_drVec(double vec[3], double factor)
{
    for (int i=0; i<3; i++)
        drVec[i] = vec[i] * factor;
}

#endif //NNP_DISTANCE_H
