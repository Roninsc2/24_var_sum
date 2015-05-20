#ifndef COORD
#define COORD
#include <cmath>

class TField
{

public:
    TField();
    void calculatedCoord(double dt, long long i);

private:
    void calculateSpeed_Up(double err[3], double errSpeedUp[3], double* errorA, int k);
    double TwoSum(double a, double b, double* error1, bool isNull);
    void Compensation(double  tempResult, double  input, double  *result, double  *error);


public:

    struct elem {
        double p[3];
    };
    struct elem speed;
    struct elem coord;
    struct elem speedUP;
    struct elem result;
    //2360591.5 == 1
    //5056.879 == 2
    //100 == 3
    //500 == 4
    //100000 == 5
    const double T = 100;
    const double g = K/(R*R);
    //6367444.7 == 1
    //384400000 == 2
    const double R = 6367444.7;
    const double K = 39.856e+13;

};

#endif // COORD

