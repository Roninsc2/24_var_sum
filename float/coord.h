#ifndef COORD
#define COORD
#include <cmath>

class TField
{

public:
    TField();
    void calculatedCoord(float dt, long long i);

private:
    void calculateSpeed_Up(float err[3], float errSpeedUp[3], float* errorA, int k);
    float TwoSum(float a, float b, float* error1, bool isNull);
    void Compensation(float  tempResult, float  input, float  *result, float  *error);


public:

    struct elem {
        float p[3];
    };
    struct elem speed;
    struct elem coord;
    struct elem speedUP;
    struct elem result;
    //2360591.5 == 1 //655.719861
    //5056.879 == 2 // 1.40468861
    //100 == 3
    //500 == 4
    //100000 == 5
    const float T = 1.40468861;
    const float g = K/(R*R);
    //6367444.7 == 1
    //384400 == 2
    const float R = 6367.4447 ;
    //
    const float K = 5.1625595801952e+12 ;

};

#endif // COORD
