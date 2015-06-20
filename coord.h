#ifndef COORD
#define COORD
#include <cmath>
#include <vector>

class TField
{

public:
    TField();
    void calculatedCoord(float dt, long long i);

public:
    void calculateSpeed_Up(float err[3],  float(& errorSpeedUp)[3],float(& errorSpeedUp3)[3],
                           float(&errorSpeedUp2)[3], float(&errorA)[3], float(& errorA2)[3] , int k);
    float TwoSum(float a, float b, float& error1, bool isNull);
    void Split(float a, int s, float& a_hi, float& a_lo);
    float TwoProduct(float a, float b, float& err);
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
    //100 == 3 //0.027
    //500 == 4 //0.138
    //100000 == 5 //27.7
    const float T = 1.40468861;
    const float g = K/(R*R);
    //6367444.7 == 1 //6367.4447
    //384400000 == 2  //384400
    const float R = 6367.4447;
    //
    const float K = 5.1625595801952e+12 ;

};

#endif // COORD

