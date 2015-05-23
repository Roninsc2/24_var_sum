#include <iostream>
#include "coord.h"
#include <fstream>

using namespace std;

ofstream fout;

TField::TField()
{
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    //speed.p[2] = 0;
}

void TField::calculateSpeed_Up(float err[3], float errSpeedUp[3], float* errorA, int k)
{
    float input = TwoSum(err[0]*err[0], err[1]*err[1], errorA, false);
    input = TwoSum(input, 2*coord.p[0]*err[0], errorA, false);
    input = TwoSum(input, 2*coord.p[1]*err[1], errorA, false);
    input = TwoSum(input, coord.p[0]*coord.p[0], errorA, false);
    input = TwoSum(input, coord.p[1]*coord.p[1], errorA, false);
    input = TwoSum(input, *errorA, errorA, true);
    float  vectorR = sqrt(input); //+ coord.p[2]*coord.p[2]
    if (vectorR == 0) {
        speedUP.p[0] = 0;
        speedUP.p[1] = 0;
        //speedUP.p[2] = 0;
    }
    else {
        float  vectorA = -1*K/vectorR;
        vectorA /= vectorR;
        vectorA /= vectorR;
        for (int l = 0; l < 2; l++) {
             speedUP.p[l] = TwoSum(coord.p[l]*vectorA,vectorA *err[l], &errSpeedUp[l], false);
        }
        if (k == 1) {
            for (int l = 0; l < 2; l++) {
                //std::cerr << errSpeedUp[l] << std::endl;
                speedUP.p[l] = TwoSum(speedUP.p[l], errSpeedUp[l], &errSpeedUp[l], true);
            }
        }
        //speedUP.p[2] = coord.p[2] * vectorA;
    }
}

float TField::TwoSum(float a, float b, float* error1, bool isNull)
{
    float x = a + b;
    float b_virt = x - a;
    float a_virt = x - b_virt;
    float b_roundoff = b - b_virt;
    float a_roudnoff = a - a_virt;
    float y = a_roudnoff + b_roundoff;
    if (!isNull) {
        *error1 += y;
    } else {
        *error1 = y;
    }
    return x;
}

void  TField::Compensation(float  tempResult, float  input, float  *result, float  *error)
{
    float  temp0 = input - *error;
    float  temp1 = tempResult + temp0;
    *error = (temp1 - tempResult) - temp0;
    *result = temp1;
}


void  TField::calculatedCoord(float time, long long i)
{
    float error[3] = {0,0,0};
    float errorSpeedUp[3] = {0,0,0};
    float errorSpeed[3] = {0,0,0};
    float errorA = 0;

    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    float  input = 0;
    long long j = 0;
    speedUP.p[0] = 0;
    speedUP.p[1] =0;
    speedUP.p[2] = 0;
    int k = 0;
    for (j = 0; j < i; j++) {
        calculateSpeed_Up(error, errorSpeedUp, &errorA, k);
        if (k == 1) {
            for (int l = 0; l < 2; l++) {
                coord.p[l] = TwoSum(coord.p[l], error[l], &error[l], true);
            }
            for (int l = 0; l < 2; l++) {
                speed.p[l] = TwoSum(speed.p[l], errorSpeed[l], &errorSpeed[l], true);
            }
            k=0;
        }
        //компенсация суммирования для координаты

        //x
        for (int l = 0; l < 2; l++) {
            input = TwoSum(speed.p[l]*time, errorSpeed[l]*time, &error[l], false);
            input = TwoSum(speedUP.p[l]*0.5*time*time, input, &error[l], false);
            input = TwoSum(input, errorSpeedUp[l]*0.5*time*time, &error[l], false);
            coord.p[l] = TwoSum(input, coord.p[l], &error[l], false);
        }

        //speed update

        for (int l = 0; l < 2; l++) {
            input = TwoSum(speedUP.p[l] * time, errorSpeedUp[l] * time, &errorSpeed[l], false);
            speed.p[l] = TwoSum(input, speed.p[l], &errorSpeed[l], false);
        }
        k++;
    }
    for (int l = 0; l < 2; l++) {
        coord.p[l] = TwoSum(coord.p[l], error[l], &error[l], true);
    }
}
int main()
{
    fout.open("float_21.txt");
    TField field;
    float resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << field.K << std::endl;
    float resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << resultx << "\t" << resulty;
    for (long long i = 17000000; i < 20000000;i+=100000 ) {
        float time = field.T/(float)i;
        field.calculatedCoord(time, i);
        //fout << field.coord.p[0] << "\t" << field.coord.p[1] << "\t";
        //fout  << (field.coord.p[0] - resultx)/resultx << "\t" << (field.coord.p[1] - resulty)/resulty << "\t";
        fout << sqrt((-resultx + field.coord.p[0])*(-resultx + field.coord.p[0])
                     + (-resulty + field.coord.p[1])*(-resulty + field.coord.p[1]))/field.R
             << "\t";
        fout << i << std::endl;
    }
    fout.close();
    return 0;
}


/*#include <iostream>
#include "coord.h"
#include <fstream>

using namespace std;

ofstream fout;

TField::TField()
{
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    //speed.p[2] = 0;
}

void TField::calculateSpeed_Up(float errX, float errY, float errZ)
{
    float  vectorR = sqrt((coord.p[0]+errX)*(coord.p[0]+errX) + (coord.p[1]+errY) * (coord.p[1]+errY)); //+ coord.p[2]*coord.p[2]
    if(vectorR == 0)
    {
        speedUP.p[0] = 0;
        speedUP.p[1] = 0;
        //speedUP.p[2] = 0;
    }
    else
    {
        float  vectorA = -1*K/vectorR;
        vectorA /= vectorR;
        vectorA /= vectorR;
        speedUP.p[0] = coord.p[0]*vectorA + vectorA *errX;
        speedUP.p[1] = coord.p[1] * vectorA + vectorA*errY;
        //speedUP.p[2] = coord.p[2] * vectorA;
    }
}

float TField::TwoSum(float a, float b, bool re, float* error1, float* error2)
{
    float x = a + b;
    float b_virt = x - a;
    float a_virt = x - b_virt;
    float b_roundoff = b - b_virt;
    float a_roudnoff = a - a_virt;
    float y = a_roudnoff + b_roundoff;
    if(re) {
        *error1 = y;
        fout << y << "\t";
        return TwoSum(x, y, false, error1, error2);
    } else {
        fout << y << std::endl;
        *error2 += y;
        return x;
    }
}

void  TField::Compensation(float  tempResult, float  input, float  *result, float  *error)
{
    float  temp0 = input - *error;
    float  temp1 = tempResult + temp0;
    *error = (temp1 - tempResult) - temp0;
    *result = temp1;
}


void  TField::calculatedCoord(float time, long long i)
{
    float error1X = 0;
    float error1Z = 0;
    float error1Y = 0;
    float error1SpeedX = 0;
    float error1SpeedY = 0;
    float error1SpeedZ = 0;
    float error2X = 0;
    float error2Z = 0;
    float error2Y = 0;
    float error2SpeedX = 0;
    float error2SpeedY = 0;
    float error2SpeedZ = 0;
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    float  input = 0;
    long long j = 0;
    speedUP.p[0] = 0;
    speedUP.p[1] =0;
    speedUP.p[2] = 0;
    long long k = 0;
    for(j = 0; j < i; j++) {
        if (k == 100) {
            coord.p[0] += error2X;
            coord.p[1] += error2Y;
            //coord.p[2] += error2Z;

            speed.p[0] += error2SpeedX;
            speed.p[1] += error2SpeedY;
            error2X = 0;
            error2Y = 0;
            error2SpeedX = 0;
            error2SpeedY = 0;
            k=0;
        }
        calculateSpeed_Up(error2X, error2Y, error2Z);
        //компенсация суммирования для координаты

        //x
        input = speed.p[0]*time + error2SpeedX*time + speedUP.p[0]*0.5*time*time;
        coord.p[0] = TwoSum(input, coord.p[0], true, &error1X, &error2X);
        coord.p[0] = coord.p[0] + error1X;

        //y
        input = speed.p[1]*time + error2SpeedY*time + speedUP.p[1]*0.5*time*time;
        coord.p[1] = TwoSum(input, coord.p[1], true, &error1Y, &error2Y);
        coord.p[1] =coord.p[1] + error1Y;
        //z
        //input = speed.p[2]*time + speedUP.p[2]*0.5*time*time;
        //Compensation(coord.p[2], input, &coord.p[2], &errorZ);

        //speed update
        input = speedUP.p[0] * time;
        speed.p[0] = TwoSum(input, speed.p[0], true, &error1SpeedX, &error2SpeedX);
        speed.p[0] = speed.p[0] + error1SpeedX;

        input = speedUP.p[1] * time;
        speed.p[1] = TwoSum(input, speed.p[1], true, &error1SpeedY, &error2SpeedY);
        speed.p[1] = speed.p[1] +  error1SpeedY;

        //input = speedUP.p[2] * time;
        //Compensation(speed.p[2], input, &speed.p[2], &errorSpeedZ);
        k++;
    }
    coord.p[0] += error2X;
    coord.p[1] += error2Y;
    //fout << error1X << "\t" << error2X << "\n";
}
long long main()
{
    fout.open("test1.txt");
    TField field;
    float resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    float resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << resultx << "\t" << resulty;
    fout << "R = " << field.R << "\t" << "Time = " << field.T << std::endl;
    fout << "NuberOfSteps;Error" << std::endl;
    for(long long i = 500000; i < 2000000;i+=1000 ) {
        float time = field.T/(float)i;
        field.calculatedCoord(time, i);
        //fout << field.coord.p[0] << "\t" << field.coord.p[1] << "\t";
        //fout  << (field.coord.p[0] - resultx)/resultx << "\t" << (field.coord.p[1] - resulty)/resulty << "\t";
        fout << i << ";";
        //fout << "=" << sqrt((-resultx + field.coord.p[0])*(-resultx + field.coord.p[0])
                    // + (-resulty + field.coord.p[1])*(-resulty + field.coord.p[1]))/field.R
            // << std::endl;
    }
    fout.close();
    return 0;
}*/
