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

void TField::calculateSpeed_Up(double err[3], double errSpeedUp[3], double* errorA, int k)
{
    double input = TwoSum(err[0]*err[0], err[1]*err[1], errorA);
    input = TwoSum(input, 2*coord.p[0]*err[0], errorA);
    input = TwoSum(input, 2*coord.p[1]*err[1], errorA);
    input = TwoSum(input, coord.p[0]*coord.p[0], errorA);
    input = TwoSum(input, coord.p[1]*coord.p[1], errorA);
    input = TwoSum(input, *errorA, errorA);
    double  vectorR = (double )sqrt(input); //+ coord.p[2]*coord.p[2]
    if (vectorR == 0)
    {
        speedUP.p[0] = 0;
        speedUP.p[1] = 0;
        //speedUP.p[2] = 0;
    }
    else
    {
        double  vectorA = -1*K/vectorR;
        vectorA /= vectorR;
        vectorA /= vectorR;
        speedUP.p[0] = TwoSum(coord.p[0]*vectorA,vectorA *err[0], &errSpeedUp[0]);
        speedUP.p[1] = TwoSum(coord.p[1] * vectorA,vectorA*err[1], &errSpeedUp[1]);
        if ( k == 100000) {
            speedUP.p[0] = TwoSum(speedUP.p[0], errSpeedUp[0], &errSpeedUp[0]);
            speedUP.p[1] = TwoSum(speedUP.p[1], errSpeedUp[1], &errSpeedUp[1]);
        }
        //speedUP.p[2] = coord.p[2] * vectorA;
    }
}

double TField::TwoSum(double a, double b, double* error1)
{
    double x = a + b;
    double b_virt = x - a;
    double a_virt = x - b_virt;
    double b_roundoff = b - b_virt;
    double a_roudnoff = a - a_virt;
    double y = a_roudnoff + b_roundoff;
    *error1 += y;
    return x;
}

void  TField::Compensation(double  tempResult, double  input, double  *result, double  *error)
{
    double  temp0 = input - *error;
    double  temp1 = tempResult + temp0;
    *error = (temp1 - tempResult) - temp0;
    *result = temp1;
}


void  TField::calculatedCoord(double time, long long i)
{
    double error[3] = {0,0,0};
    double errorSpeedUp[3] = {0,0,0};
    double errorSpeed[3] = {0,0,0};
    double errorA = 0;

    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    double  input = 0;
    long long j = 0;
    speedUP.p[0] = 0;
    speedUP.p[1] =0;
    speedUP.p[2] = 0;
    long long k = 0;
    for(j = 0; j < i; j++) {
        calculateSpeed_Up(error, errorSpeedUp, &errorA, k);
        if(k == 100000) {
            for(long long l = 0; l < 2; l++) {
                coord.p[l] = TwoSum(coord.p[l], error[l], &error[l]);
            }
            for(long long l = 0; l < 2; l++) {
                speed.p[l] = TwoSum(speed.p[l], errorSpeed[l], &errorSpeed[l]);
            }
            k=0;
        }
        //компенсация суммирования для координаты

        //x
        for(long long l = 0; l < 2; l++) {
            input = TwoSum(speed.p[l]*time, errorSpeed[l]*time, &error[l]);
            input = TwoSum(speedUP.p[l]*0.5*time*time, input, &error[l]);
            coord.p[l] = TwoSum(input, coord.p[l], &error[l]);
        }

        //speed update

        for(long long l = 0; l < 2; l++) {
            input = speedUP.p[l] * time;
            speed.p[l] = TwoSum(input, speed.p[l], &errorSpeed[l]);
        }
        k++;
    }
    for(long long l = 0; l < 2; l++) {
        coord.p[l] = TwoSum(coord.p[l], error[l], &error[l]);
    }
}
int main()
{
    fout.open("test1.txt");
    TField field;
    double resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    double resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << resultx << "\t" << resulty;
    for(long long i = 10000000; i < 20000000;i+=100000 ) {
        double time = field.T/(double)i;
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

void TField::calculateSpeed_Up(double errX, double errY, double errZ)
{
    double  vectorR = sqrt((coord.p[0]+errX)*(coord.p[0]+errX) + (coord.p[1]+errY) * (coord.p[1]+errY)); //+ coord.p[2]*coord.p[2]
    if(vectorR == 0)
    {
        speedUP.p[0] = 0;
        speedUP.p[1] = 0;
        //speedUP.p[2] = 0;
    }
    else
    {
        double  vectorA = -1*K/vectorR;
        vectorA /= vectorR;
        vectorA /= vectorR;
        speedUP.p[0] = coord.p[0]*vectorA + vectorA *errX;
        speedUP.p[1] = coord.p[1] * vectorA + vectorA*errY;
        //speedUP.p[2] = coord.p[2] * vectorA;
    }
}

double TField::TwoSum(double a, double b, bool re, double* error1, double* error2)
{
    double x = a + b;
    double b_virt = x - a;
    double a_virt = x - b_virt;
    double b_roundoff = b - b_virt;
    double a_roudnoff = a - a_virt;
    double y = a_roudnoff + b_roundoff;
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

void  TField::Compensation(double  tempResult, double  input, double  *result, double  *error)
{
    double  temp0 = input - *error;
    double  temp1 = tempResult + temp0;
    *error = (temp1 - tempResult) - temp0;
    *result = temp1;
}


void  TField::calculatedCoord(double time, long long i)
{
    double error1X = 0;
    double error1Z = 0;
    double error1Y = 0;
    double error1SpeedX = 0;
    double error1SpeedY = 0;
    double error1SpeedZ = 0;
    double error2X = 0;
    double error2Z = 0;
    double error2Y = 0;
    double error2SpeedX = 0;
    double error2SpeedY = 0;
    double error2SpeedZ = 0;
    coord.p[0] = R;
    coord.p[1] = 0;
    coord.p[2] = 0;
    speed.p[1] = sqrt(K/R);
    speed.p[0] = 0;
    double  input = 0;
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
    double resultx = cos(field.T*(field.speed.p[1]/field.R))*field.R;
    double resulty = sin(field.T*(field.speed.p[1]/field.R))*field.R;
    std::cerr << resultx << "\t" << resulty;
    fout << "R = " << field.R << "\t" << "Time = " << field.T << std::endl;
    fout << "NuberOfSteps;Error" << std::endl;
    for(long long i = 500000; i < 2000000;i+=1000 ) {
        double time = field.T/(double)i;
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
