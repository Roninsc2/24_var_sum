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

void TField::Split(float a, int s, float& a_hi, float& a_lo)
{
    float c = (pow(2, s) + 1)*a;
    float a_big = c - a;
    a_hi = c - a_big;
    a_lo = a - a_hi;
}

float TField::TwoProduct(float a, float b, float& err)
{
    float x = a*b;
    float a_hi, a_low, b_hi, b_low;
    Split(a, 16, a_hi, a_low);
    Split(b, 16, b_hi, b_low);
    float err1, err2, err3;
    err1 = x - (a_hi*b_hi);
    err2 = err1 - (a_low*b_hi);
    err3 = err2 - (a_hi*b_low);
    err += ((a_low * b_low) - err3);
    return x;

}

void TField::calculateSpeed_Up(std::vector<float> err,  std::vector<float>& errorSpeedUp,std::vector<float>& errorSpeedUp3,
                               std::vector<float>& errorSpeedUp2, std::vector<float>& errorA, std::vector<float>& errorA2 , int k)
{
    //вычисляю ускорение. vectorR - радиус вектор, вычисляется с учитываением погрешности координаты
    //vectorR = sqrt(координата_x^2 + координата_y^2) - общий вид
    std::vector<float> input;
    input.resize(3);
    for (int l=0; l < 2; l++) {
        float temp =  TwoProduct(2,coord.p[l], errorA[l]);
        temp = TwoProduct(temp, err[l], errorA[l]);
        input[l] = TwoSum(TwoProduct(err[l], err[l], errorA[l]),
                             temp, errorA2[l], false);
        input[l] = TwoSum(input[l], TwoProduct(coord.p[l],coord.p[l], errorA[l]), errorA2[l], false);
        errorA2[l] = TwoSum(errorA[l], errorA2[l], errorA[l], true);
        input[l] = TwoSum(input[l], errorA2[l], errorA2[l], true);
        errorA2[l] = 0;
        errorA[l] = 0;
    }
    float  vectorR = sqrt(input[0] + input[1]); //+ coord.p[2]*coord.p[2]
    if (vectorR == 0) {
        speedUP.p[0] = 0;
        speedUP.p[1] = 0;
        //speedUP.p[2] = 0;
    }
    else {
        //раскадываем ускорение по координатам и учитываем погрешность каждые 10 шагов
        float  vectorA = -1*K/(vectorR*vectorR*vectorR);
        for (int l = 0; l < 2; l++) {
             speedUP.p[l] = TwoSum(TwoProduct(vectorA, err[l], errorSpeedUp2[l]),
                                   TwoProduct(vectorA,coord.p[l], errorSpeedUp2[l]), errorSpeedUp[l], false);

             float temp = TwoSum(errorSpeedUp[l], errorSpeedUp2[l], errorSpeedUp2[l], true);
             errorSpeedUp3[l] = TwoSum(errorSpeedUp3[l],temp, errorSpeedUp[l], true);

        }
    }
}
//компенсация суммирования по Шевчуку
float TField::TwoSum(float a, float b, float& error1,  bool isNull)
{
    float x = a + b;
    float b_virt = x - a;
    float a_virt = x - b_virt;
    float b_roundoff = b - b_virt;
    float a_roudnoff = a - a_virt;
    float y = a_roudnoff + b_roundoff;
    if (!isNull) {
        error1 += y;
    } else {
        error1 = y;
    }
    return x;
}
//компенсация по Кэхэну, которая не используется
void  TField::Compensation(float  tempResult, float  input, float  *result, float  *error)
{
    float  temp0 = input - *error;
    float  temp1 = tempResult + temp0;
    *error = (temp1 - tempResult) - temp0;
    *result = temp1;
}


void  TField::calculatedCoord(float time, long long i)
{
    //задаем погрешности для скорости, расстояние, ускорения, погрешности от погрешности ускорения
    //в принципе, ее можно выпилить из функции, если мы ее обнуляем.
    std::vector<float> error(3,0);
    std::vector<float> errorSpeed(3,0);
    std::vector<float> errorSpeedUp(3,0);
    std::vector<float> error2(3,0);
    std::vector<float> errorSpeed2(3,0);
    std::vector<float> errorSpeedUp2(3,0);
    std::vector<float> error3(3,0);
    std::vector<float> errorSpeed3(3,0);
    std::vector<float> errorSpeedUp3(3,0);
    std::vector<float> errorA(3,0);
    std::vector<float> errorA2(3,0);

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
    float temp = 0.0;
    for (j = 0; j < i; j++) {
        calculateSpeed_Up(error, errorSpeedUp, errorSpeedUp3, errorSpeedUp2, errorA, errorA2, k);
            //каждые 10 шагов избалвляемся от погрешности
            for (int l = 0; l < 2; l++) {
                coord.p[l] = TwoSum(error3[l],coord.p[l],  error3[l], true);
                speed.p[l] = TwoSum(errorSpeed3[l], speed.p[l], errorSpeed3[l], true);
            }
        for (int l = 0; l < 2; l++) {
             //компенсация суммирования для координаты
            temp = TwoProduct(errorSpeedUp[l], 0.5, error2[l]);
            temp = TwoProduct(temp, time, error2[l]);
            temp = TwoProduct(temp, time, error2[l]);
            input = TwoSum(TwoProduct(errorSpeed[l],time, error2[l]), temp, error[l], false);
            input = TwoSum(input, TwoProduct(speed.p[l],time, error2[l]),  error[l], false);
            temp = TwoProduct(speedUP.p[l], 0.5, error2[l]);
            temp = TwoProduct(temp, time, error2[l]);
            temp = TwoProduct(temp, time, error2[l]);
            input = TwoSum(input, temp, error[l], false);
            temp = TwoSum(error[l], error2[l], error2[l], true);
            error3[l] = TwoSum(error3[l],temp, error[l], true);
            coord.p[l] = TwoSum(input, coord.p[l], error3[l], false);
             //speed update
            input = TwoSum(TwoProduct(errorSpeedUp[l],time, errorSpeed2[l]),
                           TwoProduct(speedUP.p[l] ,time, errorSpeed2[l]), errorSpeed[l], false);
            temp = TwoSum(errorSpeed[l], errorSpeed2[l], errorSpeed2[l], true);
            errorSpeed3[l] = TwoSum(errorSpeed3[l],temp, errorSpeed[l], true);
            speed.p[l] = TwoSum(input, speed.p[l], errorSpeed3[l], false);
            speedUP.p[l] = TwoSum(errorSpeedUp3[l], speedUP.p[l], errorSpeedUp3[l], true);

        }

        k++;
    }
    for (int l = 0; l < 2; l++) {
        temp = TwoSum(error[l], error2[l], error2[l], true);
        error3[l] = TwoSum(error3[l],temp, error[l], true);
        coord.p[l] = TwoSum(error3[l], coord.p[l],error3[l], true);

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
    for (long long i = 90000000; i < 100000000;i+=100000 ) {
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
        coord.p[0] = TwoSum(input, coord.p[0], true, error1X, error2X);
        coord.p[0] = coord.p[0] + error1X;

        //y
        input = speed.p[1]*time + error2SpeedY*time + speedUP.p[1]*0.5*time*time;
        coord.p[1] = TwoSum(input, coord.p[1], true, error1Y, error2Y);
        coord.p[1] =coord.p[1] + error1Y;
        //z
        //input = speed.p[2]*time + speedUP.p[2]*0.5*time*time;
        //Compensation(coord.p[2], input, coord.p[2], errorZ);

        //speed update
        input = speedUP.p[0] * time;
        speed.p[0] = TwoSum(input, speed.p[0], true, error1SpeedX, error2SpeedX);
        speed.p[0] = speed.p[0] + error1SpeedX;

        input = speedUP.p[1] * time;
        speed.p[1] = TwoSum(input, speed.p[1], true, error1SpeedY, error2SpeedY);
        speed.p[1] = speed.p[1] +  error1SpeedY;

        //input = speedUP.p[2] * time;
        //Compensation(speed.p[2], input, speed.p[2], errorSpeedZ);
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
