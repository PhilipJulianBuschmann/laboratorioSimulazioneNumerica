#include <iostream>
#include <cmath>
#include <cstdlib>
#include "error.h"

double Error(double AV[], double AV2[], int n) {
    if (n == 0) return 0;
    else return sqrt((AV2[n]-AV[n]*AV[n])/double(n));
}

double Error2(double AV[], double AV2[], int n) {
    if (n == 0 || n == 1) return 0;
    else return sqrt((AV2[n]-AV[n]*AV[n])/(double(n)*double(n-1)));
}

double Integranda(double x) {
	return (M_PI/2)*cos(M_PI*x/2);
}

double Parabola(double x) {return -(M_PI/2)*x*x + M_PI/2;}

double Cumulativa(double x) {return -(M_PI/6)*x*x*x + (M_PI/2)*x;}

double Gaussian(double mu, double sigma, double r1, double r2) {
    r2 = r2*2*M_PI;
    double mag = sigma*sqrt(-2*log(1-r1));
    double x = mag*sin(r2) + mu;
    return x;
}
