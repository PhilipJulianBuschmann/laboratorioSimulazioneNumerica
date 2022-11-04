#include <iostream>
#include <cmath>
#include <cstdlib>
#include "error.h"

double Error(double AV[], double AV2[], int n) {
    if (n == 0) return 0;
    else return sqrt((AV2[n]-AV[n]*AV[n])/double(n));
}

double Integranda(double x) {
	return (M_PI/2)*cos(M_PI*x/2);
}

double Parabola(double x) {return -(M_PI/2)*x*x + M_PI/2;}

double Cumulativa(double x) {return -(M_PI/6)*x*x*x + (M_PI/2)*x;}
