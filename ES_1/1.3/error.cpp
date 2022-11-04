#include <iostream>
#include <cmath>
#include <cstdlib>
#include "error.h"

double Error(double * AV, double * AV2, int n) {
    if (n == 0 || n == 1) return 0;
    else return sqrt((AV2[n]-pow(AV[n], 2))/double(n*(n-1)));
}