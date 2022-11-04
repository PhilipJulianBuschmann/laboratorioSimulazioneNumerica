#include "RNGExtend.h"

double cdfExpo(double y, double lambda) {
	return -(1/lambda)*log(1-y);
}

double cdfCauc(double y, double mu, double gamma) {
	return mu + gamma*tan(M_PI*(y - 0.5));
}

