#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "error.h"

using namespace std;

double retta(double x) {return (1-x)*2;}

int main() {
	
	ifstream in("OutputRNG.txt");
	ifstream inx("RNGx.txt");
	ifstream iny("RNGy.txt");
	ofstream outM("OutputMedia.txt");
	ofstream outAR("OutputAR.txt");
	ofstream outER("OutputARerr.txt");

	int M = 500000, N = 100, k = 0, L = M/N, count;
	double integMedia[N], errintegMedia[N], ave[N], av2[N], su2_prog[N], sum = 0;
	double  integAR[N], errintegAR[N], su2_progAR[N];

    //Metodo della media con N=100 blocchi di L=10000 numeri
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < L; j++) {
            sum += Integranda(rand()/static_cast<double>(RAND_MAX));
        }
        ave[i] = sum/L;
        av2[i] = ave[i]*ave[i];
    }
	for (int i = 0; i < N; i++) {
        for (int j = 0; j < i+1; j++) {
            integMedia[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        integMedia[i] = integMedia[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        errintegMedia[i] = Error(integMedia, su2_prog, i);
        outM << integMedia[i] << " " << errintegMedia[i] << endl;
    }
	cout << "2.1 Done" << endl;
	//Metodo importance sampling con parabola 
    /*for (int i = 0; i < N; i++) {
        double sum = 0;
        count = 0;
        for (int j = 0; j < L; j++) {
			double y = rand()/static_cast<double>(RAND_MAX), x = rand()/static_cast<double>(RAND_MAX);
            if(y < Parabola(x)){ 
				sum += Integranda(x)/Parabola(x);
				count++;
			}
        }
        ave[i] = sum/count;
        av2[i] = ave[i]*ave[i];
    }*/
    //Metodo importance sampling con retta
    for (int i = 0; i < N; i++) {
        double sum = 0;
        count = 0;
        for (int j = 0; j < L; j++) {
            k = j+i*L;
            //x[k] = ((1/2)+sqrt(abs(M_PI-8*x[k]))/(2*sqrt(M_PI)))*x[k];
			double u = 1 - sqrt(-rand()/static_cast<double>(RAND_MAX)+1);
			if (retta(u) != 0) sum += Integranda(u)/retta(u);
		
        }
        ave[i] = sum/L;
        av2[i] = ave[i]*ave[i];
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i+1; j++) {
            integAR[i] += ave[j];
            su2_progAR[i] += av2[j];
        }
        integAR[i] = integAR[i]/(i+1);
        su2_progAR[i] = su2_progAR[i]/(i+1);
        errintegAR[i] = Error(integAR, su2_progAR, i);
        //cout << integAR[i] << ", " << errintegAR[i] << endl;
        outAR << integAR[i] << " " <<errintegAR[i] << endl;
    }
    cout << "2.2 Done" << endl;
	return 0;
}
