#include <iostream>
#include <fstream>
#include <cmath>
#include "error.h"

int main() {
	ifstream in("OutputRNG.txt");
	ifstream inx("RNGx.txt");
	ifstream iny("RNGy.txt");
	ofstream outM("OutputMedia.txt");
	ofstream outAR("OutputAR.txt");
	int M = 100000, N = 100, k = 0, L = M/N, count;
	double r[M], integMedia[N], errintegMedia[N], ave[N], av2[N], su2_prog[N], sum = 0;
	double x[M], y[M], integAR[N], errintegAR[N], su2_progAR[N];
	while(!in.eof() && !inx.eof() && !iny.eof() && k < M) { //Carico i numeri RNG
        in >> r[k];
        inx >> x[k];
        iny >> y[k];
        k++;
    }
    in.close();
    //Metodo della media con N=100 blocchi di L=10000 numeri
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < L; j++) {
            k = j+i*L;
            sum += Integranda(r[k]);
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
        outM << integMedia[i] << endl;
        outM << errintegMedia[i] << endl;
    }
	
	//Metodo importance sampling con parabola 
    for (int i = 0; i < N; i++) {
        double sum = 0;
        count = 0;
        for (int j = 0; j < L; j++) {
            k = j+i*L;
            if(y[k] < Parabola(x[k])){ 
				sum += Integranda(x[k])*(M_PI/3)/Parabola(x[k]);
				count++;
			}
        }
        ave[i] = sum/count;
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
        cout << integAR[i] << ", " << errintegAR[i] << endl;
        outAR << integAR[i] << endl;
        outAR << errintegAR[i] << endl;
    }
	return 0;
}
