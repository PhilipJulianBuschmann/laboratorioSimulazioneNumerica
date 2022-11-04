#include <iostream>
#include <fstream>
#include <cmath>
#include "lib.h"

using namespace std;

bool check(double r, double theta) {
	double near = 0, d = 0.05 /* distance */, l = 0.03 /* needle length */;
	for(int i = 0; i < 21; i++) { //20 is the number of lines, check which one's closer
		if(abs(r - d*i) < d/2) {
			near = abs(r - d*i); //nearest line's distance
		}
	}
	if(near < (l/2)*sin(theta)) return true; //l/2 cause needle in center
	return false;
}

int main() {
	ifstream in("OutputRNG.txt");
	ofstream out("BuffonExp.txt");
	ofstream outE("BuffonExpErr.txt");
	int M = 50000, k = 0, N = 100;
	double r[M], theta[M], ave[N], av2[N], sum_prog[N], su2_prog[N], err_prog[N], l = 0.03, d = 0.05;
	while(!in.eof() && k < M) { //Load RNG numbers
        in >> r[k] >> theta[k]; //Pointing the needle in the center => theta in [0, PI/2]
		theta[k] = theta[k]*M_PI/2;
        k++;
    }
	in.close();
	for (int i = 0; i < N; i++) {
		double sum = 0;
		for(int j = i*M/N; j < (i+1)*M/N; j++) {
			if(check(r[j], theta[j])) {
				sum++;
			};
		}
		ave[i] = sum/(M/N);
		av2[i] = ave[i]*ave[i];
	}
	for (int i = 0; i < N; i++) {
        for (int j = 0; j < i+1; j++) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] = sum_prog[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        err_prog[i] = Error(sum_prog, su2_prog, i);
		err_prog[i] = 2*l/((sum_prog[i] - err_prog[i])*d) - 2*l/((sum_prog[i] + err_prog[i])*d); //difference between extremes
		out << 2*l/(sum_prog[i]*d) << endl; 
		outE << err_prog[i] << endl; 	
    }
}