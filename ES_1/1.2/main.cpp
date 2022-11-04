#include <iostream>
#include <fstream>
#include "RNGExtend.h"

using namespace std;

int main() {
	ifstream in("OutputRNG2.txt"); //File containing RNG numbers
	ofstream outE1("OutputDataExp1.txt"); //Output file containing Exponential
	ofstream outE2("OutputDataExp2.txt");
	ofstream outE10("OutputDataExp10.txt");
	ofstream outE100("OutputDataExp100.txt");
	ofstream outC1("OutputDataCau1.txt"); //Output file containing Cauchy
	ofstream outC2("OutputDataCau2.txt");
	ofstream outC10("OutputDataCau10.txt");
	ofstream outC100("OutputDataCau100.txt");
	ofstream outS1("OutputDataSTD1.txt"); //Output file containing Standard
	ofstream outS2("OutputDataSTD2.txt");
	ofstream outS10("OutputDataSTD10.txt");
	ofstream outS100("OutputDataSTD100.txt");
	int M = 1000000;
	int k = 0;
	double * r = new double[M];
	while(!in.eof() && k < M) { //Load RNG numbers
        in >> r[k];
        k++;
    }
	in.close();
	for(int j = 0; j < M; j+= 100) {
		outE1 << cdfExpo(r[j], 1) << endl;
		outC1 << cdfCauc(r[j], 0, 1) << endl;
		outS1 << r[j] << endl;
		outE2 << 0.5*(cdfExpo(r[j], 1) + cdfExpo(r[j+1], 1)) << endl;
		outC2 << 0.5*(cdfCauc(r[j], 0, 1) + cdfCauc(r[j+1], 0, 1)) << endl;
		outS2 << 0.5*(r[j] + r[j+1]) << endl;
		double Esum10 = 0, Csum10 = 0, Ssum10 = 0;
		for(int i = 0; i < 10; i++) {
			Esum10 += 0.1*cdfExpo(r[j+i], 1);
			Csum10 += 0.1*cdfCauc(r[j+i], 0, 1);
			Ssum10 += 0.1*r[j+i];
		}
		outE10 << Esum10 << endl;
		outC10 << Csum10 << endl;
		outS10 << Ssum10 << endl;	
		double Esum100 = 0, Csum100 = 0, Ssum100 = 0;
		for(int i = 0; i < 100; i++) {
			Esum100 += 0.01*cdfExpo(r[j+i], 1);
			Csum100 += 0.01*cdfCauc(r[j+i], 0, 1);
			Ssum100 += 0.01*r[j+i];
		}
		outE100 << Esum100 << endl;
		outC100 << Csum100 << endl;
		outS100 << Ssum100 << endl;
	}
	return 0;
}