#include <iostream>
#include <fstream>
#include <cmath>
#include "error.h"

class RW { //Random Walk class with position
	public:
		double pos[3];
};

void add_direction(RW &rr, double r) { //Changes position at random to a RW
	if(r*6 < 1) rr.pos[0]++;
	else if(r*6 < 2) rr.pos[0]--;
	else if(r*6 < 3) rr.pos[1]++;
	else if(r*6 < 4) rr.pos[1]--;
	else if(r*6 < 5) rr.pos[2]++;
	else rr.pos[2]--;
}

void add_free(RW &rr, double r, double p) { //Add free position at random
	rr.pos[0] += cos((r-.5)*2*M_PI)*sin(p*M_PI);
	rr.pos[1] += sin((r-.5)*2*M_PI)*sin(p*M_PI);
	rr.pos[2] += cos(p*M_PI);
}

double distance(RW r) { //Returns distance from [0,0,0] (start)
	return sqrt(r.pos[0]*r.pos[0] + r.pos[1]*r.pos[1] + r.pos[2]*r.pos[2]);
}

int main() {
	
	ifstream in("OutputRNG.txt"); //Random numbers 
	ofstream out("Lattice.txt"); //Lattice data
	ofstream outF("Free.txt"); //Free data
	
	int M = 10000 /*# of RWs*/, N = 100 /*# of blocks*/, L = M/N, k = 0, steps = 100;
	RW RWs[M]; //Generate M Random Walks
	double sum_prog[N], su2_prog[N], err_prog[N], ave[N], av2[N], r[100000], sum_progF[N], su2_progF[N], err_progF[N];
	
	while(!in.eof() && k < 100000) { //Load RNG numbers
        in >> r[k];
        k++;
    }
    in.close();
	
	for(int f = 0; f < steps; f++) { //LATTICE
		double sumtot = 0, sumerr = 0;
		for (int i = 0; i < N-1; i++) { //Loop between the N blocks, doing one add_direction for each RW
			double sum = 0;
			for (int j = 0; j < L; j++) { //Loop between a block of RWs
				k = j+i*L+f*10;
				add_direction(RWs[j+i*L], r[k]); //Lattice direction
				sum += distance(RWs[j+i*L]); 
				//cout << RWs[k].pos[0];
				//cout << RWs[k].pos[1];
				//cout << RWs[k].pos[2] << endl;
			
			}
			//cout << sum << endl;
			ave[i] = sum/L;
			//cout << ave[i] << endl;
			av2[i] = pow(ave[i], 2);
			sumtot += ave[i];
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < i+1; j++) {
				sum_prog[i] += ave[j];
				su2_prog[i] += av2[j];
			}
			sum_prog[i] = sum_prog[i]/(i+1);
			su2_prog[i] = su2_prog[i]/(i+1);
			err_prog[i] = Error(sum_prog, su2_prog, i);
		}
		out << sum_prog[N-1] << " " << err_prog[N-1] << endl;
	}
	cout << "Lattice done" << endl;
	RW RWsF[M];
	for(int f = 0; f < steps; f++) { //FREE
		double sumtot = 0, sumerr = 0;
		for (int i = 0; i < N-1; i++) { //Loop between the N blocks, doing one add_direction for each RW
			double sum = 0;
			for (int j = 0; j < L; j++) { //Loop between a block of RWs
				k = j+i*L+f*10;
				add_free(RWsF[j+i*L], r[k], r[k+2]); //Lattice direction
				sum += distance(RWsF[j+i*L]); 
				//cout << RWsF[k].pos[0];
				//cout << RWsF[k].pos[1];
				//cout << RWsF[k].pos[2] << endl;
			
			}
			ave[i] = sum/L;
			av2[i] = pow(ave[i], 2);
			sumtot += ave[i];
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < i+1; j++) {
				sum_progF[i] += ave[j];
				su2_progF[i] += av2[j];
			}
			sum_progF[i] = sum_progF[i]/(i+1);
			su2_progF[i] = su2_progF[i]/(i+1);
			err_progF[i] = Error(sum_progF, su2_progF, i);
		}
		outF << sum_progF[N-1] << " " << err_progF[N-1] << endl;
	}
	cout << "Free done" << endl;
	return 0;
}
