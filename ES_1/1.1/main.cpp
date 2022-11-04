#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "error.h"

using namespace std;

int main() {
   	ofstream outC("OutputDataChiSq.txt"); //Output file containing chi squared
    ifstream in("OutputRNG.txt"); //File containing RNG numbers
    ofstream outM("OutputDataMean.txt"); //Output file containing mean and error
	ofstream outME("OutputDataMeanErr.txt");
	ofstream outV("OutputDataVariance.txt"); //Output file containing variance and error
	ofstream outVE("OutputDataVarianceErr.txt");
    int M = 100000; //Size of RNG array
    int N = 100; //Number of blocks
    int L = M/N; //Block size
    int k = 0;
    double * r = new double[M];
    double * ave = new double[N];
    double * av2 = new double[N];
    double * sum_prog = new double[N]; 
    double * su2_prog = new double[N];  
    double * err_prog = new double[N];  
    while(!in.eof()) { //Load RNG numbers
        in >> r[k];
        k++;
    }
    in.close();
    for (int i = 0; i < N; i++) { //0.1.1
        double sum = 0;
        for (int j = 0; j < L; j++) {
            k = j+i*L;
            sum += r[k];
        }
        ave[i] = sum/double(L);
        av2[i] = pow(ave[i], 2);
    }
    for (int i = 0; i < N; i++) { //Data blocking
        for (int j = 0; j < i+1; j++) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] = sum_prog[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        err_prog[i] = Error(sum_prog, su2_prog, i);
		outM << sum_prog[i] << endl;
		outME << err_prog[i] << endl;			
    }
    cout << "1.1 Done" << endl;
	for (int i = 0; i < N; i++) { //0.1.2
        double sum = 0;
        for (int j = 0; j < L; j++) {
            k = j+i*L;
            sum += pow((r[k] - 0.5), 2);
        }
        ave[i] = sum/double(L);
        av2[i] = pow(ave[i], 2);
    }
    for (int i = 0; i < N; i++) { //Data blocking
		sum_prog[i] = 0;
        su2_prog[i] = 0;
        for (int j = 0; j < i+1; j++) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] = sum_prog[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        err_prog[i] = Error(sum_prog, su2_prog, i);
		outV << sum_prog[i] << endl;
		outVE << err_prog[i] << endl;	
    }
    cout << "1.2 Done" << endl;
	//0.1.3
    int intervals = 100, throws = 10000, count = 0;
    double width = 1.0/100.0, sumchi = 0.0, prob = double(throws)/double(intervals);
    for(int j = 0; j < 100; j++) { //100 values of Chi
        sumchi = 0.0;
        for(int k = 0; k < intervals; k++) { //Iterate for every interval
            count = 0;
            for(int i = 0; i < throws; i++) { //Count how many numbers in range [width*i,width*(i+1)]
                double rng = rand()/static_cast<double>(RAND_MAX);
                if(rng < width*(k+1) && rng > width*k) count++;
            }
            ave[k] = double(count); //Just making count a double
            sumchi += (ave[k] - prob)*(ave[k]-prob)/prob;
        }
        av2[j] = sumchi; //Collect Xj in an array
        outC << av2[j] << endl;  //Print 
    }
    cout << "1.3 Done" << endl;

    return 0;
}
