#include "lib.h"
#include <iostream>

double * ReadFromFile(const char * filename, int N) {
	double * r = new double[N];
	int k = 0;
	ifstream in(filename);
	if(!in) {
		cout << "Non-existing file" << endl;
		exit(0);
	} else {
		for(int i = 0; i < N; i++) in >> r[i];
	}
	return r;
}

void PrintDataOnFile(const char * filename, double * data, int N) {
	ofstream out(filename);
	for(int k = 0; k < N; k++) out << data[k] << endl;
}

double Error(double AV[], double AV2[], int n) {
    if (n == 0 || n == 1) return 0;
    else return sqrt((AV2[n]-AV[n]*AV[n])/double(n));
}

double Error2(double AV[], double AV2[], int n) {
    if (n == 0 || n == 1) return 0;
    else return sqrt((AV2[n]-AV[n]*AV[n])/(double(n)*double(n-1)));
}

double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

void scambiaByRef(double &a, double &b) {
  double temp;
  a = temp;
  a = b;
  b = temp;
}

void selection_sort(double * vec , int size) {

  int imin = 0;
  for(int j = 0; j < size - 1; j++) {
    imin = j;
    for(int i = j + 1; i < size; i++) {
      if(vec[i] < vec[imin]) {
        imin = i;
      }
    }
    scambiaByRef(vec[j], vec[imin]);
 }
}
 
void DataBlocking(int N, const char* fileMEAN, const char* fileERR, double ave[], double av2[]) {
	double sum_prog[N], su2_prog[N], err_prog[N];
	ofstream mean(fileMEAN), err(fileERR);
	for (int i = 0; i < N; i++) {
		sum_prog[i] = 0;
        su2_prog[i] = 0;
        for (int j = 0; j < i+1; j++) {
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
            sum_prog[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        sum_prog[i] = sum_prog[i]/(i+1);
        su2_prog[i] = su2_prog[i]/(i+1);
        err_prog[i] = Error(sum_prog, su2_prog, i);
		//cout << sum_progCall[i] << "  " << err_progCall[i] << endl;
		mean << sum_prog[i] << endl;
		err << err_prog[i] << endl;
	}
}