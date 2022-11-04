#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>

using namespace std;

double * ReadFromFile(const char * , int );
void PrintDataOnFile(const char * , double * , int );
double Error(double [], double [], int );
double Error2(double [], double [], int );
double gaussian_box_muller();
void scambiaByRef(double &, double &) ;
void selection_sort( double * , int );
void DataBlocking(int , const char* , const char* , double [], double []);