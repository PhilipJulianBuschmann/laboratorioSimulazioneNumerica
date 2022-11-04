#ifndef ga_h
#define ga_h

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "random.h"

using namespace std;

class ga{
    
private:
    int m_N, m_type, m_pop; //Numero di citta, tipo di TSP, grandezza popolazione
    double m_p_cross, m_p_mut;
    Random m_rnd;
    vector<vector<double> > m_world; //Configurazione
    vector<vector<int> > m_population; //Popolazione

public:
    ga(int, int, int, Random, double, double); 
    ~ga(){;};

    void mutation(vector<int>&);
    void crossover();
    int selection();
    void Swap(int);
    double cost_function(vector<int>);
    vector< vector<int> > Get_m_pop();
    void check(vector<int>&); 
    void print();
    void print_world();
    void First_Gen();
    void Sort();
    void Next_Gen();
    vector<int> Get_Chrom(int);
    void print_half();
};

#endif 