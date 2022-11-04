#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "error.h"

using namespace std;

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

double monte_carlo_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T, bool call) {
  double S_adjust = S * exp(T*(r-0.5*v*v)); 
  double S_cur = 0.0;   
  double payoff_sum = 0.0;

  for (int i=0; i<num_sims; i++) {
    double gauss_bm = gaussian_box_muller();  
    S_cur = S_adjust * exp(sqrt(v*v*T)*gauss_bm);
    if(call) payoff_sum += max(S_cur - K, 0.0); 
	else payoff_sum += max(K - S_cur, 0.0);
  }
  
  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}

double indPrice(double time, double init_price, int steps, double mu, double sigma, double Z[], int f) { 
	/*double dt = time/steps;
	for(int i = 0; i < steps; i++) {
		
        init_price = init_price*exp((mu - 0.5*sigma*sigma)*dt + sigma*gaussian_box_muller()*sqrt(dt));
        return init_price;
	}*/
     //Per qualche motivo ottenevo un risultato sbagliato con questo metodo
	
	double S[steps];
	S[0] = init_price;  
	double dt = time/steps;
	for(int i = 0; i < steps; i++) {
		S[i+1] = S[i]*exp((mu - 0.5*sigma*sigma)*dt + sigma*gaussian_box_muller()*sqrt(dt));
	}
    return S[steps-1];
}

double MaxZero(double a) {
    if(a > 0.0) return a;
    return 0.0;
}



int main() {
    int M=50000, N=100, L = M/N, k = 0;
    double x[M], y[M], r[M], S0 = 100, t = 0, T = 1, K = 100, rf = 0.1, vol = 0.25;
    double aveCall[N], avePut[N], av2Call[N], av2Put[N], sum_progCall[N], su2_progCall[N], sum_progPut[N], su2_progPut[N], err_progCall[N], err_progPut[N];
    double sostave[N];
	ifstream inx("RNGx.txt");
    ifstream iny("RNGy.txt");
    ofstream out("Ga.out");
	ofstream dirP("DirettoPut.txt"), dirPE("DirettoPutErr.txt");
	ofstream dirC("DirettoCall.txt"), dirCE("DirettoCallErr.txt");
	ofstream idirP("inDirettoPut.txt"), idirPE("inDirettoPutErr.txt");
	ofstream idirC("inDirettoCall.txt"), idirCE("inDirettoCallErr.txt");
    while(!inx.eof() && !iny.eof() && k < M) { //Carico i numeri RNG
        inx >> x[k];
        iny >> y[k];
        r[k] = Gaussian(0, 1, x[k], y[k]);
        out << r[k] << endl;
        k++;
    }
    //METODO DIRETTO
    for (int i = 0; i < N; i++) {
        aveCall[i] = monte_carlo_price(L, S0, K, rf, vol, T, true);//true or false stands for call or put
        av2Call[i] = aveCall[i]*aveCall[i];
        avePut[i] = monte_carlo_price(L, S0, K, rf, vol, T, false);
        av2Put[i] = avePut[i]*avePut[i];
    }
    for (int i = 0; i < N; i++) {
		sum_progCall[i] = 0;
        su2_progCall[i] = 0;
        sum_progPut[i] = 0;
        su2_progPut[i] = 0;
        for (int j = 0; j < i+1; j++) {
            sum_progCall[i] += aveCall[j];
            su2_progCall[i] += av2Call[j];
            sum_progPut[i] += avePut[j];
            su2_progPut[i] += av2Put[j];
        }
        sum_progCall[i] = sum_progCall[i]/(i+1);
        su2_progCall[i] = su2_progCall[i]/(i+1);
        err_progCall[i] = Error(sum_progCall, su2_progCall, i);
		//cout << sum_progCall[i] << "  " << err_progCall[i] << endl;
		dirC << sum_progCall[i] << endl;
		dirCE << err_progCall[i] << endl;

        sum_progPut[i] = sum_progPut[i]/(i+1);
        su2_progPut[i] = su2_progPut[i]/(i+1);
        err_progPut[i] = Error(sum_progPut, su2_progPut, i);
		//cout << sum_progPut[i] << "  " << err_progPut[i] << endl; 
		dirP << sum_progPut[i] << endl;
		dirPE << err_progPut[i] << endl;		
    }
    
    //METODO INDIRETTO
    for (int i = 0; i < N; i++) {
        double sumCall = 0;
        double sumPut = 0;
        for (int j = 0; j < L; j++) {
            k = j+i*L;
            sumCall += exp(-rf*T)*MaxZero(indPrice(T, S0, 100, rf, vol, r, k) - K);
            sumPut += exp(-rf*T)*MaxZero(K - indPrice(T, S0, 100, rf, vol, r, k));    
        }
        aveCall[i] = sumCall/L;
        av2Call[i] = aveCall[i]*aveCall[i];
        avePut[i] = sumPut/L;
        av2Put[i] = avePut[i]*avePut[i];
    }
    //indPrice(T, S0, 100, rf, vol, r);
    for (int i = 0; i < N; i++) {
		sum_progCall[i] = 0;
        su2_progCall[i] = 0;
        sum_progPut[i] = 0;
        su2_progPut[i] = 0;
        for (int j = 0; j < i+1; j++) {
            sum_progCall[i] += aveCall[j];
            su2_progCall[i] += av2Call[j];
            sum_progPut[i] += avePut[j];
            su2_progPut[i] += av2Put[j];
        }
        sum_progCall[i] = sum_progCall[i]/(i+1);
        su2_progCall[i] = su2_progCall[i]/(i+1);
        err_progCall[i] = Error(sum_progCall, su2_progCall, i);
		//cout << sum_progCall[i] << "  " << err_progCall[i] << endl;
		idirC << sum_progCall[i] << endl;
		idirCE << err_progCall[i] << endl;

        sum_progPut[i] = sum_progPut[i]/(i+1);
        su2_progPut[i] = su2_progPut[i]/(i+1);
        err_progPut[i] = Error(sum_progPut, su2_progPut, i);
		//cout << sum_progPut[i] << "  " << err_progPut[i] << endl; 
		idirP << sum_progPut[i] << endl;
		idirPE << err_progPut[i] << endl;		
		
    }
	cout << "Done" << endl;
}