#include "G_A.h"
#include "random.h"

ga :: ga(int n_cities, int type, int population, Random rnd, double p_cross, double p_mut){ //Generazione TSP, prima generazione, sorting e print

    m_N = n_cities; //Numero di città
    m_type = type; //Tipo di TSP: SQ = 0, CIRC = 1
    m_pop = population; //Numero di cromosomi in una popolazione
    m_rnd = rnd; //Generatore rnd
    m_p_cross = p_cross;
    m_p_mut = p_mut; //Probabilita di mutazione
    
    //Generazione del TSP
    if(m_type == 0){
        for (int i = 0; i<m_N; i++) {
            double theta, phi;
            vector <double> coordinates;
            
            theta = m_rnd.Rannyu(0,2*M_PI);
            coordinates.push_back(cos(theta));
            coordinates.push_back(sin(theta));
            
            m_world.push_back(coordinates);
        }
    }
    else{
        for (int i = 0; i<m_N; i++) {
            vector <double> coordinates;
            
            double x = m_rnd.Rannyu();
            double y = m_rnd.Rannyu();
            coordinates.push_back(x);
            coordinates.push_back(y);
            
            m_world.push_back(coordinates);
        }
    }
    
    print_world();
    
    First_Gen();

    Sort();
    
    ofstream out_1;
    out_1.open("pop_ini.dat");
    for (int i = 0; i<m_pop; i++) {
        for (int j = 0; j<m_N; j++) {
            out_1 << m_population[i][j] << " ";
        }
        out_1 << endl;
    }
    out_1.close();
}
    
void ga :: First_Gen(){
    vector<int> individual;
    for (int i = 0; i<m_N; i++) {
        individual.push_back(i);
    }
    for (int i=0; i<m_pop/2; i++) {
        m_population.push_back(individual);
    }
    
    for (int i = 0; i<m_pop/2; i++) { //Swaps
        int a = m_rnd.Rannyu(1,m_population.size());
        Swap(a);
    }

    int m = m_pop - m_population.size();
    for(int i = 0; i<m; i++){ //Spostamento di regioni
        int index = m_rnd.Rannyu(0, m_population.size());
        int cut = m_rnd.Rannyu(1, m_N-1);
        vector<int> v;
        v.push_back(0);
        for (int j=cut; j<m_N; j++) {
            v.push_back(m_population[index][j]);
        }
        for (int j=1; j<cut; j++) {
            v.push_back(m_population[index][j]);
        }
        m_population.push_back(v);
    }
    for (int i=1; i<m_N/2; i++) { //Specchiamento
        int index = m_rnd.Rannyu(0, m_population.size());
        int a = m_population[index][i];
        m_population[index][i] = m_population[index][m_N - i];
        m_population[index][m_N - i] = a;
    }
}

double ga :: cost_function(vector<int> cities){ //L^2
    double fitness=0;
    for (int j=0; j<m_N; j++) {
        for (int k=0; k<2; k++) {
            fitness += pow(m_world[cities[j]][k] - m_world[cities[(j+1)%m_N]][k],2);
        }
    }
    return fitness;
}

void ga :: Sort(){
    for (int i=0; i<m_pop-1; i++) {
        int p=i;
        for (int j=i+1; j<m_pop; j++) {
            if (cost_function(m_population[j])<cost_function(m_population[i]) && cost_function(m_population[j])<cost_function(m_population[p])) {
                p=j;
            }
        }
        swap(m_population[i],m_population[p]);
    }
}

int ga :: selection(){ //Questo operatore seleziona un cromosoma, assumendo popolazione gia sortata
    int index = m_pop*pow(m_rnd.Rannyu(),2);
    return index;
}

vector<int> ga :: Get_Chrom(int index){
    return m_population[index];
}

void ga :: Next_Gen(){
    vector< vector<int> > next_gen;
    for (int i=0; i<m_pop/2; i++) {//Ricreo una popolazione di figli. Da ogni crossover vengono fuori due figli -> m_N/2
        int index[2]; //Scelgo i genitori
        index[0] = selection();
        index[1] = selection();
        while (index[0] == index[1]) {
            index[1] = selection();
        }
        vector<int> son_0 = Get_Chrom(index[0]);
        vector<int> son_1 = Get_Chrom(index[1]);
        
        if(m_rnd.Rannyu()<m_p_cross){ //Crossover
            vector<int> gene;
            int cut = m_rnd.Rannyu(1,m_N);
            int h=0;
            for (int j=cut; j<m_N; j++) {
                gene.push_back(son_0[j]);
            }
            for (int j=1; j<m_N; j++) {
                for (int k=0; k<m_N-cut; k++) {
                    if (gene[k] == son_1[j]) {
                        son_0[cut + h] = gene[k];
                        h++;
                    }
                }
            }
            h=0;
            for (int j=cut; j<m_N; j++) {
                gene[j-cut] = son_1[j];
            }
            for (int j=1; j<m_N; j++) {
                for (int k=0; k<m_N-cut; k++) {
                    if (gene[k] == son_0[j]) {
                        son_1[cut + h] = gene[k];
                        h++;
                    }
                }
            }
        }
        mutation(son_0);
        mutation(son_1);
        check(son_0);
        check(son_1);
        next_gen.push_back(son_0);
        next_gen.push_back(son_1);
    }
    m_population = next_gen;
    Sort();
}

vector< vector<int> > ga :: Get_m_pop(){
    return m_population;
}

void ga :: check(vector<int>& son){
    for (int j=1; j<m_N; j++) { //Controllo che in ogni cromosoma non si ripetano i geni
        for (int h=0; h<m_N; h++) {
            if (son[j] - son[h] == 0 && j!=h) { 
                cerr << "Ripetizioni" << endl;
            }
        }
    }
}

void ga :: print(){
    ofstream out, fit;
    fit.open("fitness.dat",ios::app);
    out.open("population.dat",ios::app);
    for (int j = 0; j<m_N; j++) {
        out << m_population[0][j] << " ";
    }
    out << 0 << endl;
    
    fit << cost_function(m_population[0]) << endl;
    
    out.close();
    fit.close();
}

void ga :: mutation(vector<int>& son){ //Mutazione generica 
    int index[2] = {m_rnd.Rannyu(1,m_N),m_rnd.Rannyu(1,m_N)}; //Due numeri casuali diversi
    while (index[0]==index[1]) {
        index[1] = m_rnd.Rannyu(1,m_N);
    }

    if (m_rnd.Rannyu()<m_p_mut) { //Semplice swap di geni
       swap(son[index[0]],son[index[1]]);
    }
    
    if (m_rnd.Rannyu()<m_p_mut) { //Inversione di una regione
        vector<int> inverse;
        for (int i=index[1]; i>index[0]; i--) {
            inverse.push_back(son[i]);
        }
        for (int i=index[0]+1; i<index[1]; i++) {
            son[i] = inverse[i-index[0]];
        }
    }

    if (m_rnd.Rannyu()<m_p_mut) { //Ciclo 
        int n = m_rnd.Rannyu(2, m_N);
        int c = m_rnd.Rannyu(1, m_N-n);
        int first = son[c];
        for (int i=0; i<n; i++) {
            son[c+i] = son[c+i+1];
            if (i==n-1) {
                son[i+c] = first;
            }
        }
    }

    if (m_rnd.Rannyu()<m_p_mut) { //Traslazione di lunghezza h
        vector<int> v;
        v.push_back(0);
        int h = m_rnd.Rannyu(1,m_N);
        for (int i=h; i<m_N; i++) {
            v.push_back(son[i]);
        }
        for (int i=1; i<h; i++) {
            v.push_back(son[i]);
        }
        for (int i=1; i<m_N; i++) {
            son[i] = v[i];
        }
    }
}

void ga :: print_world(){
    ofstream out_1;
    out_1.open("cities_coordinates.dat");
    for (int i = 0; i<m_world.size(); i++) {
        out_1 << m_world[i][0] << " " << m_world[i][1] << endl;
    }
    out_1.close();
}

void ga :: Swap(int a){ 
    int index[2] = {0,0};
    index[0] = m_rnd.Rannyu(1,m_N);
    index[1] = m_rnd.Rannyu(1,m_N);
    while (index[0] == index[1]) {
        index[1] = m_rnd.Rannyu(1,m_N);
    }
    swap(m_population[a][index[0]],m_population[a][index[1]]);
}

void ga :: print_half(){
    ofstream out;
    out.open("best_half.dat",ios::app);
    double L_2 = 0, L_2_square = 0;
    for (int i=0; i<m_pop/2; i++) {
        L_2 += cost_function(m_population[i]);
        L_2_square += L_2*L_2;
    }
    out << L_2/m_pop * 2 << " " << sqrt(abs(L_2_square*2/m_pop - pow(L_2/m_pop * 2,2))/(m_pop-1)) << endl;
    out.close();
}