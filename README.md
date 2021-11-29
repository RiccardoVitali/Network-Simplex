# Network-Simplex
Vehicle Routing and Minimum Cost Flow problems (solved with Simplex method and Network Simplex Method)


ISTRUZIONI:

Per compilare il programma:
g++ -std=c++11 main.cpp -o <file_name>

Per lanciare l'eseguibile:
./<file_name>

COMMENTI:

Il file contiene due classi: Simplex e Network_Simplex. 
La prima risolve i problemi di Vehicle routing attraverso l'uso del Simplex method (ho scelto quello con il Big-M method e non quello a due fasi per trovare la prima soluzione feasible). Quindi risolve i problemi in forma generica cioè senza il vincolo del "flow balance".
La seconda risolve i problemi di Minimum cost flow con l'uso del Network simplex method.

Nel main ci sono degli esempi di problemi (tra cui quello del test) basterà sceglierne uno per testare il programma. 
Per la classe Simplex il si dovrà inserire la matrice dei vincoli, la funzione obiettivo (da minimizzare) e i segni dei vincoli. I termini noti della matrice dei vincoli devono essere positivi e il problema deve essere feasible.

Esempio:
vectordouble> objective = {0,-5,-7};  
    const vector<char> signs = {'<','<','<'};
    vector<vector<double>> table = {   {8,2,1},
                                       {9,1,2},
                                       {5,1,1}   };  
Corrisponde al problema:
min -5(x1) -7(x2):
	2(x1) + x2 <= 8
	x1 + 2(x2) <= 9
	x1 + x2 <= 5.


Per classe Network_Simplex si dovrà inserire la matrice dei grafo e il vettore supply/demand dei nodi.

Esempio:
const vector<vector<double>> rete = {       {0,0,15,27,32},
                                            {0,0,18,30,25},
                                            {0,0,0,0,0},
                                            {0,0,0,0,0},
                                            {0,0,0,0,0}   };

const vector<int> sup_dem{100,150,-50,-120,-80};
Corrisponde al problema del test.

La soluzione, insieme ai passaggi principali degli algoritmi, verranno mostrati a schermo durante l'esecuzione.

P.S. purtroppo c'è un piccolo bug che qualche volta esce fuori dopo aver scelto casualmente una soluzione iniziale (network simplex method), se accade ciò, bisogna rilanciare il programma. Mi tocca lasciarlo poiché ho finito il tempo a disposizione.
Per qualsiasi chiarimento sul codice o sulla logica dell'algoritmo o anche solamente per un feedback, scrivetemi.
