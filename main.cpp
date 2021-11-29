/*
 - Mettersi nella cartella del file.

 - Per compilare il file usare il comando:
   g++ -o <name> main.cpp

 - Per runnare l'eseguibile usare il comando:
   ./<name>




*/




#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include <cmath>
#include <list>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stdlib.h>

using namespace std;

void print_vec(vector<double> vec){
    for(auto x : vec){
        cout << x << " ";
    }
    cout << "\n";
}
bool is_not_zero(const double &d) {
    const double epsilon = 0.0001;
    return fabs(d) > epsilon;
}

bool check_optimality(const vector<double> & v) {
    for (double x : v) {
        if (x < 0) return false;
    }
    return true;
}

pair<bool,vector<double>> gauss(vector< vector<double> > A) {
    int n = A.size();
    vector<double> x(n);
    for (int i=0; i<n; i++) {
        double max = abs(A[i][i]);
        int row = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > max) {
                max = abs(A[k][i]);
                row = k;
            }
        }
        for (int k=i; k<n+1;k++) {
            double tmp = A[row][k];
            A[row][k] = A[i][k];
            A[i][k] = tmp;
        }
        for (int k=i+1; k<n; k++) {
            if(is_not_zero(A[i][i])){
                double c = -A[k][i]/A[i][i];
                for (int j=i; j<n+1; j++) {
                    if (i==j) {
                        A[k][j] = 0;
                    } else {
                        A[k][j] += c * A[i][j];
                    }
                }
            }else{
                cout << "COLUMNS ARE LINEARLY DEPENDENT\n";
                return make_pair(false,x);
            }

        }
    }

    for (int i=n-1; i>=0; i--) {
        if(is_not_zero(A[i][i])){
            x[i] = A[i][n]/A[i][i];
            for (int k=i-1;k>=0; k--) {
                A[k][n] -= A[k][i] * x[i];
            }
        }
        else{
            cout << "COLUMNS ARE LINEARLY DEPENDENT\n";
            return make_pair(false,x);
        }

    }
    return make_pair(true,x);
}

int find_min_position(const vector<double> &x, const int &start, const int &end) {
    int min = x[start];
    int index = start;
    for (int i = start; i < end; i++) {
        if (min > x[i]) {
            min = x[i];
            index = i;
        }
    }
    return index;
}

class Simplex {
private:

    const int M = 1000000; // costante M per il big-M method
    const vector<double> objective; // funzione obiettivo
    vector<double> obj; // funzione obiettivo usata per il calcolo della soluzione
    const size_t variables = objective.size(); // numero dei coefficienti della funzione obiettivo
    const vector<char> signs; // vettore dei segni dei vincoli
    vector<vector<double>> table; // matrice dei vincoli
    const size_t constraints = signs.size(); // numeri di vincoli
    vector<int> artificial_variables; // posizione delle artificial variable
    bool unbounded = false;
    //bool unfeasible = false;
    map<int,int> variables_positions;
    //int variables_positions[];  // array che tiene le posizioni delle


public:

    //constructor
    Simplex(vector<double>  a, vector<char>  b, vector<vector<double>> & c) : objective(std::move(a)),signs(std::move(b)),table(std::move(c)){
        obj = objective;
    }  // equazione, segni dei vincoli e vincoli

    void print_table() {
        for (double & it : obj) {
            cout << it << " ";
        }
        cout << "\n";
        cout << "----------------------------------";
        for (auto & i : table) {
            cout << "\n";
            for (auto ii = 0; ii < obj.size(); ii++) {
                cout << i[ii] << " ";
            }
        }
        cout << "\n\n";
    }

    void add_slack_variables() {
        for (unsigned i = 0; i < constraints; i++) {
            for (unsigned j = 0; j < constraints; j++) {
                if (j == i && signs[i] == '<') {
                    table[i].push_back(1);
                } else if (j == i && signs[i] == '>') {
                    table[i].push_back(-1);
                } else {
                    table[i].push_back(0);
                }
            }
            obj.push_back(0);
        }
    }

    bool artificial_variable(){
        for (int i = 0; i < signs.size(); i++) {
            if (signs[i] == '>' || signs[i] == '=') {
                artificial_variables.push_back(i);
            }
        }
        return !artificial_variables.empty();
    }

    void modify_obj(){
        for (int i = 1; i < obj.size(); i++) {
            double temp=0;
            for (auto & x : table) {
                temp += x[i];
            }
            obj[i] -= temp * M;
        }
        for (int i = 0; i < artificial_variables.size(); i++) {
            obj.push_back(0);
        }
    }

    void add_artificial_variables(){
        for(int i=0;i<table.size();i++){
            for(int ii=0;ii<artificial_variables.size();ii++){
                if(i==artificial_variables[ii]){
                    table[i].push_back(1);
                }
                else{
                    table[i].push_back(0);
                }

            }
        }
    }

    int leaving_variable(const int &entry) {
        vector<int> non_zeros;
        for(int j=0;j<table.size();j++){
            if(is_not_zero(table[j][entry])){
                non_zeros.push_back(j);
            }
        }
        double max = table[non_zeros[0]][entry] / table[non_zeros[0]][0];
        int index = non_zeros[0];
        double temp;
        for (auto x : non_zeros) {
            temp = table[x][entry] / table[x][0];
            if (max < temp) {
                max = temp;
                index = x;
            }
        }
        if(non_zeros.empty()){
            cout << "colonna di tutti zeri";
            return 0;
        }
        return index;
    }

    void pivoting(const int &row, const int &col) {
        double temp = table[row][col];
        for (double & i : table[row]) {
            i = i / temp;
        }
        table[row][col] = 1;// forse si può anche togliere
        temp = obj[col];
        for (int i = 0; i < obj.size(); i++) {
            obj[i] -= temp * table[row][i];
        }
        obj[col] = 0; // forse si può anche togliere

        for (int i = 0;i < table.size(); i++) {
            if (i == row) continue;
            temp = table[i][col];
            if (is_not_zero(temp)) {   // sempre meglio essere cauti quando si lavora con i double perché 0 potrebbe non essere proprio 0 (0.00001)
                for (int ii = 0; ii < obj.size(); ii++) {
                    table[i][ii] -= temp * table[row][ii];
                }
                table[i][col] = 0;// forse si può anche togliere
            }
        }
    }

    bool check_unbounded(const int & entry){
        for(vector<double> x : table){
            if(x[entry]>0) return false;
        }
        return true;
    }
    /*
    void check_feasibility(){
        for(int val=obj.size()-artificial_variables.size();val<obj.size();val++){
            cout << "val " << val << "\n";
            for(int i=0;i<constraints;i++){
                cout << "variables_positions[i] " << variables_positions[i] << "\n";
                if(variables_positions[i]==val){
                    unfeasible = true;
                    return;
                }
            }

        }
    }
    */
    void simplex_algorithm(){
        print_table();
        add_slack_variables();
        print_table();
        if(artificial_variable()) {   // ho aggiunto delle artificial variables quindi devo fare il big-M method

            add_artificial_variables();
            modify_obj();
            print_table();
        }
        int entry;
        int out;
        while(!check_optimality(obj)){
            entry = find_min_position(obj, 1, obj.size());
            if(check_unbounded(entry)){
                unbounded = true;
                return;
            }

            out = leaving_variable(entry);
            variables_positions[out]=entry; //dove è entrata la nuova variabile
            pivoting(out,entry);
            print_table();
        }
        print_solution();
    }

    void print_solution(){
        if(unbounded){
            cout << "Solution is UNBOUNDED";
            return;
        }
        double res=0;
        map<int,double> map;
        for(int i=0;i<constraints;i++){
            if(variables_positions[i]<variables){
                res += table[i][0]*objective[variables_positions[i]];
                map.insert({variables_positions[i], table[i][0]});
            }
        }
        cout << "Solution is: " << res << "\n";
        cout << "with variables: \n";
        for(int j=1;j<variables;j++){
            auto it = map.find(j);
            if(it!=map.end()){
                cout << "x" << j << " = " << it->second << "\n";
            }
            else{
                cout << "x" << j << " = 0\n";
            }
        }
    }
};

class Network_Simplex {
private:
    struct pair_hash {
        inline std::size_t operator()(const std::pair<int,int> & v) const {
            return v.first^v.second;
        }
    };
    const vector<int> b; // vettore supply/demand
    const size_t nodes = b.size(); // numero di nodi
    const vector<vector<double>> costs; // costi del network

    vector<vector<double>> csc;//complementary slackness condition
    unordered_set<pair<int,int>, pair_hash> not_present;
    vector<double> opt;
    pair<int,int> enter;
    bool unbounded = false;

public:
    vector<list<pair<int,int>>> network; // network
    //constructor
    Network_Simplex(vector<vector<double>>  a, vector<int>  b) : b(std::move(b)), costs(std::move(a)) {}

    void print_solution(){
        cout << "SOLUTION:\n";
        for(int i=0;i<nodes;i++){
            for(auto const & x : network[i]){
                if(x.second>0){
                    cout  << i << "->" << x.first << ", weight: " << x.second << "\n";
                }

            }
        }
    }

    int get_cost(const int & x, const int & y){
        return costs[x][y];
    }

    void complementary_slackness_condition(){
        vector<vector<double>> vec(nodes-1, vector<double> (nodes));
        csc = vec;
        int row = 0;
        for(int i=0;i<nodes;i++){
            for(auto p : network[i]){
                if(p.second>0){
                    if(i!=nodes-1){
                        csc[row][i] = 1;
                    }
                    if(p.first!=nodes-1){
                        csc[row][p.first] = -1;
                    }
                    csc[row][nodes-1] = get_cost(i,p.first);
                    row++;
                }
            }
        }
        cout << "COMPLEMENTARY SLACKNESS CONDITION:\n";
        for(auto i : csc){
            for(int ii=0;ii<csc[0].size();ii++){
                cout << i[ii] << " ";
            }
            cout << "\n";
        }
        opt = gauss(csc).second;
        opt.push_back(0);
        cout << "OPTIMAL VECTOR:\n";
        for(auto x : opt){
            cout << x <<" ";
        }
        cout << "\n";
    }

    bool find_arc(const int & a, const int & b){
        for(auto const & x : network[a]){
            if(x.first==b) return true;
        }
        return false;
    }

    void not_present_creation(){
        for(int i=0;i<costs.size();i++){
            for(int ii=0;ii<costs.size();ii++){
                if(get_cost(i,ii)>0 && !find_arc(i,ii)){
                    not_present.insert(make_pair(i,ii));
                }
            }
        }
        cout << "NOT PRESENT ARCS:\n";
        for (auto i : not_present) {
            cout << i.first << " " << i.second << "\n";
        }
    }

    int find_flow(const int & x, const int & y){
        for(auto const & z : network[x]){
            if(z.first==y) return z.second;
        }
    }

    bool optimality_conditions(){
        for(auto const & x : not_present){
            const int first = x.first;
            const int second = x.second;
            if(opt[first] - opt[second] > get_cost(first,second)){
                enter = make_pair(first,second);
                cout << "ENTRY ARC:\n" << enter.first << " " << enter.second << "\n";
                return false;
            }
        }
        return true;
    }

    void add_entry_arc(){
        network[enter.first].push_back(make_pair(enter.second,1));
        network[enter.second].push_back(make_pair(enter.first,-1));
    }

    static int remove_arc(vector<list<pair<int,int>>> & net,const int & x, const int & y){
        int flow = 0;
        for(auto iter=net[x].begin();iter!=net[x].end();iter++){
            if(iter->first == y){
                flow = iter->second;
                net[x].erase(iter);
            }
        }
        for(auto iter=net[y].begin();iter!=net[y].end();iter++){
            if(iter->first == x) net[y].erase(iter);
        }
        return flow;
    }

    static void add_arc(vector<list<pair<int,int>>> & net, const int & x, const int & y, const int & flow){
        net[x].push_back(make_pair(y,flow));
        net[y].push_back(make_pair(x,-flow));
    }

    bool equal_to_enter(const int & x, const int & y){
        return x==enter.first && y==enter.second;
    }

    void remove_leaving_arc(){
        unordered_map<int,int> cycle_back;
        unordered_multimap<int,int> cycle;
        add_entry_arc();
        int flow;
        auto temp_network = network;
        for(int i=0;i<nodes;i++){
            for(auto const & x : network[i]){ // se tolgo un arco e il grafo è ancora connesso allora esso farà parte di un ciclo
                flow = remove_arc(temp_network, i, x.first);
                if(!equal_to_enter(i,x.first) && is_connected(temp_network) && flow>0) {
                    cycle.insert ( {i,x.first} );
                    cycle_back[x.first]=i;
                }
                add_arc(temp_network, i, x.first, flow);
            }
        }
        cout << "CYCLE\n";
        for(auto x : cycle){
            cout << x.first << " " << x.second << "\n";
        }
        list<triple> min_flow;
        list<pair<int,int>> right_direction;
        int start = enter.second;
        int count = 0;
        while(start != enter.first){
            if(count==9){
                cout << "\nPerfavore provare a rilanciare il programma, non sempre la scelta casuale di una base della matrice di incidenza si rivela fattibile (devo fissare questo bug)\n";
                return;
            }
            count++;
            if(cycle.find(start)==cycle.end()){
                int part = cycle_back[start];
                triple v(part,start,find_flow(part,start));
                min_flow.push_back(v);
                if(cycle.count(part)>1) {
                    auto range = cycle.equal_range(part);
                    for (auto z = range.first; z != range.second; ++z) {
                        if(z->second==start){
                            cycle.erase (z);
                            break;
                        }
                    }
                }
                start = part;
            }
            else{
                right_direction.emplace_back(make_pair(cycle.find(start)->second,start));
                start = cycle.find(start)->second;
            }
        }
        if(min_flow.empty()){
            unbounded = true;
            return;
        }
        cout << "REVERSE ARCS:\n";
        for(auto x : min_flow){
            cout << x.a << " " << x.b << " " << x.f << "\n";
        }
        double min = min_flow.begin()->f;
        int leaving_arc_a = min_flow.begin()->a;
        int leaving_arc_b = min_flow.begin()->b;
        for(auto iter = min_flow.begin()++;iter!=min_flow.end();iter++){
            if(iter->f < min){
                min = iter->f;
                leaving_arc_a = iter->a;
                leaving_arc_b = iter->b;
            }
        }
        remove_arc(network,leaving_arc_a,leaving_arc_b);
        //aggiungi il flow a right_direction e sottrai il flow a min_flow (tranne chi è stato eliminato)
        for(auto & iter : right_direction){
            for(auto & x : network[iter.first]){
                if(x.first==iter.second){
                    x.second += min;
                }
            }
            for(auto & x : network[iter.second]){
                if(x.first==iter.first){
                    x.second -= min;
                }
            }
        }
        for(auto & iter : min_flow){
            if(is_not_zero(iter.f - min)){
                for(auto & x : network[iter.a]){
                    if(x.first==iter.b){
                        int iii = x.second;
                        x.second -= min;
                    }
                }
                for(auto & x : network[iter.b]){
                    if(x.first==iter.a){
                        x.second += min;
                    }
                }
            }
        }
        for(auto & x : network[enter.first]){
            if(x.first==enter.second){
                x.second = min;
            }
        }
        for(auto & x : network[enter.second]){
            if(x.first==enter.first){
                x.second = -min;
            }
        }

        not_present.erase(make_pair(enter.first,enter.second));
        not_present.insert(make_pair(leaving_arc_a,leaving_arc_b));
        cout << "NOT PRESENT ARCS:\n";
        for (auto i : not_present) {
            cout << i.first << " " << i.second << "\n";
        }
    }


    static void DepthFirstSearch(const vector<list<pair<int,int>>> & graph, int v, vector<bool> & visited){
        visited[v] = true;
        for (auto u : graph[v])
        {
            if (!visited[u.first]) {
                DepthFirstSearch(graph, u.first, visited);
            }
        }
    }

    static bool is_connected(vector<list<pair<int,int>>> & net)
    {
        for (int i = 0; i < net.size(); i++)
        {
            vector<bool> visited(net.size());
            DepthFirstSearch(net, i, visited);
            if (find(visited.begin(), visited.end(), false) != visited.end()) {
                return false;
            }
        }
        return true;
    }

    struct triple{
        int a, b;
        double f;
        triple(int a, int b, int f): a(a), b(b), f(f){}
    };

    void network_simple_algorithm(){
        first_solution();
        not_present_creation();
        complementary_slackness_condition();
        while(!optimality_conditions()){
            remove_leaving_arc();
            complementary_slackness_condition();
        }
        print_solution();
    }

    static bool sorting (int i,int j) { return (i<j); }

    void first_solution(){
        /*
        list<pair<int,int>> temp;
        temp.push_back(make_pair<int,int>(2,10));
        list<pair<int,int>> temp1;
        temp1.push_back(make_pair<int,int>(2,4));
        list<pair<int,int>> temp2;
        temp2.push_back(make_pair<int,int>(0,-10));
        temp2.push_back(make_pair<int,int>(1,-4));
        temp2.push_back(make_pair<int,int>(3,6));
        temp2.push_back(make_pair<int,int>(4,8));
        list<pair<int,int>> temp3;
        temp3.push_back(make_pair<int,int>(2,-6));
        list<pair<int,int>> temp4;
        temp4.push_back(make_pair<int,int>(2,-8));
        network.push_back(temp);
        network.push_back(temp1);
        network.push_back(temp2);
        network.push_back(temp3);
        network.push_back(temp4);
*/

        vector<pair<int,int>> arc_node;
        vector<vector<int>> incidence_matrix;
        vector<int> row(nodes,0);
        pair<int,int> link;
        for(int i=0;i<nodes;i++){
            for(int ii=0;ii<nodes;ii++){
                if(is_not_zero(costs[i][ii])){
                    row[i]=1;
                    row[ii]=-1;
                    link = make_pair(i,ii);
                    arc_node.push_back(link);
                    row.pop_back();
                    incidence_matrix.push_back(row);
                    row.push_back(0);
                    row[i]=0;
                    row[ii]=0;
                }
            }
        }
        const int N = incidence_matrix[0].size();
        bool check = false;
        unordered_set<int> random;
        while(!check){
            vector<vector<double>> inv;
            random.clear();
            int random_number = 0;
            vector<vector<double>> check_matrix_init(N, vector<double> (N, 0));
            vector<vector<double>> check_matrix = check_matrix_init;
            for(int i=0;i<N;i++){
                random_number = find_random(random, incidence_matrix.size());
                for(int ii=0;ii<N;ii++){
                    check_matrix[ii][i] = incidence_matrix[random_number][ii];
                }
            }
            if(getInverse(check_matrix)!=check_matrix_init){
                cout << "RANDOM\n";
                for(auto e : random){
                    cout << e << " ";
                }
                cout << "\n";
                inv = getInverse(check_matrix);
                vector<double> sol;
                for(int i=0;i<N;i++){
                    double x = 0;
                    for(int ii=0;ii<N;ii++){
                        x += inv[i][ii] * b[ii];
                    }
                    sol.push_back(x);
                }
                check = true;
                for(const auto & x : sol){
                    if(x<0){
                        check = false;
                    }
                }
                cout << "check: " << check << "\n";
            }
        }
        vector<int> arcs;
        arcs.reserve(random.size());
        for(auto x : random){
            arcs.push_back(x);
        }
        sort(arcs.begin(), arcs.end(), sorting);
        vector<vector<double>> flow_equations(N, vector<double> (N+1, 0));
        int ii = 0;
        for(int i=0;i<N;i++){
            for(int j = 0;j<N;j++){
                flow_equations[i][ii] = incidence_matrix[arcs[j]][i];
                ii++;
            }
            ii=0;

            flow_equations[i][N] = b[i];
        }
        list<pair<int,int>> emptylist;
        for(int i=0;i<nodes;i++){
            network.push_back(emptylist);
        }

        vector<double> initial_flows = gauss(flow_equations).second;

        for(int i=0;i<initial_flows.size();i++){
            list<pair<int,int>> dai;
            network[arc_node[arcs[i]].first].emplace_back(make_pair(arc_node[arcs[i]].second, initial_flows[i]));
            network[arc_node[arcs[i]].second].emplace_back(make_pair(arc_node[arcs[i]].first, -initial_flows[i]));
        }
    }

    static int find_random(unordered_set<int> & random, const int & N){

        int random_number = rand() % N;
        while(random.find(random_number)!=random.end()){
            random_number = rand() % N;
        }
        random.insert(random_number);
        return random_number;
    }

    static double getDeterminant(const std::vector<std::vector<double>> & vect) {
        int dimension = vect.size();
        if(dimension == 2) {
            return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
        }
        double result = 0;
        int sign = 1;
        for(int i = 0; i < dimension; i++) {
            std::vector<std::vector<double>> subVect(dimension - 1, std::vector<double> (dimension - 1));
            for(int m = 1; m < dimension; m++) {
                int z = 0;
                for(int n = 0; n < dimension; n++) {
                    if(n != i) {
                        subVect[m-1][z] = vect[m][n];
                        z++;
                    }
                }
            }
            result = result + sign * vect[0][i] * getDeterminant(subVect);
            sign = -sign;
        }

        return result;
    }

    static std::vector<std::vector<double>> getTranspose(const std::vector<std::vector<double>>& matrix1) {
        std::vector<std::vector<double>> solution(matrix1[0].size(), std::vector<double> (matrix1.size()));
        for(size_t i = 0; i < matrix1.size(); i++) {
            for(size_t j = 0; j < matrix1[0].size(); j++) {
                solution[j][i] = matrix1[i][j];
            }
        }
        return solution;
    }

    static std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>>& vect) {
        std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));
        std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double> (vect.size() - 1));
        for(std::size_t i = 0; i < vect.size(); i++) {
            for(std::size_t j = 0; j < vect[0].size(); j++) {
                int p = 0;
                for(size_t x = 0; x < vect.size(); x++) {
                    if(x == i) {
                        continue;
                    }
                    int q = 0;
                    for(size_t y = 0; y < vect.size(); y++) {
                        if(y == j) {
                            continue;
                        }
                        subVect[p][q] = vect[x][y];
                        q++;
                    }
                    p++;
                }
                solution[i][j] = pow(-1, i + j) * getDeterminant(subVect);
            }
        }
        return solution;
    }

    static std::vector<std::vector<double>> getInverse(const std::vector<std::vector<double>>& vect) {
        if(getDeterminant(vect) == 0) {
            vector<vector<double>> check_matrix_init(vect.size(), vector<double> (vect.size(), 0));
            return check_matrix_init;
        }
        double d = 1.0/getDeterminant(vect);
        std::vector<std::vector<double>> solution(vect.size(), std::vector<double> (vect.size()));

        for(size_t i = 0; i < vect.size(); i++) {
            for(size_t j = 0; j < vect.size(); j++) {
                solution[i][j] = vect[i][j] * d;
            }
        }
        return getTranspose(getCofactor(solution));
    }

};




int main() {

    cout << "\n---> SIMPLEX METHOD <---\n\n";

    const vector<double> objective = {0,1500,2700,3200,1800,3000,2500};
    const vector<char> signs = {'>','>','>','<','<'};
    vector<vector<double>> table = {         {5000,1,0,0,1,0,0},
                                             {12000,0,1,0,0,1,0},
                                             {8000,0,0,1,0,0,1},
                                             {10000,1,1,1,0,0,0},
                                             {15000,0,0,0,1,1,1}   };

    /*
    vector<double> objective = {0,-5,-7};
    const vector<char> signs = {'<','<','<'};
    vector<vector<double>> table = {   {8,2,1},
                                       {9,1,2},
                                       {5,1,1}   };
    */


    Simplex vrp(objective,signs,table);
    vrp.simplex_algorithm();


    srand (time(NULL));
    cout << "\n---> NETWORK SIMPLEX METHOD <---\n\n";

    /*
    const vector<vector<double>> rete = {   {0,10,8,1,0},
                                            {0,0,2,0,7},
                                            {0,0,0,1,4},
                                            {0,0,0,0,12},
                                            {0,0,0,0,0}   };

    const vector<int> sup_dem{10,4,0,-6,-8};
    */
    const vector<vector<double>> rete = {   {0,0,15,27,32},
                                            {0,0,18,30,25},
                                            {0,0,0,0,0},
                                            {0,0,0,0,0},
                                            {0,0,0,0,0}   };

    const vector<int> sup_dem{100,150,-50,-120,-80};

    Network_Simplex ns(rete ,sup_dem);
    ns.network_simple_algorithm();

    return 0;
}
