//
//  main.cpp
//  Graphs
//
//  Created by Vishakh Nair on 11/17/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>

using namespace std;

// Define a high constant for Max
const int MAX = 1000000;
// Entry object to store each node
class Entry{
public:
    double x;
    double y;
    double distance;
    
    Entry(double x1, double y1){
        x = x1;
        y = y1;
    }
    
    Entry(double x1, double y1, double dist){
        x = x1;
        y = y1;
        distance = dist;
    }
    
    Entry(){
        
    }
};


// Array class substitutes for two dimensional array
class Array{
public:
    // Single dimensional array initialized with n^2 spaces so that array[a][b] is array[a * size + b]
    double *array;
    int size;
    
    // Constructor to initialize n^2 spaces to array
    Array(int n){
        size = n;
        array = new double[n * n];
        for (int i = 0; i < n * n; i++){
            array[i] = 0;
        }
    }
    
    // Function to add node to matrix
    void add(int x, int y, double distance){
        array[x * size + y] = distance;
    }
};

// Main graph class
class Graph{
public:
    // Keeps track of all nodes added
    Entry *list;
    // Spanning Tree
    Entry *MSTarray;
    // Reference to the two dimenional array
    Array *matrix;
    // distance and previous and parent used for shortest path algorithm
    vector<int> distance;
    vector<int> previous;
    int *parent;
    // Keeps trakc of number of nodes added
    int count;
    // Radio distance
    double R;
    
    // Constructor
    Graph(int n, double num){
        list = new Entry[n];
        matrix = new Array(n);
        MSTarray = new Entry[n - 1];
        parent = new int[n];
        count = 0;
        R = num;
    }
    
    // Function to find the distance between 2 points
    double findDistance(double x1, double x2, double y1, double y2){
        return sqrt(pow(x2 - x1, 2) + pow(y2 -y1, 2));
    }
    
    // Find and union find are part of Kruskal Algorithm
    int find(int x){
        while (parent[x] != x){
            x = parent[x];
        }
        return x;
    }
    
    void unionFind(int i, int j){
        int a = find(i);
        int b = find(j);
        parent[a] = b;
    }
    
    void kruskal(Array *array){
        // Initalize min to find minimum weight
        int min = 0;
        // Initialize parent array
        for (int i = 0; i < array->size; i++){
            parent[i] = i;
        }
        
        // Keep track of edges added to the tree
        int countEdges = 0;
        while(countEdges < array->size - 1){
            // Set minimum to keep track of minimum weight
            double minimum = MAX;
            // Placeholders to find correct i and j value
            int a = -1;
            int b = -1;
            for (int i = 0; i < array->size; i++){
                for (int j = 0; j < array->size; j++){
                    // Have to make sure we are not adding the same node
                    if (find(i) != find(j)){
                        // Only proceed if the path has minimum weight
                        if(array->array[i * array->size + j] < minimum){
                            // Since we initialized to 0 at the beggining, have to make sure they also aren't 0
                            if(array->array[i * array->size + j] != 0){
                                minimum = array->array[i * array->size + j];
                                a = i;
                                b = j;
                            }
                        }
                    }
                }
            }
            // Call union find algorithm from before
            unionFind(a, b);
            // Add edge to the array to later print in correct order
            MSTarray[countEdges] = Entry(a+1, b+1, findDistance(list[a].x, list[b].x, list[a].y, list[b].y));
            // Incrementation
            countEdges++;
            min = min + minimum;
        }
    }
    
    // Function to print the tree in the correct order
    void printMST(){
        // Sort the array by distances
        sort(MSTarray, MSTarray + matrix->size - 1, [](Entry& lhs, Entry& rhs) {
              return lhs.distance < rhs.distance;
        });
        
        // Find total sum to print
        double sum = 0;
        for (int i = 0; i < matrix->size - 1; i++){
            // Make sure we have 2 decimal places
            cout.precision(3);
            cout << MSTarray[i].x << " " << MSTarray[i].y << " " << MSTarray[i].distance << endl;
            sum = sum + MSTarray[i].distance;
        }
        // 2 decimal places
        cout.precision(3);
        cout << sum << endl;
        
    }
    
    // Actually adds the node by calling add() function from before
    void addNode(double x, double y){
        
        list[count] = Entry(x, y);
        
        for(int i = 0; i < count; i++){
            double comp = findDistance(list[i].x, x, list[i].y, y);
            if(comp < R){
                matrix->add(i, count, comp);
                matrix->add(count, i, comp);
            }
        }
        count++;
    }
    
    // Modifies the adjacency matrix so that we have a 0s and 1s
    void modify(){
        for (int i = 0; i < matrix->size * matrix->size; i++){
            if (matrix->array[i] != 0){
                matrix->array[i] = 1;
            }
        }
    }
    
    // Gets all the neighbors of a given node
    vector<int> getNeighbors(int node) {
        vector<int> vec;
        // Check a row of matrix and see which nodes are 1
        for (int i = 0; i < matrix->size; i++){
            if (matrix->array[node * matrix->size + i] == 1) {
                // Add nodes to vector
                vec.push_back(i);
            }
        }
        return vec;
    }
    
    // Dijkstra's algorithm for finding shortest path
    void shortestPath(int source){
        // Reset distance and previous vectors
        distance.clear();
        previous.clear();
        // Initialize distance and previous vectors
        for (int i = 0; i < matrix->size; i++){
            distance.push_back(MAX);
            previous.push_back(-1);
        }
        
        
        set<int> vertices;
        for (int i = 0; i < matrix->size; i++){
            vertices.insert(i);
        }

        distance[source] = 0;
        while (!vertices.empty()){
            int min = MAX;
            set<int>::iterator minVertex;
            for (set<int>::iterator  it = vertices.begin(); it != vertices.end(); ++it){
                if (distance[*it] < min) {
                    min = distance[*it];
                    minVertex = it;
                }
            }
            
            // If min hasn't changed, break
            if (min == INT_MAX) {
                break;
            }
            
            // Erase vertex that has just been used
            vertices.erase(minVertex);

            int currNode = *minVertex;
            vector<int> neighbors = getNeighbors(currNode);
            
            // Add 1 to current distance because weight is now 1
            for (int i = 0; i < neighbors.size(); ++i) {
                int neighbor = neighbors[i];
                int newDist = min + 1;
                if (newDist < distance[neighbor]) {
                    distance[neighbor] = newDist;
                    previous[neighbor] = currNode;
                }
            }
        }
    }
    
    // Go through every node to find longest shortest path
    int diameter() {
        // Max number
        int m = 0;
        // This for loop checks for shortest path from every node
        for (int i = 0; i < matrix->size; i++) {
            shortestPath(i);
            // This for loops checks every index of the vector for the max value
            for (int j = 0; j < matrix->size; j++) {
                if (distance[j] > m) {
                    m = distance[j];
                }
            }
            
        }
        return m;
    }
    
    // Prints the path
    void printPath(){
        shortestPath(0);
        for (int i = 1; i < matrix->size; i++){
            string path = "";
            // Start with j at i, then every time we add a new number, we go to the shortest path from j
            int j = i;
            while (j != -1){
                // Need to add 1 because index starts at 0
                path = to_string(j + 1) + " " + path;
                j = previous[j];
            }
            cout << path << distance[i] << endl;
        }
    }
    
    // Function to find number of colored vertices
    int numColors(){
        int max = 0;
        
        for (int i = 0; i < matrix->size; i++){
            int num = 0;
            
            for (int j = 0; j < matrix->size; j++){
                if (matrix->array[i * matrix->size + j] == 1){
                    num++;
                }
            }
            
            if (num > max){
                max = num;
            }
            
        }
        return max;
    }
    
};


// Main function
int main(int argc, const char * argv[]) {
    // Open file
    // This section is just to get the R value before initalizing graph and matrix
    ifstream inFile("GraphData.txt");
    string str;
    getline(inFile, str);
    int num = stoi(str);
    
    for (int i = 0; i < num; i++){
        getline(inFile, str);
    }
    getline(inFile, str);
    Graph graph(num, stod(str));
    
    // This section actually puts everything in the matrix

    ifstream inFile2("GraphData.txt");
    getline(inFile2, str);
    
    for (int i = 0; i < num; i++){
        getline(inFile2, str);
        graph.addNode(stod(str.substr(0, str.find(" "))), stod(str.substr(str.find(" ") + 1)));
    }
    
    // Call all of the above function
    graph.kruskal(graph.matrix);
    graph.printMST();
    graph.modify();
    graph.printPath();
    cout << graph.diameter() << endl;
    cout << graph.numColors() << endl;
    
    return 0;
}
