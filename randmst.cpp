//#include <stdio.h>
//#include <stdlib.h>
#include <limits>
#include <iostream>
#include <cstdlib>
//#include <vector>
//#include <iterator>
#include <utility>
#include <math.h>
#include <tuple>
#include <ctime>
using namespace std;
/*** by 50943585 and ***

/**
 * 3 ARRAYS
 * Index = Vertex name
 * 1. Array of n pairs of (x,y) float coordinates between 0 and 1
 * 2. Array of shortest distance for each of n vertices
 * 3. Array of n vertices in the tree
 **/

// 2D: Generate array of n coordinate-pairs between 0 and 1
void makegraph(int vertices, tuple<float,float> coords[]){ // n = number of vertices
    for (int i = 0; i < vertices; ++i){
        float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        coords[i] = std::make_pair (x,y);
    }
}

// 3D: Generate array of n coordinate-triples between 0 and 1
void makegraph(int vertices, tuple<float,float,float> coords[]){ // n = number of vertices
    for (int i = 0; i < vertices; ++i){
        float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        coords[i] = std::make_tuple (x, y, z);
    }
}

// 4D: Generate array of n coordinate-triples between 0 and 1
void makegraph(int vertices, tuple<float,float,float,float> coords[]){ // n = number of vertices
    for (int i = 0; i < vertices; ++i){
        float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float w = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        coords[i] = std::make_tuple (x, y, z, w);
    }
}

// 2D: Returns Euclidean distance
float sqDistance(tuple<float,float> a, tuple<float,float> b){
    return sqrt(pow(std::get<0>(a) - std::get<0>(b), 2) + pow(std::get<1>(a) - std::get<1>(b), 2));
}

// 3D: Returns Euclidean distance
float sqDistance(tuple<float,float,float> a, tuple<float,float,float> b){
    return sqrt(pow(std::get<0>(a) - std::get<0>(b), 2) + pow(std::get<1>(a) - std::get<1>(b), 2) + pow(std::get<2>(a) - std::get<2>(b), 2));
}

// 4D: Returns Euclidean distance
float sqDistance(tuple<float,float,float,float> a, tuple<float,float,float,float> b){
    return sqrt(pow(std::get<0>(a) - std::get<0>(b), 2) + pow(std::get<1>(a) - std::get<1>(b), 2) + pow(std::get<2>(a) - std::get<2>(b), 2) + pow(std::get<3>(a) - std::get<3>(b), 2));
}


// Finds and returns vertex with the smallest cost and not yet in MST
int minCost(float costs[], bool visited[], int vertices){
    float shortest = std::numeric_limits<float>::max();
    int index = 0;
    for (int i = 0; i < vertices; ++i){
        if ((!visited[i]) && (costs[i] < shortest)){
            shortest = costs[i];
            index = i;
        }
    }
    return index;
}

int main(int argc, char* argv[]){
    // ./randmst 0 numpoints numtrials dimension
    // argv[0]: mdgraph argv[1]: 0 argv[2]: vertices argv[3]: trials argv[4]: dimension
    
    // Clock
    std::clock_t    start;
    // Seed random number generator
    srand (static_cast <unsigned> (time(0)));
    
    // How many vertices in graph?
    int vertices = atoi(argv[2]);
    int trials = atoi(argv[3]);
    int dimension = atoi(argv[4]);
    
    float trialweights[trials];
    double trialtimes[trials];
    
    if (dimension == 0){
        for (int trial = 0; trial < trials; ++trial){
            start = std::clock();
            float costs[vertices];
            bool visited[vertices];
            float totalcost = 0;
            
            for (int i = 0; i < vertices; ++i){
                visited[i] = false;
                costs[i] = std::numeric_limits<float>::max();
            }
            
            visited[0] = true;
            costs[0] = 0;
            
            for (int count = 0; count < vertices; ++count){
                int u = minCost(costs, visited, vertices);
                visited[u] = true;
                totalcost += costs[u];
                for (int j = 0; j < vertices; ++j){
                    float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
                    if ((!visited[j]) && (x < costs[j]))
                        costs[j] = x;
                }
            }
            trialweights[trial] = totalcost;
        }
    }
    else{
        for (int trial = 0; trial < trials; ++trial){
            start = std::clock();
            tuple <float,float> coords2[vertices];
            tuple <float,float,float> coords3[vertices];
            tuple <float,float,float,float> coords4[vertices];
            
            if (dimension == 2){
                makegraph(vertices, coords2);
                /**
                for (int i = 0; i < vertices; ++i){
                    std::cout << "(" << std::get<0>(coords2[i]) << " " << std::get<1>(coords2[i]) << ")\n";
                }**/
            }
            else if (dimension == 3){
                makegraph(vertices, coords3);
                /**
                for (int i = 0; i < vertices; ++i){
                    std::cout << "(" << std::get<0>(coords3[i]) << " " << std::get<1>(coords3[i])  << " " << std::get<2>(coords3[i]) << ")\n";
                }**/
            }
            else{ // 4D
                makegraph(vertices,coords4);
                /**
                for (int i = 0; i < vertices; ++i){
                    std::cout << "(" << std::get<0>(coords4[i]) << " " << std::get<1>(coords4[i])  << " " << std::get<2>(coords4[i]) << " " << std::get<3>(coords4[i]) << ")\n";
                }**/
            }
    
            float totalCost = 0;
            trialtimes[trial] = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
    
            float costs[vertices];
            for (int i = 0; i < vertices; ++i){
                // Initialize costs to infinity
                costs[i] = std::numeric_limits<float>::max();
            }
            costs[0] = 0; // Starting node has no cost
            
            bool visited[vertices];
            for (int i = 0; i < vertices; ++i){
                visited[i] = false;
            }
            visited[0] = true; // Starting node is added to MST
                        // MST will have V vertices for a fully connected graph
            for (int count = 0; count < vertices; ++count){
                // pick min key vertex
                int vertex = minCost(costs, visited, vertices);
                // add it to our MST
                visited[vertex] = true;
                totalCost += costs[vertex];
                /// << "total cost: " << totalCost << "\n";
                //std::cout << "vertex: " << vertex << " " << costs[vertex] << "\n";
                // update values
                for (int v = 0; v < vertices; ++v){
                    if (dimension == 2){
                        if (!visited[v] && sqDistance(coords2[v], coords2[vertex]) < costs[v])
                            costs[v] = sqDistance(coords2[v], coords2[vertex]);
                    }
                    else if (dimension == 3){
                        if (!visited[v] && sqDistance(coords3[v], coords3[vertex]) < costs[v])
                            costs[v] = sqDistance(coords3[v], coords3[vertex]);                  
                    }
                    else{
                        if (!visited[v] && sqDistance(coords4[v], coords4[vertex]) < costs[v])
                            costs[v] = sqDistance(coords4[v], coords4[vertex]);
                    }
                }
            }
            /**
            for (int i = 0; i < vertices; ++i){
                std::cout << i << "final vertex cost: " <<costs[i] << "\n";
            }**/
            
            //std::cout << totalCost << "\n";
            trialweights[trial] = totalCost;
            trialtimes[trial] = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
        }
    }
    float avg = 0.0;  //or double for higher precision
    double avgtime = 0.0;
    for (int i = 0; i < trials; ++i)
    {
        avg += trialweights[i];
        avgtime += trialtimes[i];
    }
    avg /= trials;
    avgtime /= trials;
    //std::cout << "Average time: " << avgtime << "\n";
    std::cout << avg << " " << vertices << " " << trials << " " << dimension;
    return 0;
}