// traverses vertices reachable from s.
#ifndef H_BFS_H
#define H_BFS_H
#include <iostream>
#include <list>
#include <vector>
#include <math.h>     
#include <set>
#include <utility>
 
using namespace std;
 
// This class represents a directed graph using adjacency list representation
class Graph
{
    set< pair<int, int> > s;
    int V;    // No. of vertices
    list<int> *adj;    // Pointer to an array containing adjacency lists
    vector<double> prob;
public:
    Graph(int V);  // Constructor
    void addEdge(int v, int w); // function to add an edge to graph
    vector<double> BFS(int s);  // prints BFS traversal from a given source s
};
#endif
