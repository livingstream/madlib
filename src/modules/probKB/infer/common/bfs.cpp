#include "bfs.h"
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
    prob.resize(V, 0);
}
 
void Graph::addEdge(int v, int w)
{
    if(v > w) {
       int temp = w;
       w = v;
       v = temp;
    }
    if(s.count(make_pair(v,w)) == 0) { 
      adj[v].push_back(w); // Add w to v’s list.
      adj[w].push_back(v); // Add w to v’s list.
      s.insert(make_pair(v,w));
    } 
}
 
vector<double> Graph::BFS(int s)
{
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
        visited[i] = false;
 
    // Create a queue for BFS
    list<int> queue;
    list<int> hops;
 
    // Mark the current node as visited and enqueue it
    visited[s] = true;
    queue.push_back(s);
    hops.push_back(1);
 
    // 'i' will be used to get all adjacent vertices of a vertex
    list<int>::iterator i;
    int hop = 0; 
    while(!queue.empty())
    {
        // Dequeue a vertex from queue and print it
        s = queue.front();
        hop = hops.front();
        queue.pop_front();
        hops.pop_front();
        prob[s] = 1.0/ ceil(sqrt(hop)); 
        // Get all adjacent vertices of the dequeued vertex s
        // If a adjacent has not been visited, then mark it visited
        // and enqueue it
        for(i = adj[s].begin(); i != adj[s].end(); ++i)
        {
            if(!visited[*i])
            {
                visited[*i] = true;
                queue.push_back(*i);
                hops.push_back(hop + 1);
            }
        }
    }
    return prob;
}
