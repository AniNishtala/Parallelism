#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <stack>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <chrono>
#include <omp.h>
#include <algorithm>
#include <sstream>
using namespace std;

struct Graph {
    vector<string> cities; // List of city names
    unordered_map<string, vector<pair<string, double>>> adj; // Adjacency list with weights

    void addCity(const string& city) {
        if (find(cities.begin(), cities.end(), city) == cities.end()) {
            cities.push_back(city);
        }
    }

    void addEdge(const string& city1, const string& city2, double weight) {
        adj[city1].push_back(make_pair(city2, weight));
        adj[city2].push_back(make_pair(city1, weight)); // Assuming undirected graph
    }

    vector<string> getCities() const {
        return cities;
    }
};

// Function to clean a string by removing quotes and leading/trailing spaces
string cleanString(const string& input) {
    string result = input;
    result.erase(remove(result.begin(), result.end(), '"'), result.end()); // Remove quotes
    result.erase(result.begin(), find_if(result.begin(), result.end(), [](unsigned char ch) { return !isspace(ch); })); // Trim leading spaces
    result.erase(find_if(result.rbegin(), result.rend(), [](unsigned char ch) { return !isspace(ch); }).base(), result.end()); // Trim trailing spaces
    return result;
}

// Function to read city distances from a file and generate the graph
Graph generateGraphFromFile(const string& filename) {
    Graph graph;
    ifstream file(filename);
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        return graph;
    }

    // Read lines in format: "City1,City2,Weight"
    while (getline(file, line)) {
        if (line.empty()) continue;

        stringstream ss(line);
        string city1, city2, weightStr;
        double weight;
        
        try {
    // Parse CSV line
    getline(ss, city1, ',');
    getline(ss, city2, ',');
    getline(ss, weightStr);

    // Clean strings
    city1 = cleanString(city1);
    city2 = cleanString(city2);
    weightStr = cleanString(weightStr);
    
    // Debugging: Log cleaned values
    //cout << "City1: " << city1 << ", City2: " << city2 << ", Weight: " << weightStr << endl;

    // Remove commas from weight string
    weightStr.erase(remove(weightStr.begin(), weightStr.end(), ','), weightStr.end());
    
    // Convert weight to double
    weight = stod(weightStr);

    // Add cities and edge to the graph
    graph.addCity(city1);
    graph.addCity(city2);
    graph.addEdge(city1, city2, weight);
} catch (const invalid_argument& e) {
    cerr << "Invalid weight value in line: " << line << " (Weight: " << weightStr << ")" << endl;
} catch (const exception& e) {
    cerr << "Error processing line: " << line << " - " << e.what() << endl;
}

        
    }

    file.close();
    return graph;
}



// Sequential BFS
void BFS_sequential(Graph& graph, const string& start) {
    unordered_map<string, bool> visited;
    queue<string> q;
    visited[start] = true;
    q.push(start);

    while (!q.empty()) {
        string city = q.front();
        q.pop();

        for (const auto& neighbor : graph.adj[city]) {
            if (!visited[neighbor.first]) {
                visited[neighbor.first] = true;
                q.push(neighbor.first);
            }
        }
    }
}

// Parallel BFS
void BFS_parallel(Graph& graph, const string& start) {
    unordered_map<string, bool> visited;
    vector<string> frontier = {start};
    visited[start] = true;

    #pragma omp parallel
    {
        vector<string> local_frontier;

        while (!frontier.empty()) {
            local_frontier.clear();

            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < frontier.size(); ++i) {
                const string& city = frontier[i];
                for (const auto& neighbor : graph.adj[city]) {
                    if (!visited[neighbor.first]) {
                        #pragma omp critical
                        {
                            if (!visited[neighbor.first]) {
                                visited[neighbor.first] = true;
                                local_frontier.push_back(neighbor.first);
                            }
                        }
                    }
                }
            }

            #pragma omp single
            {
                frontier.swap(local_frontier);
            }
        }
    }
}


// Sequential DFS
void DFS_sequential(Graph& graph, const string& start) {
    unordered_map<string, bool> visited;
    stack<string> s;
    s.push(start);

    while (!s.empty()) {
        string city = s.top();
        s.pop();

        if (!visited[city]) {
            visited[city] = true;
        }

        for (const auto& neighbor : graph.adj[city]) {
            if (!visited[neighbor.first]) {
                s.push(neighbor.first);
            }
        }
    }
}

// Parallel DFS
void DFS_parallel(Graph& graph, const string& start) {
    unordered_map<string, bool> visited;

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            stack<string> frontier;
            frontier.push(start);

            while (!frontier.empty()) {
                string city = frontier.top();
                frontier.pop();

                #pragma omp critical
                {
                    if (!visited[city]) {
                        visited[city] = true;
                    }
                }

                for (const auto& neighbor : graph.adj[city]) {
                    #pragma omp task firstprivate(neighbor)
                    {
                        if (!visited[neighbor.first]) {
                            frontier.push(neighbor.first);
                        }
                    }
                }
            }
        }
    }
}

// Sequential Dijkstra's algorithm
void Dijkstra_sequential(Graph& graph, const string& start) {
    unordered_map<string, int> dist;
    unordered_map<string, bool> visited;

    // Initialize distances to infinity
    for (size_t i = 0; i < graph.cities.size(); ++i) {
        dist[graph.cities[i]] = numeric_limits<int>::max();
    }
    dist[start] = 0;

    // Priority queue to get the city with the smallest distance
    typedef pair<int, string> Pair;
    priority_queue<Pair, vector<Pair>, greater<Pair>> pq;
    pq.push(make_pair(0, start));

    while (!pq.empty()) {
        string city = pq.top().second;
        pq.pop();

        if (visited[city])
            continue;

        visited[city] = true;

        for (const auto& neighbor : graph.adj[city]) {
            if (dist[city] + neighbor.second < dist[neighbor.first]) {
                dist[neighbor.first] = dist[city] + neighbor.second;
                pq.push(make_pair(dist[neighbor.first], neighbor.first));
            }
        }
    }
}

// Parallel Dijkstra's algorithm
void Dijkstra_parallel(Graph& graph, const string& start) {
    unordered_map<string, int> dist;
    unordered_map<string, bool> visited;

    // Initialize distances to infinity
    for (const auto& city : graph.cities) {
        dist[city] = numeric_limits<int>::max();
    }
    dist[start] = 0;

    typedef pair<int, string> Pair;
    priority_queue<Pair, vector<Pair>, greater<Pair>> pq;
    
    #pragma omp parallel
    {
        // Each thread has its own private priority queue
        priority_queue<Pair, vector<Pair>, greater<Pair>> local_pq;

        #pragma omp single
        {
            local_pq.push(make_pair(0, start));
        }

        while (true) {
            string current_city;
            bool has_work = false;

            // Access the priority queue safely in a critical region
            #pragma omp critical
            {
                if (!local_pq.empty()) {
                    current_city = local_pq.top().second;
                    local_pq.pop();
                    has_work = true;
                }
            }

            if (!has_work) break;

            if (visited[current_city]) continue;

            visited[current_city] = true;

            #pragma omp parallel for
            for (size_t i = 0; i < graph.adj[current_city].size(); ++i) {
                const string& neighbor = graph.adj[current_city][i].first;
                int weight = graph.adj[current_city][i].second;

                if (!visited[neighbor]) {
                    int new_distance = dist[current_city] + weight;

                    // Update the distance if the new distance is smaller
                    #pragma omp critical
                    {
                        if (new_distance < dist[neighbor]) {
                            dist[neighbor] = new_distance;
                            local_pq.push(make_pair(new_distance, neighbor));
                        }
                    }
                }
            }
        }
    }

    /*// Output distances from the starting city
    cout << "\nDistances from " << start << ":\n";
    for (const auto& entry : dist) {
        cout << entry.first << ": "
             << (entry.second == numeric_limits<int>::max() ? "Infinity" : to_string(entry.second)) 
             << endl;
    }*/
}



int main() {
    string filename = "city_distances.txt";
    Graph graph = generateGraphFromFile(filename);

    cout << "Starting graph traversal and shortest path algorithms...\n";

    // --- BFS Sequential ---
    cout << "\nBFS (Sequential) starting from " << graph.cities[0] << "...\n";
    auto start = chrono::high_resolution_clock::now();
    BFS_sequential(graph, graph.cities[0]);
    auto end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    // --- BFS Parallel ---
    cout << "\nBFS (Parallel) starting from " << graph.cities[0] << "...\n";
    start = chrono::high_resolution_clock::now();
    BFS_parallel(graph, graph.cities[0]);
    end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    // --- DFS Sequential ---
    cout << "\nDFS (Sequential) starting from " << graph.cities[0] << "...\n";
    start = chrono::high_resolution_clock::now();
    DFS_sequential(graph, graph.cities[0]);
    end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " microsec" << endl;

    // --- DFS Parallel ---
    cout << "\nDFS (Parallel) starting from " << graph.cities[0] << "...\n";
    start = chrono::high_resolution_clock::now();
    DFS_parallel(graph, graph.cities[0]);
    end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " microsec" << endl;

    // --- Dijkstra Sequential ---
    cout << "\nDijkstra's Algorithm (Sequential) starting from " << graph.cities[0] << "...\n";
    start = chrono::high_resolution_clock::now();
    Dijkstra_sequential(graph, graph.cities[0]);
    end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    // --- Dijkstra Parallel ---
    cout << "\nDijkstra's Algorithm (Parallel) starting from " << graph.cities[0] << "...\n";
    start = chrono::high_resolution_clock::now();
    Dijkstra_parallel(graph, graph.cities[0]);
    end = chrono::high_resolution_clock::now();
    cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    return 0;
}

