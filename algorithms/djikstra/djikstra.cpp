#include "djikstra.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <chrono>  // <-- added
#include "../../TFTFGraph/Helpers/helpers.h"

using json = nlohmann::json;

std::vector<Node> dijkstra_geojson(const std::string &filename, Node origin, Node destination)
{
    auto start = std::chrono::high_resolution_clock::now(); // start timing

    // Parse GeoJSON
    std::ifstream file(filename);
    json geojson;
    file >> geojson;

    // Build graph
    std::map<Node, std::vector<std::pair<Node, double>>> graph;

    for (auto &feature : geojson["features"])
    {
        if (feature["geometry"]["type"] != "LineString")
            continue;
        auto coords = feature["geometry"]["coordinates"];

        for (size_t i = 0; i + 1 < coords.size(); ++i)
        {
            Node u{coords[i][1], coords[i][0]};
            Node v{coords[i + 1][1], coords[i + 1][0]};
            double dist = haversineNode(u, v);

            graph[u].emplace_back(v, dist);
            graph[v].emplace_back(u, dist);
        }
    }

    printGraphDetails(graph);

    // Dijkstra’s algorithm
    std::map<Node, double> dist;
    std::map<Node, Node> prev;
    std::set<Node> visited;

    auto cmp = [](const std::pair<double, Node> &a, const std::pair<double, Node> &b)
    {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<double, Node>, std::vector<std::pair<double, Node>>, decltype(cmp)> pq(cmp);

    dist[origin] = 0;
    pq.emplace(0, origin);

    while (!pq.empty())
    {
        auto [d, u] = pq.top();
        pq.pop();
        if (visited.count(u))
            continue;
        visited.insert(u);

        if (u == destination)
            break;

        for (auto &[v, w] : graph[u])
        {
            if (!dist.count(v) || dist[v] > d + w)
            {
                dist[v] = d + w;
                prev[v] = u;
                pq.emplace(dist[v], v);
            }
        }
    }

    std::vector<Node> path;
    if (!dist.count(destination))
    {
        std::cout << "No path found.\n";

        auto end = std::chrono::high_resolution_clock::now(); // end timing
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Dijkstra completed in: " << duration << " ms\n";

        return path;
    }

    for (Node at = destination; at != origin; at = prev[at])
        path.push_back(at);
    path.push_back(origin);
    std::reverse(path.begin(), path.end());

    std::cout << "Shortest distance: " << dist[destination] / 1000 << " km\n";
    std::cout << "Path:\n";
    for (auto &n : path)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[" << n.lon << ", " << n.lat << "]" << "," << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now(); // end timing
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Dijkstra completed in: " << duration << " ms\n";

    return path;
}
