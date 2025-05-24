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

std::vector<Node> dijkstra_geojson(const std::string &filename, Node origin, Node destination, std::map<Node, std::vector<std::pair<Node, double>>> graph)
{
    auto start = std::chrono::high_resolution_clock::now(); 
    std::map<Node, double> dist;
    std::map<Node, Node> prev;
    std::set<Node> visited;

    auto cmp = [](const std::pair<double, Node> &a, const std::pair<double, Node> &b)
    {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<double, Node>, std::vector<std::pair<double, Node>>, decltype(cmp)> pq(cmp);

    // Initialize with transfer points near origin
    for (const auto& [node, edges] : graph) {
        if (node.isTransferPoint) continue; // Skip other transfer points
        
        double distToStart = haversineNode(node, origin);
        if (distToStart <= 300.5) { // Same transfer range as TFTFGraph
            dist[node] = distToStart;
            prev[node] = origin;
            pq.emplace(distToStart, node);
        }
    }

    while (!pq.empty())
    {
        auto [d, u] = pq.top();
        pq.pop();

        if (visited.count(u))
            continue;
        visited.insert(u);

        // Check if we're close enough to the destination
        if (haversineNode(u, destination) <= 300.5) {
            // Found a path to destination
            std::vector<Node> path;
            for (Node at = u; at != origin; at = prev[at])
                path.push_back(at);
            path.push_back(origin);
            std::reverse(path.begin(), path.end());
            path.push_back(destination);
            return path;
        }

        for (auto &[v, w] : graph[u])
        {
            // Skip invalid transfers
            if (u.routeId != -1 && v.routeId != -1 && 
                u.routeId != v.routeId && 
                !v.isTransferPoint) {
                continue;
            }

            if (!dist.count(v) || dist[v] > d + w)
            {
                dist[v] = d + w;
                prev[v] = u;
                pq.emplace(dist[v], v);
            }
        }
    }

    return std::vector<Node>();
}
