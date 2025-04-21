#include "djikstra.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <queue>
#include <map>
#include <set>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using json = nlohmann::json;

double haversine(const Node &a, const Node &b)
{
    const double R = 6371e3;
    double lat1 = a.lat * M_PI / 180;
    double lat2 = b.lat * M_PI / 180;
    double dlat = (b.lat - a.lat) * M_PI / 180;
    double dlon = (b.lon - a.lon) * M_PI / 180;

    double h = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);

    return R * 2 * atan2(sqrt(h), sqrt(1 - h));
}

std::vector<Node> dijkstra_geojson(const std::string &filename, Node origin, Node destination)
{
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
            double dist = haversine(u, v);

            graph[u].emplace_back(v, dist);
            graph[v].emplace_back(u, dist);
        }
    }

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
        return path;
    }

    for (Node at = destination; at != origin; at = prev[at])
        path.push_back(at);
    path.push_back(origin);
    std::reverse(path.begin(), path.end());

    std::cout << "Shortest distance: " << dist[destination] / 1000 << " km\n";
    // std::cout << "Path:\n";
    // for (auto &n : path)
    //     std::cout << "(" << n.lat << ", " << n.lon << ")\n";

    return path;
}
