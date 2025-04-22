#include "astar.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include "json.hpp"

using json = nlohmann::json;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

bool Node::operator<(const Node &other) const
{
    return std::tie(lat, lon) < std::tie(other.lat, other.lon);
}

bool Node::operator==(const Node &other) const
{
    return lat == other.lat && lon == other.lon;
}

bool Node::operator!=(const Node &other) const
{
    return !(*this == other);
}

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

std::vector<Node> astar_geojson(const std::string &filename, Node start, Node goal)
{
    std::ifstream file(filename);
    json geojson;
    file >> geojson;

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

    std::map<Node, double> gScore, fScore;
    std::map<Node, Node> cameFrom;
    std::set<Node> visited;

    auto cmp = [&](const std::pair<double, Node> &a, const std::pair<double, Node> &b)
    {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<double, Node>, std::vector<std::pair<double, Node>>, decltype(cmp)> openSet(cmp);

    gScore[start] = 0;
    fScore[start] = haversine(start, goal);
    openSet.emplace(fScore[start], start);

    while (!openSet.empty())
    {
        Node current = openSet.top().second;
        openSet.pop();

        if (current == goal)
            break;
        if (visited.count(current))
            continue;
        visited.insert(current);

        for (auto &[neighbor, dist] : graph[current])
        {
            double tentative_g = gScore[current] + dist;
            if (!gScore.count(neighbor) || tentative_g < gScore[neighbor])
            {
                cameFrom[neighbor] = current;
                gScore[neighbor] = tentative_g;
                fScore[neighbor] = tentative_g + haversine(neighbor, goal);
                openSet.emplace(fScore[neighbor], neighbor);
            }
        }
    }

    std::vector<Node> path;
    if (!cameFrom.count(goal))
    {
        std::cout << "No path found.\n";
        return path;
    }

    for (Node at = goal; at != start; at = cameFrom[at])
        path.push_back(at);
    path.push_back(start);
    std::reverse(path.begin(), path.end());

    std::cout << "A* Path Distance: " << gScore[goal] / 1000 << " km\n";
    std::cout << "Path:\n";
    for (auto &n : path)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[" << n.lon << ", " << n.lat << "],\n";
    }
    return path;
}
