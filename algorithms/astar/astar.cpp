#include "astar.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <chrono> // <-- added this
#include "../json.hpp"
#include "../djikstra/djikstra.h"
#include "../../TFTFGraph/Helpers/helpers.h"

using json = nlohmann::json;

std::vector<Node> astar_geojson(const std::string &filename, Node start, Node goal)
{
    auto start_time = std::chrono::high_resolution_clock::now(); // start timing

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
            double dist = haversineNode(u, v);
            graph[u].emplace_back(v, dist);
            graph[v].emplace_back(u, dist);
        }
    }
    printGraphDetails(graph);

    // A* algorithm
    std::map<Node, double> gScore, fScore;
    std::map<Node, Node> cameFrom;
    std::set<Node> visited;

    auto cmp = [&](const std::pair<double, Node> &a, const std::pair<double, Node> &b)
    {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<double, Node>, std::vector<std::pair<double, Node>>, decltype(cmp)> openSet(cmp);

    gScore[start] = 0;
    fScore[start] = haversineNode(start, goal);
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
                fScore[neighbor] = tentative_g + haversineNode(neighbor, goal);
                openSet.emplace(fScore[neighbor], neighbor);
            }
        }
    }

    std::vector<Node> path;
    if (!cameFrom.count(goal))
    {
        std::cout << "No path found.\n";

        auto end_time = std::chrono::high_resolution_clock::now(); // end timing
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "A* completed in: " << duration << " ms\n";

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

    auto end_time = std::chrono::high_resolution_clock::now(); // end timing
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "A* completed in: " << duration << " ms\n";

    return path;
}
