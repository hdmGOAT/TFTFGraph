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

std::vector<Node> astar_geojson(const std::string &filename, Node start, Node goal, std::map<Node, std::vector<std::pair<Node, double>>> graph)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // A* algorithm
    std::map<Node, double> gScore, fScore;
    std::map<Node, Node> cameFrom;
    std::set<Node> visited;

    auto cmp = [&](const std::pair<double, Node> &a, const std::pair<double, Node> &b)
    {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<double, Node>, std::vector<std::pair<double, Node>>, decltype(cmp)> openSet(cmp);

    // Initialize with transfer points near start
    for (const auto& [node, edges] : graph) {
        if (node.isTransferPoint) continue; // Skip other transfer points
        
        double distToStart = haversineNode(node, start);
        if (distToStart <= 300.5) { // Same transfer range as TFTFGraph
            gScore[node] = distToStart;
            fScore[node] = distToStart + haversineNode(node, goal);
            openSet.emplace(fScore[node], node);
            cameFrom[node] = start;
        }
    }

    while (!openSet.empty())
    {
        Node current = openSet.top().second;
        openSet.pop();

        // Check if we're close enough to the goal
        if (haversineNode(current, goal) <= 300.5) {
            // Found a path to goal
            std::vector<Node> path;
            for (Node at = current; at != start; at = cameFrom[at])
                path.push_back(at);
            path.push_back(start);
            std::reverse(path.begin(), path.end());
            path.push_back(goal);
            return path;
        }

        if (visited.count(current))
            continue;
        visited.insert(current);

        for (auto &[neighbor, dist] : graph[current])
        {
            // Skip invalid transfers
            if (current.routeId != -1 && neighbor.routeId != -1 && 
                current.routeId != neighbor.routeId && 
                !neighbor.isTransferPoint) {
                continue;
            }

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
    return path;
}
