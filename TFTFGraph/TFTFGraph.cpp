#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <limits>
#include "TFTFGraph.h"
#include <queue>
#include <unordered_set>
#include <limits>
#include <algorithm>

float TFTFEdge::totalCost(int hour) const
    {
        float densityFactor = 1.0f;
        if (hour >= 0)
        {
            for (const auto &interval : densityByInterval)
            {
                if (interval.contains(hour))
                {
                    densityFactor = interval.jeepneyDensity;
                    break;
                }
            }
        }
    return (transferCost + fare) * densityFactor;
}

void TFTFGraph::addRoute(int routeId, const std::string &routeName)
    {
        routes[routeId] = {routeId, routeName, {}};
    }
void TFTFGraph::addEdge(int fromRoute, int toRoute, const std::string &toName,
                 float transferCost, float fare,
                 const std::vector<JeepneyDensity> &densities)
    {
        routes[fromRoute].edges.push_back({toRoute, toName, transferCost, fare, densities});
    }


void TFTFGraph::visualize(int hour) const
    {
        std::cout << "\nTFTF Graph Visualization";
        if (hour >= 0)
        {
            std::cout << " (Hour: " << hour << ")";
        }
        std::cout << "\n===========================\n";

        for (const auto &pair : routes)
        {
            const RouteNode &node = pair.second;
            std::cout << "Route " << node.routeId << ": " << node.routeName << "\n";

            if (node.edges.empty())
            {
                std::cout << "  (No outgoing transfers)\n";
                continue;
            }

            for (const auto &edge : node.edges)
            {
                std::cout << "  -> Route " << edge.destinationRoute << ": "
                          << edge.destinationRouteName << "\n";
                std::cout << "     Transfer Cost: " << edge.transferCost
                          << ", Fare: " << edge.fare;

                if (hour >= 0)
                {
                    float cost = edge.totalCost(hour);
                    std::cout << ", Total Cost: " << std::fixed << std::setprecision(2) << cost;

                    float densityFactor = 1.0f;
                    for (const auto &interval : edge.densityByInterval)
                    {
                        if (interval.contains(hour))
                        {
                            densityFactor = interval.jeepneyDensity;
                            break;
                        }
                    }
                    std::cout << " (Density Factor: " << densityFactor << ")";
                }
                std::cout << "\n";
            }
        }
    }

std::vector<int> TFTFGraph::findBestPath(int startRouteId, int endRouteId, int hour) {
    using PQNode = std::pair<float, int>;
    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<PQNode>> pq;

    std::unordered_map<int, float> costToReach;
    std::unordered_map<int, int> previous;
    std::unordered_set<int> visited;

    for (const auto& [id, _] : routes) {
        costToReach[id] = std::numeric_limits<float>::infinity();
    }

    costToReach[startRouteId] = 0.0f;
    pq.push({0.0f, startRouteId});

    while (!pq.empty()) {
        auto [currCost, currRouteId] = pq.top();
        pq.pop();

        if (visited.count(currRouteId)) continue;
        visited.insert(currRouteId);

        if (currRouteId == endRouteId) break;

        const auto& currentNode = routes.at(currRouteId);
        for (const auto& edge : currentNode.edges) {
            float edgeCost = edge.totalCost(hour);
            float newCost = currCost + edgeCost;

            if (newCost < costToReach[edge.destinationRoute]) {
                costToReach[edge.destinationRoute] = newCost;
                previous[edge.destinationRoute] = currRouteId;
                pq.push({newCost, edge.destinationRoute});
            }
        }
    }

    std::vector<int> path;
    if (costToReach[endRouteId] == std::numeric_limits<float>::infinity()) {
        return {};
    }

    for (int at = endRouteId; at != startRouteId; at = previous[at]) {
        path.push_back(at);
    }
    path.push_back(startRouteId);
    std::reverse(path.begin(), path.end());

    return path;
}
