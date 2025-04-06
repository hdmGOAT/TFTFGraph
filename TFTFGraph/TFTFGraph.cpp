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
routes[toRoute].edges.push_back({fromRoute, routes[fromRoute].routeName, transferCost, fare, densities}); // Add reverse edge
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
                    std::cout << " (Density Factor : " << densityFactor << ")";
                }
                std::cout << "\n";
            }
        }
    }

    float TFTFGraph::heuristic(int current, int target) {
        if (hopDistance.count(current) && hopDistance[current].count(target)) {
            int hops = hopDistance[current][target];
            return hops * minEdgeCost;
        }
        return 0.0f; // If no path info, fallback to 0 (neutral)
    }
    
    
    void TFTFGraph::precomputeHopDistances() {
        hopDistance.clear();
        minEdgeCost = std::numeric_limits<float>::infinity();
    
        for (const auto& [sourceId, _] : routes) {
            std::queue<int> q;
            std::unordered_map<int, int> visited;
            q.push(sourceId);
            visited[sourceId] = 0;
    
            while (!q.empty()) {
                int current = q.front();
                q.pop();
    
                for (const auto& edge : routes[current].edges) {
                    if (!visited.count(edge.destinationRoute)) {
                        visited[edge.destinationRoute] = visited[current] + 1;
                        q.push(edge.destinationRoute);
                    }
                    // Update the smallest edge cost globally
                    float cost = edge.totalCost();
                    if (cost < minEdgeCost) {
                        minEdgeCost = cost;
                    }
                }
            }
    
            hopDistance[sourceId] = visited;
        }
    
        if (minEdgeCost == std::numeric_limits<float>::infinity()) {
            minEdgeCost = 1.0f; // Fallback to prevent divide-by-zero
        }
    }

    void TFTFGraph::printHopDistances() const {
        std::cout << "\n--- Hop Distances (BFS) ---\n";
        for (const auto& [from, innerMap] : hopDistance) {
            for (const auto& [to, hops] : innerMap) {
                std::cout << "From " << from << " to " << to << ": " << hops << " hops\n";
            }
        }
    }

    


    std::vector<int> TFTFGraph::findBestPath(int startRouteId, int endRouteId, int hour)
    {
        using PQNode = std::pair<float, int>;
        std::priority_queue<PQNode, std::vector<PQNode>, std::greater<PQNode>> forwardPQ, backwardPQ;
    
        std::unordered_map<int, float> forwardCost, backwardCost;
        std::unordered_map<int, int> forwardPrev, backwardPrev;
        std::unordered_set<int> forwardVisited, backwardVisited;
    
        for (const auto &[id, _] : routes)
        {
            forwardCost[id] = std::numeric_limits<float>::infinity();
            backwardCost[id] = std::numeric_limits<float>::infinity();
        }
    
        forwardCost[startRouteId] = 0.0f;
        backwardCost[endRouteId] = 0.0f;
    
        forwardPQ.push({0.0f, startRouteId});
        backwardPQ.push({0.0f, endRouteId});
    
        int meetingPoint = -1;
        float bestPathCost = std::numeric_limits<float>::infinity();
    
        while (!forwardPQ.empty() && !backwardPQ.empty())
        {
            if (!forwardPQ.empty())
            {
                auto [currCost, currRouteId] = forwardPQ.top();
                forwardPQ.pop();
    
                if (forwardVisited.count(currRouteId))
                    continue;
                forwardVisited.insert(currRouteId);
    
                if (backwardVisited.count(currRouteId))
                {
                    meetingPoint = currRouteId;
                    bestPathCost = std::min(bestPathCost, forwardCost[currRouteId] + backwardCost[currRouteId]);
                    break;
                }
    
                for (const auto &edge : routes.at(currRouteId).edges)
                {
                    float edgeCost = edge.totalCost(hour);
                    float newCost = currCost + edgeCost;
    
                    if (newCost < forwardCost[edge.destinationRoute])
                    {
                        forwardCost[edge.destinationRoute] = newCost;
                        forwardPrev[edge.destinationRoute] = currRouteId;
                        forwardPQ.push({newCost + heuristic(edge.destinationRoute, endRouteId), edge.destinationRoute});
                    }
                }
            }
    
            if (!backwardPQ.empty())
            {
                auto [currCost, currRouteId] = backwardPQ.top();
                backwardPQ.pop();
    
                if (backwardVisited.count(currRouteId))
                    continue;
                backwardVisited.insert(currRouteId);
    
                if (forwardVisited.count(currRouteId))
                {
                    meetingPoint = currRouteId;
                    bestPathCost = std::min(bestPathCost, forwardCost[currRouteId] + backwardCost[currRouteId]);
                    break;
                }
    
                for (const auto &edge : routes.at(currRouteId).edges)
                {
                    float edgeCost = edge.totalCost(hour);
                    float newCost = currCost + edgeCost;
    
                    if (newCost < backwardCost[edge.destinationRoute])
                    {
                        backwardCost[edge.destinationRoute] = newCost;
                        backwardPrev[edge.destinationRoute] = currRouteId;
                        backwardPQ.push({newCost + heuristic(edge.destinationRoute, startRouteId), edge.destinationRoute});
                    }
                }
            }
        }
    
        if (meetingPoint == -1)
        {
            std::cout << "No path found between routes " << startRouteId << " and " << endRouteId << std::endl;
            return {};
        }
    
        std::vector<int> path;
        for (int at = meetingPoint; at != startRouteId; at = forwardPrev[at])
        {
            path.push_back(at);
        }
        path.push_back(startRouteId);
        std::reverse(path.begin(), path.end());
    
        if (!path.empty() && path.back() == meetingPoint)
            path.pop_back();
    
        for (int at = meetingPoint; at != endRouteId; at = backwardPrev[at])
        {
            path.push_back(at);
        }
        path.push_back(endRouteId);
    
        return path;
    }
    