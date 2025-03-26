#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <limits>
#include "TFTFGraph.h"

// TFTF Edge Structure
float TFTFEdge::totalCost(int hour = -1) const
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


// TFTF Graph Class
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


void TFTFGraph::visualize(int hour = -1) const
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

                    // Find the density factor for display
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

std::vector<int> TFTFGraph::findBestPath(int startRoute, int endRoute, int hour = -1)
    {
        // Dijkstra's algorithm implementation would go here
        // For this example, we'll just show a placeholder
        std::cout << "\nFinding best path from Route " << startRoute
                  << " to Route " << endRoute;
        if (hour >= 0)
        {
            std::cout << " at hour " << hour;
        }
        std::cout << "\n(Implementation of pathfinding algorithm would go here)\n";
    }
