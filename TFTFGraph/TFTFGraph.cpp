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
#include "Helpers/helpers.h"

std::vector<JeepneyDensity> averageRouteDensities(const std::vector<JeepneyDensity>& a,
    const std::vector<JeepneyDensity>& b) {
    std::vector<JeepneyDensity> result;

    std::unordered_map<int, float> mapA, mapB;

    for (const auto& d : a) {
        for (int h = d.startHour; h < d.endHour; ++h) {
            mapA[h] = d.jeepneyDensity;
        }
    }

    for (const auto& d : b) {
        for (int h = d.startHour; h < d.endHour; ++h) {
            mapB[h] = d.jeepneyDensity;
        }
    }

    // Calculate average per hour
    for (int h = 0; h < 24; ++h) {
        float dA = mapA.count(h) ? mapA[h] : 1.0f; // Default to 1.0 if missing
        float dB = mapB.count(h) ? mapB[h] : 1.0f;
        float avg = (dA + dB) / 2.0f;

        result.push_back({h, h + 1, avg});
    }

    return result;
}


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
    return (transferCost) * densityFactor;
}

void TFTFGraph::addRoute(int routeId, const std::string &routeName){
        routes[routeId] = {routeId, routeName, {}};
}

void TFTFGraph::setRoutePath(int routeId, const std::vector<Coordinate>& coordinates) {
    if (routes.find(routeId) != routes.end()) {
        routes[routeId].path = coordinates;
    } else {
        std::cerr << "Route ID " << routeId << " not found.\n";
    }
}

void TFTFGraph::addEdge(int fromRoute, int toRoute, const std::string &toName,
        float transferCost,
        const std::vector<JeepneyDensity> &densities, Coordinate entryCoord, Coordinate exitCoord) 
{
    std::vector<JeepneyDensity> avgDensity = averageRouteDensities(routes[fromRoute].densities,
                                                                    routes[toRoute].densities);

    // Add edge from -> to
    routes[fromRoute].edges.push_back({toRoute, toName, transferCost, avgDensity});
    routes[fromRoute].edges.back().entryCoord = entryCoord;
    routes[fromRoute].edges.back().exitCoord = exitCoord;
}


void TFTFGraph::setRouteDensities(int routeId, const std::vector<JeepneyDensity>& densities) {
    if (routes.find(routeId) != routes.end()) {
        routes[routeId].densities = densities;
    } else {
        std::cerr << "Route ID " << routeId << " not found.\n";
    }
}


void TFTFGraph::createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer) {
    for (const auto& [fromId, fromNode] : routes) {
        for (const auto& [toId, toNode] : routes) {
            if (fromId == toId) continue;

            for (const auto& fromCoord : fromNode.path) {
                for (const auto& toCoord : toNode.path) {
                    float dist = haversine(fromCoord, toCoord);
                    if (dist <= transferRangeMeters) {
                        bool exists = false;
                        for (const auto& edge : fromNode.edges) {
                            if (edge.destinationRoute == toId && std::abs(edge.transferCost - dist) < 1e-2f) {
                                exists = true;
                                break;
                            }
                        }
                        if (!exists) {
                            float fare = computeSegmentFare(routes[fromId].path, fromCoord, toCoord);

                            addEdge(fromId, toId, toNode.routeName, dist,  {}, fromCoord, toCoord);

                        }
                    }
                }
            }
        }
    }
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
                std::cout << "     Transfer Cost: " << edge.transferCost;

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


std::vector<TFTFEdge> TFTFGraph::findBestPath(int startRouteId, int endRouteId, int hour) {
    using PQNode = std::pair<float, int>;
    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<PQNode>> dijkstraPQ;
    
    std::unordered_map<int, float> dijkstraCost;
    std::unordered_map<int, TFTFEdge> dijkstraPrevEdge;      // edge used to reach this node
    std::unordered_map<int, int> prevRouteMap;               // previous node ID
    std::unordered_set<int> dijkstraVisited;
    
    for (const auto &[id, _] : routes) {
        dijkstraCost[id] = std::numeric_limits<float>::infinity();
    }
    
    dijkstraCost[startRouteId] = 0.0f;
    dijkstraPQ.push({0.0f, startRouteId});
    
    while (!dijkstraPQ.empty()) {
        auto [currCost, currRouteId] = dijkstraPQ.top();
        dijkstraPQ.pop();
        
        if (dijkstraVisited.count(currRouteId))
            continue;
        
        dijkstraVisited.insert(currRouteId);
        
        if (currRouteId == endRouteId)
            break;
        
        for (const auto &edge : routes.at(currRouteId).edges) {
            float edgeCost = edge.totalCost(hour);
            float newCost = currCost + edgeCost;
            
            if (newCost < dijkstraCost[edge.destinationRoute]) {
                dijkstraCost[edge.destinationRoute] = newCost;
                dijkstraPrevEdge[edge.destinationRoute] = edge;
                prevRouteMap[edge.destinationRoute] = currRouteId;
                dijkstraPQ.push({newCost, edge.destinationRoute});
            }
        }
    }

    // Reconstruct path
    std::vector<TFTFEdge> pathEdges;
    int currentRouteId = endRouteId;

    while (currentRouteId != startRouteId) {
        if (dijkstraPrevEdge.find(currentRouteId) == dijkstraPrevEdge.end()) {
            std::cout << "No valid path found.\n";
            return {};
        }

        pathEdges.push_back(dijkstraPrevEdge[currentRouteId]);
        currentRouteId = prevRouteMap[currentRouteId];
    }

    std::reverse(pathEdges.begin(), pathEdges.end());

    
    return pathEdges;
}


TFTFEdge* TFTFGraph::getEdge(int fromRoute, int toRoute) const {
    // Check if the 'fromRoute' exists in the routes map
    if (routes.find(fromRoute) != routes.end()) {
        // Look through all the edges of the 'fromRoute' route
        for (const auto& edge : routes.at(fromRoute).edges) {
            if (edge.destinationRoute == toRoute) {
                // Return the pointer to the matching edge
                return const_cast<TFTFEdge*>(&edge);  // We're modifying the returned object, hence the cast
            }
        }
    }
    return nullptr;  // Return nullptr if no matching edge is found
}
    
int TFTFGraph::findClosestRoute(const Coordinate& startCoord) {
    float closestDistance = std::numeric_limits<float>::infinity();
    int closestRouteId = -1;

    for (const auto& [routeId, routeNode] : routes) {
        // Assume that we are using the first coordinate in the route as the closest point
        float dist = haversine(routeNode.path.front(), startCoord); // Compare with the first coordinate in path
        if (dist < closestDistance) {
            closestDistance = dist;
            closestRouteId = routeId;
        }
    }

    return closestRouteId;
}

// Function to calculate route based on coordinates
std::vector<TFTFEdge> TFTFGraph::calculateRouteFromCoordinates(const Coordinate& startCoord, const Coordinate& endCoord, int hour) {
    int startRouteId = findClosestRoute(startCoord);
    int endRouteId = findClosestRoute(endCoord);

    if (startRouteId == -1 || endRouteId == -1) {
        std::cerr << "Could not find a route for the given coordinates." << std::endl;
        return {};
    }

    std::cout << "Closest route from start: " << startRouteId << "\n";
    std::cout << "Closest route from end: " << endRouteId << "\n";

    std::vector<TFTFEdge> path = findBestPath(startRouteId, endRouteId, hour);

    // Inject dummy start edge
    if (!path.empty()) {
        TFTFEdge startEdge;
        startEdge.destinationRoute = path.front().destinationRoute;  // first hop's destination
        startEdge.destinationRouteName = routes[startRouteId].routeName;
        startEdge.transferCost = 0.0f;
        startEdge.densityByInterval = {};
        startEdge.entryCoord = startCoord;
        startEdge.exitCoord = path.front().entryCoord;
        path.insert(path.begin(), startEdge);
    }

    // Inject dummy end edge
    if (!path.empty()) {
        TFTFEdge endEdge;
        endEdge.destinationRoute = -1;  // -1 indicates no further route
        endEdge.destinationRouteName = "Destination";
        endEdge.transferCost = 0.0f;
        endEdge.densityByInterval = {};
        endEdge.entryCoord = path.back().exitCoord;
        endEdge.exitCoord = endCoord;
        path.push_back(endEdge);
    }

    // Print path
    std::cout << "Best Path: ";
    for (const auto& edge : path) {
        std::cout << " -> " << edge.destinationRouteName;
    }
    std::cout << std::endl;

    // Fare (skip dummy start/end edges)
    float totalFare = 0.0f;
    calculateTotalFare(path, startCoord, endCoord);

    std::cout << "Total fare: " << std::fixed << std::setprecision(2) << totalFare << " pesos" << std::endl;

    return path;
}
