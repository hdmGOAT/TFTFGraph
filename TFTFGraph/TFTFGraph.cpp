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

bool operator==(const Coordinate& lhs, const Coordinate& rhs) {
    return lhs.latitude == rhs.latitude && lhs.longitude == rhs.longitude;
}

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

    for (int h = 0; h < 24; ++h) {
        float dA = mapA.count(h) ? mapA[h] : 1.0f; 
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
    routes[fromRoute].edges.push_back({toRoute, fromRoute,toName, transferCost, avgDensity});
    routes[fromRoute].edges.back().entryCoord = entryCoord;
    routes[fromRoute].edges.back().exitCoord = exitCoord;
}
std::vector<const RouteNode*> TFTFGraph::extractTraversedRouteNodes(const std::vector<TFTFEdge>& path) const {
    std::vector<const RouteNode*> routeNodes;

    if (path.empty()) return routeNodes;

    if (routes.count(path.front().originRoute)) {
        routeNodes.push_back(&routes.at(path.front().originRoute));
    }

    for (const auto& edge : path) {
        if (routes.count(edge.destinationRoute)) {
            routeNodes.push_back(&routes.at(edge.destinationRoute));
        }
    }

    return routeNodes;
}

double TFTFGraph::calculateTotalFare(const std::vector<TFTFEdge>& path, const Coordinate& startCoord, const Coordinate& endCoord) {
    double distance = 0.0;
    double totalFare = 0.0;

    std::vector<const RouteNode*> takenRoutes = extractTraversedRouteNodes(path);
    if (takenRoutes.empty()) return 0.0;

    Coordinate projectedStart = projectOntoPath(startCoord, takenRoutes.front()->path);
    Coordinate projectedEnd = projectOntoPath(endCoord, takenRoutes.back()->path);

    for (size_t i = 0; i < path.size(); ++i) {
        const TFTFEdge& edge = path[i];

        if (i == 0) {
            distance += getActualSegmentDistance(projectedStart, edge.exitCoord, takenRoutes[0]->path);
        } else if (i == path.size() - 1) {
            distance += getActualSegmentDistance(edge.entryCoord, projectedEnd, takenRoutes.back()->path);
        } else {
            distance += getActualSegmentDistance(edge.entryCoord, edge.exitCoord, takenRoutes[i]->path);
        }
    }

    totalFare = BASE_FARE * takenRoutes.size();
    totalFare += ((distance / 1000.0) * FARE_PER_KM);
    return totalFare;
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
    std::unordered_map<int, TFTFEdge> dijkstraPrevEdge;      
    std::unordered_map<int, int> prevRouteMap;               
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

int TFTFGraph::findClosestRoute(const Coordinate& startCoord) {
    float closestDistance = std::numeric_limits<float>::infinity();
    int closestRouteId = -1;

    for (const auto& [routeId, routeNode] : routes) {
        float dist = haversine(routeNode.path.front(), startCoord); 
        if (dist < closestDistance) {
            closestDistance = dist;
            closestRouteId = routeId;
        }
    }

    return closestRouteId;
}

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

    if (!path.empty()) {
        TFTFEdge startEdge;
        startEdge.destinationRoute = path.front().destinationRoute; 
        startEdge.destinationRouteName = routes[startRouteId].routeName;
        startEdge.transferCost = 0.0f;
        startEdge.densityByInterval = {};
        startEdge.entryCoord = startCoord;
        startEdge.exitCoord = path.front().entryCoord;
        path.insert(path.begin(), startEdge);
    }

    if (!path.empty()) {
        TFTFEdge endEdge;
        endEdge.destinationRoute = -1; 
        endEdge.destinationRouteName = "Destination";
        endEdge.transferCost = 0.0f;
        endEdge.densityByInterval = {};
        endEdge.entryCoord = path.back().exitCoord;
        endEdge.exitCoord = endCoord;
        path.push_back(endEdge);
    }

    std::cout << "Best Path: Start -> ";
    for (size_t i = 0; i < path.size(); ++i) {
        if (i > 0) {
            std::cout << " -> ";
        }
        std::cout << path[i].destinationRouteName;
    }
    std::cout << std::endl;

    float totalFare = 0.0f;
    totalFare = calculateTotalFare(path, startCoord, endCoord);

    std::cout << "Total fare: " << std::fixed << std::setprecision(2) << totalFare << " pesos" << std::endl;

    return path;
}