#define _USE_MATH_DEFINES
#include <cmath>

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

bool operator == (const Coordinate& lhs, const Coordinate& rhs) {
    return lhs.latitude == rhs.latitude && lhs.longitude == rhs.longitude;
}


float TFTFEdge::totalCost() const
{

    return (transferCost);
}

void TFTFGraph::addRoute(int routeId, const std::string& routeName) {
    routes[routeId] = {routeId, routeName, {}};
}

// sets the path (sequence of coordinates) for a given route
// TO MAKE LOOP MAKE LAST COORD === TO FIRST COORD
void TFTFGraph::setRoutePath(int routeId, const std::vector<Coordinate>& coordinates) {
    if (routes.find(routeId) != routes.end()) {
        routes[routeId].path = coordinates;
        routes[routeId].isLoop = (!coordinates.empty() && coordinates.front() == coordinates.back());

    } else {
        std::cerr << "route ID " << routeId << " not found.\n";
    }
}

// adds an edge (connection) between two routes in the graph
void TFTFGraph::addEdge(
    int fromRoute, 
    int toRoute, 
    const std::string& toName,
    float transferCost,
    Coordinate entryCoord, 
    Coordinate exitCoord) {


    int entryIndex = getClosestIndex(routes[fromRoute].path, entryCoord);
    int exitIndex = getClosestIndex(routes[toRoute].path, exitCoord);


    // add a new edge from 'fromRoute' to 'toRoute' with relevant information
    routes[fromRoute].edges.push_back({toRoute, fromRoute, toName, transferCost});

    // set the entry and exit coordinates for the newly added edge
    routes[fromRoute].edges.back().entryCoord = entryCoord;
    routes[fromRoute].edges.back().exitCoord = exitCoord;

    
    //set the entry and exit indices for the newly added edge
    routes[fromRoute].edges.back().entryIndex = entryIndex;
    routes[fromRoute].edges.back().exitIndex = exitIndex;
    
}

std::vector<const RouteNode*> TFTFGraph::extractTraversedRouteNodes(
    const std::vector<TFTFEdge>& path) const {
    
    std::vector<const RouteNode*> routeNodes;

    if (path.empty()) return routeNodes;

    // check if the origin route of the first edge exists in the routes map
    if (routes.count(path.front().originRoute)) {
        // add the origin route node of the first edge to the result
        routeNodes.push_back(&routes.at(path.front().originRoute));
    }

    // iterate through the edges in the given path
    for (const auto& edge : path) {
        // check if the destination route of the current edge exists in the routes map
        if (routes.count(edge.destinationRoute)) {
            // add the destination route node to the result
            routeNodes.push_back(&routes.at(edge.destinationRoute));
        }
    }

    return routeNodes;
}

// calculates the total fare for a given path based on distance and route nodes
double TFTFGraph::calculateTotalFare(
    const std::vector<TFTFEdge>& path, 
    const Coordinate& startCoord, 
    const Coordinate& endCoord) {

    double distance = 0.0;  
    double totalFare = 0.0; 

    // get the list of route nodes traversed by the path
    std::vector<const RouteNode*> takenRoutes = extractTraversedRouteNodes(path);

    // if no valid routes are found, return 0.0 (no fare)
    if (takenRoutes.empty()) return 0.0;

    // project the start and end coordinates onto the respective route paths
    Coordinate projectedStart = projectOntoPath(startCoord, takenRoutes.front()->path);
    Coordinate projectedEnd = projectOntoPath(endCoord, takenRoutes.back()->path);

    // loop through the edges in the path to calculate total distance
    for (size_t i = 0; i < path.size(); ++i) {
        const TFTFEdge& edge = path[i];

        // if it's the first edge, calculate distance from start to the first edge's exit
        if (i == 0) {
            int startIdx = getClosestIndex(takenRoutes[0]->path, projectedStart);
            int endIdx = getClosestIndex(takenRoutes[0]->path, edge.exitCoord);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes[0]->isLoop && startIdx > endIdx) {
                continue;
            }
            
            distance += getActualSegmentDistance(projectedStart, edge.exitCoord, takenRoutes[0]->path, takenRoutes[0]->isLoop);
        }
        // if it's the last edge, calculate distance from the last edge's entry to the end
        else if (i == path.size() - 1) {
            int startIdx = getClosestIndex(takenRoutes.back()->path, edge.entryCoord);
            int endIdx = getClosestIndex(takenRoutes.back()->path, projectedEnd);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes.back()->isLoop && startIdx > endIdx) {
                continue;
            }
            distance += getActualSegmentDistance(edge.entryCoord, projectedEnd, takenRoutes.back()->path, takenRoutes.back()->isLoop);
        }
        // for all other edges, calculate distance between entry and exit coordinates
        else {
            int startIdx = getClosestIndex(takenRoutes[i]->path, edge.entryCoord);
            int endIdx = getClosestIndex(takenRoutes[i]->path, edge.exitCoord);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes[i]->isLoop && startIdx > endIdx) {
                continue;
            }

            distance += getActualSegmentDistance(edge.entryCoord, edge.exitCoord, takenRoutes[i]->path, takenRoutes[i]->isLoop);
        }
    }

    // calculate total fare based on the number of routes taken and the total distance
    totalFare = BASE_FARE * takenRoutes.size();  // add fare for number of routes taken
    totalFare += ((distance / 1000.0) * FARE_PER_KM);  // add fare based on distance (converted to kilometers)

    return totalFare;
}

std::pair<int, int> hashCoordinate(const Coordinate& coord, float cellSizeMeters) {
    float latMeters = 111320.0f; // meters per degree latitude
    float lonMeters = 111320.0f * std::cos(coord.latitude * M_PI / 180.0); // meters per degree longitude

    int x = static_cast<int>(coord.longitude * lonMeters / cellSizeMeters);
    int y = static_cast<int>(coord.latitude * latMeters / cellSizeMeters);
    return {x, y};
}

void TFTFGraph::createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer) {
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, Coordinate>>, PairHash> spatialGrid;

    // Step 1: Populate the spatial grid
    for (auto& [routeID, route] : routes) {
        auto densePath = densifyPath(route.path, 25.0f); // Every 25 meters
        for (const auto& coord : densePath) {
            auto cell = hashCoordinate(coord, transferRangeMeters);
            spatialGrid[cell].emplace_back(routeID, coord);
        }
    }

    // Step 2: Check nearby cells for potential transfers
    for (const auto& [fromID, fromRoute] : routes) {
        for (const auto& fromCoord : fromRoute.path) {
            auto baseCell = hashCoordinate(fromCoord, transferRangeMeters);

            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    std::pair<int, int> neighborCell = {baseCell.first + dx, baseCell.second + dy};
                    auto it = spatialGrid.find(neighborCell);
                    if (it == spatialGrid.end()) continue;

                    for (const auto& [toID, toCoord] : it->second) {
                        if (toID == fromID) continue;

                        float dist = haversine(fromCoord, toCoord);
                        if (dist <= transferRangeMeters) {
                            bool exists = false;
                            for (const auto& edge : routes[fromID].edges) {
                                if (edge.destinationRoute == toID && std::abs(edge.transferCost - dist) < 1e-2f) {
                                    exists = true;
                                    break;
                                }
                            }

                            if (!exists) {
                                addEdge(fromID, toID, routes[toID].routeName, dist, fromCoord, toCoord);
                            }
                        }
                    }
                }
            }
        }
    }
}
void TFTFGraph::visualize() const
    {
        std::cout << "\nTFTF Graph Visualization";
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

                std::cout << "\n";
            }
        }
    }

std::vector<int> TFTFGraph::getNearbyRoutes(const Coordinate& coord, float maxDistanceMeters) {
    std::vector<int> nearby;

    for (const auto& [routeId, route] : routes) {
        auto densePath = densifyPath(route.path, 25.0f);
      
        for (const auto& point : densePath) {
            
            // Check if the point is within the specified distance from the coordinate
            if (haversine(coord, point) <= maxDistanceMeters) {
                
                nearby.push_back(routeId);
                break; // only need one point close enough
            }
        }
    }

    return nearby;
}


std::vector<TFTFEdge> TFTFGraph::findMinFarePath(
    int startRouteId,
    int endRouteId,
    int hour,
    int projectedStartIdx) {

    constexpr float BASE_FARE = 12.0f;
    constexpr float FARE_PER_KM = 1.5f;

    struct FareState {
        float totalFare;
        int lastCoordIndex;
    };

    struct FareNode {
        int routeId;
        int entryIndex;
        float fareSoFar;
        int lastCoordIndex;

        bool operator>(const FareNode& other) const {
            return fareSoFar > other.fareSoFar;
        }
    };

    // Best[routeId][entryIndex] = FareState
    std::unordered_map<int, std::unordered_map<int, FareState>> best;
    std::unordered_map<std::pair<int, int>, TFTFEdge, PairHash> prevEdge;
    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHash> prevRoute;

    std::priority_queue<FareNode, std::vector<FareNode>, std::greater<FareNode>> pq;

    // Start with 0 fare, unknown entry index yet
    best[startRouteId][-1] = {0.0f, -1};
    pq.push({startRouteId, -1, 0.0f, -1});

    while (!pq.empty()) {
        FareNode curr = pq.top(); pq.pop();
        int routeId = curr.routeId;
        int lastCoord = curr.lastCoordIndex;
        float fareSoFar = curr.fareSoFar;

        for (const TFTFEdge& edge : routes.at(routeId).edges) {
            int nextRoute = edge.destinationRoute;
            int entry = edge.entryIndex;
            int exit = edge.exitIndex;
            bool isLoop = routes.at(routeId).isLoop;

            // Prevent backward movement on non-loop routes
            if (!isLoop) {
                // prevent backward movement even if it's the first step (no lastCoord)
                if (entry >= exit) continue;

                // also skip if continuing and moving backwards
                if (lastCoord != -1 && entry <= lastCoord) continue;

                //If only transfers ocur before onboard on non looping routes
                if (entry < projectedStartIdx) continue;
            }

            float distance = 0.0f;
            
            if (lastCoord == -1 || nextRoute != routeId) {
                // First leg OR transfer — full fare
                distance = getSubpathDistance(routes.at(nextRoute).path, edge.entryIndex, edge.exitIndex, routes.at(nextRoute).isLoop);
            } else {
                // Continuing on same route — just pay by distance
                distance = getSubpathDistance(routes.at(routeId).path, entry, exit, isLoop);
            }

            float dynamicCost = edge.totalCost();
            float fare = BASE_FARE + (distance / 1000.0f) * FARE_PER_KM + dynamicCost;
            float nextFare = fareSoFar + fare;

            auto& bestMap = best[nextRoute];
            if (!bestMap.count(exit) || nextFare < bestMap[exit].totalFare) {
                bestMap[exit] = {nextFare, exit};
                pq.push({nextRoute, entry, nextFare, exit});
                prevEdge[{nextRoute, exit}] = edge;
                prevRoute[{nextRoute, exit}] = {routeId, lastCoord};
            }
        }
    }

    // Reconstruct the path
    std::vector<TFTFEdge> path;
    const auto& endMap = best[endRouteId];
    if (endMap.empty()) {
        std::cout << "No fare-optimal path found.\n";
        return {};
    }

    // Find best final exitIndex from the end route
    auto minIt = std::min_element(endMap.begin(), endMap.end(),
        [](const auto& a, const auto& b) {
            return a.second.totalFare < b.second.totalFare;
        });

    int currRoute = endRouteId;
    int currCoord = minIt->first;

    while (currRoute != startRouteId || currCoord != -1) {
        auto it = prevEdge.find({currRoute, currCoord});
        if (it == prevEdge.end()) break;
        const auto& edge = it->second;
        path.push_back(edge);

        auto prev = prevRoute[{currRoute, currCoord}];
        currRoute = prev.first;
        currCoord = prev.second;
    }

    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<TFTFEdge> TFTFGraph::calculateRouteFromCoordinates(
    const Coordinate& startCoord, 
    const Coordinate& endCoord, 
    int hour) {

    std::vector<int> startCandidates = getNearbyRoutes(startCoord, 300.0f);
    std::vector<int> endCandidates = getNearbyRoutes(endCoord, 300.0f);

    if (startCandidates.empty()) {
        std::cerr << "No nearby routes found for the start coordinates.\n";
        return {};
    } else if (endCandidates.empty()) {
        std::cerr << "No nearby routes found for the end coordinates.\n";
        return {};
    }

    std::vector<TFTFEdge> bestPath;
    double bestFare = std::numeric_limits<double>::infinity();
    int bestStartRouteId = -1;
    int bestEndRouteId = -1;
    bool bestPathSameRoute = false;

    for (int startId : startCandidates) {
        for (int endId : endCandidates) {
            std::vector<TFTFEdge> path;
            double fare = 0.0;

            if (startId == endId) {
                double segmentDistance = getActualSegmentDistance(startCoord, endCoord, routes.at(startId).path, routes.at(startId).isLoop);
                fare = BASE_FARE + (segmentDistance / 1000.0) * FARE_PER_KM;

                // Simulate one direct TFTFEdge for uniformity
                TFTFEdge edge;
                edge.destinationRoute = startId;
                edge.destinationRouteName = routes.at(startId).routeName;
                edge.transferCost = 0.0f;
                edge.entryCoord = projectOntoPath(startCoord, routes[startId].path);
                edge.exitCoord = projectOntoPath(endCoord, routes[endId].path);
                path.push_back(edge);
            } else {
                int projStart = getClosestIndex(routes.at(startId).path, startCoord);
                path = findMinFarePath(startId, endId, hour, projStart);
                if (path.empty()) continue;
                fare = calculateTotalFare(path, startCoord, endCoord);
            }

            if (fare < bestFare) {
                bestFare = fare;
                bestPath = path;
                bestStartRouteId = startId;
                bestEndRouteId = endId;
                bestPathSameRoute = (startId == endId);
            }
        }
    }

    if (bestPath.empty()) {
        std::cerr << "No valid route found between coordinates.\n";
        return {};
    }

if (!bestPathSameRoute) {
    // Insert start edge
    TFTFEdge startEdge;
    startEdge.destinationRoute = bestPath.front().destinationRoute;
    startEdge.destinationRouteName = routes[bestStartRouteId].routeName;
    startEdge.transferCost = 0.0f;
    startEdge.entryCoord = projectOntoPath(startCoord, routes[bestStartRouteId].path);
    startEdge.exitCoord = bestPath.front().entryCoord;
    bestPath.insert(bestPath.begin(), startEdge);
}
    // Insert end edge
    TFTFEdge endEdge;
    endEdge.destinationRoute = -1;
    endEdge.destinationRouteName = "Destination";
    endEdge.transferCost = 0.0f;
    endEdge.entryCoord = bestPath.back().exitCoord;
    endEdge.exitCoord = endCoord;
    bestPath.push_back(endEdge);



    std::cout << "\n==== Route Instructions ====\n";
    for (size_t i = 0; i < bestPath.size(); ++i) {
        const TFTFEdge& edge = bestPath[i];

        if (i == 0) {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "Start at coordinates: (" 
                      << edge.entryCoord.latitude << ", " 
                      << edge.entryCoord.longitude << ")\n";
        }
// cleant hi sh
        if (edge.destinationRoute != -1) {
            std::cout << "Take route: " << edge.destinationRouteName << "\n";
            std::cout << "  Mount at: (" 
                      << edge.entryCoord.latitude << ", " 
                      << edge.entryCoord.longitude << ")\n";
            std::cout << "  Dismount at: (" 
                      << edge.exitCoord.latitude << ", " 
                      << edge.exitCoord.longitude << ")\n";
        } else {
            std::cout << "Walk to final destination at: (" 
                      << edge.exitCoord.latitude << ", " 
                      << edge.exitCoord.longitude << ")\n";
        }
    }
    std::cout << "=============================\n";

    float totalFare = calculateTotalFare(bestPath, startCoord, endCoord);
    std::cout << "Total fare: " << std::fixed << std::setprecision(2) 
              << roundUpToNearest2_5(totalFare) << " pesos" << std::endl;

    return bestPath;
}
