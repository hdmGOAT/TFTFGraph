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

// function to compute the average jeepney density per hour between two routes
std::vector<JeepneyDensity> averageRouteDensities(
    const std::vector<JeepneyDensity>& a, 
    const std::vector<JeepneyDensity>& b) {
    
    std::vector<JeepneyDensity> result;

    // maps to hold hour-to-density values for both route a and b
    std::unordered_map<int, float> mapA, mapB;

    // fill mapA with jeepney density values for each hour range in route a
    for (const auto& d : a) {
        // assign the same density value for every hour between startHour and endHour
        for (int h = d.startHour; h < d.endHour; ++h) {
            mapA[h] = d.jeepneyDensity;
        }
    }

    // fill mapB with jeepney density values for each hour range in route b
    for (const auto& d : b) {
        for (int h = d.startHour; h < d.endHour; ++h) {
            mapB[h] = d.jeepneyDensity;
        }
    }

    // loop through each hour of the day (0 to 23)
    for (int h = 0; h < 24; ++h) {
        // use the density from mapA if available, else assume 1.0
        float dA = mapA.count(h) ? mapA[h] : 1.0f;

        float dB = mapB.count(h) ? mapB[h] : 1.0f;

        float avg = (dA + dB) / 2.0f;

        // store the average as a new JeepneyDensity for the current hour
        result.push_back({h, h + 1, avg});
    }

    return result;
}


float TFTFEdge::totalCost(int hour) const
{
    float densityFactor = 1.0f;

    if (hour >= 0)
    {
        for (const auto& interval : densityByInterval)
        {
            // if the current interval contains the hour, update the density factor
            if (interval.contains(hour))
            {
                densityFactor = interval.jeepneyDensity;
                break;
            }
        }
    }

    return (transferCost) * densityFactor;
}

void TFTFGraph::addRoute(int routeId, const std::string& routeName) {
    routes[routeId] = {routeId, routeName, {}};
}

// sets the path (sequence of coordinates) for a given route
void TFTFGraph::setRoutePath(int routeId, const std::vector<Coordinate>& coordinates) {
    if (routes.find(routeId) != routes.end()) {
        routes[routeId].path = coordinates;
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
    const std::vector<JeepneyDensity>& densities, 
    Coordinate entryCoord, 
    Coordinate exitCoord) {
    
    // calculate the average density between the two routes for each hour
    std::vector<JeepneyDensity> avgDensity = averageRouteDensities(routes[fromRoute].densities, routes[toRoute].densities);

    int entryIndex = getClosestIndex(routes[fromRoute].path, entryCoord);
    int exitIndex = getClosestIndex(routes[toRoute].path, exitCoord);


    // add a new edge from 'fromRoute' to 'toRoute' with relevant information
    routes[fromRoute].edges.push_back({toRoute, fromRoute, toName, transferCost, avgDensity});

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
            distance += getActualSegmentDistance(projectedStart, edge.exitCoord, takenRoutes[0]->path);
        }
        // if it's the last edge, calculate distance from the last edge's entry to the end
        else if (i == path.size() - 1) {
            distance += getActualSegmentDistance(edge.entryCoord, projectedEnd, takenRoutes.back()->path);
        }
        // for all other edges, calculate distance between entry and exit coordinates
        else {
            distance += getActualSegmentDistance(edge.entryCoord, edge.exitCoord, takenRoutes[i]->path);
        }
    }

    // calculate total fare based on the number of routes taken and the total distance
    totalFare = BASE_FARE * takenRoutes.size();  // add fare for number of routes taken
    totalFare += ((distance / 1000.0) * FARE_PER_KM);  // add fare based on distance (converted to kilometers)

    return totalFare;
}

void TFTFGraph::setRouteDensities(int routeId, const std::vector<JeepneyDensity>& densities) {
    if (routes.find(routeId) != routes.end()) {
        routes[routeId].densities = densities;
    } else {
        std::cerr << "Route ID " << routeId << " not found.\n";
    }
}

// creates transfer edges between routes based on the distance between coordinates
void TFTFGraph::createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer) {
    // loop through all pairs of routes to check for possible transfers
    for (const auto& [fromId, fromNode] : routes) {
        for (const auto& [toId, toNode] : routes) {
            // skip creating transfers from a route to itself
            if (fromId == toId) continue;

            // loop through the coordinates of the 'from' route node
            for (const auto& fromCoord : fromNode.path) {
                // loop through the coordinates of the 'to' route node
                for (const auto& toCoord : toNode.path) {
                    // calculate the distance between the two coordinates using the Haversine formula
                    float dist = haversine(fromCoord, toCoord);

                    // check if the distance is within the transfer range
                    if (dist <= transferRangeMeters) {
                        bool exists = false;

                        // check if an edge already exists between these two routes with the same distance
                        for (const auto& edge : fromNode.edges) {
                            if (edge.destinationRoute == toId && std::abs(edge.transferCost - dist) < 1e-2f) {
                                exists = true;  // transfer edge already exists
                                break;
                            }
                        }

                        // if no such edge exists, create a new transfer edge
                        if (!exists) {
                            addEdge(fromId, toId, toNode.routeName, dist, {}, fromCoord, toCoord);
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

std::vector<TFTFEdge> TFTFGraph::findMinTransfersPath(
    int startRouteId, 
    int endRouteId, 
    int hour) {
    
    // a helper struct to store information about the number of transfers and the total cost.
    struct TransferInfo {
        int transfers; // the number of transfers
        float cost;    // the total cost to reach the route
        int lastCoordIndex = -1; // the index of the last coordinate in the path
    };

    // maps for storing the best (minimum) transfers and cost for each route,
    // as well as previous edge and route information for path reconstruction.
    std::unordered_map<int, TransferInfo> best;
    std::unordered_map<int, TFTFEdge> prevEdge;
    std::unordered_map<int, int> prevRoute;

    // a queue for breadth-first search (BFS) to explore routes.
    std::queue<int> q;
    // a set to track routes currently in the queue, preventing reprocessing.
    std::unordered_set<int> inQueue;

    // initialize all routes with a "worst-case" transfer count and cost.
    // set each route's best transfer count to max and cost to infinity.
    for (const auto& [routeId, _] : routes) {
        best[routeId] = {std::numeric_limits<int>::max(), std::numeric_limits<float>::infinity(), -1};
    }

    // start route has 0 transfers and 0 cost.
    best[startRouteId] = {0, 0.0f};
    q.push(startRouteId);  // start processing from the start route.
    inQueue.insert(startRouteId);  // mark the start route as visited.

    // perform a breadth-first search (BFS) to find the route with the minimum transfers
    // and lowest cost to each route.
    while (!q.empty()) {
        int curr = q.front(); // get the current route.
        q.pop();              // remove the route from the queue.
        inQueue.erase(curr);  // mark it as not in the queue.

        // iterate through each edge from the current route.
        for (const auto& edge : routes.at(curr).edges) {
            int next = edge.destinationRoute;

            // Get the index of the current edge's exit coordinate
            int currIndex = edge.entryIndex;

            // Skip backward or redundant transfers
            if (currIndex <= best[curr].lastCoordIndex) {
                continue;
            }

            float edgeCost = edge.totalCost(hour);
            int nextTransfers = best[curr].transfers + 1;
            float nextCost = best[curr].cost + edgeCost;

            if (nextTransfers < best[next].transfers ||
                (nextTransfers == best[next].transfers && nextCost < best[next].cost)) {
                best[next] = {nextTransfers, nextCost, edge.exitIndex};
                prevEdge[next] = edge;
                prevRoute[next] = curr;

                if (!inQueue.count(next)) {
                    q.push(next);
                    inQueue.insert(next);
                }
            }
        }

    }

    // reconstruct the path from the end route to the start route.
    std::vector<TFTFEdge> path;
    int curr = endRouteId;  // start from the end route.

    // traverse backwards from the end route, following the `prevRoute` map.
    while (curr != startRouteId) {
        if (prevEdge.find(curr) == prevEdge.end()) {
            // if no valid path is found, return an empty path.
            std::cout << "No valid path (min transfers + low cost) found.\n";
            return {};
        }
        // add the edge from the previous route to the path.
        path.push_back(prevEdge[curr]);
        // move to the previous route.
        curr = prevRoute[curr];
    }

    // reverse the path since we built it backwards (from end to start).
    std::reverse(path.begin(), path.end());

    // return the reconstructed path.
    return path;
}

std::vector<TFTFEdge> TFTFGraph::calculateRouteFromCoordinates(
    const Coordinate& startCoord, 
    const Coordinate& endCoord, 
    int hour) {
    
    int startRouteId = findClosestRoute(startCoord);
    int endRouteId = findClosestRoute(endCoord);

    if (startRouteId == -1 || endRouteId == -1) {
        std::cerr << "Could not find a route for the given coordinates." << std::endl;
        return {};
    }

    // print the closest route IDs
    std::cout << "Closest route from start: " << startRouteId << "\n";
    std::cout << "Closest route from end: " << endRouteId << "\n";

    // find the path with the minimum transfers between the start and end route
    std::vector<TFTFEdge> path = findMinTransfersPath(startRouteId, endRouteId, hour);

    // if the path is not empty, create a start edge for the path
    if (!path.empty()) {
        TFTFEdge startEdge;
        startEdge.destinationRoute = path.front().destinationRoute; 
        startEdge.destinationRouteName = routes[startRouteId].routeName;
        startEdge.transferCost = 0.0f;
        startEdge.densityByInterval = {};
        startEdge.entryCoord = startCoord;
        startEdge.exitCoord = path.front().entryCoord;
        path.insert(path.begin(), startEdge); // insert the start edge at the beginning of the path
    }

    // if the path is not empty, create an end edge for the path
    if (!path.empty()) {
        TFTFEdge endEdge;
        endEdge.destinationRoute = -1; 
        endEdge.destinationRouteName = "Destination";
        endEdge.transferCost = 0.0f;
        endEdge.densityByInterval = {};
        endEdge.entryCoord = path.back().exitCoord;
        endEdge.exitCoord = endCoord;
        path.push_back(endEdge); // add the end edge to the path
    }

    // print the best path
    std::cout << "Best Path: Start -> ";
    for (size_t i = 0; i < path.size(); ++i) {
        if (i > 0) {
            std::cout << " -> ";
        }
        std::cout << path[i].destinationRouteName;
    }
    std::cout << std::endl;

    // calculate the total fare for the path
    float totalFare = 0.0f;
    totalFare = calculateTotalFare(path, startCoord, endCoord);

    // print the total fare
    std::cout << "Total fare: " << std::fixed << std::setprecision(2) << totalFare << " pesos" << std::endl;

    return path;
}
