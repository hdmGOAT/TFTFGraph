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

bool operator==(const Coordinate &lhs, const Coordinate &rhs)
{
    return lhs.latitude == rhs.latitude && lhs.longitude == rhs.longitude;
}

float TFTFEdge::totalCost() const
{

    return (transferCost);
}

void TFTFGraph::addRoute(int routeId, const std::string &routeName)
{
    routes[routeId] = {routeId, routeName, {}};
}

// sets the path (sequence of coordinates) for a given route
// TO MAKE LOOP MAKE LAST COORD === TO FIRST COORD
void TFTFGraph::setRoutePath(int routeId, const std::vector<Coordinate> &coordinates)
{
    if (routes.find(routeId) != routes.end())
    {
        routes[routeId].path = coordinates;
        routes[routeId].isLoop = (!coordinates.empty() && coordinates.front() == coordinates.back());
    }
    else
    {
        std::cerr << "route ID " << routeId << " not found.\n";
    }
}

// adds an edge (connection) between two routes in the graph
void TFTFGraph::addEdge(
    int fromRoute,
    int toRoute,
    const std::string &toName,
    float transferCost,
    Coordinate entryCoord,
    Coordinate exitCoord)
{
    int entryIndex = getClosestIndex(routes[fromRoute].path, entryCoord);
    int exitIndex = getClosestIndex(routes[toRoute].path, exitCoord);

    // Create a TFTFEdge explicitly
    TFTFEdge edge;
    edge.type = EdgeType::Transfer; // <== ✅ SET TYPE HERE
    edge.destinationRoute = toRoute;
    edge.originRoute = fromRoute;
    edge.destinationRouteName = toName;
    edge.transferCost = transferCost;
    edge.entryCoord = entryCoord;
    edge.exitCoord = exitCoord;
    edge.entryIndex = entryIndex;
    edge.exitIndex = exitIndex;

    // Push to the edges list
    routes[fromRoute].edges.push_back(edge);
}

std::vector<const RouteNode *> TFTFGraph::extractTraversedRouteNodes(
    const std::vector<TFTFEdge> &path) const
{

    std::vector<const RouteNode *> routeNodes;

    if (path.empty())
        return routeNodes;

    // check if the origin route of the first edge exists in the routes map
    if (routes.count(path.front().originRoute))
    {
        // add the origin route node of the first edge to the result
        routeNodes.push_back(&routes.at(path.front().originRoute));
    }

    // iterate through the edges in the given path
    for (const auto &edge : path)
    {
        // check if the destination route of the current edge exists in the routes map
        if (routes.count(edge.destinationRoute))
        {
            // add the destination route node to the result
            routeNodes.push_back(&routes.at(edge.destinationRoute));
        }
    }

    return routeNodes;
}

// calculates the total fare for a given path based on distance and route nodes
double TFTFGraph::calculateTotalFare(
    const std::vector<TFTFEdge> &path,
    const Coordinate &startCoord,
    const Coordinate &endCoord)
{

    double distance = 0.0;
    double totalFare = 0.0;

    // get the list of route nodes traversed by the path
    std::vector<const RouteNode *> takenRoutes = extractTraversedRouteNodes(path);

    // if no valid routes are found, return 0.0 (no fare)
    if (takenRoutes.empty())
        return 0.0;

    // project the start and end coordinates onto the respective route paths
    Coordinate projectedStart = projectOntoPath(startCoord, takenRoutes.front()->path);
    Coordinate projectedEnd = projectOntoPath(endCoord, takenRoutes.back()->path);

    // loop through the edges in the path to calculate total distance
    for (size_t i = 0; i < path.size(); ++i)
    {
        const TFTFEdge &edge = path[i];

        // if it's the first edge, calculate distance from start to the first edge's exit
        if (i == 0)
        {
            int startIdx = getClosestIndex(takenRoutes[0]->path, projectedStart);
            int endIdx = getClosestIndex(takenRoutes[0]->path, edge.exitCoord);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes[0]->isLoop && startIdx > endIdx)
            {
                continue;
            }

            distance += getActualSegmentDistance(projectedStart, edge.exitCoord, takenRoutes[0]->path, takenRoutes[0]->isLoop);
        }
        // if it's the last edge, calculate distance from the last edge's entry to the end
        else if (i == path.size() - 1)
        {
            int startIdx = getClosestIndex(takenRoutes.back()->path, edge.entryCoord);
            int endIdx = getClosestIndex(takenRoutes.back()->path, projectedEnd);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes.back()->isLoop && startIdx > endIdx)
            {
                continue;
            }
            distance += getActualSegmentDistance(edge.entryCoord, projectedEnd, takenRoutes.back()->path, takenRoutes.back()->isLoop);
        }
        // for all other edges, calculate distance between entry and exit coordinates
        else
        {
            int startIdx = getClosestIndex(takenRoutes[i]->path, edge.entryCoord);
            int endIdx = getClosestIndex(takenRoutes[i]->path, edge.exitCoord);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes[i]->isLoop && startIdx > endIdx)
            {
                continue;
            }

            distance += getActualSegmentDistance(edge.entryCoord, edge.exitCoord, takenRoutes[i]->path, takenRoutes[i]->isLoop);
        }
    }

    // calculate total fare based on the number of routes taken and the total distance
    totalFare = BASE_FARE * takenRoutes.size();       // add fare for number of routes taken
    totalFare += ((distance / 1000.0) * FARE_PER_KM); // add fare based on distance (converted to kilometers)

    return totalFare;
}

std::pair<int, int> hashCoordinate(const Coordinate &coord, float cellSizeMeters)
{
    float latMeters = 111320.0f;                                           // meters per degree latitude
    float lonMeters = 111320.0f * std::cos(coord.latitude * M_PI / 180.0); // meters per degree longitude

    int x = static_cast<int>(coord.longitude * lonMeters / cellSizeMeters);
    int y = static_cast<int>(coord.latitude * latMeters / cellSizeMeters);
    return {x, y};
}

void TFTFGraph::createTransfersFromCoordinates(float transferRangeMeters)
{
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, Coordinate>>, PairHash> spatialGrid;

    // Step 1: Populate the spatial grid with densified paths
    for (auto &[routeID, route] : routes)
    {
        auto densePath = densifyPath(route.path, 100.0f); // Every 100 meters
        for (const auto &coord : densePath)
        {
            auto cell = hashCoordinate(coord, transferRangeMeters);
            spatialGrid[cell].emplace_back(routeID, coord);
        }
    }

    // Step 2: Add best transfer per target route every ~500 meters
    for (const auto &[fromID, fromRoute] : routes)
    {
        std::unordered_map<int, std::pair<Coordinate, float>> lastTransferPerRoute; // toID -> (lastFromCoord, distanceSinceLastTransfer)

        for (const auto &fromCoord : fromRoute.path)
        {
            auto baseCell = hashCoordinate(fromCoord, transferRangeMeters);

            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    std::pair<int, int> neighborCell = {baseCell.first + dx, baseCell.second + dy};
                    auto it = spatialGrid.find(neighborCell);
                    if (it == spatialGrid.end())
                        continue;

                    for (const auto &[toID, toCoord] : it->second)
                    {
                        if (toID == fromID)
                            continue;

                        float dist = haversine(fromCoord, toCoord);
                        if (dist > transferRangeMeters)
                            continue;

                        // Check distance from last transfer to this route
                        auto &entry = lastTransferPerRoute[toID];
                        float fromLastDist = (entry.second > 0.0f) ? haversine(entry.first, fromCoord) : std::numeric_limits<float>::max();

                        if (fromLastDist >= 500.0f || entry.second == 0.0f) // No transfer yet or far enough
                        {
                            addEdge(fromID, toID, routes[toID].routeName, dist, fromCoord, toCoord);
                            entry = {fromCoord, dist}; // Update last transfer
                        }
                        else if (dist < entry.second) // Closer transfer candidate
                        {
                            // Optional: update best known transfer if you want to keep the best within a 500m segment
                            entry = {fromCoord, dist};
                            // You could re-add/replace the edge if desired (requires edge removal or deduplication logic)
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

std::vector<int> TFTFGraph::getNearbyRoutes(const Coordinate &coord, float maxDistanceMeters)
{
    std::vector<int> nearby;

    for (const auto &[routeId, route] : routes)
    {
        auto densePath = densifyPath(route.path, 25.0f);
        for (const auto &point : densePath)
        {

            // Check if the point is within the specified distance from the coordinate
            if (haversine(coord, point) <= maxDistanceMeters)
            {

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
    int projectedStartIdx) // still passed for compatibility; unused in this version
{
    constexpr float BASE_FARE = 12.0f;
    constexpr float FARE_PER_KM = 1.5f;

    struct FareState
    {
        float totalFare;
        int prevRouteId;
    };

    struct FareNode
    {
        int routeId;
        float fareSoFar;

        bool operator>(const FareNode &other) const
        {
            return fareSoFar > other.fareSoFar;
        }
    };

    std::unordered_map<int, FareState> best;
    std::unordered_map<int, TFTFEdge> prevEdge;
    std::unordered_map<int, int> prevRoute;

    std::priority_queue<FareNode, std::vector<FareNode>, std::greater<FareNode>> pq;

    best[startRouteId] = {0.0f, -1};
    pq.push({startRouteId, 0.0f});

    while (!pq.empty())
    {
        FareNode curr = pq.top();
        pq.pop();

        int routeId = curr.routeId;
        float fareSoFar = curr.fareSoFar;

        for (const TFTFEdge &edge : routes.at(routeId).edges)
        {
            int nextRoute = edge.destinationRoute;

            if (routeId == startRouteId && !routes.at(routeId).isLoop)
            {
                if (edge.exitIndex < projectedStartIdx)
                    continue;
            }

            float fare = 0.0f;

            // Handle transfers vs regular rides
            if (edge.type == EdgeType::Transfer)
            {
                fare = edge.transferCost; // could be 0 or flat
            }
            else
            {
                float distance = getSubpathDistance(
                    routes.at(nextRoute).path,
                    edge.entryIndex,
                    edge.exitIndex,
                    routes.at(nextRoute).isLoop);

                fare = BASE_FARE + (distance / 1000.0f) * FARE_PER_KM + edge.totalCost();
            }

            float nextFare = fareSoFar + fare;

            if (!best.count(nextRoute) || nextFare < best[nextRoute].totalFare)
            {
                best[nextRoute] = {nextFare, routeId};
                pq.push({nextRoute, nextFare});
                prevEdge[nextRoute] = edge;
                prevRoute[nextRoute] = routeId;
            }
        }
    }

    std::vector<TFTFEdge> path;
    if (!best.count(endRouteId))
    {
        std::cout << "No path found from route " << startRouteId << " to route " << endRouteId << ".\n";
        return {};
    }

    int current = endRouteId;
    while (current != startRouteId && prevEdge.count(current))
    {
        path.push_back(prevEdge[current]);
        current = prevRoute[current];
    }

    std::reverse(path.begin(), path.end());
    return path;
}

bool areCoordinatesEqual(const Coordinate &a, const Coordinate &b, double tol = 1e-6)
{
    return std::abs(a.latitude - b.latitude) < tol &&
           std::abs(a.longitude - b.longitude) < tol;
}

std::vector<TFTFEdge> TFTFGraph::calculateRouteFromCoordinates(
    const Coordinate &startCoord,
    const Coordinate &endCoord,
    int hour)
{
    std::vector<int> startCandidates = getNearbyRoutes(startCoord, 300.0f);
    std::vector<int> endCandidates = getNearbyRoutes(endCoord, 300.0f);

    if (startCandidates.empty())
    {
        std::cerr << "No nearby routes found for the start coordinates.\n";
        return {};
    }
    else if (endCandidates.empty())
    {
        std::cerr << "No nearby routes found for the end coordinates.\n";
        return {};
    }

    std::vector<TFTFEdge> bestPath;
    double bestFare = std::numeric_limits<double>::infinity();
    int bestStartRouteId = -1;
    int bestEndRouteId = -1;
    bool bestPathSameRoute = false;

    for (int startId : startCandidates)
    {
        for (int endId : endCandidates)
        {
            std::vector<TFTFEdge> path;
            double fare = 0.0;

            if (startId == endId)
            {
                double segmentDistance = getActualSegmentDistance(startCoord, endCoord, routes.at(startId).path, routes.at(startId).isLoop);
                fare = BASE_FARE + (segmentDistance / 1000.0) * FARE_PER_KM;

                // Simulate one direct TFTFEdge for uniformity
                TFTFEdge edge;
                edge.destinationRoute = startId;
                edge.destinationRouteName = routes.at(startId).routeName;
                edge.transferCost = 0.0f;
                edge.entryCoord = projectOntoPath(startCoord, routes[startId].path);
                std::cout << "Entry coord: (" << edge.entryCoord.latitude << ", " << edge.entryCoord.longitude << ")\n";
                edge.exitCoord = projectOntoPath(endCoord, routes[startId].path);
                std::cout << "Exit coord: (" << edge.exitCoord.latitude << ", " << edge.exitCoord.longitude << ")\n";
                path.push_back(edge);
            }
            else
            {
                int projStart = getClosestIndex(routes.at(startId).path, startCoord);
                path = findMinFarePath(startId, endId, hour, projStart);
                if (path.empty())
                    continue;
                fare = calculateTotalFare(path, startCoord, endCoord);
            }

            if (fare < bestFare)
            {
                bestFare = fare;
                bestPath = path;
                bestStartRouteId = startId;
                bestEndRouteId = endId;
                bestPathSameRoute = (startId == endId);
            }
        }
    }

    if (bestPath.empty())
    {
        std::cerr << "No valid route found between coordinates.\n";
        return {};
    }

    // Handling transfer edge between different routes
    if (!bestPathSameRoute)
    {
        const Coordinate projectedStart = projectOntoPath(startCoord, routes[bestStartRouteId].path);
        const Coordinate firstEntry = bestPath.front().entryCoord;

        if (!areCoordinatesEqual(projectedStart, firstEntry))
        {
            TFTFEdge startEdge;
            startEdge.destinationRoute = bestPath.front().destinationRoute;
            startEdge.destinationRouteName = routes[bestStartRouteId].routeName;
            startEdge.transferCost = 0.0f;
            startEdge.entryCoord = projectedStart;
            startEdge.exitCoord = firstEntry;
            bestPath.insert(bestPath.begin(), startEdge);
        }

        // Add transfer edge from bestStartRouteId to bestEndRouteId
        TFTFEdge transferEdge;
        transferEdge.destinationRoute = bestEndRouteId;
        transferEdge.destinationRouteName = routes[bestEndRouteId].routeName;
        transferEdge.transferCost = 0.0f;                                                // Or add your transfer cost logic
        transferEdge.entryCoord = bestPath.back().exitCoord;                             // Exit coordinate from the first route
        transferEdge.exitCoord = projectOntoPath(endCoord, routes[bestEndRouteId].path); // Entry coordinate on the second route
        bestPath.push_back(transferEdge);
    }

    // Handle the final destination edge
    const Coordinate projectedEnd = projectOntoPath(endCoord, routes[bestEndRouteId].path);
    const Coordinate lastExit = bestPath.back().exitCoord;

    if (!areCoordinatesEqual(projectedEnd, lastExit))
    {
        TFTFEdge endEdge;
        endEdge.destinationRoute = -1; // Indicates final destination, no further route
        endEdge.destinationRouteName = "Destination";
        endEdge.transferCost = 0.0f;
        endEdge.entryCoord = bestPath.back().exitCoord;
        endEdge.exitCoord = endCoord;
        bestPath.push_back(endEdge);
    }

    // Output the route instructions
    std::cout << "\n==== Route Instructions ====\n";
    std::cout << "Best Path: " << bestPath.size() << " edges\n";

    for (size_t i = 0; i < bestPath.size(); ++i)
    {
        const TFTFEdge &edge = bestPath[i];

        if (i == 0)
        {
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "Start at coordinates: (" << edge.entryCoord.latitude << ", " << edge.entryCoord.longitude << ")\n";
        }

        if (edge.destinationRoute != -1)
        {
            std::cout << "Take route: " << edge.destinationRouteName << "\n";
            std::cout << "  Mount at: (" << edge.entryCoord.latitude << ", " << edge.entryCoord.longitude << ")\n";
            std::cout << "  Dismount at: (" << edge.exitCoord.latitude << ", " << edge.exitCoord.longitude << ")\n";
        }
        else
        {
            std::cout << "Walk to final destination at: (" << edge.exitCoord.latitude << ", " << edge.exitCoord.longitude << ")\n";
        }
    }
    std::cout << "=============================\n";

    float totalFare = calculateTotalFare(bestPath, startCoord, endCoord);
    std::cout << "Total fare: " << std::fixed << std::setprecision(2)
              << roundUpToNearest2_5(totalFare) << " pesos" << std::endl;

    return bestPath;
}
