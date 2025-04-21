#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <limits>
#include <fstream>
#include "TFTFGraph.h"
#include <queue>
#include <unordered_set>
#include <limits>
#include <algorithm>
#include "Helpers/helpers.h"
#include "../json.hpp"

bool operator==(const Coordinate &lhs, const Coordinate &rhs)
{
    return lhs.latitude == rhs.latitude && lhs.longitude == rhs.longitude;
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
    int routeId1,
    int routeId2,
    TransferZone route1,
    TransferZone route2,
    float transferCost)
{

    TFTFEdge edge;
    edge.transferCost = transferCost;
    edge.transferZone1 = route1;
    edge.transferZone2 = route2;
    routes[routeId1].edges.push_back(edge);
}

std::vector<const RouteNode *> TFTFGraph::extractTraversedRouteNodes(
    const std::vector<TFTFEdge> &path) const
{

    std::vector<const RouteNode *> routeNodes;

    if (path.empty())
        return routeNodes;

    // check if the origin route of the first edge exists in the routes map
    if (routes.count(path.front().transferZone1.routeId))
    {
        // add the origin route node of the first edge to the result
        routeNodes.push_back(&routes.at(path.front().transferZone1.routeId));
    }
    else if (routes.count(path.front().transferZone1.routeId))
    {
        // add the origin route node of the first edge to the result
        routeNodes.push_back(&routes.at(path.front().transferZone1.routeId));
    }

    // iterate through the edges in the given path
    for (const auto &edge : path)
    {
        // check if the destination route of the current edge exists in the routes map
        if (routes.count(edge.transferZone2.routeId))
        {
            // add the destination route node to the result
            routeNodes.push_back(&routes.at(edge.transferZone2.routeId));
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
        TFTFEdge nextEdge = path[i];
        if (i + 1 < path.size())
        {
            nextEdge = path[i + 1];
        }

        // if it's the first edge, calculate distance from start to the first edge's exit
        if (i == 0)
        {
            int startIdx = getClosestIndex(takenRoutes[0]->path, projectedStart);
            int endIdx = getClosestIndex(takenRoutes[0]->path, edge.transferZone2.closestCoord);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes[0]->isLoop && startIdx > endIdx)
            {
                continue;
            }

            distance += getActualSegmentDistance(projectedStart, edge.transferZone2.closestCoord, takenRoutes[0]->path, takenRoutes[0]->isLoop);
        }
        // if it's the last edge, calculate distance from the last edge's entry to the end
        else if (i == path.size() - 1)
        {
            int startIdx = getClosestIndex(takenRoutes.back()->path, edge.transferZone1.closestCoord);
            int endIdx = getClosestIndex(takenRoutes.back()->path, projectedEnd);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes.back()->isLoop && startIdx > endIdx)
            {
                continue;
            }
            distance += getActualSegmentDistance(edge.transferZone1.closestCoord, projectedEnd, takenRoutes.back()->path, takenRoutes.back()->isLoop);
        }
        // for all other edges, calculate distance between entry and exit coordinates
        else
        {
            int startIdx = getClosestIndex(takenRoutes[i]->path, edge.transferZone1.closestCoord);
            int endIdx = getClosestIndex(takenRoutes[i]->path, nextEdge.transferZone2.closestCoord);

            // skip if segment would wrap around a non-loop route
            if (!takenRoutes[i]->isLoop && startIdx > endIdx)
            {
                continue;
            }

            distance += getActualSegmentDistance(edge.transferZone1.closestCoord, nextEdge.transferZone2.closestCoord, takenRoutes[i]->path, takenRoutes[i]->isLoop);
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

std::vector<RoutePathInstruction> TFTFGraph::constructRoutePathInstructions(
    const std::vector<TFTFEdge> &path) const
{
    std::vector<RoutePathInstruction> routeInstructions;
    std::vector<const RouteNode *> takenRoutes = extractTraversedRouteNodes(path);

    if (takenRoutes.empty())
        return {};

    for (size_t i = 0; i < path.size(); ++i)
    {
        const TFTFEdge &edge = path[i];
        const RouteNode *route = takenRoutes[i];
        std::vector<Coordinate> subPath;

        std::cout << "Route ID: " << route->routeId << "\n";

        if (i == 0)
        {
            // First edge: from projected start to exit
            int startIdx = getClosestIndex(route->path, edge.transferZone1.closestCoord);
            int endIxd = getClosestIndex(route->path, edge.transferZone2.closestCoord);
            subPath = getFullSegmentPath(route->path, startIdx, endIxd, route->isLoop);
        }
        else if (i == path.size() - 1)
        {
            subPath = getShortestSegmentPath(route->path, edge.transferZone1.closestCoord, edge.transferZone2.closestCoord, route->isLoop);
        }
        else
        {
            // Intermediate: entry to exit of next
            const TFTFEdge &nextEdge = path[i + 1];
            int startIdx = getClosestIndex(route->path, edge.transferZone2.closestCoord);
            int endIdx = getClosestIndex(route->path, nextEdge.transferZone1.closestCoord);
            subPath = getShortestSegmentPath(route->path, edge.transferZone2.closestCoord, nextEdge.transferZone1.closestCoord, route->isLoop);
        }

        // Avoid duplicate join point
        if (!routeInstructions.empty() && !subPath.empty() &&
            routeInstructions.back().path.back() == subPath.front())
        {
            subPath.erase(subPath.begin());
        }

        if (subPath.size() <= 1)
        {
            continue;
        }

        routeInstructions.push_back({route->routeId,
                                     route->routeName,
                                     subPath});
    }

    return routeInstructions;
}
void TFTFGraph::createTransfersFromCoordinates(float transferRangeMeters)
{
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, Coordinate>>, PairHash> spatialGrid;

    // Step 1: Populate the spatial grid with densified paths
    for (auto &[routeID, route] : routes)
    {
        auto densePath = densifyPath(route.path, 25.0f);
        for (const auto &coord : densePath)
        {
            auto cell = hashCoordinate(coord, transferRangeMeters);
            spatialGrid[cell].emplace_back(routeID, coord);
        }
    }

    // Step 2: Check possible transfer zones
    for (const auto &[fromID, fromRoute] : routes)
    {
        std::unordered_map<int, std::pair<Coordinate, float>> lastTransferPerRoute;

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

                        auto &entry = lastTransferPerRoute[toID];
                        float fromLastDist = (entry.second > 0.0f) ? haversine(entry.first, fromCoord) : std::numeric_limits<float>::max();

                        if (fromLastDist >= 500.0f || entry.second == 0.0f)
                        {
                            // Build TransferZones
                            TransferZone zone1 = {
                                .routeId = fromID,
                                .start = fromRoute.path.front(),
                                .end = fromRoute.path.back(),
                                .closestCoord = fromCoord
                            };

                            TransferZone zone2 = {
                                .routeId = toID,
                                .start = routes[toID].path.front(),
                                .end = routes[toID].path.back(),
                                .closestCoord = toCoord
                            };

                            TFTFEdge edge = {
                                .transferCost = dist,
                                .transferZone1 = zone1,
                                .transferZone2 = zone2
                            };

                            entry = {fromCoord, dist};

                           addEdge(fromID, toID, zone1, zone2, dist);
                        }
                        else if (dist < entry.second)
                        {
                            // Optional improvement logic (deduplication not shown)
                            entry = {fromCoord, dist};
                        }
                    }
                }
            }
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

std::string makeEdgeKey(int routeId, int entryIdx, int exitIdx) {
    return std::to_string(routeId) + "_" + std::to_string(entryIdx) + "_" + std::to_string(exitIdx);
}

std::vector<TFTFEdge> TFTFGraph::findMinFarePath(
    int startRouteId,
    int endRouteId,
    int projectedStartIdx)
{
    struct TransferState
    {
        float totalDistance;
        int prevRouteId;
    };

    struct TransferNode
    {
        int routeId;
        float distanceSoFar;

        bool operator>(const TransferNode &other) const
        {
            return distanceSoFar > other.distanceSoFar;
        }
    };

    std::unordered_map<int, TransferState> best;
    std::unordered_map<int, TFTFEdge> prevEdge;
    std::unordered_map<int, int> prevRoute;

    std::priority_queue<TransferNode, std::vector<TransferNode>, std::greater<TransferNode>> pq;

    best[startRouteId] = {0.0f, -1};
    pq.push({startRouteId, 0.0f});
    Coordinate currPos = routes.at(startRouteId).path[projectedStartIdx];

    while (!pq.empty())
    {
        TransferNode curr = pq.top();
        pq.pop();

        int currRouteId = curr.routeId;
        float distanceSoFar = curr.distanceSoFar;

        for (const TFTFEdge &edge : routes.at(currRouteId).edges)
        {
            int nextRouteId = edge.transferZone2.routeId;

            // Only allow forward transfer if it's the first leg and the route is non-looping
            if (currRouteId == startRouteId && !routes.at(currRouteId).isLoop)
            {
                float distToTransferCoord = getActualSegmentDistance(
                    currPos,
                    edge.transferZone1.closestCoord,
                    routes.at(currRouteId).path,
                    routes.at(currRouteId).isLoop
                );

                // Skip if transferCoord is before projectedStartIdx
                int idx = getClosestIndex(routes.at(currRouteId).path, edge.transferZone1.closestCoord);
                if (idx < projectedStartIdx)
                    continue;

                distanceSoFar += distToTransferCoord; // Add to running distance
            }

            float newDistance = distanceSoFar + edge.transferCost;

            if (!best.count(nextRouteId) || newDistance < best[nextRouteId].totalDistance)
            {
                best[nextRouteId] = {newDistance, currRouteId};
                pq.push({nextRouteId, newDistance});
                prevEdge[nextRouteId] = edge;
                prevRoute[nextRouteId] = currRouteId;
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



struct PathState {
    int routeId;
    int index;
    float totalDistance;
    int transfers;
    std::vector<TFTFEdge> path;
};



void printRoutePathInstructions(const std::vector<RoutePathInstruction> &instructions)
{
    for (const auto &instr : instructions)
    {
        std::cout << "Route ID: " << instr.routeId << "\n";
        std::cout << "Route Name: " << instr.routeName << "\n";
        std::cout << "Path Coordinates (" << instr.path.size() << " points):\n";

        for (const auto &coord : instr.path)
        {
            std::cout << std::fixed << std::setprecision(6); // or however many decimal places you want
            std::cout << "  [" << coord.longitude << ", " << coord.latitude << "],\n";
        }
        std::cout << "-----------------------------\n";
    }
}

#include "../json.hpp"
#include <random>
using json = nlohmann::json;

std::string generateRandomColor() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 255);

    std::ostringstream color;
    color << "#" << std::hex << std::setfill('0') 
          << std::setw(2) << dist(gen)
          << std::setw(2) << dist(gen)
          << std::setw(2) << dist(gen);

    return color.str();
}
void writeRoutePathInstructionsToGeoJSON(
    const std::vector<RoutePathInstruction> &instructions,
    const std::string &filename,
    const Coordinate &startCoord,
    const Coordinate &endCoord)
{
    json featureCollection;
    featureCollection["type"] = "FeatureCollection";
    featureCollection["features"] = json::array();

    for (const auto &instr : instructions)
    {
        if (instr.path.empty()) continue;

        std::string color = generateRandomColor();

        // Route LineString
        json lineFeature;
        lineFeature["type"] = "Feature";
        lineFeature["geometry"]["type"] = "LineString";

        for (const auto &coord : instr.path)
        {
            lineFeature["geometry"]["coordinates"].push_back({coord.longitude, coord.latitude, 0});
        }

        lineFeature["properties"]["name"] = instr.routeName;
        lineFeature["properties"]["route_id"] = instr.routeId;
        lineFeature["properties"]["stroke"] = color;

        featureCollection["features"].push_back(lineFeature);
    }

    // Start Point Feature
    json startFeature;
    startFeature["type"] = "Feature";
    startFeature["geometry"]["type"] = "Point";
    startFeature["geometry"]["coordinates"] = {startCoord.longitude, startCoord.latitude, 0};
    startFeature["properties"]["marker-symbol"] = "circle";
    startFeature["properties"]["marker-color"] = "#00FF00"; // green
    startFeature["properties"]["marker-size"] = "medium";
    startFeature["properties"]["point_type"] = "start";

    featureCollection["features"].push_back(startFeature);

    // End Point Feature
    json endFeature;
    endFeature["type"] = "Feature";
    endFeature["geometry"]["type"] = "Point";
    endFeature["geometry"]["coordinates"] = {endCoord.longitude, endCoord.latitude, 0};
    endFeature["properties"]["marker-symbol"] = "circle";
    endFeature["properties"]["marker-color"] = "#FF0000"; // red
    endFeature["properties"]["marker-size"] = "medium";
    endFeature["properties"]["point_type"] = "end";

    featureCollection["features"].push_back(endFeature);

    // Write to file
    std::ofstream file(filename);
    if (file.is_open())
    {
        file << std::setw(2) << featureCollection << std::endl;
        file.close();
        std::cout << "GeoJSON written to: " << filename << "\n";
    }
    else
    {
        std::cerr << "Failed to open file for writing: " << filename << "\n";
    }
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

                TransferZone startZone = {
                    .routeId = startId,
                    .start = routes[startId].path.front(),
                    .end = routes[startId].path.back(),
                    .closestCoord = projectOntoPath(startCoord, routes[startId].path)
                };

                TransferZone endZone = {
                    .routeId = endId,
                    .start = routes[endId].path.front(),
                    .end = routes[endId].path.back(),
                    .closestCoord = projectOntoPath(endCoord, routes[endId].path)
                };

                TFTFEdge edge;
                edge.transferCost = segmentDistance;
                edge.transferZone1 = startZone;
                edge.transferZone2 = endZone;
                path.push_back(edge);
            }
            else
            {
                int projStart = getClosestIndex(routes.at(startId).path, startCoord);
                path = findMinFarePath(startId, endId, projStart);
                if (path.empty()) continue;

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

    // Add artificial start edge if not same route
    if (!bestPathSameRoute)
    {
        TransferZone startZone = {
            .routeId = bestStartRouteId,
            .start = routes[bestStartRouteId].path.front(),
            .end = routes[bestStartRouteId].path.back(),
            .closestCoord = projectOntoPath(startCoord, routes[bestStartRouteId].path)
        };

        TransferZone firstZone = bestPath.front().transferZone1;

        TFTFEdge startEdge = {
            .transferCost = 0.0f,
            .transferZone1 = startZone,
            .transferZone2 = firstZone
        };

        bestPath.insert(bestPath.begin(), startEdge);
    }

    // Add artificial end edge
    TransferZone lastZone = bestPath.back().transferZone2;

    TransferZone endZone = {
        .routeId = -1,
        .start = endCoord,
        .end = endCoord,
        .closestCoord = endCoord
    };

    TFTFEdge endEdge = {
        .transferCost = 0.0f,
        .transferZone1 = lastZone,
        .transferZone2 = endZone
    };

    bestPath.push_back(endEdge);

    // ==== PRINT INSTRUCTIONS ====
    std::cout << "\n==== Route Instructions ====\n";
    for (size_t i = 0; i < bestPath.size(); ++i)
    {
        const TFTFEdge &edge = bestPath[i];

        std::cout << std::fixed << std::setprecision(6);

        if (i == 0)
        {
            std::cout << "Start at coordinates: ("
                      << edge.transferZone1.closestCoord.latitude << ", "
                      << edge.transferZone1.closestCoord.longitude << ")\n";
        }

        if (edge.transferZone2.routeId != -1)
        {
            std::cout << "Take route: " << routes[edge.transferZone2.routeId].routeName << "\n";
            std::cout << "  Mount at: ("
                      << edge.transferZone1.closestCoord.latitude << ", "
                      << edge.transferZone1.closestCoord.longitude << ")\n";

            // Dismount is at the next edge's mount, or the final point
            if (i + 1 < bestPath.size())
            {
                const TFTFEdge &nextEdge = bestPath[i + 1];
                std::cout << "  Dismount at: ("
                          << nextEdge.transferZone2.closestCoord.latitude << ", "
                          << nextEdge.transferZone1.closestCoord.longitude << ")\n";
            }
            else
            {
                std::cout << "  Dismount: (unknown – no next segment)\n";
            }
        }
        else
        {
            std::cout << "Walk to final destination at: ("
                      << edge.transferZone2.closestCoord.latitude << ", "
                      << edge.transferZone2.closestCoord.longitude << ")\n";
        }
    }

    std::cout << "=============================\n";
    std::cout << "Total fare: " << std::fixed << std::setprecision(2)
              << roundUpToNearest2_5(bestFare) << " pesos\n";
    std::cout << "=============================\n";

    std::vector<RoutePathInstruction> routeInstructions = constructRoutePathInstructions(bestPath);
    writeRoutePathInstructionsToGeoJSON(routeInstructions, "route_path.geojson", startCoord, endCoord);

    return bestPath;
}
