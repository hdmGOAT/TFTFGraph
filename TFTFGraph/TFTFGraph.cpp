#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
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

    constexpr float BASE_FARE = 12.0f;
    constexpr float FARE_PER_4KM = 1.5f;
    constexpr float FREE_KM = 4.0f;

    float distanceInKm = distance / 1000.0f;

    // Base fare applies per route taken
    totalFare = BASE_FARE * takenRoutes.size();

    // Compute extra fare if distance exceeds free 4 km
    if (distanceInKm > FREE_KM)
    {
        float extraDistance = distanceInKm - FREE_KM;

        // Round up to the next 4km block
        int extraBlocks = static_cast<int>(std::ceil(extraDistance / 4.0f));

        totalFare += extraBlocks * FARE_PER_4KM;
    }

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

        // Skip artificial edges (routeId -1)
        if (route->routeId == -1)
            continue;

        std::vector<Coordinate> subPath;

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

    // Step 1: Populate spatial grid with densified paths
    for (auto &[routeID, route] : routes)
    {
        auto densePath = densifyPath(route.path, 25.0f);
        for (const auto &coord : densePath)
        {
            auto cell = hashCoordinate(coord, transferRangeMeters);
            spatialGrid[cell].emplace_back(routeID, coord);
        }
    }

    // Step 2: Store best transfers per (fromID, toID)
    struct TransferCandidate
    {
        Coordinate fromCoord;
        Coordinate toCoord;
        float dist;
    };

    std::map<std::pair<int, int>, TransferCandidate> bestTransfers;

    // Step 3: Find closest coord pairs within range
    for (const auto &[fromID, fromRoute] : routes)
    {
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

                        auto key = std::make_pair(fromID, toID);
                        if (!bestTransfers.count(key) || dist < bestTransfers[key].dist)
                        {
                            bestTransfers[key] = {fromCoord, toCoord, dist};
                        }
                    }
                }
            }
        }
    }

    // Step 4: Create edges using closest coord pair
    for (const auto &[key, candidate] : bestTransfers)
    {
        int fromID = key.first;
        int toID = key.second;

        TransferZone zone1 = {
            .routeId = fromID,
            .start = routes[fromID].path.front(),
            .end = routes[fromID].path.back(),
            .closestCoord = candidate.fromCoord};

        TransferZone zone2 = {
            .routeId = toID,
            .start = routes[toID].path.front(),
            .end = routes[toID].path.back(),
            .closestCoord = candidate.toCoord};

        addEdge(fromID, toID, zone1, zone2, candidate.dist);
    }
}

std::vector<int> TFTFGraph::getNearbyRoutes(const Coordinate &coord, float maxDistanceMeters)
{
    // Map to store minimum distance to each route
    std::map<int, float> routeMinDistances;

    // Find minimum distance to each route
    for (const auto &[routeId, route] : routes)
    {
        auto densePath = densifyPath(route.path, 25.0f);
        float minDist = std::numeric_limits<float>::infinity();

        for (const auto &point : densePath)
        {
            float dist = haversine(coord, point);
            if (dist <= maxDistanceMeters && dist < minDist)
            {
                minDist = dist;
                routeMinDistances[routeId] = dist;
            }
        }
    }

    // Sort routes by distance
    std::vector<std::pair<int, float>> sortedRoutes;
    for (const auto &[routeId, dist] : routeMinDistances)
    {
        sortedRoutes.push_back({routeId, dist});
    }
    std::sort(sortedRoutes.begin(), sortedRoutes.end(),
              [](const auto &a, const auto &b)
              { return a.second < b.second; });

    // Take only the 3 closest routes
    std::vector<int> nearby;
    for (size_t i = 0; i < std::min(size_t(3), sortedRoutes.size()); ++i)
    {
        nearby.push_back(sortedRoutes[i].first);
    }

    return nearby;
}

std::vector<TFTFEdge> TFTFGraph::findMinFarePath(
    int startRouteId,
    int endRouteId,
    int projectedStartIdx)
{
    struct TransferNode
    {
        int routeId;
        int transfers;
        Coordinate currentPos;
        std::vector<TFTFEdge> path;
        float fareSoFar;
    };

    std::queue<TransferNode> q;
    std::unordered_map<int, int> minTransfers; // routeId -> fewest transfers
    std::vector<std::vector<TFTFEdge>> candidatePaths;

    Coordinate startPos = routes.at(startRouteId).path[projectedStartIdx];
    q.push({startRouteId, 0, startPos, {}, 0.0f});
    minTransfers[startRouteId] = 0;

    int minTransferCount = INT_MAX;

    while (!q.empty())
    {
        TransferNode curr = q.front();
        q.pop();

        if (curr.transfers > minTransferCount)
            continue; // Stop exploring once deeper than the shortest transfer path

        if (curr.routeId == endRouteId)
        {
            minTransferCount = curr.transfers;
            candidatePaths.push_back(curr.path); // Store this path
            continue;
        }

        for (const TFTFEdge &edge : routes.at(curr.routeId).edges)
        {
            int nextRouteId = edge.transferZone2.routeId;

            // Prevent backward transfers on the first leg if not a loop
            if (curr.routeId == startRouteId && !routes.at(curr.routeId).isLoop)
            {
                int transferIdx = getClosestIndex(routes.at(curr.routeId).path, edge.transferZone1.closestCoord);
                if (transferIdx < projectedStartIdx)
                    continue;
            }

            if (minTransfers.count(nextRouteId) && curr.transfers + 1 > minTransfers[nextRouteId])
                continue; // We've already reached this with fewer transfers

            // Track minimal transfers
            minTransfers[nextRouteId] = curr.transfers + 1;

            // Create new path and accumulate fare
            std::vector<TFTFEdge> newPath = curr.path;
            newPath.push_back(edge);
            float newFare = curr.fareSoFar + edge.transferCost;

            q.push({nextRouteId, curr.transfers + 1, edge.transferZone2.closestCoord, newPath, newFare});
        }
    }

    // Select the candidate path with the lowest fare
    if (candidatePaths.empty())
    {
        std::cout << "No path found from route " << startRouteId << " to route " << endRouteId << ".\n";
        return {};
    }

    std::vector<TFTFEdge> bestPath = candidatePaths[0];
    float bestFare = 0.0f;
    // for (const TFTFEdge &e : bestPath)
    //     bestFare += e.transferCost;

    // for (const auto &path : candidatePaths)
    // {
    //     float totalFare = 0.0f;
    //     for (const TFTFEdge &e : path)
    //         totalFare += e.transferCost;

    //     if (totalFare < bestFare)
    //     {
    //         bestFare = totalFare;
    //         bestPath = path;
    //     }
    // }

    return bestPath;
}

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

std::string generateRandomColor()
{
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
    json geojson = generateRoutePathGeoJSON(instructions, startCoord, endCoord);

    std::ofstream file(filename);
    if (file.is_open())
    {
        file << std::setw(2) << geojson << std::endl;
        file.close();
        std::cout << "GeoJSON written to: " << filename << "\n";
    }
    else
    {
        std::cerr << "Failed to open file for writing: " << filename << "\n";
    }
}

double TFTFGraph::calculateFareFromInstructions(const std::vector<RoutePathInstruction> &routeInstructions)
{
    constexpr float BASE_FARE = 12.0f;
    constexpr float FARE_PER_4KM = 1.5f;
    constexpr float FREE_KM = 4.0f;

    double totalFare = 0.0;

    for (const auto &instr : routeInstructions)
    {
        double routeDistance = 0.0;
        const auto &path = instr.path;

        // Calculate distance for this route
        for (size_t i = 1; i < path.size(); ++i)
        {
            routeDistance += haversine(path[i - 1], path[i]);
        }

        // Calculate fare for this route
        double routeFare = BASE_FARE;
        float distanceInKm = routeDistance / 1000.0f;

        if (distanceInKm > FREE_KM)
        {
            float extraDistance = distanceInKm - FREE_KM;
            int extraBlocks = static_cast<int>(std::ceil(extraDistance / 4.0f));
            routeFare += extraBlocks * FARE_PER_4KM;
        }

        totalFare += routeFare;
    }

    return totalFare;
}

std::vector<TFTFEdge> TFTFGraph::calculateRouteFromCoordinates(
    const Coordinate &startCoord,
    const Coordinate &endCoord)
{
    auto start = std::chrono::high_resolution_clock::now();

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
                    .closestCoord = projectOntoPath(startCoord, routes[startId].path)};

                TransferZone endZone = {
                    .routeId = endId,
                    .start = routes[endId].path.front(),
                    .end = routes[endId].path.back(),
                    .closestCoord = projectOntoPath(endCoord, routes[endId].path)};

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

    // Add artificial start edge if not same route
    if (!bestPathSameRoute)
    {
        TransferZone startZone = {
            .routeId = bestStartRouteId,
            .start = routes[bestStartRouteId].path.front(),
            .end = routes[bestStartRouteId].path.back(),
            .closestCoord = projectOntoPath(startCoord, routes[bestStartRouteId].path)};

        TransferZone firstZone = bestPath.front().transferZone1;

        TFTFEdge startEdge = {
            .transferCost = 0.0f,
            .transferZone1 = startZone,
            .transferZone2 = firstZone};

        bestPath.insert(bestPath.begin(), startEdge);
    }

    // Add artificial end edge
    TransferZone lastZone = bestPath.back().transferZone2;

    TransferZone endZone = {
        .routeId = -1,
        .start = endCoord,
        .end = endCoord,
        .closestCoord = endCoord};

    TFTFEdge endEdge = {
        .transferCost = 0.0f,
        .transferZone1 = lastZone,
        .transferZone2 = endZone};

    bestPath.push_back(endEdge);

    auto end = std::chrono::high_resolution_clock::now();

    return bestPath;
}

std::unordered_map<int, RouteNode> &TFTFGraph::getRoutes()
{
    return routes;
}

void TFTFGraph::getGraphDetails() const
{
    size_t routeNodeCount = routes.size();
    size_t edgeCount = 0;
    size_t totalCoordinates = 0;

    for (const auto &[routeId, routeNode] : routes)
    {
        edgeCount += routeNode.edges.size();
        totalCoordinates += routeNode.path.size();
    }

    std::cout << "Graph Details:\n";
    std::cout << "Total Route Nodes: " << routeNodeCount << "\n";
    std::cout << "Total Edges: " << edgeCount << "\n";
    std::cout << "Total Coordinates: " << totalCoordinates << "\n";
}

json TFTFGraph::toJson() const
{
    json j;
    j["routes"] = json::array();

    for (const auto &[routeId, route] : routes)
    {
        json routeJson;
        routeJson["id"] = route.routeId;
        routeJson["name"] = route.routeName;
        routeJson["is_loop"] = route.isLoop;

        // Serialize coordinates
        routeJson["path"] = json::array();
        for (const auto &coord : route.path)
        {
            routeJson["path"].push_back({{"lat", coord.latitude}, {"lon", coord.longitude}});
        }

        // Serialize edges
        routeJson["edges"] = json::array();
        for (const auto &edge : route.edges)
        {
            json edgeJson;
            edgeJson["transfer_cost"] = edge.transferCost;

            auto serializeZone = [](const TransferZone &z) -> json
            {
                return {
                    {"route_id", z.routeId},
                    {"start", {{"lat", z.start.latitude}, {"lon", z.start.longitude}}},
                    {"end", {{"lat", z.end.latitude}, {"lon", z.end.longitude}}},
                    {"closest", {{"lat", z.closestCoord.latitude}, {"lon", z.closestCoord.longitude}}}};
            };

            edgeJson["zone1"] = serializeZone(edge.transferZone1);
            edgeJson["zone2"] = serializeZone(edge.transferZone2);

            routeJson["edges"].push_back(edgeJson);
        }

        j["routes"].push_back(routeJson);
    }

    return j;
}

void saveGraphToDisk(const TFTFGraph &graph, const std::string &filename)
{
    std::ofstream outFile(filename);
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open file for writing: " << filename << "\n";
        return;
    }

    outFile << std::setw(2) << graph.toJson() << std::endl;
    outFile.close();
    std::cout << "Graph saved to " << filename << "\n";
}

TFTFGraph loadGraphFromDisk(const std::string &filename)
{
    TFTFGraph graph;
    std::ifstream inFile(filename);
    if (!inFile.is_open())
    {
        std::cerr << "Failed to open file: " << filename << "\n";
        return graph;
    }

    json j;
    inFile >> j;
    inFile.close();

    for (const auto &routeJson : j["routes"])
    {
        int routeId = routeJson["id"];
        std::string routeName = routeJson["name"];
        bool isLoop = routeJson["is_loop"];

        std::vector<Coordinate> path;
        for (const auto &coordJson : routeJson["path"])
        {
            path.emplace_back(Coordinate{
                coordJson["lat"].get<double>(),
                coordJson["lon"].get<double>()});
        }

        graph.addRoute(routeId, routeName);
        graph.setRoutePath(routeId, path);

        // Use getRoutes() to access and update the route
        auto &routesMap = graph.getRoutes();
        routesMap.at(routeId).isLoop = isLoop;

        // Deserialize edges
        for (const auto &edgeJson : routeJson["edges"])
        {
            TFTFEdge edge;
            edge.transferCost = edgeJson["transfer_cost"];

            auto parseZone = [](const json &z) -> TransferZone
            {
                return {
                    .routeId = z["route_id"],
                    .start = {z["start"]["lat"], z["start"]["lon"]},
                    .end = {z["end"]["lat"], z["end"]["lon"]},
                    .closestCoord = {z["closest"]["lat"], z["closest"]["lon"]}};
            };

            edge.transferZone1 = parseZone(edgeJson["zone1"]);
            edge.transferZone2 = parseZone(edgeJson["zone2"]);

            routesMap.at(routeId).edges.push_back(edge);
        }
    }

    std::cout << "Graph loaded from " << filename << "\n";
    return graph;
}
json generateRoutePathGeoJSON(
    const std::vector<RoutePathInstruction> &instructions,
    const Coordinate &startCoord,
    const Coordinate &endCoord)
{
    json featureCollection;
    featureCollection["type"] = "FeatureCollection";
    featureCollection["features"] = json::array();

    for (const auto &instr : instructions)
    {
        if (instr.path.empty())
            continue;

        json lineFeature;
        lineFeature["type"] = "Feature";
        lineFeature["geometry"]["type"] = "LineString";

        for (const auto &coord : instr.path)
        {
            lineFeature["geometry"]["coordinates"].push_back({coord.longitude, coord.latitude});
        }

        lineFeature["properties"] = {
            {"route_name", instr.routeName},
            {"route_id", instr.routeId}};

        featureCollection["features"].push_back(lineFeature);
    }

    // Start point
    featureCollection["features"].push_back({{"type", "Feature"},
                                             {"geometry", {{"type", "Point"}, {"coordinates", {startCoord.longitude, startCoord.latitude}}}},
                                             {"properties", {{"marker", "start"}}}});

    // End point
    featureCollection["features"].push_back({{"type", "Feature"},
                                             {"geometry", {{"type", "Point"}, {"coordinates", {endCoord.longitude, endCoord.latitude}}}},
                                             {"properties", {{"marker", "end"}}}});

    return featureCollection;
}
