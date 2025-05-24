#include "node.h"
#include "json.hpp"
#include <unordered_map>
#include <set>
#include <chrono>
#include <iostream>

using json = nlohmann::json;

// Hash function for coordinate grid cells
struct PairHash {
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

// Function to hash a coordinate to a grid cell
std::pair<int, int> hashCoordinate(double lat, double lon, double cellSizeMeters) {
    double latMeters = 111320.0;                                    // meters per degree latitude
    double lonMeters = 111320.0 * std::cos(lat * M_PI / 180.0);    // meters per degree longitude

    int x = static_cast<int>(lon * lonMeters / cellSizeMeters);
    int y = static_cast<int>(lat * latMeters / cellSizeMeters);
    return {x, y};
}

double haversineNode(const Node &a, const Node &b)
{
    const double R = 6371e3;
    double lat1 = a.lat * M_PI / 180;
    double lat2 = b.lat * M_PI / 180;
    double dlat = (b.lat - a.lat) * M_PI / 180;
    double dlon = (b.lon - a.lon) * M_PI / 180;

    double h = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);

    return R * 2 * atan2(sqrt(h), sqrt(1 - h));
}

void findTransferPoints(std::map<Node, std::vector<std::pair<Node, double>>> &graph,
                       const std::vector<std::vector<Node>> &routes,
                       double maxTransferDistance) {
    
    // Create spatial grid
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, Node>>, PairHash> spatialGrid;

    // Step 1: Populate spatial grid
    for (size_t routeId = 0; routeId < routes.size(); ++routeId) {
        for (const auto& node : routes[routeId]) {
            auto cell = hashCoordinate(node.lat, node.lon, maxTransferDistance);
            spatialGrid[cell].emplace_back(routeId, node);
        }
    }

    // Step 2: Store best transfers per (fromID, toID)
    struct TransferCandidate {
        Node fromNode;
        Node toNode;
        double dist;
    };
    std::map<std::pair<int, int>, TransferCandidate> bestTransfers;

    // Step 3: Find closest node pairs within range
    for (size_t fromRouteId = 0; fromRouteId < routes.size(); ++fromRouteId) {
        for (const auto& fromNode : routes[fromRouteId]) {
            auto baseCell = hashCoordinate(fromNode.lat, fromNode.lon, maxTransferDistance);

            // Check neighboring cells
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    std::pair<int, int> neighborCell = {baseCell.first + dx, baseCell.second + dy};
                    auto it = spatialGrid.find(neighborCell);
                    if (it == spatialGrid.end()) continue;

                    for (const auto& [toRouteId, toNode] : it->second) {
                        if (toRouteId == fromRouteId) continue;

                        double dist = haversineNode(fromNode, toNode);
                        if (dist > maxTransferDistance) continue;

                        auto key = std::make_pair(fromRouteId, toRouteId);
                        if (!bestTransfers.count(key) || dist < bestTransfers[key].dist) {
                            bestTransfers[key] = {fromNode, toNode, dist};
                        }
                    }
                }
            }
        }
    }

    // Step 4: Create transfer points and edges for best transfers
    for (const auto& [key, candidate] : bestTransfers) {
        // Create transfer point at fromNode's location
        Node transfer(candidate.fromNode.lat, candidate.fromNode.lon);
        
        // Add edges in both directions with transfer cost
        graph[candidate.fromNode].emplace_back(transfer, candidate.dist);
        graph[transfer].emplace_back(candidate.toNode, candidate.dist);
        
        // Add reverse transfer
        graph[candidate.toNode].emplace_back(transfer, candidate.dist);
        graph[transfer].emplace_back(candidate.fromNode, candidate.dist);
    }
}

void geojsonToNodeGraph(std::map<Node, std::vector<std::pair<Node, double>>> &graph, json file) {
    std::vector<std::vector<Node>> routes;
    int routeId = 0;

    // First pass: Create nodes and route-following edges
    for (auto &feature : file["features"]) {
        if (feature["geometry"]["type"] != "LineString")
            continue;

        auto coords = feature["geometry"]["coordinates"];
        std::vector<Node> routeNodes;

        // Create nodes for this route
        for (size_t i = 0; i < coords.size(); ++i) {
            Node node(coords[i][1], coords[i][0], routeId, i);
            routeNodes.push_back(node);

            // Create edge to next node in sequence (if not last node)
            if (i + 1 < coords.size()) {
                Node nextNode(coords[i + 1][1], coords[i + 1][0], routeId, i + 1);
                double dist = haversineNode(node, nextNode);
                graph[node].emplace_back(nextNode, dist);
            }

            // If this is a loop route (first == last coordinate), add edge from last to first
            if (i == coords.size() - 1 && coords.size() > 1) {
                if (coords[0][0] == coords[i][0] && coords[0][1] == coords[i][1]) {
                    Node firstNode(coords[0][1], coords[0][0], routeId, 0);
                    double dist = haversineNode(node, firstNode);
                    graph[node].emplace_back(firstNode, dist);
                }
            }
        }

        routes.push_back(routeNodes);
        routeId++;
    }

    // Second pass: Find and add transfer points
    findTransferPoints(graph, routes, 300.5); // Using the same transfer range as TFTFGraph
}
