#ifndef NODE_H
#define NODE_H

#include <vector>
#include <cmath>
#include <string>
#include "json.hpp"
#define M_PI 3.14159265358979323846

using json = nlohmann::json;

struct Node {
    double lat, lon;
    int routeId;        // The route this node belongs to
    size_t sequenceIdx; // Position in the route sequence
    bool isTransferPoint; // Whether this is a transfer point

    // Default constructor
    Node() : lat(0), lon(0), routeId(-1), sequenceIdx(0), isTransferPoint(false) {}

    // Constructor for regular route points
    Node(double latitude, double longitude, int route_id, size_t seq_idx) 
        : lat(latitude), lon(longitude), routeId(route_id), sequenceIdx(seq_idx), isTransferPoint(false) {}
    
    // Constructor for transfer points
    Node(double latitude, double longitude) 
        : lat(latitude), lon(longitude), routeId(-1), sequenceIdx(0), isTransferPoint(true) {}

    bool operator<(const Node &other) const {
        if (routeId != other.routeId) return routeId < other.routeId;
        if (isTransferPoint != other.isTransferPoint) return isTransferPoint < other.isTransferPoint;
        if (sequenceIdx != other.sequenceIdx) return sequenceIdx < other.sequenceIdx;
        return std::tie(lat, lon) < std::tie(other.lat, other.lon);
    }

    bool operator==(const Node &other) const {
        if (isTransferPoint && other.isTransferPoint) {
            // For transfer points, only compare coordinates within a small epsilon
            const double EPSILON = 1e-10;
            return std::abs(lat - other.lat) < EPSILON && 
                   std::abs(lon - other.lon) < EPSILON;
        }
        return routeId == other.routeId && 
               sequenceIdx == other.sequenceIdx && 
               lat == other.lat && 
               lon == other.lon;
    }

    bool operator!=(const Node &other) const {
        return !(*this == other);
    }
};

// Calculate distance between two nodes using haversine formula
double haversineNode(const Node &a, const Node &b);

// Convert GeoJSON to graph representation
void geojsonToNodeGraph(std::map<Node, std::vector<std::pair<Node, double>>> &graph, json file);

// Find transfer points between routes within specified distance
void findTransferPoints(std::map<Node, std::vector<std::pair<Node, double>>> &graph,
                       const std::vector<std::vector<Node>> &routes,
                       double maxTransferDistance);

#endif