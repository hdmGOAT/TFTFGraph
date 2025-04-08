#define _USE_MATH_DEFINES
#include <cmath>
#include "helpers.h"
#include "../TFTFGraph.h"



float haversine(const Coordinate& a, const Coordinate& b)  {
    auto deg2rad = [](double deg) {
        return deg * M_PI / 180.0;
    };

    double lat1 = deg2rad(a.latitude);
    double lon1 = deg2rad(a.longitude);
    double lat2 = deg2rad(b.latitude);
    double lon2 = deg2rad(b.longitude);

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double h = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(dlon / 2) * std::sin(dlon / 2);

    double c = 2 * std::atan2(std::sqrt(h), std::sqrt(1 - h));
    return EARTH_RADIUS_METERS * c;
}

int closestCoordinateIndex(const std::vector<Coordinate>& path, const Coordinate& target) {
    int bestIdx = 0;
    float bestDist = std::numeric_limits<float>::infinity();
    for (size_t i = 0; i < path.size(); ++i) {
        float dist = haversine(path[i], target);
        if (dist < bestDist) {
            bestDist = dist;
            bestIdx = i;
        }
    }
    return bestIdx;
}

float computeRouteDistance(const std::vector<Coordinate>& path, const Coordinate& start, const Coordinate& end) {
    if (path.empty()) return 0.0f;

    int closestStart = -1;
    int closestEnd = -1;
    float minStartDist = std::numeric_limits<float>::infinity();
    float minEndDist = std::numeric_limits<float>::infinity();

    for (int i = 0; i < path.size(); ++i) {
        float dStart = haversine(path[i], start);
        float dEnd = haversine(path[i], end);
        if (dStart < minStartDist) {
            minStartDist = dStart;
            closestStart = i;
        }
        if (dEnd < minEndDist) {
            minEndDist = dEnd;
            closestEnd = i;
        }
    }

    if (closestStart > closestEnd) std::swap(closestStart, closestEnd);

    float distance = 0.0f;
    for (int i = closestStart; i < closestEnd; ++i) {
        distance += haversine(path[i], path[i+1]);
    }

    return distance;
}

float computePathDistance(const std::vector<Coordinate>& path, Coordinate entry, Coordinate exit) {
    int start = closestCoordinateIndex(path, entry);
    int end = closestCoordinateIndex(path, exit);

    if (start > end) std::swap(start, end);

    float totalDist = 0.0f;
    for (int i = start; i < end; ++i) {
        totalDist += haversine(path[i], path[i + 1]);
    }

    return totalDist; // in meters
}

float computeSegmentFare(const std::vector<Coordinate>& path, Coordinate entry, Coordinate exit) {
    float distanceMeters = computePathDistance(path, entry, exit);
    float distanceKm = distanceMeters / 1000.0f;
    return BASE_FARE + distanceKm * FARE_PER_KM;
}

double calculateTotalFare(const std::vector<TFTFEdge>& path, const Coordinate& startCoord, const Coordinate& endCoord) {
    double totalFare = 0.0;

    // If the path is empty, return 0 fare
    if (path.empty()) {
        return totalFare;
    }

    // First segment: from the start coordinate to the exit coordinate of the first edge
    std::vector<Coordinate> edgeCoordinates;
    for (const auto& edge : path) {
        edgeCoordinates.push_back(edge.entryCoord);
        edgeCoordinates.push_back(edge.exitCoord);
    }
    totalFare += computeSegmentFare(edgeCoordinates, startCoord, path[0].exitCoord);

    // Intermediate segments: from the entry to the exit coordinate of each edge
    for (size_t i = 1; i < path.size(); ++i) {
        std::vector<Coordinate> edgeCoordinates = {path[i-1].entryCoord, path[i].exitCoord};
        totalFare += computeSegmentFare(edgeCoordinates, path[i-1].entryCoord, path[i].exitCoord);
    }

    // Last segment: from the entry coordinate of the last edge to the end coordinate
    std::vector<Coordinate> lastEdgeCoordinates = {path.back().entryCoord, path.back().exitCoord};
    totalFare += computeSegmentFare(lastEdgeCoordinates, path.back().entryCoord, endCoord);

    return totalFare;
}

