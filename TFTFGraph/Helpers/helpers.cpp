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
float getActualSegmentDistance(const Coordinate& start, const Coordinate& end, const std::vector<Coordinate>& routePath) {
    if (routePath.empty()) return 0.0f;

    // Find index of the closest point to 'start'
    auto findClosestIndex = [](const Coordinate& target, const std::vector<Coordinate>& path) {
        int closestIdx = 0;
        float minDist = std::numeric_limits<float>::max();
        for (int i = 0; i < path.size(); ++i) {
            float dist = haversine(target, path[i]);
            if (dist < minDist) {
                minDist = dist;
                closestIdx = i;
            }
        }
        return closestIdx;
    };

    int startIdx = findClosestIndex(start, routePath);
    int endIdx = findClosestIndex(end, routePath);

    if (startIdx > endIdx) std::swap(startIdx, endIdx);  // Make sure we iterate forward

    float totalDist = 0.0f;
    for (int i = startIdx; i < endIdx; ++i) {
        totalDist += haversine(routePath[i], routePath[i + 1]);
    }

    return totalDist;
}

float computeSegmentFare(const std::vector<Coordinate>& path, Coordinate entry, Coordinate exit) {
    float distanceMeters = getActualSegmentDistance(entry, exit, path);
    float distanceKm = distanceMeters / 1000.0f;
    return BASE_FARE + distanceKm * FARE_PER_KM;
}

double calculateTotalFare(const std::vector<TFTFEdge>& path, const Coordinate& startCoord, const Coordinate& endCoord) {
    double totalFare = 0.0;

    

    return totalFare;
}

