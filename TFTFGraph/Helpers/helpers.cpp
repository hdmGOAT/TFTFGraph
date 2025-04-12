#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
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

double getActualSegmentDistance(const Coordinate& start, const Coordinate& end, const std::vector<Coordinate>& path, bool isLoop) {
    int startIndex = getClosestIndex(path, start);
    int endIndex = getClosestIndex(path, end);

    double distance = 0.0;

    if (startIndex <= endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            distance += haversine(path[i], path[i + 1]);
        }
    } else if (isLoop) {
        // allow wrap around
        for (int i = startIndex; i < path.size() - 1; ++i) {
            distance += haversine(path[i], path[i + 1]);
        }
        distance += haversine(path.back(), path[0]);
        for (int i = 0; i < endIndex; ++i) {
            distance += haversine(path[i], path[i + 1]);
        }
    } else {
        // invalid segment on non-looping route
        std::cerr << "Warning: Non-loop route but segment wraps around. Returning large distance.\n";
        return 1e9;
    }

    return distance;
}




Coordinate projectOntoPath(const Coordinate& point, const std::vector<Coordinate>& path) {
    double minDistance = std::numeric_limits<double>::max();
    Coordinate closestPoint;

    for (size_t i = 0; i + 1 < path.size(); ++i) {
        const Coordinate& A = path[i];
        const Coordinate& B = path[i + 1];

        double dx = B.longitude - A.latitude;
        double dy = B.latitude - A.latitude;

        double lengthSquared = dx * dx + dy * dy;
        if (lengthSquared == 0.0) continue;  

        double t = ((point.longitude - A.longitude) * dx + (point.latitude - A.latitude) * dy) / lengthSquared;
        t = std::max(0.0, std::min(1.0, t)); 

        Coordinate projection = {
            A.latitude + t * dy,
            A.latitude + t * dx
        };

        double d = haversine(point, projection);
        if (d < minDistance) {
            minDistance = d;
            closestPoint = projection;
        }
    }

    return closestPoint;
}

int getClosestIndex(const std::vector<Coordinate>& path, const Coordinate& coord) {
    float minDist = std::numeric_limits<float>::infinity();
    int index = -1;

    for (size_t i = 0; i < path.size(); ++i) {
        float dist = haversine(coord, path[i]);
        if (dist < minDist) {
            minDist = dist;
            index = static_cast<int>(i);
        }
    }

    return index;
}