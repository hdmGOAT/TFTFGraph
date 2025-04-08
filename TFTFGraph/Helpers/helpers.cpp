#define _USE_MATH_DEFINES
#include <cmath>
#include "helpers.h"
#include "../TFTFGraph.h"

const float EARTH_RADIUS_METERS = 6371000.0f;

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

