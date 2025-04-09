#ifndef HELPERS_H
#define HELPERS_H
#include "../TFTFGraph.h"

const float BASE_FARE = 10.0f;
const float FARE_PER_KM = 2.0f;
const float EARTH_RADIUS_METERS = 6371000.0f;

float haversine(const Coordinate& a, const Coordinate& b);
int closestCoordinateIndex(const std::vector<Coordinate>& path, const Coordinate& target);
float computeRouteDistance(const std::vector<Coordinate>& path, const Coordinate& start, const Coordinate& end);
float computePathDistance(const std::vector<Coordinate>& path, Coordinate entry, Coordinate exit);

float getActualSegmentDistance(const Coordinate& start, const Coordinate& end, const std::vector<Coordinate>& routePath);

#endif