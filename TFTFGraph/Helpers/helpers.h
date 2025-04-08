#ifndef HELPERS_H
#define HELPERS_H
#include "../TFTFGraph.h"

float haversine(const Coordinate& a, const Coordinate& b);
int closestCoordinateIndex(const std::vector<Coordinate>& path, const Coordinate& target);
float computeRouteDistance(const std::vector<Coordinate>& path, const Coordinate& start, const Coordinate& end);
float computeSegmentFare(const std::vector<Coordinate>& path, Coordinate entry, Coordinate exit);
float computePathDistance(const std::vector<Coordinate>& path, Coordinate entry, Coordinate exit);
#endif