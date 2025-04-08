#ifndef HELPERS_H
#define HELPERS_H
#include "../TFTFGraph.h"

float haversine(const Coordinate& a, const Coordinate& b);
int closestCoordinateIndex(const std::vector<Coordinate>& path, const Coordinate& target);
#endif