#ifndef HELPERS_H
#define HELPERS_H
#include "../TFTFGraph.h"
#include "../../algorithms/node.h"
#include <map>
#include <vector>

const float BASE_FARE = 10.0f;
const float FARE_PER_KM = 2.0f;
const float EARTH_RADIUS_METERS = 6371000.0f;

float haversine(const Coordinate& a, const Coordinate& b);
int closestCoordinateIndex(const std::vector<Coordinate>& path, const Coordinate& target);
double getActualSegmentDistance(const Coordinate& start, const Coordinate& end, const std::vector<Coordinate>& path, bool isLoop);
Coordinate projectOntoPath(const Coordinate& point, const std::vector<Coordinate>& path) ;
int getClosestIndex(const std::vector<Coordinate>& path, const Coordinate& coord);
float getSubpathDistance(const std::vector<Coordinate>& coords, int i1, int i2, bool isLoop);
std::vector<Coordinate> densifyPath(const std::vector<Coordinate>& path, float spacingMeters) ;
double roundUpToNearest2_5(double amount);
std::vector<Coordinate> getFullSegmentPath(const std::vector<Coordinate> &path, int startIdx, int endIdx, bool isLoop);
std::vector<Coordinate> getShortestSegmentPath(
    const std::vector<Coordinate> &path,
    const Coordinate &startCoord,
    const Coordinate &endCoord,
    bool isLoop);
void printGraphDetails(const std::map<Node, std::vector<std::pair<Node, double>>> &graph);
void saveRuntimesRowToCSV(
    const Coordinate& startCoord,
    const Coordinate& endCoord,
    long long tftfDuration,
    long long dijkstraDuration,
    long long astarDuration,
    const std::string& filename,
    bool writeHeader = false);

#endif