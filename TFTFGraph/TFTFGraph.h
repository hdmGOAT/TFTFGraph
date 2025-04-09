#ifndef TFTFGRAPH_H
#define TFTFGRAPH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

struct Coordinate {
    double latitude;
    double longitude;
};

bool operator==(const Coordinate& lhs, const Coordinate& rhs);

struct JeepneyDensity {
    int startHour; 
    int endHour;   
    float jeepneyDensity; 

    bool contains(int hour) const {
        return (startHour <= endHour && hour >= startHour && hour < endHour);
    }
};


struct TFTFEdge {
    int destinationRoute; 
    int originRoute;
    std::string destinationRouteName;
    float transferCost;
    std::vector<JeepneyDensity> densityByInterval;
    Coordinate entryCoord;
    Coordinate exitCoord;

    float totalCost(int hour = -1) const;
};

struct RouteNode {
    int routeId;
    std::string routeName;
    std::vector<TFTFEdge> edges;
    std::vector<Coordinate> path;
    std::vector<JeepneyDensity> densities;
};

class TFTFGraph {
public:
    void addRoute(int id, const std::string& name);
    void addEdge(int fromRoute, int toRoute, const std::string &toName,
        float transferCost, 
        const std::vector<JeepneyDensity> &densities = {}, Coordinate entryCoord = {}, Coordinate exitCoord = {});
    void visualize(int hour = -1) const;
    std::vector<TFTFEdge> findBestPath(int startRoute, int endRoute, int hour = -1);
    void setRoutePath(int routeId, const std::vector<Coordinate>& coordinates);
    void createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer = 10.0f);
    void setRouteDensities(int routeId, const std::vector<JeepneyDensity>& densities);
    int findClosestRoute(const Coordinate& startCoord);
    std::vector<TFTFEdge> calculateRouteFromCoordinates(const Coordinate& startCoord, const Coordinate& endCoord, int hour);
    double calculateTotalFare(const std::vector<TFTFEdge>& path, const Coordinate& startCoord, const Coordinate& endCoord);
    std::vector<const RouteNode*> extractTraversedRouteNodes(const std::vector<TFTFEdge>& path) const;
private:
    std::unordered_map<int, RouteNode> routes;

};

#endif
