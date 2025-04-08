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
    std::string destinationRouteName;
    float transferCost;
    float fare;
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

struct TransferPoint {
    int fromRoute;
    Coordinate fromCoord;
    int toRoute;
    Coordinate toCoord;
    float transferDistance;
};


class TFTFGraph {
public:
    void addRoute(int id, const std::string& name);
    void addEdge(int fromRoute, int toRoute, const std::string &toName,
        float transferCost, float fare,
        const std::vector<JeepneyDensity> &densities = {}, Coordinate entryCoord = {}, Coordinate exitCoord = {});
    void visualize(int hour = -1) const;
    std::vector<int> findBestPath(int startRoute, int endRoute, int hour = -1);
    void precomputeHopDistances();
    float heuristic(int current, int target);
    void printHopDistances() const;
    void setRoutePath(int routeId, const std::vector<Coordinate>& coordinates);
    void createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer = 10.0f);
    void setRouteDensities(int routeId, const std::vector<JeepneyDensity>& densities);
    TFTFEdge* getEdge(int fromRoute, int toRoute) const;



private:
    std::unordered_map<int, RouteNode> routes;
    std::unordered_map<int, std::unordered_map<int, int>> hopDistance;
    float minEdgeCost = std::numeric_limits<float>::infinity();

};

#endif
