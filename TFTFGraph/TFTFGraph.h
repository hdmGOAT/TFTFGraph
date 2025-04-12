#ifndef TFTFGRAPH_H
#define TFTFGRAPH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

struct PairHash {
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

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
    int entryIndex = -1;
    int exitIndex = -1;


    float totalCost(int hour = -1) const;
};

struct RouteNode {
    int routeId;
    std::string routeName;
    std::vector<TFTFEdge> edges;
    std::vector<Coordinate> path;
    std::vector<JeepneyDensity> densities;
    bool isLoop = false; // indicates if the route is a loop
};

class TFTFGraph {
    public:
        void addRoute(int id, const std::string& name);
        void addEdge(int fromRoute, int toRoute, const std::string &toName,
            float transferCost, 
            const std::vector<JeepneyDensity> &densities = {}, Coordinate entryCoord = {}, Coordinate exitCoord = {});
        void visualize(int hour = -1) const;
        void setRoutePath(int routeId, const std::vector<Coordinate>& coordinates);
        void createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer = 10.0f);
        void setRouteDensities(int routeId, const std::vector<JeepneyDensity>& densities);
       std::vector<int> getNearbyRoutes(const Coordinate& coord, float maxDistanceMeters);
        std::vector<TFTFEdge> calculateRouteFromCoordinates(const Coordinate& startCoord, const Coordinate& endCoord, int hour);
        double calculateTotalFare(const std::vector<TFTFEdge>& path, const Coordinate& startCoord, const Coordinate& endCoord);
        std::vector<const RouteNode*> extractTraversedRouteNodes(const std::vector<TFTFEdge>& path) const;
            std::vector<TFTFEdge> findMinFarePath(int startRouteId, int endRouteId, int hour);
        std::vector<TFTFEdge> findMinTransfersPath(int startRouteId, int endRouteId, int hour = -1); 
    private:
        std::unordered_map<int, RouteNode> routes;

};

#endif
