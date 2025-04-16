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



struct TFTFEdge {
    int destinationRoute; 
    int originRoute;
    std::string destinationRouteName;
    float transferCost;
    Coordinate entryCoord;
    Coordinate exitCoord;
    int entryIndex = -1;
    int exitIndex = -1;


    float totalCost() const;
};

struct RouteNode {
    int routeId;
    std::string routeName;
    std::vector<TFTFEdge> edges;
    std::vector<Coordinate> path;
    bool isLoop = false; // indicates if the route is a loop
};

class TFTFGraph {
    public:
        void addRoute(int id, const std::string& name);
        void addEdge(int fromRoute, int toRoute, const std::string &toName,
            float transferCost, Coordinate entryCoord = {}, Coordinate exitCoord = {});
        void visualize() const;
        void setRoutePath(int routeId, const std::vector<Coordinate>& coordinates);
        void createTransfersFromCoordinates(float transferRangeMeters, float farePerTransfer = 10.0f);
       std::vector<int> getNearbyRoutes(const Coordinate& coord, float maxDistanceMeters);
        std::vector<TFTFEdge> calculateRouteFromCoordinates(const Coordinate& startCoord, const Coordinate& endCoord, int hour);
        double calculateTotalFare(const std::vector<TFTFEdge>& path, const Coordinate& startCoord, const Coordinate& endCoord);
        std::vector<const RouteNode*> extractTraversedRouteNodes(const std::vector<TFTFEdge>& path) const;
            std::vector<TFTFEdge> findMinFarePath(int startRouteId, int endRouteId, int hour, int projectedStartIdx);
       
    private:
        std::unordered_map<int, RouteNode> routes;

};

#endif
