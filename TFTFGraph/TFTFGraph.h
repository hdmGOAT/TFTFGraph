#ifndef TFTFGRAPH_H
#define TFTFGRAPH_H

#include <string>
#include <vector>
#include <unordered_map>

struct JeepneyDensity {
    int startHour; 
    int endHour;   
    float jeepneyDensity; 

    bool contains(int hour) const {
        return (startHour <= endHour && hour >= startHour && hour < endHour);
    }
};


struct TFTFEdge {
    int destinationRoute; //What route the edge transfers to
    std::string destinationRouteName;
    float transferCost;
    float fare;
    std::vector<JeepneyDensity> densityByInterval;

    float totalCost(int hour = -1) const;
};

struct RouteNode {
    int routeId;
    std::string routeName;
    std::vector<TFTFEdge> edges;
};

class TFTFGraph {
public:
    void addRoute(int id, const std::string& name);
    void addEdge(int fromRoute, int toRoute, const std::string &toName,
                 float transferCost, float fare,
                 const std::vector<JeepneyDensity> &densities);
    void visualize(int hour = -1) const;
    std::vector<int> findBestPath(int startRoute, int endRoute, int hour = -1);


private:
    std::unordered_map<int, RouteNode> routes;
};

#endif
