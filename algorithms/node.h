#ifndef NODE_H
#define NODE_H

#include <vector>
#include <cmath>
#include "json.hpp"
#define M_PI 3.14159265358979323846


using json = nlohmann::json;

struct Node
{
    double lat, lon;

    bool operator<(const Node &other) const
    {
        return std::tie(lat, lon) < std::tie(other.lat, other.lon);
    }

    bool operator==(const Node &other) const
    {
        return lat == other.lat && lon == other.lon;
    }

    bool operator!=(const Node &other) const
    {
        return !(*this == other);
    }
};

double haversineNode(const Node &a, const Node &b);
void geojsonToNodeGraph(std::map<Node, std::vector<std::pair<Node, double>>> &graph, json file);

#endif