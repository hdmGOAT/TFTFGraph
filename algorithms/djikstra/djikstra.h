#ifndef DJIKSTRA_H
#define DJIKSTRA_H

#include <string>
#include <vector>
#include "json.hpp"

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

double haversine(const Node &a, const Node &b);

std::vector<Node> dijkstra_geojson(const std::string &filename, Node origin, Node destination);

#endif
