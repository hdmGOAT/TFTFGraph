#ifndef NODE_H
#define NODE_H

#include <vector>
#include <cmath>

#define M_PI 3.14159265358979323846

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
#endif