#ifndef ASTAR_H
#define ASTAR_H

#include <vector>
#include <string>

struct Node
{
    double lat, lon;

    bool operator<(const Node &other) const;
    bool operator==(const Node &other) const;
    bool operator!=(const Node &other) const;
};

double haversine(const Node &a, const Node &b);
std::vector<Node> astar_geojson(const std::string &filename, Node start, Node goal);

#endif // ASTAR_H
