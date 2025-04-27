#ifndef ASTAR_H
#define ASTAR_H

#include <vector>
#include <string>
#include "../node.h"


double haversine(const Node &a, const Node &b);
std::vector<Node> astar_geojson(const std::string &filename, Node start, Node goal, std::map<Node, std::vector<std::pair<Node, double>>> graph);

#endif // ASTAR_H
