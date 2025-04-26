#ifndef DJIKSTRA_H
#define DJIKSTRA_H

#include <string>
#include <vector>
#include "../json.hpp"
#include "../node.h"



double haversine(const Node &a, const Node &b);

std::vector<Node> dijkstra_geojson(const std::string &filename, Node origin, Node destination);

#endif
