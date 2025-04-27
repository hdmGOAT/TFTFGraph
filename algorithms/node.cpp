#include "node.h"
#include "json.hpp"

using json = nlohmann::json;

double haversineNode(const Node &a, const Node &b)
{
    const double R = 6371e3;
    double lat1 = a.lat * M_PI / 180;
    double lat2 = b.lat * M_PI / 180;
    double dlat = (b.lat - a.lat) * M_PI / 180;
    double dlon = (b.lon - a.lon) * M_PI / 180;

    double h = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);

    return R * 2 * atan2(sqrt(h), sqrt(1 - h));
}

void geojsonToNodeGraph(std::map<Node, std::vector<std::pair<Node, double>>> &graph, json file){
        for (auto &feature : file["features"])
    {
        if (feature["geometry"]["type"] != "LineString")
            continue;
        auto coords = feature["geometry"]["coordinates"];

        for (size_t i = 0; i + 1 < coords.size(); ++i)
        {
            Node u{coords[i][1], coords[i][0]};
            Node v{coords[i + 1][1], coords[i + 1][0]};
            double dist = haversineNode(u, v);

            graph[u].emplace_back(v, dist);
            graph[v].emplace_back(u, dist);
        }
    }
};
