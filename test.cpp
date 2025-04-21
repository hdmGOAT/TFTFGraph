#include <iostream>
#include "./algorithms/djikstra/djikstra.h"

int main()
{
    Node origin{8.508810, 124.648270};      // Example coordinates (lat, lon)
    Node destination{8.511330, 124.624290}; // Change as needed

    std::string filename = "routes.geojson"; // Your GeoJSON file

    std::vector<Node> path = dijkstra_geojson(filename, origin, destination);

    if (path.empty())
    {
        std::cout << "No path found between origin and destination.\n";
    }
    else
    {
        std::cout << "Path found with " << path.size() << " points.\n";
    }

    return 0;
}