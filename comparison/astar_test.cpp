#include "../algorithms/astar/astar.h"

int main()
{
    Node origin{8.508810, 124.648270};
    Node destination{8.511330, 124.624290};

    std::string filename = "../routes.geojson";
    auto path = astar_geojson(filename, origin, destination);
    return 0;
}

// RUN: g++ -std=c++17 comparison/astar-test.cpp algorithms/astar/astar.h algorithms/astar/astar.cpp -o comparison/astar