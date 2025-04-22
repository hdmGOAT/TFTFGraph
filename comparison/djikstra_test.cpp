#include <iostream>
#include "../algorithms/djikstra/djikstra.h"

int main()
{
    Node Bulua_Terminal{8.511330, 124.624290};

    Node Bonbon{8.50881, 124.64827};
    Node Velez_Mogchs{8.482906, 124.646094};

    Node Kauswagan_City_Engineer{8.504775, 124.642954};
    Node USTP{8.484763, 124.655977};

    Node Camp_Evangelista{8.487358, 124.629950};

    std::string filename = "../routes.geojson";

    // std::cout << "Finding path from Kauswagan City Engineer to USTP...\n";
    // std::cout << "Finding path from Bonbon to Bulua Terminal...\n";
    std::cout << "Finding path from Camp Evangelista to USTP...\n";
    std::vector<Node> path = dijkstra_geojson(filename, Camp_Evangelista, USTP);

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

// RUN: g++ -std=c++17 comparison/djikstra_test.cpp algorithms/djikstra/djikstra.h algorithms/djikstra/djikstra.cpp -o comparison/djikstra