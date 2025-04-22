#include "../algorithms/astar/astar.h"
#include <iostream>

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
    std::vector<Node> path = astar_geojson(filename, Camp_Evangelista, USTP);
    return 0;
}

// RUN: g++ -std=c++17 comparison/astar_test.cpp algorithms/astar/astar.h algorithms/astar/astar.cpp -o comparison/astar