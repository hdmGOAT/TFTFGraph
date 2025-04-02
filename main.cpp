#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h" // Assuming TFTFGraph.h is in the same directory

void printPath(const std::vector<int>& bestPath) {
    std::cout << "Best path: ";
    for (size_t i = 0; i < bestPath.size(); ++i) {
        std::cout << bestPath[i];
        if (i < bestPath.size() - 1) {
            std::cout << " -> ";
        }
    }
    std::cout << std::endl;
}

int main() {
    TFTFGraph jeepneyNetwork;

    // Add routes
    jeepneyNetwork.addRoute(1, "A");
    jeepneyNetwork.addRoute(2, "B");
    jeepneyNetwork.addRoute(3, "C");
    jeepneyNetwork.addRoute(4, "D");
    jeepneyNetwork.addRoute(5, "E");
    jeepneyNetwork.addRoute(6, "F"); // Added route 6
    jeepneyNetwork.addRoute(7, "G"); // Added route 7
    jeepneyNetwork.addRoute(8, "H"); // Added route 8

    // Define density schedules
    std::vector<JeepneyDensity> morningRush = {{6, 9, 1.5f}};
    std::vector<JeepneyDensity> afternoonRush = {{16, 19, 1.3f}};
    std::vector<JeepneyDensity> lateNight = {{22, 5, 2.0f}};
    std::vector<JeepneyDensity> normalHours = {{9, 16, 1.0f}, {19, 22, 1.0f}};
    std::vector<JeepneyDensity> allHours;
    allHours.insert(allHours.end(), morningRush.begin(), morningRush.end());
    allHours.insert(allHours.end(), afternoonRush.begin(), afternoonRush.end());
    allHours.insert(allHours.end(), lateNight.begin(), lateNight.end());
    allHours.insert(allHours.end(), normalHours.begin(), normalHours.end());

    // Add edges (transfers) with varying densities and availability
    jeepneyNetwork.addEdge(1, 2, "A-B", 5.0f, 10.0f, allHours);
    jeepneyNetwork.addEdge(2, 3, "B-C", 3.0f, 15.0f, normalHours); // Only available during normal hours
    jeepneyNetwork.addEdge(1, 4, "A-D", 2.0f, 8.0f, morningRush); // Only available during morning rush
    jeepneyNetwork.addEdge(4, 5, "D-E", 1.0f, 10.0f, afternoonRush); // Only available during afternoon rush
    jeepneyNetwork.addEdge(5, 3, "E-C", 7.0f, 20.0f, lateNight); // Only available late night
    jeepneyNetwork.addEdge(2, 6, "B-F", 4.0f, 12.0f, allHours);
    jeepneyNetwork.addEdge(6, 7, "F-G", 6.0f, 18.0f, normalHours);
    jeepneyNetwork.addEdge(7, 8, "G-H", 2.0f, 5.0f, allHours);
    jeepneyNetwork.addEdge(8, 3, "H-C", 9.0f, 25.0f, lateNight);
    jeepneyNetwork.addEdge(1, 8, "A-H", 10.0f, 30.0f, morningRush); //direct but long.

    // Visualize the graph at different times
    jeepneyNetwork.visualize();
    jeepneyNetwork.visualize(8);
    jeepneyNetwork.visualize(14);
    jeepneyNetwork.visualize(23);

    // Test edge cases
    printPath(jeepneyNetwork.findBestPath(1, 3, 8)); // Morning Rush hour test
    printPath(jeepneyNetwork.findBestPath(1, 3, 14)); // afternoon hour test
    printPath(jeepneyNetwork.findBestPath(1, 3, 23)); // Late night test
    printPath(jeepneyNetwork.findBestPath(1, 3, 10)); // normal hours test
    printPath(jeepneyNetwork.findBestPath(1, 5, 17)); // route only in the afternoon
    printPath(jeepneyNetwork.findBestPath(1, 5, 2)); // route not available at this hour
    printPath(jeepneyNetwork.findBestPath(1, 7, 10)); // test a path with more than 3 stops.
    printPath(jeepneyNetwork.findBestPath(1, 8, 8)); // direct but long path test.

    return 0;
}