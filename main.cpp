#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h"

void printPath(const std::vector<int>& bestPath) {
    std::cout << "Best path: ";
    for (size_t i = 0; i < bestPath.size(); ++i) {
        std::cout << bestPath[i];
        if (i < bestPath.size() - 1) std::cout << " -> ";
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
    jeepneyNetwork.addRoute(6, "F");
    jeepneyNetwork.addRoute(7, "G");
    jeepneyNetwork.addRoute(8, "H");

    // Add route coordinates (paths)
    jeepneyNetwork.setRoutePath(1, {{14.5995, 120.9842}, {14.6000, 120.9850}});
    jeepneyNetwork.setRoutePath(2, {{14.6001, 120.9851}, {14.6010, 120.9860}});
    jeepneyNetwork.setRoutePath(3, {{14.6030, 120.9880}, {14.6040, 120.9890}});
    jeepneyNetwork.setRoutePath(4, {{14.5980, 120.9830}});
    jeepneyNetwork.setRoutePath(5, {{14.6050, 120.9900}});
    jeepneyNetwork.setRoutePath(6, {{14.6015, 120.9865}});
    jeepneyNetwork.setRoutePath(7, {{14.6020, 120.9870}});
    jeepneyNetwork.setRoutePath(8, {{14.5990, 120.9840}, {14.6005, 120.9852}});

    // Create transfers based on proximity (e.g., 200 meters)
    jeepneyNetwork.createTransfersFromCoordinates(200.0f);

    // Precompute hop distances and visualize
    jeepneyNetwork.precomputeHopDistances();
    jeepneyNetwork.printHopDistances();
    jeepneyNetwork.visualize();
    jeepneyNetwork.visualize(8);
    jeepneyNetwork.visualize(14);
    jeepneyNetwork.visualize(23);

    // Test paths at different times of day
    printPath(jeepneyNetwork.findBestPath(1, 3, 8));   // Morning
    printPath(jeepneyNetwork.findBestPath(1, 3, 14));  // Afternoon
    printPath(jeepneyNetwork.findBestPath(1, 3, 23));  // Late Night
    printPath(jeepneyNetwork.findBestPath(1, 3, 10));  // Normal hours
    printPath(jeepneyNetwork.findBestPath(1, 5, 17));  // Try longer path
    printPath(jeepneyNetwork.findBestPath(1, 7, 10));  // Multiple hops
    printPath(jeepneyNetwork.findBestPath(1, 8, 8));   // Direct possible

    return 0;
}
