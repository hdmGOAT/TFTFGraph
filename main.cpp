#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h"

// Function to print the best path with fares
void printPathWithFares(const std::vector<int>& bestPath, const TFTFGraph& jeepneyNetwork, int timeOfDay) {
    std::cout << "Best path (with fares): ";
    float totalFare = 0.0f;

    for (size_t i = 0; i < bestPath.size() - 1; ++i) {
        int from = bestPath[i];
        int to = bestPath[i + 1];
        float segmentFare = jeepneyNetwork.getEdge(from, to)->fare;  // Get the fare for this segment
        totalFare += segmentFare;

        std::cout << from;
        if (i < bestPath.size() - 2) std::cout << " -> ";
    }
    std::cout << "\nTotal fare: PHP " << totalFare << std::endl;
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

    // Add route densities (traffic conditions at different times of day)
    jeepneyNetwork.setRouteDensities(1, {
        {0, 6, 0.6f}, {6, 10, 1.3f}, {10, 16, 1.1f}, {16, 20, 1.7f}, {20, 24, 0.8f}
    });

    jeepneyNetwork.setRouteDensities(2, {
        {0, 6, 0.5f}, {6, 10, 1.5f}, {10, 16, 1.3f}, {16, 20, 1.9f}, {20, 24, 0.9f}
    });

    jeepneyNetwork.setRouteDensities(3, {
        {0, 6, 0.7f}, {6, 10, 1.1f}, {10, 16, 1.0f}, {16, 20, 2.2f}, {20, 24, 1.0f}
    });

    jeepneyNetwork.setRouteDensities(4, {
        {0, 6, 0.4f}, {6, 10, 1.0f}, {10, 16, 1.0f}, {16, 20, 1.2f}, {20, 24, 0.6f}
    });

    jeepneyNetwork.setRouteDensities(5, {
        {0, 6, 0.3f}, {6, 10, 1.4f}, {10, 16, 1.2f}, {16, 20, 1.6f}, {20, 24, 0.7f}
    });

    jeepneyNetwork.setRouteDensities(6, {
        {0, 6, 0.5f}, {6, 10, 1.6f}, {10, 16, 1.3f}, {16, 20, 1.8f}, {20, 24, 0.9f}
    });

    jeepneyNetwork.setRouteDensities(7, {
        {0, 6, 0.8f}, {6, 10, 1.2f}, {10, 16, 1.0f}, {16, 20, 2.0f}, {20, 24, 1.1f}
    });

    jeepneyNetwork.setRouteDensities(8, {
        {0, 6, 0.6f}, {6, 10, 1.0f}, {10, 16, 0.9f}, {16, 20, 1.4f}, {20, 24, 0.7f}
    });

    // Create transfers based on proximity (e.g., 200 meters)
    jeepneyNetwork.createTransfersFromCoordinates(200.0f);

    // Precompute hop distances and visualize
    jeepneyNetwork.precomputeHopDistances();
    jeepneyNetwork.printHopDistances();
    jeepneyNetwork.visualize();
    jeepneyNetwork.visualize(8);
    jeepneyNetwork.visualize(14);
    jeepneyNetwork.visualize(23);

    // Longer and multiple-hop path tests
    printPathWithFares(jeepneyNetwork.findBestPath(1, 5, 8), jeepneyNetwork, 8);  // Morning (Longer path)
    printPathWithFares(jeepneyNetwork.findBestPath(1, 6, 14), jeepneyNetwork, 14); // Afternoon (Multiple hops)
    printPathWithFares(jeepneyNetwork.findBestPath(1, 7, 23), jeepneyNetwork, 23); // Late night (Long journey)
    printPathWithFares(jeepneyNetwork.findBestPath(1, 8, 17), jeepneyNetwork, 17); // Peak hour (Rush hour path)
    printPathWithFares(jeepneyNetwork.findBestPath(2, 8, 20), jeepneyNetwork, 20); // Evening rush (Long path)
    printPathWithFares(jeepneyNetwork.findBestPath(4, 6, 10), jeepneyNetwork, 10); // Midday path (Simple hop)
    printPathWithFares(jeepneyNetwork.findBestPath(5, 3, 16), jeepneyNetwork, 16); // Mid-afternoon (Intermediate path)
    printPathWithFares(jeepneyNetwork.findBestPath(2, 4, 6), jeepneyNetwork, 6);  // Early morning (Shorter path)

    return 0;
}
