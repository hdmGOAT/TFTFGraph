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


// Route 1: Balanced with evening peak
jeepneyNetwork.setRouteDensities(1, {
    {0, 6, 0.6f},     // Early morning: light traffic
    {6, 10, 1.3f},    // Morning rush
    {10, 16, 1.1f},   // Midday steady
    {16, 20, 1.7f},   // Evening rush
    {20, 24, 0.8f}    // Night
});

// Route 2: Busy all day, high in afternoon/evening
jeepneyNetwork.setRouteDensities(2, {
    {0, 6, 0.5f},     // Early low
    {6, 10, 1.5f},    // Morning rush
    {10, 16, 1.3f},   // Midday high
    {16, 20, 1.9f},   // Evening peak
    {20, 24, 0.9f}    // Slight night traffic
});

// Route 3: Afternoon and evening focused
jeepneyNetwork.setRouteDensities(3, {
    {0, 6, 0.7f},     // Quiet early
    {6, 10, 1.1f},    // Moderate morning
    {10, 16, 1.0f},   // Low midday
    {16, 20, 2.2f},   // Very heavy evening rush
    {20, 24, 1.0f}    // Moderate night
});

// Route 4: Light overall, with small peak
jeepneyNetwork.setRouteDensities(4, {
    {0, 6, 0.4f},     // Low early
    {6, 10, 1.0f},    // Small bump in morning
    {10, 16, 1.0f},   // Consistent low traffic
    {16, 20, 1.2f},   // Slight rise in evening
    {20, 24, 0.6f}    // Night fade
});

// Route 5: Typical commuter-heavy pattern
jeepneyNetwork.setRouteDensities(5, {
    {0, 6, 0.3f},     // Minimal pre-dawn
    {6, 10, 1.4f},    // Busy morning
    {10, 16, 1.2f},   // Steady midday
    {16, 20, 1.6f},   // Heavy rush home
    {20, 24, 0.7f}    // Calm night
});

// Route 6: High activity, especially in AM/PM
jeepneyNetwork.setRouteDensities(6, {
    {0, 6, 0.5f},     // Low pre-dawn
    {6, 10, 1.6f},    // Heavy morning
    {10, 16, 1.3f},   // Sustained midday
    {16, 20, 1.8f},   // High evening
    {20, 24, 0.9f}    // Some night traffic
});

// Route 7: Evening-heavy route
jeepneyNetwork.setRouteDensities(7, {
    {0, 6, 0.8f},     // Surprisingly active early
    {6, 10, 1.2f},    // Normal morning
    {10, 16, 1.0f},   // Stable midday
    {16, 20, 2.0f},   // Very busy evening
    {20, 24, 1.1f}    // Stays active into the night
});

// Route 8: Moderate all-day with evening bump
jeepneyNetwork.setRouteDensities(8, {
    {0, 6, 0.6f},     // Low early
    {6, 10, 1.0f},    // Normal morning
    {10, 16, 0.9f},   // Calm midday
    {16, 20, 1.4f},   // Active evening
    {20, 24, 0.7f}    // Night wind-down
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
