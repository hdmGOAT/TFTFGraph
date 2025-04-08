#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"


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

    jeepneyNetwork.visualize();
    jeepneyNetwork.visualize(8);
    jeepneyNetwork.visualize(14);
    jeepneyNetwork.visualize(23);

    // Test case: Use specific coordinates to find best path
    Coordinate startCoord1 = {14.5995, 120.9842};  // Start coordinates (route 1)
    Coordinate endCoord1 = {14.6010, 120.9860};    // End coordinates (route 2)
    int hour = 8;  // Morning hour

    std::vector<TFTFEdge> bestPath1 = jeepneyNetwork.calculateRouteFromCoordinates(startCoord1, endCoord1, hour);

    Coordinate startCoord2 = {14.6050, 120.9900};  // Start coordinates (route 5)
    Coordinate endCoord2 = {14.6005, 120.9852};    // End coordinates (route 8)
    int hour2 = 14;  // Afternoon hour

    std::vector<TFTFEdge> bestPath2 = jeepneyNetwork.calculateRouteFromCoordinates(startCoord2, endCoord2, hour2);

    Coordinate startCoord3 = {14.6030, 120.9880};  // Start coordinates (route 3)
    Coordinate endCoord3 = {14.6020, 120.9870};    // End coordinates (route 7)
    int hour3 = 20;  // Evening hour

    std::vector<TFTFEdge> bestPath3 = jeepneyNetwork.calculateRouteFromCoordinates(startCoord3, endCoord3, hour3);

    return 0;
}
