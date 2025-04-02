#include <iostream>
#include "./TFTFGraph/TFTFGraph.h"

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

int main()
{
    TFTFGraph jeepneyNetwork;

    jeepneyNetwork.addRoute(1, "Ayala - Guadalupe");
    jeepneyNetwork.addRoute(2, "Guadalupe - Cubao");
    jeepneyNetwork.addRoute(3, "Cubao - Monumento");
    jeepneyNetwork.addRoute(4, "Guadalupe - Makati");
    jeepneyNetwork.addRoute(5, "Makati - Pasay");
    
    std::vector<JeepneyDensity> morningRush = {{6, 9, 1.5f}}; // Higher cost during rush hour
    std::vector<JeepneyDensity> afternoonRush = {{16, 19, 1.3f}};
    std::vector<JeepneyDensity> lateNight = {{22, 5, 2.0f}}; // Higher cost due to scarcity
    std::vector<JeepneyDensity> normalHours = {{9, 16, 1.0f}, {19, 22, 1.0f}};

    // Combine all hours for some edges
    std::vector<JeepneyDensity> allHours;
    allHours.insert(allHours.end(), morningRush.begin(), morningRush.end());
    allHours.insert(allHours.end(), afternoonRush.begin(), afternoonRush.end());
    allHours.insert(allHours.end(), lateNight.begin(), lateNight.end());
    allHours.insert(allHours.end(), normalHours.begin(), normalHours.end());

    // Add edges (transfers) between routes
    jeepneyNetwork.addEdge(1, 2, "Guadalupe - Cubao", 5.0f, 12.0f, allHours);
    jeepneyNetwork.addEdge(1, 4, "Guadalupe - Makati", 2.0f, 8.0f, allHours);
    jeepneyNetwork.addEdge(2, 3, "Cubao - Monumento", 3.0f, 15.0f, allHours);
    jeepneyNetwork.addEdge(4, 5, "Makati - Pasay", 1.0f, 10.0f, allHours);
    jeepneyNetwork.addEdge(5, 3, "Cubao - Monumento", 7.0f, 20.0f, lateNight); // Only available late night

    // Visualize the graph at different times
    jeepneyNetwork.visualize();   // Without time consideration
    jeepneyNetwork.visualize(8);  // Morning rush hour
    jeepneyNetwork.visualize(14); // Mid-afternoon
    jeepneyNetwork.visualize(23); // Late night

    // Demonstrate path finding (placeholder)
    printPath(jeepneyNetwork.findBestPath(1, 3, 8));

    printPath(jeepneyNetwork.findBestPath(1, 3, 2));

    return 0;
}