#include "benchmark_utils.h"
#include <fstream>
#include <random>
#include <iostream>
#include <chrono>
#include "./algorithms/astar/astar.h"
#include "./algorithms/djikstra/djikstra.h"
#include <unordered_set>

void printProgressBar(int current, int total) {
    int barWidth = 50;
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

void initializeCSV(const std::string& filename) {
    std::ofstream file(filename);
    file << "from_lat,from_lon,to_lat,to_lon,tftf_ms,dijkstra_ms,astar_ms,"
         << "found_tftf,found_dijkstra,found_astar,"
         << "tftf_routes,dijkstra_routes,astar_routes\n";
    file.close();
}

void saveTestResult(const TestResult& result, const std::string& filename) {
    std::ofstream file(filename, std::ios::app);
    file << result.from.latitude << ","
         << result.from.longitude << ","
         << result.to.latitude << ","
         << result.to.longitude << ","
         << result.tftf_ms << ","
         << result.dijkstra_ms << ","
         << result.astar_ms << ","
         << result.found_tftf << ","
         << result.found_dijkstra << ","
         << result.found_astar << ","
         << result.tftf_routes_taken << ","
         << result.dijkstra_routes_taken << ","
         << result.astar_routes_taken << "\n";
    file.close();
}

std::pair<Coordinate, Coordinate> generateTestCoordinates(
    TestCategory category,
    const std::vector<Coordinate>& allCoordinates,
    TFTFGraph& network) {
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, allCoordinates.size() - 1);

    Coordinate from, to;
    auto& routes = network.getRoutes();
    auto routeCount = routes.size();

    switch (category) {
        case SAME_ROUTE: {
            if (routeCount == 0) {
                return {allCoordinates[distr(gen)], allCoordinates[distr(gen)]};
            }
            // Pick a random route and two points from it
            auto it = routes.begin();
            std::advance(it, distr(gen) % routeCount);
            const auto& route = it->second;
            const auto& path = route.path;
            
            if (path.empty()) {
                return {allCoordinates[distr(gen)], allCoordinates[distr(gen)]};
            }
            
            std::uniform_int_distribution<> pathDistr(0, path.size() - 1);
            from = path[pathDistr(gen)];
            to = path[pathDistr(gen)];
            break;
        }
        case DIFFERENT_ROUTES: {
            if (routeCount < 2) {
                return {allCoordinates[distr(gen)], allCoordinates[distr(gen)]};
            }
            // Pick two different routes and points from each
            auto it1 = routes.begin();
            std::advance(it1, distr(gen) % routeCount);
            auto it2 = routes.begin();
            std::advance(it2, distr(gen) % routeCount);
            while (it1 == it2) {
                it2 = routes.begin();
                std::advance(it2, distr(gen) % routeCount);
            }
            
            const auto& path1 = it1->second.path;
            const auto& path2 = it2->second.path;
            
            if (path1.empty() || path2.empty()) {
                return {allCoordinates[distr(gen)], allCoordinates[distr(gen)]};
            }
            
            std::uniform_int_distribution<> pathDistr1(0, path1.size() - 1);
            std::uniform_int_distribution<> pathDistr2(0, path2.size() - 1);
            
            from = path1[pathDistr1(gen)];
            to = path2[pathDistr2(gen)];
            break;
        }

        default: {
            // Fallback to random coordinates
            from = allCoordinates[distr(gen)];
            to = allCoordinates[distr(gen)];
            break;
        }
    }

    return {from, to};
}

TestResult runSingleTest(
    TFTFGraph& network,
    std::map<Node, std::vector<std::pair<Node, double>>>& nodeGraph,
    const Coordinate& from,
    const Coordinate& to) {
    
    TestResult result;
    result.from = from;
    result.to = to;
    result.tftf_routes_taken = 0;
    result.dijkstra_routes_taken = 0;
    result.astar_routes_taken = 0;

    // Test TFTF
    auto startTFTF = std::chrono::high_resolution_clock::now();
    auto tftfPath = network.calculateRouteFromCoordinates(from, to);
    auto endTFTF = std::chrono::high_resolution_clock::now();
    result.tftf_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endTFTF - startTFTF).count();
    result.found_tftf = !tftfPath.empty();

    // Count unique routes for TFTF if path was found
    if (result.found_tftf) {
        std::unordered_set<int> uniqueRoutes;
        auto routeNodes = network.extractTraversedRouteNodes(tftfPath);
        for (const auto* node : routeNodes) {
            uniqueRoutes.insert(node->routeId);
        }
        result.tftf_routes_taken = uniqueRoutes.size();
    }

    // Create transfer nodes for A* and Dijkstra
    Node fromNode(from.latitude, from.longitude);
    Node toNode(to.latitude, to.longitude);
    
    // Test A*
    auto startAStar = std::chrono::high_resolution_clock::now();
    auto astarPath = astar_geojson("routes.geojson", fromNode, toNode, nodeGraph);
    auto endAStar = std::chrono::high_resolution_clock::now();
    result.astar_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endAStar - startAStar).count();
    result.found_astar = !astarPath.empty();

    // Count unique routes for A*
    if (result.found_astar) {
        std::unordered_set<int> uniqueRoutes;
        for (const auto& node : astarPath) {
            if (node.routeId >= 0) { // Skip transfer points which have routeId = -1
                uniqueRoutes.insert(node.routeId);
            }
        }
        result.astar_routes_taken = uniqueRoutes.size();
    }

    // Test Dijkstra
    auto startDijkstra = std::chrono::high_resolution_clock::now();
    auto dijkstraPath = dijkstra_geojson("routes.geojson", fromNode, toNode, nodeGraph);
    auto endDijkstra = std::chrono::high_resolution_clock::now();
    result.dijkstra_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endDijkstra - startDijkstra).count();
    result.found_dijkstra = !dijkstraPath.empty();

    // Count unique routes for Dijkstra
    if (result.found_dijkstra) {
        std::unordered_set<int> uniqueRoutes;
        for (const auto& node : dijkstraPath) {
            if (node.routeId >= 0) { // Skip transfer points which have routeId = -1
                uniqueRoutes.insert(node.routeId);
            }
        }
        result.dijkstra_routes_taken = uniqueRoutes.size();
    }

    return result;
}

void runTestCategory(
    TFTFGraph& network,
    std::map<Node, std::vector<std::pair<Node, double>>>& nodeGraph,
    TestCategory category,
    const std::string& filename,
    int totalTests,
    const std::vector<Coordinate>& allCoordinates) {
    
    std::cout << "\nRunning tests for category: " 
              << (category == SAME_ROUTE ? "Same Route" :
                  category == DIFFERENT_ROUTES ? "Different Routes" : "Unknown") 
              << std::endl;

    initializeCSV(filename);

    for (int i = 0; i < totalTests; ++i) {
        auto [from, to] = generateTestCoordinates(category, allCoordinates, network);
        auto result = runSingleTest(network, nodeGraph, from, to);
        saveTestResult(result, filename);
        printProgressBar(i + 1, totalTests);
    }
    std::cout << std::endl;
} 