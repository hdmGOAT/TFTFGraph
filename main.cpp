#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <map>
#include "json.hpp"
#include "./TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"
#include "./algorithms/astar/astar.h"
#include "./algorithms/djikstra/djikstra.h"
#include "benchmark_utils.h"

using json = nlohmann::json;

// ---------------- existing helpers (unchanged) ----------------
void loadRoutesFromGeoJSON(const std::string &filepath, TFTFGraph &graph) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open GeoJSON file: " << filepath << std::endl;
        return;
    }
    json geojson; file >> geojson;
    const auto &features = geojson["features"];
    int routeId = 0;

    for (const auto &feature : features) {
        if (feature["geometry"]["type"] != "LineString")
            continue;
        std::vector<Coordinate> routePath;
        const auto &coords = feature["geometry"]["coordinates"];
        for (const auto &coord : coords) {
            double lon = coord[0];
            double lat = coord[1];
            routePath.emplace_back(Coordinate{lat, lon});
        }
        std::string routeName = "Route_" + std::to_string(routeId);
        if (feature.contains("properties") && feature["properties"].contains("name")) {
            routeName = feature["properties"]["name"];
        }
        graph.addRoute(routeId, routeName);
        graph.setRoutePath(routeId, routePath);
        routeId++;
    }
}

std::vector<Coordinate> extractAllCoordinates(const std::string &filepath) {
    std::vector<Coordinate> allCoordinates;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open GeoJSON file: " << filepath << std::endl;
        return allCoordinates;
    }
    json geojson; file >> geojson;
    const auto &features = geojson["features"];
    for (const auto &feature : features) {
        if (feature["geometry"]["type"] != "LineString")
            continue;
        for (const auto &coord : feature["geometry"]["coordinates"]) {
            double lon = coord[0];
            double lat = coord[1];
            allCoordinates.emplace_back(Coordinate{lat, lon});
        }
    }
    return allCoordinates;
}
// --------------------------------------------------------------

// 🔑 New reusable benchmark runner
void runBenchmark(const std::string &geojsonFile, const std::string &label) {
    std::cout << "\n=============================\n";
    std::cout << "Processing network: " << label << "\n";
    std::cout << "=============================\n";

    // --- TFTF graph construction ---
    auto tftf_start = std::chrono::high_resolution_clock::now();
    TFTFGraph network;
    loadRoutesFromGeoJSON(geojsonFile, network);
    network.createTransfersFromCoordinates(300.5f);
    auto tftf_end = std::chrono::high_resolution_clock::now();
    auto tftf_duration = std::chrono::duration_cast<std::chrono::milliseconds>(tftf_end - tftf_start);
    std::cout << "TFTF graph done in " << tftf_duration.count() << " ms\n";

    // --- Node graph construction ---
    auto node_start = std::chrono::high_resolution_clock::now();
    std::ifstream file(geojsonFile);
    json geojson; file >> geojson;
    std::map<Node, std::vector<std::pair<Node,double>>> nodeGraph;
    geojsonToNodeGraph(nodeGraph, geojson);
    auto node_end = std::chrono::high_resolution_clock::now();
    auto node_duration = std::chrono::duration_cast<std::chrono::milliseconds>(node_end - node_start);
    std::cout << "Node graph done in " << node_duration.count() << " ms\n";

    // --- Details & tests ---
    printGraphDetails(nodeGraph);
    network.getGraphDetails();
    auto allCoords = extractAllCoordinates(geojsonFile);

    const int testsPerCategory = 1500;
    runTestCategory(network, nodeGraph, SAME_ROUTE,
                    label + "_same_route.csv", testsPerCategory, allCoords);
    runTestCategory(network, nodeGraph, DIFFERENT_ROUTES,
                    label + "_different_routes.csv", testsPerCategory, allCoords);

    std::cout << "✅ " << label << " tests completed.\n";
}

int main() {
    // 🔑 Call runner for each district/network
    runBenchmark("district1.geojson", "district1");
    runBenchmark("district2.geojson", "district2");
    runBenchmark("allRoutes.geojson", "allRoutes");

    std::cout << "\nAll networks processed successfully!\n";
    return 0;
}

// Compile as before:
// g++ -std=c++17 main.cpp TFTFGraph/TFTFGraph.cpp TFTFGraph/Helpers/helpers.cpp \
//     algorithms/astar/astar.cpp algorithms/djikstra/djikstra.cpp algorithms/node.cpp \
//     benchmark_utils.cpp -o main
