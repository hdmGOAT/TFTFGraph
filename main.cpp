#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"
#include "./algorithms/astar/astar.h"
#include "./algorithms/djikstra/djikstra.h"
#include <fstream>
#include <chrono>
#include <random>
#include "json.hpp"
#include "benchmark_utils.h"

using json = nlohmann::json;

void loadRoutesFromGeoJSON(const std::string &filepath, TFTFGraph &graph)
{
    std::ifstream file(filepath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open GeoJSON file!" << std::endl;
        return;
    }

    json geojson;
    file >> geojson;

    const auto &features = geojson["features"];
    int routeId = 0;

    for (const auto &feature : features)
    {
        if (feature["geometry"]["type"] != "LineString")
            continue;

        std::vector<Coordinate> routePath;
        const auto &coords = feature["geometry"]["coordinates"];

        for (const auto &coord : coords)
        {
            double lon = coord[0];
            double lat = coord[1];
            routePath.emplace_back(Coordinate{lat, lon});
        }

        std::string routeName = "Route_" + std::to_string(routeId);
        if (feature.contains("properties") && feature["properties"].contains("name"))
        {
            routeName = feature["properties"]["name"];
        }

        bool isLoop = routePath.front() == routePath.back();

        graph.addRoute(routeId, routeName);
        graph.setRoutePath(routeId, routePath);
        routeId++;
    }
}

std::vector<Coordinate> extractAllCoordinates(const std::string &filepath) {
    std::vector<Coordinate> allCoordinates;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open GeoJSON file!" << std::endl;
        return allCoordinates;
    }

    json geojson;
    file >> geojson;
    const auto &features = geojson["features"];

    for (const auto &feature : features) {
        if (feature["geometry"]["type"] != "LineString")
            continue;

        const auto &coords = feature["geometry"]["coordinates"];
        for (const auto &coord : coords) {
            double lon = coord[0];
            double lat = coord[1];
            allCoordinates.emplace_back(Coordinate{lat, lon});
        }
    }

    return allCoordinates;
}

int main()
{
    // Start timing TFTF graph construction
    std::cout << "Constructing TFTF graph..." << std::endl;
    auto tftf_start = std::chrono::high_resolution_clock::now();

    TFTFGraph jeepneyNetwork;
    loadRoutesFromGeoJSON("routes.geojson", jeepneyNetwork);

    float transferRange = 300.5f;
    jeepneyNetwork.createTransfersFromCoordinates(transferRange);

    auto tftf_end = std::chrono::high_resolution_clock::now();
    auto tftf_duration = std::chrono::duration_cast<std::chrono::milliseconds>(tftf_end - tftf_start);
    std::cout << "TFTF graph construction completed in " << tftf_duration.count() << " ms" << std::endl;

    // Start timing node graph construction
    std::cout << "\nConstructing node graph..." << std::endl;
    auto node_start = std::chrono::high_resolution_clock::now();

    // Parse GeoJSON for node graph
    std::ifstream file("routes.geojson");
    json geojson;
    file >> geojson;
    std::map<Node, std::vector<std::pair<Node, double>>> nodeGraph;
    geojsonToNodeGraph(nodeGraph, geojson);

    auto node_end = std::chrono::high_resolution_clock::now();
    auto node_duration = std::chrono::duration_cast<std::chrono::milliseconds>(node_end - node_start);
    std::cout << "Node graph construction completed in " << node_duration.count() << " ms" << std::endl;

    // Print initial graph details
    std::cout << "\nGraph Details:" << std::endl;
    std::cout << "-------------" << std::endl;
    printGraphDetails(nodeGraph);
    jeepneyNetwork.getGraphDetails();

    // Extract all coordinates for testing
    std::vector<Coordinate> allCoordinates = extractAllCoordinates("routes.geojson");

    const int testsPerCategory = 1500;
    
    // Run tests for each category
    runTestCategory(jeepneyNetwork, nodeGraph, SAME_ROUTE, "same_route.csv", 
                   testsPerCategory, allCoordinates);
    
    runTestCategory(jeepneyNetwork, nodeGraph, DIFFERENT_ROUTES, "different_routes.csv", 
                   testsPerCategory, allCoordinates);
    

    std::cout << "\nAll tests completed successfully!" << std::endl;
    return 0;
}

// g++ -std=c++17 main.cpp TFTFGraph/TFTFGraph.cpp TFTFGraph/Helpers/helpers.cpp algorithms/astar/astar.cpp algorithms/djikstra/djikstra.cpp algorithms/node.cpp benchmark_utils.cpp -o main
