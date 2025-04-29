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

void testRoute(TFTFGraph network, std::map<Node, std::vector<std::pair<Node, double>>> nodeGraph,Coordinate from, std::string fromName, Coordinate to, std::string toName, float transferMeters)
{

    // std::cout << "Transfer range: " << transferMeters << " meters" << std::endl;
    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "From: " << fromName << " (" << from.latitude << ", " << from.longitude << ")" << std::endl;
    // std::cout << "To: " << toName << " (" << to.latitude << ", " << to.longitude << ")" << std::endl;
    // std::cout << "--------------------------------------------------" << std::endl;


    auto startTFTF = std::chrono::high_resolution_clock::now();
    network.calculateRouteFromCoordinates(from, to);
    auto endTFTF = std::chrono::high_resolution_clock::now();
    long long tftfDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTFTF - startTFTF).count();


    Node fromNode{from.latitude, from.longitude};
    Node toNode{to.latitude, to.longitude};

    // std::cout << "TFTFGRAPH" << std::endl;
    // network.getGraphDetails();

    // std::cout << "--------------------------------------------------" << std::endl;

    // std::cout << "A STAR" << std::endl;
    auto startAStar = std::chrono::high_resolution_clock::now();
    astar_geojson("routes.geojson", fromNode, toNode, nodeGraph);
    auto endAStar = std::chrono::high_resolution_clock::now();
    long long aStarDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endAStar - startAStar).count();
    // std::cout << "--------------------------------------------------" << std::endl;
    // std::cout << "DIJKSTRA" << std::endl;
    auto startDijkstra = std::chrono::high_resolution_clock::now();
    dijkstra_geojson("routes.geojson",fromNode, toNode, nodeGraph);
    auto endDijkstra = std::chrono::high_resolution_clock::now();
    long long dijkstraDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endDijkstra - startDijkstra).count();

    static bool firstRow = true;  
    saveRuntimesRowToCSV(from, to, tftfDuration, dijkstraDuration, aStarDuration, "runtimes.csv", firstRow);
    firstRow = false;
}

int main()
{
    TFTFGraph jeepneyNetwork;
    loadRoutesFromGeoJSON("routes.geojson", jeepneyNetwork);

    float transferRange = 300.5f;

    // Create transfers
    jeepneyNetwork.createTransfersFromCoordinates(transferRange);

        // Parse GeoJSON
    std::ifstream file("routes.geojson");
    json geojson;
    file >> geojson;
    std::map<Node, std::vector<std::pair<Node, double>>> nodeGraph;
    geojsonToNodeGraph(nodeGraph, geojson);


    printGraphDetails(nodeGraph); // Print graph details
    jeepneyNetwork.getGraphDetails(); // Print TFTFGraph details

    // Extract all coordinates
    std::vector<Coordinate> allCoordinates = extractAllCoordinates("routes.geojson");

    // Random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, allCoordinates.size() - 1);

    // Do 500 random tests
    int totalTests = 5000;

    for (int i = 0; i < totalTests; ++i) {
        Coordinate from = allCoordinates[distr(gen)];
        Coordinate to = allCoordinates[distr(gen)];

        while (from.latitude == to.latitude && from.longitude == to.longitude) {
            to = allCoordinates[distr(gen)];
        }

        std::string fromName = "RandomFrom_" + std::to_string(i);
        std::string toName = "RandomTo_" + std::to_string(i);

        testRoute(jeepneyNetwork, nodeGraph, from, fromName, to, toName, transferRange);

        printProgressBar(i + 1, totalTests);  // <--- Add this
    }

    std::cout << std::endl; // After loop finishes, move to new line

    return 0;
}



// g++ -std=c++17 main.cpp TFTFGraph/TFTFGraph.cpp TFTFGraph/Helpers/helpers.cpp algorithms/astar/astar.cpp algorithms/djikstra/djikstra.cpp algorithms/node.cpp -o main
