#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "./TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"
#include "json.hpp"

using json = nlohmann::json;

// Load GeoJSON routes into the graph
void loadRoutesFromGeoJSON(const std::string &filepath, TFTFGraph &graph) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open GeoJSON file!" << std::endl;
        return;
    }

    json geojson;
    file >> geojson;

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

int main() {
    const std::string geojsonInput = "routes.geojson";
    const std::string outputFilename = "graph.json";
    const float transferRangeMeters = 300.0f;

    TFTFGraph graph;
    std::cout << "Loading routes from " << geojsonInput << "...\n";
    loadRoutesFromGeoJSON(geojsonInput, graph);

    std::cout << "Creating transfers within " << transferRangeMeters << " meters...\n";
    graph.createTransfersFromCoordinates(transferRangeMeters);

    std::cout << "Saving serialized graph to " << outputFilename << "...\n";
    saveGraphToDisk(graph, outputFilename);

    std::cout << "Done.\n";
    return 0;
}

// g++ -std=c++17 graph_export.cpp TFTFGraph/TFTFGraph.cpp TFTFGraph/Helpers/helpers.cpp -o graph_export

