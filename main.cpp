#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"
#include <fstream>

#include "json.hpp"

using json = nlohmann::json;
void loadRoutesFromGeoJSON(const std::string& filepath, TFTFGraph& graph) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open GeoJSON file!" << std::endl;
        return;
    }

    json geojson;
    file >> geojson;

    const auto& features = geojson["features"];
    int routeId = 0;

    for (const auto& feature : features) {
        if (feature["geometry"]["type"] != "LineString") continue;

        std::vector<Coordinate> routePath;
        const auto& coords = feature["geometry"]["coordinates"];

        for (const auto& coord : coords) {
            double lon = coord[0];
            double lat = coord[1];
            routePath.emplace_back(Coordinate{lat, lon});
        }

        std::string routeName = "Route_" + std::to_string(routeId);
        if (feature.contains("properties") && feature["properties"].contains("name")) {
            routeName = feature["properties"]["name"];
        }

        bool isLoop = routePath.front() == routePath.back();

        

        graph.addRoute(routeId, routeName);
        graph.setRoutePath(routeId, routePath);
        routeId++;
    }
}


int main() {
    TFTFGraph jeepneyNetwork;
    loadRoutesFromGeoJSON("routes.geojson", jeepneyNetwork);
    jeepneyNetwork.createTransfersFromCoordinates(100.0f, 10.0f); 
    return 0;
}

//RUN THIS CODE

//g++ -std=c++17 main.cpp TFTFGraph/TFTFGraph.cpp TFTFGraph/Helpers/helpers.cpp -o main