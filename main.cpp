#include <iostream>
#include <vector>
#include "./TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"
#include <fstream>

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

void testRoute(TFTFGraph network, Coordinate from, std::string fromName, Coordinate to, std::string toName, float transferMeters)
{
    network.createTransfersFromCoordinates(transferMeters);
    std::cout << "Transfer range: " << transferMeters << " meters" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "From: " << fromName << " (" << from.latitude << ", " << from.longitude << ")" << std::endl;
    std::cout << "To: " << toName << " (" << to.latitude << ", " << to.longitude << ")" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    network.calculateRouteFromCoordinates(from, to);
    std::cout << "--------------------------------------------------" << std::endl;
}

int main()
{
    TFTFGraph jeepneyNetwork;
    loadRoutesFromGeoJSON("routes.geojson", jeepneyNetwork);

    float transferRange = 300.5f;

     testRoute(jeepneyNetwork, {8.508810, 124.648270}, "Bonbon", {8.511330, 124.624290}, "Westbound Bulua Terminal", transferRange);
    // testRoute(jeepneyNetwork, {8.50881, 124.64827}, "Bonbon", {8.482906, 124.646094}, "Velez Mogchs", transferRange);
    // testRoute(jeepneyNetwork, {8.504775, 124.642954}, "Kauswagan City Engineer", {8.484763, 124.655977}, "USTP", transferRange);
    testRoute(jeepneyNetwork, {8.487358, 124.629950}, "Patag Camp Evangelista", {8.484763, 124.655977}, "USTP", transferRange);

    // jeepneyNetwork.createTransfersFromCoordinates(300.0f);
    // // Bonbon - Westbound Bulua Terminal
    // std::cout << "Bonbon - Westbound Bulua Terminal" << std::endl;
    // jeepneyNetwork.calculateRouteFromCoordinates({8.50881, 124.64827}, {8.51133, 124.62429}, 10);

    std::cout << std::endl;

    // // Bonbon - Velez Mogchs
    // std::cout << "Westbound Bulua Terminal - Velez Mogchs" << std::endl;
    // jeepneyNetwork.calculateRouteFromCoordinates({8.50881, 124.64827}, {8.482906, 124.646094}, 10);

    // std::cout << std::endl;

    // // Kauswagan City Engineer - USTP
    // std::cout << "Kauswagan City Engineer - USTP" << std::endl;
    // jeepneyNetwork.calculateRouteFromCoordinates({8.504775, 124.642954}, {8.484763, 124.655977}, 10);
    return 0;
}

// RUN THIS CODE

// g++ -std=c++17 main.cpp TFTFGraph/TFTFGraph.cpp TFTFGraph/Helpers/helpers.cpp -o main