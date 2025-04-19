#include <iostream>
#include "json.hpp"
#include "TFTFGraph/TFTFGraph.h"
#include "./TFTFGraph/Helpers/helpers.h"
#include <fstream>
#include <iomanip>

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

int main(int argc, char *argv[])
{
    // SAMPLE
    // ./algo 8.487358 124.629950 "Patag Camp Evangelista" 8.484763 124.655977 "USTP" 10 10
    if (argc != 10)
    {
        std::cerr << "Usage: program <fromLat> <fromLon> <fromName> <toLat> <toLon> <toName> <transferMeters> <hour> <routesPath>\n";
        return 1;
    }

    double fromLat = std::atof(argv[1]);
    double fromLon = std::atof(argv[2]);
    std::string fromName = argv[3];
    double toLat = std::atof(argv[4]);
    double toLon = std::atof(argv[5]);
    std::string toName = argv[6];
    float transferMeters = std::stof(argv[7]);
    int hour = std::stoi(argv[8]);
    std::string routesPath = argv[9];

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Parsed arguments:\n";
    std::cout << "From: " << fromName << " (" << fromLat << ", " << fromLon << ")\n";
    std::cout << "To: " << toName << " (" << toLat << ", " << toLon << ")\n";
    std::cout << "Transfer range: " << transferMeters << " meters\n";
    std::cout << "Hour of operation: " << hour << ":00\n";
    std::cout << "GeoJSON path: " << routesPath << "\n";

    // Create the graph and coordinates
    TFTFGraph jeepneyNetwork;
    loadRoutesFromGeoJSON(routesPath, jeepneyNetwork);
    jeepneyNetwork.createTransfersFromCoordinates(transferMeters);

    Coordinate from = {fromLat, fromLon};
    Coordinate to = {toLat, toLon};

    std::vector<TFTFEdge> path = jeepneyNetwork.calculateRouteFromCoordinates(from, to, hour);

    nlohmann::json response;

    response["From"] = {{"Origin Name", fromName}, {"Lattitude", from.latitude}, {"Longitude", from.longitude}};
    response["To"] = {{"Destination Name", toName}, {"Lattitude", to.latitude}, {"Longitude", to.longitude}};
    response["Transfer Range"] = transferMeters;

    std::cout << path.size() << " edges\n";
    nlohmann::json routeArray = nlohmann::json::array();
    for (int i = 0; i < path.size(); i++)
    {
        nlohmann::json edgeInfo = {
            {"Route ID", path[i].destinationRoute},
            {"Route Name", path[i].destinationRouteName},
            {"Transfer Cost", path[i].transferCost},
            {"Entry Coordinate", {path[i].entryCoord.latitude, path[i].entryCoord.longitude}},
            {"Exit Coordinate", {path[i].exitCoord.latitude, path[i].exitCoord.longitude}}};
        routeArray.push_back(edgeInfo);
    }
    response["Routes"] = routeArray;

    std::cout << response.dump();

    return 0;
}