#include "TFTFGraph/TFTFGraph.h"
#include "json.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include "./TFTFGraph/Helpers/helpers.h"

using json = nlohmann::json;

TFTFGraph graph;

int main()
{
    std::string line;

    // Load the graph from a JSON file
    graph = loadGraphFromDisk("data/graph.json");

    while (std::getline(std::cin, line))
    {
        try
        {
            json req = json::parse(line);
            Coordinate start = {req["start"]["lat"], req["start"]["lon"]};
            Coordinate end = {req["end"]["lat"], req["end"]["lon"]};

            // Compute route path
            std::vector<TFTFEdge> path = graph.calculateRouteFromCoordinates(start, end);

            if (path.empty())
            {
                std::cerr << "{\"error\": \"No valid path found\"}" << std::endl;
                continue;
            }

            // Generate step-by-step instructions
            std::vector<RoutePathInstruction> instructions = graph.constructRoutePathInstructions(path);

            // Compute total fare
            double fare = graph.calculateFareFromInstructions(instructions);

            // Build GeoJSON output
            json geojson = generateRoutePathGeoJSON(instructions, start, end);

            // Build response
            json response = {
                {"total_fare", std::ceil(fare)},
                {"geojson", geojson},
                {"routes", json::array()}};

            for (const auto &instr : instructions)
            {
                double routeDistance = 0.0;
                const auto &path = instr.path;

                // Calculate distance for this route
                for (size_t i = 1; i < path.size(); ++i)
                {
                    routeDistance += haversine(path[i - 1], path[i]);
                }

                // Calculate fare for this route
                double routeFare = 12.0; // Base fare
                float distanceInKm = routeDistance / 1000.0f;

                // Debug print
                std::cerr << "[DEBUG] routeId: " << instr.routeId << ", name: " << instr.routeName << ", distance_km: " << distanceInKm << std::endl;

                if (distanceInKm > 4.0f) // FREE_KM
                {
                    float extraDistance = distanceInKm - 4.0f;
                    int extraBlocks = static_cast<int>(std::ceil(extraDistance / 4.0f));
                    routeFare += extraBlocks * 1.5f; // FARE_PER_4KM
                }

                // Only include real routes (distance > 0.1 km)
                if (distanceInKm > 0.3f)
                {
                    json routeInfo = {
                        {"name", instr.routeName},
                        {"fare", std::ceil(routeFare)},
                        {"distance_km", std::round(distanceInKm * 100.0) / 100.0}};
                    response["routes"].push_back(routeInfo);
                }
            }

            // Return JSON result
            std::cout << response.dump() << std::endl;
            std::cout.flush();
        }
        catch (const std::exception &e)
        {
            std::cerr << "{\"error\": \"" << e.what() << "\"}" << std::endl;
        }
    }

    return 0;
}
