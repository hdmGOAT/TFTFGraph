#include "TFTFGraph/TFTFGraph.h"
#include "json.hpp"
#include <iostream>
#include <string>
#include <sstream>

using json = nlohmann::json;

TFTFGraph graph;

int main() {
    std::string line;

    // Load the graph from a JSON file
    graph = loadGraphFromDisk("data/graph.json");

    while (std::getline(std::cin, line)) {
        try {
            json req = json::parse(line);
            Coordinate start = {req["start"]["lat"], req["start"]["lon"]};
            Coordinate end = {req["end"]["lat"], req["end"]["lon"]};

            // Compute route path
            std::vector<TFTFEdge> path = graph.calculateRouteFromCoordinates(start, end);

            if (path.empty()) {
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
                {"fare", std::ceil(fare)},
                {"geojson", geojson},
                {"route_names", json::array()}
            };

            for (const auto& instr : instructions) {
                response["route_names"].push_back(instr.routeName);
            }

            // Return JSON result
            std::cout << response.dump() << std::endl;
            std::cout.flush();

        } catch (const std::exception& e) {
            std::cerr << "{\"error\": \"" << e.what() << "\"}" << std::endl;
        }
    }

    return 0;
}
