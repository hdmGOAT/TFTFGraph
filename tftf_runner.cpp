#include "TFTFGraph/TFTFGraph.h"
#include "json.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <set>
#include <unordered_map>
#include "./TFTFGraph/Helpers/helpers.h"
#include "./algorithms/node.h"
#include "./algorithms/astar/astar.h"
#include "./algorithms/djikstra/djikstra.h"

using json = nlohmann::json;

TFTFGraph graph;

std::vector<Coordinate> flattenTftfPath(const std::vector<RoutePathInstruction> &instructions)
{
    std::vector<Coordinate> coords;
    for (const auto &instruction : instructions)
    {
        for (const auto &coord : instruction.path)
        {
            if (!coords.empty())
            {
                const auto &last = coords.back();
                if (last.latitude == coord.latitude && last.longitude == coord.longitude)
                    continue;
            }
            coords.push_back(coord);
        }
    }
    return coords;
}

std::vector<Coordinate> nodePathToCoordinates(const std::vector<Node> &nodePath)
{
    std::vector<Coordinate> coords;
    coords.reserve(nodePath.size());
    for (const auto &node : nodePath)
    {
        coords.push_back({node.lat, node.lon});
    }
    return coords;
}

double calculatePolylineDistanceMeters(const std::vector<Coordinate> &coords)
{
    double total = 0.0;
    for (size_t i = 1; i < coords.size(); ++i)
    {
        total += haversine(coords[i - 1], coords[i]);
    }
    return total;
}

json makeLineFeature(const std::vector<Coordinate> &coords, const std::string &label, const std::string &color)
{
    json coordinates = json::array();
    for (const auto &coord : coords)
    {
        coordinates.push_back({coord.longitude, coord.latitude});
    }

    return {
        {"type", "Feature"},
        {"properties", {{"label", label}, {"color", color}}},
        {"geometry", {{"type", "LineString"}, {"coordinates", coordinates}}}};
}

json makeSegmentFeature(const RoutePathInstruction &instruction, int segmentIndex, const std::string &algorithm)
{
    json coordinates = json::array();
    for (const auto &coord : instruction.path)
    {
        coordinates.push_back({coord.longitude, coord.latitude});
    }

    return {
        {"type", "Feature"},
        {"properties", {
            {"algorithm", algorithm},
            {"label", instruction.routeName},
            {"route_id", instruction.routeId},
            {"segment_index", segmentIndex}
        }},
        {"geometry", {{"type", "LineString"}, {"coordinates", coordinates}}}};
}

json buildSegmentFeatureCollection(const std::vector<RoutePathInstruction> &instructions, const std::string &algorithm)
{
    json features = json::array();
    for (size_t i = 0; i < instructions.size(); ++i)
    {
        if (instructions[i].path.size() < 2)
            continue;
        features.push_back(makeSegmentFeature(instructions[i], static_cast<int>(i), algorithm));
    }

    return {
        {"type", "FeatureCollection"},
        {"features", features}};
}

json buildComparisonGeoJSON(
    const std::vector<RoutePathInstruction> &tftfInstructions,
    const std::vector<RoutePathInstruction> &traditionalInstructions,
    const std::string &traditionalAlgorithm)
{
    json features = json::array();

    json tftfFC = buildSegmentFeatureCollection(tftfInstructions, "tftf");
    for (const auto &feature : tftfFC["features"])
    {
        features.push_back(feature);
    }

    json traditionalFC = buildSegmentFeatureCollection(traditionalInstructions, traditionalAlgorithm);
    for (const auto &feature : traditionalFC["features"])
    {
        features.push_back(feature);
    }

    return {
        {"type", "FeatureCollection"},
        {"features", features}};
}

bool loadTraditionalGraph(
    const std::string &geojsonPath,
    std::map<Node, std::vector<std::pair<Node, double>>> &nodeGraph,
    std::unordered_map<int, std::string> &routeNames)
{
    std::ifstream file(geojsonPath);
    if (!file.is_open())
    {
        return false;
    }

    json geojson;
    file >> geojson;
    nodeGraph.clear();
    routeNames.clear();

    int routeId = 0;
    for (const auto &feature : geojson["features"])
    {
        if (feature["geometry"]["type"] != "LineString")
            continue;

        std::string routeName = "Route_" + std::to_string(routeId);
        if (feature.contains("properties") && feature["properties"].contains("name"))
        {
            routeName = feature["properties"]["name"].get<std::string>();
        }
        routeNames[routeId] = routeName;
        routeId++;
    }

    geojsonToNodeGraph(nodeGraph, geojson);
    return true;
}

std::vector<Node> runTraditionalPath(const std::string &algorithm, const Coordinate &start, const Coordinate &end, const std::map<Node, std::vector<std::pair<Node, double>>> &graph)
{
    Node origin(start.latitude, start.longitude);
    Node destination(end.latitude, end.longitude);

    if (algorithm == "dijkstra")
    {
        return dijkstra_geojson("", origin, destination, graph);
    }

    return astar_geojson("", origin, destination, graph);
}

int countTraditionalRoutes(const std::vector<Node> &path)
{
    std::set<int> uniqueRoutes;
    for (const auto &node : path)
    {
        if (node.routeId >= 0)
        {
            uniqueRoutes.insert(node.routeId);
        }
    }
    return static_cast<int>(uniqueRoutes.size());
}

json serializeCoordinates(const std::vector<Coordinate> &coords)
{
    json out = json::array();
    for (const auto &coord : coords)
    {
        out.push_back({
            {"lat", coord.latitude},
            {"lon", coord.longitude}});
    }
    return out;
}

std::vector<RoutePathInstruction> buildTraditionalInstructions(
    const std::vector<Node> &path,
    const std::unordered_map<int, std::string> &routeNames)
{
    std::vector<RoutePathInstruction> instructions;
    RoutePathInstruction current;
    bool hasOpenSegment = false;

    for (const auto &node : path)
    {
        if (node.routeId < 0)
            continue;

        Coordinate coord = {node.lat, node.lon};

        if (!hasOpenSegment || current.routeId != node.routeId)
        {
            if (hasOpenSegment && current.path.size() >= 2)
            {
                instructions.push_back(current);
            }

            std::string routeName = "Route_" + std::to_string(node.routeId);
            auto routeNameIt = routeNames.find(node.routeId);
            if (routeNameIt != routeNames.end())
            {
                routeName = routeNameIt->second;
            }

            current = {node.routeId, routeName, {coord}};
            hasOpenSegment = true;
            continue;
        }

        if (current.path.empty() ||
            current.path.back().latitude != coord.latitude ||
            current.path.back().longitude != coord.longitude)
        {
            current.path.push_back(coord);
        }
    }

    if (hasOpenSegment && current.path.size() >= 2)
    {
        instructions.push_back(current);
    }

    return instructions;
}

int main()
{
    std::string line;
    std::map<Node, std::vector<std::pair<Node, double>>> traditionalGraph;
    std::unordered_map<int, std::string> traditionalRouteNames;
    std::string loadedTraditionalGeojsonPath;
    std::string loadedTftfGraphPath = "data/graph.json";

    // Load the graph from a JSON file
    graph = loadGraphFromDisk(loadedTftfGraphPath);

    while (std::getline(std::cin, line))
    {
        try
        {
            json req = json::parse(line);
            Coordinate start = {req["start"]["lat"], req["start"]["lon"]};
            Coordinate end = {req["end"]["lat"], req["end"]["lon"]};
            bool includeTraditional = req.value("include_traditional", false) || req.value("compare", false);
            std::string traditionalAlgorithm = req.value("traditional_algorithm", "astar");
            std::string traditionalGeojsonPath = req.value("traditional_geojson", "data/geojson/allRoutes.geojson");
            bool includeComparisonGeojson = req.value("include_comparison_geojson", includeTraditional);
            std::string tftfGraphPath = req.value("tftf_graph_json", loadedTftfGraphPath);

            if (tftfGraphPath != loadedTftfGraphPath)
            {
                graph = loadGraphFromDisk(tftfGraphPath);
                loadedTftfGraphPath = tftfGraphPath;
            }

            if (graph.getRoutes().empty())
            {
                throw std::runtime_error("Loaded TFTF graph has no routes. Check tftf_graph_json path.");
            }

            if (traditionalAlgorithm != "astar" && traditionalAlgorithm != "dijkstra")
            {
                traditionalAlgorithm = "astar";
            }

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
            std::vector<Coordinate> tftfCoordinates = flattenTftfPath(instructions);
            int tftfRouteCount = static_cast<int>(instructions.size());
            double tftfDistanceMeters = calculatePolylineDistanceMeters(tftfCoordinates);

            // Build response
            json response = {
                {"total_fare", std::ceil(fare)},
                {"geojson", geojson},
                {"routes", json::array()},
                {"tftf", {
                    {"found", true},
                    {"route_count", tftfRouteCount},
                    {"distance_m", std::round(tftfDistanceMeters * 100.0) / 100.0},
                    {"geojson", makeLineFeature(tftfCoordinates, "TFTF", "#0066ff")},
                    {"route_segments", json::array()},
                    {"segment_geojson", buildSegmentFeatureCollection(instructions, "tftf")}
                }}};

            for (const auto &instruction : instructions)
            {
                response["tftf"]["route_segments"].push_back({
                    {"route_id", instruction.routeId},
                    {"route_name", instruction.routeName},
                    {"point_count", instruction.path.size()},
                    {"coordinates", serializeCoordinates(instruction.path)}
                });
            }

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

            if (includeTraditional)
            {
                if (traditionalGeojsonPath != loadedTraditionalGeojsonPath)
                {
                    if (!loadTraditionalGraph(traditionalGeojsonPath, traditionalGraph, traditionalRouteNames))
                    {
                        throw std::runtime_error("Unable to open traditional_geojson file: " + traditionalGeojsonPath);
                    }
                    loadedTraditionalGeojsonPath = traditionalGeojsonPath;
                }

                auto traditionalStart = std::chrono::high_resolution_clock::now();
                std::vector<Node> traditionalPath = runTraditionalPath(traditionalAlgorithm, start, end, traditionalGraph);
                auto traditionalEnd = std::chrono::high_resolution_clock::now();
                long long traditionalMs = std::chrono::duration_cast<std::chrono::milliseconds>(traditionalEnd - traditionalStart).count();

                std::vector<Coordinate> traditionalCoordinates = nodePathToCoordinates(traditionalPath);
                double traditionalDistanceMeters = calculatePolylineDistanceMeters(traditionalCoordinates);
                std::vector<RoutePathInstruction> traditionalInstructions = buildTraditionalInstructions(traditionalPath, traditionalRouteNames);

                response["traditional"] = {
                    {"algorithm", traditionalAlgorithm},
                    {"runtime_ms", traditionalMs},
                    {"found", !traditionalPath.empty()},
                    {"route_count", countTraditionalRoutes(traditionalPath)},
                    {"distance_m", std::round(traditionalDistanceMeters * 100.0) / 100.0},
                    {"geojson", makeLineFeature(traditionalCoordinates, traditionalAlgorithm, "#ff6600")},
                    {"route_segments", json::array()},
                    {"segment_geojson", buildSegmentFeatureCollection(traditionalInstructions, traditionalAlgorithm)}
                };

                for (const auto &instruction : traditionalInstructions)
                {
                    response["traditional"]["route_segments"].push_back({
                        {"route_id", instruction.routeId},
                        {"route_name", instruction.routeName},
                        {"point_count", instruction.path.size()},
                        {"coordinates", serializeCoordinates(instruction.path)}
                    });
                }

                if (includeComparisonGeojson)
                {
                    response["comparison_geojson"] = buildComparisonGeoJSON(instructions, traditionalInstructions, traditionalAlgorithm);
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
