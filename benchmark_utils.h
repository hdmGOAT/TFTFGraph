#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#include <string>
#include "./TFTFGraph/TFTFGraph.h"
#include <map>
#include "./algorithms/node.h"

enum TestCategory {
    SAME_ROUTE,
    DIFFERENT_ROUTES,
    END_TO_END,
    LOOP_TESTS
};

struct TestResult {
    Coordinate from;
    Coordinate to;
    long long tftf_ms;
    long long dijkstra_ms;
    long long astar_ms;
    bool found_tftf;
    bool found_dijkstra;
    bool found_astar;
    int tftf_routes_taken;    // Routes taken by TFTF
    int dijkstra_routes_taken; // Routes taken by Dijkstra
    int astar_routes_taken;    // Routes taken by A*
};

void initializeCSV(const std::string& filename);
void saveTestResult(const TestResult& result, const std::string& filename);
void runTestCategory(TFTFGraph& network, 
                    std::map<Node, std::vector<std::pair<Node, double>>>& nodeGraph,
                    TestCategory category,
                    const std::string& filename,
                    int totalTests,
                    const std::vector<Coordinate>& allCoordinates);

std::pair<Coordinate, Coordinate> generateTestCoordinates(
    TestCategory category,
    const std::vector<Coordinate>& allCoordinates,
    const TFTFGraph& network);

#endif // BENCHMARK_UTILS_H 