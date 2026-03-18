CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

BIN_DIR := bin

COMMON_SOURCES := \
	TFTFGraph/TFTFGraph.cpp \
	TFTFGraph/Helpers/helpers.cpp \
	algorithms/astar/astar.cpp \
	algorithms/djikstra/djikstra.cpp \
	algorithms/node.cpp

RUNNER_SOURCES := tftf_runner.cpp $(COMMON_SOURCES)
BENCH_SOURCES := main.cpp benchmark_utils.cpp $(COMMON_SOURCES)

RUNNER_BIN := $(BIN_DIR)/tftf_runner
BENCH_BIN := $(BIN_DIR)/benchmark

START_LAT ?= 8.50895
START_LON ?= 124.64792
END_LAT ?= 8.47804
END_LON ?= 124.64303
TRAD_ALGO ?= dijkstra
TRAD_GEOJSON ?= routes.geojson
TFTF_GRAPH_JSON ?= data/graph.json

.PHONY: all runner benchmark run compare run-benchmark clean

all: runner

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

runner: $(RUNNER_BIN)

benchmark: $(BENCH_BIN)

$(RUNNER_BIN): $(RUNNER_SOURCES) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BENCH_BIN): $(BENCH_SOURCES) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

run: $(RUNNER_BIN)
	./$(RUNNER_BIN)

compare: $(RUNNER_BIN)
	printf '{"start":{"lat":$(START_LAT),"lon":$(START_LON)},"end":{"lat":$(END_LAT),"lon":$(END_LON)},"compare":true,"traditional_algorithm":"$(TRAD_ALGO)","traditional_geojson":"$(TRAD_GEOJSON)","tftf_graph_json":"$(TFTF_GRAPH_JSON)","include_comparison_geojson":true}\n' | ./$(RUNNER_BIN)

run-benchmark: $(BENCH_BIN)
	./$(BENCH_BIN)

clean:
	rm -rf $(BIN_DIR)