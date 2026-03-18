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

.PHONY: all runner benchmark run run-benchmark clean

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

run-benchmark: $(BENCH_BIN)
	./$(BENCH_BIN)

clean:
	rm -rf $(BIN_DIR)