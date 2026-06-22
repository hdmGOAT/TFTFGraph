# TFTFGraph: Transfer-Flexible Transit Fare Graph for Jeepney Routing

Official code for our published research on transfer-aware, fare-aware routing in flexible urban jeepney systems.

`Status: Under peer review` • `Language: C++17` • `Build: Makefile` • `Algorithms: TFTF, Dijkstra, A*` • `Domain: Public Utility Jeepney (PUJ)`

Jeepney networks are flexible, overlapping, and not strictly timetable-based, so classic shortest-path methods can produce commuter-unfriendly routes with too many transfers. TFTFGraph models each route as a node and transfer opportunities as coordinate-level edges, then applies lexicographic optimization (fewest transfers first, lowest fare second) to generate more practical paths.

This repository also includes baseline methods for comparison: **Dijkstra** and **A\*** on a traditional node graph.

## Why this research

- **Problem:** Informal transit routing is poorly captured by stop-level, schedule-centric assumptions.
- **Goal:** Produce realistic jeepney paths that reduce unnecessary transfers while remaining fare-aware.
- **Approach:** Compress the network via route-as-node representation and compute transfer edges from spatial proximity.
- **Outcome:** Smaller graph representation with stable query behavior and lower average route-hops in inter-route cases.

## Core idea

TFTFGraph operationalizes three principles:

1. **Route-level abstraction:** each jeepney route is treated as a graph node.
2. **Coordinate transfer modeling:** feasible transfers are formed using distance-based coordinate checks.
3. **Lexicographic path selection:** minimize transfers first, then minimize fare among transfer-tied candidates.

## Technical details

- **Graph model:** A TFTF graph $G=(V,E)$ where each $v \in V$ is an entire jeepney route, and each $e \in E$ is a coordinate-level transfer opportunity between two routes.
- **Transfer generation:** Route polylines are densified (approximately 25 m spacing), indexed in a spatial grid, then candidate inter-route transfer pairs are searched within a transfer radius; the minimum-cost feasible pair per route-pair is retained.
- **Optimization objective:** Routing uses lexicographic minimization
	$$
	\min_P \langle T(P), F(P) \rangle
	$$
	where $T(P)$ is transfer count and $F(P)$ is total fare.
- **Fare model:** Total fare combines per-boarding base fare and distance increments, with implementation-level rules for free-distance threshold and per-block add-ons.
- **Baseline comparison:** Traditional node-level Dijkstra and A* are retained as baselines to measure runtime, path structure, and route-hop behavior.

## Key findings (current)

- **Compact full-network representation:** TFTF reports `82` route nodes and `4,093` transfer edges versus `65,244` nodes and `40,716` edges in a traditional node graph setup.
- **Stable query envelope:** TFTF runtime remains bounded and predictable in tested scenarios (reported worst case around `~195 ms`) while enforcing transfer/fare-aware constraints.
- **Lower inter-route hops:** TFTF averages about `2.1` routes in inter-route cases, while baseline shortest-path methods typically require `3–4` routes.
- **Commuter-aligned paths:** By prioritizing transfer minimization before fare tie-breaking, paths are generally less fragmented and operationally closer to real jeepney travel behavior.

## Project layout

```text
.
├── TFTFGraph/                    # Core TFTF graph logic
├── algorithms/                   # Baseline algorithms and node structures
├── comparison/                   # Simple baseline test programs
├── data/
│   ├── geojson/                  # Route datasets
│   └── graph.json                # Serialized TFTF graph
├── output/
│   ├── benchmarks/               # Benchmark CSV outputs
│   ├── routes/                   # Route GeoJSON outputs
│   └── legacy/                   # Archived generated artifacts
├── visualization/                # Side-by-side map comparison UI
├── main.cpp                      # Benchmark runner
├── tftf_runner.cpp               # JSON stdin/stdout runner
├── graph_export.cpp              # GeoJSON -> TFTF graph serialization
└── Makefile
```

## Quick start

```bash
make clean
make runner
make benchmark
```

Run comparison output generation:

```bash
make compare-save
```

Run benchmark suite:

```bash
make run-benchmark
```

## Reproducible outputs

- Comparison JSON: `visualization/output/compare_response.json`
- Benchmark CSVs: `output/benchmarks/`
- Generated route GeoJSON: `output/routes/`

## Notes

- This project is actively being prepared for publication.
- Numbers and claims should be treated as research-in-progress until final paper release.
