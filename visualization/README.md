# Visualization

This folder contains the side-by-side map viewer for comparing TFTF vs a traditional algorithm path.

## Files

- `index.html` — Leaflet viewer with:
  - left map for TFTF,
  - right map for traditional algorithm,
  - north arrow on each map,
  - scale bar on each map.
- `output/` — generated JSON outputs from the comparison runner.

## Generate comparison output

From project root:

```bash
make compare-save
```

This writes:

- `visualization/output/compare_response.json`

You can override coordinates and source files via make vars:

```bash
make compare-save \
  START_LAT=8.50089 START_LON=124.616 \
  END_LAT=8.29258 END_LON=124.52 \
  TRAD_ALGO=dijkstra \
  TRAD_GEOJSON=data/geojson/allRoutes.geojson \
  TFTF_GRAPH_JSON=data/graph.json
```

## Open the viewer

Open `visualization/index.html` in your browser, then choose a JSON file:

- `visualization/output/compare_response.json`, or
- any JSON containing `comparison_geojson`, or
- a plain GeoJSON FeatureCollection with TFTF + traditional line features.