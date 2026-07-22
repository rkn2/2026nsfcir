# Water-Surface Slope ($\beta$) Sensitivity -- Downtown Montpelier, VT
Preliminary evidence for an NSF CPS-CIR proposal on flood-duration/depth modeling in Montpelier, VT: brackets the reach-scale water-surface slope $\beta$ in `WSE(x) = 158.850 m NAVD88 + \beta \cdot x` (Root Research_v2.tex, Subtask 1.2) from **public evidence only** -- no ICEYE data used -- in advance of the 2023 ICEYE depth-field calibration that the proposal describes. All numbers below are computed by `analysis/slope_sensitivity.py` from live USGS NWIS, USGS 3DEP, and OSM Overpass/Nominatim queries; see that script for exact methods.
## Method summary
1. Same DEM clip, gauge datum, and July 2023 peak stage as `analysis/stage_to_dem.py` (WSE at the gauge = **158.850 m NAVD88**).
2. Along-channel upstream distance $x$ for every DEM cell, landmark, and OSM building centroid: OSM waterway ways tagged `waterway=river` with name containing "Winooski" (captures both the Winooski River and the North Branch Winooski River) fetched via Overpass for a padded bbox, assembled into an undirected graph (vertices snapped within 5 m to merge near-coincident endpoints at confluences/way splits), and a single-source Dijkstra shortest-path from the graph node nearest the gauge gives along-channel-network distance to every other channel point -- this correctly handles the Y-shaped Winooski / North Branch confluence. Every DEM cell/landmark/building is assigned the distance value of its nearest channel point (`scipy.spatial.cKDTree`), per the task's "nearest-channel-point distance" specification. See the script's module docstring for the full method and its approximations.
   - Channel graph: **307 nodes, 306 edges**; gauge snapped to nearest channel node **36.9 m** away; **0** of 307 nodes were not reachable from the gauge's graph component (excluded from the distance field).
   - Densified channel point cloud: **4,311 points** at ~5 m spacing.
3. $\beta$ swept from 0.0 to 1.5 m/km in 0.1 m/km steps (`WSE(x) = 158.850 + \beta \cdot x`, $x$ in km). For each $\beta$: the connected inundation mask (same flood-fill-from-river-proxy logic as `stage_to_dem.py`, generalized to a spatially-varying WSE), connected inundated area, number of OSM building centroids whose grid cell falls in the connected mask, and depth (`WSE(x_landmark) - elevation`) at the three sanity-check landmarks.
## Beta sweep
| beta (m/km) | connected area (acres) | OSM buildings inundated | VT State House depth (m) | City Hall depth (m) | Confluence Park depth (m) |
|---|---|---|---|---|---|
| 0.0 | 92.4 | 31 | -9.86 | -1.48 | -0.82 |
| 0.1 | 103.0 | 41 | -9.74 | -1.30 | -0.66 |
| 0.2 | 112.7 | 48 | -9.61 | -1.13 | -0.49 |
| 0.3 | 124.1 | 67 | -9.49 | -0.96 | -0.33 |
| 0.4 | 135.4 | 88 | -9.37 | -0.79 | -0.17 |
| 0.5 | 147.3 | 112 | -9.24 | -0.61 | -0.01 |
| 0.6 | 156.1 | 130 | -9.12 | -0.44 | +0.15 |
| 0.7 | 166.8 | 146 | -9.00 | -0.27 | +0.31 |
| 0.8 | 179.4 | 188 | -8.87 | -0.10 | +0.48 |
| 0.9 | 191.7 | 234 | -8.75 | +0.08 | +0.64 |
| 1.0 | 204.3 | 266 | -8.63 | +0.25 | +0.80 |
| 1.1 | 219.2 | 305 | -8.50 | +0.42 | +0.96 |
| 1.2 | 234.2 | 336 | -8.38 | +0.59 | +1.12 |
| 1.3 | 248.2 | 374 | -8.26 | +0.77 | +1.28 |
| 1.4 | 260.0 | 407 | -8.13 | +0.94 | +1.45 |
| 1.5 | 269.7 | 429 | -8.01 | +1.11 | +1.61 |

## Admissible beta interval
- **Admissible interval: beta in [0.9, >= 1.5 (top of swept range)] m/km.**
  - Lower bound 0.9 m/km: the smallest swept beta at which BOTH Main St / downtown-core landmarks (Montpelier City Hall AND Confluence Park, both documented flooded in July 2023) are predicted inundated (depth >= 0).
  - Upper bound: the VT State House (documented NOT flooded) remains dry across the ENTIRE swept range (0.0-1.5 m/km); its own implied threshold from the flat-surface residual (8.00 m/km) is far above the swept range, so this evidence only right-censors the interval at the top of the sweep rather than pinning a tighter upper bound.
- At beta = 0.9 m/km (lower bound): connected area **191.7 acres**, **234 OSM buildings** inundated (vs. 92.4 acres / 31 buildings at beta = 0, the flat-surface baseline).
- At beta = 1.5 m/km (top of sweep): connected area **269.7 acres**, **429 OSM buildings** inundated.

## Comparison with the flat-surface-residual inference
`analysis/stage_dem.md` inferred a missing water-surface slope of roughly **0.6-1.0 m/km** from how far City Hall and Confluence Park sat above the flat WSE, divided by their straight-line distance from the gauge. This sweep's admissible lower bound (0.9 m/km) uses along-channel (not straight-line) distance and a discrete 0.1 m/km grid, so exact agreement is not expected; the two are broadly consistent (0.9 m/km vs. the 0.6-1.0 m/km prior inference). The along-channel distances used here are longer than the straight-line distances used in the flat-surface-residual calculation (channel distance >= straight-line distance always), which mechanically pushes the beta needed to flood a given landmark downward relative to a straight-line-distance estimate -- one identifiable source of any gap between the two.

## Sanity checks
- Connected area is monotonically non-decreasing in beta. OK.
- OSM building inundation count is monotonically non-decreasing in beta. OK.
- VT State House depth is monotonically non-decreasing in beta. OK.
- Montpelier City Hall depth is monotonically non-decreasing in beta. OK.
- Confluence Park depth is monotonically non-decreasing in beta. OK.

## Caveats
- **Both rivers (Winooski main stem and North Branch) are assigned the SAME beta.** This is a single reach-scale slope parameter, not two independently calibrated slopes; the North Branch's true flood-surface slope could differ from the main stem's, especially near the confluence itself.
- **Planar-per-reach surface, not a hydraulically routed profile.** WSE(x) is linear in along-channel distance; it does not represent backwater curves, local constrictions (bridges, channel narrowing), or the confluence's actual hydraulic behavior (e.g. a true water surface is not required to have equal slope on both branches meeting at a confluence).
- **Along-channel distance is a network graph shortest-path distance from OSM waterway centerline geometry, not a traced lowest-elevation thalweg path through the DEM and not a distance computed from a validated hydrography dataset (e.g. NHD).** OSM waterway digitization quality/completeness was not independently audited beyond the diagnostics reported above (node/edge counts, unreachable-node count).
- **Signed distance is approximated as unsigned.** Dijkstra distance is always >= 0; the small sliver of the bbox downstream (west) of the gauge is therefore treated as if it were upstream, which would very slightly overstate WSE there. The gauge sits within ~37 m of the bbox's west edge (see stage_to_dem.py's bbox), so this affects a negligible fraction of the domain, but this is a documented simplification specific to this domain's geometry, not a general solution for a gauge in the interior of a bbox.
- **Only two control landmarks define the admissible interval's lower bound (City Hall, Confluence Park) and one defines the upper bound (State House).** These are the same three point-sample sanity checks used in `stage_to_dem.md`, each subject to the same point-sampling caveats described there (a few meters of horizontal offset can flip a marginal point wet/dry).
- **Inundated-building counts use every OSM-tagged `building=*` footprint in the clip**, not the proposal's 271-building curated downtown masonry dataset -- a building is counted "inundated" if its centroid's DEM grid cell falls in the connected mask, a coarse point-in-cell test, not a building-footprint overlap test.
- **Connectivity is a cheap 4-connected flood-fill from an elevation-percentile channel proxy** (same heuristic and caveats as `stage_to_dem.py`) -- not a validated hydrography-derived channel network or a true 2D flood-routing connectivity solve.
- **This sweep uses NO ICEYE data.** It is deliberately built from public evidence only, to bracket beta ahead of the proposal's planned ICEYE-based calibration (Subtask 1.2); it is not a substitute for that calibration.
- **DEM is bare-earth terrain**, USGS stage data are provisional, and the DEM/geocoding/Overpass caveats in `stage_to_dem.md` and `subbasin_sensitivity.md` apply equally here (not repeated in full).
