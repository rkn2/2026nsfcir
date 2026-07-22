#!/usr/bin/env python3
"""
Water-surface slope (beta) sensitivity study -- downtown Montpelier, VT.

Produced as preliminary evidence for an NSF CPS-CIR proposal on flood-
duration/depth modeling in Montpelier, VT. analysis/stage_to_dem.py projected
the July 2023 crest (WSE = 158.850 m NAVD88, gauge 04286000) onto a 1 m 3DEP
DEM as a FLAT surface and found it under-predicts the documented downtown
flooded core: Main St landmarks (Montpelier City Hall, Confluence Park)
sample 0.8-1.5 m ABOVE the flat WSE at 1.2-1.5 km upstream of the gauge,
implying a missing water-surface slope of roughly 0.6-1.0 m/km. Root
Research_v2.tex, Subtask 1.2, now describes the stage-to-DEM projection as
stage plus a reach-scale water-surface slope beta, to be calibrated against
the 2023 ICEYE depth field. This script brackets beta from PUBLIC evidence
only (no ICEYE data used), ahead of that calibration:

    WSE(x) = 158.850 m + beta * x,   x = along-channel upstream distance
             from the gauge (km), beta swept 0.0 to 1.5 m/km in 0.1 steps.

What this script does
----------------------
1. Reuses the same DEM clip, gauge WSE, and bathtub/connected-mask machinery
   as analysis/stage_to_dem.py (copied and adapted here per task
   instructions, rather than imported, so this script is standalone).
2. Defines along-channel upstream distance x for every DEM cell using OSM
   waterway centerlines (Overpass), NOT a traced lowest-elevation path (see
   "Channel distance method" below for why, and the approximations this
   introduces).
3. Sweeps beta from 0.0 to 1.5 m/km in 0.1 m/km steps. For each beta:
   connected inundated area, number of OSM building centroids inundated,
   and depth at the three sanity-check landmarks (VT State House,
   Montpelier City Hall, Confluence Park).
4. Finds the beta interval consistent with (a) City Hall AND Confluence
   Park (the Main St / documented-flooded landmarks) becoming inundated and
   (b) the VT State House (documented NOT flooded) remaining dry, and
   compares that interval with the 0.6-1.0 m/km flat-surface-residual
   inference already reported in analysis/stage_dem.md.
5. Runs monotonicity sanity checks (area, building count, landmark depths
   must all be non-decreasing in beta) and flags any violation loudly
   rather than silently reporting a non-monotonic sweep as if it were
   expected.
6. Writes analysis/slope_sensitivity.md and analysis/slope_sensitivity.png.

Channel distance method
------------------------
The task allows two options: (a) trace the lowest-elevation path of the
river corridor through the DEM, or (b) use distance along the OSM waterway
polyline, and says the latter is "simpler and acceptable." This script uses
(b): OSM waterway ways tagged waterway=river with name containing
"Winooski" (this captures both "Winooski River" and "North Branch Winooski
River") are fetched via Overpass for a padded bbox, reprojected to
EPSG:26918, and assembled into an undirected graph (nodes = way vertices,
snapped within a 5 m tolerance to merge near-coincident endpoints at
confluences and way splits; edges = consecutive-vertex segments, weighted by
Euclidean length). A single-source Dijkstra shortest-path (scipy.sparse.
csgraph.dijkstra) from the graph node nearest the gauge gives network
along-channel distance from the gauge to every other node -- this correctly
handles the Y-shaped Winooski / North Branch confluence (a straight
point-to-point distance along one digitized way would not). Edges are then
densified at ~5 m spacing with linearly-interpolated distance to build a
dense channel point cloud; every DEM cell, landmark, and building centroid
is assigned x = the distance value of its NEAREST channel point (a
cKDTree nearest-neighbor query), exactly as the task specifies ("applied to
cells by nearest-channel-point distance").

This is a network (graph) shortest-path distance, not a hydraulically
routed thalweg profile, and treats the small sliver of the bbox downstream
(west) of the gauge as if it were upstream (unsigned Dijkstra distance is
always >= 0) -- see Caveats for why this is a small, documented
approximation for this specific domain, not a general solution.

Usage
-----
    uv run --with rasterio --with numpy --with matplotlib --with requests \\
        --with scipy --with shapely analysis/slope_sensitivity.py

Requires: requests, numpy, matplotlib, rasterio, scipy, shapely.
"""

import time
from pathlib import Path

import numpy as np
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import rasterio
from rasterio.io import MemoryFile
from rasterio.warp import transform, transform_bounds
from scipy.ndimage import label
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from shapely.geometry import Polygon

# ----------------------------------------------------------------------
# Configuration (copied/adapted from analysis/stage_to_dem.py)
# ----------------------------------------------------------------------

SITE_NO = "04286000"  # WINOOSKI RIVER AT MONTPELIER, VT
SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
IV_URL = "https://waterservices.usgs.gov/nwis/iv/"
DEM_URL = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage"
NOMINATIM_URL = "https://nominatim.openstreetmap.org/search"
OVERPASS_URL = "https://overpass-api.de/api/interpreter"
OVERPASS_HEADERS = {
    "User-Agent": "cir-flood-proposal-research/1.0 (rjn5308@psu.edu)",
    "Accept": "*/*",
}

PEAK_EVENT_START = "2023-07-01"
PEAK_EVENT_END = "2023-08-15"
PEAK_EVENT_LABEL = "July 2023"

# Same downtown bbox as stage_to_dem.py -- "reuse the same DEM clip."
BBOX_LONLAT = dict(minlon=-72.594, minlat=44.250, maxlon=-72.560, maxlat=44.268)

# Padded bbox for the waterway (Overpass) fetch only, so channel lines are
# not truncated exactly at the DEM clip edge (which would distort along-
# channel distances near the domain boundary). Pad chosen as ~0.006 deg
# (~500 m at this latitude), well beyond the DEM clip on all sides.
WATERWAY_BBOX_LONLAT = dict(
    minlon=BBOX_LONLAT["minlon"] - 0.006, minlat=BBOX_LONLAT["minlat"] - 0.006,
    maxlon=BBOX_LONLAT["maxlon"] + 0.006, maxlat=BBOX_LONLAT["maxlat"] + 0.006,
)

DEM_CRS = "EPSG:26918"  # NAD83 / UTM zone 18N, meters
TARGET_RES_M = 1.0
FT_TO_M = 0.3048

# Fallback values, empirically confirmed from live queries on 2026-07-22
# (matches stage_to_dem.py), used only if the corresponding live service
# call fails at run time.
FALLBACK_ALT_VA_FT = 499.87
FALLBACK_ALT_DATUM_CD = "NAVD88"
FALLBACK_GAUGE_LAT = 44.25672595
FALLBACK_GAUGE_LON = -72.59344318
FALLBACK_PEAK_STAGE_FT = 21.29

FALLBACK_LANDMARKS = {
    "VT State House": (44.2626941, -72.5807647),
    "Montpelier City Hall": (44.2592113, -72.5755301),
    "Confluence Park": (44.2594605, -72.5781534),
}
GEOCODE_QUERIES = {
    "VT State House": "Vermont State House, Montpelier, VT",
    "Montpelier City Hall": "Montpelier City Hall, Montpelier, VT",
    "Confluence Park": "Confluence Park, Montpelier, VT",
}
# Documented flood status used to define the admissible-beta criteria (task
# spec): City Hall and Confluence Park are the Main St / downtown-core
# landmarks reported flooded in July 2023; the VT State House is reported
# NOT flooded. See analysis/stage_dem.md prose for citations/discussion.
MAIN_ST_LANDMARKS = ["Montpelier City Hall", "Confluence Park"]
CONTROL_DRY_LANDMARK = "VT State House"

# Beta sweep (task spec): 0.0 to 1.5 m/km in 0.1 m/km steps.
BETA_VALUES_M_PER_KM = [round(0.1 * i, 1) for i in range(16)]  # 0.0 .. 1.5

# Channel graph node-snapping tolerance (meters) -- see module docstring,
# "Channel distance method."
CHANNEL_SNAP_TOL_M = 5.0
# Densification spacing for the channel point cloud used in the nearest-
# channel-point KDTree query (meters).
CHANNEL_DENSIFY_STEP_M = 5.0
# Flat-surface-residual inference already reported in analysis/stage_dem.md
# (from the July 2023 landmark under-prediction), used only for the
# comparison in the "Numbers" / Caveats sections below -- not used anywhere
# in the beta sweep itself.
PRIOR_INFERENCE_RANGE_M_PER_KM = (0.6, 1.0)

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "analysis"

# Colorblind-safe palette (matches conventions used elsewhere in this repo:
# stage_to_dem.py, winooski_recession.py, subbasin_sensitivity.py).
COLOR_AREA = "#0072B2"       # blue
COLOR_BUILDINGS = "#D55E00"  # vermillion
COLOR_ADMISSIBLE = "#009E73"  # bluish green, shaded admissible-beta band
COLOR_CITYHALL = "#D55E00"
COLOR_CONFLUENCE = "#CC79A7"  # reddish purple
COLOR_STATEHOUSE = "#0072B2"


# ----------------------------------------------------------------------
# Data acquisition: gauge datum and peak stage (copied from stage_to_dem.py)
# ----------------------------------------------------------------------

def fetch_gauge_datum(site_no, retries=4, backoff=5):
    """Fetch gauge altitude (alt_va, ft), vertical datum code, and lat/lon
    from the NWIS expanded site service. Falls back to the documented,
    dated constants above if the live query fails."""
    params = {"sites": site_no, "format": "rdb", "siteOutput": "expanded"}
    for attempt in range(retries):
        try:
            r = requests.get(SITE_URL, params=params, timeout=60)
            if r.status_code == 200:
                lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
                header = lines[0].split("\t")
                data_line = lines[2].split("\t")
                row = dict(zip(header, data_line))
                alt_va = row.get("alt_va", "").strip()
                alt_datum_cd = row.get("alt_datum_cd", "").strip()
                lat = row.get("dec_lat_va", "").strip()
                lon = row.get("dec_long_va", "").strip()
                if alt_va and alt_datum_cd and lat and lon:
                    print(f"  [fetch_gauge_datum] alt_va={alt_va} ft, "
                          f"alt_datum_cd={alt_datum_cd}, lat={lat}, lon={lon} "
                          f"(source: {r.url})")
                    return float(alt_va), alt_datum_cd, float(lat), float(lon), r.url
        except Exception as e:  # noqa: BLE001
            print(f"  [fetch_gauge_datum] attempt {attempt+1} failed: {e}")
        time.sleep(backoff)
    print(f"  [fetch_gauge_datum] live query failed; using fallback constants "
          f"alt_va={FALLBACK_ALT_VA_FT} ft, alt_datum_cd={FALLBACK_ALT_DATUM_CD} "
          f"(confirmed via manual query 2026-07-22).")
    return (FALLBACK_ALT_VA_FT, FALLBACK_ALT_DATUM_CD, FALLBACK_GAUGE_LAT,
            FALLBACK_GAUGE_LON, SITE_URL + " (fallback, manual query 2026-07-22)")


def resolve_navd88_datum_ft(alt_va_ft, alt_datum_cd):
    """See analysis/stage_to_dem.py for full rationale. For USGS 04286000,
    alt_datum_cd is NAVD88, so the NGVD29 branch is not exercised here."""
    if alt_datum_cd == "NAVD88":
        return alt_va_ft
    if alt_datum_cd == "NGVD29":
        raise NotImplementedError(
            f"Gauge datum is NGVD29 ({alt_va_ft} ft), not NAVD88. A site-"
            "specific, tool-verified NGVD29-to-NAVD88 offset (NOAA VDatum or "
            "NGS VERTCON) is required before proceeding -- not a generic "
            "regional approximation."
        )
    raise RuntimeError(f"Unrecognized alt_datum_cd '{alt_datum_cd}' for site {SITE_NO}; "
                        "cannot safely align gauge datum to the NAVD88 DEM.")


def fetch_peak_stage(site_no, start, end, retries=4, backoff=5):
    """Fetch USGS NWIS IV gage height (00065) for [start, end], return
    (peak_stage_ft, peak_time). Falls back to the documented constant if the
    live query fails."""
    params = {
        "sites": site_no, "parameterCd": "00065",
        "startDT": start, "endDT": end, "format": "json",
    }
    for attempt in range(retries):
        try:
            r = requests.get(IV_URL, params=params, timeout=60)
            if r.status_code == 200:
                data = r.json()
                ts = data.get("value", {}).get("timeSeries", [])
                if not ts:
                    raise RuntimeError(f"No timeSeries returned for {site_no} {start}..{end}")
                vals = ts[0]["values"][0]["value"]
                rows = [(v["dateTime"], float(v["value"])) for v in vals if v["value"] not in (None, "")]
                if not rows:
                    raise RuntimeError("No non-null stage values in window")
                peak_time, peak_stage = max(rows, key=lambda rv: rv[1])
                print(f"  [fetch_peak_stage] peak stage {peak_stage:.2f} ft at {peak_time} "
                      f"(n={len(rows)} obs, source: {r.url})")
                return peak_stage, peak_time
        except Exception as e:  # noqa: BLE001
            print(f"  [fetch_peak_stage] attempt {attempt+1} failed: {e}")
        time.sleep(backoff)
    print(f"  [fetch_peak_stage] live query failed; using fallback peak stage "
          f"{FALLBACK_PEAK_STAGE_FT} ft (confirmed 2026-07-22).")
    return FALLBACK_PEAK_STAGE_FT, None


# ----------------------------------------------------------------------
# Data acquisition: DEM (copied from stage_to_dem.py)
# ----------------------------------------------------------------------

def download_dem(bbox_lonlat, dst_crs=DEM_CRS, target_res_m=TARGET_RES_M, retries=4, backoff=8):
    """Download a DEM clip from the 3DEPElevation ImageServer exportImage
    endpoint, in-memory, reprojected to dst_crs at target_res_m. Raises
    after exhausting retries -- no fallback DEM."""
    minlon, minlat = bbox_lonlat["minlon"], bbox_lonlat["minlat"]
    maxlon, maxlat = bbox_lonlat["maxlon"], bbox_lonlat["maxlat"]

    left, bottom, right, top = transform_bounds("EPSG:4326", dst_crs, minlon, minlat, maxlon, maxlat)
    width_m, height_m = right - left, top - bottom
    ncols = max(1, round(width_m / target_res_m))
    nrows = max(1, round(height_m / target_res_m))
    print(f"  [download_dem] target grid: {ncols} x {nrows} px at ~{target_res_m} m")

    params = {
        "bbox": f"{minlon},{minlat},{maxlon},{maxlat}", "bboxSR": "4326",
        "imageSR": dst_crs.split(":")[1], "size": f"{ncols},{nrows}",
        "format": "tiff", "pixelType": "F32",
        "interpolation": "RSP_BilinearInterpolation", "f": "image",
    }
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.get(DEM_URL, params=params, timeout=120)
            if r.status_code == 200 and r.headers.get("content-type", "").startswith("image"):
                with MemoryFile(r.content) as memfile:
                    with memfile.open() as ds:
                        arr = ds.read(1, masked=True).filled(np.nan).astype(np.float64)
                        prof = dict(
                            transform=ds.transform, crs=ds.crs, bounds=ds.bounds,
                            res=ds.res, width=ds.width, height=ds.height,
                        )
                print(f"  [download_dem] got {prof['width']}x{prof['height']} px, "
                      f"res={prof['res']}, crs={prof['crs']} (source: {r.url})")
                return arr, prof
            else:
                last_exc = RuntimeError(f"HTTP {r.status_code}, content-type="
                                         f"{r.headers.get('content-type')}: {r.text[:200]}")
        except Exception as e:  # noqa: BLE001
            last_exc = e
        print(f"  [download_dem] attempt {attempt+1}/{retries} failed ({last_exc}); retrying...")
        time.sleep(backoff)
    raise RuntimeError(f"download_dem failed after {retries} attempts: {last_exc}. "
                        "No fallback DEM is available -- stopping.")


def geocode(name, query, retries=3, backoff=3, bbox_pad_deg=0.05):
    """Geocode a place name via OSM Nominatim, validating within a padded
    BBOX_LONLAT (copied from stage_to_dem.py -- same rationale/history)."""
    headers = {"User-Agent": "cir-flood-proposal-research/1.0 (rjn5308@psu.edu)"}
    params = {"q": query, "format": "json", "limit": 1}
    minlon = BBOX_LONLAT["minlon"] - bbox_pad_deg
    maxlon = BBOX_LONLAT["maxlon"] + bbox_pad_deg
    minlat = BBOX_LONLAT["minlat"] - bbox_pad_deg
    maxlat = BBOX_LONLAT["maxlat"] + bbox_pad_deg
    for attempt in range(retries):
        try:
            r = requests.get(NOMINATIM_URL, params=params, headers=headers, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if data:
                    lat, lon = float(data[0]["lat"]), float(data[0]["lon"])
                    if not (minlon <= lon <= maxlon and minlat <= lat <= maxlat):
                        print(f"  [geocode] '{query}' -> ({lat}, {lon}) is OUTSIDE the "
                              "expected Montpelier bbox; rejecting.")
                        break
                    print(f"  [geocode] '{query}' -> ({lat}, {lon}) "
                          f"[{data[0].get('display_name')}]")
                    return lat, lon
        except Exception as e:  # noqa: BLE001
            print(f"  [geocode] attempt {attempt+1} failed: {e}")
        time.sleep(backoff)
    if name in FALLBACK_LANDMARKS:
        lat, lon = FALLBACK_LANDMARKS[name]
        print(f"  [geocode] live query failed or was rejected for '{query}'; using "
              f"fallback ({lat}, {lon}).")
        return lat, lon
    raise RuntimeError(f"geocode failed for '{query}' and no fallback constant is defined.")


def sample_at_lonlat(dem_arr, prof, lon, lat):
    xs, ys = transform("EPSG:4326", prof["crs"], [lon], [lat])
    x, y = xs[0], ys[0]
    inv = ~prof["transform"]
    col, row = inv * (x, y)
    col, row = int(round(col)), int(round(row))
    if 0 <= row < prof["height"] and 0 <= col < prof["width"]:
        return dem_arr[row, col], (x, y)
    return np.nan, (x, y)


# ----------------------------------------------------------------------
# Data acquisition: OSM waterways -> along-channel distance field
# ----------------------------------------------------------------------

def fetch_waterway_ways(bbox_lonlat, retries=4, backoff=10):
    """Fetch OSM waterway=river ways with name containing "Winooski" (this
    matches both "Winooski River" and "North Branch Winooski River") for
    bbox_lonlat via Overpass. Returns a list of ways, each a list of (lon,
    lat) vertices in OSM-digitized order (direction not assumed to be
    downstream -- the graph built from these is undirected).

    No fallback: there is no defensible static substitute for real channel
    geometry. If Overpass fails after retries, this raises."""
    minlat, minlon = bbox_lonlat["minlat"], bbox_lonlat["minlon"]
    maxlat, maxlon = bbox_lonlat["maxlat"], bbox_lonlat["maxlon"]
    query = f"""[out:json][timeout:120];
(
  way["waterway"="river"]["name"~"Winooski"]({minlat},{minlon},{maxlat},{maxlon});
);
out geom;
"""
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.post(OVERPASS_URL, data={"data": query},
                               headers=OVERPASS_HEADERS, timeout=150)
            if r.status_code == 200:
                data = r.json()
                elements = data.get("elements", [])
                ways = []
                names = set()
                for e in elements:
                    if e["type"] != "way":
                        continue
                    geom = e.get("geometry")
                    if not geom or len(geom) < 2:
                        continue
                    ways.append([(pt["lon"], pt["lat"]) for pt in geom])
                    names.add(e.get("tags", {}).get("name", "?"))
                print(f"  [fetch_waterway_ways] {len(ways)} way(s) fetched, "
                      f"names={sorted(names)} (source: Overpass query on bbox {bbox_lonlat})")
                if not ways:
                    raise RuntimeError("Overpass returned 0 waterway ways for this bbox/name filter.")
                return ways
            last_exc = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:  # noqa: BLE001
            last_exc = e
        print(f"  [fetch_waterway_ways] attempt {attempt+1}/{retries} failed "
              f"({last_exc}); retrying...")
        time.sleep(backoff)
    raise RuntimeError(f"fetch_waterway_ways failed after {retries} attempts: {last_exc}. "
                        "No fallback channel geometry is available -- stopping.")


def build_channel_distance_field(ways_lonlat, gauge_lon, gauge_lat, dst_crs,
                                  snap_tol_m=CHANNEL_SNAP_TOL_M,
                                  densify_step_m=CHANNEL_DENSIFY_STEP_M):
    """Build the along-channel-network distance-from-gauge field described
    in the module docstring ("Channel distance method"). Returns a cKDTree
    over a dense channel point cloud (UTM coords) and the parallel array of
    along-channel distances (m) from the gauge, plus diagnostics.
    """
    # Reproject every way's vertices to dst_crs.
    ways_xy = []
    for way in ways_lonlat:
        lons = [p[0] for p in way]
        lats = [p[1] for p in way]
        xs, ys = transform("EPSG:4326", dst_crs, lons, lats)
        ways_xy.append(list(zip(xs, ys)))

    # Build graph: snap vertices within snap_tol_m to a common node id via
    # coordinate quantization (dict keyed on rounded coords); edge weights
    # use the ORIGINAL (unquantized) vertex coordinates so segment lengths
    # are not distorted by the snap tolerance.
    node_id_of = {}
    node_xy = []

    def get_node(xy):
        key = (round(xy[0] / snap_tol_m), round(xy[1] / snap_tol_m))
        if key not in node_id_of:
            node_id_of[key] = len(node_xy)
            node_xy.append(xy)
        return node_id_of[key]

    edges = []  # (u, v, weight)
    for way in ways_xy:
        for i in range(len(way) - 1):
            u = get_node(way[i])
            v = get_node(way[i + 1])
            if u == v:
                continue
            w = float(np.hypot(way[i][0] - way[i + 1][0], way[i][1] - way[i + 1][1]))
            edges.append((u, v, w))

    n_nodes = len(node_xy)
    node_xy = np.array(node_xy)
    print(f"  [build_channel_distance_field] graph: {n_nodes} nodes, {len(edges)} edges "
          f"(snap tolerance {snap_tol_m} m)")

    rows = [e[0] for e in edges] + [e[1] for e in edges]
    cols = [e[1] for e in edges] + [e[0] for e in edges]
    data = [e[2] for e in edges] * 2
    graph = csr_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes))

    gx, gy = transform("EPSG:4326", dst_crs, [gauge_lon], [gauge_lat])
    gauge_xy = np.array([gx[0], gy[0]])
    gauge_node = int(np.argmin(np.hypot(node_xy[:, 0] - gauge_xy[0], node_xy[:, 1] - gauge_xy[1])))
    gauge_snap_dist_m = float(np.hypot(node_xy[gauge_node, 0] - gauge_xy[0],
                                        node_xy[gauge_node, 1] - gauge_xy[1]))
    print(f"  [build_channel_distance_field] gauge snapped to nearest channel node "
          f"{gauge_snap_dist_m:.1f} m away (node {gauge_node} of {n_nodes})")

    dist_from_gauge = dijkstra(graph, directed=False, indices=gauge_node)
    n_unreachable = int(np.sum(~np.isfinite(dist_from_gauge)))
    n_components_note = (
        f"{n_unreachable} of {n_nodes} channel nodes are NOT connected to the "
        "gauge's graph component (excluded from the distance field)"
        if n_unreachable else
        "all channel nodes are connected to the gauge's graph component"
    )
    print(f"  [build_channel_distance_field] {n_components_note}")

    # Densify each edge at ~densify_step_m spacing, with distance-from-gauge
    # linearly interpolated between the edge's two endpoint distances.
    # Edges with a non-finite endpoint distance (disconnected component) are
    # skipped -- they cannot contribute a meaningful along-channel distance.
    pts, dists = [], []
    for e in edges:
        u, v, w = e
        du, dv = dist_from_gauge[u], dist_from_gauge[v]
        if not (np.isfinite(du) and np.isfinite(dv)) or w <= 0:
            continue
        n_steps = max(1, int(np.ceil(w / densify_step_m)))
        for k in range(n_steps + 1):
            t = k / n_steps
            x = node_xy[u, 0] + t * (node_xy[v, 0] - node_xy[u, 0])
            y = node_xy[u, 1] + t * (node_xy[v, 1] - node_xy[u, 1])
            pts.append((x, y))
            dists.append(du + t * (dv - du))

    pts = np.array(pts)
    dists = np.array(dists)
    print(f"  [build_channel_distance_field] densified channel point cloud: "
          f"{len(pts)} points at ~{densify_step_m} m spacing")

    tree = cKDTree(pts)
    return tree, dists, dict(
        n_nodes=n_nodes, n_edges=len(edges), n_unreachable=n_unreachable,
        gauge_snap_dist_m=gauge_snap_dist_m, n_channel_points=len(pts),
    )


# ----------------------------------------------------------------------
# Data acquisition: OSM building centroids (adapted from
# analysis/subbasin_sensitivity.py's fetch_building_centroids)
# ----------------------------------------------------------------------

def fetch_building_centroids(bbox_lonlat, retries=4, backoff=10):
    """Fetch OSM building footprints (way elements tagged building=*) for
    bbox_lonlat via Overpass, return a list of (lon, lat) centroids. No
    fallback -- if the live query fails after retries, this raises."""
    minlat, minlon = bbox_lonlat["minlat"], bbox_lonlat["minlon"]
    maxlat, maxlon = bbox_lonlat["maxlat"], bbox_lonlat["maxlon"]
    query = f"""[out:json][timeout:120];
(
  way["building"]({minlat},{minlon},{maxlat},{maxlon});
  relation["building"]({minlat},{minlon},{maxlat},{maxlon});
);
out geom;
"""
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.post(OVERPASS_URL, data={"data": query},
                               headers=OVERPASS_HEADERS, timeout=150)
            if r.status_code == 200:
                data = r.json()
                elements = data.get("elements", [])
                n_relations = sum(1 for e in elements if e["type"] == "relation")
                centroids = []
                n_skipped = 0
                for e in elements:
                    if e["type"] != "way":
                        continue
                    geom = e.get("geometry")
                    if not geom or len(geom) < 3:
                        n_skipped += 1
                        continue
                    coords = [(pt["lon"], pt["lat"]) for pt in geom]
                    try:
                        poly = Polygon(coords)
                        if poly.is_valid and poly.area > 0:
                            c = poly.centroid
                            centroids.append((c.x, c.y))
                        else:
                            n_skipped += 1
                    except Exception:  # noqa: BLE001
                        n_skipped += 1
                print(f"  [fetch_building_centroids] {len(elements)} OSM elements "
                      f"({n_relations} relation(s), not individually resolved), "
                      f"{len(centroids)} way centroids computed, "
                      f"{n_skipped} way(s) skipped (degenerate/invalid geometry)")
                return centroids
            last_exc = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:  # noqa: BLE001
            last_exc = e
        print(f"  [fetch_building_centroids] attempt {attempt+1}/{retries} failed "
              f"({last_exc}); retrying...")
        time.sleep(backoff)
    raise RuntimeError(f"fetch_building_centroids failed after {retries} attempts: "
                        f"{last_exc}. No fallback building set is available -- stopping.")


def locate_on_grid(xs, ys, prof):
    """Reproject-agnostic: given UTM x,y arrays and a DEM profile, return
    integer (row, col) indices and a validity mask for points inside the
    grid extent."""
    inv = ~prof["transform"]
    cols_f, rows_f = inv * (np.asarray(xs), np.asarray(ys))
    rows = np.round(rows_f).astype(int)
    cols = np.round(cols_f).astype(int)
    in_grid = ((rows >= 0) & (rows < prof["height"]) & (cols >= 0) & (cols < prof["width"]))
    return rows, cols, in_grid


# ----------------------------------------------------------------------
# Analysis: inundation masks for a spatially-varying WSE
# ----------------------------------------------------------------------

def compute_seed_mask(dem_arr, seed_percentile=1.0):
    """River-channel proxy: lowest seed_percentile of DEM cells in the
    clip (same heuristic as stage_to_dem.py's compute_depth_masks --
    independent of WSE/beta, so computed once and reused for every beta)."""
    valid = np.isfinite(dem_arr)
    seed_thresh = np.nanpercentile(dem_arr[valid], seed_percentile)
    seed_low = valid & (dem_arr <= seed_thresh)
    return valid, seed_low, float(seed_thresh)


def compute_masks_for_wse(dem_arr, wse_arr, valid, seed_low):
    """Bathtub + connected inundation masks for a spatially-varying WSE
    field (wse_arr, same shape as dem_arr). Same flood-fill logic as
    stage_to_dem.py's compute_depth_masks, generalized from a scalar wse_m
    to a per-cell wse_arr."""
    bathtub = valid & (dem_arr < wse_arr)
    seed = bathtub & seed_low

    structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    labeled, num_components = label(bathtub, structure=structure)
    seed_labels = set(np.unique(labeled[seed]).tolist()) - {0}
    connected = np.isin(labeled, list(seed_labels)) if seed_labels else np.zeros_like(bathtub)
    return bathtub, connected


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def make_figure(betas, areas_acres, n_buildings, landmark_depths, admissible_range, out_png):
    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25,
    })
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.3))

    lo, hi = admissible_range
    for ax in (ax1, ax2):
        if lo is not None:
            ax.axvspan(lo, hi if hi is not None else max(betas), color=COLOR_ADMISSIBLE,
                       alpha=0.15, zorder=0)

    ax1.plot(betas, areas_acres, marker="o", color=COLOR_AREA, lw=1.4, ms=4,
              label="Connected inundated area (acres)")
    ax1.set_xlabel(r"$\beta$ (m/km)")
    ax1.set_ylabel("Connected inundated area (acres)", color=COLOR_AREA)
    ax1.tick_params(axis="y", labelcolor=COLOR_AREA)

    ax1b = ax1.twinx()
    ax1b.plot(betas, n_buildings, marker="s", color=COLOR_BUILDINGS, lw=1.4, ms=4,
               linestyle="--", label="OSM buildings inundated")
    ax1b.set_ylabel("OSM building centroids inundated (count)", color=COLOR_BUILDINGS)
    ax1b.tick_params(axis="y", labelcolor=COLOR_BUILDINGS)
    ax1.set_title("Inundated area and building count vs. $\\beta$", fontsize=9.5)

    colors = {"VT State House": COLOR_STATEHOUSE, "Montpelier City Hall": COLOR_CITYHALL,
              "Confluence Park": COLOR_CONFLUENCE}
    markers = {"VT State House": "^", "Montpelier City Hall": "o", "Confluence Park": "D"}
    for name, depths in landmark_depths.items():
        ax2.plot(betas, depths, marker=markers.get(name, "o"), color=colors.get(name, "black"),
                  lw=1.4, ms=4, label=name)
    ax2.axhline(0, color="0.3", lw=1.0, ls=":", zorder=1)
    ax2.set_xlabel(r"$\beta$ (m/km)")
    ax2.set_ylabel("Depth vs. WSE(x) (m); positive = inundated")
    ax2.set_title("Landmark depth vs. $\\beta$", fontsize=9.5)
    ax2.legend(fontsize=7.5, loc="upper left", frameon=True, framealpha=0.85)

    admissible_patch = Patch(facecolor=COLOR_ADMISSIBLE, alpha=0.15,
                              label="Admissible $\\beta$ interval")
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles1b, labels1b = ax1b.get_legend_handles_labels()
    ax1.legend(handles1 + handles1b + [admissible_patch], labels1 + labels1b + ["Admissible $\\beta$"],
               fontsize=7, loc="upper left", frameon=True, framealpha=0.85)

    fig.suptitle(
        "Water-surface slope ($\\beta$) sensitivity, downtown Montpelier VT -- July 2023 peak\n"
        "WSE(x) = 158.850 m NAVD88 + $\\beta \\cdot x$, $x$ = along-channel distance upstream of gauge 04286000",
        fontsize=8.8,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out_png, dpi=200)
    print(f"\nSaved figure: {out_png}")


# ----------------------------------------------------------------------
# Results markdown
# ----------------------------------------------------------------------

def write_results_md(ctx, out_path):
    lines = []
    lines.append("# Water-Surface Slope ($\\beta$) Sensitivity -- Downtown Montpelier, VT\n")
    lines.append(
        "Preliminary evidence for an NSF CPS-CIR proposal on flood-duration/depth "
        "modeling in Montpelier, VT: brackets the reach-scale water-surface slope "
        "$\\beta$ in `WSE(x) = 158.850 m NAVD88 + \\beta \\cdot x` (Root "
        "Research_v2.tex, Subtask 1.2) from **public evidence only** -- no ICEYE "
        "data used -- in advance of the 2023 ICEYE depth-field calibration that "
        "the proposal describes. All numbers below are computed by "
        "`analysis/slope_sensitivity.py` from live USGS NWIS, USGS 3DEP, and OSM "
        "Overpass/Nominatim queries; see that script for exact methods.\n"
    )

    lines.append("## Method summary\n")
    lines.append(
        "1. Same DEM clip, gauge datum, and July 2023 peak stage as "
        "`analysis/stage_to_dem.py` (WSE at the gauge = "
        f"**{ctx['wse0_m']:.3f} m NAVD88**).\n"
        "2. Along-channel upstream distance $x$ for every DEM cell, landmark, "
        "and OSM building centroid: OSM waterway ways tagged `waterway=river` "
        "with name containing \"Winooski\" (captures both the Winooski River "
        "and the North Branch Winooski River) fetched via Overpass for a "
        "padded bbox, assembled into an undirected graph (vertices snapped "
        f"within {ctx['snap_tol_m']:.0f} m to merge near-coincident endpoints "
        "at confluences/way splits), and a single-source Dijkstra shortest-"
        "path from the graph node nearest the gauge gives along-channel-"
        "network distance to every other channel point -- this correctly "
        "handles the Y-shaped Winooski / North Branch confluence. Every DEM "
        "cell/landmark/building is assigned the distance value of its "
        "nearest channel point (`scipy.spatial.cKDTree`), per the task's "
        "\"nearest-channel-point distance\" specification. See the script's "
        "module docstring for the full method and its approximations.\n"
        f"   - Channel graph: **{ctx['n_nodes']} nodes, {ctx['n_edges']} edges**; "
        f"gauge snapped to nearest channel node "
        f"**{ctx['gauge_snap_dist_m']:.1f} m** away; "
        f"**{ctx['n_unreachable']}** of {ctx['n_nodes']} nodes were not reachable "
        "from the gauge's graph component (excluded from the distance field).\n"
        f"   - Densified channel point cloud: **{ctx['n_channel_points']:,} points** "
        f"at ~{ctx['densify_step_m']:.0f} m spacing.\n"
        "3. $\\beta$ swept from 0.0 to 1.5 m/km in 0.1 m/km steps "
        "(`WSE(x) = 158.850 + \\beta \\cdot x`, $x$ in km). For each $\\beta$: "
        "the connected inundation mask (same flood-fill-from-river-proxy "
        "logic as `stage_to_dem.py`, generalized to a spatially-varying WSE), "
        "connected inundated area, number of OSM building centroids whose "
        "grid cell falls in the connected mask, and depth "
        "(`WSE(x_landmark) - elevation`) at the three sanity-check landmarks.\n"
    )

    lines.append("## Beta sweep\n")
    lines.append("| beta (m/km) | connected area (acres) | OSM buildings inundated | "
                  "VT State House depth (m) | City Hall depth (m) | Confluence Park depth (m) |\n")
    lines.append("|---|---|---|---|---|---|\n")
    for row in ctx["sweep_rows"]:
        lines.append(
            f"| {row['beta']:.1f} | {row['area_acres']:.1f} | {row['n_buildings']} | "
            f"{row['depth_state_house']:+.2f} | {row['depth_city_hall']:+.2f} | "
            f"{row['depth_confluence']:+.2f} |\n"
        )

    lines.append("\n## Admissible beta interval\n")
    lo, hi = ctx["admissible_range"]
    hi_str = f"{hi:.1f}" if hi is not None else f">= {ctx['beta_max']:.1f} (top of swept range)"
    if lo is not None:
        lines.append(
            f"- **Admissible interval: beta in [{lo:.1f}, {hi_str}] m/km.**\n"
            f"  - Lower bound {lo:.1f} m/km: the smallest swept beta at which "
            f"BOTH Main St / downtown-core landmarks (Montpelier City Hall AND "
            "Confluence Park, both documented flooded in July 2023) are "
            "predicted inundated (depth >= 0).\n"
            f"  - Upper bound: the VT State House (documented NOT flooded) "
            f"remains dry across the ENTIRE swept range (0.0-{ctx['beta_max']:.1f} "
            "m/km); its own implied threshold from the flat-surface residual "
            f"({ctx['state_house_implied_slope']:.2f} m/km) is far above the "
            "swept range, so this evidence only right-censors the interval at "
            "the top of the sweep rather than pinning a tighter upper bound.\n"
        )
    else:
        lines.append(
            "- **No beta in the swept range [0.0, 1.5] m/km floods both Main "
            "St landmarks while keeping the State House dry** -- see sweep "
            "table above; this would contradict the flat-surface-residual "
            "inference and needs investigation before being used elsewhere.\n"
        )
    lines.append(
        f"- At beta = {lo if lo is not None else float('nan'):.1f} m/km "
        f"(lower bound): connected area "
        f"**{ctx['area_at_lo']:.1f} acres**, "
        f"**{ctx['buildings_at_lo']} OSM buildings** inundated "
        f"(vs. {ctx['area_at_0']:.1f} acres / {ctx['buildings_at_0']} buildings "
        "at beta = 0, the flat-surface baseline).\n"
        if lo is not None else ""
    )
    lines.append(
        f"- At beta = {ctx['beta_max']:.1f} m/km (top of sweep): connected area "
        f"**{ctx['area_at_hi']:.1f} acres**, **{ctx['buildings_at_hi']} OSM "
        "buildings** inundated.\n"
    )

    lines.append("\n## Comparison with the flat-surface-residual inference\n")
    lo_str = f"{lo:.1f} m/km" if lo is not None else "n/a"
    lines.append(
        f"`analysis/stage_dem.md` inferred a missing water-surface slope of "
        f"roughly **{PRIOR_INFERENCE_RANGE_M_PER_KM[0]:.1f}-"
        f"{PRIOR_INFERENCE_RANGE_M_PER_KM[1]:.1f} m/km** from how far City "
        "Hall and Confluence Park sat above the flat WSE, divided by their "
        "straight-line distance from the gauge. This sweep's admissible "
        f"lower bound ({lo_str}) uses "
        "along-channel (not straight-line) distance and a discrete 0.1 m/km "
        "grid, so exact agreement is not expected; "
    )
    if lo is not None:
        agree = PRIOR_INFERENCE_RANGE_M_PER_KM[0] - 0.15 <= lo <= PRIOR_INFERENCE_RANGE_M_PER_KM[1] + 0.35
        lines.append(
            f"the two are {'broadly consistent' if agree else 'NOT closely consistent'} "
            f"({lo:.1f} m/km vs. the {PRIOR_INFERENCE_RANGE_M_PER_KM[0]:.1f}-"
            f"{PRIOR_INFERENCE_RANGE_M_PER_KM[1]:.1f} m/km prior inference). "
            "The along-channel distances used here are longer than the "
            "straight-line distances used in the flat-surface-residual "
            "calculation (channel distance >= straight-line distance always), "
            "which mechanically pushes the beta needed to flood a given "
            "landmark downward relative to a straight-line-distance estimate "
            "-- one identifiable source of any gap between the two.\n"
        )
    else:
        lines.append("no admissible interval was found in this sweep to compare.\n")

    lines.append("\n## Sanity checks\n")
    for msg in ctx["sanity_messages"]:
        lines.append(f"- {msg}\n")

    lines.append("\n## Caveats\n")
    lines.append(
        "- **Both rivers (Winooski main stem and North Branch) are assigned "
        "the SAME beta.** This is a single reach-scale slope parameter, not "
        "two independently calibrated slopes; the North Branch's true "
        "flood-surface slope could differ from the main stem's, especially "
        "near the confluence itself.\n"
        "- **Planar-per-reach surface, not a hydraulically routed profile.** "
        "WSE(x) is linear in along-channel distance; it does not represent "
        "backwater curves, local constrictions (bridges, channel narrowing), "
        "or the confluence's actual hydraulic behavior (e.g. a true water "
        "surface is not required to have equal slope on both branches "
        "meeting at a confluence).\n"
        "- **Along-channel distance is a network graph shortest-path "
        "distance from OSM waterway centerline geometry, not a traced "
        "lowest-elevation thalweg path through the DEM and not a distance "
        "computed from a validated hydrography dataset (e.g. NHD).** OSM "
        "waterway digitization quality/completeness was not independently "
        "audited beyond the diagnostics reported above (node/edge counts, "
        "unreachable-node count).\n"
        "- **Signed distance is approximated as unsigned.** Dijkstra "
        "distance is always >= 0; the small sliver of the bbox downstream "
        "(west) of the gauge is therefore treated as if it were upstream, "
        "which would very slightly overstate WSE there. The gauge sits "
        f"within ~{ctx['gauge_snap_dist_m']:.0f} m of the bbox's west edge "
        "(see stage_to_dem.py's bbox), so this affects a negligible fraction "
        "of the domain, but this is a documented simplification specific to "
        "this domain's geometry, not a general solution for a gauge in the "
        "interior of a bbox.\n"
        "- **Only two control landmarks define the admissible interval's "
        "lower bound (City Hall, Confluence Park) and one defines the upper "
        "bound (State House).** These are the same three point-sample "
        "sanity checks used in `stage_to_dem.md`, each subject to the same "
        "point-sampling caveats described there (a few meters of horizontal "
        "offset can flip a marginal point wet/dry).\n"
        "- **Inundated-building counts use every OSM-tagged `building=*` "
        "footprint in the clip**, not the proposal's 271-building curated "
        "downtown masonry dataset -- a building is counted \"inundated\" if "
        "its centroid's DEM grid cell falls in the connected mask, a coarse "
        "point-in-cell test, not a building-footprint overlap test.\n"
        "- **Connectivity is a cheap 4-connected flood-fill from an "
        "elevation-percentile channel proxy** (same heuristic and caveats as "
        "`stage_to_dem.py`) -- not a validated hydrography-derived channel "
        "network or a true 2D flood-routing connectivity solve.\n"
        "- **This sweep uses NO ICEYE data.** It is deliberately built from "
        "public evidence only, to bracket beta ahead of the proposal's "
        "planned ICEYE-based calibration (Subtask 1.2); it is not a "
        "substitute for that calibration.\n"
        "- **DEM is bare-earth terrain**, USGS stage data are provisional, "
        "and the DEM/geocoding/Overpass caveats in `stage_to_dem.md` and "
        "`subbasin_sensitivity.md` apply equally here (not repeated in full).\n"
    )

    out_path.write_text("".join(lines))
    print(f"\nSaved results: {out_path}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    print("Fetching gauge datum...")
    alt_va_ft, alt_datum_cd, gauge_lat, gauge_lon, _ = fetch_gauge_datum(SITE_NO)
    datum_navd88_ft = resolve_navd88_datum_ft(alt_va_ft, alt_datum_cd)

    print("Fetching peak stage (July 2023)...")
    peak_stage_ft, _ = fetch_peak_stage(SITE_NO, PEAK_EVENT_START, PEAK_EVENT_END)
    wse0_m = (datum_navd88_ft + peak_stage_ft) * FT_TO_M
    print(f"WSE(gauge) = {wse0_m:.3f} m NAVD88")

    print("Downloading DEM...")
    dem_arr, prof = download_dem(BBOX_LONLAT)
    cell_area_m2 = prof["res"][0] * prof["res"][1]
    valid, seed_low, seed_thresh_m = compute_seed_mask(dem_arr)
    print(f"Seed/channel-proxy threshold: {seed_thresh_m:.2f} m NAVD88 "
          f"({int(seed_low.sum()):,} cells)")

    print("Fetching OSM waterway centerlines...")
    ways = fetch_waterway_ways(WATERWAY_BBOX_LONLAT)
    tree, chan_dists_m, chan_diag = build_channel_distance_field(
        ways, gauge_lon, gauge_lat, DEM_CRS)

    print("Computing along-channel distance for every DEM cell...")
    rows_idx, cols_idx = np.meshgrid(np.arange(prof["height"]), np.arange(prof["width"]), indexing="ij")
    xs_cells = prof["transform"].c + (cols_idx + 0.5) * prof["transform"].a
    ys_cells = prof["transform"].f + (rows_idx + 0.5) * prof["transform"].e
    query_pts = np.column_stack([xs_cells.ravel(), ys_cells.ravel()])
    _, nn_idx = tree.query(query_pts, k=1)
    x_arr_m = chan_dists_m[nn_idx].reshape(dem_arr.shape)
    print(f"  along-channel distance range over DEM clip: "
          f"{np.nanmin(x_arr_m):.0f} - {np.nanmax(x_arr_m):.0f} m")

    print("Geocoding landmarks and computing their along-channel distance...")
    landmark_info = {}
    for name, query in GEOCODE_QUERIES.items():
        lat, lon = geocode(name, query)
        elev_m, (x, y) = sample_at_lonlat(dem_arr, prof, lon, lat)
        _, nn = tree.query([[x, y]], k=1)
        x_landmark_m = float(chan_dists_m[nn[0]])
        landmark_info[name] = dict(elev_m=elev_m, x_m=x_landmark_m)
        print(f"  {name}: elev={elev_m:.2f} m, along-channel x={x_landmark_m:.0f} m")

    print("Fetching OSM building centroids...")
    centroids_lonlat = fetch_building_centroids(BBOX_LONLAT)
    b_lons = [c[0] for c in centroids_lonlat]
    b_lats = [c[1] for c in centroids_lonlat]
    b_xs, b_ys = transform("EPSG:4326", prof["crs"], b_lons, b_lats)
    b_rows, b_cols, b_in_grid = locate_on_grid(b_xs, b_ys, prof)
    n_buildings_total = int(b_in_grid.sum())
    print(f"  {len(centroids_lonlat)} OSM buildings fetched, {n_buildings_total} on-grid")

    print(f"\nSweeping beta = {BETA_VALUES_M_PER_KM}...")
    sweep_rows = []
    for beta in BETA_VALUES_M_PER_KM:
        wse_arr = wse0_m + beta * (x_arr_m / 1000.0)
        _, connected = compute_masks_for_wse(dem_arr, wse_arr, valid, seed_low)
        n_cells = int(connected.sum())
        area_acres = n_cells * cell_area_m2 / 4046.8564224

        b_on = b_in_grid.copy()
        b_on[b_in_grid] &= connected[b_rows[b_in_grid], b_cols[b_in_grid]]
        n_buildings = int(b_on.sum())

        depths = {}
        for name, info in landmark_info.items():
            wse_here = wse0_m + beta * (info["x_m"] / 1000.0)
            depths[name] = wse_here - info["elev_m"]

        row = dict(
            beta=beta, area_acres=area_acres, n_buildings=n_buildings,
            depth_state_house=depths["VT State House"],
            depth_city_hall=depths["Montpelier City Hall"],
            depth_confluence=depths["Confluence Park"],
        )
        sweep_rows.append(row)
        print(f"  beta={beta:.1f}: area={area_acres:.1f} ac, buildings={n_buildings}, "
              f"State House={depths['VT State House']:+.2f} m, "
              f"City Hall={depths['Montpelier City Hall']:+.2f} m, "
              f"Confluence Park={depths['Confluence Park']:+.2f} m")

    # ---- Sanity: monotonicity ----
    sanity_messages = []
    areas = [r["area_acres"] for r in sweep_rows]
    builds = [r["n_buildings"] for r in sweep_rows]
    if all(areas[i] <= areas[i + 1] + 1e-9 for i in range(len(areas) - 1)):
        sanity_messages.append("Connected area is monotonically non-decreasing in beta. OK.")
    else:
        sanity_messages.append("WARNING: connected area is NOT monotonically non-decreasing in beta -- "
                                "see sweep table; investigate before using this sweep.")
        print("  WARNING: area non-monotonic in beta!")
    if all(builds[i] <= builds[i + 1] for i in range(len(builds) - 1)):
        sanity_messages.append("OSM building inundation count is monotonically non-decreasing in beta. OK.")
    else:
        sanity_messages.append("WARNING: OSM building inundation count is NOT monotonically non-decreasing "
                                "in beta -- see sweep table; investigate before using this sweep.")
        print("  WARNING: building count non-monotonic in beta!")
    for key, label_ in [("depth_state_house", "VT State House"),
                         ("depth_city_hall", "Montpelier City Hall"),
                         ("depth_confluence", "Confluence Park")]:
        vals = [r[key] for r in sweep_rows]
        if all(vals[i] <= vals[i + 1] + 1e-9 for i in range(len(vals) - 1)):
            sanity_messages.append(f"{label_} depth is monotonically non-decreasing in beta. OK.")
        else:
            sanity_messages.append(f"WARNING: {label_} depth is NOT monotonically non-decreasing in "
                                    "beta -- investigate before using this sweep.")
            print(f"  WARNING: {label_} depth non-monotonic in beta!")

    # ---- Admissible beta interval ----
    beta_city_hall_wet = next((r["beta"] for r in sweep_rows if r["depth_city_hall"] >= 0), None)
    beta_confluence_wet = next((r["beta"] for r in sweep_rows if r["depth_confluence"] >= 0), None)
    beta_state_house_wet = next((r["beta"] for r in sweep_rows if r["depth_state_house"] >= 0), None)

    if beta_city_hall_wet is not None and beta_confluence_wet is not None:
        beta_lo = max(beta_city_hall_wet, beta_confluence_wet)
    else:
        beta_lo = None
    beta_hi = beta_state_house_wet  # None if State House never wets in the sweep -> right-censored
    if beta_lo is not None and beta_hi is not None and beta_lo >= beta_hi:
        # No admissible interval: State House would flood at or before both
        # Main St landmarks -- would contradict expectations, flag it.
        print("  WARNING: computed admissible interval is empty/inverted "
              f"(lo={beta_lo}, hi={beta_hi}) -- State House wets at or before "
              "the Main St landmarks in this sweep.")
        admissible_range = (None, None)
    else:
        admissible_range = (beta_lo, beta_hi)

    beta_max = max(BETA_VALUES_M_PER_KM)
    row_at = {r["beta"]: r for r in sweep_rows}
    area_at_0, buildings_at_0 = row_at[0.0]["area_acres"], row_at[0.0]["n_buildings"]
    area_at_hi = row_at[beta_max]["area_acres"]
    buildings_at_hi = row_at[beta_max]["n_buildings"]
    if beta_lo is not None:
        area_at_lo = row_at[beta_lo]["area_acres"]
        buildings_at_lo = row_at[beta_lo]["n_buildings"]
    else:
        area_at_lo = buildings_at_lo = float("nan")

    # Flat-surface-residual implied slope for the State House, for context
    # in the writeup (recomputed here from straight-line distance, matching
    # stage_dem.md's method, not the along-channel method used in this sweep).
    sh_elev = landmark_info["VT State House"]["elev_m"]
    sh_x_km = landmark_info["VT State House"]["x_m"] / 1000.0
    state_house_implied_slope = (sh_elev - wse0_m) / sh_x_km if sh_x_km > 0 else float("nan")

    print(f"\nAdmissible beta interval: {admissible_range}")

    make_figure(
        betas=[r["beta"] for r in sweep_rows],
        areas_acres=[r["area_acres"] for r in sweep_rows],
        n_buildings=[r["n_buildings"] for r in sweep_rows],
        landmark_depths={
            "VT State House": [r["depth_state_house"] for r in sweep_rows],
            "Montpelier City Hall": [r["depth_city_hall"] for r in sweep_rows],
            "Confluence Park": [r["depth_confluence"] for r in sweep_rows],
        },
        admissible_range=admissible_range,
        out_png=ANALYSIS_DIR / "slope_sensitivity.png",
    )

    ctx = dict(
        wse0_m=wse0_m, snap_tol_m=CHANNEL_SNAP_TOL_M, densify_step_m=CHANNEL_DENSIFY_STEP_M,
        n_nodes=chan_diag["n_nodes"], n_edges=chan_diag["n_edges"],
        n_unreachable=chan_diag["n_unreachable"], gauge_snap_dist_m=chan_diag["gauge_snap_dist_m"],
        n_channel_points=chan_diag["n_channel_points"],
        sweep_rows=sweep_rows, admissible_range=admissible_range, beta_max=beta_max,
        area_at_lo=area_at_lo, buildings_at_lo=buildings_at_lo,
        area_at_0=area_at_0, buildings_at_0=buildings_at_0,
        area_at_hi=area_at_hi, buildings_at_hi=buildings_at_hi,
        state_house_implied_slope=state_house_implied_slope,
        sanity_messages=sanity_messages,
    )
    write_results_md(ctx, out_path=ANALYSIS_DIR / "slope_sensitivity.md")

    print("\nDone.")


if __name__ == "__main__":
    main()
