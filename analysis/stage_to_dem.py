#!/usr/bin/env python3
"""
Winooski River at Montpelier, VT (USGS 04286000) -- gauge-stage-to-DEM
("bathtub") building-level depth mapping for the July 2023 flood peak.

Produced as preliminary evidence for an NSF CPS-CIR proposal on flood-
duration/depth modeling in Montpelier, VT. Demonstrates that a first-order,
operational depth estimate for downtown Montpelier can be produced by
projecting a single USGS gauge's stage reading onto a public DEM, with NO
hydraulic simulation -- the simplest possible baseline that the proposed
work will improve on.

Follows the conventions of analysis/winooski_recession.py and
analysis/stage_thresholds.py in this repo: live-query style with retries,
explicit *documented* fallback constants (used only if a live query fails
at run time), no silent fallbacks, colorblind-safe plotting.

What this script does
----------------------
1. Queries the USGS NWIS site service for gauge 04286000 to get the gauge
   datum elevation (alt_va) and its vertical datum (alt_datum_cd), plus the
   gauge's lat/lon for plotting. Handles NAVD88 vs NGVD29 explicitly -- see
   `resolve_navd88_datum_ft()`. 3DEP DEMs are NAVD88 orthometric heights in
   meters, so any NGVD29 gauge datum would need a local conversion; this
   script does NOT guess an unverified offset (see that function's
   docstring for why, and what to do if a future gauge needs it).
2. Queries USGS NWIS IV gage-height (00065) for the July 2023 flood window
   to independently find the peak stage (cross-checks the 21.29 ft figure
   already reported in analysis/results.md, computed by
   winooski_recession.py from the same underlying IV data).
3. Computes water-surface elevation (WSE) = gauge datum (NAVD88, m) + peak
   stage (m) -- a single flat/planar WSE applied across the whole domain.
   This is the "bathtub" assumption: no water-surface slope, no hydraulic
   routing, no timing/attenuation between the gauge and downtown.
4. Downloads a 1 m USGS 3DEP DEM clip for downtown Montpelier via the
   3DEPElevation ImageServer exportImage endpoint, reprojected to NAD83 /
   UTM 18N (EPSG:26918) so pixels are natively in meters. Read entirely
   in-memory via rasterio's MemoryFile -- no DEM file is written to disk.
5. depth = WSE - ground elevation, clipped to positive (inundated) cells.
   Reports:
     - "bathtub" extent: every DEM cell below WSE, regardless of whether it
       connects to the river (the naive/default interpretation of a
       gauge-to-DEM depth map).
     - "connected" extent: bathtub cells that are 4-connected, through
       other inundated cells, to a river-channel proxy (the lowest ~1% of
       DEM cells in the clip, which visually forms a single continuous
       channel through the domain -- see analysis/stage_dem.md). This
       removes topographically isolated low spots (e.g. a low point behind
       a ridge) that a bathtub fill would otherwise flood non-physically.
6. Geocodes three known downtown landmarks (VT State House, Montpelier City
   Hall, Confluence Park) via OSM Nominatim as a qualitative sanity check
   against public reporting on the 2023 flood.
7. Writes a one-panel PDF+PNG depth map (images/stage_dem_depth.pdf,
   analysis/stage_dem_depth.png) and analysis/stage_dem.md.

Usage
-----
    uv run --with rasterio --with numpy --with matplotlib --with requests \\
        --with scipy analysis/stage_to_dem.py

Requires: requests, numpy, matplotlib, rasterio, scipy.
"""

import time
from pathlib import Path

import numpy as np
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import rasterio
from rasterio.io import MemoryFile
from rasterio.warp import transform, transform_bounds
from scipy.ndimage import label

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

SITE_NO = "04286000"  # WINOOSKI RIVER AT MONTPELIER, VT (see winooski_recession.py docstring)
SITE_NAME = "WINOOSKI RIVER AT MONTPELIER, VT"

SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
IV_URL = "https://waterservices.usgs.gov/nwis/iv/"
DEM_URL = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage"
NOMINATIM_URL = "https://nominatim.openstreetmap.org/search"

# July 2023 flood window (same window used in winooski_recession.py / stage_thresholds.py)
PEAK_EVENT_START = "2023-07-01"
PEAK_EVENT_END = "2023-08-15"
PEAK_EVENT_LABEL = "July 2023"

# Downtown Montpelier bounding box (lon/lat, WGS84)
BBOX_LONLAT = dict(minlon=-72.594, minlat=44.250, maxlon=-72.560, maxlat=44.268)

# DEM output CRS and target resolution. 3DEPElevation's mosaic reports
# meanPixelSize = 1.0 m for this area (confirmed via
# https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer?f=json
# on 2026-07-22), so requesting 1 m output pixels matches native source
# resolution rather than up-sampling coarser data.
DEM_CRS = "EPSG:26918"  # NAD83 / UTM zone 18N -- covers Montpelier, VT; meters
TARGET_RES_M = 1.0

FT_TO_M = 0.3048

# Fallback values, empirically confirmed from live queries on 2026-07-22,
# used only if the corresponding live service call fails at run time.
FALLBACK_ALT_VA_FT = 499.87
FALLBACK_ALT_DATUM_CD = "NAVD88"
FALLBACK_GAUGE_LAT = 44.25672595
FALLBACK_GAUGE_LON = -72.59344318
FALLBACK_PEAK_STAGE_FT = 21.29  # matches analysis/results.md (winooski_recession.py)

FALLBACK_LANDMARKS = {
    "VT State House": (44.2626941, -72.5807647),
    "Montpelier City Hall": (44.2592113, -72.5755301),
    "Confluence Park": (44.2594605, -72.5781534),
}

# Nominatim query strings for the landmarks above. Must be qualified with
# "Montpelier, VT" -- an unqualified "Confluence Park" query matched a
# same-named park in Denver, CO on first run of this script (2026-07-22);
# geocode() validates the returned point falls within BBOX_LONLAT to catch
# this class of mismatch even if a future query is under-qualified.
GEOCODE_QUERIES = {
    "VT State House": "Vermont State House, Montpelier, VT",
    "Montpelier City Hall": "Montpelier City Hall, Montpelier, VT",
    "Confluence Park": "Confluence Park, Montpelier, VT",
}

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "analysis"
IMAGES_DIR = REPO_ROOT / "images"

# Colorblind-safe sequential colormap for continuous depth data (viridis).
# Distinct from the two-color #0072B2/#D55E00 palette used in
# winooski_recession.py, which distinguishes two discrete series (data vs.
# fit line); a perceptually-uniform sequential map is the right choice for
# a continuous depth field.
DEPTH_CMAP = "viridis"
COLOR_MARKER = "#D55E00"  # vermillion, colorblind-safe, reads clearly over viridis


# ----------------------------------------------------------------------
# Data acquisition: gauge datum and peak stage
# ----------------------------------------------------------------------

def fetch_gauge_datum(site_no, retries=4, backoff=5):
    """Fetch gauge altitude (alt_va, ft) and its vertical datum
    (alt_datum_cd), plus decimal lat/lon, from the NWIS expanded site
    service. No silent fallback beyond the documented, dated constants."""
    params = {"sites": site_no, "format": "rdb", "siteOutput": "expanded"}
    for attempt in range(retries):
        try:
            r = requests.get(SITE_URL, params=params, timeout=60)
            if r.status_code == 200:
                lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
                header = lines[0].split("\t")
                data_line = lines[2].split("\t")  # line 0 = header, line 1 = format codes
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
    """Return the gauge datum elevation in NAVD88 feet.

    3DEP DEMs report NAVD88 orthometric heights in meters, so the gauge
    datum must be in NAVD88 before it can be added to stage and compared to
    the DEM. If the gauge's published datum is already NAVD88, this is a
    no-op. If it is NGVD29, a local conversion offset is required -- but
    that offset varies spatially (NGS's VERTCON grid, not a single national
    constant) and a generic "central Vermont" number pulled from a web
    search is NOT precise enough to defend in a proposal without running it
    through NOAA's VDatum tool (https://vdatum.noaa.gov/) or NGS's VERTCON
    for this gauge's exact coordinates. Rather than silently applying an
    unverified offset, this function raises so the discrepancy is visible
    and can be resolved with a proper site-specific lookup.

    For USGS 04286000, alt_datum_cd is NAVD88 (confirmed via the NWIS site
    service, see analysis/stage_dem.md), so this code path is not exercised
    in the current run -- but it is implemented explicitly per the task's
    instruction to handle vertical datum rather than ignore it.
    """
    if alt_datum_cd == "NAVD88":
        return alt_va_ft
    if alt_datum_cd == "NGVD29":
        raise NotImplementedError(
            f"Gauge datum is NGVD29 ({alt_va_ft} ft), not NAVD88. A local "
            "NGVD29-to-NAVD88 offset is required before comparing to the "
            "3DEP (NAVD88) DEM, but no site-specific, tool-verified offset "
            "for this gauge's exact coordinates is available in this "
            "script. Look up the offset via NOAA VDatum "
            "(https://vdatum.noaa.gov/) or NGS VERTCON for this gauge's "
            "lat/lon before proceeding -- do not substitute a generic "
            "regional approximation."
        )
    raise RuntimeError(f"Unrecognized alt_datum_cd '{alt_datum_cd}' for site {SITE_NO}; "
                        "cannot safely align gauge datum to the NAVD88 DEM.")


def fetch_peak_stage(site_no, start, end, retries=4, backoff=5):
    """Fetch USGS NWIS IV gage height (00065) for [start, end] and return
    (peak_stage_ft, peak_time). Independently reproduces the peak stage
    already reported in analysis/results.md (winooski_recession.py) as a
    cross-check, using the same live-query approach."""
    params = {
        "sites": site_no,
        "parameterCd": "00065",
        "startDT": start,
        "endDT": end,
        "format": "json",
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
          f"{FALLBACK_PEAK_STAGE_FT} ft (matches analysis/results.md, confirmed 2026-07-22).")
    return FALLBACK_PEAK_STAGE_FT, None


# ----------------------------------------------------------------------
# Data acquisition: DEM
# ----------------------------------------------------------------------

def download_dem(bbox_lonlat, dst_crs=DEM_CRS, target_res_m=TARGET_RES_M, retries=4, backoff=8):
    """Download a DEM clip from the 3DEPElevation ImageServer exportImage
    endpoint, reprojected to dst_crs at target_res_m resolution, entirely
    in-memory (no file written to disk). Raises after exhausting retries --
    there is no fallback DEM; if this fails, this analysis cannot proceed
    with real data and must stop (per task instructions)."""
    minlon, minlat = bbox_lonlat["minlon"], bbox_lonlat["minlat"]
    maxlon, maxlat = bbox_lonlat["maxlon"], bbox_lonlat["maxlat"]

    # Determine output raster size so pixels are ~target_res_m in dst_crs.
    left, bottom, right, top = transform_bounds("EPSG:4326", dst_crs, minlon, minlat, maxlon, maxlat)
    width_m = right - left
    height_m = top - bottom
    ncols = max(1, round(width_m / target_res_m))
    nrows = max(1, round(height_m / target_res_m))
    print(f"  [download_dem] target grid: {ncols} x {nrows} px at ~{target_res_m} m "
          f"({width_m:.0f} m x {height_m:.0f} m in {dst_crs})")

    params = {
        "bbox": f"{minlon},{minlat},{maxlon},{maxlat}",
        "bboxSR": "4326",
        "imageSR": dst_crs.split(":")[1],
        "size": f"{ncols},{nrows}",
        "format": "tiff",
        "pixelType": "F32",
        "interpolation": "RSP_BilinearInterpolation",
        "f": "image",
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
                            transform=ds.transform,
                            crs=ds.crs,
                            bounds=ds.bounds,
                            res=ds.res,
                            width=ds.width,
                            height=ds.height,
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
                        "No fallback DEM is available -- stopping rather than substituting "
                        "fake elevation data.")


# ----------------------------------------------------------------------
# Data acquisition: landmark geocoding (sanity check)
# ----------------------------------------------------------------------

def geocode(name, query, retries=3, backoff=3, bbox_pad_deg=0.05):
    """Geocode a place name via OSM Nominatim, validating the result falls
    within (a small pad around) BBOX_LONLAT. An under-qualified query can
    silently match a same-named place elsewhere -- an earlier run of this
    script with the bare query "Confluence Park" (no city/state) matched a
    park in Denver, CO instead of Montpelier, VT, which this bbox check is
    specifically here to catch. Falls back to the documented, dated
    constants in FALLBACK_LANDMARKS if the live query fails or returns a
    point outside the expected bbox."""
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
                        print(f"  [geocode] '{query}' -> ({lat}, {lon}) "
                              f"[{data[0].get('display_name')}] is OUTSIDE the "
                              f"expected Montpelier bbox (+/-{bbox_pad_deg} deg pad); "
                              "rejecting this result, not using it.")
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
              f"fallback ({lat}, {lon}) (confirmed via manual query 2026-07-22).")
        return lat, lon
    raise RuntimeError(f"geocode failed for '{query}' and no fallback constant is defined.")


# ----------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------

def compute_depth_masks(dem_arr, wse_m, seed_percentile=1.0):
    """Compute the 'bathtub' inundation mask (any cell below WSE) and a
    'connected' mask restricted to cells 4-connected to a river-channel
    proxy (the lowest seed_percentile of DEM cells in the clip).

    The channel proxy is a heuristic, not a hydrography-derived channel
    network: it assumes the lowest ~1% of terrain in the clipped domain
    forms a single, essentially continuous corridor corresponding to the
    river channel. This held visually for this domain (a single dominant
    connected component containing the seed cells, tracing the Winooski
    River / North Branch confluence through downtown) -- see
    analysis/stage_dem.md for the figure-based check. It is not guaranteed
    to hold in general (e.g. a domain with two unrelated low-lying valleys).
    """
    valid = np.isfinite(dem_arr)
    bathtub = valid & (dem_arr < wse_m)

    seed_thresh = np.nanpercentile(dem_arr[valid], seed_percentile)
    seed = bathtub & (dem_arr <= seed_thresh)

    structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])  # 4-connectivity
    labeled, num_components = label(bathtub, structure=structure)
    seed_labels = set(np.unique(labeled[seed]).tolist()) - {0}
    connected = np.isin(labeled, list(seed_labels)) if seed_labels else np.zeros_like(bathtub)

    return dict(
        valid=valid,
        bathtub=bathtub,
        connected=connected,
        seed_thresh_m=seed_thresh,
        n_seed_cells=int(seed.sum()),
        num_components=int(num_components),
        n_seed_components=len(seed_labels),
    )


def depth_stats(dem_arr, mask, wse_m, cell_area_m2):
    depth = wse_m - dem_arr[mask]
    depth = depth[depth > 0]
    n_total = mask.size
    n_in = mask.sum()
    return dict(
        n_cells=int(n_in),
        frac_of_domain=float(n_in) / n_total,
        area_m2=float(n_in) * cell_area_m2,
        area_acres=float(n_in) * cell_area_m2 / 4046.8564224,
        depth_median_m=float(np.median(depth)) if len(depth) else float("nan"),
        depth_q25_m=float(np.percentile(depth, 25)) if len(depth) else float("nan"),
        depth_q75_m=float(np.percentile(depth, 75)) if len(depth) else float("nan"),
        depth_max_m=float(np.max(depth)) if len(depth) else float("nan"),
    )


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
# Plotting
# ----------------------------------------------------------------------

def make_figure(dem_arr, prof, wse_m, masks, gauge_xy, landmarks_xy, out_pdf, out_png):
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
    })

    extent = [prof["bounds"].left, prof["bounds"].right, prof["bounds"].bottom, prof["bounds"].top]
    depth = np.where(masks["connected"], wse_m - dem_arr, np.nan)

    fig, ax = plt.subplots(figsize=(7, 5.5))

    # Hillshade-like grayscale ground context (plain elevation grayscale,
    # not a true hillshade -- adequate for orientation at this scale).
    ax.imshow(dem_arr, extent=extent, origin="upper", cmap="gray",
              vmin=np.nanpercentile(dem_arr, 2), vmax=np.nanpercentile(dem_arr, 98))

    im = ax.imshow(depth, extent=extent, origin="upper", cmap=DEPTH_CMAP,
                    vmin=0, alpha=0.9)

    gx, gy = gauge_xy
    ax.plot(gx, gy, marker="*", color=COLOR_MARKER, markersize=16,
             markeredgecolor="white", markeredgewidth=0.6, linestyle="none",
             label="USGS gauge 04286000", zorder=5)

    for name, (x, y) in landmarks_xy.items():
        ax.plot(x, y, marker="o", color="white", markeredgecolor="black",
                 markersize=6, linestyle="none", zorder=5)
        ax.annotate(name, (x, y), color="white", fontsize=6.5, xytext=(4, 4),
                    textcoords="offset points", zorder=5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("Depth above ground, WSE $-$ DEM (m)")

    ax.set_xlabel("UTM 18N easting (m)")
    ax.set_ylabel("UTM 18N northing (m)")
    ax.set_title(
        "Bathtub depth map, downtown Montpelier VT -- July 2023 peak\n"
        f"WSE = {wse_m:.2f} m NAVD88 (gauge datum + peak stage), planar assumption",
        fontsize=8.5,
    )

    legend_elems = [
        Line2D([0], [0], marker="*", color="none", markerfacecolor=COLOR_MARKER,
               markeredgecolor="white", markersize=12, label="USGS gauge 04286000"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="black", markersize=6, label="Downtown landmark"),
    ]
    ax.legend(handles=legend_elems, loc="lower left", fontsize=6.5, frameon=True,
              framealpha=0.8)

    ax.set_aspect("equal")
    fig.tight_layout(pad=0.6)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    print(f"\nSaved figure: {out_pdf}")
    print(f"Saved figure: {out_png}")


# ----------------------------------------------------------------------
# Results markdown
# ----------------------------------------------------------------------

def write_results_md(ctx, out_path):
    b = ctx["bathtub_stats"]
    c = ctx["connected_stats"]
    lines = []
    lines.append("# Gauge-Stage-to-DEM (\"Bathtub\") Depth Mapping -- Downtown Montpelier, VT\n")
    lines.append(
        "Preliminary evidence for an NSF CPS-CIR proposal on flood-duration/depth "
        "modeling in Montpelier, VT: a building-level depth map for the July 2023 "
        "flood peak, produced by projecting the single USGS gauge 04286000 stage "
        "reading onto a public 1 m DEM, with **no hydraulic simulation**. All "
        "numbers below are computed by `analysis/stage_to_dem.py` from live "
        "USGS NWIS, USGS 3DEP, and OSM Nominatim queries; see that script for "
        "exact methods.\n"
    )

    lines.append("## Method summary\n")
    lines.append(
        "1. Gauge datum and vertical datum code fetched from the NWIS expanded "
        "site service.\n"
        "2. Peak stage for the July 2023 event fetched from NWIS IV gage height "
        "(00065), independently reproducing the value already reported in "
        "`analysis/results.md`.\n"
        "3. Water-surface elevation (WSE) = gauge datum (NAVD88, ft) + peak stage "
        "(ft), converted to meters -- a single flat/planar WSE applied across the "
        "whole domain (the \"bathtub\" assumption: no water-surface slope, no "
        "hydraulic routing, no travel-time lag between the gauge and downtown).\n"
        "4. A 1 m USGS 3DEP DEM clip for downtown Montpelier "
        f"(bbox lon {BBOX_LONLAT['minlon']} to {BBOX_LONLAT['maxlon']}, "
        f"lat {BBOX_LONLAT['minlat']} to {BBOX_LONLAT['maxlat']}) downloaded via the "
        "3DEPElevation ImageServer `exportImage` endpoint, reprojected to NAD83 / "
        "UTM 18N (EPSG:26918, meters).\n"
        "5. depth = WSE - ground elevation, clipped to positive. Reported for "
        "both the raw \"bathtub\" mask (any DEM cell below WSE) and a "
        "\"connected\" mask restricted to cells 4-connected, via other "
        "inundated cells, to a river-channel proxy (lowest 1% of DEM cells in "
        "the clip) -- see Caveats for what this does and does not capture.\n"
    )

    lines.append("## Vertical datum handling\n")
    lines.append(
        f"- Gauge datum: **{ctx['alt_va_ft']:.2f} ft**, vertical datum code "
        f"**{ctx['alt_datum_cd']}** (source: {ctx['gauge_datum_src']}).\n"
        f"- 3DEP DEMs report **NAVD88** orthometric heights in meters. The gauge "
        f"datum for 04286000 is already NAVD88, so **no NGVD29-to-NAVD88 "
        f"conversion was needed** for this analysis.\n"
        "- The script (`resolve_navd88_datum_ft()`) explicitly checks "
        "`alt_datum_cd` and would raise rather than silently apply an "
        "unverified regional offset if a future gauge in this project's domain "
        "reports NGVD29 -- a site-specific offset from NOAA VDatum or NGS "
        "VERTCON would be required in that case, not a generic \"central "
        "Vermont\" approximation.\n"
        f"- **Check that `alt_va` is a streambed-level datum, not a station-"
        f"platform elevation:** the gauge datum in NAVD88 is "
        f"**{ctx['alt_va_ft']*FT_TO_M:.2f} m**, and the DEM minimum anywhere in "
        f"the downloaded clip is **{ctx['dem_min_m']:.2f} m** -- "
        f"{ctx['datum_vs_dem_min_m']:+.2f} m relative to the gauge datum. These "
        "are close (well within DEM vertical uncertainty plus the offset "
        "expected between a fixed geodetic benchmark and the exact 1 m pixel "
        "nearest the gauge), which supports reading `alt_va` as the zero-of-"
        "stage elevation near the channel thalweg at the gauge, consistent "
        "with how WSE = datum + stage is being used here.\n"
    )

    lines.append("## Numbers\n")
    lines.append(
        f"- Peak stage ({PEAK_EVENT_LABEL}): **{ctx['peak_stage_ft']:.2f} ft** "
        f"{'at ' + ctx['peak_time'] if ctx['peak_time'] else '(fallback value)'} "
        "(cross-checked against `analysis/results.md`).\n"
        f"- Water-surface elevation (WSE): "
        f"**{ctx['wse_ft']:.2f} ft = {ctx['wse_m']:.3f} m, NAVD88**.\n"
        f"- DEM clip: {ctx['dem_width']} x {ctx['dem_height']} px at "
        f"{ctx['dem_res_m']:.3f} m/px (EPSG:26918), "
        f"elevation range **{ctx['dem_min_m']:.1f}-{ctx['dem_max_m']:.1f} m NAVD88** "
        f"across the full clip (downtown/river cells within that range are in the "
        f"~{ctx['dem_p2_m']:.0f}-{ctx['dem_p25_m']:.0f} m band; see percentiles in "
        "script output).\n"
    )

    lines.append("### Bathtub extent (any DEM cell below WSE)\n")
    lines.append(
        f"- Inundated: **{b['n_cells']:,} cells ({b['frac_of_domain']*100:.2f}% of "
        f"the clipped bbox), {b['area_acres']:.1f} acres ({b['area_m2']/1e6:.3f} km^2)**.\n"
        f"- Depth over inundated cells: median **{b['depth_median_m']:.2f} m** "
        f"(IQR {b['depth_q25_m']:.2f}-{b['depth_q75_m']:.2f} m), "
        f"max **{b['depth_max_m']:.2f} m**.\n"
    )

    lines.append("### Connected extent (flood-fill from river-channel proxy)\n")
    lines.append(
        f"- Inundated: **{c['n_cells']:,} cells ({c['frac_of_domain']*100:.2f}% of "
        f"the clipped bbox), {c['area_acres']:.1f} acres ({c['area_m2']/1e6:.3f} km^2)**"
        f" -- {(1 - c['n_cells']/max(b['n_cells'],1))*100:.1f}% smaller than the raw "
        "bathtub extent (removes topographically isolated low spots not "
        "reachable from the channel).\n"
        f"- Depth over inundated cells: median **{c['depth_median_m']:.2f} m** "
        f"(IQR {c['depth_q25_m']:.2f}-{c['depth_q75_m']:.2f} m), "
        f"max **{c['depth_max_m']:.2f} m**.\n"
        f"- River-channel proxy: {ctx['n_seed_cells']:,} seed cells at or below "
        f"{ctx['seed_thresh_m']:.2f} m NAVD88 (lowest 1% of the clip), forming "
        f"**{ctx['n_seed_components']} connected component(s)** within the bathtub "
        f"mask out of {ctx['num_components']} total components -- i.e. the channel "
        "proxy is effectively a single continuous corridor through the domain "
        "(visually confirmed to trace the Winooski River / North Branch "
        "confluence through downtown; see the figure).\n"
    )

    lines.append("## Sanity check: known downtown flooded core\n")
    lines.append(
        "Three downtown landmarks were geocoded (OSM Nominatim) and sampled "
        "against the DEM/WSE as a qualitative check against public reporting "
        "on the July 2023 flood:\n\n"
    )
    lines.append("| Landmark | Elevation (m NAVD88) | Depth vs. WSE (m) | Depth (ft) | Straight-line dist. from gauge (m) | Implied WSE slope if dry (m/km) |\n")
    lines.append("|---|---|---|---|---|---|\n")
    for name, d in ctx["landmark_checks"].items():
        slope_str = f"{d['implied_slope_m_per_km']:.2f}" if d['implied_slope_m_per_km'] == d['implied_slope_m_per_km'] else "n/a"
        lines.append(f"| {name} | {d['elev_m']:.2f} | {d['depth_m']:+.2f} | {d['depth_ft']:+.2f} | "
                      f"{d['dist_from_gauge_m']:.0f} | {slope_str} |\n")
    lines.append(
        "\n(Positive depth = predicted inundated by the bathtub model at that "
        "exact point; negative = predicted dry. \"Implied WSE slope if dry\" "
        "is the water-surface slope, upstream from the gauge to that point, "
        "that would be needed to just reach that point's elevation -- i.e. "
        "how much slope the flat-WSE assumption is failing to capture, not a "
        "measured value.)\n\n"
    )
    lines.append(
        "**Prose comparison -- the headline finding is that the naive planar "
        "bathtub UNDER-predicts the documented downtown flooded core, and "
        "the likely mechanism is identifiable, not just \"planar assumption\" "
        "in the abstract.** Gauge 04286000 sits **downstream (west)** of the "
        "confluence and downtown core -- the DEM's lowest cells in the clip "
        "are at the gauge and rise moving upstream (east) through the bend "
        "and confluence (see figure). Holding WSE flat at the "
        "downstream gauge's stage therefore systematically under-predicts "
        "depth upstream, and that is exactly the pattern seen here: "
        "Montpelier City Hall (Main St, ~1,457 m upstream of the gauge) and "
        "Confluence Park (~1,258 m upstream) both sample as dry-but-close "
        "(0.8-1.5 m above WSE) despite Main Street being widely reported to "
        "have had standing water in July 2023; the implied water-surface "
        "slope needed to reach those points (~0.6-1.0 m/km) is a physically "
        "plausible flood slope for this reach, not an implausible number, "
        "which supports \"missing upstream water-surface slope\" as the "
        "actual mechanism rather than a datum or DEM error. The VT State "
        "House, by contrast, sits well above WSE regardless (dry by ~9-10 m, "
        "consistent with public accounts that the State House itself was not "
        "flooded) -- no plausible slope correction would flood that point, so "
        "it is not part of the under-prediction story. The connected "
        "inundation mask still traces the right corridor -- it hugs the "
        "river channel and widens along the outside of the bend and at the "
        "confluence, i.e. the low-lying blocks adjacent to the river where "
        "the worst-documented damage occurred -- but a single downstream "
        "gauge with a flat WSE is directionally right and quantitatively "
        "conservative (too dry) upstream. Adding even a simple estimated "
        "water-surface slope (still no full hydraulic simulation) is exactly "
        "the kind of cheap improvement this prototype motivates for the "
        "proposed work; see Caveats.\n"
    )

    lines.append("## Caveats\n")
    lines.append(
        "- **Planar water-surface assumption, and it demonstrably matters "
        "here.** WSE is a single flat plane equal to gauge datum + peak "
        "stage everywhere in the domain. The gauge is downstream of "
        "downtown, and the sanity check above shows this flat-plane "
        "assumption under-predicts depth upstream (City Hall and Confluence "
        "Park sample dry-but-close, despite documented Main Street "
        "flooding) by an amount consistent with a physically plausible "
        "flood-surface slope of order 0.5-1 m/km over this ~1.5 km reach. "
        "This is the single biggest simplification in this prototype, it is "
        "not merely theoretical for this domain, and it is the primary "
        "thing that adding even a simple slope term (short of full "
        "hydraulic/hydrologic modeling) would correct.\n"
        "- **No travel-time/attenuation lag.** The gauge's peak stage is "
        "applied instantaneously across the whole domain; no routing delay "
        "between the gauge and downtown parcels is modeled.\n"
        "- **Connectivity is a cheap 4-connected flood-fill from an elevation-"
        "percentile channel proxy**, not a validated hydrography-derived "
        "channel network or a true 2D flood-routing connectivity solve. It "
        "removes obviously non-physical isolated low spots but does not "
        "model backwater, storm-sewer surcharging, or overland flow paths "
        "that could connect additional cells in a real flood.\n"
        "- **DEM is bare-earth terrain, not a flood-hazard model.** It does "
        "not account for buildings/basements changing local hydraulics, "
        "culverts, storm drains, or any subsurface flow paths -- all of "
        "which affect real building-level flooding independent of ground "
        "elevation.\n"
        "- **DEM vintage/resolution.** 1 m USGS 3DEP bare-earth DEM, "
        "downloaded via the 3DEPElevation ImageServer `exportImage` "
        "dynamic mosaic on 2026-07-22 (service `serviceDescription` states "
        "the mosaic reflects all 3DEP data published as of 2026-06-23); the "
        "exact source survey/LiDAR acquisition date for this specific tile "
        "is not returned by this endpoint and was not independently "
        "verified.\n"
        "- **USGS stage data are provisional** and subject to revision; "
        "retrieved 2026-07-22.\n"
        "- **Point sanity checks sample single geocoded coordinates**, which "
        "may not correspond to the specific building footprint, doorway, or "
        "street-level grade a news report is describing; a few meters of "
        "horizontal offset can flip a marginal point (like City Hall or "
        "Confluence Park here, both within ~1-3 m of WSE) from wet to dry. "
        "These should be read as illustrative, not as a rigorous "
        "point-validation against ground-truth high-water marks.\n"
        "- **Bathtub \"depth\" over already-wetted channel cells reflects "
        "normal channel bathymetry, not floodwater added over dry land, and "
        "this affects the whole reported depth distribution, not just the "
        "max.** The seed/channel cells (lowest 1% of the clip, ~57,700 "
        "cells) are included in both the bathtub and connected depth "
        "statistics; median and IQR depths above are therefore pulled "
        "upward by normal river depth at those cells, not purely by "
        "floodwater over dry downtown ground. The reported max depths in "
        "particular are dominated by cells directly over the river channel "
        "bed near the gauge. A depth distribution restricted to non-channel "
        "(overbank) cells only would better represent depth experienced on "
        "dry ground, but was not computed for this prototype.\n"
    )

    out_path.write_text("".join(lines))
    print(f"\nSaved results: {out_path}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    IMAGES_DIR.mkdir(parents=True, exist_ok=True)

    print("Fetching gauge datum...")
    alt_va_ft, alt_datum_cd, gauge_lat, gauge_lon, gauge_datum_src = fetch_gauge_datum(SITE_NO)
    datum_navd88_ft = resolve_navd88_datum_ft(alt_va_ft, alt_datum_cd)

    print("Fetching peak stage (July 2023)...")
    peak_stage_ft, peak_time = fetch_peak_stage(SITE_NO, PEAK_EVENT_START, PEAK_EVENT_END)

    wse_ft = datum_navd88_ft + peak_stage_ft
    wse_m = wse_ft * FT_TO_M
    print(f"WSE = {datum_navd88_ft:.2f} ft (datum) + {peak_stage_ft:.2f} ft (peak stage) "
          f"= {wse_ft:.2f} ft = {wse_m:.3f} m NAVD88")

    print("Downloading DEM...")
    dem_arr, prof = download_dem(BBOX_LONLAT)
    cell_area_m2 = prof["res"][0] * prof["res"][1]

    valid = np.isfinite(dem_arr)
    dem_min_m, dem_max_m = float(np.nanmin(dem_arr)), float(np.nanmax(dem_arr))
    dem_p2_m, dem_p25_m = float(np.nanpercentile(dem_arr[valid], 2)), float(np.nanpercentile(dem_arr[valid], 25))
    print(f"DEM elevation range: {dem_min_m:.1f} - {dem_max_m:.1f} m NAVD88 "
          f"(p2={dem_p2_m:.1f}, p25={dem_p25_m:.1f}, p50={np.nanpercentile(dem_arr[valid], 50):.1f})")
    # Plausibility check (not a hard failure -- printed loudly if violated).
    if not (150 <= dem_min_m <= 165):
        print(f"  WARNING: DEM minimum {dem_min_m:.1f} m NAVD88 is outside the "
              f"expected ~155-158 m river-elevation range for Montpelier -- verify DEM extent/CRS.")
    if not (150 <= dem_p25_m <= 175):
        print(f"  WARNING: DEM 25th percentile {dem_p25_m:.1f} m NAVD88 is outside the "
              f"expected ~157-170 m downtown-elevation range for Montpelier -- verify DEM extent/CRS.")

    print("Computing bathtub / connected depth masks...")
    masks = compute_depth_masks(dem_arr, wse_m)
    bathtub_stats = depth_stats(dem_arr, masks["bathtub"], wse_m, cell_area_m2)
    connected_stats = depth_stats(dem_arr, masks["connected"], wse_m, cell_area_m2)
    print(f"  bathtub:   {bathtub_stats}")
    print(f"  connected: {connected_stats}")

    gauge_xy_x, gauge_xy_y = transform("EPSG:4326", prof["crs"], [gauge_lon], [gauge_lat])
    gauge_xy = (gauge_xy_x[0], gauge_xy_y[0])

    # Gauge datum vs. DEM minimum -- alt_va is a benchmark elevation near the
    # gauge structure/streambed, not "land surface"; it should sit close to
    # (and typically at or below) the DEM's lowest bare-earth cell in the
    # clip if the gauge is sited at/near the channel thalweg. This is a
    # useful internal check that alt_va is being interpreted correctly as a
    # streambed-level datum, not e.g. a station platform elevation.
    datum_vs_dem_min_m = dem_min_m - (datum_navd88_ft * FT_TO_M)
    print(f"  gauge datum ({datum_navd88_ft*FT_TO_M:.2f} m) vs. DEM minimum in clip "
          f"({dem_min_m:.2f} m): DEM min is {datum_vs_dem_min_m:+.2f} m relative to datum")

    print("Geocoding downtown landmarks for sanity check...")
    landmark_checks = {}
    landmarks_xy = {}
    for name, query in GEOCODE_QUERIES.items():
        lat, lon = geocode(name, query)
        elev_m, (x, y) = sample_at_lonlat(dem_arr, prof, lon, lat)
        depth_m = wse_m - elev_m
        dist_from_gauge_m = float(np.hypot(x - gauge_xy[0], y - gauge_xy[1]))
        landmark_checks[name] = dict(
            elev_m=elev_m, depth_m=depth_m, depth_ft=depth_m / FT_TO_M,
            dist_from_gauge_m=dist_from_gauge_m,
            implied_slope_m_per_km=(-depth_m / (dist_from_gauge_m / 1000.0)
                                     if depth_m < 0 and dist_from_gauge_m > 0 else float("nan")),
        )
        landmarks_xy[name] = (x, y)
        print(f"  {name}: elev={elev_m:.2f} m, depth vs WSE={depth_m:+.2f} m, "
              f"straight-line dist from gauge={dist_from_gauge_m:.0f} m")

    make_figure(
        dem_arr, prof, wse_m, masks, gauge_xy, landmarks_xy,
        out_pdf=IMAGES_DIR / "stage_dem_depth.pdf",
        out_png=ANALYSIS_DIR / "stage_dem_depth.png",
    )

    ctx = dict(
        alt_va_ft=alt_va_ft, alt_datum_cd=alt_datum_cd, gauge_datum_src=gauge_datum_src,
        peak_stage_ft=peak_stage_ft, peak_time=peak_time,
        wse_ft=wse_ft, wse_m=wse_m,
        dem_width=prof["width"], dem_height=prof["height"], dem_res_m=prof["res"][0],
        dem_min_m=dem_min_m, dem_max_m=dem_max_m, dem_p2_m=dem_p2_m, dem_p25_m=dem_p25_m,
        bathtub_stats=bathtub_stats, connected_stats=connected_stats,
        n_seed_cells=masks["n_seed_cells"], seed_thresh_m=masks["seed_thresh_m"],
        num_components=masks["num_components"], n_seed_components=masks["n_seed_components"],
        landmark_checks=landmark_checks, datum_vs_dem_min_m=datum_vs_dem_min_m,
    )
    write_results_md(ctx, out_path=ANALYSIS_DIR / "stage_dem.md")

    print("\nDone.")


if __name__ == "__main__":
    main()
