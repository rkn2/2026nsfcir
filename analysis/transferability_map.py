#!/usr/bin/env python3
"""
Transferability map: bivariate choropleth of pre-1940 housing stock
vs. FEMA flood disaster declarations, with USGS gauge overlay.
Counties high on BOTH dimensions are the framework's target population.

Data sources (all free, no API key):
  - ACS 5-year 2022 table B25034 (year structure built), county-level bulk CSV
  - FEMA OpenFEMA API: flood disaster declarations by county
  - Census Bureau cartographic boundary shapefiles (5m resolution)
  - USGS Water Services: active streamgauges with stage data

Output: images/transferability_map.pdf
"""

import io
import warnings
import zipfile
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import requests

warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parent.parent / "images" / "transferability_map.pdf"
CACHE = Path(__file__).resolve().parent / ".cache"
CACHE.mkdir(exist_ok=True)


def fetch_acs_b25034():
    cache_path = CACHE / "b25034_county.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path)

    url = (
        "https://www2.census.gov/programs-surveys/acs/summary_file/"
        "2022/table-based-SF/data/5YRData/acsdt5y2022-b25034.dat"
    )
    print("Downloading ACS B25034 data...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()

    df = pd.read_csv(io.StringIO(r.text), sep="|")
    df = df[df["GEO_ID"].str.match(r"^0500000US\d{5}$")].copy()
    df["GEOID"] = df["GEO_ID"].str[-5:]
    df["total_units"] = pd.to_numeric(df["B25034_E001"], errors="coerce")
    df["pre1940"] = pd.to_numeric(df["B25034_E011"], errors="coerce")
    df = df[["GEOID", "total_units", "pre1940"]].dropna()
    df["pct_pre1940"] = df["pre1940"] / df["total_units"] * 100
    df.to_csv(cache_path, index=False)
    return df


def fetch_flood_declarations():
    cache_path = CACHE / "flood_declarations_by_county.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path, dtype={"GEOID": str})

    raise FileNotFoundError(
        "Run the curl-based download first (flood_declarations_by_county.csv)"
    )


def fetch_county_shapes():
    cache_path = CACHE / "cb_county.gpkg"
    if cache_path.exists():
        return gpd.read_file(cache_path)

    url = (
        "https://www2.census.gov/geo/tiger/GENZ2022/shp/"
        "cb_2022_us_county_5m.zip"
    )
    print("Downloading county boundaries...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()

    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        zf.extractall(CACHE / "county_shp")
    shp = list((CACHE / "county_shp").glob("*.shp"))[0]
    gdf = gpd.read_file(shp)
    gdf.to_file(cache_path, driver="GPKG")
    return gdf


def fetch_state_shapes():
    cache_path = CACHE / "cb_state.gpkg"
    if cache_path.exists():
        return gpd.read_file(cache_path)

    url = (
        "https://www2.census.gov/geo/tiger/GENZ2022/shp/"
        "cb_2022_us_state_5m.zip"
    )
    print("Downloading state boundaries...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()

    with zipfile.ZipFile(io.BytesIO(r.content)) as zf:
        zf.extractall(CACHE / "state_shp")
    shp = list((CACHE / "state_shp").glob("*.shp"))[0]
    gdf = gpd.read_file(shp)
    gdf.to_file(cache_path, driver="GPKG")
    return gdf


def fetch_usgs_gauges():
    cache_path = CACHE / "usgs_gauges.csv"
    if cache_path.exists():
        return pd.read_csv(cache_path)

    url = "https://waterservices.usgs.gov/nwis/site/"
    print("Downloading USGS gauge locations...")
    all_rows = []
    all_states = [
        "09", "23", "25", "33", "34", "36", "42", "44", "50",
        "10", "11", "12", "13", "17", "18", "19", "20", "21",
        "22", "24", "26", "27", "28", "29", "37", "39", "45",
        "47", "51", "54", "55",
        "01", "02", "04", "05", "06", "08", "15", "16", "30",
        "31", "32", "35", "38", "40", "41", "46", "48", "49", "53", "56",
    ]
    for sc in all_states:
        try:
            r = requests.get(
                url,
                params={
                    "format": "rdb",
                    "stateCd": sc,
                    "siteType": "ST",
                    "siteStatus": "active",
                    "parameterCd": "00065",
                    "hasDataTypeCd": "iv",
                },
                timeout=60,
            )
            if r.status_code != 200:
                continue
            for line in r.text.splitlines():
                if line.startswith("USGS"):
                    parts = line.split("\t")
                    if len(parts) >= 6:
                        try:
                            all_rows.append(
                                {"lat": float(parts[4]), "lon": float(parts[5])}
                            )
                        except (ValueError, IndexError):
                            pass
        except Exception:
            continue

    df = pd.DataFrame(all_rows)
    print(f"  Found {len(df)} active gauges")
    df.to_csv(cache_path, index=False)
    return df


def make_bivariate_cmap():
    """Build a 3x3 bivariate color matrix.

    Rows (bottom→top): flood declarations (low → high)
    Cols (left→right): pre-1940 housing % (low → high)

    Low-low is light grey. High housing = amber tones.
    High flood = blue tones. High-both = deep brown/violet.
    """
    # 3x3 grid: [flood_bin][housing_bin]
    # Inspired by Joshua Stevens' bivariate palette
    colors = np.array([
        # flood low
        ["#e8e8e8", "#dfb0a0", "#be6450"],
        # flood mid
        ["#a0b8d0", "#b0909c", "#a05060"],
        # flood high
        ["#4885a8", "#706080", "#6c3040"],
    ])
    return colors


def classify_bivariate(housing_pct, flood_count):
    """Classify into 3x3 bins. Returns (housing_bin, flood_bin) each 0-2."""
    # Housing: <10%, 10-25%, >25%
    h = np.where(housing_pct < 10, 0, np.where(housing_pct < 25, 1, 2))
    # Flood: 0 declarations, 1-3, 4+
    f = np.where(flood_count < 1, 0, np.where(flood_count < 4, 1, 2))
    return h, f


def main():
    acs = fetch_acs_b25034()
    floods = fetch_flood_declarations()
    counties = fetch_county_shapes()
    states = fetch_state_shapes()
    gauges = fetch_usgs_gauges()

    # CONUS filter
    conus_fips = (
        {f"{i:02d}" for i in range(1, 57)}
        - {"02", "15", "60", "66", "69", "72", "78"}
    )
    counties["STATEFP_str"] = counties["STATEFP"].astype(str).str.zfill(2)
    counties = counties[counties["STATEFP_str"].isin(conus_fips)].copy()
    states["STATEFP_str"] = states["STATEFP"].astype(str).str.zfill(2)
    states = states[states["STATEFP_str"].isin(conus_fips)].copy()

    # Join data
    counties["GEOID"] = counties["GEOID"].astype(str).str.zfill(5)
    acs["GEOID"] = acs["GEOID"].astype(str).str.zfill(5)
    floods["GEOID"] = floods["GEOID"].astype(str).str.zfill(5)

    merged = counties.merge(acs[["GEOID", "pct_pre1940"]], on="GEOID", how="left")
    merged = merged.merge(floods, on="GEOID", how="left")
    merged["flood_declarations"] = merged["flood_declarations"].fillna(0)
    merged["pct_pre1940"] = merged["pct_pre1940"].fillna(0)

    # Bivariate classification
    h_bin, f_bin = classify_bivariate(
        merged["pct_pre1940"].values, merged["flood_declarations"].values
    )
    merged["h_bin"] = h_bin
    merged["f_bin"] = f_bin
    merged["biv_class"] = merged["f_bin"] * 3 + merged["h_bin"]

    biv_colors = make_bivariate_cmap()
    color_list = biv_colors.flatten().tolist()  # 9 colors, row-major
    merged["color"] = [color_list[c] for c in merged["biv_class"]]

    # Reproject
    merged = merged.to_crs("ESRI:102003")
    states = states.to_crs("ESRI:102003")

    if len(gauges) > 0:
        gauges = gauges[
            (gauges["lat"] > 24) & (gauges["lat"] < 50)
            & (gauges["lon"] > -125) & (gauges["lon"] < -66)
        ].copy()
        gauges_gdf = gpd.GeoDataFrame(
            gauges,
            geometry=gpd.points_from_xy(gauges["lon"], gauges["lat"]),
            crs="EPSG:4326",
        ).to_crs("ESRI:102003")

    # --- Figure: map on top, legend + stats row below ---
    fig = plt.figure(figsize=(7.5, 5.5), dpi=300)
    ax = fig.add_axes([0.02, 0.20, 0.96, 0.78])

    # Plot bivariate choropleth
    merged.plot(ax=ax, color=merged["color"], edgecolor="face", linewidth=0.08)
    states.boundary.plot(ax=ax, edgecolor="#444444", linewidth=0.3)

    # Gauges — cross markers visible on any background
    if len(gauges) > 0:
        ax.scatter(
            gauges_gdf.geometry.x, gauges_gdf.geometry.y,
            marker="+", c="#000000", s=2, linewidths=0.2,
            alpha=0.3, zorder=5,
        )

    # Montpelier
    montpelier = gpd.GeoDataFrame(
        [{"lat": 44.2601, "lon": -72.5754}],
        geometry=gpd.points_from_xy([-72.5754], [44.2601]),
        crs="EPSG:4326",
    ).to_crs("ESRI:102003")
    montpelier.plot(ax=ax, color="black", markersize=40, marker="*", zorder=10)
    ax.annotate(
        "Montpelier, VT",
        xy=montpelier.geometry.iloc[0].coords[0],
        xytext=(12, -8),
        textcoords="offset points",
        fontsize=5.5,
        fontweight="bold",
        ha="left",
    )

    ax.set_xlim(merged.total_bounds[0] - 100000, merged.total_bounds[2] + 100000)
    ax.set_ylim(merged.total_bounds[1] - 100000, merged.total_bounds[3] + 100000)
    ax.set_axis_off()

    # --- Bivariate legend + stats, all via fig.text/fig.patches ---
    # Use figure coordinates directly to avoid axes clipping issues

    both_high = merged[(merged["h_bin"] == 2) & (merged["f_bin"] >= 1)]
    target = merged[(merged["h_bin"] >= 1) & (merged["f_bin"] >= 1)]
    high_high = merged[(merged["h_bin"] == 2) & (merged["f_bin"] == 2)]
    n_gauged = len(gauges) if len(gauges) > 0 else 0
    print(f"\nBivariate stats:")
    print(f"  High housing + any flood history: {len(both_high)} counties")
    print(f"  Mid+ housing + any flood history: {len(target)} counties")
    print(f"  High-high (h=2, f=2): {len(high_high)} counties")

    # Grid position in figure coords
    gx0, gy0 = 0.10, 0.04   # bottom-left of grid
    cw, ch = 0.04, 0.035     # cell width/height in figure coords

    for fi in range(3):
        for hi in range(3):
            fig.patches.append(mpatches.FancyBboxPatch(
                (gx0 + hi * cw, gy0 + fi * ch), cw, ch,
                boxstyle="square,pad=0",
                facecolor=biv_colors[fi, hi],
                edgecolor="white", linewidth=0.8,
                transform=fig.transFigure, figure=fig,
            ))

    # X-axis bin labels
    for i, lbl in enumerate(["<10%", "10–25%", ">25%"]):
        fig.text(gx0 + (i + 0.5) * cw, gy0 - 0.008, lbl,
                 fontsize=5.5, ha="center", va="top", color="#555555")
    fig.text(gx0 + 1.5 * cw, gy0 - 0.035,
             "Pre-1940 housing stock →",
             fontsize=6, ha="center", va="top",
             fontweight="bold", color="#333333")

    # Y-axis bin labels
    for i, lbl in enumerate(["0", "1–3", "4+"]):
        fig.text(gx0 - 0.01, gy0 + (i + 0.5) * ch, lbl,
                 fontsize=5.5, ha="right", va="center", color="#555555")
    fig.text(gx0 + 1.5 * cw, gy0 + 3 * ch + 0.008,
             "Flood declarations ↑",
             fontsize=6, ha="center", va="bottom",
             fontweight="bold", color="#333333")

    # Stats text
    stats_text = (
        f"{len(target):,} counties with $\\geq$10% pre-1940 stock "
        f"and flood history;  "
        f"{len(high_high):,} with $\\geq$25% and 4+ declarations.\n"
        f"+ marks {n_gauged:,} active USGS streamgauges.  "
        f"★ Montpelier, VT."
    )
    fig.text(
        0.30, 0.07, stats_text,
        fontsize=5.5, va="center", ha="left",
        color="#333333",
    )

    fig.savefig(OUT, bbox_inches="tight", pad_inches=0.15, dpi=300)
    print(f"\nSaved to {OUT}")
    plt.close()


if __name__ == "__main__":
    main()
