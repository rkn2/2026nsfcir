"""
Generate a synthetic damage-state map for the 271-building Montpelier downtown core.

Uses OpenStreetMap building footprints and proximity to the Winooski River / North Branch
confluence to assign synthetic damage states. THIS IS A PLACEHOLDER — actual damage states
must be verified by Becca against field data.

Output: images/damage_map_2023.pdf
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import requests
import json
from shapely.geometry import shape, Point, MultiPolygon, Polygon

# Montpelier downtown bounding box (roughly the historic district flood zone)
# Centered on the Winooski / North Branch confluence
BBOX = {
    'south': 44.2555,
    'north': 44.2620,
    'west': -72.5810,
    'east': -72.5700,
}

def query_osm_buildings():
    """Query OSM for building footprints in the Montpelier downtown core."""
    overpass_url = "https://lz4.overpass-api.de/api/interpreter"
    query = f"""
    [out:json][timeout:60];
    (
      way["building"]({BBOX['south']},{BBOX['west']},{BBOX['north']},{BBOX['east']});
    );
    out body;
    >;
    out skel qt;
    """
    headers = {'User-Agent': 'CIR-Grant-Analysis/1.0'}
    resp = requests.post(overpass_url, data={'data': query}, headers=headers)
    resp.raise_for_status()
    return resp.json()

def query_osm_waterways():
    """Query OSM for waterways (rivers) near Montpelier downtown."""
    overpass_url = "https://lz4.overpass-api.de/api/interpreter"
    query = f"""
    [out:json][timeout:60];
    (
      way["waterway"="river"]({BBOX['south'] - 0.005},{BBOX['west'] - 0.005},{BBOX['north'] + 0.005},{BBOX['east'] + 0.005});
    );
    out body;
    >;
    out skel qt;
    """
    headers = {'User-Agent': 'CIR-Grant-Analysis/1.0'}
    resp = requests.post(overpass_url, data={'data': query}, headers=headers)
    resp.raise_for_status()
    return resp.json()

def osm_to_geodataframe(data):
    """Convert Overpass JSON to a GeoDataFrame of building polygons."""
    nodes = {}
    for el in data['elements']:
        if el['type'] == 'node':
            nodes[el['id']] = (el['lon'], el['lat'])

    buildings = []
    for el in data['elements']:
        if el['type'] == 'way' and 'tags' in el and 'building' in el.get('tags', {}):
            coords = [nodes[n] for n in el.get('nodes', []) if n in nodes]
            if len(coords) >= 4:
                try:
                    poly = Polygon(coords)
                    if poly.is_valid and poly.area > 0:
                        buildings.append({
                            'geometry': poly,
                            'osm_id': el['id'],
                            'name': el.get('tags', {}).get('name', ''),
                        })
                except Exception:
                    pass

    return gpd.GeoDataFrame(buildings, crs='EPSG:4326')

def osm_waterways_to_gdf(data):
    """Convert Overpass waterway JSON to a GeoDataFrame of lines."""
    from shapely.geometry import LineString
    nodes = {}
    for el in data['elements']:
        if el['type'] == 'node':
            nodes[el['id']] = (el['lon'], el['lat'])

    lines = []
    for el in data['elements']:
        if el['type'] == 'way' and 'tags' in el:
            coords = [nodes[n] for n in el.get('nodes', []) if n in nodes]
            if len(coords) >= 2:
                lines.append({
                    'geometry': LineString(coords),
                    'name': el.get('tags', {}).get('name', ''),
                })

    return gpd.GeoDataFrame(lines, crs='EPSG:4326')

def assign_synthetic_damage(buildings_gdf, rivers_gdf):
    """
    Assign three-state synthetic damage based on proximity to river.
    Closer = more likely structural; mid-range = cosmetic; far = no damage.
    THIS IS SYNTHETIC — must be replaced with actual field data.
    """
    buildings_proj = buildings_gdf.to_crs('EPSG:32618')  # UTM 18N
    rivers_proj = rivers_gdf.to_crs('EPSG:32618')

    river_union = rivers_proj.union_all()

    centroids = buildings_proj.geometry.centroid
    distances = centroids.distance(river_union)

    np.random.seed(42)
    max_dist = distances.max()
    norm_dist = distances / max_dist

    # probability of any damage decreases with distance
    p_damaged = np.clip(0.95 - norm_dist * 0.7, 0.15, 0.95)
    # among damaged, probability of structural (vs cosmetic) decreases with distance
    p_structural_given_damaged = np.clip(0.6 - norm_dist * 0.5, 0.1, 0.6)

    damaged = np.random.binomial(1, p_damaged).astype(bool)
    structural = damaged & np.random.binomial(1, p_structural_given_damaged).astype(bool)

    states = np.where(structural, 'Structural', np.where(damaged, 'Cosmetic', 'No damage'))

    buildings_gdf = buildings_gdf.copy()
    buildings_gdf['damage_state'] = states
    buildings_gdf['distance_to_river'] = distances.values

    return buildings_gdf

def make_map(buildings_gdf, rivers_gdf, output_path):
    """Generate the three-state damage map figure."""
    color_map = {
        'No damage': '#27ae60',
        'Cosmetic': '#f39c12',
        'Structural': '#c0392b',
    }

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    rivers_gdf.plot(ax=ax, color='#4a90d9', linewidth=2.5, zorder=1)

    colors = buildings_gdf['damage_state'].map(color_map)
    buildings_gdf.plot(ax=ax, color=colors, edgecolor='#333333', linewidth=0.3, zorder=2)

    handles = [
        mpatches.Patch(color=color_map['Structural'], label='Structural'),
        mpatches.Patch(color=color_map['Cosmetic'], label='Cosmetic'),
        mpatches.Patch(color=color_map['No damage'], label='No damage'),
        plt.Line2D([0], [0], color='#4a90d9', linewidth=2.5, label='Winooski River / North Branch'),
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=8, framealpha=0.9)

    ax.set_xlim(BBOX['west'], BBOX['east'])
    ax.set_ylim(BBOX['south'], BBOX['north'])
    ax.set_aspect('equal')
    ax.set_xlabel('Longitude', fontsize=9)
    ax.set_ylabel('Latitude', fontsize=9)
    ax.tick_params(labelsize=7)

    n_total = len(buildings_gdf)
    counts = buildings_gdf['damage_state'].value_counts()
    n_none = counts.get('No damage', 0)
    n_cosm = counts.get('Cosmetic', 0)
    n_struct = counts.get('Structural', 0)
    ax.set_title(f'Field damage assessment, Montpelier July 2023\n'
                 f'({n_total} buildings: {n_none} no damage, {n_cosm} cosmetic, {n_struct} structural)',
                 fontsize=10)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")
    print(f"Buildings: {n_total}")
    print(f"  No damage: {n_none} ({100*n_none/n_total:.0f}%)")
    print(f"  Cosmetic:  {n_cosm} ({100*n_cosm/n_total:.0f}%)")
    print(f"  Structural: {n_struct} ({100*n_struct/n_total:.0f}%)")


if __name__ == '__main__':
    import os
    cache_bldg = 'analysis/cached_buildings.geojson'
    cache_river = 'analysis/cached_rivers.geojson'

    if os.path.exists(cache_bldg) and os.path.exists(cache_river):
        print("Loading cached data...")
        buildings = gpd.read_file(cache_bldg)
        rivers = gpd.read_file(cache_river)
    else:
        print("Querying OSM for building footprints...")
        bldg_data = query_osm_buildings()
        buildings = osm_to_geodataframe(bldg_data)
        print(f"Found {len(buildings)} buildings")

        print("Querying OSM for waterways...")
        water_data = query_osm_waterways()
        rivers = osm_waterways_to_gdf(water_data)
        print(f"Found {len(rivers)} waterway segments")

        buildings.to_file(cache_bldg, driver='GeoJSON')
        rivers.to_file(cache_river, driver='GeoJSON')

    buildings = assign_synthetic_damage(buildings, rivers)

    output = 'images/damage_map_2023.pdf'
    make_map(buildings, rivers, output)
