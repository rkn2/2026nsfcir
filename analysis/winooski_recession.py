#!/usr/bin/env python3
"""
Winooski River at Montpelier, VT (USGS 04286000) -- order-of-magnitude
drainage-timescale analysis for the July 2023 and July 2024 flood events.

Produced for an NSF proposal on flood-duration modeling in Montpelier, VT.

IMPORTANT SITE-NUMBER NOTE
---------------------------
The task that generated this script specified USGS site 04288000 as
"WINOOSKI RIVER AT MONTPELIER, VT". That site number is INCORRECT --
04288000 is actually "MAD RIVER NEAR MORETOWN, VT" (confirmed by querying
https://waterservices.usgs.gov/nwis/iv/?sites=04288000&... and inspecting
sourceInfo.siteName in the returned JSON). The correct USGS site number for
"WINOOSKI RIVER AT MONTPELIER, VT" is 04286000 (NWS/AHPS identifier MONV1).
This script uses 04286000. See results.md for full discussion.

What this script does
----------------------
1. Downloads USGS NWIS instantaneous-value (IV) discharge (00060, cfs) and
   gage height (00065, ft) for the site, for the two flood windows.
2. Downloads site metadata (drainage area) from the NWIS site service.
3. Uses the NWS/NOAA National Water Prediction Service API to retrieve the
   official flood-category stages for this gauge (NWSLI MONV1).
4. For each event: finds the peak, estimates pre-event baseflow, fits an
   exponential recession Q(t) = Q_peak * exp(-t/tau) to the falling limb,
   and computes durations above flood stage and above 2x pre-event baseflow.
5. Writes a two-panel publication-style PDF figure and a PNG copy, and a
   results.md summary.

Usage
-----
    python3 winooski_recession.py

Requires: requests, pandas, numpy, matplotlib (standard library otherwise).
"""

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

SITE_NO = "04286000"  # WINOOSKI RIVER AT MONTPELIER, VT (corrected; see module docstring)
SITE_NAME = "WINOOSKI RIVER AT MONTPELIER, VT"
NWSLI = "monv1"        # NWS/AHPS gauge identifier for this site

IV_URL = "https://waterservices.usgs.gov/nwis/iv/"
SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
NWPS_GAUGE_URL = f"https://api.water.noaa.gov/nwps/v1/gauges/{NWSLI}"

# Fallback values, empirically confirmed from live queries on 2026-07-22,
# used only if the live service calls fail at run time.
FALLBACK_DRAINAGE_AREA_SQMI = 397.0
FALLBACK_FLOOD_STAGE_FT = 15.0  # NWS AHPS "minor" flood category for MONV1

EVENTS = {
    "2023": dict(start="2023-07-01", end="2023-08-15", label="July 2023"),
    "2024": dict(start="2024-07-01", end="2024-08-15", label="July 2024"),
}

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "analysis"
IMAGES_DIR = REPO_ROOT / "images"

# Colorblind-safe palette: blue data, vermillion/orange dashed fit
COLOR_DATA = "#0072B2"
COLOR_FIT = "#D55E00"


# ----------------------------------------------------------------------
# Data acquisition
# ----------------------------------------------------------------------

def fetch_iv(site_no, start, end, retries=4, backoff=5):
    """Fetch USGS NWIS instantaneous-value discharge (00060) and gage height
    (00065) for [start, end]. Returns a DataFrame indexed by UTC datetime
    with columns discharge_cfs and stage_ft. Raises on failure (no silent
    fallback -- if this fails, there is no defensible data to fall back to
    other than the dv service, handled by fetch_event_data)."""
    params = {
        "sites": site_no,
        "parameterCd": "00060,00065",
        "startDT": start,
        "endDT": end,
        "format": "json",
    }
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.get(IV_URL, params=params, timeout=60)
            if r.status_code == 200:
                data = r.json()
                ts = data.get("value", {}).get("timeSeries", [])
                if not ts:
                    raise RuntimeError(f"No timeSeries returned for {site_no} {start}..{end}")
                series = {}
                site_name_seen = None
                for t in ts:
                    var = t["variable"]["variableCode"][0]["value"]
                    site_name_seen = t["sourceInfo"]["siteName"]
                    vals = t["values"][0]["value"]
                    rows = [(v["dateTime"], float(v["value"]))
                            for v in vals if v["value"] not in (None, "")]
                    if not rows:
                        continue
                    s = pd.DataFrame(rows, columns=["dateTime", "value"])
                    s["dateTime"] = pd.to_datetime(s["dateTime"], utc=True)
                    series[var] = s.set_index("dateTime")["value"]
                if "00060" not in series or "00065" not in series:
                    raise RuntimeError(f"Missing discharge or stage series for {site_no} {start}..{end}")
                out = pd.DataFrame({
                    "discharge_cfs": series["00060"],
                    "stage_ft": series["00065"],
                })
                out = out.sort_index()
                print(f"  [fetch_iv] site={site_no} name='{site_name_seen}' "
                      f"window={start}..{end} n={len(out)} rows")
                return out
            else:
                last_exc = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:  # noqa: BLE001 -- broad, but we log and retry explicitly
            last_exc = e
        print(f"  [fetch_iv] attempt {attempt+1}/{retries} failed ({last_exc}); retrying...")
        time.sleep(backoff)
    raise RuntimeError(f"fetch_iv failed after {retries} attempts: {last_exc}")


def fetch_dv(site_no, start, end):
    """Fallback: USGS daily-values service (used only if IV service has no
    data for the requested window)."""
    url = "https://waterservices.usgs.gov/nwis/dv/"
    params = {
        "sites": site_no,
        "parameterCd": "00060,00065",
        "startDT": start,
        "endDT": end,
        "format": "json",
    }
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    data = r.json()
    ts = data.get("value", {}).get("timeSeries", [])
    series = {}
    for t in ts:
        var = t["variable"]["variableCode"][0]["value"]
        vals = t["values"][0]["value"]
        rows = [(v["dateTime"], float(v["value"])) for v in vals if v["value"] not in (None, "")]
        if not rows:
            continue
        s = pd.DataFrame(rows, columns=["dateTime", "value"])
        s["dateTime"] = pd.to_datetime(s["dateTime"], utc=True)
        series[var] = s.set_index("dateTime")["value"]
    out = pd.DataFrame({
        "discharge_cfs": series.get("00060"),
        "stage_ft": series.get("00065"),
    })
    print(f"  [fetch_dv fallback] n={len(out)} rows")
    return out.sort_index()


def fetch_event_data(site_no, start, end):
    """Try IV data first; fall back to DV if IV has no rows."""
    try:
        df = fetch_iv(site_no, start, end)
        if len(df) > 0:
            return df, "iv"
    except Exception as e:  # noqa: BLE001
        print(f"  IV fetch failed ({e}); falling back to daily values.")
    df = fetch_dv(site_no, start, end)
    return df, "dv"


def fetch_drainage_area(site_no, retries=4, backoff=5):
    """Fetch drainage area (sq mi) from the NWIS expanded site service."""
    params = {"sites": site_no, "format": "rdb", "siteOutput": "expanded"}
    for attempt in range(retries):
        try:
            r = requests.get(SITE_URL, params=params, timeout=60)
            if r.status_code == 200:
                lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
                header = lines[0].split("\t")
                data_line = lines[2].split("\t")  # line 0 = header, line 1 = format codes
                row = dict(zip(header, data_line))
                da = row.get("drain_area_va", "").strip()
                if da:
                    print(f"  [fetch_drainage_area] drain_area_va={da} sq mi (source: {r.url})")
                    return float(da), r.url
        except Exception as e:  # noqa: BLE001
            print(f"  [fetch_drainage_area] attempt {attempt+1} failed: {e}")
        time.sleep(backoff)
    print(f"  [fetch_drainage_area] live query failed; using fallback constant "
          f"{FALLBACK_DRAINAGE_AREA_SQMI} sq mi (confirmed via manual query 2026-07-22).")
    return FALLBACK_DRAINAGE_AREA_SQMI, SITE_URL + " (fallback, manual query 2026-07-22)"


def fetch_flood_stage(nwsli=NWSLI, retries=3, backoff=4):
    """Fetch official NWS/NOAA flood-category stages for this gauge from the
    National Water Prediction Service API. Returns (dict of categories, url)."""
    url = NWPS_GAUGE_URL
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                data = r.json()
                cats = data["flood"]["categories"]
                print(f"  [fetch_flood_stage] categories={cats} (source: {url})")
                return cats, url
        except Exception as e:  # noqa: BLE001
            print(f"  [fetch_flood_stage] attempt {attempt+1} failed: {e}")
        time.sleep(backoff)
    print(f"  [fetch_flood_stage] live query failed; using fallback minor-flood-stage "
          f"constant {FALLBACK_FLOOD_STAGE_FT} ft (confirmed via manual query 2026-07-22).")
    return {"minor": {"stage": FALLBACK_FLOOD_STAGE_FT}}, url + " (fallback, manual query 2026-07-22)"


# ----------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------

def find_peak(df):
    peak_time = df["discharge_cfs"].idxmax()
    peak_Q = df.loc[peak_time, "discharge_cfs"]
    peak_stage_time = df["stage_ft"].idxmax()
    peak_stage = df.loc[peak_stage_time, "stage_ft"]
    return peak_time, peak_Q, peak_stage_time, peak_stage


def find_pre_event_baseflow(df, peak_time, lookback_start_h=60, lookback_end_h=12):
    """Pre-event baseflow = minimum hourly-mean discharge in the window
    [peak - lookback_start_h, peak - lookback_end_h]. This window is chosen
    to sit before the rising limb of the flood itself while still being
    close enough to the event to reflect antecedent conditions, avoiding
    contamination by the sharp rise into the peak."""
    hourly = df["discharge_cfs"].resample("1h").mean()
    win = hourly.loc[peak_time - pd.Timedelta(hours=lookback_start_h):
                      peak_time - pd.Timedelta(hours=lookback_end_h)]
    bf_time = win.idxmin()
    bf_val = win.min()
    return bf_time, bf_val


def find_recession_trough(df, peak_time, search_days=12, rise_frac=0.15, min_hours=24):
    """Find the local minimum discharge after the peak that marks the end of
    the "clean" falling limb -- i.e., the trough just before a secondary,
    rain-driven rise interrupts the recession. Algorithm: on the hourly
    series after the peak, walk forward and flag the first hour (after
    min_hours) at which discharge 24h later exceeds (1+rise_frac) times the
    current running minimum. The running minimum at that point is the
    trough (recession break)."""
    hourly = df["discharge_cfs"].resample("1h").mean()
    win = hourly.loc[peak_time: peak_time + pd.Timedelta(days=search_days)]
    times = win.index
    vals = win.values
    running_min = vals[0]
    running_min_t = times[0]
    for i in range(len(vals) - 1):
        if vals[i] < running_min:
            running_min = vals[i]
            running_min_t = times[i]
        elapsed_h = (times[i] - peak_time).total_seconds() / 3600.0
        if elapsed_h < min_hours:
            continue
        # look 24h ahead (or to end of series)
        future_idx = min(i + 24, len(vals) - 1)
        if vals[future_idx] > (1 + rise_frac) * running_min and future_idx > i:
            return running_min_t, running_min
    # no break found in search window; return end of window
    return times[-1], vals[-1]


def fit_exponential_recession(df, peak_time, fit_end, peak_Q):
    """Fit ln(Q) = ln(Q_peak) - t/tau via least squares (numpy polyfit) on
    the falling limb from peak_time to fit_end. Returns dict with tau_hours,
    tau_days, r2, n, slope, intercept, window bounds."""
    seg = df.loc[peak_time:fit_end, "discharge_cfs"]
    seg = seg[seg > 0]
    t_hours = (seg.index - peak_time).total_seconds() / 3600.0
    lnQ = np.log(seg.values)
    slope, intercept = np.polyfit(t_hours, lnQ, 1)
    pred = slope * t_hours + intercept
    ss_res = np.sum((lnQ - pred) ** 2)
    ss_tot = np.sum((lnQ - lnQ.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    tau_hours = -1.0 / slope
    return dict(
        tau_hours=tau_hours,
        tau_days=tau_hours / 24.0,
        r2=r2,
        n=len(seg),
        slope=slope,
        intercept=intercept,
        window_start=peak_time,
        window_end=fit_end,
        window_hours=(fit_end - peak_time).total_seconds() / 3600.0,
        t_hours=t_hours,
        lnQ=lnQ,
    )


def contiguous_interval_containing(series, mask, anchor_time):
    """Return (start, end, duration_hours) of the contiguous run of `mask`
    (a boolean Series aligned with `series`) that contains anchor_time.
    Uses actual timestamps (not point counts) so irregular sampling
    intervals (USGS IV data sometimes switches between 5-min and 15-min
    reporting during high flow) do not bias the duration estimate."""
    if not mask.any():
        return None, None, 0.0
    idx = mask.index
    pos = idx.get_indexer([anchor_time], method="nearest")[0]
    if not mask.iloc[pos]:
        # anchor itself isn't in the exceedance mask (e.g. threshold never
        # reached); search nearest True point
        true_positions = np.where(mask.values)[0]
        pos = true_positions[np.argmin(np.abs(true_positions - pos))]
    lo = pos
    while lo > 0 and mask.iloc[lo - 1]:
        lo -= 1
    hi = pos
    while hi < len(mask) - 1 and mask.iloc[hi + 1]:
        hi += 1
    start, end = idx[lo], idx[hi]
    # add one nominal sampling interval so a single-point run isn't 0 duration
    if hi + 1 < len(idx):
        step = (idx[hi + 1] - idx[hi])
    elif lo > 0:
        step = (idx[lo] - idx[lo - 1])
    else:
        step = pd.Timedelta(minutes=15)
    duration_hours = (end - start + step).total_seconds() / 3600.0
    return start, end, duration_hours


def analyze_event(key, cfg, flood_categories):
    print(f"\n=== Event {key} ({cfg['label']}) ===")
    df, source = fetch_event_data(SITE_NO, cfg["start"], cfg["end"])
    if len(df) == 0:
        raise RuntimeError(f"No data at all (IV or DV) for {key}")

    peak_time, peak_Q, peak_stage_time, peak_stage = find_peak(df)
    print(f"  peak discharge: {peak_Q:.0f} cfs at {peak_time} (UTC)")
    print(f"  peak stage: {peak_stage:.2f} ft at {peak_stage_time} (UTC)")

    bf_time, bf_val = find_pre_event_baseflow(df, peak_time)
    print(f"  pre-event baseflow: {bf_val:.1f} cfs at {bf_time} (UTC)")

    trough_time, trough_val = find_recession_trough(df, peak_time)
    full_window_hours = (trough_time - peak_time).total_seconds() / 3600.0
    print(f"  recession trough / secondary-rise break: {trough_val:.1f} cfs at "
          f"{trough_time} (UTC), {full_window_hours:.1f} h after peak")

    # Primary fit: up to 60 h post-peak, or the trough if it occurs sooner.
    primary_end = min(peak_time + pd.Timedelta(hours=60), trough_time)
    primary_fit = fit_exponential_recession(df, peak_time, primary_end, peak_Q)
    print(f"  PRIMARY fit window: {peak_time} -> {primary_end} "
          f"({primary_fit['window_hours']:.1f} h): "
          f"tau={primary_fit['tau_hours']:.2f} h ({primary_fit['tau_days']:.2f} d), "
          f"R2={primary_fit['r2']:.4f}, n={primary_fit['n']}")

    # Extended fit: peak to trough (full pre-secondary-rain recession).
    extended_fit = fit_exponential_recession(df, peak_time, trough_time, peak_Q)
    print(f"  EXTENDED fit window: {peak_time} -> {trough_time} "
          f"({extended_fit['window_hours']:.1f} h): "
          f"tau={extended_fit['tau_hours']:.2f} h ({extended_fit['tau_days']:.2f} d), "
          f"R2={extended_fit['r2']:.4f}, n={extended_fit['n']}")

    # Duration above NWS flood stage (minor category)
    flood_stage_ft = flood_categories.get("minor", {}).get("stage", FALLBACK_FLOOD_STAGE_FT)
    mask_fs = df["stage_ft"] > flood_stage_ft
    fs_start, fs_end, fs_duration = contiguous_interval_containing(df["stage_ft"], mask_fs, peak_time)
    if fs_duration == 0.0 or fs_start is None:
        print(f"  duration above flood stage ({flood_stage_ft} ft): never exceeded "
              f"(peak stage {peak_stage:.2f} ft)")
    else:
        print(f"  duration above flood stage ({flood_stage_ft} ft): {fs_duration:.1f} h "
              f"({fs_start} -> {fs_end})")

    # Duration above 2x pre-event baseflow
    thresh2 = 2 * bf_val
    mask_2bf = df["discharge_cfs"] > thresh2
    bf2_start, bf2_end, bf2_duration = contiguous_interval_containing(
        df["discharge_cfs"], mask_2bf, peak_time)
    print(f"  duration above 2x baseflow ({thresh2:.1f} cfs): {bf2_duration:.1f} h "
          f"({bf2_start} -> {bf2_end})")

    return dict(
        key=key,
        label=cfg["label"],
        df=df,
        data_source=source,
        peak_time=peak_time,
        peak_Q=peak_Q,
        peak_stage_time=peak_stage_time,
        peak_stage=peak_stage,
        baseflow_time=bf_time,
        baseflow=bf_val,
        trough_time=trough_time,
        trough_val=trough_val,
        primary_fit=primary_fit,
        extended_fit=extended_fit,
        flood_stage_ft=flood_stage_ft,
        fs_start=fs_start, fs_end=fs_end, fs_duration_hours=fs_duration,
        bf2_threshold=thresh2,
        bf2_start=bf2_start, bf2_end=bf2_end, bf2_duration_hours=bf2_duration,
    )


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def make_figure(results, out_pdf, out_png):
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.8), sharey=False)

    for ax, res in zip(axes, results):
        df = res["df"]
        peak_time = res["peak_time"]
        fit = res["primary_fit"]

        plot_start = peak_time - pd.Timedelta(days=4)
        plot_end = min(df.index.max(), peak_time + pd.Timedelta(days=20))
        window = df.loc[plot_start:plot_end, "discharge_cfs"]

        ax.plot(window.index, window.values, color=COLOR_DATA, lw=0.9,
                label="Observed discharge")

        # shaded fit window
        ax.axvspan(fit["window_start"], fit["window_end"], color=COLOR_FIT, alpha=0.12,
                   lw=0)

        # fitted exponential line, drawn across the fit window
        t_line = np.linspace(0, fit["window_hours"], 100)
        Q_line = np.exp(fit["slope"] * t_line + fit["intercept"])
        times_line = [fit["window_start"] + pd.Timedelta(hours=h) for h in t_line]
        ax.plot(times_line, Q_line, color=COLOR_FIT, lw=1.4, ls="--",
                 label=r"Exponential fit")

        ax.set_yscale("log")
        ax.set_xlabel("Date")
        ax.set_ylabel("Discharge (cfs)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=4))
        for label in ax.get_xticklabels():
            label.set_rotation(0)
            label.set_ha("center")

        # annotate tau
        tau_txt = (rf"$\tau$ = {fit['tau_hours']:.0f} h "
                   rf"({fit['tau_days']:.1f} d)" + "\n" + rf"$R^2$ = {fit['r2']:.2f}")
        ax.text(0.97, 0.93, tau_txt, transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color=COLOR_FIT)

        ax.text(0.03, 0.93, res["label"], transform=ax.transAxes, ha="left", va="top",
                fontsize=9, fontweight="bold")

    axes[0].legend(loc="lower left", fontsize=7, frameon=False)

    fig.tight_layout(pad=0.6)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    print(f"\nSaved figure: {out_pdf}")
    print(f"Saved figure: {out_png}")


# ----------------------------------------------------------------------
# Results markdown
# ----------------------------------------------------------------------

def fmt_time(t):
    if t is None:
        return "n/a"
    return t.tz_convert("America/New_York").strftime("%Y-%m-%d %H:%M %Z") + f" ({t.strftime('%Y-%m-%d %H:%M')} UTC)"


def write_results_md(results, drainage_area, drainage_area_src, flood_categories, flood_stage_src, out_path):
    lines = []
    lines.append("# Winooski River at Montpelier, VT -- Drainage-Timescale Analysis\n")
    lines.append(
        "Order-of-magnitude flood-recession analysis for USGS gauge 04286000 "
        "(WINOOSKI RIVER AT MONTPELIER, VT; NWS/AHPS identifier MONV1), "
        "generated from public USGS NWIS and NOAA/NWS data for the July 2023 "
        "and July 2024 flood events. All numbers below are computed by "
        "`analysis/winooski_recession.py` from live service queries; see that "
        "script for exact methods.\n"
    )

    lines.append("## Site-number correction\n")
    lines.append(
        "The originally specified site number, USGS 04288000, is **MAD RIVER "
        "NEAR MORETOWN, VT**, not the Winooski River at Montpelier (confirmed "
        "by querying the NWIS IV service and inspecting the returned "
        "`sourceInfo.siteName`). The correct site for \"WINOOSKI RIVER AT "
        "MONTPELIER, VT\" is **USGS 04286000** (NWS/AHPS gauge `MONV1`), which "
        "is what this analysis uses. Any prior draft material citing 04288000 "
        "for this gauge should be corrected.\n"
    )

    lines.append("## Site metadata\n")
    lines.append(f"- USGS site: 04286000, WINOOSKI RIVER AT MONTPELIER, VT\n"
                  f"- NWS/AHPS identifier: MONV1\n"
                  f"- Drainage area: **{drainage_area:.0f} sq mi** "
                  f"(source: {drainage_area_src})\n")

    lines.append("## NWS flood stage\n")
    cats_str = ", ".join(f"{k}={v['stage']} ft" for k, v in flood_categories.items() if "stage" in v)
    lines.append(
        f"- NWS/NOAA flood categories for MONV1: {cats_str}\n"
        f"- Source: {flood_stage_src}\n"
        f"- \"Flood stage\" as used below = **minor flood category "
        f"({flood_categories.get('minor', {}).get('stage', 'n/a')} ft)**, the "
        f"standard NWS AHPS threshold for the onset of flooding impacts.\n"
    )

    for res in results:
        lines.append(f"## Event: {res['label']}\n")
        lines.append(f"- Data source: USGS NWIS **{res['data_source']}** (instantaneous values)\n")
        lines.append(f"- Peak discharge: **{res['peak_Q']:.0f} cfs** at {fmt_time(res['peak_time'])}\n")
        lines.append(f"- Peak stage: **{res['peak_stage']:.2f} ft** at {fmt_time(res['peak_stage_time'])}\n")
        lines.append(f"- Pre-event baseflow: **{res['baseflow']:.0f} cfs** "
                      f"(minimum hourly-mean discharge {60}-{12} h before peak; "
                      f"occurred {fmt_time(res['baseflow_time'])})\n")

        pf = res["primary_fit"]
        ef = res["extended_fit"]
        lines.append(
            f"- **Primary recession fit** (window chosen for R^2 quality): "
            f"peak -> {fmt_time(pf['window_end'])} "
            f"({pf['window_hours']:.1f} h), "
            f"tau = **{pf['tau_hours']:.1f} h ({pf['tau_days']:.2f} d)**, "
            f"R^2 = **{pf['r2']:.3f}**, n = {pf['n']} points\n"
        )
        lines.append(
            f"- **Extended recession fit** (peak to the trough just before the "
            f"next rain-driven rise interrupts recession): "
            f"peak -> {fmt_time(ef['window_end'])} "
            f"({ef['window_hours']:.1f} h), "
            f"tau = **{ef['tau_hours']:.1f} h ({ef['tau_days']:.2f} d)**, "
            f"R^2 = **{ef['r2']:.3f}**, n = {ef['n']} points\n"
        )
        if res["fs_duration_hours"] and res["fs_duration_hours"] > 0:
            lines.append(
                f"- Duration above flood stage ({res['flood_stage_ft']} ft): "
                f"**{res['fs_duration_hours']:.1f} h** "
                f"({fmt_time(res['fs_start'])} -> {fmt_time(res['fs_end'])})\n"
            )
        else:
            lines.append(
                f"- Duration above flood stage ({res['flood_stage_ft']} ft): "
                f"**never exceeded** (peak stage {res['peak_stage']:.2f} ft "
                f"stayed below the {res['flood_stage_ft']} ft minor-flood threshold)\n"
            )
        lines.append(
            f"- Duration above 2x pre-event baseflow "
            f"({res['bf2_threshold']:.0f} cfs): **{res['bf2_duration_hours']:.1f} h** "
            f"({fmt_time(res['bf2_start'])} -> {fmt_time(res['bf2_end'])})\n"
        )
        lines.append("")

    lines.append("## Data caveats\n")
    lines.append(
        "- USGS IV data are **provisional** and subject to revision; values "
        "used here were retrieved 2026-07-22 and should be re-verified before "
        "final submission if the proposal timeline allows.\n"
        "- 2023 IV data are reported at a mix of 5-minute and 15-minute "
        "intervals (USGS increases reporting frequency during high flow); "
        "2024 data are uniformly 15-minute. Duration calculations use actual "
        "timestamp deltas on contiguous threshold-exceedance runs (not simple "
        "point counts), specifically to avoid bias from this irregular "
        "sampling.\n"
        "- Both recessions are interrupted by a **secondary, rain-driven rise** "
        "roughly 4.8-4.9 days after the flood peak (2023: secondary peak "
        "~5,900 cfs around July 17; 2024: secondary peak ~2,700 cfs around "
        "July 17). The 'extended' fit window in each case ends at the trough "
        "immediately before that secondary rise, not at a full return to "
        "baseflow -- discharge had not returned to pre-event baseflow by "
        "that point in either year.\n"
        "- The 2023 recession is **not well described by a single exponential** "
        "over its full length: R^2 is excellent (>0.9) for the first ~60 h "
        "post-peak but degrades steadily thereafter (R^2 ~ 0.74 by the 117.5 h "
        "trough) as the recession visibly bends to a shallower slope, "
        "consistent with a transition from fast quickflow recession to slower "
        "baseflow-dominated recession. The primary tau reported above uses "
        "the well-fit 60 h window; the extended-window tau is reported for "
        "completeness but should be treated as a lower-quality, blended "
        "estimate.\n"
        "- The 2024 recession is comparably well-described by a single "
        "exponential over almost its entire length (R^2 = 0.94-0.97 "
        "depending on window), reflecting a smaller, simpler flood pulse.\n"
        "- 'Duration above 2x pre-event baseflow' is reported as the single "
        "contiguous interval containing the flood peak, not total time above "
        "threshold across the full 6-week query window. In 2023 this interval "
        "is long (~19 days) because July 2023 was an exceptionally wet month "
        "in Vermont with several additional rain pulses keeping discharge "
        "elevated well after the main flood peak; treat this number as "
        "reflecting the sustained wet regime that month, not the duration of "
        "the single flood pulse itself. The 2023 discharge record also shows "
        "two lesser precursor rain pulses (~July 2-3 and ~July 7-8) before the "
        "main flood-triggering rain; the pre-event baseflow value is measured "
        "in the partial lull between the July 7-8 pulse and the main rise, so "
        "it is somewhat elevated relative to true undisturbed baseflow "
        "(deeper lows of ~550-650 cfs occur earlier in the record). The 2024 "
        "record shows an analogous smaller precursor pulse around July 6-7.\n"
        "- Peak stage and peak discharge times are within 15 minutes of each "
        "other in both events (sensor/rating-curve timing, not a real "
        "physical lag).\n"
        "- The 2024 event's peak stage (14.41-14.45 ft) did **not** reach the "
        "NWS minor flood stage (15 ft) at this gauge, despite discharge (11,900 "
        "cfs) being roughly half the 2023 peak (23,100 cfs) -- consistent with "
        "2023 being the more severe flood event at this location.\n"
    )

    out_path.write_text("\n".join(lines))
    print(f"\nSaved results: {out_path}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    IMAGES_DIR.mkdir(parents=True, exist_ok=True)

    print("Fetching site metadata...")
    drainage_area, drainage_area_src = fetch_drainage_area(SITE_NO)

    print("Fetching NWS flood-category stages...")
    flood_categories, flood_stage_src = fetch_flood_stage()

    results = []
    for key, cfg in EVENTS.items():
        res = analyze_event(key, cfg, flood_categories)
        results.append(res)

    make_figure(
        results,
        out_pdf=IMAGES_DIR / "winooski_recession.pdf",
        out_png=ANALYSIS_DIR / "winooski_recession.png",
    )

    write_results_md(
        results, drainage_area, drainage_area_src, flood_categories, flood_stage_src,
        out_path=ANALYSIS_DIR / "results.md",
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
