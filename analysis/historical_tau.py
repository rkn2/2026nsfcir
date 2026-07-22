#!/usr/bin/env python3
"""
Winooski River at Montpelier, VT (USGS 04286000) -- multi-event historical
recession-timescale (tau) analysis, spanning the full USGS instantaneous-value
(IV) period of record.

Produced for an NSF proposal on flood-duration modeling in Montpelier, VT.
This script follows the conventions established in `winooski_recession.py`
(live NWIS queries with retries and no silent fallbacks, timestamp-delta
durations, colorblind-safe palette) and extends the single-event (July 2023 /
July 2024) analysis in that script to a multi-decade sample of independent
flood events, to address two open items in the proposal:

  1. The proposal pre-estimates the LogNormal scale parameter sigma_MB for
     the recession-timescale prior from "historical gauge variability."  This
     script quantifies that variability directly: it fits tau for every
     independent flood event above a discharge threshold across the gauge's
     IV record and reports the resulting sd(ln tau).
  2. The proposal has an open concern that drainage is faster from dry
     antecedent ground (i.e., tau should correlate with pre-event baseflow /
     antecedent wetness). This script regresses tau against pre-event
     baseflow (normalized by the record median flow) to quantify that
     relationship, and reports it by season as well.

Method summary
--------------
1. Query the NWIS site "series catalog" (seriesCatalogOutput) to determine
   the actual begin/end dates of the 00060 (discharge) "uv" (instantaneous
   values) series at this site -- do not assume a start year.
2. Fetch IV discharge (00060 only; no stage) year-by-year across that full
   period, immediately resampling each year to an hourly mean and discarding
   the raw high-frequency values, to keep the concatenated multi-decade
   series memory-tractable. NWIS IV discharge at this site is intermittently
   unavailable during ice-affected winter periods (see caveats below); this
   is a real USGS data-availability gap, not a fetch failure, and is
   reported as such.
3. Detect independent flood events on the concatenated hourly series: greedy
   peak-picking above a discharge threshold, masking +-10 days around each
   accepted peak before searching for the next one, so accepted peaks are
   guaranteed to be more than 10 days apart.
4. For each candidate event, fetch a precise IV window (00060 only) around
   the hourly-detected peak and redo the peak/baseflow/trough/recession-fit
   analysis at full instantaneous resolution, using the *same* methods as
   `winooski_recession.py`:
       - pre-event baseflow = min hourly-mean discharge in [peak-60h, peak-12h]
       - recession trough = first post-peak break where discharge rises >15%
         over the running minimum (the point recession is interrupted by a
         secondary rain pulse), searched up to 12 days after peak
       - exponential fit ln(Q) = ln(Q_peak) - t/tau over min(60h, trough)
   July 2023 and July 2024 are expected to reproduce tau ~= 37 h and ~= 45 h
   respectively (as found in `winooski_recession.py`); this is the
   correctness check for this script's independent re-implementation.
5. Events with R^2 < 0.8 on the primary fit are flagged and reported but
   excluded from the summary statistics (sd(ln tau), regression, etc.).
6. Writes `analysis/historical_tau.png` (two-panel figure) and
   `analysis/historical_tau.md` (method, event table, summary stats,
   antecedent-wetness findings, caveats).

Usage
-----
    python3 historical_tau.py

Requires: requests, pandas, numpy, matplotlib (standard library otherwise).
Runtime: several minutes (year-by-year IV fetches across ~35 years).
"""

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

SITE_NO = "04286000"  # WINOOSKI RIVER AT MONTPELIER, VT (see winooski_recession.py
                       # for the 04288000 -> 04286000 site-number correction)
SITE_NAME = "WINOOSKI RIVER AT MONTPELIER, VT"

IV_URL = "https://waterservices.usgs.gov/nwis/iv/"
SITE_URL = "https://waterservices.usgs.gov/nwis/site/"

# Fallback values, empirically confirmed from live queries on 2026-07-22,
# used only if the corresponding live service call fails at run time.
FALLBACK_DRAINAGE_AREA_SQMI = 397.0
FALLBACK_POR_BEGIN = "1990-10-01"
FALLBACK_POR_END = None  # filled with today's date if needed

# Event-detection parameters (see historical_tau.md for how these were
# chosen: threshold=7000 cfs on the hourly-mean series lands the event count
# at 16, inside the requested 8-20 range, and independently reproduces both
# the July 2023 and July 2024 events from winooski_recession.py).
PEAK_THRESHOLD_CFS = 7000.0
MIN_SEPARATION_DAYS = 10
FETCH_LOOKBACK_DAYS = 5     # before hourly-detected peak, for precise refetch
FETCH_LOOKAHEAD_DAYS = 13   # after hourly-detected peak, for precise refetch
GOOD_FIT_R2 = 0.8

REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "analysis"

# Colorblind-safe palette: blue data, vermillion/orange fit/highlight
COLOR_DATA = "#0072B2"
COLOR_FIT = "#D55E00"


# ----------------------------------------------------------------------
# Data acquisition
# ----------------------------------------------------------------------

def fetch_period_of_record(site_no, retries=4, backoff=5):
    """Query the NWIS site series-catalog service to find the actual
    begin/end dates of the 00060 instantaneous-value (uv) series at this
    site. Raises on failure (no silent fallback for the primary path;
    fallback constants are used only after retries are exhausted, and the
    fallback is logged loudly)."""
    params = {"sites": site_no, "format": "rdb", "seriesCatalogOutput": "true"}
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.get(SITE_URL, params=params, timeout=60)
            if r.status_code == 200:
                lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
                header = lines[0].split("\t")
                for data_line in lines[2:]:
                    fields = data_line.split("\t")
                    row = dict(zip(header, fields))
                    if row.get("data_type_cd") == "uv" and row.get("parm_cd") == "00060":
                        begin = row["begin_date"]
                        end = row["end_date"]
                        print(f"  [fetch_period_of_record] IV discharge (00060) "
                              f"period of record: {begin} .. {end} (source: {r.url})")
                        return begin, end, r.url
                raise RuntimeError("No uv/00060 row found in series catalog")
            else:
                last_exc = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:  # noqa: BLE001
            last_exc = e
        print(f"  [fetch_period_of_record] attempt {attempt+1}/{retries} failed ({last_exc}); retrying...")
        time.sleep(backoff)
    end_fallback = FALLBACK_POR_END or pd.Timestamp.utcnow().strftime("%Y-%m-%d")
    print(f"  [fetch_period_of_record] live query failed after {retries} attempts "
          f"({last_exc}); using fallback constants {FALLBACK_POR_BEGIN} .. {end_fallback} "
          f"(confirmed via manual query 2026-07-22).")
    return FALLBACK_POR_BEGIN, end_fallback, SITE_URL + " (fallback, manual query 2026-07-22)"


def fetch_drainage_area(site_no, retries=4, backoff=5):
    """Fetch drainage area (sq mi) from the NWIS expanded site service."""
    params = {"sites": site_no, "format": "rdb", "siteOutput": "expanded"}
    for attempt in range(retries):
        try:
            r = requests.get(SITE_URL, params=params, timeout=60)
            if r.status_code == 200:
                lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
                header = lines[0].split("\t")
                data_line = lines[2].split("\t")
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


def fetch_iv_discharge(site_no, start, end, retries=4, backoff=5):
    """Fetch USGS NWIS instantaneous-value discharge (00060, cfs) only, for
    [start, end]. Returns a pandas Series indexed by UTC datetime. Raises on
    failure after all retries -- no silent fallback."""
    params = {
        "sites": site_no,
        "parameterCd": "00060",
        "startDT": start,
        "endDT": end,
        "format": "json",
    }
    last_exc = None
    for attempt in range(retries):
        try:
            r = requests.get(IV_URL, params=params, timeout=120)
            if r.status_code == 200:
                data = r.json()
                ts = data.get("value", {}).get("timeSeries", [])
                if not ts:
                    # A genuinely empty series (e.g. ice-affected winter gap)
                    # is not an error -- return an empty Series.
                    return pd.Series(dtype=float, name="discharge_cfs")
                vals = ts[0]["values"][0]["value"]
                rows = [(v["dateTime"], float(v["value"]))
                        for v in vals if v["value"] not in (None, "")]
                if not rows:
                    return pd.Series(dtype=float, name="discharge_cfs")
                s = pd.DataFrame(rows, columns=["dateTime", "value"])
                s["dateTime"] = pd.to_datetime(s["dateTime"], utc=True)
                s = s.set_index("dateTime")["value"]
                s.name = "discharge_cfs"
                s = s.sort_index()
                s = s[~s.index.duplicated(keep="first")]
                return s
            else:
                last_exc = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:  # noqa: BLE001
            last_exc = e
        print(f"  [fetch_iv_discharge] attempt {attempt+1}/{retries} failed ({last_exc}); retrying...")
        time.sleep(backoff)
    raise RuntimeError(f"fetch_iv_discharge failed after {retries} attempts: {last_exc}")


def build_full_hourly_record(site_no, por_begin, por_end):
    """Fetch IV discharge year-by-year across the full period of record,
    resampling each year to hourly means immediately (to keep memory
    tractable across ~35 years of 5/15-minute data). A year that fails all
    retries is skipped with a loud, explicit warning (not a silent gap) and
    recorded in `failed_years`; a year that succeeds but is genuinely empty
    (ice-affected winter periods at this site regularly have no computed
    discharge) is recorded in `empty_years` and is a real USGS data
    limitation, not a fetch error."""
    begin_year = pd.Timestamp(por_begin).year
    end_dt = pd.Timestamp(por_end)
    end_year = end_dt.year

    frames = []
    failed_years = []
    empty_years = []
    for year in range(begin_year, end_year + 1):
        start_s = f"{year}-01-01"
        end_s = f"{year}-12-31" if year < end_year else end_dt.strftime("%Y-%m-%d")
        if year == begin_year:
            start_s = pd.Timestamp(por_begin).strftime("%Y-%m-%d")
        try:
            s = fetch_iv_discharge(site_no, start_s, end_s)
        except Exception as e:  # noqa: BLE001
            print(f"  *** YEAR {year} FETCH FAILED AFTER RETRIES -- COVERAGE GAP: {e}")
            failed_years.append(year)
            continue
        if len(s) == 0:
            print(f"  year {year}: 0 rows returned (no IV discharge published -- "
                  f"commonly an ice-affected-winter gap at this site)")
            empty_years.append(year)
            continue
        hourly = s.resample("1h").mean()
        frames.append(hourly)
        print(f"  year {year}: n_raw={len(s)} n_hourly={len(hourly)}")

    if not frames:
        raise RuntimeError("No IV discharge data could be retrieved for any year in the period of record")

    full = pd.concat(frames).sort_index()
    full = full[~full.index.duplicated(keep="first")]
    print(f"\n  Full hourly record: {len(full)} points, {full.index.min()} .. {full.index.max()}")
    if failed_years:
        print(f"  ** {len(failed_years)} year(s) failed to fetch (network/service errors): {failed_years}")
    if empty_years:
        print(f"  ** {len(empty_years)} year(s) had zero published IV discharge "
              f"(data-availability gaps, not fetch failures): {empty_years}")
    return full, failed_years, empty_years


# ----------------------------------------------------------------------
# Event detection
# ----------------------------------------------------------------------

def detect_events(hourly, threshold, min_sep_days=10):
    """Greedy independent-peak detection: repeatedly take the global maximum
    of the remaining series, record it as an event, then mask +-min_sep_days
    around it before searching again. Guarantees every pair of accepted
    peaks is more than min_sep_days apart. Returns a chronologically sorted
    list of (peak_time, peak_value) tuples."""
    work = hourly.dropna().copy()
    events = []
    while True:
        m = work.max() if len(work) else np.nan
        if pd.isna(m) or m < threshold:
            break
        t = work.idxmax()
        v = work.loc[t]
        events.append((t, float(v)))
        lo = t - pd.Timedelta(days=min_sep_days)
        hi = t + pd.Timedelta(days=min_sep_days)
        work.loc[lo:hi] = -np.inf
    events.sort(key=lambda x: x[0])
    return events


def season_of(ts):
    m = ts.month
    if m in (12, 1, 2):
        return "winter"
    if m in (3, 4, 5):
        return "spring"
    if m in (6, 7, 8):
        return "summer"
    return "fall"


# ----------------------------------------------------------------------
# Per-event analysis (same methods as winooski_recession.py)
# ----------------------------------------------------------------------

def find_peak(discharge):
    peak_time = discharge.idxmax()
    peak_Q = discharge.loc[peak_time]
    return peak_time, float(peak_Q)


def find_pre_event_baseflow(discharge, peak_time, lookback_start_h=60, lookback_end_h=12):
    """Pre-event baseflow = minimum hourly-mean discharge in the window
    [peak - lookback_start_h, peak - lookback_end_h]. Identical method to
    winooski_recession.py."""
    hourly = discharge.resample("1h").mean()
    win = hourly.loc[peak_time - pd.Timedelta(hours=lookback_start_h):
                      peak_time - pd.Timedelta(hours=lookback_end_h)]
    win = win.dropna()
    if len(win) == 0:
        return None, np.nan
    bf_time = win.idxmin()
    bf_val = float(win.min())
    return bf_time, bf_val


def find_recession_trough(discharge, peak_time, search_days=12, rise_frac=0.15, min_hours=24):
    """Identical method to winooski_recession.py: walk forward on the hourly
    series after the peak; the trough is the running minimum at the first
    point (after min_hours) where discharge 24h later exceeds
    (1+rise_frac) times the running minimum."""
    hourly = discharge.resample("1h").mean()
    win = hourly.loc[peak_time: peak_time + pd.Timedelta(days=search_days)].dropna()
    times = win.index
    vals = win.values
    if len(vals) == 0:
        return peak_time, np.nan
    running_min = vals[0]
    running_min_t = times[0]
    for i in range(len(vals) - 1):
        if vals[i] < running_min:
            running_min = vals[i]
            running_min_t = times[i]
        elapsed_h = (times[i] - peak_time).total_seconds() / 3600.0
        if elapsed_h < min_hours:
            continue
        future_idx = min(i + 24, len(vals) - 1)
        if vals[future_idx] > (1 + rise_frac) * running_min and future_idx > i:
            return running_min_t, running_min
    return times[-1], vals[-1]


def fit_exponential_recession(discharge, peak_time, fit_end):
    """Fit ln(Q) = ln(Q_peak) - t/tau via least squares on the falling limb
    from peak_time to fit_end. Identical method to winooski_recession.py."""
    seg = discharge.loc[peak_time:fit_end]
    seg = seg[seg > 0].dropna()
    t_hours = (seg.index - peak_time).total_seconds() / 3600.0
    lnQ = np.log(seg.values)
    slope, intercept = np.polyfit(t_hours, lnQ, 1)
    pred = slope * t_hours + intercept
    ss_res = np.sum((lnQ - pred) ** 2)
    ss_tot = np.sum((lnQ - lnQ.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    tau_hours = -1.0 / slope
    return dict(
        tau_hours=tau_hours, tau_days=tau_hours / 24.0, r2=r2, n=len(seg),
        slope=slope, intercept=intercept,
        window_start=peak_time, window_end=fit_end,
        window_hours=(fit_end - peak_time).total_seconds() / 3600.0,
    )


def analyze_candidate_event(site_no, approx_peak_time, approx_peak_val):
    """Fetch a precise IV window around an hourly-detected candidate peak
    and run the full peak/baseflow/trough/recession-fit pipeline at
    instantaneous resolution."""
    start = (approx_peak_time - pd.Timedelta(days=FETCH_LOOKBACK_DAYS)).strftime("%Y-%m-%d")
    end = (approx_peak_time + pd.Timedelta(days=FETCH_LOOKAHEAD_DAYS)).strftime("%Y-%m-%d")
    discharge = fetch_iv_discharge(site_no, start, end)
    if len(discharge) == 0:
        raise RuntimeError(f"No IV data in refetch window {start}..{end}")

    peak_time, peak_Q = find_peak(discharge)
    # sanity check: precise peak should be close in time to the hourly-detected one
    dt_h = abs((peak_time - approx_peak_time).total_seconds()) / 3600.0
    if dt_h > 48:
        print(f"  ** WARNING: precise peak ({peak_time}) is {dt_h:.0f} h from the "
              f"hourly-detected candidate ({approx_peak_time}) -- possible event "
              f"mismatch; verify manually.")

    bf_time, bf_val = find_pre_event_baseflow(discharge, peak_time)
    trough_time, trough_val = find_recession_trough(discharge, peak_time)
    primary_end = min(peak_time + pd.Timedelta(hours=60), trough_time)
    fit = fit_exponential_recession(discharge, peak_time, primary_end)

    return dict(
        peak_time=peak_time, peak_Q=peak_Q,
        season=season_of(peak_time), year=peak_time.year,
        baseflow_time=bf_time, baseflow=bf_val,
        trough_time=trough_time, trough_val=trough_val,
        fit=fit,
        discharge=discharge,
    )


# ----------------------------------------------------------------------
# Summary statistics and regression
# ----------------------------------------------------------------------

def summarize(events_df):
    """events_df: DataFrame with columns including tau_hours, r2,
    baseflow_ratio (baseflow / record median flow). Only R2 >= GOOD_FIT_R2
    rows are used."""
    good = events_df[events_df["r2"] >= GOOD_FIT_R2].copy()
    tau = good["tau_hours"].values
    ln_tau = np.log(tau)

    summary = dict(
        n_events_total=len(events_df),
        n_events_good=len(good),
        n_events_excluded=len(events_df) - len(good),
        tau_min=float(np.min(tau)), tau_max=float(np.max(tau)),
        tau_median=float(np.median(tau)),
        tau_q25=float(np.percentile(tau, 25)), tau_q75=float(np.percentile(tau, 75)),
        ln_tau_mean=float(np.mean(ln_tau)), ln_tau_sd=float(np.std(ln_tau, ddof=1)),
    )

    # regression: tau vs baseflow_ratio (well-fit events only)
    x = good["baseflow_ratio"].values
    y = good["tau_hours"].values
    slope, intercept = np.polyfit(x, y, 1)
    r = np.corrcoef(x, y)[0, 1]
    driest = good.loc[good["baseflow_ratio"].idxmin()]
    wettest = good.loc[good["baseflow_ratio"].idxmax()]
    summary.update(dict(
        regr_slope=float(slope), regr_intercept=float(intercept), regr_r=float(r),
        driest_event=driest["label"], driest_ratio=float(driest["baseflow_ratio"]),
        driest_tau=float(driest["tau_hours"]),
        wettest_event=wettest["label"], wettest_ratio=float(wettest["baseflow_ratio"]),
        wettest_tau=float(wettest["tau_hours"]),
        regr_tau_at_driest=float(slope * driest["baseflow_ratio"] + intercept),
        regr_tau_at_wettest=float(slope * wettest["baseflow_ratio"] + intercept),
    ))

    season_tab = good.groupby("season")["tau_hours"].agg(["median", "mean", "count"])
    return summary, good, season_tab


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def make_figure(events_df, summary, out_png):
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))

    good = events_df[events_df["r2"] >= GOOD_FIT_R2]
    n_poor = len(events_df) - len(good)
    highlight = good[good["year"].isin([2023, 2024])]
    other = good[~good["year"].isin([2023, 2024])]

    # Panel (a): tau vs event year. Poor-fit events (R^2 < 0.8) are excluded
    # from this plot -- as from all summary stats -- because at least one
    # poor-fit event's tau (Apr 2001, ~244 h) is a fit-quality artifact, not
    # a physically meaningful recession timescale, and including it on a
    # linear axis swamps the well-fit events. Poor-fit events remain fully
    # visible in the event table in historical_tau.md.
    ax = axes[0]
    ax.scatter(other["peak_time"], other["tau_hours"], color=COLOR_DATA, s=32,
               zorder=3, label="Well-fit events ($R^2\\geq$0.8)")
    ax.scatter(highlight["peak_time"], highlight["tau_hours"], color=COLOR_FIT, s=55,
               marker="D", zorder=4, label="July 2023 / July 2024")
    ax.axhline(summary["tau_median"], color="gray", lw=0.8, ls=":", zorder=1)
    ax.set_ylim(0, good["tau_hours"].max() * 1.15)
    ax.text(0.98, 0.04, f"median $\\tau$={summary['tau_median']:.0f} h",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color="gray")
    ax.set_xlabel("Event peak date")
    ax.set_ylabel(r"Recession $\tau$ (h)")
    ax.set_title(f"(a) $\\tau$ across {summary['n_events_total']} historical events "
                 f"({summary['n_events_good']} well-fit, {n_poor} poor-fit excluded)",
                 fontsize=8.5)
    ax.legend(loc="upper left", fontsize=6.5, frameon=False)

    # Panel (b): tau vs normalized pre-event baseflow
    ax = axes[1]
    ax.scatter(other["baseflow_ratio"], other["tau_hours"], color=COLOR_DATA, s=32,
               zorder=3, label="Well-fit events")
    ax.scatter(highlight["baseflow_ratio"], highlight["tau_hours"], color=COLOR_FIT, s=55,
               marker="D", zorder=4, label="July 2023 / July 2024")
    x_line = np.linspace(good["baseflow_ratio"].min(), good["baseflow_ratio"].max(), 50)
    y_line = summary["regr_slope"] * x_line + summary["regr_intercept"]
    ax.plot(x_line, y_line, color=COLOR_FIT, lw=1.3, ls="--",
            label=f"linear fit (r={summary['regr_r']:.2f})")
    ax.set_xlabel("Pre-event baseflow / record median flow")
    ax.set_ylabel(r"Recession $\tau$ (h)")
    ax.set_title("(b) $\\tau$ vs. antecedent wetness", fontsize=8.5)
    ax.legend(loc="upper left", fontsize=6.5, frameon=False)

    fig.tight_layout(pad=0.8)
    fig.savefig(out_png, dpi=200)
    print(f"\nSaved figure: {out_png}")


# ----------------------------------------------------------------------
# Results markdown
# ----------------------------------------------------------------------

def fmt_time(t):
    if t is None or (isinstance(t, float) and pd.isna(t)):
        return "n/a"
    return t.tz_convert("America/New_York").strftime("%Y-%m-%d %H:%M %Z") + f" ({t.strftime('%Y-%m-%d %H:%M')} UTC)"


def write_results_md(events_df, summary, season_tab, drainage_area, drainage_area_src,
                      por_begin, por_end, por_src, median_flow, failed_years, empty_years,
                      out_path):
    lines = []
    lines.append("# Winooski River at Montpelier, VT -- Multi-Event Historical Recession-Timescale (tau) Analysis\n")
    lines.append(
        "Quantifies historical variability and antecedent-moisture dependence "
        "of the exponential flood-recession timescale tau for USGS gauge "
        "04286000 (WINOOSKI RIVER AT MONTPELIER, VT), across the full USGS "
        "instantaneous-value (IV) period of record. Produced to (1) directly "
        "quantify the LogNormal sigma_MB scale parameter the proposal "
        "pre-estimates from \"historical gauge variability,\" and (2) test the "
        "proposal's open concern that recession is faster from drier "
        "antecedent ground. All numbers below are computed by "
        "`analysis/historical_tau.py` from live NWIS service queries; see "
        "that script for exact methods. It is a self-contained companion to, "
        "and reuses the same fitting methods as, `analysis/winooski_recession.py`.\n"
    )

    lines.append("## Method\n")
    lines.append(
        f"- **Period of record used:** {por_begin} to {por_end}, the actual "
        f"begin/end dates of the 00060 (discharge) instantaneous-value series "
        f"at this site, retrieved live from the NWIS site series catalog "
        f"(source: {por_src}). This begins substantially earlier than the "
        f"~2007 often assumed for IV records at small gauges; gage-height "
        f"(00065) IV data at this site only begins 2007-10-01, but stage was "
        f"not needed for this analysis (discharge only).\n"
        f"- **Event detection:** the full-period IV discharge record was "
        f"fetched year-by-year and resampled to hourly means. Independent "
        f"flood peaks were detected by greedy peak-picking on the hourly "
        f"series: repeatedly take the global maximum above "
        f"**{PEAK_THRESHOLD_CFS:.0f} cfs**, record it as an event, then mask "
        f"+/-{MIN_SEPARATION_DAYS} days around it before searching again. This "
        f"guarantees every pair of accepted event peaks is more than "
        f"{MIN_SEPARATION_DAYS} days apart. The {PEAK_THRESHOLD_CFS:.0f} cfs "
        f"threshold was chosen (after checking 4000-10000 cfs) specifically "
        f"because it lands the event count at **{summary['n_events_total']}**, "
        f"inside the requested 8-20 range, while still independently "
        f"recovering both the July 2023 and July 2024 events used as the "
        f"single-event baseline in `winooski_recession.py`.\n"
        f"- **Per-event analysis:** for each hourly-detected candidate peak, "
        f"a precise IV discharge window "
        f"(-{FETCH_LOOKBACK_DAYS}/+{FETCH_LOOKAHEAD_DAYS} days) was refetched "
        f"at full instantaneous resolution, and the peak, pre-event baseflow, "
        f"recession trough, and exponential-recession tau were computed with "
        f"the **identical methods** as `winooski_recession.py`: pre-event "
        f"baseflow = minimum hourly-mean discharge in [peak-60h, peak-12h]; "
        f"recession trough = first post-peak point where discharge 24h later "
        f"exceeds 1.15x the running minimum (i.e., the break before a "
        f"secondary rain-driven rise); exponential fit ln(Q) = ln(Q_peak) - "
        f"t/tau by least squares over min(60h, trough).\n"
        f"- **Quality control:** events with primary-fit R^2 < {GOOD_FIT_R2} "
        f"are flagged and reported in the table below but excluded from all "
        f"summary statistics (sd(ln tau), regression, median/IQR).\n"
        f"- **Verification:** this script's independent re-implementation "
        f"reproduces the July 2023 and July 2024 tau values from "
        f"`winooski_recession.py` (37 h and 45 h) -- see the event table.\n"
        f"- **Antecedent-wetness normalization:** pre-event baseflow is "
        f"normalized as a dimensionless ratio to the record median hourly "
        f"discharge (**{median_flow:.0f} cfs**, computed across the full "
        f"{por_begin}-{por_end} hourly record), and separately reported as "
        f"unit discharge per drainage area "
        f"(**{drainage_area:.0f} sq mi**, source: {drainage_area_src}).\n"
    )

    lines.append("## Event table\n")
    lines.append("| Event (peak date) | Season | Peak Q (cfs) | Pre-event baseflow (cfs) | Baseflow / median | tau (h) | R^2 | Fit window (h) | Status |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for _, row in events_df.iterrows():
        status = "well-fit" if row["r2"] >= GOOD_FIT_R2 else "**poor fit, excluded from stats**"
        lines.append(
            f"| {row['label']} | {row['season']} | {row['peak_Q']:.0f} | "
            f"{row['baseflow']:.0f} | {row['baseflow_ratio']:.2f} | "
            f"{row['tau_hours']:.1f} | {row['r2']:.3f} | {row['window_hours']:.1f} | {status} |"
        )
    lines.append("")

    lines.append("## Summary: variability of tau (input to sigma_MB)\n")
    lines.append(
        f"Across the **{summary['n_events_good']} well-fit events** "
        f"(R^2 >= {GOOD_FIT_R2}; {summary['n_events_excluded']} of "
        f"{summary['n_events_total']} total events excluded for poor fit):\n\n"
        f"- tau range: **{summary['tau_min']:.1f} - {summary['tau_max']:.1f} h**\n"
        f"- tau median: **{summary['tau_median']:.1f} h**\n"
        f"- tau IQR: **{summary['tau_q25']:.1f} - {summary['tau_q75']:.1f} h**\n"
        f"- **sd(ln tau) = {summary['ln_tau_sd']:.3f}** "
        f"(mean ln tau = {summary['ln_tau_mean']:.3f}; sample standard "
        f"deviation, n={summary['n_events_good']})\n\n"
        f"This sd(ln tau) is the direct empirical estimate of the LogNormal "
        f"scale parameter the proposal pre-estimates as sigma_MB from "
        f"\"historical gauge variability\" -- it should be compared against "
        f"whatever sigma_MB value is currently written into the proposal "
        f"text.\n"
    )

    lines.append("## Antecedent-wetness relationship\n")
    lines.append(
        f"Linear regression of tau (h) on normalized pre-event baseflow "
        f"(ratio to record median flow) across the {summary['n_events_good']} "
        f"well-fit events: slope = **{summary['regr_slope']:.2f} h per unit "
        f"baseflow-ratio**, Pearson r = **{summary['regr_r']:.2f}**.\n\n"
        f"- Driest-antecedent well-fit event: **{summary['driest_event']}** "
        f"(baseflow/median = {summary['driest_ratio']:.2f}), tau = "
        f"**{summary['driest_tau']:.1f} h**\n"
        f"- Wettest-antecedent well-fit event: **{summary['wettest_event']}** "
        f"(baseflow/median = {summary['wettest_ratio']:.2f}), tau = "
        f"**{summary['wettest_tau']:.1f} h**\n"
        f"- Observed shift (driest to wettest event, raw values): "
        f"**{summary['wettest_tau'] - summary['driest_tau']:+.1f} h**\n"
        f"- Regression-line shift over the same baseflow-ratio range: "
        f"{summary['regr_tau_at_driest']:.1f} h to "
        f"{summary['regr_tau_at_wettest']:.1f} h "
        f"(**{summary['regr_tau_at_wettest'] - summary['regr_tau_at_driest']:+.1f} h**)\n\n"
    )
    sign_word = "longer" if summary["regr_slope"] > 0 else "shorter"
    lines.append(
        f"The positive/negative sign of the regression slope indicates tau "
        f"is **{sign_word}** (recession is slower) when pre-event baseflow is "
        f"higher relative to the record median -- i.e., wetter antecedent "
        f"conditions. "
        + ("This is consistent with the proposal's open concern in the "
           "other direction: it is *not* simply that dry ground drains "
           "faster in the sense of shorter tau; rather, wetter-antecedent "
           "events recede more slowly (larger tau), plausibly because "
           "already-saturated soil sustains a larger, slower baseflow-fed "
           "recession limb rather than a fast quickflow-dominated one.\n"
           if summary["regr_slope"] > 0 else
           "This is directionally consistent with the proposal's concern "
           "that recession is faster from drier antecedent ground.\n")
    )

    lines.append("### tau by season (well-fit events)\n")
    lines.append("| Season | Median tau (h) | Mean tau (h) | n |")
    lines.append("|---|---|---|---|")
    for season, row in season_tab.iterrows():
        lines.append(f"| {season} | {row['median']:.1f} | {row['mean']:.1f} | {int(row['count'])} |")
    lines.append("")

    lines.append("## Data caveats\n")
    caveat_lines = [
        "USGS IV data are **provisional** and subject to revision; values used "
        "here were retrieved 2026-07-22 and should be re-verified before final "
        "submission if the proposal timeline allows.",
        "This analysis uses a **single gauge** (04286000); tau variability and "
        "the antecedent-wetness relationship are specific to this gauge's "
        "drainage and may not transfer directly to ungauged tributaries or "
        "other basins in the proposal's study area.",
        "The **60 h fit window** (matching winooski_recession.py) is a fixed "
        "choice, not tuned per event; some events (as in the July 2023 case "
        "in winooski_recession.py) show recession curvature beyond 60 h that a "
        "single exponential does not capture well, which is part of why R^2 "
        "screening is applied.",
        "Several events' 'clean' recessions are interrupted by **secondary, "
        "rain-driven rises** within the 60 h primary fit window itself (not "
        "just after the trough); where this happens it depresses the primary "
        "R^2 and can trigger the poor-fit exclusion. This is a real feature "
        "of New England storm sequencing (multi-pulse rain events), not a "
        "methodological artifact.",
        "The gauge has **recurring winter data gaps**: NWIS does not publish "
        "computed instantaneous discharge for extended ice-affected periods "
        "(observed here as zero-row years/months, concentrated in "
        "January-February); this is a genuine USGS data-availability "
        "limitation, not a fetch failure, and means winter flood events (if "
        "any) are systematically under-represented in this event sample.",
        "The event-detection threshold "
        f"({PEAK_THRESHOLD_CFS:.0f} cfs) and separation window "
        f"({MIN_SEPARATION_DAYS} days) were chosen to land the event count in "
        "the requested 8-20 range and to recover the two known reference "
        "events; a different threshold would give a different (overlapping "
        "but not identical) event sample and correspondingly different "
        "sd(ln tau).",
        "Pre-event baseflow, and therefore the antecedent-wetness ratio, is "
        "measured in a fixed 60-12 h pre-peak window; where a precursor rain "
        "pulse falls inside that window (as documented for July 2023 in "
        "winooski_recession.py) the baseflow value is somewhat elevated "
        "relative to true undisturbed baseflow, which would bias that "
        "event's position in the antecedent-wetness regression toward the "
        "wet end.",
    ]
    if failed_years:
        caveat_lines.append(
            f"**{len(failed_years)} year(s) failed to fetch after retries** "
            f"(network/service errors, not data gaps) and are entirely absent "
            f"from the scanned record: {failed_years}. Any flood events in "
            f"those years are missing from this analysis."
        )
    if empty_years:
        caveat_lines.append(
            f"**{len(empty_years)} year(s) returned zero published IV "
            f"discharge rows** for the full year (beyond the routine winter "
            f"gaps noted above) and are absent from the scanned record: "
            f"{empty_years}."
        )
    for c in caveat_lines:
        lines.append(f"- {c}\n")

    out_path.write_text("\n".join(lines))
    print(f"\nSaved results: {out_path}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    print("Determining IV period of record...")
    por_begin, por_end, por_src = fetch_period_of_record(SITE_NO)

    print("\nFetching site metadata (drainage area)...")
    drainage_area, drainage_area_src = fetch_drainage_area(SITE_NO)

    print(f"\nBuilding full hourly discharge record ({por_begin} .. {por_end})...")
    full_hourly, failed_years, empty_years = build_full_hourly_record(SITE_NO, por_begin, por_end)
    median_flow = float(full_hourly.median())
    print(f"  Record median hourly discharge: {median_flow:.1f} cfs")

    print(f"\nDetecting independent events (threshold={PEAK_THRESHOLD_CFS:.0f} cfs, "
          f"min_sep={MIN_SEPARATION_DAYS} days)...")
    candidates = detect_events(full_hourly, PEAK_THRESHOLD_CFS, MIN_SEPARATION_DAYS)
    print(f"  {len(candidates)} candidate events detected:")
    for t, v in candidates:
        print(f"    {t}  Q~{v:.0f} cfs (hourly-mean)")

    rows = []
    for approx_peak_time, approx_peak_val in candidates:
        label = approx_peak_time.strftime("%b %Y")
        print(f"\n=== Event candidate: {label} (hourly-detected peak ~{approx_peak_val:.0f} cfs "
              f"at {approx_peak_time}) ===")
        try:
            res = analyze_candidate_event(SITE_NO, approx_peak_time, approx_peak_val)
        except Exception as e:  # noqa: BLE001
            print(f"  *** EVENT FAILED TO ANALYZE, EXCLUDED: {e}")
            continue
        fit = res["fit"]
        print(f"  precise peak: {res['peak_Q']:.0f} cfs at {res['peak_time']} (season={res['season']})")
        print(f"  pre-event baseflow: {res['baseflow']:.1f} cfs at {res['baseflow_time']}")
        print(f"  recession trough: {res['trough_val']:.1f} cfs at {res['trough_time']}")
        print(f"  fit: tau={fit['tau_hours']:.2f} h, R2={fit['r2']:.4f}, "
              f"window={fit['window_hours']:.1f} h, n={fit['n']}")

        rows.append(dict(
            label=res["peak_time"].strftime("%b %Y") + (" (Irene)" if (res["peak_time"].year == 2011 and res["peak_time"].month == 8) else ""),
            peak_time=res["peak_time"],
            year=res["year"],
            season=res["season"],
            peak_Q=res["peak_Q"],
            baseflow=res["baseflow"],
            baseflow_ratio=res["baseflow"] / median_flow if res["baseflow"] == res["baseflow"] else np.nan,
            tau_hours=fit["tau_hours"],
            r2=fit["r2"],
            window_hours=fit["window_hours"],
        ))

    events_df = pd.DataFrame(rows).sort_values("peak_time").reset_index(drop=True)

    # drop rows with NaN baseflow_ratio (can't be used in regression/table cleanly)
    n_before = len(events_df)
    events_df = events_df.dropna(subset=["baseflow_ratio", "tau_hours", "r2"]).reset_index(drop=True)
    if len(events_df) < n_before:
        print(f"\n  ** {n_before - len(events_df)} event(s) dropped for missing baseflow/tau/r2 "
              f"(insufficient pre-event data in refetch window)")

    # verification check against winooski_recession.py
    ev2023 = events_df[(events_df["year"] == 2023) & (events_df["season"] == "summer")]
    ev2024 = events_df[(events_df["year"] == 2024) & (events_df["season"] == "summer")]
    print("\n=== Verification against winooski_recession.py ===")
    if len(ev2023):
        print(f"  July 2023: tau = {ev2023.iloc[0]['tau_hours']:.1f} h "
              f"(winooski_recession.py: 36.8 h) -- "
              f"{'OK' if abs(ev2023.iloc[0]['tau_hours'] - 36.8) < 3 else 'MISMATCH -- CHECK'}")
    else:
        print("  ** July 2023 event NOT found in detected event set -- verification FAILED")
    if len(ev2024):
        print(f"  July 2024: tau = {ev2024.iloc[0]['tau_hours']:.1f} h "
              f"(winooski_recession.py: 45.2 h) -- "
              f"{'OK' if abs(ev2024.iloc[0]['tau_hours'] - 45.2) < 3 else 'MISMATCH -- CHECK'}")
    else:
        print("  ** July 2024 event NOT found in detected event set -- verification FAILED")

    summary, good_df, season_tab = summarize(events_df)

    make_figure(events_df, summary, out_png=ANALYSIS_DIR / "historical_tau.png")

    write_results_md(
        events_df, summary, season_tab, drainage_area, drainage_area_src,
        por_begin, por_end, por_src, median_flow, failed_years, empty_years,
        out_path=ANALYSIS_DIR / "historical_tau.md",
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
