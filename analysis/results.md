# Winooski River at Montpelier, VT -- Drainage-Timescale Analysis

Order-of-magnitude flood-recession analysis for USGS gauge 04286000 (WINOOSKI RIVER AT MONTPELIER, VT; NWS/AHPS identifier MONV1), generated from public USGS NWIS and NOAA/NWS data for the July 2023 and July 2024 flood events. All numbers below are computed by `analysis/winooski_recession.py` from live service queries; see that script for exact methods.

## Site-number correction

The originally specified site number, USGS 04288000, is **MAD RIVER NEAR MORETOWN, VT**, not the Winooski River at Montpelier (confirmed by querying the NWIS IV service and inspecting the returned `sourceInfo.siteName`). The correct site for "WINOOSKI RIVER AT MONTPELIER, VT" is **USGS 04286000** (NWS/AHPS gauge `MONV1`), which is what this analysis uses. Any prior draft material citing 04288000 for this gauge should be corrected.

## Site metadata

- USGS site: 04286000, WINOOSKI RIVER AT MONTPELIER, VT
- NWS/AHPS identifier: MONV1
- Drainage area: **397 sq mi** (source: https://waterservices.usgs.gov/nwis/site/?sites=04286000&format=rdb&siteOutput=expanded)

## NWS flood stage

- NWS/NOAA flood categories for MONV1: major=17.5 ft, moderate=16 ft, minor=15 ft, action=11 ft
- Source: https://api.water.noaa.gov/nwps/v1/gauges/monv1
- "Flood stage" as used below = **minor flood category (15 ft)**, the standard NWS AHPS threshold for the onset of flooding impacts.

## Event: July 2023

- Data source: USGS NWIS **iv** (instantaneous values)

- Peak discharge: **23100 cfs** at 2023-07-11 08:30 EDT (2023-07-11 12:30 UTC)

- Peak stage: **21.29 ft** at 2023-07-11 08:30 EDT (2023-07-11 12:30 UTC)

- Pre-event baseflow: **826 cfs** (minimum hourly-mean discharge 60-12 h before peak; occurred 2023-07-09 16:00 EDT (2023-07-09 20:00 UTC))

- **Primary recession fit** (window chosen for R^2 quality): peak -> 2023-07-13 20:30 EDT (2023-07-14 00:30 UTC) (60.0 h), tau = **36.8 h (1.53 d)**, R^2 = **0.910**, n = 721 points

- **Extended recession fit** (peak to the trough just before the next rain-driven rise interrupts recession): peak -> 2023-07-15 22:00 EDT (2023-07-16 02:00 UTC) (109.5 h), tau = **79.1 h (3.30 d)**, R^2 = **0.733**, n = 1315 points

- Duration above flood stage (15 ft): **28.2 h** (2023-07-10 16:30 EDT (2023-07-10 20:30 UTC) -> 2023-07-11 20:35 EDT (2023-07-12 00:35 UTC))

- Duration above 2x pre-event baseflow (1652 cfs): **466.2 h** (2023-07-10 02:00 EDT (2023-07-10 06:00 UTC) -> 2023-07-29 12:00 EDT (2023-07-29 16:00 UTC))


## Event: July 2024

- Data source: USGS NWIS **iv** (instantaneous values)

- Peak discharge: **11900 cfs** at 2024-07-11 03:00 EDT (2024-07-11 07:00 UTC)

- Peak stage: **14.45 ft** at 2024-07-11 03:15 EDT (2024-07-11 07:15 UTC)

- Pre-event baseflow: **466 cfs** (minimum hourly-mean discharge 60-12 h before peak; occurred 2024-07-09 15:00 EDT (2024-07-09 19:00 UTC))

- **Primary recession fit** (window chosen for R^2 quality): peak -> 2024-07-13 15:00 EDT (2024-07-13 19:00 UTC) (60.0 h), tau = **45.2 h (1.88 d)**, R^2 = **0.963**, n = 241 points

- **Extended recession fit** (peak to the trough just before the next rain-driven rise interrupts recession): peak -> 2024-07-15 23:00 EDT (2024-07-16 03:00 UTC) (116.0 h), tau = **67.4 h (2.81 d)**, R^2 = **0.938**, n = 465 points

- Duration above flood stage (15 ft): **never exceeded** (peak stage 14.45 ft stayed below the 15 ft minor-flood threshold)

- Duration above 2x pre-event baseflow (933 cfs): **242.5 h** (2024-07-10 19:00 EDT (2024-07-10 23:00 UTC) -> 2024-07-20 21:15 EDT (2024-07-21 01:15 UTC))


## Data caveats

- USGS IV data are **provisional** and subject to revision; values used here were retrieved 2026-07-22 and should be re-verified before final submission if the proposal timeline allows.
- 2023 IV data are reported at a mix of 5-minute and 15-minute intervals (USGS increases reporting frequency during high flow); 2024 data are uniformly 15-minute. Duration calculations use actual timestamp deltas on contiguous threshold-exceedance runs (not simple point counts), specifically to avoid bias from this irregular sampling.
- Both recessions are interrupted by a **secondary, rain-driven rise** roughly 4.8-4.9 days after the flood peak (2023: secondary peak ~5,900 cfs around July 17; 2024: secondary peak ~2,700 cfs around July 17). The 'extended' fit window in each case ends at the trough immediately before that secondary rise, not at a full return to baseflow -- discharge had not returned to pre-event baseflow by that point in either year.
- The 2023 recession is **not well described by a single exponential** over its full length: R^2 is excellent (>0.9) for the first ~60 h post-peak but degrades steadily thereafter (R^2 ~ 0.74 by the 117.5 h trough) as the recession visibly bends to a shallower slope, consistent with a transition from fast quickflow recession to slower baseflow-dominated recession. The primary tau reported above uses the well-fit 60 h window; the extended-window tau is reported for completeness but should be treated as a lower-quality, blended estimate.
- The 2024 recession is comparably well-described by a single exponential over almost its entire length (R^2 = 0.94-0.97 depending on window), reflecting a smaller, simpler flood pulse.
- 'Duration above 2x pre-event baseflow' is reported as the single contiguous interval containing the flood peak, not total time above threshold across the full 6-week query window. In 2023 this interval is long (~19 days) because July 2023 was an exceptionally wet month in Vermont with several additional rain pulses keeping discharge elevated well after the main flood peak; treat this number as reflecting the sustained wet regime that month, not the duration of the single flood pulse itself. The 2023 discharge record also shows two lesser precursor rain pulses (~July 2-3 and ~July 7-8) before the main flood-triggering rain; the pre-event baseflow value is measured in the partial lull between the July 7-8 pulse and the main rise, so it is somewhat elevated relative to true undisturbed baseflow (deeper lows of ~550-650 cfs occur earlier in the record). The 2024 record shows an analogous smaller precursor pulse around July 6-7.
- Peak stage and peak discharge times are within 15 minutes of each other in both events (sensor/rating-curve timing, not a real physical lag).
- The 2024 event's peak stage (14.41-14.45 ft) did **not** reach the NWS minor flood stage (15 ft) at this gauge, despite discharge (11,900 cfs) being roughly half the 2023 peak (23,100 cfs) -- consistent with 2023 being the more severe flood event at this location.
