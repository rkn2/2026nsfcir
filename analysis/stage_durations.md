# Winooski River at Montpelier, VT -- Stage-Threshold Duration Analysis

Duration of contiguous stage-threshold exceedance for USGS gauge 04286000 (WINOOSKI RIVER AT MONTPELIER, VT; NWS/AHPS identifier MONV1), generated from public USGS NWIS instantaneous-value (IV) gage-height (00065) data for the July 2023 and July 2024 flood events. All numbers below are computed by `analysis/stage_thresholds.py` from live NWIS IV queries; see that script for exact methods. Companion analysis: `analysis/winooski_recession.py` / `analysis/results.md` (discharge-based recession timescales for the same two events).

## Method

For each event and each stage threshold, two durations are reported:

1. **Peak-interval duration** -- the length of the single contiguous run of IV observations exceeding the threshold that contains the flood peak, computed from actual observation timestamp deltas (not point counts), identical logic to `contiguous_interval_containing()` in `winooski_recession.py`. This is the duration of *the flood pulse itself* staying above that stage.
2. **Total time above threshold** -- the sum of timestamp-delta durations across *all* contiguous exceedance runs in the full 6-week query window (2023-07-01..2023-08-15 or 2024-07-01..2024-08-15), which can exceed the peak-interval duration if the river re-exceeds the threshold during subsequent, separate rain pulses.

## Event: July 2023

- Peak stage: **21.29 ft** at 2023-07-11 08:30 EDT (2023-07-11 12:30 UTC)


| Threshold (ft) | NWS category | Peak-interval duration (h) | Peak interval start | Peak interval end | Total time above threshold in window (h) |
|---:|---|---:|---|---|---:|
| 11.0 | NWS action stage | 47.3 | 2023-07-10 10:00 EDT (2023-07-10 14:00 UTC) | 2023-07-12 09:15 EDT (2023-07-12 13:15 UTC) | 47.4 |
| 12.0 | — | 41.8 | 2023-07-10 10:45 EDT (2023-07-10 14:45 UTC) | 2023-07-12 04:30 EDT (2023-07-12 08:30 UTC) | 41.8 |
| 13.0 | — | 37.2 | 2023-07-10 12:00 EDT (2023-07-10 16:00 UTC) | 2023-07-12 01:10 EDT (2023-07-12 05:10 UTC) | 37.2 |
| 14.0 | — | 32.8 | 2023-07-10 14:00 EDT (2023-07-10 18:00 UTC) | 2023-07-11 22:40 EDT (2023-07-12 02:40 UTC) | 32.8 |
| 15.0 | NWS minor flood stage | 28.2 | 2023-07-10 16:30 EDT (2023-07-10 20:30 UTC) | 2023-07-11 20:35 EDT (2023-07-12 00:35 UTC) | 28.2 |
| 16.0 | NWS moderate flood stage | 25.6 | 2023-07-10 17:15 EDT (2023-07-10 21:15 UTC) | 2023-07-11 18:45 EDT (2023-07-11 22:45 UTC) | 25.6 |
| 17.5 | NWS major flood stage | 22.4 | 2023-07-10 18:00 EDT (2023-07-10 22:00 UTC) | 2023-07-11 16:20 EDT (2023-07-11 20:20 UTC) | 22.4 |

## Event: July 2024

- Peak stage: **14.45 ft** at 2024-07-11 03:15 EDT (2024-07-11 07:15 UTC)


| Threshold (ft) | NWS category | Peak-interval duration (h) | Peak interval start | Peak interval end | Total time above threshold in window (h) |
|---:|---|---:|---|---|---:|
| 11.0 | NWS action stage | 14.5 | 2024-07-10 22:15 EDT (2024-07-11 02:15 UTC) | 2024-07-11 12:30 EDT (2024-07-11 16:30 UTC) | 14.5 |
| 12.0 | — | 9.0 | 2024-07-10 23:00 EDT (2024-07-11 03:00 UTC) | 2024-07-11 07:45 EDT (2024-07-11 11:45 UTC) | 9.0 |
| 13.0 | — | 5.2 | 2024-07-11 01:15 EDT (2024-07-11 05:15 UTC) | 2024-07-11 06:15 EDT (2024-07-11 10:15 UTC) | 5.2 |
| 14.0 | — | 2.5 | 2024-07-11 02:15 EDT (2024-07-11 06:15 UTC) | 2024-07-11 04:30 EDT (2024-07-11 08:30 UTC) | 2.5 |
| 15.0 | NWS minor flood stage | never exceeded | — | — | 0.0 |
| 16.0 | NWS moderate flood stage | never exceeded | — | — | 0.0 |
| 17.5 | NWS major flood stage | never exceeded | — | — | 0.0 |

## Data caveats

- USGS IV data are **provisional** and subject to revision; values used here were retrieved 2026-07-22 and should be re-verified before final submission if the proposal timeline allows.
- 2023 IV data are reported at a mix of 5-minute and 15-minute intervals (USGS increases reporting frequency during high flow); 2024 data are uniformly 15-minute. Duration calculations use actual timestamp deltas on contiguous threshold-exceedance runs (not simple point counts), specifically to avoid bias from this irregular sampling.
- The July 2023 event has a **secondary, rain-driven discharge rise roughly 5 days after the flood peak** (see `results.md`, secondary discharge peak ~5,900 cfs around July 17, 2023). In the computed stage record, however, that secondary pulse does **not** re-cross any of the thresholds analyzed here (11-17.5 ft): total time above threshold across the full window equals the peak-interval duration for every threshold in 2023, to within 0.1 h (a single brief, isolated blip at the 11 ft threshold). So for the thresholds in this table, the peak-interval duration is effectively the total flood-duration number for the event; this should not be assumed to hold at lower stage thresholds not analyzed here.
- 'Never exceeded' for a given threshold means the event's peak stage did not reach that threshold at this gauge; this occurs for some thresholds in the smaller July 2024 event.
- All timestamps are USGS-reported instantaneous values at the gauge (river stage), not a direct measurement of standing water in downtown Montpelier streets or building basements -- see the separate qualitative evidence on downtown standing-water persistence for that claim.
