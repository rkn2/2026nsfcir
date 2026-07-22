# Winooski River at Montpelier, VT -- Multi-Event Historical Recession-Timescale (tau) Analysis

Quantifies historical variability and antecedent-moisture dependence of the exponential flood-recession timescale tau for USGS gauge 04286000 (WINOOSKI RIVER AT MONTPELIER, VT), across the full USGS instantaneous-value (IV) period of record. Produced to (1) directly quantify the LogNormal sigma_MB scale parameter the proposal pre-estimates from "historical gauge variability," and (2) test the proposal's open concern that recession is faster from drier antecedent ground. All numbers below are computed by `analysis/historical_tau.py` from live NWIS service queries; see that script for exact methods. It is a self-contained companion to, and reuses the same fitting methods as, `analysis/winooski_recession.py`.

## Method

- **Period of record used:** 1990-10-01 to 2026-07-22, the actual begin/end dates of the 00060 (discharge) instantaneous-value series at this site, retrieved live from the NWIS site series catalog (source: https://waterservices.usgs.gov/nwis/site/?sites=04286000&format=rdb&seriesCatalogOutput=true). This begins substantially earlier than the ~2007 often assumed for IV records at small gauges; gage-height (00065) IV data at this site only begins 2007-10-01, but stage was not needed for this analysis (discharge only).
- **Event detection:** the full-period IV discharge record was fetched year-by-year and resampled to hourly means. Independent flood peaks were detected by greedy peak-picking on the hourly series: repeatedly take the global maximum above **7000 cfs**, record it as an event, then mask +/-10 days around it before searching again. This guarantees every pair of accepted event peaks is more than 10 days apart. The 7000 cfs threshold was chosen (after checking 4000-10000 cfs) specifically because it lands the event count at **17**, inside the requested 8-20 range, while still independently recovering both the July 2023 and July 2024 events used as the single-event baseline in `winooski_recession.py`.
- **Per-event analysis:** for each hourly-detected candidate peak, a precise IV discharge window (-5/+13 days) was refetched at full instantaneous resolution, and the peak, pre-event baseflow, recession trough, and exponential-recession tau were computed with the **identical methods** as `winooski_recession.py`: pre-event baseflow = minimum hourly-mean discharge in [peak-60h, peak-12h]; recession trough = first post-peak point where discharge 24h later exceeds 1.15x the running minimum (i.e., the break before a secondary rain-driven rise); exponential fit ln(Q) = ln(Q_peak) - t/tau by least squares over min(60h, trough).
- **Quality control:** events with primary-fit R^2 < 0.8 are flagged and reported in the table below but excluded from all summary statistics (sd(ln tau), regression, median/IQR).
- **Verification:** this script's independent re-implementation reproduces the July 2023 and July 2024 tau values from `winooski_recession.py` (37 h and 45 h) -- see the event table.
- **Antecedent-wetness normalization:** pre-event baseflow is normalized as a dimensionless ratio to the record median hourly discharge (**451 cfs**, computed across the full 1990-10-01-2026-07-22 hourly record), and separately reported as unit discharge per drainage area (**397 sq mi**, source: https://waterservices.usgs.gov/nwis/site/?sites=04286000&format=rdb&siteOutput=expanded).

## Event table

| Event (peak date) | Season | Peak Q (cfs) | Pre-event baseflow (cfs) | Baseflow / median | tau (h) | R^2 | Fit window (h) | Status |
|---|---|---|---|---|---|---|---|---|
| Oct 1990 | fall | 7190 | 581 | 1.29 | 50.7 | 0.948 | 60.0 | well-fit |
| Apr 1994 | spring | 8320 | 4132 | 9.16 | 69.5 | 0.937 | 60.0 | well-fit |
| Dec 2000 | winter | 9570 | 313 | 0.69 | 43.9 | 0.877 | 60.0 | well-fit |
| Apr 2001 | spring | 7570 | 1910 | 4.24 | 243.8 | 0.430 | 60.0 | **poor fit, excluded from stats** |
| Apr 2005 | spring | 7110 | 2440 | 5.41 | 77.2 | 0.809 | 60.0 | well-fit |
| Jan 2006 | winter | 7130 | 1250 | 2.77 | 45.0 | 0.903 | 42.5 | well-fit |
| Apr 2011 | spring | 7490 | 1450 | 3.22 | 96.4 | 0.954 | 60.0 | well-fit |
| May 2011 | spring | 13100 | 848 | 1.88 | 50.6 | 0.746 | 45.8 | **poor fit, excluded from stats** |
| Aug 2011 (Irene) | summer | 14600 | 180 | 0.40 | 36.3 | 0.890 | 60.0 | well-fit |
| Apr 2014 | spring | 7080 | 3970 | 8.80 | 72.0 | 0.917 | 60.0 | well-fit |
| Feb 2016 | winter | 7740 | 464 | 1.03 | 47.5 | 0.901 | 60.0 | well-fit |
| Jul 2017 | summer | 7940 | 781 | 1.73 | 65.2 | 0.892 | 60.0 | well-fit |
| Apr 2019 | spring | 7400 | 1760 | 3.90 | 72.2 | 0.972 | 60.0 | well-fit |
| Jul 2023 | summer | 23100 | 826 | 1.83 | 36.8 | 0.910 | 60.0 | well-fit |
| Dec 2023 | winter | 14200 | 979 | 2.17 | 41.3 | 0.897 | 60.0 | well-fit |
| Mar 2024 | spring | 7250 | 1200 | 2.66 | 61.3 | 0.925 | 60.0 | well-fit |
| Jul 2024 | summer | 11900 | 466 | 1.03 | 45.2 | 0.963 | 60.0 | well-fit |

## Summary: variability of tau (input to sigma_MB)

Across the **15 well-fit events** (R^2 >= 0.8; 2 of 17 total events excluded for poor fit):

- tau range: **36.3 - 96.4 h**
- tau median: **50.7 h**
- tau IQR: **44.4 - 70.7 h**
- **sd(ln tau) = 0.298** (mean ln tau = 4.007; sample standard deviation, n=15)

This sd(ln tau) is the direct empirical estimate of the LogNormal scale parameter the proposal pre-estimates as sigma_MB from "historical gauge variability" -- it should be compared against whatever sigma_MB value is currently written into the proposal text.

## Antecedent-wetness relationship

Linear regression of tau (h) on normalized pre-event baseflow (ratio to record median flow) across the 15 well-fit events: slope = **3.80 h per unit baseflow-ratio**, Pearson r = **0.59**.

- Driest-antecedent well-fit event: **Aug 2011 (Irene)** (baseflow/median = 0.40), tau = **36.3 h**
- Wettest-antecedent well-fit event: **Apr 1994** (baseflow/median = 9.16), tau = **69.5 h**
- Observed shift (driest to wettest event, raw values): **+33.2 h**
- Regression-line shift over the same baseflow-ratio range: 47.2 h to 80.5 h (**+33.3 h**)


The positive/negative sign of the regression slope indicates tau is **longer** (recession is slower) when pre-event baseflow is higher relative to the record median -- i.e., wetter antecedent conditions. This is consistent with the proposal's open concern in the other direction: it is *not* simply that dry ground drains faster in the sense of shorter tau; rather, wetter-antecedent events recede more slowly (larger tau), plausibly because already-saturated soil sustains a larger, slower baseflow-fed recession limb rather than a fast quickflow-dominated one.

### tau by season (well-fit events)

| Season | Median tau (h) | Mean tau (h) | n |
|---|---|---|---|
| fall | 50.7 | 50.7 | 1 |
| spring | 72.1 | 74.8 | 6 |
| summer | 41.0 | 45.9 | 4 |
| winter | 44.4 | 44.4 | 4 |

## Data caveats

- USGS IV data are **provisional** and subject to revision; values used here were retrieved 2026-07-22 and should be re-verified before final submission if the proposal timeline allows.

- This analysis uses a **single gauge** (04286000); tau variability and the antecedent-wetness relationship are specific to this gauge's drainage and may not transfer directly to ungauged tributaries or other basins in the proposal's study area.

- The **60 h fit window** (matching winooski_recession.py) is a fixed choice, not tuned per event; some events (as in the July 2023 case in winooski_recession.py) show recession curvature beyond 60 h that a single exponential does not capture well, which is part of why R^2 screening is applied.

- Several events' 'clean' recessions are interrupted by **secondary, rain-driven rises** within the 60 h primary fit window itself (not just after the trough); where this happens it depresses the primary R^2 and can trigger the poor-fit exclusion. This is a real feature of New England storm sequencing (multi-pulse rain events), not a methodological artifact.

- The gauge has **recurring winter data gaps**: NWIS does not publish computed instantaneous discharge for extended ice-affected periods (observed here as zero-row years/months, concentrated in January-February); this is a genuine USGS data-availability limitation, not a fetch failure, and means winter flood events (if any) are systematically under-represented in this event sample.

- The event-detection threshold (7000 cfs) and separation window (10 days) were chosen to land the event count in the requested 8-20 range and to recover the two known reference events; a different threshold would give a different (overlapping but not identical) event sample and correspondingly different sd(ln tau).

- Pre-event baseflow, and therefore the antecedent-wetness ratio, is measured in a fixed 60-12 h pre-peak window; where a precursor rain pulse falls inside that window (as documented for July 2023 in winooski_recession.py) the baseflow value is somewhat elevated relative to true undisturbed baseflow, which would bias that event's position in the antecedent-wetness regression toward the wet end.
