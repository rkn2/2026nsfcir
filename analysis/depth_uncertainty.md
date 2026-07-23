# Depth (beta) vs. Duration Contribution to Fragility-Surface Probability Uncertainty

Monte Carlo variance-decomposition study answering `mathAdversarialReview.md` item 6 ("Depth carries no propagated uncertainty, unlike duration -- an asymmetry that undercuts the UQ framing"): with beta's calibrated uncertainty propagated into $d_i$ inside the SAME Monte Carlo already used for the duration posterior $T_i$, what share of $\text{Var}[P(DS\ge ds)]$ (and $\text{Var}[V]$, the expected damage ratio) is depth-driven versus duration-driven, and how does that split move with SAR revisit cadence? Script: `analysis/depth_uncertainty.py`. Run via `uv run --with numpy --with scipy --with matplotlib --with pandas analysis/depth_uncertainty.py`.

## Assumptions

All portfolio parameters are **anchored, not fitted** -- synthetic but structurally consistent with `ranking_sensitivity.py` (N=271, K=15 sub-basins) and `cadence_identifiability.py` (bracket mechanism, cluster tau_k convention). Fragility parameters (gamma, zeta, lambda_0) are the baseline variant from `ranking_sensitivity.py`, reused as-is for cross-study consistency.

| Parameter | Value | Note |
|---|---|---|
| N buildings | 271 | matches the 271-building 2023/2024 dataset |
| K clusters (sub-basins) | 15 | Dirichlet($\alpha$=3.0) sizes, realized: [18, 11, 23, 44, 10, 7, 15, 9, 31, 8, 33, 7, 15, 13, 27] |
| True duration center | log(96 h) | cluster log-sd 0.3, sigma_MB 0.3 |
| Along-channel distance $x_i$ | Uniform[0.0, 3.0] km | plausible downtown-reach range |
| Nominal depth at $\beta$=0.8 m/km, no DEM noise | Normal(1.1 m, 0.7 m) | realized range [0.00, 3.13] m, median 0.95 m |
| $\beta$ narrow range | Triangular(0.6, 0.8, 1.0) m/km | archived flat-surface-residual range, `stage_dem.md` |
| $\beta$ wide range | Triangular(0.4, 0.8, 1.2) m/km | robustness variant |
| DEM/elevation noise (secondary) | Normal(0, 0.15 m) | reported separately, not in primary decomposition |
| $\gamma$, $\zeta$, $\lambda_{0,ds}$ | -0.4, 0.5, 24h/72h/240h | anchored on `ranking_sensitivity.py` baseline |

Depth relation: $d_i(\beta,\epsilon_z) = \max(0,\, ND_i + (\beta-0.8)\cdot x_i - \epsilon_z)$, algebraically equivalent to the proposal's $d_i=\max(0,(s_g+\beta x_i)-z_i)$ once $z_i$ is defined relative to the $\beta$=0.8 reference. $\beta$ is drawn ONCE per MC replicate and shared across the whole portfolio (a single reach-scale calibration parameter), not independently per building; DEM noise $\epsilon_z$, where included, is independent per building per replicate.

## SAR brackets by cadence -- and what actually sets the duration posterior width

The duration posterior is the LogNormal(mu_i, sigma_MB=0.3) prior TRUNCATED to the observed bracket [L_i, U_i]. Its width therefore tracks whichever is tighter: the bracket (dense cadence) or the prior itself (very sparse cadence, where the bracket barely constrains anything and the posterior relaxes back toward sigma_MB). The **median posterior sd(log T)** column below makes this explicit -- it, not the bracket step in hours, is the number that actually drives the duration-vs-depth variance split reported below.

| Cadence | schedule | median bracket width, excl. censored (h) | right-censored fraction | median posterior sd(log T) |
|---|---|---|---|---|
| dense (6 h) | every 6h | 6.0 | 0.0% | 0.0186 |
| mid (24 h) | every 24h | 24.0 | 0.0% | 0.0806 |
| sparse (weekly) | every 168h | 168.0 | 0.0% | 0.2815 |
| 2023 on-file (single scene, 5d) | single scene @ 120h | 120.0 | 25.8% | 0.2438 |

(sigma_MB prior = 0.3, the ceiling posterior sd approaches as the bracket stops constraining anything.) The 6h bracket pins the duration posterior to sd(log T)~0.019 -- more than an order of magnitude tighter than sigma_MB (~16x) -- while the single-scene 2023 on-file cadence leaves it close to the sigma_MB=0.3 prior itself (sd(log T)~0.244), because a single overpass ~5 days post-peak right-censors most buildings and barely constrains the rest. This -- not any floor on sigma_MB at dense cadence -- is the mechanism behind the headline split below: dense cadence squeezes duration's contribution to a tiny fraction of the depth channel's beta-driven spread; sparse cadence lets duration widen back out toward its full prior and overtake depth.

## Headline: portfolio-median and worst-decile depth share of Var[V]

| Cadence | $\beta$ range | median $S_d$(V) | p90 $S_d$(V) | median $S_T$(V) | boundary-uncertain share |
|---|---|---|---|---|---|
| dense (6 h) | narrow 0.6 1.0 | 85.8% | 99.6% | 13.6% | 7.0% |
| dense (6 h) | wide 0.4 1.2 | 96.7% | 99.9% | 3.3% | 18.1% |
| mid (24 h) | narrow 0.6 1.0 | 30.8% | 93.2% | 70.2% | 7.0% |
| mid (24 h) | wide 0.4 1.2 | 67.7% | 98.1% | 32.7% | 18.1% |
| sparse (weekly) | narrow 0.6 1.0 | 2.7% | 58.8% | 97.2% | 7.0% |
| sparse (weekly) | wide 0.4 1.2 | 11.7% | 86.2% | 89.4% | 18.1% |
| 2023 on-file (single scene, 5d) | narrow 0.6 1.0 | 4.1% | 62.8% | 96.4% | 7.0% |
| 2023 on-file (single scene, 5d) | wide 0.4 1.2 | 16.7% | 87.1% | 85.0% | 18.1% |

## Per-damage-state breakdown (narrow $\beta$ range)

| Cadence | ds | median $S_d$ | p90 $S_d$ |
|---|---|---|---|
| dense (6 h) | 1 | 84.1% | 99.9% |
| dense (6 h) | 2 | 85.8% | 99.0% |
| dense (6 h) | 3 | 81.8% | 98.0% |
| mid (24 h) | 1 | 24.1% | 96.6% |
| mid (24 h) | 2 | 31.1% | 86.0% |
| mid (24 h) | 3 | 24.6% | 72.1% |
| sparse (weekly) | 1 | 0.2% | 90.5% |
| sparse (weekly) | 2 | 2.7% | 41.9% |
| sparse (weekly) | 3 | 1.5% | 10.8% |
| 2023 on-file (single scene, 5d) | 1 | 0.5% | 89.9% |
| 2023 on-file (single scene, 5d) | 2 | 3.4% | 50.4% |
| 2023 on-file (single scene, 5d) | 3 | 2.8% | 16.0% |

Per-damage-state shares track the aggregate V-based share closely (same underlying T, d draws), confirming the headline number isn't an artifact of aggregating across damage states.

## Diagnostic: depth share correlates with along-channel distance

Spearman correlation between $x_i$ (along-channel distance) and $S_d$(V), dense cadence / narrow $\beta$ range: **0.678**. This is the expected mechanism, not a coincidence: $d_i(\beta) - d_i(\beta_{center}) = (\beta-\beta_{center})\cdot x_i$, so a building at $x_i\approx0$ (near the gauge) is almost insensitive to $\beta$ uncertainty regardless of its depth, while a building far upstream accumulates the full slope uncertainty over its distance. Depth-driven variance is therefore concentrated in the upstream half of the domain, not spread uniformly.

## Secondary: DEM/elevation noise vs. beta as depth-error sub-channels

The primary decomposition above lumps all depth uncertainty into $S_d$; this section asks how much of THAT is beta versus DEM/elevation noise ($\sigma_z$=0.15 m per building), for the narrow $\beta$ range. Three depth draws, same T posterior: beta-only (primary decomposition), DEM-only ($\beta$ frozen at 0.8 m/km, only $\epsilon_z$ varies), and beta+DEM combined.

| Cadence | median $S_d$(V), beta only | median $S_d$(V), DEM only | median $S_d$(V), beta+DEM |
|---|---|---|---|
| dense (6 h) | 85.8% | 93.0% | 96.0% |
| sparse (weekly) | 2.7% | 7.2% | 10.7% |

DEM noise alone is **not negligible relative to beta -- if anything it is slightly larger** in this portfolio: at dense cadence it reaches 93.0% median share of Var[V] on its own, versus beta-only's 85.8%, and the combined beta+DEM share exceeds either alone. The two sources are comparable in magnitude, not one dominating the other, so treating DEM as a minor addition on top of beta would be wrong at this depth scale (sigma_z=0.15 m against typical beta-driven depth swings of a similar order). Mechanistically this tracks the depth relation: beta's contribution to a building's depth spread scales with $x_i$ (zero near the gauge, largest upstream, see the diagnostic above), while DEM noise is present at every $x_i$ including near the gauge where beta barely matters -- so which sub-channel dominates for a given building depends on its position, not a single reach-wide answer; both sources should be carried into any operational implementation, not just beta.

## Flood/no-flood boundary

- **narrow 0.6 1.0 $\beta$ range**: 19 of 271 buildings (7.0%) have P(wet) strictly between 5% and 95% -- their flood/no-flood status flips somewhere inside the admissible $\beta$ range.
- **wide 0.4 1.2 $\beta$ range**: 49 of 271 buildings (18.1%) have P(wet) strictly between 5% and 95% -- their flood/no-flood status flips somewhere inside the admissible $\beta$ range.

These boundary buildings are exactly where treating depth as exact is most misleading: a point estimate silently picks one side of a coin flip, discarding the fact that a plausible $\beta$ draw would zero out the loss entirely.

## Sanity checks

**1 violation(s)** beyond a 0.03 tolerance (flagged, not hidden):

- beta_variant=wide_0.4_1.2: median S_d(V) rose from 0.117 (sparse_weekly_168h) to 0.167 (onfile_2023) as cadence sparsened

**Mechanism, confirmed not a bug.** The monotonicity check assumes cadence sparsity ranks strictly as dense < mid < sparse-weekly < onfile-2023, but sparsity-by-scene-count is not the same as informativeness-about-duration -- exactly the lesson `cadence_identifiability.md` draws from its own `requested_2023` cadence. Here, the repeating weekly schedule (overpasses every 168 h) gives nearly EVERY building a bracket close to 168 h wide regardless of its true duration, while the single on-file scene at 120 h gives the majority of buildings (true duration below the ~96 h population median, so below the 120 h scene) a TIGHTER bracket [0, 120 h] -- only the ~26% of buildings with true duration above 120 h get right-censored (uninformative). The portfolio-median posterior sd(log T) is therefore slightly TIGHTER under the single on-file scene (0.244) than under repeating weekly revisits (0.282; see table above) -- one scene well-placed relative to the population's typical drying time can beat several scenes spaced too widely to bracket most buildings tightly. This is a real, if second-order, feature of the two schedules being compared, not a monotonicity bug in the decomposition.

## Verdict

The depth-vs-duration split is **sharply cadence-dependent**, and the two ends of that range matter for two different parts of the proposal: for the ACTUAL 2023 calibration data (single scene ~5 days post-peak, `onfile_2023` row), depth uncertainty is **small** (4.1% median share of Var[V] at the archived 0.6-1.0 m/km range, 16.7% at the wider 0.4-1.2 m/km range) -- the single on-file scene leaves duration's posterior so wide (sd(log T) close to the sigma_MB=0.3 prior itself; see table above) that it swamps depth's contribution for most of the portfolio, though the worst decile is not negligible (62.8% / 87.1% p90). For the PROSPECTIVE dense SAR cadence the proposal is arguing for (6 h, `dense_6h` row), the picture reverses: depth becomes the **dominant** channel (85.8% median / 99.6% p90 at the narrow range, 96.7% median / 99.9% p90 at the wide range), because the tight 6 h bracket pins the duration posterior more than an order of magnitude tighter than sigma_MB (sd(log T)~0.019 vs. sigma_MB=0.3, ~16x), leaving depth's beta-driven spread as the larger remaining error source. The intermediate weekly cadence sits between the two (2.7% / 11.7% median, narrow/wide). **The mechanism is the duration posterior's width relative to sigma_MB, not a floor duration hits at dense cadence: the bracket sets the posterior width when it is tighter than the prior, and the prior (sigma_MB) is the ceiling the posterior relaxes toward as the bracket stops constraining anything.**


Separately, 7.0% of buildings (narrow $\beta$ range) to 18.1% (wide range) have a flood/no-flood status that is genuinely uncertain across the admissible $\beta$ interval, independent of cadence (this is a pure depth/beta question) -- for these buildings specifically, the point-depth assumption is not a minor approximation but silently resolves a coin flip.


**Recommendation for the proposal sentence**: state the split, not a single share -- propagating $\beta$'s calibrated uncertainty into $d_i$ is negligible for interpreting the 2023 calibration itself (~4.1% of Var[V], single on-file scene) but material, likely dominant, for the operational/prospective use case the framework targets (~85.8%-96.7% of Var[V] under dense revisit). This is exactly the fix `mathAdversarialReview.md` item 6 asks for (propagate beta into the same Monte Carlo used for duration) rather than the alternative (assert depth error is small) -- the alternative would have been wrong for the regime the proposal is actually pitching.

## Caveats

- **Synthetic portfolio.** $x_i$, nominal depth, and duration distributions are plausibly-shaped draws, not fits to the 2023/2024 Montpelier data; re-running against the real 271-building dataset would sharpen (and could shift) the shares reported here.
- **Freeze-at-the-mean Sobol approximation, not the full double-loop estimator.** $S_T$ and $S_d$ are computed by freezing the OTHER variable at its Monte Carlo mean, not by the nested-loop conditional-expectation estimator that defines true first-order Sobol indices; this is a standard, cheap surrogate but $S_T+S_d$ is not guaranteed to equal 1 (the gap is the interaction term, which can be negative). Treat the reported shares as directional ("which channel dominates, roughly by how much"), not a certified variance budget.
- **Beta is a single shared draw per MC replicate, elevation z_i is not carried explicitly.** The depth relation $d_i=\max(0,ND_i+(\beta-\beta_{center})x_i-\epsilon_z)$ is algebraically exact given how $ND_i$ was constructed, but this study does not separately validate that the archived $\beta$ range (0.6-1.0 m/km, from `stage_dem.md`'s flat-surface-residual inference) is the right one to use -- `slope_sensitivity.md`'s later, more careful along-channel-network sweep narrows the admissible interval to [0.9, >=1.5] m/km, which barely overlaps the range used here. If the true admissible range sits closer to slope_sensitivity's estimate, both the depth shares AND the boundary-uncertain share reported here would need re-running with that range, not just the 0.4-1.2 m/km robustness variant tested.
- **Nominal depth (ND_i) is drawn independent of $x_i$.** A real reach's ground elevation profile is not independent of along-channel position; this simplification could over- or understate the depth/x_i correlation that drives the diagnostic section above.
- **sigma_MB fixed at 0.3, not re-estimated per cadence.** This study asks how the FIXED-sigma_MB duration posterior's width (via the bracket alone) compares to depth uncertainty; it does not fold in `cadence_identifiability.py`'s separate finding that sigma_MB itself is harder to identify under sparse cadence -- combining both effects would likely widen the duration channel further under sparse cadence, strengthening (not weakening) this study's sparse-cadence conclusion that duration dominates there.
- **DEM noise tested only for the narrow beta range, two cadences.** The secondary DEM analysis is a spot check, not a full sweep; the finding that DEM and beta are comparable-magnitude depth-error sub-channels (neither clearly dominant) should not be assumed to hold at every point on the grid, and was itself only checked for sigma_z=0.15 m -- a different assumed DEM vertical-error scale would shift this balance.
- **Fragility parameters anchored, not fitted**, per `ranking_sensitivity.py`'s own caveat (repeated here): $\lambda_{0,ds}$, $\gamma$, $\zeta$ are plausibility-anchored, not calibrated to observations.

