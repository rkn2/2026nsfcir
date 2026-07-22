# SAR Cadence Identifiability of sigma_MB (Subtask 1.3)

Monte Carlo study asking, in advance of Christelle's real 2023 SAR scene list, how identifiable the mass-balance error variance $\sigma_{MB}$ is by censored maximum likelihood across plausible acquisition cadences -- including the sparse cadence actually on file (single scene ~5 days post-peak, July 16 analog, then nothing useful). Produced for `adversarialReview.md` item 5. Script: `analysis/cadence_identifiability.py`. Run via `uv run --with numpy --with scipy --with matplotlib --with pandas analysis/cadence_identifiability.py`.

## Design

- **Population** (fixed once, seeded, same construction as `elicitation_power.py`): N=271 buildings across K=15 sub-basins, unequal cluster sizes from a Dirichlet(alpha=3) draw, log(tau_k) per sub-basin drawn once around log(96 h) (sd 0.30 on the log scale). No community deltas: $\mu_i = \log(\tau_k)$, isolating SAR-driven identifiability of $\sigma_{MB}$ from the Subtask 2 community fusion.
  - Realized sub-basin sizes: [18, 11, 23, 44, 10, 7, 15, 9, 31, 8, 33, 7, 15, 13, 27]
  - Realized tau_k (hours): [104.3, 51.1, 104.7, 85.2, 93.9, 92.5, 90.8, 59.0, 84.3, 75.1, 98.0, 211.3, 67.5, 122.6, 201.5]
- **Truth** (redrawn per MC rep): $\sigma_{MB} \in \{0.3, 0.5, 0.8\}$; $\log T_i = \log(\tau_k) + \sigma_{MB}\,\mathcal{N}(0,1)$.
- **Cadence scenarios** (hours post-peak, 336 h / 14-day window):
  - `uniform_6h` / `uniform_12h` / `uniform_24h` / `uniform_72h`: uniform revisit at the named spacing, random phase in $[0,\text{step})$ per rep.
  - `requested_2023`: Christelle's originally requested set -- overpasses at 0, 6, 12, 18, 30 h then one more at 13 days (312 h); fixed schedule.
  - `pessimistic_onfile`: a single overpass at 120 h (5 days, the July 16 2023 scene analog); fixed. Every building is then either dry by 120 h (bracket $[0,120]$) or still wet (right-censored, $L=120$, $U=\infty$). October 9 and December 20 are not modeled as additional observation times here: at ~90 and ~155 days post-peak they are far outside the recession timescale ($\tau\sim96$ h) and would land at $U=\infty$ or already-dry for essentially every building, adding no information about $\sigma_{MB}$ beyond what the July 16 scene already gives -- this is the concrete mechanism behind the adversarial-review concern, not just its label.
- **Estimation**: maximize the interval-censored log-likelihood over a global mu offset $c$ (mu_i(c) = log(tau_k) + c; true $c=0$) and $\sigma_{MB}$, by L-BFGS-B on $(c, \log\sigma_{MB})$ with multi-start (sigma inits 0.2/0.5/1.0, plus the prior center for the penalized fit), keeping the best objective value. **Search bounds**: $c\in[-1.0,1.0]$, $\sigma_{MB}\in[0.05,3.0]$.
  - **Plain ML**: no prior term.
  - **Penalized ML**: adds a lognormal prior on $\sigma_{MB}$, $\log\sigma_{MB}\sim\mathcal{N}(\log(1.3\,\sigma_{true}), 0.3^2)$ -- i.e. the 'gauge-based pre-estimate' is deliberately simulated 30% high of the true value being tested in that cell, to be honest that the pre-estimate itself carries error, not just a convenient correct anchor.
- **400 MC reps per (cadence, $\sigma_{MB}$) cell**, 6 cadences x 3 sigmas = 18 cells. Per cell, per estimator: bias and sampling sd of $\hat\sigma_{MB}$, 90% interval width (95th minus 5th percentile across reps), and the fraction of reps where the fit pins $\hat\sigma_{MB}$ within 0.001 (log scale) of a search bound.
- Runtime: 87s total for the full sweep.

## Headline table: 90% interval width of sigma_MB_hat

| Cadence | sigma_MB=0.3 (ML / pen) | sigma_MB=0.5 (ML / pen) | sigma_MB=0.8 (ML / pen) |
|---|---|---|---|
| uniform 6 h | 0.041 / 0.040 | 0.073 / 0.071 | 0.115 / 0.113 |
| uniform 12 h | 0.045 / 0.044 | 0.073 / 0.072 | 0.121 / 0.119 |
| uniform 24 h | 0.048 / 0.047 | 0.079 / 0.077 | 0.131 / 0.129 |
| uniform 72 h | 0.064 / 0.061 | 0.089 / 0.087 | 0.148 / 0.142 |
| requested 2023 | 0.298 / 0.111 | 0.115 / 0.105 | 0.143 / 0.140 |
| pessimistic on-file | 0.128 / 0.106 | 0.235 / 0.203 | 0.575 / 0.418 |

## Full results table

| Cadence | sigma_true | prior center | ML bias | ML sd | ML width90 | ML frac@bound | pen bias | pen sd | pen width90 | pen frac@bound |
|---|---|---|---|---|---|---|---|---|---|---|
| uniform 6 h | 0.3 | 0.390 | 0.000 | 0.012 | 0.041 | 0% | 0.002 | 0.012 | 0.040 | 0% |
| uniform 6 h | 0.5 | 0.650 | 0.000 | 0.022 | 0.073 | 0% | 0.003 | 0.022 | 0.071 | 0% |
| uniform 6 h | 0.8 | 1.040 | -0.002 | 0.036 | 0.115 | 0% | 0.003 | 0.035 | 0.113 | 0% |
| uniform 12 h | 0.3 | 0.390 | -0.001 | 0.014 | 0.045 | 0% | 0.000 | 0.013 | 0.044 | 0% |
| uniform 12 h | 0.5 | 0.650 | -0.000 | 0.022 | 0.073 | 0% | 0.002 | 0.022 | 0.072 | 0% |
| uniform 12 h | 0.8 | 1.040 | -0.002 | 0.036 | 0.121 | 0% | 0.003 | 0.035 | 0.119 | 0% |
| uniform 24 h | 0.3 | 0.390 | 0.000 | 0.014 | 0.048 | 0% | 0.002 | 0.014 | 0.047 | 0% |
| uniform 24 h | 0.5 | 0.650 | -0.001 | 0.023 | 0.079 | 0% | 0.002 | 0.023 | 0.077 | 0% |
| uniform 24 h | 0.8 | 1.040 | 0.002 | 0.040 | 0.131 | 0% | 0.007 | 0.039 | 0.129 | 0% |
| uniform 72 h | 0.3 | 0.390 | -0.002 | 0.020 | 0.064 | 0% | 0.002 | 0.019 | 0.061 | 0% |
| uniform 72 h | 0.5 | 0.650 | -0.001 | 0.028 | 0.089 | 0% | 0.003 | 0.027 | 0.087 | 0% |
| uniform 72 h | 0.8 | 1.040 | -0.000 | 0.047 | 0.148 | 0% | 0.008 | 0.046 | 0.142 | 0% |
| requested 2023 | 0.3 | 0.390 | -0.084 | 0.103 | 0.298 | 10% | 0.000 | 0.031 | 0.111 | 0% |
| requested 2023 | 0.5 | 0.650 | -0.008 | 0.038 | 0.115 | 0% | 0.000 | 0.034 | 0.105 | 0% |
| requested 2023 | 0.8 | 1.040 | -0.001 | 0.044 | 0.143 | 0% | 0.006 | 0.043 | 0.140 | 0% |
| pessimistic on-file | 0.3 | 0.390 | -0.010 | 0.041 | 0.128 | 0% | 0.006 | 0.034 | 0.106 | 0% |
| pessimistic on-file | 0.5 | 0.650 | -0.001 | 0.074 | 0.235 | 0% | 0.024 | 0.063 | 0.203 | 0% |
| pessimistic on-file | 0.8 | 1.040 | 0.023 | 0.181 | 0.575 | 0% | 0.077 | 0.127 | 0.418 | 0% |

## Sanity checks

**1 violation(s)** beyond a 0.02 width-unit tolerance:

- [sigma_monotonicity, plain] `requested_2023`: width90 went from 0.298 at sigma_MB=0.3 to 0.115 at sigma_MB=0.5 (should not decrease).


**Mechanism, confirmed not a bug.** For each flagged cell, this was checked against an independent fine-grid profile-likelihood scan over sigma_MB in [0.05, 1.0] (with the mu offset re-optimized at every grid point): the scan agrees with the optimizer's answer, confirming a genuine global maximum at the search bound for those replicates, not a multi-start failure. The affected cell (`requested_2023`, sigma_MB=0.3) has 10% of replicates pinned at the lower search bound (sigma_MB_hat=0.05). The cause is structural: at sigma_MB=0.3 and tau~96 h, essentially no building's true duration falls inside the 0-30 h window covered by the five early `requested_2023` overpasses, so those five scenes collapse to a single redundant 'still wet' observation for nearly every building, leaving effectively one wide bracket ([30 h, 312 h]) per building -- similar in kind to the single-overpass `pessimistic_onfile` scenario. Occasionally the realized sample noise makes a near-zero-sigma explanation genuinely the likelihood maximum given that single bracket, producing a fat lower tail that inflates width90 past the sigma_MB=0.5 cell. Practical reading: a cadence that looks dense by scene count is not informative unless its scene timing is matched to the true recession timescale; five early scenes bunched well before a building typically dries buys little over one well-timed scene. The penalized estimator is unaffected (0% at bound, monotonic width in sigma_MB), which is itself informative: this is exactly the regime the regularizing prior in Research_v2.tex Subtask 1.3 is there for.

## Verdict

**Plain ML identifies sigma_MB cleanly at all four uniform cadences (6-72 h, 0% bound-pinning everywhere) and, surprisingly, also under the pessimistic single-scene on-file cadence, which degrades continuously (wide but never bound-pinned) rather than collapsing; the one cadence that does show genuine boundary degeneracy at low true sigma_MB is the requested-2023 schedule, because its five early scenes are timed too early relative to the true recession scale to add information beyond the single late scene. The regularizing prior narrows every cell (never wider than plain ML, per the sanity check) at a bias cost tied to the assumed +30% pre-estimate error, largest in the cells where ML is weakest.**

At the requested-2023 cadence (5 dense early overpasses + one at 13 days), plain-ML 90% width is 0.298 at sigma_MB=0.3 (with 10% of reps pinned at the lower search bound -- see the sanity-check mechanism note above) but drops to 0.115-0.143 at sigma_MB=0.5-0.8 with 0% bound-pinning, comparable to the 24 h uniform benchmark. The lesson is not 'requested-2023 is bad' but that scene count alone does not guarantee identifiability: five scenes bunched well before the population's typical drying time buy almost nothing over the one later scene that actually straddles it, and this cadence's usefulness is more sigma-dependent than any of the uniform benchmarks.

Under the pessimistic on-file cadence (single scene at 5 days), plain-ML 90% width grows from 0.128 at sigma_MB=0.3 to 0.235 at sigma_MB=0.5 to 0.575 at sigma_MB=0.8 -- at the largest true sigma_MB the 90% width (0.575) approaches the true value itself, i.e. the single overpass leaves sigma_MB very weakly constrained -- but bound-pinning stays at 0% across all three true sigmas: the likelihood degrades continuously rather than collapsing to a numerical wall. This is a real, if second-order, correction to the a priori concern: a single overpass positioned near the population's typical drying time (tau~96 h) still carries some information via cross-sub-basin heterogeneity in tau_k, so the estimate stays wide but well-behaved. The penalized estimator narrows this further to a 90% width of 0.106-0.418 at the cost of a bias of 0.006 to 0.077 (systematically high, tracking the +30% pre-estimate misspecification used here) -- this is what 'degrades gracefully' concretely means: the posterior does not blow up or hit a numerical wall, it widens smoothly as cadence sparsens and falls back toward the (imperfect) gauge prior when penalized, and the framework's own text ('the detectability statements reported in Subtask 2.3 make explicit what the acquired cadence can and cannot resolve') is the correct hedge, not an overclaim, provided the pessimistic-cadence bias is disclosed rather than the width alone.

**Risk-register consequence**: if the real July 16 / Oct 9 / Dec 20 scene list is what governs the 2023 calibration, the reported sigma_MB should be presented as the penalized (prior-regularized) estimate with its bias caveat stated explicitly -- not because plain ML numerically fails on this cadence (it does not hit a search bound), but because a single scene near the recession timescale leaves the plain-ML 90% interval wide enough (comparable in width to sigma_MB itself at the sigma_MB=0.8 end of the grid) that the regularized estimate is materially more useful for the framework's downstream detectability statements, and its bias should be stated alongside its width rather than reporting width alone.

