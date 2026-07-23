# Sensitivity of Calibrated sigma_MB to a Mis-Scaled Regularizer Anchor (Subtask 1.3)

Monte Carlo study quantifying how much the censored-ML calibration of $\sigma_{MB}$ (the Subtask 1.3 mass-balance building-scale error variance, estimated on the 271-building 2023 SAR brackets) moves if the regularizing prior's anchor is mis-scaled. Produced in response to `mathAdversarialReview.md`'s cross-cutting critique: the anchor 0.30 is $sd(\ln \hat\tau)$ computed across **15 different historical flood events** (`analysis/historical_tau.py`) -- an event-to-event, watershed-scale quantity -- while $\sigma_{MB}$ is a **cross-sectional**, building-to-building quantity within a single event. "Nothing connects these two quantities' magnitudes"; the fix demotes 0.30 from a "pre-estimate" to an "order-of-magnitude anchor." This script asks what that demotion costs: if the anchor is off by 2x, how far does $\hat\sigma_{MB}$ move, and under which SAR cadences does the anchor end up doing most of the work? Script: `analysis/regularizer_misscale.py`. Run via `uv run --with numpy --with scipy --with matplotlib --with pandas analysis/regularizer_misscale.py`. Reuses the population construction, interval-censored likelihood, and censored-ML fitting machinery from `analysis/cadence_identifiability.py` (code copied in, not imported, to keep this study self-contained).

## Design

- **Population**: identical construction to `cadence_identifiability.py` (N=271 buildings, K=15 sub-basins, Dirichlet(alpha=3) cluster sizes, log(tau_k) ~ N(log(96h), 0.30^2) per sub-basin). This 0.30 is sub-basin heterogeneity in the *simulated population*, unrelated to the *regularizer-anchor* value 0.30 under test -- the coincidence in magnitude is not exploited anywhere in the estimation code.
- **True cross-sectional $\sigma_{MB} \in \{0.2, 0.3, 0.5\}$** (redrawn per rep): the anchor may happen to be right, too big, or too small relative to whichever of these is actually true -- unknown at submission time.
- **Regularizer anchor $\in \{0.15, 0.3, 0.6\}$** (half / nominal / double of the value currently in the proposal), swept *independently* of true $\sigma_{MB}$, since the anchor is fixed by the historical-gauge analysis regardless of the unknown true cross-sectional value.
- **Prior strength (prior sd on log $\sigma_{MB}$) $\in \{0.15, 0.3, 0.6\}$** (strong / nominal -- matching `cadence_identifiability.py`'s fixed 0.30 -- / weak).
- **Cadences**, reused from `cadence_identifiability.py`: `uniform_6h` (dense), `uniform_24h` (realistic moderate), `pessimistic_onfile` (single scene at 120 h post-peak, the July 16 2023 on-file analog).
- **Paired design**: for each of the 9 (cadence, sigma_true) base cells, 250 datasets are drawn once and the *same* dataset is fit under all 3x3=9 (anchor, prior strength) combinations plus one anchor-free plain-ML fit -- isolating the anchor's marginal effect from Monte Carlo sampling noise.
- **Anchor pass-through fraction**: paired log-log sensitivity of $\hat\sigma_{MB}$ to the anchor, $d(\log\hat\sigma_{MB})/d(\log\text{anchor})$, estimated between the half and double anchor settings (log-spaced, ratio 4, per replicate, then averaged). 1 = the anchor moves the estimate 1-for-1 (data contributes nothing beyond the prior, i.e. the anchor **dominates**); 0 = the estimate is anchor-invariant (data fully identifies $\sigma_{MB}$).
- Runtime: 148s total.

## Headline: anchor pass-through fraction (nominal prior strength, sd_log=0.30)

| Cadence | sigma_true=0.2 | sigma_true=0.3 | sigma_true=0.5 |
|---|---|---|---|
| dense (uniform 6 h) | 0.02 (sd 0.00) | 0.02 (sd 0.00) | 0.02 (sd 0.00) |
| realistic moderate | 0.03 (sd 0.00) | 0.02 (sd 0.00) | 0.02 (sd 0.00) |
| pessimistic on-file | 0.16 (sd 0.01) | 0.17 (sd 0.01) | 0.19 (sd 0.02) |

## Headline: shift in mean sigma_MB_hat when anchor doubles from nominal (0.30 -> 0.60), nominal prior strength

| Cadence | sigma_true=0.2 | sigma_true=0.3 | sigma_true=0.5 |
|---|---|---|---|
| dense (uniform 6 h) | +0.003 | +0.004 | +0.007 |
| realistic moderate | +0.004 | +0.005 | +0.008 |
| pessimistic on-file | +0.027 | +0.040 | +0.073 |

## Full grid: bias / RMSE by (cadence, sigma_true, anchor, prior strength)

| Cadence | sigma_true | anchor | prior strength | bias | RMSE | sd | frac@bound | plain-ML bias | plain-ML RMSE |
|---|---|---|---|---|---|---|---|---|---|
| dense (uniform 6 h) | 0.2 | 0.15 | strong (0.15) | -0.005 | 0.009 | 0.008 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.15 | nominal (0.3) | -0.002 | 0.009 | 0.008 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.15 | weak (0.6) | -0.001 | 0.009 | 0.009 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.3 | strong (0.15) | 0.006 | 0.010 | 0.008 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.3 | nominal (0.3) | 0.001 | 0.009 | 0.009 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.3 | weak (0.6) | -0.001 | 0.009 | 0.009 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.6 | strong (0.15) | 0.018 | 0.020 | 0.009 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.6 | nominal (0.3) | 0.004 | 0.009 | 0.009 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.2 | 0.6 | weak (0.6) | 0.000 | 0.009 | 0.009 | 0% | -0.001 | 0.009 |
| dense (uniform 6 h) | 0.3 | 0.15 | strong (0.15) | -0.016 | 0.019 | 0.011 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.15 | nominal (0.3) | -0.005 | 0.013 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.15 | weak (0.6) | -0.002 | 0.012 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.3 | strong (0.15) | -0.001 | 0.011 | 0.011 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.3 | nominal (0.3) | -0.001 | 0.012 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.3 | weak (0.6) | -0.001 | 0.012 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.6 | strong (0.15) | 0.016 | 0.020 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.6 | nominal (0.3) | 0.003 | 0.012 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.3 | 0.6 | weak (0.6) | -0.000 | 0.012 | 0.012 | 0% | -0.001 | 0.012 |
| dense (uniform 6 h) | 0.5 | 0.15 | strong (0.15) | -0.045 | 0.047 | 0.016 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.15 | nominal (0.3) | -0.015 | 0.024 | 0.018 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.15 | weak (0.6) | -0.007 | 0.020 | 0.019 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.3 | strong (0.15) | -0.022 | 0.028 | 0.017 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.3 | nominal (0.3) | -0.009 | 0.020 | 0.019 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.3 | weak (0.6) | -0.005 | 0.020 | 0.019 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.6 | strong (0.15) | 0.004 | 0.018 | 0.018 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.6 | nominal (0.3) | -0.002 | 0.019 | 0.019 | 0% | -0.003 | 0.020 |
| dense (uniform 6 h) | 0.5 | 0.6 | weak (0.6) | -0.003 | 0.019 | 0.019 | 0% | -0.003 | 0.020 |
| realistic moderate | 0.2 | 0.15 | strong (0.15) | -0.006 | 0.011 | 0.009 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.15 | nominal (0.3) | -0.002 | 0.010 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.15 | weak (0.6) | -0.001 | 0.010 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.3 | strong (0.15) | 0.007 | 0.012 | 0.009 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.3 | nominal (0.3) | 0.001 | 0.010 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.3 | weak (0.6) | -0.000 | 0.010 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.6 | strong (0.15) | 0.023 | 0.025 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.6 | nominal (0.3) | 0.005 | 0.011 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.2 | 0.6 | weak (0.6) | 0.001 | 0.010 | 0.010 | 0% | -0.001 | 0.010 |
| realistic moderate | 0.3 | 0.15 | strong (0.15) | -0.017 | 0.021 | 0.013 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.15 | nominal (0.3) | -0.004 | 0.015 | 0.015 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.15 | weak (0.6) | -0.001 | 0.015 | 0.015 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.3 | strong (0.15) | 0.000 | 0.014 | 0.014 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.3 | nominal (0.3) | 0.000 | 0.015 | 0.015 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.3 | weak (0.6) | 0.000 | 0.015 | 0.015 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.6 | strong (0.15) | 0.020 | 0.025 | 0.014 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.6 | nominal (0.3) | 0.005 | 0.016 | 0.015 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.3 | 0.6 | weak (0.6) | 0.002 | 0.015 | 0.015 | 0% | 0.000 | 0.015 |
| realistic moderate | 0.5 | 0.15 | strong (0.15) | -0.047 | 0.051 | 0.019 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.15 | nominal (0.3) | -0.015 | 0.026 | 0.021 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.15 | weak (0.6) | -0.005 | 0.023 | 0.022 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.3 | strong (0.15) | -0.023 | 0.030 | 0.020 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.3 | nominal (0.3) | -0.008 | 0.023 | 0.022 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.3 | weak (0.6) | -0.003 | 0.022 | 0.022 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.6 | strong (0.15) | 0.006 | 0.022 | 0.021 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.6 | nominal (0.3) | 0.000 | 0.022 | 0.022 | 0% | -0.002 | 0.022 |
| realistic moderate | 0.5 | 0.6 | weak (0.6) | -0.001 | 0.022 | 0.022 | 0% | -0.002 | 0.022 |
| pessimistic on-file | 0.2 | 0.15 | strong (0.15) | -0.023 | 0.027 | 0.014 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.15 | nominal (0.3) | -0.010 | 0.023 | 0.021 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.15 | weak (0.6) | -0.005 | 0.025 | 0.024 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.3 | strong (0.15) | 0.038 | 0.041 | 0.017 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.3 | nominal (0.3) | 0.010 | 0.025 | 0.023 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.3 | weak (0.6) | 0.001 | 0.025 | 0.025 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.6 | strong (0.15) | 0.154 | 0.156 | 0.023 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.6 | nominal (0.3) | 0.038 | 0.046 | 0.027 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.2 | 0.6 | weak (0.6) | 0.007 | 0.027 | 0.026 | 0% | -0.003 | 0.026 |
| pessimistic on-file | 0.3 | 0.15 | strong (0.15) | -0.069 | 0.071 | 0.018 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.15 | nominal (0.3) | -0.029 | 0.041 | 0.029 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.15 | weak (0.6) | -0.007 | 0.037 | 0.036 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.3 | strong (0.15) | 0.000 | 0.022 | 0.022 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.3 | nominal (0.3) | 0.001 | 0.033 | 0.033 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.3 | weak (0.6) | 0.002 | 0.038 | 0.038 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.6 | strong (0.15) | 0.133 | 0.136 | 0.031 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.6 | nominal (0.3) | 0.042 | 0.057 | 0.040 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.3 | 0.6 | weak (0.6) | 0.012 | 0.042 | 0.040 | 0% | 0.002 | 0.040 |
| pessimistic on-file | 0.5 | 0.15 | strong (0.15) | -0.188 | 0.190 | 0.025 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.15 | nominal (0.3) | -0.096 | 0.106 | 0.046 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.15 | weak (0.6) | -0.030 | 0.076 | 0.070 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.3 | strong (0.15) | -0.105 | 0.109 | 0.031 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.3 | nominal (0.3) | -0.046 | 0.072 | 0.055 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.3 | weak (0.6) | -0.009 | 0.077 | 0.077 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.6 | strong (0.15) | 0.056 | 0.071 | 0.043 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.6 | nominal (0.3) | 0.027 | 0.077 | 0.072 | 0% | 0.012 | 0.093 |
| pessimistic on-file | 0.5 | 0.6 | weak (0.6) | 0.016 | 0.087 | 0.086 | 0% | 0.012 | 0.093 |

## Full pass-through / shift table

| Cadence | sigma_true | prior strength | pass-through (mean, sd) | shift double (0.30->0.60) | shift half (0.30->0.15) |
|---|---|---|---|---|---|
| dense (uniform 6 h) | 0.2 | strong (0.15) | 0.08 (0.00) | +0.013 (0.000) | -0.011 (0.000) |
| dense (uniform 6 h) | 0.2 | nominal (0.3) | 0.02 (0.00) | +0.003 (0.000) | -0.003 (0.000) |
| dense (uniform 6 h) | 0.2 | weak (0.6) | 0.01 (0.00) | +0.001 (0.000) | -0.001 (0.000) |
| dense (uniform 6 h) | 0.3 | strong (0.15) | 0.08 (0.00) | +0.017 (0.001) | -0.015 (0.000) |
| dense (uniform 6 h) | 0.3 | nominal (0.3) | 0.02 (0.00) | +0.004 (0.000) | -0.004 (0.000) |
| dense (uniform 6 h) | 0.3 | weak (0.6) | 0.01 (0.00) | +0.001 (0.000) | -0.001 (0.000) |
| dense (uniform 6 h) | 0.5 | strong (0.15) | 0.07 (0.00) | +0.026 (0.001) | -0.023 (0.001) |
| dense (uniform 6 h) | 0.5 | nominal (0.3) | 0.02 (0.00) | +0.007 (0.000) | -0.007 (0.000) |
| dense (uniform 6 h) | 0.5 | weak (0.6) | 0.01 (0.00) | +0.002 (0.000) | -0.002 (0.000) |
| realistic moderate | 0.2 | strong (0.15) | 0.10 (0.00) | +0.016 (0.000) | -0.014 (0.000) |
| realistic moderate | 0.2 | nominal (0.3) | 0.03 (0.00) | +0.004 (0.000) | -0.004 (0.000) |
| realistic moderate | 0.2 | weak (0.6) | 0.01 (0.00) | +0.001 (0.000) | -0.001 (0.000) |
| realistic moderate | 0.3 | strong (0.15) | 0.09 (0.00) | +0.020 (0.001) | -0.017 (0.001) |
| realistic moderate | 0.3 | nominal (0.3) | 0.02 (0.00) | +0.005 (0.000) | -0.005 (0.000) |
| realistic moderate | 0.3 | weak (0.6) | 0.01 (0.00) | +0.001 (0.000) | -0.001 (0.000) |
| realistic moderate | 0.5 | strong (0.15) | 0.08 (0.00) | +0.029 (0.001) | -0.025 (0.001) |
| realistic moderate | 0.5 | nominal (0.3) | 0.02 (0.00) | +0.008 (0.000) | -0.008 (0.000) |
| realistic moderate | 0.5 | weak (0.6) | 0.01 (0.00) | +0.002 (0.000) | -0.002 (0.000) |
| pessimistic on-file | 0.2 | strong (0.15) | 0.50 (0.02) | +0.117 (0.008) | -0.060 (0.004) |
| pessimistic on-file | 0.2 | nominal (0.3) | 0.16 (0.01) | +0.027 (0.004) | -0.021 (0.003) |
| pessimistic on-file | 0.2 | weak (0.6) | 0.04 (0.00) | +0.006 (0.001) | -0.006 (0.001) |
| pessimistic on-file | 0.3 | strong (0.15) | 0.45 (0.02) | +0.132 (0.010) | -0.069 (0.005) |
| pessimistic on-file | 0.3 | nominal (0.3) | 0.17 (0.01) | +0.040 (0.007) | -0.030 (0.004) |
| pessimistic on-file | 0.3 | weak (0.6) | 0.05 (0.00) | +0.010 (0.002) | -0.009 (0.002) |
| pessimistic on-file | 0.5 | strong (0.15) | 0.42 (0.01) | +0.161 (0.013) | -0.084 (0.007) |
| pessimistic on-file | 0.5 | nominal (0.3) | 0.19 (0.02) | +0.073 (0.017) | -0.050 (0.010) |
| pessimistic on-file | 0.5 | weak (0.6) | 0.07 (0.01) | +0.025 (0.009) | -0.021 (0.007) |

## Sanity checks

No violations beyond a 0.03-unit tolerance: the anchor pass-through fraction increases monotonically (within MC noise) both as the prior gets stronger at fixed cadence, and as cadence sparsens at fixed prior strength, in every cell checked.

## Verdict

**Worst case at the nominal prior strength (sd_log=0.30, matching `cadence_identifiability.py`)**: a mis-scaled anchor of 0.15 produces $\hat\sigma_{MB}$ bias of -0.096 (true $\sigma_{MB}$=0.5, `pessimistic_onfile` cadence, RMSE 0.106).
**Worst case across the entire grid** (any cadence, any true sigma_MB, any prior strength): a mis-scaled anchor of 0.15 combined with a **strong** prior (sd_log=0.15) produces $\hat\sigma_{MB}$ bias of -0.188 (true $\sigma_{MB}$=0.5, `pessimistic_onfile` cadence, RMSE 0.190). Both worst cases occur, as expected, at the pessimistic single-scene cadence, where the data alone barely constrain sigma_MB and the regularizer does most of the work; the grid-wide worst case is materially larger than the nominal-strength one because a strong prior amplifies exactly the anchor error this study is stress-testing.

**Pass-through fraction (nominal prior strength, averaged over true sigma_MB)**: dense 6 h cadence 0.02, realistic-moderate 24 h cadence 0.02, pessimistic on-file cadence 0.17. Under a **strong** prior (sd_log=0.15) at the pessimistic cadence, pass-through peaks at 0.50 (true sigma_MB=0.2) -- substantial, but still short of the anchor moving the estimate fully 1-for-1; even a single overpass carries some identifying information via cross-sub-basin heterogeneity in tau_k (consistent with `cadence_identifiability.py`'s finding that this cadence degrades continuously rather than collapsing). In plain terms: under the dense and realistic-moderate cadences, doubling or halving the anchor barely moves the calibrated estimate at any prior strength tested -- the 271-building bracket data dominate. Under the pessimistic on-file cadence, a sizeable (10-50%, depending on prior strength) fraction of any anchor mis-scaling passes through to sigma_MB_hat: the data alone cannot fully correct for a bad anchor, though they are not powerless either.

**2x mis-scale, realistic-moderate cadence** (`uniform_24h`, the cadence most representative of an actual SAR tasking plan): doubling the anchor from 0.30 to 0.60 shifts mean $\hat\sigma_{MB}$ by +0.004 (true sigma_MB=0.2), +0.005 (true sigma_MB=0.3), +0.008 (true sigma_MB=0.5) -- small in absolute terms, consistent with a pass-through fraction well under 1 at this cadence.

**2x mis-scale, pessimistic on-file cadence**: doubling the anchor from 0.30 to 0.60 shifts mean $\hat\sigma_{MB}$ by +0.027 (true sigma_MB=0.2), +0.040 (true sigma_MB=0.3), +0.073 (true sigma_MB=0.5), and halving it (0.30 to 0.15) shifts it by -0.021 (true sigma_MB=0.2), -0.030 (true sigma_MB=0.3), -0.050 (true sigma_MB=0.5). Under this cadence the anchor is not a minor detail: if the July 16-analog single-scene cadence is what actually governs the 2023 calibration, the reported sigma_MB is, to a large extent, a statement about the anchor, not about the 271-building brackets.

**Honest bottom line for the proposal text**: under the two cadences that resemble a real, reasonably tasked SAR acquisition (dense and realistic-moderate), the 271-building calibration is close to anchor-robust at every prior strength tested -- a 2x mis-scaling of the demoted 0.30 anchor moves the calibrated sigma_MB by well under 0.05 in absolute terms (pass-through fraction ~0.02-0.03 at nominal prior strength, and no higher than ~0.10 even under the strong-prior setting). Under the pessimistic single-scene cadence that matches the cadence actually on file for the 2023 event, the picture is worse but not a total collapse: at nominal prior strength roughly 15-19% of any anchor mis-scaling passes through to sigma_MB_hat, rising to 40-50% if the prior is set strong. The data are not powerless even here -- pass-through never reaches 1 -- but they are not dominant either. This is exactly the situation the 'order-of-magnitude anchor' framing is meant to flag rather than hide -- the proposal should state plainly that under the on-file cadence, the reported sigma_MB leans meaningfully on the anchor's assumed order of magnitude, not on an assumption that the anchor is precisely correct, and that the prior should not be set stronger than the nominal sd_log=0.30 used elsewhere in this framework.

## Caveat

The pessimistic on-file cadence's residual identifying power (pass-through well below 1 even for a single overpass) comes from the realized sub-basin $\tau_k$ values straddling the 120 h overpass time -- i.e. some sub-basins are already dry and some are still wet at 120 h, which is informative about $\sigma_{MB}$ only because the population's $\tau_k$ spread happens to bracket that one observation time. This is a property of this study's realized population draw (seed 20260723), not a guarantee that any single-scene cadence is informative regardless of when it lands relative to the true recession timescale; a scene timed well outside the $\tau_k$ range would carry less information than reported here.

