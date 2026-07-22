# Retrofit-Priority Ranking Sensitivity to Duration Estimation Error

Monte Carlo sensitivity study answering the inline comment after the Technical KPI paragraph in `Research_v2.tex` ("Project timeline and evaluation"): what building-level duration estimation error RMSE is "decision-adequate", i.e. at or below the level at which ranking uncertainty starts to affect retrofit priority decisions? Script: `analysis/ranking_sensitivity.py`. Run via `uv run --with numpy --with scipy --with matplotlib --with pandas analysis/ranking_sensitivity.py`.

## Assumptions

All portfolio and fragility parameters below are **anchored, not fitted** -- they set the scale of a synthetic-but-structurally-realistic 271-building portfolio so the threshold question is answerable pre-submission. They are not claims about the actual Montpelier building stock.

| Parameter | Value | Note |
|---|---|---|
| N buildings | 271 | matches the 271-building 2023/2024 dataset |
| K clusters (sub-basins) | 15 | Dirichlet($\alpha$=3.0) sizes, realized: [18, 11, 23, 44, 10, 7, 15, 9, 31, 8, 33, 7, 15, 13, 27] |
| True duration center | log(96 h) | cluster log-sd 0.4, building log-sd 0.5 |
| Realized median true duration | 93.9 h | used as the operating point for the hours conversion |
| Peak depth range | [0.2, 2.5] m | Gaussian-copula correlated with duration, target Spearman 0.5, realized 0.490 |
| Replacement value | LogNormal(median \$1.5M, log-sd 0.8) | independent of T, d; heavy right tail |
| Damage states / ratios | ds1=0.10, ds2=0.35, ds3=0.80 | fraction of replacement value |
| $\lambda_{0,ds}$ (median capacity at d=1m) | 24h / 72h / 240h (ds1/ds2/ds3) | anchors from physical component discussion, Subtask 1.4 |
| $\gamma_{ds}$ (baseline) | -0.4 | deeper water lowers duration threshold; robustness variant -0.2 |
| $\zeta_{ds}$ (baseline) | 0.5 | uniform across ds (no fragility-curve crossing); robustness variant 0.7 |

Expected loss: $EL_i = \sum_{ds} P(DS=ds \mid T_i,d_i) \cdot \text{ratio}_{ds} \cdot L_i$, evaluated at true $T_i, d_i$ for the true ranking, and at $\hat T_i = T_i \cdot e^{\varepsilon}$, $\varepsilon\sim N(0,\sigma_{est}^2)$, depth held at truth, for the perturbed ranking. This isolates the duration channel; depth and replacement-value estimation error are out of scope.

## Diagnostic: what drives the ranking

Portfolio-fixed spread: sd(log replacement value) = 0.790, sd(log true vulnerability V) = 0.823 (baseline fragility; ratio 0.96). The expected-loss ranking is **shaped by comparable contributions from both channels**: replacement value's spread is neither negligible nor dominant relative to the fragility-driven vulnerability spread, so duration error has real (not confound-swamped) leverage over the ranking -- this is consistent with the threshold below sitting at a moderate, not extreme, error tolerance. Because $L_i$ is drawn independent of $T_i, d_i$ here, this decomposition would look different (more loss-driven) if replacement value and duration were themselves correlated in the real portfolio (e.g. larger buildings sited in lower-lying areas), which this synthetic design does not test.


## Variant: baseline ($\gamma$=-0.4, $\zeta$=0.5)

| $\sigma_{est}$ (log) | RMSE (h), median-scale | RMSE (h), portfolio-wide | mean top-20 overlap | mean Spearman | P(top-10 changes >3) |
|---|---|---|---|---|---|
| 0.05 | 4.7 | 7.6 | 0.930 | 0.999 | 0.000 |
| 0.10 | 9.5 | 15.3 | 0.917 | 0.995 | 0.000 |
| 0.15 | 14.4 | 23.2 | 0.906 | 0.989 | 0.000 |
| 0.20 | 19.4 | 31.0 | 0.889 | 0.981 | 0.004 |
| 0.30 | 30.4 | 49.6 | 0.846 | 0.960 | 0.062 |
| 0.40 | 43.1 | 69.9 | 0.802 | 0.932 | 0.222 |
| 0.50 | 58.1 | 94.3 | 0.768 | 0.900 | 0.404 |
| 0.60 | 76.3 | 123.4 | 0.724 | 0.865 | 0.540 |
| 0.70 | 98.9 | 161.8 | 0.689 | 0.828 | 0.656 |
| 0.80 | 127.4 | 201.5 | 0.666 | 0.793 | 0.730 |
| 0.90 | 164.1 | 261.1 | 0.637 | 0.757 | 0.840 |
| 1.00 | 211.9 | 331.2 | 0.615 | 0.723 | 0.856 |
| 1.10 | 275.1 | 427.0 | 0.587 | 0.691 | 0.888 |
| 1.20 | 360.1 | 560.2 | 0.583 | 0.659 | 0.908 |

**Top-20 overlap >= 0.8 threshold: $\sigma_{est}$ = 0.406 log-units** (~44 h RMSE at the portfolio median duration T_med=94h; ~71 h RMSE computed directly over all 271 buildings -- larger, because portfolio-wide RMSE is a quadratic-mean-like statistic pulled up by above-median durations, while the median-scale figure evaluates the same sigma_est at one fixed operating point).

**Spearman >= 0.9 threshold: $\sigma_{est}$ = 0.500 log-units** (~58 h RMSE at the portfolio median duration T_med=94h; ~94 h RMSE computed directly over all 271 buildings -- larger, because portfolio-wide RMSE is a quadratic-mean-like statistic pulled up by above-median durations, while the median-scale figure evaluates the same sigma_est at one fixed operating point).


## Variant: gamma_-0.2 ($\gamma$=-0.2, $\zeta$=0.5)

| $\sigma_{est}$ (log) | RMSE (h), median-scale | RMSE (h), portfolio-wide | mean top-20 overlap | mean Spearman | P(top-10 changes >3) |
|---|---|---|---|---|---|
| 0.05 | 4.7 | 7.6 | 0.970 | 0.998 | 0.000 |
| 0.10 | 9.5 | 15.4 | 0.945 | 0.995 | 0.000 |
| 0.15 | 14.4 | 23.3 | 0.923 | 0.989 | 0.006 |
| 0.20 | 19.4 | 31.7 | 0.901 | 0.980 | 0.012 |
| 0.30 | 30.4 | 49.5 | 0.854 | 0.958 | 0.106 |
| 0.40 | 43.1 | 69.6 | 0.800 | 0.930 | 0.252 |
| 0.50 | 58.1 | 92.6 | 0.761 | 0.897 | 0.374 |
| 0.60 | 76.3 | 120.4 | 0.719 | 0.860 | 0.490 |
| 0.70 | 98.9 | 159.6 | 0.685 | 0.823 | 0.670 |
| 0.80 | 127.4 | 206.7 | 0.653 | 0.785 | 0.728 |
| 0.90 | 164.1 | 275.3 | 0.627 | 0.750 | 0.772 |
| 1.00 | 211.9 | 356.2 | 0.603 | 0.714 | 0.826 |
| 1.10 | 275.1 | 446.6 | 0.589 | 0.682 | 0.878 |
| 1.20 | 360.1 | 537.2 | 0.569 | 0.653 | 0.850 |

**Top-20 overlap >= 0.8 threshold: $\sigma_{est}$ = 0.400 log-units** (~43 h RMSE at the portfolio median duration T_med=94h; ~70 h RMSE computed directly over all 271 buildings -- larger, because portfolio-wide RMSE is a quadratic-mean-like statistic pulled up by above-median durations, while the median-scale figure evaluates the same sigma_est at one fixed operating point).

**Spearman >= 0.9 threshold: $\sigma_{est}$ = 0.491 log-units** (~57 h RMSE at the portfolio median duration T_med=94h; ~91 h RMSE computed directly over all 271 buildings -- larger, because portfolio-wide RMSE is a quadratic-mean-like statistic pulled up by above-median durations, while the median-scale figure evaluates the same sigma_est at one fixed operating point).


## Variant: zeta_0.7 ($\gamma$=-0.4, $\zeta$=0.7)

| $\sigma_{est}$ (log) | RMSE (h), median-scale | RMSE (h), portfolio-wide | mean top-20 overlap | mean Spearman | P(top-10 changes >3) |
|---|---|---|---|---|---|
| 0.05 | 4.7 | 7.6 | 0.947 | 0.999 | 0.000 |
| 0.10 | 9.5 | 15.3 | 0.929 | 0.995 | 0.000 |
| 0.15 | 14.4 | 23.1 | 0.915 | 0.990 | 0.000 |
| 0.20 | 19.4 | 31.2 | 0.898 | 0.983 | 0.008 |
| 0.30 | 30.4 | 48.8 | 0.857 | 0.963 | 0.052 |
| 0.40 | 43.1 | 69.7 | 0.812 | 0.938 | 0.120 |
| 0.50 | 58.1 | 93.2 | 0.779 | 0.907 | 0.284 |
| 0.60 | 76.3 | 123.9 | 0.734 | 0.876 | 0.380 |
| 0.70 | 98.9 | 157.2 | 0.705 | 0.840 | 0.562 |
| 0.80 | 127.4 | 206.4 | 0.675 | 0.807 | 0.650 |
| 0.90 | 164.1 | 262.5 | 0.644 | 0.773 | 0.724 |
| 1.00 | 211.9 | 333.6 | 0.633 | 0.740 | 0.690 |
| 1.10 | 275.1 | 445.3 | 0.605 | 0.707 | 0.786 |
| 1.20 | 360.1 | 549.2 | 0.590 | 0.680 | 0.832 |

**Top-20 overlap >= 0.8 threshold: $\sigma_{est}$ = 0.437 log-units** (~48 h RMSE at the portfolio median duration T_med=94h; ~78 h RMSE computed directly over all 271 buildings -- larger, because portfolio-wide RMSE is a quadratic-mean-like statistic pulled up by above-median durations, while the median-scale figure evaluates the same sigma_est at one fixed operating point).

**Spearman >= 0.9 threshold: $\sigma_{est}$ = 0.523 log-units** (~62 h RMSE at the portfolio median duration T_med=94h; ~100 h RMSE computed directly over all 271 buildings -- larger, because portfolio-wide RMSE is a quadratic-mean-like statistic pulled up by above-median durations, while the median-scale figure evaluates the same sigma_est at one fixed operating point).

## Monotonicity check

No violations beyond MC-noise tolerance: mean top-20 overlap and mean Spearman correlation decrease monotonically, and P(top-10 changes>3) increases monotonically, with $\sigma_{est}$, in every variant.

## Robustness across fragility assumptions

**The transferable result is the log-unit $\sigma_{est}$ threshold, not either hours figure** -- the hours numbers depend on which RMSE the KPI is actually measured against (see the per-variant tables above), so both are reported here but $\sigma_{est}$ is the number that should anchor the KPI wording, with a hours figure computed the same way the real validation pipeline (Subtask 1.5) will compute it.

Top-20-overlap threshold ranges **$\sigma_{est}$ 0.400-0.437 log-units** (portfolio-wide RMSE 70-78 h) across baseline, $\gamma$=-0.2, and $\zeta$=0.7 variants (one-at-a-time).

Spearman threshold ranges **$\sigma_{est}$ 0.491-0.523 log-units** (portfolio-wide RMSE 91-100 h) across the same variants.

The threshold moves with the fragility parameters but stays within a similar order of magnitude across the swept range, i.e. it is not knife-edge on the anchored $\gamma_{ds}$/$\zeta_{ds}$ choice.

## Recommendation

The stricter of the two baseline metrics is **top-20 overlap >= 0.8**, at $\sigma_{est}$ = 0.406 log-units (~71 h portfolio-wide RMSE, ~44 h at the median-duration operating point). **Recommend stating the Technical KPI threshold in log-units ($\sigma_{est} \lesssim 0.40$, baseline fragility, robustness range 0.40-0.44), with a companion hours figure of approximately 71 h computed the same way Subtask 1.5's validation pipeline will compute portfolio RMSE, rounded down for conservatism**, and citing this study for the derivation.

## Caveats

- **Synthetic portfolio.** Cluster sizes, duration/depth/value distributions are plausibly-shaped draws, not fits to the 2023 Montpelier data; re-running against the real 271-building dataset once available would sharpen (and could shift) the threshold.
- **Duration-channel-only.** Depth and replacement value are held at truth throughout; a joint sensitivity study across all three error sources would likely tighten the tolerable duration error further, since errors could compound. This study answers 'how much duration error alone can the ranking absorb', not 'how much total estimation error'.
- **Anchored, not fitted, fragility parameters.** $\lambda_{0,ds}$, $\gamma_{ds}$, $\zeta_{ds}$ are plausibility-anchored per the physical component described in Subtask 1.4, not calibrated to observations; the robustness sweep (Section above) shows the threshold order of magnitude survives one-at-a-time perturbation of $\gamma$ and $\zeta$, but a jointly mis-anchored surface is not ruled out.
- **Stylized error model.** $\hat{\ln T} = \ln T + N(0,\sigma_{est}^2)$ is i.i.d. multiplicative lognormal noise with no bias, no depth-dependence, and no spatial correlation; a real interval-censored SAR-based duration estimator's error structure (Subtask 1.3) is unlikely to match this exactly, so $\sigma_{est}$ should be read as a stylized dial, not a validated error model.
- **Replacement value drawn independent of duration/depth.** The diagnostic section above shows replacement value and fragility-driven vulnerability contribute comparably to ranking variance in this synthetic portfolio ($L_i$ independent of $T_i, d_i$ by construction); if larger/higher-value buildings in the real Montpelier stock are systematically sited in lower-lying, longer-duration areas, the real ranking could be more loss-driven than this study models, which would loosen the real threshold relative to what is reported here -- this study does not test that correlation.
- **Threshold is for the log-predictive/ranking channel only** -- it says nothing about the separate 90% credible-interval-coverage KPI target (0.85-0.95) in the same paragraph, which is a calibration criterion, not a ranking-stability criterion.

