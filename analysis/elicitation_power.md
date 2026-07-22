# Elicitation-Rate Power Study for Subtask 2.3

Monte Carlo power study asking whether the Subtask 2.1 KPI (nonzero, mechanism-attributed $\delta_i^{\mathrm{comm}}$ adjustments for at least 20% of the 271 buildings) is enough for the Subtask 2.3 community-vs-GIS-only comparison to reliably detect a difference, and what the minimum useful elicitation rate is. Script: `analysis/elicitation_power.py`. Run via `uv run --with numpy --with scipy --with matplotlib --with pandas analysis/elicitation_power.py`.

## Design

- **Population** (fixed once, seeded): N=271 buildings across K=15 sub-basins, unequal cluster sizes from a Dirichlet(alpha=3) draw, log(tau_k) per sub-basin drawn once around log(96 h) (sd 0.30 on the log scale).
  - Realized sub-basin sizes: [18, 11, 23, 44, 10, 7, 15, 9, 31, 8, 33, 7, 15, 13, 27]
  - Realized tau_k (hours): [104.3, 51.1, 104.7, 85.2, 93.9, 92.5, 90.8, 59.0, 84.3, 75.1, 98.0, 211.3, 67.5, 122.6, 201.5]
- **Truth** (redrawn per MC rep): fraction $p_{true}$ of buildings carry a real anomaly $\delta_{true}=\pm\log 2$ (70% positive / obstruction-slowed); $\log T_i = \log\tau_k + \delta_{true,i} + \mathcal{N}(0,\sigma_{MB}^2)$.
- **Elicitation** (redrawn per MC rep): a fraction $p_{elicit}$ of buildings get a nonzero stated adjustment. Of those, fraction $q$ are correct (drawn preferentially from the true-anomaly pool, $\delta_{stated}=\delta_{true}$), fraction $1-q$ are spurious ($\delta_{stated}=\pm\log 2$, unbiased random sign, placed on a building with no real anomaly). If $q\cdot p_{elicit}\cdot N$ exceeds the true-anomaly pool ($p_{true}\cdot N$), correct calls are capped at the pool and the remainder falls back to spurious.
- **Observation**: one shared SAR acquisition schedule per MC rep, overpasses every $w$ hours from a random phase over a 336-hour (14-day) window; brackets $[L_i,U_i]$ from consecutive overpasses, right-censored past the window.
- **Scoring**: per-building interval-censored log predictive score under community vs. GIS-only priors, same pre-calibrated $\sigma_{MB}$ for both. Test statistic: mean(score$_{comm}$ - score$_{GIS}$) over 271 buildings. Inference: cluster bootstrap over the 15 sub-basins, 1000 resamples, detect = 95% percentile CI excludes 0 in favor of community.
- **Sweep**: $p_{elicit}\in\{0.05, 0.1, 0.2, 0.3, 0.5\}$, $w\in\{6, 24, 72, 168\}$ hours, $\sigma_{MB}\in\{0.3, 0.5, 0.8\}$, $q\in\{0.5, 0.8\}$, $p_{true}=0.25$ (120 cells), 500 MC reps/cell, 1000 bootstrap reps/MC rep.
- **Central point** used for the headline figure and the minimum-$p_{elicit}$ answer: $p_{elicit}=0.2$, $w=72$ h (a moderate, roughly 3-day acquisition cadence -- not the densest 6-h case), $\sigma_{MB}=0.5$ (middle of the grid), $q=0.8$ (elicitation mostly mechanism-correct).
- **Sensitivity**: $p_{true}=0.1$ re-run at the central point (fewer real anomalies to find caps how much 'correct' elicitation is even possible).
- **Assumption caveat**: sigma_MB is held fixed and shared between the truth-generating process and both scoring priors ('fine for a power study' per task spec); a real analysis recalibrates sigma_MB per condition. Building-level scores are floored at -50 log-score units purely to keep cluster-bootstrap means finite; see `frac_capped` diagnostic column for how often the correct-elicitation pool was exhausted.

## Central-point result

**Scope note**: 'detectable' here means the *aggregate* cluster-robust mean-score-difference test specified in the task design (Subtask 2.3's own KPI text asks only for a win on *at least one* disaggregated sub-group, which is a weaker bar than beating GIS-only on average across all 271 buildings -- so the power numbers below should be read as a power study of the aggregate comparison, arguably conservative relative to the literal sub-group KPI, not as the power of the Subtask 2.1 20% headcount KPI itself, which this study instead informs indirectly by asking whether 20% elicited buildings makes the aggregate comparison practically resolvable).

At the central assumptions ($p_{elicit}$=20%, $w$=72 h, $\sigma_{MB}$=0.5, $q$=0.8, $p_{true}$=0.25): **power = 70%** (mean score-difference 0.0896 log-score units, community favoring GIS-only when positive). This number is anchored to the specific choice $w$=72h; power at $p_{elicit}$=20%, $\sigma_{MB}$=0.5, $q$=0.8 ranges 62%-77% across the full cadence grid (w=6 to 168h) -- treat 70% as a mid-cadence point estimate, not a single robust number.

Sensitivity at $p_{true}$=0.1 (same central $p_{elicit}/w/\sigma/q$): **power = 3%** (frac. of reps where correct-elicitation pool was capped: 0%). Fewer real anomalies leaves less for correct elicitation to find, even at fixed $q$; this bounds how much the 20% KPI can achieve if true anomalies are rarer than assumed centrally.

**Minimum $p_{elicit}$ achieving >=80% power at central $w/\sigma_{MB}/q$: 30%.**

## Central-slice power table (w and sigma_MB, at q=0.8)


**sigma_MB = 0.3**

| p_elicit | w=6h | w=24h | w=72h | w=168h |
|---|---|---|---|---|
| 5% | 41% | 40% | 36% | 25% |
| 10% | 75% | 71% | 68% | 48% |
| 20% | 95% | 96% | 90% | 74% |
| 30% | 99% | 99% | 99% | 89% |
| 50% | 0% | 0% | 2% | 12% |

**sigma_MB = 0.5**

| p_elicit | w=6h | w=24h | w=72h | w=168h |
|---|---|---|---|---|
| 5% | 26% | 25% | 24% | 23% |
| 10% | 53% | 55% | 48% | 38% |
| 20% | 77% | 73% | 70% | 62% |
| 30% | 91% | 91% | 87% | 78% |
| 50% | 1% | 2% | 2% | 4% |

**sigma_MB = 0.8**

| p_elicit | w=6h | w=24h | w=72h | w=168h |
|---|---|---|---|---|
| 5% | 17% | 20% | 19% | 13% |
| 10% | 25% | 31% | 30% | 23% |
| 20% | 47% | 43% | 43% | 37% |
| 30% | 62% | 64% | 56% | 50% |
| 50% | 3% | 3% | 3% | 3% |

## Full main-grid power table

| p_elicit | w (h) | sigma_MB | q | power | mean score diff | frac capped |
|---|---|---|---|---|---|---|
| 5% | 6 | 0.3 | 0.5 | 2% | -0.0019 | 0% |
| 10% | 6 | 0.3 | 0.5 | 3% | 0.0043 | 0% |
| 20% | 6 | 0.3 | 0.5 | 2% | -0.0055 | 0% |
| 30% | 6 | 0.3 | 0.5 | 1% | -0.0250 | 0% |
| 50% | 6 | 0.3 | 0.5 | 0% | -0.0147 | 0% |
| 5% | 24 | 0.3 | 0.5 | 2% | -0.0022 | 0% |
| 10% | 24 | 0.3 | 0.5 | 3% | 0.0066 | 0% |
| 20% | 24 | 0.3 | 0.5 | 2% | -0.0079 | 0% |
| 30% | 24 | 0.3 | 0.5 | 1% | -0.0175 | 0% |
| 50% | 24 | 0.3 | 0.5 | 0% | -0.0158 | 0% |
| 5% | 72 | 0.3 | 0.5 | 3% | 0.0013 | 0% |
| 10% | 72 | 0.3 | 0.5 | 3% | 0.0027 | 0% |
| 20% | 72 | 0.3 | 0.5 | 2% | 0.0041 | 0% |
| 30% | 72 | 0.3 | 0.5 | 1% | -0.0043 | 0% |
| 50% | 72 | 0.3 | 0.5 | 2% | 0.0232 | 0% |
| 5% | 168 | 0.3 | 0.5 | 5% | 0.0017 | 0% |
| 10% | 168 | 0.3 | 0.5 | 6% | 0.0069 | 0% |
| 20% | 168 | 0.3 | 0.5 | 6% | 0.0043 | 0% |
| 30% | 168 | 0.3 | 0.5 | 6% | 0.0066 | 0% |
| 50% | 168 | 0.3 | 0.5 | 13% | 0.0384 | 0% |
| 5% | 6 | 0.5 | 0.5 | 4% | 0.0008 | 0% |
| 10% | 6 | 0.5 | 0.5 | 2% | 0.0011 | 0% |
| 20% | 6 | 0.5 | 0.5 | 2% | -0.0024 | 0% |
| 30% | 6 | 0.5 | 0.5 | 2% | -0.0073 | 0% |
| 50% | 6 | 0.5 | 0.5 | 2% | -0.0074 | 0% |
| 5% | 24 | 0.5 | 0.5 | 3% | -0.0005 | 0% |
| 10% | 24 | 0.5 | 0.5 | 3% | 0.0010 | 0% |
| 20% | 24 | 0.5 | 0.5 | 1% | -0.0020 | 0% |
| 30% | 24 | 0.5 | 0.5 | 3% | -0.0037 | 0% |
| 50% | 24 | 0.5 | 0.5 | 1% | -0.0026 | 0% |
| 5% | 72 | 0.5 | 0.5 | 3% | -0.0009 | 0% |
| 10% | 72 | 0.5 | 0.5 | 4% | 0.0045 | 0% |
| 20% | 72 | 0.5 | 0.5 | 3% | 0.0028 | 0% |
| 30% | 72 | 0.5 | 0.5 | 3% | -0.0028 | 0% |
| 50% | 72 | 0.5 | 0.5 | 2% | -0.0048 | 0% |
| 5% | 168 | 0.5 | 0.5 | 6% | 0.0017 | 0% |
| 10% | 168 | 0.5 | 0.5 | 5% | 0.0020 | 0% |
| 20% | 168 | 0.5 | 0.5 | 5% | 0.0023 | 0% |
| 30% | 168 | 0.5 | 0.5 | 4% | 0.0009 | 0% |
| 50% | 168 | 0.5 | 0.5 | 6% | 0.0142 | 0% |
| 5% | 6 | 0.8 | 0.5 | 5% | 0.0001 | 0% |
| 10% | 6 | 0.8 | 0.5 | 4% | 0.0003 | 0% |
| 20% | 6 | 0.8 | 0.5 | 3% | -0.0012 | 0% |
| 30% | 6 | 0.8 | 0.5 | 2% | -0.0014 | 0% |
| 50% | 6 | 0.8 | 0.5 | 2% | -0.0043 | 0% |
| 5% | 24 | 0.8 | 0.5 | 5% | 0.0000 | 0% |
| 10% | 24 | 0.8 | 0.5 | 5% | -0.0003 | 0% |
| 20% | 24 | 0.8 | 0.5 | 3% | -0.0013 | 0% |
| 30% | 24 | 0.8 | 0.5 | 3% | -0.0057 | 0% |
| 50% | 24 | 0.8 | 0.5 | 3% | -0.0032 | 0% |
| 5% | 72 | 0.8 | 0.5 | 5% | 0.0005 | 0% |
| 10% | 72 | 0.8 | 0.5 | 5% | 0.0012 | 0% |
| 20% | 72 | 0.8 | 0.5 | 5% | 0.0010 | 0% |
| 30% | 72 | 0.8 | 0.5 | 2% | -0.0028 | 0% |
| 50% | 72 | 0.8 | 0.5 | 3% | -0.0057 | 0% |
| 5% | 168 | 0.8 | 0.5 | 4% | 0.0003 | 0% |
| 10% | 168 | 0.8 | 0.5 | 5% | 0.0016 | 0% |
| 20% | 168 | 0.8 | 0.5 | 5% | 0.0010 | 0% |
| 30% | 168 | 0.8 | 0.5 | 4% | -0.0015 | 0% |
| 50% | 168 | 0.8 | 0.5 | 4% | -0.0004 | 0% |
| 5% | 6 | 0.3 | 0.8 | 41% | 0.0747 | 0% |
| 10% | 6 | 0.3 | 0.8 | 75% | 0.1580 | 0% |
| 20% | 6 | 0.3 | 0.8 | 95% | 0.3023 | 0% |
| 30% | 6 | 0.3 | 0.8 | 99% | 0.4588 | 0% |
| 50% | 6 | 0.3 | 0.8 | 0% | -0.0189 | 0% |
| 5% | 24 | 0.3 | 0.8 | 40% | 0.0715 | 0% |
| 10% | 24 | 0.3 | 0.8 | 71% | 0.1463 | 0% |
| 20% | 24 | 0.3 | 0.8 | 96% | 0.2861 | 0% |
| 30% | 24 | 0.3 | 0.8 | 99% | 0.4356 | 0% |
| 50% | 24 | 0.3 | 0.8 | 0% | -0.0185 | 0% |
| 5% | 72 | 0.3 | 0.8 | 36% | 0.0559 | 0% |
| 10% | 72 | 0.3 | 0.8 | 68% | 0.1202 | 0% |
| 20% | 72 | 0.3 | 0.8 | 90% | 0.2280 | 0% |
| 30% | 72 | 0.3 | 0.8 | 99% | 0.3459 | 0% |
| 50% | 72 | 0.3 | 0.8 | 2% | 0.0234 | 0% |
| 5% | 168 | 0.3 | 0.8 | 25% | 0.0326 | 0% |
| 10% | 168 | 0.3 | 0.8 | 48% | 0.0704 | 0% |
| 20% | 168 | 0.3 | 0.8 | 74% | 0.1365 | 0% |
| 30% | 168 | 0.3 | 0.8 | 89% | 0.2127 | 0% |
| 50% | 168 | 0.3 | 0.8 | 12% | 0.0333 | 0% |
| 5% | 6 | 0.5 | 0.8 | 26% | 0.0278 | 0% |
| 10% | 6 | 0.5 | 0.8 | 53% | 0.0585 | 0% |
| 20% | 6 | 0.5 | 0.8 | 77% | 0.1093 | 0% |
| 30% | 6 | 0.5 | 0.8 | 91% | 0.1662 | 0% |
| 50% | 6 | 0.5 | 0.8 | 1% | -0.0050 | 0% |
| 5% | 24 | 0.5 | 0.8 | 25% | 0.0256 | 0% |
| 10% | 24 | 0.5 | 0.8 | 55% | 0.0576 | 0% |
| 20% | 24 | 0.5 | 0.8 | 73% | 0.1040 | 0% |
| 30% | 24 | 0.5 | 0.8 | 91% | 0.1613 | 0% |
| 50% | 24 | 0.5 | 0.8 | 2% | -0.0144 | 0% |
| 5% | 72 | 0.5 | 0.8 | 24% | 0.0233 | 0% |
| 10% | 72 | 0.5 | 0.8 | 48% | 0.0484 | 0% |
| 20% | 72 | 0.5 | 0.8 | 70% | 0.0896 | 0% |
| 30% | 72 | 0.5 | 0.8 | 87% | 0.1370 | 0% |
| 50% | 72 | 0.5 | 0.8 | 2% | -0.0026 | 0% |
| 5% | 168 | 0.5 | 0.8 | 23% | 0.0170 | 0% |
| 10% | 168 | 0.5 | 0.8 | 38% | 0.0338 | 0% |
| 20% | 168 | 0.5 | 0.8 | 62% | 0.0673 | 0% |
| 30% | 168 | 0.5 | 0.8 | 78% | 0.1032 | 0% |
| 50% | 168 | 0.5 | 0.8 | 4% | 0.0087 | 0% |
| 5% | 6 | 0.8 | 0.8 | 17% | 0.0102 | 0% |
| 10% | 6 | 0.8 | 0.8 | 25% | 0.0209 | 0% |
| 20% | 6 | 0.8 | 0.8 | 47% | 0.0426 | 0% |
| 30% | 6 | 0.8 | 0.8 | 62% | 0.0629 | 0% |
| 50% | 6 | 0.8 | 0.8 | 3% | -0.0034 | 0% |
| 5% | 24 | 0.8 | 0.8 | 20% | 0.0108 | 0% |
| 10% | 24 | 0.8 | 0.8 | 31% | 0.0216 | 0% |
| 20% | 24 | 0.8 | 0.8 | 43% | 0.0402 | 0% |
| 30% | 24 | 0.8 | 0.8 | 64% | 0.0643 | 0% |
| 50% | 24 | 0.8 | 0.8 | 3% | -0.0031 | 0% |
| 5% | 72 | 0.8 | 0.8 | 19% | 0.0101 | 0% |
| 10% | 72 | 0.8 | 0.8 | 30% | 0.0195 | 0% |
| 20% | 72 | 0.8 | 0.8 | 43% | 0.0379 | 0% |
| 30% | 72 | 0.8 | 0.8 | 56% | 0.0567 | 0% |
| 50% | 72 | 0.8 | 0.8 | 3% | -0.0032 | 0% |
| 5% | 168 | 0.8 | 0.8 | 13% | 0.0068 | 0% |
| 10% | 168 | 0.8 | 0.8 | 23% | 0.0153 | 0% |
| 20% | 168 | 0.8 | 0.8 | 37% | 0.0301 | 0% |
| 30% | 168 | 0.8 | 0.8 | 50% | 0.0446 | 0% |
| 50% | 168 | 0.8 | 0.8 | 3% | 0.0002 | 0% |

## Vacuity check

Question: can the 20% KPI be satisfied 'on paper' by eliciting a lot of adjustments that are mostly wrong (low $q$), and still look like it produces a detectable effect? Central $w$ and $\sigma_{MB}$, $p_{true}$=0.25.

| q | p_elicit | power | mean score diff |
|---|---|---|---|
| 0.0 | 20% | 0% | -0.1574 |
| 0.0 | 30% | 0% | -0.2328 |
| 0.0 | 50% | 0% | -0.3928 |
| 0.1 | 20% | 0% | -0.1281 |
| 0.1 | 30% | 0% | -0.1874 |
| 0.1 | 50% | 0% | -0.3148 |
| 0.2 | 20% | 0% | -0.0917 |
| 0.2 | 30% | 0% | -0.1400 |
| 0.2 | 50% | 0% | -0.2352 |
| 0.8 | 20% | 70% | 0.0896 |
| 0.8 | 30% | 87% | 0.1370 |
| 0.8 | 50% | 2% | -0.0026 |

At $q$=0 (all elicited adjustments spurious), power stays at or near the nominal false-positive rate regardless of how large $p_{elicit}$ is, and the mean score difference is negative (community-informed prior is worse than GIS-only) -- confirming that volume of elicitation without mechanism accuracy cannot satisfy the *comparison*, even if it satisfies a literal count-based reading of the KPI. This is why the KPI's 'mechanism-attributed' qualifier matters: a count of nonzero deltas alone is gameable, but nonzero *and* mechanism-attributed (i.e., high-$q$) is not.

## Monotonicity check

**15 raw violation(s)** beyond a 5-percentage-point tolerance (n=500 reps/cell). All are explained by a single structural mechanism in the elicitation model, not by a bug or by chance -- derived below and confirmed against the observed cells.

**Mechanism.** Correct elicited calls (delta_stated = delta_true, magnitude log 2) and spurious elicited calls (delta_stated = +/-log 2 on a no-anomaly building) push the community score in equal-and-opposite directions, because both are the same magnitude away from the GIS-only prior and truth sits at one or the other. So the net score effect is driven by the *imbalance* (n_correct - n_spurious), not by q alone. Given n_elicit = round(p_elicit * N) and n_correct = min(round(q * n_elicit), n_true):

- **Uncapped** (q * p_elicit * N <= n_true): n_correct - n_spurious = n_elicit * (2q - 1). This is exactly zero at **q = 0.5** -- which is why the entire q=0.5 row of the main grid sits at near-null power (2-6%, essentially the false-positive rate) *for every p_elicit*, not just large ones: at q=0.5 correct and spurious calls are always tied 50/50 by construction.
- **Capped** (q * p_elicit * N > n_true, i.e. p_elicit > n_true/(q*N)): n_correct is pinned at n_true while n_spurious keeps growing with p_elicit, so the imbalance shrinks linearly and crosses zero at **p_elicit = 2 * p_true = 0.502** (using the realized p_true = 0.2509 at n_true=68), *regardless of q*, and goes negative (community reliably worse than GIS-only) beyond that.

This reproduces the data: all 12 p_elicit-axis violations are the q=0.8 rows crossing p_elicit=0.3 -> 0.5, i.e. crossing the capped-breakeven point at 0.502 (0.8 * 0.5 * N = 108 > n_true=68, so the cap engages and forces 68 correct + 68 spurious calls -- an exact tie). The remaining 3 w_hours-axis violations are noise around this same near-zero-effect point (p_elicit=0.5, q in {0.5, 0.8}) or a single small (5.8pp) jitter at p_elicit=0.1, sigma=0.8 -- not a real reversal of the cadence effect.


**Consequence for the sensitivity check**: at p_true=0.10 the capped-breakeven point is at p_elicit = 0.199, which lands almost exactly on the central p_elicit=0.20 -- this is *why* the p_true=0.10 sensitivity run collapses to near-null power at the same p_elicit where the p_true=0.25 central run still gets 70%. The 20% KPI's usefulness is therefore not an absolute volume target -- it is only informative relative to how many real drainage anomalies actually exist to be found.

## Verdict

**Central-assumption power at the literal 20% KPI is roughly 62-77% depending on cadence (70% at the mid-cadence w=72h reference point), short of a conventional 80% detection standard; 30% elicited is the minimum that clears 80% power at the central cadence/sigma/quality, and only within a bounded window.** The KPI as written (a fixed 20% headcount) is defensible only jointly with its own 'mechanism-attributed' qualifier and only if real anomaly prevalence is on the order of 25% or higher; it is not defensible as a bare volume target.

Three regimes where the 20% commitment fails or needs a caveat in the text:

1. **Sparse cadence / high sigma_MB**: at w=168h (weekly-equivalent revisit) or sigma_MB=0.8 (poorly-calibrated mass balance), power at p_elicit=20%, q=0.8 falls to 62% and 43% respectively (vs. 70% central) -- below 80% even before considering the anomaly-prevalence issue below.
2. **Low real anomaly prevalence**: if true drainage-anomaly prevalence is closer to 10% than 25% (plausible -- this is a PI judgment call, not a measured quantity), 20% elicited sits almost exactly at the point where forced padding with as-many spurious as correct calls cancels the signal (3% power in the sensitivity run). This is the single biggest risk to the KPI as written.
3. **Volume without quality (the vacuity check)**: at q<=0.5 (elicitation only as likely to be right as wrong), power stays at the nominal false-positive rate and the mean score difference is *negative* (community-informed prior actively worse than GIS-only) for every p_elicit tested, including 50%. A headcount-based KPI, taken alone, could be satisfied by low-quality volume that makes the system worse; the KPI text's existing 'mechanism-attributed' and 'no-knowledge outcomes recorded as zero' language is exactly the right guard against this, and this study is evidence that guard is load-bearing, not decorative -- it should not be loosened.

**Recommendation**: either (a) keep 20% but add a sentence noting the comparison is powered under an assumed real-anomaly prevalence and cadence, with 80% power requiring closer to 30% under central assumptions and failing at sparse cadence regardless of elicitation rate, or (b) reframe the KPI around elicitation quality relative to findable anomalies. Caution on wording (b): 'at least half of calls correct' is NOT a safe threshold -- q=0.5 is precisely this study's zero-power dead zone (Table above, q=0.5 row: 2-6% power at every p_elicit tested), because at exactly 50/50 correct-vs-spurious the two cancel by construction. A quality floor has to clear q meaningfully above 0.5 (this study's q=0.8 case is where the positive results come from) to do any work; 'a majority correct' is too weak a bar and should not be the wording used if (b) is adopted.

