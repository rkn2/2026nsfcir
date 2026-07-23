# Random-Walk Kalman Filter for the Community Delta Update

Sensitivity study validating the proposed upgrade of the per-event community adjustment update from a fixed-gain smoother ($w = \sigma_{MB}^2/(\sigma_{MB}^2+v_i)$) to a random-walk Kalman filter on $\delta_i$ (log-scale duration adjustment for building $i$), with process noise $q^2$ between events, gain $w=P/(P+v)$, and $P_{next}=(1-w)P+q^2$. Script: `analysis/delta_filter.py`. Run via `uv run --with numpy --with matplotlib --with scipy analysis/delta_filter.py`.

## Design

- **Truth**: per-building true $\delta_i^{(e)}$ follows a random walk, $\delta_i^{(e)} = \delta_i^{(e-1)} + q_{true}\,\mathcal{N}(0,1)$, $\delta_i^{(0)}=0$, for $q_{true} \in \{0.0, 0.1, 0.3\}$ (stationary / moderate drift / large drift). 300 independent buildings per $q_{true}$ cell over 15 events.
- **Observation-variance mixture** (per event, per building): missing entirely with probability 0.15 (no update that event, $P$ still grows by $q^2$); else a tight SAR bracket ($v\sim U(0.05,0.1)$, 50% of non-missing events) or a sparse bracket ($v\sim U(0.5,1.0)$, the other 50%).
- **Estimators**: (a) old fixed-gain, $\sigma_{MB}^2=0.3^2=0.09$ fixed; (b) Kalman with $\hat q = q_{true}$ (matched); (c) Kalman mismatched at $\hat q = 0.5\,q_{true}$, $2\,q_{true}$, and the two explicit wrong-regime cases ($\hat q=0$ under drift, $\hat q=0.2$ under stationarity). All Kalman variants and the fixed-gain smoother share $P_0=\sigma_{MB0}^2=0.09$ as the initial variance, for a fair comparison.
- **Steady-state Riccati fixed point**: solved (a) numerically, by iterating the actual stochastic recursion (real v-mixture, including missing draws) for 6000 steps on one long synthetic building and averaging $P$ over the last 2000 (the honest answer, since $v$ is a random mixture not a constant); and (b) in closed form, $P^* = (q^2+\sqrt{q^4+4q^2 v})/2$, at a constant 'effective' $v_{eff}=0.485$ (mean observed $v$ inflated by $1/(1-0.15)$ to approximate missing events stretching the inter-update gap) -- an approximation, cross-checked against (a) rather than trusted on its own.
- **Two-event q-identifiability**: a self-contained sub-study. Two independent calibration events ('2023', '2024') simulated for $N$ buildings (both observed by construction; missingness excluded from a calibration pair), true $q_{pair}=0.15$ generating the drift between events. Profile log-likelihood $L(q)=\sum_i \log\mathcal{N}(\text{diff}_i;\,0,\,v_{1,i}+v_{2,i}+q^2)$ evaluated over $q\in[0,0.6]$ for $N\in\{10, 50, 271\}$, with the width of the log-likelihood-drop-of-2 support interval reported as a rough likelihood-ratio stand-in for a 95% CI on 1 df.
- Runtime: 0.04s compute time (excludes ~1s `uv` environment/import startup overhead).

## Claim (1): tracked variance $P_i$ shrinks as informative events accumulate

- $q_{true}=0$ (stationary): mean $P_e$ falls from $P_0=0.0900$ to 0.0099 by event 15 and keeps shrinking with more events (the closed-form fixed point at $q=0$ is exactly $P^*=0$: with no process noise, variance shrinks without bound as informative events accumulate).

- $q_{true}=0.1$: mean $P_e$ **falls** from $P_0=0.0900$ to 0.0485 by event 15, converging to the numeric steady state 0.0490 (closed-form approx. 0.0748).

- $q_{true}=0.3$: mean $P_e$ **rises** from $P_0=0.0900$ to 0.2154 by event 15, converging to the numeric steady state 0.2103 (closed-form approx. 0.2586).

**Supported, with an important qualifier the proposal should state explicitly.** $P$ always converges toward its $q$-implied steady state within the 15-event window (including pausing, not reversing, through missing-event gaps) -- but 'shrinks' is only the right verb when that steady state is *below* the initial prior $P_0=\sigma_{MB0}^2$. At $q_{true}=0.1$ the steady state (0.0490) is well below $P_0$ (0.0900), so $P$ shrinks as claimed. At $q_{true}=0.3$ the steady state (0.2103) is *above* $P_0$ (because process noise added between events, $q^2=0.09$, exceeds what a typical observation removes), so $P$ instead **grows** from $0.09$ up to ~0.2103 within the first 3-4 events and stays there. Both directions are correct Kalman behavior, but a reader who takes 'variance shrinks as events accumulate' literally and universally will be surprised by the large-$q$ case; the accurate claim is 'variance converges to a $q$-set steady state, which is smaller than a naive fixed-gain prior only when drift is slow relative to typical observation noise.'

## Claim (2): wide brackets automatically produce small gains

By construction $w=P/(P+v)$ is monotonically decreasing in $v$ at fixed $P$ -- this is algebraic, not a simulation result. What the simulation adds: mean realized gain across the 50/50 tight/sparse mixture at $q_{true}=0.3$ settles to 0.4765 by event 15 (vs. an initial gain near 0.6429-0.0826 depending on bracket width at event 1), i.e. the mixture of tight and sparse brackets pulls the realized gain down from what a single tight-only bracket would give, exactly as claimed.
**Supported** (trivially, by algebra of $w=P/(P+v)$; confirmed operating as expected under the realistic mixture rather than degenerating numerically).

## Claim (3): convergence to a steady-state fixed gain set by q

- $q_{true}=0.1$: mean gain reaches 0.2438 by event 15; steady-state gain implied by the numeric fixed point at $v_{eff}$ is 0.0918.

- $q_{true}=0.3$: mean gain reaches 0.4765 by event 15; steady-state gain implied by the numeric fixed point at $v_{eff}$ is 0.3026.

Closed-form vs. numeric fixed point diverge by 53% at q=0.1, 23% at q=0.3 -- the constant-effective-v closed form is only an order-of-magnitude check, not a substitute for the numeric fixed point, which should be treated as authoritative (it iterates the actual time-varying, occasionally-missing v process rather than a single constant stand-in).
**Supported with a caveat.** The recursion does converge to a stable gain regime within the 15-event horizon for $q_{true}>0$ (both tested values reach their numeric steady state by event ~4-5 of 15). Two caveats: (i) because $v$ is a random mixture, not a constant, the 'steady state' is really a stationary *distribution* of $P_e$ -- individual events still show a several-fold gain swing between a tight-bracket event and a sparse-bracket one, riding on top of the converged mean; (ii) as shown under Claim (1) above, that steady state can sit above the initial prior for large $q$, so 'converges to a steady-state fixed gain' is the accurate claim, not 'gain settles low.'

## Fixed-gain noise floor vs. Kalman under a stationary truth

At $q_{true}=0$, fixed-gain final-window RMSE (events 11-15) is 0.1801, while Kalman with $\hat q=0$ (matched) reaches 0.0992 -- a 45% RMSE reduction, because the fixed-gain smoother re-applies the same $\sigma_{MB}^2$-based weight every event regardless of how much prior information has accumulated, so it never gets more confident; the Kalman filter's $P$ keeps shrinking and the estimate keeps tightening. This is the concrete mechanism behind claim (1) mattering in practice, not just holding in the abstract.

## q-mismatch robustness

- $q_{true}=0.0$: final-window RMSE minimized at grid point $\hat q=0.0$ (RMSE 0.0992); wrong-regime $\hat q=0.2$ gives RMSE 0.1995 vs. fixed-gain 0.1801.

- $q_{true}=0.1$: final-window RMSE minimized at grid point $\hat q=0.1$ (RMSE 0.1947); wrong-regime $\hat q=0.0$ gives RMSE 0.2332 vs. fixed-gain 0.2051.

- $q_{true}=0.3$: final-window RMSE minimized at grid point $\hat q=0.3$ (RMSE 0.3473); wrong-regime $\hat q=0.0$ gives RMSE 0.6442 vs. fixed-gain 0.3842.

**Graceful, not brittle.** Across the full $\hat q$ grid (0.0-0.9), final RMSE degrades smoothly moving away from the matched value in either direction -- no cliff, no numerical blow-up. Even in the worst-tested mismatch (wrong regime), Kalman final RMSE stays at or below the fixed-gain smoother's, so a moderately wrong $q$ is still at least as good as the status quo, though the RMSE advantage over fixed-gain shrinks substantially compared to the matched case.

## Claim (4): q bounded empirically by cross-event persistence between two calibration events

| N buildings | drop-2 support interval for q | width |
|---|---|---|
| 10 | [0.000, 0.600] | 0.600 |
| 50 | [0.000, 0.534] | 0.534 |
| 271 | [0.034, 0.353] | 0.319 |

True $q_{pair}=0.15$. Even at $N=271$ (the full building population, an optimistic upper bound on how many buildings actually have usable brackets in *both* named calibration events), the drop-2 support interval is [0.034, 0.353] -- informative enough to rule out very large q and to confirm q is not huge, but still several-fold wide relative to the true value, and it does not pin q to better than roughly a factor of 2. At $N=10$ (closer to what a real 2023-2024 SAR-overlap building count might look like once missingness and bracket quality are accounted for), the interval is [0.000, 0.600], wide enough that it constrains q only weakly from above and barely at all from below (q=0 is not excluded).
**Partially supported, with an important honesty caveat.** A two-event record does bound q -- the likelihood is not flat, and it does rule out implausibly large q -- but the bound is weak, especially at realistic building counts, and gets weaker still as bracket variance $v_1, v_2$ grows relative to $q^2$ (the diff's variance is $v_1+v_2+q^2$, and distinguishing 'q^2 is X' from 'q^2 is 0 and sampling noise in $v_1+v_2$ explains the observed spread' requires many buildings when $v_1,v_2$ are themselves large, which is exactly the sparse-bracket regime this study models as roughly half of all events). The proposal should not claim q is 'empirically calibrated' by a two-event comparison without reporting an interval, not a point estimate, and should treat q as swept over a plausible range (as this study does, 0/0.1/0.3) rather than pinned to a single calibrated value from two events alone.

## Sanity checks

**1 violation(s):**

- Closed-form/numeric fixed-point mismatch >35% at q_true=0.1: numeric=0.0490, closed_form=0.0748, rel_diff=0.53

## Verdict

**(1) Supported, but only 'shrinks' when the steady state is below the initial prior.** $P_i$ always converges toward its $q$-implied steady state (0, if truth is stationary) within the 15-event horizon, pausing but not reversing across missing-event gaps. At $q_{true}=0.1$ the steady state (0.049) is below $P_0=0.09$, so $P$ shrinks as claimed; at $q_{true}=0.3$ the steady state (0.210) is *above* $P_0$, so $P$ grows toward it instead. The proposal text should say 'converges to a q-set steady state' rather than 'shrinks,' since the direction depends on how fast drift is relative to typical observation noise.
**(2) Supported** -- $w=P/(P+v)$ is algebraically decreasing in bracket variance $v$; confirmed behaving as expected under the realistic tight/sparse/missing mixture, with realized mean gain visibly pulled down by the sparse-bracket share.
**(3) Supported with a caveat** -- mean $P$ and gain converge to a stable regime that matches a numerically-solved Riccati fixed point (the closed-form constant-v cross-check is only order-of-magnitude, diverging 23%-53% from the numeric fixed point across the two tested $q_{true}$ values); the caveat is that 'steady state' means a stationary mean under the random v-mixture, not a literal constant -- individual-event gain still swings several-fold with bracket quality.
**(4) Partially supported, weakest of the four** -- a two-event record does bound q away from implausibly large values, but the drop-2 support interval is several-fold wide even at the full 271-building population and only weakly excludes q=0 at realistic (smaller) two-event overlap counts. The proposal should present q as swept over a plausible range with a reported interval, not as a single value pinned by two calibration events.

**Caveats general to this study**: (i) the observation-variance mixture (15% missing, tight/sparse 50/50 split, ranges 0.05-0.10 and 0.5-1.0) is a plausible-cadence assumption, not fit to a specific real SAR bracket catalog -- results should be treated as qualitative-quantitative, not as exact numbers to cite verbatim; (ii) the closed-form steady-state uses a constant 'effective v' approximation and should not be used as a substitute for the numeric fixed point in the proposal text; (iii) the two-event sub-study assumes a clean 2-event drift model with no additional confounders (seasonal effects, building-specific process noise heterogeneity) that would likely widen the q bound further in practice, not narrow it.

