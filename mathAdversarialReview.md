# Adversarial review of the mathematical/statistical framework (2026-07-23)

Written as a skeptical technical panelist (statistician/hydrologist) reading only
the equations and their surrounding derivations, not the narrative framing.
Line numbers refer to `Research_v3.tex` at commit 3193fb4. Ranked by how much
each issue threatens the core scientific claims, not by how easy the fix is.

## The one theme underneath most of these: sigma_MB^2 is doing four different jobs

`sigma_MB^2` is introduced once (Eq. \ref{eq:prior}, line 751) as "the
building-scale error variance of the GIS-parameterized mass balance" — a
**cross-sectional** quantity: how much individual buildings within one event
deviate from their sub-basin's predicted timescale. It is then used as:

1. The pre-estimate/regularizing prior, set from `sd(ln tau_hat) = 0.30`
   computed across **15 different historical flood events** at the
   **watershed** scale (line 766-767; confirmed in
   `analysis/historical_tau.md` line 45, which calls this "the direct
   empirical estimate" of sigma_MB with no caveat). That's a **temporal,
   event-to-event** variance of one aggregate parameter, not a
   **cross-sectional, building-to-building** residual variance within an
   event. Nothing connects these two quantities' magnitudes — they're driven
   by different physical mechanisms (storm-to-storm antecedent moisture and
   size vs. building-to-building foundation/elevation heterogeneity). Using
   one as a "pre-estimate" for the other is a numerically convenient
   borrowing, not a derived relationship.
2. The MLE-calibrated building-scale error variance from the 271-building
   2023 brackets (line 861-863) — the quantity it was actually defined to be.
3. The vessel that "absorbs" the epistemic uncertainty of the community
   adjustment `delta_i^comm`, which is otherwise treated as a fixed constant
   with no distribution of its own (line 1105-1107).
4. The "prior variance" in the per-event precision-weighted update of
   `delta_i^comm` (Eq. \ref{eq:deltaupdate}, line 1274-1280) — see item 4
   below.

Four conceptually distinct sources of uncertainty are carried by one scalar,
each overloading justified individually and in isolation, but never addressed
as a compounding modeling choice. **Fix**: name and estimate these
separately — `sigma_event^2` (event-to-event, from the 15-flood historical
analysis, used only to sanity-check `sigma_MB` is in a plausible range, not
as its prior), `sigma_MB^2` (cross-sectional, estimated only from the
271-building brackets, unregularized or regularized from a genuinely
building-scale source), and an explicit `v_delta,i` for the community
adjustment's own uncertainty if it's going to be used in a precision-weighted
update.

## Major issues

### 1. The governing equation of Subtask 1.2 doesn't solve its own ODE
Eq. \ref{eq:tau} (line 722-726) states the linear-reservoir ODE
`dS_k/dt = I_k(t) - S_k/tau_k`, then defines `tau_k_hat = S_k^max / C_k`.
These are two different physical models. The ODE's `tau_k` is the reservoir's
e-folding time constant — the coefficient in a first-order linear feedback
that produces **exponential** recession (`S(t) = S_0 e^{-t/tau}`), which is
exactly what the preliminary recession analysis empirically confirms (tau =
37h/45h, R^2 = 0.91/0.96). But "storage capacity divided by conveyance
capacity" is the drain time of a **constant-outflow bucket** (`dS/dt = I - C_k`,
a zero-order model), which empties **linearly** in finite time
`S_max/C_k` and does not produce exponential decay at all. As written, the
formula under the ODE does not parameterize the ODE above it — they're two
different governing equations, and only the first one matches the empirical
result the proposal leans on. (The text already flags this line as a
placeholder for Maggie, so this may be known — but as it stands in the
document it's presented as the finished equation, and a reviewer with a
hydrology background will catch it in one read.)
**Fix**: either define `tau_k` directly as the fitted e-folding recession
constant (consistent with what's actually validated), or replace the
storage/conveyance ratio with an actual per-unit-storage outflow rate
constant `k_k = C_k / S_k` (not `S_k^max`), so that `1/tau_k = k_k` is
dimensionally and structurally the same object as the ODE's decay coefficient.

### 2. The expected-loss integral (Eq. \ref{eq:loss}, line 1187-1192) doesn't parse against the rest of the framework
```
E[loss | a] = int_0^inf sum_i P(DS_i >= ds | T, d_i, a) * L_i * p(T | a) dT
```
Two problems:
- **T has no building subscript**, but the entire cyber layer (Eq. \ref{eq:prior},
  \ref{eq:posterior}) is explicit that each building has its own duration
  distribution with sub-basin-specific `mu_i`. Subtask 1.3 also explicitly
  states "posteriors are computed independently per building" (line 884). As
  written, this integral either assumes every building shares one common
  realized duration `T` (physically wrong — different sub-basins have
  different `tau_k`), or it's shorthand that should read
  `sum_i int P(DS_i >= ds | T_i, d_i, a) L_i p(T_i | a) dT_i`, with the
  integral nested inside the sum, one per building. As currently written it
  doesn't parse under the framework's own stated independence assumption.
- **`ds` is a free index, never bound.** `P(DS_i >= ds | ...)` needs either a
  fixed reference damage state (undefined which one) or, more correctly for
  an expected-loss calculation, a sum over damage states weighted by a
  state-specific consequence function. As written, `L_i` (replacement value)
  is multiplied by an *exceedance* probability at an unspecified threshold,
  which is neither "expected damage cost" (that needs
  `sum_ds [P(DS_i=ds) - P(DS_i=ds+1)] * Loss_i(ds)`, standard in PEER/PBEE
  loss-estimation practice) nor a clearly defined single-threshold metric.
**Fix**: rewrite as a per-building sum over damage states with an explicit
state-dependent consequence function, e.g.
`E[loss|a] = sum_i sum_ds [P(DS_i=ds|a) - P(DS_i=ds+1|a)] * Loss_i(ds)`,
with `P(DS_i=ds)` obtained from successive differences of Eq. \ref{eq:fragility}
and each `P(DS_i=ds)` itself computed by Monte Carlo integration over
building i's own posterior `p(T_i|a)`.

### 3. The per-event delta update (Eq. \ref{eq:deltaupdate}) is a fixed-gain exponential smoother mislabeled as a Bayesian filter
The update is structurally a precision-weighted combination:
`w_i^(e) = sigma_MB^2 / (sigma_MB^2 + v_i^(e))`. The direction is right (more
weight on new evidence when the new evidence — `v_i^(e)`, the event's
posterior variance — is small). But for this to be a genuine sequential
Bayesian update, the "prior variance" term needs to be the actual **tracked
uncertainty of the current `delta_i^comm` estimate**, which should *shrink*
as events accumulate (the way a Kalman filter's error covariance shrinks with
more observations). Here the same static `sigma_MB^2` is reused as the prior
term at every event, forever — so the filter never converges toward
certainty even after many events; every new event gets essentially the same
weight it would have gotten as the very first update. That's the behavior of
a fixed-gain exponential smoother (which is a legitimate, simpler choice for
a genuinely non-stationary quantity — e.g., if drainage infrastructure
changes over time and old information should decay), but it is not what
"precision-weighted Bayesian update" claims to be, and the two have different
statistical properties (a real Bayesian filter's estimate becomes
increasingly stable over time; this one has a persistent noise floor).
**Fix**: either (a) track a real running variance for `delta_i^comm` and
update both the mean and the variance recursively each event (true Kalman
form), or (b) keep the fixed-gain form but describe it honestly as an
exponentially-weighted recency-biased update and justify the fixed step size
on its own terms (e.g., "local drainage conditions are non-stationary, so a
fixed discount on old evidence is appropriate" — a defensible claim, just a
different one than what's currently written).

## Moderate issues

### 4. The fragility surface has no monotonicity constraint across damage states
Eq. \ref{eq:fragility} is fit independently per damage state `ds`
(`lambda_{0,ds}`, `gamma_ds`, `zeta_ds` all carry a `ds` subscript with no
stated relationship across states). Nothing enforces
`P(DS>=1) >= P(DS>=2) >= ... >= P(DS>=D)` at every `(T,d)`, which is required
for the exceedance probabilities to correspond to a valid ordinal damage-state
distribution. Independently-fit per-state threshold curves are a known
pitfall in the fragility literature for exactly this reason — curves can
cross, producing `P(DS>=3) > P(DS>=2)` for some `(T,d)`, which is not
interpretable as a probability distribution over damage states.
**Fix**: use an ordinal structure (e.g., proportional-odds / cumulative
probit with a shared `zeta` and monotonically ordered `lambda_{0,ds}`) so the
ordering is guaranteed by construction rather than hoped for post hoc.

### 5. The physical-anchoring vs. empirical-refinement estimation procedure for lambda_0, gamma is never made precise
The text says the physical component "provides initial estimates" and the
empirical component "refines these parameters from 2023 observations" (line
931-934), but never states the actual estimation machinery. This matters a
lot for the identification defense in the very next paragraph (line 943-951):
if `lambda_{0,ds}` and `gamma_ds` are literally fixed at their physical-model
values, the "empirical component" is really just fitting `zeta_ds`
(dispersion), and the collinearity-mitigation argument holds — but then the
downstream claim "duration adds statistically significant explanatory power
beyond depth alone" (line 943-945) is not really a data-driven significance
test, since the duration effect's size was set by outside material-science
data, not estimated from the 271-building sample. If instead the parameters
are genuinely re-fit ("refined") from the correlated 2023 data, the VIF
concern the paragraph exists to mitigate comes right back. Both readings are
internally consistent on their own; the text currently supports either
reading interchangeably, which means neither claim (the identification
defense or the significance claim) is actually pinned down.
**Fix**: state the estimation procedure explicitly — e.g., Bayesian
estimation with informative priors on `lambda_{0,ds}`, `gamma_ds` centered at
the physical-model point estimates, with prior variance set from the
material-curve literature's specimen-to-specimen scatter, updated by the 2023
likelihood. That's a real, well-defined machinery that makes both claims
simultaneously coherent (shrinkage toward the physical anchor when data is
weak, genuine data-driven refinement when the empirical signal is strong).

### 6. Depth carries no propagated uncertainty, unlike duration — an asymmetry that undercuts the UQ framing
Duration's full posterior gets Monte-Carlo-propagated through Eq.
\ref{eq:fragility} (line 954-957) and into the loss integral. Depth, from
Eq. \ref{eq:stagedem}, is a deterministic point value, even though the text
elsewhere states the slope `beta` has a "building-level error" bound by the
2023 ICEYE calibration (line 797-799) — that calibrated error is never
propagated anywhere; `d_i` is simply plugged into Eq. \ref{eq:fragility} and
Eq. \ref{eq:loss} as if exact. Given the proposal now explicitly frames
itself as a "Bayesian uncertainty-quantification framework" (goals box), this
is a visible half-implementation: one of the two fragility-surface covariates
carries full uncertainty quantification, the other doesn't, with no stated
justification for the asymmetry.
**Fix**: either propagate `beta`'s calibrated uncertainty through Eq.
\ref{eq:stagedem} into a distribution over `d_i` and include it in the same
Monte Carlo integration already being done for `T_i`, or explicitly state
that depth uncertainty is judged small relative to duration uncertainty and
say why (ideally with a number from the calibration).

### 7. The soft misclassification-error likelihood is derived but never used
Eq. \ref{eq:softbracket} (line 852-856) is introduced as a generalization of
the hard bracket, with a stated correct limiting behavior. But the actual
`sigma_MB^2` calibration (line 861-863) explicitly uses "the 271
interval-censored observations" via the **hard** bracket likelihood, not the
soft one. The more careful model is derived and then set aside; nothing in
the technical sections uses it. A reviewer will ask why it's there.
**Fix**: either use it (replace the hard-bracket MLE calibration with the
soft one, which is strictly more correct given Subtask 1.1 characterizes
per-scene classification error) or drop the equation and fold the
classification-error caveat into a sentence instead.

### 8. The MCDA weighted sum never specifies commensurability
Eq. \ref{eq:mcda} (line 1203-1204) combines expected loss reduction (a dollar
or probability-weighted-dollar quantity), retrofit cost (dollars), "historic
preservation value," and "community equity" (both undefined scales) via a
linear weighted sum. Weighted-sum MCDA requires the `V_j` to be normalized
onto commensurate scales before combining, or the weights `w_j` don't
actually encode the intended trade-off preferences (a criterion measured in
raw dollars will dominate one measured on a 1-5 preservation scale regardless
of the nominal weight). This is never addressed.
**Fix**: specify a normalization step (min-max or z-score standardization of
each `V_j`) before the weighted sum, and note this is a compensatory method
(a large gain on one criterion can offset a large loss on another) if that's
the intended behavior — or note explicitly if it isn't and an outranking
method is more appropriate for combining dollar and non-dollar criteria.

## Minor items
- The depth model uses a single global reach-scale slope `beta`, while the
  duration model is explicitly disaggregated to ten-to-twenty sub-basins.
  The two covariates feeding the same fragility surface have very different
  spatial resolutions, with no discussion of why that's acceptable.
- Eq. \ref{eq:lps} (interval-censored log predictive score) is the one
  equation in the technical sections that is unambiguously correct as
  written — standard proper scoring rule, correctly degenerates to the
  survival-function score under right-censoring. Worth knowing what's
  *not* broken, not just what is.
- The wild cluster bootstrap claim for cluster-robust inference with one
  dominant sub-basin cluster (line 952) is defensible in the econometrics
  literature (Cameron, Gelbach & Miller 2008 show it performs reasonably
  down to very few clusters) but is asserted without a citation; worth
  adding one given how load-bearing this inference machinery is for the
  "duration beats depth-only" claim.

## What I'd prioritize before submission
Items 1-3 are the ones a technically sophisticated reviewer is most likely to
catch and that most directly undermine stated claims (the governing equation
not matching its own empirical validation; the loss function not being
computable as written; the "Bayesian" update not having the properties
claimed for it). Items 4-6 matter most for the RQ1/RQ2 identification
claims specifically. Items 7-8 are lower risk but cheap to fix or cheap to
cut.
