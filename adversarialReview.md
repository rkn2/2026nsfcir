# Adversarial review of Research_v2.tex (2026-07-22, post PI pass)

> Status 2026-07-22 evening: all items worked through; see the matching session
> section in FINDINGS.md for per-item resolutions. Items 4 and 5 remain on the
> risk register (blocked on the 2024 damage distribution and Christelle's scene
> list); everything else is resolved in the narrative.

Written as a skeptical CPS-CIR panelist would, ranked by damage potential. Line numbers refer to Research_v2.tex at commit f99f492.

## Major concerns (any one could sink it)

### 1. The cyber layer predicts one of the two inputs its own fragility surface requires
The fragility surface (Eq. 5) is bivariate in duration and depth. The mass balance model produces duration only. ICEYE supplies depth for calibration, but in planning mode and pre-storm mode, where does building-level depth come from? Eq. (6) integrates over T with d_i held fixed, and the text never says fixed from what (line ~1030). If the answer is "the 2023 ICEYE depth pattern," every operational assessment silently assumes the next flood spatially replays 2023, which contradicts the transfer claims. If the answer is a stage-to-DEM mapping from the gauge, that is a real modeling component missing from the subtask structure.

### 2. Does the system predict whether a building floods, or only how long, given that it floods?
T_i conditionality is never stated. Pre-storm mode promises building-level inundation risk maps, which is an extent prediction a lumped mass balance cannot produce. Related: "each building-scale sub-basin" (line ~668) is undefined; if sub-basins approach building scale it is per-building hydraulics by another name. The earlier draft's explicit scoping sentence (the error model characterizes duration error, not flood/no-flood prediction) was lost in the rewrite and needs to come back, with the approximate sub-basin count.

### 3. The new preliminary results contradict the scope conditions
The transfer scope condition requires inundation persisting approximately five or more days (line ~509), but the recession analysis reports 28 hours above flood stage in 2023 and tau of 37 hours. A reviewer with the gauge record, or the proposal's own figure, can argue the pilot community fails its own scope condition. Defense: building-level standing water outlasts river stage by days. Either the threshold changes or the text must distinguish river-stage duration from standing-water duration explicitly.

### 4. The 2024 held-out validation may be nearly signal-free
The 2024 peak stayed below the gauge's minor flood stage. If downtown 2024 flooding was marginal, most of the 271 buildings have near-zero duration and no damage, and the KPI "bivariate beats depth-only by a statistically significant margin on 2024" is unpowered by construction. Need the 2024 damage-state distribution in the field-campaign subsection; if heavily skewed toward no damage, the KPI language must change before submission.

### 5. Sigma_MB may be weakly identified by the very data that calibrates it
If the actual 2023 acquisition is sparse (the scene list on file is July 16, October 9, December 20), brackets are so wide that the censored ML for sigma_MB is nearly flat and both the RQ1 error model and RQ2 detectability collapse to "the data cannot tell." The framework degrades gracefully; the scientific claims do not. Blocked on Christelle's file list (FINDINGS item 4). Single biggest factual risk in the proposal.

## Moderate concerns

### 6. Fragility calibration contradicts the sensing model
Subtask 1.4 fits fragility with point-valued "ICEYE-measured duration" covariates (line ~791) while Subtask 1.3 says duration is only ever a bracket. Fix: state that fragility estimation treats duration as interval-censored (data augmentation or Monte Carlo over the bracket).

### 7. Physical component anchor is asserted, never argued
"The physical component provides initial estimates of lambda_0, gamma, zeta" needs one worked sentence showing how material curves bound a damage-state threshold. Notation: lambda_0 and gamma carry no ds subscript while zeta_ds does (line ~799); as written all damage states share one median surface.

### 8. The 80 percent elicitation KPI contradicts the elicitation protocol
Protocol defaults unknown buildings to delta = 0 (line ~947); KPI requires elicited values for at least 80 percent of 271 buildings. If defaults count the KPI is vacuous; if not it is implausible. Restate in terms of nonzero mechanism-attributed adjustments at a defensible rate.

### 9. CPS purists will ask where the tight bidirectional coupling is
Loop closes once per flood event through human decisions. Defense exists (human-in-the-loop actuation is the connected-communities reading; CIR invites partner-discipline science) but is never made. Add a paragraph naming loop timescales (hours pre-storm, per-event updating, years retrofit actuation) as deliberate design, and state that the RQs emerged from the post-2023 Montpelier engagement (CIR scores whether partners shaped the questions).

### 10. "First empirically derived duration-dependent fragility surface" invites counterexample hunting
Schwarz and Maiwald EDAC and Milanesi are now cited by this same proposal. The claim survives only via its qualifiers. Keep every qualifier every time or soften to "the first for this building stock."

## Minor items and nits
- Goals box says community priors elicited on duration AND depth; delta_comm adjusts duration only (line ~70).
- "log-scale likelihood in Subtask 1.3" (line ~706) is stale pre-rewrite wording; likelihood is now a bracket probability.
- Right censoring unhandled in text: building still wet at final scene has unbounded U_i; one sentence on survival machinery covers it.
- Typos in SAR paragraph (line ~518): "wavelenght," "in details," "a NDA" should be "an NDA," SAR defined twice.
- "compared on a comparable probabilistic basis" (line ~133); use "same probabilistic basis."
- Known incomplete (not review findings): field-damage subsection, Napolitano prior support, BESURE trailing sentence, Gantt placeholder, red placeholder figures, undefined decision-adequacy RMSE threshold (sensitivity analysis is the defensible route, now cheap once 2023 duration data is real).

## Fix order
Items 1, 2, 3 are narrative fixes needing no new information. Item 6 is a two-sentence fix. Items 7 to 10 and the nits are quick. Items 4 and 5 are blocked on the 2024 damage distribution and Christelle's scene list; keep on the risk register until those land.
