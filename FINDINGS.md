# CIR proposal — hypothesis check findings (2026-07-21)

Read-through of root-level `Research_v2.tex` (canonical file as of 2026-07-21; `latex/Research_v2.tex` is a stale duplicate, see project memory). Findings below are the argument-vulnerability check requested before the collaborator push. Status tags: OPEN (needs a decision/fix) vs NOTED (already flagged inline in the tex by Becca, just surfaced here).

## Genuine hypothesis vulnerabilities

### 1. OPEN — resolved diagnosis, blocked on item 4 (2026-07-21, two advisor passes): the model choice depends on ICEYE's real duration resolution, not on Busse or Wauthier answering isolated questions
- σ² ("residual local drainage uncertainty," line 615) sits on the prior mean. σ²_error (lines 640–654) is defined as mass-balance prediction error (scatter of sub-basin-level τ̂_k vs. each building's ICEYE duration) but is placed in the **likelihood** (Eq. at 642–645). Prediction error belongs on the prior side; the likelihood variance should be ICEYE's own measurement noise, which doesn't currently exist as a term.
- **Per-building mass balance is not the fix and shouldn't be pursued.** A drainage timescale is inherently a lumped sub-basin quantity; making τ̂ per-building would require building-scale hydraulic simulation, exactly what this proposal is designed to avoid, and it would gut RQ2 (δ_comm exists specifically to capture the per-building deviations the sub-basin model can't see). Sub-basin τ̂_k plus per-building δ_comm is the design, not a limitation. Confirmed with Busse: she doesn't know whether σ² and σ²_error are distinct sources, but that's not a fact to look up, it's a modeling choice.
- **The real fork is two models, decided by whether ICEYE can be treated as ground truth:**
  - Model A (ICEYE exact): collapse to one prior variance calibrated on the 2023 scatter, no likelihood term. Simple, but kills the stated "Bayesian fusion" contribution, since there's nothing to fuse if ICEYE is truth.
  - Model B (ICEYE noisy): σ²_error moves to the prior as mass-balance error; likelihood variance becomes ICEYE's own measurement error. Keeps the fusion framing, needs a real error number.
- **The deciding fact is ICEYE's actual duration resolution, which is the same question already sitting in Christelle's inbox (see item 4).** ICEYE-measured duration error is dominated by revisit spacing, not instrument precision, since duration is inferred from wet/dry across successive overpasses. Even Christelle's originally-requested set (4x July 10, 1x July 11, 1x July 23) is dense on day one then has a ~12-day gap; for slow-draining pre-code masonry, "when did it finish draining" is bracketed to that gap, nowhere near 6 hours. This makes treating ICEYE as exact ground truth hard to defend, pointing toward Model B with a real, sampling-driven error term rather than a small instrument-noise one.
- **Do not send the originally-drafted Busse/Wauthier questions** — both are superseded. The Wauthier measurement-error question is subsumed by the scene-cadence question already asked in the item 4 email. The Busse temporal-vs-spatial question is moot; it's a modeling choice, not something she can answer.
- **Status**: blocked on Christelle/Young's reply (item 4). Once the real scene cadence for both events is known, that tells us whether Model A or B is defensible and how to write the equations. Don't touch the equations before then.

### 2. OPEN — CONFIRMED: collinearity admission undercuts the "duration adds explanatory power beyond depth" claim (776–780 vs. 802–804)
- Text admits duration and depth are collinear in the 2023 data (776–778).
- The core physical-layer contribution rests on showing duration adds *statistically significant* explanatory power beyond depth alone (802–803), tested on 271 buildings in one downtown core.
- Collinear predictors + geographically clustered (non-independent) buildings make this a harder identification problem than currently written. This is the load-bearing claim for the "first empirically derived duration-dependent fragility surface" outcome — a reviewer will ask how the study is powered to separate a collinear effect.
- **Needs a decision**: add a sentence acknowledging the identification strategy (e.g., variance inflation factor check, or explaining why physical-component functional form breaks the collinearity before empirical fitting).

### 3. OPEN — JUDGMENT CALL: RQ2 test may partly measure recall of the model's own answer (936–945 vs. 1059–1074)
- δ_comm is elicited by showing partners the 2023 ICEYE map and asking them to characterize where the model is "visibly wrong" (936–940) — i.e., partners see ground truth before giving their "independent" input.
- Subtask 2.3 then scores community-informed vs. GIS-only priors against ICEYE for *both* events. For 2023 this is close to fitting to the answer; for 2024 it's a better test (persistence of drainage anomalies year to year) but the proposal doesn't distinguish "independently informative community knowledge" from "accurate recall of the error pattern the research team showed them."
- **Needs a decision**: consider whether the Year 2 cross-reference step can be reframed/strengthened to address this, or whether the RQ2 framing needs a caveat about what's actually being tested.

### 4. OPEN — PENDING REPLY: ICEYE scene dates on hand don't match Christelle's originally requested set (2026-07-21)
- The scene filenames Becca has for the 2023 event are dated July 16, October 9, and December 20, 2023 (one acquisition each from ICEYE satellites X14, X26, X24).
- Christelle's original request (May 14 and June 11 emails to NASA CSDA) described a different, tighter set actually available for the July 2023 flood: "5x Strip, 1x SLEA. 4x 7/10/2023, 1x 7/11/2023, 1x 7/23/2023."
- Becca confirmed the 6-hour interval maps shown by ICEYE/Christelle on the original project call were for this same Montpelier 2023 event, so the discrepancy is real, not a mixed-up event.
- As of the most recent email in the thread (today, 2026-07-21), NASA CSDA told Christelle the scenes are now searchable/downloadable via the Earthdata explorer, pending her signing up for an account, meaning the correct tight-interval data may not have been downloaded yet even though the proposal describes 271-building ICEYE processing as complete.
- NASA CSDA access link provided by Jordan Bell (NASA MSFC): https://csdap.earthdata.nasa.gov/explore/?bbox=-73.27135%2C43.51648%2C-72.07289%2C44.73129&date=2023-06-01%2F2023-12-31&productType=SAR&collection=iceye&itemTypes=SLC%2CGRD
- **Status**: Becca sent Christelle and Young (Christelle's postdoc) an email 2026-07-21 asking for the full July 2023 file list and confirmation of which dataset is actually in hand. Awaiting reply before touching the Preliminary Results section, the "6-hour temporal resolution" claim, or the variance model in item 1, which is blocked on this same reply.

## Already flagged inline by Becca (NOTED, not new)
- Antecedent soil moisture / AMC handling for pre-storm mode — flagged to Maggie as a submission-blocking risk if unaddressed (707–721).
- Decision-adequacy RMSE threshold undefined — flagged as needing a sensitivity analysis before submission (1294–1297).
- Single-event calibration limitation — acknowledged directly in the text itself (821–826).

## Drafting incompleteness (status, not argument issues)
- 7 placeholder citations: `placeholder_nfip`, `placeholder_acs`, `placeholder_fema_bca`, `placeholder_nist_crhp`, `placeholder_fema_hmgp`, `placeholder_scc_trust`, `placeholder_scc_twotier`.
- ~10 inline TODO comments addressed to becca/maggie/christelle for SCC citations and technical detail.
- "Results from Prior NSF Support" section is a fully empty stub for all three PIs — required NSF section, can't ship blank.
- 3 TODO figures (CPS architecture diagram, Montpelier pre-code building photo, Gantt chart) plus the ICEYE zoom-mismatch issue already noted (578–586).
- "xxx PhD students" placeholder in Broader Impacts (1213, 1236).

## Next
Awaiting Christelle/Young's reply on item 4, which also unblocks item 1's variance model. Still need decisions on items 2 and 3. More of Christelle's line-edit comments may still be incoming.
