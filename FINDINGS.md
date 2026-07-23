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
- **Update (2026-07-22), per conversation with Christelle**: confirms the Model B direction — her guidance is to validate against "whatever the best thing we can buy is," i.e., ICEYE now, framed as best-available rather than exact ground truth. She also wants the sensing/validation layer generalized from "ICEYE" specifically to a **suite of SAR data** (Sentinel-1, ICEYE, Capella, or any other SAR source as available/affordable), fused at whatever level a given community or event has, not a fixed ICEYE-exclusive pipeline. This changes the framing, not the Model A/B decision itself (still blocked on item 4 for the real variance numbers), and touches:
  - Preliminary Results / Subtask 1.1 (lines 499–592, esp. 585–592), which currently frames ICEYE as *the* commercial SAR source — should be reframed as one member of a SAR suite.
  - The Model B likelihood term, which should be source-general (each SAR product's own measurement/sampling error) rather than ICEYE-specific, so the fusion framework can accept whichever SAR product is on hand.
  - The "transferable minimum-pathway protocol requiring no ICEYE purchase" line (278), which should be tightened to name the suite explicitly (Sentinel-1 free baseline + any available commercial SAR) instead of reading as "ICEYE or nothing."

### 2. IN PROGRESS (2026-07-22): collinearity admission undercuts the "duration adds explanatory power beyond depth" claim — text drafted, needs Napolitano sign-off
- Text admits duration and depth are collinear in the 271-building empirical sample (Subtask 1.4, ~line 691).
- The core physical-layer contribution rests on showing duration adds *statistically significant* explanatory power beyond depth alone (~line 717), tested on 271 buildings in one downtown core.
- Collinear predictors + geographically clustered (non-independent) buildings make this a harder identification problem than currently written. This is the load-bearing claim for the "first empirically derived duration-dependent fragility surface" outcome — a reviewer will ask how the study is powered to separate a collinear effect.
- **Done**: inserted a sentence (after ~line 719) explaining the identification strategy already implicit in the two-component design — $\lambda_0$/$\gamma$ are anchored to the physical component's material-science functional form, not fit solely from the correlated observational data, so a depth-only fit to correlated data can't masquerade as evidence duration doesn't matter.
- **Still open**: Becca wasn't sure the specific statistical commitments were right, so rather than inserting them as settled, left an inline `% becca/napolitano` comment asking her to confirm (a) whether a VIF check on the 2023 duration-depth relationship is the right diagnostic to report, and (b) whether the bivariate-vs-depth-only comparison should use cluster-robust standard errors given the 271 buildings are geographically clustered in one downtown core, not independently sampled. Needs Napolitano's sign-off before submission — she leads Subtask 1.4.
- **Related, not yet touched**: Subtask 1.5 (~line 738) separately claims the design "maximizes available power by using all 271 buildings as independent assessment units" for the SAR-validation coverage test — same non-independence issue, different test. Left alone; flag if the same clustering caveat should apply there too.

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

## Session 2026-07-22 (PI pass, afternoon): items 1, 2, 3 resolved in the narrative

### 1. RESOLVED as interval-censored Model B, approved by Becca 2026-07-22
- Checked the CPS-CIR track description in the NSF 25-543 program webinar slides before deciding. Key realization: the fusion the CIR track evaluates is mass-balance prior with community priors, not prior with ICEYE. ICEYE is the arbiter, so the old worry that Model A "kills the fusion contribution" was the wrong lens; the model choice is purely about scientific defensibility, and Model A (ICEYE exact) fails that regardless.
- Implemented: SAR duration is now an interval-censored observation. Each overpass gives wet/dry only, so true duration lies in a bracket [L_i, U_i] between last-wet and first-dry overpass. Likelihood is the prior mass on the bracket; posterior is the prior truncated to it. Mass-balance error moved to the prior side as sigma_MB (replaces both old sigma^2 and sigma^2_error), pre-estimated from gauge variability, calibrated on the 2023 brackets by censored maximum likelihood.
- This formulation is cadence-agnostic, so the equations no longer wait on item 4. Christelle's scene list now only determines the bracket widths (and how strong the preliminary-results claims can be), not the model structure. The "6-hour temporal resolution" claims were hedged to "as fine as six hours where the acquisition sequence is dense" pending her reply.
- Downstream updates: Subtask 1.1 outputs are brackets plus per-scene classification error; 1.5 and 2.3 metrics moved to interval-censored log predictive score and bracket-consistency coverage (also resolves the old CRPS comment); KPIs reworded to match; Eq (6) conditioning cleaned up; Subtask 2.4 interval-width sentence now points at sigma_MB.

### 2. RESOLVED, Napolitano sign-off given 2026-07-22
- Becca approved both: VIF diagnostic on the 2023 duration-depth data, and cluster-robust standard errors at the drainage sub-basin level for the bivariate-vs-depth-only comparison. Written into Subtask 1.4; inline comment removed.
- The related Subtask 1.5 "271 independent assessment units" claim was fixed the same way: buildings are now assessment units with cluster-robust uncertainty at sub-basin level, consistent with 1.4.

### 3. RESOLVED via mechanism attribution and an explicit test hierarchy
- Subtask 2.1: each elicited adjustment is now recorded with its stated physical mechanism, not just magnitude.
- Subtask 2.3: the 2023 comparison is now explicitly labeled a diagnostic of elicitation mechanics (partners saw the 2023 map, so it cannot separate knowledge from recall); the 2024 held-out event carries the substantive RQ2 test via cross-event persistence, evaluated separately for mechanism-attributed vs unexplained adjustments. Also added a detectability statement: adjustments smaller than the bracket width for a building cannot be adjudicated by that event's cadence, and the reported characterization says so per sub-group.

### Citations added (2026-07-22, verified batch 1)
- Root `biblio.bib` created (copy of latex/biblio.bib plus 23 verified entries) and is now canonical; root Research_v2.tex resolves against it. All entries verified against publisher/agency records: HEC-RAS/LISFLOOD-FP, HAZUS 7.0, FEMA P-58, JRC depth-damage, Gerl 2016, Kelman & Spence 2004, Schwarz & Maiwald 2008, Milanesi 2018, Hall & Hoff, Franzoni 2015, Sathiparan 2018, StreamStats, TR-55, NHDPlus, O'Hagan 2006, Garthwaite 2005, Sendai, Twigg 2015, Goodchild 2007, Assumpcao 2018, Kim 2019 (AMC, left as a pointer for Maggie).
- Still pending from second research batch: NFIP/ACS/FEMA BCA/NIST CRPG/FEMA HMA placeholders, the FEMA P-2055 vs P-2078 question, and the NSF SCC papers (trust, two-tier outputs, closest-to-problem knowledge).

### Citations completed (2026-07-22, verified batch 2) and two catches
- All 7 placeholder citations now resolved with verified sources: NFIP count (CRS IF10988), pre-1940 stock (ACS B25034 plus the VHFA Vermont Housing Needs Assessment, which gives 25.7 percent of Vermont stock pre-1939), FEMA BCA Toolkit, NIST SP 1190 vols I and II, FEMA HMA Guide 2023, FEMA P-312 and P-259, and small-municipality capacity literature (Smith 2013, Horney 2012, Smith 2026).
- **Catch 1: FEMA P-2078 is a seismic response-spectra document**, not a flood mitigation handbook. Becca's inline question about citing it is resolved: dropped, replaced with P-312 and P-259.
- **Catch 2: the "NSF SCC-funded projects have documented technology fatigue" claim had no citable SCC paper behind it.** The Purdue nutrient SCC award has only agronomy outputs, the Maine STEM project is AISL-funded, the Champaign gun-violence SCC project has no peer-reviewed publication yet, and the on-topic Gardezi trust paper is FW-HTF-funded. Fixed honestly: the trust sentence now says "community-engaged flood technology research" and cites the Lafayette Parish SCC project (#2125472, Skilton 2022 and Habib 2023). The two-tier claim IS cleanly SCC-supported (Habib 2023 explicitly recommends two classes of tools) and keeps its SCC framing. The closest-to-problem claim cites Skilton 2022 plus the same team's 2026 EPSCoR follow-on, with the funding distinction noted in an inline comment.

### Mass balance preliminary results added (2026-07-22), plus a factual correction
- New subsection "Preliminary results: Winooski watershed drainage timescale" with real numbers from USGS instantaneous records: July 2023 peak 23,100 cfs (stage 21.3 ft), recession timescale tau = 37 h (R2 0.91); July 2024 peak 11,900 cfs, tau = 45 h (R2 0.96). Figure at images/winooski_recession.pdf; reproducible script and full results in analysis/.
- **Factual correction: the proposal's gauge number was wrong.** 04288000 is the Mad River near Moretown; the Winooski River at Montpelier gauge is 04286000 (verified against the NWIS site record). Fixed in Subtask 1.2 with an inline note.
- **Flag for Becca (not a text change): the 2024 peak stage (14.45 ft) never reached the gauge's NWS minor flood stage (15 ft).** The narrative says Montpelier "sustained significant flood damage again in 2024." The prelim subsection words this carefully (low-lying parcels took water while the gauge stayed below minor flood stage), but worth confirming the 2024 damage narrative against the field data before submission, since a reviewer can pull the same gauge record.
- Data caveat: USGS discharge values are provisional; both recessions were interrupted by secondary rain pulses about 5 days post-peak, so fits use the first 60 hours.

### Deliberately left per Becca (2026-07-22)
- "Results from Prior NSF Support" stub for Napolitano: Becca will handle.
- "xxx PhD students" placeholders: budget not settled, left as placeholders.

## Next (as of 2026-07-22 PI pass; supersedes the morning list)
- **Item 4**: still awaiting Christelle/Young's reply, but it no longer blocks the equations (interval-censored model is cadence-agnostic). Her file list determines bracket widths and how firm the preliminary-results cadence language can be.
- **Items 1, 2, 3, 5**: resolved and worked into the narrative (see session notes above).
- **Remaining drafting gaps**: second citation batch (federal/SCC), Winooski mass-balance preliminary results subsection, field damage campaign description (needs Becca's 2023 campaign details), Prior NSF Support for Napolitano (Becca), PhD student counts (Becca), figures (CPS diagram, Montpelier photo, Gantt, ICEYE panels from company).

**Scoping note (2026-07-22):** `latex/DMP2.tex`, `latex/mentoringPlan.tex`, and `latex/summary.tex` (plus the stale `latex/Research.tex` and `latex/Research_v2.tex`) still describe the older SAR-coherence-time-series architecture and haven't been reconciled to the canonical mass-balance-plus-ICEYE design in root `Research_v2.tex`. Becca confirmed this is intentional — those get reconciled in a later pass. Current focus is narrative (root `Research_v2.tex`) only; don't touch the companion docs until asked.

### Done (2026-07-22): sensing-layer text reframed from ICEYE-only to a SAR suite
Per item 1's update above, reworded the descriptive/framing text throughout `Research_v2.tex` from "ICEYE" to "a SAR suite (ICEYE commercial and Sentinel-1 open SAR, extensible to any available source)": goals box (68), CPS architecture figure caption (92), RO1 summary table row (229–233, 263), minimum-pathway language (278, 1020, 1056), sensing-layer footnote (283–286), Preliminary Results subsection heading (499), RO1 overview paragraph (573–581), Subtask 1.1 title and body (585–593), Key outcomes (759), and the transferable-protocol interval-width sentence (1004). Also resolved two of Christelle's inline TODO comments (former lines 592 and 1054) that were asking this exact question.
**Deliberately left untouched**: the Subtask 1.3/1.4 equations and notation (`T_i^{\mathrm{ICEYE}}$`, lines ~639–729) and the specific 271-building empirical descriptions (245–246, 503–505, 686, 700, 813, 834, 950) tied to the actual NDA'd ICEYE dataset in hand. Whether Sentinel-1 gets fused into that *same* likelihood term (a joint multi-source error model) or is used only descriptively/for broader flood-extent context (ICEYE remains the sole building-level duration/depth source, since Sentinel-1's ~10m resolution may not resolve individual downtown buildings) is a modeling decision for Wauthier, not something to guess at in the text — flagging as a new open item.

### 5. RESOLVED (2026-07-22): Sentinel-1 stays descriptive-only, does not join the Subtask 1.3 likelihood
- Decision: Becca's call, and independently confirmed by Christelle the same day — "Sentinel-1 images are available, but the spatial resolution is not high enough to examine the details. They could still provide some rough information." The ICEYE (<1m) vs. Sentinel-1 (~10m) resolution mismatch is too large to fuse into one building-level likelihood; downtown-core buildings are frequently smaller than or comparable to a single Sentinel-1 pixel, so a joint/multi-source error model would require either down-weighting Sentinel-1 to near-irrelevance or aggregating to block level, a different physical quantity than the per-building $\hat{\tau}$ the model targets. Not worth inventing that methodology.
- This also directly answers the open TODO Becca left at (former) line 1021 about whether Sentinel-1 is useful for community workshop flood mapping ("i thought the resolution wasn't good enough") — confirmed: usable for rough/qualitative context, not building-level detail.
- **Text updated accordingly** (reverting most of the 2026-07-22 "SAR suite" reframing back to ICEYE-specific, since only ICEYE feeds the quantitative building-level likelihood/fragility calibration): goals box (68), architecture figure caption (92–93), RO1 table row (229–234, 264), sensing-layer footnote (283–287), RO1 overview (575–586), Subtask 1.1 body (588–594), Key outcomes (761), transferable-protocol interval-width sentence (1006). Sentinel-1 is now described consistently as free broader flood-extent context / cross-check on ICEYE building classifications, not a likelihood input.
- **Left unchanged**: the "no commercial SAR purchase" language for new communities (279, transferable-protocol section) — that claim was never about Sentinel-1 joining the likelihood, it's about new communities inheriting Montpelier-fit parameters rather than needing to buy their own calibration data of any kind, so it holds regardless of this resolution.
- **Practical-use implication (flagged for awareness, not a text change)**: this means a future community with only free Sentinel-1 and no ICEYE-grade SAR does not get a cheaper path to a tighter local calibration under the current design — they still inherit the wider, Montpelier-derived default credible intervals, same as a community with zero SAR access. The minimum-pathway paragraph already implies this but doesn't say it outright; worth a sentence if a reviewer is expected to ask.

## Session 2026-07-22 (evening): adversarial review worked through as PI

All ten concerns plus the nits in `adversarialReview.md` are now addressed in the narrative except where blocked on external data. Becca approved the three decision points before work started: distinguish standing-water from river-stage duration (keep the 5-day condition, redefine it on standing water), add a stage-to-DEM mapping as the operational depth source, and hedge the 2024 KPI rather than wait for field data.

### Resolved in the narrative
- **Review 1 (depth source)**: new explicit cyber-layer component in Subtask 1.2: gauge crest stage (observed / scenario / NWS forecast) projected onto 3DEP terrain along a reach-scale water-surface slope, giving building-level depth d_i and inundation extent in every operational mode; the 2023 ICEYE depth field calibrates the slope. Eq (6) and pre-storm mode now state where d_i comes from.
- **Review 2 (extent vs duration)**: scoping sentence restored (T_i conditional on inundation; error model characterizes duration error, not flood/no-flood, which the stage projection carries). Sub-basins defined: ten to twenty across the flood-exposed core, tens of buildings each, lumping deliberate.
- **Review 3 (scope condition)**: 5-or-more-day condition explicitly redefined on building-level standing water. Backed by new stage-threshold analysis (analysis/stage_thresholds.py, stage_durations.md): 2023 gauge above minor flood stage 28.2 h, above action stage 47.3 h; and three verified news citations showing downtown building interiors pumped for 4-6 days (vermontpublic2023stateflood, vtdigger2023schools, vermontpublic2023recede, added to biblio.bib with verified quotes).
- **Review 4 (2024 KPI power)**: KPI and Subtask 1.5 now claim significance under 2023 cross-validation, with the 2024 held-out comparison stated relative to the damage variation that smaller event produced. Still on the risk register pending the 2024 field damage distribution.
- **Review 5 (sigma_MB identifiability)**: gauge-based pre-estimate now described as a regularizing prior; wide brackets pull the estimate toward it instead of leaving it unidentified. Factual risk still open pending Christelle's scene list (item 4).
- **Review 6**: fragility fitting now integrates over interval-censored duration brackets (Monte Carlo over the truncated posterior), consistent with 1.3.
- **Review 7**: lambda_0/gamma now carry ds subscripts throughout; worked sentence added showing how saturation-strength curves anchor lambda_{0,ds}, gamma_{ds}, zeta_{ds}.
- **Review 8**: elicitation KPI restated: sessions cover all 271 buildings, nonzero mechanism-attributed adjustments for at least 20% of buildings, zero defaults recorded as explicit no-knowledge outcomes (20% is a PI judgment call, Becca can adjust).
- **Review 9**: new paragraph in the research plan naming the three loop timescales (hours / per-event / years) as deliberate human-in-the-loop design and noting the RQs emerged from the post-2023 Montpelier engagement.
- **Review 10**: audited all four "first" claim sites; every one carries the pre-code load-bearing masonry qualifiers. No change needed.
- **Nits**: goals box now duration-only priors (plus the depth projection named separately); right-censoring sentence added; stale "log-scale likelihood" wording fixed; wavelenght/in details/a NDA typos fixed; SAR defined once; "comparable probabilistic basis" fixed.

### New preliminary results (both reproducible, uv run per script headers)
- **Stage-threshold durations** (analysis/stage_thresholds.py): 2023 47.3 h above action stage down to 22.4 h above major; 2024 14.5 h above action, never reached minor. Secondary rain pulse never re-crossed any threshold, so peak-interval = total duration.
- **Stage-to-DEM prototype** (analysis/stage_to_dem.py, stage_dem.md, images/stage_dem_depth.pdf): gauge datum 499.87 ft NAVD88 (live-queried); 2023 crest flat-WSE projection gives ~92 connected acres tracing the river corridor but underpredicts the documented downtown core by 0.8-1.5 m at Main Street landmarks, consistent with a ~0.6-1.0 m/km water-surface slope. Written into the narrative as motivation for the ICEYE-calibrated slope parameter, not hidden. Figure exists but left out of the tex for space (inline % becca comment marks it).

### Still blocked / for Becca
- Review 4 residual: 2024 field damage-state distribution (Becca's campaign data) before the KPI hedge can be firmed either way.
- Review 5 residual / item 4: Christelle/Young scene list still pending; determines bracket widths and prelim cadence language.
- 20% elicitation KPI rate and the ten-to-twenty sub-basin count are PI-judgment numbers inserted today; Becca should sanity-check both (sub-basin count is also flagged to Maggie implicitly via her Subtask 1.2 comment block).
- Known-incomplete list unchanged: field damage subsection, Prior NSF stub (Napolitano), PhD counts, CPS/Gantt/photo figures, RMSE threshold sensitivity analysis.

## Session 2026-07-22 (late): sensitivity studies on the two PI-judgment numbers

Becca asked for pythonic sensitivity studies on the two judgment calls inserted during the adversarial-review pass. Both ran as seeded, reproducible scripts in analysis/ (run commands in each script header; note `uv run --with ...` is required, system pythons lack packages, and subbasin needs numpy<2 for pysheds).

### Sub-basin count (analysis/subbasin_sensitivity.py, .md, .png)
- D8 delineation on the 3DEP DEM swept from K=3 to K=106 sub-basins; 176 OSM-proxy flood-exposed buildings (order-of-magnitude proxy for the curated 271).
- Raw groupwise R2 is mechanically inflated as K grows, so the study uses a permutation-null-corrected gap (verified against the exact analytic E[R2]=(k-1)/(n-1) to within 0.004).
- Verdict: **ten-to-twenty holds**; the null-corrected gap is flat (~0.65-0.69) across K=3-22 and declines steadily to 0.32 by K=106, so the claimed range sits inside the information plateau and finer partitioning adds parameters, not signal.
- Correction adopted in the tex: "each containing tens of buildings" was wrong; counts are strongly skewed (one ~60-building main-corridor sub-basin, many single-digit tributary basins). Subtask 1.2 now says so and cites the archived sensitivity analysis.
- Caveats: terrain-depth proxy (not a drainage-timescale simulation); OSM building set; 2m DEM (1m returned HTTP 500 at the enlarged bbox); permutation null is non-contiguous, so the gap partly reflects generic spatial autocorrelation (does not affect the granularity claim).

### Elicitation KPI rate (analysis/elicitation_power.py, .md, .png, 3 csv grids)
- Monte Carlo power study of the Subtask 2.3 community-vs-GIS comparison: 271 buildings, 15 skewed clusters, interval-censored scoring, cluster bootstrap; sweeps p_elicit x bracket width x sigma_MB x elicitation quality q.
- Verdict: **20% is defensible but not lavish**: ~70% power centrally (62-77% across cadence); 30% clears 80%. The aggregate test used is stricter than the proposal's literal "at least one sub-group" KPI, so these are conservative.
- The vacuity check validates the KPI restatement: at q<=0.5 (elicited volume mostly wrong), power sits at the false-positive floor at every rate up to 50%, and the community prior is on average worse than GIS-only. Mechanism attribution is load-bearing, exactly as the restated KPI claims.
- Structural fragilities worth knowing: power collapses when p_elicit approaches 2x the true anomaly prevalence (correct and spurious calls balance out), and sparse cadence (weekly brackets) or sigma_MB=0.8 pulls 20% down to 43-62%.
- Tex updated: Subtask 2.3 now cites the archived power study (70-80% power at 20-30% mechanism-attributed rates, cadence dependence, volume-without-quality contributes nothing).

### For Becca
- If she wants more headroom on the KPI, moving 20% to 25-30% buys real power under the central assumptions; left at 20% since the sub-group KPI it feeds is less strict than the simulated aggregate test.
- The one-dominant-sub-basin skew also means cluster-robust inference at sub-basin level leans on ~15 clusters with one heavy one; consistent with the 1.4/1.5 text as written, but worth a look from Napolitano's stats side before submission.

## Session 2026-07-22 (night): equations pass + four more sensitivity studies

### Equations added (CPS-community formalization pass)
All seven existing equations converted to \label/\eqref (references were hardcoded numbers and would have silently broken). Five new numbered equations: linear-reservoir tau_k definition (eq:tau, flagged % maggie as a placeholder anchor consistent with the exponential recessions), stage-to-DEM depth projection (eq:stagedem, d_i = max{0, s_g + beta x_i - z_i}), classification-error-softened bracket likelihood (eq:softbracket, recovers the hard bracket as eps -> 0), interval-censored log predictive score (eq:lps), and the precision-weighted per-event delta update (eq:deltaupdate, the actuation law of the per-event loop).

### Four sensitivity studies (all in analysis/, seeded/reproducible, run commands in script headers; agents hit the 6pm session limit mid-batch and were resumed after reset)
- **historical_tau**: 15 well-fit floods, 1990-2026 IV record (IV goes back to 1990 at this gauge, not 2007). tau 36-96 h, median 51, sd(ln tau) = 0.30 = the sigma_MB pre-estimate, now cited in Subtask 1.2. Antecedent: wetter pre-event baseflow correlates with SLOWER recession (r=0.59, ~+33 h dry-to-wet), direction opposite the AMC worry as stated; pre-event baseflow is an operational AMC proxy from the same gauge. Written into prelim results; % maggie AMC comment annotated. 2023/2024 reproduce exactly.
- **cadence_identifiability**: sigma_MB stays identified at every cadence tested including a single usable post-peak scene (estimation uncertainty grows continuously, no collapse, 0% bound-pinning in the pessimistic scenario); regularizer narrows every cell at a quantified bias cost. One honest non-monotonicity: the "requested 2023" schedule at sigma=0.3 has genuine boundary-MLE degeneracy because scenes clustered at 0-30 h never straddle multi-day drying times. That insight (spacing beats count) is now a design sentence in Subtask 1.1. De-risks adversarial review item 5 ahead of Christelle's scene list.
- **ranking_sensitivity**: top-20 retrofit overlap stays >= 0.8 up to ~0.40-0.44 log-units duration error (70-78 h portfolio-wide RMSE) across fragility variants; Spearman >= 0.9 to ~0.5. The decision-adequacy KPI now carries this number (stricter overlap criterion) and Becca's inline "is it worth it?" comment is resolved. Caveats: synthetic portfolio, duration channel only; redo on real 2023 ICEYE durations when in hand.
- **slope_sensitivity**: admissible reach slope bracketed at >= 0.9 m/km from documented 2023 outcomes (City Hall + Confluence Park flood, State House dry across sweep; State House only floods near 8 m/km so evidence sets a floor, not a ceiling); 234 OSM buildings inundated at 0.9 m/km, same order as the curated 271. Consistent with the 0.6-1.0 m/km flat-residual inference; written into prelim results.

### For Becca / open
- Wild cluster bootstrap committed in Subtask 1.4 text (consequence of the dominant-cluster skew), inline comment asks Napolitano to confirm.
- eq:tau linear-reservoir form and the AMC formulation are still Maggie's calls; both flagged inline with the new empirical numbers attached.
- Unchanged blockers: 2024 damage distribution, Christelle scene list, Prior NSF stub, PhD counts, figures.

## Session 2026-07-23: Research_v3.tex branched, title decided, FY2028-memo framing pass

Branched `Research_v3.tex` from the frozen `Research_v2.tex` (commit ff6895b) as the new working file. Also resolved an Overleaf/repo sync scare: the freshly-downloaded Overleaf copy (`Research_vold.tex`, untracked) turned out to match the repo exactly at commit 9cb7b76, the point right before the 2026-07-22 editing session -- Overleaf was simply never updated after that session, nothing was lost or overwritten. Needs pushing back up to Overleaf when convenient.

Read `Science-A-New-Golden-Age.md` (OSTP/OMB FY2028 R&D priorities memo) and layered in three small, page-neutral framing edits tying the existing methodology to named memo priorities (not adding new claims, just naming what the framework already does): "Bayesian uncertainty-quantification framework" in the goals box (line 69), "state-estimation problem" opening Subtask 1.3 (line 808), and "dynamics, estimation, and control" reframing the closed-loop-timescales paragraph (line 658). Explicitly steered away from claiming "AI" anywhere -- the memo itself warns reviewers to screen out AI-labeling without justification, and this framework is Bayesian statistics, not ML.

**Title decided**: CPS-CIR: Closing the Cyber-Physical Loop between Community Knowledge and Hydrological Modeling for Flood-Resilient Historic Communities. Drafted in the style of funded S&CC/CIR titles (checked via web search against real awards, including #2125472, the Lafayette Parish project already cited elsewhere in this tex for its two-tier design and tech-skepticism findings). Recorded as a comment at the top of Research_v3.tex; needs entering on the actual Research.gov submission form since it isn't part of the typeset Project Description.

### Still open
- Page limit: body is 22 pages against the 15-page target, unaddressed since the last tightening pass.
- Prior-NSF stub for PI Napolitano still empty.
- Placeholder figures: CPS architecture diagram, Gantt chart.
- Blocked on others: 2024 field damage distribution, Christelle/Young ICEYE scene list.
- latex/summary.tex is stale relative to Research_v3.tex (still describes a two-site Montpelier+Marshall design and Capella X-band SAR that no longer match the current framework) -- needs a rewrite pass before submission, not touched this session.

## Session 2026-07-23 (afternoon): math adversarial review fixed end to end

All 8 items plus the sigma_MB-overloading umbrella in `mathAdversarialReview.md` are resolved in Research_v3.tex, verified by an independent adversarial re-review whose four residuals were also fixed. Full per-item log now lives at the bottom of `mathAdversarialReview.md`. Becca's two modeling calls going in: random-walk Kalman for the per-event delta update (tracked P_i, process noise q^2), ordered probit for the fragility surface (shared gamma/zeta, ordered lambda_0,ds).

Three new reproducible studies in analysis/ (uv run commands in headers): `delta_filter` (validates the Kalman claims; forced honest softening of "P shrinks" and of q-from-two-events), `regularizer_misscale` (sigma_event anchor mis-scaling moves calibrated sigma_MB <0.01 at realistic cadence), `depth_uncertainty` (depth error ~4% of damage-probability uncertainty at 2023 cadence but dominant at dense 6h cadence; 7-18% of buildings flip flood/no-flood across the beta range — justified propagating depth everywhere).

Also: cameron2008bootstrap added to root biblio.bib (verified); historical_tau.md sigma_MB conflation corrected; both % maggie eq:tau comments updated to reflect the committed fractional-storage closure (still her call).

**New for Becca / watch items**
- PDF is now 28 pages (was 26); the fixes cost ~1.5 body pages. The 15-page cut pass is still pending and now has to absorb this too.
- q (Kalman step variance) is committed as a swept design hyperparameter, only coarsely bounded by a 2023-2024 pair; if Napolitano prefers a fixed q, the delta_filter study shows mismatch is benign (no cliff, always at/below fixed-gain RMSE).
- Blockers unchanged: Christelle/Young scene list, 2024 field damage distribution, Prior NSF stub, PhD counts, CPS/Gantt figures, summary.tex rewrite.
