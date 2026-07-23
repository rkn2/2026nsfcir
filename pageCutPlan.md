# Page-cut plan for Research_v3.tex (2026-07-23)

Body is ~23.3 pages against the 15-page NSF limit, and ~0.5 page must be reserved for Napolitano's Prior NSF Support text (currently a stub), so the real cut target is ~8.5-8.8 pages (~360 rendered lines at ~42 lines/page). Three scouting passes (Opus on the Research plan, Sonnet on front matter and prelim/back matter) produced verbatim-anchored candidates; each scout reported honestly rather than padding to target. Full candidate lists with exact anchors and replacement text live in the session transcripts; this file records the plan of record.

## Committed cut set (safe, ~6.2 pages)

| Territory | Content | Est. lines | Est. pages |
|---|---|---|---|
| Research plan Tier 1 (pure redundancy: RO2 prose mini-TOC, subtask enumerations, mode definitions restated, scoring metrics restated from 1.5, tier descriptions, parallel-run narration) | 12 cuts | ~62 | 1.5 |
| Research plan Tier 2 (narrative compression: three-timescale framing, workshop year narration, MC prose after eq:loss, null-findings philosophy, transfer-interval prose, Key outcomes blocks) | 16 cuts | ~78 | 1.85 |
| Front matter (four-time repeated "no framework connects drainage to buildings" claim, uncited CPS-literature paragraph, per-subsection gap restatements that Gap summary re-restates, federal-frameworks list compression) | 16 cuts | ~37 | 0.9 |
| Prelim/back prose (partner-site retelling of Status quo, SAR background, Winooski narrative connective tissue with all numbers kept, Broader impacts boilerplate, KPI compression with all thresholds kept, triplicated "gauge ID + GIS + one workshop" list) | ~31 cuts | ~59 | 1.4 |
| Figures: drop the fig:precode placeholder wrapfigure (img.jpg #2, "pre-code masonry buildings" photo slot) | 1 | ~13 | 0.3 |
| Figures: drop the fig:timeline Gantt figure (img.jpg #3, placeholder); the compressed timeline prose carries the schedule; NSF does not require a Gantt figure | 1 | ~15 | 0.35 |
| **Subtotal** | | **~264** | **~6.3** |

Protections honored throughout: all numbered equations, all % co-PI comments, every archived-study number (0.30, 37/45 h, 70-80% power, ~4% depth share, 0.4 log-units/70 h RMSE, stage-threshold hours), the sigma_event/sigma_MB distinction, ordered-probit no-crossing statement, Kalman recursion, the 2024 gauge hedge, verified citations, KPI thresholds.

## Decision levers (Becca's call, ~1.5-2.1 pages available)

1. **Research-plan Tier 3 (~0.65 pp):** compressions that touch reviewer-facing material: the single-reach-beta defense, the graceful-degradation study sentence, the depth-share block, the independent-posteriors caveat, the two Subtask 1.5 honesty caveats. Compress-only variants exist for each; protected clauses stay verbatim.
2. **Table 1 (~0.35-0.7 pp):** the gaps/objectives/approaches/outcomes summary table spans ~1 page and is by design redundant with the prose. Options: compress the Gaps column cells (~0.35 pp) or restructure to a tighter 3-column table (~0.7 pp).
3. **vt.pdf photo (~0.12 pp):** the only real event photograph (Montpelier during the 2023 flood). Human-interest value for panelists; small saving if dropped.
4. **Goals box + contributions block (~0.3-0.4 pp):** the opening tcolorbox (lines 62-89) and the "three scientific contributions" paragraph overlap Table 1's columns; both could compress.

## Decisions (Becca, 2026-07-23)

- Tier 3: **skipped entirely**; all reviewer-facing defenses and honesty caveats stay at full length.
- Table 1: **compress cells only** (~0.35 pp), keep the 4-column structure.
- Goals box + contributions block: **compress** (~0.3-0.4 pp).
- vt.pdf flood photo: **keep**.
- Expected landing point with these choices: ~16-16.5 pages; the last stretch to 15.0 is a Becca judgment pass or a later revisit of the skipped levers.

## Arithmetic

Committed 6.3 + all four levers ~2.0 = ~8.3 pp, which reaches 15.0 pages only if everything is taken and estimates hold. Rendered-line estimates historically overpredict savings (page breaks are quantized), so expect to land at 15.5-16 and need a final manual pass. Execution order: apply committed set, recompile, measure, then apply approved levers largest-first until the measurement says 15.0, keeping a running log so any lever can be reverted.

## Execution protocol

Batches by territory (research plan, front, prelim/back), compile and page-check after each batch, adversarial spot-review of the compressed research-plan sections at the end (the compressions must not have broken the math-review fixes), then commit and push.

## Execution results (2026-07-23, same day)

Executed: committed safe set (all three territories, ~58 cuts applied by Sonnet executor agents from the scout lists, orchestrator-applied goals box, contributions paragraph, Table 1 cell compression, and both placeholder-figure removals). An Opus diff review of all ~60 edits found exactly one defect (the deleted "Three observations anchor the research plan." lead-in left the First/Second/Third enumeration cold); restored. All equations, cross-references, % comments, protected numbers, and math-review-hardened content verified intact.

**Measured outcome: body 23.3 -> 20.0 pages (references start p21).** The scouts' rendered-line estimates overpredicted by ~2x (estimated ~6.3 pp from the safe set, realized ~3.3 pp), consistent with the project's earlier lesson that page breaks are quantized and prose reflows absorb savings.

**Remaining gap: 5 pages (plus ~0.5 for the Prior NSF stub when filled).** Unused levers, in order of likely yield:
1. Tier 3 compressions (skipped by Becca this round): ~0.4 pp realistic.
2. Deeper content cuts only the PI can make: dropping or table-izing whole passages (candidate: the 2.2 worked dollar example, the stage-to-DEM prototype narrative, one of the two preliminary-results figures, further Table 1 shrink to 3 columns).
3. Layout: the summary table spans ~2 pages rendered; converting to a half-page compact table or moving detail to subtask headers is the single largest remaining object.
4. Accept a Becca judgment pass on scientific content (which subtask descriptions can lose their second paragraph).
