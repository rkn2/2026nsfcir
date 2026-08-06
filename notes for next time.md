# CIR (NSF CPS-CIR) — resume notes

**Deadline: 2026-09-01 (hard).** Becca is PI + primary writer.
Working file: `Research_v3.tex`. Compile with `~/.local/bin/tectonic Research_v3.tex` (no brew/sudo/Docker on this Mac).

## Current status (2026-08-06)
- Body = **18 pages** against the **15-page** Project Description limit → **3 pages over**.
  (Gantt chart added ~1 page; body was at 17 before adding it.)
- References start on p19. Compiles clean.
- **Full framing sweep is DONE** — §1 through §7 all checked against the closed-loop framing,
  as of commit `d733989`.

## PICK UP HERE NEXT TIME

### 1. Cut for length (~3 pp needed)
This is the main remaining task. Cut levers ranked by size:
1. **§2.2 worked dollar example** (~0.5–0.75 pp) — biggest self-contained object, low risk.
2. **Gantt chart** — could be compressed (smaller font, fewer labels) or cut if space is tight.
3. **Subtask description trims** — judgment pass across 1.2, 1.3, 2.1, 2.2 (longest subtasks).
4. **One of the two preliminary-results figures** (recession analysis) (~0.3–0.4 pp).
5. **Tier 3 compressions** (reviewer-facing defenses/honesty caveats) (~0.4 pp).

### 2. Remaining TODOs from Becca
- **Damage state counts** — how many of the 271 buildings are no damage / cosmetic / structural?
- **BESURE details** — how many undergrads per year, recruitment from underrepresented groups, REU supplement?
- **Postdoc and PhD advising** — who advises each, any co-advising?
- **Expert panel** — "APT DRI" for decision-adequacy validation in Subtask 1.5; confirm full name.
- **Transferable community count** — quantify using NFIP enrollment, ACS pre-1940 housing, USGS gauge
  coverage (inputs are in biblio.bib; cross-tabulation still needed).
- **`Research_vold.tex`** — untrack or leave? Not resolved.

### 3. Co-PI items (visible in PDF as yellow highlights)
**Christelle (Wauthier):**
- Confirm ICEYE figures are OK to use or provide replacements.
- Confirm 2023/2024 scene cadence.
- Add ICEYE processing pipeline description and depth/duration measurement flags (Subtask 1.1).

**Maggie (Busse):**
- Confirm or replace τ̂_k formulation (linear-reservoir, fractional-storage closure).
- Provide sub-basin discretization (approximate size/count).
- Address antecedent moisture condition (AMC) handling for pre-storm mode.
- Confirm whether this uses an existing Winooski model or requires new development.
- Add prior NSF support if any beyond CBET-2400672.

### 4. After cuts
- Read-through / hypothesis check (outstanding since 2026-06-18).
- Push to collaborators (Wauthier, Busse) for review — the yellow highlights mark their items.
- Fill/confirm any remaining Prior NSF Support text.

### Framing reference paragraph
The intellectual-merit paragraph (RQ subsection, ~line 177) is the **canonical statement of the
framing**. It says: the closed loop is THE contribution; the three advances (computational /
physical / community knowledge) are subordinate to it. Masonry appears only as two-sided
physical coupling. Anything written later should match it.

## Structural changes made this session (2026-08-06)

### Framing sweep (full document)
- **§2 opener** — short capability-chain setup, no longer redundant with gap summary.
- **§2.3 community knowledge gap** — added tacit-process-knowledge framing.
- **§3 recession analysis** — cut redundant baseflow paragraph, compressed stage-to-DEM.
- **§3 field damage assessment** — three ordinal states (no damage / cosmetic / structural),
  pre-2024 timing confirmed, each building paired with ICEYE measurements.
- **§4 research plan intro** — bridge sentence; three-timescale closing as one loop at three speeds.
- **RO1 intro** — rewritten around four capabilities with community knowledge handoff to RO2.
- **RO2 box** — added substitute-versus-complement language to match RQ2.
- **Subtask 2.4/2.5 split** — 2.4 is Closed-Loop Validation, 2.5 is Transferable Protocol.
- **Subtask 1.2** — retitled "Mass Balance and Depth Projection" with framing sentence.
- **Subtask 1.4** — material-science prose trimmed; wild cluster bootstrap → Bayesian model comparison;
  three-state ordinal taxonomy integrated; building archetypes added.
- **RO1 Key Outcomes** — community knowledge folded into outcome 1.
- **RO2 Subtask 2.5** — substitute/complement connection to 2.3 findings added.
- **§5 broader impacts** — first-ness claim killed; tacit-knowledge framing threaded through
  §2.3, §5 intro, and §5.1 (Golden Age Ch.IV parallel, absorbed not cited).
- **§5.2 Education** — 1 postdoc + 1 PhD, BESURE sentence completed.
- **§6 timeline** — Gantt chart added (pgfgantt). Stale cluster-bootstrap KPI fixed.
- **§7 prior support** — co-PI comment formatting cleaned up.
- **Co-PI comments** — all Christelle/Maggie comments ALL CAPS with yellow PDF highlights.
- **Lara citation** added (Maalouf & Napolitano 2025, IJHS).
- **README.md** created explaining file organization.
