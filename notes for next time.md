# CIR (NSF CPS-CIR) — resume notes

**Deadline: 2026-09-01 (hard).** Becca is PI + primary writer.
Working file: `Research_v3.tex`. Compile with `~/.local/bin/tectonic Research_v3.tex` (no brew/sudo/Docker on this Mac).

## Current status (2026-08-06)
- Body = **18 pages** against the **15-page** Project Description limit → **~3 pages still over**.
  (Added synthetic damage map figure which cost ~1 page; body was at 17 before the figure.)
- References start on p19. Compiles clean.
- **Framing sweep through §2, §3, and §4 intro/RO1 subtasks is DONE** as of commit `0f98f5a`.
- **Still unswept: RO2 subtasks (2.1–2.5), §5 (broader impacts), §6 (timeline — mostly done),
  §7 (prior support).**

## PICK UP HERE NEXT TIME

### 1. RO2 subtask sweep (highest priority)
Sweep Subtasks 2.1–2.5 against the closed-loop framing. Two flags already planted:
- **Substitute-vs-complement** must flow through all subtasks, not just 2.3 where it currently
  lands (see `% FRAMING SWEEP TODO` comment above the RO2 subtasks).
- **Subtask 2.2, Year 2** (~line 921): "cross-references community-observed damage against 2023
  field damage states" — against the three-state ordinal taxonomy (no damage / cosmetic /
  structural), check whether this cross-reference is informative or near-trivial.

### 2. §5 broader impacts
**Known leftover:** ~line 1120, "This framework is the first to place community-scale drainage
investment and…" — a first-ness novelty claim the merit paragraph deliberately abandoned.

### 3. Remaining cut levers (~3 pp needed)
The damage map figure added ~1 page. Updated lever list:
1. **§2.2 worked dollar example** (~0.5–0.75 pp) — biggest self-contained object, low risk.
2. **Damage map figure** — consider sizing down or combining with ICEYE panel (~0.3–0.5 pp).
3. **One of the two preliminary-results figures** (recession analysis) (~0.3–0.4 pp).
4. **Subtask description trims** — judgment pass, many small cuts (~0.5–1 pp).
5. **Tier 3 compressions** (reviewer-facing defenses/honesty caveats) (~0.4 pp).
Stage-to-DEM prototype was already compressed in this session.

### 4. Damage state data (from Becca)
- Three ordinal states confirmed: no damage, cosmetic, structural.
- Synthetic map is a PLACEHOLDER (yellow-highlighted in PDF). Replace with actual field data
  before submission. Script: `analysis/synthetic_damage_map.py`.
- Still needed: actual counts per damage state for the 271-building pre-code masonry subset.
- Building archetypes: predominantly load-bearing brick URM, small number mixed brick-and-stone.

### 5. Co-PI items (visible in PDF as yellow highlights)
**Christelle (Wauthier):**
- Confirm ICEYE figures are OK to use or provide replacements.
- Confirm 2023/2024 scene cadence.
- Add ICEYE processing pipeline description and depth/duration measurement flags (Subtask 1.1).

**Maggie (Busse):**
- Confirm or replace τ̂_k formulation (linear-reservoir, fractional-storage closure).
- Provide sub-basin discretization (approximate size/count).
- Address antecedent moisture condition (AMC) handling for pre-storm mode.
- Confirm whether this uses an existing Winooski model or requires new development.

### Framing reference paragraph
The intellectual-merit paragraph (RQ subsection, ~line 177) is the **canonical statement of the
framing**. It says: the closed loop is THE contribution; the three advances (computational /
physical / community knowledge) are subordinate to it. Masonry appears only as two-sided
physical coupling. Anything written later should match it.

## Structural changes made this session (2026-08-06)
- **§2 opener rewritten** — now a short capability-chain setup ("Closing the loop requires four
  capabilities… Each has a mature literature. None addresses the interface with the next."). No
  longer redundant with the gap summary. Approved by Becca.
- **§3 recession analysis tightened** — cut redundant baseflow paragraph, compressed stage-to-DEM
  prototype from 11 to 4 lines.
- **Field damage assessment filled in** — three ordinal states (no damage / cosmetic / structural),
  pre-2024 timing confirmed, synthetic damage map figure added as placeholder.
- **§4 research plan intro** — bridge sentence connecting §2 gaps to research plan; three-timescale
  closing rewritten as one loop at three speeds ("at every timescale" not "all three loops").
- **RO1 intro paragraph** — rewritten around the four capabilities from §2, with explicit community
  knowledge handoff to RO2 via δ_i^comm.
- **RO2 box** — added substitute-versus-complement language to match RQ2 fully.
- **Subtask 2.4 split into 2.4 + 2.5** — 2.4 is Closed-Loop Validation (the reentry mechanism),
  2.5 is Transferable Protocol. RO2 now has five subtasks. All cross-references updated. Adoption
  KPIs updated to include 2.5.
- **Subtask 1.2 retitled** — "GIS-Parameterizable Mass Balance and Depth Projection" with framing
  sentence naming both components. Addresses reviewer concern that the stage-to-DEM model was
  hidden inside a "mass balance" subtask.
- **Subtask 1.4 trimmed** — material-science prose compressed so masonry reads as one mechanism
  serving the loop. Wild cluster bootstrap replaced with Bayesian model comparison (interval-censored
  log predictive score), consistent with the rest of the framework. Three-state ordinal damage
  taxonomy integrated into the fragility equation description. Building archetypes added.
- **RO1 Key Outcomes** — community knowledge slot folded into outcome 1, kept at three outputs.
- **Co-PI comments** — all Christelle/Maggie comments converted to ALL CAPS with visible
  yellow-highlighted notes in the rendered PDF.
- **Lara citation added** (Maalouf & Napolitano 2025, IJHS).
- **README.md created** explaining file organization.

## Open decisions for Becca
- **`Research_vold.tex`** — untrack it (`git rm --cached`) or leave it? (Not resolved.)
- **Damage map figure** — keep at current size or shrink to save space?
- **Expert panel for decision-adequacy validation** — Becca said "APT DRI"; confirm full name and
  composition for Subtask 1.5 text.

## After framing sweep + cuts
- Read-through / hypothesis check (outstanding since 2026-06-18).
- Push to collaborators (Wauthier, Busse) for review — the yellow highlights mark their items.
- Fill/confirm any remaining Prior NSF Support text (Wauthier CAREER text is in; confirm
  Napolitano/Busse as needed).
