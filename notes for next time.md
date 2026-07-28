# CIR (NSF CPS-CIR) — resume notes

**Deadline: 2026-09-01 (hard).** Becca is PI + primary writer.
Working file: `Research_v3.tex`. Compile with `~/.local/bin/tectonic Research_v3.tex` (no brew/sudo/Docker on this Mac).

## Current status (2026-07-28)
- Body = **18 pages** against the **15-page** Project Description limit → **~3 pages still over**. (References start on p19; they do not count toward the 15.)
- Compiles clean (exit 0, no undefined refs).
- Length is the one big open item. Framing/argument work is in good shape.

## Done in this session (2026-07-28)
- **De-layered the CPS framing** (commit `9b53c2f`). Removed the rigid "four coupled layers" taxonomy (physical/sensing/cyber/decision) that was fighting the three-contribution merit argument. The solicitation does not ask for a layer taxonomy — it was self-imposed scaffolding.
  - Goals box: enumeration → loop-story prose (sense → compute → community closes the loop).
  - Merit subsection: three named layers → three plain contributions (**computational / physical / community knowledge**).
  - Removed the `fig:cps_arch` four-layer diagram (`img.jpg`).
  - Swept ~15 downstream "X layer" mentions (Table row labels + body shorthand) to plain phrasing. Left "GIS layers" / "national data layers" alone (different meaning) and the two `%` PI comments.
- **Dropped Table 1** (the gaps/objectives/approaches/outcomes summary, commit `dc6ea0b`). It was ~1 full page and fully redundant with the prose. Fixed all four `\ref{tab:research_summary}` cross-references. Saved ~1 page (19 → 18).

## Remaining cut levers (largest first, ~3 pp needed)
1. **§2.2 worked dollar example** — full worked allocation example; table-ize or trim (~0.5–0.75 pp). **Suggested next** — biggest self-contained object, low risk.
2. **One of the two preliminary-results figures** (recession analysis + one other) (~0.3–0.4 pp).
3. **Stage-to-DEM prototype narrative** — compress to equation + one sentence (~0.3 pp).
4. **Second-paragraph trims on subtask descriptions** — judgment pass, many small cuts (~0.5–1 pp).
5. **Tier 3 compressions** skipped last round (reviewer-facing defenses/honesty caveats) (~0.4 pp).

Note: page savings are quantized — rendered-line estimates historically overpredict ~2x because reflow absorbs cuts. Measure after each batch (`mdls -name kMDItemNumberOfPages Research_v3.pdf` + find where "References:" starts via pdftotext).

## Open decisions for Becca
- **`Research_vold.tex`** — a 1191-line backup got swept into git via `git add -A` in commit `9b53c2f`. Untrack it (`git rm --cached`) or leave it? (Not resolved.)
- **Which cut lever next** — recommended the §2.2 dollar example.

## After cuts
- Read-through / hypothesis check (was outstanding since 2026-06-18).
- Push to collaborators (Wauthier, Busse) for review.
- Fill/confirm any remaining Prior NSF Support text (Wauthier CAREER text is in; confirm Napolitano/Busse as needed).
