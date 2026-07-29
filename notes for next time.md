# CIR (NSF CPS-CIR) — resume notes

**Deadline: 2026-09-01 (hard).** Becca is PI + primary writer.
Working file: `Research_v3.tex`. Compile with `~/.local/bin/tectonic Research_v3.tex` (no brew/sudo/Docker on this Mac).

## Current status (2026-07-28)
- Body = **18 pages** against the **15-page** Project Description limit → **~3 pages still over**. (References start on p19; they do not count toward the 15.)
- Compiles clean (exit 0, no undefined refs).
- **Framing/argument work is DONE** as of commit `d6f0ff8`. Length is the one remaining big item.

## PICK UP HERE NEXT TIME (2026-07-29)
**The framing sweep was NOT actually complete** — the 2026-07-28 sweep was wording-level and never
reached §2. §2 (Current state-of-the-art) has now been rethesised to the loop (see below).
**Still unswept against the closed-loop framing: §3 (partner site/preliminary results, ~348–500),
§4 (research plan, ~501–1080 — the big one), §5–7 (broader impacts, timeline, prior support).**

**Known leftover already spotted, fix when sweeping §5:** line ~1120, "This framework is the first
to place community-scale drainage investment and…" — a first-ness novelty claim of exactly the kind
the merit paragraph deliberately abandoned (same species as the RO1 fix in item 3 below).

**Verification method that works** (use it on §3/§4, a read-through will not catch this): build a
mapping of each merit-paragraph claim → which paragraph of the target section establishes its gap
or delivers it. Empty cells are the findings. Also grep for residue from *both* old framings:
`grep -n 'layer\|epistemic\|first \(to\|for\)\|four \(coupled\|layers\|domains\)\|co-equal\|novel'`.

Length is still ~3 pp over; cut levers below are untouched.

### Framing reference paragraph
The intellectual-merit paragraph (RQ subsection, ~line 177) is the **canonical statement of the
framing** — anything written later should match it, not drift from it. It says: the closed loop is
THE contribution; the three advances (computational / physical / community knowledge) are
subordinate to it, not co-equal novelties. Masonry appears only as two-sided physical coupling
(sensing flood signatures; drainage/hardening interventions whose effects reenter future flood
behavior) — the fragility-surface novelty clause was deliberately removed from this paragraph.

## Done in this session (2026-07-29)
- **Carried the closed-loop framing into §2 Current state-of-the-art.** The mapping check found the
  structural hole the wording sweep missed: the words *loop / closes / reenter* appeared **zero
  times** in §2, so the central contribution had no state-of-the-art support at all. §2's own thesis
  was still "integrating these four domains creates two scientific tests" — the old co-equal-novelty
  structure. Four fixes, **net length delta exactly 0** (body 18 pp before and after, References
  still at top of p19):
  1. **§2 opener rethesised** — "four domains → two tests" → "four domains each supply one piece of
     the loop; none couples them; the loop has never been evaluated end to end." Also fixes the
     dangling "these four domains" (the antecedent came *after* the sentence). Ran ~2 lines shorter,
     which paid for fix 3.
  2. **Gap summary rethesised** — was landing on two co-equal gaps (no model + community knowledge
     untested); now lands on "the pieces exist; the loop does not," enumerates where each domain
     stops at the interface with the next, and points forward into the research plan.
  3. **Post-event updating gap added** (§2.4) — the mode that actually *closes* the loop had no gap
     anywhere. Careful wording needed: post-event *planning* does exist, so the claim is that
     practice doesn't *feed back* (assessments reissued from static inputs, nothing improves as a
     community accumulates events).
  4. **Two-sided physical coupling labeled** (§2.4) — the intervention side was present but
     unlabeled; added the reentry clause ("nothing carries a completed retrofit back into the flood
     behavior a town should expect next time") to match the merit paragraph's physical advance.
  - **Deliberately left alone:** §2.2 "no duration-dependent fragility surface exists" — a gap
    section is *supposed* to say no X exists; the RO1 fix dropped a first-ness **headline in the
    merit argument**, it did not make gap statements illegitimate. §2.3 close already matches
    advance three. Transferability scope ("any gauged municipality" here vs "any domain with…" in
    the merit ¶) — different altitudes, not a contradiction.

## Done in previous session (2026-07-28)
- **Consistency sweep — carried the closed-loop framing through the rest of the document** (commit `d6f0ff8`). Four real mismatches fixed, net-neutral on length:
  1. **Goals box close** — "The defining intellectual contribution is treating community knowledge as a structural calibration input…" flatly contradicted the merit paragraph. Now leads on the loop, with RQ2's measured-against-a-baseline claim kept but subordinate. Also killed "epistemic contribution."
  2. **Status quo close** (~line 159) — terminated the section on community knowledge alone; now points forward to closing the loop.
  3. **RO1 Key outcomes** (~line 858) — dropped the "the first for the older load-bearing masonry" novelty claim, since the merit paragraph deliberately abandoned that first-ness claim. Fragility surfaces stay as an RO1 outcome, just not a headline novelty.
  4. **Transformative impact** (~line 1123) — "free geospatial data" → "free national data" so the generalizable-principle sentence matches its twin in the merit paragraph.
  - Also de-layered the stale "PHYSICAL layer" wording in the Subtask 1.4 `%` comment. **Trim TODO left in place** — that's a cut task.
  - **Checked and deliberately left alone:** goals box body 68–78 (already IS the loop story; masonry there is a subordinate appositive that makes the physical coupling concrete); RO1 box and Subtask 1.4 prose (the deemphasis was about the *merit argument*, not erasing masonry from the research plan — Subtask 1.4 genuinely is fragility surface development); "these four domains" at line 208 (a count of the four literature subsections, unrelated to the old four-*layer* taxonomy — verified).
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
