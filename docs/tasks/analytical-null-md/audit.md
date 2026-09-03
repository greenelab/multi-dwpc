# Analytical null integration — drift audit

**Task** `analytical-null-md` · **Step** audit (plan.md Task 5)
**Pass** first pass — fresh reader (the auditor wrote none of this work)
**HEAD audited** `662838bbd58595beddf5764a448d2819b26ad00d`
("Correct verification prose: warm-up asymmetry and negative-control fallback")
**Base for the reverse pass** `106cb7236d9472c1db617a1735f31590d69389b3`
(confirmed: `git merge-base HEAD upstream/main` == 106cb72)
**As of** 2026-09-01.

Method: `docs/tasks/analytical-null-md/design.md` was read as 42 numbered,
checkable claims (C1–C42) and each was checked against the landed tree with
named evidence. The branch diff (`git diff 106cb72..HEAD`, 29 files) was then
walked file by file and each change tied back to a design clause or a
`decisions.md` entry. `plan.md` and `verification.md` were read as secondary
claim sources; departures found in them are carried in the findings list.

Verdicts: **aligned** (claim and tree agree) · **recorded** (they differ and a
dated `decisions.md` entry records the difference) · **drift** (they differ
silently). One claim (C34, `summary.ipynb`) is marked **deferred**: the design's
own step order places audit *before* summarise, so there is nothing yet to check
and it is not counted as drift.

**Counts** — 42 claims audited: **27 aligned · 2 recorded · 12 drift · 1
deferred**. 18 findings (F1–F18), all documentary or trivial in severity; none
requires a code change to make a correct document true.

---

## Forward pass — design claims against the tree

| # | Claim (short) | Evidence | Verdict | Note |
|---|---|---|---|---|
| C1 | Branch `analytical-null-md` cut from `upstream/main`, 106cb72 | `git merge-base HEAD upstream/main` → 106cb72 | aligned | |
| C2 | Steps: promote, adapt, wire, verify, audit, summarise | commits 99985b2 / e53018e / 9ae627c / c47c7ec; plan.md Tasks 1–6 | aligned | audit is this document; summarise pending |
| C3 | Validated by `docs/tasks/capacity-hurdle-adaptive-null/` on `fix/random-null-stratified-srswor`, evidence `output/tier0_capacity_hurdle/`, commit 3b0bdee | `git ls-tree 3b0bdee docs/tasks/`, `… output/` — both paths present | aligned | |
| C4 | Today's `query_metapath_z` scores against `b` random same-size subsets of the whole gene universe | `git show 106cb72:src/multi_dwpc_query.py` (`rng.choice(gene_pool, …)` × `b`) | aligned | |
| C5 | MC null replaced by exact resampling moments stratified by leave-target-out capacity; same z every time; no `b` to tune | `src/analytical_null.py:88-137`; `src/multi_dwpc_query.py:145`; verification.md "Determinism … True" | aligned | |
| C6 | Latency, measured at verify, is disk-dominated and roughly unchanged | `tables/timing.csv` medians 103.217 s / 85.302 s (recomputed); verification.md §Timing | recorded | decisions.md 2026-09-01 "task 4, measured" replaces the original order-of-magnitude expectation |
| C7 | Reads `data/dwpc_cache/` as `dwpc_<metapath>_d0.50.npz`, raw scale; `scripts/prewarm_dwpc_cache.py` is the warm-up path; a missing matrix is computed on demand, read-only | `src/dwpc_direct.py:499-524`; `scripts/prewarm_dwpc_cache.py` exists; `data/dwpc_cache/` holds 52 `_d0.5.npz` + 406 legacy `_d0.50.npz` | **drift** | F1 (filename), F2 (read-only) |
| C8 | Per-metapath `raw_mean` from `data/metapath-dwpc-stats.tsv` | `src/dwpc_direct.py:89` (`cache_path = data_dir / "metapath-dwpc-stats.tsv"`); `src/analytical_null.py:114` | aligned | |
| C9 | Node tables map Entrez IDs and the target identifier to matrix positions | `src/multi_dwpc_query.py:_gene_index_map` / `_target_position`; `hetmat.get_nodes` | aligned | |
| C10 | `T = mean(arcsinh(DWPC/raw_mean))` — the transformed scale `get_dwpc_for_pairs(..., transform=True)` reports | `src/dwpc_direct.py:180-198` (`np.arcsinh(raw/mean)`), `:684` (`transform: bool = True`), `:726-728`; `src/analytical_null.py:113-115` | aligned | old call omitted `transform=`; the default supplied it, so the scale is unchanged |
| C11 | Why the null is stratified (rationale paragraph) | design.md:37-39 | **drift** | F3 — the paragraph does not parse |
| C12 | `c_g = sum_{t'≠t} DWPC_m(g,t')` — raw scale, same damping, excludes the tested column | `src/capacity.py:20-24`; `src/analytical_null.py:111,117` (one `compute_dwpc_matrix_csc` at `DEFAULT_DAMPING`); tests/test_capacity.py::test_capacity_is_rowsum_minus_target_column | aligned | |
| C13 | Zero-capacity genes form one stratum | `src/hurdle_adaptive_bins.py:37-39`; tests/test_hurdle_adaptive_bins.py::test_zero_keys_form_exclusive_hurdle_stratum | **drift** | F4 — the justifying sentence does not parse |
| C14 | Positive-capacity genes grouped into strata of at least 50, never splitting equal capacity | `src/hurdle_adaptive_bins.py:45-68`; tests …::test_equal_keys_never_split, ::test_min_stratum_size_met_except_none | **drift** | F6 — `:57` allows one under-sized stratum |
| C15 | Null set = same-stratum draw, w/o replacement, user's genes excluded; mean/sd exact via `exact_resampling_moments`; z and p derived locally | `src/pool_assembly.py:41-49`; `src/analytical_null.py:56,68-70`; `src/hetnetex_md_import.py:22,24`; tests/test_hetnetex_md_import.py::test_srswor_moments_match_exhaustive_enumeration | aligned | |
| C16 | Calibrated: upper tail exactly 0.0500, CI [0.0449, 0.0554] | `3b0bdee:docs/tasks/capacity-hurdle-adaptive-null/verification.md:148,152-153` (0.04494026, 0.05544988, 0.050000000000000044) | aligned | inherited claim; the cited row checks out to 4 dp |
| C17 | Deficient stratum merges into its lower-capacity neighbour; every merge reported on the result | `src/hurdle_adaptive_bins.py:95-108`; `src/analytical_null.py:122,136`; tests …::test_lowest_stratum_merges_upward, ::test_cascading_merge | **drift** | F7 — the lowest stratum merges *upward* |
| C18 | Zero-variance null → NaN z and p on a kept row, matching `std == 0` today | `src/analytical_null.py:58-66`; tests/test_query_metapath_z.py::test_zero_variance_metapath_row_is_nan_not_dropped; verification.md "6 of 52 … NaN" (recomputed: 6) | aligned | |
| C19 | An unresolvable metapath is skipped | `src/multi_dwpc_query.py:145-150`; tests/test_query_metapath_z.py::test_unresolvable_metapath_is_skipped_not_fatal | **drift** | F5 — clause ends "skipped, as 2026-09-01" |
| C20 | `analytical_gene_set_z(hetmat, metapath, source_idx, target_pos, *, min_stratum_size=50)` → `real_mean, null_mean, null_std, z, p_value, n_active_strata, merges` | `src/analytical_null.py:77-95` (`GeneSetZResult` fields match one-for-one) | aligned | `n_active_strata` defined in decisions.md 2026-09-01; implemented at `:127` as `count > 0` |
| C21 | `query_metapath_z` keeps its signature; `b`/`seed` accepted, ignored, warn; output gains `p_value`; scores stay transformed; ranked by `effect_size_z` | `src/multi_dwpc_query.py:94-95,118-125,152-160,164`; tests/test_query_metapath_z.py::test_b_and_seed_deprecated_and_inert, ::test_frame_columns_and_ranking | **drift** | F8 — parameter defaults/types changed; the rest is exact |
| C22 | `app.py` stops plumbing `b`/`seed`, names the analytical null in progress text, no other UI change | `git diff 106cb72..HEAD -- app.py` — 2 lines: dropped `DEFAULT_B` import, `app.py:463` spinner text | aligned | |
| C23 | `pyproject.toml` pins `hetnetex-md @ …tghosh30/HetNetEx-MD.git@26f5ba85…`; no submodule | `pyproject.toml:25`; `direct_url.json` in the env reports commit_id 26f5ba85…; `.gitmodules` unchanged from base and lists only `connectivity-search-backend` | aligned | |
| C24 | Every figure drawn from a table committed in this folder; controls beside real data on shared axes; figures land in `figures/`, drawn by the verification run | `verify_query.py:213,251,336,357` (tables) and `:395,412,441,455` (figures) — one script, one pass | aligned | figures 1–2 use the in-memory frame that is written to `per_metapath_comparison.csv` in the same run |
| C25 | adapt *positive*: planted enrichment → z > 1.65; a pool planted with known mu/sigma recovered within float tolerance | tests/test_analytical_gene_set_z.py::test_planted_enrichment_yields_high_z (aligned); ::test_null_moments_match_manual_pipeline | **drift** | F10 — the moments test re-derives with the same kernel |
| C26 | adapt *negative*: stratum-matched draw ≈ 0 and standard-normal shaped; zero-variance column → NaN row | tests/test_analytical_gene_set_z.py::test_stratum_matched_null_draws_are_standard_normal_shaped (200 draws, \|mean\|<0.15, \|sd−1\|<0.15), ::test_zero_variance_column_yields_nan_without_raising | aligned | |
| C27 | adapt *figures*: old-vs-new z scatter (MC at 3 seeds vs the single analytical z) and the seed-to-seed MC spread as bars | `figures/old_vs_new_z_scatter.png`, `figures/mc_seed_spread.png` (both inspected); `tables/per_metapath_comparison.csv` | aligned | spread figure shows bar + dot on shared axes, as declared |
| C28 | wire *positive*: the query call returns the ranked frame with `p_value` | tests/test_query_metapath_z.py::test_frame_columns_and_ranking (asserts the 7-column set incl. `p_value`) | aligned | |
| C29 | wire *negative*: passing `b`/`seed` warns and changes nothing (two calls, identical frames) | tests/test_query_metapath_z.py::test_b_and_seed_deprecated_and_inert (`pdt.assert_frame_equal` ×3) | aligned | |
| C30 | verify *positive*: a real query **in the running streamlit app** returns a ranked table; run twice, byte-identical | verification.md §"Positive control" — two direct `query_metapath_z` calls; §Streamlit records AppTest load + headless-serve 200, run-query flow not driven | **drift** | F11 — real departure, described in verification.md but not in decisions.md |
| C31 | verify *negative*: unknown gene IDs fail with the existing clear error, not a crash | verification.md §"Negative control" — `ValueError: None of the provided gene IDs were found in Hetionet`; `src/multi_dwpc_query.py:139-140` | **drift** | F11 — behaviour correct; ran via direct call, not the app |
| C32 | verify *figures*: wall-clock old-vs-new from a committed timing table; rank-agreement scatter with Spearman rho printed | `figures/timing_comparison.png` (inspected: 2 panels, 0.83x / 0.1x); `figures/rank_agreement.png`; `tables/timing.csv`, `tables/warm_metapath_timing.csv`, `tables/rank_agreement.csv` | aligned | rho recomputed from the committed table: 0.2605, p = 0.08771, n = 44 — matches verification.md |
| C33 | audit: design read as numbered claims, forward and reverse, verdicts aligned/recorded/drift, artifacts checksummed | this document | aligned | |
| C34 | summarise: `summary.ipynb`, committed without outputs — gate 2 | absent from the tree | deferred | Task 6; the design places audit before summarise, so not drift |
| C35 | Deterministic: identical queries return identical z; `b`/`seed` inert | verification.md `new_df_1.equals(new_df_2) → True`; tests/test_query_metapath_z.py::test_b_and_seed_deprecated_and_inert | aligned | |
| C36 | Cost no longer buys noise: new 103 s vs old 85 s over 52 metapaths; warm-matrix 70 ms vs 8.6 ms; ~16 % relative sd error at b = 20; validated 213x / 2,160x vs B = 1,000 / 10,000 | `tables/timing.csv` medians 103.217/85.302 ✓; `tables/warm_metapath_timing.csv` GiGiGpBP medians 69.98/8.59 ms ✓; 1/sqrt(2·19) = 0.162 ✓; `3b0bdee:…/verification.md:307-308` = 214.03x / 2144.82x ✗ | **drift** | F12 — three of four sub-claims recompute exactly; the speedup pair is misquoted |
| C37 | Five modules carried file-for-file from `src/tier0/`, tests included; hardening = merge-path + capacity-key tests, `min_stratum_size <= 0` guard, damping defined once, no dead parameters | `diff` of each module against `3b0bdee:src/tier0/…` shows only import-path edits + the three itemised deviations (commit 99985b2 lists them); `src/dwpc_direct.py:431`; `src/hurdle_adaptive_bins.py:27-28`; tests …::test_forced_deficient_stratum_reports_merges, ::test_one_matrix_load_per_evaluation; no unused parameter found in the five modules | aligned | `src/pool_assembly.py` is byte-identical to `3b0bdee:src/tier0/_pool_assembly.py` |
| C38 | The 5-row promotion mapping table | verified row by row via `git show 3b0bdee:<path>` + `diff` | aligned | |
| C39 | The example query's analytical z correlates with the Monte-Carlo z | recomputed from `tables/rank_agreement.csv`: Spearman rho 0.2605, p = 0.0877 (n = 44); Pearson on z 0.881, driven by 3 high-z rows | **drift** | F13 — bullet not revisited when Expected result was rewritten |
| C40 | Seed-to-seed MC spread visible at default `b`; the analytical value has none | `figures/mc_seed_spread.png`; `tables/per_metapath_comparison.csv` `z_mc_std` (e.g. GpPWpGpBP 9.72, GpBPpGpBP 8.29) | aligned | |
| C41 | Latency dominated by matrix loads and roughly unchanged; the null share no longer scales with any `b` | verification.md §Timing; `src/multi_dwpc_query.py` has no `b` in the compute path | recorded | rewritten per decisions.md 2026-09-01 "task 4, measured" |
| C42 | Out of scope: `scripts/permutation_null_datasets.py`, `scripts/random_null_datasets.py` | absent from `git diff 106cb72..HEAD --stat` | aligned | |

---

## Reverse pass — landed changes against a motivating clause

Every one of the 29 changed files traces to a clause or a decisions entry. One
*behaviour* inside a motivated file does not.

| Landed file / group | Motivating clause | Verdict |
|---|---|---|
| `.gitignore` (+6) | decisions.md 2026-09-01 — "task's docs land … via a scoped negation of the repository's blanket `docs/` gitignore" | recorded |
| `app.py` (−1 import, spinner text) | design Interface bullet 3 (C22) | aligned |
| `pyproject.toml` (+1) | design Interface bullet 4 (C23); decisions.md 2026-09-01 (pinned direct reference, not a submodule) | aligned |
| `src/dwpc_direct.py` (+8: `DEFAULT_DAMPING`, default swap) | plan Task 1 hardening item 2 | aligned (two stray blank lines at `:428-430` are cosmetic) |
| `src/hurdle_adaptive_bins.py` (new) | design Promoted apparatus row 1; plan Task 1 hardening item 1 | aligned |
| `src/pool_assembly.py` (new) | design Promoted apparatus row 2 | aligned |
| `src/capacity.py` (new) | design Promoted apparatus row 3; hardening item 2 | aligned |
| `src/hetnetex_md_import.py` (new) | design Promoted apparatus row 4; docstring rewritten submodule → pyproject pin per decisions.md 2026-09-01 | recorded |
| `src/analytical_null.py` (new) | design Promoted apparatus row 5 + Interface bullet 1; decisions.md 2026-09-01 (task 2, capacity path) | aligned |
| `src/multi_dwpc_query.py` (modified) | design Interface bullet 2; plan Task 3 | **drift** — F9, the new `if not rows: raise ValueError` at `:162-163` has no clause behind it |
| `tests/test_{hurdle_adaptive_bins,capacity,hetnetex_md_import,analytical_null}.py` | plan Task 1 "and their test files under `tests/`"; commit 99985b2 itemises the two beyond the explicit list | recorded (F15 notes the missing `test_pool_assembly.py`) |
| `tests/test_analytical_gene_set_z.py` (new, 7 tests) | plan Task 2 "Expected test outcomes" — all seven land | aligned (see F10) |
| `tests/test_query_metapath_z.py` (new, 6 tests) | plan Task 3 "Expected test outcomes" | aligned except the one test backing F9 |
| `docs/tasks/analytical-null-md/{design,plan,decisions,verification}.md` | design §"Figures, controls and the reviewer's path"; plan Tasks 4–5 | aligned |
| `docs/tasks/analytical-null-md/verify_query.py` (483 lines) | plan Task 4 — "the comparison script"; old implementation from `git show upstream/main:…` (`verify_query.py:106-121`) | aligned |
| `docs/tasks/analytical-null-md/tables/*.csv` (4) | design adapt + verify figure clauses; plan Task 4 | aligned |
| `docs/tasks/analytical-null-md/figures/*.png` (4) | design's four declared figures; decisions.md 2026-09-01 (task 4) corrects plan.md's "three" | recorded |

Untracked-but-present: the `data` symlink (`git status` → `?? data`). Not
committed, correctly — but see F18: it is **not** in fact ignored.

---

## Artifact checksums

Every artifact read during this audit, `md5 -r`:

```
16a20da20f1270074f62ba037528ba25  docs/tasks/analytical-null-md/tables/per_metapath_comparison.csv
cb1a66f7d33ff78f1ed95148b383bcd4  docs/tasks/analytical-null-md/tables/rank_agreement.csv
299f73d8e59a28666484c189be0fc4f1  docs/tasks/analytical-null-md/tables/timing.csv
1d3a686a5c1dec8c61f88d0aa48412c5  docs/tasks/analytical-null-md/tables/warm_metapath_timing.csv
69ca30e965376cd0a34ac9e252033c77  docs/tasks/analytical-null-md/figures/mc_seed_spread.png
f7eddbff3946dbc16e7161bdea96d235  docs/tasks/analytical-null-md/figures/old_vs_new_z_scatter.png
76429c4f2476241eb7251ebe776f3a8c  docs/tasks/analytical-null-md/figures/rank_agreement.png
118bef83b9c0a4258293515b8541e0bd  docs/tasks/analytical-null-md/figures/timing_comparison.png
```

Recomputations run against these exact files (conda env `multi_dwpc`), all
matching the committed prose: timing medians 103.217 s / 85.302 s; warm medians
GpBP 12.600 / 7.458 ms and GiGiGpBP 69.980 / 8.590 ms; 6 of 52 analytical NaN
and 10 of 52 any-MC-seed NaN, the former a strict subset of the latter; Spearman
rho 0.2605 (p = 0.08771, n = 44); 52 − 8 = 44 rank rows. `pytest -q
--collect-only` reports 65 tests across 8 files, matching verification.md's
"65 passed".

---

## Findings

Severity key: **documentary** = the document is wrong about a tree that is
right; **trivial** = imprecision a reader would not be misled by but the fix
wave should still take.

**F1 — `dwpc_<metapath>_d0.50.npz` names files the query path does not read.**
`src/dwpc_direct.py:506` builds the cache name with `repr(float(damping))`, i.e.
`dwpc_<metapath>_d0.5.npz`. `data/dwpc_cache/` holds 52 `_d0.5.npz` files (the
`G -> BP` set this query reads, mtime 2026-05-29) beside 406 legacy `_d0.50.npz`
files (mtime 2026-04-23, written before the naming change in b054f76). *Fix
(document):* design.md:22 — `d0.50` → `d0.5`. Severity documentary.

**F2 — "read-only" / "No file under `data/` is written" is asserted, not
enforced.** design.md:23 says a missing matrix is "computed on demand,
read-only" and plan.md:26-27 makes it a global invariant. `HetMat.__init__`
defaults `write_disk_cache=True` and `_save_to_disk` (`src/dwpc_direct.py:518-524`)
writes the `.npz` on a miss; `load_metapath_stats` (`:89-93`) downloads the stats
TSV into `data/` when absent. Neither `app.py:55` nor `verify_query.py:157`
passes `write_disk_cache=False`. No write actually occurred in the verify run —
all 52 caches predate it — so the tree behaved as claimed by luck of a warm
cache, not by construction. *Fix (document):* reword both to "reads the existing
cache; a miss falls through to HetMat's existing compute-and-write behaviour,
unchanged by this task", since changing `HetMat`'s default is out of this task's
scope. Severity documentary (latent).

**F3 — the stratification rationale does not parse.** design.md:37-39:
"Therefore, what is "a random gene set".  Stratifying adjust for
degree-confounding in the random set by …". Sentence fragment plus
subject-verb disagreement. *Fix (document):* "Therefore the question is what
counts as 'a random gene set'. Stratifying adjusts for degree confounding by
comparing each of the user's genes only against genes of similar connectivity
…". Severity documentary.

**F4 — the zero-capacity justification does not parse.** design.md:48-49: "Zero
capacity means no these genes score zero for every target." The tree's own
statement of the same fact is correct: `src/capacity.py:9-10`, "a zero row means
zero DWPC to every target of the metapath". *Fix (document):* "Zero capacity
means these genes score zero for every other target of the metapath, so they
form an exact exchangeability class." Severity documentary.

**F5 — dangling date in the degenerate-cases clause.** design.md:66: "A metapath
whose matrix cannot be resolved is skipped, as 2026-09-01." *Fix (document):*
drop ", as 2026-09-01" (the section already carries the doc-level **As of**
date at design.md:6). Severity documentary.

**F6 — "strata of at least 50 genes each" has an implemented exception.**
`src/hurdle_adaptive_bins.py:57` merges a trailing under-filled stratum into its
predecessor *only* when `len(boundaries) > 1`; when the whole positive-capacity
set is smaller than `min_stratum_size` it stands alone below the minimum. The
module docstring (`:24-25`) and tests/test_hurdle_adaptive_bins.py::test_min_stratum_size_met_except_none
both record this; the design does not. *Fix (document):* design.md:50-52 — append
"(a single positive stratum smaller than the minimum is kept rather than
dropped)". Severity documentary.

**F7 — "merges into its lower-capacity neighbour" is wrong for the lowest
stratum.** `src/hurdle_adaptive_bins.py:102`: `into = deficient - 1 if deficient
> 0 else deficient + 1`; tests/test_hurdle_adaptive_bins.py::test_lowest_stratum_merges_upward
pins the upward case. The same incomplete wording was copied into the production
docstring at `src/multi_dwpc_query.py:105`. *Fix (document):* design.md:60-62
and the `query_metapath_z` docstring — "merges into its lower-capacity
neighbour (the higher-capacity one when it is already the lowest stratum)".
Severity documentary.

**F8 — "`query_metapath_z` keeps its signature" overstates.** Parameter names and
order are unchanged and every existing call still type-checks, but the
annotations and defaults did change: `b: int = DEFAULT_B` → `b: int | None =
None` and `seed: int = 42` → `seed: int | None = None`
(`src/multi_dwpc_query.py:94-95`). The change is load-bearing — `None` is what
distinguishes "not passed" from "passed", and so gates the DeprecationWarning at
`:118`. *Fix (document):* design.md:75 — "keeps its call signature; `b`/`seed`
now default to `None` so that passing them can be detected and warned on".
Severity trivial.

**F9 — landed behaviour with no motivating clause: the all-skipped
`ValueError`.** `src/multi_dwpc_query.py:162-163` raises
`ValueError("No metapaths could be scored with the analytical null")` when every
metapath was skipped, pinned by
tests/test_query_metapath_z.py::test_all_metapaths_unresolvable_raises_value_error.
No design clause, plan interface line, or decisions entry describes it; the
design's degenerate-cases paragraph stops at "is skipped". The behaviour is
sensible (an empty frame would break the ranked-table contract downstream), so
the tree is right and the documents are short. *Fix (document):* one clause in
design.md's "Degenerate cases" — "if no metapath resolves, the query raises
rather than returning an empty frame" — or a dated decisions.md entry.
Severity documentary.

**F10 — the declared "planted mu and sigma" control is a same-kernel
re-derivation.** design.md:100-101 declares "a pool planted with a known mu and
sigma is recovered by the moments within float tolerance". The landed adapt-step
test, tests/test_analytical_gene_set_z.py::test_null_moments_match_manual_pipeline,
rebuilds the partition by hand and then calls the *same*
`exact_resampling_moments` (`:93`), so it verifies the adapter's wiring, not the
moments themselves; a wrong kernel would pass it. The independent ground truth
does exist, one step earlier:
tests/test_hetnetex_md_import.py::test_srswor_moments_match_exhaustive_enumeration
checks the moments against exhaustive enumeration. *Fix (document):*
design.md:100-101 — split the control in two, naming the enumeration test as the
moments check (promote step) and the manual-pipeline test as the adapter-wiring
check (adapt step). Severity documentary.

**F11 — the verify step's controls were not run through the running app, and the
departure is not in the decisions ledger.** design.md:118-121 declares
*positive*: "a real query in the running streamlit app returns a ranked table;
run twice, byte-identical" and *negative*: an unknown-ID gene set "fails with the
existing clear error"; plan.md:109-110 repeats "app launched in the env, example
query executed". What landed: an `AppTest` headless render, a
`streamlit run --server.headless` HTTP-200 smoke, and two direct
`query_metapath_z` calls for determinism and for the unknown-ID error
(verification.md §Streamlit, §Positive control, §Negative control). The reason —
`app.py:55-58` preloads all 52 `G -> BP` matrices (~30 GB decompressed) on the
first "Run query" click, which this machine cannot absorb — is documented at
length in verification.md and is a good reason; verification.md even calls it
"per the task's documented fallback". But `decisions.md` is this task's
departure ledger and carries no entry for it, so by the audit's own rule this
is drift rather than a recorded departure. The evidence itself is sound:
`query_metapath_z` is the identical code path the app calls at `app.py:464`.
*Fix (document):* add a dated `decisions.md` entry under a "verification.md"
heading recording the memory-bounded fallback and what it substitutes for, and
add the same one-line pointer to plan.md Task 4. Severity documentary.

**F12 — "213x / 2,160x" misquotes the validation's measured speedups.** The
source of truth,
`3b0bdee:docs/tasks/capacity-hurdle-adaptive-null/verification.md:307-308`,
reports `214.031893` at B = 1,000 and `2144.821029` at B = 10,000, and its prose
at `:315-316` rounds these to "214x" and "2,145x". Neither "213" nor "2,160"
appears anywhere in that branch's design.md or verification.md. The wrong pair
appears three times on this branch: design.md:152, decisions.md:76,
verification.md:318. *Fix (document):* replace all three with "214x /
2,145x". Severity documentary.

**F13 — "The example query's analytical z correlates with the Monte-Carlo z" is
not what was measured.** Recomputed from `tables/rank_agreement.csv`: Spearman
rho = 0.2605, p = 0.0877 over n = 44 — weak, and not significant at alpha = 0.05.
(Pearson on the raw z is 0.881, but that is carried almost entirely by the three
high-z metapaths, which is exactly what `figures/old_vs_new_z_scatter.png`
shows.) verification.md §"Rank agreement" states this honestly and frames the
reordering as the desired behaviour; the design's Expected result bullet was
simply not revisited when decisions.md 2026-09-01 rewrote the latency bullet
beside it. *Fix (document):* design.md:175 — "The example query's analytical z
agrees with the Monte-Carlo z at the top of the ranking and reorders the middle
substantially (measured Spearman rho 0.26, verify step)". Severity documentary.

**F14 — plan.md Task 4's expected outcome still carries the superseded latency
expectation.** plan.md:113-114: "latency drops by an order of magnitude or more;
analytical z generally at or below MC z". The first is the exact expectation
decisions.md 2026-09-01 records as wrong and corrects in the design; the second
holds for only 28 of 44 comparable rows (recomputed). The same decisions entry
did correct plan.md's "three declared figures", so the plan was in scope for
correction and this line was missed. A pre-registered expectation is worth
preserving rather than silently rewriting. *Fix (document):* leave the text and
append "(superseded — see decisions.md, 2026-09-01, task 4, measured)".
Severity documentary.

**F15 — plan.md Task 1 promises a test file `src/pool_assembly.py` never had.**
plan.md:32-35 says the five modules arrive "and their test files under
`tests/`". No `tests/test_pool_assembly.py` landed, and none exists at 3b0bdee
either (the validation branch covered `_pool_assembly` through
`tests/tier0/test_pool_construction.py`, which was not promoted). Its coverage
here is indirect, through tests/test_analytical_gene_set_z.py, which drives
`pools_from_bins` on every call. *Fix (document):* plan.md:32-35 — "and the test
files that exist for them (`pool_assembly` is covered indirectly via the adapter
tests)". Severity trivial.

**F16 — "10-30+ z units" overstates one of the three named metapaths.**
verification.md:226-227 names `GpPWpGpBP`, `GpBPpGpBP`, `GpMFpGpBP` as showing
"MC seed ranges spanning 10-30+ z units". Recomputed from
`tables/per_metapath_comparison.csv`: 19.34, 15.88 and **4.75** respectively.
*Fix (document):* "spanning 5-20 z units", or drop `GpMFpGpBP` from the list.
Severity documentary.

**F17 — "`pip show hetnetex-md` reports the pinned commit's version" is not what
pip prints.** plan.md:55-56. `pip show hetnetex-md` reports `Version: 0.1.0`; the
commit lives in `…/hetnetex_md-0.1.0.dist-info/direct_url.json`, which was
checked and does carry `26f5ba85bd9b9886eafebcbc544879e531f5262e`. The pin is
correct; the stated check would not have proved it. *Fix (document):* plan.md —
name `direct_url.json` (or `pip list --format=freeze`) as the check. Severity
trivial.

**F18 — the `data` symlink is not actually ignored.** verification.md:29-32 says
the symlink is "intentionally **not** committed — `data/` is gitignored". It is
not: `.gitignore:50` is `/data/`, whose trailing slash matches a *directory*,
while `data` here is a symlink, which git treats as a file. `git check-ignore -v
data` exits 1 and `git status --short` reports `?? data`. Nothing is committed
today, but a broad `git add` would stage a machine-specific absolute-path
symlink. *Fix:* correct verification.md's sentence, and — the one place a
one-character tree change is warranted — consider `/data` (no trailing slash) in
`.gitignore`. Severity documentary.

---

## What the audit did not find

- No landed file lacks a motivating clause; the only unmotivated item is one
  behaviour inside a motivated file (F9).
- The five promoted modules are byte-faithful to `3b0bdee:src/tier0/` apart from
  import paths and the three hardening items commit 99985b2 itemises — checked
  by `diff` on all five, not by reading.
- Every number quoted from a committed table in verification.md and design.md
  recomputes exactly, with the single exception of the inherited speedup pair
  (F12) and one span description (F16).
- The import policy holds: `hetnetex_md` is imported in exactly one place
  (`src/hetnetex_md_import.py:22`), only `exact_resampling_moments`, with a
  negative test guarding the surface.
- No test predating this branch was edited: `tests/` held only
  `test_dwpc_direct.py` and `test_dwpc_validation.py` at 106cb72, and neither
  appears in the branch diff — so plan.md Task 3's "any other existing-test edit
  is a finding" is satisfied vacuously and correctly.

---

## Re-audit — fresh reader, after the fix wave (044132e)

**Pass** re-audit, fresh reader (wrote none of the original work, the first
audit, or the fix wave). **Verified against** `git show 044132e` (6 files:
`.gitignore`, `docs/tasks/analytical-null-md/{decisions,design,plan,verification}.md`,
`src/multi_dwpc_query.py`) and the fix map at
`.superpowers/sdd/plan/task-5-fixwave-report.md`. Each finding below was
checked by reading the actual current document text at the cited location,
not the diff alone.

| Finding | Verdict | Evidence |
|---|---|---|
| F1 | fixed | design.md:22-25 now reads "`dwpc_<metapath>_d0.5.npz`, raw scale — the name the read path itself builds via `repr(float(damping))`", and separately notes the 406 legacy `_d0.50.npz` files "are not the read path's names." |
| F2 | fixed | design.md:25-27: "a matrix missing from the cache is read through HetMat's existing cache semantics, which fall through to compute-and-write on a miss, unchanged by this task." plan.md:26-28 matches: "write a cache entry under `data/` on a miss; this task's runs wrote nothing because the cache was warm throughout." |
| F3 | fixed | design.md:41-45 now parses: "Therefore the question is what counts as 'a random gene set'. Stratifying adjusts for degree confounding by comparing each of the user's genes only against genes of similar connectivity, …". |
| F4 | fixed | design.md:54-56: "Zero capacity means no paths into the metapath's target layer, so these genes score zero for every target." Parses; matches `src/capacity.py`'s own statement in substance. |
| F5 | fixed | design.md:75: "resolved is skipped, as before this change." No dangling date; grammatical. (Ruling differs from the audit's own suggested rewrite — drop the clause entirely — but resolves the defect the finding named.) |
| F6 | fixed | design.md:57-60: "grouped into strata of at least 50 genes each where the gene count allows; when the positive-capacity genes number fewer than 50 they form a single stratum, …" — states the exception the tree implements. |
| F7 | fixed | design.md:68-71: "a deficient stratum merges into its lower-capacity neighbour, and the lowest positive stratum merges upward." `src/multi_dwpc_query.py:104-105` docstring: "…lower-capacity neighbour, and the lowest stratum merged upward)". Diff of that file is one hunk inside the function's docstring only (5 lines changed, all prose); `python -m pytest tests -q` still reports 65 passed, 4 subtests passed — behavior unchanged. |
| F8 | fixed | design.md:86-91: "`query_metapath_z` keeps its parameters. … their defaults change to `None` so that passing them can be told apart from not passing them, which is what gates the warning; …". |
| F9 | fixed | design.md:75-77: "When every metapath is skipped, the existing `ValueError` ('No metapaths could be scored with the analytical null') contract is preserved." |
| F10 | fixed | design.md:110-116: "the adapter's moments match a direct kernel invocation on the same pools (a wiring check) — the independent enumeration check of the kernel itself lives in the promoted wrapper's tests (`tests/test_hetnetex_md_import.py`)." Splits the control as the finding asked. |
| F11 | recorded | decisions.md:97-110, new `## verification.md` heading, dated 2026-09-01: records both verify controls ran via direct `query_metapath_z` calls (not the live UI) because `app.py`'s `get_hetmat()` preload exceeds this machine's memory, and that verification.md documented the fallback but the ledger previously did not. The residual departure (controls not run through the live app) still stands — this is the recorded case, not a tree change. |
| F12 | fixed | All three sites corrected: design.md:166 "validated 214x / 2,145x against B = 1,000 / 10,000"; decisions.md:76 "The validated 214x / 2,145x speedup is real"; verification.md:322 "**214x / 2,145x per-row**". No remaining `213x`/`2,160x`/`2160x` in these three files. |
| F13 | fixed | design.md:190-191: "correlates weakly and positively with the Monte-Carlo z (Spearman rho 0.26, p = 0.088, n = 44)." decisions.md:87-95 adds the dated entry explaining the pre-run bullet expected stronger agreement than measured. |
| F14 | fixed | plan.md:118-120: "latency is disk-dominated and roughly unchanged" replaces "latency drops by an order of magnitude or more" — aligned with the corrected design rather than annotated as superseded (a different remediation than the finding's suggested fix, but the document and tree now agree, which is what "fixed" requires). |
| F15 | fixed | plan.md:35-38: "`src/pool_assembly.py` has no dedicated test file; its behaviour is covered indirectly, through the adapter (`tests/test_analytical_gene_set_z.py`) and binning (`tests/test_hurdle_adaptive_bins.py`) tests". |
| F16 | fixed | verification.md:229-231: "spanning roughly 5-20 z units (min-max: 19.34, 15.88, 4.75 respectively)" — matches the recomputed per-metapath figures from `tables/per_metapath_comparison.csv`. |
| F17 | fixed | plan.md:58-60: "the pin is verified via the installed package's `direct_url.json` (`pip show hetnetex-md` reports only the package version, not the commit)". |
| F18 | fixed | `.gitignore:50-53` adds a bare `data` line beside the existing `/data/`; `git check-ignore -v data` now reports `.gitignore:53:data\tdata` (exit 0), and `git status --short` is clean (no `?? data`). verification.md:30-35 rewrites the false "gitignored" claim into an accurate account of the gap and its closure. |

**Test suite:** `conda activate multi_dwpc && python -m pytest tests -q` →
65 passed, 4 subtests passed — unchanged from the fix map's own report, and
consistent with the tree-touching fixes (F7, F18) being non-behavioral.
`git status --short` on the worktree is clean before this audit's own commit.

**Tally:** 17 fixed, 1 recorded (F11), 0 still drift.

**Closing verdict: PASSES.** Every finding F1-F18 is either fixed (document
and tree now agree) or recorded (a dated `decisions.md` entry covers the
residual difference); none remains as silent drift.

### Notes (out of scope for this re-audit)

- plan.md Task 2's "Expected test outcomes" (line ~80) still says "planted
  mu/sigma recovered within float tolerance" — the same imprecision F10
  flagged in design.md, but plan.md was not in F10's fix location and the fix
  map did not touch it. Not a new finding under this re-audit's charter.

## Audit delta — B-sweep addition (6e33629), fresh reader

Fresh reader with respect to this surface: wrote none of the B-sweep
evidence, the design doc, or any prior audit pass. Read against
`git diff ee71b84..6e33629` (7 files) and the committed tree at 6e33629.
Constraint checked first: `design.md` does not appear in the commit's file
list (`decisions.md`, two `figures/*.png`, two `tables/*.csv`,
`verification.md`, `verify_query.py` — seven files, no `design.md`), so the
declared-figure list at design.md "Figures, controls and the reviewer's
path" (lines 99-124) is untouched and does not name `b_sweep_agreement.png`
or `b_sweep_tradeoff.png` — the no-new-hypothesis constraint holds.

| Artifact | Verdict | Evidence |
|---|---|---|
| `tables/b_sweep_per_metapath.csv` | aligned | 624 data rows (52 metapaths x 4 `B` x 3 seeds), columns `metapath,B,seed,z_mc,z_analytical,z_oldnull_b10k`. `z_analytical` is a map from the already-committed `per_metapath_comparison.csv` (`finalize_b_sweep`'s `z_analytical_by_mp`), not a second computation — the NaN pattern confirms it: 6 unique metapaths carry NaN `z_analytical` in the B-sweep table, exactly the "6/52" NaN count already recorded for `per_metapath_comparison.csv` in this audit's own summary table. Spot-recomputed Spearman rho directly from the raw rows at `B=10000, seed=42` (dropping NaNs, n=46): `0.9995065764076316`, which equals `b_sweep_summary.csv`'s `rho_seed42` for that row (`0.999507`) to displayed precision — the summary derives from this table, not independently. |
| `tables/b_sweep_summary.csv` | aligned | 7 rows (4 `same_null_mc` + 1 `analytical` + 1 `old_null_b10k`, one header). All four claimed same-null rho_mean values reproduce: 0.975194/0.992154/0.997944/0.999548 round to the verification.md table's stated 0.9752/0.9922/0.9979/0.9995. Old-null b10k rho_seed42 = 0.282921, rounds to the stated 0.283/0.2829. Wall-time entries match verbatim: analytical total_wall_time_s = 0.0063935830112313 (6.4ms stated), same-null B=10000 = 173.080257 (173.08s stated), old_null_b10k = 222.279652 (222.28s stated). |
| `figures/b_sweep_agreement.png` | aligned | Named against `b_sweep_summary.csv` by content, not just filename prefix: `finalize_b_sweep()` builds this figure directly from the in-memory `summary` DataFrame that is written to `b_sweep_summary.csv` in the same function call (`mc_summary`/`old_row` slices of `summary`), so table and figure are drawn from one shared object, not regenerated independently. File present, 97,304 bytes, matches the diff's stated size. |
| `figures/b_sweep_tradeoff.png` | aligned | Same construction as above — built from the same `summary` DataFrame (`mc_summary`, `analytical_row`, `old_row`) in the same `finalize_b_sweep()` call that writes `b_sweep_summary.csv`. File present, 60,171 bytes, matches the diff's stated size. |
| verification.md, "Agreement and cost versus B" section | aligned | All eight `b-sweep-*` commands present verbatim, both in this section and in "Reproducing" below it. Statistical method as described: same-null MC rebuilds the S1 partition via the public modules (`hurdle_adaptive_bins`, `pools_from_bins`, `merge_deficient_strata`) in the identical sequence `analytical_gene_set_z` uses internally (`src/analytical_null.py:88-137`, confirmed by direct comparison of the two code bodies), draws SRSWOR (argpartition-of-random-keys, without replacement) from `pools`/`counts` via `numpy.random.default_rng(seed)`, and computes `z_mc = (observed - mu) / sigma` from the empirical mean/std of `B` set-means (`T = totals / K`); `z_analytical` is read back from the committed table, never recomputed by a second method (`finalize_b_sweep` maps it in from `per_metapath_comparison.csv`; `run_analytical_reference_timing` calls `analytical_null` only for isolated wall-clock timing, its return value is discarded, not written into the per-metapath rows). Both stated caveats present, quoted: (a) "`old_null_b10k`'s row uses a single seed (42), so no seed-to-seed spread is reported for it (`n/a` in the table above) ... this matches what the task asked for ('one pass') but means the old-null point is not directly seed-range-comparable to the same-null curve"; (b) "The same-null MC wall-times and the analytical reference wall-time are measured on the identical basis (cached S1 partitions, no disk I/O, one full query across 52 metapaths) and are directly comparable; the old-null-b10k wall-time additionally includes matrix disk I/O ... so part of its 222s is not 'null generation' in the same narrow sense — it is still the correct number to report". No smoothing found: the section states the residual gap plainly (same-null MC at B=10,000 does not reach rho=1 exactly) and calls the old-null b=10,000 point "the single worst point on the whole plot" rather than softening it. |
| decisions.md entry | recorded | Dated `2026-09-02` entry names all five new artifacts, states the addition was made "after the audit passed, at Lucas's direction," frames it as Gate-2 presentation evidence for something "the design already implies" rather than a new hypothesis, and states explicitly that it "is not added to `design.md`" and the design's declared-figure list "is not retroactively edited to include these two figures." This is the entry that makes the whole addition an intentional, recorded departure from the audited baseline rather than silent drift. |

**Checksums (md5):**

```
f12211e54276b92494828d9bdb62bb23  tables/b_sweep_per_metapath.csv
97a683e4f37c4e93539c496e32058d94  tables/b_sweep_summary.csv
c5d1871981817f67e10f60cf88adbb73  figures/b_sweep_agreement.png
4869d8fb4f0ff62bf1fd3c7df94ac200  figures/b_sweep_tradeoff.png
```

**Closing verdict: PASSES.** All six new-artifact rows are aligned except
the decisions.md entry itself, which is recorded (by construction — it is
the record) rather than aligned against some prior baseline; nothing in this
delta is silent drift. design.md is untouched, confirmed from the commit's
file list, so the no-new-hypothesis constraint holds.
