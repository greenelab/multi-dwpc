# Laptop null bundle: audit

**As of** 2026-09-13. design.md read as numbered claims against the landed tree
(branch `laptop-null-bundle`), forward (design -> tree) and reverse (tree ->
design). Verdicts: **aligned** (tree matches), **recorded** (tree departs and
decisions.md says why), **drift** (tree departs silently; fixed in the fix wave),
**outcome** (a measured expectation, reported in verification.md, not a tree claim).

## Forward: design claims

| # | claim (design.md) | evidence | verdict |
|---|---|---|---|
| 1 | The app no longer preloads the G->BP matrices | `app.py` `get_query_resources` loads `NullBundle` + `QueryDwpc`; no HetMat import | aligned |
| 2 | Bundle is Parquet + manifest, built by one Slurm array task per metapath | `hpc/submit_null_bundle.sh` (array 0-51 + finalize); `src/null_bundle.py` | aligned |
| 3 | `strata.parquet` columns: metapath, target_position, stratum, capacity_min, capacity_max, n_genes, score_mean, score_centered_ss | parquet schema of the built bundle | aligned |
| 4 | `strata.parquet` sorted by `(metapath, target_position)` | sorted by target, then metapath, then stratum | recorded |
| 5 | About 83M strata rows | 97,790,661 (`tables/bundle_manifest.json`) | recorded |
| 6 | `row_sums.parquet`: metapath, gene_position, row_sum, 1.1M rows | 1,089,140 rows | aligned |
| 7 | Manifest: schema version, damping 0.5, min_stratum_size 50, SHA-256 of matrices, node tables, edge files, stats, git commit | manifest keys; also records `metagraph.json`, raw means, zero tolerance | aligned |
| 8 | Strata built exactly as today (hurdle_adaptive_bins on leave-target-out capacity; arcsinh scores) | `build_strata`; `scripts/build_null_bundle.py` | aligned |
| 9 | `src/query_dwpc.py` restricts hetmatpy products to query rows; no_repeats, short_repeat, BABA | `QueryDwpc.rows`; `tests/test_query_dwpc.py` | aligned |
| 10 | `src/null_bundle.py` reads the target's 52 slices and row sums after checking hashes | `NullBundle.__init__`, `strata_for_target` | aligned |
| 11 | `src/summary_null.py`: capacity placement, pool removal, merge as `merge_deficient_strata`, exact SRSWOR moments, z and p | `summary_gene_set_z`; tests against the matrix adapter | aligned |
| 12 | `query_metapath_z` and `app.py` use these instead of HetMat | `query_metapath_z(..., bundle, query_dwpc)` | aligned |
| 13 | "Intermediates and paths ... are unchanged" | `query_intermediates_and_paths` behaviour unchanged, but its `hetmat` argument was removed | drift -> recorded |
| 14 | HetMat, the matrix-based adapter and the pipeline stay for Alpine | `src/dwpc_direct.py`, `src/analytical_null.py` unchanged in role; used by build/validate | aligned |
| 15 | Rule: |DWPC| <= 1e-15 is 0 in the bundle and on the laptop | `DWPC_ZERO_TOL`; build script; `QueryDwpc.target_values` | aligned |
| 16 | Rule: pool SS <= 64 eps (stratum SS + removed squares) is 0 | eps x (64 + n) x scale after validation run 1 | recorded |
| 17 | Rule: missing/mismatched bundle errors with regeneration instructions, no fallback | `BundleMismatchError`, `REGENERATE_HINT`; app shows the error and stops | aligned |
| 18 | port controls (positive and negative) | `tests/test_summary_null.py`, `test_query_dwpc.py`, `test_null_bundle.py`; 179 pass | aligned |
| 19 | build scripts and checks (n_genes sum 20,945; no negative capacity) | `scripts/build_null_bundle.py` raises on both; all tasks completed | aligned |
| 20 | validate: 500 targets x random sizes 5-500 and GO-annotated sets | 496 targets; one annotated set (up to 200 genes) for the 360 targets with >= 5 annotated genes | drift -> recorded |
| 21 | validate pass criteria, shuffled-strata negative control, two figures | `validate_verdict.json` (its `pass` flag omits NaN agreement, which verification.md reports as not met for 23 rows); control uses capacities of random other genes; figures present | recorded |
| 22 | verify: app with 52 metapaths via Playwright; worked example ranked table identical on rerun | app rerun identity not checked in the browser; checked in-process | drift -> recorded |
| 23 | verify negative control and figures from committed tables | `verify_app_unknown_genes.png`; `verify_memory_and_latency.png` from `tables/verify_*.csv` | aligned |
| 24 | Expected result | `verification.md` table; 18-gene first query 1.29 s vs <= 1 s | outcome |
| 25 | Out of scope (Postgres loader, non-BP, pipeline nulls) | none implemented | aligned |

## Reverse: tree changes not described in design.md

| # | change | verdict |
|---|---|---|
| R1 | Gene IDs treated as a set in `query_metapath_z` | recorded |
| R2 | `GeneSetZResult` moved to `src/summary_null.py` | recorded |
| R3 | `CAPACITY_REL_TOL` rule for query-gene capacity | recorded |
| R4 | `scripts/summarize_null_bundle_validation.py`, `hpc/submit_validate_null_bundle.sh`, `scripts/verify_null_bundle_app.py` | aligned (the design's validate/verify steps) |
| R5 | `.claude/hooks/alpine-guard.sh` relaxed (outside the repository) | recorded |
| R6 | Zero-variance tolerance adds a pool-mean rounding term, n x (64 eps x mean)^2 | recorded (re-audit) |
| R7 | `assign_strata` sends a capacity in a gap between two ranges to the nearer range | recorded (re-audit) |
| R8 | `query_metapath_z` rejects non-BP targets and takes its default metapaths from the bundle | recorded (re-audit) |
| R9 | Port positive control compares against `analytical_gene_set_z` (which calls `exact_resampling_moments`), not the kernel directly | recorded (re-audit) |
| R10 | README "Run the web tool prototype" section describes the bundle | aligned |

## Fix wave

Drifts 13, 20 and 22 recorded in decisions.md (2026-09-13 audit entries). No
code change was needed: each is a documentation gap, not a behaviour defect.

## Re-audit (fresh reader)

A fresh reader checked every row, every number in verification.md, and the tree
for silent drift; tests passed (179). Findings and dispositions:

- Rows 20, 21, 24: verdicts corrected above (annotated sets exist for 360 of
  496 targets; row 21 is a departure; row 24 is an outcome).
- verification.md: the 17.3 GB peak and 3.07 GB size had no committed source
  (now `tables/alpine_jobs.csv`, `tables/bundle_files.csv`); run 1 had no
  logs cited (Alpine log paths now named; logs are not committed); the app check waited for the tab label, not 52
  rows, and the browser negative control hits the input parser's error (wording
  corrected); NaN agreement is a design criterion that 23 rows do not meet
  (now stated as not met).
- Silent drift R6-R9 recorded in decisions.md.

All numbers the reader checked matched the committed tables and logs.
