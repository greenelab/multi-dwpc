# Laptop null bundle: verification

**As of** 2026-09-13. Measured outcomes against design.md's declared controls and
expected result. Tables in `tables/`, figures in `figures/`. Job logs are not
committed; they stay on Alpine `/projects/$USER/repositories/multi-dwpc-laptop-null-bundle/hpc/logs/null_bundle/`.

## port (local unit tests)

179 tests pass (`pytest tests`).

- *positive*: the summary null matches `analytical_gene_set_z` on 25 random
  columns with ties, hurdle zeros, merges and k >= N within 1e-12 relative
  (`tests/test_summary_null.py`); a planted enrichment gives z > 1.65;
  `QueryDwpc` rows match the cached hetmatpy matrices for GpBP, GbCbGpBP,
  GcGpBP, G<rGcGpBP and GpBPpGpBP within 1e-12 (`tests/test_query_dwpc.py`).
- *negative*: stratum-matched null draws have |mean z| < 0.15 and
  |sd - 1| < 0.15; a zero-variance pool gives NaN while 1e-6 relative spread
  stays finite; a large constant pool gives NaN like the reference; a changed
  data file is refused (`tests/test_null_bundle.py`); a gene with no edges
  gets zero DWPC; the query path imports neither `hetnetex_md` nor
  `scipy.stats`.

## build (Alpine, jobs 32526771 / 32526772)

- 52 array tasks completed (longest 7 min, max RSS 17.3 GB of the 32 GB requested;
  `tables/alpine_jobs.csv`); finalize wrote 52 metapaths, 97,790,661 strata rows
  (`tables/bundle_manifest.json`), 3.07 GB (`tables/bundle_files.csv`).
- *positive*: every task checked that stratum `n_genes` sum to 20,945 for all
  11,381 targets; no task raised.
- *negative*: no task raised on negative capacity; 578 residues were zeroed,
  all in GpBPpGpBP (Alpine log `build_32526771_34.out`).

## validate (Alpine, run 2: jobs 32527987 / 32527988)

`tables/validate_verdict.json`: 147,680 rows over 52 metapaths, 496 targets,
random gene sets of 5-500, and an annotated set for the 360 targets with at
least 5 annotated genes.

| criterion | expected | measured |
|---|---|---|
| stratum mismatches | 0 | **0** |
| max z relative difference, finite rows | <= 1e-9 | **7.1e-10** (GeAuGpBP, random 500) |
| merges identical to the reference | (implied) | all rows |
| NaN agreement | agree; rows inside the zero-variance band tabled | **not met for 23 rows**, all inside the band (reference null std <= 1.8e-16), tabled |
| shuffled-strata control | disagrees | median abs(dz) **1.46** |

`validate_verdict.json`'s `pass` flag checks only the first two criteria.

All 23 NaN disagreements (`tables/validate_nan_disagreements.csv`) have a
reference null std <= 1.8e-16 and |z_ref| <= 1.5: true zero-variance nulls that
the reference turns into a finite z through rounding. The summary path returns
NaN for them, as design.md's zero-variance rule intends.

Run 1 (jobs 32527631 / 32527632; Alpine logs `validate_32527631_*.out`,
`validate_summary_32527632.out`) met the two numeric criteria but had 24
NaN disagreements, one in the other direction (GpBP, target 1669: reference NaN,
summary z = 1.3e10; Alpine log `validate_32527631_33.out`). The tolerance was corrected
(decisions.md) and run 2 has none in that direction.
Figures: `validate_z_scatter.png` (laptop path on the diagonal, control on
shared axes), `validate_z_rel_diff_hist.png`.

## verify (this laptop: 24 GB, Apple Silicon)

`scripts/verify_null_bundle_app.py`; `tables/verify_query_timing.csv`,
`tables/verify_app.csv`, `tables/verify_app_memory_trace.csv`.

- *positive*: the Streamlit app, driven by Playwright, showed the ranked
  metapath table for the worked example and, after a rerun, for 200 random
  genes (`figures/verify_app_worked_example.png`, `verify_app_random_200.png`;
  the script waits for the tab label and for the app to go idle). The same
  queries in-process score all 52 metapaths and return identical frames on
  repeat (`tables/verify_query_timing.csv`).
- *negative*: in the app, symbols absent from Hetionet showed the existing
  input error, "No valid gene symbols or Entrez IDs parsed from input"
  (`figures/verify_app_unknown_genes.png`); in-process, Entrez IDs absent from
  Hetionet raise the existing "None of the provided gene IDs" error.

## Expected result against measured

| design.md expectation | measured | |
|---|---|---|
| full 52-metapath query on a 24 GB laptop, peak RSS <= 1.5 GB | streamlit server peak **1.19 GB** (previous app's preload: 33.6 GB) | met |
| query <= 1 s for 18 genes | **1.29 s** first query, **0.32 s** repeat | first query not met |
| query <= 3 s for 200 genes | **1.67 s** first and repeat | met |
| z matches the matrix adapter within validate tolerances | see validate | met |

The 18-gene first query misses the 1 s expectation: the first call loads and
degree-weights the edge matrices for all metapaths (~1 s); later queries reuse
them. In the app the page completes in 5.3-5.6 s per query, which includes the
unchanged path enumeration and subgraph drawing.
