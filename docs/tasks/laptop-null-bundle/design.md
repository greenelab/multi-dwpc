# Laptop null bundle: run the web query without the DWPC matrices

**Task** `laptop-null-bundle` · **Status** OPEN · **Branch** `laptop-null-bundle` (cut from `upstream/main` 69237b5; its first commit carries the `fix/app-null-csc-memory` changes)
**Steps** port, build, wire, validate, verify, audit, summarise.
**As of** 2026-09-13. Decision lineage in [`decisions.md`](decisions.md).

## What this task establishes

The Streamlit app preloads every G->BP DWPC matrix: 52 metapaths, 33.6 GB in
memory, more than a 24 GB laptop has. The query does not need the matrices. Per
metapath it needs the query genes' DWPC to the target and, for the null, a
summary of every stratum of the target's column. This task precomputes those
summaries on Alpine, computes the query genes' DWPC from the edge files at query
time, and returns the same z as the exact matrix-based adapter.

## How it works

**Alpine builds a bundle** (`null_bundle/`, Parquet + manifest) from the existing
matrices, one Slurm array task per metapath:

- `strata.parquet`, one row per stratum per (metapath, target), sorted by
  `(metapath, target_position)`: `metapath`, `target_position`, `stratum`,
  `capacity_min`, `capacity_max`, `n_genes`, `score_mean`, `score_centered_ss`.
  About 83M rows (median 141 strata per metapath-target, max 355).
- `row_sums.parquet`: `metapath`, `gene_position`, `row_sum` (1.1M rows).
- `manifest.json`: schema version, damping 0.5, `min_stratum_size` 50, SHA-256
  of the matrices, node tables, edge files and `metapath-dwpc-stats.tsv`, git
  commit.

Strata are built exactly as today (`hurdle_adaptive_bins` on leave-target-out
capacity, scores `arcsinh(raw / raw_mean)`). The same tables load into Postgres
for production.

**The laptop answers a query** without HetMat:

1. `src/query_dwpc.py` computes the query genes' DWPC to the target from the
   degree-weighted edge matrices (~52 MB), restricting hetmatpy's products to
   those rows (categories present: no_repeats, short_repeat, BABA).
2. `src/null_bundle.py` reads the 52 stratum slices for the target and the row
   sums, after checking the manifest hashes against local files.
3. `src/summary_null.py` places each query gene in its stratum by capacity
   (`row_sum - dwpc`), removes the query genes from the pools, merges deficient
   strata as `merge_deficient_strata` does, and computes the exact SRSWOR mean
   and variance from the summaries, then z and p.
4. `query_metapath_z` and `app.py` use these instead of HetMat. Intermediates and
   paths already use only the edge files and are unchanged.

HetMat, the matrix-based `analytical_gene_set_z`, and the pipeline stay for
Alpine; the matrix-based adapter is the validation reference.

## Rules

- DWPC values with absolute value at most 1e-15 are 0, in the bundle and on the
  laptop. These are cancellation residues (~1e-18, some negative) where the true
  value is 0; real values are ~1e-5.
- A pool's centered sum of squares at most 64 eps (stratum SS + removed squares)
  is 0: z and p are NaN on a kept row.
- A missing or mismatched bundle is an error with regeneration instructions; no
  fallback to matrices.

## Steps, controls, figures

- **port**: the three laptop modules, unit-tested.
  - *positive*: summary moments equal `exact_resampling_moments` on random pools,
    including merges and k >= N (<= 1e-12); planted enrichment gives z > 1.65;
    `query_dwpc` matches cached matrices for all three categories.
  - *negative*: stratum-matched null draws score standard-normal; a zero-variance
    pool gives NaN while a pool with 1e-6 relative spread stays finite; a
    hash-mismatched bundle is refused; a gene with no edges gets zeros.
- **build**: `scripts/build_null_bundle.py`, `hpc/submit_null_bundle.sh`,
  `scripts/finalize_null_bundle.py`.
  - *positive*: per (metapath, target), stratum `n_genes` sum to 20,945; row
    counts match the manifest.
  - *negative*: no negative capacity after the zero rule.
- **wire**: `query_metapath_z`, `query_intermediates_and_paths`, `app.py`.
- **validate** (Alpine): `scripts/validate_null_bundle.py`, 52 metapaths x 500
  targets (worked example plus a degree-stratified sample) x gene sets of size 5,
  18, 50, 200, 500 (random and GO-annotated), against the matrix-based adapter.
  - *pass*: 0 stratum mismatches; max z relative difference <= 1e-9 on finite
    rows; NaN agreement, with rows inside the zero-variance band tabled.
  - *negative control*: the same comparison with shuffled strata must disagree.
  - *figures*: summary z against reference z with the shuffled control on shared
    axes; histogram of log10 relative difference.
- **verify** (laptop): the Streamlit app with all 52 metapaths, driven by
  Playwright.
  - *positive*: the worked example returns a ranked table, identical on rerun.
  - *negative*: genes absent from Hetionet give the existing error.
  - *figures*: peak RSS and query time for 18 and 200 genes, from a committed
    table.
- **audit** and **summarise**: as in `analytical-null-md`; `summary.ipynb` is
  gate 2.

## Expected result

- Full 52-metapath query on a 24 GB laptop with peak RSS <= 1.5 GB.
- Query time <= 1 s for 18 genes and <= 3 s for 200 genes.
- z matches the matrix-based adapter within the validate tolerances.

## Out of scope

The Postgres loader (DBMI), non-BP targets and non-gene sources, the batch
pipeline's nulls.
