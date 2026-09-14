# Laptop null bundle: decision lineage

**As of** 2026-09-13. Dated entries, oldest first.

## design.md

- **2026-09-13**: Task opened. Testing the full Streamlit app after the
  `fix/app-null-csc-memory` work showed the preload itself cannot fit: all 52
  G->BP matrices are 2.80B nonzeros, 33.6 GB as CSR, on a 24 GB laptop. The
  memory fix removed a 2x regression but not the architecture's floor.
- **2026-09-13**: The bundle doubles as the production feed (user choice), so it
  is tabular and Postgres-loadable rather than a laptop-only file format.
- **2026-09-13**: Agreement target float64, ~1e-12 relative, not bit-identical
  (user choice). This allows summary statistics in place of full columns.
- **2026-09-13**: Design B (per-stratum summaries + on-the-fly query-gene DWPC)
  chosen over design A (on-disk column store, 15-34 GB, not Postgres-friendly),
  gated by a spike.
- **2026-09-13**: Spike (throwaway, not committed) passed.
  - Row-restricted sparse DWPC from edges against the cached matrices, 52
    metapaths, 238 genes: identical zero patterns; max relative difference
    1.6e-14 except 11 of 2.5M GpBPpGpBP entries that are ~1e-18 cancellation
    residues of true zeros (absolute difference <= 6e-17). 225 ms per query for
    18 genes, 1.6 s for 200 genes, ~52 MB of edge matrices.
  - hetmatpy's `dense_threshold=0` densifies adjacency (3.5 GB for G x G); the
    laptop path must keep it sparse (`dense_threshold=1`).
  - Summary-statistic z against `analytical_gene_set_z`, 52 metapaths x 5 targets
    x 3 gene sets (780): 0 stratum mismatches, max z relative difference 3.5e-13,
    one NaN disagreement (a zero-variance pool that the summary computed as
    ~1e-30), which the zero-variance rule addresses.
  - Stratum count is higher than first estimated: median 141 per
    metapath-target, max 355, about 83M rows.
- **2026-09-13**: DWPC values with absolute value <= 1e-15 are treated as 0.
  The GpBPpGpBP matrix holds negative residues (~-3e-19) that could make a
  capacity negative, which `hurdle_adaptive_bins` rejects.
- **2026-09-13**: `min_stratum_size` fixed at 50 by the bundle; a different
  value means regenerating it.

## implementation

- **2026-09-13**: Gate 1 approved ("approved ... Autonomously implement. Do not
  overcomplicate things."). The stale first-pass HetNetEX edits in the main
  checkout were reverted, and no HetNetEX-MD PR is needed: the laptop query
  path no longer imports `hetnetex_md` or `scipy.stats`.
- **2026-09-13**: The query treats gene IDs as a set; a duplicated ID counts
  once. The matrix-based adapter counted duplicates twice in stratum counts
  while excluding them once from pools.
- **2026-09-13**: `GeneSetZResult` moved to `src/summary_null.py` so the query
  path does not import `src.analytical_null` (and through it `hetnetex_md`).
- **2026-09-13**: Measured correction (mechanism). design.md says
  `strata.parquet` is sorted by `(metapath, target_position)`. It is sorted by
  `target_position`, then manifest metapath order, then stratum, so one
  target's 52 slices are contiguous and read in ~4 ms. Arrow cannot sort a
  dictionary column, so finalize relies on a stable sort by target over parts
  concatenated in manifest order.
- **2026-09-13**: Added rule. A query gene's capacity within 1e-10 of its row
  sum is 0 (`CAPACITY_REL_TOL`): on-the-fly DWPC rounding must not move a
  capacity-0 gene out of the hurdle stratum.
- **2026-09-13**: Known limit, measured. The exact binning can put genes whose
  capacities differ only in the last bit (same off-target row, different
  target value) in different strata; no on-the-fly capacity can reproduce
  that ordering. Unit tests use exactly tied rows; the validate step counts
  stratum mismatches on real data.
- **2026-09-13**: Measured correction (mechanism). Alpine's checkout had no
  G->BP `_d0.5` matrices; the build computed all 52 into
  `/scratch/alpine/$USER/multi-dwpc/null_bundle/dwpc_cache` (22 GB). Alpine's
  37 data files (metagraph, stats, node tables, edges) are byte-identical to
  the laptop's, so the bundle passes the laptop's hash check.
- **2026-09-13**: Measured correction (mechanism). The bundle has 97,790,661
  strata rows (design.md: about 83M), 3.07 GB. GpBPpGpBP was the only
  metapath with residues (578 zeroed).
- **2026-09-13**: Measured correction (mechanism). Validation run 1 (jobs
  32527631/32527632) met both pass criteria but had 24 NaN disagreements: 23
  where the reference null std was ~1e-17 (a true zero-variance null given a
  finite z by rounding) and the summary returned NaN, and 1 (GpBP, target
  1669, 164 annotated genes) where the summary kept 148 eps x scale of
  summation error from a 6,176-gene stratum. The zero-variance tolerance
  changed from 64 eps x scale to eps x (64 + n) x scale; validation reran.
- **2026-09-13**: The user granted full Alpine permission mid-task; the
  project's Alpine guard hook now only blocks credential access.

## audit

- **2026-09-13**: Audit drift 13. design.md says intermediates and paths are
  unchanged. Their behaviour is, but `query_intermediates_and_paths` no longer
  takes `hetmat`: it read node positions through HetMat only, and now reads the
  node tables under `repo_root/data`.
- **2026-09-13**: Audit drift 20. The validate step scored 496 targets, not
  500: annotated-count deciles collapse to 9 bins where counts tie, each bin
  draws 56 targets, and bins include both edges, so a target on a shared edge
  can be drawn twice; de-duplication (with the worked example added) leaves
  496. Annotated gene sets are one per target (its annotated
  genes, up to 200, sampled when more) rather than one per size; random sets
  cover sizes 5, 18, 50, 200, 500. Targets with fewer than 5 annotated genes
  get no annotated set, so 360 of the 496 targets have one. The shuffled-strata
  control gives the query genes the capacities of random other genes.
- **2026-09-13**: Audit drift 22. Rerun identity was checked in-process
  (`scripts/verify_null_bundle_app.py`, two `query_metapath_z` calls per gene
  set, frames equal), not by comparing two browser renders: Streamlit draws the
  table on a canvas that Playwright cannot read as text.

## re-audit

- **2026-09-13**: R6. The zero-variance tolerance has a second term,
  n x (64 eps x pool mean)^2, for rounding of the pool mean; design.md states
  only the sum-of-squares term.
- **2026-09-13**: R7. `assign_strata` places a query capacity that rounding
  puts between two strata's ranges in the nearer range (plan.md; not in
  design.md).
- **2026-09-13**: R8. `query_metapath_z` raises `ValueError` for a target type
  other than BP (the bundle is G->BP only) and defaults to the bundle's
  metapaths instead of discovering them from `metapath-dwpc-stats.tsv`.
- **2026-09-13**: R9. The port step's positive control compares the summary null
  with `analytical_gene_set_z` on the same columns, which computes the moments
  with `exact_resampling_moments`, rather than calling the kernel directly: the
  comparison then also covers strata, pools and merges.

## publication

- **2026-09-14**: The null bundle is published on Zenodo as
  [10.5281/zenodo.22752562](https://doi.org/10.5281/zenodo.22752562) (three files
  uploaded unzipped: `strata.parquet`, `row_sums.parquet`, `manifest.json`), so
  users without Alpine access can download it; the README downloads it into
  `data/null_bundle/`. `manifest.json` records build commit `0940ac6`, a
  pre-squash commit; the build and read code are unchanged in the squashed
  commits (`cb2adcd`, `2130fe1`), and the manifest was left as built so its
  checksum matches the upload.

## review

- **2026-09-14**: Copilot review of PR #11, four findings, all verified and
  fixed in one follow-up commit (no change to results; the published bundle
  still loads):
  - The build could read DWPC matrices cached from other data: the cache is
    looked up by file name only, while finalize hashes the current `data/`.
    `--cache-dir` is now namespaced by a fingerprint of the data-file hashes,
    each part records that fingerprint, and finalize refuses parts built from
    other data.
  - `NullBundle` did not check the manifest's `schema_version`, `damping`,
    `min_stratum_size` or `dwpc_zero_tol`; it now refuses a bundle whose
    settings differ from the query code's.
  - A bundle missing `strata.parquet` or `row_sums.parquet` raised an uncaught
    `FileNotFoundError`; it now raises `BundleMismatchError`, which the app
    shows. The error now points to the Zenodo download as well as the rebuild.
  - The drill-down (`query_intermediates_and_paths`) enumerated a repeated gene
    ID twice; gene IDs are de-duplicated there and in the app's input parser.
