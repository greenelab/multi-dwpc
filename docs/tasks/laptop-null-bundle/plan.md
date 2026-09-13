# Laptop null bundle: plan

**As of** 2026-09-13. Interfaces, invariants and ordered steps; no embedded source.

## 1. port: `src/summary_null.py` (pure math)

- `DWPC_ZERO_TOL = 1e-15`; `zero_residues(values) -> ndarray` sets |x| <= tol to 0.
- `StratumTable` (arrays per stratum, ascending capacity): `capacity_min`,
  `capacity_max`, `n_genes`, `score_mean`, `score_centered_ss`.
- `build_strata(capacity, scores, min_stratum_size=50) -> StratumTable`: bins with
  `hurdle_adaptive_bins`, sums with `np.bincount`. Used by the Alpine build.
- `assign_strata(table, capacity) -> ndarray`: capacity 0 goes to the hurdle
  stratum when one exists; otherwise the positive stratum whose range contains it,
  or the nearest range when rounding lands in a gap.
- `summary_gene_set_z(table, scores, capacities) -> GeneSetZResult`: pool removal,
  deficient-stratum merge (same order and `merges` as `merge_deficient_strata`),
  SRSWOR mean/variance, zero-variance rule, z, p via `math.erfc`.
- `GeneSetZResult` moves here; `src/analytical_null.py` imports it.
- Tests: moments equal `exact_resampling_moments` on random pools, with merges and
  k >= N (<= 1e-12); merges equal `merge_deficient_strata`; assignment reproduces
  `hurdle_adaptive_bins` for every gene; planted enrichment z > 1.65;
  stratum-matched null draws standard-normal; zero-variance pool NaN, 1e-6 spread
  finite.

## 2. port: `src/query_dwpc.py`

- `QueryDwpc(data_dir)`; `target_values(metapath, gene_positions, target_position)
  -> ndarray`; `UnsupportedMetapathError`.
- Sparse degree-weighted adjacency per metaedge via hetmatpy
  (`dense_threshold=1`), row-restricted products for no_repeats, short_repeat
  (no head segment), BABA (cached diagonals); zero rule applied.
- Tests (skip without `data/dwpc_cache`): matches cached GpBP, GiGpBP, GpBPpGpBP
  rows within 1e-12 relative (plus the residue case); zero-edge gene gives zeros;
  unsupported category raises.

## 3. port: `src/null_bundle.py`

- Writer: `write_part(parts_dir, index, metapath, strata_frame, row_sums,
  matrix_sha256)`; `finalize(bundle_dir, data_dir)` concatenates parts, sorts
  strata by `(target_position, metapath, stratum)` with small row groups, writes
  `row_sums.parquet` and `manifest.json` (hashes of node TSVs, edge files,
  `metapath-dwpc-stats.tsv`; matrix hashes from parts; damping, min stratum size,
  schema version, git commit, row counts).
- Reader: `NullBundle(bundle_dir, data_dir)` verifies hashes;
  `metapaths`, `row_sums(metapath)`, `strata_for_target(target_position)`.
- Tests: tiny synthetic bundle round-trips; a changed local file is refused.

## 4. build (Alpine)

- `scripts/build_null_bundle.py --index I --parts-dir D` (one metapath: zero rule,
  row sums, per-target `build_strata`), `scripts/finalize_null_bundle.py`,
  `hpc/submit_null_bundle.sh` (array 0-51, then finalize).
- Checks: `n_genes` per (metapath, target) sums to 20,945; no negative capacity.

## 5. wire

- `query_metapath_z(gene_ids, target_id, *, bundle, query_dwpc, metapaths=None)`
  and `query_intermediates_and_paths` read node positions from node TSVs; no HetMat.
- `app.py`: cached `NullBundle(data/null_bundle)` and `QueryDwpc`; no matrix preload.

## 6. validate (Alpine), 7. verify (laptop), 8. audit, 9. summarise

As in `design.md`. Tables in `tables/`, figures in `figures/`, `verification.md`,
`audit.md`, `summary.ipynb` (no outputs).
