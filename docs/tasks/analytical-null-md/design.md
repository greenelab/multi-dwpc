# MD analytical null integration (query path) — design

**As of** 2026-09-01.
**Branch** `analytical-null-md`, cut from `upstream/main` (106cb72).
**Status** OPEN — PR 1 of the two-PR plan (code first; the PROCESS/SPEC/STATE
docs system and the Tier 0 write-up are PR 2).

## Context

The web-query null in `src/multi_dwpc_query.py` is a Monte-Carlo draw of `b`
same-size gene subsets from the full, unstratified gene universe — stochastic,
degree-blind, and the slow step of every query. The validation task
(`docs/tasks/capacity-hurdle-adaptive-null/` on branch
`fix/random-null-stratified-srswor`, evidence committed under
`output/tier0_capacity_hurdle/` through commit 3b0bdee) established the
replacement: HetNetEX-MD's exact resampling moments over a
capacity-keyed hurdle+adaptive stratification (strategy S1,
`capacity_hurdle_adaptive`). All spec criteria passed — analytical/MC
concordance at B = 10,000, negative-control calibration at 0.0500
(CI [0.0449, 0.0554]), planted-signal recovery on every row, and a measured
213x/2,160x per-row speedup at B = 1,000/10,000. This task wires that
validated null into the query path and the streamlit app.

The earlier configuration-model integration (5 never-pushed commits on
`pr-pipeline-updates-and-webtool-reqs`) is superseded: its adapter shape
(`analytical_gene_set_z`, deprecation of `b`/`seed`, added `p_value` column)
is reused here; its backend and its HetNetEX submodule are not.

## The statistic and null (binding)

Identical to what the validation validated, applied at query time:

- Observed, per metapath m and target t: `T = mean_{g in S} arcsinh(DWPC_m(g, t) / raw_mean)`
  over the user's mapped gene set S — the same transformed scale
  `query_metapath_z` reports today (`get_dwpc_for_pairs(..., transform=True)`).
- Null: stratified SRSWOR over the S1 partition — key
  `c_g = sum_{t' != t} DWPC_m(g, t')` on the RAW scale (leave-target-out row
  sums), hurdle stratum at exactly 0, value-respecting adaptive bins with
  `min_stratum_size = 50`, deficient strata merged into the lower-key
  neighbor and logged. The partition depends only on (metapath, target,
  graph); the gene set enters via self-exclusion and the logged merges only.
- Moments via `hetnetex_md.exact_resampling_moments` exclusively;
  `z = (T - mu) / sigma`, `p = Phi_bar(z)` derived locally; the library's
  p-values are never surfaced.

## Dependency

`pyproject.toml` gains a pinned direct reference —
`hetnetex-md @ git+https://github.com/tghosh30/HetNetEx-MD.git@26f5ba85bd9b9886eafebcbc544879e531f5262e`
— the canonical upstream at the exact commit the validation ran (the
`lagillenwater` fork is currently identical and remains the contribution
vehicle only). No submodule.

## Code layout

Promoted verbatim (with provenance noted in the commit message) from the
validation branch's `src/tier0/` namespace into production `src/`, tests
included:

| New file | From | Contents |
|---|---|---|
| `src/hurdle_adaptive_bins.py` | `src/tier0/hurdle_adaptive_bins.py` | `hurdle_adaptive_bins`, `merge_deficient_strata` |
| `src/pool_assembly.py` | `src/tier0/_pool_assembly.py` | `pools_from_bins` (self-exclusion) |
| `src/capacity.py` | `src/tier0/capacity.py` | `leave_target_out_capacity`, `CapacityProvider` |
| `src/hetnetex_md_import.py` | `src/tier0/hetnetex_md_import.py` | single monitored import point (only `exact_resampling_moments`) |
| `src/analytical_null.py` | `src/tier0/analytical_null.py` + new adapter | `analytical_null(scores, pools, counts, observed)` wrapper, plus the query adapter below |

New query adapter in `src/analytical_null.py`:

```python
def analytical_gene_set_z(
    hetmat, metapath: str, source_idx: np.ndarray, target_pos: int,
    *, min_stratum_size: int = 50,
) -> GeneSetZResult  # fields: real_mean, null_mean, null_std, z, p_value, n_active_strata, merges
```

It pulls the metapath's transformed score column for the target, builds the
S1 partition from raw capacity keys, applies the deficiency merge, and
scores via `analytical_null`. Zero-variance nulls yield NaN z/p (row kept,
matching upstream's `std == 0` behavior). A missing DWPC matrix is computed
on demand by HetMat (read-only, `write_disk_cache=False` semantics preserved);
`scripts/prewarm_dwpc_cache.py` (already on main) remains the warm-up path.

Integration-hardening items carried from the validation's final review, done
here because the code is now production surface: unit tests for the
fallback-merge path and the concrete capacity key wiring; the
`min_stratum_size <= 0` guard; centralize the damping constant
(`PIPELINE_DAMPING`) in one importable place; no dead parameters.

## Query and app changes

- `src/multi_dwpc_query.py::query_metapath_z`: same signature; `b`/`seed`
  become deprecated-and-ignored (DeprecationWarning, kept for call
  compatibility); backend swapped to `analytical_gene_set_z`; output gains
  `p_value`; `real_mean_score`/`null_mean_score` stay on the transformed
  scale (continuity with upstream's columns); rows still ranked by
  `effect_size_z`.
- `app.py`: drop the `b`/`seed` plumbing into the query call, adjust the
  spinner text to name the analytical null. No other UI change.

## Behavior changes a PR reviewer must know

1. **Deterministic**: identical queries return identical z's; no `b` to tune.
2. **The null got honest, so z's shrink.** With capacity conditioning, the
   fraction of features clearing z >= 1.65 fell 85% -> 47% -> 14.7% across
   the null's three generations on the validation rows
   (`output/tier0_capacity_hurdle/pass_rates.csv` on the validation branch).
   Ranking still works identically; any downstream consumer thresholding on
   z will select fewer metapaths, by design — the old passes were inflated
   by degree confounding.
3. **Speed**: per-metapath null cost drops from `b` DWPC resamples to
   sub-millisecond moments; end-to-end query latency is measured before and
   after in this task's verification and recorded there.

## Testing and verification

- Promoted module tests run as-is under `tests/` (renamed imports only).
- New adapter tests (stub hetmat with small dense matrices): z matches a
  direct moments computation; hurdle exactness for zero-capacity genes;
  zero-variance -> NaN row; deprecation warning fires; `p_value` present;
  merge-log populated when a stratum is forced deficient.
- Full suite green locally before every push.
- Streamlit verification: run the app in the `multi_dwpc` env against the
  local cache, execute a real query end-to-end, record commands, sanity of
  ranked output, and old-vs-new wall-clock in
  `docs/tasks/analytical-null-md/verification.md`.

## Out of scope

- Batch pipeline Monte-Carlo nulls (`scripts/permutation_null_datasets.py`,
  `scripts/random_null_datasets.py`) — unchanged.
- The docs system, Tier 0 write-up, and evidence promotion (PR 2).
- Porting the old branch's `benchmark_null_methods.py` (measured the dropped
  configuration-model backend); the validation's committed runtime benchmark
  plus this task's latency measurement carry the speed claim.
- Any change to `query_intermediates_and_paths` beyond what the
  `query_metapath_z` swap forces.

## Decisions

- **2026-09-01** — initial design, from the approved 2026-08-28 brainstorm
  (approach A: fresh branch off upstream/main, pip pin not submodule) plus
  the validation verdict (S1 ships; integration checklist folded in). The
  transformed-scale statistic is retained (upstream continuity and it is the
  validated object); the configuration-model commit's raw-scale switch is
  not carried over.
