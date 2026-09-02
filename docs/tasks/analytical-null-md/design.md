# Analytical null integration — the query path

**Task** `analytical-null-md` · **Status** OPEN · **Branch** `analytical-null-md` (cut from `upstream/main`, 106cb72)
**Steps** promote, adapt, wire, verify, audit, summarise.
**Validated by** [`docs/tasks/capacity-hurdle-adaptive-null/`](https://github.com/lagillenwater/multi-dwpc/tree/fix/random-null-stratified-srswor/docs/tasks/capacity-hurdle-adaptive-null) on branch `fix/random-null-stratified-srswor`, evidence at `output/tier0_capacity_hurdle/` (commit 3b0bdee).
**As of** 2026-09-01. Decision lineage in [`decisions.md`](decisions.md).

## What this task establishes

The web query's null is updated and validated. The current `query_metapath_z` scores a
user's gene set against `b` random same-size gene subsets drawn from the whole
gene universe — stochastic, blind to degree, and the slowest step of every query.
This task replaces that monte carlo null with the exact resampling moments of the same null,
stratified by leave-target-out capacity, so a query returns the same z every
time, with no `b` to tune and no sampling error to hide in (latency, measured
at the verify step, is disk-dominated and roughly unchanged).

## Data

No new data. The query path reads what the repository already
serves: the per-metapath DWPC matrices under `data/dwpc_cache/`
(`dwpc_<metapath>_d0.50.npz`, raw scale; `scripts/prewarm_dwpc_cache.py` is the
warm-up path and a missing matrix is computed on demand, read-only), the
per-metapath `raw_mean` from `data/metapath-dwpc-stats.tsv`, and the node
tables that map the user's Entrez IDs and the target identifier to matrix
positions. The user's gene set arrives at query time and is the only input
that varies.

## What is computed

**The statistic.** For metapath m and target t, `T = mean over the mapped gene
set of arcsinh(DWPC_m(gene, t) / raw_mean)` — the transformed scale
`query_metapath_z` reports today (`get_dwpc_for_pairs(..., transform=True)`).

**The null, and why it is stratified.** The null must answer: is this gene
set's connectivity to *this* target higher than a random gene set's would be?
Therefore, what is "a random gene set".  Stratifying adjust for degree-confounding in the random set by comparing each of the user's genes only against genes of similar
connectivity, so the z-score reports whether these genes reach a
particular target more than genes with the same overall connectivity. 

**How the strata are built.** Per metapath and target:

1. Each gene gets a **capacity**: its total raw DWPC to every *other* target
   of the metapath (`c_g = sum over t' != t of DWPC_m(g, t')` — the matrix
   row sum minus the tested column). This is the gene's generic reach,
   computed with the same damping as the statistic, and it deliberately
   excludes the tested target.
2. Genes with capacity of zero form one stratum. Zero capacity means no
   these genes score zero for every target.
3. Genes with positive capacity are grouped into
   strata of at least 50 genes each, never splitting genes with equal
   capacity across strata. 
4. A null gene set replaces each of the user's genes with a gene drawn from
   the same stratum, without replacement, the user's own genes excluded. The
   mean and standard deviation of `T` over *all* such draws are computed
   exactly by `hetnetex_md.exact_resampling_moments` — no `b`, no seed, no sampling error.
   `z = (T - mu) / sigma` and `p = Phi_bar(z)` are derived locally. This z is calibrated on the validation rows, stratum-matched random sets
   scored z standard-normal with an upper tail of exactly 0.0500
   (CI [0.0449, 0.0554]).
5. Excluding the user's genes can leave a stratum with fewer candidates than
   it must supply; such a stratum merges into its lower-capacity neighbour,
   and every merge is reported on the result.

**Degenerate cases.** A zero-variance null yields a NaN z and p on a kept row,
matching the current `std == 0` behaviour. A metapath whose matrix cannot be
resolved is skipped, as 2026-09-01.

## Interface

- `src/analytical_null.py` gains the query adapter:
  `analytical_gene_set_z(hetmat, metapath, source_idx, target_pos, *,
  min_stratum_size=50)` returning `real_mean, null_mean, null_std, z, p_value,
  n_active_strata, merges` — the shape the superseded configuration-model line
  proved out, backed now by the validated machinery.
- `query_metapath_z` keeps its signature. `b` and `seed` are accepted,
  ignored, and warn (DeprecationWarning); the output gains `p_value`;
  `real_mean_score` and `null_mean_score` stay on the transformed scale;
  rows are still ranked by `effect_size_z`.
- `app.py` stops plumbing `b`/`seed` into the query and names the analytical
  null in its progress text. No other UI change.
- `pyproject.toml` pins the dependency:
  `hetnetex-md @ git+https://github.com/tghosh30/HetNetEx-MD.git@26f5ba85bd9b9886eafebcbc544879e531f5262e`
  — the canonical upstream at the exact commit the validation ran. No
  submodule; the `lagillenwater` fork stays the contribution vehicle only.

## Figures, controls and the reviewer's path

Every step declares a positive control that plants a known answer through the
shipped code, a negative control that must return null, and the figures a
reviewer sees. Both figure rules from the validated practice hold: every
figure is drawn from a table committed in this task's folder, so its value is
recomputable; and where a control exists it appears beside the real data on
shared axes. Figures land in `docs/tasks/analytical-null-md/figures/`, drawn
by the verification run, and the reviewer meets them in `summary.ipynb` in
this order.

- **adapt** — the query backend.
  - *positive*: on a stub hetmat with a planted enrichment (the gene set
    holds the target column's top scorers), the real adapter returns z above
    the 1.65 convention; a pool planted with a known mu and sigma is
    recovered by the moments within float tolerance.
  - *negative*: a stratum-matched draw from the null scores near zero, and
    over many draws the z distribution is standard-normal shaped (the
    validation's calibration result, re-checked in miniature on the stub);
    a zero-variance column yields the NaN row.
  - *figures*: old-against-new z scatter for one real query — the Monte-Carlo
    z at three seeds against the single analytical z per metapath, drawn from
    a committed per-metapath table; the same table's seed-to-seed MC spread
    as horizontal bars, which is the determinism claim made visible.
- **wire** — the app.
  - *positive*: the app's query call returns the ranked frame with `p_value`
    present.
  - *negative*: passing `b`/`seed` warns and changes nothing (two calls,
    identical frames).
  - *figures*: none of its own — the wiring is one call site; its evidence is
    the verify step's end-to-end run.
- **verify** — the app runs, and the change is measured.
  - *positive*: a real query in the running streamlit app returns a ranked
    table; run twice, byte-identical.
  - *negative*: a gene set of IDs absent from Hetionet fails with the
    existing clear error, not a crash.
  - *figures*: wall-clock per query, old null (at its default `b`) against
    new, from a committed timing table; rank agreement between the old and
    new orderings for the example query (scatter of ranks with Spearman rho
    printed), which shows the reviewer how much the honest null reorders, not
    just that it does.
- **audit** — the design read as numbered claims against the landed tree,
  forward and reverse, per the fm-pm-evaluator audit standard adapted to this
  repository: every claim verified with named evidence, departures classified
  aligned / recorded / drift, a fix wave, then a fresh-reader re-audit. The
  audit reads this folder's tables and figures before they are pushed and
  records what it read. It passes before the summary is reviewed and before
  the PR opens.
- **summarise** — `summary.ipynb`, committed without outputs: the hypothesis,
  each step's figures beside the tables they are drawn from, the behaviour
  changes a reviewer must accept (below), the conclusions. **This is gate 2**:
  the human review happens here, and the pull request description opens by
  sending the reader to it.

## Behaviour changes the summary must present

1. **Deterministic.** Identical queries return identical z; `b`/`seed` are
   inert.
2. **Cost no longer buys noise.** Measured on the example query
   (verify step, 2026-09-01): end-to-end latency is disk-dominated and
   roughly unchanged (new 103 s vs old 85 s over 52 metapaths), and the
   warm-matrix null step is slower per metapath (70 ms vs 8.6 ms) — at the
   app's default `b = 20` the old null was fast because it was imprecise: a
   null sd estimated from 20 draws carries roughly 16% relative error, which
   every z inherits. The analytical moments price out at the Monte-Carlo
   size needed for comparable precision (validated 213x / 2,160x against
   B = 1,000 / 10,000), and the `b` knob is gone. The design's original
   order-of-magnitude latency expectation was written against the wrong
   yardstick and is corrected here (decisions.md, 2026-09-01).

## Promoted apparatus

Carried file-for-file from the validation branch's `src/tier0/` namespace,
tests included; nothing is re-typed, and the integration-hardening items from
the validation's final review are done here because this is now production
surface: unit tests for the fallback-merge path and the concrete capacity-key
wiring, a `min_stratum_size <= 0` guard, the damping constant defined once
and imported, no dead parameters.

| Production path | From (validation branch) | Role |
|---|---|---|
| `src/hurdle_adaptive_bins.py` | `src/tier0/hurdle_adaptive_bins.py` | hurdle + value-respecting bins; deficient-stratum merge |
| `src/pool_assembly.py` | `src/tier0/_pool_assembly.py` | per-stratum pools with self-exclusion |
| `src/capacity.py` | `src/tier0/capacity.py` | leave-target-out capacity, cached per metapath |
| `src/hetnetex_md_import.py` | `src/tier0/hetnetex_md_import.py` | the single, policy-carrying import point |
| `src/analytical_null.py` | `src/tier0/analytical_null.py` + new adapter | moments-to-z wrapper; `analytical_gene_set_z` |

## Expected result

- The example query's analytical z correlates with the Monte-Carlo z.
- Seed-to-seed MC spread is visible at default `b`; the analytical value has
  none.
- Query latency is dominated by matrix loads and roughly unchanged; the
  null-computation share no longer scales with any `b`.

## Out of scope

The batch pipeline's Monte-Carlo nulls (`scripts/permutation_null_datasets.py`,
`scripts/random_null_datasets.py`). 
