# Analytical null integration — implementation plan

**Task** `analytical-null-md` · **Plan for** [`design.md`](design.md) (approved 2026-09-01, a9f7dc6)
**As of** 2026-09-01. Interfaces, invariants and expected outcomes only — code
exists once, in the repository; this plan never embeds it.

Execution: subagent-driven, one implementer and one reviewer per task, ledger
in the plan workspace; the run proceeds through the audit without pausing and
stops at the summary (gate 2). Every python/pytest call runs inside the
`multi_dwpc` conda env; staging is explicit by path; pushes go to `origin`
only.

## Global invariants

- The five promoted modules come file-for-file from branch
  `fix/random-null-stratified-srswor` (`src/tier0/` namespace, commit
  3b0bdee); changes beyond import paths are limited to the three hardening
  items named in Task 1 and are listed in the commit message.
- Only `exact_resampling_moments` is imported from `hetnetex_md`, only via
  `src/hetnetex_md_import.py`; z and p are derived locally everywhere.
- The partition (capacity key, hurdle, adaptive bins at
  `min_stratum_size = 50`, deficiency merge) is a function of (metapath,
  target, graph); the queried gene set enters via self-exclusion and reported
  merges only.
- The query statistic stays on the transformed scale; capacity stays raw.
- No file under `data/` is written; DWPC matrices are read with the existing
  read-only cache semantics.
- Full test suite green before every push.

## Task 1 — promote

**Brings in.** `src/hurdle_adaptive_bins.py`, `src/pool_assembly.py`,
`src/capacity.py`, `src/hetnetex_md_import.py`, `src/analytical_null.py`
(wrapper only at this stage) and their test files under `tests/`, from the
validation branch via `git show <branch>:<path>`; the `pyproject.toml`
dependency pin from the design's Interface section; the pin installed into
the local env.

**Interfaces (unchanged from the validation branch).**
`hurdle_adaptive_bins(keys, min_stratum_size=50) -> ndarray` ·
`merge_deficient_strata(pools, counts) -> (pools, counts, merges)` ·
`pools_from_bins(bin_of_row, real_row_idx, n_bins) -> (pools, counts)` ·
`CapacityProvider(hetmat, damping).capacity(metapath, target_position) -> ndarray` ·
`leave_target_out_capacity(matrix_csc, target_position) -> ndarray` ·
`analytical_null(scores, pools, counts, observed) -> AnalyticalNullResult(mean, var, std, z, p, n_pool, k_total)`.

**Hardening (the only permitted deviations).** (1) `hurdle_adaptive_bins`
raises `ValueError` on `min_stratum_size <= 0`, with a test. (2) The damping
constant is defined once — the promoted modules import the default from
`src/dwpc_direct.py` instead of redefining `PIPELINE_DAMPING`. (3) A negative
test asserts the import-policy module exposes `exact_resampling_moments` and
nothing else from the library.

**Expected outcome.** Promoted suites pass with only import-path edits;
`pip show hetnetex-md` reports the pinned commit's version; full suite green.

## Task 2 — adapt

**Adds to `src/analytical_null.py`.**
`analytical_gene_set_z(hetmat, metapath, source_idx, target_pos, *, min_stratum_size=50) -> GeneSetZResult`
with fields `real_mean, null_mean, null_std, z, p_value, n_active_strata,
merges`. Internals per the design's "How the strata are built": transformed
score column for the statistic, raw capacity via `CapacityProvider`,
partition via `hurdle_adaptive_bins` + `pools_from_bins` +
`merge_deficient_strata`, moments via `analytical_null`.

**Invariants.** One matrix load serves both the score column and the capacity
key. NaN z/p (never an exception) on a zero-variance null. `merges` reports
`(from_stratum, into_stratum)` pairs.

**Expected test outcomes (stub hetmat, small dense matrices).** Planted
enrichment (gene set = the target column's top scorers) yields z > 1.65;
planted mu/sigma recovered within float tolerance; ~200 stratum-matched null
draws give z with |mean| < 0.15 and std within 0.15 of 1; zero-variance
column yields the NaN result; a forced-deficient stratum yields a non-empty
`merges`; `min_stratum_size=0` raises.

## Task 3 — wire

**Modifies.** `src/multi_dwpc_query.py::query_metapath_z` — same signature,
backend swapped to `analytical_gene_set_z`; `b`/`seed` accepted, ignored,
DeprecationWarning; output frame gains `p_value`; ranking unchanged
(`effect_size_z` descending). `app.py` — stops passing `b`/`seed`, spinner
text names the analytical null (the current text interpolates `DEFAULT_B`).

**Expected test outcomes.** Warning fires exactly when `b` or `seed` is
passed, and two calls with and without them return identical frames; frame
carries the design's column set including `p_value`; unresolvable metapaths
are skipped not raised; the existing query-path tests pass with at most
column-set updates (any other existing-test edit is a finding, not a fix).

## Task 4 — verify

**Produces** `docs/tasks/analytical-null-md/tables/` and `figures/`, and
`verification.md` with the exact commands and real output.

- Comparison run on the app's built-in example query (`EXAMPLE_BP_NAME` and
  its example gene set in `app.py`): per-metapath table with the old
  Monte-Carlo z at seeds {42, 43, 44} (`b = DEFAULT_B = 20`) and the
  analytical z. The old implementation is materialized from git history
  (`git show upstream/main:src/multi_dwpc_query.py` loaded as a temporary
  module by the comparison script) — never a maintained second copy.
- Timing table: wall-clock for the same query, old at default `b` against
  new, medians over 3 runs each.
- The three declared figures (design: adapt and verify steps), each drawn
  from its committed table.
- Streamlit: app launched in the env, example query executed; run twice —
  byte-identical results; a gene set of unknown IDs produces the existing
  clear error.

**Expected outcome.** Determinism holds; analytical z generally at or below
MC z; latency drops by an order of magnitude or more; rank agreement reported
as measured, whatever it is.

## Task 5 — audit

`audit.md` per the fm-pm-evaluator standard adapted here: the design read as
numbered claims, forward (every claim reaches the tree) and reverse (every
landed change traces to a claim), each verdicted aligned / recorded / drift
with named evidence; one fix wave for drift; a fresh-reader re-audit confirms.
The auditor and re-auditor are fresh with respect to the work. Passes before
Task 6.

## Task 6 — summarise, then stop

`summary.ipynb`, committed without outputs: the hypothesis, each step's
figures beside their tables in design order, the three behaviour changes,
conclusions. Push, then **stop at gate 2** — Lucas reviews the summary before
the pull request opens. The PR description will open by pointing at it.
