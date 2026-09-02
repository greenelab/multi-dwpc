# Analytical null integration — verify step

**Task** `analytical-null-md` · **Step** verify (plan.md Task 4)
**As of** 2026-09-01. Machine: local dev machine (macOS, 25.7 GB physical RAM,
already carrying ~21 GB from other processes and ~7 GB/8 GB swap in use at
the time of this run — see "Memory-bounded HetMat" below).

All commands below ran in the `multi_dwpc` conda env
(`conda activate multi_dwpc && ...`), from the repo root
`/Users/gillenlu/Repositories/multi-dwpc/.worktrees/analytical-null-md`, on
branch `analytical-null-md`.

## Example query

The app's built-in worked example (`app.py`, `EXAMPLE_GENES` /
`EXAMPLE_BP_ID` / `EXAMPLE_BP_NAME`):

- Target: `GO:0006244`, "pyrimidine nucleotide catabolic process"
- Genes: 18 symbols (`DUT, ENTPD4, MBD4, NEIL1, NEIL2, NT5C, NT5M, NTHL1,
  OGG1, SMUG1, TDG, UNG, DPYD, DPYS, TYMP, UPB1, UPP1, UPP2`) — all 18
  resolved to Entrez IDs against the repo's Gene node table (verify_query.py
  log: "18 genes resolved from 18 tokens").

`HetMat` is constructed the same way `app.py`'s `get_hetmat()` does:
`HetMat(data_dir=DATA_DIR)` with `DATA_DIR = <repo root>/data`. This
worktree has no `data/` of its own (it's gitignored); a symlink
`data -> /Users/gillenlu/Repositories/multi-dwpc/data` was created so the
worktree reaches the same read-only DWPC cache the main checkout uses (same
pattern already used by the sibling `fix-random-null-stratified-srswor`
worktree). The symlink is intentionally **not** committed. At the time of
this run `.gitignore`'s `/data/` entry matched only a directory, not the
symlink (`git check-ignore -v data` exited 1, `git status --short` showed
`?? data`), so it was untracked local machine setup rather than actually
ignored; the audit fix wave that followed added a bare `data` line to
`.gitignore` to close that gap.

## Memory-bounded HetMat (methodology note, read this first)

`HetMat`'s in-memory matrix cache (`_dwpc_cache` / `_dwpc_cache_csc`) never
evicts. There are 52 valid `G -> BP` metapaths in this repo's `dwpc_cache`
(discovered via `discover_source_target_metapaths`); measured directly by
loading each and inspecting `.data.nbytes + .indices.nbytes + .indptr.nbytes`,
their decompressed CSR footprint totals roughly 30+ GB (the single largest,
`GeAeGpBP`, is ~2.46 GB decompressed from a 1.75 GB npz file; total on-disk
size across all 52 is 23.1 GB). `app.py`'s `get_hetmat()` explicitly
preloads **all** of them at first launch ("Preload all G -> BP matrices so
first query skips disk I/O") — a pre-existing characteristic of the app,
unrelated to this task's null swap.

This dev machine had ~2.2 GB physical memory free and ~1.4 GB swap free at
the time of this run (`top -l 1`, `sysctl vm.swapusage`). A first attempt at
the comparison run with a plain `HetMat` was killed by the OS (`exit 137`,
SIGKILL) after accumulating matrices across many metapaths without eviction.

`verify_query.py` therefore wraps `HetMat` in `BoundedHetMat`, defined in the
script: an LRU-1 cache keyed on metapath name. Both `query_metapath_z`
implementations process metapaths one at a time internally (finish one
fully, including all `b=20` MC resamples for the old implementation, before
moving to the next), so this only changes *when* a matrix is read from disk,
never a numeric result. `BoundedHetMat` also drops the redundant CSR copy
once a CSC view has been built for the current metapath (analytical only
needs CSC). Peak observed RSS for the main run was ~2.1 GB (`ps aux`,
sampled during the run), which is what let it complete on this machine
without further OOM kills. **This is a verification-script-only
workaround for this machine's memory headroom — it does not change
`src/dwpc_direct.py`, `src/analytical_null.py`, `src/multi_dwpc_query.py`,
or `app.py`.**

The practical consequence: the end-to-end timing figure below is dominated
by disk I/O (every metapath's matrix is read fresh on every timed run, for
both implementations), not by the null-computation cost alone. A
supplementary warm-matrix measurement (below) isolates the latter.

## Command

```
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py
```

The old implementation is materialized from git history inside the script
(`git show upstream/main:src/multi_dwpc_query.py`, written to a temp file,
imported via `importlib`, and deleted at the end of the run) — never a
hand-copied second implementation.

Wall-clock: the full run (52-metapath comparison across 3 MC seeds + 2
analytical determinism calls + 1 warm-up + 3+3 end-to-end timing runs, all
under `BoundedHetMat`'s forced per-metapath disk reload) took **~24.5
minutes** on this machine (process start to `tables/timing.csv` written:
19:49 -> 20:13:30, file mtimes). The two-metapath warm-matrix supplement
(below) added under a minute more. This matches the "expect minutes not
hours" guidance; it is longer than a properly warmed deployment would take
because of the memory-bounded eviction described above.

## Positive control: ranked table + run-twice determinism

`verify_query.py` calls `query_metapath_z` (the production function,
`src/multi_dwpc_query.py`) directly, twice, against the same warm `HetMat`:

```
new_df_1 = query_metapath_z(gene_ids, EXAMPLE_BP_ID, hetmat=hetmat)
new_df_2 = query_metapath_z(gene_ids, EXAMPLE_BP_ID, hetmat=hetmat)
new_df_1.equals(new_df_2)  # -> True
```

Real output (from the run log):

```
Determinism (direct function, run1 vs run2 byte-identical frames): True
```

Both calls return a 52-row, ranked (`effect_size_z` descending) DataFrame
with columns `metapath, real_mean_score, null_mean_score, null_std_score,
diff, effect_size_z, p_value` — the frame shape `query_metapath_z`'s
docstring promises.

## Negative control: unknown gene IDs

```
query_metapath_z([999999901, 999999902], EXAMPLE_BP_ID, hetmat=hetmat)
```

Real output:

```
Negative control (unknown gene IDs): ValueError: None of the provided gene IDs were found in Hetionet
```

A clean `ValueError`, not a crash — matches the existing behavior
(`source_idx.size == 0` guard in `query_metapath_z`).

Like the positive-control determinism check above, this ran via a direct
`query_metapath_z` call rather than through the live UI, for the same
`get_hetmat()` preload memory reason discussed under "Streamlit" below
(driving it through the app would require clicking "Run query", which
triggers the ~30+ GB matrix preload this machine cannot safely absorb).

## Streamlit: AppTest smoke + fallback determinism check (explicit record)

`streamlit==1.56.0` is installed, and `streamlit.testing.v1.AppTest` is
available:

```
conda activate multi_dwpc && python -c "from streamlit.testing.v1 import AppTest; print('AppTest available')"
# -> AppTest available
```

**AppTest was used** to load `app.py` headless and confirm it renders
without exception, with the example query prefilled by default (matching
`app.py`'s `st.session_state` initialization):

```python
from streamlit.testing.v1 import AppTest
at = AppTest.from_file("app.py", default_timeout=120)
at.run()
print("exception:", list(at.exception))
print("title:", at.title[0].value)
print("target_name_choice:", at.session_state["target_name_choice"])
print("genes_raw (first 40):", at.session_state["genes_raw"][:40])
print("buttons:", [b.label for b in at.button])
```

Real output:

```
exception: []
title: Multi-DWPC query tool
target_name_choice: pyrimidine nucleotide catabolic process
genes_raw (first 40): DUT
ENTPD4
MBD4
NEIL1
NEIL2
NT
buttons: ['Example query', 'Random BP', 'Random genes', 'Run query']
```

**Fallback triggered for the actual query execution** (recording this
explicitly, per the task contract): `get_hetmat()` is only called inside
`app.py`'s `if run:` block, i.e. only after the "Run query" button is
clicked — and that call preloads all 52 `G -> BP` matrices (see "Memory-
bounded HetMat" above, ~30+ GB decompressed). This machine had ~2-3 GB of
headroom at the time of testing, so driving `AppTest` through an actual
`at.button[3].click().run()` (clicking "Run query") would very likely
reproduce the same SIGKILL already observed once in this session for a
comparably-sized load, on a machine other people are actively using for
other work — not a risk worth taking to prove a code path that
`query_metapath_z` (identical code, called directly) already exercises
under a memory-safe wrapper above. **App.py "resists" the full AppTest
flow here for this specific memory-related reason**, so per the task's
documented fallback:

1. Server-serves smoke check — `streamlit run app.py --server.headless true`
   launched briefly, curl'd, killed:

   ```
   conda activate multi_dwpc
   nohup streamlit run app.py --server.headless true --server.port 8765 > _streamlit_serve.log 2>&1 &
   curl -s -o /dev/null -w "%{http_code}" http://localhost:8765   # polled until 200
   kill $SERVER_PID
   ```

   Real output: `HTTP 200`; server log confirmed
   `You can now view your Streamlit app in your browser. Local URL:
   http://localhost:8765`, followed by `Stopping...` after the kill. Process
   confirmed gone via `ps -p $SERVER_PID` (no output). No server left
   running.

2. Determinism check — via the direct `query_metapath_z` calls documented
   under "Positive control" above (`new_df_1.equals(new_df_2) == True`).
   This is exactly the fallback's "run the query function directly twice
   for the determinism check."

## Per-metapath comparison table and figures

Table: [`tables/per_metapath_comparison.csv`](tables/per_metapath_comparison.csv)
(52 rows: `metapath, z_analytical, z_mc_seed42, z_mc_seed43, z_mc_seed44,
z_mc_mean, z_mc_std, z_mc_min, z_mc_max`).

- **Figure 1** — [`figures/old_vs_new_z_scatter.png`](figures/old_vs_new_z_scatter.png):
  analytical z (x) vs MC z at each of the 3 seeds (y, one color per seed),
  with `y = x` for reference. The MC seeds sit systematically above the
  `y = x` line for the highest-z metapaths (MC over-states enrichment
  relative to the degree-stratified analytical null there) and scatter both
  sides of it near z = 0.
- **Figure 2** — [`figures/mc_seed_spread.png`](figures/mc_seed_spread.png):
  per-metapath MC seed range (min-max across seeds 42/43/44, horizontal
  bars) against the single deterministic analytical z (dot) — the
  determinism claim made visible. Several metapaths (`GpPWpGpBP`,
  `GpBPpGpBP`, `GpMFpGpBP`) show MC seed ranges spanning roughly 5-20 z units
  (min-max: 19.34, 15.88, 4.75 respectively) at `b = 20`; the analytical
  value has none by construction.

NaN rows: **6 of 52** analytical z values are NaN
(`GbCuGpBP, GdCbGpBP, GdCuGpBP, GpBP, GuCbGpBP, GuCuGpBP`) — all 6 are
zero-variance-null rows, matching the documented `std == 0 -> NaN z/p`
behavior (a kept row, not an exception). **10 of 52** rows have a NaN in at
least one MC seed; the 6 analytical NaNs are a strict subset of these 10 —
i.e. every metapath the analytical null calls zero-variance, the MC null
also called zero-variance in every seed it hit NaN on, a useful cross-check.
The other 4 (`GbCbGpBP, GiGpBP, GcGpBP, Gr>GpBP`) are NaN in *some* MC seeds
but not all — the MC null's own zero-variance/non-zero-variance classification
is itself seed-dependent at `b = 20` (20 random draws can, by chance, all
land on the same value), which the deterministic analytical null does not
suffer from. `query_metapath_z`'s output frame does not surface
`n_active_strata`/`merges` (those live on `GeneSetZResult`, one level below
the adapter query sees, consistent with what `app.py` itself displays), so
merge counts for this run are not reported here.

## Timing

Table: [`tables/timing.csv`](tables/timing.csv) (end-to-end, 52 metapaths,
under `BoundedHetMat`) and
[`tables/warm_metapath_timing.csv`](tables/warm_metapath_timing.csv)
(warm-matrix, single metapath, null-computation only). Figure:
[`figures/timing_comparison.png`](figures/timing_comparison.png) (two
panels, one per table).

Warm-up is not symmetric between the two timing tables, and it's worth being
precise about it. For the **end-to-end** table, `verify_query.py` issues one
explicit warm-up call (`new_query_metapath_z(...)`, not itself timed)
immediately before the analytical timed loop only; the Monte-Carlo timed
loop that follows has no dedicated warm-up call directly before it — its
cache warmth is incidental, inherited from the earlier seeds 42/43/44
per-metapath comparison runs (same `hetmat`, same 52 metapaths, run just
before the timing section). For the **warm-matrix** table, each of the two
metapaths gets one explicit warm-up call for *each* implementation
(`old_mod.query_metapath_z(...)` then `new_query_metapath_z(...)`) before
its own 5-run timed loop, so that comparison is warm-up-symmetric. This
asymmetry in the end-to-end table is immaterial to the qualitative finding
below (analytical slower than MC at `b=20`): both loops there read from a
`hetmat` that had already touched every one of the 52 metapaths at least
once by the time either timed loop starts, and the disk-I/O-free
warm-matrix table — which is symmetric — corroborates the same direction and
a *larger* relative gap, not a smaller one. Medians below are of 3 runs
(end-to-end) / 5 runs (warm-matrix).

*Provenance note*: the warm-matrix block (both metapaths) was added to
`verify_query.py` after the end-to-end run above had already completed and
written `tables/per_metapath_comparison.csv` / `tables/timing.csv` /
`figures/{old_vs_new_z_scatter,mc_seed_spread}.png`. Rather than re-run the
full ~24.5-minute, 52-metapath comparison a second time to pick up the new
block (which does not touch or depend on the earlier stages), the added
code was executed directly against the same repo state to produce
`tables/warm_metapath_timing.csv` and regenerate
`figures/timing_comparison.png`. The logic is identical to what is now
committed in `verify_query.py::main()`; running the full script end-to-end
from a clean state reproduces both tables and all four figures in one pass.

| Comparison | Analytical (new) | Monte-Carlo b=20 (old) | "Speedup" (old/new) |
|---|---|---|---|
| End-to-end, 52 metapaths (disk-bound, `BoundedHetMat`) | 103.22 s | 85.30 s | **0.83x — new is slower** |
| Warm matrix, `GpBP` (single hop, smallest, ~3.7 MB) | 12.60 ms | 7.46 ms | 0.59x — new is slower |
| Warm matrix, `GiGiGpBP` (3 hops, ~1035 MB) | 69.98 ms | 8.59 ms | 0.12x — new is ~8x slower |

Raw timed runs (seconds):
`new_times = [108.037, 103.217, 102.706]`,
`old_times = [84.933, 85.302, 86.146]` (end-to-end);
`GpBP: new=[12.371,12.951,12.636,12.600,12.371]ms old=[7.552,7.488,7.396,7.443,7.458]ms`;
`GiGiGpBP: new=[69.654,70.031,69.980,69.708,70.022]ms old=[8.590,8.342,8.775,8.755,8.203]ms`.

### This contradicts the design doc's expected "order of magnitude or more" — an honest, surprising finding

At the app's actual `b = 20` (`DEFAULT_B`), the analytical null is
**measurably slower**, not faster, than the Monte-Carlo null it replaces,
end-to-end and in the isolated warm-matrix comparison, and the gap widens
with metapath size:

- **End-to-end** is dominated by disk I/O forced by `BoundedHetMat` (see
  the methodology note above) — this number does not isolate the
  algorithmic cost and should not be read as a null-computation comparison.
- The **warm-matrix** measurements do isolate it (matrix already resident,
  no disk I/O in the timed loop), and they still show `analytical` slower:
  ~1.7x slower on the smallest metapath, ~8.1x slower on a mid-sized one.
  The gap scaling with metapath size points to the cause: `analytical_gene_set_z`
  recomputes `leave_target_out_capacity`'s full-matrix row sum
  (`matrix.sum(axis=1)`, O(nnz)) on **every call**, with no cross-call
  caching (the design's own decisions.md entry for Task 2 notes this choice
  explicitly — one load per evaluation, not routed through
  `CapacityProvider`'s cross-call row-sum cache). The old MC null's cost is
  `O(b * n_genes)` sparse lookups against an already-loaded matrix —
  independent of matrix size and, at `b = 20`, cheap.
- The design doc's validated **214x / 2,145x per-row** speedup
  (`docs/tasks/analytical-null-md/design.md`, "Behaviour changes") was
  measured against `B = 1,000` / `B = 10,000` Monte Carlo draws, where MC's
  cost scales linearly with `B` and eventually dominates the analytical
  null's largely-`B`-independent cost. At the *much* smaller `b = 20` the
  live app actually uses, that crossover has evidently not been reached for
  at least these two metapaths — this task did not locate the break-even
  `B`.

This is reported as measured, per the task's instruction to report timing
"whatever it is." It is a genuine concern for the audit step: the design's
"faster" behavior-change claim, as written, does not hold at the app's
actual default in this measurement, and should be revisited (either
re-scoped to explicitly compare against large-`B` MC, where the recorded
Tier-0 validation number applies, or `analytical_gene_set_z`'s per-call
capacity recomputation should be addressed if end-user query latency at
`b = 20` parity is the intended comparison).

## Rank agreement

Table: [`tables/rank_agreement.csv`](tables/rank_agreement.csv). Figure:
[`figures/rank_agreement.png`](figures/rank_agreement.png).

```
Spearman rho (analytical rank vs MC seed=42 rank), n=44: rho=0.2605, p=0.08771
```

`n = 44` (52 rows minus the 8 rows where either side is NaN). A rho of 0.26
(not significant at alpha=0.05, p=0.088) is weak-to-moderate agreement — the
degree-stratified analytical null reorders the ranking substantially
relative to the naive whole-gene-universe MC null, which is the expected
and desired behavior the design doc calls out ("shows the reviewer how much
the honest null reorders, not just that it does"): the top MC-ranked
metapaths (`GpBPpGpBP`, `GpPWpGpBP`, `GpMFpGpBP`) are also analytically
high, but mid-table ordering diverges considerably (visible as the broad
scatter around, not tight against, the rank_agreement.png `y = x` line).

## Agreement and cost versus B

Added 2026-09-02, after the audit, as Gate-2 presentation evidence at
Lucas's direction (see `decisions.md`'s `## verification.md` entry): does
same-null agreement with the analytical value increase with `B`, and what
is the noise-versus-cost tradeoff of getting there by raising `B` instead of
using the closed form? This is presentation framing for evidence already
implied by the design (the analytical null *is* the `B -> infinity` limit
of the same stratified null); it is not a new hypothesis, and `design.md`'s
declared-figure list is not retroactively edited.

**Method.** For each of the 52 metapaths in the committed
`per_metapath_comparison.csv`, the S1 partition (transformed scores, raw
leave-target-out capacity, `hurdle_adaptive_bins`, `pools_from_bins` on the
source rows, `merge_deficient_strata`) was rebuilt from the public modules
— the same sequence `analytical_gene_set_z` uses internally, stopping short
of the final `analytical_null` moments call. `z_analytical` itself is
**not** recomputed; it is read back from the committed table. A
**same-null** Monte Carlo estimate is then drawn directly against that
partition: for `B` in `{20, 100, 1,000, 10,000}` x seed in `{42, 43, 44}`,
`B` independent stratified sets are drawn (`counts[s]` from `pools[s]`,
without replacement, vectorized per stratum via `numpy.random.default_rng(seed)`
and an argpartition-of-random-keys SRSWOR trick, chunked over `B` to bound
memory for large strata), giving `B` set-means of the transformed scores;
`z_mc = (observed - mean) / std` over those `B` draws. A quick correctness
check against the exact closed form before running the full sweep (metapath
`GiGiGpBP`, `B=2,000` and `B=20,000`, seed 42) gave `z_mc = 4.51` and
`4.46` against the exact `z = 4.4721` — converging, as it must, since the
same-null MC and `exact_resampling_moments` estimate the same quantity by
different means.

Separately, the **old-null reference** is the git-history *unstratified*
implementation (materialized the same way as everywhere else in this file)
at `b = 10,000`, seed 42, one real pass against the live `HetMat` (disk-
bound, like the earlier end-to-end timing) — not derived from the cached
partition, so its wall-time and its answer both carry full fidelity to the
actual old code path, including whatever inefficiency is in it.

**Commands** (chunked, each its own foreground invocation per the
instruction to never background this — the `B=10,000` same-null chunk and
the old-null chunk both ran a few minutes, matching the stated expectation):

```
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-build
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 20
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 100
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 1000
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 10000
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-oldnull
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-analytical-ref
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-finalize
```

Real wall-clock, this machine: `b-sweep-build` (rebuild all 52 partitions
from the live `HetMat`, disk-bound) 108.12s; `b-sweep-run 20` 0.34s total (3
seeds); `100` 1.62s; `1000` 16.05s; `10000` 173.08s (57-58s/seed);
`b-sweep-oldnull` 222.28s (one seed, one pass); `b-sweep-analytical-ref`
(closed form, all 52 metapaths, 5 repeats, no disk I/O) median 6.39ms.

Tables:
[`tables/b_sweep_per_metapath.csv`](tables/b_sweep_per_metapath.csv) (624
rows: 52 metapaths x 4 `B` x 3 seeds, columns `metapath, B, seed, z_mc,
z_analytical, z_oldnull_b10k`) and
[`tables/b_sweep_summary.csv`](tables/b_sweep_summary.csv) (one row per `B`
plus one `analytical` row and one `old_null_b10k` row). Figures:
[`figures/b_sweep_agreement.png`](figures/b_sweep_agreement.png) (rho and
median absolute error vs `B`, log x, seed range shaded, analytical as the
`B -> infinity` reference line, old-null-b10k plotted as a separate marker)
and [`figures/b_sweep_tradeoff.png`](figures/b_sweep_tradeoff.png) (null
wall-time, log y, vs agreement — one point per `B`, plus analytical and
old-null-b10k).

**Real numbers, exactly as measured** (`tables/b_sweep_summary.csv`):

| `B` | rho (mean over 3 seeds) | rho range | median &#124;z_mc − z_analytical&#124; | total wall-time (3 seeds) |
|---|---|---|---|---|
| 20 | 0.9752 | 0.0057 | 0.198 | 0.339 s |
| 100 | 0.9922 | 0.0056 | 0.078 | 1.620 s |
| 1,000 | 0.9979 | 0.0007 | 0.026 | 16.05 s |
| 10,000 | 0.9995 | 0.0004 | 0.007 | 173.08 s |
| analytical (closed form) | 1.0000 (exact) | 0 | 0 | 0.0064 s |
| old-null, b=10,000 (1 seed) | 0.2829 | n/a (1 seed) | 0.660 | 222.28 s |

**(1) Agreement increases with `B`, monotonically, and clearly converges
toward the analytical value**: rho climbs 0.975 -> 0.992 -> 0.998 -> 0.9995
as `B` goes 20 -> 100 -> 1,000 -> 10,000, and median absolute z-error falls
0.198 -> 0.078 -> 0.026 -> 0.007 over the same span — both heading toward
the analytical value's own row (rho = 1, error = 0) exactly as the theory
requires, since the same-null MC and the closed form estimate the identical
quantity (`exact_resampling_moments`' docstring: "the `B -> infinity` limit
of the stratified... null"). Even the *smallest* `B = 20` already gets rho
= 0.975 — the stratification itself, not a large `B`, is doing most of the
agreement work; raising `B` from there buys diminishing, but real,
additional convergence. Reported as measured: same-null MC at `B = 10,000`
does not reach rho = 1 exactly (rho = 0.9995, not 1.0) — finite `B` leaves
a small residual gap, honestly present in the number.

**(2) The old-null series answers "why not just raise b" directly, and the
answer is no** — raising `b` on the *old, unstratified* null does not
converge to the analytical answer at all. At `b = 10,000` (500x the app's
own `b = 20` default, and the same `B` as the same-null MC's best-agreement
point), the old null's rho against the analytical z is **0.283** — *worse*
than the same-null MC gets at its cheapest setting (`B = 20`, rho = 0.975),
and close to the old null's own `b = 20` rank-agreement figure measured
earlier in this document (0.2605, "Rank agreement" section above) — i.e.
raising `b` 500x on the old null moved its agreement with the analytical
answer from 0.2605 to 0.283, essentially nowhere. This is the expected
result, not a surprise: the old null draws uniformly from the *whole* gene
universe, ignoring degree/capacity entirely, so it is consistently
estimating a **different, degree-blind quantity** — more `b` reduces *its
own* sampling noise around *that* quantity, not the distance between that
quantity and the degree-stratified analytical answer. The `b_sweep_agreement.png`
figure makes this visually explicit: the old-null-b10k marker sits far
below the same-null MC's whole convergence trajectory, not on it.

**(3) The cost side of the tradeoff, in one number each**: the analytical
closed form scores all 52 metapaths for one query in **6.4 ms**
(`b-sweep-analytical-ref`, median of 5, no disk I/O — the same inputs
`analytical_gene_set_z` would use, timed in isolation). Same-null MC at
`B = 10,000` (the setting that gets closest to the analytical answer) costs
**173 s** per 3-seed sweep — roughly **9,000x** the closed form's cost for
a result that is still not exactly it (rho = 0.9995, not 1). The old null
at `b = 10,000` costs **222 s**, worse than same-null MC at the same `B`
(disk I/O plus the old implementation's own Python-loop-bound draw
generation — 10,000 sequential `rng.choice` calls, then a further 10,000-
iteration Python loop per metapath for `get_dwpc_for_pairs`), **for a rho
of 0.283**: it is the single worst point on the whole plot, both slowest
and least accurate — genuinely worse than doing nothing extra, i.e. worse
than staying at the app's own `b = 20` default (rho 0.2605, 100x cheaper).
`b_sweep_tradeoff.png` puts this in one panel exactly as asked: the old
null sits top-left (expensive, wrong), the same-null MC curve runs from
lower-left toward upper-right as `B` grows (cheaper-but-noisier to
costlier-but-closer), and the analytical point sits alone in the
bottom-right corner (cheapest and exact) — nothing on the same-null MC
curve, however large `B` gets, reaches that corner without also paying the
`O(B)` cost the closed form skips entirely.

**Honest caveats**: (a) `old_null_b10k`'s row uses a single seed (42), so
no seed-to-seed spread is reported for it (`n/a` in the table above) —
unlike the same-null MC rows, which average over 3 seeds; this matches
what the task asked for ("one pass") but means the old-null point is not
directly seed-range-comparable to the same-null curve. (b) The same-null MC
wall-times and the analytical reference wall-time are measured on the
identical basis (cached S1 partitions, no disk I/O, one full query across
52 metapaths) and are directly comparable; the old-null-b10k wall-time
additionally includes matrix disk I/O (like the original end-to-end timing
section above), so part of its 222s is not "null generation" in the same
narrow sense — it is still the correct number to report, since it is what
the real old code path actually costs end-to-end, and the task asked for
this reference to be run "genuinely" rather than approximated from the
cached scores.

## Full test suite

```
conda activate multi_dwpc && python -m pytest -q
```

```
.................................................................    [100%]
65 passed, 4 subtests passed in 1.10s
```

Green, all 65 collected tests across `tests/test_*.py` (8 files).

## Summary of headline numbers

| Item | Value |
|---|---|
| n metapaths scored (new) | 52 |
| n metapaths scored (old, seed=42) | 52 |
| n NaN (analytical) | 6 / 52 |
| n NaN (any MC seed) | 10 / 52 |
| Determinism (direct call, run1 vs run2) | identical (`.equals() == True`) |
| Negative control | `ValueError`, not a crash |
| AppTest | available, used for headless load + example-query-prefill smoke check; full run-query flow not exercised via AppTest (memory; fallback used, recorded above) |
| Streamlit headless serve | HTTP 200, served and killed cleanly |
| Median latency, end-to-end (disk-bound) | new 103.22 s vs old 85.30 s (new **slower**, 0.83x) |
| Median latency, warm matrix (`GiGiGpBP`) | new 69.98 ms vs old 8.59 ms (new **slower**, 0.12x) |
| Spearman rho (analytical vs MC seed=42 rank) | 0.2605 (p=0.088, n=44) |
| B-sweep: rho, same-null MC at B=20/100/1,000/10,000 | 0.9752 / 0.9922 / 0.9979 / 0.9995 (mean over 3 seeds; monotonically increasing) |
| B-sweep: rho, old (unstratified) null at b=10,000 | 0.2829 (single seed; does not converge) |
| B-sweep: wall-time, analytical (52 metapaths, 1 query) vs same-null MC B=10,000 (3 seeds) vs old-null b=10,000 (1 seed) | 6.4 ms vs 173.08 s (~9,000x) vs 222.28 s (~34,700x) |

## Reproducing

```
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py
```

The B-sweep (Gate-2 evidence, added after the audit — see "Agreement and
cost versus B" above) is a separate, chunked set of stages on the same
script, each its own foreground invocation:

```
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-build
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 20
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 100
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 1000
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-run 10000
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-oldnull
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-analytical-ref
conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py b-sweep-finalize
```

`b-sweep-build` regenerates `_b_sweep_artifacts.pkl` (a transient, uncommitted
cache of the 52 metapaths' S1 partitions); each `b-sweep-run`/`b-sweep-oldnull`/
`b-sweep-analytical-ref` stage writes partial CSVs to a transient, uncommitted
`_b_sweep_partials/` directory (so an interrupted sweep survives — resume by
re-running only the missing `B` stages); `b-sweep-finalize` merges everything
into the two committed tables and two figures, then deletes both transient
working locations.

regenerates `tables/*.csv` and `figures/*.png` from scratch (~25 minutes on
a memory-constrained machine under `BoundedHetMat`; faster on a machine that
can hold all 52 metapaths resident, since disk I/O would then only happen
once per metapath rather than on every timed call).
