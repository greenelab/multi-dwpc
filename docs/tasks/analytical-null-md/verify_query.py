"""Verify step for analytical-null-md: measured comparison, figures, tables.

Compares the new analytical query null (src.multi_dwpc_query.query_metapath_z,
backed by analytical_gene_set_z) against the old Monte-Carlo null, materialized
from git history (upstream/main:src/multi_dwpc_query.py, loaded as a temporary
module -- never a maintained second copy of the old code).

Run with: conda activate multi_dwpc && python docs/tasks/analytical-null-md/verify_query.py

Produces:
  docs/tasks/analytical-null-md/tables/per_metapath_comparison.csv
  docs/tasks/analytical-null-md/tables/timing.csv
  docs/tasks/analytical-null-md/figures/old_vs_new_z_scatter.png
  docs/tasks/analytical-null-md/figures/mc_seed_spread.png
  docs/tasks/analytical-null-md/figures/timing_comparison.png
  docs/tasks/analytical-null-md/figures/rank_agreement.png

Also runs the streamlit determinism/negative-control checks and prints their
results to stdout for verification.md to record.
"""

from __future__ import annotations

import functools
import gc
import importlib.util
import pickle
import subprocess
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

print = functools.partial(print, flush=True)  # noqa: A001 -- see memory note below

REPO_ROOT = Path(__file__).resolve().parents[3]
TASK_DIR = Path(__file__).resolve().parent
TABLES_DIR = TASK_DIR / "tables"
FIGURES_DIR = TASK_DIR / "figures"
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(REPO_ROOT))

DATA_DIR = REPO_ROOT / "data"

# ---------------------------------------------------------------------------
# B-sweep (Gate-2 presentation evidence, added 2026-09-02 at Lucas's request
# -- see decisions.md's "## verification.md" entry). Transient working
# files, cleaned up by finalize_b_sweep(); not committed.
# ---------------------------------------------------------------------------
B_SWEEP_ARTIFACTS_PATH = TASK_DIR / "_b_sweep_artifacts.pkl"
B_SWEEP_PARTIAL_DIR = TASK_DIR / "_b_sweep_partials"
B_SWEEP_SEEDS = [42, 43, 44]
B_SWEEP_VALUES = [20, 100, 1000, 10000]

# ---------------------------------------------------------------------------
# Memory note (see verification.md, "Memory-bounded HetMat" section):
# HetMat's in-memory matrix cache (_dwpc_cache / _dwpc_cache_csc) never
# evicts. There are 52 G->BP metapaths in this repo's dwpc_cache; measured
# directly, their decompressed CSR footprint totals ~32 GB (largest single
# matrix, GeAeGpBP, is ~2.46 GB decompressed). This dev machine has 25.7 GB
# physical RAM and was already using ~21 GB / had ~1.3 GB free swap at the
# time of this run, so holding all 52 metapaths resident (as a bare HetMat
# would, and as app.py's get_hetmat() explicitly does for production) is
# infeasible here and produced a SIGKILL (exit 137) on the first attempt.
# BoundedHetMat below keeps only the single metapath currently being scored
# resident (LRU-1 across both the plain query_metapath_z callers -- old and
# new -- which both process metapaths one at a time internally), and drops
# the redundant CSR copy once a CSC view has been built for it. This is a
# verification-script-only workaround for this machine's memory headroom; it
# does not change the app or library code under test.
# ---------------------------------------------------------------------------


def make_bounded_hetmat_cls(base_cls):
    class BoundedHetMat(base_cls):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._resident_metapath = None

        def _evict_if_new(self, metapath: str) -> None:
            if metapath != self._resident_metapath:
                self._dwpc_cache.clear()
                self._dwpc_cache_csc.clear()
                gc.collect()
                self._resident_metapath = metapath

        def compute_dwpc_matrix(self, metapath, damping=None):
            self._evict_if_new(metapath)
            return super().compute_dwpc_matrix(metapath, damping=damping)

        def compute_dwpc_matrix_csc(self, metapath, damping=None):
            self._evict_if_new(metapath)
            result = super().compute_dwpc_matrix_csc(metapath, damping=damping)
            damp = damping if damping is not None else self.damping
            # analytical_gene_set_z only needs the CSC view; drop the
            # redundant CSR copy compute_dwpc_matrix_csc leaves behind.
            self._dwpc_cache.pop((metapath, damp), None)
            gc.collect()
            return result

    return BoundedHetMat

# ---------------------------------------------------------------------------
# Materialize the old query module from git history (never hand-copied).
# ---------------------------------------------------------------------------


def load_old_query_module() -> "module":
    refs = ["upstream/main", "origin/main", "main"]
    last_err = None
    for ref in refs:
        try:
            old_src = subprocess.run(
                ["git", "show", f"{ref}:src/multi_dwpc_query.py"],
                cwd=REPO_ROOT,
                check=True,
                capture_output=True,
                text=True,
            ).stdout
            break
        except subprocess.CalledProcessError as err:
            last_err = err
    else:
        raise RuntimeError(
            "Could not load base src/multi_dwpc_query.py via git; tried: "
            + ", ".join(refs)
        ) from last_err

    tmp_path = TASK_DIR / "_old_multi_dwpc_query_snapshot.py"
    tmp_path.write_text(old_src)
    spec = importlib.util.spec_from_file_location("_old_multi_dwpc_query", tmp_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, tmp_path


# ---------------------------------------------------------------------------
# The example query (from app.py: EXAMPLE_GENES / EXAMPLE_BP_ID / EXAMPLE_BP_NAME).
# ---------------------------------------------------------------------------

EXAMPLE_GENES_RAW = (
    "DUT\nENTPD4\nMBD4\nNEIL1\nNEIL2\nNT5C\nNT5M\nNTHL1\nOGG1\nSMUG1\nTDG\nUNG\n"
    "DPYD\nDPYS\nTYMP\nUPB1\nUPP1\nUPP2"
)
EXAMPLE_BP_ID = "GO:0006244"
EXAMPLE_BP_NAME = "pyrimidine nucleotide catabolic process"


def parse_gene_input(raw: str, hetmat) -> list[int]:
    from src.path_enumeration import NODE_TYPE_NAMES

    genes = hetmat.get_nodes(NODE_TYPE_NAMES["G"])
    symbol_to_entrez = dict(zip(genes["name"].astype(str), genes["identifier"].astype(int)))
    entrez_set = set(genes["identifier"].astype(int))
    out: list[int] = []
    for tok in raw.replace(",", "\n").split():
        t = tok.strip()
        if not t:
            continue
        if t.lstrip("-").isdigit() and int(t) in entrez_set:
            out.append(int(t))
        elif t in symbol_to_entrez:
            out.append(symbol_to_entrez[t])
    return out


def main() -> None:
    from src.dwpc_direct import HetMat

    BoundedHetMat = make_bounded_hetmat_cls(HetMat)
    print(f"Loading HetMat (memory-bounded, LRU-1 metapath) from {DATA_DIR} ...")
    hetmat = BoundedHetMat(data_dir=DATA_DIR)

    gene_ids = parse_gene_input(EXAMPLE_GENES_RAW, hetmat)
    print(f"Example query: target={EXAMPLE_BP_NAME!r} ({EXAMPLE_BP_ID}), "
          f"{len(gene_ids)} genes resolved from {len(EXAMPLE_GENES_RAW.split())} tokens")

    from src.multi_dwpc_query import query_metapath_z as new_query_metapath_z

    old_mod, old_tmp_path = load_old_query_module()

    # ------------------------------------------------------------------
    # New (analytical) run -- once, then a repeat for the determinism check.
    # ------------------------------------------------------------------
    print("Running new analytical query (run 1) ...")
    t0 = time.perf_counter()
    new_df_1 = new_query_metapath_z(gene_ids, EXAMPLE_BP_ID, hetmat=hetmat)
    new_elapsed_1 = time.perf_counter() - t0

    print("Running new analytical query (run 2, determinism check) ...")
    t0 = time.perf_counter()
    new_df_2 = new_query_metapath_z(gene_ids, EXAMPLE_BP_ID, hetmat=hetmat)
    new_elapsed_2 = time.perf_counter() - t0

    determinism_ok = new_df_1.equals(new_df_2)
    print(f"Determinism (direct function, run1 vs run2 byte-identical frames): {determinism_ok}")

    # ------------------------------------------------------------------
    # Old (Monte-Carlo) runs -- seeds 42/43/44, b=20 (DEFAULT_B).
    # ------------------------------------------------------------------
    old_seeds = [42, 43, 44]
    old_dfs = {}
    old_elapsed = {}
    for seed in old_seeds:
        print(f"Running old Monte-Carlo query (seed={seed}, b={old_mod.DEFAULT_B}) ...")
        t0 = time.perf_counter()
        old_dfs[seed] = old_mod.query_metapath_z(
            gene_ids, EXAMPLE_BP_ID, b=old_mod.DEFAULT_B, seed=seed, hetmat=hetmat
        )
        old_elapsed[seed] = time.perf_counter() - t0

    # ------------------------------------------------------------------
    # Per-metapath comparison table.
    # ------------------------------------------------------------------
    new_z = new_df_1.set_index("metapath")["effect_size_z"]
    comparison = pd.DataFrame({"metapath": new_z.index, "z_analytical": new_z.values})
    for seed in old_seeds:
        s = old_dfs[seed].set_index("metapath")["effect_size_z"]
        comparison[f"z_mc_seed{seed}"] = comparison["metapath"].map(s)

    mc_cols = [f"z_mc_seed{s}" for s in old_seeds]
    comparison["z_mc_mean"] = comparison[mc_cols].mean(axis=1)
    comparison["z_mc_std"] = comparison[mc_cols].std(axis=1, ddof=1)
    comparison["z_mc_min"] = comparison[mc_cols].min(axis=1)
    comparison["z_mc_max"] = comparison[mc_cols].max(axis=1)

    comparison_path = TABLES_DIR / "per_metapath_comparison.csv"
    comparison.to_csv(comparison_path, index=False)
    print(f"Wrote {comparison_path} ({len(comparison)} rows)")

    n_nan_analytical = comparison["z_analytical"].isna().sum()
    n_nan_mc = comparison[mc_cols].isna().any(axis=1).sum()
    print(f"NaN rows: analytical={n_nan_analytical}, any-MC-seed={n_nan_mc}")

    # ------------------------------------------------------------------
    # Timing table -- median of 3 runs each, warm cache (one warm-up run first).
    # ------------------------------------------------------------------
    print("Warm-up run (new, not timed for the table) ...")
    _ = new_query_metapath_z(gene_ids, EXAMPLE_BP_ID, hetmat=hetmat)

    n_timing_runs = 3
    new_times = []
    for i in range(n_timing_runs):
        t0 = time.perf_counter()
        _ = new_query_metapath_z(gene_ids, EXAMPLE_BP_ID, hetmat=hetmat)
        new_times.append(time.perf_counter() - t0)
    print(f"New (analytical) timed runs: {new_times}")

    old_times = []
    for i in range(n_timing_runs):
        t0 = time.perf_counter()
        _ = old_mod.query_metapath_z(
            gene_ids, EXAMPLE_BP_ID, b=old_mod.DEFAULT_B, seed=42, hetmat=hetmat
        )
        old_times.append(time.perf_counter() - t0)
    print(f"Old (Monte-Carlo, seed=42) timed runs: {old_times}")

    timing = pd.DataFrame(
        {
            "implementation": ["analytical"] * n_timing_runs + ["monte_carlo_b20"] * n_timing_runs,
            "run_index": list(range(n_timing_runs)) * 2,
            "wall_clock_s": new_times + old_times,
        }
    )
    timing_path = TABLES_DIR / "timing.csv"
    timing.to_csv(timing_path, index=False)
    print(f"Wrote {timing_path}")

    median_new = float(np.median(new_times))
    median_old = float(np.median(old_times))
    speedup = median_old / median_new if median_new > 0 else float("inf")
    print(f"Median wall-clock: new={median_new:.4f}s old={median_old:.4f}s speedup={speedup:.1f}x")

    # ------------------------------------------------------------------
    # Supplementary: warm-matrix, single-metapath null-computation timing.
    #
    # The end-to-end timing above uses BoundedHetMat (see the memory note at
    # the top of this file), which evicts every metapath's matrix from
    # memory as soon as the next one is touched -- necessary so a 52-
    # metapath query fits this machine's RAM, but it means EVERY metapath's
    # matrix is read from disk and deserialized on EVERY one of the 6 timed
    # end-to-end runs, for both implementations. That disk/deserialization
    # cost dominates the end-to-end number and can hide (or invert) the
    # actual null-computation speedup. This block isolates the
    # null-computation cost alone for two real, already-scored metapaths of
    # very different size (see per_metapath table + dwpc_cache file sizes):
    # GpBP (single hop, ~3.7 MB on disk -- the smallest) and GiGiGpBP (three
    # hops, ~1035 MB on disk -- a large one). Each is loaded ONCE into an
    # unbounded HetMat and kept resident (no eviction), then both
    # implementations are timed directly against that warm, in-memory
    # matrix -- the apples-to-apples comparison the design doc's validated
    # 213x/2160x per-row figures refer to (measured there against
    # B = 1,000 / 10,000 Monte Carlo, not the app's b = 20 default used
    # here -- see verification.md for why that distinction matters).
    # ------------------------------------------------------------------
    warm_metapaths = ["GpBP", "GiGiGpBP"]
    n_warm_runs = 5
    warm_rows = []
    warm_medians: dict[str, tuple[float, float]] = {}
    warm_times_by_metapath: dict[str, tuple[list[float], list[float]]] = {}
    for warm_metapath in warm_metapaths:
        print(f"Warm-matrix single-metapath timing (metapath={warm_metapath}) ...")
        warm_hetmat = HetMat(data_dir=DATA_DIR, write_disk_cache=False)
        # Warm-up: load the matrix once for each implementation's access
        # pattern (this is a plain HetMat, so csr and csc simply stay
        # resident together -- no eviction).
        _ = old_mod.query_metapath_z(
            gene_ids, EXAMPLE_BP_ID, b=old_mod.DEFAULT_B, seed=42,
            metapaths=[warm_metapath], hetmat=warm_hetmat,
        )
        _ = new_query_metapath_z(gene_ids, EXAMPLE_BP_ID, metapaths=[warm_metapath], hetmat=warm_hetmat)

        warm_new_times = []
        for _ in range(n_warm_runs):
            t0 = time.perf_counter()
            _ = new_query_metapath_z(
                gene_ids, EXAMPLE_BP_ID, metapaths=[warm_metapath], hetmat=warm_hetmat
            )
            warm_new_times.append(time.perf_counter() - t0)
        warm_old_times = []
        for _ in range(n_warm_runs):
            t0 = time.perf_counter()
            _ = old_mod.query_metapath_z(
                gene_ids, EXAMPLE_BP_ID, b=old_mod.DEFAULT_B, seed=42,
                metapaths=[warm_metapath], hetmat=warm_hetmat,
            )
            warm_old_times.append(time.perf_counter() - t0)
        print(f"Warm-matrix new (analytical) timed runs ({warm_metapath}): {warm_new_times}")
        print(f"Warm-matrix old (Monte-Carlo, seed=42) timed runs ({warm_metapath}): {warm_old_times}")

        for i, t in enumerate(warm_new_times):
            warm_rows.append({"metapath": warm_metapath, "implementation": "analytical",
                               "run_index": i, "wall_clock_s": t})
        for i, t in enumerate(warm_old_times):
            warm_rows.append({"metapath": warm_metapath, "implementation": "monte_carlo_b20",
                               "run_index": i, "wall_clock_s": t})

        median_warm_new = float(np.median(warm_new_times))
        median_warm_old = float(np.median(warm_old_times))
        warm_speedup = median_warm_old / median_warm_new if median_warm_new > 0 else float("inf")
        print(f"Warm-matrix median wall-clock ({warm_metapath}): "
              f"new={median_warm_new*1000:.3f}ms old={median_warm_old*1000:.3f}ms "
              f"speedup={warm_speedup:.2f}x")
        warm_medians[warm_metapath] = (median_warm_new, median_warm_old)
        warm_times_by_metapath[warm_metapath] = (warm_new_times, warm_old_times)
        del warm_hetmat
        gc.collect()

    warm_timing = pd.DataFrame(warm_rows)
    warm_timing_path = TABLES_DIR / "warm_metapath_timing.csv"
    warm_timing.to_csv(warm_timing_path, index=False)
    print(f"Wrote {warm_timing_path}")

    # The large metapath is the more informative single number for the
    # headline (closer to a "typical" scored metapath than the trivial
    # single-hop case); keep both for the figure and verification.md.
    median_warm_new, median_warm_old = warm_medians["GiGiGpBP"]
    warm_speedup = median_warm_old / median_warm_new if median_warm_new > 0 else float("inf")
    warm_metapath = "GiGiGpBP"

    # ------------------------------------------------------------------
    # Rank agreement (Spearman rho): new z ranks vs old (seed=42) z ranks.
    # ------------------------------------------------------------------
    rank_df = comparison.dropna(subset=["z_analytical", "z_mc_seed42"]).copy()
    rank_df["rank_analytical"] = rank_df["z_analytical"].rank(ascending=False)
    rank_df["rank_mc_seed42"] = rank_df["z_mc_seed42"].rank(ascending=False)
    rho, pval = spearmanr(rank_df["rank_analytical"], rank_df["rank_mc_seed42"])
    print(f"Spearman rho (analytical rank vs MC seed=42 rank), n={len(rank_df)}: "
          f"rho={rho:.4f}, p={pval:.4g}")

    rank_path = TABLES_DIR / "rank_agreement.csv"
    rank_df.to_csv(rank_path, index=False)
    print(f"Wrote {rank_path}")

    # ------------------------------------------------------------------
    # Negative control: unknown gene IDs -> ValueError, not a crash.
    # ------------------------------------------------------------------
    unknown_ids = [999999901, 999999902]
    try:
        new_query_metapath_z(unknown_ids, EXAMPLE_BP_ID, hetmat=hetmat)
        neg_control_result = "NO ERROR RAISED (unexpected)"
    except ValueError as exc:
        neg_control_result = f"ValueError: {exc}"
    except Exception as exc:  # noqa: BLE001
        neg_control_result = f"UNEXPECTED {type(exc).__name__}: {exc}"
    print(f"Negative control (unknown gene IDs): {neg_control_result}")

    # ------------------------------------------------------------------
    # Figures.
    # ------------------------------------------------------------------

    # 1. Old-vs-new z scatter with the MC seed values.
    fig, ax = plt.subplots(figsize=(6, 6))
    for seed, color in zip(old_seeds, ["#74c476", "#fd8d3c", "#9ecae1"]):
        ax.scatter(
            comparison["z_analytical"], comparison[f"z_mc_seed{seed}"],
            label=f"MC seed {seed}", alpha=0.7, color=color, s=30,
        )
    lims = [
        np.nanmin(comparison[["z_analytical"] + mc_cols].values),
        np.nanmax(comparison[["z_analytical"] + mc_cols].values),
    ]
    ax.plot(lims, lims, "k--", linewidth=1, label="y = x")
    ax.set_xlabel("Analytical z")
    ax.set_ylabel("Monte-Carlo z (b=20)")
    ax.set_title("Old (MC) vs new (analytical) z, per metapath\nexample query: "
                  f"{EXAMPLE_BP_NAME}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "old_vs_new_z_scatter.png", dpi=150)
    plt.close(fig)

    # 2. MC seed-to-seed spread as horizontal bars.
    plot_df = comparison.dropna(subset=mc_cols).sort_values("z_mc_mean").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7, max(4, 0.3 * len(plot_df))))
    y = np.arange(len(plot_df))
    ax.hlines(y, plot_df["z_mc_min"], plot_df["z_mc_max"], color="#9ecae1", linewidth=4,
              label="MC seed range (min-max, seeds 42/43/44)")
    ax.scatter(plot_df["z_analytical"], y, color="#fb6a4a", zorder=3, s=25,
               label="Analytical z (deterministic)")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["metapath"], fontsize=6)
    ax.set_xlabel("z")
    ax.set_title("Monte-Carlo seed-to-seed spread vs deterministic analytical z")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "mc_seed_spread.png", dpi=150)
    plt.close(fig)

    # 3. Timing comparison -- two panels: end-to-end (52 metapaths, disk-
    # bound under the memory-bounded HetMat) and warm-matrix single-metapath
    # null-computation-only (see the supplementary timing block above).
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
    impls = ["analytical", "monte_carlo_b20"]

    medians = [median_new, median_old]
    all_pts = [new_times, old_times]
    ax1.bar(impls, medians, color=["#74c476", "#fb6a4a"], alpha=0.7)
    for i, pts in enumerate(all_pts):
        ax1.scatter([i] * len(pts), pts, color="black", zorder=3, s=20)
    ax1.set_ylabel("Wall-clock seconds")
    ax1.set_title(f"End-to-end, 52 metapaths\n(disk-bound; speedup = {speedup:.2f}x)")

    warm_new_times, warm_old_times = warm_times_by_metapath[warm_metapath]
    warm_medians_ms = [median_warm_new * 1000, median_warm_old * 1000]
    warm_all_pts_ms = [[t * 1000 for t in warm_new_times], [t * 1000 for t in warm_old_times]]
    ax2.bar(impls, warm_medians_ms, color=["#74c476", "#fb6a4a"], alpha=0.7)
    for i, pts in enumerate(warm_all_pts_ms):
        ax2.scatter([i] * len(pts), pts, color="black", zorder=3, s=20)
    ax2.set_ylabel("Wall-clock milliseconds")
    ax2.set_title(f"Warm matrix, 1 metapath ({warm_metapath})\n"
                  f"(null-computation only; speedup = {warm_speedup:.1f}x)")

    fig.suptitle("Query latency: analytical vs Monte-Carlo (b=20)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "timing_comparison.png", dpi=150)
    plt.close(fig)

    # 4. Rank-agreement scatter with Spearman rho printed.
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(rank_df["rank_analytical"], rank_df["rank_mc_seed42"], alpha=0.7, s=30,
               color="#9ecae1")
    lims = [1, max(rank_df["rank_analytical"].max(), rank_df["rank_mc_seed42"].max())]
    ax.plot(lims, lims, "k--", linewidth=1)
    ax.set_xlabel("Rank (analytical z, 1 = highest)")
    ax.set_ylabel("Rank (MC z, seed=42, 1 = highest)")
    ax.set_title(f"Rank agreement: analytical vs MC (seed=42)\nSpearman rho = {rho:.3f} "
                 f"(n={len(rank_df)}, p={pval:.3g})")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "rank_agreement.png", dpi=150)
    plt.close(fig)

    print("Figures written to", FIGURES_DIR)

    # Clean up the materialized old-module snapshot (kept out of the commit;
    # it is a transient artifact of this run, not a maintained second copy).
    old_tmp_path.unlink(missing_ok=True)

    # ------------------------------------------------------------------
    # Summary for verification.md.
    # ------------------------------------------------------------------
    print("\n=== SUMMARY ===")
    print(f"n_metapaths_scored (new) = {len(new_df_1)}")
    print(f"n_metapaths_scored (old, seed=42) = {len(old_dfs[42])}")
    print(f"n_nan_analytical = {n_nan_analytical}")
    print(f"determinism_ok (direct function call) = {determinism_ok}")
    print(f"median_latency_new_s (end-to-end, 52 metapaths, disk-bound) = {median_new:.4f}")
    print(f"median_latency_old_s (end-to-end, 52 metapaths, disk-bound) = {median_old:.4f}")
    print(f"speedup (end-to-end, disk-bound) = {speedup:.2f}x")
    print(f"median_latency_new_ms (warm matrix, {warm_metapath}) = {median_warm_new*1000:.4f}")
    print(f"median_latency_old_ms (warm matrix, {warm_metapath}) = {median_warm_old*1000:.4f}")
    print(f"speedup (warm matrix, null-computation only) = {warm_speedup:.1f}x")
    print(f"spearman_rho = {rho:.4f} (p={pval:.4g}, n={len(rank_df)})")
    print(f"negative_control = {neg_control_result}")


def _draw_stratum_sums(rng, pool_scores: np.ndarray, k: int, B: int,
                        max_chunk_bytes: int = 150_000_000) -> np.ndarray:
    """B independent SRSWOR(k, pool) draws' per-draw score sums, vectorized
    over B in memory-bounded chunks (a large stratum's pool can hold
    thousands of genes; a full (B, n) random-key matrix at B=10,000 would
    risk this machine's tight RAM headroom -- see verify_query.py's memory
    note). Chunk size is capped so one chunk's key matrix stays under
    ``max_chunk_bytes``; the draw itself is the standard "argpartition of
    per-element random keys" trick, exact SRSWOR, no approximation.
    """
    n = len(pool_scores)
    if k == 0:
        return np.zeros(B)
    if k == n:
        # Every draw takes the whole pool: deterministic, no need to
        # generate random keys.
        return np.full(B, pool_scores.sum())
    chunk_B = max(1, min(B, max_chunk_bytes // max(n * 8, 1)))
    out = np.empty(B)
    start = 0
    while start < B:
        cb = min(chunk_B, B - start)
        keys = rng.random((cb, n))
        idx = np.argpartition(keys, k - 1, axis=1)[:, :k]
        out[start:start + cb] = pool_scores[idx].sum(axis=1)
        start += cb
    return out


def build_b_sweep_artifacts() -> None:
    """Stage 1 of the B-sweep: rebuild the S1 partition -- transformed
    scores, raw leave-target-out capacity, hurdle_adaptive_bins,
    pools_from_bins(source rows), merge_deficient_strata -- for every
    metapath in the committed per-metapath table, using only the public S1
    modules (matches analytical_gene_set_z's own body; z_analytical itself
    is read back from the committed table at finalize time, never
    recomputed a second way).

    Caches the small per-metapath products (scores/pools/counts/observed --
    kilobytes, not the multi-GB sparse matrices) to disk, so every later
    B-sweep stage (run once per B, each its own foreground invocation) can
    run without ever touching HetMat's disk-backed cache again.
    """
    from src.dwpc_direct import HetMat, DEFAULT_DAMPING, get_dwpc_raw_mean, transform_dwpc
    from src.capacity import leave_target_out_capacity
    from src.hurdle_adaptive_bins import hurdle_adaptive_bins, merge_deficient_strata
    from src.pool_assembly import pools_from_bins
    from src.multi_dwpc_query import _gene_index_map, _target_position

    BoundedHetMat = make_bounded_hetmat_cls(HetMat)
    hetmat = BoundedHetMat(data_dir=DATA_DIR)
    gene_ids = parse_gene_input(EXAMPLE_GENES_RAW, hetmat)
    gmap = _gene_index_map(hetmat)
    source_idx = np.array([gmap[g] for g in gene_ids if g in gmap], dtype=np.int64)
    target_pos = _target_position(hetmat, "BP", EXAMPLE_BP_ID)

    comparison = pd.read_csv(TABLES_DIR / "per_metapath_comparison.csv")
    metapaths = comparison["metapath"].tolist()

    artifacts = {}
    t0 = time.perf_counter()
    for mp in metapaths:
        matrix = hetmat.compute_dwpc_matrix_csc(mp, damping=DEFAULT_DAMPING)
        raw_target_col = np.asarray(matrix[:, target_pos].todense()).ravel()
        raw_mean = get_dwpc_raw_mean(hetmat.metapath_stats, mp)
        scores = transform_dwpc(raw_target_col, raw_mean)
        capacity = leave_target_out_capacity(matrix, target_pos)
        bins = hurdle_adaptive_bins(capacity, min_stratum_size=50)
        n_bins = int(bins.max()) + 1
        pools, counts = pools_from_bins(bins, source_idx, n_bins)
        pools, counts, merges = merge_deficient_strata(pools, counts)
        observed = float(scores[source_idx].mean())
        artifacts[mp] = {
            "scores": scores, "pools": pools, "counts": counts,
            "observed": observed, "K_total": int(sum(counts)), "merges": merges,
        }
        print(f"  built S1 partition: {mp} (K={sum(counts)}, n_strata={len(pools)})", flush=True)
    build_elapsed = time.perf_counter() - t0
    print(f"Stage 1 (build S1 partitions, {len(metapaths)} metapaths): {build_elapsed:.2f}s")

    B_SWEEP_PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "gene_ids": gene_ids, "source_idx": source_idx, "target_pos": target_pos,
        "metapaths": metapaths, "artifacts": artifacts, "build_elapsed_s": build_elapsed,
    }
    with open(B_SWEEP_ARTIFACTS_PATH, "wb") as f:
        pickle.dump(payload, f)
    print(f"Wrote {B_SWEEP_ARTIFACTS_PATH}")


def run_b_sweep_for_B(B: int) -> None:
    """Stage 2 (run once per B, foreground): SAME-NULL Monte Carlo at this
    B, for each of the 3 seeds, over every metapath's cached S1 partition
    (no HetMat, no disk I/O -- pure null-computation timing). Writes a
    partial CSV pair so an interrupted sweep survives (each B stands alone).
    """
    if not B_SWEEP_ARTIFACTS_PATH.exists():
        raise SystemExit(
            f"{B_SWEEP_ARTIFACTS_PATH} not found -- run 'b-sweep-build' first."
        )
    with open(B_SWEEP_ARTIFACTS_PATH, "rb") as f:
        payload = pickle.load(f)
    metapaths = payload["metapaths"]
    artifacts = payload["artifacts"]

    rows = []
    timing_rows = []
    for seed in B_SWEEP_SEEDS:
        rng = np.random.default_rng(seed)
        t0 = time.perf_counter()
        for mp in metapaths:
            a = artifacts[mp]
            scores, pools, counts = a["scores"], a["pools"], a["counts"]
            observed, K = a["observed"], a["K_total"]
            if K <= 0:
                z_mc = float("nan")
            else:
                totals = np.zeros(B, dtype=np.float64)
                for pool, k in zip(pools, counts):
                    if k == 0:
                        continue
                    totals += _draw_stratum_sums(rng, scores[pool], k, B)
                T = totals / K
                mu = float(T.mean())
                sigma = float(T.std(ddof=1))
                z_mc = (observed - mu) / sigma if sigma > 0 else float("nan")
            rows.append({"metapath": mp, "B": B, "seed": seed, "z_mc": z_mc})
        elapsed = time.perf_counter() - t0
        timing_rows.append({"B": B, "seed": seed, "wall_time_s": elapsed,
                             "n_metapaths": len(metapaths)})
        print(f"B={B} seed={seed}: {elapsed:.3f}s over {len(metapaths)} metapaths", flush=True)

    B_SWEEP_PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(B_SWEEP_PARTIAL_DIR / f"rows_B{B}.csv", index=False)
    pd.DataFrame(timing_rows).to_csv(B_SWEEP_PARTIAL_DIR / f"timing_B{B}.csv", index=False)
    print(f"Wrote partials for B={B} to {B_SWEEP_PARTIAL_DIR}")


def run_old_null_b10k() -> None:
    """Stage 3 (foreground, one pass): the git-history UNSTRATIFIED
    implementation at b=10,000, seed=42 -- the "why not just raise b on the
    old null" reference. Genuinely executes the materialized old module
    against the real HetMat (BoundedHetMat, disk-bound, like the original
    end-to-end timing) rather than approximating it from the cached S1
    scores, so this number carries the same disk-I/O-inclusive honesty as
    the earlier old-vs-new end-to-end comparison -- no new assumption is
    introduced just because it's convenient.
    """
    from src.dwpc_direct import HetMat

    BoundedHetMat = make_bounded_hetmat_cls(HetMat)
    hetmat = BoundedHetMat(data_dir=DATA_DIR)
    gene_ids = parse_gene_input(EXAMPLE_GENES_RAW, hetmat)
    old_mod, old_tmp_path = load_old_query_module()

    print("Running old-null (unstratified) b=10,000, seed=42, one pass ...", flush=True)
    t0 = time.perf_counter()
    df = old_mod.query_metapath_z(gene_ids, EXAMPLE_BP_ID, b=10_000, seed=42, hetmat=hetmat)
    elapsed = time.perf_counter() - t0
    print(f"old-null b=10,000 seed=42: {elapsed:.3f}s over {len(df)} metapaths", flush=True)

    B_SWEEP_PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    df[["metapath", "effect_size_z"]].rename(
        columns={"effect_size_z": "z_oldnull_b10k"}
    ).to_csv(B_SWEEP_PARTIAL_DIR / "oldnull_b10k_rows.csv", index=False)
    pd.DataFrame([{"wall_time_s": elapsed, "n_metapaths": len(df)}]).to_csv(
        B_SWEEP_PARTIAL_DIR / "oldnull_b10k_timing.csv", index=False
    )
    old_tmp_path.unlink(missing_ok=True)
    print(f"Wrote old-null-b10k partials to {B_SWEEP_PARTIAL_DIR}")


def run_analytical_reference_timing() -> None:
    """Stage: measure the analytical null's own total wall-time across all
    52 metapaths for one query, from the cached S1 partitions (no HetMat,
    no disk I/O -- the exact same inputs analytical_gene_set_z would have
    used, just the closed-form moments step timed in isolation). This is
    the number the B-sweep tradeoff figure and summary table need on the
    same "one full query, 52 metapaths" basis the same-null-MC wall-times
    are measured on; the earlier warm_metapath_timing.csv (kept, still
    reported in verification.md) only sampled 2 of the 52 metapaths
    individually and isn't on that basis.
    """
    from src.analytical_null import analytical_null

    if not B_SWEEP_ARTIFACTS_PATH.exists():
        raise SystemExit(
            f"{B_SWEEP_ARTIFACTS_PATH} not found -- run 'b-sweep-build' first."
        )
    with open(B_SWEEP_ARTIFACTS_PATH, "rb") as f:
        payload = pickle.load(f)
    metapaths = payload["metapaths"]
    artifacts = payload["artifacts"]

    n_repeats = 5
    times = []
    for _ in range(n_repeats):
        t0 = time.perf_counter()
        for mp in metapaths:
            a = artifacts[mp]
            analytical_null(a["scores"], a["pools"], a["counts"], a["observed"])
        times.append(time.perf_counter() - t0)
    print(f"analytical, {len(metapaths)} metapaths, {n_repeats} repeats: {times}", flush=True)

    B_SWEEP_PARTIAL_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"run_index": range(n_repeats), "wall_time_s": times,
                  "n_metapaths": len(metapaths)}).to_csv(
        B_SWEEP_PARTIAL_DIR / "analytical_reference_timing.csv", index=False
    )
    print(f"Wrote analytical reference timing to {B_SWEEP_PARTIAL_DIR}")


def finalize_b_sweep() -> None:
    """Stage 4: merge all partials into the two committed tables and draw
    the two declared figures. Requires build + all 4 B's + oldnull-b10k +
    analytical-reference-timing partials to already exist.
    """
    from scipy.stats import spearmanr as _spearmanr

    comparison = pd.read_csv(TABLES_DIR / "per_metapath_comparison.csv")
    z_analytical_by_mp = comparison.set_index("metapath")["z_analytical"]

    missing = [B for B in B_SWEEP_VALUES if not (B_SWEEP_PARTIAL_DIR / f"rows_B{B}.csv").exists()]
    if missing or not (B_SWEEP_PARTIAL_DIR / "oldnull_b10k_rows.csv").exists():
        raise SystemExit(f"Missing B-sweep partials in {B_SWEEP_PARTIAL_DIR}: "
                          f"B={missing}, oldnull_b10k={not (B_SWEEP_PARTIAL_DIR / 'oldnull_b10k_rows.csv').exists()}")

    # ------------------------------------------------------------------
    # b_sweep_per_metapath.csv
    # ------------------------------------------------------------------
    row_frames = [pd.read_csv(B_SWEEP_PARTIAL_DIR / f"rows_B{B}.csv") for B in B_SWEEP_VALUES]
    per_mp = pd.concat(row_frames, ignore_index=True)
    per_mp["z_analytical"] = per_mp["metapath"].map(z_analytical_by_mp)

    oldnull_rows = pd.read_csv(B_SWEEP_PARTIAL_DIR / "oldnull_b10k_rows.csv")
    oldnull_by_mp = oldnull_rows.set_index("metapath")["z_oldnull_b10k"]
    per_mp["z_oldnull_b10k"] = per_mp["metapath"].map(oldnull_by_mp)

    per_mp_path = TABLES_DIR / "b_sweep_per_metapath.csv"
    per_mp.to_csv(per_mp_path, index=False)
    print(f"Wrote {per_mp_path} ({len(per_mp)} rows)")

    # ------------------------------------------------------------------
    # b_sweep_summary.csv
    # ------------------------------------------------------------------
    timing_frames = [pd.read_csv(B_SWEEP_PARTIAL_DIR / f"timing_B{B}.csv") for B in B_SWEEP_VALUES]
    timing = pd.concat(timing_frames, ignore_index=True)

    summary_rows = []
    for B in B_SWEEP_VALUES:
        sub = per_mp[per_mp["B"] == B]
        rhos = {}
        for seed in B_SWEEP_SEEDS:
            s = sub[sub["seed"] == seed].dropna(subset=["z_mc", "z_analytical"])
            rho_val = _spearmanr(s["z_mc"], s["z_analytical"])[0] if len(s) >= 2 else float("nan")
            rhos[seed] = rho_val
        rho_vals = list(rhos.values())
        dz = (sub["z_mc"] - sub["z_analytical"]).abs().dropna()
        pivot = sub.pivot(index="metapath", columns="seed", values="z_mc")
        seed_spread = (pivot.max(axis=1) - pivot.min(axis=1)).dropna()
        total_wt = timing.loc[timing["B"] == B, "wall_time_s"].sum()
        summary_rows.append({
            "row_type": "same_null_mc", "B": B,
            "rho_seed42": rhos[42], "rho_seed43": rhos[43], "rho_seed44": rhos[44],
            "rho_mean": float(np.mean(rho_vals)), "rho_range": float(np.ptp(rho_vals)),
            "median_abs_dz": float(dz.median()) if len(dz) else float("nan"),
            "median_seed_spread": float(seed_spread.median()) if len(seed_spread) else float("nan"),
            "total_wall_time_s": float(total_wt), "n_metapaths": len(sub["metapath"].unique()),
        })

    # analytical reference row
    ref_timing = pd.read_csv(B_SWEEP_PARTIAL_DIR / "analytical_reference_timing.csv")
    summary_rows.append({
        "row_type": "analytical", "B": float("nan"),
        "rho_seed42": 1.0, "rho_seed43": 1.0, "rho_seed44": 1.0,
        "rho_mean": 1.0, "rho_range": 0.0,
        "median_abs_dz": 0.0, "median_seed_spread": 0.0,
        "total_wall_time_s": float(ref_timing["wall_time_s"].median()),
        "n_metapaths": int(ref_timing["n_metapaths"].iloc[0]),
    })

    # old-null-b10k reference row (single seed=42 only -- no cross-seed spread)
    oldnull_timing = pd.read_csv(B_SWEEP_PARTIAL_DIR / "oldnull_b10k_timing.csv")
    old_sub = per_mp[["metapath", "z_oldnull_b10k", "z_analytical"]].drop_duplicates("metapath").dropna()
    rho_old = _spearmanr(old_sub["z_oldnull_b10k"], old_sub["z_analytical"])[0] if len(old_sub) >= 2 else float("nan")
    dz_old = (old_sub["z_oldnull_b10k"] - old_sub["z_analytical"]).abs()
    summary_rows.append({
        "row_type": "old_null_b10k", "B": 10_000,
        "rho_seed42": rho_old, "rho_seed43": float("nan"), "rho_seed44": float("nan"),
        "rho_mean": rho_old, "rho_range": float("nan"),
        "median_abs_dz": float(dz_old.median()) if len(dz_old) else float("nan"),
        "median_seed_spread": float("nan"),
        "total_wall_time_s": float(oldnull_timing["wall_time_s"].iloc[0]),
        "n_metapaths": int(oldnull_timing["n_metapaths"].iloc[0]),
    })

    summary = pd.DataFrame(summary_rows)
    summary_path = TABLES_DIR / "b_sweep_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"Wrote {summary_path}")
    print(summary.to_string(index=False))

    # ------------------------------------------------------------------
    # Figure: b_sweep_agreement.png
    # ------------------------------------------------------------------
    mc_summary = summary[summary["row_type"] == "same_null_mc"].sort_values("B")
    old_row = summary[summary["row_type"] == "old_null_b10k"].iloc[0]

    fig, (axr, axd) = plt.subplots(1, 2, figsize=(11, 4.5))

    axr.plot(mc_summary["B"], mc_summary["rho_mean"], "o-", color="#74c476",
             label="same-null MC (mean over 3 seeds)")
    axr.fill_between(mc_summary["B"],
                      mc_summary["rho_mean"] - mc_summary["rho_range"] / 2,
                      mc_summary["rho_mean"] + mc_summary["rho_range"] / 2,
                      color="#74c476", alpha=0.2, label="seed range")
    axr.axhline(1.0, color="#3182bd", linestyle="--", linewidth=1.5,
                label="analytical (B -> infinity, exact)")
    axr.scatter([old_row["B"] * 1.6], [old_row["rho_mean"]], color="#fb6a4a",
                marker="X", s=120, zorder=5,
                label="old (unstratified) null, b=10,000")
    axr.set_xscale("log")
    axr.margins(x=0.15, y=0.12)
    axr.set_xlabel("B (same-null MC draws)")
    axr.set_ylabel("Spearman rho vs analytical z")
    axr.set_title("Rank agreement vs B")
    axr.legend(fontsize=7, loc="lower left")

    axd.plot(mc_summary["B"], mc_summary["median_abs_dz"], "o-", color="#74c476")
    axd.fill_between(mc_summary["B"], 0,
                      mc_summary["median_seed_spread"],
                      color="#9ecae1", alpha=0.3, label="median seed spread")
    axd.axhline(0.0, color="#3182bd", linestyle="--", linewidth=1.5,
                label="analytical (B -> infinity, exact)")
    axd.scatter([old_row["B"] * 1.6], [old_row["median_abs_dz"]], color="#fb6a4a",
                marker="X", s=120, zorder=5, label="old null, b=10,000")
    axd.set_xscale("log")
    axd.margins(x=0.15, y=0.12)
    axd.set_xlabel("B (same-null MC draws)")
    axd.set_ylabel("median |z_mc - z_analytical|")
    axd.set_title("Absolute error vs B")
    axd.legend(fontsize=7, loc="center right")

    fig.suptitle(f"Agreement with the analytical null vs B, example query "
                 f"({EXAMPLE_BP_NAME})")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "b_sweep_agreement.png", dpi=150)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure: b_sweep_tradeoff.png
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.plot(mc_summary["rho_mean"], mc_summary["total_wall_time_s"], "o-",
            color="#74c476", zorder=3)
    for _, r in mc_summary.iterrows():
        ax.annotate(f"B={int(r['B'])}", (r["rho_mean"], r["total_wall_time_s"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=8)

    analytical_row = summary[summary["row_type"] == "analytical"].iloc[0]
    ax.scatter([analytical_row["rho_mean"]], [analytical_row["total_wall_time_s"]],
               color="#3182bd", marker="*", s=250, zorder=4, label="analytical")
    ax.annotate("analytical", (analytical_row["rho_mean"], analytical_row["total_wall_time_s"]),
                textcoords="offset points", xytext=(8, -4), fontsize=8)

    ax.scatter([old_row["rho_mean"]], [old_row["total_wall_time_s"]],
               color="#fb6a4a", marker="X", s=180, zorder=4, label="old null, b=10,000")
    ax.annotate("old null\nb=10,000", (old_row["rho_mean"], old_row["total_wall_time_s"]),
                textcoords="offset points", xytext=(8, 0), fontsize=8)

    ax.set_yscale("log")
    ax.set_xlabel("Spearman rho vs analytical z (agreement)")
    ax.set_ylabel("Total null wall-time, one query / 52 metapaths (s, log scale)")
    ax.set_title("Cost vs agreement: same-null MC (by B), old null, and analytical")
    ax.legend(loc="center left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "b_sweep_tradeoff.png", dpi=150)
    plt.close(fig)

    print("Wrote figures/b_sweep_agreement.png and figures/b_sweep_tradeoff.png")

    # ------------------------------------------------------------------
    # Clean up transient working files (not committed).
    # ------------------------------------------------------------------
    import shutil
    shutil.rmtree(B_SWEEP_PARTIAL_DIR, ignore_errors=True)
    B_SWEEP_ARTIFACTS_PATH.unlink(missing_ok=True)
    print("Cleaned up transient B-sweep working files.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "b-sweep-build":
        build_b_sweep_artifacts()
    elif len(sys.argv) > 1 and sys.argv[1] == "b-sweep-run":
        run_b_sweep_for_B(int(sys.argv[2]))
    elif len(sys.argv) > 1 and sys.argv[1] == "b-sweep-oldnull":
        run_old_null_b10k()
    elif len(sys.argv) > 1 and sys.argv[1] == "b-sweep-analytical-ref":
        run_analytical_reference_timing()
    elif len(sys.argv) > 1 and sys.argv[1] == "b-sweep-finalize":
        finalize_b_sweep()
    else:
        main()
