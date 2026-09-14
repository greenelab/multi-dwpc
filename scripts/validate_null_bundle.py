#!/usr/bin/env python3
"""Validate the null bundle against the matrix-based adapter (one metapath per task).

For metapath number ``--index`` (manifest order) and a fixed sample of targets
(the worked example plus a sample stratified by annotated-gene count), score
gene sets two ways:

- reference: ``analytical_gene_set_z`` on the full DWPC matrix (zero rule applied)
- laptop path: ``QueryDwpc`` (edges) + bundle strata + ``summary_gene_set_z``

and, as a negative control, the laptop path with the query genes given the
capacities of random other genes (wrong strata). Gene sets: random sets of
5, 18, 50, 200 and 500 genes, and each target's annotated genes (up to 200).
Writes ``validate_<index>.parquet`` under ``--out-dir``.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.analytical_null import analytical_gene_set_z  # noqa: E402
from src.dwpc_direct import DEFAULT_DAMPING, HetMat, transform_dwpc  # noqa: E402
from src.hurdle_adaptive_bins import hurdle_adaptive_bins  # noqa: E402
from src.multi_dwpc_query import _target_position  # noqa: E402
from src.null_bundle import MIN_STRATUM_SIZE, NullBundle, data_fingerprint  # noqa: E402
from src.query_dwpc import QueryDwpc  # noqa: E402
from src.sparse_column import dense_column  # noqa: E402
from src.summary_null import DWPC_ZERO_TOL, assign_strata, query_capacity, summary_gene_set_z  # noqa: E402

RANDOM_SIZES = [5, 18, 50, 200, 500]
ANNOTATED_MAX = 200
N_TARGETS = 500
WORKED_EXAMPLE_TARGET = "GO:0006244"
SEED = 20260913


def sample_targets(data_dir: Path, gpbp, n_targets: int = N_TARGETS) -> np.ndarray:
    """Worked example plus targets stratified by annotated-gene count (10 bins)."""
    rng = np.random.default_rng(SEED)
    counts = np.diff(gpbp.indptr)
    bins = np.unique(np.quantile(counts, np.linspace(0, 1, 11)))
    groups = [np.flatnonzero((counts >= lo) & (counts <= hi)) for lo, hi in zip(bins[:-1], bins[1:])]
    per_group = (n_targets - 1) // len(groups) + 1
    sampled = np.concatenate([rng.choice(g, size=min(per_group, g.size), replace=False) for g in groups])
    example = _target_position(data_dir, "BP", WORKED_EXAMPLE_TARGET)
    return np.unique(np.concatenate([[example], rng.permutation(sampled)[: n_targets - 1]]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-targets", type=int, default=N_TARGETS)
    args = parser.parse_args()

    bundle = NullBundle(args.bundle_dir, args.data_dir)
    query_dwpc = QueryDwpc(args.data_dir)
    cache_dir = args.cache_dir / data_fingerprint(args.data_dir) if args.cache_dir else None  # as the build does
    hetmat = HetMat(data_dir=args.data_dir, cache_dir=cache_dir, write_disk_cache=False)
    metapath = bundle.metapaths[args.index]
    start = time.perf_counter()

    matrix = hetmat.compute_dwpc_matrix(metapath, damping=DEFAULT_DAMPING).tocsr()
    matrix.data[np.abs(matrix.data) <= DWPC_ZERO_TOL] = 0.0
    matrix.eliminate_zeros()
    hetmat._dwpc_cache[(metapath, DEFAULT_DAMPING)] = matrix  # reference sees the zero rule too
    row_sums = hetmat.get_dwpc_row_sums(metapath)
    n_genes = matrix.shape[0]

    gpbp = query_dwpc._weight(query_dwpc._metagraph.metapath_from_abbrev("GpBP")[0]).tocsc()
    targets = sample_targets(args.data_dir, gpbp, args.n_targets)
    rng = np.random.default_rng(SEED)
    random_sets = {f"random_{s}": np.sort(rng.choice(n_genes, size=s, replace=False)) for s in RANDOM_SIZES}
    union = np.unique(np.concatenate(list(random_sets.values())))
    union_rows = query_dwpc.rows(metapath, union)

    records = []
    for target in targets:
        strata = bundle.strata_for_target(int(target))[metapath]
        column = dense_column(matrix, int(target))
        true_bins = hurdle_adaptive_bins(row_sums - column, min_stratum_size=MIN_STRATUM_SIZE)
        annotated = gpbp.indices[gpbp.indptr[target]:gpbp.indptr[target + 1]]
        gene_sets = dict(random_sets)
        if annotated.size >= 5:
            gene_sets["annotated"] = np.sort(rng.choice(annotated, size=min(annotated.size, ANNOTATED_MAX),
                                                        replace=False))
        for name, genes in gene_sets.items():
            if name == "annotated":
                values = query_dwpc.target_values(metapath, genes, int(target))
            else:
                values = np.asarray(union_rows[np.searchsorted(union, genes)][:, [int(target)]].todense()).ravel()
                values[np.abs(values) <= DWPC_ZERO_TOL] = 0.0
            scores = transform_dwpc(values, bundle.raw_mean(metapath))
            capacity = query_capacity(bundle.row_sums(metapath)[genes], values)
            summary = summary_gene_set_z(strata, scores, capacity)
            shuffled_capacity = (row_sums - column)[rng.choice(n_genes, size=genes.size, replace=False)]
            shuffled = summary_gene_set_z(strata, scores, shuffled_capacity)
            reference = analytical_gene_set_z(hetmat, metapath, genes, int(target))
            records.append({
                "metapath": metapath, "target_position": int(target), "gene_set": name,
                "n_genes": int(genes.size),
                "strata_mismatches": int(np.count_nonzero(assign_strata(strata, capacity) != true_bins[genes])),
                "max_value_rel_diff": float(np.max(np.abs(values - column[genes])
                                                   / np.maximum(np.abs(column[genes]), 1e-300))),
                "z_ref": reference.z, "z_summary": summary.z, "z_shuffled": shuffled.z,
                "null_mean_ref": reference.null_mean, "null_mean_summary": summary.null_mean,
                "null_std_ref": reference.null_std, "null_std_summary": summary.null_std,
                "p_ref": reference.p_value, "p_summary": summary.p_value,
                "merges_equal": reference.merges == summary.merges,
            })

    out = pd.DataFrame(records)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out_dir / f"validate_{args.index:02d}.parquet", index=False)
    finite = np.isfinite(out.z_ref) & np.isfinite(out.z_summary)
    rel = (out.z_summary - out.z_ref).abs() / out.z_ref.abs().clip(lower=1e-300)
    print(f"[{args.index}] {metapath}: {len(out)} rows, strata mismatches {out.strata_mismatches.sum()}, "
          f"NaN disagreements {(np.isnan(out.z_ref) != np.isnan(out.z_summary)).sum()}, "
          f"max z rel diff {rel[finite].max():.2e}, {time.perf_counter() - start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
