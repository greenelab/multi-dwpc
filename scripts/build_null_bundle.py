#!/usr/bin/env python3
"""Build one metapath's part of the null bundle (one Slurm array task).

For G->BP metapath number ``--index`` (in ``discover_source_target_metapaths``
order): load or compute its DWPC matrix, set residues of true zeros to 0, and
summarise every target column into strata. Writes ``strata_<index>.parquet``,
``row_sums_<index>.parquet`` and ``part_<index>.json`` under ``--parts-dir``.

    python scripts/build_null_bundle.py --index 0 --parts-dir /scratch/.../parts
    python scripts/build_null_bundle.py --list-metapaths
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.dwpc_direct import DEFAULT_DAMPING, HetMat, get_dwpc_raw_mean, transform_dwpc  # noqa: E402
from src.multi_dwpc_query import discover_source_target_metapaths  # noqa: E402
from src.null_bundle import MIN_STRATUM_SIZE, sha256, write_part  # noqa: E402
from src.summary_null import DWPC_ZERO_TOL, build_strata  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--index", type=int)
    parser.add_argument("--parts-dir", type=Path)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    parser.add_argument("--cache-dir", type=Path, default=None,
                        help="DWPC matrix cache (default: <data-dir>/dwpc_cache)")
    parser.add_argument("--list-metapaths", action="store_true")
    args = parser.parse_args()

    hetmat = HetMat(data_dir=args.data_dir, cache_dir=args.cache_dir)
    metapaths = discover_source_target_metapaths(hetmat, "G", "BP")
    if args.list_metapaths:
        print("\n".join(metapaths))
        return

    metapath = metapaths[args.index]
    start = time.perf_counter()
    matrix = hetmat.compute_dwpc_matrix(metapath, damping=DEFAULT_DAMPING).tocsr()
    n_residues = int(np.count_nonzero(np.abs(matrix.data) <= DWPC_ZERO_TOL))
    matrix.data[np.abs(matrix.data) <= DWPC_ZERO_TOL] = 0.0
    matrix.eliminate_zeros()
    print(f"[{args.index}] {metapath}: {matrix.nnz:,} nonzeros, {n_residues} residues zeroed, "
          f"loaded in {time.perf_counter() - start:.0f}s", flush=True)

    row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    raw_mean = get_dwpc_raw_mean(hetmat.metapath_stats, metapath)
    columns = matrix.tocsc()
    n_genes, n_targets = columns.shape
    rows = {name: [] for name in ["target_position", "stratum", "capacity_min", "capacity_max",
                                  "n_genes", "score_mean", "score_centered_ss"]}
    for target in range(n_targets):
        column = np.zeros(n_genes)
        lo, hi = columns.indptr[target], columns.indptr[target + 1]
        column[columns.indices[lo:hi]] = columns.data[lo:hi]
        capacity = row_sums - column
        if capacity.min() < 0:
            raise ValueError(f"{metapath} target {target}: negative capacity {capacity.min()}")
        table = build_strata(capacity, transform_dwpc(column, raw_mean), MIN_STRATUM_SIZE)
        if table.n_genes.sum() != n_genes:
            raise ValueError(f"{metapath} target {target}: strata cover {table.n_genes.sum()} genes")
        n_strata = table.n_genes.size
        rows["target_position"].append(np.full(n_strata, target))
        rows["stratum"].append(np.arange(n_strata))
        for name in ["capacity_min", "capacity_max", "n_genes", "score_mean", "score_centered_ss"]:
            rows[name].append(getattr(table, name))

    strata = {name: np.concatenate(values) for name, values in rows.items()}
    write_part(args.parts_dir, args.index, metapath, strata, row_sums, raw_mean,
               sha256(hetmat._get_cache_path(metapath, DEFAULT_DAMPING)))
    print(f"[{args.index}] {metapath}: {strata['stratum'].size:,} strata rows over {n_targets} targets "
          f"in {time.perf_counter() - start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
