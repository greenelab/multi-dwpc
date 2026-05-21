#!/usr/bin/env python3
"""Pre-warm the on-disk DWPC matrix cache for web-tool queries.

Modes:

- local parallel (default):
    python scripts/prewarm_dwpc_cache.py                      # G -> BP, 4 workers
    python scripts/prewarm_dwpc_cache.py --target-type D --workers 8

- list mode (for HPC array driver):
    python scripts/prewarm_dwpc_cache.py --list-metapaths

- single-metapath mode (for HPC array task):
    python scripts/prewarm_dwpc_cache.py --single-metapath GpPWpGpBP
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.dwpc_direct import HetMat  # noqa: E402
from src.multi_dwpc_query import discover_source_target_metapaths  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-type", default="G", help="Source metanode abbrev (default: G)")
    p.add_argument("--target-type", default="BP", help="Target metanode abbrev (default: BP)")
    p.add_argument("--data-dir", default=str(REPO_ROOT / "data"))
    p.add_argument("--workers", type=int, default=4, help="Parallel workers (default: 4)")
    p.add_argument("--damping", type=float, default=0.5)
    p.add_argument("--list-metapaths", action="store_true",
                   help="Print discovered metapaths (one per line) and exit.")
    p.add_argument("--single-metapath", default=None,
                   help="Compute and cache exactly this metapath, then exit (for HPC array tasks).")
    args = p.parse_args()

    hetmat = HetMat(data_dir=args.data_dir, damping=args.damping)

    if args.single_metapath:
        t0 = time.perf_counter()
        hetmat.compute_dwpc_matrix(args.single_metapath, args.damping)
        print(f"Computed {args.single_metapath} in {time.perf_counter() - t0:.1f}s")
        return

    metapaths = discover_source_target_metapaths(hetmat, args.source_type, args.target_type)

    if args.list_metapaths:
        for mp in metapaths:
            print(mp)
        return

    print(f"Discovered {len(metapaths)} {args.source_type} -> {args.target_type} metapaths")
    for mp in metapaths:
        print(f"  {mp}")

    t0 = time.perf_counter()
    hetmat.precompute_matrices(
        metapaths=metapaths,
        damping=args.damping,
        n_workers=args.workers,
        show_progress=True,
        cache_in_memory=False,
    )
    print(f"Done in {time.perf_counter() - t0:.1f}s. Cache dir: {hetmat.cache_dir}")


if __name__ == "__main__":
    main()
