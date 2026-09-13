#!/usr/bin/env python3
"""Merge the per-metapath parts into the null bundle and write its manifest.

    python scripts/finalize_null_bundle.py --parts-dir /scratch/.../parts --bundle-dir data/null_bundle
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from src.dwpc_direct import DEFAULT_DAMPING  # noqa: E402
from src.null_bundle import finalize  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--parts-dir", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, default=REPO_ROOT / "data" / "null_bundle")
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    parser.add_argument("--expected-parts", type=int, default=None,
                        help="Fail unless exactly this many parts exist (52 for G->BP)")
    args = parser.parse_args()
    manifest = finalize(args.parts_dir, args.bundle_dir, args.data_dir, damping=DEFAULT_DAMPING,
                        expected_parts=args.expected_parts)
    print(f"Wrote {args.bundle_dir}: {len(manifest['metapaths'])} metapaths, "
          f"{manifest['n_strata_rows']:,} strata rows")


if __name__ == "__main__":
    main()
