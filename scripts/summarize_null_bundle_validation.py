#!/usr/bin/env python3
"""Combine per-metapath validation parts into committed tables, figures and a verdict.

    python scripts/summarize_null_bundle_validation.py --in-dir /scratch/.../validate \
        --task-dir docs/tasks/laptop-null-bundle

Pass criteria (design.md, validate step): 0 stratum mismatches; max z relative
difference <= 1e-9 on rows finite in both; NaN disagreements are counted and
tabled. The shuffled-strata control is expected to disagree.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

Z_REL_TOL = 1e-9


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in-dir", type=Path, required=True)
    parser.add_argument("--task-dir", type=Path, required=True)
    args = parser.parse_args()
    tables, figures = args.task_dir / "tables", args.task_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    rows = pd.concat([pd.read_parquet(p) for p in sorted(args.in_dir.glob("validate_*.parquet"))],
                     ignore_index=True)
    finite = np.isfinite(rows.z_ref) & np.isfinite(rows.z_summary)
    rows["z_rel_diff"] = np.where(finite, (rows.z_summary - rows.z_ref).abs() / rows.z_ref.abs().clip(lower=1e-300),
                                  np.nan)
    rows["nan_disagree"] = np.isnan(rows.z_ref) != np.isnan(rows.z_summary)
    rows.to_parquet(tables / "validate_rows.parquet", index=False, compression="zstd")

    by_metapath = rows.groupby("metapath").agg(
        rows=("z_ref", "size"), finite_rows=("z_rel_diff", "count"),
        strata_mismatches=("strata_mismatches", "sum"), nan_disagreements=("nan_disagree", "sum"),
        max_z_rel_diff=("z_rel_diff", "max"), max_value_rel_diff=("max_value_rel_diff", "max"),
        merges_equal=("merges_equal", "all"),
    ).reset_index()
    by_metapath.to_csv(tables / "validate_by_metapath.csv", index=False)
    rows.loc[rows.nan_disagree].to_csv(tables / "validate_nan_disagreements.csv", index=False)

    verdict = {
        "rows": int(len(rows)),
        "metapaths": int(rows.metapath.nunique()),
        "targets": int(rows.target_position.nunique()),
        "strata_mismatches": int(rows.strata_mismatches.sum()),
        "max_z_rel_diff": float(rows.z_rel_diff.max()),
        "nan_disagreements": int(rows.nan_disagree.sum()),
        "merges_equal_all": bool(rows.merges_equal.all()),
        "shuffled_control_median_abs_dz": float((rows.z_shuffled - rows.z_ref).abs().median()),
    }
    verdict["pass"] = verdict["strata_mismatches"] == 0 and verdict["max_z_rel_diff"] <= Z_REL_TOL
    (tables / "validate_verdict.json").write_text(json.dumps(verdict, indent=2))

    both = rows.loc[finite]
    control = rows.loc[np.isfinite(rows.z_ref) & np.isfinite(rows.z_shuffled)]
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(control.z_ref, control.z_shuffled, s=3, alpha=0.25, color="#bbbbbb",
               label="control: shuffled strata", rasterized=True)
    ax.scatter(both.z_ref, both.z_summary, s=3, alpha=0.5, color="#1f77b4", label="laptop path", rasterized=True)
    lim = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lim, lim, color="black", lw=0.8)
    ax.set_xscale("symlog"); ax.set_yscale("symlog")
    ax.set_xlabel("z, matrix-based adapter (reference)")
    ax.set_ylabel("z, bundle + on-the-fly DWPC")
    ax.set_title(f"{verdict['rows']:,} rows, {verdict['metapaths']} metapaths")
    ax.legend(loc="upper left", markerscale=4)
    fig.tight_layout()
    fig.savefig(figures / "validate_z_scatter.png", dpi=150)

    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    exponent = np.log10(both.z_rel_diff.clip(lower=1e-17))
    ax.hist(exponent, bins=np.arange(-17.25, -8, 0.5), color="#1f77b4")
    ax.axvline(np.log10(Z_REL_TOL), color="black", ls="--", lw=0.8, label="pass threshold 1e-9")
    ax.set_xlabel("log10 relative z difference (exact agreement plotted at -17)")
    ax.set_ylabel("rows")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "validate_z_rel_diff_hist.png", dpi=150)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
