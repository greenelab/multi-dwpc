#!/usr/bin/env python3
"""Validate direct DWPC concordance against historical API outputs with plots."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

if Path.cwd().name == "scripts":
    REPO_ROOT = Path("..").resolve()
else:
    REPO_ROOT = Path.cwd()

sys.path.insert(0, str(REPO_ROOT / "src"))

from dwpc_direct import HetMat  # noqa: E402
from dwpc_validation import (  # noqa: E402
    allocate_samples_across_metapaths,
    list_api_metapaths,
    load_api_results,
    sample_metapath_concordance,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--api-results",
        default="output/dwpc_com/all_GO_positive_growth/results/res_all_GO_positive_growth_2016_real.csv",
    )
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--output-dir", default="output/dwpc_validation/all_GO_positive_growth")
    parser.add_argument("--n-samples-total", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--metapath-api", default=None)
    parser.add_argument("--list-metapaths", action="store_true")
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Skip HetMat sampling; read concordance_samples.csv + concordance_summary.csv from --output-dir and re-render plots only.",
    )
    return parser.parse_args()


def _safe_slug(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in str(text))


def _plot_metapath_scatter(samples_df: pd.DataFrame, out_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.scatter(
        samples_df["api_dwpc"].astype(float),
        samples_df["direct_dwpc"].astype(float),
        s=14,
        alpha=0.45,
        edgecolors="none",
    )
    combined = pd.concat(
        [samples_df["api_dwpc"].astype(float), samples_df["direct_dwpc"].astype(float)],
        ignore_index=True,
    )
    lower = float(combined.min())
    upper = float(combined.max())
    ax.plot([lower, upper], [lower, upper], linestyle="--", linewidth=1.0, color="#444444")
    ax.set_xlabel("API DWPC")
    ax.set_ylabel("Direct DWPC")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_overall_scatter_density(samples_df: pd.DataFrame, out_path: Path) -> None:
    x = samples_df["api_dwpc"].astype(float).to_numpy()
    y = samples_df["direct_dwpc"].astype(float).to_numpy()
    bins = 45
    x_edges = np.linspace(float(np.min(x)), float(np.max(x)), bins + 1)
    y_edges = np.linspace(float(np.min(y)), float(np.max(y)), bins + 1)
    hist, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])

    x_idx = np.clip(np.digitize(x, x_edges) - 1, 0, bins - 1)
    y_idx = np.clip(np.digitize(y, y_edges) - 1, 0, bins - 1)
    counts = hist[x_idx, y_idx]
    order = np.argsort(counts)
    x = x[order]
    y = y[order]
    counts = counts[order]

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    scatter = ax.scatter(
        x,
        y,
        c=counts,
        s=18,
        cmap="viridis",
        alpha=1.0,
        linewidths=0,
    )
    combined = pd.concat(
        [samples_df["api_dwpc"].astype(float), samples_df["direct_dwpc"].astype(float)],
        ignore_index=True,
    )
    lower = float(combined.min())
    upper = float(combined.max())
    ax.plot([lower, upper], [lower, upper], linestyle="--", linewidth=1.0, color="#444444")
    ax.set_xlabel("API DWPC")
    ax.set_ylabel("Direct DWPC")
    ax.set_title("Direct vs API DWPC concordance")
    ax.grid(alpha=0.15)
    fig.colorbar(scatter, ax=ax, label="Local point density")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_abs_diff_hist(samples_df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.hist(samples_df["abs_diff"].astype(float), bins=50, color="#4c78a8", alpha=0.8)
    ax.set_xlabel("Absolute difference |direct - api|")
    ax.set_ylabel("Count")
    ax.set_title("DWPC absolute difference distribution")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    output_dir = REPO_ROOT / args.output_dir

    if args.plot_only:
        samples_path = output_dir / "concordance_samples.csv"
        summary_path = output_dir / "concordance_summary.csv"
        for p in (samples_path, summary_path):
            if not p.exists():
                raise FileNotFoundError(
                    f"--plot-only requires {p} to exist. Run without --plot-only first."
                )
        all_samples_df = pd.read_csv(samples_path)
        summary_df = pd.read_csv(summary_path)
        print(f"--plot-only: loaded {samples_path} ({len(all_samples_df)} rows) and {summary_path} ({len(summary_df)} rows)")

        # Per-metapath scatters: one per metapath_api group in the samples CSV.
        if "metapath_api" in all_samples_df.columns:
            for metapath_api, group_df in all_samples_df.groupby("metapath_api"):
                metapath_slug = _safe_slug(str(metapath_api))
                _plot_metapath_scatter(
                    group_df,
                    output_dir / f"scatter_{metapath_slug}.pdf",
                    title=f"{metapath_api} direct vs API",
                )
        _plot_overall_scatter_density(all_samples_df, output_dir / "concordance_overall_scatter.pdf")
        print(f"Saved: {output_dir / 'concordance_overall_scatter.pdf'}")
        _plot_abs_diff_hist(all_samples_df, output_dir / "concordance_abs_diff_hist.pdf")
        print(f"Saved: {output_dir / 'concordance_abs_diff_hist.pdf'}")
        return 0

    api_results_path = REPO_ROOT / args.api_results
    data_dir = REPO_ROOT / args.data_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    api_df = load_api_results(api_results_path)
    metapaths = list_api_metapaths(api_df)
    if args.list_metapaths:
        for metapath in metapaths:
            print(metapath)
        return 0

    if args.metapath_api:
        metapaths = [str(args.metapath_api)]
        sample_plan = {str(args.metapath_api): int(args.n_samples_total)}
    else:
        sample_plan = allocate_samples_across_metapaths(
            api_df,
            total_samples=int(args.n_samples_total),
            metapaths=metapaths,
        )

    hetmat = HetMat(data_dir, damping=0.5, use_disk_cache=True, write_disk_cache=False)
    all_samples = []
    summary_rows = []
    for idx, metapath_api in enumerate(sorted(sample_plan)):
        samples_df, summary = sample_metapath_concordance(
            hetmat,
            api_df,
            data_dir=data_dir,
            metapath_api=metapath_api,
            n_samples=int(sample_plan[metapath_api]),
            seed=int(args.seed) + idx,
        )
        summary_rows.append(summary)
        if samples_df.empty:
            continue
        all_samples.append(samples_df)
        metapath_slug = _safe_slug(summary["metapath_api"])
        _plot_metapath_scatter(
            samples_df,
            output_dir / f"scatter_{metapath_slug}.pdf",
            title=f"{summary['metapath_api']} direct vs API",
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(output_dir / "concordance_summary.csv", index=False)
    print(f"Saved: {output_dir / 'concordance_summary.csv'}")
    if sample_plan:
        print(f"Sample plan: {sample_plan}")

    if not all_samples:
        print("No valid DWPC concordance samples were produced.")
        return 0

    all_samples_df = pd.concat(all_samples, ignore_index=True)
    all_samples_df.to_csv(output_dir / "concordance_samples.csv", index=False)
    print(f"Saved: {output_dir / 'concordance_samples.csv'}")

    _plot_overall_scatter_density(all_samples_df, output_dir / "concordance_overall_scatter.pdf")
    print(f"Saved: {output_dir / 'concordance_overall_scatter.pdf'}")
    _plot_abs_diff_hist(all_samples_df, output_dir / "concordance_abs_diff_hist.pdf")
    print(f"Saved: {output_dir / 'concordance_abs_diff_hist.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
