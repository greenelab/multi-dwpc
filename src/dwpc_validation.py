"""Helpers for validating direct DWPC values against historical API outputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from dwpc_direct import HetMat, reverse_metapath_abbrev
from result_normalization import load_neo4j_mappings


NODE_ABBREVIATIONS = ["BP", "CC", "MF", "PW", "SE", "PC", "G", "A", "D", "C", "S"]


def load_api_results(api_results_path: Path | str) -> pd.DataFrame:
    """Load an API result CSV and validate the required columns."""
    api_results_path = Path(api_results_path)
    if not api_results_path.exists():
        raise FileNotFoundError(f"API results not found: {api_results_path}")
    df = pd.read_csv(api_results_path)
    required = {"neo4j_source_id", "neo4j_target_id", "metapath_abbreviation", "dwpc"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"API results {api_results_path} missing required columns: {sorted(missing)}"
        )
    return df


def list_api_metapaths(api_df: pd.DataFrame) -> list[str]:
    """Return sorted API metapath abbreviations present in a result frame."""
    return sorted(api_df["metapath_abbreviation"].dropna().astype(str).unique().tolist())


def metapath_length(metapath: str) -> int:
    """Infer edge length from the metapath abbreviation by counting metanodes."""
    text = str(metapath)
    node_tokens: list[str] = []
    idx = 0
    node_abbrevs = sorted(NODE_ABBREVIATIONS, key=len, reverse=True)
    while idx < len(text):
        matched = False
        for token in node_abbrevs:
            if text.startswith(token, idx):
                node_tokens.append(token)
                idx += len(token)
                matched = True
                break
        if not matched:
            idx += 1
    return max(0, len(node_tokens) - 1)


def allocate_samples_across_metapaths(
    api_df: pd.DataFrame,
    *,
    total_samples: int,
    metapaths: list[str] | None = None,
    allowed_lengths: tuple[int, ...] = (2, 3, 4),
) -> dict[str, int]:
    """
    Allocate a total sample budget across metapaths while covering multiple lengths.

    The allocation is stratified by inferred metapath length, then distributed
    approximately evenly within each length bucket.
    """
    selected = metapaths or list_api_metapaths(api_df)
    if not selected:
        return {}

    total_samples = max(1, int(total_samples))
    allowed_lengths = tuple(int(length) for length in allowed_lengths)
    length_to_metapaths: dict[int, list[str]] = {}
    available_counts: dict[str, int] = {}
    for metapath in selected:
        length = metapath_length(str(metapath))
        if allowed_lengths and length not in allowed_lengths:
            continue
        count = int((api_df["metapath_abbreviation"].astype(str) == str(metapath)).sum())
        if count <= 0:
            continue
        available_counts[str(metapath)] = count
        length_to_metapaths.setdefault(length, []).append(str(metapath))

    if not length_to_metapaths:
        return {}

    allocations = {metapath: 0 for metapath in available_counts}
    lengths = sorted(length_to_metapaths)

    if total_samples < len(available_counts):
        selected_subset: list[str] = []
        per_length_queues = {
            length: sorted(length_to_metapaths[length])
            for length in lengths
        }
        while len(selected_subset) < total_samples:
            progressed = False
            for length in lengths:
                queue = per_length_queues[length]
                if not queue:
                    continue
                selected_subset.append(queue.pop(0))
                progressed = True
                if len(selected_subset) >= total_samples:
                    break
            if not progressed:
                break
        return {metapath: 1 for metapath in selected_subset}

    base_per_length = total_samples // len(lengths)
    extra_length_budget = total_samples % len(lengths)
    length_budgets = {
        length: base_per_length + (1 if idx < extra_length_budget else 0)
        for idx, length in enumerate(lengths)
    }

    for length in lengths:
        bucket = sorted(length_to_metapaths[length])
        if not bucket:
            continue
        budget = min(length_budgets[length], sum(available_counts[mp] for mp in bucket))
        if budget <= 0:
            continue

        base_per_metapath = budget // len(bucket)
        extra_bucket_budget = budget % len(bucket)

        for idx, metapath in enumerate(bucket):
            allocation = base_per_metapath + (1 if idx < extra_bucket_budget else 0)
            allocations[metapath] = min(allocation, available_counts[metapath])

    return {metapath: n for metapath, n in allocations.items() if n > 0}


def _build_index_maps(hetmat: HetMat) -> tuple[dict, dict]:
    genes = hetmat.get_nodes("Gene")
    bps = hetmat.get_nodes("Biological Process")
    entrez_to_idx = dict(zip(genes["identifier"], genes["position"]))
    go_to_idx = dict(zip(bps["identifier"], bps["position"]))
    return entrez_to_idx, go_to_idx


def sample_metapath_concordance(
    hetmat: HetMat,
    api_df: pd.DataFrame,
    *,
    data_dir: Path | str,
    metapath_api: str,
    n_samples: int = 1000,
    seed: int = 42,
) -> tuple[pd.DataFrame, dict]:
    """Sample API rows for one metapath and compare direct vs API DWPC values."""
    api_subset = api_df[api_df["metapath_abbreviation"].astype(str) == str(metapath_api)].copy()
    if api_subset.empty:
        return pd.DataFrame(), {
            "metapath_api": str(metapath_api),
            "status": "skipped",
            "reason": "no api rows",
        }

    neo4j_to_entrez, neo4j_to_go = load_neo4j_mappings(Path(data_dir))
    api_subset["go_id"] = api_subset["neo4j_source_id"].map(neo4j_to_go)
    api_subset["entrez_gene_id"] = api_subset["neo4j_target_id"].map(neo4j_to_entrez)
    api_subset = api_subset.dropna(subset=["go_id", "entrez_gene_id"]).copy()
    if api_subset.empty:
        return pd.DataFrame(), {
            "metapath_api": str(metapath_api),
            "status": "skipped",
            "reason": "no mapped rows",
        }

    entrez_to_idx, go_to_idx = _build_index_maps(hetmat)
    api_subset["source_idx"] = api_subset["entrez_gene_id"].map(entrez_to_idx)
    api_subset["target_idx"] = api_subset["go_id"].map(go_to_idx)
    api_subset = api_subset.dropna(subset=["source_idx", "target_idx"]).copy()
    if api_subset.empty:
        return pd.DataFrame(), {
            "metapath_api": str(metapath_api),
            "status": "skipped",
            "reason": "no valid indexed rows",
        }

    api_subset["source_idx"] = api_subset["source_idx"].astype(int)
    api_subset["target_idx"] = api_subset["target_idx"].astype(int)
    metapath_direct = reverse_metapath_abbrev(str(metapath_api))
    api_subset["direct_dwpc"] = hetmat.get_dwpc_for_pairs(
        metapath_direct,
        api_subset["source_idx"].to_numpy(),
        api_subset["target_idx"].to_numpy(),
        transform=True,
    )
    api_subset["api_dwpc"] = api_subset["dwpc"].astype(float)
    api_subset["diff"] = api_subset["direct_dwpc"] - api_subset["api_dwpc"]
    api_subset["abs_diff"] = api_subset["diff"].abs()
    api_subset["metapath_api"] = str(metapath_api)

    requested_pairs = int(n_samples)
    sample_size = min(requested_pairs, len(api_subset))
    rng = np.random.RandomState(seed)
    shuffled = api_subset.iloc[rng.permutation(len(api_subset))].reset_index(drop=True)

    sample_df = shuffled.iloc[:sample_size].copy()
    while (
        sample_size < len(shuffled)
        and (
            sample_df["api_dwpc"].nunique(dropna=True) < 2
            or sample_df["direct_dwpc"].nunique(dropna=True) < 2
        )
    ):
        sample_size = min(len(shuffled), max(sample_size + requested_pairs, sample_size * 2))
        sample_df = shuffled.iloc[:sample_size].copy()

    correlation_defined = (
        sample_df["api_dwpc"].nunique(dropna=True) >= 2
        and sample_df["direct_dwpc"].nunique(dropna=True) >= 2
    )
    correlation = sample_df["direct_dwpc"].corr(sample_df["api_dwpc"]) if correlation_defined else np.nan
    exact_matches = int((sample_df["abs_diff"] < 1e-5).sum())
    summary = {
        "metapath_api": str(metapath_api),
        "metapath_length": int(metapath_length(str(metapath_api))),
        "status": "validated",
        "n_pairs_requested": requested_pairs,
        "n_pairs": int(len(sample_df)),
        "correlation_defined": bool(correlation_defined),
        "exact_matches": exact_matches,
        "pct_exact": 100.0 * exact_matches / float(len(sample_df)),
        "correlation": float(correlation) if pd.notna(correlation) else np.nan,
        "mean_abs_diff": float(sample_df["abs_diff"].mean()),
        "median_abs_diff": float(sample_df["abs_diff"].median()),
        "max_abs_diff": float(sample_df["abs_diff"].max()),
    }
    return sample_df, summary
