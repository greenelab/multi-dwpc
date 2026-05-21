"""Shared intermediate sharing computation for year and LV analyses.

Functions here are used by ``scripts/pipeline/intermediate_sharing.py`` (which
dispatches on ``--analysis-type {year,lv}``).  Analysis-specific logic
(metapath selection strategy, sharing-stat definitions, DWPC numpy loading)
stays in that script.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.dwpc_direct import reverse_metapath_abbrev
from src.path_enumeration import (
    EdgeLoader,
    NodeMaps,
    enumerate_paths,
    parse_metapath,
    select_paths,
)


def load_runs_at_b(
    runs_path: Path,
    b: int,
    group_cols: list[str],
) -> pd.DataFrame:
    """Load all_runs_long.csv filtered to permutation null at the given B.

    Averages effect-size across seeds within each group defined by
    *group_cols*.  Year analysis passes
    ``["go_id", "metapath"]`` (+ optional ``"year"``); LV analysis
    passes ``["lv_id", "target_id", "target_name", "node_type", "metapath"]``.
    """
    runs_df = pd.read_csv(runs_path)
    runs_df = runs_df[(runs_df["b"] == b) & (runs_df["control"] == "permuted")].copy()
    if runs_df.empty:
        raise ValueError(f"No permuted rows found for b={b} in {runs_path}")

    if "effect_size_z" in runs_df.columns:
        z_col = "effect_size_z"
    elif "effect_size_d" in runs_df.columns:
        z_col = "effect_size_d"
    elif "d" in runs_df.columns:
        z_col = "d"
    else:
        raise ValueError(
            f"{runs_path} has neither 'effect_size_z', 'effect_size_d', nor 'd' column; "
            f"columns present: {sorted(runs_df.columns.tolist())}"
        )

    present_cols = [c for c in group_cols if c in runs_df.columns]
    return runs_df.groupby(present_cols, as_index=False).agg(
        effect_size_z=(z_col, "mean"),
        diff_perm=("diff", "mean"),
    )


def compute_dwpc_thresholds(
    dwpc_lookup: dict[tuple, float],
    percentile: float,
) -> dict[tuple[str, str], float]:
    """Compute DWPC percentile thresholds per (gene_set_id, metapath).

    Works for both year keys ``(go_id, metapath, gene_id)`` and LV keys
    ``(lv_id, metapath, gene_id)`` -- groups by the first two tuple elements.
    """
    grouped: dict[tuple[str, str], list[float]] = {}
    for key, dwpc in dwpc_lookup.items():
        group_key = (key[0], key[1])
        grouped.setdefault(group_key, []).append(dwpc)

    return {
        k: float(np.percentile(v, percentile))
        for k, v in grouped.items()
        if v
    }


def compute_dwpc_z_stats(
    dwpc_lookup: dict[tuple, float],
) -> dict[tuple[str, str], tuple[float, float]]:
    """Compute (mean, std) of DWPC per (gene_set_id, metapath) over gene universe.

    Used for gene filtering via z-score: keep gene if
    ``(gene_dwpc - mean) / std >= z_min``. Groups may have n>=2; smaller
    groups are omitted (std undefined).
    """
    grouped: dict[tuple[str, str], list[float]] = {}
    for key, dwpc in dwpc_lookup.items():
        group_key = (key[0], key[1])
        grouped.setdefault(group_key, []).append(dwpc)

    stats: dict[tuple[str, str], tuple[float, float]] = {}
    for k, v in grouped.items():
        if len(v) < 2:
            continue
        arr = np.asarray(v, dtype=float)
        std = float(arr.std(ddof=1))
        if std <= 0:
            continue
        stats[k] = (float(arr.mean()), std)
    return stats


def load_dwpc_from_numpy(
    output_dirs: list[Path] | Path,
    analysis_type: str = "lv",
    gene_set_id_override: str | None = None,
) -> dict[tuple, float]:
    """Load per-gene DWPC values from numpy arrays.

    Returns dict mapping ``(gene_set_id, metapath, gene_id) -> dwpc``.
    Works for both LV and year analyses by selecting the appropriate ID
    column from the feature manifest.

    Parameters
    ----------
    output_dirs
        One or more experiment output directories containing
        ``gene_feature_scores.npy``, ``gene_ids.npy``, and
        ``feature_manifest.csv``.
    analysis_type
        ``"lv"`` reads ``lv_id`` from manifest; ``"year"`` reads ``go_id``.
    gene_set_id_override
        If set, use this as the gene-set ID for all entries instead of
        reading from the manifest.  Useful when processing a single LV.
    """
    if isinstance(output_dirs, Path):
        output_dirs = [output_dirs]

    id_col = "lv_id" if analysis_type == "lv" else "go_id"
    dwpc_lookup: dict[tuple, float] = {}

    for out_dir in output_dirs:
        scores_path = out_dir / "gene_feature_scores.npy"
        genes_path = out_dir / "gene_ids.npy"
        manifest_path = out_dir / "feature_manifest.csv"

        if not all(p.exists() for p in [scores_path, genes_path, manifest_path]):
            continue

        scores = np.load(scores_path)
        gene_ids = np.load(genes_path)
        manifest = pd.read_csv(manifest_path)

        fallback_id = gene_set_id_override or (
            out_dir.name if id_col not in manifest.columns else None
        )

        for _, row in manifest.iterrows():
            feature_idx = int(row["feature_idx"])
            metapath = row["metapath"]
            gs_id = gene_set_id_override or row.get(id_col, fallback_id)

            feature_scores = scores[:, feature_idx]
            nonzero = feature_scores > 0
            for gene_idx in np.flatnonzero(nonzero):
                dwpc_lookup[(gs_id, metapath, int(gene_ids[gene_idx]))] = float(
                    feature_scores[gene_idx]
                )

    return dwpc_lookup


def enumerate_gene_intermediates(
    genes: list[int],
    target_pos: int,
    metapath: str,
    edge_loader: EdgeLoader,
    maps: NodeMaps,
    *,
    path_top_k: int = 100,
    degree_d: float = 0.5,
    gene_set_id: str | None = None,
    dwpc_lookup: dict[tuple, float] | None = None,
    dwpc_thresholds: dict[tuple[str, str], float] | None = None,
    dwpc_z_stats: dict[tuple[str, str], tuple[float, float]] | None = None,
    dwpc_z_min: float = 1.65,
    path_z_min: float | None = None,
    path_enumeration_cap: int | None = None,
    record_paths: list[dict] | None = None,
    record_extra: dict | None = None,
    debug: bool = False,
) -> dict[int, set[str]]:
    """Enumerate paths for *genes* and return ``{gene_id: set_of_intermediate_ids}``.

    The metapath is in Gene-first orientation.  Internally it is reversed for
    enumeration (``enumerate_paths`` walks source -> target, and source here is
    the non-Gene end).

    Parameters
    ----------
    gene_set_id
        Identifier used as the first key element when looking up DWPC
        thresholds (``go_id`` for year, ``lv_id`` for LV).
    record_extra
        Extra key-value pairs merged into every ``record_paths`` entry.
        Year analysis passes ``{"year": "2016"}``; LV passes
        ``{"lv_id": "LV246"}``.
    """
    reversed_mp = reverse_metapath_abbrev(metapath)
    nodes, edges = parse_metapath(reversed_mp)

    if debug:
        print(f"      Metapath {metapath} -> reversed {reversed_mp}: nodes={nodes}, edges={edges}")

    gene_intermediates: dict[int, set[str]] = {}
    genes_found = 0
    genes_with_paths = 0
    genes_filtered_by_dwpc = 0

    gene_id_map = maps.id_to_pos.get("G", {})

    # Gene filter: z-based (mean+std of universe DWPC) overrides percentile path.
    dwpc_threshold = 0.0
    dwpc_mean: float | None = None
    dwpc_std: float | None = None
    if dwpc_z_stats is not None and gene_set_id is not None:
        stats = dwpc_z_stats.get((gene_set_id, metapath))
        if stats is not None:
            dwpc_mean, dwpc_std = stats
            if debug:
                print(
                    f"      DWPC z-filter for {gene_set_id}/{metapath}: "
                    f"mean={dwpc_mean:.4f}, std={dwpc_std:.4f}, z_min={dwpc_z_min}"
                )
    elif dwpc_thresholds is not None and gene_set_id is not None:
        dwpc_threshold = dwpc_thresholds.get((gene_set_id, metapath), 0.0)
        if debug and dwpc_threshold > 0:
            print(f"      DWPC threshold for {gene_set_id}/{metapath}: {dwpc_threshold:.4f}")

    # First pass: enumerate paths per surviving gene. Scores pooled for optional
    # path-z threshold computed across all genes for this (gene_set, metapath).
    enumeration_cap = path_enumeration_cap if path_enumeration_cap is not None else path_top_k
    per_gene_paths: list[tuple[int, int, list[tuple[float, list[int]]]]] = []
    for gene_id in genes:
        if dwpc_lookup is not None and gene_set_id is not None:
            gene_dwpc = dwpc_lookup.get((gene_set_id, metapath, gene_id), 0.0)
            if dwpc_mean is not None and dwpc_std is not None and dwpc_std > 0:
                z = (gene_dwpc - dwpc_mean) / dwpc_std
                if z < dwpc_z_min:
                    genes_filtered_by_dwpc += 1
                    continue
            elif dwpc_threshold > 0 and gene_dwpc < dwpc_threshold:
                genes_filtered_by_dwpc += 1
                continue

        gene_pos = (
            gene_id_map.get(gene_id)
            or gene_id_map.get(str(gene_id))
            or gene_id_map.get(int(gene_id) if isinstance(gene_id, str) else gene_id)
        )
        if gene_pos is None:
            continue
        genes_found += 1

        try:
            paths = enumerate_paths(
                target_pos, gene_pos, nodes, edges, edge_loader,
                top_k=enumeration_cap, degree_d=degree_d,
            )
        except Exception as exc:
            if debug:
                print(f"      Exception for gene {gene_id}: {exc}")
            continue

        if not paths:
            continue
        genes_with_paths += 1
        per_gene_paths.append((gene_id, gene_pos, paths))

    # Optional path-z threshold: pool scores across all genes for this metapath.
    path_score_cutoff: float | None = None
    if path_z_min is not None and per_gene_paths:
        all_scores = np.asarray(
            [score for _, _, paths in per_gene_paths for score, _ in paths],
            dtype=float,
        )
        if all_scores.size >= 2:
            pool_std = float(all_scores.std(ddof=1))
            if pool_std > 0:
                pool_mean = float(all_scores.mean())
                path_score_cutoff = pool_mean + path_z_min * pool_std
                if debug:
                    print(
                        f"      Path-z filter for {gene_set_id}/{metapath}: "
                        f"pool_mean={pool_mean:.4g}, pool_std={pool_std:.4g}, "
                        f"z_min={path_z_min}, cutoff={path_score_cutoff:.4g}, "
                        f"n_paths={all_scores.size}"
                    )

    for gene_id, gene_pos, paths in per_gene_paths:
        if path_score_cutoff is not None:
            selected = [p for p in paths if p[0] >= path_score_cutoff]
        else:
            selected = select_paths(
                paths,
                selection_method="effective_number",
                top_paths=path_top_k,
                path_cumulative_frac=None,
                path_min_count=1,
                path_max_count=None,
            )
        if not selected:
            continue

        intermediates: set[str] = set()
        for path_rank, (score, pos_path) in enumerate(selected, start=1):
            path_node_ids: list[str | None] = []
            for i, (node_type, pos) in enumerate(zip(nodes, pos_path)):
                node_id = maps.pos_to_id.get(node_type, {}).get(int(pos))
                qualified = f"{node_type}:{node_id}" if node_id is not None else None
                path_node_ids.append(qualified)
                if 0 < i < len(nodes) - 1 and qualified is not None:
                    intermediates.add(qualified)

            if record_paths is not None:
                reversed_ids = list(reversed(path_node_ids))
                reversed_types = list(reversed(nodes))
                record: dict = {
                    "metapath": metapath,
                    "gene_id": gene_id,
                    "path_rank": int(path_rank),
                    "path_score": float(score),
                }
                if record_extra:
                    record.update(record_extra)
                for hop_idx, (nt, qid) in enumerate(zip(reversed_types, reversed_ids)):
                    record[f"hop_{hop_idx}_type"] = nt
                    record[f"hop_{hop_idx}_id"] = qid
                record_paths.append(record)

        if intermediates:
            gene_intermediates[gene_id] = intermediates

    if debug:
        print(
            f"      genes_found={genes_found}, genes_with_paths={genes_with_paths}, "
            f"intermediates={len(gene_intermediates)}, filtered_by_dwpc={genes_filtered_by_dwpc}"
        )

    return gene_intermediates, genes_filtered_by_dwpc


def compute_intermediate_coverage(
    gene_intermediates: dict[int, set[str]],
    node_name_maps: dict[str, dict[str, str]] | None = None,
) -> tuple[dict, list[dict]]:
    """Compute intermediate coverage statistics (superset of year + LV metrics)."""
    n_genes = len(gene_intermediates)
    if n_genes == 0:
        return {
            "n_shared_intermediates_2plus": 0,
            "n_shared_intermediates_quarter": 0,
            "n_shared_intermediates_majority": 0,
            "n_shared_intermediates_all": 0,
            "pct_intermediates_shared_2plus": 0.0,
            "pct_intermediates_shared_quarter": 0.0,
            "pct_intermediates_shared_majority": 0.0,
            "pct_intermediates_shared_all": 0.0,
            "n_intermediates_cover_50pct": None,
            "n_intermediates_cover_80pct": None,
            "top1_intermediate_coverage": 0.0,
            "top5_intermediate_coverage": 0.0,
        }, []

    intermediate_gene_counts: dict[str, set[int]] = {}
    for gene_id, ints in gene_intermediates.items():
        for int_id in ints:
            intermediate_gene_counts.setdefault(int_id, set()).add(gene_id)

    intermediate_stats: list[dict] = []
    for int_id, genes_using in intermediate_gene_counts.items():
        intermediate_name = None
        if node_name_maps and ":" in int_id:
            node_type, node_identifier = int_id.split(":", 1)
            if node_type in node_name_maps:
                intermediate_name = node_name_maps[node_type].get(node_identifier)

        intermediate_stats.append({
            "intermediate_id": int_id,
            "intermediate_name": intermediate_name,
            "n_genes_using": len(genes_using),
            "pct_genes_using": len(genes_using) / n_genes * 100,
            "genes_using": sorted(genes_using),
        })

    intermediate_stats.sort(key=lambda x: x["n_genes_using"], reverse=True)

    n_total = len(intermediate_gene_counts)
    n_shared_2plus = sum(1 for g in intermediate_gene_counts.values() if len(g) >= 2)
    n_shared_quarter = sum(1 for g in intermediate_gene_counts.values() if len(g) > n_genes / 4)
    n_shared_majority = sum(1 for g in intermediate_gene_counts.values() if len(g) > n_genes / 2)
    n_shared_all = sum(1 for g in intermediate_gene_counts.values() if len(g) == n_genes)

    genes_covered: set[int] = set()
    n_for_50pct = None
    n_for_80pct = None
    for i, stat in enumerate(intermediate_stats, 1):
        genes_covered.update(stat["genes_using"])
        pct = len(genes_covered) / n_genes * 100
        if n_for_50pct is None and pct >= 50:
            n_for_50pct = i
        if n_for_80pct is None and pct >= 80:
            n_for_80pct = i
        if pct >= 100:
            break

    top1_coverage = intermediate_stats[0]["pct_genes_using"] if intermediate_stats else None
    if len(intermediate_stats) >= 5:
        top5_genes: set[int] = set()
        for stat in intermediate_stats[:5]:
            top5_genes.update(stat["genes_using"])
        top5_coverage = len(top5_genes) / n_genes * 100
    else:
        top5_coverage = None

    coverage_stats = {
        "n_shared_intermediates_2plus": n_shared_2plus,
        "n_shared_intermediates_quarter": n_shared_quarter,
        "n_shared_intermediates_majority": n_shared_majority,
        "n_shared_intermediates_all": n_shared_all,
        "pct_intermediates_shared_2plus": n_shared_2plus / n_total * 100 if n_total > 0 else None,
        "pct_intermediates_shared_quarter": n_shared_quarter / n_total * 100 if n_total > 0 else None,
        "pct_intermediates_shared_majority": n_shared_majority / n_total * 100 if n_total > 0 else None,
        "pct_intermediates_shared_all": n_shared_all / n_total * 100 if n_total > 0 else None,
        "n_intermediates_cover_50pct": n_for_50pct,
        "n_intermediates_cover_80pct": n_for_80pct,
        "top1_intermediate_coverage": top1_coverage,
        "top5_intermediate_coverage": top5_coverage,
    }

    return coverage_stats, intermediate_stats
