"""Query logic behind the API; mirrors the Streamlit app in ``app.py``."""

from __future__ import annotations

from functools import lru_cache
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
from django.conf import settings

from src.multi_dwpc_query import (
    DEFAULT_PATH_Z_MIN,
    query_intermediates_and_paths,
    query_metapath_z,
)
from src.null_bundle import NullBundle
from src.path_enumeration import load_node_names
from src.query_dwpc import QueryDwpc

POOL_Z_MIN = 1.65  # matches HPC end-to-end --effect-size-threshold default
DRILLDOWN_PATH_TOP_K = 500
OVERALL_TOP_PATHS = 40  # matches plot_lv_top_paths_subgraph top_n_paths default
OVERALL_TOP_SHARED_PER_HOP = 15  # matches plot_lv_top_shared_subgraph top_n_per_hop default
HEATMAP_TOP_INTERMEDIATES = 30
RANDOM_GENE_COUNT = 20


@lru_cache(maxsize=1)
def get_query_resources() -> tuple[NullBundle, QueryDwpc]:
    """Precomputed stratum summaries plus on-the-fly DWPC from edge files."""
    data_dir = settings.MULTI_DWPC_DATA_DIR
    return NullBundle(data_dir / "null_bundle", data_dir), QueryDwpc(data_dir)


@lru_cache(maxsize=1)
def load_bp_nodes() -> pd.DataFrame:
    df = pd.read_csv(settings.MULTI_DWPC_DATA_DIR / "nodes" / "Biological Process.tsv", sep="\t")
    return df.sort_values("name").reset_index(drop=True)


@lru_cache(maxsize=1)
def load_gene_nodes() -> pd.DataFrame:
    return pd.read_csv(settings.MULTI_DWPC_DATA_DIR / "nodes" / "Gene.tsv", sep="\t")


@lru_cache(maxsize=1)
def load_name_maps() -> dict[str, dict[str, str]]:
    return load_node_names(settings.MULTI_DWPC_REPO_ROOT)


@lru_cache(maxsize=1)
def symbol_to_entrez_map() -> dict[str, int]:
    genes = load_gene_nodes()
    return dict(zip(genes["name"].astype(str), genes["identifier"].astype(int)))


@lru_cache(maxsize=1)
def valid_entrez_ids() -> frozenset[int]:
    return frozenset(load_gene_nodes()["identifier"].astype(int).tolist())


def parse_gene_input(raw: str) -> list[int]:
    symbol_to_entrez = symbol_to_entrez_map()
    entrez_set = valid_entrez_ids()
    out: list[int] = []
    for tok in raw.replace(",", "\n").split():
        if tok.lstrip("-").isdigit() and int(tok) in entrez_set:
            out.append(int(tok))
        elif tok in symbol_to_entrez:
            out.append(symbol_to_entrez[tok])
    gene_ids = list(dict.fromkeys(out))  # a gene set: repeated symbols or IDs count once
    if not gene_ids:
        raise ValueError("No valid gene symbols or Entrez IDs parsed from input.")
    return gene_ids


def _df_records(df: pd.DataFrame) -> list[dict[str, Any]]:
    if df.empty:
        return []
    # object dtype so null survives in float columns; NaN is not valid JSON
    clean = df.astype(object).where(pd.notnull(df), None)
    return clean.to_dict(orient="records")


def _qualified_to_name(qid: Any, name_maps: dict[str, dict[str, str]]) -> Any:
    if not isinstance(qid, str) or ":" not in qid:
        return qid
    node_type, _, node_id = qid.partition(":")
    return name_maps.get(node_type, {}).get(node_id, qid)


def _hop_id_cols(paths_df: pd.DataFrame) -> list[str]:
    return sorted(c for c in paths_df.columns if c.startswith("hop_") and c.endswith("_id"))


def _paths_with_names(paths_df: pd.DataFrame) -> pd.DataFrame:
    if paths_df.empty:
        return paths_df
    name_maps = load_name_maps()
    gene_names = name_maps.get("G", {})
    out = paths_df.copy()
    out.insert(
        out.columns.get_loc("gene_id") + 1,
        "gene",
        out["gene_id"].astype(str).map(lambda g: gene_names.get(g, g)),
    )
    for col in _hop_id_cols(paths_df):
        out[f"{col[: -len('_id')]}_name"] = out[col].map(
            lambda qid: _qualified_to_name(qid, name_maps)
        )
    return out


def _sharing_with_names(sharing_df: pd.DataFrame) -> pd.DataFrame:
    if sharing_df.empty:
        return sharing_df
    name_maps = load_name_maps()
    gene_names = name_maps.get("G", {})
    out = sharing_df.copy()
    out.insert(
        1,
        "intermediate",
        out["intermediate_id"].map(lambda qid: _qualified_to_name(qid, name_maps)),
    )
    out["genes"] = out["gene_ids"].map(
        lambda ids: ", ".join(gene_names.get(g, g) for g in str(ids).split(","))
    )
    return out


def _per_row_intermediates(paths_df: pd.DataFrame) -> list[list[str]]:
    """For each row, the intermediate-hop IDs (excludes source and the row's own target)."""
    hop_id_cols = _hop_id_cols(paths_df)
    out: list[list[str]] = []
    for _, row in paths_df.iterrows():
        vals = [row[c] for c in hop_id_cols]
        last = next((i for i in range(len(vals) - 1, -1, -1) if isinstance(vals[i], str)), -1)
        out.append([vals[i] for i in range(1, last) if isinstance(vals[i], str)])
    return out


def filter_top_shared_pooled(
    paths_df: pd.DataFrame, top_n_per_hop: int = OVERALL_TOP_SHARED_PER_HOP
) -> pd.DataFrame:
    """Per-hop top-N intermediates by gene coverage. Rows from metapaths with zero
    intermediate hops pass trivially (no hops to filter)."""
    if paths_df.empty:
        return paths_df
    hop_id_cols = _hop_id_cols(paths_df)

    per_row: list[tuple[list[int], list[str]]] = []
    for _, row in paths_df.iterrows():
        vals = [row[c] if c in row.index else None for c in hop_id_cols]
        last = next((i for i in range(len(vals) - 1, -1, -1) if isinstance(vals[i], str)), -1)
        if last < 1:
            per_row.append(([], []))
            continue
        hops = [i for i in range(1, last) if isinstance(vals[i], str)]
        per_row.append((hops, [vals[i] for i in hops]))

    per_hop_genes: dict[int, dict[str, set[int]]] = {}
    for (hops, iids), gene_id in zip(per_row, paths_df["gene_id"]):
        for h, iid in zip(hops, iids):
            per_hop_genes.setdefault(h, {}).setdefault(iid, set()).add(int(gene_id))

    allowed: dict[int, set[str]] = {}
    for h, m in per_hop_genes.items():
        ranked = sorted(m.items(), key=lambda kv: len(kv[1]), reverse=True)
        allowed[h] = {iid for iid, _ in ranked[:top_n_per_hop]}

    keep = [
        all(iid in allowed.get(h, set()) for h, iid in zip(hops, iids))
        for hops, iids in per_row
    ]
    return paths_df[keep].reset_index(drop=True)


def filter_top_paths_stratified(
    paths_df: pd.DataFrame, top_n: int = OVERALL_TOP_PATHS
) -> pd.DataFrame:
    if paths_df.empty or "path_score" not in paths_df.columns:
        return paths_df
    if "metapath" in paths_df.columns and paths_df["metapath"].nunique() > 1:
        per = max(1, top_n // max(paths_df["metapath"].nunique(), 1))
        return (
            paths_df.sort_values("path_score", ascending=False)
            .groupby("metapath", group_keys=False)
            .head(per)
            .reset_index(drop=True)
        )
    return paths_df.nlargest(top_n, "path_score").reset_index(drop=True)


def filter_top_paths(paths_df: pd.DataFrame, top_n: int = OVERALL_TOP_PATHS) -> pd.DataFrame:
    if paths_df.empty or "path_score" not in paths_df.columns:
        return paths_df
    return paths_df.nlargest(top_n, "path_score").reset_index(drop=True)


def build_heatmap(paths_df: pd.DataFrame) -> dict[str, Any] | None:
    """Binary matrix: intermediates (rows) x genes (cols), 1 if gene uses intermediate."""
    if paths_df.empty:
        return None
    per_row = _per_row_intermediates(paths_df)
    flat = [(g, iid) for g, ints in zip(paths_df["gene_id"], per_row) for iid in ints]
    if not flat:
        return None
    long = pd.DataFrame(flat, columns=["gene_id", "intermediate_id"])
    coverage = long.groupby("intermediate_id")["gene_id"].nunique().sort_values(ascending=False)
    top_ids = coverage.head(HEATMAP_TOP_INTERMEDIATES).index.tolist()
    long = long[long["intermediate_id"].isin(top_ids)]
    pivot = (
        long.assign(value=1)
        .pivot_table(
            index="intermediate_id",
            columns="gene_id",
            values="value",
            aggfunc="max",
            fill_value=0,
        )
        .reindex(top_ids)
    )
    name_maps = load_name_maps()
    gene_names = name_maps.get("G", {})
    return {
        "intermediates": [
            {"id": qid, "name": _qualified_to_name(qid, name_maps)} for qid in pivot.index
        ],
        "genes": [
            {"id": int(g), "name": gene_names.get(str(g), str(g))} for g in pivot.columns
        ],
        "cells": pivot.to_numpy().astype(int).tolist(),
    }


def build_subgraph(
    paths_df: pd.DataFrame,
    title: str,
    total_genes: int,
    weight_by: str,
) -> dict[str, Any] | None:
    """Layered graph of surviving subpaths. Columns = hop positions; nodes stacked within column.

    Handles mixed-length paths: each row contributes only its non-null hops. Source
    column is left-aligned; target column is right-aligned at max_hops-1.

    ``weight_by``: ``"coverage"`` (share of genes through the node) or ``"score"``
    (summed path_score through the node, relative to the largest intermediate).
    ``nodes`` are in display order within each column. ``scale`` is the relative
    node size, floored at 0.2 so tiny values remain visible. A node's ``id`` includes
    its column because the same node can appear in more than one.
    """
    if paths_df.empty:
        return None
    name_maps = load_name_maps()
    hop_id_cols = _hop_id_cols(paths_df)
    hop_type_cols = [c[: -len("_id")] + "_type" for c in hop_id_cols]
    n_hops = len(hop_id_cols)

    column_nodes: list[dict[str, str]] = [{} for _ in range(n_hops)]
    column_genes: list[dict[str, set[int]]] = [{} for _ in range(n_hops)]
    column_scores: list[dict[str, float]] = [{} for _ in range(n_hops)]
    edges: set[tuple[int, str, int, str]] = set()
    for _, row in paths_df.iterrows():
        ids = [row[c] if c in row else None for c in hop_id_cols]
        types = [row[c] if c in row else None for c in hop_type_cols]
        gene_id = row.get("gene_id") if "gene_id" in row else None
        score = float(row.get("path_score", 0.0) or 0.0)
        last_col = next(
            (i for i in range(n_hops - 1, -1, -1) if isinstance(ids[i], str)), -1,
        )
        if last_col < 1:
            continue
        col_map = {i: i for i in range(last_col)}
        col_map[last_col] = n_hops - 1
        prev = None
        for local_i in range(last_col + 1):
            qid = ids[local_i]
            if not isinstance(qid, str):
                prev = None
                continue
            dst_col = col_map[local_i]
            column_nodes[dst_col][qid] = types[local_i] if isinstance(types[local_i], str) else "?"
            if gene_id is not None and not pd.isna(gene_id):
                column_genes[dst_col].setdefault(qid, set()).add(int(gene_id))
            column_scores[dst_col][qid] = column_scores[dst_col].get(qid, 0.0) + score
            if prev is not None:
                edges.add((prev[0], prev[1], dst_col, qid))
            prev = (dst_col, qid)

    # Global denominator for score weights: max summed score across intermediate nodes.
    max_score_any_node = 0.0
    for ci in range(1, n_hops - 1):
        for qid in column_nodes[ci]:
            max_score_any_node = max(max_score_any_node, column_scores[ci].get(qid, 0.0))

    nodes: list[dict[str, Any]] = []
    for ci, col in enumerate(column_nodes):
        def coverage(q: str, ci: int = ci) -> int:
            return len(column_genes[ci].get(q, set()))

        def score_weight(q: str, ci: int = ci) -> float:
            return column_scores[ci].get(q, 0.0)

        sort_key = score_weight if weight_by == "score" else coverage
        ordered = sorted(col, key=lambda q: (-sort_key(q), _qualified_to_name(q, name_maps)))
        for qid in ordered:
            name = _qualified_to_name(qid, name_maps)
            if 0 < ci < n_hops - 1:
                if weight_by == "score" and max_score_any_node > 0:
                    frac = score_weight(qid) / max_score_any_node
                    label = f"{name} ({100 * frac:.0f}% score)"
                elif weight_by == "coverage" and total_genes:
                    frac = coverage(qid) / total_genes
                    label = f"{name} ({100 * frac:.0f}% genes)"
                else:
                    frac = 1.0
                    label = name
                scale = max(0.2, min(1.0, frac))
            else:
                scale = 1.0
                label = name
            nodes.append(
                {
                    "id": f"{ci}|{qid}",
                    "node_id": qid,
                    "name": name,
                    "type": col[qid],
                    "column": ci,
                    "label": label,
                    "scale": scale,
                }
            )

    return {
        "title": title,
        "columns": [sorted(set(col.values()))[0] if col else "" for col in column_nodes],
        "nodes": nodes,
        "edges": [
            {"source": f"{src_col}|{src}", "target": f"{dst_col}|{dst}"}
            for src_col, src, dst_col, dst in sorted(edges)
        ],
    }


def _drilldown(
    gene_ids: list[int], target_id: str, metapath: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    return query_intermediates_and_paths(
        gene_ids=gene_ids,
        target_id=target_id,
        metapath=metapath,
        repo_root=settings.MULTI_DWPC_REPO_ROOT,
        path_top_k=DRILLDOWN_PATH_TOP_K,
        path_z_min=DEFAULT_PATH_Z_MIN,
    )


def _metapath_z(gene_ids: list[int], target_id: str) -> pd.DataFrame:
    bundle, query_dwpc = get_query_resources()
    return query_metapath_z(gene_ids, target_id, bundle=bundle, query_dwpc=query_dwpc)


def run_metapath_ranking(payload: dict[str, Any]) -> dict[str, Any]:
    query_start = perf_counter()
    target_id = payload["target_id"]
    gene_ids = parse_gene_input(payload["genes_raw"])
    z_df = _metapath_z(gene_ids, target_id)
    return {
        "target_id": target_id,
        "gene_ids": gene_ids,
        "metapath_ranking": _df_records(z_df),
        "elapsed_s": perf_counter() - query_start,
    }


def run_metapath_drilldown(payload: dict[str, Any]) -> dict[str, Any]:
    target_id = payload["target_id"]
    gene_ids = parse_gene_input(payload["genes_raw"])
    metapath = payload["metapath"]
    if metapath not in get_query_resources()[0].metapaths:
        raise ValueError(f"Unknown metapath: {metapath}")

    sharing_df, paths_df, diagnostics = _drilldown(gene_ids, target_id, metapath)
    top_df = filter_top_paths(paths_df)
    return {
        "target_id": target_id,
        "gene_ids": gene_ids,
        "metapath": metapath,
        "path_z_min": DEFAULT_PATH_Z_MIN,
        "intermediate_sharing": _df_records(_sharing_with_names(sharing_df)),
        "heatmap": build_heatmap(paths_df),
        "paths": _df_records(_paths_with_names(paths_df)),
        "subgraph": build_subgraph(
            top_df,
            f"{metapath}: top {len(top_df)} paths by score",
            len(gene_ids),
            "score",
        ),
        "diagnostics": diagnostics,
    }


def _pool_across_metapaths(
    gene_ids: list[int], target_id: str, z_df: pd.DataFrame
) -> tuple[pd.DataFrame, list[str], float]:
    """Pool at POOL_Z_MIN; fall back to z > 0 only if nothing clears."""
    z_used = POOL_Z_MIN
    pool_mps = z_df.loc[z_df["effect_size_z"] >= POOL_Z_MIN, "metapath"].tolist()
    if not pool_mps:
        z_used = 0.0
        pool_mps = z_df.loc[z_df["effect_size_z"] > 0, "metapath"].tolist()

    frames: list[pd.DataFrame] = []
    for mp in pool_mps:
        _, mp_paths, _ = _drilldown(gene_ids, target_id, mp)
        if not mp_paths.empty:
            frames.append(mp_paths)
    if not frames:
        return pd.DataFrame(), pool_mps, z_used
    return pd.concat(frames, ignore_index=True), pool_mps, z_used


def run_overall_subgraphs(payload: dict[str, Any]) -> dict[str, Any]:
    target_id = payload["target_id"]
    gene_ids = parse_gene_input(payload["genes_raw"])
    z_df = _metapath_z(gene_ids, target_id)
    pooled_df, pooled_mps, z_used = _pool_across_metapaths(gene_ids, target_id, z_df)

    top_df = filter_top_paths_stratified(pooled_df)
    return {
        "target_id": target_id,
        "gene_ids": gene_ids,
        "pooled_metapaths": pooled_mps,
        "pool_z_min": POOL_Z_MIN,
        "pool_z_used": z_used,
        "top_shared_subgraph": build_subgraph(
            filter_top_shared_pooled(pooled_df),
            f"Top shared intermediates; {len(pooled_mps)} metapaths",
            len(gene_ids),
            "coverage",
        ),
        "top_paths_subgraph": build_subgraph(
            top_df,
            f"Top {len(top_df)} paths by score; {len(pooled_mps)} metapaths",
            len(gene_ids),
            "score",
        ),
    }


def list_biological_processes(q: str = "", limit: int = 100) -> list[dict[str, str]]:
    df = load_bp_nodes()
    if q:
        q_lower = q.lower()
        names = df["name"].astype(str).str.lower()
        mask = (
            names.str.contains(q_lower, regex=False)
            | df["identifier"].astype(str).str.lower().str.contains(q_lower, regex=False)
        )
        # Best matches first, so an exact name is never cut off by the limit.
        rank = names.ne(q_lower).astype(int) + (~names.str.startswith(q_lower)).astype(int)
        df = df.assign(_rank=rank).loc[mask].sort_values("_rank", kind="stable")
    return _df_records(df[["identifier", "name"]].head(limit))


def random_biological_process() -> dict[str, str]:
    row = load_bp_nodes().sample(n=1).iloc[0]
    return {"identifier": str(row["identifier"]), "name": str(row["name"])}


def random_genes(size: int = RANDOM_GENE_COUNT) -> list[str]:
    rng = np.random.default_rng()
    symbols = rng.choice(load_gene_nodes()["name"].astype(str).values, size=size, replace=False)
    return sorted(symbols.tolist())
