from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
import pandas as pd
from django.conf import settings

from src.dwpc_direct import HetMat
from src.multi_dwpc_query import (
    discover_source_target_metapaths,
    query_intermediates_and_paths,
    query_metapath_z,
)


def _df_records(df: pd.DataFrame) -> list[dict[str, Any]]:
    if df.empty:
        return []
    clean = df.where(pd.notnull(df), None)
    return clean.to_dict(orient="records")


@lru_cache(maxsize=1)
def get_hetmat() -> HetMat:
    hetmat = HetMat(data_dir=settings.MULTI_DWPC_DATA_DIR)
    if settings.MULTI_DWPC_PRELOAD_MATRICES:
        for mp in discover_source_target_metapaths(hetmat, "G", "BP"):
            hetmat.compute_dwpc_matrix(mp)
    return hetmat


@lru_cache(maxsize=1)
def load_bp_nodes() -> pd.DataFrame:
    df = pd.read_csv(settings.MULTI_DWPC_DATA_DIR / "nodes" / "Biological Process.tsv", sep="\t")
    return df.sort_values("name").reset_index(drop=True)


@lru_cache(maxsize=1)
def load_gene_nodes() -> pd.DataFrame:
    return pd.read_csv(settings.MULTI_DWPC_DATA_DIR / "nodes" / "Gene.tsv", sep="\t")


@lru_cache(maxsize=1)
def bp_name_to_id_map() -> dict[str, str]:
    df = load_bp_nodes()
    return dict(zip(df["name"].astype(str), df["identifier"].astype(str)))


@lru_cache(maxsize=1)
def symbol_to_entrez_map() -> dict[str, int]:
    genes = load_gene_nodes()
    return dict(zip(genes["name"].astype(str), genes["identifier"].astype(int)))


@lru_cache(maxsize=1)
def valid_entrez_ids() -> set[int]:
    genes = load_gene_nodes()
    return set(genes["identifier"].astype(int).tolist())


def resolve_target_id(payload: dict[str, Any]) -> str:
    target_id = str(payload.get("target_id") or "").strip()
    if target_id:
        return target_id
    target_name = str(payload.get("target_name") or "").strip()
    if not target_name:
        raise ValueError("target_id or target_name is required")
    mapped = bp_name_to_id_map().get(target_name)
    if not mapped:
        raise ValueError(f"Unknown target_name: {target_name}")
    return mapped


def parse_gene_payload(payload: dict[str, Any]) -> list[int]:
    symbol_to_entrez = symbol_to_entrez_map()
    entrez_set = valid_entrez_ids()

    out: list[int] = []

    for gid in payload.get("gene_ids") or []:
        val = int(gid)
        if val in entrez_set:
            out.append(val)

    for symbol in payload.get("gene_symbols") or []:
        sym = str(symbol).strip()
        if sym in symbol_to_entrez:
            out.append(symbol_to_entrez[sym])

    raw = str(payload.get("genes_raw") or "")
    if raw:
        for tok in raw.replace(",", "\n").split():
            t = tok.strip()
            if not t:
                continue
            if t.lstrip("-").isdigit():
                val = int(t)
                if val in entrez_set:
                    out.append(val)
            elif t in symbol_to_entrez:
                out.append(symbol_to_entrez[t])

    seen: set[int] = set()
    deduped: list[int] = []
    for gid in out:
        if gid not in seen:
            deduped.append(gid)
            seen.add(gid)

    if not deduped:
        raise ValueError("No valid gene symbols or Entrez IDs parsed from input")
    return deduped


def filter_top_shared_pooled(paths_df: pd.DataFrame, top_n_per_hop: int = 15) -> pd.DataFrame:
    if paths_df.empty:
        return paths_df
    hop_id_cols = sorted(c for c in paths_df.columns if c.startswith("hop_") and c.endswith("_id"))

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


def filter_top_paths_stratified(paths_df: pd.DataFrame, top_n: int = 40) -> pd.DataFrame:
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


def run_metapath_ranking(payload: dict[str, Any]) -> dict[str, Any]:
    target_id = resolve_target_id(payload)
    gene_ids = parse_gene_payload(payload)

    z_df = query_metapath_z(
        gene_ids,
        target_id,
        b=int(payload.get("b", 20)),
        seed=int(payload.get("seed", 42)),
        hetmat=get_hetmat(),
    )

    return {
        "target_id": target_id,
        "n_genes": len(gene_ids),
        "gene_ids": gene_ids,
        "metapath_ranking": _df_records(z_df),
    }


def run_metapath_drilldown(payload: dict[str, Any]) -> dict[str, Any]:
    target_id = resolve_target_id(payload)
    gene_ids = parse_gene_payload(payload)
    metapath = str(payload["metapath"])

    sharing_df, paths_df, diagnostics = query_intermediates_and_paths(
        gene_ids=gene_ids,
        target_id=target_id,
        metapath=metapath,
        repo_root=settings.MULTI_DWPC_REPO_ROOT,
        hetmat=get_hetmat(),
        path_top_k=int(payload.get("path_top_k", 500)),
        path_z_min=float(payload.get("path_z_min", 1.65)),
        debug=bool(payload.get("debug", False)),
    )

    return {
        "target_id": target_id,
        "n_genes": len(gene_ids),
        "gene_ids": gene_ids,
        "metapath": metapath,
        "intermediate_sharing": _df_records(sharing_df),
        "paths": _df_records(paths_df),
        "diagnostics": diagnostics,
    }


def _pool_across_metapaths(
    gene_ids: list[int],
    target_id: str,
    z_df: pd.DataFrame,
    *,
    pool_z_min: float,
    path_top_k: int,
    path_z_min: float,
) -> tuple[pd.DataFrame, list[str], float]:
    z_used = float(pool_z_min)
    pool_mps = z_df.loc[z_df["effect_size_z"] >= z_used, "metapath"].tolist()
    if not pool_mps:
        z_used = 0.0
        pool_mps = z_df.loc[z_df["effect_size_z"] > 0, "metapath"].tolist()

    frames: list[pd.DataFrame] = []
    for mp in pool_mps:
        _, mp_paths, _ = query_intermediates_and_paths(
            gene_ids=gene_ids,
            target_id=target_id,
            metapath=mp,
            repo_root=settings.MULTI_DWPC_REPO_ROOT,
            hetmat=get_hetmat(),
            path_top_k=path_top_k,
            path_z_min=path_z_min,
            debug=False,
        )
        if not mp_paths.empty:
            frames.append(mp_paths)

    if not frames:
        return pd.DataFrame(), pool_mps, z_used
    return pd.concat(frames, ignore_index=True), pool_mps, z_used


def run_overall_subgraphs(payload: dict[str, Any]) -> dict[str, Any]:
    target_id = resolve_target_id(payload)
    gene_ids = parse_gene_payload(payload)

    z_df = query_metapath_z(
        gene_ids,
        target_id,
        b=int(payload.get("b", 20)),
        seed=int(payload.get("seed", 42)),
        hetmat=get_hetmat(),
    )

    pooled_df, pooled_mps, z_used = _pool_across_metapaths(
        gene_ids,
        target_id,
        z_df,
        pool_z_min=float(payload.get("pool_z_min", 1.65)),
        path_top_k=int(payload.get("path_top_k", 500)),
        path_z_min=float(payload.get("path_z_min", 1.65)),
    )

    shared_pool = filter_top_shared_pooled(
        pooled_df,
        top_n_per_hop=int(payload.get("top_n_shared_per_hop", 15)),
    )
    top_paths = filter_top_paths_stratified(
        pooled_df,
        top_n=int(payload.get("top_n_paths", 40)),
    )

    return {
        "target_id": target_id,
        "n_genes": len(gene_ids),
        "gene_ids": gene_ids,
        "pooled_metapaths": pooled_mps,
        "pool_z_used": z_used,
        "metapath_ranking": _df_records(z_df),
        "pooled_paths": _df_records(pooled_df),
        "top_shared_paths": _df_records(shared_pool),
        "top_paths_by_score": _df_records(top_paths),
    }


def list_biological_processes(q: str = "", limit: int = 100) -> list[dict[str, str]]:
    df = load_bp_nodes()
    if q:
        q_lower = q.lower()
        mask = (
            df["name"].astype(str).str.lower().str.contains(q_lower)
            | df["identifier"].astype(str).str.lower().str.contains(q_lower)
        )
        df = df.loc[mask]
    return _df_records(df[["identifier", "name"]].head(limit))


def list_genes(q: str = "", limit: int = 100) -> list[dict[str, str]]:
    df = load_gene_nodes()
    if q:
        q_lower = q.lower()
        mask = (
            df["name"].astype(str).str.lower().str.contains(q_lower)
            | df["identifier"].astype(str).str.lower().str.contains(q_lower)
        )
        df = df.loc[mask]
    result_df = df[["identifier", "name"]].head(limit).copy()
    result_df["identifier"] = result_df["identifier"].astype(str)
    return _df_records(result_df)
