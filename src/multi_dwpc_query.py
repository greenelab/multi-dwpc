"""Web-tool adapter: user gene list + target -> metapath z, intermediates, subpaths.

See ``docs/web_tool_plan_2026-04-23.md`` for scope. MVP is Gene -> BP only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.dwpc_direct import HetMat, reverse_metapath_abbrev
from src.intermediate_sharing import enumerate_gene_intermediates
from src.path_enumeration import (
    EdgeLoader, NODE_FILES, NODE_TYPE_NAMES, load_node_maps, parse_metapath,
)

USER_QUERY_ID = "user_query"
DEFAULT_B = 20
DEFAULT_PATH_Z_MIN = 1.65


def _gene_index_map(hetmat: HetMat) -> dict[int, int]:
    nodes = hetmat.get_nodes(NODE_TYPE_NAMES["G"])
    return dict(zip(nodes["identifier"].astype(int), nodes["position"].astype(int)))


def _valid_metapaths(hetmat: HetMat, candidates: list[str]) -> list[str]:
    """Drop metapaths whose edge abbreviations aren't in the metagraph."""
    mg = hetmat._hetmatpy.metagraph
    valid = []
    for mp in candidates:
        try:
            mg.metapath_from_abbrev(mp)
        except KeyError:
            continue
        valid.append(mp)
    return valid


def discover_source_target_metapaths(
    hetmat: HetMat, source_type: str, target_type: str
) -> list[str]:
    """All metapaths from metapath-dwpc-stats.tsv between the two types,
    oriented source-first. Validated against the hetmat metagraph."""
    stats = hetmat.metapath_stats
    node_abbrevs = sorted(NODE_FILES.keys(), key=len, reverse=True)

    def tokenize(mp: str) -> list[str]:
        s = mp.replace("<", "").replace(">", "")
        out: list[str] = []
        i = 0
        while i < len(s):
            for ab in node_abbrevs:
                if s.startswith(ab, i):
                    out.append(ab)
                    i += len(ab)
                    break
            else:
                i += 1
        return out

    oriented: list[str] = []
    for mp in stats["metapath"].astype(str).tolist():
        toks = tokenize(mp)
        if not toks:
            continue
        if toks[0] == source_type and toks[-1] == target_type:
            oriented.append(mp)
        elif toks[0] == target_type and toks[-1] == source_type:
            oriented.append(reverse_metapath_abbrev(mp))
    return _valid_metapaths(hetmat, sorted(set(oriented)))


def _target_position(hetmat: HetMat, target_type: str, target_id: str) -> int:
    nodes = hetmat.get_nodes(NODE_TYPE_NAMES[target_type])
    hit = nodes.loc[nodes["identifier"].astype(str) == str(target_id), "position"]
    if len(hit) == 0:
        raise ValueError(f"Target {target_id!r} not found among {target_type} nodes")
    return int(hit.iloc[0])


def query_metapath_z(
    gene_ids: list[int],
    target_id: str,
    target_type: str = "BP",
    *,
    b: int = DEFAULT_B,
    seed: int = 42,
    metapaths: list[str] | None = None,
    hetmat: HetMat,
) -> pd.DataFrame:
    """Return per-metapath z-score for the user's gene set against ``target_id``.

    Null model: ``b`` random gene subsets of the same size drawn from the full
    gene universe. z = (real - null_mean) / null_std, matching
    ``build_b_seed_runs``. Differs from the paper's degree-preserving XSWAP
    null; see docs/web_tool_plan_2026-04-23.md.
    """
    if metapaths is None:
        metapaths = discover_source_target_metapaths(hetmat, "G", target_type)
    if not metapaths:
        raise ValueError("No valid metapaths for this source/target type")

    gene_idx_map = _gene_index_map(hetmat)
    source_idx = np.array(
        [gene_idx_map[g] for g in gene_ids if g in gene_idx_map], dtype=np.int64
    )
    if source_idx.size == 0:
        raise ValueError("None of the provided gene IDs were found in Hetionet")

    gene_pool = np.fromiter(gene_idx_map.values(), dtype=np.int64)
    target_pos = _target_position(hetmat, target_type, target_id)
    target_idx_arr = np.full(source_idx.size, target_pos, dtype=np.int64)

    rng = np.random.default_rng(seed)
    perm_idx = np.stack(
        [rng.choice(gene_pool, size=source_idx.size, replace=False) for _ in range(b)]
    )

    rows = []
    for mp in metapaths:
        real_mean = float(hetmat.get_dwpc_for_pairs(mp, source_idx, target_idx_arr).mean())
        null_means = np.array(
            [
                hetmat.get_dwpc_for_pairs(
                    mp, perm_idx[i], np.full(source_idx.size, target_pos, dtype=np.int64)
                ).mean()
                for i in range(b)
            ]
        )
        null_mean = float(null_means.mean())
        null_std = float(null_means.std(ddof=1))
        z = (real_mean - null_mean) / null_std if null_std > 0 else np.nan
        rows.append(
            {
                "metapath": mp,
                "real_mean_score": real_mean,
                "null_mean_score": null_mean,
                "null_std_score": null_std,
                "diff": real_mean - null_mean,
                "effect_size_z": z,
            }
        )
    return pd.DataFrame(rows).sort_values("effect_size_z", ascending=False, ignore_index=True)


def query_intermediates_and_paths(
    gene_ids: list[int],
    target_id: str,
    metapath: str,
    *,
    repo_root: Path,
    hetmat: HetMat,
    target_type: str = "BP",
    path_top_k: int = 100,
    path_z_min: float = DEFAULT_PATH_Z_MIN,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return ``(intermediates_df, paths_df, diagnostics)`` for a single metapath.

    ``intermediates_df`` rows: intermediate_id, n_genes_sharing, gene_ids.
    ``paths_df`` rows: one per surviving (path_z >= path_z_min) subpath.
    ``diagnostics``: counts / pool stats for debugging empty results.

    The pipeline's DWPC gene-level z-filter is intentionally not applied here:
    the user's gene list is the query, not a pool to screen. Only the path-z
    filter is retained.
    """
    node_types = sorted({"G", target_type, *_metapath_nodes(metapath)})
    maps = load_node_maps(repo_root, node_types)
    edge_loader = EdgeLoader(repo_root / "data" / "edges")

    gene_idx_map = _gene_index_map(hetmat)
    valid_gene_ids = [int(g) for g in gene_ids if int(g) in gene_idx_map]
    target_pos = _target_position(hetmat, target_type, target_id)

    # Manual pool-score diagnostic: call enumerate_paths directly for each gene.
    from src.path_enumeration import enumerate_paths, parse_metapath
    from src.dwpc_direct import reverse_metapath_abbrev
    reversed_mp = reverse_metapath_abbrev(metapath)
    nodes, edges = parse_metapath(reversed_mp)
    gene_id_map = maps.id_to_pos.get("G", {})
    pool_scores: list[float] = []
    genes_with_paths = 0
    for gid in valid_gene_ids:
        gpos = gene_id_map.get(gid) or gene_id_map.get(str(gid))
        if gpos is None:
            continue
        try:
            paths = enumerate_paths(
                target_pos, gpos, nodes, edges, edge_loader,
                top_k=path_top_k, degree_d=0.5,
            )
        except Exception:
            continue
        if paths:
            genes_with_paths += 1
            pool_scores.extend(s for s, _ in paths)

    pool_mean = float(np.mean(pool_scores)) if len(pool_scores) >= 2 else 0.0
    pool_std = float(np.std(pool_scores, ddof=1)) if len(pool_scores) >= 2 else 0.0
    path_cutoff = pool_mean + path_z_min * pool_std if pool_std > 0 else None
    n_surviving = sum(1 for s in pool_scores if path_cutoff is None or s >= path_cutoff)

    record_paths: list[dict] = []
    intermediates, _filtered = enumerate_gene_intermediates(
        genes=valid_gene_ids,
        target_pos=target_pos,
        metapath=metapath,
        edge_loader=edge_loader,
        maps=maps,
        path_top_k=path_top_k,
        gene_set_id=USER_QUERY_ID,
        dwpc_lookup=None,
        dwpc_z_stats=None,
        path_z_min=path_z_min,
        record_paths=record_paths,
    )

    diagnostics = {
        "n_input_genes": len(gene_ids),
        "n_genes_in_hetmat": len(valid_gene_ids),
        "n_genes_with_paths": genes_with_paths,
        "n_paths_total": len(pool_scores),
        "path_pool_mean": pool_mean,
        "path_pool_std": pool_std,
        "path_score_cutoff": path_cutoff,
        "n_paths_surviving_cutoff": n_surviving,
        "n_records_actually_produced": len(record_paths),
    }

    sharing_rows = []
    inv: dict[str, list[int]] = {}
    for gene_id, intermediate_set in intermediates.items():
        for intermediate in intermediate_set:
            inv.setdefault(intermediate, []).append(int(gene_id))
    for intermediate, genes in sorted(inv.items(), key=lambda kv: -len(kv[1])):
        sharing_rows.append(
            {
                "intermediate_id": intermediate,
                "n_genes_sharing": len(genes),
                "gene_ids": ",".join(str(g) for g in sorted(genes)),
            }
        )
    sharing_df = pd.DataFrame(sharing_rows)
    paths_df = pd.DataFrame(record_paths)
    return sharing_df, paths_df, diagnostics


def _metapath_nodes(metapath: str) -> list[str]:
    from src.path_enumeration import parse_metapath

    nodes, _ = parse_metapath(metapath)
    return nodes
