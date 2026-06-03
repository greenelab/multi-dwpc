"""Graph traversal utilities for hetionet path enumeration.

Provides sparse-matrix-based path enumeration through the hetionet knowledge
graph.  Extracted from the former ``extract_top_paths_local.py`` script and
the duplicate copies that lived in ``src/lv_subgraphs.py``.
"""

from __future__ import annotations

import heapq
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

from src.dwpc_direct import reverse_metapath_abbrev  # noqa: F401 — re-exported

NODE_FILES: dict[str, str] = {
    "G": "Gene.tsv",
    "BP": "Biological Process.tsv",
    "CC": "Cellular Component.tsv",
    "MF": "Molecular Function.tsv",
    "PW": "Pathway.tsv",
    "A": "Anatomy.tsv",
    "D": "Disease.tsv",
    "C": "Compound.tsv",
    "S": "Symptom.tsv",
    "SE": "Side Effect.tsv",
    "PC": "Pharmacologic Class.tsv",
}

NODE_TYPE_NAMES: dict[str, str] = {
    "G": "Gene",
    "A": "Anatomy",
    "BP": "Biological Process",
    "CC": "Cellular Component",
    "C": "Compound",
    "D": "Disease",
    "MF": "Molecular Function",
    "PC": "Pharmacologic Class",
    "PW": "Pathway",
    "SE": "Side Effect",
    "S": "Symptom",
}


def parse_metapath(metapath: str) -> tuple[list[str], list[str]]:
    """Parse a metapath abbreviation into node-type and edge-token lists."""
    node_abbrevs = sorted(NODE_FILES.keys(), key=len, reverse=True)
    nodes: list[str] = []
    edges: list[str] = []
    i = 0
    edge_token = ""
    while i < len(metapath):
        matched = False
        for ab in node_abbrevs:
            if metapath.startswith(ab, i):
                if nodes:
                    edges.append(edge_token)
                    edge_token = ""
                nodes.append(ab)
                i += len(ab)
                matched = True
                break
        if not matched:
            edge_token += metapath[i]
            i += 1
    return nodes, edges


@dataclass
class NodeMaps:
    id_to_pos: Dict[str, Dict[str, int]]
    pos_to_id: Dict[str, Dict[int, str]]
    id_to_name: Dict[str, Dict[str, str]]


def load_node_maps(repo_root: Path, node_types: list[str]) -> NodeMaps:
    """Load position / ID / name mappings for the requested node types."""
    id_to_pos: dict[str, dict[str, int]] = {}
    pos_to_id: dict[str, dict[int, str]] = {}
    id_to_name: dict[str, dict[str, str]] = {}
    for node_type in node_types:
        filename = NODE_FILES[node_type]
        df = pd.read_csv(repo_root / "data" / "nodes" / filename, sep="\t")
        id_to_pos[node_type] = dict(zip(df["identifier"], df["position"]))
        pos_to_id[node_type] = dict(zip(df["position"], df["identifier"]))
        id_to_name[node_type] = dict(zip(df["identifier"], df["name"]))
    return NodeMaps(id_to_pos, pos_to_id, id_to_name)


def load_node_names(repo_root: Path) -> dict[str, dict[str, str]]:
    """Load node ID-to-name mappings for all known node types."""
    nodes_dir = repo_root / "data" / "nodes"
    name_maps: dict[str, dict[str, str]] = {}
    for abbrev, full_name in NODE_TYPE_NAMES.items():
        node_file = nodes_dir / f"{full_name}.tsv"
        if not node_file.exists():
            continue
        df = pd.read_csv(node_file, sep="\t")
        name_maps[abbrev] = dict(zip(df["identifier"].astype(str), df["name"]))
    return name_maps


class EdgeLoader:
    """Loads and caches sparse adjacency matrices from ``data/edges/``."""

    def __init__(self, edges_dir: Path):
        self.edges_dir = edges_dir
        self.edge_files = {p.name for p in edges_dir.glob("*.sparse.npz")}
        self.cache: dict[tuple, sparse.csr_matrix] = {}
        self.deg_cache: dict[tuple, np.ndarray] = {}

    def _canonical_edge_token(self, edge_token: str) -> Tuple[str, bool]:
        reverse = "<" in edge_token and ">" not in edge_token
        edge_base = edge_token.replace("<", "").replace(">", "")
        if "<" in edge_token or ">" in edge_token:
            file_edge = f"{edge_base}>"
        else:
            file_edge = edge_base
        return file_edge, reverse

    def load(self, src: str, edge_token: str, dst: str) -> sparse.csr_matrix:
        key = (src, edge_token, dst)
        if key in self.cache:
            return self.cache[key]

        file_edge, reverse = self._canonical_edge_token(edge_token)
        name_direct = f"{src}{file_edge}{dst}.sparse.npz"
        name_reverse = f"{dst}{file_edge}{src}.sparse.npz"

        if name_direct in self.edge_files:
            mat = sparse.load_npz(self.edges_dir / name_direct)
            if reverse:
                mat = mat.T
        elif name_reverse in self.edge_files:
            mat = sparse.load_npz(self.edges_dir / name_reverse).T
            if reverse:
                mat = mat.T
        else:
            raise FileNotFoundError(f"No edge matrix for {src}{edge_token}{dst}")

        mat = mat.tocsr()
        self.cache[key] = mat
        return mat

    def degree(self, src: str, edge_token: str, dst: str) -> np.ndarray:
        key = (src, edge_token, dst)
        if key in self.deg_cache:
            return self.deg_cache[key]
        mat = self.load(src, edge_token, dst)
        deg = np.asarray(mat.getnnz(axis=1)).astype(float)
        self.deg_cache[key] = deg
        return deg


def enumerate_paths(
    source_pos: int,
    target_pos: int,
    nodes: list[str],
    edges: list[str],
    edge_loader: EdgeLoader,
    top_k: int,
    degree_d: float = 0.5,
) -> list[tuple[float, list[int]]]:
    """Enumerate top-k weighted paths for metapaths up to length 3."""
    num_edges = len(edges)
    if num_edges == 1:
        A0 = edge_loader.load(nodes[0], edges[0], nodes[1])
        if A0[source_pos, target_pos] == 0:
            return []
        return [(1.0, [source_pos, target_pos])]

    A0 = edge_loader.load(nodes[0], edges[0], nodes[1])
    neighbors1 = A0.getrow(source_pos).indices

    if num_edges == 2:
        A1 = edge_loader.load(nodes[1], edges[1], nodes[2])
        deg1 = edge_loader.degree(nodes[1], edges[1], nodes[2])
        targets = set(A1[:, target_pos].nonzero()[0])
        heap: list[tuple[float, list[int]]] = []
        for n1 in neighbors1:
            if n1 in targets:
                score = 1.0 / math.sqrt(max(deg1[n1], 1.0)) ** degree_d
                heapq.heappush(heap, (score, [source_pos, n1, target_pos]))
        return sorted(heap, key=lambda x: x[0], reverse=True)[:top_k]

    if num_edges == 3:
        A1 = edge_loader.load(nodes[1], edges[1], nodes[2])
        A2 = edge_loader.load(nodes[2], edges[2], nodes[3])
        deg1 = edge_loader.degree(nodes[1], edges[1], nodes[2])
        deg2 = edge_loader.degree(nodes[2], edges[2], nodes[3])
        targets = set(A2[:, target_pos].nonzero()[0])

        heap = []
        for n1 in neighbors1:
            n2s = A1.getrow(n1).indices
            if len(n2s) == 0:
                continue
            for n2 in n2s:
                if n2 in targets:
                    score = 1.0 / (
                        (math.sqrt(max(deg1[n1], 1.0)) ** degree_d)
                        * (math.sqrt(max(deg2[n2], 1.0)) ** degree_d)
                    )
                    if len(heap) < top_k:
                        heapq.heappush(heap, (score, [source_pos, n1, n2, target_pos]))
                    elif score > heap[0][0]:
                        heapq.heapreplace(heap, (score, [source_pos, n1, n2, target_pos]))
        return sorted(heap, key=lambda x: x[0], reverse=True)

    raise ValueError(f"Metapath length {num_edges} not supported for enumeration.")


def select_paths(
    paths: list[tuple[float, list[int]]],
    *,
    selection_method: str,
    top_paths: int,
    path_cumulative_frac: float | None,
    path_min_count: int,
    path_max_count: int | None,
) -> list[tuple[float, list[int]]]:
    """Select paths by effective number, cumulative fraction, or top-N."""
    if selection_method == "effective_number":
        vals = np.array([max(float(s), 0.0) for s, _ in paths], dtype=float)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if vals.size == 0:
            k = 1
        else:
            weights = vals / vals.sum()
            entropy = float(-(weights * np.log(weights)).sum())
            k = int(np.ceil(np.exp(entropy)))
        k = max(int(path_min_count), k)
        if path_max_count is not None:
            k = min(k, int(path_max_count))
        return paths[:k]

    if path_cumulative_frac is None:
        return paths[:top_paths]

    if not paths:
        return []
    total = sum(max(float(s), 0.0) for s, _ in paths)
    cumulative = 0.0
    selected: list[tuple[float, list[int]]] = []
    for item in paths:
        selected.append(item)
        cumulative += max(float(item[0]), 0.0)
        enough_by_frac = total <= 0 or cumulative / total >= float(path_cumulative_frac)
        enough_by_min = len(selected) >= int(path_min_count)
        enough_by_max = path_max_count is not None and len(selected) >= int(path_max_count)
        if enough_by_max or (enough_by_min and enough_by_frac):
            break
    return selected
