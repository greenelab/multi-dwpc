"""Streamlit MVP for the multi-DWPC web tool.

Run with: ``streamlit run app.py`` (from the repo root, in the multi_dwpc env).

See ``docs/web_tool_plan_2026-04-23.md`` for scope.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from src.dwpc_direct import HetMat
from src.multi_dwpc_query import (
    DEFAULT_B,
    DEFAULT_PATH_Z_MIN,
    discover_source_target_metapaths,
    query_intermediates_and_paths,
    query_metapath_z,
)
from src.path_enumeration import NODE_TYPE_NAMES, load_node_names

TYPE_COLORS = {
    "G": "#74c476",
    "BP": "#fd8d3c",
    "CC": "#fdae6b",
    "MF": "#fdd0a2",
    "PW": "#9ecae1",
    "C": "#c994c7",
    "D": "#fb6a4a",
    "A": "#bcbddc",
}

REPO_ROOT = Path(__file__).resolve().parent
DATA_DIR = REPO_ROOT / "data"
SAMPLE_SIZE = 20
SAMPLE_SEED = 7

# Worked example: pyrimidine nucleotide catabolic process gene set (from year
# 2016-vs-2024 pipeline output) + target BP. Used to pre-populate the form.
EXAMPLE_GENES = (
    "DUT\nENTPD4\nMBD4\nNEIL1\nNEIL2\nNT5C\nNT5M\nNTHL1\nOGG1\nSMUG1\nTDG\nUNG\n"
    "DPYD\nDPYS\nTYMP\nUPB1\nUPP1\nUPP2"
)
EXAMPLE_BP_ID = "GO:0006244"
EXAMPLE_BP_NAME = "pyrimidine nucleotide catabolic process"


@st.cache_resource(show_spinner="Preloading DWPC matrices (first launch only)...")
def get_hetmat() -> HetMat:
    hetmat = HetMat(data_dir=DATA_DIR)
    # Preload all G -> BP matrices so first query skips disk I/O.
    for mp in discover_source_target_metapaths(hetmat, "G", "BP"):
        hetmat.compute_dwpc_matrix(mp)
    return hetmat


@st.cache_data
def load_bp_nodes() -> pd.DataFrame:
    df = pd.read_csv(DATA_DIR / "nodes" / "Biological Process.tsv", sep="\t")
    return df.sort_values("name").reset_index(drop=True)


@st.cache_data
def load_gene_nodes() -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / "nodes" / "Gene.tsv", sep="\t")


@st.cache_data
def load_name_maps() -> dict[str, dict[str, str]]:
    return load_node_names(REPO_ROOT)


def random_gene_text(size: int = SAMPLE_SIZE) -> str:
    rng = np.random.default_rng()
    symbols = rng.choice(load_gene_nodes()["name"].values, size=size, replace=False)
    return "\n".join(sorted(symbols))


def random_bp_name() -> str:
    rng = np.random.default_rng()
    names = load_bp_nodes()["name"].values
    return str(rng.choice(names))


def parse_gene_input(raw: str) -> list[int]:
    genes = load_gene_nodes()
    symbol_to_entrez = dict(zip(genes["name"].astype(str), genes["identifier"].astype(int)))
    entrez_set = set(genes["identifier"].astype(int))
    out: list[int] = []
    for tok in raw.replace(",", "\n").split():
        t = tok.strip()
        if not t:
            continue
        if t.lstrip("-").isdigit() and int(t) in entrez_set:
            out.append(int(t))
        elif t in symbol_to_entrez:
            out.append(symbol_to_entrez[t])
    return out


def _qualified_to_name(qid, name_maps: dict[str, dict[str, str]]) -> str:
    if not isinstance(qid, str) or ":" not in qid:
        return qid
    node_type, _, node_id = qid.partition(":")
    return name_maps.get(node_type, {}).get(node_id, qid)


def build_subpath_figure(
    paths_df: pd.DataFrame,
    name_maps: dict[str, dict[str, str]],
    title: str,
    total_genes: int | None = None,
    weight_by: str = "coverage",
) -> plt.Figure:
    """Layered diagram of surviving subpaths. Columns = hop positions; nodes stacked within column.

    Handles mixed-length paths: each row contributes only its non-null hops. Source
    column is left-aligned; target column is right-aligned at max_hops-1.

    ``weight_by``: either ``"coverage"`` (node size reflects gene coverage) or
    ``"score"`` (node size reflects summed path_score through that node). Keeps
    Top-Shared vs Top-Paths plots visually distinct even when they share a node
    set.
    """
    hop_id_cols = sorted(c for c in paths_df.columns if c.startswith("hop_") and c.endswith("_id"))
    hop_type_cols = [c[: -len("_id")] + "_type" for c in hop_id_cols]
    n_hops = len(hop_id_cols)
    if total_genes is None:
        total_genes = paths_df["gene_id"].nunique() if "gene_id" in paths_df else 0

    column_nodes: list[dict[str, str]] = [{} for _ in range(n_hops)]
    column_genes: list[dict[str, set]] = [dict() for _ in range(n_hops)]
    column_scores: list[dict[str, float]] = [dict() for _ in range(n_hops)]
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
        mp_hops = last_col
        col_map = {i: i for i in range(mp_hops)}
        col_map[mp_hops] = n_hops - 1
        prev = None
        for local_i in range(mp_hops + 1):
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

    pos: dict[tuple[int, str], tuple[float, float]] = {}
    half_widths: dict[tuple[int, str], float] = {}
    max_col = max((len(c) for c in column_nodes), default=1)
    fig, ax = plt.subplots(figsize=(2.6 * n_hops, 0.55 * max_col + 2))

    # Global denominator for score-based sizing: max summed score across all intermediate nodes.
    max_score_any_node = 0.0
    for ci in range(1, n_hops - 1):
        for qid in column_nodes[ci]:
            max_score_any_node = max(max_score_any_node, column_scores[ci].get(qid, 0.0))

    for ci, col in enumerate(column_nodes):
        def _coverage(q):
            return len(column_genes[ci].get(q, set()))
        def _score_weight(q):
            return column_scores[ci].get(q, 0.0)
        sort_key = _score_weight if weight_by == "score" else _coverage
        sorted_ids = sorted(col, key=lambda q: (-sort_key(q), _qualified_to_name(q, name_maps)))
        n_in_col = len(sorted_ids)
        for ri, qid in enumerate(sorted_ids):
            y = (n_in_col - 1) / 2 - ri
            pos[(ci, qid)] = (ci, y)
            color = TYPE_COLORS.get(col[qid], "#d9d9d9")
            name = _qualified_to_name(qid, name_maps)
            # Sizing and label: both driven by ``weight_by`` so the visual matches
            # the number shown. Floor at 0.2 so tiny values remain visible.
            if 0 < ci < n_hops - 1:
                if weight_by == "score" and max_score_any_node > 0:
                    frac = _score_weight(qid) / max_score_any_node
                    label = f"{name} ({100 * frac:.0f}% score)"
                elif weight_by == "coverage" and total_genes:
                    frac = _coverage(qid) / total_genes
                    label = f"{name} ({100 * frac:.0f}% genes)"
                else:
                    frac = 1.0
                    label = name
                scale = max(0.2, min(1.0, frac))
            else:
                scale = 1.0
                label = name
            box_w = 0.9 * scale
            box_h = 0.44 * scale
            half_widths[(ci, qid)] = box_w / 2
            ax.add_patch(mpatches.FancyBboxPatch(
                (ci - box_w / 2, y - box_h / 2), box_w, box_h,
                boxstyle="round,pad=0.02", linewidth=0.5,
                edgecolor="black", facecolor=color,
            ))
            ax.text(ci, y, label, ha="center", va="center", fontsize=7)

    for (ci, cid, cj, jid) in edges:
        x0, y0 = pos[(ci, cid)]
        x1, y1 = pos[(cj, jid)]
        hw_src = half_widths.get((ci, cid), 0.45)
        hw_dst = half_widths.get((cj, jid), 0.45)
        ax.annotate(
            "", xy=(x1 - hw_dst, y1), xytext=(x0 + hw_src, y0),
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.5, alpha=0.5),
        )

    ax.set_xlim(-0.8, n_hops - 0.2)
    ax.set_ylim(-max_col / 2 - 0.7, max_col / 2 + 0.7)
    ax.set_xticks(range(n_hops))
    ax.set_xticklabels(
        [sorted(set(col.values()))[0] if col else "" for col in column_nodes]
    )
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title(title)

    types_present = {t for col in column_nodes for t in col.values() if t}
    if types_present:
        handles = [
            mpatches.Patch(facecolor=TYPE_COLORS.get(t, "#d9d9d9"), edgecolor="black",
                           label=NODE_TYPE_NAMES.get(t, t))
            for t in sorted(types_present)
        ]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=8, frameon=False, title="Node type")
    fig.tight_layout()
    return fig


POOL_Z_MIN = 1.65  # matches HPC end-to-end --effect-size-threshold default
DRILLDOWN_PATH_TOP_K = 500
OVERALL_TOP_PATHS = 40  # matches plot_lv_top_paths_subgraph top_n_paths default
OVERALL_TOP_SHARED_PER_HOP = 15  # matches plot_lv_top_shared_subgraph top_n_per_hop default


def _per_row_intermediates(paths_df: pd.DataFrame) -> list[list[str]]:
    """For each row, return the intermediate-hop qualified IDs (excludes source and the row's own target)."""
    hop_id_cols = sorted(c for c in paths_df.columns if c.startswith("hop_") and c.endswith("_id"))
    out: list[list[str]] = []
    for _, row in paths_df.iterrows():
        vals = [row[c] for c in hop_id_cols]
        last = next((i for i in range(len(vals) - 1, -1, -1) if isinstance(vals[i], str)), -1)
        out.append([vals[i] for i in range(1, last) if isinstance(vals[i], str)])
    return out


def filter_top_shared_pooled(
    paths_df: pd.DataFrame, top_n_per_hop: int = OVERALL_TOP_SHARED_PER_HOP
) -> pd.DataFrame:
    """Per-hop top-N intermediates by gene coverage, matching
    ``plot_lv_top_shared_subgraph``. Rows from metapaths with zero intermediate
    hops pass trivially (no hops to filter)."""
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

    per_hop_genes: dict[int, dict[str, set]] = {}
    for (hops, iids), gene_id in zip(per_row, paths_df["gene_id"]):
        for h, iid in zip(hops, iids):
            per_hop_genes.setdefault(h, {}).setdefault(iid, set()).add(int(gene_id))

    allowed: dict[int, set] = {}
    for h, m in per_hop_genes.items():
        ranked = sorted(m.items(), key=lambda kv: len(kv[1]), reverse=True)
        allowed[h] = {iid for iid, _ in ranked[:top_n_per_hop]}

    keep = [
        all(iid in allowed.get(h, set()) for h, iid in zip(hops, iids))
        for hops, iids in per_row
    ]
    return paths_df[keep].reset_index(drop=True)


def filter_top_paths_stratified(paths_df: pd.DataFrame, top_n: int = OVERALL_TOP_PATHS) -> pd.DataFrame:
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


def filter_top_shared(paths_df: pd.DataFrame, top_n_per_hop: int = 15) -> pd.DataFrame:
    """Keep paths where each intermediate hop is in that hop's top-N by gene coverage."""
    hop_id_cols = sorted(c for c in paths_df.columns if c.startswith("hop_") and c.endswith("_id"))
    intermediate_hops = hop_id_cols[1:-1]
    if not intermediate_hops or paths_df.empty:
        return paths_df
    allowed: dict[str, set] = {}
    for col in intermediate_hops:
        coverage = paths_df.groupby(col)["gene_id"].nunique().sort_values(ascending=False)
        allowed[col] = set(coverage.head(top_n_per_hop).index)
    mask = np.ones(len(paths_df), dtype=bool)
    for col, ids in allowed.items():
        mask &= paths_df[col].isin(ids).to_numpy()
    return paths_df.loc[mask].reset_index(drop=True)


def filter_top_paths(paths_df: pd.DataFrame, top_n: int = 40) -> pd.DataFrame:
    if paths_df.empty or "path_score" not in paths_df.columns:
        return paths_df
    return paths_df.nlargest(top_n, "path_score").reset_index(drop=True)


def build_intermediate_heatmap(
    paths_df: pd.DataFrame,
    name_maps: dict[str, dict[str, str]],
    metapath: str,
    top_n_intermediates: int = 30,
) -> plt.Figure | None:
    """Binary heatmap: intermediates (rows) x genes (cols), 1 if gene uses intermediate."""
    if paths_df.empty:
        return None
    per_row = _per_row_intermediates(paths_df)
    flat = [(g, iid) for g, ints in zip(paths_df["gene_id"], per_row) for iid in ints]
    if not flat:
        return None
    long = pd.DataFrame(flat, columns=["gene_id", "intermediate_id"])
    coverage = long.groupby("intermediate_id")["gene_id"].nunique().sort_values(ascending=False)
    top_ids = coverage.head(top_n_intermediates).index.tolist()
    long = long[long["intermediate_id"].isin(top_ids)]
    pivot = (
        long.assign(value=1)
        .pivot_table(index="intermediate_id", columns="gene_id", values="value", aggfunc="max", fill_value=0)
        .reindex(top_ids)
    )
    gene_names = name_maps.get("G", {})
    col_labels = [gene_names.get(str(g), str(g)) for g in pivot.columns]
    row_labels = [_qualified_to_name(q, name_maps) for q in pivot.index]
    fig, ax = plt.subplots(figsize=(max(4, 0.35 * len(col_labels) + 2), max(3, 0.28 * len(row_labels) + 1)))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Intermediate")
    ax.set_title(f"Intermediate sharing: {metapath}")
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, ticks=[0, 1])
    cbar.ax.set_yticklabels(["absent", "uses\nintermediate"], fontsize=7)
    fig.tight_layout()
    return fig


def rename_ids_to_names(df: pd.DataFrame, name_maps: dict[str, dict[str, str]]) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    if "intermediate_id" in out.columns:
        out.insert(0, "intermediate", out["intermediate_id"].map(lambda s: _qualified_to_name(s, name_maps)))
        out = out.drop(columns=["intermediate_id"])
    if "gene_ids" in out.columns:
        gene_names = name_maps.get("G", {})
        out["genes"] = out["gene_ids"].map(
            lambda s: ", ".join(gene_names.get(g, g) for g in str(s).split(","))
        )
        out = out.drop(columns=["gene_ids"])
    if "gene_id" in out.columns:
        gene_names = name_maps.get("G", {})
        out.insert(0, "gene", out["gene_id"].astype(str).map(lambda g: gene_names.get(g, g)))
        out = out.drop(columns=["gene_id"])
    hop_cols = sorted(c for c in out.columns if c.startswith("hop_") and c.endswith("_id"))
    for col in hop_cols:
        prefix = col[: -len("_id")]
        out[prefix] = out[col].map(lambda s: _qualified_to_name(s, name_maps))
        drop = [col]
        if f"{prefix}_type" in out.columns:
            drop.append(f"{prefix}_type")
        out = out.drop(columns=drop)
    return out


st.set_page_config(page_title="Multi-DWPC", layout="wide")
st.title("Multi-DWPC query tool")
st.caption("Gene list -> target -> metapath z, intermediate sharing, surviving subpaths")

bp_df = load_bp_nodes()
name_maps = load_name_maps()
bp_name_to_id = dict(zip(bp_df["name"], bp_df["identifier"]))

if "genes_raw" not in st.session_state:
    st.session_state["genes_raw"] = EXAMPLE_GENES
if "target_name_choice" not in st.session_state:
    st.session_state["target_name_choice"] = EXAMPLE_BP_NAME

bp_options = bp_df["name"].tolist()

with st.sidebar:
    st.header("Query")
    target_name = st.selectbox(
        "Target Biological Process", bp_options,
        placeholder="Search...", key="target_name_choice",
    )
    target_id = bp_name_to_id.get(target_name, "") if target_name else ""
    def _fill_example_query():
        st.session_state["target_name_choice"] = EXAMPLE_BP_NAME
        st.session_state["genes_raw"] = EXAMPLE_GENES

    def _fill_random_bp():
        st.session_state["target_name_choice"] = random_bp_name()

    def _fill_random_genes():
        st.session_state["genes_raw"] = random_gene_text()

    col_a, col_b, col_c = st.columns(3)
    with col_a:
        st.button("Example query", on_click=_fill_example_query,
                  help="Pyrimidine nucleotide catabolism worked example")
    with col_b:
        st.button("Random BP", on_click=_fill_random_bp)
    with col_c:
        st.button("Random genes", on_click=_fill_random_genes)
    genes_raw = st.text_area(
        "Gene symbols (one per line or comma-separated; Entrez IDs also accepted)",
        height=200,
        key="genes_raw",
    )
    run = st.button("Run query", type="primary", disabled=not (target_id and genes_raw.strip()))

if run:
    gene_ids = parse_gene_input(genes_raw)
    if not gene_ids:
        st.error("No valid gene symbols or Entrez IDs parsed from input.")
        st.stop()
    hetmat = get_hetmat()
    with st.spinner(f"Computing metapath z for {len(gene_ids)} genes (b={DEFAULT_B})"):
        z_df = query_metapath_z(gene_ids, target_id, hetmat=hetmat)
    st.session_state["z_df"] = z_df
    st.session_state["gene_ids"] = gene_ids
    st.session_state["target_id"] = target_id
    st.session_state["_drilldown_cache"] = {}

z_df = st.session_state.get("z_df")
if z_df is None:
    st.info("Enter a target and gene list in the sidebar, then run the query.")
    st.stop()

with st.sidebar:
    st.header("Drill-down metapath")
    metapath_choice = st.selectbox("Metapath", z_df["metapath"].tolist(), key="mp_choice")

tab_z, tab_int, tab_overall, tab_metapath = st.tabs(
    ["Metapath ranking", "Intermediate sharing", "Overall subgraphs", "Metapath subgraphs"]
)

with tab_z:
    st.dataframe(z_df, use_container_width=True)


def _load_drilldown(gene_ids, target_id, metapath):
    key = (tuple(gene_ids), target_id, metapath)
    cache = st.session_state.setdefault("_drilldown_cache", {})
    if key not in cache:
        hetmat = get_hetmat()
        cache[key] = query_intermediates_and_paths(
            gene_ids=gene_ids,
            target_id=target_id,
            metapath=metapath,
            repo_root=REPO_ROOT,
            hetmat=hetmat,
            path_top_k=DRILLDOWN_PATH_TOP_K,
            path_z_min=DEFAULT_PATH_Z_MIN,
        )
    return cache[key]


with st.spinner(f"Enumerating paths for {metapath_choice}"):
    sharing_df, paths_df, diagnostics = _load_drilldown(
        st.session_state["gene_ids"], st.session_state["target_id"], metapath_choice,
    )

with tab_int:
    if sharing_df.empty or paths_df.empty:
        st.warning("No intermediates survived the path-z filter.")
    else:
        heatmap = build_intermediate_heatmap(paths_df, name_maps, metapath_choice)
        if heatmap is None:
            st.info("No intermediate hops for this metapath.")
        else:
            st.pyplot(heatmap)
        with st.expander("Intermediate sharing table"):
            st.dataframe(rename_ids_to_names(sharing_df, name_maps), use_container_width=True)


def _pool_across_metapaths(gene_ids, target_id, z_df) -> tuple[pd.DataFrame, list[str], float]:
    """Pool at POOL_Z_MIN (1.65); fall back to z > 0 only if nothing clears. Returns
    (pooled_df, metapath_list, z_used)."""
    z_used = POOL_Z_MIN
    pool_mps = z_df.loc[z_df["effect_size_z"] >= POOL_Z_MIN, "metapath"].tolist()
    if not pool_mps:
        z_used = 0.0
        pool_mps = z_df.loc[z_df["effect_size_z"] > 0, "metapath"].tolist()
    frames: list[pd.DataFrame] = []
    for mp in pool_mps:
        _, mp_paths, _ = _load_drilldown(gene_ids, target_id, mp)
        if not mp_paths.empty:
            frames.append(mp_paths)
    if not frames:
        return pd.DataFrame(), pool_mps, z_used
    return pd.concat(frames, ignore_index=True), pool_mps, z_used


with tab_overall:
    with st.spinner("Pooling paths across significant metapaths"):
        pooled_df, pooled_mps, z_used = _pool_across_metapaths(
            st.session_state["gene_ids"], st.session_state["target_id"], z_df
        )
    if pooled_df.empty:
        st.warning(
            f"No paths across metapaths with z >= {POOL_Z_MIN:.2f} (nor z > 0)."
        )
    else:
        n_genes = len(st.session_state["gene_ids"])
        fallback = " (fallback: nothing cleared z>=1.65)" if z_used < POOL_Z_MIN else ""
        st.caption(
            f"Pooled over {len(pooled_mps)} metapaths at z >= {z_used:.2f}{fallback}: "
            + ", ".join(pooled_mps)
        )
        st.subheader("Top shared intermediates")
        st.caption("Node size = number of genes using that intermediate.")
        shared_pool = filter_top_shared_pooled(pooled_df)
        if shared_pool.empty:
            st.info("No paths survived the top-shared filter.")
        else:
            st.pyplot(build_subpath_figure(
                shared_pool, name_maps,
                f"Top shared intermediates; {len(pooled_mps)} metapaths",
                total_genes=n_genes, weight_by="coverage",
            ))

        st.subheader("Top paths by score")
        st.caption("Node size = summed path score through that intermediate.")
        top_pool = filter_top_paths_stratified(pooled_df)
        if top_pool.empty:
            st.info("No paths have path_score populated.")
        else:
            st.pyplot(build_subpath_figure(
                top_pool, name_maps,
                f"Top {len(top_pool)} paths by score; {len(pooled_mps)} metapaths",
                total_genes=n_genes, weight_by="score",
            ))

with tab_metapath:
    if paths_df.empty:
        st.warning(f"No paths survived the path-z filter (>= {DEFAULT_PATH_Z_MIN:.2f}).")
        with st.expander("Diagnostics"):
            st.json(diagnostics)
    else:
        n_genes = len(st.session_state["gene_ids"])
        st.subheader(f"Top paths by score: {metapath_choice}")
        st.caption(
            "Node size = summed path score through that intermediate. "
            "Shared-intermediate information for this metapath is shown as a heatmap "
            "in the *Intermediate sharing* tab."
        )
        top_df = filter_top_paths(paths_df)
        if top_df.empty:
            st.info("No paths have path_score populated.")
        else:
            st.pyplot(build_subpath_figure(
                top_df, name_maps,
                f"{metapath_choice}: top {len(top_df)} paths by score",
                total_genes=n_genes, weight_by="score",
            ))

        with st.expander("Table of surviving subpaths"):
            st.dataframe(rename_ids_to_names(paths_df, name_maps), use_container_width=True)
