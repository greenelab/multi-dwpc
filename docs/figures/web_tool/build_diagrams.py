"""Generate editable SVG flowcharts for the multi-DWPC web tool requirements.

SVG keeps text live (svg.fonttype=none -> <text> elements, not paths) so
diagrams remain editable in Adobe Illustrator while embedding inline in the
markdown requirements doc. Run from repo root:

    python docs/figures/web_tool/build_diagrams.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.fonttype"] = "none"
# Explicit font stack so the SVG carries `font-family: Helvetica, Arial, ...`
# instead of matplotlib's default DejaVu Sans, which Illustrator substitutes.
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "sans-serif"]

OUT_DIR = Path(__file__).resolve().parent

REUSE = "#bcd9f0"      # blue: reuse from connectivity-search
NEW = "#fdd0a2"        # orange: new for multi-dwpc
DATA = "#d9d9d9"       # grey: data store
USER = "#c7e9c0"       # green: user-facing
EDGE = "#525252"


def _box(ax, x, y, w, h, label, color, fontsize=9, weight="normal"):
    ax.add_patch(mpatches.FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.04", linewidth=0.7,
        edgecolor="black", facecolor=color,
    ))
    ax.text(x, y, label, ha="center", va="center",
            fontsize=fontsize, fontweight=weight, wrap=True)


def _arrow(ax, x0, y0, x1, y1, label=None, style="->", lw=0.8, ls="-"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle=style, color=EDGE,
                                lw=lw, linestyle=ls))
    if label:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.08, label,
                ha="center", va="bottom", fontsize=7, color=EDGE,
                style="italic")


def architecture_diagram() -> None:
    fig, ax = plt.subplots(figsize=(11, 7.5))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 7.5)
    ax.set_axis_off()

    # Layer headings
    for y, label in [(7.0, "Browser"), (5.6, "API"),
                     (3.6, "Compute"), (1.2, "Data")]:
        ax.text(0.15, y, label, ha="left", va="center",
                fontsize=10, fontweight="bold", color="#333")

    # Browser layer
    _box(ax, 5.5, 7.0, 4.6, 0.7,
         "React frontend (fork of connectivity-search-frontend)",
         REUSE, fontsize=9, weight="bold")

    # API layer
    _box(ax, 2.0, 5.6, 2.6, 0.9, "Reused endpoints\nNode search\nMetanode lookup", REUSE, fontsize=8)
    _box(ax, 5.5, 5.6, 3.0, 0.9, "New endpoint\nPOST /v1/gene-set/query", NEW, weight="bold", fontsize=8)
    _box(ax, 9.0, 5.6, 3.0, 0.9, "New endpoint\nPOST /v1/gene-set/paths", NEW, weight="bold", fontsize=8)

    # Compute layer
    _box(ax, 2.0, 3.6, 2.6, 1.0,
         "Django ORM\nNode / Metanode\nMetapath lookup", REUSE, fontsize=8)
    _box(ax, 5.5, 3.9, 3.0, 0.8,
         "query_metapath_z\n(on-demand DWPC + null)", NEW, weight="bold", fontsize=8)
    _box(ax, 5.5, 2.9, 3.0, 0.8,
         "query_intermediates_and_paths\n(Cypher + scoring layer)", NEW, weight="bold", fontsize=8)

    # Data layer
    _box(ax, 2.0, 1.0, 2.6, 1.2,
         "Postgres\nNode, Metanode,\nMetapath\n(no PathCount/DGP)",
         DATA, fontsize=8)
    _box(ax, 5.5, 1.0, 3.0, 1.2,
         "HetMat data dir\nDWPC sparse matrices\n(.npz, .tsv)\n[mounted volume]",
         DATA, weight="bold", fontsize=8)
    _box(ax, 9.0, 1.0, 3.0, 1.2,
         "Neo4j\nHetionet graph\nfor Cypher path\ntraversal",
         REUSE, fontsize=8)

    # Vertical wires (browser -> API -> compute -> data)
    _arrow(ax, 5.5, 6.65, 5.5, 6.05)
    _arrow(ax, 5.5, 5.15, 5.5, 4.30)        # POST /query -> query_metapath_z
    _arrow(ax, 5.5, 3.50, 5.5, 3.30)        # query_metapath_z -> query_intermediates_and_paths
    _arrow(ax, 5.5, 2.50, 5.5, 1.62)        # compute (5.5) -> HetMat
    _arrow(ax, 2.0, 5.15, 2.0, 4.10)        # reused API -> Django ORM
    _arrow(ax, 2.0, 3.10, 2.0, 1.62)        # Django ORM -> Postgres

    # Diagonal: POST /v1/gene-set/paths -> query_intermediates_and_paths
    _arrow(ax, 9.0, 5.15, 7.0, 3.30)

    # Cross-link: query_intermediates_and_paths -> Neo4j (Cypher)
    _arrow(ax, 7.0, 2.9, 9.0, 1.62, label="Cypher")

    # Cross-link: ORM -> compute query helpers (node id lookup)
    _arrow(ax, 3.3, 3.6, 4.0, 3.6, ls=":")
    ax.text(3.65, 3.45, "node id\nlookup", fontsize=7,
            ha="center", va="top", style="italic", color=EDGE)

    # Legend
    legend = [
        mpatches.Patch(facecolor=REUSE, edgecolor="black",
                       label="Reused from connectivity-search"),
        mpatches.Patch(facecolor=NEW, edgecolor="black",
                       label="New for multi-DWPC"),
        mpatches.Patch(facecolor=DATA, edgecolor="black",
                       label="Data store"),
    ]
    ax.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.02, 0.02),
              fontsize=8, frameon=False)

    ax.text(5.5, 7.45, "Multi-DWPC web tool: system architecture",
            ha="center", va="top", fontsize=12, fontweight="bold")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "architecture.svg", bbox_inches="tight")
    plt.close(fig)


def query_flow_diagram() -> None:
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 8.5)
    ax.set_axis_off()

    ax.text(5.5, 8.3, "Multi-DWPC web tool: query flow",
            ha="center", va="top", fontsize=12, fontweight="bold")

    # Stage 1: input
    _box(ax, 5.5, 7.5, 5.2, 0.6,
         "User submits: gene list (symbols/Entrez) + target node",
         USER, fontsize=8, weight="bold")

    # Stage 2: validation (reuse)
    _box(ax, 5.5, 6.6, 5.2, 0.6,
         "Validate genes against Node table; resolve target via /v1/nodes",
         REUSE, fontsize=8)

    # Stage 3a: enumerate metapaths
    _box(ax, 5.5, 5.7, 5.2, 0.7,
         "Enumerate G->target metapaths\n(discover_source_target_metapaths)",
         NEW, fontsize=8)

    # Stage 4: real + null DWPC (parallel)
    _box(ax, 2.7, 4.5, 4.0, 0.9,
         "Real DWPC\nHetMat.get_dwpc_for_pairs\n(vectorized lookup, b=1)",
         NEW, fontsize=8)
    _box(ax, 8.3, 4.5, 4.0, 0.9,
         "Null DWPC\nb=20 random gene subsets\nsame target",
         NEW, fontsize=8)

    # Stage 5: effect size
    _box(ax, 5.5, 3.3, 5.2, 0.7,
         "z = (real - null_mean) / null_std  per metapath",
         NEW, fontsize=8)

    # Stage 6: result table -> user choice
    _box(ax, 5.5, 2.4, 5.2, 0.6,
         "Return ranked metapath table -> user selects one",
         USER, fontsize=8)

    # Stage 7a: intermediates
    _box(ax, 2.7, 1.2, 4.0, 0.9,
         "Intermediate sharing\nenumerate_gene_intermediates\n(per-gene top-K, path_z >= 1.65)",
         NEW, fontsize=8)
    # Stage 7b: subpaths
    _box(ax, 8.3, 1.2, 4.0, 0.9,
         "Surviving subpaths\n(same call as intermediates,\nrecord_paths flag)",
         NEW, fontsize=8)

    # Stage 8: visualization
    _box(ax, 5.5, 0.25, 5.2, 0.5,
         "Frontend: ranking table, sharing heatmap, subpath subgraph",
         USER, fontsize=8)

    # Arrows: top-down central spine for stages 1->2->3
    _arrow(ax, 5.5, 7.20, 5.5, 6.93)        # input -> validation
    _arrow(ax, 5.5, 6.30, 5.5, 6.05)        # validation -> enumerate
    # Stage 3 -> Stage 4 (diagonals to two parallel boxes)
    _arrow(ax, 4.0, 5.35, 3.0, 4.96)        # enumerate (left corner) -> real DWPC top
    _arrow(ax, 7.0, 5.35, 8.0, 4.96)        # enumerate (right corner) -> null DWPC top
    # Stage 4 -> Stage 5 (diagonals back to center)
    _arrow(ax, 3.5, 4.05, 4.0, 3.66)        # real DWPC -> effect size (left)
    _arrow(ax, 7.5, 4.05, 7.0, 3.66)        # null DWPC -> effect size (right)
    # Stage 5 -> Stage 6 (vertical)
    _arrow(ax, 5.5, 2.94, 5.5, 2.71)
    # Stage 6 -> Stage 7 (diagonals to two parallel boxes)
    _arrow(ax, 4.0, 2.09, 3.0, 1.66)        # ranked table -> intermediates
    _arrow(ax, 7.0, 2.09, 8.0, 1.66)        # ranked table -> subpaths
    # Stage 7 -> Stage 8 (diagonals back to center)
    _arrow(ax, 3.5, 0.74, 4.0, 0.51)        # intermediates -> viz
    _arrow(ax, 7.5, 0.74, 7.0, 0.51)        # subpaths -> viz

    # Side notes: cache key composition (matches §5 in doc)
    ax.text(0.2, 4.5, "Cache key:\nsorted source_ids,\ntarget_id, b, seed,\nmetapath_list",
            fontsize=7, ha="left", va="center", style="italic", color=EDGE)
    ax.text(0.2, 1.3, "Cache key:\nmetapath, path_top_k,\npath_z_min",
            fontsize=7, ha="left", va="center", style="italic", color=EDGE)

    legend = [
        mpatches.Patch(facecolor=USER, edgecolor="black", label="User-facing"),
        mpatches.Patch(facecolor=REUSE, edgecolor="black",
                       label="Reused from connectivity-search"),
        mpatches.Patch(facecolor=NEW, edgecolor="black",
                       label="New for multi-DWPC"),
    ]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.0, 0.02),
              fontsize=7, frameon=False, ncol=3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "query_flow.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    architecture_diagram()
    query_flow_diagram()
    print(f"Wrote: {OUT_DIR / 'architecture.svg'}")
    print(f"Wrote: {OUT_DIR / 'query_flow.svg'}")
