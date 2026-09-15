"""DWPC for a few source genes, computed from the edge files at query time.

The full G->BP DWPC matrices are 33.6 GB; a query needs only its own genes'
rows. hetmatpy computes DWPC as products of degree-weighted adjacency matrices
with diagonal corrections for repeated metanodes. Row r of each product
depends only on row r of its first factor, and removing the diagonal becomes
zeroing entry (r, gene r), so the same computation restricted to the query
genes' rows gives their DWPC. Agreement with hetmatpy's full (dense) matrices
is within ~1e-14 relative; entries below ``DWPC_ZERO_TOL`` are set to 0.

Supported hetmatpy categories are the ones G->BP metapaths use: no_repeats,
short_repeat without a head segment, and BABA without a head segment.
"""

from __future__ import annotations

import threading
from pathlib import Path

import numpy as np
from hetmatpy.degree_weight import _degree_weight, categorize, get_segments, remove_diag
from hetmatpy.hetmat import HetMat as HetMatPy
from hetmatpy.matrix import metaedge_to_adjacency_matrix
from scipy import sparse

from src.dwpc_direct import DEFAULT_DAMPING
from src.summary_null import zero_residues


class UnsupportedMetapathError(ValueError):
    """The metapath is not in the metagraph or uses an unsupported category."""


def _zero_row_diagonal(rows: sparse.csr_matrix, genes: np.ndarray) -> sparse.csr_matrix:
    """``remove_diag`` for the rows of a square matrix taken at ``genes``."""
    r = np.arange(genes.size)
    diagonal = np.asarray(rows[r, genes]).ravel()
    return rows - sparse.csr_matrix((diagonal, (r, genes)), shape=rows.shape)


class QueryDwpc:
    """Row-restricted DWPC from ``data_dir`` edges (hetmat layout)."""

    def __init__(self, data_dir: Path, damping: float = DEFAULT_DAMPING):
        self._graph = HetMatPy(Path(data_dir))
        self._metagraph = self._graph.metagraph
        self.damping = damping
        self._weighted: dict[str, sparse.csr_matrix] = {}
        self._baba_diagonals: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self._lock = threading.Lock()

    def target_values(self, metapath: str, gene_positions: np.ndarray, target_position: int) -> np.ndarray:
        """Raw DWPC from each gene to the target, residues set to 0."""
        genes = np.asarray(gene_positions, dtype=np.int64)
        rows = self.rows(metapath, genes)
        return zero_residues(np.asarray(rows[:, [int(target_position)]].todense()).ravel())

    def rows(self, metapath: str, genes: np.ndarray) -> sparse.csr_matrix:
        """DWPC rows (genes x all targets) for one metapath."""
        try:
            mp = self._metagraph.metapath_from_abbrev(metapath)
        except KeyError as exc:
            raise UnsupportedMetapathError(f"{metapath}: not in the metagraph") from exc
        category = categorize(mp)
        if category == "no_repeats":
            return self._product(mp, genes)
        if category == "short_repeat":
            return self._short_repeat_rows(mp, genes)
        if category == "BABA":
            return self._baba_rows(mp, genes)
        raise UnsupportedMetapathError(f"{metapath}: category {category} is not supported")

    def _weight(self, metaedge) -> sparse.csr_matrix:
        key = str(metaedge)
        with self._lock:
            if key not in self._weighted:
                # dense_threshold=1 keeps the adjacency sparse (0 would densify it).
                _, _, adjacency = metaedge_to_adjacency_matrix(
                    self._graph, metaedge, dtype=np.float64, dense_threshold=1
                )
                weighted = _degree_weight(sparse.csc_matrix(adjacency), self.damping, dtype=np.float64)
                self._weighted[key] = sparse.csr_matrix(weighted)
            return self._weighted[key]

    def _product(self, segment, genes: np.ndarray | None = None) -> sparse.csr_matrix:
        """Degree-weighted walk product of a segment, optionally for rows ``genes``."""
        if categorize(segment) != "no_repeats":
            raise UnsupportedMetapathError(f"{segment}: segment with repeats is not supported")
        product = None
        for metaedge in segment:
            weight = self._weight(metaedge)
            if product is None:
                product = weight if genes is None else weight[genes]
            else:
                product = product @ weight
        return sparse.csr_matrix(product)

    def _short_repeat_rows(self, mp, genes: np.ndarray) -> sparse.csr_matrix:
        head = tail = repeat = None
        for i, segment in enumerate(get_segments(self._metagraph, mp)):
            if segment.source() == segment.target():
                repeat = segment
            elif i == 0:
                head = segment
            else:
                tail = segment
        if head is not None:
            raise UnsupportedMetapathError(f"{mp}: short_repeat with a head segment is not supported")
        node = repeat.source()
        repeat_at = [i for i, v in enumerate(repeat.get_nodes()) if v == node]
        rows = None
        for metaedge in repeat[: repeat_at[1]]:
            rows = self._weight(metaedge)[genes] if rows is None else rows @ self._weight(metaedge)
        rows = _zero_row_diagonal(sparse.csr_matrix(rows), genes)
        if len(repeat_at) == 3:
            inner = None
            for metaedge in repeat[repeat_at[1]:]:
                inner = self._weight(metaedge) if inner is None else inner @ self._weight(metaedge)
            inner = sparse.csr_matrix(remove_diag(inner, dtype=np.float64))
            rows = _zero_row_diagonal(rows @ inner, genes)
        if tail is not None:
            rows = rows @ self._product(tail)
        return sparse.csr_matrix(rows)

    def _baba_rows(self, mp, genes: np.ndarray) -> sparse.csr_matrix:
        segments = get_segments(self._metagraph, mp)
        seg_axb = seg_cda = seg_bed = None
        for i, segment in enumerate(segments[:-2]):
            if segment.source() == segments[i + 2].source() and not seg_axb:
                seg_axb, seg_bya, seg_azb = segment, segments[i + 1], segments[i + 2]
                seg_cda = segments[0] if i == 1 else None
                seg_bed = segments[-1] if segments[-1] != seg_azb else None
        if seg_axb is None or seg_cda is not None:
            raise UnsupportedMetapathError(f"{mp}: BABA with a head segment is not supported")
        axb, bya, azb = self._product(seg_axb), self._product(seg_bya), self._product(seg_azb)
        key = str(mp)
        with self._lock:
            if key not in self._baba_diagonals:
                self._baba_diagonals[key] = ((axb @ bya).diagonal(), (bya @ azb).diagonal())
            diag_ab, diag_ba = self._baba_diagonals[key]
        axb_rows = axb[genes]
        correction_a = sparse.diags(diag_ab[genes]) @ azb[genes]
        correction_b = axb_rows @ sparse.diags(diag_ba)
        correction_c = axb_rows.multiply(sparse.csr_matrix(bya.T)[genes]).multiply(azb[genes])
        rows = axb_rows @ bya @ azb - correction_a - correction_b + correction_c
        if seg_axb.source == seg_azb.target:  # mirrors hetmatpy's (bound-method) check
            rows = _zero_row_diagonal(sparse.csr_matrix(rows), genes)
        if seg_bed is not None:
            rows = rows @ self._product(seg_bed)
        return sparse.csr_matrix(rows)
