"""Shared pool/count assembly for Tier 0 stratified-SRSWOR pool-construction
strategies (promiscuity and metaedge-degree).

Both strategies stratify the full gene universe into `n_bins` bins by some
per-row degree/promiscuity key, then build one candidate pool per bin
(excluding the LV's own real genes from their own candidate pools) and one
count per bin (how many real genes landed in that bin). Extracted here so
that behavior which must stay identical across both strategies -- self-
exclusion logic in particular -- can't silently drift when only one of them
is edited.
"""

from __future__ import annotations

import numpy as np


def pools_from_bins(
    bin_of_row: np.ndarray, real_row_idx: np.ndarray, n_bins: int
) -> tuple[list[np.ndarray], list[int]]:
    """Build per-bin candidate pools and real-gene counts.

    Parameters
    ----------
    bin_of_row : np.ndarray
        Bin id (0..n_bins-1) for every row in the gene universe, in
        gene_ids order.
    real_row_idx : np.ndarray
        Row indices (into the same gene universe) of the LV's own real
        genes for this row -- excluded from every candidate pool.
    n_bins : int
        Number of bins.

    Returns
    -------
    (pools, counts)
        `pools[b]` is the array of candidate row indices in bin `b` with
        `real_row_idx` removed. `counts[b]` is how many real genes fall in
        bin `b`.
    """
    bin_of_row = np.asarray(bin_of_row)
    real_row_idx = np.asarray(real_row_idx)
    real_bins = bin_of_row[real_row_idx]

    # One stable sort groups rows by bin while keeping ascending row order
    # within each bin, so every pool is sliced out of a single array instead
    # of rescanning the whole universe once per bin.
    is_real = np.zeros(bin_of_row.shape[0], dtype=bool)
    is_real[real_row_idx[real_row_idx >= 0]] = True
    candidates = np.flatnonzero(~is_real)
    candidate_bins = bin_of_row[candidates]
    order = np.argsort(candidate_bins, kind="stable")
    candidates = candidates[order]
    bin_edges = np.searchsorted(candidate_bins[order], np.arange(n_bins + 1), side="left")
    pools = [candidates[bin_edges[b]:bin_edges[b + 1]] for b in range(n_bins)]

    real_edges = np.searchsorted(np.sort(real_bins), np.arange(n_bins + 1), side="left")
    counts = np.diff(real_edges).tolist()

    return pools, counts
