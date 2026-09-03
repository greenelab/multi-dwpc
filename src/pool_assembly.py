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
    real_bins = bin_of_row[real_row_idx]

    pools: list[np.ndarray] = []
    counts: list[int] = []
    for b in range(n_bins):
        candidate_rows = np.flatnonzero(bin_of_row == b)
        candidate_rows = candidate_rows[~np.isin(candidate_rows, real_row_idx)]
        pools.append(candidate_rows)
        counts.append(int((real_bins == b).sum()))

    return pools, counts
