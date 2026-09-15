"""Dense column extraction from a CSR matrix without scanning every entry.

``matrix[:, j]`` on a CSR matrix walks all stored entries (67 ms on the
largest G->BP DWPC matrix, 205M nonzeros). The web query needs one target
column per metapath, and keeping a CSC copy to make that cheap doubles
resident memory. For CSR with sorted, duplicate-free indices, column ``j``
can instead be read with one binary search per row, run for all rows at once
(~2 ms on the same matrix, independent of the number of nonzeros).
"""

from __future__ import annotations

import numpy as np
from scipy import sparse


def dense_column(matrix: sparse.spmatrix, position: int) -> np.ndarray:
    """Return column ``position`` of ``matrix`` as a dense 1-D array.

    Values and dtype are identical to
    ``np.asarray(matrix[:, position].todense()).ravel()``. Canonical CSR input
    takes the per-row binary-search path; any other format, or CSR with
    unsorted or duplicate indices, falls back to that slice.
    """
    n_rows, n_cols = matrix.shape
    if not -n_cols <= position < n_cols:
        raise IndexError(f"column index {position} out of range for {n_cols} columns")
    if position < 0:
        position += n_cols
    if not (sparse.isspmatrix_csr(matrix) and matrix.has_canonical_format):
        return np.asarray(matrix[:, position].todense()).ravel()

    column = np.zeros(n_rows, dtype=matrix.dtype)
    if matrix.nnz == 0:
        return column

    indices = matrix.indices
    lo = matrix.indptr[:-1].astype(np.int64)
    hi = matrix.indptr[1:].astype(np.int64)
    row_end = hi.copy()
    # Lower-bound bisection of every row's sorted column indices at once:
    # afterwards lo[r] is the first entry of row r with column >= position.
    active = lo < hi
    while active.any():
        mid = (lo + hi) // 2
        go_right = active & (indices[np.minimum(mid, matrix.nnz - 1)] < position)
        lo = np.where(go_right, mid + 1, lo)
        hi = np.where(active & ~go_right, mid, hi)
        active = lo < hi

    found = lo < row_end
    found[found] = indices[lo[found]] == position
    column[found] = matrix.data[lo[found]]
    return column
