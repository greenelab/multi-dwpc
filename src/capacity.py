"""Leave-target-out metapath capacity: the stratification key for the
primary hurdle+adaptive strategy (design doc, strategy S1).

Capacity is computed on the RAW DWPC scale -- the per-entry arcsinh
transform used for observed scores is nonlinear, so summing transformed
entries would not be a walk-capacity quantity. Excluding the feature's own
target column keeps the key a generic "reach along this metapath"
covariate rather than the outcome; capacity 0 is an exact exchangeability
class (a zero row means zero DWPC to every target of the metapath).
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from src.dwpc_direct import DEFAULT_DAMPING
from src.sparse_column import dense_column


def leave_target_out_capacity(
    matrix: sparse.spmatrix,
    target_position: int,
    target_col: np.ndarray | None = None,
    row_sums: np.ndarray | None = None,
) -> np.ndarray:
    """Row sums of the raw DWPC matrix minus the target column's entries.

    Accepts any scipy sparse format (CSR from the web app's cache, CSC from
    ``CapacityProvider``). Pass ``target_col`` (the dense raw target column)
    when the caller already has it, to avoid reading the column twice, and
    ``row_sums`` (e.g. from ``HetMat.get_dwpc_row_sums``) to avoid rescanning
    the matrix.
    """
    if row_sums is None:
        row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    if target_col is None:
        target_col = dense_column(matrix, target_position)
    return row_sums - target_col


class CapacityProvider:
    """Per-metapath capacity vectors from a HetMat's cached DWPC matrices.

    The row-sum vector is cached per metapath (the expensive part is the
    matrix load); the per-feature target-column subtraction is done per
    call. ``hetmat`` must expose ``compute_dwpc_matrix_csc(metapath,
    damping)`` returning raw DWPC with gene rows -- the same source
    ``precompute_gene_feature_scores`` builds the observed scores from, so
    key and statistic describe the same matrix.
    """

    def __init__(self, hetmat, damping: float = DEFAULT_DAMPING):
        self._hetmat = hetmat
        self._damping = damping
        self._rowsums: dict[str, np.ndarray] = {}

    def capacity(self, metapath: str, target_position: int) -> np.ndarray:
        matrix = self._hetmat.compute_dwpc_matrix_csc(metapath, damping=self._damping)
        if metapath not in self._rowsums:
            self._rowsums[metapath] = np.asarray(matrix.sum(axis=1)).ravel()
        target_col = np.asarray(matrix[:, target_position].todense()).ravel()
        return self._rowsums[metapath] - target_col
