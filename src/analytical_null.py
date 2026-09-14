"""Thin project-owned wrapper around the pinned `exact_resampling_moments`
kernel (imported via :mod:`src.hetnetex_md_import`). Derives z/p locally via
the normal approximation -- never surfaces HetNetEX-MD's own
p_edgeworth/p_normal/exact_median_pvalue, per the validation spec's
tail-calibration finding (those are anti-conservative,
1.21x-13.5x excess in the tail).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from src.capacity import leave_target_out_capacity
from src.dwpc_direct import DEFAULT_DAMPING, get_dwpc_raw_mean, transform_dwpc
from src.hetnetex_md_import import exact_resampling_moments
from src.hurdle_adaptive_bins import hurdle_adaptive_bins, merge_deficient_strata
from src.pool_assembly import pools_from_bins
from src.sparse_column import dense_column
from src.summary_null import GeneSetZResult


@dataclass
class AnalyticalNullResult:
    mean: float
    var: float
    std: float
    z: float
    p: float
    n_pool: int
    k_total: int


def analytical_null(
    scores: np.ndarray,
    pools: Sequence[np.ndarray],
    counts: Sequence[int],
    observed: float,
) -> AnalyticalNullResult:
    scores = np.asarray(scores, dtype=float)
    if not np.all(np.isfinite(scores)):
        raise ValueError("scores must be finite")

    k_total = int(sum(counts))
    if k_total <= 0:
        raise ValueError("K (sum of counts) must be > 0")

    n_pool = 0
    for pool, k in zip(pools, counts):
        n_pool += len(pool)
        if k > len(pool):
            raise ValueError(
                f"stratum count {k} exceeds pool size {len(pool)}"
            )

    mean, var, _mu3 = exact_resampling_moments(scores, pools, counts)

    if not np.isfinite(var) or var <= 0.0:
        # Zero-variance null: undefined, not +/-inf. Maps to NaN so callers
        # must handle it explicitly rather than silently treating it as an
        # extreme-but-finite result.
        return AnalyticalNullResult(
            mean=mean, var=var, std=float("nan"),
            z=float("nan"), p=float("nan"),
            n_pool=n_pool, k_total=k_total,
        )

    std = float(np.sqrt(var))
    z = (observed - mean) / std
    # Normal survival function via math.erfc (equal to scipy's norm.sf to
    # ~1e-13 relative, tail included) so this module, loaded at web-app
    # startup, does not pull in scipy.stats (~50 MB resident).
    p = 0.5 * math.erfc(z / math.sqrt(2.0))
    return AnalyticalNullResult(
        mean=mean, var=var, std=std, z=float(z), p=p,
        n_pool=n_pool, k_total=k_total,
    )


def analytical_gene_set_z(
    hetmat,
    metapath: str,
    source_idx: np.ndarray,
    target_pos: int,
    *,
    min_stratum_size: int = 50,
) -> GeneSetZResult:
    """Query adapter: gene-set enrichment z-score against the analytical,
    capacity-stratified null for one metapath/target pair (design doc,
    "What is computed").

    The statistic stays on the transformed scale (arcsinh(raw / raw_mean),
    mean over the mapped gene set); the stratification key -- leave-target-
    out capacity -- stays raw. One matrix load serves both: this function
    calls ``hetmat.compute_dwpc_matrix`` exactly once and applies the
    pure ``leave_target_out_capacity`` to that same matrix directly, rather
    than routing the capacity through ``CapacityProvider`` (which would
    issue its own call into the hetmat and depend on the hetmat's internal
    cache to avoid a second load).

    The matrix is requested in CSR, the format the web app preloads and
    ``get_dwpc_for_pairs`` reads. ``compute_dwpc_matrix_csc`` would keep a
    second, CSC copy of every queried metapath in the hetmat's cache,
    doubling resident DWPC memory for a column slice that CSR also serves.
    The target column is read with :func:`src.sparse_column.dense_column`
    (per-row binary search) rather than a CSR column slice, which would scan
    every stored entry. The target-independent row sums come from the
    hetmat's per-metapath cache (``get_dwpc_row_sums``) for the same reason.
    """
    source_idx = np.asarray(source_idx, dtype=np.int64)

    matrix = hetmat.compute_dwpc_matrix(metapath, damping=DEFAULT_DAMPING)

    raw_target_col = dense_column(matrix, target_pos)
    raw_mean = get_dwpc_raw_mean(hetmat.metapath_stats, metapath)
    scores = transform_dwpc(raw_target_col, raw_mean)

    capacity = leave_target_out_capacity(
        matrix,
        target_pos,
        target_col=raw_target_col,
        row_sums=hetmat.get_dwpc_row_sums(metapath, damping=DEFAULT_DAMPING),
    )

    bins = hurdle_adaptive_bins(capacity, min_stratum_size=min_stratum_size)
    n_bins = int(bins.max()) + 1
    pools, counts = pools_from_bins(bins, source_idx, n_bins)
    pools, counts, merges = merge_deficient_strata(pools, counts)

    observed = float(scores[source_idx].mean())
    null_result = analytical_null(scores, pools, counts, observed)

    n_active_strata = sum(1 for c in counts if c > 0)

    return GeneSetZResult(
        real_mean=observed,
        null_mean=null_result.mean,
        null_std=null_result.std,
        z=null_result.z,
        p_value=null_result.p,
        n_active_strata=n_active_strata,
        merges=merges,
    )
