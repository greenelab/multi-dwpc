"""Capacity-stratified null computed from per-stratum summaries.

The exact null (``src.analytical_null.analytical_gene_set_z``) needs the target
column of the DWPC matrix for every gene. Its mean and variance depend on each
stratum's pool only through the pool's size, mean and centered sum of squares,
so those can be precomputed once per (metapath, target) on the full gene
universe (``build_strata``) and adjusted at query time by removing the query
genes (``summary_gene_set_z``). No matrix is needed at query time.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from src.hurdle_adaptive_bins import hurdle_adaptive_bins

# DWPC entries at or below this magnitude are floating-point cancellation
# residues of true zeros (~1e-18, some negative); real values are ~1e-5.
DWPC_ZERO_TOL = 1e-15

# A query gene's leave-target-out capacity within this fraction of its row sum
# is zero: rounding in the on-the-fly DWPC must not move a capacity-0 gene out
# of the hurdle stratum.
CAPACITY_REL_TOL = 1e-10

_EPS = np.finfo(float).eps


@dataclass
class GeneSetZResult:
    real_mean: float
    null_mean: float
    null_std: float
    z: float
    p_value: float
    n_active_strata: int
    merges: list[tuple[int, int]]


@dataclass
class StratumTable:
    """Per-stratum summaries for one (metapath, target), ascending capacity."""

    capacity_min: np.ndarray
    capacity_max: np.ndarray
    n_genes: np.ndarray
    score_mean: np.ndarray
    score_centered_ss: np.ndarray


def zero_residues(values: np.ndarray) -> np.ndarray:
    """Return a copy with entries of magnitude <= ``DWPC_ZERO_TOL`` set to 0."""
    values = np.array(values, dtype=float)
    values[np.abs(values) <= DWPC_ZERO_TOL] = 0.0
    return values


def query_capacity(row_sums: np.ndarray, target_values: np.ndarray) -> np.ndarray:
    """Leave-target-out capacity of query genes, with rounding residue set to 0."""
    capacity = np.asarray(row_sums, dtype=float) - np.asarray(target_values, dtype=float)
    capacity[np.abs(capacity) <= CAPACITY_REL_TOL * np.abs(row_sums)] = 0.0
    return capacity


def build_strata(capacity: np.ndarray, scores: np.ndarray, min_stratum_size: int = 50) -> StratumTable:
    """Summarise every stratum of one target column over all genes."""
    capacity = np.asarray(capacity, dtype=float)
    scores = np.asarray(scores, dtype=float)
    bins = hurdle_adaptive_bins(capacity, min_stratum_size=min_stratum_size)
    n_strata = int(bins.max()) + 1
    n_genes = np.bincount(bins, minlength=n_strata)
    score_mean = np.bincount(bins, weights=scores, minlength=n_strata) / n_genes
    score_centered_ss = np.bincount(bins, weights=(scores - score_mean[bins]) ** 2, minlength=n_strata)
    capacity_min = np.full(n_strata, np.inf)
    capacity_max = np.full(n_strata, -np.inf)
    np.minimum.at(capacity_min, bins, capacity)
    np.maximum.at(capacity_max, bins, capacity)
    return StratumTable(capacity_min, capacity_max, n_genes, score_mean, score_centered_ss)


def assign_strata(table: StratumTable, capacity: np.ndarray) -> np.ndarray:
    """Stratum index of each capacity value.

    Capacity 0 belongs to the hurdle stratum (index 0) when the table has one.
    A positive capacity belongs to the stratum whose range contains it; when
    rounding lands it between two ranges, to the nearer one.
    """
    capacity = np.asarray(capacity, dtype=float)
    has_hurdle = bool(table.capacity_max[0] == 0.0)
    first = 1 if has_hurdle else 0
    lo = table.capacity_min[first:]
    hi = table.capacity_max[first:]
    if lo.size == 0:
        return np.zeros(capacity.size, dtype=int)
    idx = np.clip(np.searchsorted(lo, capacity, side="right") - 1, 0, lo.size - 1)
    nxt = np.minimum(idx + 1, lo.size - 1)
    in_gap_nearer_next = (capacity > hi[idx]) & ((lo[nxt] - capacity) < (capacity - hi[idx]))
    strata = np.where(in_gap_nearer_next, nxt, idx) + first
    if has_hurdle:
        strata = np.where(capacity == 0.0, 0, strata)
    return strata


def _zero_small_ss(ss: np.ndarray, scale: np.ndarray, n: np.ndarray, mean: np.ndarray) -> np.ndarray:
    """Treat a pool's centered sum of squares as 0 when it is within rounding error.

    The stratum sum of squares accumulates ~eps per summed gene, and removing
    the query genes cancels terms of size ``scale``, so the bound grows with the
    pool size ``n``; the second term covers rounding of the pool mean.
    """
    tol = _EPS * (64.0 + n) * scale + n * (64.0 * _EPS * mean) ** 2
    return np.where(ss <= tol, 0.0, ss)


def summary_gene_set_z(table: StratumTable, scores: np.ndarray, capacities: np.ndarray) -> GeneSetZResult:
    """z of the query genes' mean score against the stratified SRSWOR null.

    ``scores`` and ``capacities`` describe the distinct query genes. Matches
    ``analytical_gene_set_z`` on the full matrix: pools exclude the query
    genes, deficient strata merge as in ``merge_deficient_strata``, and the
    null mean and variance are the exact SRSWOR moments.
    """
    scores = np.asarray(scores, dtype=float)
    strata = assign_strata(table, capacities)
    n_strata = table.n_genes.size

    counts = np.bincount(strata, minlength=n_strata)
    d = scores - table.score_mean[strata]
    removed_d = np.bincount(strata, weights=d, minlength=n_strata)
    removed_sq = np.bincount(strata, weights=d ** 2, minlength=n_strata)
    n = (table.n_genes - counts).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(n > 0, table.score_mean - removed_d / n, 0.0)
        ss = np.where(n > 0, table.score_centered_ss - removed_sq - removed_d ** 2 / n, 0.0)
    scale = table.score_centered_ss + removed_sq
    ss = _zero_small_ss(np.maximum(ss, 0.0), scale, n, mean)

    if counts.sum() > n.sum():
        raise ValueError(
            "total real-gene count exceeds total candidate pool; "
            "no partition of these pools is feasible"
        )

    # Deficient-stratum merge, same order and bookkeeping as merge_deficient_strata.
    pools = [list(p) for p in zip(n, mean, ss, scale)]
    counts = counts.tolist()
    orig_idx = list(range(n_strata))
    merges: list[tuple[int, int]] = []
    while True:
        deficient = next((i for i, (p, c) in enumerate(zip(pools, counts)) if c > 0 and c > p[0]), None)
        if deficient is None:
            break
        into = deficient - 1 if deficient > 0 else deficient + 1
        merges.append((orig_idx[deficient], orig_idx[into]))
        lo, hi = sorted((deficient, into))
        (n1, m1, ss1, sc1), (n2, m2, ss2, sc2) = pools[lo], pools[hi]
        n12 = n1 + n2
        m12 = (n1 * m1 + n2 * m2) / n12 if n12 else 0.0
        between = n1 * (m1 - m12) ** 2 + n2 * (m2 - m12) ** 2
        pools[lo] = [n12, m12, ss1 + ss2 + between, sc1 + sc2 + between]
        counts[lo] += counts[hi]
        orig_idx[lo] = orig_idx[into]
        del pools[hi], counts[hi], orig_idx[hi]

    n, mean, ss, scale = (np.array(col, dtype=float) for col in zip(*pools))
    ss = _zero_small_ss(ss, scale, n, mean)
    k = np.array(counts, dtype=float)
    active = k > 0
    n, mean, ss, k = n[active], mean[active], ss[active], k[active]
    weight = k / k.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        var_s = np.where(k >= n, 0.0, (ss / n / k) * (n / (n - 1.0)) * (1.0 - k / n))
    null_mean = float(np.sum(weight * mean))
    null_var = float(np.sum(weight ** 2 * var_s))

    real_mean = float(scores.mean())
    if not np.isfinite(null_var) or null_var <= 0.0:
        z = p = std = float("nan")
    else:
        std = math.sqrt(null_var)
        z = (real_mean - null_mean) / std
        p = 0.5 * math.erfc(z / math.sqrt(2.0))
    return GeneSetZResult(
        real_mean=real_mean,
        null_mean=null_mean,
        null_std=std,
        z=float(z),
        p_value=float(p),
        n_active_strata=int(active.sum()),
        merges=merges,
    )
