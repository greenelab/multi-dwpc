"""Hurdle + value-respecting adaptive binning for null stratification keys.

Replaces ``_assign_rank_bins`` for the hurdle+adaptive strategies: rank
binning with positional tie-breaking can split equal keys into different
strata (the Entrez-ID-block artifact behind the vacuous Tier 0 promiscuity
result -- see docs/tasks/capacity-hurdle-adaptive-null/design.md). Here
strata are contiguous in key value, equal keys always share a stratum, and
key == 0 (an exact exchangeability class for both first-hop degree and
leave-target-out capacity) is isolated as its own hurdle stratum.
"""

from __future__ import annotations

import numpy as np


def hurdle_adaptive_bins(keys: np.ndarray, min_stratum_size: int = 50) -> np.ndarray:
    """Return a dense bin id (0..M-1) per element, ascending in key value.

    Bin ids ascend with key value. If any key is exactly 0, those elements
    form their own stratum (the lowest bin) regardless of its size. Positive
    keys are grouped greedily in ascending distinct-value order, closing a
    stratum once it holds >= min_stratum_size elements; a trailing stratum
    below the minimum is merged into its predecessor (or stands alone when
    it is the only positive stratum).
    """
    if min_stratum_size <= 0:
        raise ValueError("min_stratum_size must be > 0")

    keys = np.asarray(keys, dtype=float)
    if keys.ndim != 1:
        raise ValueError("keys must be 1-D")
    if (keys < 0).any():
        raise ValueError("keys must be non-negative")

    bins = np.empty(len(keys), dtype=int)
    has_hurdle = bool((keys == 0).any())
    if has_hurdle:
        bins[keys == 0] = 0

    pos_mask = keys > 0
    if not pos_mask.any():
        return bins

    values, counts = np.unique(keys[pos_mask], return_counts=True)
    # Greedy pass over distinct values in ascending order: a stratum closes
    # once it reaches min_stratum_size. Equal keys enter together, so they
    # can never straddle a boundary.
    boundaries = []  # index into `values` where each stratum starts
    acc = 0
    for i, c in enumerate(counts):
        if acc == 0:
            boundaries.append(i)
        acc += c
        if acc >= min_stratum_size:
            acc = 0
    if acc > 0 and len(boundaries) > 1:
        # Trailing under-filled stratum: merge into its predecessor.
        boundaries.pop()

    value_bin = np.empty(len(values), dtype=int)
    first = 1 if has_hurdle else 0
    for b, start in enumerate(boundaries):
        end = boundaries[b + 1] if b + 1 < len(boundaries) else len(values)
        value_bin[start:end] = first + b

    idx = np.searchsorted(values, keys[pos_mask])
    bins[pos_mask] = value_bin[idx]
    return bins


def merge_deficient_strata(
    pools: list[np.ndarray], counts: list[int]
) -> tuple[list[np.ndarray], list[int], list[tuple[int, int]]]:
    """Deterministic feasibility fallback (design doc: binning contract).

    A stratum whose candidate pool (post self-exclusion) is smaller than its
    real-gene count cannot be drawn from; it is merged into its lower-key
    neighbor (next-higher for the lowest stratum), repeatedly, until every
    stratum with count > 0 satisfies count <= pool size. Merge events are
    returned in original stratum indices so callers can log them -- part of
    the null's definition, applied identically to analytical and MC arms.
    A stratum with count == 0 is never deficient (nothing is drawn from it).
    """
    if sum(int(c) for c in counts) > sum(len(p) for p in pools):
        raise ValueError(
            "total real-gene count exceeds total candidate pool; "
            "no partition of these pools is feasible"
        )
    pools = [np.asarray(p) for p in pools]
    counts = [int(c) for c in counts]
    orig_idx = list(range(len(pools)))
    merges: list[tuple[int, int]] = []

    while True:
        deficient = next(
            (i for i, (p, c) in enumerate(zip(pools, counts)) if c > 0 and c > len(p)),
            None,
        )
        if deficient is None:
            return pools, counts, merges
        into = deficient - 1 if deficient > 0 else deficient + 1
        merges.append((orig_idx[deficient], orig_idx[into]))
        lo, hi = sorted((deficient, into))
        pools[lo] = np.concatenate([pools[lo], pools[hi]])
        counts[lo] += counts[hi]
        orig_idx[lo] = orig_idx[into]
        del pools[hi], counts[hi], orig_idx[hi]
