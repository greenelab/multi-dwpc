"""Thin project-owned wrapper around the vendored exact_resampling_moments
kernel. Derives z/p locally via the normal approximation -- never surfaces
HetNetEX-MD's own p_edgeworth/p_normal/exact_median_pvalue, per the
validation spec's tail-calibration finding (those are anti-conservative,
1.21x-13.5x excess in the tail).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.stats import norm

from src.hetnetex_md_import import exact_resampling_moments


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
    p = float(norm.sf(z))
    return AnalyticalNullResult(
        mean=mean, var=var, std=std, z=float(z), p=p,
        n_pool=n_pool, k_total=k_total,
    )
