import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest

from src.pool_assembly import pools_from_bins


def _reference_pools_from_bins(bin_of_row, real_row_idx, n_bins):
    """The original per-bin loop, kept as the behavioral reference."""
    real_bins = bin_of_row[real_row_idx]
    pools, counts = [], []
    for b in range(n_bins):
        candidate_rows = np.flatnonzero(bin_of_row == b)
        candidate_rows = candidate_rows[~np.isin(candidate_rows, real_row_idx)]
        pools.append(candidate_rows)
        counts.append(int((real_bins == b).sum()))
    return pools, counts


def _assert_same(got, want):
    got_pools, got_counts = got
    want_pools, want_counts = want
    assert got_counts == want_counts
    assert all(type(c) is int for c in got_counts)
    assert len(got_pools) == len(want_pools)
    for g, w in zip(got_pools, want_pools):
        assert g.dtype == w.dtype
        np.testing.assert_array_equal(g, w)


@pytest.mark.parametrize("seed", range(20))
def test_matches_reference_on_random_partitions(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(1, 3000))
    n_bins = int(rng.integers(1, 300))
    # Some bins left empty, some rows binned beyond n_bins (excluded from all pools).
    bin_of_row = rng.integers(0, n_bins + 3, size=n)
    k = int(rng.integers(0, min(n, 60) + 1))
    real_row_idx = rng.choice(n, size=k, replace=False)
    _assert_same(
        pools_from_bins(bin_of_row, real_row_idx, n_bins),
        _reference_pools_from_bins(bin_of_row, real_row_idx, n_bins),
    )


def test_matches_reference_with_duplicate_and_negative_real_indices():
    bin_of_row = np.array([0, 1, 1, 2, 0, 2, 1])
    real_row_idx = np.array([1, 1, -1, 4])
    _assert_same(
        pools_from_bins(bin_of_row, real_row_idx, 3),
        _reference_pools_from_bins(bin_of_row, real_row_idx, 3),
    )


def test_real_genes_excluded_and_counted():
    bin_of_row = np.array([0, 0, 1, 1, 1])
    pools, counts = pools_from_bins(bin_of_row, np.array([1, 3]), 2)
    np.testing.assert_array_equal(pools[0], [0])
    np.testing.assert_array_equal(pools[1], [2, 4])
    assert counts == [1, 1]
