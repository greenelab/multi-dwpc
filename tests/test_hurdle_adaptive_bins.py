import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest

from src.hurdle_adaptive_bins import hurdle_adaptive_bins


def test_zero_keys_form_exclusive_hurdle_stratum():
    keys = np.array([0, 0, 0, 1, 1, 2, 3, 5, 8, 8])
    bins = hurdle_adaptive_bins(keys, min_stratum_size=3)
    assert set(bins[keys == 0]) == {0}
    assert 0 not in set(bins[keys > 0])


def test_equal_keys_never_split():
    rng = np.random.default_rng(0)
    keys = rng.integers(0, 6, size=200)
    bins = hurdle_adaptive_bins(keys, min_stratum_size=7)
    for v in np.unique(keys):
        assert len(set(bins[keys == v])) == 1, f"key {v} split across strata"


def test_min_stratum_size_met_except_none():
    rng = np.random.default_rng(1)
    keys = np.concatenate([np.zeros(50), rng.integers(1, 40, size=300)])
    bins = hurdle_adaptive_bins(keys, min_stratum_size=25)
    sizes = np.bincount(bins)
    assert (sizes[1:] >= 25).all(), f"non-hurdle stratum below min size: {sizes}"


def test_bins_ascend_with_key_value():
    keys = np.array([0.0, 0.5, 0.5, 1.2, 1.2, 3.0, 3.0, 9.9, 9.9, 9.9])
    bins = hurdle_adaptive_bins(keys, min_stratum_size=2)
    order = np.argsort(keys)
    assert (np.diff(bins[order]) >= 0).all()


def test_all_zero_keys_single_stratum():
    bins = hurdle_adaptive_bins(np.zeros(10), min_stratum_size=5)
    assert set(bins) == {0}


def test_no_zero_keys_no_empty_hurdle():
    keys = np.array([1.0, 1.0, 2.0, 2.0, 3.0, 3.0])
    bins = hurdle_adaptive_bins(keys, min_stratum_size=2)
    sizes = np.bincount(bins)
    assert (sizes > 0).all()


def test_continuous_keys_work():
    rng = np.random.default_rng(2)
    keys = np.concatenate([np.zeros(20), rng.gamma(2.0, 3.0, size=500)])
    bins = hurdle_adaptive_bins(keys, min_stratum_size=50)
    sizes = np.bincount(bins)
    assert sizes[0] == 20
    assert (sizes[1:] >= 50).all()


def test_negative_keys_rejected():
    with pytest.raises(ValueError):
        hurdle_adaptive_bins(np.array([-1.0, 0.0, 1.0]), min_stratum_size=2)


def test_non_positive_min_stratum_size_rejected():
    # Hardening deviation (Task 1): min_stratum_size <= 0 raises ValueError
    # rather than silently producing degenerate binning.
    with pytest.raises(ValueError):
        hurdle_adaptive_bins(np.array([0.0, 1.0, 2.0]), min_stratum_size=0)
    with pytest.raises(ValueError):
        hurdle_adaptive_bins(np.array([0.0, 1.0, 2.0]), min_stratum_size=-5)


from src.hurdle_adaptive_bins import merge_deficient_strata


def _pools(*sizes):
    start = 0
    out = []
    for s in sizes:
        out.append(np.arange(start, start + s))
        start += s
    return out


def test_no_merge_when_feasible():
    pools, counts, merges = merge_deficient_strata(_pools(5, 5, 5), [2, 3, 1])
    assert merges == []
    assert [len(p) for p in pools] == [5, 5, 5]
    assert counts == [2, 3, 1]


def test_deficient_stratum_merges_into_lower_neighbor():
    # stratum 2 needs 4 but has 2 candidates
    pools, counts, merges = merge_deficient_strata(_pools(5, 5, 2), [1, 1, 4])
    assert merges == [(2, 1)]
    assert len(pools) == 2
    assert counts == [1, 5]
    assert len(pools[1]) == 7


def test_lowest_stratum_merges_upward():
    pools, counts, merges = merge_deficient_strata(_pools(1, 6), [3, 1])
    assert merges == [(0, 1)]
    assert counts == [4]
    assert len(pools[0]) == 7


def test_cascading_merge():
    # After (2 -> 1), stratum 1 holds pool 3+2=5 vs count 2+4=6: still short,
    # so it merges into 0.
    pools, counts, merges = merge_deficient_strata(_pools(9, 3, 2), [1, 2, 4])
    assert merges == [(2, 1), (1, 0)]
    assert counts == [7]
    assert len(pools[0]) == 14


def test_globally_infeasible_raises():
    with pytest.raises(ValueError):
        merge_deficient_strata(_pools(2, 2), [3, 3])


def test_zero_count_strata_untouched():
    pools, counts, merges = merge_deficient_strata(_pools(4, 0, 4), [2, 0, 2])
    assert merges == []
    assert counts == [2, 0, 2]
