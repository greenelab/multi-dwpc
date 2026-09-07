import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest

from src.hetnetex_md_import import exact_resampling_moments

TOL = 1e-12


@pytest.mark.parametrize("N,k", [(9, 3), (10, 4), (12, 5), (11, 2), (8, 4), (7, 1)])
def test_srswor_moments_match_exhaustive_enumeration(N, k):
    """Enumerate all C(N,k) samples; compare mean, variance, third moment.

    Adapted from HetNetEX-MD's own
    tests/test_exactness.py::test_srswor_moments_match_exhaustive_enumeration,
    run here against the actual installed package (not a copy) to confirm it
    behaves correctly against our own numpy/scipy versions.
    """
    rng = np.random.default_rng(N * 100 + k)
    x = rng.gamma(2.0, 1.0, size=N)
    pool = np.arange(N)

    means = np.array([x[list(c)].mean() for c in itertools.combinations(range(N), k)])
    brute_mean = means.mean()
    brute_var = ((means - brute_mean) ** 2).mean()
    brute_mu3 = ((means - brute_mean) ** 3).mean()

    mean, var, mu3 = exact_resampling_moments(x, [pool], [k])

    assert abs(mean - brute_mean) < TOL
    assert abs(var - brute_var) < TOL * max(1.0, abs(brute_var))
    assert abs(mu3 - brute_mu3) < 1e-10 * max(1.0, abs(brute_mu3))


@pytest.mark.parametrize(
    "N1,k1,N2,k2", [(7, 3, 6, 2), (8, 3, 7, 4), (6, 2, 6, 3), (5, 2, 5, 2)]
)
def test_srswor_moments_stratified(N1, k1, N2, k2):
    """Two independent strata: means add, variances and third cumulants add."""
    N = N1 + N2
    rng = np.random.default_rng(N1 * 1000 + k1 * 100 + N2 * 10 + k2)
    x = rng.gamma(2.0, 1.0, size=N)
    p1, p2 = np.arange(N1), np.arange(N1, N)

    vals = [
        x[list(c1) + list(c2)].mean()
        for c1 in itertools.combinations(p1, k1)
        for c2 in itertools.combinations(p2, k2)
    ]
    vals = np.asarray(vals)
    bm = vals.mean()

    mean, var, mu3 = exact_resampling_moments(x, [p1, p2], [k1, k2])
    assert abs(mean - bm) < TOL
    assert abs(var - ((vals - bm) ** 2).mean()) < 1e-11
    assert abs(mu3 - ((vals - bm) ** 3).mean()) < 1e-11


def test_zero_count_stratum_contributes_nothing():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    mean, var, mu3 = exact_resampling_moments(x, [np.arange(4), np.arange(4)], [2, 0])
    mean2, var2, mu32 = exact_resampling_moments(x, [np.arange(4)], [2])
    assert mean == pytest.approx(mean2)
    assert var == pytest.approx(var2)
    assert mu3 == pytest.approx(mu32)


def test_forbidden_functions_are_not_reexported():
    """Guards against accidentally widening the import surface later.

    Hardening deviation (Task 1, item 3): the import-policy module exposes
    only exact_resampling_moments, never edgeworth_upper_tail,
    exact_median_pvalue, aggregate_network_null, or network_null_moments.
    """
    import src.hetnetex_md_import as wrapper

    assert hasattr(wrapper, "exact_resampling_moments")
    for name in ("edgeworth_upper_tail", "exact_median_pvalue", "aggregate_network_null", "network_null_moments"):
        assert not hasattr(wrapper, name)
