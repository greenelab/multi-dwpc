import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.analytical_null import analytical_gene_set_z
from src.dwpc_direct import DEFAULT_DAMPING, transform_dwpc
from src.hurdle_adaptive_bins import hurdle_adaptive_bins
from src.summary_null import (
    assign_strata,
    build_strata,
    query_capacity,
    summary_gene_set_z,
    zero_residues,
)

METAPATH = "GaDlA"
RAW_MEAN = 3.0


class _StubHetMat:
    """Full-matrix reference: column 0 is the target, the rest set capacity."""

    def __init__(self, dense):
        self._matrix = sparse.csr_matrix(dense)
        self.metapath_stats = pd.DataFrame({"metapath": [METAPATH], "dwpc_raw_mean": [RAW_MEAN]})

    def compute_dwpc_matrix(self, metapath, damping=DEFAULT_DAMPING):
        return self._matrix

    def get_dwpc_row_sums(self, metapath, damping=DEFAULT_DAMPING):
        return np.asarray(self._matrix.sum(axis=1)).ravel()


def _random_matrix(rng, n_genes, zero_frac=0.2, tie_frac=0.3, n_other=3):
    other = rng.exponential(scale=1.0, size=(n_genes, n_other))
    other[rng.random(n_genes) < zero_frac] = 0.0
    target = rng.exponential(scale=1.0, size=n_genes)
    target[rng.random(n_genes) < 0.5] = 0.0
    dense = np.column_stack([target, other])
    # Identical rows give exactly tied capacities, as structurally identical genes do.
    dense[rng.random(n_genes) < tie_frac] = dense[0]
    return dense


def _summary_z(dense, genes, min_stratum_size):
    target = dense[:, 0]
    row_sums = _StubHetMat(dense).get_dwpc_row_sums(METAPATH)
    scores = transform_dwpc(target, RAW_MEAN)
    table = build_strata(row_sums - target, scores, min_stratum_size)
    return summary_gene_set_z(table, scores[genes], query_capacity(row_sums[genes], target[genes]))


def _assert_matches(summary, reference):
    assert math.isnan(summary.z) == math.isnan(reference.z)
    for field in ("real_mean", "null_mean", "null_std", "z", "p_value"):
        got, want = getattr(summary, field), getattr(reference, field)
        if math.isnan(want):
            assert math.isnan(got), field
        else:
            assert got == pytest.approx(want, rel=1e-12, abs=1e-300), field
    assert summary.merges == reference.merges
    assert summary.n_active_strata == reference.n_active_strata


@pytest.mark.parametrize("seed", range(25))
def test_matches_exact_adapter_on_random_columns(seed):
    rng = np.random.default_rng(seed)
    n_genes = int(rng.integers(40, 600))
    min_stratum_size = int(rng.choice([5, 20, 50]))
    dense = _random_matrix(rng, n_genes)
    genes = rng.choice(n_genes, size=int(rng.integers(1, min(n_genes // 3, 80))), replace=False)
    reference = analytical_gene_set_z(_StubHetMat(dense), METAPATH, genes, 0, min_stratum_size=min_stratum_size)
    _assert_matches(_summary_z(dense, genes, min_stratum_size), reference)


def test_deficient_strata_merge_like_the_exact_adapter():
    # Most of the low-capacity stratum is in the query, forcing merges.
    n = 12
    dense = np.column_stack([np.arange(1.0, n + 1.0), [1.0] * 6 + [2.0] * 6])
    genes = np.array([0, 1, 2, 3, 4])
    reference = analytical_gene_set_z(_StubHetMat(dense), METAPATH, genes, 0, min_stratum_size=1)
    summary = _summary_z(dense, genes, min_stratum_size=1)
    assert reference.merges
    _assert_matches(summary, reference)


def test_stratum_fully_drawn_has_zero_variance_contribution():
    # k == N in a stratum: SRSWOR of the whole pool has no variance.
    dense = np.column_stack([np.arange(1.0, 9.0), [1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]])
    genes = np.array([0, 2])
    reference = analytical_gene_set_z(_StubHetMat(dense), METAPATH, genes, 0, min_stratum_size=2)
    _assert_matches(_summary_z(dense, genes, min_stratum_size=2), reference)


def test_assignment_reproduces_hurdle_adaptive_bins():
    rng = np.random.default_rng(7)
    dense = _random_matrix(rng, 2000)
    capacity = _StubHetMat(dense).get_dwpc_row_sums(METAPATH) - dense[:, 0]
    table = build_strata(capacity, dense[:, 0], min_stratum_size=50)
    np.testing.assert_array_equal(assign_strata(table, capacity), hurdle_adaptive_bins(capacity, 50))


def test_rounding_near_boundaries_and_zero_keeps_the_stratum():
    rng = np.random.default_rng(8)
    dense = _random_matrix(rng, 2000)
    row_sums = _StubHetMat(dense).get_dwpc_row_sums(METAPATH)
    capacity = row_sums - dense[:, 0]
    table = build_strata(capacity, dense[:, 0], min_stratum_size=50)
    jitter = 1 + rng.choice([-1, 1], size=capacity.size) * 1e-13
    noisy = query_capacity(row_sums, dense[:, 0] * jitter)
    np.testing.assert_array_equal(assign_strata(table, noisy), hurdle_adaptive_bins(capacity, 50))


def test_planted_enrichment_is_detected():
    rng = np.random.default_rng(42)
    n = 400
    dense = np.column_stack([rng.exponential(1.0, n), rng.exponential(2.0, n) + 0.1])
    genes = rng.choice(n, size=20, replace=False)
    dense[genes, 0] += 50.0
    assert _summary_z(dense, genes, min_stratum_size=20).z > 1.65


def test_stratum_matched_null_draws_are_standard_normal():
    rng = np.random.default_rng(7)
    n = 500
    dense = np.column_stack([rng.exponential(1.0, n), rng.exponential(2.0, n) + 0.1])
    capacity = dense[:, 1]
    bins = hurdle_adaptive_bins(capacity, 20)
    template = rng.choice(n, size=15, replace=False)
    z = []
    for _ in range(300):
        # Same stratum composition as the template, genes drawn at random.
        genes = np.concatenate(
            [rng.choice(np.flatnonzero(bins == b), size=int(c), replace=False)
             for b, c in zip(*np.unique(bins[template], return_counts=True))]
        )
        z.append(_summary_z(dense, genes, min_stratum_size=20).z)
    z = np.array(z)
    assert abs(z.mean()) < 0.15
    assert abs(z.std() - 1.0) < 0.15


def test_zero_variance_pool_is_nan_but_tiny_spread_is_finite():
    n = 40
    genes = np.array([0, 10, 20, 30])  # one gene per stratum, pools keep 9 genes
    flat = np.column_stack([np.full(n, 2.0), np.linspace(0.1, 5.0, n)])
    reference = analytical_gene_set_z(_StubHetMat(flat), METAPATH, genes, 0, min_stratum_size=10)
    assert math.isnan(reference.z)
    assert math.isnan(_summary_z(flat, genes, min_stratum_size=10).z)
    spread = flat.copy()
    spread[:, 0] *= 1 + 1e-6 * np.random.default_rng(0).standard_normal(n)
    assert math.isfinite(_summary_z(spread, genes, min_stratum_size=10).z)


def test_zero_residues():
    np.testing.assert_array_equal(zero_residues(np.array([1e-18, -3e-19, 1e-15, 2e-15, 1e-5])),
                                  [0.0, 0.0, 0.0, 2e-15, 1e-5])


@pytest.mark.parametrize("seed", range(10))
def test_large_constant_pool_is_nan_like_the_exact_adapter(seed):
    # One big stratum whose pool is all zeros once the query genes (which carry
    # every nonzero score) are removed: the null variance is exactly 0.
    # Summation error in the stratum's sum of squares grows with its size, so
    # the zero-variance rule must scale with it.
    rng = np.random.default_rng(seed)
    n = 8000
    genes = np.sort(rng.choice(n, size=int(rng.integers(2, 6)), replace=False))
    target = np.zeros(n)
    target[genes] = rng.uniform(1.0, 20.0, genes.size)
    dense = np.column_stack([target, np.full(n, 1.0)])
    reference = analytical_gene_set_z(_StubHetMat(dense), METAPATH, genes, 0, min_stratum_size=5000)
    assert math.isnan(reference.z)
    assert math.isnan(_summary_z(dense, genes, min_stratum_size=5000).z)
