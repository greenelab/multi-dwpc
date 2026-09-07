import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.analytical_null import analytical_gene_set_z, analytical_null
from src.dwpc_direct import DEFAULT_DAMPING, transform_dwpc
from src.hetnetex_md_import import exact_resampling_moments
from src.hurdle_adaptive_bins import hurdle_adaptive_bins, merge_deficient_strata
from src.pool_assembly import pools_from_bins

METAPATH = "GaDlA"
RAW_MEAN = 3.0


class _StubHetMat:
    """Minimal hetmat double: one metapath, one dense raw DWPC matrix.

    Mirrors tests/test_capacity.py's _StubHetMat pattern, extended with
    ``metapath_stats`` so ``get_dwpc_raw_mean`` resolves.
    """

    def __init__(self, dense_matrix: np.ndarray, raw_mean: float = RAW_MEAN):
        self._matrix = sparse.csc_matrix(dense_matrix)
        self.metapath_stats = pd.DataFrame(
            {"metapath": [METAPATH], "dwpc_raw_mean": [raw_mean]}
        )
        self.calls = []

    def compute_dwpc_matrix_csc(self, metapath, damping=DEFAULT_DAMPING):
        self.calls.append(metapath)
        return self._matrix


def test_one_matrix_load_per_evaluation():
    rng = np.random.default_rng(1)
    n = 60
    dense = rng.exponential(scale=1.0, size=(n, 2))
    hetmat = _StubHetMat(dense)
    source_idx = np.arange(5)
    analytical_gene_set_z(hetmat, METAPATH, source_idx, target_pos=0, min_stratum_size=10)
    assert hetmat.calls == [METAPATH]


def test_planted_enrichment_yields_high_z():
    rng = np.random.default_rng(42)
    n = 400
    target0 = rng.exponential(scale=1.0, size=n)
    target1 = rng.exponential(scale=2.0, size=n) + 0.1  # capacity key, always > 0
    dense = np.column_stack([target0, target1])

    # Plant enrichment: boost a chosen gene set's target-0 (tested) values,
    # leaving their capacity (derived only from target1) untouched.
    source_idx = rng.choice(n, size=20, replace=False)
    dense[source_idx, 0] += 50.0

    hetmat = _StubHetMat(dense)
    result = analytical_gene_set_z(
        hetmat, METAPATH, source_idx, target_pos=0, min_stratum_size=20
    )
    assert result.z > 1.65


def test_null_moments_match_manual_pipeline():
    # Single stratum (all-equal positive capacity, no hurdle split) so the
    # exact moments can be reproduced by hand-driving the same primitives
    # the adapter uses internally.
    n = 12
    target0 = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0])
    target1 = np.full(n, 5.0)  # constant capacity -> one stratum
    dense = np.column_stack([target0, target1])
    hetmat = _StubHetMat(dense)

    source_idx = np.array([0, 1, 2, 3])
    result = analytical_gene_set_z(
        hetmat, METAPATH, source_idx, target_pos=0, min_stratum_size=1
    )

    # Manual pipeline, independent of the adapter's internals.
    matrix = sparse.csc_matrix(dense)
    raw_target_col = np.asarray(matrix[:, 0].todense()).ravel()
    scores = transform_dwpc(raw_target_col, RAW_MEAN)
    capacity = np.asarray(matrix.sum(axis=1)).ravel() - raw_target_col
    bins = hurdle_adaptive_bins(capacity, min_stratum_size=1)
    pools, counts = pools_from_bins(bins, source_idx, int(bins.max()) + 1)
    pools, counts, merges = merge_deficient_strata(pools, counts)
    mean, var, _mu3 = exact_resampling_moments(scores, pools, counts)

    assert result.null_mean == pytest.approx(mean, rel=1e-9, abs=1e-12)
    assert result.null_std == pytest.approx(math.sqrt(var), rel=1e-9, abs=1e-12)
    assert result.merges == merges


def test_stratum_matched_null_draws_are_standard_normal_shaped():
    rng = np.random.default_rng(7)
    n = 500
    target0 = rng.exponential(scale=1.0, size=n)
    target1 = rng.exponential(scale=2.0, size=n) + 0.1
    dense = np.column_stack([target0, target1])
    hetmat = _StubHetMat(dense)

    source_idx = rng.choice(n, size=15, replace=False)
    result = analytical_gene_set_z(
        hetmat, METAPATH, source_idx, target_pos=0, min_stratum_size=20
    )

    # Rebuild the same partition and scores to draw stratum-matched null
    # gene sets and confirm they score standard-normal under (null_mean,
    # null_std) from the adapter.
    matrix = sparse.csc_matrix(dense)
    raw_target_col = np.asarray(matrix[:, 0].todense()).ravel()
    scores = transform_dwpc(raw_target_col, RAW_MEAN)
    capacity = np.asarray(matrix.sum(axis=1)).ravel() - raw_target_col
    bins = hurdle_adaptive_bins(capacity, min_stratum_size=20)
    pools, counts = pools_from_bins(bins, source_idx, int(bins.max()) + 1)
    pools, counts, merges = merge_deficient_strata(pools, counts)

    n_draws = 200
    z_draws = np.empty(n_draws)
    for i in range(n_draws):
        drawn = [
            rng.choice(pool, size=k, replace=False)
            for pool, k in zip(pools, counts)
            if k > 0
        ]
        drawn_idx = np.concatenate(drawn)
        t = scores[drawn_idx].mean()
        z_draws[i] = (t - result.null_mean) / result.null_std

    assert abs(z_draws.mean()) < 0.15
    assert abs(z_draws.std() - 1.0) < 0.15


def test_zero_variance_column_yields_nan_without_raising():
    n = 40
    target0 = np.full(n, 2.0)  # constant -> zero variance under any draw
    target1 = np.linspace(0.1, 5.0, n)
    dense = np.column_stack([target0, target1])
    hetmat = _StubHetMat(dense)

    source_idx = np.array([0, 1, 2, 3, 4])
    result = analytical_gene_set_z(
        hetmat, METAPATH, source_idx, target_pos=0, min_stratum_size=5
    )
    assert math.isnan(result.z)
    assert math.isnan(result.p_value)
    assert math.isfinite(result.real_mean)
    assert math.isfinite(result.null_mean)


def test_forced_deficient_stratum_reports_merges():
    # Two positive-capacity strata (values 1 and 2, min_stratum_size small
    # enough to keep them separate); nearly all of the low-capacity
    # stratum's genes are placed in the queried gene set, leaving too few
    # candidates behind to satisfy that stratum's count.
    n = 10
    target1 = np.array([1.0] * 5 + [2.0] * 5)  # capacity key
    target0 = np.arange(1.0, n + 1.0)
    dense = np.column_stack([target0, target1])
    hetmat = _StubHetMat(dense)

    source_idx = np.array([0, 1, 2, 3])  # 4 of the 5 capacity==1 genes
    result = analytical_gene_set_z(
        hetmat, METAPATH, source_idx, target_pos=0, min_stratum_size=1
    )
    assert len(result.merges) > 0


def test_min_stratum_size_zero_raises():
    n = 20
    dense = np.column_stack(
        [np.arange(1.0, n + 1.0), np.arange(1.0, n + 1.0)]
    )
    hetmat = _StubHetMat(dense)
    with pytest.raises(ValueError):
        analytical_gene_set_z(
            hetmat, METAPATH, np.array([0, 1]), target_pos=0, min_stratum_size=0
        )
