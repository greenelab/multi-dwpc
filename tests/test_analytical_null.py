import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest

from src.analytical_null import analytical_null


def test_normal_case_z_and_p():
    rng = np.random.default_rng(0)
    scores = rng.normal(size=200)
    pools = [np.arange(200)]
    counts = [20]
    result = analytical_null(scores, pools, counts, observed=scores[:20].mean() + 5.0)
    assert math.isfinite(result.z)
    assert result.z > 0
    assert 0.0 <= result.p <= 1.0


def test_zero_variance_maps_to_nan_not_inf():
    # k == N: SRSWOR with no randomness left -> var == 0.
    scores = np.array([1.0, 2.0, 3.0])
    result = analytical_null(scores, [np.arange(3)], [3], observed=10.0)
    assert math.isnan(result.z)
    assert math.isnan(result.p)


def test_k_total_zero_raises():
    scores = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="K"):
        analytical_null(scores, [np.arange(3)], [0], observed=1.0)


def test_non_finite_scores_raises():
    scores = np.array([1.0, np.nan, 3.0])
    with pytest.raises(ValueError, match="finite"):
        analytical_null(scores, [np.arange(3)], [2], observed=1.0)


def test_oversized_k_raises():
    scores = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="exceeds pool"):
        analytical_null(scores, [np.arange(3)], [5], observed=1.0)
