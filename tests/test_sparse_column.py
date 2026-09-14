import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest
from scipy import sparse

from src.sparse_column import dense_column


def _reference(matrix, position):
    return np.asarray(matrix[:, position].todense()).ravel()


def _assert_identical(got, want):
    assert got.dtype == want.dtype
    assert got.shape == want.shape
    np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize("seed", range(15))
def test_matches_column_slice_on_random_csr(seed):
    rng = np.random.default_rng(seed)
    n_rows, n_cols = int(rng.integers(1, 400)), int(rng.integers(1, 300))
    density = float(rng.choice([0.0, 0.001, 0.05, 0.4, 1.0]))
    m = sparse.random(n_rows, n_cols, density=density, format="csr", random_state=seed)
    assert m.has_canonical_format
    for position in {0, n_cols - 1, int(rng.integers(0, n_cols)), -1}:
        _assert_identical(dense_column(m, position), _reference(m, position))


def test_empty_rows_and_empty_matrix():
    m = sparse.csr_matrix(np.array([[0.0, 2.0, 0.0], [0.0, 0.0, 0.0], [5.0, 0.0, 7.0]]))
    for position in range(3):
        _assert_identical(dense_column(m, position), _reference(m, position))
    empty = sparse.csr_matrix((4, 3), dtype=np.float64)
    _assert_identical(dense_column(empty, 1), _reference(empty, 1))


def test_preserves_integer_and_float32_dtypes():
    dense = np.array([[1, 0, 3], [0, 4, 0]])
    for dtype in (np.int32, np.float32):
        m = sparse.csr_matrix(dense.astype(dtype))
        _assert_identical(dense_column(m, 2), _reference(m, 2))


def test_non_canonical_csr_falls_back_to_slice():
    # Unsorted indices with a duplicate entry (summed by the slice).
    m = sparse.csr_matrix(
        (np.array([1.0, 2.0, 3.0]), np.array([2, 0, 2]), np.array([0, 3])), shape=(1, 3)
    )
    assert not m.has_canonical_format
    _assert_identical(dense_column(m, 2), _reference(m, 2))


def test_csc_input_supported():
    m = sparse.random(50, 20, density=0.2, format="csc", random_state=3)
    _assert_identical(dense_column(m, 7), _reference(m, 7))


def test_out_of_range_position_raises():
    m = sparse.random(5, 4, density=0.5, format="csr", random_state=0)
    for position in (4, -5):
        with pytest.raises(IndexError):
            dense_column(m, position)
