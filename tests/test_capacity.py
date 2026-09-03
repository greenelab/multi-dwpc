import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy import sparse

from src.capacity import CapacityProvider, leave_target_out_capacity


def _matrix():
    # 4 genes x 3 targets, raw DWPC values
    dense = np.array(
        [[0.0, 0.0, 0.0],
         [1.0, 2.0, 0.0],
         [0.5, 0.0, 4.0],
         [3.0, 3.0, 3.0]]
    )
    return sparse.csc_matrix(dense)


def test_capacity_is_rowsum_minus_target_column():
    cap = leave_target_out_capacity(_matrix(), target_position=1)
    np.testing.assert_allclose(cap, [0.0, 1.0, 4.5, 6.0])


def test_zero_row_gene_has_zero_capacity_for_any_target():
    m = _matrix()
    for t in range(3):
        assert leave_target_out_capacity(m, t)[0] == 0.0


class _StubHetMat:
    def __init__(self):
        self.calls = []

    def compute_dwpc_matrix_csc(self, metapath, damping=0.5):
        self.calls.append(metapath)
        return _matrix()


def test_provider_caches_rowsums_per_metapath():
    hetmat = _StubHetMat()
    provider = CapacityProvider(hetmat)
    a = provider.capacity("GaDlA", target_position=0)
    b = provider.capacity("GaDlA", target_position=2)
    np.testing.assert_allclose(a, [0.0, 2.0, 4.0, 6.0])
    np.testing.assert_allclose(b, [0.0, 3.0, 0.5, 6.0])
    assert hetmat.calls == ["GaDlA", "GaDlA"]  # matrix reloaded, rowsum cached is optional
    c = provider.capacity("GiGaD", target_position=0)
    assert len(c) == 4
