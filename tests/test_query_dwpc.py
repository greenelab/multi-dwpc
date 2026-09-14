import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pytest
from scipy import sparse

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
CACHE = DATA_DIR / "dwpc_cache"

pytestmark = pytest.mark.skipif(not (DATA_DIR / "edges").exists(), reason="requires bundled data/")


@pytest.fixture(scope="module")
def query_dwpc():
    from src.query_dwpc import QueryDwpc

    return QueryDwpc(DATA_DIR)


@pytest.fixture(scope="module")
def genes():
    rng = np.random.default_rng(11)
    return np.unique(np.concatenate([rng.choice(20945, 150, replace=False), [0, 20944]]))


# One metapath per supported shape: no repeats, two-G repeat with an insert,
# adjacent G-G repeat, three-G repeat, and BABA (1 GB cached matrix).
@pytest.mark.parametrize("metapath", ["GpBP", "GbCbGpBP", "GcGpBP", "G<rGcGpBP", "GpBPpGpBP"])
def test_rows_match_cached_hetmatpy_matrix(query_dwpc, genes, metapath):
    path = CACHE / f"dwpc_{metapath}_d0.5.npz"
    if not path.exists():
        pytest.skip(f"{path.name} not cached")
    reference = sparse.load_npz(path).tocsr()[genes].toarray()
    got = query_dwpc.rows(metapath, genes).toarray()
    big = np.abs(reference) > 1e-15
    np.testing.assert_array_equal(big, np.abs(got) > 1e-15)
    np.testing.assert_allclose(got[big], reference[big], rtol=1e-12, atol=0)
    # Everything else is a cancellation residue of a true zero.
    assert np.abs(got[~big]).max(initial=0.0) <= 1e-15


def test_target_values_apply_the_zero_rule(query_dwpc, genes):
    values = query_dwpc.target_values("GpBPpGpBP", genes, 7)
    assert not np.any((values != 0) & (np.abs(values) <= 1e-15))


def test_gene_without_edges_has_zero_dwpc(query_dwpc):
    degree = np.asarray(query_dwpc.rows("GpBP", np.arange(20945)).sum(axis=1)).ravel()
    isolated = np.flatnonzero(degree == 0)[:3]
    assert isolated.size
    assert np.all(query_dwpc.target_values("GcGpBP", isolated, 0) == 0)


def test_unsupported_metapaths_raise(query_dwpc):
    from src.query_dwpc import UnsupportedMetapathError

    with pytest.raises(UnsupportedMetapathError):
        query_dwpc.rows("GiGbCrC", np.array([0]))  # category 'disjoint'
    with pytest.raises(UnsupportedMetapathError):
        query_dwpc.rows("GxBP", np.array([0]))  # not in the metagraph
