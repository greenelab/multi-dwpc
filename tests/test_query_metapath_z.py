import math
import subprocess
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from scipy import sparse

from src.analytical_null import analytical_gene_set_z
from src.dwpc_direct import transform_dwpc
from src.multi_dwpc_query import query_metapath_z
from src.query_dwpc import UnsupportedMetapathError
from src.summary_null import build_strata

GOOD_MP = "GaDlA"
FLAT_MP = "GaDlB"
BAD_MP = "GbAlG"
RAW_MEAN = 3.0
N = 60
REPO_ROOT = Path(__file__).resolve().parents[1]


def _matrices():
    rng = np.random.default_rng(3)
    target0 = rng.exponential(scale=1.0, size=N)
    target1 = rng.exponential(scale=2.0, size=N) + 0.1  # capacity key, > 0
    source_positions = rng.choice(N, size=10, replace=False)
    enriched = target0.copy()
    enriched[source_positions] += 20.0  # planted enrichment
    return source_positions, {
        GOOD_MP: sparse.csr_matrix(np.column_stack([enriched, target1])),
        FLAT_MP: sparse.csr_matrix(np.column_stack([np.full(N, 2.0), target1])),  # zero variance
    }


class _StubBundle:
    """Bundle built from small matrices exactly as build_null_bundle.py would."""

    def __init__(self, data_dir, matrices, metapaths):
        self.data_dir = data_dir
        self.metapaths = metapaths
        self._matrices = matrices

    def raw_mean(self, metapath):
        return RAW_MEAN

    def row_sums(self, metapath):
        return np.asarray(self._matrices[metapath].sum(axis=1)).ravel()

    def strata_for_target(self, target_position):
        tables = {}
        for metapath, matrix in self._matrices.items():
            column = matrix[:, target_position].toarray().ravel()
            capacity = self.row_sums(metapath) - column
            tables[metapath] = build_strata(capacity, transform_dwpc(column, RAW_MEAN), 50)
        return tables


class _StubQueryDwpc:
    def __init__(self, matrices):
        self._matrices = matrices

    def target_values(self, metapath, genes, target_position):
        if metapath not in self._matrices:
            raise UnsupportedMetapathError(metapath)
        return self._matrices[metapath][genes, target_position].toarray().ravel()


class _ReferenceHetMat:
    def __init__(self, matrices):
        self._matrices = matrices
        self.metapath_stats = pd.DataFrame({"metapath": list(matrices), "dwpc_raw_mean": RAW_MEAN})

    def compute_dwpc_matrix(self, metapath, damping=None):
        return self._matrices[metapath]

    def get_dwpc_row_sums(self, metapath, damping=None):
        return np.asarray(self._matrices[metapath].sum(axis=1)).ravel()


@pytest.fixture
def setup(tmp_path):
    nodes = tmp_path / "nodes"
    nodes.mkdir()
    pd.DataFrame({"identifier": np.arange(1, N + 1), "name": [f"g{i}" for i in range(N)],
                  "position": np.arange(N)}).to_csv(nodes / "Gene.tsv", sep="\t", index=False)
    pd.DataFrame({"identifier": ["GO:1"], "name": ["bp"], "position": [0]}).to_csv(
        nodes / "Biological Process.tsv", sep="\t", index=False)
    source_positions, matrices = _matrices()
    kwargs = dict(bundle=_StubBundle(tmp_path, matrices, [GOOD_MP, FLAT_MP]),
                  query_dwpc=_StubQueryDwpc(matrices))
    return (source_positions + 1).tolist(), matrices, kwargs


EXPECTED_COLUMNS = ["metapath", "real_mean_score", "null_mean_score", "null_std_score", "diff",
                    "effect_size_z", "p_value"]


def test_matches_the_matrix_based_adapter(setup):
    gene_ids, matrices, kwargs = setup
    df = query_metapath_z(gene_ids, "GO:1", **kwargs).set_index("metapath")
    positions = np.array(sorted(g - 1 for g in gene_ids))
    for metapath in matrices:
        reference = analytical_gene_set_z(_ReferenceHetMat(matrices), metapath, positions, 0)
        row = df.loc[metapath]
        for column, value in [("real_mean_score", reference.real_mean), ("null_mean_score", reference.null_mean),
                              ("effect_size_z", reference.z)]:
            if math.isnan(value):
                assert math.isnan(row[column])
            else:
                assert row[column] == pytest.approx(value, rel=1e-12)


def test_frame_columns_and_ranking(setup):
    gene_ids, _, kwargs = setup
    df = query_metapath_z(gene_ids, "GO:1", **kwargs)
    assert list(df.columns) == EXPECTED_COLUMNS
    assert df.iloc[0]["metapath"] == GOOD_MP
    assert df.iloc[0]["effect_size_z"] > 1.65


def test_zero_variance_metapath_row_is_nan_not_dropped(setup):
    gene_ids, _, kwargs = setup
    row = query_metapath_z(gene_ids, "GO:1", **kwargs).set_index("metapath").loc[FLAT_MP]
    assert math.isnan(row["effect_size_z"])
    assert math.isnan(row["p_value"])


def test_duplicate_gene_ids_count_once(setup):
    gene_ids, _, kwargs = setup
    pdt.assert_frame_equal(query_metapath_z(gene_ids + gene_ids[:3], "GO:1", **kwargs),
                           query_metapath_z(gene_ids, "GO:1", **kwargs))


def test_unsupported_metapath_is_skipped_not_fatal(setup):
    gene_ids, _, kwargs = setup
    df = query_metapath_z(gene_ids, "GO:1", metapaths=[GOOD_MP, BAD_MP], **kwargs)
    assert df["metapath"].tolist() == [GOOD_MP]


def test_all_metapaths_unsupported_raises_value_error(setup):
    gene_ids, _, kwargs = setup
    with pytest.raises(ValueError):
        query_metapath_z(gene_ids, "GO:1", metapaths=[BAD_MP], **kwargs)


def test_unknown_genes_raise_value_error(setup):
    _, _, kwargs = setup
    with pytest.raises(ValueError, match="None of the provided gene IDs"):
        query_metapath_z([999999], "GO:1", **kwargs)


def test_b_and_seed_deprecated_and_inert(setup):
    gene_ids, _, kwargs = setup
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        plain = query_metapath_z(gene_ids, "GO:1", **kwargs)
    with pytest.warns(DeprecationWarning):
        with_b = query_metapath_z(gene_ids, "GO:1", b=5, seed=1, **kwargs)
    pdt.assert_frame_equal(plain, with_b)


def test_query_path_does_not_load_heavy_null_dependencies():
    # hetnetex_md and scipy.stats (~50 MB resident) belong to the matrix-based
    # reference only; the web query path must not import them.
    code = ("import sys; sys.path.insert(0, '.'); import src.multi_dwpc_query; "
            "print('hetnetex_md' in sys.modules, 'scipy.stats' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    assert out.stdout.split()[-2:] == ["False", "False"]
