import math
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest
from scipy import sparse

from src.multi_dwpc_query import query_metapath_z

GOOD_MP = "GaDlA"
FLAT_MP = "GaDlB"
BAD_MP = "GbAlG"


class _StubHetMat:
    """Minimal hetmat double covering the calls query_metapath_z makes when
    ``metapaths`` is supplied explicitly (discovery is bypassed): node lookup
    for genes/target, and ``analytical_gene_set_z``'s matrix + stats access.
    """

    def __init__(self):
        rng = np.random.default_rng(3)
        n = 60
        target0 = rng.exponential(scale=1.0, size=n)
        target1 = rng.exponential(scale=2.0, size=n) + 0.1  # capacity key, > 0

        self.source_positions = rng.choice(n, size=10, replace=False)
        enriched = target0.copy()
        enriched[self.source_positions] += 20.0  # planted enrichment

        self._matrices = {
            GOOD_MP: sparse.csc_matrix(np.column_stack([enriched, target1])),
            FLAT_MP: sparse.csc_matrix(
                np.column_stack([np.full(n, 2.0), target1])
            ),  # zero-variance target column
        }
        self.metapath_stats = pd.DataFrame(
            {"metapath": [GOOD_MP, FLAT_MP], "dwpc_raw_mean": [3.0, 3.0]}
        )
        gene_identifiers = np.arange(1, n + 1)
        self._gene_nodes = pd.DataFrame(
            {"identifier": gene_identifiers, "position": np.arange(n)}
        )
        self._target_nodes = pd.DataFrame({"identifier": ["GO:1"], "position": [0]})

    @property
    def gene_ids(self) -> list[int]:
        return (self.source_positions + 1).tolist()

    def get_nodes(self, name):
        if name == "Gene":
            return self._gene_nodes
        if name == "Biological Process":
            return self._target_nodes
        raise KeyError(name)

    def compute_dwpc_matrix_csc(self, metapath, damping=None):
        if metapath not in self._matrices:
            raise KeyError(metapath)
        return self._matrices[metapath]


EXPECTED_COLUMNS = [
    "metapath",
    "real_mean_score",
    "null_mean_score",
    "null_std_score",
    "diff",
    "effect_size_z",
    "p_value",
]


def test_frame_columns_and_ranking():
    hetmat = _StubHetMat()
    df = query_metapath_z(
        hetmat.gene_ids,
        "GO:1",
        metapaths=[GOOD_MP, FLAT_MP],
        hetmat=hetmat,
    )
    assert list(df.columns) == EXPECTED_COLUMNS
    # Finite rows sorted effect_size_z descending; NaN rows sort last (pandas default).
    finite = df["effect_size_z"].dropna()
    assert list(finite) == sorted(finite, reverse=True)
    assert math.isnan(df.iloc[-1]["effect_size_z"])
    assert df.iloc[-1]["metapath"] == FLAT_MP


def test_zero_variance_metapath_row_is_nan_not_dropped():
    hetmat = _StubHetMat()
    df = query_metapath_z(
        hetmat.gene_ids,
        "GO:1",
        metapaths=[GOOD_MP, FLAT_MP],
        hetmat=hetmat,
    )
    flat_row = df[df["metapath"] == FLAT_MP].iloc[0]
    assert math.isnan(flat_row["effect_size_z"])
    assert math.isnan(flat_row["p_value"])


def test_unresolvable_metapath_is_skipped_not_fatal():
    hetmat = _StubHetMat()
    df = query_metapath_z(
        hetmat.gene_ids,
        "GO:1",
        metapaths=[GOOD_MP, BAD_MP],
        hetmat=hetmat,
    )
    assert list(df["metapath"]) == [GOOD_MP]


def test_all_metapaths_unresolvable_raises_value_error():
    hetmat = _StubHetMat()
    with pytest.raises(ValueError, match="No metapaths could be scored"):
        query_metapath_z(
            hetmat.gene_ids,
            "GO:1",
            metapaths=[BAD_MP],
            hetmat=hetmat,
        )


def test_b_and_seed_deprecated_and_inert():
    hetmat = _StubHetMat()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df_plain = query_metapath_z(
            hetmat.gene_ids, "GO:1", metapaths=[GOOD_MP, FLAT_MP], hetmat=hetmat
        )
        assert not any(issubclass(w.category, DeprecationWarning) for w in caught)

    with pytest.warns(DeprecationWarning):
        df_b = query_metapath_z(
            hetmat.gene_ids, "GO:1", metapaths=[GOOD_MP, FLAT_MP], hetmat=hetmat, b=5
        )
    with pytest.warns(DeprecationWarning):
        df_seed = query_metapath_z(
            hetmat.gene_ids, "GO:1", metapaths=[GOOD_MP, FLAT_MP], hetmat=hetmat, seed=1
        )
    with pytest.warns(DeprecationWarning):
        df_both = query_metapath_z(
            hetmat.gene_ids,
            "GO:1",
            metapaths=[GOOD_MP, FLAT_MP],
            hetmat=hetmat,
            b=5,
            seed=1,
        )

    pdt.assert_frame_equal(df_plain, df_b)
    pdt.assert_frame_equal(df_plain, df_seed)
    pdt.assert_frame_equal(df_plain, df_both)


def test_planted_enrichment_yields_high_z():
    hetmat = _StubHetMat()
    df = query_metapath_z(
        hetmat.gene_ids,
        "GO:1",
        metapaths=[GOOD_MP],
        hetmat=hetmat,
    )
    assert df.iloc[0]["effect_size_z"] > 1.65
