import math
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

fake_dwpc_direct = types.ModuleType("dwpc_direct")
fake_dwpc_direct.HetMat = object
fake_dwpc_direct.reverse_metapath_abbrev = lambda metapath: metapath
sys.modules.setdefault("dwpc_direct", fake_dwpc_direct)

fake_result_normalization = types.ModuleType("result_normalization")
fake_result_normalization.load_neo4j_mappings = lambda data_dir: ({}, {})
sys.modules.setdefault("result_normalization", fake_result_normalization)

import dwpc_validation  # noqa: E402



class FakeRandomState:
    def __init__(self, seed):
        self.seed = seed

    def permutation(self, n):
        return list(range(n))


class FakeHetMat:
    def __init__(self, values):
        self.values = values

    def get_dwpc_for_pairs(self, metapath, source_idx, target_idx, transform=True):
        return self.values[: len(source_idx)]


class DwpcValidationTests(unittest.TestCase):
    def test_allocate_samples_stratifies_by_length_and_ignores_length_one(self):
        api_df = pd.DataFrame(
            {
                "metapath_abbreviation": (
                    ["BPpG"] * 20
                    + ["BPpG<rG"] * 20
                    + ["BPpGcG"] * 20
                    + ["BPpG<rG<rG"] * 20
                    + ["BPpGcGcG"] * 20
                )
            }
        )

        allocation = dwpc_validation.allocate_samples_across_metapaths(
            api_df,
            total_samples=9,
        )

        self.assertNotIn("BPpG", allocation)
        self.assertEqual(sum(allocation.values()), 9)
        self.assertEqual(
            allocation,
            {
                "BPpG<rG": 3,
                "BPpGcG": 2,
                "BPpG<rG<rG": 2,
                "BPpGcGcG": 2,
            },
        )

    def test_sample_metapath_concordance_expands_when_initial_sample_is_constant(self):
        api_df = pd.DataFrame(
            {
                "neo4j_source_id": ["go:1"] * 5,
                "neo4j_target_id": [f"gene:{i}" for i in range(5)],
                "metapath_abbreviation": ["BPpG<rG"] * 5,
                "dwpc": [0.0, 0.0, 0.0, 1.0, 2.0],
            }
        )

        with (
            mock.patch.object(
                dwpc_validation,
                "load_neo4j_mappings",
                return_value=(
                    {f"gene:{i}": 100 + i for i in range(5)},
                    {"go:1": "GO:1"},
                ),
            ),
            mock.patch.object(
                dwpc_validation,
                "_build_index_maps",
                return_value=(
                    {100 + i: i for i in range(5)},
                    {"GO:1": 0},
                ),
            ),
            mock.patch.object(dwpc_validation.np.random, "RandomState", FakeRandomState),
        ):
            samples, summary = dwpc_validation.sample_metapath_concordance(
                FakeHetMat([0.0, 0.0, 0.0, 1.0, 2.0]),
                api_df,
                data_dir="unused",
                metapath_api="BPpG<rG",
                n_samples=2,
                seed=7,
            )

        self.assertEqual(summary["n_pairs_requested"], 2)
        self.assertGreater(summary["n_pairs"], 2)
        self.assertTrue(summary["correlation_defined"])
        self.assertAlmostEqual(summary["correlation"], 1.0)
        self.assertEqual(len(samples), summary["n_pairs"])

    def test_sample_metapath_concordance_marks_constant_samples_as_undefined(self):
        api_df = pd.DataFrame(
            {
                "neo4j_source_id": ["go:1"] * 4,
                "neo4j_target_id": [f"gene:{i}" for i in range(4)],
                "metapath_abbreviation": ["BPpG<rG"] * 4,
                "dwpc": [0.0, 0.0, 0.0, 0.0],
            }
        )

        with (
            mock.patch.object(
                dwpc_validation,
                "load_neo4j_mappings",
                return_value=(
                    {f"gene:{i}": 100 + i for i in range(4)},
                    {"go:1": "GO:1"},
                ),
            ),
            mock.patch.object(
                dwpc_validation,
                "_build_index_maps",
                return_value=(
                    {100 + i: i for i in range(4)},
                    {"GO:1": 0},
                ),
            ),
            mock.patch.object(dwpc_validation.np.random, "RandomState", FakeRandomState),
        ):
            _, summary = dwpc_validation.sample_metapath_concordance(
                FakeHetMat([0.0, 0.0, 0.0, 0.0]),
                api_df,
                data_dir="unused",
                metapath_api="BPpG<rG",
                n_samples=2,
                seed=7,
            )

        self.assertFalse(summary["correlation_defined"])
        self.assertTrue(math.isnan(summary["correlation"]))


if __name__ == "__main__":
    unittest.main()
