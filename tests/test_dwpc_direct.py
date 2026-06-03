"""Unit tests for the pure helpers in ``src/dwpc_direct.py``.

These cover the metapath string handling, disk-cache key construction, and
node-to-matrix-index mapping. They use a small synthetic metagraph so they do
not depend on the bundled ``data/metagraph.json`` and run without any matrix
files on disk.
"""

import sys
import types
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.dwpc_direct import (  # noqa: E402
    HetMat,
    create_node_index_mapping,
    parse_metapath,
    reverse_metapath_abbrev,
)


# A minimal metagraph with a multi-character node abbreviation ("BP") so the
# parser's shared-node advancement is exercised, not just single-char nodes.
METAGRAPH = {
    "metanode_kinds": ["Gene", "Biological Process", "Pathway"],
    "kind_to_abbrev": {
        "Gene": "G",
        "Biological Process": "BP",
        "Pathway": "PW",
        "participates": "p",
        "interacts": "i",
    },
    "metaedge_tuples": [
        ["Gene", "Biological Process", "participates", "both"],
        ["Gene", "Gene", "interacts", "both"],
        ["Gene", "Pathway", "participates", "both"],
    ],
}


class ReverseMetapathAbbrevTests(unittest.TestCase):
    def test_simple_reverse(self):
        self.assertEqual(reverse_metapath_abbrev("GpBP"), "BPpG")

    def test_multi_edge_reverse(self):
        self.assertEqual(reverse_metapath_abbrev("GiGpBP"), "BPpGiG")

    def test_multi_character_nodes_round_trip(self):
        # Reversing twice returns the original metapath.
        self.assertEqual(reverse_metapath_abbrev(reverse_metapath_abbrev("GpBP")), "GpBP")

    def test_direction_arrows_round_trip(self):
        # Reversing twice restores the original, including directed edges.
        self.assertEqual(reverse_metapath_abbrev(reverse_metapath_abbrev("Gr>G")), "Gr>G")

    def test_unrecognized_character_raises(self):
        with self.assertRaises(ValueError):
            reverse_metapath_abbrev("Gp$BP")


class ParseMetapathTests(unittest.TestCase):
    def test_single_edge(self):
        edges = parse_metapath("GpBP", METAGRAPH)
        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0]["source"], "Gene")
        self.assertEqual(edges[0]["target"], "Biological Process")

    def test_two_edges_through_repeated_gene(self):
        edges = parse_metapath("GiGpBP", METAGRAPH)
        self.assertEqual([e["abbrev"] for e in edges], ["GiG", "GpBP"])

    def test_multi_character_node_advancement(self):
        # "BP" is two characters; the parser must resume at the shared node,
        # not one character past the matched edge.
        edges = parse_metapath("GpBPpGpBP", METAGRAPH)
        self.assertEqual(len(edges), 3)

    def test_lone_trailing_node_terminates(self):
        # A trailing metanode is not an error; it just ends the walk.
        self.assertEqual(len(parse_metapath("GpBP", METAGRAPH)), 1)

    def test_unparseable_metapath_raises(self):
        for bad in ("Zz", "GpBPX", "Gq", "GpB"):
            with self.subTest(metapath=bad):
                with self.assertRaises(ValueError):
                    parse_metapath(bad, METAGRAPH)


class CachePathTests(unittest.TestCase):
    def _path(self, damping):
        # _get_cache_path only reads self.cache_dir, so a stub avoids the heavy
        # HetMat.__init__ (which would load matrices and metagraph from disk).
        stub = types.SimpleNamespace(cache_dir=Path("/tmp/dwpc_cache"))
        return HetMat._get_cache_path(stub, "GpBP", damping)

    def test_distinct_damping_values_do_not_alias(self):
        # The old "{damping:.2f}" key mapped 0.33 and 0.333 to the same file.
        self.assertNotEqual(self._path(0.33), self._path(0.333))

    def test_key_round_trips_value(self):
        self.assertEqual(self._path(0.5).name, "dwpc_GpBP_d0.5.npz")
        self.assertEqual(self._path(0.333).name, "dwpc_GpBP_d0.333.npz")


class FakeHetMat:
    """Stand-in exposing only ``get_nodes`` for create_node_index_mapping."""

    def __init__(self, nodes_by_type):
        self._nodes_by_type = nodes_by_type

    def get_nodes(self, node_type):
        return self._nodes_by_type[node_type]


class CreateNodeIndexMappingTests(unittest.TestCase):
    def test_prefers_position_column_over_dataframe_index(self):
        # positions deliberately differ from the default 0..n-1 index so the
        # test fails if .index were used instead of "position".
        genes = pd.DataFrame({"position": [10, 11, 12], "identifier": ["a", "b", "c"]})
        bps = pd.DataFrame({"position": [20, 21], "identifier": ["x", "y"]})
        hetmat = FakeHetMat({"Gene": genes, "Biological Process": bps})

        df = pd.DataFrame({"gene": ["a", "c"], "bp": ["y", "x"]})
        out = create_node_index_mapping(
            hetmat, df, "Gene", "Biological Process", "gene", "bp"
        )

        self.assertEqual(list(out["source_idx"]), [10, 12])
        self.assertEqual(list(out["target_idx"]), [21, 20])

    def test_falls_back_to_index_when_no_position_column(self):
        genes = pd.DataFrame({"identifier": ["a", "b", "c"]})
        bps = pd.DataFrame({"identifier": ["x", "y"]})
        hetmat = FakeHetMat({"Gene": genes, "Biological Process": bps})

        df = pd.DataFrame({"gene": ["a", "c"], "bp": ["y", "x"]})
        out = create_node_index_mapping(
            hetmat, df, "Gene", "Biological Process", "gene", "bp"
        )

        self.assertEqual(list(out["source_idx"]), [0, 2])
        self.assertEqual(list(out["target_idx"]), [1, 0])


if __name__ == "__main__":
    unittest.main()
