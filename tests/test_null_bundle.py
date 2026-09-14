import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import json

import numpy as np
import pytest

from src.null_bundle import BundleMismatchError, NullBundle, data_fingerprint, finalize, write_part
from src.summary_null import build_strata

METAPATHS = ["GbCbGpBP", "G<rGpBP"]  # includes a '<' to exercise names
N_GENES, N_TARGETS = 300, 4


def _fake_data_dir(root: Path) -> Path:
    data = root / "data"
    (data / "nodes").mkdir(parents=True)
    (data / "edges").mkdir()
    (data / "metagraph.json").write_text("{}")
    (data / "metapath-dwpc-stats.tsv").write_text("metapath\tdwpc_raw_mean\n")
    (data / "nodes" / "Gene.tsv").write_text("identifier\tname\tposition\n")
    (data / "edges" / "GpBP.sparse.npz").write_bytes(b"edges")
    return data


def _build(root: Path, fingerprint: str | None = None):
    data = _fake_data_dir(root)
    fingerprint = fingerprint or data_fingerprint(data)
    rng = np.random.default_rng(0)
    expected = {}
    for index, metapath in enumerate(METAPATHS):
        row_sums = rng.exponential(2.0, N_GENES)
        rows = {k: [] for k in ["target_position", "stratum", "capacity_min", "capacity_max",
                                "n_genes", "score_mean", "score_centered_ss"]}
        for target in range(N_TARGETS):
            column = rng.exponential(1.0, N_GENES) * (rng.random(N_GENES) < 0.5)
            table = build_strata(np.maximum(row_sums - column, 0.0), column, min_stratum_size=40)
            expected[(metapath, target)] = table
            n = table.n_genes.size
            rows["target_position"].append(np.full(n, target))
            rows["stratum"].append(np.arange(n))
            for name in ["capacity_min", "capacity_max", "n_genes", "score_mean", "score_centered_ss"]:
                rows[name].append(getattr(table, name))
        write_part(root / "parts", index, metapath, {k: np.concatenate(v) for k, v in rows.items()},
                   row_sums, raw_mean=1.5 + index, matrix_sha256="abc", data_fingerprint=fingerprint)
        expected[(metapath, "row_sums")] = row_sums
    finalize(root / "parts", root / "bundle", data, damping=0.5, expected_parts=len(METAPATHS))
    return data, expected


def test_round_trip(tmp_path):
    data, expected = _build(tmp_path)
    bundle = NullBundle(tmp_path / "bundle", data)
    assert bundle.metapaths == METAPATHS
    assert bundle.raw_mean("G<rGpBP") == 2.5
    for metapath in METAPATHS:
        np.testing.assert_array_equal(bundle.row_sums(metapath), expected[(metapath, "row_sums")])
    for target in range(N_TARGETS):
        tables = bundle.strata_for_target(target)
        assert sorted(tables) == sorted(METAPATHS)
        for metapath, table in tables.items():
            want = expected[(metapath, target)]
            for name in ["capacity_min", "capacity_max", "n_genes", "score_mean", "score_centered_ss"]:
                np.testing.assert_array_equal(getattr(table, name), getattr(want, name))


def test_changed_data_file_is_refused(tmp_path):
    data, _ = _build(tmp_path)
    (data / "edges" / "GpBP.sparse.npz").write_bytes(b"different edges")
    with pytest.raises(BundleMismatchError, match="edges/GpBP.sparse.npz"):
        NullBundle(tmp_path / "bundle", data)


def test_missing_bundle_is_refused(tmp_path):
    with pytest.raises(BundleMismatchError, match="No null bundle"):
        NullBundle(tmp_path / "nowhere", tmp_path)


def test_data_fingerprint_tracks_data_files(tmp_path):
    data = _fake_data_dir(tmp_path)
    first = data_fingerprint(data)
    assert data_fingerprint(data) == first
    (data / "edges" / "GpBP.sparse.npz").write_bytes(b"other edges")
    assert data_fingerprint(data) != first


def test_parts_built_from_other_data_are_refused_at_finalize(tmp_path):
    # A matrix cache left over from another dataset must not end up in a bundle
    # that is then stamped with the current data's hashes.
    with pytest.raises(ValueError, match="different data"):
        _build(tmp_path, fingerprint="0" * 16)


@pytest.mark.parametrize("key, value", [("schema_version", 2), ("damping", 0.4),
                                        ("min_stratum_size", 25), ("dwpc_zero_tol", 1e-12)])
def test_incompatible_settings_are_refused(tmp_path, key, value):
    data, _ = _build(tmp_path)
    manifest_path = tmp_path / "bundle" / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest[key] = value
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(BundleMismatchError, match=key):
        NullBundle(tmp_path / "bundle", data)


@pytest.mark.parametrize("missing", ["strata.parquet", "row_sums.parquet"])
def test_incomplete_bundle_is_refused(tmp_path, missing):
    data, _ = _build(tmp_path)
    (tmp_path / "bundle" / missing).unlink()
    with pytest.raises(BundleMismatchError, match=missing):
        NullBundle(tmp_path / "bundle", data)
