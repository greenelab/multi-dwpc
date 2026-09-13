"""The null bundle: precomputed stratum summaries the web query reads.

Layout of ``bundle_dir``::

    strata.parquet     one row per stratum per (metapath, target), sorted by
                       target_position, then manifest metapath order, then stratum
    row_sums.parquet   per-gene DWPC row sums per metapath
    manifest.json      settings, per-metapath raw_mean and matrix hash, and the
                       SHA-256 of every data file the query path reads

Alpine writes one part per metapath (``write_part``) and ``finalize`` merges
them. ``NullBundle`` refuses a bundle whose data-file hashes do not match the
local ``data/``.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from src.summary_null import DWPC_ZERO_TOL, StratumTable

SCHEMA_VERSION = 1
MIN_STRATUM_SIZE = 50
STRATA_COLUMNS = ["capacity_min", "capacity_max", "n_genes", "score_mean", "score_centered_ss"]
REGENERATE_HINT = "Regenerate it on Alpine with hpc/submit_null_bundle.sh."


class BundleMismatchError(RuntimeError):
    """The bundle is missing or was built from different data files."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def data_files(data_dir: Path) -> list[Path]:
    """Data files the query path reads, relative to ``data_dir``."""
    data_dir = Path(data_dir)
    files = [data_dir / "metagraph.json", data_dir / "metapath-dwpc-stats.tsv"]
    files += sorted((data_dir / "nodes").glob("*.tsv"))
    files += sorted((data_dir / "edges").glob("*.sparse.npz"))
    return [f.relative_to(data_dir) for f in files]


def write_part(parts_dir: Path, index: int, metapath: str, strata: dict[str, np.ndarray],
               row_sums: np.ndarray, raw_mean: float, matrix_sha256: str) -> None:
    """Write one metapath's strata rows, row sums and metadata."""
    parts_dir = Path(parts_dir)
    parts_dir.mkdir(parents=True, exist_ok=True)
    n_rows = len(strata["target_position"])
    table = pa.table({
        "metapath": pa.array([metapath] * n_rows).dictionary_encode(),
        "target_position": pa.array(strata["target_position"], pa.int32()),
        "stratum": pa.array(strata["stratum"], pa.int16()),
        "capacity_min": pa.array(strata["capacity_min"], pa.float64()),
        "capacity_max": pa.array(strata["capacity_max"], pa.float64()),
        "n_genes": pa.array(strata["n_genes"], pa.int32()),
        "score_mean": pa.array(strata["score_mean"], pa.float64()),
        "score_centered_ss": pa.array(strata["score_centered_ss"], pa.float64()),
    })
    pq.write_table(table, parts_dir / f"strata_{index:02d}.parquet")
    pq.write_table(pa.table({
        "metapath": pa.array([metapath] * row_sums.size).dictionary_encode(),
        "gene_position": pa.array(np.arange(row_sums.size), pa.int32()),
        "row_sum": pa.array(row_sums, pa.float64()),
    }), parts_dir / f"row_sums_{index:02d}.parquet")
    meta = {"metapath": metapath, "raw_mean": float(raw_mean), "matrix_sha256": matrix_sha256,
            "n_strata_rows": n_rows}
    (parts_dir / f"part_{index:02d}.json").write_text(json.dumps(meta, indent=2))


def finalize(parts_dir: Path, bundle_dir: Path, data_dir: Path, *, damping: float,
             expected_parts: int | None = None) -> dict:
    """Merge every part in ``parts_dir`` into ``bundle_dir`` and write the manifest."""
    parts_dir, bundle_dir, data_dir = Path(parts_dir), Path(bundle_dir), Path(data_dir)
    indices = sorted(int(p.stem.split("_")[1]) for p in parts_dir.glob("part_*.json"))
    if expected_parts is not None and len(indices) != expected_parts:
        raise ValueError(f"{parts_dir} has {len(indices)} parts, expected {expected_parts}")
    metas = [json.loads((parts_dir / f"part_{i:02d}.json").read_text()) for i in indices]
    bundle_dir.mkdir(parents=True, exist_ok=True)

    strata = pa.concat_tables(
        [pq.read_table(parts_dir / f"strata_{i:02d}.parquet") for i in indices],
        promote_options="permissive",
    ).unify_dictionaries().combine_chunks()
    # Parts are concatenated in manifest order, each sorted by (target, stratum);
    # a stable sort on target gives (target_position, metapath, stratum).
    order = pc.sort_indices(strata, sort_keys=[("target_position", "ascending")])
    pq.write_table(strata.take(order), bundle_dir / "strata.parquet", row_group_size=16384,
                   compression="zstd")
    row_sums = pa.concat_tables(
        [pq.read_table(parts_dir / f"row_sums_{i:02d}.parquet") for i in indices],
        promote_options="permissive",
    )
    pq.write_table(row_sums, bundle_dir / "row_sums.parquet", compression="zstd")

    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                check=True, cwd=Path(__file__).resolve().parent).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": commit,
        "damping": damping,
        "min_stratum_size": MIN_STRATUM_SIZE,
        "dwpc_zero_tol": DWPC_ZERO_TOL,
        "n_strata_rows": strata.num_rows,
        "metapaths": metas,
        "data_files": {str(f): sha256(data_dir / f) for f in data_files(data_dir)},
    }
    (bundle_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


class NullBundle:
    """Read-only view of a finalized bundle, checked against local data files."""

    def __init__(self, bundle_dir: Path, data_dir: Path):
        self.bundle_dir = Path(bundle_dir)
        self.data_dir = Path(data_dir)
        manifest_path = self.bundle_dir / "manifest.json"
        if not manifest_path.exists():
            raise BundleMismatchError(f"No null bundle at {self.bundle_dir}. {REGENERATE_HINT}")
        self.manifest = json.loads(manifest_path.read_text())
        mismatched = [
            name for name, digest in self.manifest["data_files"].items()
            if not (Path(data_dir) / name).exists() or sha256(Path(data_dir) / name) != digest
        ]
        if mismatched:
            raise BundleMismatchError(
                f"Null bundle {self.bundle_dir} was built from different data files "
                f"({', '.join(mismatched[:5])}{' ...' if len(mismatched) > 5 else ''}). {REGENERATE_HINT}"
            )
        self.metapaths = [m["metapath"] for m in self.manifest["metapaths"]]
        self._raw_mean = {m["metapath"]: m["raw_mean"] for m in self.manifest["metapaths"]}
        table = pq.read_table(self.bundle_dir / "row_sums.parquet")
        metapath = table.column("metapath").to_numpy().astype(str)
        position = table.column("gene_position").to_numpy()
        values = table.column("row_sum").to_numpy()
        self._row_sums = {}
        for mp in self.metapaths:
            mask = metapath == mp
            sums = np.zeros(int(position[mask].max()) + 1)
            sums[position[mask]] = values[mask]
            self._row_sums[mp] = sums

    def raw_mean(self, metapath: str) -> float:
        return self._raw_mean[metapath]

    def row_sums(self, metapath: str) -> np.ndarray:
        return self._row_sums[metapath]

    def strata_for_target(self, target_position: int) -> dict[str, StratumTable]:
        """Stratum tables of every metapath for one target."""
        table = pq.read_table(self.bundle_dir / "strata.parquet",
                              filters=[("target_position", "==", int(target_position))])
        metapath = table.column("metapath").to_numpy().astype(str)
        columns = {name: table.column(name).to_numpy() for name in STRATA_COLUMNS}
        names, starts = np.unique(metapath, return_index=True)
        bounds = np.append(np.sort(starts), metapath.size)
        by_start = {int(start): str(name) for start, name in zip(starts, names)}
        return {
            by_start[start]: StratumTable(**{name: col[start:end] for name, col in columns.items()})
            for start, end in zip(bounds[:-1].tolist(), bounds[1:].tolist())
        }
