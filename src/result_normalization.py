"""Normalize direct and API DWPC outputs into a shared result schema."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


YEAR_DATASET_RE = re.compile(
    r"^(?P<base>.+)_(?P<year>\d{4})_(?P<kind>real|perm|random)(?:_(?P<rep>\d+))?$"
)
YEAR_SHARED_COLUMNS = [
    "domain",
    "name",
    "year",
    "control",
    "replicate",
    "score_source",
    "go_id",
    "entrez_gene_id",
    "metapath",
    "dwpc",
]
YEAR_SUMMARY_COLUMNS = [
    "domain",
    "name",
    "control",
    "replicate",
    "year",
    "go_id",
    "metapath",
    "mean_score",
    "n_pairs",
]


def parse_year_dataset_name(dataset_name: str) -> dict:
    """Parse a year dataset name into standardized metadata."""
    match = YEAR_DATASET_RE.match(str(dataset_name))
    if not match:
        raise ValueError(f"Unrecognized year dataset name: {dataset_name}")
    kind = match.group("kind")
    control_map = {"real": "real", "perm": "permuted", "random": "random"}
    rep = int(match.group("rep")) if match.group("rep") else 0
    return {
        "base_name": str(match.group("base")),
        "year": int(match.group("year")),
        "control": control_map[kind],
        "replicate": int(rep),
        "name": str(dataset_name),
    }


def load_neo4j_mappings(data_dir: Path) -> tuple[dict, dict]:
    """Load mappings from Neo4j node ids to external GO and Entrez identifiers."""
    gene_map = pd.read_csv(Path(data_dir) / "neo4j_gene_mapping.csv")
    bp_map = pd.read_csv(Path(data_dir) / "neo4j_bp_mapping.csv")
    neo4j_to_entrez = dict(zip(gene_map["neo4j_id"], gene_map["identifier"]))
    neo4j_to_go = dict(zip(bp_map["neo4j_id"], bp_map["identifier"]))
    return neo4j_to_entrez, neo4j_to_go


def normalize_direct_year_result(df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """Normalize a direct-DWPC year result frame into the shared schema."""
    meta = parse_year_dataset_name(dataset_name)
    required = {"go_id", "entrez_gene_id", "metapath", "dwpc"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Direct year result {dataset_name} missing required columns: {sorted(missing)}"
        )

    out = df.copy()
    out["metapath"] = out["metapath"].astype(str)
    out.insert(0, "score_source", "direct")
    out.insert(0, "replicate", meta["replicate"])
    out.insert(0, "control", meta["control"])
    out.insert(0, "year", meta["year"])
    out.insert(0, "name", meta["name"])
    out.insert(0, "domain", "year")
    return out[YEAR_SHARED_COLUMNS].copy()


def normalize_api_year_result(
    df: pd.DataFrame,
    dataset_name: str,
    neo4j_to_entrez: dict,
    neo4j_to_go: dict,
) -> pd.DataFrame:
    """Normalize an API year result frame into the shared schema."""
    meta = parse_year_dataset_name(dataset_name)
    required = {"neo4j_source_id", "neo4j_target_id", "dwpc"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"API year result {dataset_name} missing required columns: {sorted(missing)}"
        )

    out = df.copy()
    if "metapath" not in out.columns:
        if "metapath_abbreviation" not in out.columns:
            raise ValueError(
                f"API year result {dataset_name} must include `metapath` or `metapath_abbreviation`."
            )
        out = out.rename(columns={"metapath_abbreviation": "metapath"})

    out["go_id"] = out["neo4j_source_id"].map(neo4j_to_go)
    target_col = "neo4j_pseudo_target_id" if "neo4j_pseudo_target_id" in out.columns else "neo4j_target_id"
    out["entrez_gene_id"] = out[target_col].map(neo4j_to_entrez)
    out = out.dropna(subset=["go_id", "entrez_gene_id", "metapath", "dwpc"]).copy()

    out.insert(0, "score_source", "api")
    out.insert(0, "replicate", meta["replicate"])
    out.insert(0, "control", meta["control"])
    out.insert(0, "year", meta["year"])
    out.insert(0, "name", meta["name"])
    out.insert(0, "domain", "year")
    return out[YEAR_SHARED_COLUMNS].copy()


def _dataset_name_from_result_file(path: Path, score_source: str) -> str:
    stem = path.stem
    if str(score_source) == "direct":
        if stem.startswith("dwpc_"):
            return stem.replace("dwpc_", "", 1)
    elif str(score_source) == "api":
        if stem.startswith("res_"):
            return stem.replace("res_", "", 1)
    else:
        raise ValueError(f"Unsupported score_source: {score_source}")
    return stem


def discover_year_result_files(results_dir: Path | str, score_source: str) -> list[Path]:
    """Discover year-domain result files for a score source."""
    results_dir = Path(results_dir)
    if str(score_source) == "direct":
        files = sorted(results_dir.glob("dwpc_*.csv"))
    elif str(score_source) == "api":
        files = sorted(results_dir.glob("res_*.csv"))
    else:
        raise ValueError(f"Unsupported score_source: {score_source}")

    valid_files: list[Path] = []
    for path in files:
        dataset_name = _dataset_name_from_result_file(path, str(score_source))
        try:
            parse_year_dataset_name(dataset_name)
        except ValueError:
            continue
        valid_files.append(path)
    return valid_files


def load_normalized_year_results(
    results_dir: Path | str,
    score_source: str,
    data_dir: Path | str | None = None,
) -> pd.DataFrame:
    """
    Load and normalize all year result files for a score source.

    This is intentionally source-only normalization: no API calls and no
    downstream replicate analysis assumptions.
    """
    paths = discover_year_result_files(results_dir, score_source=score_source)
    if not paths:
        raise FileNotFoundError(
            f"No parseable year result files found under {Path(results_dir)} for source={score_source}"
        )

    neo4j_to_entrez: dict | None = None
    neo4j_to_go: dict | None = None
    if str(score_source) == "api":
        if data_dir is None:
            raise ValueError("`data_dir` is required to normalize API year results.")
        neo4j_to_entrez, neo4j_to_go = load_neo4j_mappings(Path(data_dir))

    frames: list[pd.DataFrame] = []
    for path in paths:
        dataset_name = _dataset_name_from_result_file(path, str(score_source))
        raw_df = pd.read_csv(path)
        if str(score_source) == "api":
            if "neo4j_target_id" not in raw_df.columns and "neo4j_pseudo_target_id" in raw_df.columns:
                raw_df = raw_df.copy()
                raw_df["neo4j_target_id"] = raw_df["neo4j_pseudo_target_id"]
            norm_df = normalize_api_year_result(
                raw_df,
                dataset_name=dataset_name,
                neo4j_to_entrez=neo4j_to_entrez or {},
                neo4j_to_go=neo4j_to_go or {},
            )
        else:
            norm_df = normalize_direct_year_result(raw_df, dataset_name=dataset_name)
        if not norm_df.empty:
            frames.append(norm_df)

    if not frames:
        return pd.DataFrame(columns=YEAR_SHARED_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def summarize_normalized_year_results(normalized_df: pd.DataFrame) -> pd.DataFrame:
    """Convert normalized year rows into replicate summary rows used downstream."""
    if normalized_df.empty:
        return pd.DataFrame(columns=YEAR_SUMMARY_COLUMNS)

    required = {"domain", "name", "control", "replicate", "year", "go_id", "metapath", "dwpc"}
    missing = required - set(normalized_df.columns)
    if missing:
        raise ValueError(f"Normalized year dataframe missing required columns: {sorted(missing)}")

    summary_df = (
        normalized_df.groupby(
            ["domain", "name", "control", "replicate", "year", "go_id", "metapath"],
            as_index=False,
        )["dwpc"]
        .agg(mean_score="mean", n_pairs="size")
        .sort_values(["year", "control", "replicate", "name", "go_id", "metapath"])
        .reset_index(drop=True)
    )
    return summary_df[YEAR_SUMMARY_COLUMNS].copy()
