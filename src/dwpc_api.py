"""
DWPC computation via Hetionet API.

This module provides functions for computing Degree-Weighted Path Counts (DWPC)
using the connectivity-search-backend Docker API. It includes health checking,
batch processing with checkpointing, and validation utilities.
"""

import asyncio
import io
import json
import logging
import os
import re
import shutil
import time
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import httpx
import pandas as pd
from tqdm import tqdm

from hetionet_utils_patched.udf import async_get_complete_metapaths


@contextmanager
def quiet_http_logs(names=("httpx", "urllib3", "chardet"), level=logging.WARNING):
    """
    Context manager to temporarily suppress HTTP client logging.

    Parameters
    ----------
    names : tuple of str
        Logger names to suppress
    level : int
        Logging level to set during suppression

    Yields
    ------
    None
    """
    prev = {}
    try:
        for n in names:
            lg = logging.getLogger(n)
            prev[n] = lg.level
            lg.setLevel(level)
        yield
    finally:
        for n, lvl in prev.items():
            logging.getLogger(n).setLevel(lvl)


@contextmanager
def capture_console():
    """
    Context manager to capture stdout and stderr.

    Yields
    ------
    io.StringIO
        Buffer containing captured output
    """
    buf = io.StringIO()
    with redirect_stdout(buf), redirect_stderr(buf):
        yield buf


def _now_tag() -> str:
    """Generate a timestamp string for file naming."""
    return time.strftime("%Y%m%d_%H%M%S")


async def check_neo4j_health(max_wait_minutes: int = 10) -> bool:
    """
    Check if Neo4j API is healthy.

    Parameters
    ----------
    max_wait_minutes : int
        Maximum time to wait for Neo4j to become healthy

    Returns
    -------
    bool
        True if Neo4j is healthy, False otherwise
    """
    api_url = os.environ.get('CONNECTIVITY_SEARCH_API', 'http://localhost:8015')
    health_url = f"{api_url}/v1/nodes/"

    print(f"Checking Neo4j health at {health_url}...")

    async with httpx.AsyncClient(timeout=5.0) as client:
        start_time = time.time()
        max_wait_seconds = max_wait_minutes * 60
        attempt = 0

        while True:
            attempt += 1
            try:
                response = await client.get(health_url)
                if response.status_code == 200:
                    print("  Neo4j is healthy!")
                    return True
            except (httpx.ConnectError, httpx.TimeoutException):
                pass

            elapsed = time.time() - start_time
            if elapsed > max_wait_seconds:
                print(f"  Neo4j health check timeout after {max_wait_minutes} minutes")
                return False

            if attempt % 6 == 0:
                print(f"  Waiting for Neo4j... ({int(elapsed)}s elapsed)")

            await asyncio.sleep(10)


def _build_unique_pairs(
    df_pairs: pd.DataFrame,
    col_source: str,
    col_target: str,
) -> List[Tuple[int, int, int]]:
    """
    Build list of unique source-target pairs from dataframe.

    Parameters
    ----------
    df_pairs : pd.DataFrame
        DataFrame containing source and target columns
    col_source : str
        Name of source column
    col_target : str
        Name of target column

    Returns
    -------
    List[Tuple[int, int, int]]
        List of (pair_index, source_id, target_id) tuples
    """
    pairs_df = (
        df_pairs[[col_source, col_target]]
        .dropna()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    return [(i + 1, int(r[col_source]), int(r[col_target]))
            for i, r in pairs_df.iterrows()]


def _prepare_out_dir(base_out_dir: Path, group: str, clear_group: bool) -> Path:
    """
    Prepare output directory for a processing group.

    Parameters
    ----------
    base_out_dir : Path
        Base output directory
    group : str
        Group name (subdirectory)
    clear_group : bool
        Whether to delete existing group directory

    Returns
    -------
    Path
        Path to the group directory
    """
    group_dir = base_out_dir / group
    if clear_group and group_dir.exists():
        print(f"Deleting existing group folder: {group_dir}")
        shutil.rmtree(group_dir)
    group_dir.mkdir(parents=True, exist_ok=True)
    return group_dir


def _result_path(group_dir: Path, pair_idx: int,
                 source_id: int, target_id: int) -> Path:
    """Generate parquet result file path."""
    return group_dir / f"s{source_id}_t{target_id}_pair{pair_idx:04d}.parquet"


def _meta_path(group_dir: Path, pair_idx: int,
               source_id: int, target_id: int) -> Path:
    """Generate metadata JSON file path."""
    return group_dir / f"s{source_id}_t{target_id}_pair{pair_idx:04d}.json"


def _get_existing_pairs(group_dir: Path) -> Set[Tuple[int, int]]:
    """
    Get set of already processed pairs for checkpointing.

    Parameters
    ----------
    group_dir : Path
        Directory containing parquet files

    Returns
    -------
    Set[Tuple[int, int]]
        Set of (source_id, target_id) tuples already processed
    """
    existing = set()
    if not group_dir.exists():
        return existing

    for parquet_file in group_dir.glob("s*_t*_pair*.parquet"):
        try:
            parts = parquet_file.stem.split('_')
            source_id = int(parts[0][1:])
            target_id = int(parts[1][1:])
            existing.add((source_id, target_id))
        except (ValueError, IndexError):
            continue

    return existing


def check_early_validation(summaries: List[Dict[str, Any]],
                           threshold: int = 10) -> None:
    """
    Check if early pairs are returning valid data.

    Raises RuntimeError if all first N pairs return 0 rows.

    Parameters
    ----------
    summaries : List[Dict]
        List of summary dictionaries from completed pairs
    threshold : int
        Number of pairs to check (default 10)

    Raises
    ------
    RuntimeError
        If all first threshold pairs returned 0 rows
    """
    if len(summaries) < threshold:
        return

    first_n = summaries[:threshold]
    zero_count = sum(1 for s in first_n if s.get('rows', 0) == 0)

    if zero_count == threshold:
        print(f"\n{'='*80}")
        print("EARLY VALIDATION FAILED")
        print(f"{'='*80}")
        print(f"All first {threshold} pairs returned 0 rows.")
        print("This indicates a systematic failure (likely invalid Neo4j IDs).")
        print("\nSample failed pair:")
        sample = first_n[0]
        print(f"  Source ID: {sample.get('neo4j_source_id', 'N/A')}")
        print(f"  Target ID: {sample.get('neo4j_target_id', 'N/A')}")
        print(f"  Failed metapaths: {len(sample.get('failed_metapaths', []))}")
        print("\nRecommendation: Stop processing and investigate Neo4j ID mappings.")
        print(f"{'='*80}\n")
        raise RuntimeError(
            f"Early validation failed: All first {threshold} pairs returned 0 rows"
        )


def check_periodic_progress(summaries: List[Dict[str, Any]],
                            check_interval: int = 100,
                            failure_threshold: float = 0.5) -> None:
    """
    Periodically check if recent pairs are failing.

    Prints warning if more than threshold of recent pairs have 0 rows.

    Parameters
    ----------
    summaries : List[Dict]
        List of summary dictionaries from completed pairs
    check_interval : int
        Check every N pairs (default 100)
    failure_threshold : float
        Fraction of failures to trigger warning (default 0.5)
    """
    if len(summaries) % check_interval != 0 or len(summaries) < check_interval:
        return

    recent = summaries[-check_interval:]
    zero_count = sum(1 for s in recent if s.get('rows', 0) == 0)
    failure_rate = zero_count / len(recent)

    if failure_rate > failure_threshold:
        print(f"\n{'='*80}")
        print("PERIODIC PROGRESS WARNING")
        print(f"{'='*80}")
        print(f"Recent {check_interval} pairs: {zero_count} returned 0 rows "
              f"({failure_rate:.1%} failure rate)")
        print(f"This exceeds the {failure_threshold:.0%} threshold.")
        print("\nRecommendation: Consider stopping and investigating before continuing.")
        print(f"Total pairs processed so far: {len(summaries)}")
        print(f"{'='*80}\n")


def validate_parquet_files(group_dir: Path, expected_pairs: int,
                           min_completion_rate: float = 0.2
                           ) -> Tuple[int, int, float]:
    """
    Validate parquet files before merging.

    Raises RuntimeError if completion rate is below threshold.

    Parameters
    ----------
    group_dir : Path
        Directory containing parquet files
    expected_pairs : int
        Expected number of GO-gene pairs
    min_completion_rate : float
        Minimum acceptable completion rate (default 0.2 = 20%)

    Returns
    -------
    Tuple[int, int, float]
        (n_valid_files, n_empty_files, completion_rate)

    Raises
    ------
    FileNotFoundError
        If no metadata files found
    RuntimeError
        If completion rate is below minimum threshold
    """
    metadata_files = list(group_dir.glob("*.json"))

    if not metadata_files:
        raise FileNotFoundError(f"No metadata files found in {group_dir}")

    n_valid = 0
    n_empty = 0

    for meta_file in metadata_files:
        try:
            with open(meta_file, 'r') as f:
                meta = json.load(f)
                if meta.get('rows', 0) > 0:
                    n_valid += 1
                else:
                    n_empty += 1
        except Exception:
            continue

    total_files = n_valid + n_empty
    completion_rate = n_valid / total_files if total_files > 0 else 0

    if completion_rate < min_completion_rate:
        print(f"\n{'='*80}")
        print("FINAL VALIDATION FAILED")
        print(f"{'='*80}")
        print("Parquet file validation:")
        print(f"  Valid files (rows > 0): {n_valid}")
        print(f"  Empty files (rows = 0): {n_empty}")
        print(f"  Completion rate: {completion_rate:.1%}")
        print(f"  Minimum required: {min_completion_rate:.0%}")
        print("\nThis indicates systematic failure during DWPC computation.")
        print("Do not merge these files into CSV.")
        print("\nRecommendation: Investigate failed pairs before proceeding.")
        print(f"{'='*80}\n")
        raise RuntimeError(
            f"Dataset has insufficient valid data: "
            f"{completion_rate:.1%} < {min_completion_rate:.0%}"
        )

    return n_valid, n_empty, completion_rate


async def _one_call(
    task_id: int,
    pair_idx: int,
    source_id: int,
    target_id: int,
    source_col_name: str,
    target_col_name: str,
    group_dir: Path,
    client: httpx.AsyncClient,
    sem: asyncio.Semaphore,
    compression: str,
    retries: int,
    backoff_first: float,
) -> Dict[str, Any]:
    """
    Fetch DWPC for a single GO-gene pair using shared client and semaphore.

    Parameters
    ----------
    task_id : int
        Unique task identifier
    pair_idx : int
        Pair index in original dataframe
    source_id : int
        Source node ID (GO term)
    target_id : int
        Target node ID (Gene)
    source_col_name : str
        Name of source column
    target_col_name : str
        Name of target column
    group_dir : Path
        Output directory for parquet files
    client : httpx.AsyncClient
        Shared HTTP client with connection pooling
    sem : asyncio.Semaphore
        Shared semaphore controlling global concurrency
    compression : str
        Parquet compression algorithm
    retries : int
        Number of retry attempts
    backoff_first : float
        Initial backoff delay in seconds

    Returns
    -------
    Dict[str, Any]
        Summary dictionary with results metadata
    """
    start_perf = time.perf_counter()
    summary: Dict[str, Any] = {
        "task_id": task_id,
        "pair_idx": pair_idx,
        source_col_name: source_id,
        target_col_name: target_id,
        "rows": None,
        "total_sec": None,
        "status": "ok",
        "file": None,
        "meta_file": None,
        "error": None,
        "failed_metapaths": [],
        "group_dir": str(group_dir),
    }

    tries = max(1, int(retries))
    delay = max(0.0, float(backoff_first))

    failed_metapaths: List[str] = []
    captured_tail: List[str] = []

    for attempt in range(1, tries + 1):
        try:
            with quiet_http_logs(), capture_console() as cap:
                res = await async_get_complete_metapaths(
                    source_id, target_id, client, sem
                )

            text_dump = cap.getvalue()
            if text_dump:
                failed_metapaths = re.findall(r"metapath='([^']+)'", text_dump)
                captured_tail = text_dump.splitlines()[-50:]

            df = pd.DataFrame(res) if res is not None else pd.DataFrame()
            df.insert(0, "record_idx", range(1, len(df) + 1))
            df.insert(1, "pair_idx", pair_idx)
            df.insert(2, source_col_name, source_id)
            df.insert(3, target_col_name, target_id)
            df.insert(4, "total_sec", None)

            out_fp = _result_path(group_dir, pair_idx, source_id, target_id)
            tmp_fp = out_fp.with_suffix(".parquet.tmp")
            df.to_parquet(tmp_fp, engine="pyarrow", index=False,
                          compression=compression)
            os.replace(tmp_fp, out_fp)

            total_sec = time.perf_counter() - start_perf

            if len(df) > 0:
                df.loc[:, "total_sec"] = total_sec
                tmp_fp2 = out_fp.with_suffix(".parquet.tmp2")
                df.to_parquet(tmp_fp2, engine="pyarrow", index=False,
                              compression=compression)
                os.replace(tmp_fp2, out_fp)

            meta = {
                "task_id": task_id,
                "pair_idx": pair_idx,
                source_col_name: source_id,
                target_col_name: target_id,
                "rows": int(len(df)),
                "total_sec": float(total_sec),
                "parquet_file": str(out_fp),
                "compression": compression,
                "failed_metapaths": failed_metapaths,
            }
            if failed_metapaths:
                meta["captured_log_tail"] = captured_tail

            meta_fp = _meta_path(group_dir, pair_idx, source_id, target_id)
            with open(meta_fp, "w") as f:
                json.dump(meta, f, indent=2)

            summary.update({
                "rows": len(df),
                "total_sec": total_sec,
                "file": str(out_fp),
                "meta_file": str(meta_fp),
                "failed_metapaths": failed_metapaths,
                "status": "partial" if failed_metapaths else "ok",
            })
            break

        except Exception as e:
            if attempt == tries:
                summary["status"] = "error"
                summary["error"] = str(e)
                summary["total_sec"] = time.perf_counter() - start_perf
            else:
                await asyncio.sleep(delay)
                delay = min(delay * 2, 60.0) if delay > 0 else 10.0

    return summary


async def run_metapaths_for_df(
    df_pairs: pd.DataFrame,
    *,
    col_source: str,
    col_target: str,
    base_out_dir: Path = Path("../output/dwpc_com/metapaths_parquet"),
    group: Optional[str] = None,
    clear_group: bool = False,
    max_concurrency: int = 120,
    compression: str = "zstd",
    retries: int = 10,
    backoff_first: float = 10.0,
    show_progress: bool = True,
    check_health: bool = True,
    enable_early_validation: bool = True,
    enable_periodic_monitoring: bool = True,
    early_check_threshold: int = 10,
    periodic_check_interval: int = 100,
    periodic_failure_threshold: float = 0.5,
) -> pd.DataFrame:
    """
    Process all GO-gene pairs and fetch DWPC via Het.io API.

    Creates a single shared httpx client and semaphore to avoid nested
    concurrency bottleneck. All HTTP requests are controlled by one global
    semaphore for optimal throughput.

    Parameters
    ----------
    df_pairs : pd.DataFrame
        DataFrame with GO-gene pairs
    col_source : str
        Source column name (GO term neo4j_id)
    col_target : str
        Target column name (Gene neo4j_id)
    base_out_dir : Path
        Base output directory
    group : str, optional
        Group name for this run
    clear_group : bool
        Whether to delete existing group folder
    max_concurrency : int
        Maximum concurrent HTTP requests (default 120)
    compression : str
        Parquet compression algorithm
    retries : int
        Number of retry attempts per pair
    backoff_first : float
        Initial backoff delay
    show_progress : bool
        Show progress bar
    check_health : bool
        Check Neo4j health before starting
    enable_early_validation : bool
        Enable early validation check
    enable_periodic_monitoring : bool
        Enable periodic progress monitoring
    early_check_threshold : int
        Number of pairs for early validation
    periodic_check_interval : int
        Check interval for periodic monitoring
    periodic_failure_threshold : float
        Failure rate threshold for warnings

    Returns
    -------
    pd.DataFrame
        Summary dataframe with processing results
    """
    if check_health:
        is_healthy = await check_neo4j_health(max_wait_minutes=10)
        if not is_healthy:
            raise RuntimeError("Neo4j is not healthy, cannot proceed")

    group = group or f"run_{_now_tag()}"
    group_dir = _prepare_out_dir(base_out_dir, group, clear_group)

    pairs = _build_unique_pairs(df_pairs, col_source, col_target)
    if not pairs:
        print("No valid pairs found.")
        return pd.DataFrame()

    existing_pairs = _get_existing_pairs(group_dir)

    if existing_pairs:
        print(f"CHECKPOINT: Found {len(existing_pairs)} existing pair results")
        pairs_to_process = [
            (idx, s, t) for idx, s, t in pairs if (s, t) not in existing_pairs
        ]
        skipped_count = len(pairs) - len(pairs_to_process)
        print(f"  Skipping {skipped_count} completed pairs")
        print(f"  Processing {len(pairs_to_process)} remaining pairs")
    else:
        pairs_to_process = pairs

    if not pairs_to_process:
        print("All pairs already processed!")
        return pd.DataFrame()

    limits = httpx.Limits(
        max_keepalive_connections=200,
        max_connections=250,
        keepalive_expiry=60.0
    )

    timeout = httpx.Timeout(
        timeout=90.0,
        connect=10.0,
        read=70.0,
        write=10.0
    )

    print("Connection pooling configuration:")
    print(f"  Max connections: {limits.max_connections}")
    print(f"  Keepalive connections: {limits.max_keepalive_connections}")
    print(f"  Global concurrency: {max_concurrency}")
    print(f"  Timeout: {90.0}s")

    async with httpx.AsyncClient(limits=limits, timeout=timeout) as shared_client:
        sem = asyncio.Semaphore(max_concurrency)

        coroutines = []
        task_id = 0
        for pair_idx, s, t in pairs_to_process:
            task_id += 1
            coroutines.append(
                _one_call(
                    task_id, pair_idx, s, t,
                    source_col_name=col_source,
                    target_col_name=col_target,
                    group_dir=group_dir,
                    client=shared_client,
                    sem=sem,
                    compression=compression,
                    retries=retries,
                    backoff_first=backoff_first,
                )
            )

        start = time.perf_counter()

        if show_progress:
            summaries = []
            pbar_desc = f"Processing [{group}]"
            with tqdm(total=len(coroutines), desc=pbar_desc, unit="task") as pbar:
                chunk_size = max_concurrency
                for i in range(0, len(coroutines), chunk_size):
                    chunk = coroutines[i:i+chunk_size]
                    chunk_results = await asyncio.gather(
                        *chunk, return_exceptions=False
                    )
                    summaries.extend(chunk_results)
                    pbar.update(len(chunk))

                    if enable_early_validation:
                        try:
                            check_early_validation(
                                summaries, threshold=early_check_threshold
                            )
                        except RuntimeError:
                            print("Stopping processing due to early validation failure")
                            raise

                    if enable_periodic_monitoring:
                        check_periodic_progress(
                            summaries,
                            check_interval=periodic_check_interval,
                            failure_threshold=periodic_failure_threshold
                        )
        else:
            summaries = await asyncio.gather(*coroutines, return_exceptions=False)

    total_sec = time.perf_counter() - start
    print(f"\n[{group}] Finished {len(coroutines)} unique pairs in {total_sec:.2f}s")

    summary_df = pd.DataFrame(summaries)
    cols_front = [
        "task_id", "pair_idx", col_source, col_target, "rows", "total_sec",
        "status", "file", "meta_file", "error", "failed_metapaths", "group_dir"
    ]
    cols_front = [c for c in cols_front if c in summary_df.columns]
    return summary_df[cols_front]
