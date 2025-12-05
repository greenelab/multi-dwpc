"""
Creates UDFs for various Het.io API calls.

PATCHED VERSION: Added URL encoding for metapath parameters to fix failures
with directional metapaths containing < and > characters.
"""

import asyncio
import json
import os
from urllib.parse import quote

import httpx
import requests

CONNECTIVITY_SEARCH_API = os.environ.get(
    "CONNECTIVITY_SEARCH_API", "https://search-api.het.io"
)


def get_paths_json(source: int, target: int, metapath: str) -> str:
    """
    Fetch the full Het.io "paths" JSON blob for the given triple
    and return it as a raw JSON string.
    """
    url = (
        f"{CONNECTIVITY_SEARCH_API}/v1/paths/"
        f"source/{source}/target/{target}/metapath/{quote(metapath, safe='')}/"
        "?format=json"
    )
    resp = requests.get(url)
    resp.raise_for_status()
    raw_paths = resp.json()["paths"]

    EXPECTED_KEYS = [
        "metapath",
        "node_ids",
        "rel_ids",
        "percent_of_DWPC",
        "PC",
        "DWPC",
        "score",
        "PDP",
    ]

    return json.dumps([{k: p.get(k) for k in EXPECTED_KEYS} for p in raw_paths])


async def fetch_metapaths(
    client: httpx.AsyncClient,
    source: int,
    target: int,
    semaphore: asyncio.Semaphore,
) -> list:
    """
    Fetch list of available metapaths for a source-target pair.

    Parameters
    ----------
    client : httpx.AsyncClient
        Shared HTTP client with connection pooling
    source : int
        Source node ID
    target : int
        Target node ID
    semaphore : asyncio.Semaphore
        Shared semaphore to control concurrency

    Returns
    -------
    list
        List of metapath IDs
    """
    url = (
        f"{CONNECTIVITY_SEARCH_API}/v1/metapaths/"
        f"source/{source}/target/{target}/"
        "?format=json&complete=true"
    )
    async with semaphore:
        resp = await client.get(url, timeout=30.0)
        resp.raise_for_status()
        data = resp.json()["path_counts"]
    return [item["metapath_id"] for item in data]


async def fetch_path_count(
    client: httpx.AsyncClient,
    source: int,
    target: int,
    metapath: str,
    semaphore: asyncio.Semaphore,
) -> dict:
    """
    Fetch path count information for a specific metapath.

    Parameters
    ----------
    client : httpx.AsyncClient
        Shared HTTP client with connection pooling
    source : int
        Source node ID
    target : int
        Target node ID
    metapath : str
        Metapath abbreviation
    semaphore : asyncio.Semaphore
        Shared semaphore to control concurrency

    Returns
    -------
    dict
        Path count information with DWPC statistics
    """
    url = (
        f"{CONNECTIVITY_SEARCH_API}/v1/paths/"
        f"source/{source}/target/{target}/metapath/{quote(metapath, safe='')}/"
        "?format=json"
    )
    async with semaphore:
        resp = await client.get(url, timeout=30.0)
        resp.raise_for_status()
        raw = resp.json()["path_count_info"]

    EXPECTED_KEYS = [
        "source",
        "target",
        "metapath_abbreviation",
        "path_count",
        "adjusted_p_value",
        "p_value",
        "dwpc",
        "dgp_source_degree",
        "dgp_target_degree",
        "dgp_n_dwpcs",
        "dgp_n_nonzero_dwpcs",
        "dgp_nonzero_mean",
        "dgp_nonzero_sd",
    ]

    return {k: raw.get(k) for k in EXPECTED_KEYS}


async def async_get_complete_metapaths(
    source: int,
    target: int,
    client: httpx.AsyncClient,
    semaphore: asyncio.Semaphore,
) -> list:
    """
    Fetch all complete metapaths concurrently via asyncio+httpx.

    OPTIMIZED VERSION: Uses shared httpx client and semaphore to avoid
    nested concurrency bottleneck. All HTTP requests are controlled by
    a single global semaphore for optimal throughput.

    Parameters
    ----------
    source : int
        Source node ID
    target : int
        Target node ID
    client : httpx.AsyncClient
        Shared HTTP client with connection pooling
    semaphore : asyncio.Semaphore
        Shared semaphore controlling all HTTP requests globally

    Returns
    -------
    list
        List of path count dictionaries with DWPC statistics
    """
    metapaths = await fetch_metapaths(client, source, target, semaphore)

    async def worker(metapath_id: str) -> dict:
        return await fetch_path_count(client, source, target, metapath_id, semaphore)

    tasks = [asyncio.create_task(worker(m)) for m in metapaths]
    results = await asyncio.gather(*tasks, return_exceptions=True)

    output, errors = [], []
    for mid, res in zip(metapaths, results):
        if isinstance(res, Exception):
            errors.append((mid, res))
        else:
            output.append(res)

    if errors:
        for mid, exc in errors:
            print(f"[ERROR] metapath={mid!r} → {exc!r}")

    return output


def get_complete_metapaths(source: str, target: str) -> list:
    return asyncio.run(async_get_complete_metapaths(source, target))
