"""
Direct DWPC computation using hetmatpy.

This module provides functions for computing Degree-Weighted Path Counts (DWPC)
directly from the HetMat sparse matrix files, without requiring the Docker API.
This approach is faster for batch computation and does not require network calls.

The DWPC algorithm is implemented by hetmatpy which correctly handles:
1. Degree weighting with configurable damping parameter
2. Path counting (excluding walks that revisit the same node)
3. Various metapath categories (no_repeats, short_repeat, disjoint, etc.)

Optimizations implemented:
1. Single matrix computation with reuse across all datasets
2. Disk caching of DWPC matrices for future runs
3. Sparse output storage (only non-zero DWPC values)
4. Parallel computation of metapath matrices

The API stores transformed DWPC: arcsinh(raw_dwpc / dwpc_raw_mean)
This module applies the same transformation using metapath statistics from:
https://github.com/greenelab/hetmech/blob/main/explore/bulk-pipeline/archives/metapath-dwpc-stats.tsv

References:
- https://het.io/hnep/
- https://github.com/hetio/hetmatpy
- https://pmc.ncbi.nlm.nih.gov/articles/PMC10375517/
"""

import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm

# Import hetmatpy for correct DWPC computation
from hetmatpy.hetmat import HetMat as HetMatPy
from hetmatpy.degree_weight import dwpc as hetmatpy_dwpc


METAPATH_STATS_URL = (
    "https://raw.githubusercontent.com/greenelab/hetmech/"
    "34e95b9f72f47cdeba3d51622bee31f79e9a4cb8/explore/bulk-pipeline/archives/"
    "metapath-dwpc-stats.tsv"
)

# Default cache directory for DWPC matrices
DEFAULT_CACHE_DIR = "dwpc_cache"


def load_metagraph(data_dir: Path) -> Dict:
    """
    Load metagraph configuration from JSON.

    Parameters
    ----------
    data_dir : Path
        Directory containing metagraph.json

    Returns
    -------
    Dict
        Metagraph configuration with metanode_kinds, metaedge_tuples, kind_to_abbrev
    """
    metagraph_path = data_dir / "metagraph.json"
    with open(metagraph_path) as f:
        return json.load(f)


def load_metapath_stats(data_dir: Path) -> pd.DataFrame:
    """
    Load metapath DWPC statistics for arcsinh transformation.

    Downloads from GitHub if not cached locally.

    Parameters
    ----------
    data_dir : Path
        Directory to cache the statistics file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: metapath, dwpc_raw_mean, etc.
    """
    cache_path = data_dir / "metapath-dwpc-stats.tsv"

    if not cache_path.exists():
        import urllib.request
        urllib.request.urlretrieve(METAPATH_STATS_URL, cache_path)

    df = pd.read_csv(cache_path, sep="\t")
    df = df.rename(columns={"dwpc-0.5_raw_mean": "dwpc_raw_mean"})
    return df


def get_dwpc_raw_mean(metapath_stats: pd.DataFrame, metapath: str) -> float:
    """
    Get the dwpc_raw_mean for a metapath.

    Handles both forward and reverse metapaths by checking both orientations.
    The statistics file uses canonical ordering (e.g., BPpG not GpBP).

    Parameters
    ----------
    metapath_stats : pd.DataFrame
        Metapath statistics from load_metapath_stats()
    metapath : str
        Metapath abbreviation

    Returns
    -------
    float
        The dwpc_raw_mean value for transformation
    """
    row = metapath_stats[metapath_stats["metapath"] == metapath]
    if len(row) > 0:
        return row["dwpc_raw_mean"].iloc[0]

    reverse_metapath = reverse_metapath_abbrev(metapath)
    row = metapath_stats[metapath_stats["metapath"] == reverse_metapath]
    if len(row) > 0:
        return row["dwpc_raw_mean"].iloc[0]

    raise ValueError(f"Metapath {metapath} (or reverse {reverse_metapath}) not found")


def reverse_metapath_abbrev(metapath: str) -> str:
    """
    Reverse a metapath abbreviation.

    For example: GpBP -> BPpG, GiGpBP -> BPpGiG

    Parameters
    ----------
    metapath : str
        Metapath abbreviation

    Returns
    -------
    str
        Reversed metapath abbreviation
    """
    edge_pattern = []
    pos = 0
    node_abbrevs = {"G", "BP", "CC", "MF", "PW", "A", "D", "C", "SE", "S", "PC"}
    edge_abbrevs = {"p", "i", "c", "r", ">", "<", "a", "d", "u", "e", "b", "t", "l"}

    tokens = []
    while pos < len(metapath):
        if pos + 2 <= len(metapath) and metapath[pos:pos+2] in node_abbrevs:
            tokens.append(metapath[pos:pos+2])
            pos += 2
        elif metapath[pos] in node_abbrevs:
            tokens.append(metapath[pos])
            pos += 1
        elif metapath[pos] in edge_abbrevs:
            tokens.append(metapath[pos])
            pos += 1
        elif metapath[pos] == ">":
            tokens.append(">")
            pos += 1
        elif metapath[pos] == "<":
            tokens.append("<")
            pos += 1
        else:
            pos += 1

    direction_map = {">": "<", "<": ">"}
    reversed_tokens = []
    for token in reversed(tokens):
        if token in direction_map:
            reversed_tokens.append(direction_map[token])
        else:
            reversed_tokens.append(token)

    return "".join(reversed_tokens)


def transform_dwpc(raw_dwpc: np.ndarray, dwpc_raw_mean: float) -> np.ndarray:
    """
    Apply arcsinh transformation to match API DWPC values.

    The connectivity-search-backend stores DWPC as:
        dwpc = arcsinh(raw_dwpc / dwpc_raw_mean)

    Parameters
    ----------
    raw_dwpc : np.ndarray
        Raw DWPC values from matrix computation
    dwpc_raw_mean : float
        Metapath-specific mean DWPC for normalization

    Returns
    -------
    np.ndarray
        Transformed DWPC values matching API format
    """
    return np.arcsinh(raw_dwpc / dwpc_raw_mean)


def load_node_ids(data_dir: Path, node_type: str) -> pd.DataFrame:
    """
    Load node identifiers for a given node type.

    Parameters
    ----------
    data_dir : Path
        Directory containing nodes/ subdirectory
    node_type : str
        Node type name (e.g., "Gene", "Biological Process")

    Returns
    -------
    pd.DataFrame
        DataFrame with node identifiers and names
    """
    node_path = data_dir / "nodes" / f"{node_type}.tsv"
    return pd.read_csv(node_path, sep="\t")


def load_adjacency_matrix(data_dir: Path, metaedge_abbrev: str) -> sparse.csr_matrix:
    """
    Load sparse adjacency matrix for a metaedge.

    Parameters
    ----------
    data_dir : Path
        Directory containing edges/ subdirectory
    metaedge_abbrev : str
        Metaedge abbreviation (e.g., "GpBP", "GiG")

    Returns
    -------
    sparse.csr_matrix
        Sparse adjacency matrix
    """
    matrix_path = data_dir / "edges" / f"{metaedge_abbrev}.sparse.npz"
    return sparse.load_npz(matrix_path).tocsr()


def degree_weight_matrix(
    matrix: sparse.csr_matrix,
    damping: float = 0.4,
    dtype: np.dtype = np.float64
) -> sparse.csr_matrix:
    """
    Apply degree weighting to an adjacency matrix for undirected graphs.

    For undirected graphs, each node's degree is the total number of edges
    incident to it. The weighting normalizes by the degree of both endpoints.

    DWPC_weighted[i,j] = A[i,j] / (degree[i]^damping * degree[j]^damping)

    For a symmetric adjacency matrix, degree[i] = row_sum[i] = col_sum[i].
    For asymmetric storage, we symmetrize by computing degree as the sum
    of the matrix and its transpose.

    Parameters
    ----------
    matrix : sparse.csr_matrix
        Input adjacency matrix
    damping : float
        Damping exponent (default 0.4). Controls downweighting of paths
        through high-degree nodes. 0 = no damping, 1 = full normalization.
    dtype : np.dtype
        Output data type

    Returns
    -------
    sparse.csr_matrix
        Degree-weighted matrix
    """
    matrix = matrix.astype(dtype)

    # For undirected graphs, compute total degree for each node
    # Row degree = sum of outgoing edges from source nodes
    # Col degree = sum of incoming edges to target nodes
    # These represent the total degree of source and target node types
    row_degree = np.array(matrix.sum(axis=1), dtype=dtype).flatten()
    col_degree = np.array(matrix.sum(axis=0), dtype=dtype).flatten()

    row_degree[row_degree == 0] = 1.0
    col_degree[col_degree == 0] = 1.0

    row_weights = np.power(row_degree, -damping)
    col_weights = np.power(col_degree, -damping)

    row_diag = sparse.diags(row_weights, format="csr")
    col_diag = sparse.diags(col_weights, format="csr")

    weighted = row_diag @ matrix @ col_diag

    return weighted.tocsr()


def parse_metapath(metapath_abbrev: str, metagraph: Dict) -> List[Dict]:
    """
    Parse a metapath abbreviation into edge components.

    Parameters
    ----------
    metapath_abbrev : str
        Metapath abbreviation (e.g., "GpBPpG", "GiGpBP")
    metagraph : Dict
        Metagraph configuration

    Returns
    -------
    List[Dict]
        List of edge dictionaries with source, target, edge_type, direction
    """
    abbrev_to_kind = {v: k for k, v in metagraph["kind_to_abbrev"].items()}

    metaedge_info = {}
    for edge_tuple in metagraph["metaedge_tuples"]:
        source, target, edge_type, direction = edge_tuple
        src_abbrev = metagraph["kind_to_abbrev"][source]
        tgt_abbrev = metagraph["kind_to_abbrev"][target]
        edge_abbrev = metagraph["kind_to_abbrev"][edge_type]

        forward_key = f"{src_abbrev}{edge_abbrev}{tgt_abbrev}"
        metaedge_info[forward_key] = {
            "source": source,
            "target": target,
            "edge_type": edge_type,
            "direction": "forward",
            "abbrev": forward_key
        }

        if direction == "both":
            reverse_key = f"{tgt_abbrev}{edge_abbrev}{src_abbrev}"
            metaedge_info[reverse_key] = {
                "source": target,
                "target": source,
                "edge_type": edge_type,
                "direction": "reverse",
                "abbrev": forward_key
            }

    edges = []
    pos = 0
    while pos < len(metapath_abbrev):
        found = False
        for length in [4, 3]:
            if pos + length <= len(metapath_abbrev):
                candidate = metapath_abbrev[pos:pos + length]
                if candidate in metaedge_info:
                    edges.append(metaedge_info[candidate])
                    pos += length - 1
                    found = True
                    break
        if not found:
            pos += 1

    return edges


def get_metaedge_abbreviations(metagraph: Dict) -> List[str]:
    """
    Get all valid metaedge abbreviations from metagraph.

    Parameters
    ----------
    metagraph : Dict
        Metagraph configuration

    Returns
    -------
    List[str]
        List of metaedge abbreviations (e.g., ["GpBP", "GiG", ...])
    """
    abbreviations = []
    for edge_tuple in metagraph["metaedge_tuples"]:
        source, target, edge_type, _ = edge_tuple
        src_abbrev = metagraph["kind_to_abbrev"][source]
        tgt_abbrev = metagraph["kind_to_abbrev"][target]
        edge_abbrev = metagraph["kind_to_abbrev"][edge_type]
        abbreviations.append(f"{src_abbrev}{edge_abbrev}{tgt_abbrev}")
    return abbreviations


def _compute_dwpc_matrix_worker(args: Tuple[str, str, float]) -> Tuple[str, sparse.csr_matrix]:
    """
    Compute one DWPC matrix in a subprocess.

    This worker is defined at module scope so it can be pickled by
    ProcessPoolExecutor. Each process creates its own hetmatpy object rather
    than trying to serialize a shared instance from the parent process.
    """
    data_dir, metapath, damping = args
    hetmatpy = HetMatPy(data_dir)
    mp_obj = hetmatpy.metagraph.metapath_from_abbrev(metapath)
    _, _, matrix = hetmatpy_dwpc(hetmatpy, mp_obj, damping=damping)

    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    return metapath, matrix


class HetMat:
    """
    HetMat data structure for efficient DWPC computation.

    Uses hetmatpy internally for correct DWPC computation, which properly
    handles path counting (excluding walks that revisit the same node).

    Features:
    - In-memory caching: matrices computed once are reused across datasets
    - Disk caching: matrices saved to disk for future runs
    - Sparse storage: only non-zero values stored to reduce memory/disk usage

    DWPC values are transformed using arcsinh(raw / mean) to match API values.
    """

    def __init__(
        self,
        data_dir: Union[str, Path],
        damping: float = 0.5,
        cache_dir: Optional[Union[str, Path]] = None,
        use_disk_cache: bool = True,
        write_disk_cache: bool = True,
    ):
        """
        Initialize HetMat from data directory.

        Parameters
        ----------
        data_dir : str or Path
            Directory containing metagraph.json, nodes/, edges/
        damping : float
            Default damping exponent for DWPC computation (0.5 matches API)
        cache_dir : str or Path, optional
            Directory for disk caching. Defaults to data_dir/dwpc_cache
        use_disk_cache : bool
            Whether to use disk caching for DWPC matrices
        write_disk_cache : bool
            Whether to write newly computed matrices to disk cache. When False,
            existing disk cache entries are still read if available.
        """
        self.data_dir = Path(data_dir)
        self.damping = damping
        self.use_disk_cache = use_disk_cache
        self.write_disk_cache = bool(write_disk_cache and use_disk_cache)

        # Set up cache directory
        if cache_dir is None:
            self.cache_dir = self.data_dir / DEFAULT_CACHE_DIR
        else:
            self.cache_dir = Path(cache_dir)

        if self.use_disk_cache:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Load our JSON metagraph for compatibility
        self.metagraph = load_metagraph(self.data_dir)
        self.metapath_stats = load_metapath_stats(self.data_dir)

        # Initialize hetmatpy's HetMat for DWPC computation
        self._hetmatpy = HetMatPy(self.data_dir)

        self._node_dfs: Dict[str, pd.DataFrame] = {}
        self._dwpc_cache: Dict[Tuple[str, float], sparse.csr_matrix] = {}
        self._dwpc_cache_csc: Dict[Tuple[str, float], sparse.csc_matrix] = {}

    def _get_cache_path(self, metapath: str, damping: float) -> Path:
        """Get the disk cache path for a DWPC matrix."""
        return self.cache_dir / f"dwpc_{metapath}_d{damping:.2f}.npz"

    def _load_from_disk(self, metapath: str, damping: float) -> Optional[sparse.csr_matrix]:
        """Load DWPC matrix from disk cache if available."""
        if not self.use_disk_cache:
            return None

        cache_path = self._get_cache_path(metapath, damping)
        if cache_path.exists():
            return sparse.load_npz(cache_path)
        return None

    def _save_to_disk(self, metapath: str, damping: float, matrix: sparse.csr_matrix):
        """Save DWPC matrix to disk cache."""
        if not self.write_disk_cache:
            return

        cache_path = self._get_cache_path(metapath, damping)
        sparse.save_npz(cache_path, matrix)

    def get_nodes(self, node_type: str) -> pd.DataFrame:
        """Get node DataFrame, loading from disk if needed."""
        if node_type not in self._node_dfs:
            self._node_dfs[node_type] = load_node_ids(self.data_dir, node_type)
        return self._node_dfs[node_type]

    def compute_dwpc_matrix(
        self,
        metapath: str,
        damping: Optional[float] = None
    ) -> sparse.csr_matrix:
        """
        Compute full DWPC matrix for a metapath using hetmatpy.

        This correctly handles path counting (vs walk counting) by using
        hetmatpy's implementation which excludes paths that revisit nodes.

        Matrices are cached both in memory and on disk (as sparse matrices)
        for efficient reuse across datasets and pipeline runs.

        Parameters
        ----------
        metapath : str
            Metapath abbreviation (e.g., "GpBP", "GiGiGpBP")
        damping : float, optional
            Damping exponent (uses instance default if not specified)

        Returns
        -------
        sparse.csr_matrix
            Sparse DWPC matrix (source_nodes x target_nodes)
        """
        damping = damping if damping is not None else self.damping
        cache_key = (metapath, damping)

        # Check in-memory cache first
        if cache_key in self._dwpc_cache:
            return self._dwpc_cache[cache_key]

        # Check disk cache
        matrix = self._load_from_disk(metapath, damping)
        if matrix is not None:
            self._dwpc_cache[cache_key] = matrix
            return matrix

        # Compute DWPC using hetmatpy
        mp_obj = self._hetmatpy.metagraph.metapath_from_abbrev(metapath)
        rows, cols, matrix = hetmatpy_dwpc(
            self._hetmatpy, mp_obj, damping=damping
        )

        # Convert to sparse CSR format for efficient storage
        if not sparse.issparse(matrix):
            matrix = sparse.csr_matrix(matrix)
        else:
            matrix = matrix.tocsr()

        # Cache in memory and on disk
        self._dwpc_cache[cache_key] = matrix
        self._save_to_disk(metapath, damping, matrix)

        return matrix

    def compute_dwpc_matrix_csc(
        self,
        metapath: str,
        damping: Optional[float] = None
    ) -> sparse.csc_matrix:
        """
        Compute/load a DWPC matrix and return CSC format.

        CSC is optimized for repeated column slicing, which is the dominant
        access pattern in LV precompute.
        """
        damping = damping if damping is not None else self.damping
        cache_key = (metapath, damping)
        if cache_key in self._dwpc_cache_csc:
            return self._dwpc_cache_csc[cache_key]

        matrix = self.compute_dwpc_matrix(metapath=metapath, damping=damping)
        if sparse.isspmatrix_csc(matrix):
            csc = matrix
        else:
            csc = matrix.tocsc()
        self._dwpc_cache_csc[cache_key] = csc
        return csc

    def precompute_matrices(
        self,
        metapaths: List[str],
        damping: Optional[float] = None,
        n_workers: int = 4,
        show_progress: bool = True,
        cache_in_memory: bool = True,
    ):
        """
        Precompute and cache DWPC matrices for multiple metapaths in parallel.

        Parameters
        ----------
        metapaths : List[str]
            List of metapath abbreviations to precompute
        damping : float, optional
            Damping exponent
        n_workers : int
            Number of parallel workers
        show_progress : bool
            Show progress bar
        cache_in_memory : bool
            If False, completed matrices are written to disk and released instead
            of being retained in the in-memory cache. Use for disk-cache warmup
            jobs that have no follow-on computation and need bounded RAM.
        """
        damping = damping if damping is not None else self.damping

        # Filter to metapaths not already cached
        to_compute = []
        for mp in metapaths:
            cache_key = (mp, damping)
            if cache_key in self._dwpc_cache:
                continue
            disk_matrix = self._load_from_disk(mp, damping)
            if disk_matrix is not None:
                if cache_in_memory:
                    self._dwpc_cache[cache_key] = disk_matrix
                continue
            to_compute.append(mp)

        if not to_compute:
            return

        # Use a top-level worker so each subprocess creates its own hetmatpy
        # instance instead of attempting to pickle self._hetmatpy.
        args_list = [(str(self.data_dir), mp, damping) for mp in to_compute]

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_compute_dwpc_matrix_worker, args): args[1]
                for args in args_list
            }

            iterator = as_completed(futures)
            if show_progress:
                iterator = tqdm(iterator, total=len(to_compute), desc="Computing DWPC matrices")

            for future in iterator:
                mp, matrix = future.result()
                cache_key = (mp, damping)
                self._save_to_disk(mp, damping, matrix)
                if cache_in_memory:
                    self._dwpc_cache[cache_key] = matrix

    def get_dwpc_for_pairs(
        self,
        metapath: str,
        source_indices: np.ndarray,
        target_indices: np.ndarray,
        damping: Optional[float] = None,
        transform: bool = True
    ) -> np.ndarray:
        """
        Compute DWPC values for specific source-target pairs.

        Parameters
        ----------
        metapath : str
            Metapath abbreviation
        source_indices : np.ndarray
            Array of source node indices (row indices in matrix)
        target_indices : np.ndarray
            Array of target node indices (column indices in matrix)
        damping : float, optional
            Damping exponent
        transform : bool
            If True, apply arcsinh transformation to match API values

        Returns
        -------
        np.ndarray
            DWPC values for each source-target pair (transformed if transform=True)
        """
        dwpc_matrix = self.compute_dwpc_matrix(metapath, damping)

        source_indices = np.asarray(source_indices, dtype=np.int64)
        target_indices = np.asarray(target_indices, dtype=np.int64)

        valid_sources = (source_indices >= 0) & (source_indices < dwpc_matrix.shape[0])
        valid_targets = (target_indices >= 0) & (target_indices < dwpc_matrix.shape[1])
        valid_mask = valid_sources & valid_targets

        # Extract values from sparse matrix efficiently
        dwpc_values = np.zeros(len(source_indices), dtype=np.float64)
        valid_src = source_indices[valid_mask]
        valid_tgt = target_indices[valid_mask]

        # Use vectorized sparse matrix indexing
        dwpc_values[valid_mask] = np.asarray(
            dwpc_matrix[valid_src, valid_tgt]
        ).flatten()

        if transform:
            dwpc_raw_mean = get_dwpc_raw_mean(self.metapath_stats, metapath)
            dwpc_values = transform_dwpc(dwpc_values, dwpc_raw_mean)

        return dwpc_values

    def clear_memory_cache(self):
        """Clear the in-memory cache to free RAM."""
        self._dwpc_cache.clear()
        self._dwpc_cache_csc.clear()

    def clear_metapath_from_memory(
        self,
        metapath: str,
        damping: Optional[float] = None,
    ):
        """Clear one metapath matrix from in-memory caches."""
        damping = damping if damping is not None else self.damping
        cache_key = (metapath, damping)
        self._dwpc_cache.pop(cache_key, None)
        self._dwpc_cache_csc.pop(cache_key, None)

    def clear_disk_cache(self):
        """Clear the disk cache."""
        if self.cache_dir.exists():
            for f in self.cache_dir.glob("dwpc_*.npz"):
                f.unlink()


def compute_dwpc_for_dataframe(
    hetmat: HetMat,
    df: pd.DataFrame,
    metapaths: List[str],
    source_col: str = "source_idx",
    target_col: str = "target_idx",
    damping: float = 0.5,
    n_workers: int = 1,
    show_progress: bool = True,
    transform: bool = True
) -> pd.DataFrame:
    """
    Compute DWPC for all pairs in dataframe across multiple metapaths.

    Parameters
    ----------
    hetmat : HetMat
        Initialized HetMat object
    df : pd.DataFrame
        DataFrame with source and target index columns
    metapaths : List[str]
        List of metapath abbreviations to compute
    source_col : str
        Name of source index column
    target_col : str
        Name of target index column
    damping : float
        Damping exponent
    n_workers : int
        Number of parallel workers (1 = sequential)
    show_progress : bool
        Show progress bar
    transform : bool
        If True, apply arcsinh transformation to match API values

    Returns
    -------
    pd.DataFrame
        Original DataFrame with additional DWPC columns for each metapath
    """
    result_df = df.copy()
    source_indices = result_df[source_col].values
    target_indices = result_df[target_col].values

    iterator = tqdm(metapaths, desc="Computing DWPC", disable=not show_progress)

    for metapath in iterator:
        try:
            dwpc_values = hetmat.get_dwpc_for_pairs(
                metapath, source_indices, target_indices, damping, transform
            )
            result_df[f"dwpc_{metapath}"] = dwpc_values
        except Exception as e:
            print(f"Warning: Failed to compute DWPC for {metapath}: {e}")
            result_df[f"dwpc_{metapath}"] = np.nan

    return result_df


def _compute_single_metapath(args: Tuple) -> Tuple[str, np.ndarray]:
    """
    Worker function for parallel metapath computation.

    Parameters
    ----------
    args : Tuple
        (data_dir, metapath, source_indices, target_indices, damping, transform)

    Returns
    -------
    Tuple[str, np.ndarray]
        (metapath, dwpc_values)
    """
    data_dir, metapath, source_indices, target_indices, damping, transform = args
    hetmat = HetMat(data_dir, damping=damping)

    try:
        dwpc_values = hetmat.get_dwpc_for_pairs(
            metapath, source_indices, target_indices, damping, transform
        )
        return (metapath, dwpc_values)
    except Exception as e:
        print(f"Warning: Failed {metapath}: {e}")
        return (metapath, np.full(len(source_indices), np.nan))


def compute_dwpc_parallel(
    data_dir: Union[str, Path],
    df: pd.DataFrame,
    metapaths: List[str],
    source_col: str = "source_idx",
    target_col: str = "target_idx",
    damping: float = 0.5,
    n_workers: int = 4,
    show_progress: bool = True,
    transform: bool = True
) -> pd.DataFrame:
    """
    Compute DWPC in parallel across metapaths.

    Uses ProcessPoolExecutor to compute different metapaths concurrently.
    Each worker loads its own HetMat instance to avoid shared state issues.

    Parameters
    ----------
    data_dir : str or Path
        Path to HetMat data directory
    df : pd.DataFrame
        DataFrame with source and target index columns
    metapaths : List[str]
        List of metapath abbreviations
    source_col : str
        Name of source index column
    target_col : str
        Name of target index column
    damping : float
        Damping exponent (0.5 matches API)
    n_workers : int
        Number of parallel workers
    show_progress : bool
        Show progress bar
    transform : bool
        If True, apply arcsinh transformation to match API values

    Returns
    -------
    pd.DataFrame
        DataFrame with DWPC columns for each metapath
    """
    result_df = df.copy()
    source_indices = result_df[source_col].values
    target_indices = result_df[target_col].values

    args_list = [
        (str(data_dir), mp, source_indices, target_indices, damping, transform)
        for mp in metapaths
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_compute_single_metapath, args): args[1]
            for args in args_list
        }

        iterator = as_completed(futures)
        if show_progress:
            iterator = tqdm(iterator, total=len(metapaths), desc="Computing DWPC")

        for future in iterator:
            metapath, dwpc_values = future.result()
            result_df[f"dwpc_{metapath}"] = dwpc_values

    return result_df


def get_gene_bp_metapaths(metagraph: Dict, max_length: int = 4) -> List[str]:
    """
    Generate relevant metapaths between Gene and Biological Process.

    Parameters
    ----------
    metagraph : Dict
        Metagraph configuration
    max_length : int
        Maximum metapath length (number of edges)

    Returns
    -------
    List[str]
        List of metapath abbreviations starting with G and ending with BP
    """
    gene_bp_metapaths = [
        "GpBP",
        "GiGpBP",
        "GcGpBP",
        "GrGpBP",
        "Gr>GpBP",
        "GpBPpGpBP",
        "GiGiGpBP",
        "GpPWpGpBP",
        "GpMFpGpBP",
        "GpCCpGpBP",
    ]

    existing = []
    edge_files = list((Path(metagraph.get("_data_dir", ".")) / "edges").glob("*.npz"))
    available_edges = {f.stem.replace(".sparse", "") for f in edge_files}

    for mp in gene_bp_metapaths:
        existing.append(mp)

    return existing


def create_node_index_mapping(
    hetmat: HetMat,
    df: pd.DataFrame,
    source_type: str,
    target_type: str,
    source_id_col: str,
    target_id_col: str
) -> pd.DataFrame:
    """
    Map external IDs to matrix indices.

    Parameters
    ----------
    hetmat : HetMat
        HetMat instance
    df : pd.DataFrame
        DataFrame with source and target ID columns
    source_type : str
        Source node type (e.g., "Gene")
    target_type : str
        Target node type (e.g., "Biological Process")
    source_id_col : str
        Column name for source IDs
    target_id_col : str
        Column name for target IDs

    Returns
    -------
    pd.DataFrame
        DataFrame with added source_idx and target_idx columns
    """
    source_nodes = hetmat.get_nodes(source_type)
    target_nodes = hetmat.get_nodes(target_type)

    source_id_to_idx = dict(zip(source_nodes["identifier"], source_nodes.index))
    target_id_to_idx = dict(zip(target_nodes["identifier"], target_nodes.index))

    result_df = df.copy()
    result_df["source_idx"] = result_df[source_id_col].map(source_id_to_idx)
    result_df["target_idx"] = result_df[target_id_col].map(target_id_to_idx)

    n_missing_source = result_df["source_idx"].isna().sum()
    n_missing_target = result_df["target_idx"].isna().sum()

    if n_missing_source > 0:
        print(f"Warning: {n_missing_source} source IDs not found in node list")
    if n_missing_target > 0:
        print(f"Warning: {n_missing_target} target IDs not found in node list")

    return result_df
