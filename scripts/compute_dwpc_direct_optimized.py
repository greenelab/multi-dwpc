# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Compute DWPC Directly via Matrix Multiplication
#
# This script computes Degree-Weighted Path Counts (DWPC) directly from
# the HetMat sparse matrices, without requiring the Docker API.
#
# ## Inputs
# - `data/` directory containing HetMat sparse matrices
# - `output/intermediate/hetio_bppg_all_GO_positive_growth*_filtered.csv` files
# - `output/permutations/` and `output/random_samples/` directories
#
# ## Outputs
# - `output/dwpc_direct/all_GO_positive_growth/` directory with DWPC results
#
# ## Advantages over API approach
# - No Docker dependency
# - Faster computation via vectorized operations
# - Parallel processing across metapaths
# - Full control over damping parameter

# %%
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import pandas as pd

# Setup repo root
if Path.cwd().name == "scripts":
    repo_root = Path("..").resolve()
else:
    repo_root = Path.cwd()

sys.path.insert(0, str(repo_root / "src"))

from dwpc_direct import (
    HetMat,
    reverse_metapath_abbrev,
)

print(f"Repository root: {repo_root}")
print("DWPC direct computation module loaded")

# %% [markdown]
# ## Configuration

def process_dataset(config, hetmat, metapaths, gene_id_to_idx, bp_id_to_idx):
    """
    Process a single dataset: load data, compute DWPC, save results.

    Uses vectorized operations for efficient computation. DWPC matrices
    are already precomputed and cached, so this just extracts values.

    Parameters
    ----------
    config : dict
        Dataset configuration
    hetmat : HetMat
        HetMat instance with precomputed matrices
    metapaths : list
        List of metapath strings
    gene_id_to_idx : dict
        Gene ID to matrix index mapping
    bp_id_to_idx : dict
        BP ID to matrix index mapping

    Returns
    -------
    dict
        Processing summary
    """
    name = config["name"]
    path = config["path"]

    if not path.exists():
        return {"name": name, "status": "skipped", "reason": "file not found"}

    # Load data
    df = pd.read_csv(path)
    n_pairs = len(df)

    # Map IDs to indices (vectorized)
    df["source_idx"] = df[config["gene_col"]].map(gene_id_to_idx)
    df["target_idx"] = df[config["go_col"]].map(bp_id_to_idx)

    # Check for unmapped IDs
    n_unmapped_genes = df["source_idx"].isna().sum()
    n_unmapped_bps = df["target_idx"].isna().sum()

    if n_unmapped_genes > 0 or n_unmapped_bps > 0:
        print(f"  Warning: {n_unmapped_genes} unmapped genes, {n_unmapped_bps} unmapped BPs")

    # Filter to mapped pairs only
    df_mapped = df.dropna(subset=["source_idx", "target_idx"]).copy()
    df_mapped["source_idx"] = df_mapped["source_idx"].astype(int)
    df_mapped["target_idx"] = df_mapped["target_idx"].astype(int)

    # Extract arrays for vectorized lookup
    source_indices = df_mapped["source_idx"].values
    target_indices = df_mapped["target_idx"].values
    gene_ids = df_mapped[config["gene_col"]].values
    go_ids = df_mapped[config["go_col"]].values

    # Compute DWPC for each metapath (vectorized extraction from cached matrices)
    results_list = []
    for metapath in metapaths:
        try:
            dwpc_values = hetmat.get_dwpc_for_pairs(
                metapath, source_indices, target_indices
            )

            # Build result arrays vectorized (no row-by-row iteration)
            mp_results = pd.DataFrame({
                "entrez_gene_id": gene_ids,
                "go_id": go_ids,
                "source_idx": source_indices,
                "target_idx": target_indices,
                "metapath": metapath,
                "dwpc": dwpc_values
            })
            results_list.append(mp_results)

        except (KeyError, TypeError, ValueError) as e:
            print(f"  Error computing {metapath}: {e}")
            continue

    # Concatenate all results
    results_df = pd.concat(results_list, ignore_index=True) if results_list else pd.DataFrame()

    return {
        "name": name,
        "status": "completed",
        "n_input_pairs": n_pairs,
        "n_mapped_pairs": len(df_mapped),
        "n_results": len(results_df),
        "results_df": results_df
    }


def main():
    # %% [markdown]
    # ## Configuration

    # %%
    # Configuration
    DATA_DIR = repo_root / "data"
    OUTPUT_DIR = repo_root / "output" / "dwpc_direct" / "all_GO_positive_growth"
    INTERMEDIATE_DIR = repo_root / "output" / "intermediate"

    # DWPC parameters (0.5 matches API, arcsinh transformation applied)
    DAMPING = 0.5

    # Parallel processing.
    # The worker pool can be unstable on some Debian + Python 3.13 VMs when
    # NumPy/OpenBLAS is forked. Keep the default serialized to avoid crashes on
    # memory-constrained hosts.
    N_WORKERS = 1

    # Create output directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "results").mkdir(exist_ok=True)

    print(f"Data directory: {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Damping exponent: {DAMPING}")
    print(f"Parallel workers: {N_WORKERS}")

    # %% [markdown]
    # ## Load HetMat Data

    # %%
    # Initialize HetMat with disk caching enabled
    hetmat = HetMat(DATA_DIR, damping=DAMPING, use_disk_cache=True)

    # Load node information
    gene_nodes = hetmat.get_nodes("Gene")
    bp_nodes = hetmat.get_nodes("Biological Process")

    print(f"Loaded {len(gene_nodes)} genes")
    print(f"Loaded {len(bp_nodes)} biological processes")
    print(f"Disk cache directory: {hetmat.cache_dir}")

    # Show sample nodes
    print("\nSample genes:")
    print(gene_nodes.head())

    print("\nSample biological processes:")
    print(bp_nodes.head())

    # %% [markdown]
    # ## Define Metapaths
    #
    # We compute DWPC for metapaths connecting genes to biological processes.
    # These are the same metapaths used in the connectivity-search-backend.

    # %%
    # Define metapaths to compute
    # Use all BP->G metapaths (length <= 3) from the metapath stats file,
    # excluding the direct BPpG path, then reverse to G->BP for direct computation.
    metapath_stats = pd.read_csv(DATA_DIR / "metapath-dwpc-stats.tsv", sep="\t")
    bp_to_g = metapath_stats[
        metapath_stats["metapath"].str.startswith("BP")
        & metapath_stats["metapath"].str.endswith("G")
    ]
    bp_to_g = bp_to_g[bp_to_g["metapath"] != "BPpG"]
    METAPATHS = [
        reverse_metapath_abbrev(mp) for mp in sorted(bp_to_g["metapath"].unique())
    ]

    print(f"Computing {len(METAPATHS)} metapaths:")
    for mp in METAPATHS:
        print(f"  - {mp}")

    # %% [markdown]
    # ## Precompute DWPC Matrices
    #
    # Precompute all DWPC matrices in parallel before processing datasets.
    # This ensures each matrix is computed only once and reused across all datasets.
    # Matrices are cached to disk for future pipeline runs.

    # %%
    # Precompute all DWPC matrices (parallel, with disk caching)
    print("\n" + "=" * 80)
    print("PRECOMPUTING DWPC MATRICES")
    print("=" * 80)
    print(f"Computing {len(METAPATHS)} matrices with {N_WORKERS} workers...")
    print("(Cached matrices will be loaded from disk if available)")

    precompute_start = time.perf_counter()
    hetmat.precompute_matrices(METAPATHS, n_workers=N_WORKERS, show_progress=True)
    precompute_time = time.perf_counter() - precompute_start

    print(f"\nPrecomputation complete: {precompute_time:.1f}s")

    # Show cache status
    cached_files = list(hetmat.cache_dir.glob("dwpc_*.npz"))
    print(f"Cached matrices on disk: {len(cached_files)}")
    for f in sorted(cached_files):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name}: {size_mb:.2f} MB")

    # %% [markdown]
    # ## Dataset Configuration

    # %%
    # Define datasets to process
    # Using all_GO_positive_growth (Jaccard-filtered all GO positive growth terms)
    datasets_config = [
        # 2016 Real
        {
            "name": "all_GO_positive_growth_2016_real",
            "path": INTERMEDIATE_DIR / "hetio_bppg_all_GO_positive_growth_filtered.csv",
            "gene_col": "entrez_gene_id",
            "go_col": "go_id",
            "year": 2016,
            "type": "real"
        },
        # 2024 Real
        {
            "name": "all_GO_positive_growth_2024_real",
            "path": INTERMEDIATE_DIR / "hetio_bppg_all_GO_positive_growth_2024_filtered.csv",
            "gene_col": "entrez_gene_id",
            "go_col": "go_id",
            "year": 2024,
            "type": "real"
        },
    ]

    # Add permutation datasets if they exist
    perm_dir_2016 = repo_root / "output" / "permutations" / "all_GO_positive_growth_2016"
    perm_dir_2024 = repo_root / "output" / "permutations" / "all_GO_positive_growth_2024"

    for perm_file in sorted(perm_dir_2016.glob("perm_*.csv")) if perm_dir_2016.exists() else []:
        perm_num = perm_file.stem.split("_")[1]
        datasets_config.append({
            "name": f"all_GO_positive_growth_2016_perm_{perm_num}",
            "path": perm_file,
            "gene_col": "entrez_gene_id",
            "go_col": "go_id",
            "year": 2016,
            "type": "permuted"
        })

    for perm_file in sorted(perm_dir_2024.glob("perm_*.csv")) if perm_dir_2024.exists() else []:
        perm_num = perm_file.stem.split("_")[1]
        datasets_config.append({
            "name": f"all_GO_positive_growth_2024_perm_{perm_num}",
            "path": perm_file,
            "gene_col": "entrez_gene_id",
            "go_col": "go_id",
            "year": 2024,
            "type": "permuted"
        })

    # Add random datasets if they exist
    random_dir_2016 = repo_root / "output" / "random_samples" / "all_GO_positive_growth_2016"
    random_dir_2024 = repo_root / "output" / "random_samples" / "all_GO_positive_growth_2024"

    for random_file in sorted(random_dir_2016.glob("random_*.csv")) if random_dir_2016.exists() else []:
        random_num = random_file.stem.split("_")[1]
        datasets_config.append({
            "name": f"all_GO_positive_growth_2016_random_{random_num}",
            "path": random_file,
            "gene_col": "entrez_gene_id",
            "go_col": "go_id",
            "year": 2016,
            "type": "random"
        })

    for random_file in sorted(random_dir_2024.glob("random_*.csv")) if random_dir_2024.exists() else []:
        random_num = random_file.stem.split("_")[1]
        datasets_config.append({
            "name": f"all_GO_positive_growth_2024_random_{random_num}",
            "path": random_file,
            "gene_col": "entrez_gene_id",
            "go_col": "go_id",
            "year": 2024,
            "type": "random"
        })

    print(f"\nConfigured {len(datasets_config)} datasets:")
    for ds in datasets_config:
        exists = ds["path"].exists()
        status = "EXISTS" if exists else "MISSING"
        print(f"  {ds['name']}: {status}")

    # %% [markdown]
    # ## Create Node Index Mappings

    # %%
    # Create ID to index mappings
    # Prefer the explicit 'position' column when available; otherwise, fall back to the DataFrame index
    gene_positions = gene_nodes["position"] if "position" in gene_nodes.columns else gene_nodes.index
    bp_positions = bp_nodes["position"] if "position" in bp_nodes.columns else bp_nodes.index
    gene_id_to_idx = dict(zip(gene_nodes["identifier"], gene_positions))
    bp_id_to_idx = dict(zip(bp_nodes["identifier"], bp_positions))

    print(f"Gene ID mapping: {len(gene_id_to_idx)} entries")
    print(f"BP ID mapping: {len(bp_id_to_idx)} entries")

    # Sample mappings
    print("\nSample gene mappings (entrez_id -> matrix_index):")
    sample_genes = list(gene_id_to_idx.items())[:5]
    for gid, idx in sample_genes:
        print(f"  {gid} -> {idx}")

    print("\nSample BP mappings (go_id -> matrix_index):")
    sample_bps = list(bp_id_to_idx.items())[:5]
    for bpid, idx in sample_bps:
        print(f"  {bpid} -> {idx}")

    # %% [markdown]
    # ## Process All Datasets

    # %%
    # Process all datasets
    print("=" * 80)
    print("PROCESSING DATASETS")
    print("=" * 80)

    batch_start = time.perf_counter()
    all_summaries = []

    for i, config in enumerate(datasets_config, 1):
        print(f"\n[{i}/{len(datasets_config)}] Processing {config['name']}...")

        output_path = OUTPUT_DIR / "results" / f"dwpc_{config['name']}.csv"

        # Check if already processed
        if output_path.exists():
            existing_df = pd.read_csv(output_path)
            print(f"  Already processed: {len(existing_df)} rows")
            all_summaries.append({
                "name": config["name"],
                "status": "skipped (exists)",
                "n_results": len(existing_df)
            })
            continue

        start_time = time.perf_counter()

        summary = process_dataset(
            config, hetmat, METAPATHS, gene_id_to_idx, bp_id_to_idx
        )

        if summary["status"] == "completed":
            # Save results
            summary["results_df"].to_csv(output_path, index=False)
            print(f"  Saved: {output_path}")
            print(f"  Rows: {summary['n_results']}")

            elapsed = time.perf_counter() - start_time
            print(f"  Time: {elapsed:.1f}s")

            summary["time_seconds"] = elapsed
            del summary["results_df"]  # Don't keep in memory

        all_summaries.append(summary)

    batch_time = time.perf_counter() - batch_start
    print(f"\n{'=' * 80}")
    print(f"BATCH COMPLETE: {batch_time / 60:.1f} minutes")
    print(f"{'=' * 80}")

    # %% [markdown]
    # ## Summary

    # %%
    # Print summary
    summary_df = pd.DataFrame(all_summaries)
    print("\nProcessing Summary:")
    print(summary_df.to_string(index=False))

    # Count results
    completed = sum(1 for s in all_summaries if s.get("status") == "completed")
    skipped_exists = sum(1 for s in all_summaries if "exists" in str(s.get("status", "")))
    skipped_missing = sum(1 for s in all_summaries if s.get("status") == "skipped")

    print("\nResults:")
    print(f"  Completed: {completed}")
    print(f"  Skipped (already exists): {skipped_exists}")
    print(f"  Skipped (file missing): {skipped_missing}")
    print(f"  Total: {len(datasets_config)}")

    # %% [markdown]
    # ## Verify Output Files

    # %%
    # List output files
    print("\nOutput files:")
    for f in sorted((OUTPUT_DIR / "results").glob("*.csv")):
        size_mb = f.stat().st_size / (1024 * 1024)
        df = pd.read_csv(f)
        print(f"  {f.name}: {len(df):,} rows ({size_mb:.2f} MB)")

    # %%
    # Sample output
    sample_file = OUTPUT_DIR / "results" / "dwpc_all_GO_positive_growth_2016_real.csv"
    if sample_file.exists():
        sample_df = pd.read_csv(sample_file)
        print("\nSample output (first 10 rows):")
        print(sample_df.head(10))

        print("\nDWPC statistics by metapath:")
        stats = sample_df.groupby("metapath")["dwpc"].agg(["mean", "std", "min", "max", "count"])
        print(stats)


if __name__ == "__main__":
    main()

