# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] papermill={"duration": 0.001869, "end_time": "2025-12-05T20:29:36.322002", "exception": false, "start_time": "2025-12-05T20:29:36.320133", "status": "completed"}
# # Data Loading
#
# Download or Load 2016 Hetionet BP-Gene associations and 2024 GO annotations.
#
# ## Description
# This notebook loads data from Hetionet v1.0 (2016) and GO annotations
# from 2024. It identifies common GO terms between both years and prepares
# the data for downstream processing.
#
# Filters applied:
# - Biological Process domain only
# - Genes present in Hetionet
# - GO terms with 2-1000 genes (Hetionet criterion)
# - Common genes between both years for fair comparison
#
# ## Data Sources
# - **Hetionet v1.0 (2016)**: https://github.com/dhimmel/hetionet
# - **GO Annotations (2024)**: https://github.com/NegarJanani/gene-ontology
# - **GO Ontology (2016)**: http://release.geneontology.org/2016-02-01/ontology/go-basic.obo
# - **GO Ontology (2024)**: http://purl.obolibrary.org/obo/go/go-basic.obo
#
# ## Inputs
# - `data/nodes/Gene.tsv`
# - `data/nodes/Biological Process.tsv`
# - `data/edges/GpBP.sparse.npz`
#
# ## Outputs
# - `output/intermediate/hetio_bppg_2016.csv`
# - `output/intermediate/upd_go_bp_2024.csv`
# - `output/intermediate/common_go_terms.csv`

# %% [markdown] papermill={"duration": 0.00381, "end_time": "2025-12-05T20:29:36.330528", "exception": false, "start_time": "2025-12-05T20:29:36.326718", "status": "completed"}
# ### Read existing Biological Process (BP) GO terms in Hetionet (2016)
# This section loads the Hetionet 2016 BP-Gene associations from a CSV file, renames columns for clarity, and displays the first few rows. It then calculates and prints:
# The number of unique GO terms in Hetio (2016)
# The number of unique genes associated with these GO terms

# %% papermill={"duration": 2.108204, "end_time": "2025-12-05T20:29:38.442376", "exception": false, "start_time": "2025-12-05T20:29:36.334172", "status": "completed"}
import gc
import urllib.request

import pandas as pd
import scipy.sparse
from pathlib import Path
import hetnetpy.readwrite
import hetmatpy.hetmat

print("*** testing what's going on here", flush=True)

# Constants
HETIONET_GENE_MIN = 2
HETIONET_GENE_MAX = 1000
HETIONET_URL = (
    "https://github.com/dhimmel/hetionet/raw/"
    "76550e6c93fbe92124edc71725e8c7dd4ca8b1f5/hetnet/json/hetionet-v1.0.json.bz2"
)
GO_2024_URL = (
    "https://raw.githubusercontent.com/NegarJanani/gene-ontology/"
    "refs/heads/gh-pages/annotations/taxid_9606/"
    "GO_annotations-9606-inferred-allev.tsv"
)

# Expected ranges based on Hetionet v1.0 statistics
EXPECTED_EDGE_RANGE = (500_000, 1_500_000)
MIN_EXPECTED_GENES = 14_000
MIN_EXPECTED_GO_TERMS = 5_000

# Setup directories relative to repo root
# Works whether notebook is run from repo root or notebooks/ subdirectory
if Path.cwd().name == "notebooks":
    repo_root = Path("..").resolve()
else:
    repo_root = Path.cwd()

data_dir = repo_root / "data"
output_dir = repo_root / "output" / "intermediate"
data_dir.mkdir(exist_ok=True)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Data directory: {data_dir.resolve()}")
print(f"Output directory: {output_dir.resolve()}")

# %% papermill={"duration": 0.00773, "end_time": "2025-12-05T20:29:38.451624", "exception": false, "start_time": "2025-12-05T20:29:38.443894", "status": "completed"}
hetmat_path = data_dir
gpbp_edge_file = hetmat_path / "edges" / "GpBP.sparse.npz"

# Check if GpBP edge file already exists
if gpbp_edge_file.exists():
    print(f"GpBP edge file exists at {gpbp_edge_file}")
else:
    print("GpBP edge file not found, creating hetmat...")
    print(f"Downloading graph from {HETIONET_URL}")
    graph = hetnetpy.readwrite.read_graph(HETIONET_URL)
    print("Graph loaded, creating hetmat matrices...")
    hetmat = hetmatpy.hetmat.hetmat_from_graph(graph, str(hetmat_path))
    print(f"Hetmat created and saved to {hetmat_path}")

# %% papermill={"duration": 0.112978, "end_time": "2025-12-05T20:29:41.640242", "exception": false, "start_time": "2025-12-05T20:29:41.527264", "status": "completed"}
# Load sparse adjacency matrix directly
gpbp_matrix_file = data_dir / "edges" / "GpBP.sparse.npz"
adj_matrix = scipy.sparse.load_npz(gpbp_matrix_file)

# Load node mappings
gene_nodes = pd.read_csv(data_dir / "nodes" / "Gene.tsv", sep="\t")
bp_nodes = pd.read_csv(data_dir / "nodes" / "Biological Process.tsv", sep="\t")

# Convert sparse matrix to edge list
rows, cols = adj_matrix.nonzero()

# Map matrix indices to identifiers
hetio_BPpG_df = pd.DataFrame({
    'entrez_gene_id': gene_nodes.loc[rows, 'identifier'].values,
    'metaedge': 'GpBP',
    'go_id': bp_nodes.loc[cols, 'identifier'].values
})

print(f"Loaded {len(hetio_BPpG_df)} edges")
print(f"Matrix shape: {adj_matrix.shape}")
hetio_BPpG_df.head()

# %% papermill={"duration": 0.040313, "end_time": "2025-12-05T20:29:41.823463", "exception": false, "start_time": "2025-12-05T20:29:41.783150", "status": "completed"}
# Validate loaded data against expected Hetionet v1.0 statistics
edge_count = len(hetio_BPpG_df)
assert EXPECTED_EDGE_RANGE[0] <= edge_count <= EXPECTED_EDGE_RANGE[1], \
    f"Edge count {edge_count} outside expected range {EXPECTED_EDGE_RANGE}"
print(f"PASSED: Edge count ({edge_count}) in expected range")

assert hetio_BPpG_df.isnull().sum().sum() == 0, "Found null values in data"
print("PASSED: No null values")

expected_cols = {'entrez_gene_id', 'metaedge', 'go_id'}
assert set(hetio_BPpG_df.columns) == expected_cols, "Unexpected columns"
print("PASSED: Correct columns")

gene_count = hetio_BPpG_df['entrez_gene_id'].nunique()
go_count = hetio_BPpG_df['go_id'].nunique()
assert gene_count >= MIN_EXPECTED_GENES, f"Gene count {gene_count} too low"
assert go_count >= MIN_EXPECTED_GO_TERMS, f"GO term count {go_count} too low"
print(f"PASSED: Unique genes ({gene_count}) and GO terms ({go_count})")

print("\nAll tests PASSED")

# %% [markdown] papermill={"duration": 0.003343, "end_time": "2025-12-05T20:29:41.830511", "exception": false, "start_time": "2025-12-05T20:29:41.827168", "status": "completed"}
# ### Count the number of genes associated with each GO (BP) term
# This section counts how many genes are associated with each GO term in Hetio (2016) and displays the resulting DataFrame.

# %% papermill={"duration": 0.020824, "end_time": "2025-12-05T20:29:41.854693", "exception": false, "start_time": "2025-12-05T20:29:41.833869", "status": "completed"}
hetio_bp_g_freq_df = hetio_BPpG_df['go_id'].value_counts().reset_index()
hetio_bp_g_freq_df.columns = ['go_id', 'no_of_genes_in_hetio_GO_2016']

print(f"hetio_bp_g_freq_df shape: {hetio_bp_g_freq_df.shape}")
print(hetio_bp_g_freq_df.head())

# %% [markdown]
# ### Read updated Biological Process (BP) GO terms (2024)
# The annotation file can expand substantially when gene_ids/gene_symbols are
# exploded. To keep peak memory bounded, download it to disk once and process
# it in chunks, filtering each chunk before retaining it.

# Small lookup sets needed while streaming the 2024 annotations.
hetio_genes = set(hetio_BPpG_df['entrez_gene_id'].unique())
hetio_go_terms = set(hetio_BPpG_df['go_id'].unique())

# The sparse matrix and node tables are no longer needed after constructing
# hetio_BPpG_df. Release them before loading/expanding the 2024 annotations.
del adj_matrix, rows, cols, gene_nodes, bp_nodes
gc.collect()
print("Released Hetionet matrix/node intermediates", flush=True)

# Cache the remote annotation file locally. Streaming a local file keeps the
# network layer out of pandas' parsing path and makes subsequent runs cheaper.
go_2024_file = data_dir / "GO_annotations-9606-inferred-allev.tsv"
if not go_2024_file.exists():
    print(f"Downloading GO 2024 annotations to {go_2024_file}...", flush=True)
    try:
        urllib.request.urlretrieve(GO_2024_URL, go_2024_file)
    except Exception as e:
        raise RuntimeError(f"Failed to fetch GO 2024 annotations: {e}") from e
    print("Download complete", flush=True)
else:
    print(f"Using cached GO 2024 annotations: {go_2024_file}", flush=True)

# Only these columns are used below. A conservative chunk size keeps the
# temporary list/object allocations created by str.split/explode bounded.
GO_CHUNK_SIZE = 25_000
GO_USECOLS = ['go_id', 'go_name', 'go_domain', 'gene_ids', 'gene_symbols']

filtered_chunks = []
bp_pairs_before_hetio_gene_filter = 0
bp_genes_before_hetio_gene_filter = set()
raw_rows_processed = 0

try:
    reader = pd.read_csv(
        go_2024_file,
        sep='\t',
        usecols=GO_USECOLS,
        dtype='string',
        chunksize=GO_CHUNK_SIZE,
    )

    for chunk_number, chunk in enumerate(reader, start=1):
        raw_rows_processed += len(chunk)

        # Filter before explode: this is the main memory-saving step.
        chunk = chunk[chunk['go_domain'] == 'biological_process'].copy()
        chunk = chunk[chunk['go_id'].isin(hetio_go_terms)].copy()

        if chunk.empty:
            del chunk
            continue

        chunk['gene_id'] = chunk['gene_ids'].str.split('|')
        chunk['gene_symbol'] = chunk['gene_symbols'].str.split('|')
        chunk = chunk.explode(['gene_id', 'gene_symbol'], ignore_index=True)

        # Clean IDs without making a second full-size raw dataframe.
        chunk['gene_id'] = pd.to_numeric(
            chunk['gene_id'].str.strip(), errors='coerce'
        )
        chunk['gene_symbol'] = chunk['gene_symbol'].str.strip()
        chunk = chunk.dropna(subset=['gene_id'])
        chunk['gene_id'] = chunk['gene_id'].astype('int64')

        bp_pairs_before_hetio_gene_filter += len(chunk)
        bp_genes_before_hetio_gene_filter.update(chunk['gene_id'].unique())

        # Keep only genes that can survive the downstream Hetionet filter.
        chunk = chunk[chunk['gene_id'].isin(hetio_genes)]
        if not chunk.empty:
            filtered_chunks.append(
                chunk[['go_id', 'go_name', 'gene_id', 'gene_symbol']].rename(
                    columns={'gene_id': 'entrez_gene_id'}
                )
            )

        if chunk_number % 10 == 0:
            print(
                f"Processed {raw_rows_processed:,} raw GO annotation rows...",
                flush=True,
            )

        del chunk
        gc.collect()
except Exception as e:
    raise RuntimeError(f"Failed while processing GO 2024 annotations: {e}") from e

if filtered_chunks:
    upd_go_bp_2024_df = pd.concat(filtered_chunks, ignore_index=True)
else:
    upd_go_bp_2024_df = pd.DataFrame(
        columns=['go_id', 'go_name', 'entrez_gene_id', 'gene_symbol']
    )

del filtered_chunks, reader
gc.collect()

print(f"Processed {raw_rows_processed:,} raw GO annotation rows")
print(
    f"Biological Process pairs in Hetionet terms before gene filtering: "
    f"{bp_pairs_before_hetio_gene_filter:,}"
)
genes_before = len(bp_genes_before_hetio_gene_filter)
genes_after = upd_go_bp_2024_df['entrez_gene_id'].nunique()
print(f"After filtering to Hetionet genes: {len(upd_go_bp_2024_df)} pairs "
      f"({genes_after}/{genes_before} genes)")
del bp_genes_before_hetio_gene_filter
gc.collect()

# Filter 3: Terms with 2-1000 genes (Hetionet criterion)
upd_go_bp_2024_freq_df = (
    upd_go_bp_2024_df
    .groupby('go_id')
    .size()
    .reset_index(name='no_of_genes_in_GO_2024')
)
upd_go_bp_2024_freq_df = upd_go_bp_2024_freq_df[
    (upd_go_bp_2024_freq_df['no_of_genes_in_GO_2024'] >= HETIONET_GENE_MIN) & 
    (upd_go_bp_2024_freq_df['no_of_genes_in_GO_2024'] <= HETIONET_GENE_MAX)
]
upd_go_bp_2024_df = upd_go_bp_2024_df[
    upd_go_bp_2024_df['go_id'].isin(upd_go_bp_2024_freq_df['go_id'])
]
print(f"After filtering to {HETIONET_GENE_MIN}-{HETIONET_GENE_MAX} genes: "
      f"{len(upd_go_bp_2024_df)} pairs ({len(upd_go_bp_2024_freq_df)} terms)")

# Filter 4: Terms present in Hetionet
upd_go_bp_2024_df = upd_go_bp_2024_df[
    upd_go_bp_2024_df['go_id'].isin(hetio_go_terms)
]
upd_go_bp_2024_freq_df = upd_go_bp_2024_freq_df[
    upd_go_bp_2024_freq_df['go_id'].isin(hetio_go_terms)
]
print(f"After filtering to Hetionet terms: {len(upd_go_bp_2024_df)} pairs "
      f"({len(upd_go_bp_2024_freq_df)} terms)")

upd_go_bp_2024_df.head()

# %% [markdown] papermill={"duration": 0.007123, "end_time": "2025-12-05T20:29:44.413331", "exception": false, "start_time": "2025-12-05T20:29:44.406208", "status": "completed"}
# ### As hetionet considered GO terms with Biological processes with 2–1000 annotated genes were included. Same filter were applied for GO 2024  https://github.com/dhimmel/gene-ontology/issues/9

# %% papermill={"duration": 0.091813, "end_time": "2025-12-05T20:29:44.529038", "exception": false, "start_time": "2025-12-05T20:29:44.437225", "status": "completed"}
# Filter both datasets to common gene universe for fair comparison
genes_2024 = set(upd_go_bp_2024_df['entrez_gene_id'].unique())
genes_2016 = set(hetio_BPpG_df['entrez_gene_id'].unique())
common_genes = genes_2016 & genes_2024

print(f"Genes in 2016: {len(genes_2016)}")
print(f"Genes in 2024: {len(genes_2024)}")
print(f"Common genes: {len(common_genes)}")
print(f"Genes only in 2016: {len(genes_2016 - genes_2024)}")
print(f"Genes only in 2024: {len(genes_2024 - genes_2016)}")

# Filter both datasets to common genes
hetio_BPpG_df = hetio_BPpG_df[hetio_BPpG_df['entrez_gene_id'].isin(common_genes)]
upd_go_bp_2024_df = upd_go_bp_2024_df[
    upd_go_bp_2024_df['entrez_gene_id'].isin(common_genes)
]

print(f"\nAfter filtering to common genes:")
print(f"Hetionet 2016: {len(hetio_BPpG_df)} pairs")
print(f"2024: {len(upd_go_bp_2024_df)} pairs")

# Recalculate frequencies with common gene universe
hetio_bp_g_freq_df = (
    hetio_BPpG_df
    .groupby('go_id')
    .size()
    .reset_index(name='no_of_genes_in_hetio_GO_2016')
)
upd_go_bp_2024_freq_df = (
    upd_go_bp_2024_df
    .groupby('go_id')
    .size()
    .reset_index(name='no_of_genes_in_GO_2024')
)

# Reapply 2-1000 filter
upd_go_bp_2024_freq_df = upd_go_bp_2024_freq_df[
    (upd_go_bp_2024_freq_df['no_of_genes_in_GO_2024'] >= HETIONET_GENE_MIN) & 
    (upd_go_bp_2024_freq_df['no_of_genes_in_GO_2024'] <= HETIONET_GENE_MAX)
]
upd_go_bp_2024_df = upd_go_bp_2024_df[
    upd_go_bp_2024_df['go_id'].isin(upd_go_bp_2024_freq_df['go_id'])
]

print(f"After reapplying {HETIONET_GENE_MIN}-{HETIONET_GENE_MAX} gene filter: "
      f"{len(upd_go_bp_2024_freq_df)} terms in 2024")

# %% papermill={"duration": 0.010995, "end_time": "2025-12-05T20:29:44.543889", "exception": false, "start_time": "2025-12-05T20:29:44.532894", "status": "completed"}
# Merge 2016 and 2024 frequencies for common GO terms
common_terms_df = pd.merge(
    hetio_bp_g_freq_df,
    upd_go_bp_2024_freq_df,
    on='go_id',
    how='inner'
)
common_terms_df

# %% papermill={"duration": 0.665853, "end_time": "2025-12-05T20:29:45.220565", "exception": false, "start_time": "2025-12-05T20:29:44.554712", "status": "completed"}
# Save intermediate outputs
print('Saving intermediate outputs...')

# Save 2016 BP-Gene associations
hetio_BPpG_df.to_csv(output_dir / 'hetio_bppg_2016.csv', index=False)
print(f'Saved hetio_bppg_2016.csv: {len(hetio_BPpG_df)} rows')

# Save 2024 BP-Gene associations
upd_go_bp_2024_df.to_csv(output_dir / 'upd_go_bp_2024.csv', index=False)
print(f'Saved upd_go_bp_2024.csv: {len(upd_go_bp_2024_df)} rows')

# Save common GO terms with both 2016 and 2024 counts
common_terms_df.to_csv(output_dir / 'common_go_terms.csv', index=False)
print(f'Saved common_go_terms.csv: {len(common_terms_df)} GO terms')

print('\n Data loading complete!')
