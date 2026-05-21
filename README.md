# multi-dwpc
This repository contains scripts testing aggregate degree weighted path count (DWPC) statistics. 

# Motivation
Pathway search methods focus on explaining the relationship between two entities, such as a gene and a disease. However, most diseases are the result of more complex interactions between multiple genes and molecular features. For example, many cancers result from mutations across multiple genes. Trisomy 21 is a condition based on the upregulation of the ~225 genes on HSA21, and single-gene searches do not capture the collective mechanistic evidence. Sets of coexpressed or interacting genes thus better capture the cause or effect of disease.  Methods of gene set enrichment do not clearly represent relationships across data modalities that are better captured within a heterogeneous knowledge graph. 

Here we present the scripts for a proof of concept experiment that considers whether the aggregated pathway statistics (mean/median of DWPC) between GO terms and gene annotations differ between the 2015 and 2024 additions. The hypothesis is that the connectivity of hetionet (i.e., the shared connections of genes and pathways to other nodes such as compound, anatomy, and disease) between the 2015 and 2024 GO annotations will result in similar aggregated degree-adjusted pathway scores. Similar scores for the updated GO terms demonstrate the potential of aggregated pathway scores for identifying shared mechanisms of genes not previously annotated to a single shared node (e.g., GO pathway or disease).

# Previous work
Several of the scripts build on previous work from the Greene Lab and [hetionet project](https://het.io/), including the [connectivity-search-analyses](https://github.com/greenelab/connectivity-search-analyses) repository, [hetio/hetnetpy](https://github.com/hetio/hetnetpy), and [hetio/hetmatpy](https://github.com/hetio/hetmatpy). 

# Getting Started

## Clone the repository

```bash
git clone https://github.com/greenelab/multi-dwpc.git
cd multi-dwpc
```

## Create the environment

```bash
conda env create -f env/environment.yml
conda activate multi_dwpc
pip install -e ".[dev]"
```

## Prepare the Data

```bash
# Run individual data-prep steps (no top-level pipeline runner; chain
# the steps you need)
poe load-data
poe filter-change
poe go-hierarchy-analysis
poe filter-jaccard

# Generate year null datasets
poe gen-permutation
poe gen-random

# Compute DWPC matrices for the prepared data
poe compute-dwpc-direct
```

For end-to-end runs, use the unified pipeline runners (or the HPC submit scripts):

```bash
# Year track end-to-end (B-selection -> intermediate sharing -> consumables)
python3 scripts/pipeline/run_year_pipeline.py --output-dir output/year_full_analysis

# LV track end-to-end
python3 scripts/pipeline/run_lv_pipeline.py --output-dir output/lv_full_analysis
```

### Available tasks

Run `poe --help` to see all available tasks:

| Task | Description |
|------|-------------|
| **Data prep** | |
| `load-data` | Load Hetionet v1.0 (2016) and GO annotations (2024) |
| `filter-change` | Percent-change + IQR filtering (all_GO_positive_growth) |
| `go-hierarchy-analysis` | GO hierarchy metrics + parents_GO_postive_growth |
| `filter-jaccard` | Jaccard filtering (all_GO, plus parents_GO when available) |
| `gen-permutation` | Generate permutation null datasets (year cohorts) |
| `gen-random` | Generate promiscuity-matched random null datasets (year cohorts) |
| `pipeline-null` | Run both null dataset generators in sequence |
| **DWPC compute / validation** | |
| `compute-dwpc-direct` | Compute DWPC via direct sparse-matrix multiplication |
| `warmup-dwpc-cache` | Precompute and cache year DWPC metapath matrices only |
| `validate-dwpc-concordance` | Compare direct DWPC values against historical API results |
| `lookup-dwpc-api` | Lookup DWPC via the Docker API stack |
| `lookup-dwpc-api-with-docker` | Start Docker stack, wait for API, then run DWPC |
| `test-dwpc-accuracy` | Validate direct DWPC computation against API |
| `benchmark-dwpc` | Benchmark direct vs API DWPC computation |
| **Year experiments + plots** | |
| `year-null-variance-exp` / `year-null-variance-plots` | Year null-variance analysis + plots |
| `year-rank-stability-exp` / `year-rank-stability-plots` | Year metapath rank-stability analysis + plots |
| `year-rank-seed-comparisons` | Per-seed metapath rank scatter against a reference seed |
| `year-go-term-support` | Build year GO-term and global metapath support tables |
| `year-snapshot-rank-similarity` | Year 2016-vs-2024 metapath rank similarity at fixed support |
| `year-api-compare-to-mean-dwpc`, `year-api-sampled-statistics-analysis` | Year API analyses |
| **LV preparation + experiments + plots** | |
| `lv-prepare-exp` | Prepare shared LV inputs for explicit null-replicate experiments |
| `lv-gen-permutation` / `lv-gen-random` | Generate explicit LV permuted / random replicate artifacts |
| `lv-compute-replicate-summaries` | Compute one LV metapath-summary artifact per replicate |
| `lv-null-variance-exp` / `lv-null-variance-plots` | LV null-variance analysis + plots |
| `lv-rank-stability-exp` / `lv-rank-stability-plots` | LV rank-stability analysis + plots |
| `lv-rank-seed-comparisons` | Per-seed metapath rank scatter against a reference seed |
| **LV multi-DWPC pipeline (staged orchestrator)** | |
| `lv-top-genes` / `lv-target-sets` / `lv-real-dwpc` | Per-stage entrypoints |
| `lv-precompute-scores` / `lv-nulls` / `lv-stats` | Per-stage entrypoints |
| `lv-top-subgraphs` / `lv-plot-subgraphs` | Per-stage entrypoints |
| `lv-pipeline-fast` / `lv-pipeline-fast-smoke` | Full / smoke runs |
| **Maintenance** | |
| `clean` | Remove generated data directories |

The null-variance / rank-stability / rank-seed-comparison experiments and
plots are dispatched on `--analysis-type {year,lv}` behind the scenes; the
`year-*` and `lv-*` task names are thin wrappers that pre-fill the flag.

Note: `filter-jaccard` includes parents_GO_postive_growth only when run with
`python scripts/data_prep/jaccard_similarity_and_filtering.py --include-parents`.

### Pipeline scripts

Layout (under `scripts/`):

#### `data_prep/` — load, filter, generate nulls, compute DWPC
- `load_data.py` — Loads Hetionet v1.0 (2016) and GO annotations (2024), filters to common genes and GO terms
- `percent_change_and_filtering.py` — Percent-change + IQR filtering for all_GO_positive_growth
- `go_hierarchy_analysis.py` — GO hierarchy metrics and parents_GO_postive_growth generation
- `jaccard_similarity_and_filtering.py` — Jaccard filtering for all_GO (and parents_GO when available)
- `build_null_datasets.py --method {permutation,random}` — Permutation- or promiscuity-matched random null datasets (year cohorts)
- `compute_dwpc_direct.py` — Direct DWPC computation via sparse-matrix multiplication
- `build_year_top_genes.py` — Build the per-(GO, gene) table consumed by `intermediate_sharing.py --analysis-type year`
- `download_lv_loadings.py` — Download MultiPLIER Z matrix
- `go_hierarchy_analysis.py` — GO hierarchy metric computation

#### `pipeline/` — analysis stages chained by the runners
- `intermediate_sharing.py --analysis-type {year,lv}` — Computes intermediate sharing per metapath
- `select_optimal_b.py` — Joint stabilization-curve B selection from null-variance + rank-stability outputs
- `select_top_go_terms.py` — Pick top-N GO terms by intermediate-sharing signal (year only)
- `generate_global_summary.py` — Roll up per-metapath sharing into a global summary
- `generate_gene_table.py` — Per-gene connectivity + DWPC table (consumable for downstream visualization)
- `run_year_pipeline.py` / `run_lv_pipeline.py` — End-to-end orchestrators that chain the above stages

#### `experiments/` — null-variance and rank-stability experiments
- `null_variance_experiment.py --analysis-type {year,lv}` — Null variance across B/seed
- `rank_stability_experiment.py --analysis-type {year,lv}` — Metapath rank stability across B/seed
- `lv_prepare_experiment.py`, `lv_generate_control_replicates.py`, `lv_compute_replicate_summaries.py` — LV experiment setup steps
- `lv_multidwpc_analysis.py` — Staged LV multi-DWPC orchestrator (run with `--stage <name>`)

#### `visualization/` — plots
- `plot_null_variance_results.py --analysis-type {year,lv}` — Null-variance plots
- `plot_rank_stability_results.py --analysis-type {year,lv}` — Rank-stability plots
- `plot_rank_seed_comparisons.py --analysis-type {year,lv}` — Seed-comparison rank scatter
- `plot_year_intermediate_sharing.py` / `plot_lv_intermediate_sharing.py` — Intermediate-sharing plots (domain-specific, not unified)
- `plot_metapath_subgraphs.py` — Subgraph rendering for top metapaths
- `plot_year_effective_metapath_selection.py` — Year-only metapath-selection diagnostic
- `year_go_term_support.py`, `year_snapshot_rank_similarity.py` — Year-only diagnostics
- `plot_dwpc_benchmark.py` — DWPC benchmarking plots

#### `api/` — API stack utilities
- `wait_for_api.py`, `lookup_dwpc_api.py`, `test_dwpc_accuracy.py`, `validate_dwpc_concordance.py`, `benchmark_dwpc_methods.py`, `benchmark_dwpc_permutation_scaling.py`, `year_api_compare_statistics_to_mean_dwpc.py`, `year_api_sampled_statistics_analysis.py`

### Dataset naming

- `all_GO_positive_growth`: all GO terms with positive growth after IQR filtering
- `parents_GO_postive_growth`: parents of leaf terms within the same filtered set

### LV loadings input

The LV pipeline expects an external loadings file. Use the curated MultiPLIER
Z matrix from:

- `greenelab/phenoplier`
- `data/input/multiplier/multiplier_model_z.tsv.gz`

Download this file into `data/lv_loadings/`:

```bash
poe download-lv-loadings
```

Then set `LV_LOADINGS_PATH` when running LV pipeline tasks:

```bash
export LV_LOADINGS_PATH=data/lv_loadings/multiplier_model_z.tsv.gz
```

## Run the DWPC computation

There are two methods for computing Degree-Weighted Path Counts (DWPC):

### Option A: Direct computation 

Computes DWPC directly from the HetMat sparse matrices using hetmatpy. This method is significantly faster and does not require Docker.

Before the larger direct-DWPC run, validate concordance against the historical
API outputs:

```bash
poe validate-dwpc-concordance
```

This writes a concordance summary table, sampled comparison rows, per-metapath
scatter plots, and overall scatter/error plots under
`output/dwpc_validation/all_GO_positive_growth/`.

```bash
poe compute-dwpc-direct
```

**First run:** Computes and caches all DWPC matrices 

**Subsequent runs:** Loads cached matrices from disk

**Testing accuracy:**

Validate that direct computation matches the API gold-standard values:

```bash
poe test-dwpc-accuracy
```

**Benchmarking:**

Compare direct computation vs API performance:

```bash
poe benchmark-dwpc
```

### Option B: API-based computation

Computes DWPC via the Hetionet API. This requires running the connectivity-search-backend Docker stack.

**1. Start the Docker stack:**

```bash
cd connectivity-search-backend
./run_stack.sh
```

This will:
- Download the postgres database dump (~5 GB) on first run
- Set up the `.env` file with required secrets
- Start the postgres, neo4j, and API containers

The initial database load takes approximately 30 minutes for postgres and 10 minutes for neo4j.

**2. Verify the API is running:**

Wait for the containers to become healthy, then test the API:

```bash
curl http://localhost:8015/v1/nodes/
```

**3. Run the API-based computation:**

```bash
poe lookup-dwpc-api
```

If you want to start the Docker stack and wait for the API in one step:

```bash
poe lookup-dwpc-api-with-docker
```

**4. Stop the Docker stack when finished:**

```bash
cd connectivity-search-backend
docker compose down
```

### Performance comparison

Direct computation is 1,000-6,000x faster than API lookups after matrices are cached:

| Method | Time per pair | Notes |
|--------|---------------|-------|
| Direct (cached) | 0.002-0.01 ms | After initial matrix computation |
| API | 12-15 ms | Network overhead per request |

`Time per pair` is computed as:

\[
\text{time\_per\_pair\_ms} = \frac{\text{elapsed\_seconds}}{\text{n\_pairs}} \times 1000
\]

This metric is per gene-BP pair for each `(method, metapath, sample_size)` benchmark run.

Generate runtime comparison plots from `output/benchmark_results.csv`:

```bash
poe benchmark-dwpc-plots
```

Outputs are written to `output/benchmark_plots/`, including:
- `time_per_pair_by_metapath.png`
- `time_per_pair_overall.png`
- `speedup_api_over_direct.png`
- `benchmark_summary_by_method.csv`
- `benchmark_summary_by_metapath.csv`
- `benchmark_speedup_summary.csv`

### Year experiment pipeline

Run the staged year comparison pipeline (2016 vs 2024) with null generation,
direct DWPC computation, year-vs-control effect-size summaries, and
metapath rank-stability analyses by chaining the per-stage poe tasks (or
use the end-to-end runner / HPC submit script):

```bash
# Local end-to-end (uses run_year_pipeline.py to orchestrate stages)
python3 scripts/pipeline/run_year_pipeline.py --output-dir output/year_full_analysis

# Or assemble manually from the per-stage poe tasks
poe gen-permutation
poe gen-random
poe compute-dwpc-direct
poe year-null-variance-exp
poe year-rank-stability-exp
```

### Unified year experiment workflow

The year workflow follows the same explicit-replicate pattern as LV, but on the
`GO term x gene` bipartite graph.

Generate null replicate artifacts:

```bash
N_PERMUTATIONS=20 poe gen-permutation
N_RANDOM_SAMPLES=20 poe gen-random
poe warmup-dwpc-cache
poe compute-dwpc-direct
```

`compute-dwpc-direct` now writes both pair-level DWPC outputs and one
metapath-summary artifact per replicate under:

- `output/dwpc_direct/all_GO_positive_growth/results/`
- `output/dwpc_direct/all_GO_positive_growth/replicate_summaries/`
- `output/dwpc_direct/all_GO_positive_growth/replicate_manifest.csv`

For HPC, the intended order is:

1. generate the year control datasets
2. warm the shared DWPC cache across metapaths
3. build the dataset manifest
4. submit the dataset array
5. run year variance and rank-stability

Recommended warmup sequence:

```bash
mkdir -p output/dwpc_direct/all_GO_positive_growth
python scripts/compute_dwpc_direct.py --list-metapaths > output/dwpc_direct/all_GO_positive_growth/metapath_manifest.txt
sbatch --array=0-$(($(wc -l < output/dwpc_direct/all_GO_positive_growth/metapath_manifest.txt)-1)) hpc/year_dwpc_cache_warmup_array.sbatch
```

The dataset array assumes a warm cache and uses read-only cache access.

Run the downstream analyses, using the same explicit-replicate design as the LV
workflow:

```bash
poe year-null-variance-exp
poe year-rank-stability-exp
```

Use `20` null replicates by default for both year and LV so the workflows stay parallel. Increase that explicitly for final production runs if you need tighter stability estimates.

Outputs are written under:

- `output/year_rank_stability_exp/year_rank_stability_experiment/all_runs_long.csv`
- `output/year_rank_stability_exp/year_rank_stability_experiment/metapath_rank_table.csv`
- `output/year_rank_stability_exp/year_rank_stability_experiment/pairwise_metrics.csv`
- `output/year_rank_stability_exp/year_rank_stability_experiment/go_term_stability_summary.csv`
- `output/year_rank_stability_exp/year_rank_stability_experiment/overall_stability_summary.csv`
- `output/year_rank_stability_exp/year_rank_stability_experiment/spearman_overall_by_group.png`
- `output/year_rank_stability_exp/year_rank_stability_experiment/topk_jaccard_overall_by_group.png`

### Unified LV experiment workflow

The LV experiment now mirrors the year experiment:

1. prepare shared inputs
2. generate explicit null replicate artifacts
3. compute one metapath-summary artifact per replicate
4. run variance analysis
5. run rank-stability analysis

Prepare once:

```bash
export LV_LOADINGS_PATH=data/lv_loadings/multiplier_model_z.tsv.gz
poe lv-prepare-exp
```

For HPC, the recommended path mirrors the year workflow:

1. prepare LV metadata
2. warm the shared DWPC cache across LV metapaths
3. finalize LV precompute outputs
4. generate null replicate artifacts
5. summarize artifacts
6. run variance and rank-stability aggregates

Generate explicit LV null artifacts:

```bash
LV_N_REPLICATES=20 poe lv-gen-permutation
LV_N_REPLICATES=20 poe lv-gen-random
```

Compute one summary file per replicate artifact:

```bash
poe lv-compute-replicate-summaries
```

Run the two downstream analyses:

```bash
poe lv-null-variance-exp
poe lv-rank-stability-exp
```

All-metapath preparation:

```bash
LV_OUTPUT_DIR=output/lv_experiment_all_metapaths poe lv-prepare-exp-all-metapaths
```

Key outputs:

- `output/lv_experiment/replicate_manifest.csv`
- `output/lv_experiment/replicate_artifacts/`
- `output/lv_experiment/replicate_summaries/`
- `output/lv_experiment/lv_null_variance_experiment/feature_variance_summary.csv`
- `output/lv_experiment/lv_rank_stability_experiment/overall_stability_summary.csv`



## Run the web tool prototype

A Streamlit app (`app.py`) lets you enter a gene list and a target node and see the metapath ranking, intermediate sharing, and surviving sub-paths inline.

```bash
streamlit run app.py
```

On first launch the app preloads all G→BP DWPC matrices (~90 seconds against a warm `data/dwpc_cache/`; substantially longer on a cold cache while the matrices are computed and saved). Subsequent queries reuse the cached matrices.

# AI Assistance
This project utilized the AI assistant Claude, developed by Anthropic, during the development process. Its assistance included generating initial code snippets and improving documentation. All AI-generated content was reviewed, tested, and validated by human developers.
