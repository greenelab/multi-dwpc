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
```

## Prepare the Data

```bash
papermill notebooks/1.1_data_loading.ipynb notebooks/1.1_data_loading.ipynb && \
papermill notebooks/1.2_percent_change_and_filtering.ipynb notebooks/1.2_percent_change_and_filtering.ipynb && \
papermill notebooks/1.3_jaccard_similarity_and_filtering.ipynb notebooks/1.3_jaccard_similarity_and_filtering.ipynb
```

Pipeline notebooks:

1. **1.1_data_loading.ipynb** - Loads Hetionet v1.0 (2016) and GO annotations (2024), filters to common genes and GO terms
2. **1.2_percent_change_and_filtering.ipynb** - Filters GO ontology terms by positive change between 2024 and 2016, GO terms in the IQR of positive change, and GO terms that are the immediate parents of leaf terms
3. **1.3_jaccard_similarity_and_filtering.ipynb** - Calculates Jaccard similarity between GO terms, filters overlapping terms (>0.1 threshold), and generates Neo4j ID mappings from the Hetionet API

### Papermill Parameters

Some notebooks support papermill parameters for customizing execution:

**1.3_jaccard_similarity_and_filtering.ipynb**:
- `FORCE_RECOMPUTE` (default: `False`) - Force recomputation of Jaccard similarity matrices and heatmaps, ignoring cached results

Example usage:
```bash
papermill notebooks/1.3_jaccard_similarity_and_filtering.ipynb \
    notebooks/1.3_jaccard_similarity_and_filtering.ipynb \
    -p FORCE_RECOMPUTE True
```

## Run the DWPC computation 

Notebook 2 computes Degree-Weighted Path Counts (DWPC) via the Hetionet API. This requires running the connectivity-search-backend Docker stack.

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

**3. Run notebook 2:**

```bash
cd ..
papermill notebooks/2_compute_hetio_stat_via_docker.ipynb \
    notebooks/2_compute_hetio_stat_via_docker.ipynb
```

**4. Stop the Docker stack when finished:**

```bash
cd connectivity-search-backend
docker compose down
```

# AI Assistance
This project utilized the AI assistant Claude, developed by Anthropic, during the development process. Its assistance included generating initial code snippets and improving documentation. All AI-generated content was reviewed, tested, and validated by human developers.