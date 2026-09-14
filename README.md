# multi-dwpc

Multi-DWPC asks how a *set* of genes connects to a target node in Hetionet,
beyond what each gene's own connectivity predicts, and which intermediate nodes
carry that signal. It aggregates degree-weighted path counts (DWPC) over the
gene set and scores them against a degree-matched null.

Single-pair pathway search explains how one gene relates to one disease or
process, but most phenotypes involve many genes at once (for example, the ~225
genes upregulated in trisomy 21). Gene set enrichment does not show how those
genes connect across data types; a heterogeneous knowledge graph does.

This work builds on the Greene Lab and [Hetionet](https://het.io/) projects,
including [connectivity-search-analyses](https://github.com/greenelab/connectivity-search-analyses),
[hetnetpy](https://github.com/hetio/hetnetpy) and [hetmatpy](https://github.com/hetio/hetmatpy).

## Setup

```bash
git clone https://github.com/greenelab/multi-dwpc.git
cd multi-dwpc
conda env create -f env/environment.yml
conda activate multi_dwpc
pip install -e ".[dev]"
```

## Run the web app

The Streamlit app takes a gene list and a target Biological Process and shows
the metapath ranking, intermediate sharing and the top paths. It runs on a
laptop (~1.2 GB of memory).

**1. Get the data.** The app reads everything from `data/` in the repository
root:

```
data/
├── metagraph.json            poe load-data
├── nodes/                    poe load-data (node tables)
├── edges/                    poe load-data (edge matrices)
├── metapath-dwpc-stats.tsv   downloaded below
└── null_bundle/              downloaded below from Zenodo (3.1 GB)
    ├── manifest.json
    ├── row_sums.parquet
    └── strata.parquet
```

```bash
poe load-data    # Hetionet v1.0 -> data/metagraph.json, data/nodes/, data/edges/
curl -L -o data/metapath-dwpc-stats.tsv \
  https://raw.githubusercontent.com/greenelab/hetmech/34e95b9f72f47cdeba3d51622bee31f79e9a4cb8/explore/bulk-pipeline/archives/metapath-dwpc-stats.tsv

# Null bundle: https://doi.org/10.5281/zenodo.22752562
mkdir -p data/null_bundle
for f in manifest.json row_sums.parquet strata.parquet; do
  curl -L -o "data/null_bundle/$f" "https://zenodo.org/records/22752562/files/$f?download=1"
done
```

Compare the files' MD5 checksums with the ones listed on the Zenodo record.
The null bundle holds precomputed null summaries for every Gene -> Biological
Process metapath and target; the app refuses it if `data/` was built from
different Hetionet files (its manifest records their hashes). To build the
bundle yourself, see [hpc/README.md](hpc/README.md).

**2. Start the app.**

```bash
streamlit run app.py
```

Open http://localhost:8501, click **Run query** (a worked example is
prefilled), or use **Random genes** / **Random BP**, or paste your own genes.
Pick a metapath in the sidebar for intermediate sharing and subgraphs. A query
takes ~0.3-1.7 s.

**3. Test.**

```bash
python -m pytest tests -q
```

## More documentation

- [scripts/README.md](scripts/README.md): data preparation, null datasets and
  direct DWPC computation for the year analysis.
- [hpc/README.md](hpc/README.md): Alpine jobs that build and validate the null
  bundle and warm the DWPC cache.
- [docs/tasks/laptop-null-bundle/](docs/tasks/laptop-null-bundle/): how the app
  computes its null without DWPC matrices, with validation and measurements.
- [docs/web_tool_requirements.md](docs/web_tool_requirements.md): requirements
  for the production web tool.

## AI assistance

This project used the AI assistant Claude, developed by Anthropic, to generate
initial code and improve documentation. All AI-generated content was reviewed,
tested and validated by human developers.
