# scripts

Data preparation and DWPC computation for the year analysis: do GO terms'
gene annotations added between 2016 (Hetionet v1.0) and 2024 show the same
aggregated DWPC as the 2016 annotations? Run from the repository root in the
`multi_dwpc` environment. Outputs go to `output/`.

The web app's Alpine scripts (`build_null_bundle.py`, `finalize_null_bundle.py`,
`validate_null_bundle.py`, `summarize_null_bundle_validation.py`,
`prewarm_dwpc_cache.py`) are described in [../hpc/README.md](../hpc/README.md);
`verify_null_bundle_app.py` measures the app (see
[../docs/tasks/laptop-null-bundle/](../docs/tasks/laptop-null-bundle/)).

## Pipeline

```bash
poe load-data            # Hetionet v1.0 hetmat + 2016/2024 GO annotations
poe filter-change        # percent change + IQR filtering
poe filter-jaccard       # Jaccard filtering of overlapping GO terms
poe gen-permutation      # permuted GO-label null datasets
poe gen-random           # promiscuity-matched random null datasets
poe compute-dwpc-direct  # DWPC for real and null datasets
```

`poe pipeline-production` runs the first five steps in order (logs under
`output/logs/`); `poe pipeline-null` runs the two null generators.
`python scripts/pipeline_publication.py` adds the GO hierarchy analysis.

| Script | poe task | What it does |
|---|---|---|
| `load_data.py` | `load-data` | Downloads Hetionet v1.0 and GO annotations; builds the hetmat under `data/`; keeps BP terms with 2-1000 genes and genes common to both years |
| `percent_change_and_filtering.py` | `filter-change` | Percent change in gene counts 2016 -> 2024 and IQR filtering (`all_GO_positive_growth`) |
| `jaccard_similarity_and_filtering.py` | `filter-jaccard` | Removes GO terms that overlap by Jaccard similarity |
| `permutation_null_datasets.py` | `gen-permutation` | Permuted GO-gene associations for the null |
| `random_null_datasets.py` | `gen-random` | Promiscuity-matched random gene sets for the null |
| `compute_dwpc_direct.py` | `compute-dwpc-direct` | DWPC from the hetmat sparse matrices (hetmatpy); matrices cached in `data/dwpc_cache/` |
| `pipeline_production.py` | `pipeline-production` | Runs load, filter and null steps in order |
| `pipeline_publication.py` | | Same, plus GO hierarchy analysis |

## Datasets

- `all_GO_positive_growth`: GO terms with positive growth after IQR filtering.
- `parents_GO_postive_growth`: parents of leaf terms within that set (publication pipeline).

## DWPC computation

`compute_dwpc_direct.py` computes DWPC by sparse matrix multiplication with
hetmatpy, transformed as `arcsinh(raw / raw_mean)` using
`data/metapath-dwpc-stats.tsv`. The first run computes and caches each metapath
matrix; later runs load the cache. Agreement with the Hetionet API values is
tested in `tests/test_dwpc_validation.py` (`src/dwpc_validation.py`).

The year and LV null-variance and rank-stability experiments are on the
`refactor/year-lv-unify` branch and will be documented here when it merges.
