# hpc

Slurm jobs for CU Boulder Alpine. Submit from the repository root on a login
node; logs go to `hpc/logs/`. The null-bundle scripts default to account
`amc-general`, partition `acpu`, QOS `cpu-normal` and conda from
`/curc/sw/anaconda3/2023.09`; override with `ACCOUNT`, `PARTITION`, `QOS` and
`CONDA_SH`.

## Build the null bundle (web app)

The web app scores gene sets against precomputed per-stratum null summaries for
every G->BP metapath and target, so it needs no DWPC matrices at query time.

```bash
DATA_DIR=$PWD/data bash hpc/submit_null_bundle.sh
```

- 52 array tasks (`scripts/build_null_bundle.py --index i`), one per metapath:
  load or compute the DWPC matrix into `CACHE_DIR`, then summarise every target
  column into strata. Up to ~7 min and ~17 GB each.
- A finalize job (`scripts/finalize_null_bundle.py`) runs after all tasks
  succeed and writes `$DATA_DIR/null_bundle/`: `strata.parquet` (~98M rows,
  3.1 GB), `row_sums.parquet` and `manifest.json` (settings and SHA-256 of the
  data files).

| Variable | Default |
|---|---|
| `DATA_DIR` | `<repo>/data` |
| `WORK_DIR` | `/scratch/alpine/$USER/multi-dwpc/null_bundle` |
| `CACHE_DIR` | `$WORK_DIR/dwpc_cache` (~22 GB of DWPC matrices; scratch is purged periodically) |
| `PARTS_DIR` | `$WORK_DIR/parts` |
| `BUNDLE_DIR` | `$DATA_DIR/null_bundle` |

Copy the bundle to the machine that runs the app (use the data-transfer node):

```bash
rsync -a alpine-dtn:<repo on Alpine>/data/null_bundle/ data/null_bundle/
```

The `data/` files on both machines must be identical; the app checks the
manifest hashes. The published bundle is on Zenodo
([10.5281/zenodo.22752562](https://doi.org/10.5281/zenodo.22752562)); upload a
rebuilt bundle as a new version of that record and update the record ID in the
main README.

## Validate the null bundle

```bash
DATA_DIR=$PWD/data bash hpc/submit_validate_null_bundle.sh
AFTER_JOB=<finalize job id> DATA_DIR=$PWD/data bash hpc/submit_validate_null_bundle.sh  # wait for a build
```

52 array tasks (`scripts/validate_null_bundle.py`) score 500 targets and random
and annotated gene sets both through the bundle and through the full-matrix
adapter; a summary job (`scripts/summarize_null_bundle_validation.py`) writes
tables, figures and a pass/fail verdict to `docs/tasks/laptop-null-bundle/`.

## Warm the DWPC cache

```bash
PARTITION=acpu QOS=cpu-normal bash hpc/submit_prewarm_dwpc.sh              # G -> BP
SOURCE_TYPE=G TARGET_TYPE=D PARTITION=acpu QOS=cpu-normal bash hpc/submit_prewarm_dwpc.sh
```

One task per metapath (`scripts/prewarm_dwpc_cache.py --single-metapath`)
computes DWPC matrices into `data/dwpc_cache/`. This older script's partition
and QOS defaults (`amilan`, `normal`) no longer exist on Alpine, so pass them;
it also sets `MEM` and `TIME_LIMIT`, uses `module load anaconda`, and passes no
account.
